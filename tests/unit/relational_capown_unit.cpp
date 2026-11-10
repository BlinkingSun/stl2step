// Synthetic cases for the lateral-facet ownership law of committed cylinders
// (D-gid-8.2, releaseCapTriangles): a cap fan a cylinder region owns (every
// vertex on the rim arc, facet normal parallel to the axis) moves whole to the
// one adjacent tau-certified plane whose vertex certificate covers it; an
// uncertified plane, a vertex off the plane, no plane, two qualifying planes,
// a region made only of caps, or a lateral partner that misses its own tau
// certificate leave it where it is. The cylinder keeps its fit; lateral
// triangles never move. Generic quarter-round geometry only; no
// training mesh.
// SPDX-License-Identifier: MIT

#include "refit.hpp"
#include "refit_internal.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;
    double q = 1e-6;

    int add(const gp_XYZ& p) {
        pts.push_back(p);
        return (int)pts.size() - 1;
    }
    int addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
        return (int)(triFlat.size() / 3) - 1;
    }
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2], triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]), hi = std::max(v[s], v[s + 1]);
                auto it = ids.find({lo, hi});
                int id;
                if (it == ids.end()) {
                    id = (int)edges.size();
                    ids.emplace(std::make_pair(lo, hi), id);
                    edges.emplace_back(lo, hi);
                } else {
                    id = it->second;
                }
                triEdges[t][(std::size_t)s] = id;
                if (v[s] < v[s + 1]) triDirs[t] = (uint8_t)(triDirs[t] | (uint8_t)(1u << s));
            }
        }
        compVtx.resize(pts.size());
        for (std::size_t i = 0; i < pts.size(); ++i) compVtx[i] = (int)i;
    }
    MeshView view() const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int (*)[3]>(triFlat.data());
        mv.compTris = compTris.data();
        mv.compVtx = compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = triEdges.data();
        mv.triDirs = triDirs.data();
        mv.nTri = triFlat.size() / 3;
        mv.nVtx = pts.size();
        mv.nEdge = edges.size();
        mv.diag = 100.0;
        mv.quantFloor = q;
        return mv;
    }
};

// A quarter-round wall (axis +Z through the origin, radius kR, z in [0, kH]) as
// kSeg sectors of two triangles; its top rim a[0..kSeg]. The top cap's circular
// segment between the arc and the chord a[0]-a[kSeg] is a fan from a[0]
// (every vertex on the arc); the plane triangle (a[0], a[kSeg], centre) lies
// across that chord.
constexpr int kSeg = 6;
constexpr double kR = 2.0, kH = 3.0;

struct Part {
    Synth mesh;
    std::vector<int> wall, fan;
    int planeTri = -1;
    int floorTri = -1;  // across the wall's first bottom edge (a lateral partner)
    SegmentWork work;
};

static Part quarter(double liftFanVertex = 0.0, bool withPlane = true, bool withFloor = false) {
    Part p;
    std::vector<int> a, b;
    for (int k = 0; k <= kSeg; ++k) {
        const double th = 0.5 * kPi * (double)k / (double)kSeg;
        const double lift = (k == 2) ? liftFanVertex : 0.0;
        a.push_back(p.mesh.add(gp_XYZ(kR * std::cos(th), kR * std::sin(th), kH + lift)));
        b.push_back(p.mesh.add(gp_XYZ(kR * std::cos(th), kR * std::sin(th), 0.0)));
    }
    for (int k = 0; k < kSeg; ++k) {
        p.wall.push_back(p.mesh.addTri(b[k], b[k + 1], a[k + 1]));
        p.wall.push_back(p.mesh.addTri(b[k], a[k + 1], a[k]));
    }
    for (int k = 1; k < kSeg; ++k) p.fan.push_back(p.mesh.addTri(a[0], a[k], a[k + 1]));
    const int c = p.mesh.add(gp_XYZ(0.0, 0.0, kH));
    if (withPlane) p.planeTri = p.mesh.addTri(a[0], c, a[kSeg]);
    if (withFloor) p.floorTri = p.mesh.addTri(b[1], b[0], p.mesh.add(gp_XYZ(0.0, 0.0, 0.0)));
    p.mesh.finish();
    return p;
}

static Region cylinderOf(const std::vector<int>& tris) {
    Region r;
    r.type = SurfType::Cylinder;
    r.origin = Origin::FilletStrip;
    r.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.radius = kR;
    r.tris = tris;
    std::sort(r.tris.begin(), r.tris.end());
    r.maxVertexDev = 0.0;
    return r;
}

static Region planeOf(const std::vector<int>& tris, double maxDev, double z = kH) {
    Region r;
    r.type = SurfType::Plane;
    r.origin = Origin::PlaneGrow;
    r.ax = gp_Ax3(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.tris = tris;
    std::sort(r.tris.begin(), r.tris.end());
    r.maxVertexDev = maxDev;
    return r;
}

static DerivedTols tols() {
    DerivedTols t;
    t.thetaSharp = 30.0 * kPi / 180.0;  // SegmentParams::thetaSharpDeg default
    return t;
}

static const Region* byType(const SegmentWork& w, SurfType st) {
    for (const Region& r : w.accepted)
        if (r.type == st) return &r;
    return nullptr;
}

static std::vector<int> joined(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> v = a;
    v.insert(v.end(), b.begin(), b.end());
    std::sort(v.begin(), v.end());
    return v;
}

int main() {
    const DerivedTols tol = tols();
    {  // the fan moves whole to the certified cap plane; the cylinder keeps its fit
        Part p = quarter();
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        expect(releaseCapTriangles(mv, tol, p.work), "move: stage returns true");
        const Region* cy = byType(p.work, SurfType::Cylinder);
        const Region* pl = byType(p.work, SurfType::Plane);
        std::vector<int> wall = p.wall;
        std::sort(wall.begin(), wall.end());
        expect(cy && cy->tris == wall, "move: cylinder keeps exactly its lateral triangles");
        expect(pl && pl->tris == joined({p.planeTri}, p.fan), "move: plane owns the whole fan");
        expect(cy && cy->radius == kR && cy->ax.Location().Distance(gp_Pnt(0, 0, 0)) == 0.0,
               "move: cylinder fit untouched");
        expect(pl && pl->maxVertexDev <= tauOf(mv), "move: plane stays tau-certified");
    }
    {  // an uncertified plane takes nothing
        Part p = quarter();
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 100.0 * tauOf(mv)));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        expect(cy && cy->tris == joined(p.wall, p.fan), "uncertified plane: fan stays");
    }
    {  // one fan vertex 3 tau off the plane (still on the cylinder): the set is not split
        Part p = quarter(3.0 * 2e-6);
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        const Region* pl = byType(p.work, SurfType::Plane);
        expect(cy && cy->tris == joined(p.wall, p.fan), "vertex-off: whole fan stays");
        expect(pl && pl->tris.size() == 1, "vertex-off: plane takes no subset");
    }
    {  // no committed plane across the chord: nothing moves
        Part p = quarter(0.0, true);
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        expect(cy && cy->tris == joined(p.wall, p.fan), "no plane: fan stays");
    }
    {  // two certified planes qualify for one set: refused
        Part p = quarter();
        // a second plane region across the fan's first spoke: split the fan so
        // its triangle at a[0]-a[1]-a[2] belongs to a coplanar plane region
        const MeshView mv = p.mesh.view();
        std::vector<int> fanRest(p.fan.begin() + 1, p.fan.end());
        p.work.accepted.push_back(cylinderOf(joined(p.wall, fanRest)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        p.work.accepted.push_back(planeOf({p.fan.front()}, 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        expect(cy && cy->tris == joined(p.wall, fanRest), "two planes: ambiguous set stays");
    }
    {  // a region made only of caps is not a wall to release from
        Part p = quarter();
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(p.fan));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        std::vector<int> fan = p.fan;
        std::sort(fan.begin(), fan.end());
        expect(cy && cy->tris == fan, "all-cap: region unchanged");
    }
    {  // a wall without caps is untouched
        Part p = quarter();
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(p.wall));
        p.work.accepted.push_back(planeOf(joined({p.planeTri}, p.fan), 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        std::vector<int> wall = p.wall;
        std::sort(wall.begin(), wall.end());
        expect(cy && cy->tris == wall, "lateral only: nothing moves");
    }
    {  // a lateral partner that misses tau: the cylinder's caps stay (no partial move)
        Part p = quarter(0.0, true, true);
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        p.work.accepted.push_back(planeOf({p.floorTri}, 100.0 * tauOf(mv), 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        expect(cy && cy->tris == joined(p.wall, p.fan), "uncertified partner: fan stays");
    }
    {  // the same partner certified: the fan moves
        Part p = quarter(0.0, true, true);
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        p.work.accepted.push_back(planeOf({p.floorTri}, 0.0, 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        std::vector<int> wall = p.wall;
        std::sort(wall.begin(), wall.end());
        expect(cy && cy->tris == wall, "certified partner: fan moves");
    }
    {  // q <= 0 admits nothing
        Part p = quarter();
        p.mesh.q = 0.0;
        const MeshView mv = p.mesh.view();
        p.work.accepted.push_back(cylinderOf(joined(p.wall, p.fan)));
        p.work.accepted.push_back(planeOf({p.planeTri}, 0.0));
        releaseCapTriangles(mv, tol, p.work);
        const Region* cy = byType(p.work, SurfType::Cylinder);
        expect(cy && cy->tris == joined(p.wall, p.fan), "q=0: nothing moves");
    }
    std::fprintf(stderr, "%s relational_capown_unit (%d failed)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
