// D-fc-6.1: one-cycle closed360 demotion, and dissolving one region that
// fails classification while a sibling region stays. Synthetic census only.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <array>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

#if !defined(_WIN32)
#include <pthread.h>
#include <sys/resource.h>
#endif

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;

    int add(double x, double y, double z) {
        pts.push_back(gp_XYZ(x, y, z));
        return (int)pts.size() - 1;
    }
    int addTri(int a, int b, int c) {
        const int t = (int)triFlat.size() / 3;
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
        return t;
    }
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2],
                              triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]);
                const int hi = std::max(v[s], v[s + 1]);
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
        mv.diag = 20.0;
        mv.quantFloor = 1e-6;
        return mv;
    }
};

static gp_Ax3 zAx(double z) {
    return gp_Ax3(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
}

static Region cylRegion(std::vector<int> tris, bool closed) {
    Region r;
    r.id = 0;
    r.type = SurfType::Cylinder;
    r.origin = Origin::CylGrow;
    r.ax = zAx(0);
    r.radius = 1.0;
    r.uMin = 0.0;
    r.uMax = 6.283185307179586;
    r.vMin = 0.0;
    r.vMax = 1.0;
    r.closed360 = closed;
    r.outwardNormal = true;
    r.tris = std::move(tris);
    return r;
}

static Region planeRegion(std::vector<int> tris) {
    Region r;
    r.id = 0;
    r.type = SurfType::Plane;
    r.origin = Origin::PlaneGrow;
    r.ax = gp_Ax3(gp_Pnt(10, 0, 5), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.closed360 = false;
    r.outwardNormal = true;
    r.tris = std::move(tris);
    return r;
}

// Three quads of a unit cylinder, gap between the last and the first station.
// The boundary is one cycle. Returns the six triangle ids.
static std::vector<int> addOpenWall(Synth& s) {
    const int b0 = s.add(1, 0, 0), b1 = s.add(0, 1, 0), b2 = s.add(-1, 0, 0), b3 = s.add(0, -1, 0);
    const int t0 = s.add(1, 0, 1), t1 = s.add(0, 1, 1), t2 = s.add(-1, 0, 1), t3 = s.add(0, -1, 1);
    std::vector<int> tris;
    auto quad = [&](int a, int b, int c, int d) {
        tris.push_back(s.addTri(a, b, c));
        tris.push_back(s.addTri(a, c, d));
    };
    quad(b0, b1, t1, t0);
    quad(b1, b2, t2, t1);
    quad(b2, b3, t3, t2);
    return tris;
}

static std::vector<int> addClosedWall(Synth& s) {
    const int b0 = s.add(1, 0, 0), b1 = s.add(0, 1, 0), b2 = s.add(-1, 0, 0), b3 = s.add(0, -1, 0);
    const int t0 = s.add(1, 0, 1), t1 = s.add(0, 1, 1), t2 = s.add(-1, 0, 1), t3 = s.add(0, -1, 1);
    std::vector<int> tris;
    auto quad = [&](int a, int b, int c, int d) {
        tris.push_back(s.addTri(a, b, c));
        tris.push_back(s.addTri(a, c, d));
    };
    quad(b0, b1, t1, t0);
    quad(b1, b2, t2, t1);
    quad(b2, b3, t3, t2);
    quad(b3, b0, t0, t3);
    return tris;
}

static std::vector<int> addPlane(Synth& s) {
    const int p0 = s.add(10, 0, 5), p1 = s.add(12, 0, 5), p2 = s.add(12, 2, 5), p3 = s.add(10, 2, 5);
    return {s.addTri(p0, p1, p2), s.addTri(p0, p2, p3)};
}

// Flat annulus: outer square and inner square, both at z = 0, so two boundary
// cycles share one axial mean. A closed360 cylinder cannot pick two caps.
static std::vector<int> addAnnulus(Synth& s) {
    const int o0 = s.add(0, 0, 0), o1 = s.add(4, 0, 0), o2 = s.add(4, 4, 0), o3 = s.add(0, 4, 0);
    const int i0 = s.add(1, 1, 0), i1 = s.add(3, 1, 0), i2 = s.add(3, 3, 0), i3 = s.add(1, 3, 0);
    std::vector<int> tris;
    auto quad = [&](int a, int b, int c, int d) {
        tris.push_back(s.addTri(a, b, c));
        tris.push_back(s.addTri(a, c, d));
    };
    quad(o0, o1, i1, i0);
    quad(o1, o2, i2, i1);
    quad(o2, o3, i3, i2);
    quad(o3, o0, i0, i3);
    return tris;
}

static const Region* findType(const RegionSet& rs, SurfType ty) {
    for (const Region& r : rs.regions)
        if (r.type == ty) return &r;
    return nullptr;
}

static int countRole(const Region& r, LoopRole role) {
    int n = 0;
    for (const Loop& lp : r.loops)
        if (lp.role == role) ++n;
    return n;
}

static bool run(Synth& s, std::vector<Region> accepted, RegionSet& out) {
    s.finish();
    SegmentWork work;
    work.accepted = std::move(accepted);
    SegmentParams p{};
    DerivedTols tol{};
    tol.epsPlane = 0.02;
    tol.epsMesh = 1e-4;
    return buildTopologyD(s.view(), p, tol, work, out);
}

static void caseDemoteOneCycle() {
    Synth s;
    const std::vector<int> wall = addOpenWall(s);
    const std::vector<int> plane = addPlane(s);
    RegionSet out;
    const bool ok = run(s, {cylRegion(wall, true), planeRegion(plane)}, out);
    expect(ok, "demote: stage returns true");
    const Region* cy = findType(out, SurfType::Cylinder);
    const Region* pl = findType(out, SurfType::Plane);
    expect(cy != nullptr && pl != nullptr, "demote: cylinder and sibling plane both kept");
    if (!cy || !pl) return;
    expect(!cy->closed360, "demote: one cycle clears closed360");
    expect(cy->loops.size() == 1 && countRole(*cy, LoopRole::Outer) == 1 &&
               countRole(*cy, LoopRole::CapLow) == 0 && countRole(*cy, LoopRole::CapHigh) == 0,
           "demote: partial arm assigns one Outer");
    expect(pl->loops.size() == 1 && countRole(*pl, LoopRole::Outer) == 1,
           "demote: sibling plane keeps one Outer");
    bool wallClaimed = true;
    for (int t : wall)
        if (t < 0 || (size_t)t >= out.triRegion.size() || out.triRegion[(size_t)t] < 0)
            wallClaimed = false;
    expect(wallClaimed, "demote: wall triangles stay accepted");
}

static void caseTwoCapsStayClosed() {
    Synth s;
    const std::vector<int> wall = addClosedWall(s);
    RegionSet out;
    const bool ok = run(s, {cylRegion(wall, true)}, out);
    expect(ok, "two-cap: stage returns true");
    const Region* cy = findType(out, SurfType::Cylinder);
    expect(cy != nullptr, "two-cap: region kept");
    if (!cy) return;
    expect(cy->closed360, "two-cap: closed360 stays set");
    expect(countRole(*cy, LoopRole::CapLow) == 1 && countRole(*cy, LoopRole::CapHigh) == 1 &&
               countRole(*cy, LoopRole::Outer) == 0,
           "two-cap: CapLow and CapHigh, no Outer");
}

static void caseTorusMouthStaysClosed() {
    Synth s;
    const std::vector<int> wall = addOpenWall(s);
    Region mouth = cylRegion(wall, true);
    mouth.type = SurfType::Torus;
    mouth.origin = Origin::TorusBlend;
    mouth.radius = 3.0;
    mouth.radius2 = 1.0;
    mouth.vMin = -0.5;
    mouth.vMax = 0.5;
    RegionSet out;
    const bool ok = run(s, {mouth}, out);
    expect(ok, "torus-mouth: stage returns true");
    const Region* ty = findType(out, SurfType::Torus);
    expect(ty != nullptr && ty->closed360, "torus-mouth: one cycle stays closed360");
    if (!ty) return;
    expect(ty->loops.size() == 1 && countRole(*ty, LoopRole::Outer) == 1,
           "torus-mouth: existing one-loop arm assigns Outer");
}

static void caseDissolveOneRegion() {
    Synth s;
    const std::vector<int> ring = addAnnulus(s);
    const std::vector<int> plane = addPlane(s);
    RegionSet out;
    const bool ok = run(s, {cylRegion(ring, true), planeRegion(plane)}, out);
    expect(ok, "dissolve: stage returns true");
    expect(findType(out, SurfType::Cylinder) == nullptr, "dissolve: failed region leaves the set");
    const Region* pl = findType(out, SurfType::Plane);
    expect(pl != nullptr && pl->loops.size() == 1 && countRole(*pl, LoopRole::Outer) == 1,
           "dissolve: sibling plane still classifies");
    bool islanded = !ring.empty();
    for (int t : ring) {
        if (t < 0 || (size_t)t >= out.triRegion.size() || (size_t)t >= out.triIsland.size() ||
            out.triRegion[(size_t)t] >= 0 || out.triIsland[(size_t)t] < 0)
            islanded = false;
    }
    expect(islanded, "dissolve: failed triangles are a facet island");
    bool planeKept = !plane.empty();
    for (int t : plane) {
        if (t < 0 || (size_t)t >= out.triRegion.size() || out.triRegion[(size_t)t] < 0)
            planeKept = false;
    }
    expect(planeKept, "dissolve: sibling triangles stay accepted");
}

static void caseOverlapStillFails() {
    Synth s;
    const std::vector<int> plane = addPlane(s);
    RegionSet out;
    const bool ok = run(s, {planeRegion(plane), planeRegion(plane)}, out);
    expect(!ok, "overlap: structural failure still aborts the stage");
}

// K flat annuli, each a closed360 cylinder whose two boundary cycles share one
// axial mean (iLow == iHigh), plus one sibling plane. Every annulus dissolves.
// The retry used to recurse once per failure and keep every RegionSet, which
// overflowed a 512 KB worker stack near this K. One live attempt must finish
// on that stack, and the process high-water must stay far below the recursive
// measurement (hundreds of MB at this K on this shape).
struct StressOut {
    bool ok = false;
    bool planeKept = false;
    int regions = -1;
    int islandTris = -1;
    int nIslands = -1;
    int thrown = 0;
};

static void runManyFailing(int k, StressOut* o) {
    try {
        Synth s;
        std::vector<Region> acc;
        std::vector<int> ring;
        acc.reserve((std::size_t)k + 1);
        ring.reserve((std::size_t)k * 8);
        for (int i = 0; i < k; ++i) {
            const std::vector<int> tris = addAnnulus(s);
            ring.insert(ring.end(), tris.begin(), tris.end());
            acc.push_back(cylRegion(tris, true));
        }
        const std::vector<int> plane = addPlane(s);
        acc.push_back(planeRegion(plane));
        RegionSet out;
        o->ok = run(s, std::move(acc), out);
        o->regions = (int)out.regions.size();
        const Region* pl = findType(out, SurfType::Plane);
        o->planeKept = pl != nullptr && pl->loops.size() == 1 &&
                       countRole(*pl, LoopRole::Outer) == 1 &&
                       findType(out, SurfType::Cylinder) == nullptr;
        int nIsland = 0;
        for (int t : ring) {
            if (t >= 0 && (std::size_t)t < out.triRegion.size() &&
                (std::size_t)t < out.triIsland.size() && out.triRegion[(std::size_t)t] < 0 &&
                out.triIsland[(std::size_t)t] >= 0)
                ++nIsland;
        }
        o->islandTris = nIsland;
        o->nIslands = out.nIslands;
        bool planeClaimed = !plane.empty();
        for (int t : plane) {
            if (t < 0 || (std::size_t)t >= out.triRegion.size() || out.triRegion[(std::size_t)t] < 0)
                planeClaimed = false;
        }
        if (!planeClaimed) o->planeKept = false;
    } catch (...) {
        o->thrown = 1;
        o->ok = false;
    }
}

#if !defined(_WIN32)
static void* stressTrampoline(void* arg) {
    auto* pack = static_cast<std::pair<int, StressOut*>*>(arg);
    runManyFailing(pack->first, pack->second);
    return nullptr;
}

static long rssHighWaterBytes() {
    struct rusage ru {};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return -1;
#if defined(__APPLE__)
    return (long)ru.ru_maxrss;
#else
    return (long)ru.ru_maxrss * 1024L;
#endif
}
#endif

static void caseManyFailingRegionsBounded() {
    constexpr int K = 300;
    StressOut out;
#if !defined(_WIN32)
    // Same stack the audit overflowed at this K (macOS std::thread default).
    constexpr std::size_t kStackBytes = 512 * 1024;
    const long rss0 = rssHighWaterBytes();
    pthread_attr_t attr;
    if (pthread_attr_init(&attr) != 0) {
        expect(false, "stress: pthread_attr_init");
        return;
    }
    const int stackRc = pthread_attr_setstacksize(&attr, kStackBytes);
    pthread_t th;
    std::pair<int, StressOut*> pack{K, &out};
    const int createRc = (stackRc == 0) ? pthread_create(&th, &attr, stressTrampoline, &pack) : stackRc;
    pthread_attr_destroy(&attr);
    if (createRc != 0) {
        expect(false, "stress: small-stack thread did not start");
        return;
    }
    pthread_join(th, nullptr);
    const long rss1 = rssHighWaterBytes();
    const long delta = (rss0 >= 0 && rss1 >= rss0) ? (rss1 - rss0) : -1;
    // One attempt's working set on this mesh. The recursive form measured
    // 658 MB at K=400 on the same annulus shape; 128 MiB is the ceiling
    // that still fails that blow-up and passes a single live attempt.
    constexpr long kRssCeiling = 128L * 1024L * 1024L;
    std::fprintf(stderr,
                 "INFO stress K=%d stack=%zu rss0=%ld rss1=%ld rssDelta=%ld ceiling=%ld islands=%d islandTris=%d\n",
                 K, kStackBytes, rss0, rss1, delta, kRssCeiling, out.nIslands, out.islandTris);
    expect(delta >= 0 && delta < kRssCeiling, "stress: resident set stays within one attempt");
#else
    runManyFailing(K, &out);
#endif
    expect(out.thrown == 0 && out.ok, "stress: stage returns true");
    expect(out.planeKept && out.regions == 1, "stress: sibling plane kept, failing regions leave");
    expect(out.islandTris == 8 * K, "stress: every failing triangle is a facet island");
}

int main() {
    caseDemoteOneCycle();
    caseTwoCapsStayClosed();
    caseTorusMouthStaysClosed();
    caseDissolveOneRegion();
    caseOverlapStillFails();
    caseManyFailingRegionsBounded();
    std::fprintf(stderr, "%s (%d failed)\n", gFail ? "FAILED" : "OK", gFail);
    return gFail ? 1 : 0;
}
