// D-fc-6.1(d). hostSagitta and the F1 island predicate. A sliver whose vertices
// already belong to one closed360 cylinder is absorbed; a triangle whose sagitta
// exceeds that host's own is refused; two hosts are refused. Synthetic census.
// SPDX-License-Identifier: MIT

#include "refit_island_absorb.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846;

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

struct Wall {
    Synth s;
    int b0, b1, b2, b3, t0, t1, t2, t3;
    std::vector<int> tris;
};

static Wall makeWall() {
    Wall w;
    w.b0 = w.s.add(1, 0, 0);
    w.b1 = w.s.add(0, 1, 0);
    w.b2 = w.s.add(-1, 0, 0);
    w.b3 = w.s.add(0, -1, 0);
    w.t0 = w.s.add(1, 0, 1);
    w.t1 = w.s.add(0, 1, 1);
    w.t2 = w.s.add(-1, 0, 1);
    w.t3 = w.s.add(0, -1, 1);
    auto quad = [&](int a, int b, int c, int d) {
        w.tris.push_back(w.s.addTri(a, b, c));
        w.tris.push_back(w.s.addTri(a, c, d));
    };
    quad(w.b0, w.b1, w.t1, w.t0);
    quad(w.b1, w.b2, w.t2, w.t1);
    quad(w.b2, w.b3, w.t3, w.t2);
    quad(w.b3, w.b0, w.t0, w.t3);
    return w;
}

static Region cylOf(std::vector<int> tris) {
    Region r;
    r.id = 0;
    r.type = SurfType::Cylinder;
    r.origin = Origin::CylGrow;
    r.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.radius = 1.0;
    r.uMin = 0.0;
    r.uMax = 2.0 * kPi;
    r.vMin = 0.0;
    r.vMax = 1.0;
    r.closed360 = true;
    r.outwardNormal = true;
    r.tris = std::move(tris);
    return r;
}

int main() {
    const double lateral = std::sin(30.0 * kPi / 180.0);
    Wall w = makeWall();
    const int sliver = w.s.addTri(w.b0, w.t0, w.b1);
    const int flat = w.s.addTri(w.b0, w.b2, w.t0);
    const int strangerPt = w.s.add(2, 0, 0);
    const int stranger = w.s.addTri(w.b0, w.t0, strangerPt);
    w.s.finish();
    const MeshView mv = w.s.view();
    const Region host = cylOf(w.tris);

    const double wallSag = 1.0 - std::sqrt(5.0 / 9.0);
    const double got = hostSagitta(mv, host);
    expect(std::fabs(got - wallSag) < 1e-12, "hostSagitta: max chord sagitta of the host triangles");

    double sag = -1.0;
    expect(f1WouldAbsorb(mv, host, {sliver}, lateral, &sag),
           "sliver: vertices are the host's and the sagitta is within hostSagitta");
    expect(sag <= got && sag > 0.0, "sliver: reported sag does not exceed hostSagitta");

    int one = -1;
    double oneSag = 0, oneHost = 0;
    expect(uniqueF1Host(mv, {host}, {sliver}, lateral, one, oneSag, oneHost) && one == 0,
           "sliver: exactly one host");

    expect(!f1WouldAbsorb(mv, host, {flat}, lateral, nullptr),
           "flat: sagitta above hostSagitta is refused");
    expect(!f1WouldAbsorb(mv, host, {stranger}, lateral, nullptr),
           "stranger: a vertex the host does not own is refused");

    Region other = host;
    other.id = 1;
    int two = -1;
    double twoSag = 0, twoHost = 0;
    expect(!uniqueF1Host(mv, {host, other}, {sliver}, lateral, two, twoSag, twoHost),
           "two hosts: refused");
    expect(f1WouldAbsorb(mv, host, {sliver}, lateral, nullptr) &&
               f1WouldAbsorb(mv, other, {sliver}, lateral, nullptr),
           "two hosts: each host alone would take the sliver");

    // A horizontal triangle of host vertices, with the host sagitta raised by an
    // inset triangle so the sag gate would pass. The lateral gate still refuses.
    {
        Wall h = makeWall();
        const int axis = h.s.add(0, 0, 0);
        std::vector<int> tris = h.tris;
        tris.push_back(h.s.addTri(h.b0, h.b2, axis));
        const int cap = h.s.addTri(h.b0, h.b1, h.b2);
        h.s.finish();
        const Region raised = cylOf(tris);
        expect(hostSagitta(h.s.view(), raised) + 1e-12 >= (2.0 / 3.0),
               "lateral: the raised host sagitta covers the flat triangle");
        expect(!f1WouldAbsorb(h.s.view(), raised, {cap}, lateral, nullptr),
               "lateral: |n·axis| above sin(thetaSharp) is refused");
        expect(f1WouldAbsorb(h.s.view(), raised, {cap}, 1.0, nullptr),
               "lateral: the same flat passes when the lateral gate is open");
    }

    if (gFail) std::fprintf(stderr, "%d FAIL\n", gFail);
    else std::fprintf(stderr, "ALL PASS\n");
    return gFail ? 1 : 0;
}
