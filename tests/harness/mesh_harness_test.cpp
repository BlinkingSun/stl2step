// ctest fixture checks: degenerate-triangle drop, two-body split, cube truth.
//
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;

static int gFails = 0;

static void check(bool cond, const char* msg) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", msg);
        gFails++;
    }
}

static void checkLocalIds(const stl2step::harness::HarnessComponent& c,
                          const char* label) {
    const int nVtx = (int)c.compVtx.size();
    const int nEdge = (int)c.compEdges.size();
    const int nTri = (int)c.compTris.size();
    check(nVtx > 0, "component has vertices");
    check(nEdge > 0, "component has edges");
    check(nTri > 0, "component has triangles");
    check((int)c.triEdges.size() == nTri, "triEdges size matches tris");
    check((int)c.triDirs.size() == nTri, "triDirs size matches tris");

    int minV = nVtx, maxV = -1;
    for (const auto& e : c.compEdges) {
        check(e.first >= 0 && e.first < nVtx, "edge vLo in local range");
        check(e.second >= 0 && e.second < nVtx, "edge vHi in local range");
        if (e.first < minV) minV = e.first;
        if (e.second < minV) minV = e.second;
        if (e.first > maxV) maxV = e.first;
        if (e.second > maxV) maxV = e.second;
    }
    check(minV == 0, "local vertex ids start at 0");

    for (const auto& te : c.triEdges) {
        for (int s = 0; s < 3; s++) {
            check(te[s] >= 0 && te[s] < nEdge, "triEdges id in local range");
        }
    }
    for (int i = 1; i < nTri; i++) {
        check(c.compTris[i - 1] < c.compTris[i],
              "compTris stored in ascending global triangle id");
    }
    (void)label;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: %s <degenerate.stl> <two_cubes.stl> [cube.stl]\n",
                     argv[0]);
        return 2;
    }
    const char* degPath = argv[1];
    const char* twoPath = argv[2];
    const char* cubePath = argc > 3 ? argv[3] : nullptr;

    {
        HarnessMesh m;
        std::string err;
        check(loadMesh(degPath, 1.0, 0.0, 0.0, m, err), "load degenerate.stl");
        if (!err.empty() && m.tris.empty()) {
            std::fprintf(stderr, "  err: %s\n", err.c_str());
        }
        check(m.comps.size() == 1, "degenerate: 1 component");
        check(m.tris.size() == 12, "degenerate: extra triangle dropped -> 12 tris");
        check(m.pts.size() == 8, "degenerate: 8 welded vertices");
        if (!m.comps.empty()) {
            check(m.comps[0].compTris.size() == 12, "degenerate: 12 local tris");
            check(m.comps[0].compVtx.size() == 8, "degenerate: 8 local vtx");
            check(m.comps[0].compEdges.size() == 18, "degenerate: 18 local edges");
            check(m.comps[0].clean(), "degenerate: clean");
            checkLocalIds(m.comps[0], "degenerate");
        }
    }

    {
        HarnessMesh m;
        std::string err;
        check(loadMesh(twoPath, 1.0, 0.0, 0.0, m, err), "load two_cubes.stl");
        if (!err.empty() && m.tris.empty()) {
            std::fprintf(stderr, "  err: %s\n", err.c_str());
        }
        check(m.comps.size() == 2, "two_cubes: 2 components");
        check(m.tris.size() == 24, "two_cubes: 24 global tris");
        check(m.pts.size() == 16, "two_cubes: 16 welded vertices");
        if (m.comps.size() == 2) {
            check(m.comps[0].root < m.comps[1].root, "two_cubes: comps ordered by ascending root");
            for (int i = 0; i < 2; i++) {
                check(m.comps[i].compTris.size() == 12, "two_cubes: 12 tris per body");
                check(m.comps[i].compVtx.size() == 8, "two_cubes: 8 vtx per body");
                check(m.comps[i].compEdges.size() == 18, "two_cubes: 18 edges per body");
                check(m.comps[i].clean(), "two_cubes: each body clean");
                checkLocalIds(m.comps[i], "two_cubes");
            }
            // Independent local id spaces: both map local 0, and the global
            // points they name are distinct.
            check(m.comps[0].compVtx[0] != m.comps[1].compVtx[0],
                  "two_cubes: independent global points at local vtx 0");
        }
    }

    if (cubePath) {
        HarnessMesh m;
        std::string err;
        check(loadMesh(cubePath, 1.0, 0.0, 0.0, m, err), "load cube.stl");
        check(m.comps.size() == 1, "cube: 1 component");
        check(m.tris.size() == 12, "cube: 12 triangles");
        check(m.pts.size() == 8, "cube: 8 vertices");
        if (!m.comps.empty()) {
            check(m.comps[0].compEdges.size() == 18, "cube: 18 edges");
            check(m.comps[0].open == 0, "cube: open==0");
            check(m.comps[0].conflict == 0, "cube: conflict==0");
            check(m.comps[0].nonManifold == 0, "cube: nonManifold==0");
            check(m.comps[0].clean(), "cube: clean");
            const double vol = std::fabs(m.comps[0].vol);
            const double rel = std::fabs(vol - 1000.0) / 1000.0;
            check(rel < 1e-9, "cube: |volume| == 1000 within 1e-9 relative");
            checkLocalIds(m.comps[0], "cube");
        }
    }

    if (gFails) {
        std::fprintf(stderr, "%d assertion(s) failed\n", gFails);
        return 1;
    }
    std::printf("ok\n");
    return 0;
}
