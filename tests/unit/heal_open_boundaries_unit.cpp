// healOpenBoundaries (D-fc-5.10) on synthetic meshes:
//  - an open unit box plus a detached lid whose rim is a jittered duplicate
//    (inside stitchTol) and whose front edge carries a T-junction closes to
//    one closed, consistently wound mesh; through convert() it is one solid on
//    the TrueForm path and stays two open shells on the Verbatim path;
//  - a closed mesh is a byte no-op;
//  - a weld that would make an edge of valence 3 is refused, and only that
//    cluster keeps its input (a healable cluster beside it still heals);
//  - a lid offset by more than stitchTol stays open;
//  - a T-vertex 3·tau off the edge is not split.
// tau = 2·quantFloor (D-130-12); stitchTol = min(sewTol, sMin / 2).
// SPDX-License-Identifier: MIT

#include "mesh_heal.hpp"
#include "stl2step/stl2step.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>

using stl2step::detail::HealStats;
using stl2step::detail::healOpenBoundaries;
using Tris = std::vector<std::array<int, 3>>;
using Pts = std::vector<std::array<double, 3>>;

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (ok) { std::fprintf(stderr, "PASS %s\n", name); return; }
    std::fprintf(stderr, "FAIL %s\n", name);
    ++gFail;
}

struct EdgeHealth { int boundary = 0; int nonManifold = 0; int conflict = 0; };

static EdgeHealth health(const Tris& tris) {
    struct Rec { int cnt = 0; int fwd = 0; };
    std::map<std::pair<int, int>, Rec> edges;
    for (const auto& t : tris)
        for (int s = 0; s < 3; ++s) {
            const int a = t[s], b = t[(s + 1) % 3];
            Rec& e = edges[{ std::min(a, b), std::max(a, b) }];
            e.cnt++;
            if (a < b) e.fwd++;
        }
    EdgeHealth h;
    for (const auto& kv : edges) {
        if (kv.second.cnt == 1) h.boundary++;
        else if (kv.second.cnt > 2) h.nonManifold++;
        else if (kv.second.fwd != 1) h.conflict++;
    }
    return h;
}

// Unit box without its top (vertices 0..7, top rim 4..7 at z = 1) and a lid
// (vertices 8..11 = rim copies offset by `lidOff`, vertex 12 on the lid's
// front edge 8-9 offset by `tOff` along -y). Lid winding matches the box.
static void openBoxWithLid(Pts& pts, Tris& tris, const std::array<double, 3>& lidOff,
                           double tOff, const std::array<double, 3>& origin = { 0, 0, 0 }) {
    const int base = (int)pts.size();
    const double box[8][3] = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
        { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 },
    };
    for (const auto& p : box)
        pts.push_back({ origin[0] + p[0], origin[1] + p[1], origin[2] + p[2] });
    for (int i = 4; i < 8; ++i)
        pts.push_back({ origin[0] + box[i][0] + lidOff[0], origin[1] + box[i][1] + lidOff[1],
                        origin[2] + box[i][2] + lidOff[2] });
    pts.push_back({ origin[0] + 0.5, origin[1] - tOff, origin[2] + 1.0 });
    const int side[10][3] = {
        { 0, 2, 1 }, { 0, 3, 2 },              // bottom
        { 0, 1, 5 }, { 0, 5, 4 },              // y = 0
        { 3, 7, 6 }, { 3, 6, 2 },              // y = 1
        { 0, 4, 7 }, { 0, 7, 3 },              // x = 0
        { 1, 2, 6 }, { 1, 6, 5 },              // x = 1
    };
    for (const auto& t : side) tris.push_back({ base + t[0], base + t[1], base + t[2] });
    const int lid[3][3] = { { 8, 12, 11 }, { 12, 9, 10 }, { 12, 10, 11 } };
    for (const auto& t : lid) tris.push_back({ base + t[0], base + t[1], base + t[2] });
}

static void closedBox(Pts& pts, Tris& tris) {
    pts = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
            { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
    tris = { { 0, 2, 1 }, { 0, 3, 2 }, { 4, 5, 6 }, { 4, 6, 7 }, { 0, 1, 5 }, { 0, 5, 4 },
             { 3, 7, 6 }, { 3, 6, 2 }, { 0, 4, 7 }, { 0, 7, 3 }, { 1, 2, 6 }, { 1, 6, 5 } };
}

// Three triangles hinged on one jittered edge: every weld makes it valence 3.
static void threeFlaps(Pts& pts, Tris& tris, const std::array<double, 3>& origin, double jit) {
    const int b = (int)pts.size();
    auto P = [&](double x, double y, double z) {
        pts.push_back({ origin[0] + x, origin[1] + y, origin[2] + z });
    };
    P(0, 0, 0); P(1, 0, 0); P(0.5, 1, 0);                  // b+0..2
    P(jit, 0, 0); P(1 + jit, 0, 0); P(0.5, -1, 0);         // b+3..5
    P(0, jit, 0); P(1, jit, 0); P(0.5, 0, 1);              // b+6..8
    tris.push_back({ b + 0, b + 1, b + 2 });
    tris.push_back({ b + 4, b + 3, b + 5 });
    tris.push_back({ b + 7, b + 6, b + 8 });
}

static void putU32(std::string& b, uint32_t v) {
    for (int i = 0; i < 4; ++i) b.push_back((char)((v >> (8 * i)) & 0xFFu));
}
static void putF32(std::string& b, float v) {
    uint32_t u;
    std::memcpy(&u, &v, 4);
    putU32(b, u);
}

static bool writeStl(const std::string& path, const Tris& tris, const Pts& pts) {
    std::string b(80, '\0');
    putU32(b, (uint32_t)tris.size());
    for (const auto& t : tris) {
        putF32(b, 0.f); putF32(b, 0.f); putF32(b, 0.f);
        for (int k = 0; k < 3; ++k) {
            const auto& p = pts[(size_t)t[k]];
            putF32(b, (float)p[0]); putF32(b, (float)p[1]); putF32(b, (float)p[2]);
        }
        b.push_back('\0'); b.push_back('\0');
    }
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = std::fwrite(b.data(), 1, b.size(), f) == b.size();
    std::fclose(f);
    return ok;
}

int main(int argc, char** argv) {
    const std::string scratch = argc > 1 ? argv[1] : std::string(".");
    const double q = 1e-7;
    const double tau = 2.0 * q;
    const double sewTol = 1e-3;   // the repair path's tolerance, as the caller passes it

    // ---- 1. jittered lid + T-junction closes ---------------------------------
    {
        Pts pts;
        Tris tris;
        openBoxWithLid(pts, tris, { 3 * tau, -2 * tau, tau }, 0.4 * tau);
        const EdgeHealth before = health(tris);
        check(before.boundary == 9 && before.nonManifold == 0 && before.conflict == 0,
              "lid: fixture has 9 open edges, no non-manifold, no conflict");
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        const EdgeHealth after = health(tris);
        check(st.boundaryBefore == 9 && st.boundaryAfter == 0, "lid: 9 open edges -> 0");
        check(after.boundary == 0 && after.nonManifold == 0 && after.conflict == 0,
              "lid: closed 2-manifold, consistent winding");
        check(st.welds == 4, "lid: 4 welds (lid rim onto the box rim)");
        check(st.splits == 1, "lid: 1 vertex-on-edge split (the T-junction)");
        check(st.collapsed == 0, "lid: nothing collapsed");
        check(tris.size() == 14, "lid: 13 + 1 split = 14 triangles");
        check(st.clustersHealed == 1 && st.clustersRefused == 0, "lid: one cluster healed");
        check(std::fabs(st.sMin - 0.5) <= 4 * tau, "lid: sMin = shortest edge > tau (lid half-edge, 0.5 less jitter)");
        check(st.stitchTol == sewTol, "lid: stitchTol = min(sewTol, sMin / 2) = sewTol");
        check(st.weldMax > tau && st.weldMax <= st.stitchTol, "lid: welds beyond tau, within stitchTol");
        check(st.splitMax <= tau, "lid: split vertex within tau of its edge");
        bool repsAreBox = true;
        for (const auto& t : tris)
            for (int v : t)
                if (v >= 8 && v <= 11) repsAreBox = false;
        check(repsAreBox, "lid: representative = lowest id (the box's own rim points)");
        Tris again = tris;
        const HealStats st2 = healOpenBoundaries(again, pts, q, sewTol, false);
        check(again == tris && st2.boundaryBefore == 0, "lid: second pass is a no-op");
    }

    // ---- 2. closed mesh: byte no-op ------------------------------------------
    {
        Pts pts;
        Tris tris;
        closedBox(pts, tris);
        const Tris in = tris;
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        check(tris == in, "closed: triangles byte-identical");
        check(st.boundaryBefore == 0 && st.clustersHealed == 0 && st.welds == 0 && st.splits == 0,
              "closed: no open edge, no work");
    }

    // ---- 3. valence-3 weld refused, per cluster ------------------------------
    {
        Pts pts;
        Tris tris;
        threeFlaps(pts, tris, { 0, 0, 0 }, 3 * tau);
        const Tris in = tris;
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        check(tris == in, "valence-3: refused, input unchanged");
        check(st.clustersHealed == 0 && st.clustersRefused == 1, "valence-3: cluster refused");
        check(st.boundaryAfter == st.boundaryBefore, "valence-3: open edges unchanged");

        Pts pts2;
        Tris tris2;
        threeFlaps(pts2, tris2, { 0, 0, 0 }, 3 * tau);
        const size_t nFlaps = tris2.size();
        openBoxWithLid(pts2, tris2, { 3 * tau, -2 * tau, tau }, 0.4 * tau, { 10, 0, 0 });
        const Tris flapsIn(tris2.begin(), tris2.begin() + (long)nFlaps);
        const HealStats st2 = healOpenBoundaries(tris2, pts2, q, sewTol, false);
        check(st2.clustersHealed == 1 && st2.clustersRefused == 1,
              "per-cluster: box heals, flaps refused");
        check(Tris(tris2.begin(), tris2.begin() + (long)nFlaps) == flapsIn,
              "per-cluster: refused flaps keep their input triangles");
        Tris boxOnly(tris2.begin() + (long)nFlaps, tris2.end());
        const EdgeHealth hb = health(boxOnly);
        check(hb.boundary == 0 && hb.nonManifold == 0 && hb.conflict == 0,
              "per-cluster: the box beside them closes");
        check(st2.boundaryAfter == 9, "per-cluster: only the flaps' 9 open edges remain");
    }

    // ---- 4. a crack wider than stitchTol stays open --------------------------
    {
        Pts pts;
        Tris tris;
        openBoxWithLid(pts, tris, { 0, 0, 5 * sewTol }, 0.0);
        const Tris in = tris;
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        check(tris == in, "wide crack: input unchanged");
        check(st.welds == 0 && st.splits == 0 && st.clustersHealed == 0, "wide crack: no weld, no split");
        check(st.boundaryAfter == 9, "wide crack: 9 open edges remain");
    }

    // ---- 5. a T-vertex 3·tau off the edge is not split -----------------------
    {
        Pts pts;
        Tris tris;
        openBoxWithLid(pts, tris, { 3 * tau, -2 * tau, tau }, 3 * tau);
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        check(st.splits == 0, "3·tau T-vertex: no split");
        check(st.welds == 4 && st.boundaryAfter == 3,
              "3·tau T-vertex: rim welds still close 6 of 9 open edges");
    }

    // ---- 6. quantFloor <= 0 does nothing -------------------------------------
    {
        Pts pts;
        Tris tris;
        openBoxWithLid(pts, tris, { 3 * tau, -2 * tau, tau }, 0.0);
        const Tris in = tris;
        const HealStats st = healOpenBoundaries(tris, pts, 0.0, sewTol, false);
        check(tris == in && st.boundaryBefore == 0, "q <= 0: nothing");
    }

    // ---- 7. convert(): TrueForm closes the lid, Verbatim does not -------------
    {
        // Float32 file: rim jitter 4e-6 is above the file's tau (~2e-7 at |x| = 1)
        // and below its auto sewTol (diag·1e-5 ~ 1.7e-5).
        Pts pts;
        Tris tris;
        openBoxWithLid(pts, tris, { 4e-6, -4e-6, 4e-6 }, 0.0);
        const std::string stl = scratch + "/heal_open_lid.stl";
        check(writeStl(stl, tris, pts), "wrote synthetic STL");

        stl2step::Options on;
        on.input = stl;
        on.output = scratch + "/heal_open_lid_on.step";
        on.threads = 1;
        on.smooth = true;
        const stl2step::Result r = stl2step::convert(on);
        if (!r.ok) std::fprintf(stderr, "  error: %s\n", r.error.c_str());
        check(r.ok && r.solids == 1 && r.openShells == 0 && r.watertight && r.components == 1,
              "convert TrueForm: one watertight solid");
        check(r.smoothSkippedComponents == 0, "convert TrueForm: component not skipped as dirty");

        stl2step::Options off = on;
        off.output = scratch + "/heal_open_lid_off.step";
        off.smooth = false;
        const stl2step::Result v = stl2step::convert(off);
        check(v.ok && v.components == 2 && !v.watertight, "convert Verbatim: unchanged (no heal)");
    }

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "heal_open_boundaries_unit ok\n");
    return 0;
}
