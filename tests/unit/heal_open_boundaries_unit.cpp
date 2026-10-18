// healOpenBoundaries (D-fc-5.10) on synthetic meshes:
//  - an open unit box plus a detached lid whose rim is a jittered duplicate
//    (inside stitchTol) and whose front edge carries a T-junction closes to
//    one closed, consistently wound mesh; through convert() it is one solid on
//    the TrueForm path and stays two open shells on the Verbatim path;
//  - a closed mesh is a byte no-op;
//  - a weld that would make an edge of valence 3 is refused, and only that
//    cluster keeps its input (a healable cluster beside it still heals);
//  - a lid offset by more than stitchTol stays open;
//  - a T-vertex 3·tau off the edge is not split;
//  - an unwelded sheet whose boundary vertices all share one x (the case an
//    x-sorted sweep handles in quadratic time) gets every weld and T split;
//  - a non-finite vertex coordinate or edge length refuses the heal, tris
//    untouched, even when a disjoint cluster would otherwise close (D-fc-5.19);
//  - the owning apex of an open edge is discarded alone, so a later vertex on
//    that edge still splits it (D-fc-5.19).
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
#include <limits>
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

static bool hasEdge(const Tris& tris, int a, int b) {
    for (const auto& t : tris)
        for (int s = 0; s < 3; ++s) {
            const int u = t[s], v = t[(s + 1) % 3];
            if ((u == a && v == b) || (u == b && v == a)) return true;
        }
    return false;
}

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

    // ---- 8. unwelded sheet in one x slab: every weld and every T split -------
    // An open sheet on x = 0 (lower half n x n/2 cells, upper half with the y
    // step halved, so the seam carries n T-vertices) where every triangle owns
    // its vertex copies, jittered by up to tau/8 in y and z. All boundary
    // vertices share one x, so an x-window holds every one of them; the grid
    // index must still find each weld partner and each T-vertex.
    {
        const int n = 32;
        const double h = 1.0 / n;
        uint32_t lcg = 12345u;
        auto jitter = [&]() {
            lcg = lcg * 1664525u + 1013904223u;
            return (((double)(lcg >> 8) / (double)(1u << 24)) * 2.0 - 1.0) * (tau / 8);
        };
        Pts pts;
        Tris tris;
        auto addTri = [&](std::array<double, 2> a, std::array<double, 2> b, std::array<double, 2> c) {
            const int base = (int)pts.size();
            for (const auto& p : { a, b, c }) pts.push_back({ 0.0, p[0] + jitter(), p[1] + jitter() });
            tris.push_back({ base, base + 1, base + 2 });
        };
        auto addGrid = [&](int ni, double hy, int j0, int j1) {
            for (int i = 0; i < ni; ++i)
                for (int j = j0; j < j1; ++j) {
                    const std::array<double, 2> a{ i * hy, j * h }, b{ (i + 1) * hy, j * h },
                        c{ (i + 1) * hy, (j + 1) * h }, d{ i * hy, (j + 1) * h };
                    addTri(a, c, b);
                    addTri(a, d, c);
                }
        };
        addGrid(n, h, 0, n / 2);
        addGrid(2 * n, h / 2, n / 2, n);
        const int distinct = (n + 1) * (n / 2 + 1) + (2 * n + 1) * (n / 2 + 1) - (n + 1);
        const HealStats st = healOpenBoundaries(tris, pts, q, sewTol, false);
        const EdgeHealth after = health(tris);
        check(st.boundaryBefore == 3 * 3 * n * n, "slab: every edge of the unwelded sheet is open");
        check(st.welds == 3 * 3 * n * n - distinct, "slab: every copy welds onto its lowest-id twin");
        check(st.splits == n, "slab: every seam T-vertex splits its edge");
        check(st.clustersHealed == 1 && st.clustersRefused == 0, "slab: one stitched cluster heals");
        check(after.boundary == 5 * n && st.boundaryAfter == 5 * n,
              "slab: only the sheet's outer rim stays open (no capping)");
        check(after.nonManifold == 0 && after.conflict == 0, "slab: 2-manifold, consistent winding");
    }

    // ---- 9. non-finite input refuses the heal (D-fc-5.19 B1) -----------------
    // The review mesh: an open box whose lid is jittered by 3·tau (it welds and
    // closes on its own) plus a disjoint flap. A finite flap leaves the lid
    // free to heal. A NaN or infinite coordinate, or a finite coordinate whose
    // edge length overflows, refuses the whole heal before any weld.
    {
        const double qB = 1e-4;
        const double tauB = 2.0 * qB;
        const double sewB = 1e-3;
        Pts base;
        Tris baseTris;
        openBoxWithLid(base, baseTris, { 3 * tauB, -2 * tauB, tauB }, 0.4 * tauB);
        {
            Pts pts = base;
            Tris tris = baseTris;
            const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
            const EdgeHealth h = health(tris);
            check(st.clustersHealed == 1 && h.boundary == 0 && h.nonManifold == 0 && h.conflict == 0,
                  "B1 control: jittered lid closes");
        }
        {
            Pts pts = base;
            Tris tris = baseTris;
            const int b = (int)pts.size();
            pts.push_back({ 0, 10, 0 });
            pts.push_back({ 1, 10, 0 });
            pts.push_back({ 0.5, 10, 1 });
            tris.push_back({ b, b + 1, b + 2 });
            const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
            bool flapKept = false;
            for (const auto& t : tris)
                if (t[0] == b && t[1] == b + 1 && t[2] == b + 2) flapKept = true;
            check(st.clustersHealed == 1 && st.clustersRefused == 1 && flapKept,
                  "B1 finite flap: lid heals, disjoint flap keeps its triangle");
        }
        {
            Pts pts = base;
            Tris tris = baseTris;
            const int b = (int)pts.size();
            pts.push_back({ 0, 10, 0 });
            pts.push_back({ 1, 10, 0 });
            pts.push_back({ std::numeric_limits<double>::quiet_NaN(), 10, 0 });
            tris.push_back({ b, b + 1, b + 2 });
            const Tris in = tris;
            const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
            check(tris == in && st.clustersHealed == 0 && st.welds == 0 && st.splits == 0,
                  "B1 NaN coordinate: heal refused, lid not welded");
        }
        {
            Pts pts = base;
            Tris tris = baseTris;
            const int b = (int)pts.size();
            pts.push_back({ 0, 10, 0 });
            pts.push_back({ 1, 10, 0 });
            pts.push_back({ std::numeric_limits<double>::infinity(), 10, 0 });
            tris.push_back({ b, b + 1, b + 2 });
            const Tris in = tris;
            const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
            check(tris == in && st.clustersHealed == 0 && st.welds == 0 && st.splits == 0,
                  "B1 +Inf coordinate: heal refused");
        }
        {
            Pts pts = base;
            Tris tris = baseTris;
            const int b = (int)pts.size();
            const double big = std::numeric_limits<double>::max();
            pts.push_back({ 0, 10, 0 });
            pts.push_back({ big, 10, 0 });
            pts.push_back({ 0, 11, 0 });
            tris.push_back({ b, b + 1, b + 2 });
            const Tris in = tris;
            const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
            check(tris == in && st.clustersHealed == 0 && st.welds == 0 && st.splits == 0,
                  "B1 non-finite edge length: heal refused");
        }
    }

    // ---- 10. owning apex does not mask a later on-edge split (D-fc-5.19 B2) --
    // A(0,0,0) B(1,0,0) C(0.25, 1e-4, 0): C is the owning apex and the minimum-t
    // candidate on AB (t = 0.25, height 1e-4 <= tau). D(0.75, 0, 0) lies on AB
    // at t = 0.75. Triangle B-D-E owns DB, so inserting D closes that edge and
    // the certificate keeps the split. E-D-F is the disjoint-from-AB flap.
    {
        const double qB = 1e-4;
        const double tauB = 2.0 * qB;
        const double sewB = 1e-3;
        const int A = 0, B = 1, C = 2, D = 3, E = 4, F = 5;
        Pts pts = {
            { 0, 0, 0 }, { 1, 0, 0 }, { 0.25, 1e-4, 0 },
            { 0.75, 0, 0 }, { 0.75, 0.4, 0 }, { 0.75, 0.4, 0.3 },
        };
        Tris tris = { { A, B, C }, { B, D, E }, { E, D, F } };
        check(hasEdge(tris, A, B), "B2: fixture has edge AB");
        const HealStats st = healOpenBoundaries(tris, pts, qB, sewB, false);
        const EdgeHealth h = health(tris);
        check(st.splits == 1 && st.welds == 0, "B2: D is the one split; nothing welds");
        check(st.clustersHealed == 1 && st.clustersRefused == 0, "B2: the cluster heals");
        check(st.boundaryBefore == 7 && st.boundaryAfter == 6, "B2: open edges 7 -> 6");
        check(h.boundary == 6 && h.nonManifold == 0 && h.conflict == 0,
              "B2: manifold, consistent winding");
        check(hasEdge(tris, A, D) && hasEdge(tris, D, B) && !hasEdge(tris, A, B),
              "B2: D splits AB; the apex did not abandon the edge");
        bool distinct = true;
        for (const auto& t : tris)
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) distinct = false;
        check(distinct && tris.size() == 4, "B2: one added triangle, none degenerate");
        check(st.splitMax <= tauB, "B2: the split lies within tau of its edge");
    }

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "heal_open_boundaries_unit ok\n");
    return 0;
}
