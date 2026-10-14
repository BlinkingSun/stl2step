// Closed box with one zero-area T-junction. After absorb: same triangle
// count, 0 boundary edges, 0 zero-area triangles. A vertex 3·tau off the
// opposite edge is not absorbed. quantFloor <= 0 absorbs nothing.
// tau = 2·quantFloor (D-130-12). The engine predicate is the one in
// detail::triangleAreaDegenerate (mag < 1e-12 or mag^2 < l2^2 * 1e-20).
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace stl2step {
namespace detail {
int absorbZeroAreaTJunctions(std::vector<std::array<int, 3>>& tris,
                             const std::vector<std::array<double, 3>>& pts,
                             double quantFloor, bool diag);
}
}

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (ok) { std::fprintf(stderr, "PASS %s\n", name); return; }
    std::fprintf(stderr, "FAIL %s\n", name);
    ++gFail;
}

static bool areaDegenerate(const std::array<double, 3>& A,
                           const std::array<double, 3>& B,
                           const std::array<double, 3>& C) {
    const double bx = B[0] - A[0], by = B[1] - A[1], bz = B[2] - A[2];
    const double cx = C[0] - A[0], cy = C[1] - A[1], cz = C[2] - A[2];
    const double nx = by * cz - bz * cy;
    const double ny = bz * cx - bx * cz;
    const double nz = bx * cy - by * cx;
    const double mag = std::sqrt(nx * nx + ny * ny + nz * nz);
    const double lAB = bx * bx + by * by + bz * bz;
    const double lAC = cx * cx + cy * cy + cz * cz;
    const double l2 = lAB > lAC ? lAB : lAC;
    return mag < 1e-12 || mag * mag < l2 * l2 * 1e-20;
}

struct EdgeHealth { int boundary = 0; int nonManifold = 0; int conflict = 0; int zeroArea = 0; };

static uint64_t ekey(int a, int b) {
    return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                 : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
}

static EdgeHealth health(const std::vector<std::array<int, 3>>& tris,
                         const std::vector<std::array<double, 3>>& pts) {
    EdgeHealth h;
    struct Rec { int cnt = 0; int fwd = 0; };
    std::vector<std::pair<uint64_t, Rec>> edges;
    auto find = [&](uint64_t k) -> Rec& {
        for (auto& e : edges) if (e.first == k) return e.second;
        edges.push_back({ k, {} });
        return edges.back().second;
    };
    for (const auto& t : tris) {
        if (areaDegenerate(pts[(size_t)t[0]], pts[(size_t)t[1]], pts[(size_t)t[2]]))
            h.zeroArea++;
        const int v[4] = { t[0], t[1], t[2], t[0] };
        for (int s = 0; s < 3; ++s) {
            Rec& e = find(ekey(v[s], v[s + 1]));
            e.cnt++;
            if (v[s] < v[s + 1]) e.fwd++;
        }
    }
    for (const auto& e : edges) {
        if (e.second.cnt == 1) h.boundary++;
        else if (e.second.cnt > 2) h.nonManifold++;
        else if (e.second.fwd != 1) h.conflict++;
    }
    return h;
}

// Unit cube. Top diagonal 4-6 carries midpoint vertex 8 and a zero-area
// triangle; the other side of that diagonal stays one triangle (the neighbour).
static std::vector<std::array<double, 3>> boxPoints(bool offEdge, double tau) {
    std::vector<std::array<double, 3>> pts = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
        {0.5, 0.5, 1},
    };
    if (offEdge) {
        const double off = 3.0 * tau;
        const double inv = 1.0 / std::sqrt(2.0);
        pts[8][0] += off * inv;
        pts[8][1] -= off * inv;
    }
    return pts;
}

static std::vector<std::array<int, 3>> boxTris() {
    return {
        {0, 2, 1}, {0, 3, 2},
        {4, 5, 6},
        {4, 6, 8},
        {4, 8, 7},
        {8, 6, 7},
        {0, 1, 5}, {0, 5, 4},
        {3, 7, 6}, {3, 6, 2},
        {0, 4, 7}, {0, 7, 3},
        {1, 2, 6}, {1, 6, 5},
    };
}

static void putU32(std::string& b, uint32_t v) {
    for (int i = 0; i < 4; ++i) b.push_back((char)((v >> (8 * i)) & 0xFFu));
}
static void putF32(std::string& b, float v) {
    uint32_t u;
    std::memcpy(&u, &v, 4);
    putU32(b, u);
}

static bool writeStl(const std::string& path,
                     const std::vector<std::array<int, 3>>& tris,
                     const std::vector<std::array<double, 3>>& pts) {
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

    // q is a fixture small enough that a vertex 3·tau off a  √2 mm edge is
    // still area-degenerate (mag < 1e-12), so the miss is the tau gate.
    const double q = 1e-14;
    const double tau = 2.0 * q;

    auto ptsOn = boxPoints(false, tau);
    auto trisOn = boxTris();
    const int nBefore = (int)trisOn.size();
    const EdgeHealth before = health(trisOn, ptsOn);
    check(before.boundary == 0 && before.nonManifold == 0 && before.conflict == 0,
          "fixture: closed box, consistent winding");
    check(before.zeroArea == 1, "fixture: exactly one zero-area triangle");

    const int nOn = stl2step::detail::absorbZeroAreaTJunctions(trisOn, ptsOn, q, false);
    const EdgeHealth after = health(trisOn, ptsOn);
    check(nOn == 1, "on-edge: one absorb");
    check((int)trisOn.size() == nBefore, "on-edge: triangle count unchanged");
    check(after.boundary == 0, "on-edge: 0 boundary edges");
    check(after.nonManifold == 0, "on-edge: 0 non-manifold edges");
    check(after.conflict == 0, "on-edge: winding preserved");
    check(after.zeroArea == 0, "on-edge: 0 zero-area triangles");

    auto trisAgain = trisOn;
    const int nAgain = stl2step::detail::absorbZeroAreaTJunctions(trisAgain, ptsOn, q, false);
    check(nAgain == 0, "on-edge: second pass absorbs nothing");

    auto ptsOff = boxPoints(true, tau);
    auto trisOff = boxTris();
    const EdgeHealth offBefore = health(trisOff, ptsOff);
    check(offBefore.zeroArea == 1, "3·tau: triangle is still zero-area");
    {
        const auto& P = ptsOff[8];
        const auto& A = ptsOff[4];
        const auto& B = ptsOff[6];
        const double abx = B[0] - A[0], aby = B[1] - A[1], abz = B[2] - A[2];
        const double ab2 = abx * abx + aby * aby + abz * abz;
        const double param = ((P[0] - A[0]) * abx + (P[1] - A[1]) * aby + (P[2] - A[2]) * abz) / ab2;
        const double dx = P[0] - (A[0] + abx * param);
        const double dy = P[1] - (A[1] + aby * param);
        const double dz = P[2] - (A[2] + abz * param);
        const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        check(param > 0.0 && param < 1.0, "3·tau: projection is on the open segment");
        // |coord| = 1, offset = 3·tau = 6e-14. Adding that onto 0.5 rounds by
        // ~2e-17 mm (measured). Still 3·tau, not 1·tau and not 4·tau.
        check(std::fabs(dist / (3.0 * tau) - 1.0) < 1e-3, "3·tau: distance is 3·tau");
        check(dist > tau, "3·tau: distance is outside tau");
    }
    const int nOff = stl2step::detail::absorbZeroAreaTJunctions(trisOff, ptsOff, q, false);
    check(nOff == 0, "3·tau: nothing absorbed");
    check((int)trisOff.size() == nBefore, "3·tau: triangle count unchanged");
    check(health(trisOff, ptsOff).zeroArea == 1, "3·tau: zero-area triangle remains");

    auto trisQ = boxTris();
    const int nQ = stl2step::detail::absorbZeroAreaTJunctions(trisQ, ptsOn, 0.0, false);
    check(nQ == 0, "q <= 0: absorb nothing");
    check(health(trisQ, ptsOn).zeroArea == 1, "q <= 0: zero-area triangle remains");

    const std::string stl = scratch + "/tjunction_box.stl";
    const std::string step = scratch + "/tjunction_box.step";
    auto ptsFile = boxPoints(false, tau);
    auto trisFile = boxTris();
    check(writeStl(stl, trisFile, ptsFile), "wrote synthetic STL");
    stl2step::Options opt;
    opt.input = stl;
    opt.output = step;
    opt.threads = 1;
    opt.smooth = false;
    opt.verify = true;
    const stl2step::Result r = stl2step::convert(opt);
    check(r.ok, "convert ok");
    if (!r.ok) std::fprintf(stderr, "  error: %s\n", r.error.c_str());
    check(r.solids == 1, "convert: 1 solid");
    check(r.openShells == 0, "convert: 0 open shells");
    check(r.watertight, "convert: watertight");
    const std::string js = r.toJson();
    check(js.find("\"meshTJunctionsAbsorbed\":1") != std::string::npos,
          "RESULT meshTJunctionsAbsorbed=1");

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "tjunction_absorb_test ok\n");
    return 0;
}
