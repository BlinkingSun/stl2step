// tess-law unit cases (SPEC-tess-law). Synthetic meshes only.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

using stl2step::refit::DerivedTols;
using stl2step::refit::MeshView;
using stl2step::refit::Origin;
using stl2step::refit::ProvClaim;
using stl2step::refit::Provisional;
using stl2step::refit::Region;
using stl2step::refit::SegmentParams;
using stl2step::refit::SegmentWork;
using stl2step::refit::SurfType;
using stl2step::refit::refuteNgonByLaw;

static const double kPi = 3.14159265358979323846264338327950288;
static const double kTwoPi = 2.0 * kPi;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    if (ok) {
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    }
}

// One component, local id == global id. Tubes are open prism walls: two
// triangles per side, vertices exactly on the cylinder.
struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris;
    std::vector<int> compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::uint8_t> triDirs;

    struct Tube {
        int firstTri = 0;
        int nTri = 0;
        int firstVert = 0;
        int N = 0;
        double R = 0;
        double z0 = 0, z1 = 0;
        double cx = 0, cy = 0;
    };

    Tube addTube(double R, int N, double z0, double z1, double cx, double cy) {
        Tube tb;
        tb.firstVert = static_cast<int>(pts.size());
        tb.firstTri = static_cast<int>(triFlat.size() / 3);
        tb.N = N;
        tb.R = R;
        tb.z0 = z0;
        tb.z1 = z1;
        tb.cx = cx;
        tb.cy = cy;
        for (int i = 0; i < N; ++i) {
            const double a = kTwoPi * static_cast<double>(i) / static_cast<double>(N);
            pts.push_back(gp_XYZ(cx + R * std::cos(a), cy + R * std::sin(a), z0));
        }
        for (int i = 0; i < N; ++i) {
            const double a = kTwoPi * static_cast<double>(i) / static_cast<double>(N);
            pts.push_back(gp_XYZ(cx + R * std::cos(a), cy + R * std::sin(a), z1));
        }
        for (int i = 0; i < N; ++i) {
            const int j = (i + 1) % N;
            const int b0 = tb.firstVert + i;
            const int b1 = tb.firstVert + j;
            const int t0 = tb.firstVert + N + i;
            const int t1 = tb.firstVert + N + j;
            triFlat.push_back(b0);
            triFlat.push_back(b1);
            triFlat.push_back(t1);
            triFlat.push_back(b0);
            triFlat.push_back(t1);
            triFlat.push_back(t0);
        }
        tb.nTri = N * 2;
        return tb;
    }

    // Open strip. angA / angB are increasing azimuths (radians) of the two
    // boundary rails at z0 and z1. Interior edges are shared; the end caps are
    // the only other boundary edges.
    Tube addStrip(double R, const std::vector<double>& angA, const std::vector<double>& angB,
                  double z0, double z1, double cx, double cy) {
        Tube tb;
        tb.firstVert = static_cast<int>(pts.size());
        tb.firstTri = static_cast<int>(triFlat.size() / 3);
        tb.N = 0;
        tb.R = R;
        tb.z0 = z0;
        tb.z1 = z1;
        tb.cx = cx;
        tb.cy = cy;
        const int nA = static_cast<int>(angA.size());
        const int nB = static_cast<int>(angB.size());
        const int baseA = tb.firstVert;
        for (double a : angA)
            pts.push_back(gp_XYZ(cx + R * std::cos(a), cy + R * std::sin(a), z0));
        const int baseB = static_cast<int>(pts.size());
        for (double a : angB)
            pts.push_back(gp_XYZ(cx + R * std::cos(a), cy + R * std::sin(a), z1));
        int i = 0, j = 0;
        auto tri = [&](int a, int b, int c) {
            triFlat.push_back(a);
            triFlat.push_back(b);
            triFlat.push_back(c);
        };
        while (i < nA - 1 || j < nB - 1) {
            if (i == nA - 1) {
                tri(baseA + i, baseB + j, baseB + j + 1);
                ++j;
            } else if (j == nB - 1) {
                tri(baseA + i, baseA + i + 1, baseB + j);
                ++i;
            } else if (angA[static_cast<size_t>(i + 1)] <= angB[static_cast<size_t>(j + 1)]) {
                tri(baseA + i, baseA + i + 1, baseB + j);
                ++i;
            } else {
                tri(baseA + i, baseB + j, baseB + j + 1);
                ++j;
            }
        }
        tb.nTri = static_cast<int>(triFlat.size() / 3) - tb.firstTri;
        return tb;
    }

    // Move one bottom vertex radially by `delta` millimetres (breaks τ).
    void bump(const Tube& tb, int corner, double delta) {
        gp_XYZ& p = pts[static_cast<size_t>(tb.firstVert + corner)];
        const double x = p.X() - tb.cx;
        const double y = p.Y() - tb.cy;
        const double m = std::hypot(x, y);
        if (!(m > 0.0)) return;
        p.SetX(tb.cx + x * (1.0 + delta / m));
        p.SetY(tb.cy + y * (1.0 + delta / m));
    }

    MeshView view(double quantFloor) {
        const int nT = static_cast<int>(triFlat.size() / 3);
        compTris.resize(static_cast<size_t>(nT));
        for (int i = 0; i < nT; ++i) compTris[static_cast<size_t>(i)] = i;
        compVtx.resize(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) compVtx[i] = static_cast<int>(i);
        edges.clear();
        triDirs.assign(static_cast<size_t>(nT), 0);
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int(*)[3]>(triFlat.data());
        mv.compTris = compTris.data();
        mv.compVtx = compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = nullptr;
        mv.triDirs = triDirs.data();
        mv.nTri = static_cast<size_t>(nT);
        mv.nVtx = pts.size();
        mv.nEdge = 0;
        mv.diag = 40.0;
        mv.weldTol = 0.0;
        mv.sewTol = 0.0;
        mv.quantFloor = quantFloor;
        return mv;
    }
};

static Region makeCyl(const Synth::Tube& tb, Origin origin, bool closed) {
    Region r;
    r.type = SurfType::Cylinder;
    r.origin = origin;
    r.radius = tb.R;
    r.nSides = tb.N;
    r.closed360 = closed;
    r.ax = gp_Ax3(gp_Pnt(tb.cx, tb.cy, 0.5 * (tb.z0 + tb.z1)), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.vMin = tb.z0;
    r.vMax = tb.z1;
    r.tris.reserve(static_cast<size_t>(tb.nTri));
    for (int i = 0; i < tb.nTri; ++i) r.tris.push_back(tb.firstTri + i);
    return r;
}

// Production claimNgonWallsA: ConsumedCylinder, triangles cleared, A2 plane kept.
static void addClearedWalls(SegmentWork& work, const Synth::Tube& tb) {
    for (int i = 0; i < tb.N; ++i) {
        const double a0 = kTwoPi * static_cast<double>(i) / static_cast<double>(tb.N);
        const double a1 = kTwoPi * static_cast<double>(i + 1) / static_cast<double>(tb.N);
        const double am = 0.5 * (a0 + a1);
        Provisional p;
        p.claim = ProvClaim::ConsumedCylinder;
        p.maxVertexDev = 0.0;
        p.area = tb.R * std::sin(kPi / static_cast<double>(tb.N)) * (tb.z1 - tb.z0);
        const gp_XYZ loc(tb.cx + tb.R * std::cos(a0), tb.cy + tb.R * std::sin(a0), tb.z0);
        p.plane = gp_Ax3(gp_Pnt(loc.X(), loc.Y(), loc.Z()),
                         gp_Dir(std::cos(am), std::sin(am), 0.0));
        work.provisionals.push_back(p);
    }
}

struct Outcome {
    int ngon6 = 0;
    int ngon8 = 0;
    int ngon24 = 0;
    int ngon64 = 0;
    int ngon126 = 0;
    int cylGrow = 0;
    int unclaimedTris = 0;
    int consumedWalls = 0;
    std::string sig;
};

static Outcome run(Synth& synth, SegmentWork& work, int threads, double q) {
    MeshView mv = synth.view(q);
    SegmentParams p;
    p.requestedThreads = threads;
    p.maxRefineIters = 16;
    DerivedTols tol;
    refuteNgonByLaw(mv, p, tol, work);
    Outcome o;
    std::vector<std::pair<int, int>> regs;
    for (const Region& r : work.accepted) {
        const int n0 = r.tris.empty() ? -1 : *std::min_element(r.tris.begin(), r.tris.end());
        regs.push_back({n0, static_cast<int>(r.origin) * 1000 + r.nSides});
        if (r.origin == Origin::NgonWall && r.nSides == 6) ++o.ngon6;
        if (r.origin == Origin::NgonWall && r.nSides == 8) ++o.ngon8;
        if (r.origin == Origin::NgonWall && r.nSides == 24) ++o.ngon24;
        if (r.origin == Origin::NgonWall && r.nSides == 64) ++o.ngon64;
        if (r.origin == Origin::NgonWall && r.nSides == 126) ++o.ngon126;
        if (r.origin == Origin::CylGrow) ++o.cylGrow;
    }
    std::sort(regs.begin(), regs.end());
    for (const Provisional& pr : work.provisionals) {
        if (pr.claim == ProvClaim::Unclaimed) {
            o.unclaimedTris += static_cast<int>(pr.tris.size());
            std::vector<int> ts = pr.tris;
            std::sort(ts.begin(), ts.end());
            o.sig += "U";
            for (int t : ts) o.sig += "," + std::to_string(t);
            o.sig += ";";
        } else if (pr.claim == ProvClaim::ConsumedCylinder) {
            ++o.consumedWalls;
        }
    }
    for (auto& rg : regs)
        o.sig += "R" + std::to_string(rg.first) + ":" + std::to_string(rg.second) + ";";
    return o;
}

static void testHexRefutedByFineRound() {
    Synth s;
    const auto hex = s.addTube(8.66, 6, 0.0, 5.0, 0.0, 0.0);
    const auto fine = s.addTube(1.0, 48, 0.0, 5.0, 20.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
    w.accepted.push_back(makeCyl(fine, Origin::CylGrow, true));
    addClearedWalls(w, hex);
    addClearedWalls(w, fine);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon6 == 0, "case1 hex refuted");
    expect(o.cylGrow == 1, "case1 fine round kept");
    expect(o.unclaimedTris == hex.nTri, "case1 hex walls unclaimed");
}

static void testOctagonNoSmallerCurve() {
    Synth s;
    const auto oct = s.addTube(5.0, 8, 0.0, 4.0, 0.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(oct, Origin::NgonWall, true));
    addClearedWalls(w, oct);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon8 == 1, "case2 8-gon kept");
    expect(o.unclaimedTris == 0, "case2 walls stay claimed");
}

static void testHexRefutedBy64() {
    Synth s;
    const auto fine = s.addTube(5.0, 64, 0.0, 4.0, 0.0, 0.0);
    const auto hex = s.addTube(8.66, 6, 0.0, 4.0, 30.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(fine, Origin::NgonWall, true));
    w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
    addClearedWalls(w, fine);
    addClearedWalls(w, hex);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon64 == 1, "case3 64-gon kept");
    expect(o.ngon6 == 0, "case3 hex refuted");
    expect(o.unclaimedTris == hex.nTri, "case3 hex walls unclaimed");
}

static void testHexKeptWhenRefsAreLarger() {
    Synth s;
    const auto hex = s.addTube(2.0, 6, 0.0, 3.0, 0.0, 0.0);
    const auto big = s.addTube(10.0, 64, 0.0, 3.0, 40.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
    w.accepted.push_back(makeCyl(big, Origin::CylGrow, true));
    addClearedWalls(w, hex);
    addClearedWalls(w, big);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon6 == 1, "case4 hex kept (no smaller reference)");
    expect(o.cylGrow == 1, "case4 larger round kept");
    expect(o.unclaimedTris == 0, "case4 nothing released");
}

static void testBadRefitIgnored() {
    Synth s;
    const auto hex = s.addTube(8.66, 6, 0.0, 5.0, 0.0, 0.0);
    const auto bad = s.addTube(1.0, 16, 0.0, 5.0, 25.0, 0.0);
    s.bump(bad, 0, 1.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
    w.accepted.push_back(makeCyl(bad, Origin::CylGrow, true));
    addClearedWalls(w, hex);
    addClearedWalls(w, bad);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon6 == 1, "case5 hex kept (reference fails tau)");
    expect(o.cylGrow == 1, "case5 bad reference still present");
    expect(o.unclaimedTris == 0, "case5 nothing released");
}

static void testThreadsAgree() {
    auto once = [](int threads) {
        Synth s;
        const auto hex = s.addTube(8.66, 6, 0.0, 5.0, 0.0, 0.0);
        const auto fine = s.addTube(1.0, 48, 0.0, 5.0, 20.0, 0.0);
        SegmentWork w;
        w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
        w.accepted.push_back(makeCyl(fine, Origin::CylGrow, true));
        addClearedWalls(w, hex);
        addClearedWalls(w, fine);
        return run(s, w, threads, 1e-4);
    };
    const Outcome a = once(1);
    const Outcome b = once(8);
    expect(a.sig == b.sig, "case6 threads 1 and 8 agree");
    expect(a.ngon6 == 0 && a.unclaimedTris == 12, "case6 still refutes the hex");
}

static std::vector<double> cumulateDeg(const std::vector<double>& stepsDeg) {
    std::vector<double> ang;
    ang.push_back(0.0);
    double a = 0.0;
    for (double s : stepsDeg) {
        a += s * kPi / 180.0;
        ang.push_back(a);
    }
    return ang;
}

// Case 7. Interleaved quarter-turn: one rail is 6×15°, the other is
// 5.625, 5.625, then 7×11.25°. The union of those azimuths is a 11.25° step
// with k=13 (bound under 15°). The equal rail's bound is 18°, so the 24-gon
// at R2.5 stays. The old union reading refutes this case.
static void testInterleavedStripKeeps24() {
    Synth s;
    const std::vector<double> railA = cumulateDeg(std::vector<double>(6, 15.0));
    std::vector<double> bSteps{5.625, 5.625};
    bSteps.insert(bSteps.end(), 7, 11.25);
    const std::vector<double> railB = cumulateDeg(bSteps);
    const auto strip = s.addStrip(1.0, railA, railB, 0.0, 6.0, 0.0, 0.0);
    const auto gon = s.addTube(2.5, 24, 0.0, 4.0, 30.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(strip, Origin::CylGrow, false));
    w.accepted.push_back(makeCyl(gon, Origin::NgonWall, true));
    addClearedWalls(w, gon);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon24 == 1, "case7 24-gon kept (interleaved rails)");
    expect(o.cylGrow == 1, "case7 strip kept");
    expect(o.unclaimedTris == 0, "case7 nothing released");
}

// Case 8. Both rails graded 15, 15, 15, 13.5, 10.935, 9.8415, 5.3618, 5.3618°.
// No equal-step sample, so the octagon at R7 stays.
static void testGradedRailKeepsOctagon() {
    const std::vector<double> steps{15.0, 15.0, 15.0, 13.5, 10.935, 9.8415, 5.3618, 5.3618};
    const std::vector<double> ang = cumulateDeg(steps);
    Synth s;
    const auto strip = s.addStrip(1.0, ang, ang, 0.0, 6.0, 0.0, 0.0);
    const auto oct = s.addTube(7.0, 8, 0.0, 4.0, 30.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(strip, Origin::CylGrow, false));
    w.accepted.push_back(makeCyl(oct, Origin::NgonWall, true));
    addClearedWalls(w, oct);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon8 == 1, "case8 octagon kept (no law sample)");
    expect(o.cylGrow == 1, "case8 graded strip kept");
    expect(o.unclaimedTris == 0, "case8 nothing released");
}

// Case 9. R1.5 sample is 3×8.3336° (U 12.5°). A witness quarter at R2 steps
// 6×15°, which exceeds that bound, so the 24-gon at R2.5 stays. Without the
// witness the same sample refutes it.
static Outcome guardCase(bool withWitness) {
    Synth s;
    const std::vector<double> sample = cumulateDeg({8.3336, 8.3336, 8.3336});
    const std::vector<double> wit = cumulateDeg(std::vector<double>(6, 15.0));
    const auto strip = s.addStrip(1.5, sample, sample, 0.0, 6.0, 0.0, 0.0);
    Synth::Tube witness;
    if (withWitness) witness = s.addStrip(2.0, wit, wit, 0.0, 6.0, 20.0, 0.0);
    const auto gon = s.addTube(2.5, 24, 0.0, 4.0, 40.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(strip, Origin::CylGrow, false));
    if (withWitness) w.accepted.push_back(makeCyl(witness, Origin::CylGrow, false));
    w.accepted.push_back(makeCyl(gon, Origin::NgonWall, true));
    addClearedWalls(w, gon);
    return run(s, w, 1, 1e-4);
}

static void testGuardRejectsSample() {
    const Outcome kept = guardCase(true);
    expect(kept.ngon24 == 1, "case9 24-gon kept when the R2 witness contradicts U");
    expect(kept.cylGrow == 2, "case9 both witness strips kept");
    expect(kept.unclaimedTris == 0, "case9 nothing released");
    const Outcome dropped = guardCase(false);
    expect(dropped.ngon24 == 0, "case9 24-gon refuted without the R2 witness");
    expect(dropped.cylGrow == 1, "case9 sample strip kept");
    expect(dropped.unclaimedTris == 48, "case9 walls unclaimed");
}

// Case 10. One rail is 32×2.5944° then 2.1928°; the other is 12.812° down by
// ×0.9 per step while the step stays at least 0.559°. Neither is equal-step,
// so the 126-gon at R16 stays.
static void testRemainderAndGradedKeep126() {
    std::vector<double> rem(32, 2.5944);
    rem.push_back(2.1928);
    std::vector<double> graded;
    for (double s = 12.812; s + 1e-12 >= 0.559; s *= 0.9) graded.push_back(s);
    const std::vector<double> angA = cumulateDeg(rem);
    const std::vector<double> angB = cumulateDeg(graded);
    Synth s;
    const auto strip = s.addStrip(2.0, angA, angB, 0.0, 6.0, 0.0, 0.0);
    const auto gon = s.addTube(16.0, 126, 0.0, 4.0, 80.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(strip, Origin::CylGrow, false));
    w.accepted.push_back(makeCyl(gon, Origin::NgonWall, true));
    addClearedWalls(w, gon);
    const Outcome o = run(s, w, 1, 1e-4);
    expect(o.ngon126 == 1, "case10 126-gon kept");
    expect(o.cylGrow == 1, "case10 remainder strip kept");
    expect(o.unclaimedTris == 0, "case10 nothing released");
}

static void testTauZeroChangesNothing() {
    Synth s;
    const auto hex = s.addTube(8.66, 6, 0.0, 5.0, 0.0, 0.0);
    const auto fine = s.addTube(1.0, 48, 0.0, 5.0, 20.0, 0.0);
    SegmentWork w;
    w.accepted.push_back(makeCyl(hex, Origin::NgonWall, true));
    w.accepted.push_back(makeCyl(fine, Origin::CylGrow, true));
    addClearedWalls(w, hex);
    addClearedWalls(w, fine);
    const Outcome o = run(s, w, 1, 0.0);
    expect(o.ngon6 == 1 && o.cylGrow == 1 && o.unclaimedTris == 0, "tau 0 changes nothing");
}

int main() {
    setenv("STL2STEP_DIAG_TESSLAW", "1", 1);
    testHexRefutedByFineRound();
    testOctagonNoSmallerCurve();
    testHexRefutedBy64();
    testHexKeptWhenRefsAreLarger();
    testBadRefitIgnored();
    testThreadsAgree();
    testInterleavedStripKeeps24();
    testGradedRailKeepsOctagon();
    testGuardRejectsSample();
    testRemainderAndGradedKeep126();
    testTauZeroChangesNothing();
    if (gFail) {
        std::fprintf(stderr, "tesslaw_unit %d FAILED\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "tesslaw_unit ok\n");
    return 0;
}
