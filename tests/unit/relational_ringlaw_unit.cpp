// D-gid-8.4 unit cases: the tessellation-law decision for a closed ring band
// that folds at or above thetaCylLo (tessLawSamples / tessLawBandStep /
// tessLawDecide). Synthetic meshes only.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"
#include "refit_relational.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

using stl2step::refit::MeshView;
using stl2step::refit::Origin;
using stl2step::refit::Region;
using stl2step::refit::SegmentParams;
using stl2step::refit::SegmentWork;
using stl2step::refit::SurfType;
using stl2step::refit::TessLawDecision;
using stl2step::refit::TessLawSample;
using stl2step::refit::TessLawVerdict;
using stl2step::refit::tessLawBandStep;
using stl2step::refit::tessLawDecide;
using stl2step::refit::tessLawSamples;

static const double kPi = 3.14159265358979323846264338327950288;
static const double kTwoPi = 2.0 * kPi;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static bool near(double a, double b, double tol) { return std::abs(a - b) <= tol; }

// One component, local id == global id. A ring is a closed wall between a
// bottom loop of nA and a top loop of nB vertices, both equal-step from
// azimuth 0, every vertex exactly on the cylinder.
struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris;
    std::vector<int> compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::uint8_t> triDirs;

    struct Ring {
        int firstTri = 0;
        int nTri = 0;
        double R = 0;
        double z0 = 0, z1 = 0;
        double cx = 0;
    };

    Ring addRing(double R, int nA, int nB, double z0, double z1, double cx) {
        Ring rg;
        rg.firstTri = static_cast<int>(triFlat.size() / 3);
        rg.R = R;
        rg.z0 = z0;
        rg.z1 = z1;
        rg.cx = cx;
        const int baseA = static_cast<int>(pts.size());
        for (int i = 0; i < nA; ++i) {
            const double a = kTwoPi * i / nA;
            pts.push_back(gp_XYZ(cx + R * std::cos(a), R * std::sin(a), z0));
        }
        const int baseB = static_cast<int>(pts.size());
        for (int j = 0; j < nB; ++j) {
            const double a = kTwoPi * j / nB;
            pts.push_back(gp_XYZ(cx + R * std::cos(a), R * std::sin(a), z1));
        }
        auto tri = [&](int a, int b, int c) {
            triFlat.push_back(a);
            triFlat.push_back(b);
            triFlat.push_back(c);
        };
        int i = 0, j = 0;
        while (i < nA || j < nB) {
            const double ai1 = kTwoPi * (i + 1) / nA;
            const double bj1 = kTwoPi * (j + 1) / nB;
            if (j == nB || (i < nA && ai1 <= bj1)) {
                tri(baseA + i % nA, baseA + (i + 1) % nA, baseB + j % nB);
                ++i;
            } else {
                tri(baseA + i % nA, baseB + (j + 1) % nB, baseB + j % nB);
                ++j;
            }
        }
        rg.nTri = static_cast<int>(triFlat.size() / 3) - rg.firstTri;
        return rg;
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

static std::vector<int> trisOf(const Synth::Ring& rg) {
    std::vector<int> t;
    for (int i = 0; i < rg.nTri; ++i) t.push_back(rg.firstTri + i);
    return t;
}

static Region makeCyl(const Synth::Ring& rg, Origin origin, int nSides) {
    Region r;
    r.type = SurfType::Cylinder;
    r.origin = origin;
    r.radius = rg.R;
    r.nSides = nSides;
    r.closed360 = true;
    r.ax = gp_Ax3(gp_Pnt(rg.cx, 0.0, 0.5 * (rg.z0 + rg.z1)), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.vMin = rg.z0;
    r.vMax = rg.z1;
    r.tris = trisOf(rg);
    return r;
}

static TessLawSample sample(int minTri, double r, int k, double maxStep) {
    TessLawSample s;
    s.minTri = minTri;
    s.r = r;
    s.k = k;
    s.maxStep = maxStep;
    s.U = maxStep * k / (k - 1);
    return s;
}

// Decide: pure function cases.
static void testDecide() {
    const double tau = 2e-4;
    const std::vector<TessLawSample> none;
    expect(tessLawDecide(none, 5.0, 0.1, tau).verdict == TessLawVerdict::NoRef,
           "decide no samples -> no-ref");
    // Only a sample of a larger circle: it bounds nothing at R (D-fc-2.1 r <= R + tau).
    const std::vector<TessLawSample> larger{sample(7, 6.0, 48, kTwoPi / 48)};
    expect(tessLawDecide(larger, 5.0, kTwoPi / 6, tau).verdict == TessLawVerdict::NoRef,
           "decide larger-only sample -> no-ref");
    // A circle of R + tau is still a sample (boundary inclusive).
    const std::vector<TessLawSample> edge{sample(7, 5.0 + tau, 48, kTwoPi / 48)};
    expect(tessLawDecide(edge, 5.0, kTwoPi / 6, tau).verdict == TessLawVerdict::Refute,
           "decide r = R + tau sample counts");
    // Fine round R1 in 48 steps: a hexagon at R7 is refuted, a 64-step ring kept.
    const std::vector<TessLawSample> fine{sample(3, 1.0, 48, kTwoPi / 48)};
    const TessLawDecision hex = tessLawDecide(fine, 7.0, kTwoPi / 6, tau);
    expect(hex.verdict == TessLawVerdict::Refute && hex.refRid == 3 && hex.refK == 48,
           "decide hexagon refuted by a fine round");
    const TessLawDecision c64 = tessLawDecide(fine, 7.0, kTwoPi / 64, tau);
    expect(c64.verdict == TessLawVerdict::Keep && near(c64.bound, kTwoPi / 47, 1e-12),
           "decide 64-step ring kept");
    // step == U refutes (the D-fc-2.1 >=).
    const double U = kTwoPi / 47;
    expect(tessLawDecide(fine, 7.0, U, tau).verdict == TessLawVerdict::Refute,
           "decide step == U refutes");
    // A hit beats a miss; among hits the smallest U, then the smallest region id.
    const std::vector<TessLawSample> mixed{sample(9, 1.0, 4, kPi / 2), sample(5, 1.0, 16, kTwoPi / 16),
                                           sample(2, 1.0, 16, kTwoPi / 16)};
    const TessLawDecision m = tessLawDecide(mixed, 3.0, kTwoPi / 12, tau);
    expect(m.verdict == TessLawVerdict::Refute && m.refRid == 2 && m.refK == 16,
           "decide hit beats miss, tie -> smaller id");
    // A miss reports the tightest bound.
    const TessLawDecision k = tessLawDecide(mixed, 3.0, kTwoPi / 32, tau);
    expect(k.verdict == TessLawVerdict::Keep && k.refRid == 2 && near(k.bound, kTwoPi / 15, 1e-12),
           "decide keep reports tightest bound");
    // tau <= 0 or a non-positive step: no verdict.
    expect(tessLawDecide(fine, 7.0, kTwoPi / 6, 0.0).verdict == TessLawVerdict::NoRef,
           "decide tau 0 -> no-ref");
    expect(tessLawDecide(fine, 7.0, 0.0, tau).verdict == TessLawVerdict::NoRef,
           "decide step 0 -> no-ref");
}

// Band step: the coarsest boundary rail of a closed ring band.
static void testBandStep() {
    const double q = 1e-4, tau = 2 * q;
    Synth s;
    const auto r16 = s.addRing(3.0, 16, 16, 0.0, 10.0, 0.0);
    const auto r1218 = s.addRing(3.0, 12, 18, 0.0, 10.0, 10.0);
    const auto hex = s.addRing(7.0, 6, 6, 0.0, 5.0, 40.0);
    const MeshView mv = s.view(q);
    double step = 0;
    int k = 0;
    const gp_Dir z(0, 0, 1);
    bool ok = tessLawBandStep(mv, trisOf(r16), z, gp_Pnt(0, 0, 0), 3.0, tau, step, k);
    expect(ok && k == 16 && near(step, kTwoPi / 16, 1e-9), "band step 16-ring");
    ok = tessLawBandStep(mv, trisOf(r1218), z, gp_Pnt(10, 0, 0), 3.0, tau, step, k);
    expect(ok && k == 12 && near(step, kTwoPi / 12, 1e-9), "band step 12/18 ring -> coarser rim");
    ok = tessLawBandStep(mv, trisOf(hex), z, gp_Pnt(40, 0, 0), 7.0, tau, step, k);
    expect(ok && k == 6 && near(step, kTwoPi / 6, 1e-9), "band step hexagon");
    expect(!tessLawBandStep(mv, {}, z, gp_Pnt(0, 0, 0), 3.0, tau, step, k), "band step empty");
    expect(!tessLawBandStep(mv, trisOf(r16), z, gp_Pnt(0, 0, 0), 3.0, 0.0, step, k),
           "band step tau 0");
}

// Samples: certified non-N-gon references only, admissible rails only.
static void testSamples() {
    const double q = 1e-4, tau = 2 * q;
    Synth s;
    const auto fine = s.addRing(1.0, 48, 48, 0.0, 5.0, 0.0);
    const auto ngon = s.addRing(5.0, 16, 16, 0.0, 5.0, 20.0);
    const auto bent = s.addRing(2.0, 32, 32, 0.0, 5.0, 40.0);
    // Bump one vertex of `bent` off its circle by 50 tau: its refit fails tau.
    const int v0 = s.triFlat[static_cast<size_t>(3 * bent.firstTri)];
    gp_XYZ& p = s.pts[static_cast<size_t>(v0)];
    const double rx = p.X() - 40.0, ry = p.Y();
    const double m = std::hypot(rx, ry);
    p.SetX(40.0 + rx * (1.0 + 50 * tau / m));
    p.SetY(ry * (1.0 + 50 * tau / m));
    const MeshView mv = s.view(q);
    SegmentParams prm;
    prm.maxRefineIters = 16;
    SegmentWork w;
    w.accepted.push_back(makeCyl(ngon, Origin::NgonWall, 16));
    w.accepted.push_back(makeCyl(fine, Origin::CylGrow, 48));
    w.accepted.push_back(makeCyl(bent, Origin::CylGrow, 32));
    const std::vector<TessLawSample> ss = tessLawSamples(mv, prm, tau, w);
    bool onlyFine = !ss.empty();
    for (const TessLawSample& x : ss)
        if (x.minTri != fine.firstTri || x.k != 48 || !near(x.r, 1.0, tau)) onlyFine = false;
    expect(onlyFine && ss.size() == 2, "samples: two rails of the certified round only");
    expect(!ss.empty() && near(ss[0].U, kTwoPi / 47, 1e-9), "samples: U = maxStep k/(k-1)");
    expect(tessLawSamples(mv, prm, 0.0, w).empty(), "samples: tau 0 -> none");

    // End to end: the hexagon is refuted, a 12/18 ring at R3 kept, no sample -> no-ref.
    double step = 0;
    int k = 0;
    Synth t;
    const auto hex = t.addRing(7.0, 6, 6, 0.0, 5.0, 60.0);
    const auto mixed = t.addRing(3.0, 12, 18, 0.0, 10.0, 80.0);
    const MeshView mt = t.view(q);
    tessLawBandStep(mt, trisOf(hex), gp_Dir(0, 0, 1), gp_Pnt(60, 0, 0), 7.0, tau, step, k);
    expect(tessLawDecide(ss, 7.0, step, tau).verdict == TessLawVerdict::Refute,
           "e2e hexagon refuted");
    tessLawBandStep(mt, trisOf(mixed), gp_Dir(0, 0, 1), gp_Pnt(80, 0, 0), 3.0, tau, step, k);
    const std::vector<TessLawSample> coarse{sample(1, 2.0, 3, kTwoPi / 12)};
    expect(tessLawDecide(coarse, 3.0, step, tau).verdict == TessLawVerdict::Keep,
           "e2e 12/18 ring kept by a coarser-law sample");
    expect(tessLawDecide(ss, 3.0, step, tau).verdict == TessLawVerdict::Refute,
           "e2e 12/18 ring refuted by a 48-step round");
    expect(tessLawDecide({}, 3.0, step, tau).verdict == TessLawVerdict::NoRef,
           "e2e no sample -> no-ref");
}

int main() {
    testDecide();
    testBandStep();
    testSamples();
    std::fprintf(stderr, "%s relational_ringlaw_unit (%d failures)\n", gFail ? "FAIL" : "OK",
                 gFail);
    return gFail ? 1 : 0;
}
