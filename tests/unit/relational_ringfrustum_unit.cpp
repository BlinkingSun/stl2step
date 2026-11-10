// D-gid-8.4. The late ring-frustum closure: between two certified
// coaxial rims the wall closes as the cylinder or cone it certifies as. The
// profile line (profileLine), then claimRingFrustums on a synthetic chamfered
// through-hole: bore + two 45 deg chamfers whose rims carry different vertex
// counts, a planar annulus, a bent wall, a certified bore that donates its
// frame, an uncertified member robbed, a chart held whole by one region and a
// frustum with no donor.
// Synthetic geometry only.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_ringwall.hpp"

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

#include <algorithm>
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

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kDeg = kPi / 180.0;

// ---- profileLine (pure) -----------------------------------------------------------

static void testProfileLine() {
    // A 45 deg frustum profile: rho = 1 + (z - 1) for z in [1, 1.3].
    std::vector<double> rho{1.0, 1.0, 1.3, 1.3}, z{1.0, 1.0, 1.3, 1.3};
    ProfileLine L = profileLine(rho, z);
    expect(L.ok && L.maxDev < 1e-12, "line: frustum profile is a line");
    expect(std::abs(L.rhoAt(1.0) - 1.0) < 1e-12 && std::abs(L.rhoAt(1.3) - 1.3) < 1e-12,
           "line: radius at both rims");
    // A vertical line (a cylinder) is a line too; its radius is constant.
    ProfileLine C = profileLine({2.0, 2.0, 2.0}, {0.0, 0.5, 1.0});
    expect(C.ok && std::abs(C.rhoAt(0.25) - 2.0) < 1e-12, "line: cylinder profile");
    // Constant z: an annulus, not a frustum.
    expect(!profileLine({1.0, 1.3}, {0.0, 0.0}).ok, "line: constant z refused");
    expect(!profileLine({1.0}, {0.0}).ok, "line: one point refused");
    expect(!profileLine({1.0, 2.0}, {0.0}).ok, "line: size mismatch refused");
    // A middle point 1e-3 off the chord shows as a perpendicular deviation.
    ProfileLine B = profileLine({1.0, 1.15 + 1e-3, 1.3}, {1.0, 1.15, 1.3});
    expect(B.ok && B.maxDev > 1e-4 && B.maxDev < 1e-3, "line: bent profile deviates");
    expect(std::abs(profileLineDist(L, 1.15, 1.15)) < 1e-12, "line: on-line point distance 0");
}

// ---- synthetic mesh -----------------------------------------------------------------

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
        mv.diag = 100.0;
        mv.weldTol = 0.0;
        mv.sewTol = 0.0;
        mv.quantFloor = q;
        return mv;
    }
};

// A ring of n vertices at radius R, height z, azimuth offset th0 (axis +Z).
static std::vector<int> ring(Synth& s, int n, double R, double z, double th0) {
    std::vector<int> v;
    for (int i = 0; i < n; ++i) {
        const double th = th0 + 2.0 * kPi * i / n;
        v.push_back(s.add(gp_XYZ(R * std::cos(th), R * std::sin(th), z)));
    }
    return v;
}

static double azOf(const Synth& s, int v) {
    const gp_XYZ& p = s.pts[(std::size_t)v];
    double a = std::atan2(p.Y(), p.X());
    if (a < 0.0) a += 2.0 * kPi;
    return a;
}

// The strip between two closed rings, by azimuth merge (ring counts may differ).
static std::vector<int> strip(Synth& s, const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> out;
    const int na = (int)a.size(), nb = (int)b.size();
    int i = 0, j = 0;
    auto aNext = [&](int k) { return azOf(s, a[(std::size_t)(k % na)]) + (k >= na ? 2 * kPi : 0); };
    auto bNext = [&](int k) { return azOf(s, b[(std::size_t)(k % nb)]) + (k >= nb ? 2 * kPi : 0); };
    while (i < na || j < nb) {
        const bool stepA = j >= nb || (i < na && aNext(i + 1) <= bNext(j + 1));
        if (stepA) {
            out.push_back(s.addTri(a[(std::size_t)(i % na)], a[(std::size_t)((i + 1) % na)],
                                   b[(std::size_t)(j % nb)]));
            ++i;
        } else {
            out.push_back(s.addTri(a[(std::size_t)(i % na)], b[(std::size_t)((j + 1) % nb)],
                                   b[(std::size_t)(j % nb)]));
            ++j;
        }
    }
    return out;
}

static void addChart(SegmentWork& w, const std::vector<int>& tris, int chart) {
    for (int t : tris) {
        if ((int)w.triChart.size() <= t) w.triChart.resize((std::size_t)t + 1, -1);
        w.triChart[(std::size_t)t] = chart;
    }
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.chartId = chart;
    pr.tris = tris;
    std::sort(pr.tris.begin(), pr.tris.end());
    w.provisionals.push_back(pr);
}

static void params(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kDeg;
    tol.thetaCylLo = p.thetaCylLoDeg * kDeg;
}

static void runStage(const Synth& s, SegmentWork& w) {
    const MeshView mv = s.view();
    SegmentParams p;
    DerivedTols tol;
    params(p, tol);
    claimRingFrustums(mv, p, tol, w);
}

// Bore R = 1 between z = +-1 on 40-vertex rims; 45 deg chamfers to R = 1.3 on
// 24-vertex rims at z = +-1.3. Charts: chamfer 0, bore 1, chamfer 2.
constexpr double kRb = 1.0, kC = 0.3, kH = 1.0;
constexpr int kNb = 40, kNa = 24;

struct Hole {
    Synth s;
    SegmentWork w;
    std::vector<int> top, bore, bottom;
};

static void buildHole(Hole& h, double bendTop = 0.0) {
    Synth& s = h.s;
    const std::vector<int> A = ring(s, kNa, kRb + kC, kH + kC, 0.0);
    const std::vector<int> B = ring(s, kNb, kRb, kH, 0.0);
    const std::vector<int> C = ring(s, kNb, kRb, -kH, 0.0);
    const std::vector<int> D = ring(s, kNa, kRb + kC, -kH - kC, 0.0);
    if (bendTop != 0.0) {
        // A middle row on the top chamfer, pushed off its line: a round, not a cone.
        const std::vector<int> M =
            ring(s, kNa, kRb + 0.5 * kC + bendTop, kH + 0.5 * kC, 0.5 * kPi / kNa);
        h.top = strip(s, A, M);
        const std::vector<int> t2 = strip(s, M, B);
        h.top.insert(h.top.end(), t2.begin(), t2.end());
    } else {
        h.top = strip(s, A, B);
    }
    h.bore = strip(s, B, C);
    h.bottom = strip(s, C, D);
    s.finish();
}

static void chartsOf(Hole& h) {
    addChart(h.w, h.top, 0);
    addChart(h.w, h.bore, 1);
    addChart(h.w, h.bottom, 2);
}

static int count(const SegmentWork& w, SurfType t) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.type == t) ++n;
    return n;
}

static int cylIndex(const SegmentWork& w) {
    for (int i = 0; i < (int)w.accepted.size(); ++i)
        if (w.accepted[(std::size_t)i].type == SurfType::Cylinder) return i;
    return -1;
}

static bool conesOnBore(const SegmentWork& w, int ci, double tau, bool exactFrame) {
    if (ci < 0) return false;
    const Region& cy = w.accepted[(std::size_t)ci];
    int n = 0;
    for (const Region& r : w.accepted) {
        if (r.type != SurfType::Cone) continue;
        ++n;
        if (r.origin != Origin::ChamferCone || !r.closed360) return false;
        if (!(std::abs(r.vMin - kRb) <= tau) || !(std::abs(r.vMax - (kRb + kC)) <= tau)) return false;
        if (!(std::abs(std::abs(r.maxVertexSnap) - kC) <= tau)) return false;
        if (!(r.maxVertexDev <= tau)) return false;
        if (r.filletNbrA != ci) return false;
        if (std::abs(r.ax.Direction().Dot(cy.ax.Direction())) < 1.0 - 1e-12) return false;
        if (exactFrame && r.ax.XDirection().Dot(cy.ax.XDirection()) < 1.0 - 1e-12) return false;
    }
    return n == 2;
}

static bool allConsumed(const SegmentWork& w) {
    for (const Provisional& p : w.provisionals)
        if (p.claim == ProvClaim::Unclaimed && !p.tris.empty()) return false;
    return true;
}

// Free hole: the bore closes on the ring certificate, both chamfers close as cones
// on the bore's frame.
static void caseHoleAccept() {
    Hole h;
    buildHole(h);
    chartsOf(h);
    const double tau = 2 * h.s.q;
    runStage(h.s, h.w);
    const int ci = cylIndex(h.w);
    expect(h.w.accepted.size() == 3 && count(h.w, SurfType::Cylinder) == 1 &&
               count(h.w, SurfType::Cone) == 2,
           "hole: one cylinder, two cones");
    expect(ci >= 0 && h.w.accepted[(std::size_t)ci].ringCertified &&
               h.w.accepted[(std::size_t)ci].origin == Origin::CylGrow &&
               h.w.accepted[(std::size_t)ci].closed360 &&
               std::abs(h.w.accepted[(std::size_t)ci].radius - kRb) <= tau,
           "hole: bore is a ring-certified closed cylinder");
    expect(conesOnBore(h.w, ci, tau, true), "hole: chamfers are cones on the bore frame");
    expect(allConsumed(h.w), "hole: no provisional left unclaimed");
}

static bool sameXYZ(const gp_XYZ& a, const gp_XYZ& b) {
    return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z();
}

// Determinism: the same input closes to the same fields.
static void caseHoleDeterministic() {
    Hole a, b;
    buildHole(a);
    chartsOf(a);
    buildHole(b);
    chartsOf(b);
    runStage(a.s, a.w);
    runStage(b.s, b.w);
    bool same = a.w.accepted.size() == b.w.accepted.size();
    for (std::size_t i = 0; same && i < a.w.accepted.size(); ++i) {
        const Region& x = a.w.accepted[i];
        const Region& y = b.w.accepted[i];
        same = x.type == y.type && x.tris == y.tris && x.radius == y.radius &&
               x.vMin == y.vMin && x.vMax == y.vMax && x.dVolPredicted == y.dVolPredicted &&
               sameXYZ(x.ax.Location().XYZ(), y.ax.Location().XYZ()) &&
               sameXYZ(x.ax.Direction().XYZ(), y.ax.Direction().XYZ()) &&
               sameXYZ(x.ax.XDirection().XYZ(), y.ax.XDirection().XYZ());
    }
    expect(same, "hole: two runs identical");
}

// A wall between two coplanar rims is a plane: left to A3.
static void casePlanarAnnulus() {
    Synth s;
    SegmentWork w;
    const std::vector<int> A = ring(s, kNa, 1.3, 0.0, 0.0);
    const std::vector<int> B = ring(s, kNb, 1.0, 0.0, 0.0);
    const std::vector<int> t = strip(s, A, B);
    s.finish();
    addChart(w, t, 0);
    runStage(s, w);
    expect(w.accepted.empty() && w.provisionals[0].claim == ProvClaim::Unclaimed,
           "annulus: plane refused, provisional untouched");
}

// A chamfer whose middle row stands 50 tau off the rims' line is no cone (a round
// is not invented either); the bore and the other chamfer still close.
static void caseBentWall() {
    Hole h;
    const double tau = 2 * h.s.q;
    buildHole(h, 50 * tau);
    chartsOf(h);
    runStage(h.s, h.w);
    bool topFree = h.w.provisionals[0].claim == ProvClaim::Unclaimed;
    for (const Region& r : h.w.accepted)
        for (int t : h.top)
            if (std::find(r.tris.begin(), r.tris.end(), t) != r.tris.end()) topFree = false;
    expect(topFree, "bent: chamfer with an off-line row refused");
    expect(count(h.w, SurfType::Cylinder) == 1 && count(h.w, SurfType::Cone) == 1,
           "bent: bore and the straight chamfer close");
}

// A certified cylinder already on the bore is never robbed; it donates its frame.
static void caseCertifiedBoreDonates() {
    Hole h;
    buildHole(h);
    const double tau = 2 * h.s.q;
    addChart(h.w, h.top, 0);
    for (int t : h.bore) {
        if ((int)h.w.triChart.size() <= t) h.w.triChart.resize((std::size_t)t + 1, -1);
        h.w.triChart[(std::size_t)t] = 1;
    }
    Region cy;
    cy.type = SurfType::Cylinder;
    cy.origin = Origin::CylGrow;
    cy.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(0, 1, 0));
    cy.radius = kRb;
    cy.closed360 = true;
    cy.tris = h.bore;
    cy.maxVertexDev = 0.0;
    h.w.accepted.push_back(cy);
    addChart(h.w, h.bottom, 2);
    runStage(h.s, h.w);
    const Region& b = h.w.accepted[0];
    expect(b.tris == h.bore && !b.ringCertified && sameXYZ(b.ax.XDirection().XYZ(), gp_XYZ(0, 1, 0)),
           "certified: bore region untouched");
    expect(conesOnBore(h.w, 0, tau, true), "certified: cones take the bore frame");
}

// An uncertified CylGrow member inside the free bore chart is robbed (D-fc-4.3).
static void caseUncertifiedMemberRobbed() {
    Hole h;
    buildHole(h);
    const double tau = 2 * h.s.q;
    chartsOf(h);
    Region chip;
    chip.type = SurfType::Cylinder;
    chip.origin = Origin::CylGrow;
    chip.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    chip.radius = kRb + 10 * tau;
    chip.tris = {h.bore[0], h.bore[1], h.bore[2], h.bore[3]};
    chip.maxVertexDev = 10 * tau;
    h.w.accepted.push_back(chip);
    runStage(h.s, h.w);
    const int ci = cylIndex(h.w);
    expect(h.w.accepted.size() == 3 && ci >= 0 &&
               h.w.accepted[(std::size_t)ci].tris.size() == h.bore.size() &&
               h.w.accepted[(std::size_t)ci].ringCertified,
           "robbed: the chip dissolves into the closed bore");
    expect(conesOnBore(h.w, ci, tau, true), "robbed: cones name the closed bore");
}

// A chart held whole by one uncertified region is not entered (RQ-H6). That
// region is no certified hole wall, closed or not, so the chamfers have no donor
// and stay as they are (D-gid-2.5 / D-gid-3.6).
static void caseHeldWholeNotEntered() {
    for (int closed = 0; closed < 2; ++closed) {
        Hole h;
        buildHole(h);
        const double tau = 2 * h.s.q;
        addChart(h.w, h.top, 0);
        for (int t : h.bore) {
            if ((int)h.w.triChart.size() <= t) h.w.triChart.resize((std::size_t)t + 1, -1);
            h.w.triChart[(std::size_t)t] = 1;
        }
        Region cy;
        cy.type = SurfType::Cylinder;
        cy.origin = Origin::CylGrow;
        cy.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        cy.radius = kRb + 10 * tau;
        cy.closed360 = closed == 1;
        cy.tris = h.bore;
        cy.maxVertexDev = 10 * tau;
        h.w.accepted.push_back(cy);
        addChart(h.w, h.bottom, 2);
        runStage(h.s, h.w);
        const Region& b = h.w.accepted[0];
        expect(h.w.accepted.size() == 1 && b.tris == h.bore && !b.ringCertified &&
                   b.maxVertexDev == 10 * tau,
               closed ? "held: closed uncertified bore untouched, no cone"
                      : "held: open uncertified bore untouched, no cone");
        expect(h.w.provisionals[0].claim == ProvClaim::Unclaimed &&
                   h.w.provisionals[1].claim == ProvClaim::Unclaimed,
               closed ? "held: chamfers stay free (closed)" : "held: chamfers stay free (open)");
    }
}

// A frustum between two rims with no cylinder beside it (both rims meet other
// walls) has no donor: refused, provisional untouched.
static void caseLoneFrustumNoDonor() {
    Synth s;
    SegmentWork w;
    const std::vector<int> A = ring(s, kNa, kRb + kC, kH + kC, 0.0);
    const std::vector<int> B = ring(s, kNb, kRb, kH, 0.0);
    const std::vector<int> t = strip(s, A, B);
    s.finish();
    addChart(w, t, 0);
    runStage(s, w);
    expect(w.accepted.empty() && w.provisionals[0].claim == ProvClaim::Unclaimed,
           "lone: frustum without a donor refused");
}

// No quantum, no certificate: nothing is admitted.
static void caseNoQuantum() {
    Hole h;
    h.s.q = 0.0;
    buildHole(h);
    chartsOf(h);
    runStage(h.s, h.w);
    expect(h.w.accepted.empty(), "q0: nothing admitted");
}

int main() {
    testProfileLine();
    caseHoleAccept();
    caseHoleDeterministic();
    casePlanarAnnulus();
    caseBentWall();
    caseCertifiedBoreDonates();
    caseUncertifiedMemberRobbed();
    caseHeldWholeNotEntered();
    caseLoneFrustumNoDonor();
    caseNoQuantum();
    if (gFail) std::fprintf(stderr, "relational_ringfrustum_unit: %d FAIL\n", gFail);
    return gFail ? 1 : 0;
}
