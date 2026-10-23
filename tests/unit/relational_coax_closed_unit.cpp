// Synthetic coaxial-wall cases (D-fc-6.8 RQ-1 closed walls, RQ-3 tangent parallel
// cylinder anchors, D-fc-6.11 fills). No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_coax_closed.hpp"
#include "refit_relational.hpp"
#include "refit_relational_arc.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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
static constexpr double kQ = 1e-6;
static constexpr double kTau = 2 * kQ;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;
    double quantFloor = kQ;

    int add(double x, double y, double z) {
        pts.emplace_back(x, y, z);
        return (int)pts.size() - 1;
    }
    int addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
        return (int)(triFlat.size() / 3) - 1;
    }
    // Two triangles between columns (a0, a1) and (b0, b1) (index 0: z low).
    void quad(int a0, int a1, int b0, int b1, std::vector<int>& out) {
        out.push_back(addTri(a0, b0, a1));
        out.push_back(addTri(a1, b0, b1));
    }
    // The same column with a T-junction vertex m on its low chord (a0, b0).
    void quadT(int a0, int a1, int b0, int b1, int m, std::vector<int>& out) {
        out.push_back(addTri(a0, m, a1));
        out.push_back(addTri(m, b0, b1));
        out.push_back(addTri(m, b1, a1));
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
        mv.quantFloor = quantFloor;
        return mv;
    }
};

static void params(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kPi / 180.0;
    tol.thetaCylLo = 5.0 * kPi / 180.0;  // B1's seed floor (refit.hpp thetaCylLoDeg)
}

// Columns (bottom, top) of an arc of radius R about (cx, cy), azimuths a0..a1
// in n steps, z in [0, h]. `bump` moves column `bumpAt` radially by that much.
static std::vector<std::array<int, 2>> arcCols(Synth& s, double cx, double cy, double R,
                                               double a0, double a1, int n, double h,
                                               int bumpAt = -1, double bump = 0) {
    std::vector<std::array<int, 2>> c;
    for (int i = 0; i <= n; ++i) {
        const double a = a0 + (a1 - a0) * (double)i / (double)n;
        const double r = R + (i == bumpAt ? bump : 0.0);
        const double x = cx + r * std::cos(a), y = cy + r * std::sin(a);
        c.push_back({{s.add(x, y, 0.0), s.add(x, y, h)}});
    }
    return c;
}

static std::array<gp_XYZ, 3> tri3(const gp_XYZ& a, const gp_XYZ& b, const gp_XYZ& c) {
    return std::array<gp_XYZ, 3>{{a, b, c}};
}

static void casePredicates() {
    const gp_XYZ o(0, 0, 0), a(0, 0, 1);
    const double R = 10, step = kPi / 36;
    const gp_XYZ p0(R, 0, 0), p1(R * std::cos(step), R * std::sin(step), 0), p0h(R, 0, 5);
    const double sag = R * (1 - std::cos(step / 2));
    // The chord bound of one column triangle is its rim chord's sagitta.
    const double s = coaxChordBound({tri3(p0, p1, p0h)}, o, a, R);
    expect(std::fabs(s - sag) <= 1e-12, "chord bound = rim sagitta");
    // A T-junction vertex a third along the chord deviates less than s.
    const gp_XYZ m = p0 + (p1 - p0) / 3.0;
    expect(coaxCoverAdmits(tri3(p0, m, p0h), o, a, R, kTau, s, 0, 5), "cover: T-junction column");
    expect(!coaxCoverAdmits(tri3(p0, m, p0h), o, a, R, 0.0, s, 0, 5), "cover: tau <= 0 admits nothing");
    expect(!coaxCoverAdmits(tri3(p0, m, p1), o, a, R, kTau, s, 0, 5), "cover: one level (ear)");
    const gp_XYZ mOff = m * (1 - 2 * sag / R);
    expect(!coaxCoverAdmits(tri3(p0, mOff, p0h), o, a, R, kTau, s, 0, 5), "cover: vertex beyond s");
    expect(!coaxCoverAdmits(tri3(m, m + gp_XYZ(0, 0, 5), mOff), o, a, R, kTau, s, 0, 5),
           "cover: no vertex within tau");
    expect(!coaxCoverAdmits(tri3(p0, m, p0h), o, a, R, kTau, s, 0, 4), "cover: beyond I's extent");
    // Angular extent: the gap [1, 2) excludes 1.5, keeps 0.5 and 2 - slack.
    expect(!coaxInArc(1.5, 1.0, 2.0, 1e-9) && coaxInArc(0.5, 1.0, 2.0, 1e-9) &&
               coaxInArc(2.0, 1.0, 2.0, 1e-9) && coaxInArc(1.5 - 2 * kPi, 1.0, 2.0, 1e-9) == false,
           "arc: largest-gap complement");
    // Fold scope: strictly under the seed floor.
    expect(coaxFoldInScope(4.9 * kPi / 180, 5 * kPi / 180) &&
               !coaxFoldInScope(5 * kPi / 180, 5 * kPi / 180),
           "fold scope: under thetaCylLo only");
    expect(std::fabs(coaxFold(gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0)) - kPi / 2) <= 1e-12, "fold: 90 deg");
    // Closed wall: D2, one component, two rims winding once (sign-free).
    expect(coaxClosedWall(true, 1, {1, -1}) && coaxClosedWall(true, 1, {-1, -1}),
           "closed wall: two rims");
    expect(!coaxClosedWall(false, 1, {1, -1}), "closed wall: D2 open");
    expect(!coaxClosedWall(true, 2, {1, -1}), "closed wall: two components");
    expect(!coaxClosedWall(true, 1, {1, -1, 0}), "closed wall: a window is not committed here");
    expect(!coaxClosedWall(true, 1, {0}), "closed wall: a disk");
    expect(coaxLiftRoutes(false) && !coaxLiftRoutes(true), "route: an annulus preimage is F1's");
    // Fill: two edges shared with I + Cover, the third along the top rim (z = 5),
    // every vertex within s but none within tau.
    {
        const gp_XYZ q0 = p0 + (p1 - p0) * 0.4, q1 = p0h + (p1 + gp_XYZ(0, 0, 5) - p0h) * 0.3,
                     q2 = p0h + (p1 + gp_XYZ(0, 0, 5) - p0h) * 0.6;
        const std::array<gp_XYZ, 3> t = tri3(q0, q2, q1);
        const std::array<std::array<gp_XYZ, 2>, 3> e{{{{q0, q2}}, {{q2, q1}}, {{q1, q0}}}};
        expect(coaxFillAdmits(t, e, {{true, false, true}}, o, a, R, kTau, s, 0, 5),
               "fill: enclosed by I + Cover and the rim");
        expect(!coaxFillAdmits(t, e, {{true, false, false}}, o, a, R, kTau, s, 0, 5),
               "fill: an open side off the rim is an extension");
        expect(!coaxFillAdmits(t, e, {{true, false, true}}, o, a, R, 0.0, s, 0, 5),
               "fill: tau <= 0 admits nothing");
        expect(!coaxFillAdmits(t, e, {{true, false, true}}, o, a, R, kTau, 0.5 * s, 0, 5),
               "fill: a vertex beyond s");
        const std::array<gp_XYZ, 3> ear = tri3(p0h, q1, q2);
        const std::array<std::array<gp_XYZ, 2>, 3> ee{{{{p0h, q1}}, {{q1, q2}}, {{q2, p0h}}}};
        expect(!coaxFillAdmits(ear, ee, {{true, true, false}}, o, a, R, kTau, s, 0, 5),
               "fill: a one-level ear is never a fill");
    }
    // Tangent parallel cylinder: internal (R - r) and external (R + r).
    const double L = 6, tolA = kTau, rA = 5;
    const double rWallIn = 24, rWallOut = 14;
    CoaxTangent t = coaxTangentCylinder(o, a, rWallIn, gp_XYZ(rWallIn - rA, 0, 0), a, rA, tolA, kTau, L);
    expect(t.ok && t.internal && t.tanDist <= kTau, "tangent: internal R - r");
    t = coaxTangentCylinder(o, a, rWallOut, gp_XYZ(0, rWallOut + rA, 0), a, rA, tolA, kTau, L);
    expect(t.ok && !t.internal, "tangent: external R + r");
    // A fillet tangent to the inscribed polygon's facet sits a sagitta inside.
    const double facetHalf = kPi / 48;
    const double sagFacet = rWallIn * (1 - std::cos(facetHalf));
    t = coaxTangentCylinder(o, a, rWallIn, gp_XYZ(rWallIn - rA - sagFacet, 0, 0), a, rA, tolA, kTau, L);
    expect(!t.ok, "tangent: facet-tangent (polygon) refused");
    t = coaxTangentCylinder(o, a, rWallIn, gp_XYZ(rWallIn - rA, 0, 0), gp_XYZ(1e-4, 0, 1), rA, tolA, kTau, L);
    expect(!t.ok, "tangent: tilted axis refused");
    t = coaxTangentCylinder(o, a, rA, gp_XYZ(0, 0, 0), a, rA, tolA, kTau, L);
    expect(!t.ok, "tangent: a coaxial cylinder is no anchor");
    t = coaxTangentCylinder(o, a, rWallIn, gp_XYZ(rWallIn - rA, 0, 0), a, rA, tolA, 0.0, L);
    expect(!t.ok, "tangent: tau <= 0 admits nothing");
    expect(coaxContact(gp_XYZ(rWallIn, 0, 2), o, a, rWallIn, gp_XYZ(rWallIn - rA, 0, 0), a, rA, tolA, kTau) &&
               !coaxContact(gp_XYZ(rWallIn - 1, 0, 2), o, a, rWallIn, gp_XYZ(rWallIn - rA, 0, 0), a, rA,
                            tolA, kTau),
           "contact: on both surfaces");
}

static void addDonor(Synth& s, SegmentWork& w, std::vector<int>& donorTris, double h) {
    const auto dc = arcCols(s, 0, 0, 10, 0, kPi / 3, 12, h);
    for (std::size_t i = 0; i + 1 < dc.size(); ++i)
        s.quad(dc[i][0], dc[i][1], dc[i + 1][0], dc[i + 1][1], donorTris);
    (void)w;
}

static void acceptDonor(SegmentWork& w, const std::vector<int>& donorTris) {
    Region d;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = 10;
    d.tris = donorTris;
    w.accepted.push_back(d);
}

struct WallOut {
    int n = 0;          // accepted R regions
    bool closed = false;
    std::size_t nTri = 0;
    double radius = 0;
};

// D-fc-6.8 RQ-1 end to end: a coaxial family (an R10 donor arc) and a closed
// R wall of n columns about the same axis. `tAt` columns carry a T-junction
// vertex a third along their low chord (the chord triangles split the
// tau-preimage); `bumpAt` moves one column radially by `bump`. `family` false
// drops the donor. `holeAt` triangulates one column through chord vertices on
// both rims so that its middle triangle has no vertex within tau (D-fc-6.11 (a));
// `earAt` adds a cap-plane ear on three consecutive top-rim vertices.
static WallOut closedWall(double R, int n, const std::vector<int>& tAt, int bumpAt, double bump,
                          bool family, double q = kQ, int holeAt = -1, int earAt = -1) {
    Synth s;
    s.quantFloor = q;
    SegmentWork w;
    const double h = 20;
    std::vector<int> donorTris, wallTris;
    if (family) addDonor(s, w, donorTris, h);
    auto bc = arcCols(s, 0, 0, R, 0, 2 * kPi, n, h, bumpAt, bump);
    bc.back() = bc.front();
    for (int i = 0; i + 1 < (int)bc.size(); ++i) {
        if (i == holeAt) {
            const int a0 = bc[(std::size_t)i][0], a1 = bc[(std::size_t)i][1];
            const int b0 = bc[(std::size_t)i + 1][0], b1 = bc[(std::size_t)i + 1][1];
            auto on = [&](int u, int v, double t) {
                const gp_XYZ p = s.pts[(std::size_t)u] + (s.pts[(std::size_t)v] - s.pts[(std::size_t)u]) * t;
                return s.add(p.X(), p.Y(), p.Z());
            };
            const int m0 = on(a0, b0, 0.4), m1 = on(a1, b1, 0.3), m2 = on(a1, b1, 0.6);
            wallTris.push_back(s.addTri(a0, m0, a1));
            wallTris.push_back(s.addTri(m0, m1, a1));
            wallTris.push_back(s.addTri(m0, m2, m1));  // no vertex within tau
            wallTris.push_back(s.addTri(m0, b1, m2));
            wallTris.push_back(s.addTri(m0, b0, b1));
        } else if (std::find(tAt.begin(), tAt.end(), i) != tAt.end()) {
            const gp_XYZ a0 = s.pts[(std::size_t)bc[(std::size_t)i][0]];
            const gp_XYZ b0 = s.pts[(std::size_t)bc[(std::size_t)i + 1][0]];
            const gp_XYZ m = a0 + (b0 - a0) / 3.0;
            const int mi = s.add(m.X(), m.Y(), m.Z());
            s.quadT(bc[(std::size_t)i][0], bc[(std::size_t)i][1], bc[(std::size_t)i + 1][0],
                    bc[(std::size_t)i + 1][1], mi, wallTris);
        } else {
            s.quad(bc[(std::size_t)i][0], bc[(std::size_t)i][1], bc[(std::size_t)i + 1][0],
                   bc[(std::size_t)i + 1][1], wallTris);
        }
    }
    std::vector<int> capTris;
    if (earAt >= 0)
        capTris.push_back(s.addTri(bc[(std::size_t)earAt + 2][1], bc[(std::size_t)earAt + 1][1],
                                   bc[(std::size_t)earAt][1]));
    s.finish();
    if (family) acceptDonor(w, donorTris);
    if (!capTris.empty()) {
        Provisional cap;
        cap.claim = ProvClaim::Unclaimed;
        cap.tris = capTris;
        cap.plane = gp_Ax3(gp_Pnt(0, 0, 20), gp_Dir(0, 0, 1));
        w.provisionals.push_back(cap);
    }
    Provisional band;
    band.claim = ProvClaim::Unclaimed;
    band.tris = wallTris;
    band.plane = gp_Ax3(gp_Pnt(R, 0, 0), gp_Dir(1, 0, 0));
    band.maxVertexDev = 1.0;  // a curved band is no plane anchor
    w.provisionals.push_back(band);
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, w, RelRevolveMode::Line);
    WallOut o;
    for (const Region& r : w.accepted)
        if (r.type == SurfType::Cylinder && std::fabs(r.radius - R) <= kTau) {
            ++o.n;
            o.closed = r.closed360;
            o.nTri = r.tris.size();
            o.radius = r.radius;
        }
    return o;
}

static void caseClosedWall() {
    // 120 columns (3-degree folds, under B1's seed floor), T-junctions on the low
    // rim of four columns: the tau-preimage is four pieces, the cover joins them.
    const std::vector<int> tAt{10, 40, 70, 100};
    WallOut o = closedWall(8, 120, tAt, -1, 0, true);
    expect(o.n == 1 && o.closed && o.nTri == 120 * 2 + 4, "closed wall: ONE closed cylinder over I + cover");
    // The same wall with one column off the family axis by more than s (twice the
    // rim sagitta): no closed wall, and no passing subset of it ships.
    const double sag = 8 * (1 - std::cos(kPi / 120));
    o = closedWall(8, 120, tAt, 45, 2 * sag, true);
    expect(o.n == 0, "closed wall: a column beyond s refuses the band whole");
    // A clean wall (no chord triangle) is the F1 ring path's: deferred unchanged.
    o = closedWall(8, 120, {}, -1, 0, true);
    expect(o.n == 0, "closed wall: an annulus preimage stays closed-deferred");
    // An inscribed N-gon (concyclic corners) with no family: nothing.
    o = closedWall(8, 6, {}, -1, 0, false);
    expect(o.n == 0, "N-gon without a family refused");
    // An octagon on a family, its preimage split by T-junctions: 45-degree folds
    // are B1 / NgonWall / tess-law's (F1 scope), refused.
    o = closedWall(8, 8, {2, 6}, -1, 0, true);
    expect(o.n == 0, "N-gon on a family refused by the fold scope");
    // D-fc-6.11 (a): a column whose middle triangle has no vertex within tau (every
    // vertex a chord vertex within s) is a hole of the closure; it is filled and the
    // wall commits whole. A cap-plane ear on the rim (one level) is never filled.
    o = closedWall(8, 120, tAt, -1, 0, true, kQ, 55);
    expect(o.n == 1 && o.closed && o.nTri == 120 * 2 + 4 + 3, "closed wall: an enclosed chord triangle is filled");
    o = closedWall(8, 120, tAt, -1, 0, true, kQ, 55, 20);
    expect(o.n == 1 && o.closed && o.nTri == 120 * 2 + 4 + 3, "closed wall: a rim ear is never a fill");
    // q <= 0 admits nothing.
    o = closedWall(8, 120, tAt, -1, 0, true, 0.0);
    expect(o.n == 0, "closed wall: q <= 0 admits nothing");
}

// D-fc-6.8 RQ-3 end to end: a partial wall (30 degrees, 10 columns) on the
// family, with no tangent plane. `anchor`: an accepted small cylinder whose axis
// stands at offset (R - r + shift) on the wall's 0-degree column, sharing that
// column's vertices.
static int partialWall(bool anchor, double shift) {
    Synth s;
    SegmentWork w;
    const double h = 20, R = 14, r = 5;
    const int nCols = 10;
    std::vector<int> donorTris, wallTris, anchorTris;
    addDonor(s, w, donorTris, h);
    const auto bc = arcCols(s, 0, 0, R, 0, kPi / 6, nCols, h);
    for (std::size_t i = 0; i + 1 < bc.size(); ++i)
        s.quad(bc[i][0], bc[i][1], bc[i + 1][0], bc[i + 1][1], wallTris);
    const double cx = R - r + shift;
    if (anchor) {
        // A short arc from the wall's 0-degree column (azimuth 0 about its centre)
        // a quarter turn round, internally tangent to the wall there.
        auto ac = arcCols(s, cx, 0, r, 0, -kPi / 2, 6, h);
        ac.front() = bc.front();
        for (std::size_t i = 0; i + 1 < ac.size(); ++i)
            s.quad(ac[i + 1][0], ac[i + 1][1], ac[i][0], ac[i][1], anchorTris);
    }
    s.finish();
    acceptDonor(w, donorTris);
    if (anchor) {
        Region a;
        a.type = SurfType::Cylinder;
        a.origin = Origin::CylGrow;
        a.ax = gp_Ax3(gp_Pnt(cx, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        a.radius = r;
        a.maxVertexDev = 0;
        a.tris = anchorTris;
        w.accepted.push_back(a);
    }
    Provisional band;
    band.claim = ProvClaim::Unclaimed;
    band.tris = wallTris;
    band.plane = gp_Ax3(gp_Pnt(R, 0, 0), gp_Dir(1, 0, 0));
    band.maxVertexDev = 1.0;
    w.provisionals.push_back(band);
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, w, RelRevolveMode::Line);
    int n = 0;
    for (const Region& reg : w.accepted)
        if (reg.type == SurfType::Cylinder && reg.origin == Origin::FilletStrip &&
            std::fabs(reg.radius - R) <= kTau && reg.tris.size() == wallTris.size() &&
            !reg.closed360)
            ++n;
    return n;
}

static void casePartialWall() {
    expect(partialWall(true, 0.0) == 1, "partial wall anchored by a tangent parallel cylinder");
    expect(partialWall(false, 0.0) == 0, "partial wall with no anchor refused (no-anchor)");
    // A cylinder tangent to the 0-degree facet of the inscribed polygon instead of
    // the circle: a sagitta off, no anchor.
    const double sag = 14 * (1 - std::cos(kPi / 6 / (2 * 10)));
    expect(partialWall(true, -sag) == 0, "partial wall: facet-tangent cylinder is no anchor");
}

int main() {
    casePredicates();
    caseClosedWall();
    casePartialWall();
    std::fprintf(stderr, "relational_coax_closed_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS",
                 gFail);
    return gFail ? 1 : 0;
}
