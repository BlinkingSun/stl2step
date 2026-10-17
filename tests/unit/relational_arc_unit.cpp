// Synthetic arc-relation cases (D-fc-5.1 coaxial family, D-fc-5.2 section
// fillets). No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_relational.hpp"
#include "refit_relational_arc.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
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
        mv.quantFloor = kQ;
        return mv;
    }
};

static void params(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kPi / 180.0;
}

// Columns (bottom, top) of an arc of radius R about (cx, cy), azimuths a0..a1
// in n steps, z in [0, h].
static std::vector<std::array<int, 2>> arcCols(Synth& s, double cx, double cy, double R,
                                               double a0, double a1, int n, double h) {
    std::vector<std::array<int, 2>> c;
    for (int i = 0; i <= n; ++i) {
        const double a = a0 + (a1 - a0) * (double)i / (double)n;
        const double x = cx + R * std::cos(a), y = cy + R * std::sin(a);
        c.push_back({{s.add(x, y, 0.0), s.add(x, y, h)}});
    }
    return c;
}

static std::vector<gp_XYZ> arcPts(double R, double a0, double a1, int n, double h, double ox = 0,
                                  double oy = 0) {
    std::vector<gp_XYZ> v;
    for (int i = 0; i <= n; ++i) {
        const double a = a0 + (a1 - a0) * (double)i / (double)n;
        v.emplace_back(ox + R * std::cos(a), oy + R * std::sin(a), 0.0);
        v.emplace_back(ox + R * std::cos(a), oy + R * std::sin(a), h);
    }
    return v;
}

static void casePredicates() {
    // sectionLineResid: collinear in section vs an arc.
    std::vector<gp_XYZ> line{{0, 0, 0}, {1, 2, 5}, {2, 4, 0}, {3, 6, 9}};
    expect(sectionLineResid(line, gp_XYZ(0, 0, 1)) <= 1e-12, "line resid collinear");
    const std::vector<gp_XYZ> arc = arcPts(10, 0, kPi / 3, 6, 5);
    expect(sectionLineResid(arc, gp_XYZ(0, 0, 1)) > kTau, "line resid arc");
    // coaxialSameRadius.
    JointAxisFit f;
    f.radius = {10.0, 8.0};
    expect(!coaxialSameRadius(f, 1, kTau), "same radius: 8 vs 10");
    f.radius = {10.0, 10.0 + 0.5 * kTau};
    expect(coaxialSameRadius(f, 1, kTau), "same radius within tau");
    // sectionContactSlack = tanDistTol / r.
    expect(std::fabs(sectionContactSlack(3 * kTau, kTau, 2.0) - 2 * kTau) <= 1e-18, "contact slack");
    // sectionFilletCentre, line|circle: circle (0,0) R10 inside, line y = 5
    // below it; r = 1 puts the centre at y = 4, |p| = 9.
    SectionAnchor C, L;
    C.circle = true;
    C.R = 10;
    C.side = -1;
    L.nx = 0;
    L.ny = 1;
    L.d = 5;
    L.side = -1;
    double px = 0, py = 0;
    const bool lc = sectionFilletCentre(C, L, 1.0, 8, 4, px, py);
    expect(lc && std::fabs(py - 4.0) <= 1e-12 && std::fabs(px - std::sqrt(65.0)) <= 1e-12,
           "fillet centre line|circle");
    // circle|circle external: (0,0) R5 and (20,0) R5, r = 6.
    SectionAnchor A, B;
    A.circle = B.circle = true;
    A.R = B.R = 5;
    B.cx = 20;
    const bool cc = sectionFilletCentre(A, B, 6.0, 10, 5, px, py);
    expect(cc && std::fabs(px - 10.0) <= 1e-12 && std::fabs(py - std::sqrt(21.0)) <= 1e-12,
           "fillet centre circle|circle");
    SectionAnchor L2 = L;
    expect(!sectionFilletCentre(L, L2, 1.0, 0, 0, px, py), "plane|plane is not a section fillet");
    // sectionConcentric.
    SectionAnchor in = A, out = A;
    in.R = 10.0;
    out.R = 12.5;
    in.tol = out.tol = kTau;
    expect(sectionConcentric(in, out), "concentric");
    out.cx = 1e-3;
    expect(!sectionConcentric(in, out), "not concentric");
    // D-fc-5.17 (a). equalStepTol = 4τ/r. A rail needs minStations(kParamsCircle).
    const double est = equalStepTol(kTau, 2.0);
    expect(std::fabs(est - 2.0 * kTau) <= 1e-18, "equal step tol 4τ/r");
    expect(equalStepTol(kTau, 0.0) == 0.0 && equalStepTol(0.0, 2.0) == 0.0, "equal step tol refuses");
    // est + est - est is est, so the bound comparison does not depend on a
    // decimal that is not a binary fraction.
    expect(equalStepRail(minStations(kParamsCircle), est, est + est, kTau, 2.0),
           "equal step rail at the bound");
    expect(!equalStepRail(minStations(kParamsCircle) - 1, est, est, kTau, 2.0),
           "equal step rail under stations");
    expect(!equalStepRail(minStations(kParamsCircle), est, est + est + est, kTau, 2.0),
           "equal step rail over the bound");
    expect(sectionTangentLine(0.0, -2.0, 2.0, kTau, kTau, 4.0), "section tangent line");
    expect(!sectionTangentLine(1.0, -2.0, 2.0, kTau, kTau, 4.0), "section line not parallel");
    expect(!sectionTangentLine(0.0, -3.0, 2.0, kTau, kTau, 4.0), "section line not one radius off");
    int side = 0;
    expect(sectionTangentCircle(0.0, 12.0, 10.0, 2.0, kTau, kTau, 4.0, side) && side == 1,
           "section tangent circle external");
    expect(sectionTangentCircle(0.0, 8.0, 10.0, 2.0, kTau, kTau, 4.0, side) && side == -1,
           "section tangent circle internal");
    expect(!sectionTangentCircle(1.0, 12.0, 10.0, 2.0, kTau, kTau, 4.0, side),
           "section circle axes not parallel");
    double nx = 0, ny = 0, d = 0;
    sectionTangentPlane(0.0, 0.0, 2.0, kPi / 2.0, nx, ny, d);
    expect(std::fabs(nx) <= 1e-12 && std::fabs(ny - 1.0) <= 1e-12 && std::fabs(d - 2.0) <= 1e-12,
           "section tangent plane");
}

static void caseCoaxialCert() {
    const gp_XYZ o(0, 0, 0), a(0, 0, 1);
    const std::vector<std::vector<gp_XYZ>> fam{arcPts(10, 0, kPi / 3, 12, 20)};
    CoaxialCert c = certifyCoaxialBand(fam, arcPts(8, 0, kPi / 3, 12, 20), o, a, kTau, 8);
    expect(c.ok && !c.refit && std::strcmp(c.why, "accept") == 0, "coaxial accept on family axis");
    expect(c.fit.radius.size() == 2 && std::fabs(c.fit.radius.back() - 8.0) <= kTau,
           "coaxial radius");
    c = certifyCoaxialBand(fam, arcPts(8, 0, kPi / 3, 2, 20), o, a, kTau, 8);
    expect(!c.ok && std::strcmp(c.why, "underdetermined") == 0, "coaxial 3 stations underdetermined");
    // A huge circle over a short span: its sagitta is under tau, a line.
    c = certifyCoaxialBand(fam, arcPts(1e5, 0, 1e-6, 4, 20), o, a, kTau, 8);
    expect(!c.ok && std::strcmp(c.why, "line") == 0, "coaxial flat strip is a line");
    std::vector<gp_XYZ> off = arcPts(8, 0, kPi / 3, 12, 20);
    off[7] = off[7] * (1.0 + 3 * kTau / 8.0);
    c = certifyCoaxialBand(fam, off, o, a, kTau, 8);
    expect(!c.ok && std::strcmp(c.why, "vertex-off") == 0, "coaxial vertex off tau refused");
    // The family axis seed is 1.5 tau off the true axis: only the joint least
    // squares over donor plus band certifies, on the true axis.
    const std::vector<std::vector<gp_XYZ>> famFull{arcPts(10, 0, 1.9 * kPi, 40, 20)};
    c = certifyCoaxialBand(famFull, arcPts(8, 0, kPi / 3, 12, 20), gp_XYZ(1.5 * kTau, 0, 0), a,
                           kTau, 16);
    expect(c.ok && c.refit, "coaxial joint least squares");
    expect(std::hypot(c.fit.origin.X(), c.fit.origin.Y()) <= kTau, "coaxial joint axis");
}

// D-fc-5.1 end to end: a certified R10 arc (donor), a coaxial R8 arc band and,
// when `plane`, the plane x = 8 tangent to the band at its 0-degree end.
static int lineModeCase(bool plane, bool closed, std::string& why) {
    Synth s;
    SegmentWork w;
    const double h = 20;
    std::vector<int> donorTris, bandTris, planeTris;
    const auto dc = arcCols(s, 0, 0, 10, 0, kPi / 3, 12, h);
    for (std::size_t i = 0; i + 1 < dc.size(); ++i)
        s.quad(dc[i][0], dc[i][1], dc[i + 1][0], dc[i + 1][1], donorTris);
    const double span = closed ? 2 * kPi : kPi / 3;
    const int nb = closed ? 72 : 12;
    auto bc = arcCols(s, 0, 0, 8, 0, span, nb, h);
    if (closed) bc.back() = bc.front();
    for (std::size_t i = 0; i + 1 < bc.size(); ++i)
        s.quad(bc[i][0], bc[i][1], bc[i + 1][0], bc[i + 1][1], bandTris);
    if (plane) {
        const int p0 = s.add(8, -10, 0), p1 = s.add(8, -10, h);
        s.quad(p0, p1, bc[0][0], bc[0][1], planeTris);
    }
    s.finish();
    Region d;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = 10;
    d.tris = donorTris;
    w.accepted.push_back(d);
    Provisional band;
    band.claim = ProvClaim::Unclaimed;
    band.tris = bandTris;
    band.plane = gp_Ax3(gp_Pnt(8, 0, 0), gp_Dir(1, 0, 0));
    band.maxVertexDev = 1.0;  // a curved band is no plane anchor
    w.provisionals.push_back(band);
    if (plane) {
        Provisional pl;
        pl.claim = ProvClaim::Unclaimed;
        pl.tris = planeTris;
        pl.plane = gp_Ax3(gp_Pnt(8, -5, 10), gp_Dir(1, 0, 0));
        pl.area = 10 * h;
        pl.maxVertexDev = 0;
        w.provisionals.push_back(pl);
    }
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, w, RelRevolveMode::Line);
    int n8 = 0;
    for (const Region& r : w.accepted)
        if (r.type == SurfType::Cylinder && r.origin == Origin::FilletStrip &&
            std::fabs(r.radius - 8.0) <= kTau && r.tris.size() == bandTris.size() &&
            r.maxVertexDev <= kTau)
            ++n8;
    why = n8 ? "accept" : "refused";
    return n8;
}

static void caseLineMode() {
    std::string why;
    expect(lineModeCase(true, false, why) == 1, "coaxial band anchored by its tangent plane");
    expect(lineModeCase(false, false, why) == 0, "coaxial band with no plane anchor refused");
    expect(lineModeCase(false, true, why) == 0, "closed coaxial band deferred");
}

// D-fc-5.2 end to end: circle anchor R10 about the origin, plane x = 3, and
// the r = 1 fillet tangent to both (centre (4, sqrt(105))). An accepted CylGrow
// atom holds all but the fillet's last strip at the circle; that strip is
// Unclaimed. `complete`: the atom holds every strip.
static void caseSectionFillet(bool complete) {
    Synth s;
    SegmentWork w;
    const double h = 10, r = 1, cy = std::sqrt(105.0);
    const double phL = kPi, phC = std::atan2(10.0 / 11.0 * cy - cy, 10.0 / 11.0 * 4.0 - 4.0) + 2 * kPi;
    const auto fc = arcCols(s, 4, cy, r, phL, phC, 10, h);
    std::vector<int> atomTris, stripTris, circTris, planeTris;
    for (std::size_t i = 0; i + 1 < fc.size(); ++i) {
        std::vector<int>& dst = (!complete && i + 2 == fc.size()) ? stripTris : atomTris;
        s.quad(fc[i][0], fc[i][1], fc[i + 1][0], fc[i + 1][1], dst);
    }
    // Circle anchor from its contact (fillet's last column) on to 120 degrees.
    const double a0 = std::atan2(10.0 / 11.0 * cy, 10.0 / 11.0 * 4.0);
    auto cc = arcCols(s, 0, 0, 10, a0, 2 * kPi / 3, 12, h);
    cc.front() = fc.back();
    for (std::size_t i = 0; i + 1 < cc.size(); ++i)
        s.quad(cc[i][0], cc[i][1], cc[i + 1][0], cc[i + 1][1], circTris);
    // Plane x = 3 from the line contact (fillet's first column) up to y = 20.
    const int p0 = s.add(3, 20, 0), p1 = s.add(3, 20, h);
    s.quad(fc.front()[0], fc.front()[1], p0, p1, planeTris);
    s.finish();
    Region circ;
    circ.type = SurfType::Cylinder;
    circ.origin = Origin::FilletStrip;
    circ.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    circ.radius = 10;
    circ.tris = circTris;
    w.accepted.push_back(circ);
    Region atom;
    atom.type = SurfType::Cylinder;
    atom.origin = Origin::CylGrow;
    atom.ax = gp_Ax3(gp_Pnt(4, cy, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    atom.radius = r;
    atom.tris = atomTris;
    w.accepted.push_back(atom);
    Provisional pl;
    pl.claim = ProvClaim::Unclaimed;
    pl.tris = planeTris;
    pl.plane = gp_Ax3(gp_Pnt(3, 15, 5), gp_Dir(-1, 0, 0));
    pl.area = 10 * h;
    w.provisionals.push_back(pl);
    if (!stripTris.empty()) {
        Provisional st;
        st.claim = ProvClaim::Unclaimed;
        st.tris = stripTris;
        st.plane = gp_Ax3(gp_Pnt(3.5, cy - 0.9, 5), gp_Dir(0.6, -0.8, 0));
        st.maxVertexDev = 1.0;
        w.provisionals.push_back(st);
    }
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    SectionCircle c;
    c.origin = gp_XYZ(0, 0, 0);
    c.axis = gp_XYZ(0, 0, 1);
    c.radius = 10;
    c.tol = kTau;
    c.length = h;
    c.acc = 0;
    claimSectionFillets(mv, sp, tol, w, {c});
    const Region& R = w.accepted[1];
    if (complete) {
        expect(R.origin == Origin::CylGrow && R.tris.size() == atomTris.size(),
               "section fillet: complete atom unchanged (leave-one-out underdetermined)");
        return;
    }
    expect(R.origin == Origin::FilletStrip, "section fillet claimed (line|circle)");
    expect(R.tris.size() == atomTris.size() + stripTris.size(), "section fillet absorbs the strip");
    expect(std::fabs(R.radius - r) <= kTau && R.maxVertexDev <= kTau, "section fillet r within tau");
    const gp_XYZ loc = R.ax.Location().XYZ();
    expect(std::hypot(loc.X() - 4.0, loc.Y() - cy) <= kTau, "section fillet centre on the offsets");
    expect(w.provisionals.size() >= 2 && w.provisionals[1].claim == ProvClaim::ConsumedFillet,
           "section fillet strip provisional consumed");
}

static void casePlaneAnchor() {
    // A wall whose A2 plane absorbed a tilted facet: the peel certifies the
    // wall on its maximum-area facet's plane (D-fc-4.1).
    Synth s;
    const int a = s.add(0, 0, 0), b = s.add(10, 0, 0), c = s.add(10, 0, 10), d = s.add(0, 0, 10);
    const int e = s.add(11, 0.5, 0), f = s.add(11, 0.5, 10);
    std::vector<int> tris{s.addTri(a, b, c), s.addTri(a, c, d), s.addTri(b, e, f), s.addTri(b, f, c)};
    s.finish();
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = tris;
    pr.plane = gp_Ax3(gp_Pnt(5, 0.1, 5), gp_Dir(0, 1, 0));
    pr.area = 110;
    pr.maxVertexDev = 0.4;
    PlaneAnchor pa;
    expect(planeAnchorOf(mv, sp, tol, pr, pa) && pa.certified && pa.body.size() == 2 &&
               pa.dev <= kTau,
           "plane anchor peel certifies the wall");
    pr.claim = ProvClaim::ConsumedFillet;
    expect(!planeAnchorOf(mv, sp, tol, pr, pa), "plane anchor needs an Unclaimed provisional");
}

static void caseOsbGuards() {
    // D-fc-5.5 F-1: past the fitted contact OR the station generator, tau slack.
    expect(osbPastContact(-2 * kTau, 0.5 * kTau, 1.0, kTau), "osb slab: past the station generator");
    expect(osbPastContact(0.0, -2 * kTau, 1.0, kTau), "osb slab: past the fitted contact");
    expect(!osbPastContact(-2 * kTau, -2 * kTau, 1.0, kTau), "osb slab: behind both generators");
    expect(osbPastContact(2 * kTau, 2 * kTau, -1.0, kTau) == false, "osb slab: band side from sg");
    // D-fc-5.5 F-2: stations left without any one source. An extruded arc of 5
    // columns (stations) at 2 levels, one quad (2 triangles) between columns.
    auto quad = [](double R, double a0, double a1, double h) {
        const gp_XYZ p0(R * std::cos(a0), R * std::sin(a0), 0), p1(R * std::cos(a0), R * std::sin(a0), h);
        const gp_XYZ q0(R * std::cos(a1), R * std::sin(a1), 0), q1(R * std::cos(a1), R * std::sin(a1), h);
        return std::vector<std::array<gp_XYZ, 3>>{{p0, q0, p1}, {p1, q0, q1}};
    };
    const gp_XYZ ax(0, 0, 1);
    std::vector<std::array<gp_XYZ, 3>> tp;
    std::vector<int> ts;
    for (int i = 0; i < 4; ++i) {
        for (const auto& t : quad(10, 0.1 * i, 0.1 * (i + 1), 5)) {
            tp.push_back(t);
            ts.push_back(i < 2 ? 7 : 9);
        }
    }
    expect(leaveOneSourceOutStations(tp, ts, ax, kTau) == 3, "leave-one-source-out: 3 stations left");
    std::vector<int> one(ts.size(), 4);
    expect(leaveOneSourceOutStations(tp, one, ax, kTau) == 0, "leave-one-source-out: one source leaves 0");
    // I-starved: two quads from two sources, each alone 2 stations (a line).
    std::vector<std::array<gp_XYZ, 3>> tp2;
    std::vector<int> ts2;
    for (int i : {0, 2})
        for (const auto& t : quad(10, 0.1 * i, 0.1 * (i + 1), 5)) {
            tp2.push_back(t);
            ts2.push_back(i);
        }
    expect(leaveOneSourceOutStations(tp2, ts2, ax, kTau) == 2,
           "leave-one-source-out: two one-quad sources underdetermined");
    // The band holding an anchor's certified triangle.
    std::vector<char> inBand(6, 0);
    inBand[4] = 1;
    expect(bandHoldsAny({1, 4}, inBand) && !bandHoldsAny({0, 1, 2}, inBand) && !bandHoldsAny({9}, inBand),
           "band holds an anchor body triangle");
}

// D-fc-5.5 partial cones. Family: one R10 arc (0..60 deg, z 0..20). Cone:
// rho = 10 - (z - 20) over z 20..21 (45 deg), the same azimuths.
static std::vector<gp_XYZ> coneRingPts(double a0, double a1, int n, std::vector<double> zs,
                                       double r0 = 10, double z0 = 20) {
    std::vector<gp_XYZ> v;
    for (int i = 0; i <= n; ++i) {
        const double a = a0 + (a1 - a0) * (double)i / (double)n;
        for (double z : zs) {
            const double r = r0 - (z - z0);
            v.emplace_back(r * std::cos(a), r * std::sin(a), z);
        }
    }
    return v;
}

static void casePartialCone() {
    const gp_XYZ o(0, 0, 0), a(0, 0, 1);
    const std::vector<std::vector<gp_XYZ>> fam{arcPts(10, 0, kPi / 3, 12, 20)};
    PartialConeCert c = certifyPartialCone(fam, coneRingPts(0, kPi / 3, 12, {20, 20.5, 21}), o, a, kTau);
    expect(c.ok && std::strcmp(c.why, "accept") == 0, "partial cone accept on family axis");
    expect(std::fabs(c.k + 1.0) <= 1e-9 && std::fabs(c.rhoAt(20) - 10) <= kTau &&
               std::fabs(c.rhoAt(21) - 9) <= kTau,
           "partial cone profile line rho = 30 - z");
    c = certifyPartialCone(fam, coneRingPts(0, kPi / 3, 12, {20, 21}), o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "underdetermined") == 0, "partial cone 2 profile stations underdetermined");
    // Two rings whose vertices scatter 1.7 tau across the profile line: four
    // tau-merged profile points, but two stations along the line.
    std::vector<gp_XYZ> noisy = coneRingPts(0, kPi / 3, 12, {20, 21});
    for (std::size_t i = 0; i < noisy.size(); ++i) {
        const double dn = ((i / 2) % 2 ? 0.85 : -0.85) * kTau / std::sqrt(2.0);
        noisy[i] += gp_XYZ(dn * std::cos(kPi / 36 * (double)(i / 2)), dn * std::sin(kPi / 36 * (double)(i / 2)), dn);
    }
    c = certifyPartialCone(fam, noisy, o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "underdetermined") == 0 && c.stations == 2,
           "partial cone: scatter across the line is not a station");
    c = certifyPartialCone(fam, coneRingPts(0, kPi / 3, 2, {20, 20.5, 21}), o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "underdetermined") == 0, "partial cone 3 azimuth stations underdetermined");
    // 1e-3 rad at R ~10: columns 2.5e-3 mm apart, sagitta 1.25e-6 mm < tau.
    c = certifyPartialCone(fam, coneRingPts(0, 1e-3, 4, {20, 20.5, 21}), o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "line") == 0, "partial cone planar strip is a line");
    std::vector<gp_XYZ> off = coneRingPts(0, kPi / 3, 12, {20, 20.5, 21});
    off[10] = off[10] * (1.0 + 3 * kTau / off[10].Modulus());
    c = certifyPartialCone(fam, off, o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "vertex-off") == 0, "partial cone vertex off tau refused");
    std::vector<std::vector<gp_XYZ>> famOff = fam;
    famOff[0][5] = famOff[0][5] * (1.0 + 3 * kTau / 10.0);
    c = certifyPartialCone(famOff, coneRingPts(0, kPi / 3, 12, {20, 20.5, 21}), o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "vertex-off") == 0, "partial cone family vertex off tau refused");
    std::vector<gp_XYZ> cyl8;
    for (int i = 0; i <= 12; ++i)
        for (double z : {0.0, 10.0, 20.0})
            cyl8.emplace_back(8 * std::cos(kPi / 36 * i), 8 * std::sin(kPi / 36 * i), z);
    c = certifyPartialCone(fam, cyl8, o, a, kTau);
    expect(!c.ok && std::strcmp(c.why, "same-surface") == 0, "partial cone of no radius span is a cylinder");

    // Builder predicates.
    Synth s0;
    s0.add(0, 0, 0);
    s0.finish();
    const MeshView mv0 = s0.view();
    Region cone;
    cone.type = SurfType::Cone;
    cone.origin = Origin::ChamferCone;
    cone.ax = gp_Ax3(gp_Pnt(0, 0, 21), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    cone.radius = 9;
    cone.vMin = 9;
    cone.vMax = 10;
    cone.maxVertexSnap = -1;
    cone.closed360 = false;
    Region cyl;
    cyl.type = SurfType::Cylinder;
    cyl.origin = Origin::FilletStrip;
    cyl.ax = gp_Ax3(gp_Pnt(0, 0, 10), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    cyl.radius = 10 + kTau;
    expect(partialConeRimCoaxial(cyl, cone, mv0), "partial cone rim on a coaxial cylinder");
    Region cylR = cyl;
    cylR.radius = 10.5;
    expect(!partialConeRimCoaxial(cylR, cone, mv0), "partial cone rim: other radius");
    Region cylS = cyl;
    cylS.ax = gp_Ax3(gp_Pnt(3 * kTau, 0, 10), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    expect(!partialConeRimCoaxial(cylS, cone, mv0), "partial cone rim: shifted axis");
    Region cylT = cyl;
    cylT.ax = gp_Ax3(gp_Pnt(0, 0, 10), gp_Dir(0.001, 0, 1), gp_Dir(1, 0, -0.001));
    expect(!partialConeRimCoaxial(cylT, cone, mv0), "partial cone rim: tilted axis");
    Region closedCone = cone;
    closedCone.closed360 = true;
    expect(partialConeSide(&cyl, &cone) && partialConeSide(&cone, nullptr) &&
               !partialConeSide(&cyl, &closedCone) && !partialConeSide(&cyl, nullptr),
           "partial cone side");

    // Line mode end to end: donor R10 arc band, a 45 deg roof above it.
    Synth s;
    SegmentWork w;
    const double h = 20;
    std::vector<int> donorTris, coneTris;
    const auto dc = arcCols(s, 0, 0, 10, 0, kPi / 3, 12, h);
    for (std::size_t i = 0; i + 1 < dc.size(); ++i)
        s.quad(dc[i][0], dc[i][1], dc[i + 1][0], dc[i + 1][1], donorTris);
    std::vector<std::array<int, 3>> rings;
    for (int i = 0; i <= 12; ++i) {
        const double ang = kPi / 3 * (double)i / 12.0;
        rings.push_back({{dc[(std::size_t)i][1], s.add(9.5 * std::cos(ang), 9.5 * std::sin(ang), 20.5),
                          s.add(9 * std::cos(ang), 9 * std::sin(ang), 21)}});
    }
    for (int i = 0; i < 12; ++i)
        for (int k = 0; k < 2; ++k)
            s.quad(rings[(std::size_t)i][(std::size_t)k], rings[(std::size_t)i][(std::size_t)k + 1],
                   rings[(std::size_t)i + 1][(std::size_t)k], rings[(std::size_t)i + 1][(std::size_t)k + 1],
                   coneTris);
    s.finish();
    Region d;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = 10;
    d.tris = donorTris;
    w.accepted.push_back(d);
    Provisional roof;
    roof.claim = ProvClaim::Unclaimed;
    roof.tris = coneTris;
    roof.plane = gp_Ax3(gp_Pnt(9.5, 0, 20.5), gp_Dir(0.7071, 0, 0.7071));
    roof.maxVertexDev = 1.0;
    w.provisionals.push_back(roof);
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, w, RelRevolveMode::Line);
    int nCone = 0;
    for (const Region& r : w.accepted) {
        if (r.type != SurfType::Cone || r.origin != Origin::ChamferCone) continue;
        ++nCone;
        expect(!r.closed360 && std::fabs(r.radius - 9.0) <= kTau && std::fabs(r.vMax - 10.0) <= kTau &&
                   std::fabs((r.uMax - r.uMin) - kPi / 3) <= 1e-6 && r.tris.size() == coneTris.size() &&
                   r.maxVertexDev <= kTau,
               "partial cone region: Rlo 9, Rhi 10, span 60 deg, every roof triangle");
    }
    expect(nCone == 1, "partial cone claimed in Line mode");
}

// audit-fx-arch F-1: geometric refusals and the preimage closure.
static void caseGeomFailed() {
    expect(lineChainGeometric("vertex-off") && lineChainGeometric("line"), "geometric: vertex-off, line");
    expect(!lineChainGeometric("underdetermined") && !lineChainGeometric("same-surface") &&
               !lineChainGeometric("partial-cone") && !lineChainGeometric(nullptr),
           "not geometric: underdetermined, same-surface, partial-cone");
    std::vector<char> failed(5, 0), own(5, 0);
    failed[2] = failed[3] = 1;
    own[3] = 1;
    expect(preimageAdmits(failed, own, 1) && !preimageAdmits(failed, own, 2) &&
               preimageAdmits(failed, own, 3) && !preimageAdmits(failed, own, -1),
           "closure never takes a failed chain's triangle outside its own seed");
}

static int nFillet(const SegmentWork& w) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.origin == Origin::FilletStrip) ++n;
    return n;
}

// D-fc-5.17 (a). Quarter circle, six equal 15° chords, one certified plane at
// the 0° end, and a shatter-class strip in the tangent plane at the far end.
// No torus atom. areaMin sits between the strip and the anchor.
static void caseEqualStepOneAnchor() {
    Synth s;
    SegmentWork w;
    const double h = 4, R = 2;
    std::vector<int> arcTris, planeTris, stripTris;
    const auto fc = arcCols(s, 0, 0, R, 0, kPi / 2, 6, h);
    for (std::size_t i = 0; i + 1 < fc.size(); ++i)
        s.quad(fc[i][0], fc[i][1], fc[i + 1][0], fc[i + 1][1], arcTris);
    const int p0 = s.add(R, -6, 0), p1 = s.add(R, -6, h);
    s.quad(fc.front()[0], fc.front()[1], p0, p1, planeTris);
    const int q0 = s.add(-0.4, R, 0), q1 = s.add(-0.4, R, h);
    s.quad(fc.back()[0], fc.back()[1], q0, q1, stripTris);
    s.finish();
    Provisional arc;
    arc.claim = ProvClaim::Unclaimed;
    arc.tris = arcTris;
    arc.maxVertexDev = 1.0;
    arc.plane = gp_Ax3(gp_Pnt(R, 0, h / 2), gp_Dir(1, 0, 0));
    w.provisionals.push_back(arc);
    Provisional pl;
    pl.claim = ProvClaim::Unclaimed;
    pl.tris = planeTris;
    pl.area = 6 * h;
    pl.maxVertexDev = 0;
    pl.plane = gp_Ax3(gp_Pnt(R, -3, h / 2), gp_Dir(1, 0, 0));
    w.provisionals.push_back(pl);
    Provisional st;
    st.claim = ProvClaim::Unclaimed;
    st.tris = stripTris;
    st.area = 0.4 * h;
    st.maxVertexDev = 0;
    st.plane = gp_Ax3(gp_Pnt(-0.2, R, h / 2), gp_Dir(0, 1, 0));
    w.provisionals.push_back(st);
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    tol.epsPlane = 0.05;  // areaMin = 5: the strip is shatter, the anchor is not
    claimSectionFillets(mv, sp, tol, w, {});
    int got = 0;
    for (const Region& r : w.accepted) {
        if (r.origin != Origin::FilletStrip) continue;
        ++got;
        const gp_XYZ loc = r.ax.Location().XYZ();
        expect(std::fabs(r.radius - R) <= kTau && r.maxVertexDev <= kTau, "one anchor: r within tau");
        expect(std::hypot(loc.X(), loc.Y()) <= kTau, "one anchor: centre");
        expect(r.tris.size() == arcTris.size(), "one anchor: the arc, not the strip");
        expect((r.filletNbrA == 2 || r.filletNbrB == 2) && (r.filletNbrA == 1 || r.filletNbrB == 1),
               "one anchor: names the plane and the far strip");
    }
    expect(got == 1, "one anchor: one section fillet");
    expect(w.provisionals.size() >= 3 && w.provisionals[1].claim == ProvClaim::Unclaimed &&
               w.provisionals[2].claim == ProvClaim::Unclaimed,
           "one anchor: plane and strip stay unclaimed");
    expect(!a3WouldCommit(mv, tol, w.provisionals[2], false) &&
               a3WouldCommit(mv, tol, w.provisionals[2], true),
           "one anchor: the strip commits only as a fillet flank");
}

// Graded chords on a true circle fail the equal-step law. No anchor at all
// fails for lack of a certificate. Neither claims a fillet.
static void caseEqualStepRefusals() {
    auto graded = []() {
        Synth s;
        SegmentWork w;
        const double h = 4, R = 2;
        std::vector<double> az{0};
        for (double deg : {8.0, 10.0, 12.0, 14.0, 16.0, 18.0})
            az.push_back(az.back() + deg * kPi / 180.0);
        std::vector<std::array<int, 2>> cols;
        for (double a : az) {
            cols.push_back({{s.add(R * std::cos(a), R * std::sin(a), 0.0),
                             s.add(R * std::cos(a), R * std::sin(a), h)}});
        }
        std::vector<int> arcTris, planeTris;
        for (std::size_t i = 0; i + 1 < cols.size(); ++i)
            s.quad(cols[i][0], cols[i][1], cols[i + 1][0], cols[i + 1][1], arcTris);
        const int p0 = s.add(R, -6, 0), p1 = s.add(R, -6, h);
        s.quad(cols.front()[0], cols.front()[1], p0, p1, planeTris);
        s.finish();
        Provisional arc;
        arc.claim = ProvClaim::Unclaimed;
        arc.tris = arcTris;
        arc.maxVertexDev = 1.0;
        w.provisionals.push_back(arc);
        Provisional pl;
        pl.claim = ProvClaim::Unclaimed;
        pl.tris = planeTris;
        pl.area = 6 * h;
        pl.maxVertexDev = 0;
        pl.plane = gp_Ax3(gp_Pnt(R, -3, h / 2), gp_Dir(1, 0, 0));
        w.provisionals.push_back(pl);
        return std::make_pair(std::move(s), std::move(w));
    };
    {
        auto sw = graded();
        const MeshView mv = sw.first.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol);
        tol.epsPlane = 0.05;
        claimSectionFillets(mv, sp, tol, sw.second, {});
        expect(nFillet(sw.second) == 0, "graded steps: no section fillet");
    }
    {
        Synth s;
        SegmentWork w;
        const double h = 4, R = 2;
        std::vector<int> arcTris;
        const auto fc = arcCols(s, 0, 0, R, 0, kPi / 2, 6, h);
        for (std::size_t i = 0; i + 1 < fc.size(); ++i)
            s.quad(fc[i][0], fc[i][1], fc[i + 1][0], fc[i + 1][1], arcTris);
        s.finish();
        Provisional arc;
        arc.claim = ProvClaim::Unclaimed;
        arc.tris = arcTris;
        arc.maxVertexDev = 1.0;
        w.provisionals.push_back(arc);
        const MeshView mv = s.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol);
        claimSectionFillets(mv, sp, tol, w, {});
        expect(nFillet(w) == 0, "equal steps with no certified anchor: no section fillet");
    }
}

// D-fc-5.2 circle|circle with no atom. Internal to R10, external to R3.
static void caseEqualStepCircleCircle() {
    Synth s;
    SegmentWork w;
    const double h = 4, r = 2;
    std::vector<int> arcTris, t10, t3;
    const auto fc = arcCols(s, 8, 0, r, 0, kPi / 2, 6, h);
    for (std::size_t i = 0; i + 1 < fc.size(); ++i)
        s.quad(fc[i][0], fc[i][1], fc[i + 1][0], fc[i + 1][1], arcTris);
    auto c10 = arcCols(s, 0, 0, 10, 0, kPi / 3, 6, h);
    c10.front() = fc.front();
    for (std::size_t i = 0; i + 1 < c10.size(); ++i)
        s.quad(c10[i][0], c10[i][1], c10[i + 1][0], c10[i + 1][1], t10);
    auto c3 = arcCols(s, 8, 5, 3, -kPi / 2, -kPi / 2 + kPi / 3, 6, h);
    c3.front() = fc.back();
    for (std::size_t i = 0; i + 1 < c3.size(); ++i)
        s.quad(c3[i][0], c3[i][1], c3[i + 1][0], c3[i + 1][1], t3);
    s.finish();
    Region a10;
    a10.type = SurfType::Cylinder;
    a10.origin = Origin::FilletStrip;
    a10.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    a10.radius = 10;
    a10.tris = t10;
    w.accepted.push_back(a10);
    Region a3;
    a3.type = SurfType::Cylinder;
    a3.origin = Origin::FilletStrip;
    a3.ax = gp_Ax3(gp_Pnt(8, 5, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    a3.radius = 3;
    a3.tris = t3;
    w.accepted.push_back(a3);
    Provisional arc;
    arc.claim = ProvClaim::Unclaimed;
    arc.tris = arcTris;
    arc.maxVertexDev = 1.0;
    w.provisionals.push_back(arc);
    SectionCircle s10, s3;
    s10.origin = gp_XYZ(0, 0, 0);
    s10.axis = gp_XYZ(0, 0, 1);
    s10.radius = 10;
    s10.tol = kTau;
    s10.length = h;
    s10.acc = 0;
    s3.origin = gp_XYZ(8, 5, 0);
    s3.axis = gp_XYZ(0, 0, 1);
    s3.radius = 3;
    s3.tol = kTau;
    s3.length = h;
    s3.acc = 1;
    const MeshView mv = s.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimSectionFillets(mv, sp, tol, w, {s10, s3});
    int got = 0;
    for (int i = 2; i < (int)w.accepted.size(); ++i) {
        const Region& r = w.accepted[(std::size_t)i];
        if (r.origin != Origin::FilletStrip || r.tris.size() != arcTris.size()) continue;
        ++got;
        const gp_XYZ loc = r.ax.Location().XYZ();
        expect(std::fabs(r.radius - 2.0) <= kTau && r.maxVertexDev <= kTau, "circle|circle: r within tau");
        expect(std::hypot(loc.X() - 8.0, loc.Y()) <= kTau, "circle|circle: centre");
    }
    expect(got == 1, "circle|circle: one section fillet, no atom");
}

int main() {
    caseGeomFailed();
    casePredicates();
    caseOsbGuards();
    casePartialCone();
    caseCoaxialCert();
    caseLineMode();
    caseSectionFillet(false);
    caseSectionFillet(true);
    casePlaneAnchor();
    caseEqualStepOneAnchor();
    caseEqualStepRefusals();
    caseEqualStepCircleCircle();
    std::fprintf(stderr, "relational_arc_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
