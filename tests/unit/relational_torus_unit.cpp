// Synthetic torus-relation cases (D-fc-4.6 F1/F2/F4). No training mesh.
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "posix_compat.hpp"
#include "refit_relational.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gp_Ax1.hxx>
#include <gp_Lin.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;
    double q = 1e-6;

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
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2], triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]), hi = std::max(v[s], v[s + 1]);
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
        mv.quantFloor = q;
        return mv;
    }
};

static std::vector<int> ringOf(Synth& s, double R, double z, int N, int cols) {
    std::vector<int> ix;
    for (int i = 0; i < cols; ++i) {
        const double th = 2.0 * kPi * (double)i / (double)N;
        ix.push_back(s.add(R * std::cos(th), R * std::sin(th), z));
    }
    return ix;
}

// Rows of quads between consecutive rings; `sectors` quads per row.
static void row(Synth& s, const std::vector<int>& a, const std::vector<int>& b, int sectors, bool wrap,
                std::vector<int>& tris) {
    const int n = (int)a.size();
    for (int i = 0; i < sectors; ++i) {
        const int i1 = wrap ? (i + 1) % n : i + 1;
        tris.push_back(s.addTri(a[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i]));
        tris.push_back(s.addTri(b[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i1]));
    }
}

struct Part {
    Synth mesh;
    SegmentWork work;
    std::vector<std::vector<int>> rows;  // triangles of each torus row
    std::vector<int> donorTris;
};

// A closed cylinder donor (R = Rmaj + r at the bottom of the profile) and a
// quarter torus of nPhi rows above it, N columns, fully closed.
static Part ringTorus(int N, int nPhi) {
    Part p;
    const double Rmaj = 12.0, r = 2.0, z0 = 8.0;
    std::vector<std::vector<int>> rings;
    for (int j = 0; j <= nPhi; ++j) {
        const double phi = kPi + 0.5 * kPi * (double)j / (double)nPhi;
        rings.push_back(ringOf(p.mesh, Rmaj + r * std::cos(phi), z0 + r * std::sin(phi), N, N));
    }
    const double rhoC = Rmaj - r;
    auto lower = ringOf(p.mesh, rhoC, z0 - 6.0, N, N);
    row(p.mesh, lower, rings[0], N, true, p.donorTris);
    for (int j = 0; j < nPhi; ++j) {
        p.rows.emplace_back();
        row(p.mesh, rings[(std::size_t)j], rings[(std::size_t)j + 1], N, true, p.rows.back());
    }
    p.mesh.finish();
    Region d;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, z0 - 3.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = rhoC;
    d.closed360 = true;
    d.tris = p.donorTris;
    p.work.accepted.push_back(d);
    return p;
}

static Region atomOf(const std::vector<int>& tris) {
    Region a;
    a.type = SurfType::Cylinder;
    a.origin = Origin::CylGrow;
    a.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0));
    a.radius = 2.0;
    a.closed360 = false;
    a.tris = tris;
    std::sort(a.tris.begin(), a.tris.end());
    return a;
}

static void params(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kPi / 180.0;
}

static int countOrigin(const SegmentWork& w, Origin o, bool closed) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.origin == o && r.closed360 == closed) ++n;
    return n;
}

static const Region* torusOf(const SegmentWork& w) {
    for (const Region& r : w.accepted)
        if (r.origin == Origin::TorusBlend) return &r;
    return nullptr;
}

// F1/F2: rows 0..3 free (5 stations), rows 4..5 one atom wholly on the torus.
// The band is the whole chain; the atom joins whole and leaves the accepted set.
static void caseAtomJoins() {
    Part p = ringTorus(24, 6);
    Provisional pr;
    for (int j = 0; j < 4; ++j) pr.tris.insert(pr.tris.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
    std::sort(pr.tris.begin(), pr.tris.end());
    p.work.provisionals.push_back(pr);
    std::vector<int> at;
    for (int j = 4; j < 6; ++j) at.insert(at.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
    p.work.accepted.push_back(atomOf(at));
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    expect(claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus), "atom call");
    const Region* t = torusOf(p.work);
    expect(t != nullptr, "atom: torus claimed");
    if (!t) return;
    expect(t->tris.size() == 6u * 48u, "atom: band is the whole chain");
    expect(countOrigin(p.work, Origin::CylGrow, false) == 0, "atom: absorbed region left accepted");
    expect(t->torusPending >= 0 && p.work.torusPending[(std::size_t)t->torusPending].stashedCylinders.size() == 1,
           "atom: stashed for U-3");
    expect(std::fabs(t->radius - 12.0) <= tauOf(mv) && std::fabs(t->radius2 - 2.0) <= tauOf(mv), "atom: radii");
}

// F2: an atom with one vertex off the torus leaves whole; it keeps its triangles.
static void caseAtomOffSurface() {
    Part p = ringTorus(24, 6);
    Provisional pr;
    for (int j = 0; j < 4; ++j) pr.tris.insert(pr.tris.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
    std::sort(pr.tris.begin(), pr.tris.end());
    p.work.provisionals.push_back(pr);
    std::vector<int> at(p.rows[4].begin(), p.rows[4].end());
    at.insert(at.end(), p.rows[5].begin(), p.rows[5].end());
    // Lift one top vertex of the atom by 50 tau.
    const int gt = at.back();
    const int v = p.mesh.triFlat[3 * (std::size_t)gt + 2];
    p.mesh.pts[(std::size_t)v].SetZ(p.mesh.pts[(std::size_t)v].Z() + 100.0 * p.mesh.q);
    p.work.accepted.push_back(atomOf(at));
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const Region* t = torusOf(p.work);
    expect(t != nullptr, "off: free rows still certify");
    if (t) expect(t->tris.size() == 4u * 48u, "off: band keeps only the free rows");
    expect(countOrigin(p.work, Origin::CylGrow, false) == 1, "off: atom kept whole");
}

// F2 leave-one-region-out: free rows alone cannot determine the circle, so an
// atom that carries the rest of the stations does not join (rows 0..2 free =
// 4 stations, below certifyTorus's 5).
static void caseLeaveOneOut() {
    Part p = ringTorus(24, 6);
    Provisional pr;
    pr.tris = p.rows[0];
    std::sort(pr.tris.begin(), pr.tris.end());
    p.work.provisionals.push_back(pr);
    std::vector<int> at;
    for (int j = 1; j < 6; ++j) at.insert(at.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
    p.work.accepted.push_back(atomOf(at));
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    expect(torusOf(p.work) == nullptr, "loo: one atom cannot carry the certificate");
    expect(countOrigin(p.work, Origin::CylGrow, false) == 1, "loo: atom kept");
}

static void caseLeaveOneOutPure() {
    const double tau = 2e-6;
    std::vector<std::array<ProfPt, 3>> tp;
    std::vector<int> ta;
    int lv = 0;
    for (int st = 0; st < 8; ++st) {
        const double phi = kPi + 0.5 * kPi * (double)st / 7.0;
        for (int c = 0; c < 3; ++c) {
            std::array<ProfPt, 3> a;
            for (int k = 0; k < 3; ++k) {
                a[(std::size_t)k].rho = 12.0 + 2.0 * std::cos(phi);
                a[(std::size_t)k].z = 2.0 * std::sin(phi);
                a[(std::size_t)k].psi = 0.1 * (double)(c * 3 + k);
                a[(std::size_t)k].lv = lv++;
            }
            tp.push_back(a);
            ta.push_back(st < 4 ? 0 : 1);
        }
    }
    expect(!leaveOneOutDeterminate(tp, ta, 2, tau), "loo pure: half bands are each 4 stations");
    for (int& a : ta) a = (a == 0) ? 0 : -1;
    expect(leaveOneOutDeterminate(tp, ta, 1, tau) == false, "loo pure: 4 free stations");
    std::vector<int> small(ta.size(), -1);
    small[0] = 0;
    expect(leaveOneOutDeterminate(tp, small, 1, tau), "loo pure: a small atom leaves 8 stations");
    expect(leaveOneOutDeterminate(tp, small, 0, tau), "loo pure: no atom");
}

// F4: rails. A full ring is closed; a half ring is not.
static void caseRails() {
    for (int half = 0; half < 2; ++half) {
        Synth s;
        const int N = 24;
        const int cols = half ? N / 2 + 1 : N;
        std::vector<std::vector<int>> rings;
        for (int j = 0; j <= 4; ++j) {
            const double phi = kPi + 0.5 * kPi * (double)j / 4.0;
            rings.push_back(ringOf(s, 12.0 + 2.0 * std::cos(phi), 2.0 * std::sin(phi), N, cols));
        }
        std::vector<int> tris;
        for (int j = 0; j < 4; ++j)
            row(s, rings[(std::size_t)j], rings[(std::size_t)j + 1], half ? cols - 1 : N, !half, tris);
        s.finish();
        MeshView mv = s.view();
        Region t;
        t.type = SurfType::Torus;
        t.origin = Origin::TorusBlend;
        t.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        t.radius = 12.0;
        t.radius2 = 2.0;
        t.tris = tris;
        const TorusRails rl = torusRails(mv, t);
        expect(rl.full == !half, half ? "rails: half ring is partial" : "rails: full ring is closed");
        // reclaimRelational(Torus) normalises closure by the rails.
        SegmentWork w;
        t.closed360 = false;
        w.accepted.push_back(t);
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol);
        reclaimRelational(mv, sp, tol, w, RelRevolveMode::Torus);
        expect(w.accepted[0].closed360 == !half, half ? "reclaim: half stays open" : "reclaim: full ring closed");
    }
}

// F2 in reclaimRelational(Torus): an atom beside a certified torus, within
// tau, joins it; one off by 50 tau does not.
static void caseReclaimAbsorb() {
    Part p = ringTorus(24, 6);
    Region t;
    t.type = SurfType::Torus;
    t.origin = Origin::TorusBlend;
    t.ax = gp_Ax3(gp_Pnt(0, 0, 8.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    t.radius = 12.0;
    t.radius2 = 2.0;
    t.closed360 = true;
    t.uMax = 2.0 * kPi;
    for (int j = 0; j < 4; ++j) t.tris.insert(t.tris.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
    std::sort(t.tris.begin(), t.tris.end());
    t.nSides = 24;
    p.work.accepted.push_back(t);
    std::vector<int> near(p.rows[4].begin(), p.rows[4].end());
    p.work.accepted.push_back(atomOf(near));
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    reclaimRelational(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const Region* tr = torusOf(p.work);
    expect(tr && tr->tris.size() == 5u * 48u, "reclaim: adjacent atom absorbed whole");
    expect(countOrigin(p.work, Origin::CylGrow, false) == 0, "reclaim: atom removed");
    expect(tr && tr->maxVertexDev <= tauOf(mv), "reclaim: union within tau");
}

// F8 axis line: a closed torus touching a ChamferCone whose axis differs by
// rounding takes the cone's line (circle-on-torus is exactly coaxial only).
static void caseAxisLine() {
    Part p = ringTorus(24, 6);
    p.work.accepted[0].origin = Origin::ChamferCone;
    p.work.accepted[0].type = SurfType::Cone;
    p.work.accepted[0].ax = gp_Ax3(gp_Pnt(1e-9, 0, 5.0), gp_Dir(1e-12, 0, 1), gp_Dir(1, 0, 0));
    Region t;
    t.type = SurfType::Torus;
    t.origin = Origin::TorusBlend;
    t.ax = gp_Ax3(gp_Pnt(0, 0, 8.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    t.radius = 12.0;
    t.radius2 = 2.0;
    t.closed360 = true;
    t.uMax = 2.0 * kPi;
    for (const auto& rw : p.rows) t.tris.insert(t.tris.end(), rw.begin(), rw.end());
    std::sort(t.tris.begin(), t.tris.end());
    p.work.accepted.push_back(t);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol);
    reclaimRelational(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const Region* tr = torusOf(p.work);
    expect(tr != nullptr, "axis: torus kept");
    if (!tr) return;
    const gp_Ax1 ca(p.work.accepted[0].ax.Location(), p.work.accepted[0].ax.Direction());
    const gp_Lin cl(ca);
    expect(cl.Distance(tr->ax.Location()) < 1e-12 &&
               tr->ax.Direction().IsParallel(ca.Direction(), 1e-15),
           "axis: torus on the cone's line");
    expect(std::fabs(tr->ax.Location().Z() - 8.0) < 1e-9, "axis: centre height kept");
}

static void casePredicates() {
    Region a = atomOf({0, 1});
    expect(torusAtom(a), "atom predicate: non-closed CylGrow");
    a.closed360 = true;
    expect(!torusAtom(a), "atom predicate: closed is a donor, not an atom");
    a.closed360 = false;
    a.origin = Origin::FilletStrip;
    expect(!torusAtom(a), "atom predicate: FilletStrip is not an atom");
    Region t;
    t.type = SurfType::Torus;
    t.origin = Origin::TorusBlend;
    t.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    t.radius = 10.0;
    t.radius2 = 2.0;
    expect(std::fabs(torusSurfaceDev(t, gp_XYZ(12, 0, 0))) < 1e-12, "dev: outer equator");
    expect(std::fabs(torusSurfaceDev(t, gp_XYZ(0, 10, 2)) ) < 1e-12, "dev: top");
    expect(std::fabs(torusSurfaceDev(t, gp_XYZ(12.5, 0, 0)) - 0.5) < 1e-12, "dev: 0.5 out");
}

// D-fc-4.8 (b): points of a partial cylinder of radius R about (origin, axis),
// azimuth [u0, u1], heights [h0, h1], each pushed radially by a deterministic
// +-amp pattern (amp < tau keeps it certified).
static std::vector<gp_XYZ> cylPts(const gp_XYZ& o, const gp_XYZ& axisIn, double R, double u0,
                                  double u1, double h0, double h1, double amp) {
    gp_XYZ a = axisIn / axisIn.Modulus();
    const gp_XYZ tmp = std::fabs(a.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    gp_XYZ x = a.Crossed(tmp);
    x /= x.Modulus();
    const gp_XYZ y = a.Crossed(x);
    std::vector<gp_XYZ> out;
    int k = 0;
    for (int i = 0; i <= 8; ++i)
        for (int j = 0; j <= 2; ++j, ++k) {
            const double u = u0 + (u1 - u0) * i / 8.0, h = h0 + (h1 - h0) * j / 2.0;
            const double r = R + ((k % 2) ? amp : -amp);
            out.push_back(o + a * h + x * (r * std::cos(u)) + y * (r * std::sin(u)));
        }
    return out;
}

static void caseJointAxis() {
    const double tau = 1e-4;
    const gp_XYZ o(10, -5, 3), ax(0, 1, 0);
    std::vector<std::vector<gp_XYZ>> g = {cylPts(o, ax, 92.5, 0.2, 0.9, -3, 0, 0.3 * tau),
                                          cylPts(o, ax, 72.5, 1.0, 1.6, -1, 2, 0.3 * tau),
                                          cylPts(o, ax, 25.0, 0.0, 6.0, 0, 3, 0.3 * tau)};
    const JointAxisFit about = jointAxisAbout(g, o, ax);
    expect(about.radius.size() == 3 && std::fabs(about.radius[0] - 92.5) < 0.3 * tau &&
               std::fabs(about.radius[2] - 25.0) < 0.3 * tau,
           "joint about: one radius per donor");
    expect(jointAxisCertifies(about, tau), "joint about: true axis certifies");
    // A seed off by 2 tau in position and 1e-5 rad in tilt: the joint fit recovers it.
    JointAxisFit fit;
    const bool ok = jointAxisFit(g, o + gp_XYZ(2 * tau, 0, 0), gp_XYZ(1e-5, 1, 0), 16, fit);
    const JointAxisFit seed = jointAxisAbout(g, o + gp_XYZ(2 * tau, 0, 0), gp_XYZ(1e-5, 1, 0));
    expect(!jointAxisCertifies(seed, tau), "joint fit: the perturbed seed does not certify");
    expect(ok && jointAxisCertifies(fit, tau), "joint fit: refit certifies every donor");
    expect(ok && fit.axis.Y() > 0.0 && std::fabs(fit.axis.X()) < 1e-6, "joint fit: axis recovered");
    // One donor off the axis by 5 tau cannot certify jointly.
    std::vector<std::vector<gp_XYZ>> bad = g;
    bad.push_back(cylPts(o + gp_XYZ(5 * tau, 0, 0), ax, 40.0, 0.0, 6.2, 0, 3, 0.0));
    JointAxisFit fb;
    jointAxisFit(bad, o, ax, 16, fb);
    expect(!jointAxisCertifies(fb, tau), "joint fit: an off-axis donor breaks the certificate");
    expect(!jointAxisCertifies(about, 0.0), "joint certificate: tau 0 admits nothing");
}

static void caseCoaxialFamilies() {
    const double tau = 1e-4;
    const gp_XYZ o(0, 0, 0), ax(0, 1, 0);
    std::vector<AxisDonor> d(4);
    auto setD = [&](int i, const gp_XYZ& oi, const gp_XYZ& ai, double R, double u0, double u1) {
        d[(std::size_t)i].pts = cylPts(oi, ai, R, u0, u1, 0, 3, 0.2 * tau);
        d[(std::size_t)i].origin = oi;
        d[(std::size_t)i].axis = ai / ai.Modulus();
        d[(std::size_t)i].radius = R;
        d[(std::size_t)i].length = 3.0;
    };
    // Donor 0 leads with its own frame 1.5 tau off (its short arc still certifies
    // on it); donor 1's half circle does not, the joint LS does. Donor 2 is a
    // parallel bore elsewhere. Donor 3 is coaxial.
    setD(0, o + gp_XYZ(1.5 * tau, 0, 0), ax, 92.5, 0.2, 0.9);
    d[0].pts = cylPts(o, ax, 92.5, 0.2, 0.9, 0, 3, 0.2 * tau);
    setD(1, o, ax, 72.5, 0.0, 3.0);
    setD(2, o + gp_XYZ(30, 0, 0), ax, 5.0, 0.0, 6.2);
    setD(3, o, ax, 100.0, 4.0, 4.5);
    const std::vector<CoaxialFamily> f = coaxialFamilies(d, {0, 1, 2, 3}, tau, 16);
    expect(f.size() == 2, "families: coaxial three and one bore");
    expect(f.size() == 2 && f[0].members == std::vector<int>({0, 1, 3}) && f[0].refit,
           "families: leader 0 takes 1 and 3 by joint refit");
    expect(f.size() == 2 && f[1].members == std::vector<int>({2}) && !f[1].refit,
           "families: the parallel bore stays alone");
    expect(f.size() == 2 && jointAxisCertifies(f[0].fit, tau),
           "families: the family axis certifies");
    // A leader whose own frame certifies everyone keeps it (no refit).
    std::vector<AxisDonor> e = {d[1], d[3]};
    const std::vector<CoaxialFamily> g = coaxialFamilies(e, {0, 1}, tau, 16);
    expect(g.size() == 1 && !g[0].refit && g[0].members.size() == 2,
           "families: a certifying leader frame is kept");
    expect(coaxialFamilies(d, {0, 1, 2, 3}, 0.0, 16).empty(), "families: tau 0 forms none");
}

static void caseChainSuperset() {
    const double tau = 1e-4;
    // Profile points of a torus Rmaj 20, r 3 over 7 stations.
    std::vector<ProfPt> good;
    for (int i = 0; i < 7; ++i) {
        const double v = 0.2 + 0.2 * i;
        for (int k = 0; k < 3; ++k) {
            ProfPt p;
            p.rho = 20.0 + 3.0 * std::cos(v);
            p.z = 3.0 * std::sin(v);
            p.psi = 0.3 * k;
            good.push_back(p);
        }
    }
    std::vector<ProfPt> bad = good;
    bad[4].rho += 5.0 * tau;
    expect(!torusChainSupersetFails(good, 60, 30, tau),
           "superset: a certified chain set keeps the band");
    expect(torusChainSupersetFails(bad, 60, 30, tau),
           "superset: a refused chain set refuses its subset");
    expect(!torusChainSupersetFails(bad, 30, 30, tau),
           "superset: the band that is the whole set is not a subset");
}

// D-fc-4.9: points of a partial torus (Rmaj, r) about (origin, axis), azimuth
// [u0, u1] in nu steps, tube angle [v0, v1] in nv steps, each pushed along the
// tube normal by a deterministic +-amp pattern.
static std::vector<gp_XYZ> torusPts(const gp_XYZ& o, const gp_XYZ& axisIn, double Rmaj, double r,
                                    double u0, double u1, int nu, double v0, double v1, int nv,
                                    double amp) {
    gp_XYZ a = axisIn / axisIn.Modulus();
    const gp_XYZ tmp = std::fabs(a.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    gp_XYZ x = a.Crossed(tmp);
    x /= x.Modulus();
    const gp_XYZ y = a.Crossed(x);
    std::vector<gp_XYZ> out;
    int k = 0;
    for (int i = 0; i <= nu; ++i)
        for (int j = 0; j <= nv; ++j, ++k) {
            const double u = u0 + (u1 - u0) * i / nu, v = v0 + (v1 - v0) * j / nv;
            const double rr = r + ((k % 2) ? amp : -amp);
            const gp_XYZ radial = x * std::cos(u) + y * std::sin(u);
            out.push_back(o + radial * (Rmaj + rr * std::cos(v)) + a * (rr * std::sin(v)));
        }
    return out;
}

static void caseJointTorus() {
    const double tau = 1e-4;
    const gp_XYZ o(10, -5, 3), ax(0, 1, 0);
    // Short donor arcs about the true axis (the arc ribs and their walls) and a
    // partial fillet torus R 89.5 / r 3 on the same axis, 1 rad of azimuth.
    const std::vector<std::vector<gp_XYZ>> donors = {
        cylPts(o, ax, 92.5, 0.2, 0.25, 0, 3, 0.3 * tau),
        cylPts(o, ax, 72.5, 1.0, 1.05, -1, 2, 0.3 * tau),
        cylPts(o, ax, 100.0, 4.0, 4.04, 0, 3, 0.3 * tau)};
    const std::vector<gp_XYZ> chain =
        torusPts(o, ax, 89.5, 3.0, 0.1, 1.1, 20, 0.5 * kPi, kPi, 6, 0.3 * tau);
    // The seed axis is tilted by 4e-6 rad (0.36 mm at rho 89.5 ... 3.6 tau): the
    // chain does not certify on it, the donors' short arcs do not notice.
    const gp_XYZ seedAx(4e-6, 1, 0);
    const gp_XYZ xHint(1, 0, 0);
    const TorusCert onSeed = certifyTorus(profileAbout(chain, o, seedAx, xHint), tau);
    expect(!onSeed.ok, "joint torus: the chain does not certify on the tilted seed axis");
    expect(jointAxisCertifies(jointAxisAbout(donors, o, seedAx), tau),
           "joint torus: the donors alone certify the tilted seed axis");
    const JointTorusCert jc = jointTorusCertify(donors, chain, o, seedAx, xHint, tau, 16);
    expect(jc.fitted && jc.ok, "joint torus: donors + chain certify on the joint axis");
    expect(jc.ok && std::fabs(jc.chain.Rmaj - 89.5) < tau && std::fabs(jc.chain.Rmin - 3.0) < tau,
           "joint torus: the chain's torus radii are recovered");
    expect(jc.ok && std::fabs(jc.family.axis.X()) < 5e-7 && jc.family.axis.Y() > 0.0,
           "joint torus: the joint axis is the true axis, in the seed's hemisphere");
    expect(jc.ok && std::fabs(jc.tilt - 4e-6) < 5e-7, "joint torus: tilt reported vs the seed");
    expect(jc.family.radius.size() == 3 && jc.family.maxDev.size() == 3,
           "joint torus: one radius and deviation per donor");
    // The same torus on an axis tilted 1e-4 rad from the donors' certifies alone
    // (a band-only refit would take it); with the donors in the fit it cannot.
    const gp_XYZ tiltAx(1e-4, 1, 0);
    const std::vector<gp_XYZ> tilted =
        torusPts(o, tiltAx, 89.5, 3.0, 0.1, 1.1, 20, 0.5 * kPi, kPi, 6, 0.3 * tau);
    expect(certifyTorus(profileAbout(tilted, o, tiltAx, xHint), tau).ok,
           "joint torus: the tilted chain certifies on its own axis");
    const JointTorusCert jt = jointTorusCertify(donors, tilted, o, ax, xHint, tau, 16);
    expect(!jt.ok, "joint torus: never a band-only refit (donors pin the axis)");
    // A chain whose tube radius drifts by 4 tau is not one torus on any axis.
    std::vector<gp_XYZ> drift =
        torusPts(o, ax, 89.5, 3.0, 0.1, 0.6, 10, 0.5 * kPi, kPi, 6, 0.0);
    const std::vector<gp_XYZ> drift2 =
        torusPts(o, ax, 89.5, 3.0 + 4 * tau, 0.6, 1.1, 10, 0.5 * kPi, kPi, 6, 0.0);
    drift.insert(drift.end(), drift2.begin(), drift2.end());
    expect(!jointTorusCertify(donors, drift, o, ax, xHint, tau, 16).ok,
           "joint torus: a variable-radius chain is refused");
    expect(!jointTorusCertify(donors, chain, o, seedAx, xHint, 0.0, 16).ok,
           "joint torus: tau 0 admits nothing");
    expect(!jointTorusCertify({}, chain, o, seedAx, xHint, tau, 16).ok,
           "joint torus: no donor, no family certificate");
}

static void caseProfileAbout() {
    const gp_XYZ o(1, 2, 3), ax(0, 0, 2), xHint(1, 1, 0);
    const std::vector<ProfPt> p = profileAbout({gp_XYZ(1, 2, 3) + gp_XYZ(0, -2, 5)}, o, ax, xHint);
    const double expPsi = 2.0 * kPi - 0.75 * kPi;  // -135 deg about (1,1,0)/sqrt2, in [0, 2pi)
    expect(p.size() == 1 && std::fabs(p[0].z - 5.0) < 1e-12 && std::fabs(p[0].rho - 2.0) < 1e-12 &&
               std::fabs(p[0].psi - expPsi) < 1e-12,
           "profileAbout: rho, z, psi in the scan frame");
}

static void caseTorusSides() {
    expect(torusSidesOf(48, 2.0 * kPi, true) == 48, "torusSides: a closed band keeps its cols");
    expect(torusSidesOf(21, 1.0, false) == (int)std::lround(2.0 * kPi * 20.0),
           "torusSides: a partial band's steps over its own span");
    expect(torusSidesOf(1, 0.5, false) == 1, "torusSides: one column is one side");
    expect(torusSidesOf(9, 0.0, false) == 9, "torusSides: no span keeps cols");
}

static void caseWindowOverlap() {
    std::vector<char> failed(10, 0);
    failed[4] = failed[5] = 1;
    expect(torusWindowOverlapsRefused({3, 5, 7}, failed),
           "F1: a window sharing a triangle with a refused chain set is refused");
    expect(!torusWindowOverlapsRefused({0, 1, 2}, failed),
           "F1: a disjoint window is not refused");
    expect(!torusWindowOverlapsRefused({0, 12}, failed), "F1: out-of-range ids are ignored");
}

// D-fc-4.11 (F1 refined): which certifyTorus refusals are geometric.
static void caseRefusalGeometric() {
    const double tau = 1e-4;
    const gp_XYZ o(0, 0, 0), ax(0, 0, 1), xHint(1, 0, 0);
    const std::vector<gp_XYZ> good =
        torusPts(o, ax, 20.0, 3.0, 0.0, 1.0, 12, 0.5 * kPi, kPi, 6, 0.0);
    const TorusCert ok = certifyTorus(profileAbout(good, o, ax, xHint), tau);
    expect(ok.ok && !torusRefusalGeometric(ok), "F1 refined: a certified band is not a refusal");
    // Two tube stations: below certifyTorus's levels / profPts -> over-determination.
    const std::vector<gp_XYZ> thin =
        torusPts(o, ax, 20.0, 3.0, 0.0, 1.0, 12, 0.5 * kPi, 0.6 * kPi, 1, 0.0);
    const TorusCert under = certifyTorus(profileAbout(thin, o, ax, xHint), tau);
    expect(!under.ok && !under.determinate && !torusRefusalGeometric(under),
           "F1 refined: an under-determined band is not a geometric refusal");
    // One vertex 10 tau off the tube: vertex-off on determinate evidence.
    std::vector<gp_XYZ> off = good;
    off[off.size() / 2] += gp_XYZ(0, 0, 10.0 * tau);
    const TorusCert voff = certifyTorus(profileAbout(off, o, ax, xHint), tau);
    expect(!voff.ok && voff.determinate && torusRefusalGeometric(voff),
           "F1 refined: a vertex off tau is a geometric refusal");
    // A cone frustum's profile is a line: line-first on determinate evidence.
    std::vector<gp_XYZ> cone;
    for (int i = 0; i <= 12; ++i)
        for (int j = 0; j <= 6; ++j) {
            const double u = (double)i / 12.0, rr = 10.0 + 0.5 * j;
            cone.emplace_back(rr * std::cos(u), rr * std::sin(u), 0.5 * j);
        }
    const TorusCert line = certifyTorus(profileAbout(cone, o, ax, xHint), tau);
    expect(!line.ok && line.determinate && torusRefusalGeometric(line),
           "F1 refined: line-first on determinate evidence is a geometric refusal");
    TorusCert sing;
    sing.determinate = true;
    sing.why = "fit-singular";
    expect(!torusRefusalGeometric(sing),
           "F1 refined: a singular circle fit is not a geometric refusal");
}

// D-fc-4.11 (F1 refined): the chain-level verdict -- the joint certificate decides
// when it ran on a finite axis, the family certificate otherwise.
static void caseChainFailsGeometric() {
    const double tau = 1e-4;
    const gp_XYZ o(10, -5, 3), ax(0, 1, 0), xHint(1, 0, 0);
    const std::vector<std::vector<gp_XYZ>> donors = {
        cylPts(o, ax, 92.5, 0.2, 0.25, 0, 3, 0.3 * tau),
        cylPts(o, ax, 72.5, 1.0, 1.05, -1, 2, 0.3 * tau)};
    const std::vector<gp_XYZ> chain =
        torusPts(o, ax, 89.5, 3.0, 0.1, 1.1, 20, 0.5 * kPi, kPi, 6, 0.3 * tau);
    const gp_XYZ seedAx(4e-6, 1, 0);
    const TorusCert onSeed = certifyTorus(profileAbout(chain, o, seedAx, xHint), tau);
    const TorusCert onTrue = certifyTorus(profileAbout(chain, o, ax, xHint), tau);
    expect(onTrue.ok && !torusChainFailsGeometric(onTrue, nullptr, tau),
           "F1 refined: a chain the family axis certifies does not fail");
    expect(!onSeed.ok && onSeed.determinate && torusChainFailsGeometric(onSeed, nullptr, tau),
           "F1 refined: vertex-off on the family axis fails when no joint fit ran");
    const JointTorusCert jc = jointTorusCertify(donors, chain, o, seedAx, xHint, tau, 16);
    expect(jc.ok && !torusChainFailsGeometric(onSeed, &jc, tau),
           "F1 refined: a chain that certifies jointly does not fail (it may still be "
           "refused for over-determination, and then consumes nothing)");
    // Under-determined on the joint axis: the deviations of the donors about an
    // axis the chain cannot determine refuse nothing ...
    JointTorusCert und = jc;
    und.ok = false;
    und.chain = certifyTorus(
        profileAbout(torusPts(o, ax, 89.5, 3.0, 0.1, 1.1, 20, 0.5 * kPi, 0.6 * kPi, 1, 0.0), o, ax,
                     xHint),
        tau);
    for (double& m : und.family.maxDev) m = 1e3 * tau;
    const TorusCert thinSeed = certifyTorus(
        profileAbout(torusPts(o, ax, 89.5, 3.0, 0.1, 1.1, 20, 0.5 * kPi, 0.6 * kPi, 1, 0.0), o, ax,
                     xHint),
        tau);
    expect(!und.chain.determinate && !thinSeed.determinate &&
               !torusChainFailsGeometric(thinSeed, &und, tau),
           "F1 refined: a chain under-determined on both axes does not fail");
    // ... but D-fc-4.12: a refusal on a determinate family-axis certificate is
    // geometric even when the joint fit is under-determined.
    expect(!und.chain.determinate && torusChainFailsGeometric(onSeed, &und, tau),
           "D-fc-4.12: vertex-off on the family axis fails though the joint fit is "
           "under-determined");
    // Determinate on the joint axis but a donor vertex off tau: geometric.
    JointTorusCert don = jc;
    don.ok = false;
    don.family.maxDev.front() = 2.0 * tau;
    expect(torusChainFailsGeometric(onSeed, &don, tau),
           "F1 refined: a donor vertex off tau under the joint axis fails");
    // A tube radius drifting by 4 tau is vertex-off on the joint axis: geometric.
    std::vector<gp_XYZ> drift = torusPts(o, ax, 89.5, 3.0, 0.1, 0.6, 10, 0.5 * kPi, kPi, 6, 0.0);
    const std::vector<gp_XYZ> drift2 =
        torusPts(o, ax, 89.5, 3.0 + 4 * tau, 0.6, 1.1, 10, 0.5 * kPi, kPi, 6, 0.0);
    drift.insert(drift.end(), drift2.begin(), drift2.end());
    const JointTorusCert jd = jointTorusCertify(donors, drift, o, seedAx, xHint, tau, 16);
    const TorusCert dSeed = certifyTorus(profileAbout(drift, o, seedAx, xHint), tau);
    expect(jd.fitted && !jd.ok && jd.chain.determinate && torusChainFailsGeometric(dSeed, &jd, tau),
           "F1 refined: a variable-radius chain fails geometrically on the joint axis");
    JointTorusCert unfit;
    expect(torusChainFailsGeometric(onSeed, &unfit, tau),
           "F1 refined: without a finite joint axis the family verdict decides");
}

// D-fc-4.12 (audit-fx-tori-r3 F-1): a ring chain that bridges two surfaces --
// the base-plane annulus between a rho 72.5 rim and the fillet's contact ring at
// rho 89.5, plus the fillet's next row -- has three profile stations, so no
// axis over-determines it, yet on a seed axis tilted by the donors' residual the
// stations split and certifyTorus refuses it vertex-off on determinate counts.
// The joint fit cannot determine it; under the strict rule it is a geometric
// failure (it consumes its edges and is a failing superset).
static void caseBridgeChain() {
    const double tau = 1e-4;
    const gp_XYZ o(10, -5, 3), ax(0, 1, 0), xHint(1, 0, 0);
    const std::vector<std::vector<gp_XYZ>> donors = {
        cylPts(o, ax, 92.5, 0.2, 0.25, 0, 3, 0.3 * tau),
        cylPts(o, ax, 72.5, 1.0, 1.05, -1, 2, 0.3 * tau)};
    const gp_XYZ y = ax, x = ax.Crossed(gp_XYZ(1, 0, 0)) / ax.Crossed(gp_XYZ(1, 0, 0)).Modulus();
    const gp_XYZ e2 = ax.Crossed(x);
    std::vector<gp_XYZ> bridge;
    const double vNext = 75.0 * kPi / 180.0;
    const std::array<std::pair<double, double>, 3> stations = {
        std::make_pair(72.5, 4.5), std::make_pair(89.5, 4.5),
        std::make_pair(89.5 + 3.0 * std::cos(vNext), 1.5 + 3.0 * std::sin(vNext))};
    for (const auto& st : stations)
        for (int i = 0; i <= 20; ++i) {
            const double u = 0.1 + 1.0 * i / 20.0;
            bridge.push_back(o + y * st.second + x * (st.first * std::cos(u)) +
                             e2 * (st.first * std::sin(u)));
        }
    const TorusCert onTrue = certifyTorus(profileAbout(bridge, o, ax, xHint), tau);
    expect(!onTrue.ok && !onTrue.determinate,
           "bridge: three stations are under-determined on the true axis");
    const gp_XYZ seedAx(8e-6, 1, 0);
    const TorusCert onSeed = certifyTorus(profileAbout(bridge, o, seedAx, xHint), tau);
    expect(!onSeed.ok && onSeed.determinate && torusRefusalGeometric(onSeed),
           "bridge: on the tilted family axis the split stations refuse vertex-off, "
           "determinately");
    const JointTorusCert jb = jointTorusCertify(donors, bridge, o, seedAx, xHint, tau, 16);
    expect(jb.fitted && !jb.ok, "bridge: the joint certificate refuses it");
    expect(torusChainFailsGeometric(onSeed, &jb, tau),
           "D-fc-4.12: the bridge chain fails geometrically (consumes, failing superset)");
    JointTorusCert jbu = jb;
    jbu.chain = onTrue;
    expect(!jbu.chain.determinate && torusChainFailsGeometric(onSeed, &jbu, tau),
           "D-fc-4.12: ... also when the joint fit finds it under-determined");
    expect(!torusChainFailsGeometric(onTrue, &jbu, tau),
           "D-fc-4.12: under-determined on the family axis and jointly is not geometric");
}

// D-fc-4.11 (a): a band holding a withheld triangle is gathered but not emitted;
// its triangles and its atom keep their pre-T state.
static void caseWithheld() {
    for (int variant = 0; variant < 2; ++variant) {
        Part p = ringTorus(24, 6);
        Provisional pr;
        for (int j = 0; j < 4; ++j)
            pr.tris.insert(pr.tris.end(), p.rows[(std::size_t)j].begin(),
                           p.rows[(std::size_t)j].end());
        std::sort(pr.tris.begin(), pr.tris.end());
        p.work.provisionals.push_back(pr);
        std::vector<int> at;
        for (int j = 4; j < 6; ++j)
            at.insert(at.end(), p.rows[(std::size_t)j].begin(), p.rows[(std::size_t)j].end());
        p.work.accepted.push_back(atomOf(at));
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol);
        // variant 0: a triangle of the band; variant 1: a donor triangle (not a band).
        sp.torusWithheldTris = {variant == 0 ? p.rows[2][5] : p.donorTris[3]};
        expect(claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus), "withheld call");
        if (variant == 0) {
            expect(torusOf(p.work) == nullptr, "withheld: the band is not emitted");
            expect(countOrigin(p.work, Origin::CylGrow, false) == 1,
                   "withheld: its atom stays accepted");
            expect(p.work.provisionals[0].claim == ProvClaim::Unclaimed,
                   "withheld: its free triangles stay unclaimed (pre-T)");
            expect(!p.work.torusAdmitted && p.work.torusPending.empty(),
                   "withheld: no torus admitted, no pending pack");
        } else {
            expect(torusOf(p.work) != nullptr, "withheld: a band without the triangle still ships");
        }
    }
}

static bool sameRegionSet(const RegionSet& a, const RegionSet& b) {
    if (a.regions.size() != b.regions.size() || a.chains.size() != b.chains.size() ||
        a.triRegion != b.triRegion || a.triIsland != b.triIsland ||
        a.torusRevertValid != b.torusRevertValid ||
        a.torusRevert.regions.size() != b.torusRevert.regions.size() ||
        a.torusRevert.triRegion != b.torusRevert.triRegion)
        return false;
    for (std::size_t i = 0; i < a.regions.size(); ++i) {
        const Region& x = a.regions[i];
        const Region& y = b.regions[i];
        if (x.type != y.type || x.origin != y.origin || x.tris != y.tris || x.radius != y.radius ||
            x.radius2 != y.radius2 || x.closed360 != y.closed360 ||
            x.dVolPredicted != y.dVolPredicted)
            return false;
        // Bit-for-bit (gp_Dir::IsEqual measures an angle, never exactly 0).
        const gp_XYZ lx = x.ax.Location().XYZ(), ly = y.ax.Location().XYZ();
        const gp_XYZ dx = x.ax.Direction().XYZ(), dy = y.ax.Direction().XYZ();
        if (lx.X() != ly.X() || lx.Y() != ly.Y() || lx.Z() != ly.Z() || dx.X() != dy.X() ||
            dx.Y() != dy.Y() || dx.Z() != dy.Z())
            return false;
    }
    for (std::size_t i = 0; i < a.chains.size(); ++i)
        if (a.chains[i].regA != b.chains[i].regA || a.chains[i].regB != b.chains[i].regB ||
            a.chains[i].meshVerts != b.chains[i].meshVerts)
            return false;
    return true;
}

// D-fc-4.11 (a): resegmentTorus (the stages from T on, from segment()'s pre-T
// state) equals a full segment() with the same withheld triangles. On the corpus
// pickup the two R12 hub tori are withheld; the two R89.5 arc tori stay.
static void caseResegment(const char* stl) {
    stl2step::harness::HarnessMesh mesh;
    std::string err;
    if (!stl2step::harness::loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        expect(false, "reseg: corpus mesh loads");
        return;
    }
    MeshView mv{};
    stl2step::harness::toMeshView(mesh, mesh.comps[0], mv);
    SegmentParams sp;
    RegionSet probe;
    expect(segment(mv, sp, probe, nullptr), "reseg: segment");
    expect(probe.torusRevertValid && probe.preTorusWork != nullptr,
           "reseg: a torus was admitted and the pre-T state is kept");
    std::vector<int> held;
    int nArc = 0;
    for (const Region& r : probe.regions) {
        if (r.origin != Origin::TorusBlend) continue;
        if (r.radius < 20.0) held.insert(held.end(), r.tris.begin(), r.tris.end());
        else ++nArc;
    }
    std::sort(held.begin(), held.end());
    expect(!held.empty() && nArc >= 1, "reseg: hub and arc tori claimed");
    SegmentParams hp = sp;
    hp.torusWithheldTris = held;
    RegionSet full, tail;
    expect(segment(mv, hp, full, nullptr), "reseg: segment with withheld tori");
    expect(resegmentTorus(mv, hp, probe, tail, nullptr), "reseg: resegmentTorus");
    expect(sameRegionSet(full, tail), "reseg: resegmentTorus == segment, region by region");
    int nTorus = 0, nHeldTorus = 0;
    for (const Region& r : tail.regions) {
        if (r.origin != Origin::TorusBlend) continue;
        ++nTorus;
        for (int t : r.tris)
            if (std::binary_search(held.begin(), held.end(), t)) {
                ++nHeldTorus;
                break;
            }
    }
    expect(nTorus == nArc && nHeldTorus == 0, "reseg: only the withheld tori are gone");
    RegionSet none;
    expect(!resegmentTorus(mv, hp, RegionSet{}, none, nullptr),
           "reseg: no pre-T state, no resegment");
}

int main(int argc, char** argv) {
    if (argc > 1) caseResegment(argv[1]);
    caseRefusalGeometric();
    caseChainFailsGeometric();
    caseBridgeChain();
    caseWithheld();
    caseWindowOverlap();
    caseJointTorus();
    caseProfileAbout();
    caseTorusSides();
    caseJointAxis();
    caseCoaxialFamilies();
    caseChainSuperset();
    casePredicates();
    caseLeaveOneOutPure();
    caseRails();
    caseReclaimAbsorb();
    caseAxisLine();
    caseAtomJoins();
    caseAtomOffSurface();
    caseLeaveOneOut();
    if (gFail) std::fprintf(stderr, "relational_torus_unit: %d failure(s)\n", gFail);
    return gFail ? 1 : 0;
}
