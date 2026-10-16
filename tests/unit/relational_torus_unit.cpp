// Synthetic torus-relation cases (D-fc-4.6 F1/F2/F4). No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_relational.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
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

int main() {
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
