// Synthetic relational-revolve cases (D-fc-1.6, D-fc-1.7). No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#define rr_dup _dup
#define rr_dup2 _dup2
#define rr_fileno _fileno
#define rr_close _close
#else
// dup/dup2/close are libc. Declared here so this file does not include <unistd.h>.
// fileno comes from <cstdio>.
extern "C" {
int dup(int);
int dup2(int, int);
int close(int);
}
#define rr_dup dup
#define rr_dup2 dup2
#define rr_fileno fileno
#define rr_close close
#endif

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    if (ok) {
        std::fprintf(stderr, "PASS %s\n", name);
        return;
    }
    std::fprintf(stderr, "FAIL %s\n", name);
    ++gFail;
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
    void addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
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

static void band(Synth& s, const std::vector<int>& a, const std::vector<int>& b, int sectors,
                 std::vector<int>& tris) {
    const int N = (int)a.size();
    for (int i = 0; i < sectors; ++i) {
        const int i1 = (i + 1) % N;
        const int t0 = (int)(s.triFlat.size() / 3);
        s.addTri(a[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i]);
        s.addTri(b[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i1]);
        tris.push_back(t0);
        tris.push_back(t0 + 1);
    }
}

static std::vector<int> ring(Synth& s, double R, double z, int N) {
    std::vector<int> ix;
    for (int i = 0; i < N; ++i) {
        const double th = kPi * 2.0 * (double)i / (double)N;
        ix.push_back(s.add(R * std::cos(th), R * std::sin(th), z));
    }
    return ix;
}

struct Part {
    Synth mesh;
    SegmentWork work;
    std::vector<int> feature;
};

static Part frustum(int N, int stations, int sectors, double Rlo, double zLo, double Rhi,
                    double zHi, double cylH, bool perturbDonor) {
    Part p;
    const double tau = 2.0 * p.mesh.q;
    std::vector<std::vector<int>> rings;
    rings.push_back(ring(p.mesh, Rlo, zLo - cylH, N));
    rings.push_back(ring(p.mesh, Rlo, zLo, N));
    for (int s = 1; s < stations; ++s) {
        const double t = (double)s / (double)(stations - 1);
        const double R = Rlo + t * (Rhi - Rlo);
        const double z = zLo + t * (zHi - zLo);
        rings.push_back(ring(p.mesh, R, z, N));
    }
    if (perturbDonor) {
        for (std::size_t i = 0; i < rings[0].size(); i += 2) {
            gp_XYZ& v = p.mesh.pts[(std::size_t)rings[0][(std::size_t)i]];
            const double m = std::hypot(v.X(), v.Y());
            if (m > 0.0) {
                v.SetX(v.X() * (1.0 + 50.0 * tau / m));
                v.SetY(v.Y() * (1.0 + 50.0 * tau / m));
            }
        }
    }
    std::vector<int> donorTris;
    band(p.mesh, rings[0], rings[1], N, donorTris);
    for (int s = 1; s + 1 < (int)rings.size(); ++s)
        band(p.mesh, rings[(std::size_t)s], rings[(std::size_t)s + 1], sectors, p.feature);
    p.mesh.finish();
    Region d;
    d.id = 0;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, zLo - 0.5 * cylH), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = Rlo;
    d.closed360 = true;
    d.tris = donorTris;
    d.maxVertexDev = 0;
    p.work.accepted.push_back(d);
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = p.feature;
    p.work.provisionals.push_back(pr);
    return p;
}

static void params(SegmentParams& p, DerivedTols& tol, int threads) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    p.requestedThreads = threads;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kPi / 180.0;
}

static int countType(const SegmentWork& w, SurfType t, Origin o) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.type == t && r.origin == o) ++n;
    return n;
}

static const Region* findNew(const SegmentWork& w, SurfType t, Origin o) {
    for (const Region& r : w.accepted)
        if (r.type == t && r.origin == o) return &r;
    return nullptr;
}

struct Grab {
    int saved = -1;
    bool on = false;
    Grab() {
        ::setenv("STL2STEP_DIAG_REL", "1", 1);
        saved = rr_dup(rr_fileno(stderr));
        FILE* f = std::fopen("rel-revolve-diag.txt", "w");
        if (!f) return;
        rr_dup2(rr_fileno(f), rr_fileno(stderr));
        std::fclose(f);
        on = true;
    }
    std::string stop() {
        std::fflush(stderr);
        if (on && saved >= 0) rr_dup2(saved, rr_fileno(stderr));
        if (saved >= 0) rr_close(saved);
        saved = -1;
        std::ifstream in("rel-revolve-diag.txt");
        std::stringstream ss;
        ss << in.rdbuf();
        std::remove("rel-revolve-diag.txt");
        return ss.str();
    }
};

static bool hasWhy(const std::string& diag, const char* why) {
    const std::string key = std::string("why=") + why;
    return diag.find(key) != std::string::npos;
}

static void checkRim(const Region& cone, const Synth& mesh, const std::vector<int>& feat,
                     double tau, const char* name) {
    const double h = cone.maxVertexSnap;
    const double Rlo = cone.radius;
    const double Rhi = cone.vMax;
    const gp_XYZ loc = cone.ax.Location().XYZ();
    const gp_XYZ ax = cone.ax.Direction().XYZ();
    double worst = 0;
    std::vector<int> seen;
    for (int t : feat) {
        for (int k = 0; k < 3; ++k) seen.push_back(mesh.triFlat[(std::size_t)(3 * t + k)]);
    }
    std::sort(seen.begin(), seen.end());
    seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
    for (int g : seen) {
        const gp_XYZ p = mesh.pts[(std::size_t)g];
        const gp_XYZ d = p - loc;
        const double v = d.Dot(ax);
        const double rho = (d - ax * v).Modulus();
        const double t = (std::fabs(h) > 0.0) ? v / h : 0.0;
        const double R = Rlo + t * (Rhi - Rlo);
        worst = std::max(worst, std::fabs(rho - R));
        worst = std::max(worst, std::fabs(v - t * h));
    }
    std::fprintf(stderr, "rim %s worst=%.3g tau=%.3g Rlo=%.6g Rhi=%.6g h=%.6g\n", name, worst,
                 tau, Rlo, Rhi, h);
    expect(worst <= tau, name);
}

static void caseClosed45() {
    Part p = frustum(16, 2, 16, 10.0, 0.0, 14.0, 4.0, 4.0, false);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    expect(claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line), "45 call");
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr, "45 cone claimed");
    if (!c) return;
    expect(c->closed360, "45 closed360");
    expect(std::fabs(c->radius - 10.0) <= tauOf(mv), "45 Rlo");
    expect(std::fabs(c->vMax - 14.0) <= tauOf(mv), "45 Rhi");
    expect(std::fabs(std::fabs(c->maxVertexSnap) - 4.0) <= tauOf(mv), "45 height");
    const double slope = (c->vMax - c->radius) / std::fabs(c->maxVertexSnap);
    const double ang = std::atan(slope) * 180.0 / kPi;
    const double tolDeg = (tauOf(mv) / 4.0) * 180.0 / kPi;
    std::fprintf(stderr, "45 halfAng=%.6f tolDeg=%.6g\n", ang, tolDeg);
    expect(std::fabs(ang - 45.0) <= std::max(tolDeg, 1e-6), "45 half-angle");
    expect(p.work.provisionals[0].claim == ProvClaim::ConsumedCylinder, "45 provisional");
    checkRim(*c, p.mesh, p.feature, tauOf(mv), "45 rims within tau");
    expect(countType(p.work, SurfType::Torus, Origin::TorusBlend) == 0, "45 no torus");
}

static void caseFunnel() {
    const double slope = std::tan(21.1247 * kPi / 180.0);
    Part p = frustum(24, 6, 24, 12.0, 0.0, 12.0 + slope * 10.0, 10.0, 6.0, false);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr, "21 cone claimed");
    if (!c) return;
    const double got = (c->vMax - c->radius) / std::fabs(c->maxVertexSnap);
    const double ang = std::atan(got) * 180.0 / kPi;
    const double tolDeg = (tauOf(mv) / 10.0) * 180.0 / kPi;
    std::fprintf(stderr, "21 halfAng=%.6f tolDeg=%.6g nTri=%d\n", ang, tolDeg, (int)c->tris.size());
    expect(c->closed360, "21 closed");
    expect(std::fabs(ang - 21.1247) <= std::max(tolDeg, 1e-4), "21 half-angle");
    expect(c->maxVertexDev <= tauOf(mv), "21 vertices within tau");
}

static Part torusPart() {
    Part p;
    const int N = 16;
    const int nPhi = 6;
    const double Rmaj = 12.0, r = 2.0, z0 = 8.0;
    const double phi0 = kPi;
    const double phi1 = 0.5 * kPi;
    std::vector<std::vector<int>> rings;
    for (int j = 0; j <= nPhi; ++j) {
        const double t = (double)j / (double)nPhi;
        const double phi = phi0 + t * (phi1 - phi0);
        const double rho = Rmaj + r * std::cos(phi);
        const double z = z0 + r * std::sin(phi);
        rings.push_back(ring(p.mesh, rho, z, N));
    }
    const double rhoC = Rmaj + r * std::cos(phi0);
    const double zC = z0 + r * std::sin(phi0);
    auto lower = ring(p.mesh, rhoC, zC - 6.0, N);
    std::vector<int> donorTris;
    band(p.mesh, lower, rings[0], N, donorTris);
    for (int j = 0; j < nPhi; ++j)
        band(p.mesh, rings[(std::size_t)j], rings[(std::size_t)j + 1], N, p.feature);
    p.mesh.finish();
    Region d;
    d.id = 0;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, zC - 3.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = rhoC;
    d.closed360 = true;
    d.tris = donorTris;
    p.work.accepted.push_back(d);
    Provisional pr;
    pr.tris = p.feature;
    p.work.provisionals.push_back(pr);
    return p;
}

static void caseQuarterTorus() {
    Part line = torusPart();
    MeshView mv = line.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    claimRelRevolve(mv, sp, tol, line.work, RelRevolveMode::Line);
    expect(countType(line.work, SurfType::Cone, Origin::ChamferCone) == 0, "torus not a cone");
    Part tor = torusPart();
    MeshView mv2 = tor.mesh.view();
    claimRelRevolve(mv2, sp, tol, tor.work, RelRevolveMode::Torus);
    const Region* t = findNew(tor.work, SurfType::Torus, Origin::TorusBlend);
    expect(t != nullptr, "torus claimed in Torus mode");
    if (!t) return;
    std::fprintf(stderr, "torus Rmaj=%.6g Rmin=%.6g dev=%.3g tau=%.3g pending=%d admitted=%d\n",
                 t->radius, t->radius2, t->maxVertexDev, tauOf(mv2), t->torusPending,
                 tor.work.torusAdmitted ? 1 : 0);
    expect(std::fabs(t->radius - 12.0) <= tauOf(mv2), "torus Rmaj");
    expect(std::fabs(t->radius2 - 2.0) <= tauOf(mv2), "torus Rmin");
    expect(t->maxVertexDev <= tauOf(mv2), "torus vertices");
    expect(tor.work.torusAdmitted, "torusAdmitted");
    expect(t->torusPending >= 0, "torusPending index");
    expect(tor.work.provisionals[0].claim == ProvClaim::PendingTorus, "PendingTorus");
}

static void caseTwoLevelNotTorus() {
    Part p = frustum(16, 2, 16, 8.0, 0.0, 11.0, 3.0, 3.0, false);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 1, "two-level is a cone");
    Part q = frustum(16, 2, 16, 8.0, 0.0, 11.0, 3.0, 3.0, false);
    MeshView mv2 = q.mesh.view();
    claimRelRevolve(mv2, sp, tol, q.work, RelRevolveMode::Torus);
    expect(countType(q.work, SurfType::Torus, Origin::TorusBlend) == 0, "two-level is not a torus");
}

static Part lathe(double Rmaj, double Rmin, double phi0, double phi1, int nPhi, int N,
                  double cylR, double cylZ0, double cylZ1) {
    Part p;
    std::vector<std::vector<int>> rings;
    for (int j = 0; j <= nPhi; ++j) {
        const double t = (double)j / (double)nPhi;
        const double phi = phi0 + t * (phi1 - phi0);
        const double rho = Rmaj + Rmin * std::cos(phi);
        const double z = Rmin * std::sin(phi);
        rings.push_back(ring(p.mesh, rho, z, N));
    }
    int join = 0;
    double best = 1e300;
    for (int j = 0; j <= nPhi; ++j) {
        const double t = (double)j / (double)nPhi;
        const double phi = phi0 + t * (phi1 - phi0);
        const double rho = Rmaj + Rmin * std::cos(phi);
        const double d = std::fabs(rho - cylR);
        if (d < best) {
            best = d;
            join = j;
        }
    }
    auto lower = ring(p.mesh, cylR, cylZ0, N);
    std::vector<int> donorTris;
    band(p.mesh, lower, rings[(std::size_t)join], N, donorTris);
    for (int j = 0; j < nPhi; ++j) {
        if (j == join || j + 1 == join) continue;
        band(p.mesh, rings[(std::size_t)j], rings[(std::size_t)j + 1], N, p.feature);
    }
    // The join ring is the cylinder rim. Bands on either side that touch it
    // stay in the feature so the skin is edge-adjacent to the donor.
    if (join > 0)
        band(p.mesh, rings[(std::size_t)join - 1], rings[(std::size_t)join], N, p.feature);
    if (join < nPhi)
        band(p.mesh, rings[(std::size_t)join], rings[(std::size_t)join + 1], N, p.feature);
    (void)cylZ1;
    p.mesh.finish();
    Region d;
    d.id = 0;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, 0.5 * (cylZ0 + cylZ1)), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = cylR;
    d.closed360 = true;
    d.tris = donorTris;
    p.work.accepted.push_back(d);
    Provisional pr;
    pr.tris = p.feature;
    p.work.provisionals.push_back(pr);
    return p;
}

static void caseSphere() {
    // Profile circle centred on the axis: Rmaj = 0, Rmin = 5.
    const double phi0 = 0.4;
    const double rho0 = 5.0 * std::cos(phi0);
    const double zJoin = 5.0 * std::sin(phi0);
    Part p = lathe(0.0, 5.0, phi0, 1.2, 5, 16, rho0, zJoin - 4.0, zJoin);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const std::string diag = g.stop();
    expect(countType(p.work, SurfType::Torus, Origin::TorusBlend) == 0, "sphere not claimed");
    expect(hasWhy(diag, "sphere"), "sphere why");
    if (!hasWhy(diag, "sphere")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
}

static void caseSpindle() {
    Part p = lathe(3.0, 5.0, 0.0, 0.8, 5, 16, 8.0, 0.0, 4.0);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const std::string diag = g.stop();
    expect(countType(p.work, SurfType::Torus, Origin::TorusBlend) == 0, "spindle not claimed");
    expect(hasWhy(diag, "spindle"), "spindle why");
    if (!hasWhy(diag, "spindle")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
}

// D-fc-5.17. A quarter whose fitted major radius sits half a tau inside the
// minor radius is a spindle to certifyTorus. The donor cylinder is the outer
// equator, so the horn that ships has Rmaj = r exactly.
static void caseHorn() {
    const double tau = 2.0 * 1e-6;
    const double r = 2.0;
    const double Rmaj = r - 0.5 * tau;
    Part p = lathe(Rmaj, r, 0.0, 0.5 * kPi, 5, 16, Rmaj + r, -4.0, 0.0);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    expect(tauOf(mv) == tau, "horn tau is 2q");
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const Region* t = findNew(p.work, SurfType::Torus, Origin::TorusBlend);
    expect(t != nullptr, "horn claimed");
    if (!t) return;
    std::fprintf(stderr, "horn Rmaj=%.10g Rmin=%.10g dev=%.3g tau=%.3g\n", t->radius, t->radius2,
                 t->maxVertexDev, tauOf(mv));
    expect(t->radius == t->radius2, "horn Rmaj = r");
    expect(t->radius > tauOf(mv), "horn radius above tau");
    expect(t->maxVertexDev <= tauOf(mv), "horn vertices within tau");
}

static void casePartial() {
    Part p = frustum(16, 2, 4, 10.0, 0.0, 14.0, 4.0, 4.0, false);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const std::string diag = g.stop();
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0, "partial not claimed");
    expect(hasWhy(diag, "partial-cone"), "partial-cone why");
    if (!hasWhy(diag, "partial-cone")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
}

static void caseBadDonor() {
    Part p = frustum(16, 2, 16, 10.0, 0.0, 14.0, 4.0, 4.0, true);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const std::string diag = g.stop();
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0, "bad donor claims nothing");
    expect(hasWhy(diag, "no-donor"), "no-donor why");
    if (!hasWhy(diag, "no-donor")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
}

static void caseWedge() {
    Part p;
    const int N = 8;
    auto bot = ring(p.mesh, 5.0, 0.0, N);
    auto top = ring(p.mesh, 5.0, 6.0, N);
    std::vector<int> donorTris, cap;
    band(p.mesh, bot, top, N, donorTris);
    // Elliptical cap shares each rim edge with the donor wall and nothing else.
    // A center fan on that edge used to fill both adjacency slots, so the cap
    // was never scored.
    std::vector<int> ell;
    for (int i = 0; i < N; ++i) {
        const double th = kPi * 2.0 * (double)i / (double)N;
        const double rho = (4.0 * 7.0) / std::sqrt(std::pow(7.0 * std::cos(th), 2.0) +
                                                   std::pow(4.0 * std::sin(th), 2.0));
        // Stagger z so no three ellipse vertices share (ρ, z). An axis-aligned
        // ellipse repeats ρ four times; those vertices would be a ring.
        ell.push_back(p.mesh.add(rho * std::cos(th), rho * std::sin(th), 6.4 + 0.05 * (double)i));
    }
    band(p.mesh, top, ell, N, cap);
    p.mesh.finish();
    Region d;
    d.id = 0;
    d.type = SurfType::Cylinder;
    d.origin = Origin::NgonWall;
    d.ax = gp_Ax3(gp_Pnt(0, 0, 3), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = 5.0;
    d.closed360 = true;
    d.nSides = 8;
    d.tris = donorTris;
    p.work.accepted.push_back(d);
    Provisional pr;
    pr.tris = cap;
    p.work.provisionals.push_back(pr);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Torus);
    const std::string diag = g.stop();
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0, "wedge claims no cone");
    expect(countType(p.work, SurfType::Torus, Origin::TorusBlend) == 0, "wedge claims no torus");
    expect(hasWhy(diag, "vertex-off"), "wedge cap scored");
    expect(diag.find("nTri=16") != std::string::npos, "wedge cap triangle count");
    if (!hasWhy(diag, "vertex-off")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
}

static void addCyl(Part& p, std::vector<int> tris, double R, double zMid) {
    Region d;
    d.id = 0;
    d.type = SurfType::Cylinder;
    d.origin = Origin::CylGrow;
    d.ax = gp_Ax3(gp_Pnt(0, 0, zMid), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    d.radius = R;
    d.closed360 = true;
    d.tris = std::move(tris);
    d.maxVertexDev = 0;
    p.work.accepted.push_back(std::move(d));
}

static void caseProfileRow() {
    const char* names[] = {"unclaimed", "pending-torus", "consumed-cylinder"};
    const ProvClaim claims[] = {ProvClaim::Unclaimed, ProvClaim::PendingTorus,
                                ProvClaim::ConsumedCylinder};
    for (int run = 0; run < 3; ++run) {
        Part p;
        const int N = 16;
        const int nRing = 6;
        const double Rmaj = 20.0, rad = 3.0;
        const double phi0 = kPi, phi1 = 0.5 * kPi;
        std::vector<std::vector<int>> rings;
        for (int j = 0; j < nRing; ++j) {
            const double t = (double)j / (double)(nRing - 1);
            const double phi = phi0 + t * (phi1 - phi0);
            rings.push_back(ring(p.mesh, Rmaj + rad * std::cos(phi), rad * std::sin(phi), N));
        }
        const double rho0 = Rmaj + rad * std::cos(phi0);
        const double z0 = rad * std::sin(phi0);
        auto lower = ring(p.mesh, rho0, z0 - 5.0, N);
        std::vector<int> donorTris, mid, nbr;
        band(p.mesh, lower, rings[0], N, donorTris);
        for (int j = 0; j + 1 < nRing; ++j) {
            std::vector<int>& dest = (j == 2) ? mid : nbr;
            band(p.mesh, rings[(std::size_t)j], rings[(std::size_t)j + 1], N, dest);
        }
        p.mesh.finish();
        addCyl(p, donorTris, rho0, z0 - 2.5);
        Provisional a;
        a.claim = ProvClaim::Unclaimed;
        a.tris = mid;
        p.work.provisionals.push_back(a);
        Provisional b;
        b.claim = claims[run];
        b.tris = nbr;
        p.work.provisionals.push_back(b);
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol, 1);
        Grab g;
        claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
        const std::string diag = g.stop();
        const std::string label = std::string("profile-row ") + names[run];
        expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0, label.c_str());
        expect(hasWhy(diag, "profile-row"), label.c_str());
        if (!hasWhy(diag, "profile-row")) std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
    }
}

static void caseCreaseFrustum() {
    Part p;
    const int N = 16;
    auto b0 = ring(p.mesh, 10.0, -6.0, N);
    auto b1 = ring(p.mesh, 10.0, -3.0, N);
    auto b2 = ring(p.mesh, 10.0, 0.0, N);
    auto t0 = ring(p.mesh, 14.0, 4.0, N);
    auto t1 = ring(p.mesh, 14.0, 7.0, N);
    auto t2 = ring(p.mesh, 14.0, 10.0, N);
    std::vector<int> lower, fr, upper;
    band(p.mesh, b0, b1, N, lower);
    band(p.mesh, b1, b2, N, lower);
    band(p.mesh, b2, t0, N, fr);
    band(p.mesh, t0, t1, N, upper);
    band(p.mesh, t1, t2, N, upper);
    p.mesh.finish();
    addCyl(p, lower, 10.0, -3.0);
    addCyl(p, upper, 14.0, 7.0);
    Provisional pr;
    pr.tris = fr;
    p.work.provisionals.push_back(pr);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const std::string diag = g.stop();
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr, "crease cone");
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 1, "crease one cone");
    if (!c) {
        std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
        return;
    }
    expect(c->closed360, "crease closed");
    expect(std::fabs(c->radius - 10.0) <= tauOf(mv), "crease Rlo");
    expect(std::fabs(c->vMax - 14.0) <= tauOf(mv), "crease Rhi");
    const double ang = std::atan((c->vMax - c->radius) / std::fabs(c->maxVertexSnap)) * 180.0 / kPi;
    expect(std::fabs(ang - 45.0) <= 1e-4, "crease 45");
    expect(!hasWhy(diag, "profile-row"), "crease is not a profile row");
}

static void caseGapFunnel() {
    Part p;
    const int N = 16;
    auto d0 = ring(p.mesh, 10.0, -4.0, N);
    auto d1 = ring(p.mesh, 10.0, 0.0, N);
    auto g1 = ring(p.mesh, 12.0, 2.0, N);
    auto f1 = ring(p.mesh, 18.0, 10.0, N);
    std::vector<int> donor, gap, fun;
    band(p.mesh, d0, d1, N, donor);
    band(p.mesh, d1, g1, N, gap);
    band(p.mesh, g1, f1, N, fun);
    p.mesh.finish();
    addCyl(p, donor, 10.0, -2.0);
    Provisional g;
    g.claim = ProvClaim::ConsumedCylinder;
    g.tris = gap;
    p.work.provisionals.push_back(g);
    Provisional u;
    u.tris = fun;
    p.work.provisionals.push_back(u);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr, "gap funnel cone");
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 1, "gap funnel one");
    if (!c) return;
    expect(c->closed360, "gap funnel closed");
    expect((int)c->tris.size() == (int)fun.size(), "gap funnel triangles");
    expect(c->maxVertexDev <= tauOf(mv), "gap funnel within tau");
}

static void caseBlendBesideCone() {
    Part p;
    const int N = 16;
    auto c0 = ring(p.mesh, 10.0, -4.0, N);
    auto c1 = ring(p.mesh, 10.0, 0.0, N);
    auto c2 = ring(p.mesh, 14.0, 4.0, N);
    auto at = [&](double deg) {
        const double phi = deg * kPi / 180.0;
        return ring(p.mesh, 16.0 + 2.0 * std::cos(phi), 4.0 + 2.0 * std::sin(phi), N);
    };
    // Six stations on the circle. Every blend pair then has two further rings
    // on the circle on at least one side, so each pair is a profile row.
    auto t1 = at(150.0);
    auto t2 = at(120.0);
    auto t3 = at(90.0);
    auto t4 = at(60.0);
    auto t5 = at(30.0);
    std::vector<int> donor, cone, blend;
    band(p.mesh, c0, c1, N, donor);
    band(p.mesh, c1, c2, N, cone);
    band(p.mesh, c2, t1, N, blend);
    band(p.mesh, t1, t2, N, blend);
    band(p.mesh, t2, t3, N, blend);
    band(p.mesh, t3, t4, N, blend);
    band(p.mesh, t4, t5, N, blend);
    p.mesh.finish();
    addCyl(p, donor, 10.0, -2.0);
    Provisional pr;
    pr.tris = cone;
    pr.tris.insert(pr.tris.end(), blend.begin(), blend.end());
    p.work.provisionals.push_back(pr);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    Grab g;
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    const std::string diag = g.stop();
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr, "blend-side cone");
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 1, "blend-side one cone");
    expect(hasWhy(diag, "profile-row"), "blend edge row refused");
    if (!c) {
        std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
        return;
    }
    expect((int)c->tris.size() == (int)cone.size(), "cone does not take the blend row");
    expect(c->maxVertexDev <= tauOf(mv), "blend-side cone within tau");
}

static void caseThreeRing() {
    Part p = frustum(12, 3, 12, 8.0, 0.0, 14.0, 6.0, 4.0, false);
    MeshView mv = p.mesh.view();
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
    expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 1, "three-ring one cone");
    const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
    expect(c != nullptr && (int)c->tris.size() == (int)p.feature.size(), "three-ring takes every station");
    if (c) expect(c->closed360, "three-ring closed");
}

static void snapshot(const Region& c, double& Rlo, double& Rhi, double& h, int& n, double& vol) {
    Rlo = c.radius;
    Rhi = c.vMax;
    h = c.maxVertexSnap;
    n = (int)c.tris.size();
    vol = c.dVolPredicted;
}

static void caseThreads() {
    auto run = [](int threads, double& Rlo, double& Rhi, double& h, int& n, double& vol) {
        Part p = frustum(16, 2, 16, 10.0, 0.0, 14.0, 4.0, 4.0, false);
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol, threads);
        claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
        const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
        expect(c != nullptr, threads == 1 ? "threads1 cone" : "threadsN cone");
        if (!c) return;
        snapshot(*c, Rlo, Rhi, h, n, vol);
    };
    double a0, a1, a2, b0, b1, b2, va, vb;
    int na = 0, nb = 0;
    run(1, a0, a1, a2, na, va);
    run(8, b0, b1, b2, nb, vb);
    expect(a0 == b0 && a1 == b1 && a2 == b2 && na == nb && va == vb, "threads 1 and N agree");
}

int main() {
    caseClosed45();
    caseFunnel();
    caseQuarterTorus();
    caseTwoLevelNotTorus();
    caseSphere();
    caseSpindle();
    caseHorn();
    casePartial();
    caseBadDonor();
    caseWedge();
    caseProfileRow();
    caseCreaseFrustum();
    caseGapFunnel();
    caseBlendBesideCone();
    caseThreeRing();
    caseThreads();
    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_revolve_unit ok\n");
    return 0;
}
