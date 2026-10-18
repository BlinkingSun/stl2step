// Synthetic relational-revolve cases (D-fc-1.6, D-fc-1.7). No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_cone_math.hpp"
#include "refit_relational.hpp"
#include "refit_band_volume.hpp"
#include "refit_relational_arc.hpp"

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

// D-fc-5.7 O-3: revolvedChordTerm summed over a band is the exact chord volume.
// A ring: n vertices at axial u, radius rho, first azimuth phi.
struct BandRing {
    int n;
    double u, rho, phi;
};

static gp_XYZ ringPt(const BandRing& r, int k, const gp_XYZ& o, const gp_XYZ& ax,
                     const gp_XYZ& e1, const gp_XYZ& e2) {
    const double a = r.phi + 2.0 * kPi * (double)k / (double)r.n;
    return o + ax * r.u + e1 * (r.rho * std::cos(a)) + e2 * (r.rho * std::sin(a));
}

// Closed band between consecutive rings, stitched by azimuth (rings may carry
// different vertex counts and offsets, as an exporter's rows do). Triangles are
// oriented away from the axis.
static std::vector<std::array<gp_XYZ, 3>> bandTris(const std::vector<BandRing>& rings,
                                                   const gp_XYZ& o, const gp_XYZ& ax,
                                                   const gp_XYZ& e1, const gp_XYZ& e2) {
    std::vector<std::array<gp_XYZ, 3>> out;
    for (std::size_t j = 0; j + 1 < rings.size(); ++j) {
        const BandRing& A = rings[j];
        const BandRing& B = rings[j + 1];
        auto angA = [&](int i) { return A.phi + 2.0 * kPi * (double)i / (double)A.n; };
        auto angB = [&](int i) { return B.phi + 2.0 * kPi * (double)i / (double)B.n; };
        int ia = 0, ib = 0;
        while (ia < A.n || ib < B.n) {
            const bool stepA = ib >= B.n || (ia < A.n && angA(ia + 1) <= angB(ib + 1));
            std::array<gp_XYZ, 3> t;
            if (stepA) {
                t = {ringPt(A, ia, o, ax, e1, e2), ringPt(A, ia + 1, o, ax, e1, e2),
                     ringPt(B, ib, o, ax, e1, e2)};
                ++ia;
            } else {
                t = {ringPt(A, ia, o, ax, e1, e2), ringPt(B, ib + 1, o, ax, e1, e2),
                     ringPt(B, ib, o, ax, e1, e2)};
                ++ib;
            }
            const gp_XYZ c = (t[0] + t[1] + t[2]) / 3.0 - o;
            const gp_XYZ n = (t[1] - t[0]).Crossed(t[2] - t[0]);
            if ((c - ax * c.Dot(ax)).Dot(n) < 0.0) std::swap(t[1], t[2]);
            out.push_back(t);
        }
    }
    return out;
}

static double bandSum(const std::vector<std::array<gp_XYZ, 3>>& tris, const gp_XYZ& o,
                      const gp_XYZ& ax, const gp_XYZ& e1, const gp_XYZ& e2, double rho0,
                      double slope) {
    double v = 0;
    for (const auto& t : tris) v += revolvedChordTerm(t[0], t[1], t[2], o, ax, e1, e2, rho0, slope);
    return v;
}

// Independent reference: frustum volume minus the polyhedron the band's facets
// enclose with its two rim polygons (divergence theorem over a closed surface).
static double bandReference(const std::vector<BandRing>& rings,
                            const std::vector<std::array<gp_XYZ, 3>>& tris, const gp_XYZ& o,
                            const gp_XYZ& ax, const gp_XYZ& e1, const gp_XYZ& e2) {
    double six = 0;
    auto add = [&](const gp_XYZ& a, const gp_XYZ& b, const gp_XYZ& c) {
        six += a.Dot(b.Crossed(c));
    };
    for (const auto& t : tris) add(t[0], t[1], t[2]);
    const BandRing& lo = rings.front();
    const BandRing& hi = rings.back();
    const gp_XYZ cLo = o + ax * lo.u, cHi = o + ax * hi.u;
    for (int k = 0; k < lo.n; ++k)
        add(cLo, ringPt(lo, k + 1, o, ax, e1, e2), ringPt(lo, k, o, ax, e1, e2));
    for (int k = 0; k < hi.n; ++k)
        add(cHi, ringPt(hi, k, o, ax, e1, e2), ringPt(hi, k + 1, o, ax, e1, e2));
    const double h = hi.u - lo.u;
    const double frustum =
        kPi * h * (lo.rho * lo.rho + lo.rho * hi.rho + hi.rho * hi.rho) / 3.0;
    return frustum - six / 6.0;
}

static bool relClose(double a, double b, double rel) {
    return std::fabs(a - b) <= rel * std::max(std::fabs(a), std::fabs(b));
}

static void caseBandVolume() {
    const gp_XYZ o(3.0, -2.0, 7.0);
    const gp_XYZ ax = gp_XYZ(1.0, 2.0, 2.0) / 3.0;
    gp_XYZ e1 = gp_XYZ(2.0, -1.0, 0.0);
    e1 /= e1.Modulus();
    const gp_XYZ e2 = ax.Crossed(e1);
    // (1) Regular 24-gon prism: dVolCylinderSector's closed form.
    {
        const int N = 24;
        const double R = 5.0, h = 3.0;
        std::vector<BandRing> rings = {{N, 0.0, R, 0.0}, {N, h, R, 0.0}};
        const auto tris = bandTris(rings, o, ax, e1, e2);
        const double area = (double)N * 2.0 * R * std::sin(kPi / N) * h;
        const double v = bandSum(tris, o, ax, e1, e2, R, 0.0);
        expect(relClose(v, dVolCylinderSector(area, R, N, true), 1e-9),
               "band volume: regular prism = dVolCylinderSector");
    }
    // (2) Regular 36-gon frustum, aligned columns: coneFrustumChordVolume.
    {
        const int N = 36;
        const double R1 = 4.0, R2 = 7.0, h = 5.0;
        std::vector<BandRing> rings = {{N, 0.0, R1, 0.0}, {N, h, R2, 0.0}};
        const auto tris = bandTris(rings, o, ax, e1, e2);
        const double v = bandSum(tris, o, ax, e1, e2, R1, (R2 - R1) / h);
        expect(relClose(v, coneFrustumChordVolume(R1, R2, h, N), 1e-9),
               "band volume: regular frustum = coneFrustumChordVolume");
    }
    // (3) Rims of 2N, staggered interior rows of N (the exporter's S26 cone
    // pattern): exact against the independent polyhedron reference; one-N
    // densities miss it either way.
    const int N = 45;
    const double R1 = 11.5, R2 = 24.0, h = 19.25, k = (R2 - R1) / h;
    std::vector<BandRing> rings;
    rings.push_back({2 * N, 0.0, R1, 0.0});
    for (int j = 1; j < 7; ++j) {
        const double u = h * j / 7.0;
        rings.push_back({N, u, R1 + k * u, (j % 2) ? kPi / N : 0.0});
    }
    rings.push_back({2 * N, h, R2, 0.0});
    const auto tris = bandTris(rings, o, ax, e1, e2);
    const double ref = bandReference(rings, tris, o, ax, e1, e2);
    const double v = bandSum(tris, o, ax, e1, e2, R1, k);
    expect(ref > 0.0 && relClose(v, ref, 1e-9), "band volume: staggered rows exact");
    double area = 0, dRim = 0, dRow = 0;
    for (const auto& t : tris) {
        const gp_XYZ c = (t[0] + t[1] + t[2]) / 3.0 - o;
        const double a = 0.5 * (t[1] - t[0]).Crossed(t[2] - t[0]).Modulus();
        const double rho = R1 + k * c.Dot(ax);
        area += a;
        dRim += a * coneFrustumDVolDensity(rho, R1, R2, h, 2 * N);
        dRow += a * coneFrustumDVolDensity(rho, R1, R2, h, N);
    }
    expect(dRim < 0.5 * ref && dRow > ref, "band volume: no single side count fits");
    // (4) Frame and origin invariance: another e1, origin moved along the axis.
    {
        const double phi = 0.7;
        const gp_XYZ f1 = e1 * std::cos(phi) + e2 * std::sin(phi);
        const gp_XYZ f2 = ax.Crossed(f1);
        const double du = -4.25;
        const gp_XYZ o2 = o + ax * du;
        const double v2 = bandSum(tris, o2, ax, f1, f2, R1 + k * du, k);
        expect(relClose(v, v2, 1e-9), "band volume: frame and origin invariant");
    }
    // (5) A bore's winding (facets toward the axis) gives the same magnitude.
    {
        auto rev = tris;
        for (auto& t : rev) std::swap(t[1], t[2]);
        expect(relClose(bandSum(rev, o, ax, e1, e2, R1, k), v, 1e-12),
               "band volume: winding independent");
    }
    (void)area;
}

// ---- D-fc-5.11 F1 / F2 / F4 (fx-ringdonor) --------------------------------------

static std::vector<gp_XYZ> circlePts(double R, double z, int N, double cx = 0, double cy = 0) {
    std::vector<gp_XYZ> v;
    for (int i = 0; i < N; ++i) {
        const double th = kPi * 2.0 * (double)i / (double)N;
        v.emplace_back(cx + R * std::cos(th), cy + R * std::sin(th), z);
    }
    return v;
}

// F1 ring certificate: planar and circular within tau, minStations(kParamsCircle).
static void caseMeshRing() {
    const double tau = 2e-6;
    MeshRing r;
    expect(certifyMeshRing(circlePts(5.0, 1.0, 24, 2.0, -3.0), tau, r), "ring: circle accepted");
    expect(std::fabs(r.R - 5.0) <= tau && std::fabs(std::fabs(r.n.Z()) - 1.0) <= 1e-9 &&
               (r.c - gp_XYZ(2.0, -3.0, 1.0)).Modulus() <= tau && r.stations == 24,
           "ring: centre, normal, radius, stations");
    expect(!certifyMeshRing(circlePts(5.0, 0.0, 3), tau, r), "ring: 3 points under-determined");
    std::vector<gp_XYZ> bump = circlePts(5.0, 0.0, 24);
    bump[5].SetZ(10.0 * tau);
    expect(!certifyMeshRing(bump, tau, r), "ring: off-plane vertex refused");
    std::vector<gp_XYZ> off = circlePts(5.0, 0.0, 24);
    off[7] = off[7] * (1.0 + 10.0 * tau / 5.0);
    expect(!certifyMeshRing(off, tau, r), "ring: off-circle vertex refused");
    expect(!certifyMeshRing(circlePts(0.5 * tau, 0.0, 24), tau, r), "ring: radius <= tau refused");
    expect(!certifyMeshRing(circlePts(5.0, 0.0, 24), 0.0, r), "ring: tau 0 admits nothing");
}

// F1 helper: the closed boundary loops of a triangle set.
static void caseTriSetLoops() {
    Synth s;
    const int N = 16;
    auto a = ring(s, 5.0, 0.0, N);
    auto b = ring(s, 5.0, 2.0, N);
    std::vector<int> all, half;
    band(s, a, b, N, all);
    s.finish();
    for (int t = 0; t < N; ++t) half.push_back(all[(std::size_t)t]);
    const MeshView mv = s.view();
    const auto loops = triSetLoops(mv, all);
    expect(loops.size() == 2 && loops[0].size() == (std::size_t)N &&
               loops[1].size() == (std::size_t)N,
           "loops: closed wall has its two rims");
    expect(triSetLoops(mv, half).size() == 1, "loops: half wall is a disk (one loop)");
}

// A planar annulus chart at z=0 (Ri..Ro) and, below it, a closed cone band from
// the annulus outer ring (Ro at z=0) to (Rc at z=-h). No region donor.
struct RingPart {
    Synth mesh;
    SegmentWork work;
    std::vector<int> annulus, cone;
};

static RingPart ringCone(int N, double Ri, double Ro, double Rc, double h, bool everyOther,
                         bool annulusChart) {
    RingPart p;
    auto in = ring(p.mesh, Ri, 0.0, N);
    auto out = ring(p.mesh, Ro, 0.0, N);
    auto low = ring(p.mesh, Rc, -h, N);
    band(p.mesh, in, out, N, p.annulus);
    std::vector<int> cone;
    band(p.mesh, out, low, N, cone);
    for (std::size_t i = 0; i < cone.size(); i += 2)
        if (!everyOther || (i / 2) % 2 == 0) {
            p.cone.push_back(cone[i]);
            p.cone.push_back(cone[i + 1]);
        }
    p.mesh.finish();
    Provisional ann;
    ann.claim = ProvClaim::Unclaimed;
    ann.tris = p.annulus;
    if (annulusChart) {
        p.work.provisionals.push_back(ann);
    } else {
        // The annulus split into two charts: neither has two loops.
        Provisional a1 = ann, a2 = ann;
        a1.tris.assign(p.annulus.begin(), p.annulus.begin() + (long)N);
        a2.tris.assign(p.annulus.begin() + (long)N, p.annulus.end());
        p.work.provisionals.push_back(a1);
        p.work.provisionals.push_back(a2);
    }
    Provisional cp;
    cp.claim = ProvClaim::Unclaimed;
    cp.tris = p.cone;
    p.work.provisionals.push_back(cp);
    return p;
}

static void caseAnnulusRings() {
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, false, true);
        const MeshView mv = p.mesh.view();
        const auto a = annulusRingCharts(mv, p.work, 2.0 * p.mesh.q);
        expect(a.size() == 1 && a[0].chart == 0, "annulus: one chart");
        if (a.size() == 1) {
            const double r0 = a[0].ring[0].R, r1 = a[0].ring[1].R;
            expect(std::fabs(std::min(r0, r1) - 5.0) < 1e-6 &&
                       std::fabs(std::max(r0, r1) - 8.0) < 1e-6,
                   "annulus: both rims certified");
        }
    }
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, false, false);
        const MeshView mv = p.mesh.view();
        expect(annulusRingCharts(mv, p.work, 2.0 * p.mesh.q).empty(),
               "annulus: split chart (one loop each) is no ring donor");
    }
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, false, true);
        Region owner;
        owner.type = SurfType::Plane;
        owner.tris = {p.annulus.front()};
        p.work.accepted.push_back(owner);
        const MeshView mv = p.mesh.view();
        expect(annulusRingCharts(mv, p.work, 2.0 * p.mesh.q).empty(),
               "annulus: chart holding an accepted triangle refused");
    }
}

// D-fc-5.17 F-2: an N-gon plate outline concentric with a round bore is no ring.
static void caseMeshRingTessLaw() {
    const double tau = 2e-6;
    MeshRing outer, inner;
    expect(certifyMeshRing(circlePts(10.0, 0.0, 4), tau, outer), "tesslaw: square certifies");
    expect(certifyMeshRing(circlePts(4.0, 0.0, 48), tau, inner), "tesslaw: bore certifies");
    outer.verts.assign(4, 0);
    inner.verts.assign(48, 0);
    Synth s;
    SegmentWork work;
    SegmentParams sp;
    std::vector<const MeshRing*> peers = {&inner};
    expect(meshRingRefutedByTessLaw(s.view(), sp, tau, work, outer, peers),
           "tesslaw: square refuted by round bore witness");
    peers = {&outer};
    expect(!meshRingRefutedByTessLaw(s.view(), sp, tau, work, inner, peers),
           "tesslaw: round bore not refuted by square");
}

// D-fc-5.20 B-1. Two flat annuli and no cylinder rail, so the only F-2 witness
// is the other chart. The washer (48 stations, radii 5 and 6) refutes the
// square frame (4 stations, radii 10 and 20). Both provisional orders must
// keep the washer and refuse the square: accepting the washer must not empty
// its rings before the square is tested.
static void pushUnclaimed(SegmentWork& work, const std::vector<int>& tris) {
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = tris;
    work.provisionals.push_back(std::move(pr));
}

static void addAnnulus(Synth& s, double Ri, double Ro, double z, int N, std::vector<int>& tris) {
    const std::vector<int> inn = ring(s, Ri, z, N);
    const std::vector<int> out = ring(s, Ro, z, N);
    band(s, inn, out, N, tris);
}

static void caseAnnulusWitnessOrder() {
    auto run = [](bool washerFirst) {
        Synth s;
        std::vector<int> washerTris, squareTris;
        auto addWasher = [&]() { addAnnulus(s, 5.0, 6.0, 0.0, 48, washerTris); };
        auto addSquare = [&]() { addAnnulus(s, 10.0, 20.0, 5.0, 4, squareTris); };
        if (washerFirst) {
            addWasher();
            addSquare();
        } else {
            addSquare();
            addWasher();
        }
        s.finish();
        SegmentWork work;
        if (washerFirst) {
            pushUnclaimed(work, washerTris);
            pushUnclaimed(work, squareTris);
        } else {
            pushUnclaimed(work, squareTris);
            pushUnclaimed(work, washerTris);
        }
        return annulusRingCharts(s.view(), work, 2.0 * s.q);
    };
    for (int ord = 0; ord < 2; ++ord) {
        const bool washerFirst = ord == 0;
        const auto charts = run(washerFirst);
        const char* tag = washerFirst ? "washer-first" : "square-first";
        char buf[160];
        std::snprintf(buf, sizeof(buf), "witness order (%s): only the washer chart", tag);
        bool onlyWasher = charts.size() == 1;
        if (onlyWasher) {
            const double lo = std::min(charts[0].ring[0].R, charts[0].ring[1].R);
            const double hi = std::max(charts[0].ring[0].R, charts[0].ring[1].R);
            const int wantChart = washerFirst ? 0 : 1;
            onlyWasher = std::fabs(lo - 5.0) < 1e-6 && std::fabs(hi - 6.0) < 1e-6 &&
                         charts[0].chart == wantChart;
        }
        expect(onlyWasher, buf);
        if (!onlyWasher) {
            std::fprintf(stderr, "  charts=%zu\n", charts.size());
            for (const AnnulusRings& a : charts)
                std::fprintf(stderr, "  chart=%d R=%.6f/%.6f stations=%d/%d\n", a.chart, a.ring[0].R,
                             a.ring[1].R, a.ring[0].stations, a.ring[1].stations);
        }
    }
}

// F1 end to end: the annulus's rings donate the cone's axis; scattered columns
// (every other sector) pass the column rule but are no wall.
static void caseRingDonorCone() {
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, false, true);
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol, 1);
        Grab g;
        claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
        const std::string diag = g.stop();
        const Region* c = findNew(p.work, SurfType::Cone, Origin::ChamferCone);
        expect(c != nullptr && c->closed360 && c->tris.size() == p.cone.size(),
               "ring donor: closed cone claimed");
        expect(diag.find("DIAG_REL_RINGDONOR") != std::string::npos, "ring donor: diag line");
        if (c) {
            expect(std::fabs(c->radius - 8.0) < 1e-6 && std::fabs(c->vMax - 9.0) < 1e-6 &&
                       std::fabs(c->ax.Direction().Z()) > 1.0 - 1e-9,
                   "ring donor: cone on the ring axis");
        } else {
            std::fprintf(stderr, "DIAG\n%s\n", diag.c_str());
        }
    }
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, false, false);
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol, 1);
        Grab g;
        claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
        g.stop();
        expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0,
               "ring donor: no annulus, no donor, no cone");
    }
    {
        RingPart p = ringCone(32, 5.0, 8.0, 9.0, 3.0, true, true);
        MeshView mv = p.mesh.view();
        SegmentParams sp;
        DerivedTols tol;
        params(sp, tol, 1);
        Grab g;
        claimRelRevolve(mv, sp, tol, p.work, RelRevolveMode::Line);
        const std::string diag = g.stop();
        expect(countType(p.work, SurfType::Cone, Origin::ChamferCone) == 0,
               "ring donor: scattered columns not claimed");
        expect(hasWhy(diag, "domain-not-annulus"), "ring donor: domain-not-annulus why");
    }
}

// F4: the apex. A cone may stand on its apex rim (radius 0); a vertex on the axis
// is no azimuth station of a partial cone.
static void caseConeApex() {
    gp_Cone cone;
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    expect(coneFromRims(ax, 0.0, 1.0, 1.0, cone) && cone.RefRadius() == 0.0 &&
               std::fabs(cone.SemiAngle() - 0.25 * kPi) < 1e-12,
           "apex: cone from an apex rim");
    expect(!coneFromRims(ax, -1e-9, 1.0, 1.0, cone), "apex: negative rim radius refused");
    expect(!coneFromRims(ax, 0.0, 0.0, 1.0, cone), "apex: both rims at the axis refused");
    const double tau = 2e-6;
    std::vector<std::vector<gp_XYZ>> fam{circlePts(3.0, 0.0, 48), circlePts(3.0, 2.0, 48)};
    std::vector<gp_XYZ> band;
    for (int st = 0; st <= 4; ++st) {
        const double z = 0.25 * st, R = 1.0 - z;
        if (R <= 0.0) continue;
        for (int i = 0; i <= 8; ++i) {
            const double th = 0.5 * kPi * (double)i / 8.0;
            band.emplace_back(R * std::cos(th), R * std::sin(th), z);
        }
    }
    const PartialConeCert a = certifyPartialCone(fam, band, gp_XYZ(0, 0, 0), gp_XYZ(0, 0, 1), tau);
    band.emplace_back(0.0, 0.0, 1.0);  // the apex
    const PartialConeCert b = certifyPartialCone(fam, band, gp_XYZ(0, 0, 0), gp_XYZ(0, 0, 1), tau);
    expect(a.ok && b.ok, "apex: partial cone certifies with its apex");
    expect(a.cols == b.cols && b.stations == a.stations + 1,
           "apex: apex is a line station, not an azimuth column");
}

// F2: a same-radius coaxial band merges into its member by the D-fc-4.4 closure.
struct MergePart {
    Synth mesh;
    SegmentWork work;
    std::vector<int> member, band;
};

static MergePart halfWall(double bandR, double bandZ0, bool closedMember) {
    MergePart p;
    const int N = 32;
    auto a = ring(p.mesh, 10.0, 0.0, N);
    auto b = ring(p.mesh, 10.0, 4.0, N);
    std::vector<int> wall;
    band(p.mesh, a, b, 16, wall);
    for (std::size_t i = 0; i < 16; ++i) p.member.push_back(wall[i]);
    if (bandR == 10.0 && bandZ0 == 0.0) {
        for (std::size_t i = 16; i < wall.size(); ++i) p.band.push_back(wall[i]);
    } else {
        auto c = ring(p.mesh, bandR, bandZ0, N);
        auto d = ring(p.mesh, bandR, bandZ0 + 4.0, N);
        std::vector<int> other;
        band(p.mesh, c, d, 16, other);
        for (std::size_t i = 16; i < other.size(); ++i) p.band.push_back(other[i]);
    }
    p.mesh.finish();
    Region r;
    r.id = 0;
    r.type = SurfType::Cylinder;
    r.origin = Origin::CylGrow;
    r.ax = gp_Ax3(gp_Pnt(0, 0, 2.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.radius = 10.0;
    r.closed360 = closedMember;
    r.tris = p.member;
    p.work.accepted.push_back(r);
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = p.band;
    p.work.provisionals.push_back(pr);
    return p;
}

static void caseSameSurfaceMerge() {
    SegmentParams sp;
    DerivedTols tol;
    params(sp, tol, 1);
    {
        MergePart p = halfWall(10.0, 0.0, false);
        const MeshView mv = p.mesh.view();
        const bool ok = mergeSameSurfaceBand(mv, sp, tol, p.work, 0, p.band, gp_XYZ(0, 0, 0),
                                             gp_XYZ(0, 0, 1), 10.0, 10.0);
        expect(ok && p.work.accepted.size() == 1 &&
                   p.work.accepted[0].tris.size() == p.member.size() + p.band.size() &&
                   std::fabs(p.work.accepted[0].radius - 10.0) < 1e-6 &&
                   p.work.provisionals[0].claim != ProvClaim::Unclaimed,
               "merge: band joins its member");
    }
    {
        const double tau = 2e-6;
        MergePart p = halfWall(10.0 + 10.0 * tau, 0.0, false);
        const MeshView mv = p.mesh.view();
        expect(!mergeSameSurfaceBand(mv, sp, tol, p.work, 0, p.band, gp_XYZ(0, 0, 0),
                                     gp_XYZ(0, 0, 1), 10.0, 10.0 + 10.0 * tau) &&
                   p.work.accepted[0].tris.size() == p.member.size(),
               "merge: off-radius band refused, member unchanged");
    }
    {
        MergePart p = halfWall(10.0, 8.0, false);
        const MeshView mv = p.mesh.view();
        expect(!mergeSameSurfaceBand(mv, sp, tol, p.work, 0, p.band, gp_XYZ(0, 0, 0),
                                     gp_XYZ(0, 0, 1), 10.0, 10.0),
               "merge: stacked (no axial overlap) refused");
    }
    {
        MergePart p = halfWall(10.0, 0.0, true);
        const MeshView mv = p.mesh.view();
        expect(!mergeSameSurfaceBand(mv, sp, tol, p.work, 0, p.band, gp_XYZ(0, 0, 0),
                                     gp_XYZ(0, 0, 1), 10.0, 10.0),
               "merge: closed member refused");
    }
    {
        MergePart p = halfWall(10.0, 0.0, false);
        p.work.accepted[0].origin = Origin::FilletStrip;
        const MeshView mv = p.mesh.view();
        expect(!mergeSameSurfaceBand(mv, sp, tol, p.work, 0, p.band, gp_XYZ(0, 0, 0),
                                     gp_XYZ(0, 0, 1), 10.0, 10.0) &&
                   p.work.accepted[0].tris.size() == p.member.size(),
               "merge: FilletStrip member refused (D-fc-4.4 class is CylGrow)");
    }
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
    caseBandVolume();
    caseMeshRing();
    caseMeshRingTessLaw();
    caseAnnulusWitnessOrder();
    caseTriSetLoops();
    caseAnnulusRings();
    caseRingDonorCone();
    caseConeApex();
    caseSameSurfaceMerge();
    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_revolve_unit ok\n");
    return 0;
}
