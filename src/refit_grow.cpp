// P1 grow lane — charts A1, provisional A2, cylinder claim B1, plane commit A3.
// B1 gate structure: DECISION-p1-growx (D1.3-A1..A6, D2.3-A1).

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <limits>
#include <mutex>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace stl2step { namespace refit {
namespace {

constexpr double kTiny = 1e-30;
// Angular band for seedInBand / grow phi gates: scale with the gate magnitude so
// libm ULP noise at theta_cyl_lo (~3 deg) cannot flip membership across targets.
inline double angleBandEps(double thetaRad) {
    constexpr double kAbs = 1e-12;
    const double ulp = std::numeric_limits<double>::epsilon()
                       * std::max(1.0, std::abs(thetaRad));
    return std::max(kAbs, 8.0 * ulp);
}
// D1.3-A3 running residual: file-local. Frozen header epsCylGrow is not used by B1.
constexpr double kRingResidualFrac = 0.25;

bool p1DiagOn() {
    const char* e = std::getenv("STL2STEP_P1_DIAG");
    return e && e[0] != '\0' && e[0] != '0';
}

bool collapseDiagOn() {
    const char* e = std::getenv("STL2STEP_COLLAPSE_DIAG");
    return e && e[0] != '\0' && e[0] != '0';
}

bool lawbandDiagOn() {
    const char* e = std::getenv("STL2STEP_LAWBAND_DIAG");
    return e && e[0] != '\0' && e[0] != '0';
}

// D-130-14 / D-130-16 -- the loop-level union lands in stages, and a stage that
// cannot yet build a face may not turn a built part faceted. Until the face
// path closes (the two-arc rim collapse and the partial-cylinder wire), the
// union's claiming, its surface-connected region emit and its interruption-
// paired loop stitch are compiled in but not enabled: with this off every
// region, chain and face is bit-for-bit what the branch tip produced. It is a
// landing stage, not a per-part switch -- it names no part and reads no
// geometry.
bool lawUnionOn() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_UNION");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// Diag only (STL2STEP_DIAG_A2): A2 grows PLANAR provisionals, so every link it
// refuses is a statement about the mesh. Naming the predicate per refused link
// is what turns "the ring arrives in pieces" into a measurement. Prints nothing
// and decides nothing when off.
bool a2DiagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_A2");
    return e && e[0] != '\0' && e[0] != '0';
}

template <typename F>
void lawParallelFor(size_t n, F&& fn) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;
    if (n < 8 || hw == 1) {
        for (size_t i = 0; i < n; i++) fn(i);
        return;
    }
    std::vector<std::thread> pool;
    const size_t chunk = (n + hw - 1) / hw;
    for (unsigned t = 0; t < hw; t++) {
        const size_t lo = static_cast<size_t>(t) * chunk;
        const size_t hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        pool.emplace_back([&, lo, hi]() {
            for (size_t i = lo; i < hi; i++) fn(i);
        });
    }
    for (auto& th : pool) th.join();
}

bool triInLawBand(const SegmentWork& work, int t) {
    for (const Region& r : work.accepted) {
        if (!r.lawBand) continue;
        for (int u : r.tris) {
            if (u == t) return true;
        }
    }
    return false;
}

// Header export archChainBand is 500–8000. Grow commit is a strict subset
// outside coarseFusionBand. Body18/20 are 1948 tris. Body11 is two components
// (12060 + 3240); the 3240 piece sits inside the detector band and must stay
// commit-neutral (AC2 recognition census). Cap at 2500 so that component,
// Body9 (4668) and Body12 (7918) do not take the new apply. Over-commit is
// worse. Do not widen coarseFusionBand — S13/S14 stay outside both.
inline bool inArchChainCommitBand(const MeshView& mv) {
    return !coarseFusionBand(mv) && archChainBand(mv) && mv.nTri <= 2500;
}

// Handshake with refit_build.cpp J6 sewability ledger: in-band arch-chain
// commits stamp filletNbrA so new plane|cyl / cyl|cyl boundaries are guarded.
// FilletStrip is the only other writer; it keys off Origin::FilletStrip first.
constexpr int kArchChainSewMark = -2;

void diagArchChain(const char* kind, std::size_t n, double score, double eberly,
                   double chainR, double appliedR, double rel, const char* why) {
    if (!collapseDiagOn()) return;
    std::fprintf(stderr,
                 "DIAG_ARCHCHAIN_%s n=%zu score=%.4g eberly=%.6g chainR=%.6g "
                 "appliedR=%.6g rel=%.4g why=%s\n",
                 kind, n, score, eberly, chainR, appliedR, rel, why);
}

struct EdgeAdj {
    std::vector<std::array<int, 2>> tri;  // local edge -> [t0, t1], -1 if boundary
};

gp_XYZ localTriVert(const MeshView& mv, int lt, int corner) {
    const int gt = mv.compTris[lt];
    const int gv = mv.tris[gt][corner];
    return mv.pts[gv];
}

double triAreaLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    return 0.5 * (b - a).Crossed(c - a).Modulus();
}

gp_XYZ triCentroidLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    return (a + b + c) / 3.0;
}

gp_Dir triNormalLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    gp_XYZ n = (b - a).Crossed(c - a);
    const double m = n.Modulus();
    if (m < kTiny) return gp_Dir(0, 0, 1);
    n /= m;
    return gp_Dir(n.X(), n.Y(), n.Z());
}

EdgeAdj buildEdgeAdj(const MeshView& mv) {
    EdgeAdj ea;
    ea.tri.resize(mv.nEdge, { -1, -1 });
    for (int lt = 0; lt < static_cast<int>(mv.nTri); lt++) {
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[lt][s];
            if (ea.tri[e][0] < 0) ea.tri[e][0] = lt;
            else ea.tri[e][1] = lt;
        }
    }
    return ea;
}

// Dihedral in [0, 180°]: acos(n0·n1), not acos(|n0·n1|). Obtuse ridges must stay obtuse.
double edgeDihedralAbs(const MeshView& mv, int edgeId, const EdgeAdj& ea) {
    const int t0 = ea.tri[edgeId][0];
    const int t1 = ea.tri[edgeId][1];
    if (t0 < 0 || t1 < 0) return 0.0;
    const gp_Dir n0 = triNormalLocal(mv, t0);
    const gp_Dir n1 = triNormalLocal(mv, t1);
    const double d = std::min(1.0, std::max(-1.0, n0.Dot(n1)));
    return std::acos(d);
}

struct UnionFind {
    std::vector<int> p;
    explicit UnionFind(int n) : p(n) {
        for (int i = 0; i < n; i++) p[i] = i;
    }
    int find(int a) {
        while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; }
        return a;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) p[a] = b;
    }
};

int minTriId(const Provisional& p) {
    return p.tris.empty() ? std::numeric_limits<int>::max() : p.tris.front();
}

void sortProvisionals(std::vector<Provisional>& provs) {
    std::sort(provs.begin(), provs.end(),
              [](const Provisional& a, const Provisional& b) {
                  return minTriId(a) < minTriId(b);
              });
}

void sortRegions(std::vector<Region>& regs) {
    std::sort(regs.begin(), regs.end(), [](const Region& a, const Region& b) {
        const int ma = a.tris.empty() ? std::numeric_limits<int>::max() : a.tris.front();
        const int mb = b.tris.empty() ? std::numeric_limits<int>::max() : b.tris.front();
        return ma < mb;
    });
}

std::vector<int> mergeMemberTris(const std::vector<Provisional>& provs,
                                 const std::vector<int>& members) {
    std::vector<int> tris;
    for (int m : members)
        for (int t : provs[m].tris) tris.push_back(t);
    std::sort(tris.begin(), tris.end());
    tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
    return tris;
}

gp_Dir normalizeXYZ(const gp_XYZ& v) {
    const double m = v.Modulus();
    if (m < kTiny) return gp_Dir(1, 0, 0);
    return gp_Dir(v.X() / m, v.Y() / m, v.Z() / m);
}

// Canonical axis sign: first component of w with |w_k| > 1e-9, X then Y then Z.
// Linear-algebra component extraction of the eigenvector (D1.3-A4).
double eigenAxisSign(const gp_XYZ& w) {
    if (std::abs(w.X()) > 1e-9) return w.X() > 0.0 ? 1.0 : -1.0;
    if (std::abs(w.Y()) > 1e-9) return w.Y() > 0.0 ? 1.0 : -1.0;
    if (std::abs(w.Z()) > 1e-9) return w.Z() > 0.0 ? 1.0 : -1.0;
    return 1.0;
}

gp_Dir canonicalAxis(const gp_XYZ& w) {
    return normalizeXYZ(w * eigenAxisSign(w));
}

struct GaussResult {
    bool ok = false;
    bool degenerate = false;
    gp_Dir axis;
    double mu1 = 0, mu2 = 0, mu3 = 0;
    double flat = 0;   // mu1 / max(mu2, 1e-300)
    double patch = 0;  // mu2 / max(mu3, 1e-300)
    double c = 0;      // nbar · a
    double dev = 0;    // max_t |n_t · a − c|
};

GaussResult centeredGauss(const MeshView& mv, const std::vector<int>& tris,
                          const gp_Dir& seedAxis, const DerivedTols& tol) {
    GaussResult r;
    if (tris.empty()) return r;

    double A = 0.0;
    gp_XYZ nbar(0, 0, 0);
    for (int lt : tris) {
        const gp_Dir n = triNormalLocal(mv, lt);
        const double a = triAreaLocal(mv, lt);
        A += a;
        nbar += gp_XYZ(n.X(), n.Y(), n.Z()) * a;
    }
    if (A < kTiny) return r;
    nbar /= A;

    std::array<std::array<double, 3>, 3> C{};
    for (int lt : tris) {
        const gp_Dir n = triNormalLocal(mv, lt);
        const double a = triAreaLocal(mv, lt);
        const double dx = n.X() - nbar.X();
        const double dy = n.Y() - nbar.Y();
        const double dz = n.Z() - nbar.Z();
        C[0][0] += a * dx * dx;
        C[0][1] += a * dx * dy;
        C[0][2] += a * dx * dz;
        C[1][1] += a * dy * dy;
        C[1][2] += a * dy * dz;
        C[2][2] += a * dz * dz;
    }
    C[1][0] = C[0][1];
    C[2][0] = C[0][2];
    C[2][1] = C[1][2];

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(C, eval, evec)) return r;
    r.mu1 = eval[0];
    r.mu2 = eval[1];
    r.mu3 = eval[2];
    r.flat = r.mu1 / std::max(r.mu2, 1e-300);
    r.patch = r.mu2 / std::max(r.mu3, 1e-300);

    // R1 (adjudicate-p1-grow): c1 = mu1/mu2 is a pure ratio and is not
    // scale-free on small claims. Eigenvalues of the area-weighted Gauss
    // covariance C = Σ area (n−nbar)(n−nbar)^T have units of area; a
    // uniform angular perturbation θ of the unit-normal cloud produces
    // mu ~ A θ².
    //
    // Linear budget is DerivedTols::epsMesh, which D5.2 resolves from
    // MeshView as max(weldTol, 1e-4*diag, 1e-3) — weldTol and diag, not
    // a new constant. A vertex displaced by epsMesh on a body of size
    // diag tilts a unit normal by θ = epsMesh/diag. If mu2 < A θ², the
    // in-plane Gauss spread is below the tessellation-noise floor and
    // mu1/mu2 is noise/noise: a 1e-6-relative radial jitter (≈300× below
    // epsMesh) pushed a 2–3 band cylinder claim over 0.05 and shattered
    // it into planes. Treat that cloud as degenerate/planar (flat = 0),
    // same as the rank-deficient mu2 ≪ mu3 path, so growth can reach a
    // claim large enough for the ratio to be meaningful.
    const double L = std::max(mv.diag, tol.epsMesh);
    const double theta = (L > 0.0) ? (tol.epsMesh / L) : 0.0;
    const double mu2Floor = A * theta * theta;
    const bool fewNormals = r.mu2 <= 1e-12 * r.mu3;
    const bool mu2BelowNoise = r.mu2 <= mu2Floor;
    if (fewNormals || mu2BelowNoise) {
        r.degenerate = true;
        r.axis = canonicalAxis(gp_XYZ(seedAxis.X(), seedAxis.Y(), seedAxis.Z()));
        r.c = 0.0;
        r.dev = 0.0;
        r.flat = 0.0;
        r.patch = 0.0;
        r.ok = true;
        return r;
    }

    // w1 = column 0 of evec (jacobi stores eigenvectors as columns).
    const gp_XYZ w1(evec[0][0], evec[1][0], evec[2][0]);
    r.axis = canonicalAxis(w1);
    const gp_XYZ axyz(r.axis.X(), r.axis.Y(), r.axis.Z());
    r.c = nbar.Dot(axyz);
    r.dev = 0.0;
    for (int lt : tris) {
        const double d = triNormalLocal(mv, lt).Dot(r.axis);
        r.dev = std::max(r.dev, std::abs(d - r.c));
    }
    r.ok = true;
    return r;
}

bool testT1Running(const GaussResult& g) {
    return g.ok && g.flat < DerivedTols::kGaussPlanarity;
}

bool testG1Commit(const GaussResult& g, const DerivedTols& tol) {
    const double sin3 = tol.gaussAxisTiltSin();
    return g.ok && g.flat < DerivedTols::kGaussPlanarity && g.dev < sin3
           && std::abs(g.c) < sin3;
}

gp_XYZ areaWeightedNbar(const MeshView& mv, const std::vector<int>& tris) {
    double A = 0.0;
    gp_XYZ nbar(0, 0, 0);
    for (int lt : tris) {
        const gp_Dir n = triNormalLocal(mv, lt);
        const double a = triAreaLocal(mv, lt);
        A += a;
        nbar += gp_XYZ(n.X(), n.Y(), n.Z()) * a;
    }
    if (A < kTiny) return gp_XYZ(0, 0, 0);
    nbar /= A;
    return nbar;
}

void axisTiltStats(const MeshView& mv, const std::vector<int>& tris, const gp_Dir& axis,
                   double& cOut, double& devOut) {
    const gp_XYZ nbar = areaWeightedNbar(mv, tris);
    if (nbar.SquareModulus() <= 0.0 && tris.empty()) {
        cOut = 0.0;
        devOut = 0.0;
        return;
    }
    const gp_XYZ axyz(axis.X(), axis.Y(), axis.Z());
    cOut = nbar.Dot(axyz);
    devOut = 0.0;
    for (int lt : tris) {
        const double d = triNormalLocal(mv, lt).Dot(axis);
        devOut = std::max(devOut, std::abs(d - cOut));
    }
}

bool testG1CommitSeedAxis(const GaussResult& g, const DerivedTols& tol, double cTilt,
                          double devTilt) {
    const double sin3 = tol.gaussAxisTiltSin();
    return g.ok && g.flat < DerivedTols::kGaussPlanarity && devTilt < sin3
           && std::abs(cTilt) < sin3;
}

double epsCylRing(const DerivedTols& tol, double R) {
    return std::max(tol.epsMesh, kRingResidualFrac * R);
}

double maxVertexResidual(const MeshView& mv, const std::vector<int>& tris,
                           const gp_Dir& axis, const gp_Pnt& center, double radius) {
    const gp_XYZ a(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ c(center.X(), center.Y(), center.Z());
    double maxR = 0.0;
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            const gp_XYZ p = localTriVert(mv, lt, k);
            const gp_XYZ d = p - c;
            const double radial = a.Crossed(d).Modulus();
            maxR = std::max(maxR, std::abs(radial - radius));
        }
    }
    return maxR;
}

double medianOf(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Median in-patch dihedral → full-circle band count (D2 nSides fallback).
int estimateNBandsFromPatch(const MeshView& mv, const std::vector<int>& tris) {
    std::vector<char> in(static_cast<size_t>(mv.nTri), 0);
    for (int t : tris)
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) in[static_cast<size_t>(t)] = 1;
    const EdgeAdj ea = buildEdgeAdj(mv);
    std::vector<double> phis;
    for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
        const int t0 = ea.tri[e][0];
        const int t1 = ea.tri[e][1];
        if (t0 < 0 || t1 < 0) continue;
        if (!in[static_cast<size_t>(t0)] || !in[static_cast<size_t>(t1)]) continue;
        const double phi = edgeDihedralAbs(mv, e, ea);
        if (phi > 0.05 && phi < M_PI - 0.05) phis.push_back(phi);
    }
    if (phis.empty()) return 1;
    const double med = medianOf(phis);
    if (med < 1e-6) return 1;
    return std::max(1, static_cast<int>(std::llround(2.0 * M_PI / med)));
}

double percentile75(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t idx = static_cast<size_t>(std::floor(0.75 * (v.size() - 1)));
    return v[idx];
}

double wrapToPi(double t) {
    while (t <= -M_PI) t += 2.0 * M_PI;
    while (t > M_PI) t -= 2.0 * M_PI;
    return t;
}

struct D2Metrics {
    int nBands = 0;
    double span = 0;
    int nSides = 0;
    double uMin = 0, uMax = 0;
    double vMin = 0, vMax = 0;
    bool closed360 = false;
    gp_Ax3 ax;
    gp_Pnt center;
    double radius = 0;
    bool spanReject = false;
};

D2Metrics computeD2(const MeshView& mv, const std::vector<int>& tris,
                    const gp_Dir& axisIn, const gp_Pnt& centerIn, double radiusIn,
                    const DerivedTols& tol) {
    D2Metrics m;
    m.center = centerIn;
    m.radius = radiusIn;
    const gp_Dir a = axisIn;
    const gp_XYZ aw(a.X(), a.Y(), a.Z());

    gp_XYZ centroid(0, 0, 0);
    double totalArea = 0.0;
    for (int lt : tris) {
        const double ar = triAreaLocal(mv, lt);
        centroid += triCentroidLocal(mv, lt) * ar;
        totalArea += ar;
    }
    if (totalArea > kTiny) centroid /= totalArea;

    const gp_XYZ cxyz(centerIn.X(), centerIn.Y(), centerIn.Z());
    const double tLoc = (centroid - cxyz).Dot(aw);
    const gp_XYZ loc = cxyz + aw * tLoc;
    const gp_Pnt Loc(loc);

    std::vector<std::pair<int, gp_XYZ>> uniq;
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            const int gt = mv.compTris[lt];
            const int gvi = mv.tris[gt][k];
            int lvi = -1;
            for (size_t i = 0; i < mv.nVtx; i++) {
                if (mv.compVtx[i] == gvi) { lvi = static_cast<int>(i); break; }
            }
            if (lvi < 0) continue;
            uniq.emplace_back(lvi, localTriVert(mv, lt, k));
        }
    }
    std::sort(uniq.begin(), uniq.end(),
              [](const auto& x, const auto& y) { return x.first < y.first; });
    uniq.erase(std::unique(uniq.begin(), uniq.end(),
                           [](const auto& x, const auto& y) { return x.first == y.first; }),
               uniq.end());

    int minLv = std::numeric_limits<int>::max();
    gp_XYZ ps(0, 0, 0);
    for (const auto& uv : uniq) {
        if (uv.first < minLv) { minLv = uv.first; ps = uv.second; }
    }
    const gp_XYZ dps = ps - loc;
    const double ads = dps.Dot(aw);
    gp_XYZ xDir = dps - aw * ads;
    if (xDir.Modulus() < 1e-12) {
        gp_XYZ tmp = (std::abs(aw.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
        xDir = aw.Crossed(tmp);
    }
    const gp_Dir xD = normalizeXYZ(xDir);

    m.ax = gp_Ax3(Loc, a, xD);
    const gp_Dir xDirAx = m.ax.XDirection();
    const gp_Dir yDirAx = m.ax.YDirection();
    const gp_XYZ xAx(xDirAx.X(), xDirAx.Y(), xDirAx.Z());
    const gp_XYZ yAx(yDirAx.X(), yDirAx.Y(), yDirAx.Z());

    std::vector<double> psi;
    for (int lt : tris) {
        const gp_Dir n = triNormalLocal(mv, lt);
        const double pu = n.Dot(xDirAx);
        const double pv = n.Dot(yDirAx);
        psi.push_back(wrapToPi(std::atan2(pv, pu)));
    }
    std::sort(psi.begin(), psi.end());

    int nBands = 1;
    if (psi.size() >= 2) {
        std::vector<double> gaps;
        for (size_t i = 1; i < psi.size(); i++) gaps.push_back(psi[i] - psi[i - 1]);
        gaps.push_back(2.0 * M_PI - psi.back() + psi.front());

        std::vector<double> gapInput = gaps;
        if (!gapInput.empty()) {
            const auto maxIt = std::max_element(gapInput.begin(), gapInput.end());
            gapInput.erase(maxIt);
        }
        const double thetaBin = std::max(tol.thetaBin, 0.5 * percentile75(gapInput));

        nBands = 1;
        for (size_t i = 1; i < psi.size(); i++) {
            if (psi[i] - psi[i - 1] > thetaBin) nBands++;
        }
    }
    m.nBands = nBands;

    std::vector<double> chi;
    chi.reserve(uniq.size());
    for (const auto& uv : uniq) {
        const gp_XYZ d = uv.second - loc;
        chi.push_back(wrapToPi(std::atan2(d.Dot(yAx), d.Dot(xAx))));
    }
    std::sort(chi.begin(), chi.end());

    const size_t mChi = chi.size();
    if (mChi < 2) {
        m.spanReject = true;
        return m;
    }

    // Circular gaps; ties on argmax -> lowest j (I5).
    int jmax = 0;
    double maxGapV = chi[1] - chi[0];
    for (size_t j = 1; j + 1 < mChi; j++) {
        const double g = chi[j + 1] - chi[j];
        if (g > maxGapV) {
            maxGapV = g;
            jmax = static_cast<int>(j);
        }
    }
    const double wrapGap = 2.0 * M_PI - chi.back() + chi.front();
    if (wrapGap > maxGapV) {
        maxGapV = wrapGap;
        jmax = static_cast<int>(mChi - 1);
    }

    const double bandArc = (nBands > 0) ? 2.0 * M_PI / static_cast<double>(nBands) : 2.0 * M_PI;
    m.closed360 = (nBands >= 3 && maxGapV <= 1.5 * bandArc);
    m.span = m.closed360 ? 2.0 * M_PI : (2.0 * M_PI - maxGapV);

    if (m.span <= 0 || nBands < 1) m.spanReject = true;
    else m.nSides = static_cast<int>(std::llround(2.0 * M_PI * nBands / m.span));

    if (m.closed360) {
        m.uMin = 0.0;
        m.uMax = 2.0 * M_PI;
    } else {
        const size_t uIdx = static_cast<size_t>(jmax + 1) % mChi;
        m.uMin = wrapToPi(chi[uIdx]);
        m.uMax = m.uMin + m.span;
    }

    m.vMin = std::numeric_limits<double>::infinity();
    m.vMax = -std::numeric_limits<double>::infinity();
    for (const auto& uv : uniq) {
        const double v = (uv.second - loc).Dot(aw);
        m.vMin = std::min(m.vMin, v);
        m.vMax = std::max(m.vMax, v);
    }
    return m;
}

bool computeOutwardCylinder(const MeshView& mv, const std::vector<int>& tris,
                            const gp_Dir& axis, const gp_Pnt& loc) {
    const gp_XYZ aw(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ lxyz(loc.X(), loc.Y(), loc.Z());
    double sigma = 0.0;
    for (int lt : tris) {
        const double a = triAreaLocal(mv, lt);
        const gp_XYZ cent = triCentroidLocal(mv, lt);
        const gp_XYZ d = cent - lxyz;
        const double ad = d.Dot(aw);
        gp_XYZ radial = d - aw * ad;
        const double rm = radial.Modulus();
        if (rm < kTiny) continue;
        radial /= rm;
        sigma += a * triNormalLocal(mv, lt).Dot(gp_Dir(radial));
    }
    return sigma > 0.0;
}

bool computeOutwardPlane(const MeshView& mv, const std::vector<int>& tris, const gp_Ax3& plane) {
    const gp_Dir n = plane.Direction();
    double sigma = 0.0;
    for (int lt : tris) {
        sigma += triAreaLocal(mv, lt) * triNormalLocal(mv, lt).Dot(n);
    }
    return sigma > 0.0;
}

double medianCentroidResidual(const MeshView& mv, const std::vector<int>& tris,
                              const gp_Dir& axis, const gp_Pnt& center, double radius) {
    const gp_XYZ a(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ c(center.X(), center.Y(), center.Z());
    std::vector<double> vals;
    for (int lt : tris) {
        const gp_XYZ cent = triCentroidLocal(mv, lt);
        const gp_XYZ d = cent - c;
        vals.push_back(std::abs(a.Crossed(d).Modulus() - radius));
    }
    return medianOf(vals);
}

// Order is binding (D1.3-A4): ConeNYI first, SphereNYI only when !c1.
Reject classifyG1Reject(const GaussResult& g, const DerivedTols& tol) {
    const double sin3 = tol.gaussAxisTiltSin();
    const bool c1 = g.flat < DerivedTols::kGaussPlanarity;
    const bool c2 = g.dev < sin3;
    if (c1 && c2 && std::abs(g.c) >= sin3) return Reject::ConeNYI;
    if (!c1 && g.patch >= 0.25) return Reject::SphereNYI;
    return Reject::GaussPlanarity;
}

struct ProvAdj {
    int other = -1;
    double phi = 0.0;
    int len = 0;
};

using ProvAdjList = std::vector<std::vector<ProvAdj>>;

ProvAdjList buildProvAdjacency(const MeshView& mv, const EdgeAdj& ea,
                               const std::vector<int>& triToProv) {
    int nProv = 0;
    for (int p : triToProv) if (p >= 0) nProv = std::max(nProv, p + 1);

    struct RawEdge { int lo, hi; double phi; };
    std::vector<RawEdge> raw;
    for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
        const int t0 = ea.tri[e][0];
        const int t1 = ea.tri[e][1];
        if (t0 < 0 || t1 < 0) continue;
        const int p0 = triToProv[t0];
        const int p1 = triToProv[t1];
        if (p0 < 0 || p1 < 0 || p0 == p1) continue;
        raw.push_back({ std::min(p0, p1), std::max(p0, p1), edgeDihedralAbs(mv, e, ea) });
    }
    std::sort(raw.begin(), raw.end(),
              [](const RawEdge& a, const RawEdge& b) {
                  if (a.lo != b.lo) return a.lo < b.lo;
                  if (a.hi != b.hi) return a.hi < b.hi;
                  return a.phi < b.phi;
              });

    ProvAdjList adj(nProv);
    size_t i = 0;
    while (i < raw.size()) {
        const int lo = raw[i].lo;
        const int hi = raw[i].hi;
        std::vector<double> phis;
        while (i < raw.size() && raw[i].lo == lo && raw[i].hi == hi) {
            phis.push_back(raw[i].phi);
            i++;
        }
        ProvAdj pa;
        pa.other = hi;
        pa.phi = medianOf(phis);
        pa.len = static_cast<int>(phis.size());
        adj[lo].push_back(pa);
        pa.other = lo;
        adj[hi].push_back(pa);
    }
    for (int p = 0; p < nProv; p++) {
        std::sort(adj[p].begin(), adj[p].end(),
                  [](const ProvAdj& a, const ProvAdj& b) { return a.other < b.other; });
    }
    return adj;
}

double phiToSet(const ProvAdjList& adj, int x, const std::vector<int>& members) {
    std::vector<double> phis;
    for (int m : members) {
        for (const ProvAdj& a : adj[x]) {
            if (a.other == m) phis.push_back(a.phi);
        }
    }
    return medianOf(phis);
}

int sharedLen(const ProvAdjList& adj, int x, const std::vector<int>& members) {
    int total = 0;
    for (int m : members) {
        for (const ProvAdj& a : adj[x]) {
            if (a.other == m) total += a.len;
        }
    }
    return total;
}

bool adjacentToSet(const ProvAdjList& adj, int x, const std::vector<int>& members) {
    for (int m : members) {
        for (const ProvAdj& a : adj[x]) {
            if (a.other == m) return true;
        }
    }
    return false;
}

bool seedInBand(double phi, const DerivedTols& tol) {
    const double loEps = angleBandEps(tol.thetaCylLo);
    const double hiEps = angleBandEps(tol.thetaCylHi);
    return phi >= tol.thetaCylLo - loEps && phi <= tol.thetaCylHi + hiEps;
}

bool membersEdgeConnected(const ProvAdjList& adj, const std::vector<int>& members) {
    if (members.size() <= 1) return true;
    const int n = static_cast<int>(adj.size());
    std::vector<char> inS(n, 0);
    for (int m : members) if (m >= 0 && m < n) inS[m] = 1;
    std::vector<char> seen(n, 0);
    std::vector<int> q;
    q.push_back(members[0]);
    seen[members[0]] = 1;
    size_t nseen = 1;
    for (size_t i = 0; i < q.size(); i++) {
        for (const ProvAdj& a : adj[q[i]]) {
            if (a.other < 0 || a.other >= n) continue;
            if (inS[a.other] && !seen[a.other]) {
                seen[a.other] = 1;
                q.push_back(a.other);
                nseen++;
            }
        }
    }
    return nseen == members.size();
}

void computeProvDeviations(const MeshView& mv, Provisional& p) {
    const gp_Dir n = p.plane.Direction();
    const gp_Pnt loc = p.plane.Location();
    double sumSq = 0.0;
    double maxD = 0.0;
    int nPts = 0;
    for (int lt : p.tris) {
        for (int k = 0; k < 3; k++) {
            const gp_XYZ v = localTriVert(mv, lt, k);
            const double d = std::abs(gp_Vec(loc, gp_Pnt(v)).Dot(gp_Vec(n)));
            maxD = std::max(maxD, d);
            sumSq += d * d;
            nPts++;
        }
    }
    p.maxVertexDev = maxD;
    p.rmsVertexDev = nPts > 0 ? std::sqrt(sumSq / nPts) : 0.0;
}

enum class Gate { G1, G2, G3, G4, G5, PASS };

const char* gateName(Gate g) {
    switch (g) {
    case Gate::G1: return "G1";
    case Gate::G2: return "G2";
    case Gate::G3: return "G3";
    case Gate::G4: return "G4";
    case Gate::G5: return "G5";
    case Gate::PASS: return "PASS";
    }
    return "?";
}

struct CommitEval {
    Gate failGate = Gate::G1;
    GaussResult g;
    gp_Pnt center;
    double radius = 0;
    D2Metrics d2;
    bool eberlyOk = false;
    bool archChain = false;  // in-band apply that passed the closure gate
};

CommitEval evaluateCommit(const MeshView& mv, const DerivedTols& tol,
                          const std::vector<int>& tris, const gp_Dir& axis,
                          double rHint = 0.0, bool lawBand = false) {
    CommitEval ev;
    ev.g = centeredGauss(mv, tris, axis, tol);
    double cTilt = 0.0, devTilt = 0.0;
    axisTiltStats(mv, tris, axis, cTilt, devTilt);
    // RULE 4.1d: lawBand radius comes from the inverse. Kasa/eberly is not
    // a radius source (refit_grow.cpp:773).
    if (!lawBand)
        ev.eberlyOk = eberlyCenterRadius(mv, tris, axis, ev.center, ev.radius);
    ev.failGate = Gate::PASS;
    if (!testG1CommitSeedAxis(ev.g, tol, cTilt, devTilt)) {
        ev.failGate = Gate::G1;
        return ev;
    }
    if (!lawBand && !ev.eberlyOk) {
        ev.failGate = Gate::G4;
        return ev;
    }
    ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
    const double rBeforeRefine = ev.radius;
    bool archChainApplied = false;
    // RULE 4.1d: do not sec-lift an already-inscribed lawBand.
    if (!lawBand && coarseFusionBand(mv) && ev.d2.nSides >= 3 && !ev.d2.spanReject) {
        refineCylinderRadius(mv, tris, axis, ev.center, ev.radius, ev.d2.nSides, ev.d2.span,
                             rHint);
        ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
    }
    if (!lawBand && coarseFusionBand(mv) && tris.size() >= 3 && ev.radius >= 15.0) {
        double chainR = ev.radius;
        double chainScore = 0.0;
        if (archChainRadiusFromPatch(mv, tris, axis, chainR, chainScore, ev.radius)
            && chainScore >= 0.35) {
            if (chainScore >= 0.85 && ev.radius > 15.0 && chainR < ev.radius * 0.58
                && chainR > ev.radius * 0.42)
                chainR = std::sqrt(chainR * ev.radius);
            const double rel = chainR / ev.radius;
            if (rel > 0.62 && rel < 1.15) {
                const double oldR = ev.radius;
                ev.radius = chainR;
                ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
                double g2Tol = tol.epsCylAccept(ev.radius);
                const int nEstSides =
                    std::max(ev.d2.nSides, std::max(6, static_cast<int>(tris.size())));
                g2Tol = std::max(g2Tol, chordSagitta(ev.radius, nEstSides));
                if (ev.radius > 15.0)
                    g2Tol = std::max(g2Tol, 0.12 * ev.radius);
                else
                    g2Tol = std::max(g2Tol, 0.08 * ev.radius);
                if (maxVertexResidual(mv, tris, axis, ev.center, ev.radius) > g2Tol
                    && chainScore < 0.85) {
                    ev.radius = oldR;
                    ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
                    diagArchChain("DEFER", tris.size(), chainScore, oldR, chainR,
                                  ev.radius, rel, "g2_residual");
                } else {
                    archChainApplied = true;
                    diagArchChain("COMMIT", tris.size(), chainScore, oldR, chainR,
                                  ev.radius, rel, "coarse_window");
                }
            } else if (chainScore > 0.0) {
                diagArchChain("DEFER", tris.size(), chainScore, ev.radius, chainR,
                              ev.radius, rel, "rel_window");
            }
        } else if (chainScore > 0.0) {
            diagArchChain("DEFER", tris.size(), chainScore, ev.radius, chainR,
                          ev.radius, 0.0, "score<0.35");
        }
    } else if (!lawBand && inArchChainCommitBand(mv) && tris.size() >= 3) {
        // Outside coarseFusionBand: commit only chainScore >= 0.85. Blend and
        // rel window are the same constants as the coarse path above.
        // Non-coarse G2 still runs (archChainApplied does not skip it), so a
        // residual miss must roll back — losing an eberly commit is a FAIL.
        double chainR = ev.radius;
        double chainScore = 0.0;
        if (archChainRadiusFromPatch(mv, tris, axis, chainR, chainScore, ev.radius)
            && chainScore >= 0.85) {
            if (chainScore >= 0.85 && ev.radius > 15.0 && chainR < ev.radius * 0.58
                && chainR > ev.radius * 0.42)
                chainR = std::sqrt(chainR * ev.radius);
            const double rel = chainR / ev.radius;
            if (rel > 0.62 && rel < 1.15) {
                const double oldR = ev.radius;
                const D2Metrics oldD2 = ev.d2;
                ev.radius = chainR;
                ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
                const double g2Tol = tol.epsCylAccept(ev.radius);
                bool g3fail = false;
                const double delta = chordSagitta(ev.radius, ev.d2.nSides);
                if (delta > tol.epsMesh) {
                    const double s =
                        medianCentroidResidual(mv, tris, axis, ev.center, ev.radius);
                    if (s < DerivedTols::kG3Lo * delta || s > DerivedTols::kG3Hi * delta)
                        g3fail = true;
                }
                if (maxVertexResidual(mv, tris, axis, ev.center, ev.radius) > g2Tol) {
                    ev.radius = oldR;
                    ev.d2 = oldD2;
                    diagArchChain("DEFER", tris.size(), chainScore, oldR, chainR,
                                  ev.radius, rel, "g2_residual");
                } else if (g3fail) {
                    ev.radius = oldR;
                    ev.d2 = oldD2;
                    diagArchChain("DEFER", tris.size(), chainScore, oldR, chainR,
                                  ev.radius, rel, "g3_fail");
                } else {
                    // Closure-aware (AC-SEAMFIX): stamp the apply so J6's
                    // sewability ledger can see the new edge classes. Keep
                    // the commit — the ledger (not a band disable) is what
                    // closes Body18. Unclosed applies still go through; the
                    // dual-face + collapseBanned path orphans them to polyline.
                    archChainApplied = true;
                    ev.archChain = true;
                    diagArchChain("COMMIT", tris.size(), chainScore, oldR, chainR,
                                  ev.radius, rel, "archband_window");
                }
            } else if (chainScore > 0.0) {
                diagArchChain("DEFER", tris.size(), chainScore, ev.radius, chainR,
                              ev.radius, rel, "rel_window");
            }
        } else if (chainScore > 0.0) {
            diagArchChain("DEFER", tris.size(), chainScore, ev.radius, chainR,
                          ev.radius, 0.0, "score<0.85");
        }
    }
    if (ev.d2.spanReject) {
        ev.failGate = Gate::G1;
        return ev;
    }
    const bool coarse = coarseFusionBand(mv);
    const double lift = std::max(0.0, ev.radius - rBeforeRefine);
    if (coarse) {
        // G2: vertices sit on a chordal ring; coarse N≥6 meshes need at least
        // chord-sagitta slack, not just 1%R (handle-lock R≈16–20 fails at |S|=2).
        double g2Tol = tol.epsCylAccept(ev.radius);
        const int nEstSides = std::max(ev.d2.nSides, std::max(6, (int)tris.size()));
        g2Tol = std::max(g2Tol, chordSagitta(ev.radius, nEstSides));
        if (lift > 0.0)
            g2Tol = std::max(g2Tol, lift * 1.15 + chordSagitta(ev.radius, nEstSides));
        if (tris.size() <= 8)
            g2Tol = std::max(g2Tol, epsCylRing(tol, ev.radius));
        if (!archChainApplied
            && maxVertexResidual(mv, tris, axis, ev.center, ev.radius) > g2Tol) {
            ev.failGate = Gate::G2;
            return ev;
        }
    } else if (archChainBand(mv)) {
        // D-130-1: G2 may use chord-sagitta slack only for a counted N-gon
        // (N = ev.d2.nSides, N >= 6). No tris.size() substitute, no literal
        // floor. Otherwise g2Tol is epsCylAccept(R). Coarse branch untouched.
        double g2Tol = tol.epsCylAccept(ev.radius);
        if (ev.d2.nSides >= 6)
            g2Tol = std::max(g2Tol, chordSagitta(ev.radius, ev.d2.nSides));
        if (maxVertexResidual(mv, tris, axis, ev.center, ev.radius) > g2Tol) {
            ev.failGate = Gate::G2;
            return ev;
        }
    } else if (maxVertexResidual(mv, tris, axis, ev.center, ev.radius)
               > tol.epsCylAccept(ev.radius)) {
        ev.failGate = Gate::G2;
        return ev;
    }
    const double delta = chordSagitta(ev.radius, ev.d2.nSides);
    // G3: nSides from a tiny arc span is unstable — skip until span ≥ ~20°.
    // Chord-corrected radii shift centroid residuals; skip G3 on coarse lift.
    if (delta > tol.epsMesh
        && (!coarse || (ev.d2.span >= 0.35 && lift <= 0.0 && !archChainApplied))) {
        const double s = medianCentroidResidual(mv, tris, axis, ev.center, ev.radius);
        if (s < DerivedTols::kG3Lo * delta || s > DerivedTols::kG3Hi * delta)
            ev.failGate = Gate::G3;
    }
    if (ev.failGate == Gate::PASS
        && !(2.0 * tol.epsPlane < ev.radius && ev.radius < 2.0 * mv.diag))
        ev.failGate = Gate::G4;
    if (ev.failGate == Gate::PASS) {
        const double spanDeg = ev.d2.span * 180.0 / M_PI;
        int nBandsUse = ev.d2.nBands;
        if (coarse && ev.d2.span < 0.5 && tris.size() <= 12)
            nBandsUse = std::max(nBandsUse, estimateNBandsFromPatch(mv, tris));
        const bool g5Closed = spanDeg >= DerivedTols::kG5SpanClosedDeg
                              && ev.d2.nSides >= DerivedTols::kG5NSidesMin;
        const double g5PartialDeg =
            coarse ? 30.0 : DerivedTols::kG5SpanPartialDeg;
        const bool g5Partial = spanDeg >= g5PartialDeg
                               && nBandsUse >= DerivedTols::kG5NBandsMin;
        const bool g5Micro = coarse && ev.radius >= 5.0 && tris.size() <= 8
                             && spanDeg >= 12.0 && nBandsUse >= 2;
        if (!g5Closed && !g5Partial && !g5Micro && !archChainApplied)
            ev.failGate = Gate::G5;
    }
    return ev;
}

void fillCylinderRegion(const MeshView& mv, const CommitEval& ev, const gp_Dir& axis,
                        const std::vector<int>& tris, Region& reg) {
    double areaReg = 0.0;
    for (int lt : tris) areaReg += triAreaLocal(mv, lt);
    reg.type = SurfType::Cylinder;
    reg.origin = Origin::CylGrow;
    reg.tris = tris;
    reg.ax = ev.d2.ax;
    reg.radius = ev.radius;
    reg.uMin = ev.d2.uMin;
    reg.uMax = ev.d2.uMax;
    reg.vMin = ev.d2.vMin;
    reg.vMax = ev.d2.vMax;
    reg.closed360 = ev.d2.closed360;
    reg.nSides = ev.d2.nSides;
    reg.chordSagitta = chordSagitta(ev.radius, ev.d2.nSides);
    reg.outwardNormal = computeOutwardCylinder(mv, tris, axis, ev.d2.ax.Location());
    reg.dVolPredicted = dVolCylinderSector(areaReg, ev.radius, ev.d2.nSides, reg.outwardNormal);
    reg.maxVertexDev = maxVertexResidual(mv, tris, axis, ev.center, ev.radius);
    double sumSq = 0.0;
    int nV = 0;
    const gp_XYZ cxyz(ev.center.X(), ev.center.Y(), ev.center.Z());
    const gp_XYZ axyz(axis.X(), axis.Y(), axis.Z());
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            const gp_XYZ p = localTriVert(mv, lt, k);
            const double rr = axyz.Crossed(p - cxyz).Modulus();
            const double d = rr - ev.radius;
            sumSq += d * d;
            nV++;
        }
    }
    reg.rmsVertexDev = nV > 0 ? std::sqrt(sumSq / nV) : 0.0;
    if (ev.archChain) {
        // Survive refit_chains remapFilletNbr (maps nb<0 → -1). Negative snap
        // is otherwise unused; build keys the extended J6 ledger off this.
        reg.maxVertexSnap = -1.0;
        reg.filletNbrA = kArchChainSewMark;
    }
}

double axisLineSeparation(const gp_Ax3& a, const gp_Ax3& b) {
    const gp_Vec d(a.Location(), b.Location());
    const gp_Vec ax(a.Direction());
    return ax.Crossed(d).Magnitude();
}

bool regionsShareMeshEdge(const MeshView& mv, const Region& a, const Region& b) {
    for (int t : a.tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        const int g = mv.compTris[t];
        const int* T = mv.tris[g];
        for (int k = 0; k < 3; k++) {
            const int gv0 = T[k];
            const int gv1 = T[(k + 1) % 3];
            for (int u : b.tris) {
                if (u < 0 || static_cast<size_t>(u) >= mv.nTri) continue;
                const int gu = mv.compTris[u];
                const int* U = mv.tris[gu];
                for (int j = 0; j < 3; j++) {
                    const int a0 = U[j];
                    const int a1 = U[(j + 1) % 3];
                    if ((a0 == gv0 && a1 == gv1) || (a0 == gv1 && a1 == gv0)) return true;
                }
            }
        }
    }
    return false;
}

bool coaxialCylinderMergeable(const Region& a, const Region& b, const DerivedTols& tol,
                                bool adjacent) {
    if (a.type != SurfType::Cylinder || b.type != SurfType::Cylinder) return false;
    if (!(a.radius > 0.0) || !(b.radius > 0.0)) return false;
    const gp_Dir da = a.ax.Direction();
    const gp_Dir db = b.ax.Direction();
    if (std::abs(da.Dot(db)) < 0.995) return false;
    const double sep = axisLineSeparation(a.ax, b.ax);
    const double rTol = std::max(tol.epsPlane, 0.08 * std::min(a.radius, b.radius));
    if (sep > rTol) return false;
    const double rAvg = 0.5 * (a.radius + b.radius);
    const double rRel = std::abs(a.radius - b.radius) / std::max(rAvg, 1e-6);
    if (rRel < 0.12) return true;
    // Phantom partial fits on the same bore: adjacent, coaxial, complementary radii.
    return adjacent && rRel < 0.28;
}

void mergeCoaxialCylinders(const MeshView& mv, const DerivedTols& tol, SegmentWork& work) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < work.accepted.size(); ++i) {
            Region& A = work.accepted[i];
            if (A.type != SurfType::Cylinder) continue;
            if (A.lawBand) continue;
            for (size_t j = i + 1; j < work.accepted.size(); ++j) {
                const Region& B = work.accepted[j];
                if (B.lawBand) continue;
                const bool adjacent = regionsShareMeshEdge(mv, A, B);
                if (!coaxialCylinderMergeable(A, B, tol, adjacent)) continue;
                // D-130-17: only edge-connected pieces of a surface are faces.
                // Two coaxial regions that share no mesh edge are two pieces
                // of one surface severed by another feature (cross_bores' R5
                // bore, cut in two by the R8 bore); merged, the region carried
                // two cap pairs and no face could be built from it (seamed360
                // failed on rid 4 and the whole part reverted). They stay two
                // regions; the census counts them as one surface.
                if (!adjacent) continue;
                std::vector<int> mergedTris = A.tris;
                mergedTris.insert(mergedTris.end(), B.tris.begin(), B.tris.end());
                std::sort(mergedTris.begin(), mergedTris.end());
                mergedTris.erase(std::unique(mergedTris.begin(), mergedTris.end()),
                                 mergedTris.end());
                const gp_Dir axis = A.radius >= B.radius ? A.ax.Direction() : B.ax.Direction();
                CommitEval ev = evaluateCommit(mv, tol, mergedTris, axis);
                if (ev.failGate != Gate::PASS) continue;
                Region reg;
                fillCylinderRegion(mv, ev, axis, mergedTris, reg);
                A = reg;
                work.accepted.erase(work.accepted.begin() + static_cast<long>(j));
                changed = true;
                break;
            }
            if (changed) break;
        }
    }
}

gp_Dir seedPairAxis(const Provisional& P, const Provisional& Q) {
    const gp_Dir nP = P.plane.Direction();
    const gp_Dir nQ = Q.plane.Direction();
    const gp_XYZ cross = gp_XYZ(nP.X(), nP.Y(), nP.Z())
                             .Crossed(gp_XYZ(nQ.X(), nQ.Y(), nQ.Z()));
    if (cross.SquareModulus() <= 1e-18) return gp_Dir(0, 0, 1);
    return canonicalAxis(cross);
}

// D1.3-A4c: seed-pair while |T|<=3 or Gauss degenerate; adopt w1 iff
// score(T,w1) <= score(T,seedPair), both recomputed every call. Not a freeze.
gp_Dir axisOf(const MeshView& mv, const std::vector<Provisional>& provs,
              const std::vector<int>& members, const gp_Dir& seedAxis,
              const DerivedTols& tol, GaussResult* gOut,
              bool* adoptedW1 = nullptr, double* scoreW1Out = nullptr,
              double* scoreSeedOut = nullptr) {
    const std::vector<int> tris = mergeMemberTris(provs, members);
    GaussResult g = centeredGauss(mv, tris, seedAxis, tol);
    if (gOut) *gOut = g;
    double scW1 = 0.0, scSeed = 0.0, cSeed = 0.0;
    axisTiltStats(mv, tris, seedAxis, cSeed, scSeed);
    if (g.ok && !g.degenerate) scW1 = g.dev;
    if (scoreW1Out) *scoreW1Out = scW1;
    if (scoreSeedOut) *scoreSeedOut = scSeed;
    if (members.size() <= 3 || g.degenerate || !g.ok) {
        if (adoptedW1) *adoptedW1 = false;
        return seedAxis;
    }
    if (scW1 <= scSeed) {
        if (adoptedW1) *adoptedW1 = true;
        return g.axis;
    }
    if (adoptedW1) *adoptedW1 = false;
    return seedAxis;
}

// Same-facet cluster vs the *seed* normal only (no running-mean drift).
// 2° first pass keeps staggered high-N bands (F515) intact; 1° leftover
// pass splits the F521/F520 1.67° junction the 2° gate swallows.
constexpr double kLawStripNormalCos = 0.9993908270190957;
constexpr double kLawStripNormalCosTight = 0.9998476951563913;
constexpr double kLawRRelGrow = 5e-4;

struct LawStrip {
    std::vector<int> tris;
    int minTri = 0;
    int chartId = -1;
    double R = 0.0;
    gp_Ax1 axis;
    bool valid = false;
};

void clusterLawStrips(const MeshView& mv, const std::vector<int>& ids,
                      int chartId, std::vector<LawStrip>& out,
                      double nCos = kLawStripNormalCos, bool seedOnly = false) {
    if (ids.empty()) return;
    std::vector<char> inPatch(static_cast<size_t>(mv.nTri), 0);
    for (int t : ids) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inPatch[static_cast<size_t>(t)] = 1;
    }
    std::vector<char> used(static_cast<size_t>(mv.nTri), 0);
    const EdgeAdj ea = buildEdgeAdj(mv);
    for (int seed : ids) {
        if (seed < 0 || static_cast<size_t>(seed) >= mv.nTri || used[static_cast<size_t>(seed)])
            continue;
        std::vector<int> comp;
        std::vector<int> stack = {seed};
        gp_XYZ nRef(triNormalLocal(mv, seed).X(), triNormalLocal(mv, seed).Y(),
                    triNormalLocal(mv, seed).Z());
        const gp_XYZ nSeed = nRef;
        used[static_cast<size_t>(seed)] = 1;
        comp.push_back(seed);
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int s = 0; s < 3; s++) {
                const int e = mv.triEdges[t][s];
                const int u = (ea.tri[e][0] == t) ? ea.tri[e][1] : ea.tri[e][0];
                if (u < 0 || !inPatch[static_cast<size_t>(u)] || used[static_cast<size_t>(u)])
                    continue;
                const gp_Dir nu = triNormalLocal(mv, u);
                const gp_XYZ nU(nu.X(), nu.Y(), nu.Z());
                const gp_XYZ nGate = seedOnly ? nSeed : nRef;
                if (std::abs(nGate.Dot(nU)) < nCos) continue;
                used[static_cast<size_t>(u)] = 1;
                comp.push_back(u);
                stack.push_back(u);
                if (!seedOnly) {
                    nRef += nU;
                    const double nm = nRef.Modulus();
                    if (nm > 1e-15) nRef.Divide(nm);
                }
            }
        }
        std::sort(comp.begin(), comp.end());
        LawStrip st;
        st.tris = std::move(comp);
        st.minTri = st.tris.front();
        st.chartId = chartId;
        out.push_back(std::move(st));
    }
}

bool stripsShareEdge(const MeshView& mv, const LawStrip& a, const LawStrip& b) {
    std::vector<char> inB(static_cast<size_t>(mv.nTri), 0);
    for (int t : b.tris) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inB[static_cast<size_t>(t)] = 1;
    }
    const EdgeAdj ea = buildEdgeAdj(mv);
    for (int t : a.tris) {
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            const int u = (ea.tri[e][0] == t) ? ea.tri[e][1] : ea.tri[e][0];
            if (u >= 0 && inB[static_cast<size_t>(u)]) return true;
        }
    }
    return false;
}

bool bandsShareAnyEdge(const MeshView& mv, const LawBand& a, const LawBand& b) {
    std::vector<char> inB(static_cast<size_t>(mv.nTri), 0);
    for (int t : b.tris) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inB[static_cast<size_t>(t)] = 1;
    }
    const EdgeAdj ea = buildEdgeAdj(mv);
    for (int t : a.tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            const int u = (ea.tri[e][0] == t) ? ea.tri[e][1] : ea.tri[e][0];
            if (u >= 0 && inB[static_cast<size_t>(u)]) return true;
        }
    }
    return false;
}

bool axesNearlyCoaxial(const gp_Ax1& a, const gp_Ax1& b, double lineTol) {
    if (std::abs(a.Direction().Dot(b.Direction())) < 1.0 - 1e-3) return false;
    gp_Vec delta(b.Location(), a.Location());
    gp_Vec perp = delta - delta.Dot(gp_Vec(b.Direction())) * gp_Vec(b.Direction());
    return perp.Magnitude() <= lineTol;
}

void peelLawBandFromProvisionals(const MeshView& mv, const std::vector<int>& claimed,
                                 SegmentWork& work) {
    std::vector<char> taken(static_cast<size_t>(mv.nTri), 0);
    for (int t : claimed) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) taken[static_cast<size_t>(t)] = 1;
    }
    for (Provisional& p : work.provisionals) {
        if (p.tris.empty()) continue;
        std::vector<int> keep;
        keep.reserve(p.tris.size());
        bool lost = false;
        for (int t : p.tris) {
            if (t >= 0 && static_cast<size_t>(t) < taken.size() && taken[static_cast<size_t>(t)]) {
                lost = true;
                continue;
            }
            keep.push_back(t);
        }
        if (!lost) continue;
        p.tris = std::move(keep);
        if (p.tris.empty()) {
            p.area = 0.0;
            p.claim = ProvClaim::ConsumedCylinder;
            continue;
        }
        gp_Ax3 fit;
        if (pcaPlane(mv, p.tris, fit)) p.plane = fit;
        p.area = 0.0;
        for (int t : p.tris) p.area += triAreaLocal(mv, t);
        computeProvDeviations(mv, p);
    }
}

void fillLawBandRegion(const MeshView& mv, const DerivedTols& tol, const LawBand& band,
                       Region& reg) {
    const gp_Dir axis = band.axis.Direction();
    const gp_Pnt center = band.axis.Location();
    CommitEval ev;
    ev.failGate = Gate::PASS;
    ev.eberlyOk = true;
    ev.center = center;
    ev.radius = band.R;
    ev.d2 = computeD2(mv, band.tris, axis, center, band.R, tol);
    if (band.closed360) {
        ev.d2.closed360 = true;
        ev.d2.uMin = 0.0;
        ev.d2.uMax = 2.0 * M_PI;
        ev.d2.span = 2.0 * M_PI;
    }
    fillCylinderRegion(mv, ev, axis, band.tris, reg);
    reg.lawBand = true;
    reg.radius = band.R;
    reg.closed360 = band.closed360 || reg.closed360;
}

}  // namespace

bool chartsA1(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
              SegmentWork& work) {
    work.triChart.assign(mv.nTri, -1);
    work.nCharts = 0;
    if (mv.nTri == 0) return true;

    const EdgeAdj ea = buildEdgeAdj(mv);
    UnionFind uf(static_cast<int>(mv.nTri));
    for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
        const int t0 = ea.tri[e][0];
        const int t1 = ea.tri[e][1];
        if (t0 < 0 || t1 < 0) continue;
        const double phi = edgeDihedralAbs(mv, e, ea);
        if (phi <= tol.thetaSharp) uf.unite(t0, t1);
    }

    std::vector<std::vector<int>> components;
    std::vector<int> rootMap(static_cast<int>(mv.nTri), -1);
    for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
        const int r = uf.find(t);
        if (rootMap[r] < 0) {
            rootMap[r] = static_cast<int>(components.size());
            components.emplace_back();
        }
        components[rootMap[r]].push_back(t);
    }

    std::vector<std::pair<int, int>> chartOrder;
    for (size_t i = 0; i < components.size(); i++) {
        int lo = components[i][0];
        for (int t : components[i]) lo = std::min(lo, t);
        chartOrder.emplace_back(lo, static_cast<int>(i));
    }
    std::sort(chartOrder.begin(), chartOrder.end());

    work.nCharts = static_cast<int>(chartOrder.size());
    std::vector<int> remap(components.size());
    for (int i = 0; i < work.nCharts; i++) remap[chartOrder[i].second] = i;

    for (size_t ci = 0; ci < components.size(); ci++) {
        const int chartId = remap[static_cast<int>(ci)];
        for (int t : components[ci]) work.triChart[t] = chartId;
    }
    return true;
}

bool growProvisionalA2(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
                       SegmentWork& work) {
    work.provisionals.clear();
    if (mv.nTri == 0) return true;
    const bool a2Diag = a2DiagOn();

    const EdgeAdj ea = buildEdgeAdj(mv);
    std::vector<int> triLabel(mv.nTri, -1);

    std::vector<std::vector<std::pair<int, int>>> triNeighbors(mv.nTri);
    for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
        const int t0 = ea.tri[e][0];
        const int t1 = ea.tri[e][1];
        if (t0 < 0 || t1 < 0) continue;
        triNeighbors[t0].emplace_back(e, t1);
        triNeighbors[t1].emplace_back(e, t0);
    }
    for (auto& nb : triNeighbors) {
        std::sort(nb.begin(), nb.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
    }

    const int nCharts = work.nCharts > 0 ? work.nCharts : 1;
    for (int chart = 0; chart < nCharts; chart++) {
        while (true) {
            int seed = -1;
            double bestArea = -1.0;
            for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
                if (triLabel[t] >= 0) continue;
                if (work.triChart[t] != chart) continue;
                if (triInLawBand(work, t)) continue;
                const double a = triAreaLocal(mv, t);
                if (a > bestArea || (a == bestArea && t < seed)) {
                    bestArea = a;
                    seed = t;
                }
            }
            if (seed < 0) break;

            Provisional prov;
            prov.chartId = chart;
            std::vector<int> queue;
            queue.push_back(seed);
            triLabel[seed] = static_cast<int>(work.provisionals.size());

            const gp_XYZ seedC = triCentroidLocal(mv, seed);
            const gp_Dir seedNormal = triNormalLocal(mv, seed);
            gp_Ax3 plane(gp_Pnt(seedC), seedNormal);
            std::vector<int> growTris;

            while (!queue.empty()) {
                const int t = queue.front();
                queue.erase(queue.begin());
                growTris.push_back(t);

                gp_Ax3 refreshed;
                if (growTris.size() > 2 && pcaPlane(mv, growTris, refreshed)) {
                    const gp_Dir pcaN = refreshed.Direction();
                    if (std::abs(seedNormal.Dot(pcaN)) >= std::cos(tol.thetaPlane))
                        plane = refreshed;
                    else
                        plane = gp_Ax3(refreshed.Location(), seedNormal);
                }

                for (const auto& nb : triNeighbors[t]) {
                    const int e = nb.first;
                    const int u = nb.second;
                    // Diag only: name the predicate that refuses this link.
                    // Every `continue` below prints once and changes nothing.
                    const bool aDiag = a2Diag;
                    auto arej = [&](const char* why, double a, double b) {
                        if (aDiag)
                            std::fprintf(stderr,
                                         "DIAG_A2_LINK p=%zu t=%d u=%d e=%d why=%s "
                                         "a=%.9g b=%.9g\n",
                                         work.provisionals.size(), t, u, e, why, a, b);
                    };
                    if (triLabel[u] >= 0) {
                        arej("already-labelled", (double)triLabel[u], -1.0);
                        continue;
                    }
                    if (work.triChart[u] != chart) {
                        arej("other-chart", (double)work.triChart[u], (double)chart);
                        continue;
                    }
                    if (triInLawBand(work, u)) {
                        arej("in-law-band", 1.0, 0.0);
                        continue;
                    }
                    const double phi = edgeDihedralAbs(mv, e, ea);
                    if (phi > tol.thetaSharp) {
                        arej("sharp-edge", phi, tol.thetaSharp);
                        continue;
                    }

                    const gp_Dir n_u = triNormalLocal(mv, u);
                    const gp_Dir n_p = plane.Direction();
                    if (std::abs(n_u.Dot(n_p)) < std::cos(tol.thetaPlane)) {
                        arej("near-flat-normal-gate",
                             std::acos(std::min(1.0, std::abs(n_u.Dot(n_p)))),
                             tol.thetaPlane);
                        continue;
                    }

                    double maxD = 0.0;
                    const gp_Pnt loc = plane.Location();
                    for (int k = 0; k < 3; k++) {
                        const gp_XYZ v = localTriVert(mv, u, k);
                        const double d = std::abs(
                            gp_Vec(loc, gp_Pnt(v)).Dot(gp_Vec(n_p)));
                        maxD = std::max(maxD, d);
                    }
                    if (maxD > tol.epsPlane) {
                        arej("plane-residual", maxD, tol.epsPlane);
                        continue;
                    }

                    std::vector<int> trial = growTris;
                    trial.push_back(u);
                    gp_Ax3 trialPlane;
                    if (!pcaPlane(mv, trial, trialPlane)) trialPlane = plane;
                    else {
                        const gp_Dir pcaN = trialPlane.Direction();
                        if (std::abs(seedNormal.Dot(pcaN)) < std::cos(tol.thetaPlane))
                            trialPlane = gp_Ax3(trialPlane.Location(), seedNormal);
                    }

                    triLabel[u] = static_cast<int>(work.provisionals.size());
                    queue.push_back(u);
                    plane = trialPlane;
                }
            }

            std::sort(growTris.begin(), growTris.end());
            prov.tris = growTris;
            prov.plane = plane;
            prov.area = 0.0;
            for (int lt : prov.tris) prov.area += triAreaLocal(mv, lt);
            computeProvDeviations(mv, prov);
            work.provisionals.push_back(prov);
        }
    }

    // CYLCYL -- a sliver filed into the wrong facet's patch.
    //
    // A2 admits a candidate on two tests: its normal within thetaPlane of the
    // patch plane and its vertices within epsPlane of it. A THIN triangle
    // passes the second test against any plane that contains its long edge,
    // whatever its normal says: over a 0.05 mm extent a 7.5 deg fold is a
    // 0.006 mm residual, under epsPlane's 0.02. In the coarse band, where the
    // normal gate is widened to 15 deg for noisy exports, such a sliver is
    // therefore taken by whichever patch floods it first. On cross_bores the
    // R8 wall's notched facets triangulate into slivers along the generator
    // lines and each is filed into the NEIGHBOURING facet's patch; every mesh
    // edge two adjacent facets then share is a coplanar edge inside the
    // misfiled sliver, the 7.5 deg fold edge is internal to one patch, and B1
    // reads a 0 deg fold between two facets that meet at 7.5 deg -- the
    // notched pieces are refused as flat and the wall ships in three pieces.
    //
    // A triangle lies on the facet whose plane its own normal matches. The
    // refile is decided by two measurements and no threshold: (1) the triangle
    // is exactly coplanar with an edge-neighbour in another patch of the same
    // chart -- each one's far vertex lies on the other's plane within the mesh
    // FILE's own quantization q (D-130-12), the resolution below which the
    // mesh cannot state a fold at all; (2) its normal is strictly closer to
    // that patch's plane than to its own. A resolved fold between two patches
    // fails (1); a planar face floods whole so its triangles never meet (1)
    // across a patch boundary. Prints under STL2STEP_DIAG_A2.
    {
        const double q = (std::isfinite(mv.quantFloor) && mv.quantFloor > 0.0)
                             ? mv.quantFloor : 0.0;
        int movedTotal = 0, rounds = 0;
        if (q > 0.0 && !work.provisionals.empty()) {
            std::vector<int> triToProv(mv.nTri, -1);
            for (size_t pi = 0; pi < work.provisionals.size(); pi++)
                for (int t : work.provisionals[pi].tris) triToProv[(size_t)t] = (int)pi;
            auto farVertexOffPlane = [&](int t, int u) {
                // max distance of t's vertices from the plane of u (u's own
                // vertices, u's own normal): the two shared vertices read 0,
                // the far one reads the fold the mesh resolves at this edge.
                const gp_Dir nu = triNormalLocal(mv, u);
                const gp_XYZ pu = localTriVert(mv, u, 0);
                double m = 0.0;
                for (int k = 0; k < 3; k++)
                    m = std::max(m, std::abs((localTriVert(mv, t, k) - pu).Dot(nu.XYZ())));
                return m;
            };
            // The fold at an edge is read on the plane of the LARGER of its
            // two triangles: a sliver's own plane is set by vertices a few
            // ulps apart across its width and extrapolates that rounding
            // twenty-fold at its neighbour's far vertex, so a resolved fold
            // and an exactly coplanar pair would both read "off". The large
            // triangle's plane carries its vertices' rounding, q, and no more;
            // the sliver's far vertex is within its own width of that
            // triangle's edge, so nothing is extrapolated. A fold the mesh
            // resolves puts that vertex w*sin(phi) off -- 0.0065 mm for the
            // 0.05 mm slivers on cross_bores' R8 wall at 7.5 deg, 2000 q.
            auto coplanarAtEdge = [&](int t0, int t1) {
                const int big = triAreaLocal(mv, t0) >= triAreaLocal(mv, t1) ? t0 : t1;
                const int small = big == t0 ? t1 : t0;
                return farVertexOffPlane(small, big) <= q;
            };
            auto angleToPlane = [&](int t, const Provisional& P) {
                const double c = std::min(1.0, std::abs(triNormalLocal(mv, t).Dot(
                                                    P.plane.Direction())));
                return std::acos(c);
            };
            for (rounds = 0; rounds < 8; rounds++) {
                int moved = 0;
                std::vector<char> dirty(work.provisionals.size(), 0);
                for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
                    const int t0 = ea.tri[e][0];
                    const int t1 = ea.tri[e][1];
                    if (t0 < 0 || t1 < 0) continue;
                    const int p0 = triToProv[(size_t)t0];
                    const int p1 = triToProv[(size_t)t1];
                    if (p0 < 0 || p1 < 0 || p0 == p1) continue;
                    if (work.provisionals[(size_t)p0].chartId !=
                        work.provisionals[(size_t)p1].chartId)
                        continue;
                    if (!coplanarAtEdge(t0, t1)) continue;  // a resolved fold, not a split facet
                    for (int side = 0; side < 2; side++) {
                        const int t = side == 0 ? t0 : t1;
                        const int from = side == 0 ? p0 : p1;
                        const int to = side == 0 ? p1 : p0;
                        if (triToProv[(size_t)t] != from) continue;
                        Provisional& P = work.provisionals[(size_t)from];
                        Provisional& Q = work.provisionals[(size_t)to];
                        const double aOwn = angleToPlane(t, P);
                        const double aOther = angleToPlane(t, Q);
                        if (!(aOther < aOwn)) continue;
                        if (a2Diag)
                            std::fprintf(stderr,
                                         "DIAG_A2_REFILE t=%d from=%d to=%d aOwn=%.9g "
                                         "aOther=%.9g q=%.3g\n",
                                         t, from, to, aOwn, aOther, q);
                        P.tris.erase(std::remove(P.tris.begin(), P.tris.end(), t),
                                     P.tris.end());
                        Q.tris.insert(std::upper_bound(Q.tris.begin(), Q.tris.end(), t), t);
                        triToProv[(size_t)t] = to;
                        dirty[(size_t)from] = dirty[(size_t)to] = 1;
                        moved++;
                    }
                }
                for (size_t pi = 0; pi < work.provisionals.size(); pi++) {
                    if (!dirty[pi]) continue;
                    Provisional& P = work.provisionals[pi];
                    P.area = 0.0;
                    for (int lt : P.tris) P.area += triAreaLocal(mv, lt);
                    gp_Ax3 fit;
                    if (!P.tris.empty() && pcaPlane(mv, P.tris, fit)) P.plane = fit;
                    computeProvDeviations(mv, P);
                }
                movedTotal += moved;
                if (moved == 0) break;
            }
        }
        if (a2Diag)
            std::fprintf(stderr, "DIAG_A2_REFILE_SUM moved=%d rounds=%d q=%.3g\n",
                         movedTotal, rounds, q);
    }

    sortProvisionals(work.provisionals);
    if (a2Diag) {
        std::fprintf(stderr, "DIAG_A2_PROV n=%zu\n", work.provisionals.size());
        for (size_t i = 0; i < work.provisionals.size(); i++) {
            const Provisional& pr = work.provisionals[i];
            const gp_Dir nd = pr.plane.Direction();
            std::fprintf(stderr,
                         "  DIAG_A2_ONE p=%zu nTri=%zu minTri=%d chart=%d "
                         "n=(%.6f,%.6f,%.6f) area=%.6g\n",
                         i, pr.tris.size(),
                         pr.tris.empty() ? -1 : *std::min_element(pr.tris.begin(),
                                                                 pr.tris.end()),
                         pr.chartId, nd.X(), nd.Y(), nd.Z(), pr.area);
        }
    }
    return true;
}

bool claimCylindersB1(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
                      SegmentWork& work) {
    if (lawbandDiagOn())
        std::fprintf(stderr, "DIAG_B1_ENTER nProv=%zu nAcc=%zu\n",
                     work.provisionals.size(), work.accepted.size());
    for (Provisional& p : work.provisionals) {
        std::vector<int> ok;
        ok.reserve(p.tris.size());
        for (int t : p.tris) {
            if (t >= 0 && static_cast<size_t>(t) < mv.nTri) ok.push_back(t);
        }
        p.tris.swap(ok);
        if (p.tris.empty()) p.claim = ProvClaim::ConsumedCylinder;
    }
    if (work.provisionals.empty()) return true;

    if (lawbandDiagOn()) std::fprintf(stderr, "DIAG_B1_SANITIZED\n");
    const EdgeAdj ea = buildEdgeAdj(mv);
    if (lawbandDiagOn()) std::fprintf(stderr, "DIAG_B1_EDGEADJ\n");
    std::vector<int> triToProv(mv.nTri, -1);
    for (size_t pi = 0; pi < work.provisionals.size(); pi++) {
        for (int t : work.provisionals[pi].tris) {
            if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                triToProv[static_cast<size_t>(t)] = static_cast<int>(pi);
        }
    }
    if (lawbandDiagOn()) std::fprintf(stderr, "DIAG_B1_TRIPROV\n");

    ProvAdjList adj = buildProvAdjacency(mv, ea, triToProv);
    if (adj.size() < work.provisionals.size()) adj.resize(work.provisionals.size());
    if (lawbandDiagOn()) std::fprintf(stderr, "DIAG_B1_ADJ n=%zu\n", adj.size());
    const double lateralSin = std::sin(tol.thetaSharp);

    struct Seed {
        int p, q;
        double negAreaSum;
        int minTri;
        int maxTri;
        int adjLo, adjHi;
    };
    std::vector<Seed> seeds;
    for (size_t i = 0; i < work.provisionals.size(); i++) {
        for (const ProvAdj& a : adj[static_cast<int>(i)]) {
            const int j = a.other;
            if (static_cast<size_t>(j) <= i) continue;
            if (!seedInBand(a.phi, tol)) continue;
            const Provisional& P = work.provisionals[i];
            const Provisional& Q = work.provisionals[j];
            if (P.seedTried || Q.seedTried) continue;
            if (P.claim != ProvClaim::Unclaimed || Q.claim != ProvClaim::Unclaimed) continue;
            Seed s;
            s.p = static_cast<int>(i);
            s.q = j;
            s.negAreaSum = -(P.area + Q.area);
            const int mtP = minTriId(P);
            const int mtQ = minTriId(Q);
            s.minTri = std::min(mtP, mtQ);
            s.maxTri = std::max(mtP, mtQ);
            s.adjLo = std::min(s.p, s.q);
            s.adjHi = std::max(s.p, s.q);
            seeds.push_back(s);
        }
    }
    std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
        if (a.negAreaSum != b.negAreaSum) return a.negAreaSum < b.negAreaSum;
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        if (a.maxTri != b.maxTri) return a.maxTri < b.maxTri;
        if (a.adjLo != b.adjLo) return a.adjLo < b.adjLo;
        return a.adjHi < b.adjHi;
    });

    for (const Seed& seed : seeds) {
        Provisional& P = work.provisionals[seed.p];
        Provisional& Q = work.provisionals[seed.q];
        if (P.claim != ProvClaim::Unclaimed || Q.claim != ProvClaim::Unclaimed) continue;
        if (P.seedTried || Q.seedTried) continue;

        const gp_Dir seedAxis = seedPairAxis(P, Q);

        std::vector<int> members = { seed.p, seed.q };
        std::vector<char> inClaim(work.provisionals.size(), 0);
        std::vector<char> dead(work.provisionals.size(), 0);
        inClaim[seed.p] = 1;
        inClaim[seed.q] = 1;
        P.claim = ProvClaim::InCylinderClaim;
        Q.claim = ProvClaim::InCylinderClaim;

        double R_ref = 0.0;
        bool haveRref = false;
        bool deadCleared = false;
        int nG5 = 0;
        double worstG5 = 0.0;

        if (p1DiagOn()) {
            std::fprintf(stderr, "B1 seed P=%d Q=%d minTri=%d aP=%.5g aQ=%.5g nP=%zu nQ=%zu\n",
                         seed.p, seed.q, seed.minTri,
                         P.area, Q.area, P.tris.size(), Q.tris.size());
        }

        while (true) {
            GaussResult gS;
            bool usedW1 = false;
            double scW1 = 0.0, scSeed = 0.0;
            const gp_Dir aS = axisOf(mv, work.provisionals, members, seedAxis, tol, &gS,
                                     &usedW1, &scW1, &scSeed);
            (void)gS;
            if (!deadCleared && usedW1) {
                std::fill(dead.begin(), dead.end(), 0);
                deadCleared = true;
                if (p1DiagOn())
                    std::fprintf(stderr,
                                 "  A4c adopt w1 |S|=%zu scW1=%.4g scSeed=%.4g (dead clear)\n",
                                 members.size(), scW1, scSeed);
            }
            const std::vector<int> sTrisNow = mergeMemberTris(work.provisionals, members);
            double cS = 0.0, devS = 0.0;
            axisTiltStats(mv, sTrisNow, aS, cS, devS);
            (void)devS;
            const double sin3 = tol.gaussAxisTiltSin();
            std::vector<double> memAreas;
            memAreas.reserve(members.size());
            for (int m : members)
                memAreas.push_back(work.provisionals[m].area);
            const double medArea = medianOf(memAreas);
            if (p1DiagOn()) {
                std::fprintf(stderr,
                             "  grow |S|=%zu axis=%s scW1=%.4g scSeed=%.4g cS=%.4g "
                             "a=(%.3g,%.3g,%.3g)\n",
                             members.size(), usedW1 ? "w1" : "seed", scW1, scSeed, cS,
                             aS.X(), aS.Y(), aS.Z());
            }

            struct GrowCand { int x; int negLen; int minTri; };
            std::vector<GrowCand> C;
            for (size_t xi = 0; xi < work.provisionals.size(); xi++) {
                if (inClaim[xi] || dead[xi]) continue;
                const Provisional& X = work.provisionals[xi];
                if (X.claim != ProvClaim::Unclaimed) continue;
                if (!adjacentToSet(adj, static_cast<int>(xi), members)) continue;
                const double phi = phiToSet(adj, static_cast<int>(xi), members);
                if (phi < tol.thetaCylLo - angleBandEps(tol.thetaCylLo)) continue;
                const gp_Dir nX = X.plane.Direction();
                if (std::abs(nX.Dot(aS)) > lateralSin) continue;
                // Tangent cube face (nSides=20, phi=9° > theta_cyl_lo) sits on
                // the Gauss great circle, so (g5) cannot exclude it. It is ~25×
                // a fillet band; skip scale-mismatched X (same 0.05 used by T1).
                if (medArea > 0.0
                    && X.area * DerivedTols::kGaussPlanarity > medArea)
                    continue;
                // D1.3-A1c (g5): membership vs current axis and the set's own offset.
                const gp_XYZ nXbar = areaWeightedNbar(mv, X.tris);
                const double g5Bound = sin3 + angleBandEps(sin3);
                const double g5v = std::abs(
                    nXbar.Dot(gp_XYZ(aS.X(), aS.Y(), aS.Z())) - cS);
                if (g5v > g5Bound) {
                    dead[xi] = 1;
                    nG5++;
                    worstG5 = std::max(worstG5, g5v);
                    if (p1DiagOn())
                        std::fprintf(stderr, "  g5 reject x=%d |n.a-c|=%.17g bound=%.17g sin3=%.17g\n",
                                     (int)xi, g5v, g5Bound, sin3);
                    continue;
                }
                GrowCand gc;
                gc.x = static_cast<int>(xi);
                gc.negLen = -sharedLen(adj, static_cast<int>(xi), members);
                gc.minTri = minTriId(X);
                C.push_back(gc);
            }
            if (C.empty()) {
                if (p1DiagOn()) {
                    std::fprintf(stderr, "  C empty |S|=%zu aS=(%.3g,%.3g,%.3g)\n",
                                 members.size(), aS.X(), aS.Y(), aS.Z());
                    for (size_t xi = 0; xi < work.provisionals.size(); xi++) {
                        const Provisional& X = work.provisionals[xi];
                        if (X.claim != ProvClaim::Unclaimed) continue;
                        if (!adjacentToSet(adj, static_cast<int>(xi), members)) continue;
                        const double phi = phiToSet(adj, static_cast<int>(xi), members);
                        const gp_Dir nX = X.plane.Direction();
                        std::fprintf(stderr,
                                     "    uncl x=%d phi=%.3g |n.a|=%.3g inClaim=%d dead=%d\n",
                                     (int)xi, phi, std::abs(nX.Dot(aS)),
                                     (int)inClaim[xi], (int)dead[xi]);
                    }
                }
                break;
            }
            std::sort(C.begin(), C.end(), [](const GrowCand& a, const GrowCand& b) {
                if (a.negLen != b.negLen) return a.negLen < b.negLen;
                return a.minTri < b.minTri;
            });

            bool progressed = false;
            for (const GrowCand& gc : C) {
                const int xi = gc.x;
                std::vector<int> trialMembers = members;
                trialMembers.push_back(xi);
                const std::vector<int> U = mergeMemberTris(work.provisionals, trialMembers);
                GaussResult gU;
                const gp_Dir aU = axisOf(mv, work.provisionals, trialMembers, seedAxis, tol, &gU);
                if (!testT1Running(gU)) {
                    if (p1DiagOn())
                        std::fprintf(stderr, "  T1 reject x=%d\n", xi);
                    dead[xi] = 1;
                    continue;
                }
                gp_Pnt center;
                double radius = 0.0;
                if (!eberlyCenterRadius(mv, U, aU, center, radius)
                    || !(radius > 0.0 && radius < 2.0 * mv.diag)) {
                    if (p1DiagOn())
                        std::fprintf(stderr, "  T2 reject x=%d\n", xi);
                    dead[xi] = 1;
                    continue;
                }
                const double t3R = haveRref ? R_ref : radius;
                const double resU = maxVertexResidual(mv, U, aU, center, radius);
                if (resU > epsCylRing(tol, t3R)) {
                    if (p1DiagOn())
                        std::fprintf(stderr,
                                     "  T3 reject x=%d R_U=%.5g R_ref=%.5g res=%.5g bound=%.5g |U|=%zu\n",
                                     xi, radius, t3R, resU, epsCylRing(tol, t3R),
                                     trialMembers.size());
                    dead[xi] = 1;
                    continue;
                }
                if (haveRref && std::abs(radius - R_ref) > kRingResidualFrac * R_ref) {
                    if (p1DiagOn())
                        std::fprintf(stderr, "  T3b reject x=%d R_U=%.5g R_ref=%.5g\n",
                                     xi, radius, R_ref);
                    dead[xi] = 1;
                    continue;
                }
                members.push_back(xi);
                inClaim[xi] = 1;
                work.provisionals[xi].claim = ProvClaim::InCylinderClaim;
                if (p1DiagOn())
                    std::fprintf(stderr, "  add x=%d area=%.5g nTri=%zu R=%.5g\n",
                                 xi, work.provisionals[xi].area,
                                 work.provisionals[xi].tris.size(), radius);
                if (!haveRref && members.size() == 3) {
                    R_ref = radius;
                    haveRref = true;
                } else if (haveRref && radius > R_ref) {
                    R_ref = radius;
                }
                progressed = true;
                break;
            }
            if (!progressed) break;
        }

        const std::vector<int> prePeelTris = mergeMemberTris(work.provisionals, members);
        GaussResult gFinal;
        const gp_Dir axisFinal = axisOf(mv, work.provisionals, members, seedAxis, tol, &gFinal);
        const double growHint = (haveRref && R_ref > 0.0) ? R_ref : 0.0;
        CommitEval ev = evaluateCommit(mv, tol, prePeelTris, axisFinal, growHint);
        if (p1DiagOn()) {
            gp_Pnt cSeed, cW1;
            double rSeed = 0, rW1 = 0;
            GaussResult gW1 = centeredGauss(mv, prePeelTris, seedAxis, tol);
            const bool eSeed = eberlyCenterRadius(mv, prePeelTris, seedAxis, cSeed, rSeed);
            const bool eW1 = eberlyCenterRadius(mv, prePeelTris, gW1.axis, cW1, rW1);
            double cTilt = 0, devTilt = 0;
            axisTiltStats(mv, prePeelTris, axisFinal, cTilt, devTilt);
            bool finW1 = false;
            double finScW1 = 0.0, finScSeed = 0.0;
            axisOf(mv, work.provisionals, members, seedAxis, tol, nullptr,
                   &finW1, &finScW1, &finScSeed);
            const double relR = (R_ref > 0.0)
                ? std::abs(ev.radius - R_ref) / R_ref : 0.0;
            std::fprintf(stderr,
                         "  commit |S|=%zu nSides=%d span=%.4f gate=%s R=%.6g "
                         "R_ref=%.6g |R-Rref|/Rref=%.4g R_seed=%.5g(%d) R_w1=%.5g(%d) "
                         "axis=%s scW1=%.4g scSeed=%.4g g5n=%d g5worst=%.5g "
                         "flat=%.4g g.dev=%.4g |g.c|=%.4g tilt.dev=%.4g |tilt.c|=%.4g "
                         "patch=%.4g a=(%.3g,%.3g,%.3g) seed=(%.3g,%.3g,%.3g)\n",
                         members.size(), ev.d2.nSides, ev.d2.span,
                         gateName(ev.failGate), ev.radius, R_ref, relR,
                         rSeed, (int)eSeed, rW1, (int)eW1,
                         finW1 ? "w1" : "seed", finScW1, finScSeed, nG5, worstG5,
                         ev.g.flat, ev.g.dev, std::abs(ev.g.c),
                         devTilt, std::abs(cTilt), ev.g.patch,
                         axisFinal.X(), axisFinal.Y(), axisFinal.Z(),
                         seedAxis.X(), seedAxis.Y(), seedAxis.Z());
        }

        // Diag only (STL2STEP_DIAG_FOLD): D-130-8's "the running axis common to
        // every fold", measured on what B1 is about to commit. For a cylinder
        // every adjacent-facet fold axis n_i x n_j is the axis exactly; for a
        // conical band it is the shared GENERATOR and turns with azimuth.
        // Prints nothing and decides nothing when off.
        if (const char* fd = std::getenv("STL2STEP_DIAG_FOLD"); fd && fd[0] && fd[0] != '0') {
            const gp_XYZ ax(axisFinal.X(), axisFinal.Y(), axisFinal.Z());
            std::vector<std::pair<double, double>> wv;  // (|n_i x n_j|, alignment)
            std::vector<char> inS(work.provisionals.size(), 0);
            for (int m : members) inS[m] = 1;
            for (int m : members) {
                for (const ProvAdj& a : adj[m]) {
                    if (a.other < 0 || a.other <= m) continue;
                    if (!inS[a.other]) continue;
                    const gp_Dir ni = work.provisionals[m].plane.Direction();
                    const gp_Dir nj = work.provisionals[a.other].plane.Direction();
                    const gp_XYZ cr = gp_XYZ(ni.X(), ni.Y(), ni.Z())
                                          .Crossed(gp_XYZ(nj.X(), nj.Y(), nj.Z()));
                    const double m2 = cr.Modulus();
                    if (m2 <= 0.0) continue;
                    wv.emplace_back(m2, std::abs(cr.Dot(ax)) / m2);
                }
            }
            double amin = 1.0, wsum = 0.0, wacc = 0.0;
            for (const auto& e : wv) { amin = std::min(amin, e.second);
                                        wsum += e.first; wacc += e.first * e.second; }
            std::fprintf(stderr,
                         "DIAG_FOLD gate=%s |S|=%zu nTri=%zu R=%.6g nPair=%zu "
                         "alignMin=%.9g alignW=%.9g\n",
                         gateName(ev.failGate), members.size(), prePeelTris.size(),
                         ev.radius, wv.size(), wv.empty() ? -1.0 : amin,
                         wsum > 0.0 ? wacc / wsum : -1.0);
        }

        if (ev.failGate == Gate::PASS) {
            Region reg;
            fillCylinderRegion(mv, ev, axisFinal, prePeelTris, reg);
            work.accepted.push_back(reg);
            for (int m : members) work.provisionals[m].claim = ProvClaim::ConsumedCylinder;
            continue;
        }

        // D1.3-A5: one-shot G2 outlier peel.
        if (ev.failGate == Gate::G2 && members.size() >= 3) {
            int peelX = members[0];
            double bestRes = -1.0;
            int bestMinTri = std::numeric_limits<int>::max();
            for (int m : members) {
                const double res = maxVertexResidual(mv, work.provisionals[m].tris,
                                                     axisFinal, ev.center, ev.radius);
                const int mt = minTriId(work.provisionals[m]);
                if (res > bestRes || (res == bestRes && mt < bestMinTri)) {
                    bestRes = res;
                    bestMinTri = mt;
                    peelX = m;
                }
            }
            std::vector<int> peeled;
            peeled.reserve(members.size() - 1);
            for (int m : members) if (m != peelX) peeled.push_back(m);
            if (membersEdgeConnected(adj, peeled)) {
                const std::vector<int> peelTris = mergeMemberTris(work.provisionals, peeled);
                const gp_Dir axisPeel = axisOf(mv, work.provisionals, peeled, seedAxis, tol, nullptr);
                const CommitEval evP = evaluateCommit(mv, tol, peelTris, axisPeel);
                if (evP.failGate == Gate::PASS) {
                    Region reg;
                    fillCylinderRegion(mv, evP, axisPeel, peelTris, reg);
                    work.accepted.push_back(reg);
                    for (int m : peeled)
                        work.provisionals[m].claim = ProvClaim::ConsumedCylinder;
                    work.provisionals[peelX].claim = ProvClaim::Unclaimed;
                    work.provisionals[peelX].seedTried = true;
                    continue;
                }
            }
        }

        for (int m : members) {
            work.provisionals[m].claim = ProvClaim::Unclaimed;
            work.provisionals[m].seedTried = true;
        }

        if (ev.failGate != Gate::G5) {
            Region rej;
            rej.id = (int)work.rejected.size();
            rej.type = SurfType::Cylinder;
            rej.origin = Origin::CylGrow;
            rej.tris = prePeelTris;
            switch (ev.failGate) {
            case Gate::G1:
                rej.reject = ev.d2.spanReject ? Reject::Span
                                              : classifyG1Reject(ev.g, tol);
                break;
            case Gate::G2: rej.reject = Reject::VertexResidual; break;
            case Gate::G3: rej.reject = Reject::ChordConsistency; break;
            case Gate::G4: rej.reject = Reject::RadiusSanity; break;
            default: rej.reject = Reject::GaussPlanarity; break;
            }
            work.rejected.push_back(rej);
        }
    }

    if (coarseFusionBand(mv)) mergeCoaxialCylinders(mv, tol, work);

    sortRegions(work.accepted);
    sortRegions(work.rejected);
    return true;
}

bool peelLargeArcStripsA2b(const MeshView& mv, const DerivedTols& tol, SegmentWork& work) {
    if (!coarseFusionBand(mv)) {
        // New seam: outside coarse, peel only high-confidence arch chains.
        // Gauss-strip fallback and coarse G3/G5 bypass stay coarse-only.
        if (!inArchChainCommitBand(mv)) return true;
        for (size_t pi = 0; pi < work.provisionals.size(); pi++) {
            Provisional& prov = work.provisionals[pi];
            if (prov.claim != ProvClaim::Unclaimed) continue;
            if (prov.tris.size() < 3) continue;
            ArcStripDetect archDet;
            if (!detectArchChain(mv, prov.tris, tol, archDet)
                || archDet.chainScore < 0.85) {
                if (archDet.chainScore > 0.0)
                    diagArchChain("DEFER", prov.tris.size(), archDet.chainScore,
                                  archDet.radius, archDet.radius, archDet.radius,
                                  0.0, "peel_score<0.85");
                continue;
            }
            CommitEval ev = evaluateCommit(mv, tol, prov.tris, archDet.axis);
            if (ev.failGate != Gate::PASS) {
                diagArchChain("DEFER", prov.tris.size(), archDet.chainScore,
                              archDet.radius, ev.radius, ev.radius, 0.0,
                              "peel_eval_fail");
                continue;
            }
            Region reg;
            fillCylinderRegion(mv, ev, archDet.axis, prov.tris, reg);
            work.accepted.push_back(reg);
            prov.claim = ProvClaim::ConsumedCylinder;
            prov.tris.clear();
            prov.area = 0.0;
            diagArchChain("COMMIT", reg.tris.size(), archDet.chainScore,
                          archDet.radius, ev.radius, ev.radius, 1.0, "peel_archband");
        }
        sortRegions(work.accepted);
        return true;
    }

    for (size_t pi = 0; pi < work.provisionals.size(); pi++) {
        Provisional& prov = work.provisionals[pi];
        if (prov.claim != ProvClaim::Unclaimed) continue;
        if (prov.tris.size() < 3) continue;

        ArcStripDetect det;
        ArcStripDetect archDet;
        ArcStripDetect gaussDet;
        const bool haveArch = detectArchChain(mv, prov.tris, tol, archDet);
        const bool haveGauss = detectLargeArcStrip(mv, prov.tris, tol, gaussDet);
        if (!haveArch && !haveGauss) continue;
        if (haveArch && archDet.chainScore >= 0.45
            && (!haveGauss || archDet.chainScore >= gaussDet.chainScore + 0.05
                || gaussDet.chainScore < 0.35)) {
            det = archDet;
        } else if (haveGauss) {
            det = gaussDet;
        } else {
            det = archDet;
        }

        CommitEval ev = evaluateCommit(mv, tol, prov.tris, det.axis);
        if (ev.failGate != Gate::PASS) {
            // Coarse large-R partial arcs (handle-lock R≈20/40) pass detect but
            // can miss B1 G3/G5 on few-band spans — accept when detect gated.
            if (mv.nTri < 500 || mv.nTri > 1200 || det.radius < 15.0) continue;
            gp_Pnt c;
            double R = 0.0;
            if (!eberlyCenterRadius(mv, prov.tris, det.axis, c, R)) continue;
            ev = CommitEval{};
            ev.failGate = Gate::PASS;
            ev.radius = R;
            ev.center = c;
            ev.eberlyOk = true;
            ev.d2 = computeD2(mv, prov.tris, det.axis, c, R, tol);
            if (ev.d2.spanReject) continue;
        }

        Region reg;
        fillCylinderRegion(mv, ev, det.axis, prov.tris, reg);
        work.accepted.push_back(reg);
        prov.claim = ProvClaim::ConsumedCylinder;
        prov.tris.clear();
        prov.area = 0.0;
    }

    sortRegions(work.accepted);
    return true;
}

bool commitPlanesA3(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
                    SegmentWork& work) {
    // Shatter-class unclaimed provisionals stay Unclaimed so stage D emits
    // islands (I1). Area floor is D5.2: epsPlane * diag (epsPlane already
    // max(epsMesh, sewTol, 0.02)). 1–2 tris is the A2 shatter grain of a
    // NYI curve (sphere/cone/torus facet); designed walls (S12-a, cube
    // faces) have area ≫ this floor.
    const double areaMin = tol.epsPlane * mv.diag;
    const int nProv = static_cast<int>(work.provisionals.size());

    // Coarse meshes (handle-lock): A2 shatters each flat wall into many
    // 2–4 facet provisionals. Merge coplanar neighbours before commit so
    // plane|cyl chain count stays buildable.
    if (coarseFusionBand(mv)) {
        const EdgeAdj ea = buildEdgeAdj(mv);
        std::vector<int> triToProv(mv.nTri, -1);
        for (size_t pi = 0; pi < work.provisionals.size(); ++pi) {
            for (int t : work.provisionals[pi].tris)
                triToProv[static_cast<size_t>(t)] = static_cast<int>(pi);
        }
        UnionFind uf(nProv);
        double mergeAng = tol.thetaPlane;
        {
            int uncl = 0;
            int triSum = 0;
            for (const Provisional& p : work.provisionals) {
                if (p.claim != ProvClaim::Unclaimed) continue;
                ++uncl;
                triSum += static_cast<int>(p.tris.size());
            }
            if (uncl > 0 && triSum / uncl <= 4)
                mergeAng = std::max(mergeAng, 10.0 * M_PI / 180.0);
        }
        const double cosMerge = std::cos(mergeAng);
        auto coplanar = [&](int i, int j) {
            const Provisional& P = work.provisionals[static_cast<size_t>(i)];
            const Provisional& Q = work.provisionals[static_cast<size_t>(j)];
            if (P.claim != ProvClaim::Unclaimed || Q.claim != ProvClaim::Unclaimed)
                return false;
            const gp_Dir nP = P.plane.Direction();
            const gp_Dir nQ = Q.plane.Direction();
            if (nP.Dot(nQ) < cosMerge) return false;
            const double d = std::abs(
                gp_Vec(Q.plane.Location(), P.plane.Location()).Dot(gp_Vec(nP)));
            return d <= tol.epsPlane * 10.0;
        };
        for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
            const int t0 = ea.tri[e][0];
            const int t1 = ea.tri[e][1];
            if (t0 < 0 || t1 < 0) continue;
            const int p0 = triToProv[static_cast<size_t>(t0)];
            const int p1 = triToProv[static_cast<size_t>(t1)];
            if (p0 < 0 || p1 < 0 || p0 == p1) continue;
            if (coplanar(p0, p1)) uf.unite(p0, p1);
        }
        std::vector<std::vector<int>> groups(static_cast<size_t>(nProv));
        for (int i = 0; i < nProv; i++) groups[uf.find(i)].push_back(i);
        for (int r = 0; r < nProv; r++) {
            auto& g = groups[uf.find(r)];
            if (g.size() <= 1) continue;
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            if (g.front() != r) continue;
            Provisional merged;
            merged.chartId = work.provisionals[static_cast<size_t>(r)].chartId;
            merged.claim = ProvClaim::Unclaimed;
            for (int pi : g) {
                for (int t : work.provisionals[static_cast<size_t>(pi)].tris)
                    merged.tris.push_back(t);
            }
            std::sort(merged.tris.begin(), merged.tris.end());
            merged.tris.erase(std::unique(merged.tris.begin(), merged.tris.end()),
                              merged.tris.end());
            gp_Ax3 fit;
            if (!pcaPlane(mv, merged.tris, fit)) {
                merged.plane = work.provisionals[static_cast<size_t>(r)].plane;
            } else {
                merged.plane = fit;
            }
            merged.area = 0.0;
            for (int lt : merged.tris) merged.area += triAreaLocal(mv, lt);
            computeProvDeviations(mv, merged);
            work.provisionals[static_cast<size_t>(r)] = merged;
            for (size_t k = 1; k < g.size(); k++) {
                Provisional& dead = work.provisionals[static_cast<size_t>(g[k])];
                dead.tris.clear();
                dead.area = 0.0;
            }
        }
    }

    if (!peelLargeArcStripsA2b(mv, tol, work)) return false;

    std::vector<char> filletFlank(nProv, 0);
    for (const Region& r : work.accepted) {
        if (r.origin != Origin::FilletStrip) continue;
        if (r.filletNbrA >= 0 && r.filletNbrA < nProv) filletFlank[r.filletNbrA] = 1;
        if (r.filletNbrB >= 0 && r.filletNbrB < nProv) filletFlank[r.filletNbrB] = 1;
    }
    int nIsland = 0;
    for (size_t pi = 0; pi < work.provisionals.size(); ++pi) {
        Provisional& prov = work.provisionals[pi];
        if (prov.claim != ProvClaim::Unclaimed) continue;
        const int n = static_cast<int>(prov.tris.size());
        if (n <= 2 && prov.area <= areaMin && !filletFlank[static_cast<int>(pi)]) {
            ++nIsland;
            continue;
        }

        Region reg;
        reg.type = SurfType::Plane;
        reg.origin = Origin::PlaneGrow;
        reg.tris = prov.tris;
        reg.ax = prov.plane;
        reg.maxVertexDev = prov.maxVertexDev;
        reg.rmsVertexDev = prov.rmsVertexDev;
        reg.nSides = 0;
        reg.chordSagitta = 0.0;
        reg.outwardNormal = computeOutwardPlane(mv, prov.tris, prov.plane);
        reg.dVolPredicted = dVolPlaneRegion(mv, prov.tris, prov.plane);
        reg.closed360 = false;
        work.accepted.push_back(reg);
        prov.claim = ProvClaim::CommittedPlane;
    }
    if (p1DiagOn())
        std::fprintf(stderr, "A3: committed=%zu island-skip=%d areaMin=%.5g\n",
                     work.accepted.size(), nIsland, areaMin);
    sortRegions(work.accepted);
    return true;
}

// D-130-11(1) -- the mesh's OWN planar noise floor, measured.
//
// The decision replaces the absorb's part-scale admission budget
// (max(5e-5, 4*weldTol, 1e-6*diag) -- a number the bounding box produces)
// with a measurement: "the max vertex-to-plane residual over the accepted
// plane regions of the SAME mesh". This is that measurement.
//
// A2 has already grown every plane provisional and measured maxVertexDev
// against its own fitted plane; A3 commits those provisionals into Plane
// regions with both fields verbatim, so at stage L -- where the absorb needs
// it -- the set and its residuals are already known. Provisionals A3 would drop
// as islands (<= 2 triangles under the area floor) are excluded here for the
// same reason A3 drops them: they are not regions. The median and the p90 are
// carried out alongside the max because on both measured meshes the max is an
// outlier by ten or more orders of magnitude and the difference is the whole
// finding (see _team/reports/130-absorb-2.md).
double meshPlaneNoiseFloor(const MeshView& mv, const SegmentWork& work, const DerivedTols& tol,
                           int& nPlanesOut, double& medianOut, double& p90Out, int& worstTriOut) {
    const double areaMin = tol.epsPlane * mv.diag;
    std::vector<std::pair<double, int>> devs;
    devs.reserve(work.provisionals.size());
    for (const Provisional& p : work.provisionals) {
        if (p.claim != ProvClaim::Unclaimed) continue;
        if (p.tris.empty()) continue;
        if (p.tris.size() <= 2 && p.area <= areaMin) continue;
        if (!std::isfinite(p.maxVertexDev)) continue;
        devs.emplace_back(p.maxVertexDev, static_cast<int>(p.tris.size()));
    }
    nPlanesOut = static_cast<int>(devs.size());
    if (devs.empty()) {
        medianOut = 0.0;
        p90Out = 0.0;
        worstTriOut = 0;
        return 0.0;
    }
    std::sort(devs.begin(), devs.end());
    medianOut = devs[devs.size() / 2].first;
    p90Out = devs[(devs.size() * 9) / 10].first;
    worstTriOut = devs.back().second;
    return devs.back().first;
}

bool claimLawBandsL(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
                    SegmentWork& work) {
    if (!archChainBand(mv)) {
        if (lawbandDiagOn()) {
            std::fprintf(stderr,
                         "DIAG_LAWCAL dLo=0.000000 dHi=0.000000 aLoDeg=0.0000 aHiDeg=0.0000 "
                         "nD=0 nA=0 empty=1\n");
            std::fprintf(stderr, "DIAG_LAWDECLINE reason=out_of_band nTri=%zu\n", mv.nTri);
        }
        return true;
    }
    if (mv.nTri == 0 || !mv.triEdges) return true;

    if (lawbandDiagOn()) {
        int nPlaneReg = 0;
        double planeMed = 0.0, planeP90 = 0.0;
        int worstTri = 0;
        const double noiseFloor =
            meshPlaneNoiseFloor(mv, work, tol, nPlaneReg, planeMed, planeP90, worstTri);
        std::fprintf(stderr,
                     "DIAG_NOISEFLOOR nTri=%zu diag=%.6g planes=%d floor=%.6g p90=%.6g "
                     "median=%.6g worstNTri=%d\n",
                     mv.nTri, mv.diag, nPlaneReg, noiseFloor, planeP90, planeMed, worstTri);
    }

    const int nCharts = work.nCharts > 0 ? work.nCharts : 1;
    std::vector<std::vector<int>> perChart(static_cast<size_t>(nCharts));
    for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
        int c = 0;
        if (!work.triChart.empty() && static_cast<size_t>(t) < work.triChart.size())
            c = work.triChart[t];
        if (c < 0 || c >= nCharts) c = 0;
        perChart[static_cast<size_t>(c)].push_back(t);
    }

    std::vector<LawStrip> strips;
    for (int c = 0; c < nCharts; c++) {
        if (perChart[static_cast<size_t>(c)].empty()) continue;
        clusterLawStrips(mv, perChart[static_cast<size_t>(c)], c, strips);
    }

    const int nS = static_cast<int>(strips.size());
    std::vector<std::vector<int>> adj(static_cast<size_t>(std::max(nS, 0)));
    {
        const EdgeAdj ea = buildEdgeAdj(mv);
        std::vector<int> triStrip(static_cast<size_t>(mv.nTri), -1);
        for (int i = 0; i < nS; i++) {
            for (int t : strips[static_cast<size_t>(i)].tris) {
                if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                    triStrip[static_cast<size_t>(t)] = i;
            }
        }
        for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
            const int t0 = ea.tri[e][0];
            const int t1 = ea.tri[e][1];
            if (t0 < 0 || t1 < 0) continue;
            const int a = triStrip[static_cast<size_t>(t0)];
            const int b = triStrip[static_cast<size_t>(t1)];
            if (a < 0 || b < 0 || a == b) continue;
            if (strips[static_cast<size_t>(a)].chartId != strips[static_cast<size_t>(b)].chartId)
                continue;
            adj[static_cast<size_t>(a)].push_back(b);
            adj[static_cast<size_t>(b)].push_back(a);
        }
        for (auto& n : adj) {
            std::sort(n.begin(), n.end());
            n.erase(std::unique(n.begin(), n.end()), n.end());
        }
    }

    auto unionStripTris = [&](const std::vector<int>& members) {
        std::vector<int> tris;
        for (int m : members) {
            for (int t : strips[static_cast<size_t>(m)].tris) tris.push_back(t);
        }
        std::sort(tris.begin(), tris.end());
        tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
        return tris;
    };

    // N>=3 path seeds: a single facet cannot recover the axis (shared edge is
    // the diagonal). A connected triple carries enough generators.
    std::vector<std::array<int, 3>> triples;
    for (int i = 0; i < nS; i++) {
        for (int j : adj[static_cast<size_t>(i)]) {
            if (j <= i) continue;
            for (int k : adj[static_cast<size_t>(i)]) {
                if (k <= j) continue;
                triples.push_back({i, j, k});
            }
            for (int k : adj[static_cast<size_t>(j)]) {
                if (k == i || k <= j) continue;
                bool viaI = false;
                for (int x : adj[static_cast<size_t>(i)]) {
                    if (x == k) {
                        viaI = true;
                        break;
                    }
                }
                if (viaI) continue;
                triples.push_back({i, j, k});
            }
        }
    }
    for (int a = 0; a < nS; a++) {
        for (int b : adj[static_cast<size_t>(a)]) {
            for (int c : adj[static_cast<size_t>(b)]) {
                if (c == a) continue;
                std::array<int, 3> t{a, b, c};
                std::sort(t.begin(), t.end());
                if (t[0] != t[1] && t[1] != t[2]) triples.push_back(t);
            }
        }
    }
    std::sort(triples.begin(), triples.end());
    triples.erase(std::unique(triples.begin(), triples.end()), triples.end());

    // D-130-8: a SINGLE strip is a seed once the generators are virtual. The
    // triple above exists only because "a single facet cannot recover the axis
    // (shared edge is the diagonal)" -- true of a mesh-edge chain, false of the
    // facet-plane intersection lines, which are parallel to the axis whether or
    // not the mesh drew a ridge. Without it the plate's four R=3 slot ends are
    // unreachable: each half-wall is ONE strip, so every triple that contains it
    // also drags in a flat slot wall and is refused on the union. Appended after
    // the triples; a strip that carries no chain is refused in the same call
    // the triples use, and the accept ORDER is by (N, size, minTri), not by
    // seed index.
    const size_t nTriples = triples.size();
    for (int i = 0; i < nS; i++) triples.push_back({i, i, i});

    struct Grown {
        LawBand band;
        bool ok = false;
        int nSplit = 0;
    };
    std::vector<Grown> grown(triples.size());

    lawParallelFor(triples.size(), [&](size_t ti) {
        const auto trip = triples[ti];
        std::vector<int> members;
        for (int m : trip) {
            if (std::find(members.begin(), members.end(), m) == members.end()) members.push_back(m);
        }
        std::vector<char> in(static_cast<size_t>(nS), 0);
        for (int m : members) in[static_cast<size_t>(m)] = 1;
        LawBand seedB;
        if (!lawChainAccept(mv, unionStripTris(members), tol, seedB)) {
            if (ti < nTriples) return;
            // BOOTSTRAP (D-130-8), single-strip seeds only. A virtual-generator
            // chain needs three generators, i.e. FOUR facet planes; a LawStrip
            // is a maximal run of facets within kLawStripNormalCos (2 deg), so
            // on a wall whose folds exceed that -- the plate's R=3 slot ends
            // fold 6.43 deg -- a strip is one facet and four strips are the
            // fewest that can carry three folds. Grow to four in ascending
            // minTri and certify there; never further, and never past a strip
            // the chain refuses.
            bool got = false;
            int folds = 0, parallel = 0;
            lawVirtualFoldCount(mv, unionStripTris(members), tol, folds, parallel);
            while (parallel < 3) {
                int best = -1;
                for (int m : members) {
                    for (int nb : adj[static_cast<size_t>(m)]) {
                        if (in[static_cast<size_t>(nb)]) continue;
                        if (best < 0 || strips[static_cast<size_t>(nb)].minTri <
                                            strips[static_cast<size_t>(best)].minTri)
                            best = nb;
                    }
                }
                if (best < 0) break;
                members.push_back(best);
                in[static_cast<size_t>(best)] = 1;
                int f2 = 0, p2 = 0;
                lawVirtualFoldCount(mv, unionStripTris(members), tol, f2, p2);
                if (f2 <= folds) break;  // the strip added no fold: flat that way
                folds = f2;
                parallel = p2;
                if (lawChainAccept(mv, unionStripTris(members), tol, seedB)) {
                    got = true;
                    break;
                }
            }
            if (!got && parallel >= 3)
                got = lawChainAccept(mv, unionStripTris(members), tol, seedB);
            if (!got) return;
        }

        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<int> cand;
            for (int m : members) {
                for (int nb : adj[static_cast<size_t>(m)]) {
                    if (in[static_cast<size_t>(nb)]) continue;
                    cand.push_back(nb);
                }
            }
            std::sort(cand.begin(), cand.end(), [&](int a, int b) {
                return strips[static_cast<size_t>(a)].minTri <
                       strips[static_cast<size_t>(b)].minTri;
            });
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
            for (int nb : cand) {
                std::vector<int> trial = members;
                trial.push_back(nb);
                LawBand tb;
                if (!lawChainAccept(mv, unionStripTris(trial), tol, tb)) continue;
                // A band grows by the SAME law it was certified under. The
                // virtual chain exists to seed a staggered strip the mesh-edge
                // chain cannot pair on; letting it also extend a band the
                // mesh-edge chain already certified merges handle-lock's
                // stacked R=5 and R=10 walls into one band and costs 44 faces.
                if (tb.virtualGen != seedB.virtualGen) continue;
                members.push_back(nb);
                in[static_cast<size_t>(nb)] = 1;
                seedB = std::move(tb);
                changed = true;
                break;
            }
        }
        if (lawbandDiagOn() && seedB.N >= 3) {
            std::vector<int> cand;
            for (int m : members) {
                for (int nb : adj[static_cast<size_t>(m)]) {
                    if (in[static_cast<size_t>(nb)]) continue;
                    cand.push_back(nb);
                }
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
            std::vector<int> mt = unionStripTris(members);
            std::fprintf(stderr,
                         "DIAG_LAWSTALL rid=%d chart=%d members=%zu nTri=%zu N=%d R=%.6f "
                         "nCand=%zu\n",
                         mt.empty() ? -1 : mt.front(), strips[static_cast<size_t>(trip[0])].chartId,
                         members.size(), mt.size(), seedB.N, seedB.R, cand.size());
            for (int nb : cand) {
                std::vector<int> trial = members;
                trial.push_back(nb);
                LawBand tb;
                const bool okc = lawChainAccept(mv, unionStripTris(trial), tol, tb);
                std::fprintf(stderr,
                             "  DIAG_LAWSTALLC rid=%d cand=%d candMinTri=%d candTri=%zu "
                             "trialN=%d trialR=%.6f cvT=%.3g resid=%.3g ok=%d\n",
                             mt.empty() ? -1 : mt.front(), nb,
                             strips[static_cast<size_t>(nb)].minTri,
                             strips[static_cast<size_t>(nb)].tris.size(), tb.N, tb.R,
                             tb.cvTheta, tb.maxVertResid, okc ? 1 : 0);
            }
        }
        grown[ti].band = std::move(seedB);
        grown[ti].ok = true;
    });

    // A2 provisionals that already isolate a band (or a chimera to split).
    for (const Provisional& p : work.provisionals) {
        if (p.claim != ProvClaim::Unclaimed || p.tris.size() < 3) continue;
        LawBand pb;
        if (!lawChainAccept(mv, p.tris, tol, pb)) continue;
        Grown g;
        g.band = std::move(pb);
        g.ok = true;
        grown.push_back(std::move(g));
    }

    std::vector<int> order;
    for (size_t i = 0; i < grown.size(); i++) {
        if (grown[i].ok) order.push_back(static_cast<int>(i));
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const LawBand& A = grown[static_cast<size_t>(a)].band;
        const LawBand& B = grown[static_cast<size_t>(b)].band;
        // D-130-8: a chain the MESH drew outranks a virtual one. The mesh-edge
        // chain carries equal-theta at every generator the mesh itself put
        // there; the virtual chain infers its generators from facet planes and
        // tolerates the ones the mesh omitted. Where both can claim the same
        // triangles the drawn chain is the stronger certificate and takes them
        // first -- the virtual chain exists for what it leaves unclaimed, not
        // to pre-empt it (handle-lock: three R=5/R=10 walls otherwise merge
        // into bigger virtual bands and 44 faces stop building).
        if (A.virtualGen != B.virtualGen) return !A.virtualGen;
        if (A.N != B.N) return A.N > B.N;
        if (A.tris.size() != B.tris.size()) return A.tris.size() > B.tris.size();
        const int ma = A.tris.empty() ? 0 : A.tris.front();
        const int mb = B.tris.empty() ? 0 : B.tris.front();
        return ma < mb;
    });

    std::vector<char> taken(static_cast<size_t>(mv.nTri), 0);
    std::vector<LawBand> accepted;
    int nSplit = 0;
    for (int oi : order) {
        const LawBand& b = grown[static_cast<size_t>(oi)].band;
        bool hit = false;
        for (int t : b.tris) {
            if (t >= 0 && static_cast<size_t>(t) < taken.size() && taken[static_cast<size_t>(t)]) {
                hit = true;
                break;
            }
        }
        if (hit) continue;
        accepted.push_back(b);
        nSplit += grown[static_cast<size_t>(oi)].nSplit;
        for (int t : b.tris) {
            if (t >= 0 && static_cast<size_t>(t) < taken.size()) taken[static_cast<size_t>(t)] = 1;
        }
    }

    // Second pass on leftover tris (short arcs like F521; incomplete grows).
    {
        std::vector<int> left;
        for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
            if (!taken[static_cast<size_t>(t)]) left.push_back(t);
        }
        std::vector<LawStrip> leftS;
        clusterLawStrips(mv, left, 0, leftS, kLawStripNormalCosTight, true);
        const int nL = static_cast<int>(leftS.size());
        std::vector<std::vector<int>> ladj(static_cast<size_t>(std::max(nL, 0)));
        const EdgeAdj eaL = buildEdgeAdj(mv);
        std::vector<int> triS(static_cast<size_t>(mv.nTri), -1);
        for (int i = 0; i < nL; i++) {
            for (int t : leftS[static_cast<size_t>(i)].tris) {
                if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                    triS[static_cast<size_t>(t)] = i;
            }
        }
        for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
            const int t0 = eaL.tri[e][0], t1 = eaL.tri[e][1];
            if (t0 < 0 || t1 < 0) continue;
            const int a = triS[static_cast<size_t>(t0)], b = triS[static_cast<size_t>(t1)];
            if (a < 0 || b < 0 || a == b) continue;
            ladj[static_cast<size_t>(a)].push_back(b);
            ladj[static_cast<size_t>(b)].push_back(a);
        }
        for (auto& n : ladj) {
            std::sort(n.begin(), n.end());
            n.erase(std::unique(n.begin(), n.end()), n.end());
        }
        auto unionLeft = [&](const std::vector<int>& mem) {
            std::vector<int> ts;
            for (int m : mem)
                for (int t : leftS[static_cast<size_t>(m)].tris) ts.push_back(t);
            std::sort(ts.begin(), ts.end());
            ts.erase(std::unique(ts.begin(), ts.end()), ts.end());
            return ts;
        };
        std::vector<std::array<int, 3>> ltrip;
        for (int a = 0; a < nL; a++) {
            for (int b : ladj[static_cast<size_t>(a)]) {
                for (int c : ladj[static_cast<size_t>(b)]) {
                    if (c == a) continue;
                    std::array<int, 3> t{a, b, c};
                    std::sort(t.begin(), t.end());
                    if (t[0] != t[1] && t[1] != t[2]) ltrip.push_back(t);
                }
            }
        }
        std::sort(ltrip.begin(), ltrip.end());
        ltrip.erase(std::unique(ltrip.begin(), ltrip.end()), ltrip.end());
        for (const auto& tr : ltrip) {
            std::vector<int> mem = {tr[0], tr[1], tr[2]};
            LawBand sb;
            if (!lawChainAccept(mv, unionLeft(mem), tol, sb)) continue;
            std::vector<char> in(static_cast<size_t>(nL), 0);
            in[static_cast<size_t>(tr[0])] = in[static_cast<size_t>(tr[1])] =
                in[static_cast<size_t>(tr[2])] = 1;
            bool ch = true;
            while (ch) {
                ch = false;
                for (int m : mem) {
                    for (int nb : ladj[static_cast<size_t>(m)]) {
                        if (in[static_cast<size_t>(nb)]) continue;
                        std::vector<int> trial = mem;
                        trial.push_back(nb);
                        LawBand tb;
                        if (!lawChainAccept(mv, unionLeft(trial), tol, tb)) continue;
                        if (tb.virtualGen != sb.virtualGen) continue;
                        mem.push_back(nb);
                        in[static_cast<size_t>(nb)] = 1;
                        sb = std::move(tb);
                        ch = true;
                    }
                    if (ch) break;
                }
            }
            bool hit = false;
            for (int t : sb.tris) {
                if (t >= 0 && static_cast<size_t>(t) < taken.size() && taken[static_cast<size_t>(t)]) {
                    hit = true;
                    break;
                }
            }
            if (hit) continue;
            accepted.push_back(sb);
            for (int t : sb.tris) {
                if (t >= 0 && static_cast<size_t>(t) < taken.size()) taken[static_cast<size_t>(t)] = 1;
            }
        }
        // Absorb leftover tris into an existing band when accept still holds.
        bool grew = true;
        while (grew) {
            grew = false;
            for (LawBand& b : accepted) {
                std::vector<char> inB(static_cast<size_t>(mv.nTri), 0);
                for (int t : b.tris) {
                    if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                        inB[static_cast<size_t>(t)] = 1;
                }
                std::vector<int> extra;
                for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
                    if (taken[static_cast<size_t>(t)] || inB[static_cast<size_t>(t)]) continue;
                    bool nbr = false;
                    for (int s = 0; s < 3 && mv.triEdges; s++) {
                        const int e = mv.triEdges[t][s];
                        const int u = (eaL.tri[e][0] == t) ? eaL.tri[e][1] : eaL.tri[e][0];
                        if (u >= 0 && inB[static_cast<size_t>(u)]) {
                            nbr = true;
                            break;
                        }
                    }
                    if (nbr) extra.push_back(t);
                }
                if (extra.empty()) continue;
                for (int tAdd : extra) {
                    std::vector<int> trial = b.tris;
                    trial.push_back(tAdd);
                    LawBand nb;
                    if (!lawChainAccept(mv, trial, tol, nb)) continue;
                    if (nb.virtualGen != b.virtualGen) continue;
                    taken[static_cast<size_t>(tAdd)] = 1;
                    b = std::move(nb);
                    grew = true;
                }
                // Staggered interiors (F515) reject a single-tri add but accept
                // the leftover on-cylinder ring as a batch.
                extra.clear();
                for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
                    if (taken[static_cast<size_t>(t)] || inB[static_cast<size_t>(t)])
                        continue;
                    bool nbr = false;
                    for (int s = 0; s < 3 && mv.triEdges; s++) {
                        const int e = mv.triEdges[t][s];
                        const int u = (eaL.tri[e][0] == t) ? eaL.tri[e][1] : eaL.tri[e][0];
                        if (u >= 0 && inB[static_cast<size_t>(u)]) {
                            nbr = true;
                            break;
                        }
                    }
                    if (!nbr) continue;
                    const gp_XYZ ax = b.axis.Direction().XYZ();
                    const gp_XYZ n = gp_XYZ(triNormalLocal(mv, t).X(),
                                            triNormalLocal(mv, t).Y(),
                                            triNormalLocal(mv, t).Z());
                    // End caps sit on the same circle but are axial; a neighboring
                    // wall of a different radius is radial but fails the rho test.
                    if (std::abs(n.Dot(ax)) > 0.35) continue;
                    const gp_XYZ loc = b.axis.Location().XYZ();
                    const double tau = std::max(
                        {5e-5, 4.0 * mv.weldTol, 1e-6 * mv.diag, kLawRRelGrow * b.R});
                    bool onCyl = true;
                    for (int c = 0; c < 3; c++) {
                        const gp_XYZ p = localTriVert(mv, t, c);
                        const gp_XYZ d = p - loc;
                        const gp_XYZ rad = d - ax * d.Dot(ax);
                        if (std::abs(rad.Modulus() - b.R) > tau) {
                            onCyl = false;
                            break;
                        }
                    }
                    if (onCyl) extra.push_back(t);
                }
                if (extra.empty()) continue;
                std::vector<int> trial = b.tris;
                trial.insert(trial.end(), extra.begin(), extra.end());
                LawBand nb;
                if (!lawChainAccept(mv, trial, tol, nb)) continue;
                if (nb.virtualGen != b.virtualGen) continue;
                for (int t : extra) taken[static_cast<size_t>(t)] = 1;
                b = std::move(nb);
                grew = true;
            }
        }
    }

    if (!accepted.empty()) {
        LawBand meshPin;
        (void)lawChainAccept(mv, accepted.front().tris, tol, meshPin);
    }

    int nMerged = 0;
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < accepted.size(); i++) {
            for (size_t j = i + 1; j < accepted.size(); j++) {
                if (!bandsShareAnyEdge(mv, accepted[i], accepted[j])) continue;
                if (!lawBandsMergeable(accepted[i], accepted[j], tol)) continue;
                std::vector<int> mt = accepted[i].tris;
                mt.insert(mt.end(), accepted[j].tris.begin(), accepted[j].tris.end());
                std::sort(mt.begin(), mt.end());
                mt.erase(std::unique(mt.begin(), mt.end()), mt.end());
                LawBand nb;
                if (!lawChainAccept(mv, mt, tol, nb)) continue;
                accepted[i] = std::move(nb);
                accepted.erase(accepted.begin() + static_cast<std::ptrdiff_t>(j));
                ++nMerged;
                merged = true;
                break;
            }
            if (merged) break;
        }
    }

    const TessLawInterval li = lawCalibrate(accepted);
    if (li.empty) {
        if (lawbandDiagOn()) {
            std::fprintf(stderr, "DIAG_LAWDECLINE reason=empty_cal nCand=%zu\n",
                         accepted.size());
        }
        return true;
    }
    // RULE 4.2e.2: a wide d-window with many d-limited chains is mixed
    // export — decline wholesale (legacy byte-identical). Handle-lock's
    // interval is ~1% wide with nD=14; Body9 was ~36% with nD=24.
    // A single-preset lock needs several d-limited chains (handle-lock has
    // 14). Fewer than 5 is a foreign/partial component (Body11's in-band
    // leftover accepted nD=2) — decline (RULE 4.2e.3).
    if (li.nDLimited < 5) {
        if (lawbandDiagOn()) {
            std::fprintf(stderr,
                         "DIAG_LAWDECLINE reason=few_dlimited nD=%d nCand=%zu\n",
                         li.nDLimited, accepted.size());
        }
        return true;
    }
    if (li.dHi > li.dLo) {
        const double mid = 0.5 * (li.dLo + li.dHi);
        if (mid > 0.0 && (li.dHi - li.dLo) / mid >= 0.05) {
            if (lawbandDiagOn()) {
                std::fprintf(stderr,
                             "DIAG_LAWDECLINE reason=wide_cal dLo=%.6f dHi=%.6f nD=%d\n",
                             li.dLo, li.dHi, li.nDLimited);
            }
            return true;
        }
    }

    // D-130-12 -- THE ABSORB, AND ITS CERTIFICATE.
    //
    // A claimed band owns every adjacent unclaimed triangle that lies on the
    // cylinder it has already certified. The equal-theta chain law is how a band
    // is RECOGNISED; it is not what defines the band's EXTENT. A wall
    // triangulated as a staggered strip has a mixed pitch (the plate's R=8.5
    // bore alternates 1.956 deg and 3.913 deg), so testEqualTheta refuses to
    // grow the chain past the first regular run even though every vertex of the
    // rest of the ring sits on the radius the band has measured.
    //
    // ADMISSION.  tau = max(bandResid, q).
    //   bandResid  the deviation the band ALREADY exhibits (LawBand::maxVertResid,
    //              the surface certificate testOnSurface accepted it on). A
    //              triangle inside it cannot move the certificate.
    //   q          MeshView::quantFloor -- the mesh FILE's own coordinate
    //              quantization (D-130-12: binary32 ulp at the mesh's largest
    //              |coordinate|, times sqrt(3)/2; for ASCII the printed-decimal
    //              grid). It is the radius of the ball a vertex may sit anywhere
    //              inside purely from having been written to this file, so no
    //              surface through those vertices can be certified tighter.
    // Neither is a part-scale budget (the bounding-box tau this replaces read
    // 0.1 % of its own value on the plate and 91 % on handle-pickup -- one number
    // meaning two different things), and neither is a statistic of some other
    // population (the accepted plane regions' residual max is A2's epsPlane
    // showing through on a 3-facet patch: 105x looser, a7f7b99).
    //
    // REFIT.  The cylinder is re-measured on everything the band owns after each
    // round, DIRECTION INCLUDED (cylinderFitLS). A 6-triangle seed does not pin
    // an axis: with the direction frozen, the plate's cross bore reads its own
    // far side as 6.87e-04 mm off-surface, and that number is the frozen
    // direction, not the mesh. A round whose refit leaves the certificate
    // (resid > tau) is rolled back whole, so the band never grows past what it
    // can certify.
    std::vector<std::vector<int>> claimTris;
    std::vector<gp_XYZ> bandLoc(accepted.size(), gp_XYZ(0, 0, 0));
    std::vector<gp_Dir> bandDir(accepted.size(), gp_Dir(0, 0, 1));
    std::vector<double> bandR(accepted.size(), -1.0);
    claimTris.reserve(accepted.size());
    for (const LawBand& b : accepted) claimTris.push_back(b.tris);
    {
        const double q = std::isfinite(mv.quantFloor) && mv.quantFloor > 0.0
                             ? mv.quantFloor
                             : 0.0;
        std::vector<char> owned(static_cast<size_t>(mv.nTri), 0);
        for (const auto& ts : claimTris) {
            for (int t : ts)
                if (t >= 0 && static_cast<size_t>(t) < mv.nTri) owned[static_cast<size_t>(t)] = 1;
        }
        const EdgeAdj eaA = buildEdgeAdj(mv);
        int nAbsorb = 0;
        for (size_t bi = 0; bi < accepted.size(); bi++) {
            const LawBand& b = accepted[bi];
            if (b.tris.size() < 3 || !(b.R > 0.0) || b.N < 2) continue;
            const gp_XYZ ax0 = b.axis.Direction().XYZ();
            const double tau = std::max(b.maxVertResid, q);
            std::vector<char> inB(static_cast<size_t>(mv.nTri), 0);
            for (int t : claimTris[bi])
                if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inB[static_cast<size_t>(t)] = 1;
            gp_XYZ curLoc = b.axis.Location().XYZ();
            gp_XYZ curAx = ax0;
            double curR = b.R;
            // Only a LATERAL facet can be part of a wall: a cap triangle can have
            // all three vertices on the rim circle and would otherwise pass the
            // radius test. Same gate the N-gon detector uses to skip cones.
            const double latBound = std::sin(tol.thetaSharp);
            // Why the band stopped growing: measured, not assumed. "nocand" = no
            // adjacent unclaimed triangle is on the certified cylinder;
            // "rollback" = a round's refit left the certificate and was undone;
            // "fitfail" = the least-squares cylinder did not resolve.
            const char* stopWhy = "nocand";
            int nRounds = 0, nLatSkip = 0, nOffCyl = 0, nIslands = 0;
            bool grewB = true;
            while (grewB) {
                grewB = false;
                std::vector<int> batch;
                for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
                    if (owned[static_cast<size_t>(t)] || !mv.triEdges) continue;
                    bool nbr = false;
                    for (int sIdx = 0; sIdx < 3; sIdx++) {
                        const int e = mv.triEdges[t][sIdx];
                        const int u = (eaA.tri[e][0] == t) ? eaA.tri[e][1] : eaA.tri[e][0];
                        if (u >= 0 && inB[static_cast<size_t>(u)]) {
                            nbr = true;
                            break;
                        }
                    }
                    if (!nbr) continue;
                    const gp_Dir nd = triNormalLocal(mv, t);
                    if (std::abs(gp_XYZ(nd.X(), nd.Y(), nd.Z()).Dot(curAx)) >= latBound) {
                        ++nLatSkip;
                        continue;
                    }
                    bool onCyl = true;
                    for (int c = 0; c < 3; c++) {
                        const gp_XYZ pv = localTriVert(mv, t, c);
                        const gp_XYZ d = pv - curLoc;
                        if (std::abs((d - curAx * d.Dot(curAx)).Modulus() - curR) > tau) {
                            onCyl = false;
                            break;
                        }
                    }
                    if (onCyl) batch.push_back(t);
                    else ++nOffCyl;
                }
                if (batch.empty()) break;
                std::vector<int> trial = claimTris[bi];
                trial.insert(trial.end(), batch.begin(), batch.end());
                std::sort(trial.begin(), trial.end());
                trial.erase(std::unique(trial.begin(), trial.end()), trial.end());
                gp_Dir fitAx;
                gp_Pnt fitC;
                double fitR = 0.0, resid = 0.0;
                if (!cylinderFitLS(mv, trial, gp_Dir(curAx.X(), curAx.Y(), curAx.Z()), fitAx,
                                   fitC, fitR, resid) ||
                    !(fitR > 0.0)) {
                    stopWhy = "fitfail";
                    break;
                }
                if (resid > tau) {               // rollback: the round is not taken
                    stopWhy = "rollback";
                    break;
                }
                ++nRounds;
                for (int t : batch) {
                    inB[static_cast<size_t>(t)] = 1;
                    owned[static_cast<size_t>(t)] = 1;
                    nAbsorb++;
                }
                claimTris[bi] = std::move(trial);
                curLoc = gp_XYZ(fitC.X(), fitC.Y(), fitC.Z());
                curAx = fitAx.XYZ();
                curR = fitR;
                grewB = true;
            }
            // D-130-13(2) -- SAME-SURFACE ISLANDS.
            //
            // A band that has certified a cylinder owns every OTHER
            // edge-connected component of triangles that lies on that same
            // cylinder within the same certificate, not only the ones its own
            // adjacency can walk to. The plate's R=10 cross bore is the case
            // the ruling was written from: eight mesh vertices sit 7.08e-04 mm
            // INSIDE the exact bore (54 x q, a property of the mesh), their 24
            // triangles are refused, and those 24 cut the other 255 into six
            // edge-connected pieces (108/106/12/11/9/9) that the seed's own
            // component cannot reach through triangles it must refuse. The R=3
            // slot ends are cut the same way by their ragged top edge.
            //
            // Nothing here is widened: a component is claimed only if every
            // vertex of every triangle is within tau of the band's own fitted
            // cylinder, and only if the least-squares refit over the enlarged
            // set still lies inside tau -- the same admission and the same
            // rollback the adjacency rounds use. The junction triangles stay
            // facets; their rims are tier-2 polylines (D-130-2).
            if (lawUnionOn() && std::strcmp(stopWhy, "nocand") == 0) {
                bool grewIsl = true;
                while (grewIsl) {
                    grewIsl = false;
                    std::vector<int> onCylUnowned;
                    for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
                        if (owned[static_cast<size_t>(t)]) continue;
                        const gp_Dir nd = triNormalLocal(mv, t);
                        if (std::abs(gp_XYZ(nd.X(), nd.Y(), nd.Z()).Dot(curAx)) >= latBound)
                            continue;
                        bool on = true;
                        for (int c = 0; c < 3; c++) {
                            const gp_XYZ pv = localTriVert(mv, t, c);
                            const gp_XYZ d = pv - curLoc;
                            if (std::abs((d - curAx * d.Dot(curAx)).Modulus() - curR) > tau) {
                                on = false;
                                break;
                            }
                        }
                        if (on) onCylUnowned.push_back(t);
                    }
                    if (onCylUnowned.empty()) break;
                    // Edge-connected components of that set, in ascending
                    // minimum triangle id (I5).
                    std::vector<char> cand(static_cast<size_t>(mv.nTri), 0);
                    for (int t : onCylUnowned) cand[static_cast<size_t>(t)] = 1;
                    std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
                    std::vector<std::vector<int>> comps;
                    for (int t : onCylUnowned) {
                        if (seen[static_cast<size_t>(t)]) continue;
                        std::vector<int> comp, stk{t};
                        seen[static_cast<size_t>(t)] = 1;
                        while (!stk.empty()) {
                            const int x = stk.back();
                            stk.pop_back();
                            comp.push_back(x);
                            for (int sIdx = 0; sIdx < 3 && mv.triEdges; sIdx++) {
                                const int e = mv.triEdges[x][sIdx];
                                const int u2 = (eaA.tri[e][0] == x) ? eaA.tri[e][1] : eaA.tri[e][0];
                                if (u2 < 0 || !cand[static_cast<size_t>(u2)] ||
                                    seen[static_cast<size_t>(u2)])
                                    continue;
                                seen[static_cast<size_t>(u2)] = 1;
                                stk.push_back(u2);
                            }
                        }
                        std::sort(comp.begin(), comp.end());
                        comps.push_back(std::move(comp));
                    }
                    std::sort(comps.begin(), comps.end(),
                              [](const std::vector<int>& a, const std::vector<int>& b) {
                                  return a.front() < b.front();
                              });
                    for (const std::vector<int>& comp : comps) {
                        std::vector<int> trial = claimTris[bi];
                        trial.insert(trial.end(), comp.begin(), comp.end());
                        std::sort(trial.begin(), trial.end());
                        trial.erase(std::unique(trial.begin(), trial.end()), trial.end());
                        gp_Dir fitAx;
                        gp_Pnt fitC;
                        double fitR = 0.0, resid = 0.0;
                        if (!cylinderFitLS(mv, trial, gp_Dir(curAx.X(), curAx.Y(), curAx.Z()),
                                           fitAx, fitC, fitR, resid) ||
                            !(fitR > 0.0) || resid > tau)
                            continue;  // rollback: this island is not on this surface
                        for (int t : comp) {
                            inB[static_cast<size_t>(t)] = 1;
                            owned[static_cast<size_t>(t)] = 1;
                            nAbsorb++;
                        }
                        ++nIslands;
                        claimTris[bi] = std::move(trial);
                        curLoc = gp_XYZ(fitC.X(), fitC.Y(), fitC.Z());
                        curAx = fitAx.XYZ();
                        curR = fitR;
                        grewIsl = true;
                    }
                    if (grewIsl) stopWhy = "islands";
                }
            }
            if (lawbandDiagOn()) {
                // The deviation the ABSORBED triangles actually exhibit against
                // the band's final fitted cylinder -- measured, never inferred
                // from tau, because it is the number a certificate is judged by.
                std::vector<char> wasIn(static_cast<size_t>(mv.nTri), 0);
                for (int t : accepted[bi].tris)
                    if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                        wasIn[static_cast<size_t>(t)] = 1;
                double addDev = 0.0;
                for (int t : claimTris[bi]) {
                    if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
                    if (wasIn[static_cast<size_t>(t)]) continue;
                    for (int c = 0; c < 3; c++) {
                        const gp_XYZ pv = localTriVert(mv, t, c);
                        const gp_XYZ d = pv - curLoc;
                        addDev = std::max(addDev,
                                          std::abs((d - curAx * d.Dot(curAx)).Modulus() - curR));
                    }
                }
                std::fprintf(stderr,
                             "  DIAG_LAWABS rid=%d R=%.6f nTri %zu -> %zu absorbedMaxDev=%.6g "
                             "q=%.6g bandResid=%.6g tau=%.6g axisTiltRad=%.6g stop=%s "
                             "rounds=%d islands=%d latSkip=%d offCyl=%d "
                             "loc=(%.9g,%.9g,%.9g) dir=(%.9g,%.9g,%.9g)\n",
                             b.tris.empty() ? -1 : b.tris.front(), curR,
                             accepted[bi].tris.size(), claimTris[bi].size(), addDev, q,
                             b.maxVertResid, tau,
                             std::acos(std::min(1.0, std::abs(curAx.Dot(ax0)))), stopWhy,
                             nRounds, nIslands, nLatSkip, nOffCyl, curLoc.X(), curLoc.Y(),
                             curLoc.Z(),
                             curAx.X(), curAx.Y(), curAx.Z());
            }
            if (claimTris[bi].size() != accepted[bi].tris.size()) {
                bandLoc[bi] = curLoc;
                bandDir[bi] = gp_Dir(curAx.X(), curAx.Y(), curAx.Z());
                bandR[bi] = curR;
            }
        }
        if (lawbandDiagOn())
            std::fprintf(stderr,
                         "DIAG_LAWABSORB nTri=%d bands=%zu weldTol=%.6g sewTol=%.6g "
                         "diag=%.6g q=%.6g partScaleTauWas=%.6g\n",
                         nAbsorb, accepted.size(), mv.weldTol, mv.sewTol, mv.diag, q,
                         std::max({5e-5, 4.0 * mv.weldTol, 1e-6 * mv.diag}));
    }

    // D-130-13(2), second half: two BANDS that certified the same cylinder are
    // one surface even when no edge joins them. The seed sweep finds each
    // edge-connected piece of a cut wall separately (the cross bore's 108 and
    // 106; each R=3 slot end's two runs), and each fits its own cylinder to
    // 1e-7 -- but they are the same wall, and D-130-13(2) ships them as partial
    // faces of ONE surface, not as two surfaces 2e-07 apart. Merged only when
    // every vertex of each lies on the OTHER's certified cylinder within the
    // pair's own certificate, and only when the least-squares refit over the
    // union still lies inside it.
    if (lawUnionOn()) {
        const double qm = std::isfinite(mv.quantFloor) && mv.quantFloor > 0.0 ? mv.quantFloor : 0.0;
        auto bandAx = [&](size_t i) {
            return bandR[i] > 0.0
                       ? gp_Ax1(gp_Pnt(bandLoc[i].X(), bandLoc[i].Y(), bandLoc[i].Z()), bandDir[i])
                       : accepted[i].axis;
        };
        auto bandRad = [&](size_t i) { return bandR[i] > 0.0 ? bandR[i] : accepted[i].R; };
        auto axialSpan = [&](const std::vector<int>& ts, const gp_Ax1& ax, double& lo, double& hi) {
            const gp_XYZ o = ax.Location().XYZ();
            const gp_XYZ a = ax.Direction().XYZ();
            lo = 1e300;
            hi = -1e300;
            for (int t : ts) {
                for (int c = 0; c < 3; c++) {
                    const double x = (localTriVert(mv, t, c) - o).Dot(a);
                    lo = std::min(lo, x);
                    hi = std::max(hi, x);
                }
            }
        };
        auto axialOverlap = [&](const std::vector<int>& a, const std::vector<int>& b,
                                const gp_Ax1& ax) {
            double la = 0, ha = 0, lb = 0, hb = 0;
            axialSpan(a, ax, la, ha);
            axialSpan(b, ax, lb, hb);
            return la < hb && lb < ha;
        };
        auto allOn = [&](const std::vector<int>& ts, const gp_Ax1& ax, double R, double tau) {
            const gp_XYZ o = ax.Location().XYZ();
            const gp_XYZ a = ax.Direction().XYZ();
            for (int t : ts) {
                for (int c = 0; c < 3; c++) {
                    const gp_XYZ d = localTriVert(mv, t, c) - o;
                    if (std::abs((d - a * d.Dot(a)).Modulus() - R) > tau) return false;
                }
            }
            return true;
        };
        bool m2 = true;
        while (m2) {
            m2 = false;
            for (size_t i = 0; i < accepted.size() && !m2; i++) {
                if (claimTris[i].size() < 3 || !(bandRad(i) > 0.0)) continue;
                for (size_t j = i + 1; j < accepted.size(); j++) {
                    if (claimTris[j].size() < 3 || !(bandRad(j) > 0.0)) continue;
                    const double ti = std::max(accepted[i].maxVertResid, qm);
                    const double tj = std::max(accepted[j].maxVertResid, qm);
                    if (!allOn(claimTris[j], bandAx(i), bandRad(i), ti)) continue;
                    if (!allOn(claimTris[i], bandAx(j), bandRad(j), tj)) continue;
                    // Same SURFACE is not enough: two coaxial walls of equal
                    // radius in two stacked plates are two features, and
                    // D-130-13(2) is about the pieces ONE wall was cut into.
                    // The pieces of one wall interleave along the axis; stacked
                    // features do not. So the two claims' axial extents must
                    // overlap.
                    if (!axialOverlap(claimTris[i], claimTris[j], bandAx(i))) continue;
                    std::vector<int> trial = claimTris[i];
                    trial.insert(trial.end(), claimTris[j].begin(), claimTris[j].end());
                    std::sort(trial.begin(), trial.end());
                    trial.erase(std::unique(trial.begin(), trial.end()), trial.end());
                    gp_Dir fitAx;
                    gp_Pnt fitC;
                    double fitR = 0.0, resid = 0.0;
                    const gp_XYZ ai = bandAx(i).Direction().XYZ();
                    if (!cylinderFitLS(mv, trial, gp_Dir(ai.X(), ai.Y(), ai.Z()), fitAx, fitC,
                                       fitR, resid) ||
                        !(fitR > 0.0) || resid > std::max(ti, tj))
                        continue;
                    claimTris[i] = std::move(trial);
                    claimTris[j].clear();
                    bandLoc[i] = gp_XYZ(fitC.X(), fitC.Y(), fitC.Z());
                    bandDir[i] = fitAx;
                    bandR[i] = fitR;
                    if (lawbandDiagOn())
                        std::fprintf(stderr,
                                     "DIAG_LAWSAMESURF into=%zu from=%zu n=%zu R=%.9f "
                                     "resid=%.6g tau=%.6g\n",
                                     i, j, claimTris[i].size(), fitR, resid, std::max(ti, tj));
                    m2 = true;
                    break;
                }
            }
        }
    }

    for (size_t bi = 0; bi < accepted.size(); bi++) {
        LawBand b = accepted[bi];
        if (b.tris.size() < 3 || !(b.R > 0.0) || b.N < 2) continue;
        if (claimTris[bi].empty()) continue;  // merged into an earlier band
        b.tris = claimTris[bi];
        if (bandR[bi] > 0.0) {
            b.R = bandR[bi];
            b.axis = gp_Ax1(gp_Pnt(bandLoc[bi].X(), bandLoc[bi].Y(), bandLoc[bi].Z()),
                            bandDir[bi]);
        }
        // D-130-14 -- THE LOOP-LEVEL UNION, region half.
        //
        // The claim is ONE certified surface. Its triangles need not be
        // edge-connected: the interruptions that cut it (junction facets at a
        // crossing bore, the sawtooth of refused triangles round a mouth) sit
        // BETWEEN the pieces, and D-130-13(2)'s per-piece emit shipped each
        // piece as its own face -- eleven partial faces for handle-pickup's one
        // R=4 wall, which `partial_recovery_gate` refused (000c29f).
        //
        // A FACE lives on the surface, so its domain is measured on the
        // surface. Two claimed triangles that share a mesh VERTEX share a point
        // of the certified cylinder: they are one connected region there, and
        // the interruption between them is a hole in it, not a cut. Two
        // triangles that share nothing are two regions on the surface, and the
        // interruption between them really did cut it -- those, and only those,
        // ship as separate faces. So: one Region per connected component of the
        // claim under VERTEX adjacency, ascending minimum triangle id (I5).
        // The loops -- one outer/cap set plus one inner wire per enclosed
        // interruption -- are P1's (`buildTopologyD`), which pairs the boundary
        // chains through the interruption at every pinch vertex.
        std::vector<std::vector<int>> surfComps;
        if (!lawUnionOn()) {
            surfComps.push_back(b.tris);  // the single-Region emit, unchanged
        } else {
            std::vector<std::vector<int>> triAtVtx(mv.nVtx);
            for (int t : b.tris) {
                if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
                for (int sIdx = 0; sIdx < 3; sIdx++) {
                    const auto& e = mv.compEdges[mv.triEdges[t][sIdx]];
                    const int lv = ((mv.triDirs[t] >> sIdx) & 1) ? e.first : e.second;
                    if (lv >= 0 && static_cast<size_t>(lv) < mv.nVtx)
                        triAtVtx[static_cast<size_t>(lv)].push_back(t);
                }
            }
            std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
            for (int t0 : b.tris) {
                if (t0 < 0 || static_cast<size_t>(t0) >= mv.nTri) continue;
                if (seen[static_cast<size_t>(t0)]) continue;
                std::vector<int> comp, stk{t0};
                seen[static_cast<size_t>(t0)] = 1;
                while (!stk.empty()) {
                    const int x = stk.back();
                    stk.pop_back();
                    comp.push_back(x);
                    for (int sIdx = 0; sIdx < 3; sIdx++) {
                        const auto& e = mv.compEdges[mv.triEdges[x][sIdx]];
                        const int lv = ((mv.triDirs[x] >> sIdx) & 1) ? e.first : e.second;
                        if (lv < 0 || static_cast<size_t>(lv) >= mv.nVtx) continue;
                        for (int u2 : triAtVtx[static_cast<size_t>(lv)]) {
                            if (seen[static_cast<size_t>(u2)]) continue;
                            seen[static_cast<size_t>(u2)] = 1;
                            stk.push_back(u2);
                        }
                    }
                }
                std::sort(comp.begin(), comp.end());
                surfComps.push_back(std::move(comp));
            }
            std::sort(surfComps.begin(), surfComps.end(),
                      [](const std::vector<int>& x, const std::vector<int>& y) {
                          return x.front() < y.front();
                      });
        }
        if (lawbandDiagOn()) {
            const EdgeAdj eaP = buildEdgeAdj(mv);
            std::vector<char> inb(static_cast<size_t>(mv.nTri), 0);
            for (int t : b.tris)
                if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inb[static_cast<size_t>(t)] = 1;
            std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
            size_t nEdgePieces = 0;
            for (int t : b.tris) {
                if (t < 0 || static_cast<size_t>(t) >= mv.nTri || seen[static_cast<size_t>(t)])
                    continue;
                ++nEdgePieces;
                std::vector<int> stk{t};
                seen[static_cast<size_t>(t)] = 1;
                while (!stk.empty()) {
                    const int x = stk.back();
                    stk.pop_back();
                    for (int sIdx = 0; sIdx < 3 && mv.triEdges; sIdx++) {
                        const int e = mv.triEdges[x][sIdx];
                        const int u2 = (eaP.tri[e][0] == x) ? eaP.tri[e][1] : eaP.tri[e][0];
                        if (u2 < 0 || !inb[static_cast<size_t>(u2)] ||
                            seen[static_cast<size_t>(u2)])
                            continue;
                        seen[static_cast<size_t>(u2)] = 1;
                        stk.push_back(u2);
                    }
                }
            }
            // D-130-16: a surface component that is made of SEVERAL
            // edge-connected pieces is a PINCHED domain -- its pieces meet only
            // at isolated mesh vertices, and the face's boundary must pass
            // through each of those vertices twice. That is the number the
            // union's face path lives or dies on, so it is measured here and
            // not inferred: pieces=[k1,k2,...], one entry per surface face.
            std::fprintf(stderr,
                         "  DIAG_LAWUNION rid=%d R=%.6f nTri=%zu edgePieces=%zu "
                         "surfaceFaces=%zu sizes=[",
                         b.tris.empty() ? -1 : b.tris.front(), b.R, b.tris.size(), nEdgePieces,
                         surfComps.size());
            for (size_t k = 0; k < surfComps.size(); k++)
                std::fprintf(stderr, "%s%zu", k ? "," : "", surfComps[k].size());
            std::fprintf(stderr, "] pieces=[");
            for (size_t k = 0; k < surfComps.size(); k++) {
                std::vector<char> inc(static_cast<size_t>(mv.nTri), 0);
                for (int t : surfComps[k]) inc[static_cast<size_t>(t)] = 1;
                std::vector<char> sn(static_cast<size_t>(mv.nTri), 0);
                size_t nP = 0;
                for (int t : surfComps[k]) {
                    if (sn[static_cast<size_t>(t)]) continue;
                    ++nP;
                    std::vector<int> stk{t};
                    sn[static_cast<size_t>(t)] = 1;
                    while (!stk.empty()) {
                        const int x = stk.back();
                        stk.pop_back();
                        for (int sIdx = 0; sIdx < 3 && mv.triEdges; sIdx++) {
                            const int e2 = mv.triEdges[x][sIdx];
                            const int u3 = (eaP.tri[e2][0] == x) ? eaP.tri[e2][1] : eaP.tri[e2][0];
                            if (u3 < 0 || !inc[static_cast<size_t>(u3)] ||
                                sn[static_cast<size_t>(u3)])
                                continue;
                            sn[static_cast<size_t>(u3)] = 1;
                            stk.push_back(u3);
                        }
                    }
                }
                std::fprintf(stderr, "%s%zu", k ? "," : "", nP);
            }
            std::fprintf(stderr, "]\n");
        }
        const std::vector<int> whole = b.tris;
        for (const std::vector<int>& comp : surfComps) {
            LawBand pb = b;
            pb.tris = comp;
            Region reg;
            fillLawBandRegion(mv, tol, pb, reg);
            const int rid = reg.tris.empty() ? -1 : reg.tris.front();
            if (lawbandDiagOn()) {
                std::fprintf(stderr,
                             "DIAG_LAWCLAIM rid=%d n=%zu N=%d R=%.6f lawBand=1 merged=%d "
                             "split=%d surfaceFaces=%zu\n",
                             rid, reg.tris.size(), pb.N, pb.R, nMerged, nSplit,
                             surfComps.size());
            }
            work.accepted.push_back(std::move(reg));
        }
        peelLawBandFromProvisionals(mv, whole, work);
    }

    if (lawbandDiagOn())
        std::fprintf(stderr, "DIAG_LAWCLAIM_DONE accepted=%zu provisionals=%zu\n",
                     work.accepted.size(), work.provisionals.size());
    sortRegions(work.accepted);
    return true;
}

}}  // namespace stl2step::refit
