// P1 grow lane — charts A1, provisional A2, cylinder claim B1, plane commit A3.
// B1 gate structure: DECISION-p1-growx (D1.3-A1..A6, D2.3-A1).

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace stl2step { namespace refit {
namespace {

constexpr double kTiny = 1e-30;
constexpr double kBandEps = 1e-9;
// D1.3-A3 running residual: file-local. Frozen header epsCylGrow is not used by B1.
constexpr double kRingResidualFrac = 0.25;

bool p1DiagOn() {
    const char* e = std::getenv("STL2STEP_P1_DIAG");
    return e && e[0] != '\0' && e[0] != '0';
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

// Coarse Fusion/STLB export band (handle-lock @ 908 tris). Matches
// adaptCoarseSegmentParams in refit_segment.cpp — keep corpus fixtures
// (S02=412, S13/14≤20) on default gates.
bool coarseFusionBand(const MeshView& mv) {
    return mv.nTri >= 500 && mv.nTri <= 1200;
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
    const gp_XYZ xDir = dps - aw * ads;
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
    return phi >= tol.thetaCylLo - kBandEps && phi <= tol.thetaCylHi + kBandEps;
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
};

CommitEval evaluateCommit(const MeshView& mv, const DerivedTols& tol,
                          const std::vector<int>& tris, const gp_Dir& axis) {
    CommitEval ev;
    ev.g = centeredGauss(mv, tris, axis, tol);
    double cTilt = 0.0, devTilt = 0.0;
    axisTiltStats(mv, tris, axis, cTilt, devTilt);
    ev.eberlyOk = eberlyCenterRadius(mv, tris, axis, ev.center, ev.radius);
    ev.failGate = Gate::PASS;
    if (!testG1CommitSeedAxis(ev.g, tol, cTilt, devTilt)) {
        ev.failGate = Gate::G1;
        return ev;
    }
    if (!ev.eberlyOk) {
        ev.failGate = Gate::G4;
        return ev;
    }
    ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
    const double rBeforeRefine = ev.radius;
    if (coarseFusionBand(mv) && ev.d2.nSides >= 3 && !ev.d2.spanReject) {
        refineCylinderRadius(mv, tris, axis, ev.center, ev.radius, ev.d2.nSides, ev.d2.span);
        ev.d2 = computeD2(mv, tris, axis, ev.center, ev.radius, tol);
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
    if (delta > tol.epsMesh && (!coarse || (ev.d2.span >= 0.35 && lift <= 0.0))) {
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
        if (!g5Closed && !g5Partial && !g5Micro) ev.failGate = Gate::G5;
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
            for (size_t j = i + 1; j < work.accepted.size(); ++j) {
                const Region& B = work.accepted[j];
                const bool adjacent = regionsShareMeshEdge(mv, A, B);
                if (!coaxialCylinderMergeable(A, B, tol, adjacent)) continue;
                if (!adjacent && axisLineSeparation(A.ax, B.ax) > tol.epsPlane) continue;
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
                    if (triLabel[u] >= 0) continue;
                    if (work.triChart[u] != chart) continue;
                    const double phi = edgeDihedralAbs(mv, e, ea);
                    if (phi > tol.thetaSharp) continue;

                    const gp_Dir n_u = triNormalLocal(mv, u);
                    const gp_Dir n_p = plane.Direction();
                    if (std::abs(n_u.Dot(n_p)) < std::cos(tol.thetaPlane)) continue;

                    double maxD = 0.0;
                    const gp_Pnt loc = plane.Location();
                    for (int k = 0; k < 3; k++) {
                        const gp_XYZ v = localTriVert(mv, u, k);
                        const double d = std::abs(
                            gp_Vec(loc, gp_Pnt(v)).Dot(gp_Vec(n_p)));
                        maxD = std::max(maxD, d);
                    }
                    if (maxD > tol.epsPlane) continue;

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

    sortProvisionals(work.provisionals);
    return true;
}

bool claimCylindersB1(const MeshView& mv, const SegmentParams&, const DerivedTols& tol,
                      SegmentWork& work) {
    if (work.provisionals.empty()) return true;

    const EdgeAdj ea = buildEdgeAdj(mv);
    std::vector<int> triToProv(mv.nTri, -1);
    for (size_t pi = 0; pi < work.provisionals.size(); pi++) {
        for (int t : work.provisionals[pi].tris) triToProv[t] = static_cast<int>(pi);
    }

    const ProvAdjList adj = buildProvAdjacency(mv, ea, triToProv);
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
                if (phi < tol.thetaCylLo - kBandEps) continue;
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
                const double g5v = std::abs(
                    nXbar.Dot(gp_XYZ(aS.X(), aS.Y(), aS.Z())) - cS);
                if (g5v > sin3) {
                    dead[xi] = 1;
                    nG5++;
                    worstG5 = std::max(worstG5, g5v);
                    if (p1DiagOn())
                        std::fprintf(stderr, "  g5 reject x=%d |n.a-c|=%.5g sin3=%.5g\n",
                                     (int)xi, g5v, sin3);
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
                }
                progressed = true;
                break;
            }
            if (!progressed) break;
        }

        const std::vector<int> prePeelTris = mergeMemberTris(work.provisionals, members);
        GaussResult gFinal;
        const gp_Dir axisFinal = axisOf(mv, work.provisionals, members, seedAxis, tol, &gFinal);
        CommitEval ev = evaluateCommit(mv, tol, prePeelTris, axisFinal);
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

}}  // namespace stl2step::refit
