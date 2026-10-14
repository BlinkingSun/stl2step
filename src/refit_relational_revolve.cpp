// Relational revolve (D-fc-1.3 … 1.7, D-fc-1.10). Line mode claims closed
// cones; Torus mode claims tori. Axis comes from a donor refit and is not
// re-fit on the band. Certificate is τ = tauOf(mv). Never throws.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_profile_cert.hpp"
#include "refit_cone_math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>
#include <Standard_Failure.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

bool diagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
}

gp_XYZ triVert(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.pts[mv.tris[gt][k]];
}

double triAreaNormal(const MeshView& mv, int t, gp_XYZ& nOut) {
    const gp_XYZ a = triVert(mv, t, 0);
    const gp_XYZ b = triVert(mv, t, 1);
    const gp_XYZ c = triVert(mv, t, 2);
    const gp_XYZ n = (b - a).Crossed(c - a);
    const double mag = n.Modulus();
    if (mag > gp::Resolution()) nOut = n / mag;
    else nOut = gp_XYZ(0, 0, 1);
    return 0.5 * mag;
}

gp_XYZ triCentroid(const MeshView& mv, int t) {
    return (triVert(mv, t, 0) + triVert(mv, t, 1) + triVert(mv, t, 2)) / 3.0;
}

struct Frame {
    gp_XYZ origin;
    gp_XYZ axis;
    gp_XYZ e1;
    gp_XYZ e2;
};

Frame makeFrame(const gp_XYZ& origin, const gp_XYZ& axisIn, const gp_XYZ& xHint) {
    Frame f;
    f.origin = origin;
    f.axis = axisIn;
    const double am = f.axis.Modulus();
    if (am > gp::Resolution()) f.axis /= am;
    else f.axis = gp_XYZ(0, 0, 1);
    gp_XYZ x = xHint - f.axis * xHint.Dot(f.axis);
    if (!(x.Modulus() > gp::Resolution())) {
        const gp_XYZ tmp = (std::abs(f.axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
        x = f.axis.Crossed(tmp);
    }
    f.e1 = x / x.Modulus();
    f.e2 = f.axis.Crossed(f.e1);
    return f;
}

ProfPt profOfPoint(const gp_XYZ& p, const Frame& f, int lv) {
    ProfPt o;
    const gp_XYZ d = p - f.origin;
    o.lv = lv;
    o.z = d.Dot(f.axis);
    const double x = d.Dot(f.e1);
    const double y = d.Dot(f.e2);
    o.rho = std::hypot(x, y);
    o.psi = std::atan2(y, x);
    if (o.psi < 0.0) o.psi += kTwoPi;
    return o;
}

struct Vtx {
    int lv = -1;
    int tri = -1;
    gp_XYZ p;
    ProfPt pr;
};

struct SurfaceFit {
    bool line = false;
    bool circle = false;
    double a = 0, b = 0;          // rho = a + b z, when regressOk
    bool regressRho = false;      // dependent variable is rho (else z)
    double orthMax = 0;
    double Rmaj = 0, z0 = 0, Rmin = 0;
    double circMax = 0;
    // Regression normal matrix for rho = a + b z, or z = a + b rho.
    double n = 0, sX = 0, sXX = 0, det = 0;
    bool hatOk = false;
};

double hypot2(double x, double y) { return std::hypot(x, y); }

SurfaceFit fitSurface(const std::vector<ProfPt>& P, double tau) {
    SurfaceFit f;
    if (P.size() < 2) return f;
    f.orthMax = profileLineMaxResid(P);
    f.line = f.orthMax <= tau;
    double sZ = 0, sR = 0, sZZ = 0, sRR = 0, sZR = 0;
    const double n = static_cast<double>(P.size());
    double zMin = P[0].z, zMax = P[0].z, rMin = P[0].rho, rMax = P[0].rho;
    for (const ProfPt& p : P) {
        sZ += p.z;
        sR += p.rho;
        sZZ += p.z * p.z;
        sRR += p.rho * p.rho;
        sZR += p.z * p.rho;
        zMin = std::min(zMin, p.z);
        zMax = std::max(zMax, p.z);
        rMin = std::min(rMin, p.rho);
        rMax = std::max(rMax, p.rho);
    }
    const double zSpan = zMax - zMin;
    const double rSpan = rMax - rMin;
    // rho = a + b z when z is the long direction; otherwise z = a + b rho.
    f.regressRho = !(rSpan > zSpan && zSpan <= tau);
    if (f.regressRho && zSpan > gp::Resolution()) {
        const double det = n * sZZ - sZ * sZ;
        if (std::abs(det) > gp::Resolution()) {
            f.b = (n * sZR - sZ * sR) / det;
            f.a = (sR - f.b * sZ) / n;
            f.n = n;
            f.sX = sZ;
            f.sXX = sZZ;
            f.det = det;
            f.hatOk = true;
        }
    } else if (!f.regressRho && rSpan > gp::Resolution()) {
        const double det = n * sRR - sR * sR;
        if (std::abs(det) > gp::Resolution()) {
            f.b = (n * sZR - sR * sZ) / det;
            f.a = (sZ - f.b * sR) / n;
            f.n = n;
            f.sX = sR;
            f.sXX = sRR;
            f.det = det;
            f.hatOk = true;
        }
    } else if (zSpan <= tau) {
        // Horizontal: z = mean. Hat row of a mean is 1/n.
        f.a = sZ / n;
        f.b = 0;
        f.regressRho = false;
        f.n = n;
        f.sX = 0;
        f.sXX = 0;
        f.det = n;  // sentinel: constant model
        f.hatOk = n >= 1.0;
    } else {
        // Vertical: rho = mean.
        f.a = sR / n;
        f.b = 0;
        f.regressRho = true;
        f.n = n;
        f.sX = 0;
        f.sXX = 0;
        f.det = n;
        f.hatOk = n >= 1.0;
    }
    if (P.size() >= 3 && fitProfileCircle(P, tau, f.Rmaj, f.z0, f.Rmin)) {
        f.circle = true;
        f.circMax = 0;
        for (const ProfPt& p : P) {
            const double d = std::fabs(hypot2(p.rho - f.Rmaj, p.z - f.z0) - f.Rmin);
            f.circMax = std::max(f.circMax, d);
        }
    }
    return f;
}

double dependentOf(const SurfaceFit& f, const ProfPt& p) {
    return f.regressRho ? p.z : p.rho;
}

double predictedDependent(const SurfaceFit& f, double indep) {
    return f.a + f.b * indep;
}

double radialMiss(const SurfaceFit& f, const ProfPt& p) {
    if (!f.hatOk) return std::numeric_limits<double>::infinity();
    if (f.det == f.n && f.b == 0.0 && f.sXX == 0.0) {
        return f.regressRho ? std::fabs(p.rho - f.a) : std::fabs(p.z - f.a);
    }
    if (f.regressRho) return std::fabs(p.rho - (f.a + f.b * p.z));
    return std::fabs(p.z - (f.a + f.b * p.rho));
}

// ℓ1 of the hat row that maps fitted samples onto the predicted dependent
// variable at `p`. Constant model: every weight is 1/n, ℓ1 = 1.
double hatL1(const SurfaceFit& f, const std::vector<ProfPt>& samples, const ProfPt& p) {
    if (!f.hatOk || samples.empty()) return std::numeric_limits<double>::infinity();
    if (f.sXX == 0.0 && f.b == 0.0) return 1.0;
    const double xq = dependentOf(f, p);
    const double det = f.det;
    if (!(std::abs(det) > gp::Resolution())) return std::numeric_limits<double>::infinity();
    double l1 = 0.0;
    for (const ProfPt& s : samples) {
        const double xi = dependentOf(f, s);
        const double h = ((f.sXX - xq * f.sX) + xi * (f.n * xq - f.sX)) / det;
        l1 += std::abs(h);
    }
    return l1;
}

double extrapOf(double q, double l1) { return q * l1; }

// Perpendicular profile distance from the regression line, so growth and the
// τ certificate speak about the same surface distance.
double linePerp(const SurfaceFit& f, const ProfPt& p) {
    const double miss = radialMiss(f, p);
    if (!std::isfinite(miss)) return miss;
    const double scale = std::sqrt(1.0 + f.b * f.b);
    return miss / scale;
}

double growTolLine(const SurfaceFit& f, const std::vector<ProfPt>& samples, const ProfPt& p,
                   double tau, double q) {
    const double scale = std::sqrt(1.0 + f.b * f.b);
    const double prop = extrapOf(q, hatL1(f, samples, p)) / scale;
    return tau + prop;
}

bool solve3(const double m[3][3], const double rhs[3], double x[3]) {
    double a[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) a[i][j] = m[i][j];
        a[i][3] = rhs[i];
    }
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int row = col + 1; row < 3; ++row)
            if (std::abs(a[row][col]) > std::abs(a[piv][col])) piv = row;
        if (!(std::abs(a[piv][col]) > gp::Resolution())) return false;
        if (piv != col)
            for (int j = col; j < 4; ++j) std::swap(a[piv][j], a[col][j]);
        const double div = a[col][col];
        for (int j = col; j < 4; ++j) a[col][j] /= div;
        for (int row = 0; row < 3; ++row) {
            if (row == col) continue;
            const double w = a[row][col];
            for (int j = col; j < 4; ++j) a[row][j] -= w * a[col][j];
        }
    }
    for (int i = 0; i < 3; ++i) x[i] = a[i][3];
    return true;
}

// Geometric-circle hat row. Residual is already the profile distance.
double circleHatL1(const std::vector<ProfPt>& samples, double Rmaj, double z0, double Rmin,
                   const ProfPt& q) {
    if (samples.size() < 3) return std::numeric_limits<double>::infinity();
    double G[3][3] = {};
    std::vector<std::array<double, 3>> J;
    J.reserve(samples.size());
    for (const ProfPt& s : samples) {
        const double dx = s.rho - Rmaj, dy = s.z - z0;
        const double d = std::hypot(dx, dy);
        if (!(d > gp::Resolution())) return std::numeric_limits<double>::infinity();
        std::array<double, 3> j{{-(dx / d), -(dy / d), -1.0}};
        J.push_back(j);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) G[a][b] += j[a] * j[b];
    }
    const double dx = q.rho - Rmaj, dy = q.z - z0;
    const double d = std::hypot(dx, dy);
    if (!(d > gp::Resolution())) return std::numeric_limits<double>::infinity();
    const double jq[3] = {-(dx / d), -(dy / d), -1.0};
    double l1 = 0.0;
    for (const auto& j : J) {
        double rhs[3] = {j[0], j[1], j[2]};
        double col[3];
        if (!solve3(G, rhs, col)) return std::numeric_limits<double>::infinity();
        const double h = jq[0] * col[0] + jq[1] * col[1] + jq[2] * col[2];
        l1 += std::abs(h);
    }
    return l1;
}

double circDist(const SurfaceFit& f, const ProfPt& p) {
    if (!f.circle) return std::numeric_limits<double>::infinity();
    return std::fabs(hypot2(p.rho - f.Rmaj, p.z - f.z0) - f.Rmin);
}

struct Band {
    std::vector<int> tris;
    std::vector<char> in;
    std::vector<Vtx> verts;
    std::vector<ProfPt> pts;
    int sheet = 0;
    int seed = -1;
};

int sheetOf(const MeshView& mv, int t, const Frame& f) {
    gp_XYZ n;
    triAreaNormal(mv, t, n);
    const gp_XYZ c = triCentroid(mv, t);
    const gp_XYZ d = c - f.origin;
    gp_XYZ rad = d - f.axis * d.Dot(f.axis);
    const double rm = rad.Modulus();
    if (!(rm > gp::Resolution())) return 0;
    rad /= rm;
    const double s = n.Dot(rad);
    if (std::abs(s) <= gp::Resolution()) return 0;
    return s > 0.0 ? 1 : -1;
}

void collectBand(const MeshView& mv, const Frame& f, const std::vector<int>& g2l, Band& b) {
    b.verts.clear();
    b.pts.clear();
    std::vector<int> seen;
    for (int t : b.tris) {
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int g = mv.tris[gt][k];
            int lv = (g >= 0 && g < (int)g2l.size()) ? g2l[g] : -1;
            if (lv < 0) lv = g;
            seen.push_back(lv);
            Vtx v;
            v.lv = lv;
            v.tri = t;
            v.p = mv.pts[g];
            v.pr = profOfPoint(v.p, f, lv);
            b.verts.push_back(v);
        }
    }
    std::sort(seen.begin(), seen.end());
    seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
    std::vector<char> took(b.verts.size(), 0);
    for (int lv : seen) {
        for (std::size_t i = 0; i < b.verts.size(); ++i) {
            if (took[i]) continue;
            if (b.verts[i].lv == lv) {
                b.pts.push_back(b.verts[i].pr);
                took[i] = 1;
                break;
            }
        }
    }
    std::sort(b.tris.begin(), b.tris.end());
    b.tris.erase(std::unique(b.tris.begin(), b.tris.end()), b.tris.end());
}

struct ClosedM {
    int nBands = 0;
    double maxGap = kTwoPi;
    double span = 0;
    int nSides = 0;
    bool closed = false;
    double uMin = 0, uMax = 0;
};

double wrapPi(double t) {
    while (t <= -kPi) t += kTwoPi;
    while (t > kPi) t -= kTwoPi;
    return t;
}

ClosedM closedOf(const MeshView& mv, const Frame& f, const std::vector<int>& tris,
                 const std::vector<Vtx>& uniq, double thetaBin) {
    ClosedM m;
    if (tris.empty()) return m;
    std::vector<double> psi;
    psi.reserve(tris.size());
    for (int t : tris) {
        gp_XYZ n;
        triAreaNormal(mv, t, n);
        const double pu = n.Dot(f.e1);
        const double pv = n.Dot(f.e2);
        psi.push_back(wrapPi(std::atan2(pv, pu)));
    }
    // Same gap clustering as chamfer's nBandsAzimuth (existing column rule).
    int nBands = 0;
    if (!psi.empty()) {
        std::sort(psi.begin(), psi.end());
        const int n = (int)psi.size();
        if (n == 1) {
            nBands = 1;
        } else {
            std::vector<double> gaps((std::size_t)n);
            for (int i = 0; i < n - 1; ++i)
                gaps[(std::size_t)i] = psi[(std::size_t)i + 1] - psi[(std::size_t)i];
            gaps[(std::size_t)n - 1] = (psi[0] + kTwoPi) - psi[(std::size_t)n - 1];
            int imax = 0;
            for (int i = 1; i < n; ++i)
                if (gaps[(std::size_t)i] > gaps[(std::size_t)imax]) imax = i;
            std::vector<double> rest;
            rest.reserve((std::size_t)n - 1);
            for (int i = 0; i < n; ++i)
                if (i != imax) rest.push_back(gaps[(std::size_t)i]);
            double p75 = 0.0;
            if (!rest.empty()) {
                std::sort(rest.begin(), rest.end());
                const std::size_t idx = (rest.size() * 3u) / 4u;
                p75 = rest[idx < rest.size() ? idx : rest.size() - 1u];
            }
            const double th = std::max(thetaBin, 0.5 * p75);
            const int start = (imax + 1) % n;
            nBands = 1;
            double prev = psi[(std::size_t)start];
            for (int k = 1; k < n; ++k) {
                const int idx = (start + k) % n;
                const double cur = psi[(std::size_t)idx];
                double gap = cur - prev;
                if (gap < 0.0) gap += kTwoPi;
                if (gap > th) ++nBands;
                prev = cur;
            }
        }
    }
    m.nBands = nBands;
    std::vector<double> chi;
    for (const Vtx& v : uniq) {
        const gp_XYZ d = v.p - f.origin;
        chi.push_back(wrapPi(std::atan2(d.Dot(f.e2), d.Dot(f.e1))));
    }
    std::sort(chi.begin(), chi.end());
    chi.erase(std::unique(chi.begin(), chi.end()), chi.end());
    if (chi.size() < 2) return m;
    int jmax = static_cast<int>(chi.size()) - 1;
    double maxGap = kTwoPi - chi.back() + chi.front();
    for (std::size_t j = 0; j + 1 < chi.size(); ++j) {
        const double g = chi[j + 1] - chi[j];
        if (g > maxGap) {
            maxGap = g;
            jmax = static_cast<int>(j);
        }
    }
    m.maxGap = maxGap;
    const double bandArc = (nBands > 0) ? kTwoPi / static_cast<double>(nBands) : kTwoPi;
    m.closed = nBands >= 3 && maxGap <= 1.5 * bandArc;
    m.span = m.closed ? kTwoPi : (kTwoPi - maxGap);
    if (m.span > 0.0 && nBands >= 1)
        m.nSides = static_cast<int>(std::llround(kTwoPi * nBands / m.span));
    if (m.closed) {
        m.uMin = 0;
        m.uMax = kTwoPi;
    } else {
        const std::size_t uIdx = static_cast<std::size_t>(jmax + 1) % chi.size();
        m.uMin = wrapPi(chi[uIdx]);
        m.uMax = m.uMin + m.span;
    }
    return m;
}

bool stationsConcyclic(const std::vector<ProfPt>& P, double tau, const SurfaceFit& fit) {
    if (P.empty() || !fit.hatOk) return false;
    std::vector<ProfPt> s = P;
    std::sort(s.begin(), s.end(), [](const ProfPt& a, const ProfPt& b) { return a.z < b.z; });
    std::size_t i = 0;
    while (i < s.size()) {
        std::size_t j = i + 1;
        while (j < s.size() && s[j].z - s[i].z <= tau) ++j;
        double rMin = s[i].rho, rMax = s[i].rho;
        for (std::size_t k = i; k < j; ++k) {
            rMin = std::min(rMin, s[k].rho);
            rMax = std::max(rMax, s[k].rho);
        }
        // A station bin τ wide on a slanted generator moves ρ by |b|·τ.
        const double allow = tau * std::sqrt(1.0 + fit.b * fit.b) + tau;
        if (rMax - rMin > allow) return false;
        for (std::size_t k = i; k < j; ++k)
            if (linePerp(fit, s[k]) > tau) return false;
        i = j;
    }
    return true;
}

struct Adj {
    std::vector<std::array<int, 2>> edgeTris;
    std::vector<std::vector<int>> nbr;
};

Adj buildAdj(const MeshView& mv) {
    Adj a;
    const int nTri = (int)mv.nTri;
    const int nEdge = (int)mv.nEdge;
    a.edgeTris.assign((std::size_t)std::max(nEdge, 0), std::array<int, 2>{-1, -1});
    a.nbr.assign((std::size_t)std::max(nTri, 0), {});
    if (!mv.triEdges) return a;
    for (int t = 0; t < nTri; ++t) {
        for (int s = 0; s < 3; ++s) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || e >= nEdge) continue;
            if (a.edgeTris[(std::size_t)e][0] < 0) a.edgeTris[(std::size_t)e][0] = t;
            else if (a.edgeTris[(std::size_t)e][1] < 0) a.edgeTris[(std::size_t)e][1] = t;
        }
    }
    for (int e = 0; e < nEdge; ++e) {
        const int t0 = a.edgeTris[(std::size_t)e][0];
        const int t1 = a.edgeTris[(std::size_t)e][1];
        if (t0 < 0 || t1 < 0 || t0 == t1) continue;
        a.nbr[(std::size_t)t0].push_back(t1);
        a.nbr[(std::size_t)t1].push_back(t0);
    }
    for (auto& v : a.nbr) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return a;
}

struct Donor {
    int idx = -1;
    Frame fr;
    double radius = 0;
    double maxDev = 0;
    double rms = 0;
    double length = 0;
    int minTri = 0;
    bool updated = false;
};

bool regionIsDonorClass(const Region& r) {
    if (r.type == SurfType::Cylinder) return true;
    return r.type == SurfType::Cone && r.origin == Origin::ChamferCone;
}

double axisGap(const Frame& a, const Frame& b) {
    const double cross = a.axis.Crossed(b.axis).Modulus();
    const gp_XYZ d = b.origin - a.origin;
    const double along = d.Dot(a.axis);
    const double dist = (d - a.axis * along).Modulus();
    return dist + cross;  // radians mixed with mm — callers compare parts separately
}

bool sameAxis(const Frame& a, const Frame& b, double tau, double length) {
    const double cross = a.axis.Crossed(b.axis).Modulus();
    const double L = std::max(length, tau);
    const double ang = tanAngTol(tau, tau, L);
    if (cross > ang) return false;
    const gp_XYZ d = b.origin - a.origin;
    const double along = d.Dot(a.axis);
    const double dist = (d - a.axis * along).Modulus();
    return dist <= tau;
}

struct DiagRec {
    const char* mode = "line";
    int seed = -1;
    int nTri = 0;
    int nV = 0;
    int stations = 0;
    int levels = 0;
    int params = 0;
    const char* fit = "anchored";
    double r0 = 0, r1 = 0, halfAng = 0;
    gp_XYZ axis;
    double maxDev = 0, maxDevTau = 0, lineResTau = 0;
    double span = 0, dihedral = 0, spanTol = 0;
    double tan0 = 0, tan1 = 0;
    int straddle = 0;
    const char* why = "underdetermined";
};

void printBand(const DiagRec& d) {
    if (!diagOn()) return;
    std::fprintf(stderr,
                 "DIAG_REL_BAND stage=revolve mode=%s seed=%d nTri=%d nV=%d stations=%d "
                 "levels=%d params=%d fit=%s radii=%.8g/%.8g halfAng=%.6f "
                 "axis=%.6g,%.6g,%.6g maxDev=%.6g maxDevTau=%.4f lineResTau=%.4f "
                 "span=%.6g dihedral=%.6g spanTol=%.6g tan=%.4f,%.4f straddle=%d why=%s\n",
                 d.mode, d.seed, d.nTri, d.nV, d.stations, d.levels, d.params, d.fit, d.r0,
                 d.r1, d.halfAng, d.axis.X(), d.axis.Y(), d.axis.Z(), d.maxDev, d.maxDevTau,
                 d.lineResTau, d.span, d.dihedral, d.spanTol, d.tan0, d.tan1, d.straddle,
                 d.why);
}

void printSum(const char* mode, int donors, int claimed, int refused) {
    if (!diagOn()) return;
    std::fprintf(stderr, "DIAG_REL_SUM stage=revolve mode=%s donors=%d claimed=%d refused=%d\n",
                 mode, donors, claimed, refused);
}

struct WorkMesh {
    const MeshView* mv = nullptr;
    Adj adj;
    std::vector<int> g2l;
    std::vector<int> triProv;
    std::vector<char> freeTri;
    double tau = 0, q = 0;
    double thetaBin = 0;
    int maxIter = 8;
};

int localId(const WorkMesh& wm, int t, int k) {
    const MeshView& mv = *wm.mv;
    const int gt = mv.compTris ? mv.compTris[t] : t;
    const int g = mv.tris[gt][k];
    if (g >= 0 && g < (int)wm.g2l.size() && wm.g2l[g] >= 0) return wm.g2l[g];
    return g;
}

bool vertOkLine(const SurfaceFit& fit, const std::vector<ProfPt>& samples, const ProfPt& p,
                double tau, double q) {
    return linePerp(fit, p) <= growTolLine(fit, samples, p, tau, q);
}

std::vector<int> frontierOf(const WorkMesh& wm, const Band& b) {
    std::vector<int> out;
    const int nTri = (int)wm.mv->nTri;
    for (int t : b.tris) {
        if (t < 0 || t >= nTri) continue;
        for (int n : wm.adj.nbr[(std::size_t)t]) {
            if (n < 0 || n >= nTri) continue;
            if (b.in[(std::size_t)n]) continue;
            if (!wm.freeTri[(std::size_t)n]) continue;
            out.push_back(n);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// Maximal same-sheet free component. Returns false when it is the band itself.
bool largerSheet(const WorkMesh& wm, const Frame& fr, const Band& band, Band& ext) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> seen((std::size_t)std::max(nTri, 0), 0);
    std::vector<int> stack;
    int nBand = 0;
    for (int t : band.tris) {
        if (t < 0 || t >= nTri || seen[(std::size_t)t]) continue;
        seen[(std::size_t)t] = 1;
        stack.push_back(t);
        ++nBand;
    }
    std::vector<int> comp;
    while (!stack.empty()) {
        const int t = stack.back();
        stack.pop_back();
        comp.push_back(t);
        for (int n : wm.adj.nbr[(std::size_t)t]) {
            if (n < 0 || n >= nTri || seen[(std::size_t)n]) continue;
            if (!wm.freeTri[(std::size_t)n]) continue;
            if (sheetOf(*wm.mv, n, fr) != band.sheet) continue;
            seen[(std::size_t)n] = 1;
            stack.push_back(n);
        }
    }
    if ((int)comp.size() <= nBand) return false;
    ext = Band();
    ext.in.assign((std::size_t)std::max(nTri, 0), 0);
    ext.tris = std::move(comp);
    for (int t : ext.tris)
        if (t >= 0 && t < nTri) ext.in[(std::size_t)t] = 1;
    collectBand(*wm.mv, fr, wm.g2l, ext);
    return true;
}

// A line band that is only a subset of a certifyTorus sheet is a fillet
// chord, not a cone (D-fc-1.6.9). A sheet that is still a line is a cone.
bool lineIsCircleSubset(const WorkMesh& wm, const Frame& fr, const Band& band, double tau) {
    Band ext;
    if (!largerSheet(wm, fr, band, ext)) return false;
    if ((int)ext.pts.size() < minStations(kParamsCircle)) return false;
    if (profileLineMaxResid(ext.pts) <= tau) return false;
    return certifyTorus(ext.pts, tau).ok;
}

struct ConeCand {
    Band band;
    SurfaceFit fit;
    ClosedM cl;
    double maxP = 0;
};

bool ptsOnCircle(const std::vector<ProfPt>& pts, double Rmaj, double z0, double Rmin, double tau) {
    if (pts.empty()) return false;
    for (const ProfPt& p : pts)
        if (std::fabs(std::hypot(p.rho - Rmaj, p.z - z0) - Rmin) > tau) return false;
    return true;
}

// A 2-station line band that grows, through free neighbours, into a
// certifyTorus profile is one row of that blend (D-fc-1.6.9). A cone that
// only meets the blend does not lie on the circle, so the grow stops.
bool rowOfTorus(const WorkMesh& wm, const Frame& fr, const Band& band, double tau) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> in((std::size_t)std::max(nTri, 0), 0);
    for (int t : band.tris)
        if (t >= 0 && t < nTri) in[(std::size_t)t] = 1;
    std::vector<ProfPt> pts = band.pts;
    bool grew = false;
    for (int guard = 0; guard < nTri; ++guard) {
        std::vector<int> frnt;
        for (int t = 0; t < nTri; ++t) {
            if (!in[(std::size_t)t]) continue;
            for (int n : wm.adj.nbr[(std::size_t)t]) {
                if (n < 0 || n >= nTri || in[(std::size_t)n] || !wm.freeTri[(std::size_t)n]) continue;
                if (sheetOf(*wm.mv, n, fr) != band.sheet) continue;
                frnt.push_back(n);
            }
        }
        std::sort(frnt.begin(), frnt.end());
        frnt.erase(std::unique(frnt.begin(), frnt.end()), frnt.end());
        if (frnt.empty()) break;
        std::vector<char> infr((std::size_t)nTri, 0);
        for (int t : frnt) infr[(std::size_t)t] = 1;
        bool any = false;
        std::vector<char> seen((std::size_t)nTri, 0);
        for (int seed : frnt) {
            if (seen[(std::size_t)seed]) continue;
            std::vector<int> comp;
            std::vector<int> stack{seed};
            seen[(std::size_t)seed] = 1;
            while (!stack.empty()) {
                const int t = stack.back();
                stack.pop_back();
                comp.push_back(t);
                for (int n : wm.adj.nbr[(std::size_t)t]) {
                    if (n < 0 || n >= nTri || !infr[(std::size_t)n] || seen[(std::size_t)n]) continue;
                    seen[(std::size_t)n] = 1;
                    stack.push_back(n);
                }
            }
            std::vector<ProfPt> trial = pts;
            for (int t : comp)
                for (int k = 0; k < 3; ++k)
                    trial.push_back(profOfPoint(triVert(*wm.mv, t, k), fr, 0));
            double Rmaj = 0, z0 = 0, Rmin = 0;
            if (!fitProfileCircle(trial, tau, Rmaj, z0, Rmin)) continue;
            if (!ptsOnCircle(trial, Rmaj, z0, Rmin, tau)) continue;
            pts.swap(trial);
            for (int t : comp) in[(std::size_t)t] = 1;
            any = true;
            grew = true;
        }
        if (!any) break;
    }
    return grew && certifyTorus(pts, tau).ok;
}

void markTorusChords(const WorkMesh& wm, const Frame& fr, const std::vector<ConeCand>& cs, double tau,
                     std::vector<char>& drop) {
    drop.assign(cs.size(), 0);
    for (int i = 0; i < (int)cs.size(); ++i)
        if (rowOfTorus(wm, fr, cs[(std::size_t)i].band, tau)) drop[(std::size_t)i] = 1;
}

// A torus band that is a strict subset of a sheet certifyTorus refuses is
// the plate-mouth window (D-fc-1.6.9). The sheet that is the band is kept.
bool torusIsFailingSubset(const WorkMesh& wm, const Frame& fr, const Band& band, double tau) {
    Band ext;
    if (!largerSheet(wm, fr, band, ext)) return false;
    return !certifyTorus(ext.pts, tau).ok;
}

void growLineDrain(WorkMesh& wm, const Frame& fr, Band& band, const std::vector<int>& stamp) {
    const int nTri = (int)wm.mv->nTri;
    for (int round = 0; round < wm.maxIter; ++round) {
        collectBand(*wm.mv, fr, wm.g2l, band);
        const SurfaceFit fit = fitSurface(band.pts, wm.tau);
        if (!fit.hatOk || !fit.line) break;
        bool any = false;
        for (;;) {
            std::vector<int> add;
            for (int t : band.tris) {
                if (t < 0 || t >= nTri) continue;
                for (int n : wm.adj.nbr[(std::size_t)t]) {
                    if (n < 0 || n >= nTri || band.in[(std::size_t)n]) continue;
                    if (!wm.freeTri[(std::size_t)n] || stamp[(std::size_t)n] != 0) continue;
                    if (sheetOf(*wm.mv, n, fr) != band.sheet) continue;
                    bool ok = true;
                    for (int k = 0; k < 3 && ok; ++k) {
                        const ProfPt pp = profOfPoint(triVert(*wm.mv, n, k), fr, 0);
                        if (!vertOkLine(fit, band.pts, pp, wm.tau, wm.q)) ok = false;
                    }
                    if (ok) add.push_back(n);
                }
            }
            std::sort(add.begin(), add.end());
            add.erase(std::unique(add.begin(), add.end()), add.end());
            if (add.empty()) break;
            for (int n : add) {
                band.in[(std::size_t)n] = 1;
                band.tris.push_back(n);
            }
            any = true;
        }
        if (!any) break;
    }
    collectBand(*wm.mv, fr, wm.g2l, band);
}

// Torus growth. While the circle fit cannot yet resolve τ (extrapBound > τ),
// or the profile is still a line, take the next same-sheet layer. After the
// fit tightens, admit only vertices inside growTol. The caller certifies the
// whole grown band once (D-fc-1.6.6, D-fc-1.6.9).
void growTorusDrain(WorkMesh& wm, const Frame& fr, Band& band, const std::vector<int>& stamp) {
    const int nTri = (int)wm.mv->nTri;
    for (int guard = 0; guard < nTri; ++guard) {
        collectBand(*wm.mv, fr, wm.g2l, band);
        const double lineRes = profileLineMaxResid(band.pts);
        const bool line = lineRes <= wm.tau;
        const SurfaceFit fit = fitSurface(band.pts, wm.tau);
        const bool circ = fit.circle && fit.circMax <= wm.tau && !line;
        double probeL1 = std::numeric_limits<double>::infinity();
        int probe = -1;
        for (int t : band.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int n : wm.adj.nbr[(std::size_t)t]) {
                if (n < 0 || n >= nTri || band.in[(std::size_t)n]) continue;
                if (!wm.freeTri[(std::size_t)n] || stamp[(std::size_t)n] != 0) continue;
                if (sheetOf(*wm.mv, n, fr) != band.sheet) continue;
                probe = n;
                break;
            }
            if (probe >= 0) break;
        }
        if (probe >= 0 && circ) {
            const ProfPt pp = profOfPoint(triVert(*wm.mv, probe, 0), fr, 0);
            probeL1 = circleHatL1(band.pts, fit.Rmaj, fit.z0, fit.Rmin, pp);
        }
        const double extra = extrapOf(wm.q, probeL1);
        const bool loose = !circ || !std::isfinite(extra) || extra > wm.tau;
        std::vector<int> add;
        for (int t : band.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int n : wm.adj.nbr[(std::size_t)t]) {
                if (n < 0 || n >= nTri || band.in[(std::size_t)n]) continue;
                if (!wm.freeTri[(std::size_t)n] || stamp[(std::size_t)n] != 0) continue;
                if (sheetOf(*wm.mv, n, fr) != band.sheet) continue;
                if (!loose) {
                    bool ok = true;
                    for (int k = 0; k < 3 && ok; ++k) {
                        const ProfPt pp = profOfPoint(triVert(*wm.mv, n, k), fr, 0);
                        const double l1 = circleHatL1(band.pts, fit.Rmaj, fit.z0, fit.Rmin, pp);
                        if (circDist(fit, pp) > wm.tau + extrapOf(wm.q, l1)) ok = false;
                    }
                    if (!ok) continue;
                }
                add.push_back(n);
            }
        }
        std::sort(add.begin(), add.end());
        add.erase(std::unique(add.begin(), add.end()), add.end());
        if (add.empty()) break;
        for (int n : add) {
            band.in[(std::size_t)n] = 1;
            band.tris.push_back(n);
        }
    }
    collectBand(*wm.mv, fr, wm.g2l, band);
}

void endTangency(const WorkMesh& wm, const Frame& fr, const Band& band, const SurfaceFit& fit,
                 double& ang0, double& ang1) {
    ang0 = ang1 = 0;
    const std::vector<int> frnt = frontierOf(wm, band);
    if (frnt.empty() || band.pts.empty() || !fit.hatOk) return;
    double zMin = band.pts[0].z, zMax = band.pts[0].z;
    for (const ProfPt& p : band.pts) {
        zMin = std::min(zMin, p.z);
        zMax = std::max(zMax, p.z);
    }
    auto angleAt = [&](bool lo) {
        double best = 1e300;
        int pick = -1;
        for (int t : frnt) {
            for (int k = 0; k < 3; ++k) {
                const ProfPt p = profOfPoint(triVert(*wm.mv, t, k), fr, 0);
                const double dz = std::abs(p.z - (lo ? zMin : zMax));
                if (dz < best) {
                    best = dz;
                    pick = t;
                }
            }
        }
        if (pick < 0) return 0.0;
        std::vector<ProfPt> q;
        for (int k = 0; k < 3; ++k) q.push_back(profOfPoint(triVert(*wm.mv, pick, k), fr, 0));
        const SurfaceFit nf = fitSurface(q, wm.tau);
        if (!nf.hatOk && !fit.hatOk) return 0.0;
        const double b0 = fit.regressRho ? fit.b : ((std::abs(fit.b) > gp::Resolution()) ? 1.0 / fit.b : 0.0);
        double b1 = 0;
        if (nf.hatOk)
            b1 = nf.regressRho ? nf.b : ((std::abs(nf.b) > gp::Resolution()) ? 1.0 / nf.b : 0.0);
        const double d0 = std::hypot(b0, 1.0);
        const double d1 = std::hypot(b1, 1.0);
        double c = (b0 * b1 + 1.0) / (d0 * d1);
        c = std::max(-1.0, std::min(1.0, c));
        return std::acos(c) * 180.0 / kPi;
    };
    ang0 = angleAt(true);
    ang1 = angleAt(false);
}

int countStraddle(const SegmentWork& work, const std::vector<int>& tris, int nTri) {
    if (work.provisionals.empty()) return 0;
    std::vector<char> take((std::size_t)nTri, 0);
    for (int t : tris)
        if (t >= 0 && t < nTri) take[(std::size_t)t] = 1;
    int n = 0;
    for (const Provisional& p : work.provisionals) {
        if (p.claim != ProvClaim::Unclaimed) continue;
        int hit = 0, miss = 0;
        for (int t : p.tris) {
            if (t < 0 || t >= nTri) continue;
            if (take[(std::size_t)t]) ++hit;
            else ++miss;
        }
        if (hit > 0 && miss > 0) ++n;
    }
    return n;
}

void refitProv(const MeshView& mv, Provisional& p) {
    if (p.tris.empty()) {
        p.area = 0;
        p.maxVertexDev = 0;
        p.rmsVertexDev = 0;
        return;
    }
    std::sort(p.tris.begin(), p.tris.end());
    gp_Ax3 pl = p.plane;
    pcaPlane(mv, p.tris, pl);
    p.plane = pl;
    const gp_XYZ n = pl.Direction().XYZ();
    const gp_XYZ o = pl.Location().XYZ();
    double area = 0, maxD = 0, sum2 = 0;
    int nv = 0;
    std::vector<int> seen;
    for (int t : p.tris) {
        gp_XYZ nn;
        area += triAreaNormal(mv, t, nn);
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int g = mv.tris[gt][k];
            seen.push_back(g);
        }
    }
    std::sort(seen.begin(), seen.end());
    seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
    for (int g : seen) {
        const double d = std::abs((mv.pts[g] - o).Dot(n));
        maxD = std::max(maxD, d);
        sum2 += d * d;
        ++nv;
    }
    p.area = area;
    p.maxVertexDev = maxD;
    p.rmsVertexDev = nv ? std::sqrt(sum2 / nv) : 0;
}

void trimProvisionals(const MeshView& mv, SegmentWork& work, const std::vector<int>& taken,
                      int nTri, ProvClaim fullClaim, std::vector<int>& provTouched) {
    std::vector<char> take((std::size_t)std::max(nTri, 0), 0);
    for (int t : taken)
        if (t >= 0 && t < nTri) take[(std::size_t)t] = 1;
    const int nProv = (int)work.provisionals.size();
    for (int i = 0; i < nProv; ++i) {
        Provisional& p = work.provisionals[(std::size_t)i];
        if (p.claim != ProvClaim::Unclaimed) continue;
        std::vector<int> keep, hit;
        for (int t : p.tris) {
            if (t >= 0 && t < nTri && take[(std::size_t)t]) hit.push_back(t);
            else keep.push_back(t);
        }
        if (hit.empty()) continue;
        provTouched.push_back(i);
        if (keep.empty()) {
            p.claim = fullClaim;
            continue;
        }
        // Split the remainder into edge-connected components. The first keeps
        // this provisional; the rest are new Unclaimed provisionals.
        Adj adj = buildAdj(mv);
        std::vector<char> inKeep((std::size_t)nTri, 0);
        for (int t : keep)
            if (t >= 0 && t < nTri) inKeep[(std::size_t)t] = 1;
        std::vector<char> seen((std::size_t)nTri, 0);
        std::vector<std::vector<int>> comps;
        for (int t : keep) {
            if (t < 0 || t >= nTri || seen[(std::size_t)t]) continue;
            std::vector<int> stack{t};
            seen[(std::size_t)t] = 1;
            std::vector<int> comp;
            while (!stack.empty()) {
                const int c = stack.back();
                stack.pop_back();
                comp.push_back(c);
                for (int n : adj.nbr[(std::size_t)c]) {
                    if (n < 0 || n >= nTri || !inKeep[(std::size_t)n] || seen[(std::size_t)n])
                        continue;
                    seen[(std::size_t)n] = 1;
                    stack.push_back(n);
                }
            }
            std::sort(comp.begin(), comp.end());
            comps.push_back(std::move(comp));
        }
        std::sort(comps.begin(), comps.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
            return minTriOf(a) < minTriOf(b);
        });
        if (comps.empty()) {
            p.claim = fullClaim;
            continue;
        }
        p.tris = comps[0];
        refitProv(mv, p);
        for (std::size_t c = 1; c < comps.size(); ++c) {
            Provisional np = p;
            np.tris = comps[c];
            np.claim = ProvClaim::Unclaimed;
            refitProv(mv, np);
            work.provisionals.push_back(std::move(np));
        }
    }
}

bool emitCone(const WorkMesh& wm, SegmentWork& work, const Donor& donor, const Band& band,
              const SurfaceFit& fit, const ClosedM& cl, double& ang0, double& ang1, int& straddle) {
    const MeshView& mv = *wm.mv;
    double zMin = band.pts[0].z, zMax = band.pts[0].z;
    for (const ProfPt& p : band.pts) {
        zMin = std::min(zMin, p.z);
        zMax = std::max(zMax, p.z);
    }
    const double RzMin = fit.a + fit.b * zMin;
    const double RzMax = fit.a + fit.b * zMax;
    if (!(RzMin > wm.tau) || !(RzMax > wm.tau)) return false;
    const bool loAtMin = RzMin <= RzMax;
    const double Rlo = loAtMin ? RzMin : RzMax;
    const double Rhi = loAtMin ? RzMax : RzMin;
    const double zLo = loAtMin ? zMin : zMax;
    const double height = std::abs(zMax - zMin);
    const double hSigned = loAtMin ? height : -height;
    if (!(Rhi - Rlo > wm.tau)) return false;
    const gp_XYZ loc = donor.fr.origin + donor.fr.axis * zLo;
    gp_XYZ xDir = donor.fr.e1;
    Frame placed = makeFrame(loc, donor.fr.axis, xDir);
    double sigma = 0, area = 0;
    for (int t : band.tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
        area += a;
        gp_XYZ rad = triCentroid(mv, t) - loc;
        rad -= placed.axis * rad.Dot(placed.axis);
        const double rm = rad.Modulus();
        if (rm > gp::Resolution()) sigma += a * n.Dot(rad / rm);
    }
    const int nSides = std::max(cl.nSides, 3);
    double dVol = 0;
    for (int t : band.tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
        if (!(a > 0.0)) continue;
        const double u = (triCentroid(mv, t) - loc).Dot(placed.axis);
        const double rhoC = Rlo + (std::abs(hSigned) > gp::Resolution() ? (u / hSigned) : 0.0) * (Rhi - Rlo);
        dVol += a * coneFrustumDVolDensity(rhoC, Rlo, Rhi, hSigned, nSides);
    }
    if (!std::isfinite(dVol)) dVol = 0;
    double maxD = 0, sum2 = 0;
    for (const ProfPt& p : band.pts) {
        const double d = linePerp(fit, p);
        maxD = std::max(maxD, d);
        sum2 += d * d;
    }
    Region r;
    r.id = 0;
    r.type = SurfType::Cone;
    r.origin = Origin::ChamferCone;
    r.ax = gp_Ax3(gp_Pnt(loc.X(), loc.Y(), loc.Z()), gp_Dir(placed.axis.X(), placed.axis.Y(), placed.axis.Z()),
                  gp_Dir(placed.e1.X(), placed.e1.Y(), placed.e1.Z()));
    r.radius = Rlo;
    r.uMin = 0;
    r.uMax = kTwoPi;
    r.vMin = Rlo;
    r.vMax = Rhi;
    r.closed360 = true;
    r.outwardNormal = sigma > 0.0;
    r.tris = band.tris;
    r.maxVertexDev = maxD;
    r.rmsVertexDev = band.pts.empty() ? 0 : std::sqrt(sum2 / (double)band.pts.size());
    r.chordSagitta = chordSagitta(Rlo, nSides);
    r.nSides = nSides;
    r.dVolPredicted = (sigma > 0.0 ? 1.0 : -1.0) * dVol;
    r.maxVertexSnap = hSigned;
    r.reject = Reject::None;
    r.builtAs = BuiltAs::NotBuilt;
    r.filletNbrA = donor.idx;
    r.filletNbrB = -1;
    endTangency(wm, donor.fr, band, fit, ang0, ang1);
    straddle = countStraddle(work, band.tris, (int)mv.nTri);
    std::vector<int> touched;
    trimProvisionals(mv, work, band.tris, (int)mv.nTri, ProvClaim::ConsumedCylinder, touched);
    work.accepted.push_back(std::move(r));
    (void)area;
    return true;
}

void arcExtent(std::vector<double> a, double& lo, double& hi, double& maxGap) {
    lo = hi = maxGap = 0;
    if (a.empty()) return;
    std::sort(a.begin(), a.end());
    if (a.size() == 1) {
        lo = hi = a.front();
        maxGap = kTwoPi;
        return;
    }
    std::size_t ig = a.size() - 1;
    maxGap = (a.front() + kTwoPi) - a.back();
    for (std::size_t i = 0; i + 1 < a.size(); ++i) {
        const double g = a[i + 1] - a[i];
        if (g > maxGap) {
            maxGap = g;
            ig = i;
        }
    }
    lo = a[(ig + 1) % a.size()];
    hi = a[ig];
    if (hi < lo) hi += kTwoPi;
}

bool emitTorus(const WorkMesh& wm, SegmentWork& work, const Donor& donor, const Band& band,
               const TorusCert& cert, const ClosedM& cl, int& straddle) {
    const MeshView& mv = *wm.mv;
    const Frame& fr = donor.fr;
    int lvMin = std::numeric_limits<int>::max();
    gp_XYZ pMin(0, 0, 0);
    bool have = false;
    for (const Vtx& v : band.verts) {
        if (!have || v.lv < lvMin) {
            lvMin = v.lv;
            pMin = v.p;
            have = true;
        }
    }
    gp_XYZ xdir = pMin - fr.origin;
    xdir -= fr.axis * xdir.Dot(fr.axis);
    if (!(xdir.Modulus() > gp::Resolution())) xdir = fr.e1;
    else xdir /= xdir.Modulus();
    const gp_XYZ centre = fr.origin + fr.axis * cert.z0;
    std::vector<double> psis, vs;
    for (const ProfPt& p : band.pts) {
        psis.push_back(p.psi);
        double v = std::atan2(p.z - cert.z0, p.rho - cert.Rmaj);
        if (v < 0.0) v += kTwoPi;
        vs.push_back(v);
    }
    double uLo = 0, uHi = 0, uGap = 0, vLo = 0, vHi = 0, vGap = 0;
    arcExtent(psis, uLo, uHi, uGap);
    arcExtent(vs, vLo, vHi, vGap);
    while (vLo > kPi) vLo -= kTwoPi;
    while (vLo <= -kPi) vLo += kTwoPi;
    while (vHi > kPi) vHi -= kTwoPi;
    while (vHi <= -kPi) vHi += kTwoPi;
    if (vHi < vLo) vHi += kTwoPi;
    const bool closed = cl.closed;
    double sigma = 0;
    const int nSides = std::max(cert.cols, 1);
    const double gamma = kTwoPi / (double)nSides;
    const double sHalf = std::sin(0.5 * gamma);
    const double kSeg = (sHalf > gp::Resolution()) ? (gamma - std::sin(gamma)) / (4.0 * sHalf) : 0.0;
    double dVolSum = 0;
    for (int t : band.tris) {
        gp_XYZ n;
        const double a2 = triAreaNormal(mv, t, n);
        const gp_XYZ d = triCentroid(mv, t) - fr.origin;
        const double zc = d.Dot(fr.axis);
        gp_XYZ rad = d - fr.axis * zc;
        const double rc = rad.Modulus();
        if (!(rc > gp::Resolution())) continue;
        rad /= rc;
        gp_XYZ mhat = rad * (rc - cert.Rmaj) + fr.axis * (zc - cert.z0);
        const double mm = mhat.Modulus();
        if (mm > gp::Resolution()) sigma += a2 * n.Dot(mhat / mm);
        dVolSum += a2 * rc * kSeg * std::fabs(n.Dot(rad));
    }
    Region r;
    r.id = 0;
    r.type = SurfType::Torus;
    r.origin = Origin::TorusBlend;
    r.ax = gp_Ax3(gp_Pnt(centre.X(), centre.Y(), centre.Z()),
                  gp_Dir(fr.axis.X(), fr.axis.Y(), fr.axis.Z()),
                  gp_Dir(xdir.X(), xdir.Y(), xdir.Z()));
    r.radius = cert.Rmaj;
    r.radius2 = cert.Rmin;
    r.uMin = closed ? 0.0 : uLo;
    r.uMax = closed ? kTwoPi : uHi;
    r.vMin = vLo;
    r.vMax = vHi;
    r.closed360 = closed;
    r.outwardNormal = sigma >= 0.0;
    r.tris = band.tris;
    r.maxVertexDev = cert.maxDev;
    r.rmsVertexDev = cert.rmsDev;
    r.chordSagitta = chordSagitta(cert.Rmaj, nSides);
    r.nSides = nSides;
    r.dVolPredicted = (sigma >= 0.0 ? 1.0 : -1.0) * dVolSum;
    r.reject = Reject::None;
    r.builtAs = BuiltAs::NotBuilt;
    r.filletNbrA = donor.idx;
    r.filletNbrB = -1;
    straddle = countStraddle(work, band.tris, (int)mv.nTri);
    std::vector<int> touched;
    trimProvisionals(mv, work, band.tris, (int)mv.nTri, ProvClaim::PendingTorus, touched);
    TorusPendingReclaim pend;
    pend.provIdx = touched;
    const int packIdx = (int)work.torusPending.size();
    work.torusPending.push_back(std::move(pend));
    r.torusPending = packIdx;
    work.torusAdmitted = true;
    (void)uGap;
    (void)vGap;
    work.accepted.push_back(std::move(r));
    return true;
}

void updateDonor(SegmentWork& work, Donor& donor) {
    if (donor.updated || donor.idx < 0 || donor.idx >= (int)work.accepted.size()) return;
    Region& r = work.accepted[(std::size_t)donor.idx];
    r.ax = gp_Ax3(gp_Pnt(donor.fr.origin.X(), donor.fr.origin.Y(), donor.fr.origin.Z()),
                  gp_Dir(donor.fr.axis.X(), donor.fr.axis.Y(), donor.fr.axis.Z()),
                  gp_Dir(donor.fr.e1.X(), donor.fr.e1.Y(), donor.fr.e1.Z()));
    if (r.type == SurfType::Cylinder) r.radius = donor.radius;
    r.maxVertexDev = donor.maxDev;
    r.rmsVertexDev = donor.rms;
    donor.updated = true;
}

// The cone face reuses a neighbour cylinder's rim circle. That circle's
// parameter is the cylinder's XDirection, and the cone pcurve assumes the
// two frames share it. Rotate every coaxial cylinder that touches the band
// onto the donor frame's XDirection. Location, radius, and axis stay.
void shareSeamFrame(const WorkMesh& wm, SegmentWork& work, const Donor& donor, const Band& band) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> touch((std::size_t)std::max(nTri, 0), 0);
    for (int t : band.tris) {
        if (t < 0 || t >= nTri) continue;
        for (int n : wm.adj.nbr[(std::size_t)t])
            if (n >= 0 && n < nTri) touch[(std::size_t)n] = 1;
    }
    const gp_Dir xdir(donor.fr.e1.X(), donor.fr.e1.Y(), donor.fr.e1.Z());
    for (Region& r : work.accepted) {
        if (r.type != SurfType::Cylinder) continue;
        bool hit = false;
        for (int t : r.tris) {
            if (t >= 0 && t < nTri && touch[(std::size_t)t]) {
                hit = true;
                break;
            }
        }
        if (!hit) continue;
        Frame rf = makeFrame(r.ax.Location().XYZ(), r.ax.Direction().XYZ(), r.ax.XDirection().XYZ());
        if (!sameAxis(donor.fr, rf, wm.tau, std::max(donor.length, wm.tau))) continue;
        try {
            r.ax = gp_Ax3(r.ax.Location(), r.ax.Direction(), xdir);
        } catch (const Standard_Failure&) {
        }
    }
}

bool adjacentSameSurface(const WorkMesh& wm, SegmentWork& work, const Band& band,
                         const SurfaceFit& fit, bool circle) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<int> accHit(work.accepted.size(), 0);
    for (int t : band.tris) {
        if (t < 0 || t >= nTri) continue;
        for (int n : wm.adj.nbr[(std::size_t)t]) {
            if (n < 0 || n >= nTri || band.in[(std::size_t)n]) continue;
            for (int i = 0; i < (int)work.accepted.size(); ++i) {
                const Region& r = work.accepted[(std::size_t)i];
                if (std::binary_search(r.tris.begin(), r.tris.end(), n)) accHit[(std::size_t)i] = 1;
            }
        }
    }
    for (int i = 0; i < (int)work.accepted.size(); ++i) {
        if (!accHit[(std::size_t)i]) continue;
        const Region& r = work.accepted[(std::size_t)i];
        if (r.tris.empty()) continue;
        bool all = true;
        std::vector<int> gids;
        for (int t : r.tris) {
            for (int k = 0; k < 3; ++k) {
                const int gt = wm.mv->compTris ? wm.mv->compTris[t] : t;
                gids.push_back(wm.mv->tris[gt][k]);
            }
        }
        std::sort(gids.begin(), gids.end());
        gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
        for (int g : gids) {
            const ProfPt p = profOfPoint(wm.mv->pts[g], /*frame*/ Frame{}, 0);
            (void)p;
        }
        // Distance uses the band frame stored on the fit via profile coords.
        // Recompute below with the caller frame — this stub is replaced.
        (void)all;
        (void)fit;
        (void)circle;
        (void)r;
    }
    return false;
}

bool sameSurfaceFr(const WorkMesh& wm, const SegmentWork& work, const Frame& fr, const Band& band,
                   const SurfaceFit& fit, bool asCircle) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> bandTri((std::size_t)nTri, 0);
    for (int t : band.tris)
        if (t >= 0 && t < nTri) bandTri[(std::size_t)t] = 1;
    for (const Region& r : work.accepted) {
        if (r.tris.empty()) continue;
        bool touch = false;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int n : wm.adj.nbr[(std::size_t)t]) {
                if (n >= 0 && n < nTri && bandTri[(std::size_t)n]) touch = true;
            }
            if (touch) break;
        }
        if (!touch) continue;
        std::vector<int> gids;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int k = 0; k < 3; ++k) {
                const int gt = wm.mv->compTris ? wm.mv->compTris[t] : t;
                gids.push_back(wm.mv->tris[gt][k]);
            }
        }
        std::sort(gids.begin(), gids.end());
        gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
        if (gids.empty()) continue;
        bool all = true;
        for (int g : gids) {
            const ProfPt p = profOfPoint(wm.mv->pts[g], fr, 0);
            const double d = asCircle ? circDist(fit, p) : linePerp(fit, p);
            if (!(d <= wm.tau)) {
                all = false;
                break;
            }
        }
        if (all) return true;
    }
    return false;
}

DiagRec baseDiag(const char* mode, const Donor& donor, const Band& band, double tau) {
    DiagRec d;
    d.mode = mode;
    d.seed = band.seed;
    d.nTri = (int)band.tris.size();
    d.nV = (int)band.pts.size();
    d.fit = "anchored";
    d.axis = donor.fr.axis;
    d.maxDevTau = (tau > 0) ? d.maxDev / tau : 0;
    std::vector<double> zs, ps, rs;
    for (const ProfPt& p : band.pts) {
        zs.push_back(p.z);
        ps.push_back(p.psi);
        rs.push_back(p.rho);
    }
    d.levels = countLevels(zs, tau);
    d.stations = countColumns(ps, rs, tau);
    return d;
}

const char* torusWhy(const TorusCert& c) {
    if (c.ok) return "accept";
    if (!c.why) return "underdetermined";
    if (std::strcmp(c.why, "profile-underdetermined") == 0) return "underdetermined";
    if (std::strcmp(c.why, "fit-singular") == 0) return "underdetermined";
    if (std::strcmp(c.why, "cone-or-cylinder-or-plane") == 0) return "line";
    if (std::strcmp(c.why, "sphere") == 0) return "sphere";
    if (std::strcmp(c.why, "degenerate") == 0) return "degenerate";
    if (std::strcmp(c.why, "spindle-torus") == 0) return "spindle";
    if (std::strcmp(c.why, "vertex-off-surface") == 0) return "vertex-off";
    return "vertex-off";
}

}  // namespace

bool claimRelRevolve(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work, RelRevolveMode mode) {
    try {
        const double tau = tauOf(mv);
        const char* modeName = (mode == RelRevolveMode::Torus) ? "torus" : "line";
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) {
            printSum(modeName, 0, 0, 0);
            return true;
        }
        const int nTri = (int)mv.nTri;
        const int nVtx = (int)mv.nVtx;
        WorkMesh wm;
        wm.mv = &mv;
        wm.adj = buildAdj(mv);
        wm.tau = tau;
        wm.q = mv.quantFloor;
        wm.thetaBin = tol.thetaBin;
        wm.maxIter = p.maxRefineIters > 0 ? p.maxRefineIters : 1;
        wm.g2l.assign((std::size_t)std::max(nVtx, 1), -1);
        for (int v = 0; v < nVtx; ++v) {
            const int g = mv.compVtx ? mv.compVtx[v] : v;
            if (g >= (int)wm.g2l.size()) wm.g2l.resize((std::size_t)g + 1, -1);
            if (g >= 0) wm.g2l[(std::size_t)g] = v;
        }
        wm.freeTri.assign((std::size_t)nTri, 0);
        wm.triProv.assign((std::size_t)nTri, -1);
        for (int i = 0; i < (int)work.provisionals.size(); ++i) {
            Provisional& pr = work.provisionals[(std::size_t)i];
            if (pr.claim != ProvClaim::Unclaimed) continue;
            for (int t : pr.tris) {
                if (t < 0 || t >= nTri) continue;
                wm.freeTri[(std::size_t)t] = 1;
                wm.triProv[(std::size_t)t] = i;
            }
        }
        for (const Region& r : work.accepted) {
            for (int t : r.tris)
                if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
        }

        std::vector<Donor> donors;
        int refusedDonors = 0;
        for (int i = 0; i < (int)work.accepted.size(); ++i) {
            Region& r = work.accepted[(std::size_t)i];
            if (!regionIsDonorClass(r) || r.tris.empty()) continue;
            std::vector<int> tris = r.tris;
            std::sort(tris.begin(), tris.end());
            gp_Dir axisSeed = r.ax.Direction();
            gp_Dir axisOut = axisSeed;
            gp_Pnt center;
            double radius = 0, maxResid = 0;
            const bool fitOk = cylinderFitLS(mv, tris, axisSeed, axisOut, center, radius, maxResid,
                                             wm.maxIter);
            auto failDonor = [&]() {
                ++refusedDonors;
                if (!diagOn()) return;
                DiagRec d;
                d.mode = modeName;
                d.seed = minTriOf(r.tris);
                d.nTri = (int)r.tris.size();
                d.fit = "anchored";
                d.axis = r.ax.Direction().XYZ();
                d.r0 = r.radius;
                d.why = "no-donor";
                d.params = kParamsLine;
                printBand(d);
            };
            if (!fitOk || !(radius > tau)) {
                failDonor();
                continue;
            }
            Frame fr = makeFrame(center.XYZ(), axisOut.XYZ(), r.ax.XDirection().XYZ());
            double maxD = 0, sum2 = 0;
            int nv = 0;
            double zMin = 1e300, zMax = -1e300;
            std::vector<double> zs, ps, rs;
            std::vector<int> gids;
            for (int t : tris) {
                for (int k = 0; k < 3; ++k) {
                    const int gt = mv.compTris ? mv.compTris[t] : t;
                    gids.push_back(mv.tris[gt][k]);
                }
            }
            std::sort(gids.begin(), gids.end());
            gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
            for (int g : gids) {
                const ProfPt pp = profOfPoint(mv.pts[g], fr, 0);
                const double d = std::fabs(pp.rho - radius);
                maxD = std::max(maxD, d);
                sum2 += d * d;
                ++nv;
                zMin = std::min(zMin, pp.z);
                zMax = std::max(zMax, pp.z);
                zs.push_back(pp.z);
                ps.push_back(pp.psi);
                rs.push_back(pp.rho);
            }
            const int cols = countColumns(ps, rs, tau);
            const int levels = countLevels(zs, tau);
            const bool partialOk = cols >= 4 && levels >= 2;
            if (!(maxD <= tau) || !(r.closed360 || partialOk)) {
                failDonor();
                continue;
            }
            Donor d;
            d.idx = i;
            d.fr = fr;
            d.radius = radius;
            d.maxDev = maxD;
            d.rms = nv ? std::sqrt(sum2 / nv) : 0;
            d.length = (zMax > zMin) ? (zMax - zMin) : tau;
            d.minTri = minTriOf(r.tris);
            donors.push_back(d);
        }
        std::sort(donors.begin(), donors.end(), [](const Donor& a, const Donor& b) {
            if (a.minTri != b.minTri) return a.minTri < b.minTri;
            return a.idx < b.idx;
        });

        int claimed = 0, refused = 0;
        std::vector<int> stamp(nTri, 0);
        std::vector<char> hold((std::size_t)std::max(nTri, 0), 0);
        int stampGen = 0;

        for (Donor& donor : donors) {
            // A failed band on one donor is not a refusal on another axis.
            // Held triangles are rows of a blend already refused as cones.
            for (int t = 0; t < nTri; ++t)
                stamp[(std::size_t)t] = hold[(std::size_t)t] ? 1 : 0;
            stampGen = 1;
            std::vector<ConeCand> cands;
            // The candidate set is every free triangle reachable from the
            // donor without crossing an accepted triangle. A cone that meets
            // the donor only through an unclaimed blend is in this set; growth
            // at growTol then keeps the blend and the cone as separate bands.
            std::vector<char> pool((std::size_t)nTri, 0);
            std::vector<int> stack;
            auto seedBoundary = [&](const std::vector<int>& tris) {
                for (int t : tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int n : wm.adj.nbr[(std::size_t)t]) {
                        if (n < 0 || n >= nTri || !wm.freeTri[(std::size_t)n] || pool[(std::size_t)n])
                            continue;
                        pool[(std::size_t)n] = 1;
                        stack.push_back(n);
                    }
                }
            };
            if (donor.idx >= 0 && donor.idx < (int)work.accepted.size())
                seedBoundary(work.accepted[(std::size_t)donor.idx].tris);
            for (const Region& r : work.accepted) {
                if (r.origin != Origin::ChamferCone && r.origin != Origin::TorusBlend) continue;
                Frame rf = makeFrame(r.ax.Location().XYZ(), r.ax.Direction().XYZ(),
                                     r.ax.XDirection().XYZ());
                if (!sameAxis(donor.fr, rf, tau, std::max(donor.length, tau))) continue;
                seedBoundary(r.tris);
            }
            while (!stack.empty()) {
                const int t = stack.back();
                stack.pop_back();
                for (int n : wm.adj.nbr[(std::size_t)t]) {
                    if (n < 0 || n >= nTri || !wm.freeTri[(std::size_t)n] || pool[(std::size_t)n])
                        continue;
                    pool[(std::size_t)n] = 1;
                    stack.push_back(n);
                }
            }
            int scan = 0;
            while (true) {
                int seed = -1;
                for (int t = scan; t < nTri; ++t) {
                    if (!pool[(std::size_t)t] || !wm.freeTri[(std::size_t)t] || stamp[(std::size_t)t] != 0)
                        continue;
                    if (sheetOf(mv, t, donor.fr) == 0) continue;
                    seed = t;
                    scan = t + 1;
                    break;
                }
                if (seed < 0) break;
                const int sh = sheetOf(mv, seed, donor.fr);
                if (sh == 0) continue;

                Band band;
                band.in.assign((std::size_t)nTri, 0);
                band.sheet = sh;
                band.seed = seed;
                band.tris.push_back(seed);
                band.in[(std::size_t)seed] = 1;
                collectBand(mv, donor.fr, wm.g2l, band);

                const bool wantTorus = mode == RelRevolveMode::Torus;
                if (!wantTorus) growLineDrain(wm, donor.fr, band, stamp);
                else growTorusDrain(wm, donor.fr, band, stamp);

                collectBand(mv, donor.fr, wm.g2l, band);
                if (band.tris.empty()) continue;
                // Subset of an already-refused band: do not accept.
                const int g0 = stamp[(std::size_t)band.tris.front()];
                bool subset = g0 != 0;
                if (subset) {
                    for (int t : band.tris)
                        if (stamp[(std::size_t)t] != g0) subset = false;
                }
                if (subset) continue;

                const SurfaceFit fit = fitSurface(band.pts, tau);
                const ClosedM cl = closedOf(mv, donor.fr, band.tris, band.verts, wm.thetaBin);
                std::vector<double> zs, psis, rhos;
                for (const ProfPt& pp : band.pts) {
                    zs.push_back(pp.z);
                    psis.push_back(pp.psi);
                    rhos.push_back(pp.rho);
                }
                const int levels = countLevels(zs, tau);
                const int cols = countColumns(psis, rhos, tau);
                const int stations = countProfPts(band.pts, tau);

                auto stampFail = [&](const char* why) {
                    ++refused;
                    const int g = ++stampGen;
                    for (int t : band.tris)
                        if (t >= 0 && t < nTri && stamp[(std::size_t)t] == 0) stamp[(std::size_t)t] = g;
                    DiagRec d = baseDiag(modeName, donor, band, tau);
                    d.stations = stations;
                    d.levels = levels;
                    d.params = wantTorus ? kParamsCircle : kParamsLine;
                    d.maxDev = wantTorus ? fit.circMax : fit.orthMax;
                    d.maxDevTau = tau > 0 ? d.maxDev / tau : 0;
                    d.lineResTau = tau > 0 ? fit.orthMax / tau : 0;
                    d.span = cl.span;
                    d.spanTol = (cl.nBands > 0) ? 1.5 * kTwoPi / cl.nBands : 0;
                    d.halfAng = std::atan(std::abs(fit.b)) * 180.0 / kPi;
                    d.r0 = fit.regressRho ? fit.a : 0;
                    d.why = why;
                    if (fit.circle) {
                        d.r0 = fit.Rmaj;
                        d.r1 = fit.Rmin;
                    }
                    printBand(d);
                };

                if (!wantTorus) {
                    if (stations < 2 || !fit.hatOk) {
                        stampFail("underdetermined");
                        continue;
                    }
                    if (!(fit.orthMax <= tau)) {
                        stampFail("vertex-off");
                        continue;
                    }
                    double zMin = band.pts[0].z, zMax = zMin, rMin = band.pts[0].rho, rMax = rMin;
                    for (const ProfPt& pp : band.pts) {
                        zMin = std::min(zMin, pp.z);
                        zMax = std::max(zMax, pp.z);
                        rMin = std::min(rMin, pp.rho);
                        rMax = std::max(rMax, pp.rho);
                    }
                    const double zSpan = zMax - zMin;
                    const double rSpan = rMax - rMin;
                    if (zSpan <= tau && rSpan <= tau) {
                        stampFail("underdetermined");
                        continue;
                    }
                    if (zSpan <= tau) {
                        stampFail("line");
                        continue;
                    }
                    if (rSpan <= tau) {
                        stampFail("same-surface");
                        continue;
                    }
                    // Two stations are a cone only when each is a τ-ring.
                    // Three or more still have to be concyclic stations.
                    if (stations < minStations(kParamsLine) && stations < 2) {
                        stampFail("underdetermined");
                        continue;
                    }
                    if (stations == 2 || stations < minStations(kParamsLine)) {
                        if (stations < 2 || !stationsConcyclic(band.pts, tau, fit)) {
                            stampFail("underdetermined");
                            continue;
                        }
                    }
                    if (!stationsConcyclic(band.pts, tau, fit)) {
                        stampFail("vertex-off");
                        continue;
                    }
                    if (lineIsCircleSubset(wm, donor.fr, band, tau)) {
                        stampFail("line");
                        continue;
                    }
                    if (sameSurfaceFr(wm, work, donor.fr, band, fit, false)) {
                        stampFail("same-surface");
                        continue;
                    }
                    double maxP = 0;
                    for (const ProfPt& pp : band.pts) maxP = std::max(maxP, linePerp(fit, pp));
                    if (maxP > tau) {
                        stampFail("vertex-off");
                        continue;
                    }
                    if (!cl.closed) {
                        stampFail("partial-cone");
                        continue;
                    }
                    ConeCand cand;
                    cand.band = band;
                    cand.fit = fit;
                    cand.cl = cl;
                    cand.maxP = maxP;
                    cands.push_back(std::move(cand));
                    const int g = ++stampGen;
                    for (int t : band.tris)
                        if (t >= 0 && t < nTri) {
                            stamp[(std::size_t)t] = g;
                            hold[(std::size_t)t] = 1;
                        }
                    continue;
                }

                // Torus mode. Line class was claimed, or refused, in Line mode.
                if (fit.line && fit.orthMax <= tau) {
                    stampFail("line");
                    continue;
                }
                const TorusCert cert = certifyTorus(band.pts, tau);
                if (!cert.ok) {
                    stampFail(torusWhy(cert));
                    continue;
                }
                if (torusIsFailingSubset(wm, donor.fr, band, tau)) {
                    stampFail("vertex-off");
                    continue;
                }
                if (sameSurfaceFr(wm, work, donor.fr, band, fit, true)) {
                    stampFail("same-surface");
                    continue;
                }
                int straddle = countStraddle(work, band.tris, nTri);
                if (!emitTorus(wm, work, donor, band, cert, cl, straddle)) {
                    stampFail("degenerate");
                    continue;
                }
                // emitTorus trims before its own straddle recount; keep the
                // pre-trim count captured above.
                updateDonor(work, donor);
                for (int t : band.tris)
                    if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
                ++claimed;
                DiagRec d = baseDiag(modeName, donor, band, tau);
                d.stations = cert.profPts;
                d.levels = cert.levels;
                d.params = kParamsCircle;
                d.maxDev = cert.maxDev;
                d.maxDevTau = cert.maxDev / tau;
                d.lineResTau = cert.lineRes / tau;
                d.span = cl.span;
                d.spanTol = (cl.nBands > 0) ? 1.5 * kTwoPi / cl.nBands : 0;
                d.r0 = cert.Rmaj;
                d.r1 = cert.Rmin;
                d.halfAng = 0;
                double ang0 = 0, ang1 = 0;
                endTangency(wm, donor.fr, band, fit, ang0, ang1);
                d.tan0 = ang0;
                d.tan1 = ang1;
                d.straddle = straddle;
                d.why = "accept";
                d.nV = (int)band.pts.size();
                printBand(d);
            }
            if (!cands.empty()) {
                std::vector<char> drop;
                markTorusChords(wm, donor.fr, cands, tau, drop);
                for (int i = 0; i < (int)cands.size(); ++i) {
                    ConeCand& cand = cands[(std::size_t)i];
                    if (!drop.empty() && drop[(std::size_t)i]) {
                        ++refused;
                        DiagRec d = baseDiag(modeName, donor, cand.band, tau);
                        d.stations = countProfPts(cand.band.pts, tau);
                        d.levels = d.stations;
                        d.params = kParamsLine;
                        d.maxDev = cand.maxP;
                        d.maxDevTau = tau > 0 ? cand.maxP / tau : 0;
                        d.lineResTau = tau > 0 ? cand.fit.orthMax / tau : 0;
                        d.span = cand.cl.span;
                        d.halfAng = std::atan(std::abs(cand.fit.b)) * 180.0 / kPi;
                        d.why = "line";
                        printBand(d);
                        continue;
                    }
                    // The cone pcurve treats the rim circle's first parameter as
                    // cone-u. A vertex a hair on the negative side of X is stored
                    // as 2π − ε and the pcurve doubles it. Put X through a band
                    // vertex so that parameter is 0.
                    {
                        double best = 1e300;
                        gp_XYZ bestR(0, 0, 0);
                        bool have = false;
                        for (const Vtx& v : cand.band.verts) {
                            gp_XYZ rad = v.p - donor.fr.origin;
                            rad -= donor.fr.axis * rad.Dot(donor.fr.axis);
                            const double m = rad.Modulus();
                            if (!(m > gp::Resolution())) continue;
                            rad /= m;
                            const double miss = 1.0 - rad.Dot(donor.fr.e1);
                            if (!have || miss < best) {
                                best = miss;
                                bestR = rad;
                                have = true;
                            }
                        }
                        if (have) donor.fr.e1 = bestR;
                    }
                    double ang0 = 0, ang1 = 0;
                    int straddle = 0;
                    if (!emitCone(wm, work, donor, cand.band, cand.fit, cand.cl, ang0, ang1,
                                  straddle)) {
                        ++refused;
                        continue;
                    }
                    shareSeamFrame(wm, work, donor, cand.band);
                    updateDonor(work, donor);
                    for (int t : cand.band.tris)
                        if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
                    ++claimed;
                    DiagRec d = baseDiag(modeName, donor, cand.band, tau);
                    d.stations = countProfPts(cand.band.pts, tau);
                    d.levels = d.stations;
                    d.params = kParamsLine;
                    d.maxDev = cand.maxP;
                    d.maxDevTau = tau > 0 ? cand.maxP / tau : 0;
                    d.lineResTau = tau > 0 ? cand.fit.orthMax / tau : 0;
                    d.span = cand.cl.span;
                    d.spanTol = (cand.cl.nBands > 0) ? 1.5 * kTwoPi / cand.cl.nBands : 0;
                    d.halfAng = std::atan(std::abs(cand.fit.b)) * 180.0 / kPi;
                    d.dihedral = d.halfAng;
                    d.r0 = work.accepted.back().radius;
                    d.r1 = work.accepted.back().vMax;
                    d.tan0 = ang0;
                    d.tan1 = ang1;
                    d.straddle = straddle;
                    d.why = "accept";
                    printBand(d);
                }
            }
        }
        printSum(modeName, (int)donors.size(), claimed, refused + refusedDonors);
        (void)axisGap;
        return true;
    } catch (...) {
        return true;
    }
}

}  // namespace refit
}  // namespace stl2step
