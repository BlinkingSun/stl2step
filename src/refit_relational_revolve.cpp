// Relational revolve (D-fc-1.3 … 1.7, D-fc-1.10, D-fc-3.5). Line mode gathers
// rings about a certified donor axis and claims closed cones. A 2-station
// chain whose next rings are a profile circle is why=profile-row. Torus mode
// gathers by axis too (D-fc-4.6): a chain of concyclic rings, atoms joining
// whole. Axis comes from a donor refit and is not re-fit on the band.
// Certificate is τ = tauOf(mv). Never throws.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_profile_cert.hpp"
#include "refit_cone_math.hpp"
#include "refit_relational_arc.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <tuple>
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

// Perpendicular profile distance from the regression line, so growth and the
// τ certificate speak about the same surface distance.
double linePerp(const SurfaceFit& f, const ProfPt& p) {
    const double miss = radialMiss(f, p);
    if (!std::isfinite(miss)) return miss;
    const double scale = std::sqrt(1.0 + f.b * f.b);
    return miss / scale;
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
    for (int t : b.tris) {
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int g = mv.tris[gt][k];
            int lv = (g >= 0 && g < (int)g2l.size()) ? g2l[g] : -1;
            if (lv < 0) lv = g;
            Vtx v;
            v.lv = lv;
            v.tri = t;
            v.p = mv.pts[g];
            v.pr = profOfPoint(v.p, f, lv);
            b.verts.push_back(v);
        }
    }
    // One point per vertex, ascending lv, from its first occurrence in verts.
    std::vector<std::pair<int, std::size_t>> firstOf;
    firstOf.reserve(b.verts.size());
    for (std::size_t i = 0; i < b.verts.size(); ++i) firstOf.emplace_back(b.verts[i].lv, i);
    std::sort(firstOf.begin(), firstOf.end());
    for (std::size_t i = 0; i < firstOf.size(); ++i) {
        if (i > 0 && firstOf[i].first == firstOf[i - 1].first) continue;
        b.pts.push_back(b.verts[firstOf[i].second].pr);
    }
    std::sort(b.tris.begin(), b.tris.end());
    b.tris.erase(std::unique(b.tris.begin(), b.tris.end()), b.tris.end());
}

// The points of b.pts in their order: one per vertex, ascending lv, from its
// first occurrence in b.verts (collectBand's rule).
void bandPoints(const Band& b, std::vector<gp_XYZ>& out) {
    std::vector<std::pair<int, std::size_t>> firstOf;
    firstOf.reserve(b.verts.size());
    for (std::size_t i = 0; i < b.verts.size(); ++i) firstOf.emplace_back(b.verts[i].lv, i);
    std::sort(firstOf.begin(), firstOf.end());
    out.clear();
    for (std::size_t i = 0; i < firstOf.size(); ++i) {
        if (i > 0 && firstOf[i].first == firstOf[i - 1].first) continue;
        out.push_back(b.verts[firstOf[i].second].p);
    }
}

// certifyTorus's own count of distinct profile points (D-fc-1.6.2).
constexpr int kTorusProfPts = 5;

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
    std::vector<int> gids;     // its own certified vertices (global ids)
    std::vector<int> members;  // D-fc-4.8 (b): accepted indices of its coaxial family
    std::vector<std::vector<gp_XYZ>> famPts;  // D-fc-4.9: each member's certified vertices
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
    const char* side0 = "na";
    const char* side1 = "na";
    const char* why = "underdetermined";
};

void printBand(const DiagRec& d) {
    if (!diagOn()) return;
    std::fprintf(stderr,
                 "DIAG_REL_BAND stage=revolve mode=%s seed=%d nTri=%d nV=%d stations=%d "
                 "levels=%d params=%d fit=%s radii=%.8g/%.8g halfAng=%.6f "
                 "axis=%.6g,%.6g,%.6g maxDev=%.6g maxDevTau=%.4f lineResTau=%.4f "
                 "span=%.6g dihedral=%.6g spanTol=%.6g tan=%.4f,%.4f straddle=%d "
                 "side0=%s side1=%s why=%s\n",
                 d.mode, d.seed, d.nTri, d.nV, d.stations, d.levels, d.params, d.fit, d.r0,
                 d.r1, d.halfAng, d.axis.X(), d.axis.Y(), d.axis.Z(), d.maxDev, d.maxDevTau,
                 d.lineResTau, d.span, d.dihedral, d.spanTol, d.tan0, d.tan1, d.straddle,
                 d.side0, d.side1, d.why);
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
    // U-3 per torus (D-fc-4.11 (a)): triangles of torus faces the builder could
    // not ship (SegmentParams::torusWithheldTris), and the atoms a withheld band
    // would have absorbed (they stay accepted, in their pre-T state).
    std::vector<char> withheldTri;
    std::vector<char> atomHeld;
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

// D-fc-5.5. A partial cone (a ChamferCone that does not close) about the
// family axis (origin, axis): Location on the R_lo rim, XDirection the radial
// direction of the band's first azimuth (closedOf's uMin), so the admitted arc
// is u in [0, span] and the face never crosses the surface's period. Radii from
// the certified profile line; nSides by the partial column rule.
bool emitPartialCone(const WorkMesh& wm, SegmentWork& work, const Donor& donor, const Band& band,
                     const gp_XYZ& origin, const gp_XYZ& axis, const PartialConeCert& pc,
                     double& spanOut, int& straddle) {
    const MeshView& mv = *wm.mv;
    const Frame f0 = makeFrame(origin, axis, donor.fr.e1);
    std::vector<gp_XYZ> pts;
    bandPoints(band, pts);
    std::vector<Vtx> uniq;
    for (const gp_XYZ& p : pts) {
        Vtx v;
        v.p = p;
        v.pr = profOfPoint(p, f0, 0);
        uniq.push_back(v);
    }
    const ClosedM cl = closedOf(mv, f0, band.tris, uniq, wm.thetaBin);
    if (cl.closed || !(cl.span > 0.0)) return false;
    const double RzMin = pc.rhoAt(pc.zLo);
    const double RzMax = pc.rhoAt(pc.zHi);
    // The face must not pass its apex (vertices on the axis are never members).
    if (!(RzMin > 0.0) || !(RzMax > 0.0)) return false;
    const bool loAtMin = RzMin <= RzMax;
    const double Rlo = loAtMin ? RzMin : RzMax;
    const double Rhi = loAtMin ? RzMax : RzMin;
    const double zLo = loAtMin ? pc.zLo : pc.zHi;
    const double height = std::abs(pc.zHi - pc.zLo);
    const double hSigned = loAtMin ? height : -height;
    if (!(Rhi - Rlo > wm.tau)) return false;
    const gp_XYZ loc = f0.origin + f0.axis * zLo;
    const gp_XYZ xDir = f0.e1 * std::cos(cl.uMin) + f0.e2 * std::sin(cl.uMin);
    const Frame placed = makeFrame(loc, f0.axis, xDir);
    double sigma = 0;
    for (int t : band.tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
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
    const double kScale = std::sqrt(1.0 + pc.k * pc.k);
    double sum2 = 0;
    for (const Vtx& v : uniq) {
        const ProfPt pp = profOfPoint(v.p, f0, 0);
        const double d = std::fabs(pp.rho - pc.rhoAt(pp.z)) / kScale;
        sum2 += d * d;
    }
    Region r;
    r.id = 0;
    r.type = SurfType::Cone;
    r.origin = Origin::ChamferCone;
    r.ax = gp_Ax3(gp_Pnt(loc.X(), loc.Y(), loc.Z()),
                  gp_Dir(placed.axis.X(), placed.axis.Y(), placed.axis.Z()),
                  gp_Dir(placed.e1.X(), placed.e1.Y(), placed.e1.Z()));
    r.radius = Rlo;
    r.uMin = 0;
    r.uMax = cl.span;
    r.vMin = Rlo;
    r.vMax = Rhi;
    r.closed360 = false;
    r.outwardNormal = sigma > 0.0;
    r.tris = band.tris;
    r.maxVertexDev = pc.lineRes;
    r.rmsVertexDev = uniq.empty() ? 0 : std::sqrt(sum2 / (double)uniq.size());
    r.chordSagitta = chordSagitta(Rlo, nSides);
    r.nSides = nSides;
    r.dVolPredicted = (sigma > 0.0 ? 1.0 : -1.0) * dVol;
    r.maxVertexSnap = hSigned;
    r.reject = Reject::None;
    r.builtAs = BuiltAs::NotBuilt;
    r.filletNbrA = donor.idx;
    r.filletNbrB = -1;
    straddle = countStraddle(work, band.tris, (int)mv.nTri);
    std::vector<int> touched;
    trimProvisionals(mv, work, band.tris, (int)mv.nTri, ProvClaim::ConsumedCylinder, touched);
    work.accepted.push_back(std::move(r));
    spanOut = cl.span;
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
    // D-fc-4.6 F8: one seam azimuth per revolution axis. The torus takes the
    // donor frame's XDirection, the frame its band azimuths are measured in and
    // the one shareSeamFrame gives the coaxial cones and cylinders.
    const gp_XYZ xdir = fr.e1;
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
    // The side count takes the closure the torus ships with: closed by its
    // bands, or a full ring by the F4 rail rule reclaimTorusMode applies.
    bool fullRing = closed;
    if (!fullRing) {
        Region probe;
        probe.ax = gp_Ax3(gp_Pnt(centre.X(), centre.Y(), centre.Z()),
                          gp_Dir(fr.axis.X(), fr.axis.Y(), fr.axis.Z()),
                          gp_Dir(xdir.X(), xdir.Y(), xdir.Z()));
        probe.tris = band.tris;
        fullRing = torusRails(mv, probe).full;
    }
    const int nSides = torusSidesOf(cert.cols, uHi - uLo, fullRing);
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

// One scan per axis. Later donors share it only when the longer axis length
// makes tanAngTol strict, so a short tilted sliver cannot lock out a real axis.
struct AxisScan {
    Frame fr;
    double length = 0;
    bool scanned = false;
};

struct SideFit {
    bool have = false;
    double circDev = 0;
    double lineRes = 0;
};

void fmtSide(char* buf, int n, const SideFit& s, double tau) {
    if (!s.have || !(tau > 0.0)) {
        std::snprintf(buf, (std::size_t)n, "na");
        return;
    }
    if (!std::isfinite(s.circDev))
        std::snprintf(buf, (std::size_t)n, "circTau=na,lineTau=%.4f", s.lineRes / tau);
    else
        std::snprintf(buf, (std::size_t)n, "circTau=%.4f,lineTau=%.4f", s.circDev / tau, s.lineRes / tau);
}

bool sideIsProfileCircle(const SideFit& s, double tau) {
    return s.have && std::isfinite(s.circDev) && s.circDev <= tau && s.lineRes > tau;
}

// D-fc-5.1. A coaxial family in Line mode: its donors' certified vertices
// (D-fc-4.8 (b)) and every band claimed on it, all on one axis.
struct CoaxFamily {
    gp_XYZ origin, axis;
    std::vector<std::vector<gp_XYZ>> groups;  // donors first, then claimed bands
    std::size_t nDonors = 0;
    std::vector<int> bandRegion;  // work.accepted index of each claimed band
    bool moved = false;           // the axis changed after a band was placed on it
};

struct LineFamilies {
    std::vector<CoaxFamily> fams;
    std::vector<int> famOfAcc;  // work.accepted index -> family, -1 when not a donor
    const DerivedTols* tol = nullptr;
    const SegmentParams* sp = nullptr;
};

// D-fc-5.1 anchor ("planes define the arcs"): a certified plane (D-fc-1.4)
// parallel to the family axis over the band's length (tanAngTol), tangent to
// the band's cylinder (tanDistTol) and sharing a contact vertex with the band
// on it. A vertex test alone cannot tell an arc from an inscribed polygon
// (concyclic corners), so a band that no plane anchors is refused (no-anchor).
// Returns the plane's provisional index, -1 when none.
int coaxialTangentPlane(const WorkMesh& wm, const SegmentWork& work, const SegmentParams& sp,
                        const DerivedTols& tol, const std::vector<int>& bandTris,
                        const gp_XYZ& o, const gp_XYZ& a, double R, double& tanDist) {
    const MeshView& mv = *wm.mv;
    const int nTri = (int)mv.nTri;
    std::vector<int> bg;
    double zLo = 1e300, zHi = -1e300;
    for (int t : bandTris) {
        if (t < 0 || t >= nTri) continue;
        const int gt = mv.compTris ? mv.compTris[t] : t;
        for (int k = 0; k < 3; ++k) {
            const int g = mv.tris[gt][k];
            bg.push_back(g);
            const double z = (mv.pts[g] - o).Dot(a);
            zLo = std::min(zLo, z);
            zHi = std::max(zHi, z);
        }
    }
    std::sort(bg.begin(), bg.end());
    bg.erase(std::unique(bg.begin(), bg.end()), bg.end());
    const double L = std::max(zHi - zLo, wm.tau);
    tanDist = -1;
    for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
        PlaneAnchor pa;
        if (!planeAnchorOf(mv, sp, tol, work.provisionals[(std::size_t)pi], pa) || !pa.certified)
            continue;
        const double tolT = pa.tol;
        const gp_XYZ n = pa.n;
        if (std::fabs(n.Dot(a)) > tanAngTol(tolT, wm.tau, L)) continue;
        const double d = n.Dot(pa.p0);
        const double td = std::fabs(std::fabs(n.Dot(o) - d) - R);
        if (td > tanDistTol(tolT, wm.tau)) continue;
        bool contact = false;
        for (int t : pa.body) {
            if (t < 0 || t >= nTri) continue;
            const int gt = mv.compTris ? mv.compTris[t] : t;
            for (int k = 0; k < 3 && !contact; ++k) {
                const int g = mv.tris[gt][k];
                if (std::binary_search(bg.begin(), bg.end(), g) &&
                    std::fabs(n.Dot(mv.pts[g]) - d) <= tanDistTol(tolT, wm.tau))
                    contact = true;
            }
            if (contact) break;
        }
        if (!contact) continue;
        tanDist = td;
        return pi;
    }
    return -1;
}

// The cylinder fields of a coaxial band about (origin, axis): radius = the mean
// rho of its own vertices, frame XDirection = the radial direction of its
// first vertex, u span by the existing column rule (closedOf), v relative to
// the frame's location (the band's area centroid on the axis).
void placeCoaxial(const MeshView& mv, Region& r, const std::vector<int>& tris,
                  const std::vector<gp_XYZ>& pts, const gp_XYZ& origin, const gp_XYZ& axisIn,
                  double thetaBin) {
    gp_XYZ axis = axisIn;
    axis /= axis.Modulus();
    double sum = 0;
    for (const gp_XYZ& p : pts) {
        const gp_XYZ d = p - origin;
        sum += (d - axis * d.Dot(axis)).Modulus();
    }
    const double R = pts.empty() ? 0.0 : sum / (double)pts.size();
    gp_XYZ cent(0, 0, 0);
    double areaSum = 0;
    for (int t : tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
        areaSum += a;
        cent += triCentroid(mv, t) * a;
    }
    if (areaSum > 0.0) cent /= areaSum;
    const gp_XYZ loc = origin + axis * ((cent - origin).Dot(axis));
    gp_XYZ xHint(1, 0, 0);
    for (const gp_XYZ& p : pts) {
        gp_XYZ rad = p - loc;
        rad -= axis * rad.Dot(axis);
        if (rad.Modulus() > gp::Resolution()) {
            xHint = rad;
            break;
        }
    }
    const Frame fr = makeFrame(loc, axis, xHint);
    std::vector<Vtx> uniq;
    uniq.reserve(pts.size());
    double vMin = 1e300, vMax = -1e300, maxD = 0, sum2 = 0;
    for (const gp_XYZ& p : pts) {
        Vtx v;
        v.p = p;
        v.pr = profOfPoint(p, fr, 0);
        uniq.push_back(v);
        vMin = std::min(vMin, v.pr.z);
        vMax = std::max(vMax, v.pr.z);
        const double d = std::fabs(v.pr.rho - R);
        maxD = std::max(maxD, d);
        sum2 += d * d;
    }
    const ClosedM cl = closedOf(mv, fr, tris, uniq, thetaBin);
    double sigma = 0;
    for (int t : tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
        gp_XYZ rad = triCentroid(mv, t) - loc;
        rad -= axis * rad.Dot(axis);
        const double rm = rad.Modulus();
        if (rm > gp::Resolution()) sigma += a * n.Dot(rad / rm);
    }
    const bool outward = sigma > 0.0;
    const int nSides = std::max(cl.nSides, 1);
    r.type = SurfType::Cylinder;
    r.origin = Origin::FilletStrip;
    r.ax = gp_Ax3(gp_Pnt(loc.X(), loc.Y(), loc.Z()), gp_Dir(axis.X(), axis.Y(), axis.Z()),
                  gp_Dir(fr.e1.X(), fr.e1.Y(), fr.e1.Z()));
    r.radius = R;
    r.uMin = cl.closed ? 0.0 : cl.uMin;
    r.uMax = cl.closed ? kTwoPi : cl.uMax;
    r.vMin = (vMin > vMax) ? 0.0 : vMin;
    r.vMax = (vMin > vMax) ? 0.0 : vMax;
    r.closed360 = cl.closed;
    r.outwardNormal = outward;
    r.tris = tris;
    r.maxVertexDev = maxD;
    r.rmsVertexDev = pts.empty() ? 0.0 : std::sqrt(sum2 / (double)pts.size());
    r.nSides = nSides;
    r.chordSagitta = chordSagitta(R, nSides);
    r.dVolPredicted = dVolCylinderSector(areaSum, R, nSides, outward);
    r.maxVertexSnap = 0.0;
    r.reject = Reject::None;
    r.builtAs = BuiltAs::NotBuilt;
    r.filletNbrA = -1;
    r.filletNbrB = -1;
}

// D-fc-1.6.11 for a coaxial band: an accepted region beside the band whose
// every vertex is within tau of the band's cylinder is a second face of it.
bool sameSurfaceCyl(const WorkMesh& wm, const SegmentWork& work, const Band& band,
                    const gp_XYZ& origin, const gp_XYZ& axis, double R) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> bandTri((std::size_t)nTri, 0);
    for (int t : band.tris)
        if (t >= 0 && t < nTri) bandTri[(std::size_t)t] = 1;
    for (const Region& r : work.accepted) {
        if (r.tris.empty()) continue;
        bool touch = false;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int n : wm.adj.nbr[(std::size_t)t])
                if (n >= 0 && n < nTri && bandTri[(std::size_t)n]) touch = true;
            if (touch) break;
        }
        if (!touch) continue;
        bool all = true;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int k = 0; k < 3 && all; ++k) {
                const gp_XYZ d = triVert(*wm.mv, t, k) - origin;
                const double rho = (d - axis * d.Dot(axis)).Modulus();
                if (!(std::fabs(rho - R) <= wm.tau)) all = false;
            }
            if (!all) break;
        }
        if (all) return true;
    }
    return false;
}

// D-fc-3.5 line mode. Rings are vertices at equal (ρ, z) within τ, ≥ 3.
// A band is a maximal collinear chain of those rings. Edge adjacency to the
// donor is not required. D-fc-5.1: a chain is first tested as a cylinder on
// its donor's coaxial family (one parameter, the radius, about the family
// axis); only a chain that is not one goes on to the cone (fewest parameters).
void gatherLines(WorkMesh& wm, SegmentWork& work, Donor& donor, std::vector<AxisScan>& scans,
                 LineFamilies& lf, int& claimed, int& refused) {
    const MeshView& mv = *wm.mv;
    const int nTri = (int)mv.nTri;
    const double tau = wm.tau;
    if (!(tau > 0.0) || nTri <= 0) return;
    // A later donor shares this gather only when its axis stays inside τ over
    // the whole mesh. Donor length alone is too short: a sliver's angle
    // tolerance swings a distant ring by many τ.
    const double axisLen = std::max(donor.length, tau);
    const double reach = std::max(mv.diag, axisLen);
    for (const AxisScan& s : scans) {
        const double L = std::max(reach, std::max(s.length, tau));
        if (s.scanned && sameAxis(s.fr, donor.fr, tau, L)) return;
    }
    AxisScan scan;
    scan.fr = donor.fr;
    scan.length = donor.length;
    const Frame& fr = donor.fr;
    int maxV = 0;
    for (int t = 0; t < nTri; ++t) {
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            if (v > maxV) maxV = v;
        }
    }
    std::vector<double> rhoOf((std::size_t)maxV + 1, 0.0), zOf((std::size_t)maxV + 1, 0.0);
    std::vector<char> seenV((std::size_t)maxV + 1, 0);
    std::vector<int> allV;
    allV.reserve((std::size_t)nTri);
    for (int t = 0; t < nTri; ++t) {
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            if (v < 0 || v > maxV || seenV[(std::size_t)v]) continue;
            const ProfPt p = profOfPoint(mv.pts[v], fr, 0);
            if (!std::isfinite(p.rho) || !std::isfinite(p.z)) continue;
            seenV[(std::size_t)v] = 1;
            rhoOf[(std::size_t)v] = p.rho;
            zOf[(std::size_t)v] = p.z;
            allV.push_back(v);
        }
    }

    struct Cell {
        int parent = 0;
        long long kr = 0, kz = 0;
        double rMin = 0, rMax = 0, zMin = 0, zMax = 0;
        std::vector<int> verts;
    };
    std::vector<Cell> cells;
    std::map<std::pair<long long, long long>, int> buckets;
    for (int v : allV) {
        const double rho = rhoOf[(std::size_t)v];
        const double z = zOf[(std::size_t)v];
        const long long kr = std::llround(rho / tau);
        const long long kz = std::llround(z / tau);
        const auto key = std::make_pair(kr, kz);
        auto it = buckets.find(key);
        int id;
        if (it == buckets.end()) {
            id = (int)cells.size();
            Cell c;
            c.parent = id;
            c.kr = kr;
            c.kz = kz;
            c.rMin = c.rMax = rho;
            c.zMin = c.zMax = z;
            c.verts.push_back(v);
            cells.push_back(std::move(c));
            buckets.emplace(key, id);
        } else {
            Cell& c = cells[(std::size_t)it->second];
            c.verts.push_back(v);
            c.rMin = std::min(c.rMin, rho);
            c.rMax = std::max(c.rMax, rho);
            c.zMin = std::min(c.zMin, z);
            c.zMax = std::max(c.zMax, z);
        }
    }
    auto findCell = [&](int x) {
        while (cells[(std::size_t)x].parent != x) {
            cells[(std::size_t)x].parent = cells[(std::size_t)cells[(std::size_t)x].parent].parent;
            x = cells[(std::size_t)x].parent;
        }
        return x;
    };
    auto uniteCell = [&](int a, int b) {
        a = findCell(a);
        b = findCell(b);
        if (a == b) return false;
        const double rMin = std::min(cells[(std::size_t)a].rMin, cells[(std::size_t)b].rMin);
        const double rMax = std::max(cells[(std::size_t)a].rMax, cells[(std::size_t)b].rMax);
        const double zMin = std::min(cells[(std::size_t)a].zMin, cells[(std::size_t)b].zMin);
        const double zMax = std::max(cells[(std::size_t)a].zMax, cells[(std::size_t)b].zMax);
        if (rMax - rMin > tau || zMax - zMin > tau) return false;
        cells[(std::size_t)b].parent = a;
        cells[(std::size_t)a].rMin = rMin;
        cells[(std::size_t)a].rMax = rMax;
        cells[(std::size_t)a].zMin = zMin;
        cells[(std::size_t)a].zMax = zMax;
        return true;
    };
    if (!cells.empty()) {
        bool changed = true;
        for (int guard = 0; changed && guard < (int)cells.size() + 2; ++guard) {
            changed = false;
            for (int i = 0; i < (int)cells.size(); ++i) {
                for (int dr = -1; dr <= 1; ++dr) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        if (dr == 0 && dz == 0) continue;
                        const auto it = buckets.find(
                            std::make_pair(cells[(std::size_t)i].kr + dr, cells[(std::size_t)i].kz + dz));
                        if (it == buckets.end()) continue;
                        if (uniteCell(i, it->second)) changed = true;
                    }
                }
            }
        }
    }
    struct Ring {
        double rho = 0, z = 0;
        int minV = 0;
        std::vector<int> verts;
    };
    std::vector<Ring> rings;
    std::vector<int> ringOf((std::size_t)maxV + 1, -1);
    {
        std::vector<std::vector<int>> groups(cells.size());
        for (int i = 0; i < (int)cells.size(); ++i) {
            const int r = findCell(i);
            auto& g = groups[(std::size_t)r];
            g.insert(g.end(), cells[(std::size_t)i].verts.begin(), cells[(std::size_t)i].verts.end());
        }
        for (const auto& g : groups) {
            if ((int)g.size() < 3) continue;
            Ring ring;
            ring.minV = g.front();
            double sr = 0, sz = 0;
            for (int v : g) {
                sr += rhoOf[(std::size_t)v];
                sz += zOf[(std::size_t)v];
                ring.minV = std::min(ring.minV, v);
            }
            ring.rho = sr / (double)g.size();
            ring.z = sz / (double)g.size();
            ring.verts = g;
            rings.push_back(std::move(ring));
        }
    }
    std::sort(rings.begin(), rings.end(), [](const Ring& a, const Ring& b) {
        if (a.z != b.z) return a.z < b.z;
        if (a.rho != b.rho) return a.rho < b.rho;
        return a.minV < b.minV;
    });
    for (int i = 0; i < (int)rings.size(); ++i)
        for (int v : rings[(std::size_t)i].verts)
            if (v >= 0 && v <= maxV) ringOf[(std::size_t)v] = i;
    const int nRings = (int)rings.size();
    auto triRings = [&](int t, int out[3], int& nOut) {
        nOut = 0;
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            if (v < 0 || v > maxV) continue;
            const int r = ringOf[(std::size_t)v];
            if (r < 0) continue;
            bool dup = false;
            for (int i = 0; i < nOut; ++i)
                if (out[i] == r) dup = true;
            if (!dup && nOut < 3) out[nOut++] = r;
        }
    };
    auto vertOnRing = [&](int t) {
        int on = 0;
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            if (v >= 0 && v <= maxV && ringOf[(std::size_t)v] >= 0) ++on;
        }
        return on;
    };

    std::vector<std::vector<int>> fullAdj((std::size_t)std::max(nRings, 0));
    struct Cand {
        int t = -1;
        int sheet = 0;
        int nR = 0;
        int ids[3] = {-1, -1, -1};
    };
    std::vector<Cand> cands;
    for (int t = 0; t < nTri; ++t) {
        int ids[3], nIds = 0;
        triRings(t, ids, nIds);
        if (nIds >= 2) {
            for (int a = 0; a < nIds; ++a) {
                for (int b = a + 1; b < nIds; ++b) {
                    fullAdj[(std::size_t)ids[a]].push_back(ids[b]);
                    fullAdj[(std::size_t)ids[b]].push_back(ids[a]);
                }
            }
        }
        if (!wm.freeTri[(std::size_t)t]) continue;
        if (vertOnRing(t) < 3 || nIds < 2) continue;
        const int sh = sheetOf(mv, t, fr);
        if (sh == 0) continue;
        Cand c;
        c.t = t;
        c.sheet = sh;
        c.nR = nIds;
        for (int i = 0; i < nIds; ++i) c.ids[i] = ids[i];
        cands.push_back(c);
    }
    for (auto& nbr : fullAdj) {
        std::sort(nbr.begin(), nbr.end());
        nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
    }

    struct Edge {
        int u = 0, v = 0;
        int minTri = 0;
        bool used = false;
        int sheet = 0;
    };
    std::vector<Edge> edges;
    std::map<std::tuple<int, int, int>, int> edgeKey;
    for (const Cand& c : cands) {
        for (int a = 0; a < c.nR; ++a) {
            for (int b = a + 1; b < c.nR; ++b) {
                int u = c.ids[a], v = c.ids[b];
                if (u > v) std::swap(u, v);
                const auto key = std::make_tuple(c.sheet, u, v);
                auto it = edgeKey.find(key);
                if (it == edgeKey.end()) {
                    Edge e;
                    e.u = u;
                    e.v = v;
                    e.sheet = c.sheet;
                    e.minTri = c.t;
                    edgeKey.emplace(key, (int)edges.size());
                    edges.push_back(e);
                } else if (c.t < edges[(std::size_t)it->second].minTri) {
                    edges[(std::size_t)it->second].minTri = c.t;
                }
            }
        }
    }
    std::vector<std::vector<int>> edgeAt((std::size_t)std::max(nRings, 0));
    for (int i = 0; i < (int)edges.size(); ++i) {
        edgeAt[(std::size_t)edges[(std::size_t)i].u].push_back(i);
        edgeAt[(std::size_t)edges[(std::size_t)i].v].push_back(i);
    }
    auto chainCollinear = [&](const std::vector<int>& ids) {
        if ((int)ids.size() <= 2) return true;
        std::vector<ProfPt> pts;
        pts.reserve(ids.size());
        for (int id : ids) {
            ProfPt p;
            p.rho = rings[(std::size_t)id].rho;
            p.z = rings[(std::size_t)id].z;
            pts.push_back(p);
        }
        return profileLineMaxResid(pts) <= tau;
    };
    struct Chain {
        std::vector<int> rings;
        std::vector<int> tris;
        int sheet = 0;
        int seed = 0;
    };
    std::vector<Chain> chains;
    std::vector<int> edgeOrder(edges.size());
    for (int i = 0; i < (int)edges.size(); ++i) edgeOrder[(std::size_t)i] = i;
    std::sort(edgeOrder.begin(), edgeOrder.end(), [&](int a, int b) {
        if (edges[(std::size_t)a].minTri != edges[(std::size_t)b].minTri)
            return edges[(std::size_t)a].minTri < edges[(std::size_t)b].minTri;
        if (edges[(std::size_t)a].u != edges[(std::size_t)b].u)
            return edges[(std::size_t)a].u < edges[(std::size_t)b].u;
        return edges[(std::size_t)a].v < edges[(std::size_t)b].v;
    });
    for (int ei0 : edgeOrder) {
        if (edges[(std::size_t)ei0].used) continue;
        edges[(std::size_t)ei0].used = true;
        std::vector<int> chain{edges[(std::size_t)ei0].u, edges[(std::size_t)ei0].v};
        const int sheet = edges[(std::size_t)ei0].sheet;
        auto grow = [&](bool front) {
            while ((int)chain.size() >= 2) {
                const int cur = front ? chain.front() : chain.back();
                const int prev = front ? chain[1] : chain[(int)chain.size() - 2];
                int bestE = -1, bestN = -1, bestMin = std::numeric_limits<int>::max();
                for (int ei : edgeAt[(std::size_t)cur]) {
                    Edge& e = edges[(std::size_t)ei];
                    if (e.used || e.sheet != sheet) continue;
                    const int nxt = (e.u == cur) ? e.v : e.u;
                    if (nxt == prev) continue;
                    std::vector<int> trial = chain;
                    if (front) trial.insert(trial.begin(), nxt);
                    else trial.push_back(nxt);
                    if (!chainCollinear(trial)) continue;
                    if (e.minTri < bestMin || (e.minTri == bestMin && (bestN < 0 || nxt < bestN))) {
                        bestMin = e.minTri;
                        bestE = ei;
                        bestN = nxt;
                    }
                }
                if (bestE < 0) break;
                edges[(std::size_t)bestE].used = true;
                if (front) chain.insert(chain.begin(), bestN);
                else chain.push_back(bestN);
            }
        };
        grow(false);
        grow(true);
        std::vector<char> inR((std::size_t)nRings, 0);
        for (int id : chain)
            if (id >= 0 && id < nRings) inR[(std::size_t)id] = 1;
        for (Edge& e : edges)
            if (e.sheet == sheet && inR[(std::size_t)e.u] && inR[(std::size_t)e.v]) e.used = true;
        Chain ch;
        ch.rings = chain;
        ch.sheet = sheet;
        for (const Cand& c : cands) {
            if (c.sheet != sheet || c.nR < 2) continue;
            bool inside = true;
            for (int i = 0; i < c.nR; ++i)
                if (c.ids[i] < 0 || c.ids[i] >= nRings || !inR[(std::size_t)c.ids[i]]) inside = false;
            if (inside) ch.tris.push_back(c.t);
        }
        if (ch.tris.empty()) continue;
        std::sort(ch.tris.begin(), ch.tris.end());
        ch.tris.erase(std::unique(ch.tris.begin(), ch.tris.end()), ch.tris.end());
        ch.seed = minTriOf(ch.tris);
        chains.push_back(std::move(ch));
    }
    std::sort(chains.begin(), chains.end(), [](const Chain& a, const Chain& b) {
        if (a.seed != b.seed) return a.seed < b.seed;
        return a.sheet < b.sheet;
    });

    auto pickNext = [&](int cur, int prev, int extra) {
        std::vector<int> opts;
        if (cur < 0 || cur >= nRings) return -1;
        for (int n : fullAdj[(std::size_t)cur]) {
            if (n == prev || n == extra) continue;
            opts.push_back(n);
        }
        if (opts.empty()) return -1;
        if (opts.size() == 1) return opts[0];
        if (prev < 0 || prev >= nRings) return *std::min_element(opts.begin(), opts.end());
        const double dx = rings[(std::size_t)cur].rho - rings[(std::size_t)prev].rho;
        const double dz = rings[(std::size_t)cur].z - rings[(std::size_t)prev].z;
        int best = -1;
        double bestF = -1e300;
        for (int n : opts) {
            const double fwd = (rings[(std::size_t)n].rho - rings[(std::size_t)cur].rho) * dx +
                               (rings[(std::size_t)n].z - rings[(std::size_t)cur].z) * dz;
            if (best < 0 || fwd > bestF || (fwd == bestF && n < best)) {
                best = n;
                bestF = fwd;
            }
        }
        if (!(bestF > 0.0)) return -1;
        return best;
    };
    auto measureSide = [&](int end, int other) {
        SideFit s;
        const int n1 = pickNext(end, other, -1);
        if (n1 < 0) return s;
        const int n2 = pickNext(n1, end, other);
        if (n2 < 0) return s;
        const int seq[4] = {other, end, n1, n2};
        std::vector<ProfPt> pts;
        pts.reserve(4);
        for (int id : seq) {
            if (id < 0 || id >= nRings) return SideFit();
            ProfPt p;
            p.rho = rings[(std::size_t)id].rho;
            p.z = rings[(std::size_t)id].z;
            pts.push_back(p);
        }
        s.have = true;
        s.lineRes = profileLineMaxResid(pts);
        double Rmaj = 0, z0 = 0, Rmin = 0;
        if (!fitProfileCircle(pts, tau, Rmaj, z0, Rmin) || !std::isfinite(Rmin)) {
            s.circDev = std::numeric_limits<double>::infinity();
            return s;
        }
        s.circDev = 0;
        for (const ProfPt& p : pts)
            s.circDev = std::max(s.circDev, std::fabs(std::hypot(p.rho - Rmaj, p.z - z0) - Rmin));
        return s;
    };

    for (const Chain& ch : chains) {
        std::vector<int> live;
        for (int t : ch.tris)
            if (t >= 0 && t < nTri && wm.freeTri[(std::size_t)t]) live.push_back(t);
        if (live.empty()) continue;
        Band band;
        band.in.assign((std::size_t)nTri, 0);
        band.tris = live;
        band.sheet = ch.sheet;
        band.seed = minTriOf(live);
        for (int t : live)
            if (t >= 0 && t < nTri) band.in[(std::size_t)t] = 1;
        collectBand(mv, fr, wm.g2l, band);
        if (band.pts.empty()) continue;
        const bool two = (int)ch.rings.size() == 2;
        SideFit s0, s1;
        if (two) {
            s0 = measureSide(ch.rings[0], ch.rings[1]);
            s1 = measureSide(ch.rings[1], ch.rings[0]);
        }
        char sideBuf0[96], sideBuf1[96];
        fmtSide(sideBuf0, (int)sizeof(sideBuf0), s0, tau);
        fmtSide(sideBuf1, (int)sizeof(sideBuf1), s1, tau);
        const SurfaceFit fit = fitSurface(band.pts, tau);
        const ClosedM cl = closedOf(mv, fr, band.tris, band.verts, wm.thetaBin);
        std::vector<double> zs, psis, rhos;
        for (const ProfPt& pp : band.pts) {
            zs.push_back(pp.z);
            psis.push_back(pp.psi);
            rhos.push_back(pp.rho);
        }
        const int levels = countLevels(zs, tau);
        const int stations = countProfPts(band.pts, tau);
        double maxP = 0;
        for (const ProfPt& pp : band.pts) maxP = std::max(maxP, linePerp(fit, pp));
        double rLoS = 1e300, rHiS = -1e300;
        for (int id : ch.rings) {
            if (id < 0 || id >= nRings) continue;
            rLoS = std::min(rLoS, rings[(std::size_t)id].rho);
            rHiS = std::max(rHiS, rings[(std::size_t)id].rho);
        }
        auto refuse = [&](const char* why) {
            ++refused;
            DiagRec d = baseDiag("line", donor, band, tau);
            d.stations = stations;
            d.levels = levels;
            d.params = kParamsLine;
            d.maxDev = std::isfinite(maxP) ? maxP : 0;
            d.maxDevTau = tau > 0 ? d.maxDev / tau : 0;
            d.lineResTau = (tau > 0 && fit.hatOk) ? fit.orthMax / tau : 0;
            d.span = cl.span;
            d.spanTol = (cl.nBands > 0) ? 1.5 * kTwoPi / cl.nBands : 0;
            d.halfAng = fit.hatOk ? std::atan(std::fabs(fit.b)) * 180.0 / kPi : 0;
            d.dihedral = d.halfAng;
            d.r0 = rLoS;
            d.r1 = rHiS;
            d.side0 = sideBuf0;
            d.side1 = sideBuf1;
            d.why = why;
            printBand(d);
        };
        if (two && (sideIsProfileCircle(s0, tau) || sideIsProfileCircle(s1, tau))) {
            refuse("profile-row");
            continue;
        }
        // D-fc-5.1. A cylinder on the donor's coaxial family first: one free
        // parameter (its radius) about the family axis, certified jointly with
        // the family's certified vertices. A chain that is not one goes on to
        // the cone below unchanged.
        const int famId = (donor.idx >= 0 && donor.idx < (int)lf.famOfAcc.size())
                              ? lf.famOfAcc[(std::size_t)donor.idx]
                              : -1;
        if (famId >= 0) {
            CoaxFamily& fam = lf.fams[(std::size_t)famId];
            // The claimed set is the complete tau-preimage of the certified
            // cylinder that is edge-connected to the chain (D-fc-3.2.4: every
            // free triangle with 3 vertices within tau and >= 2 axial levels),
            // refit and re-certified to a fixed point (D-fc-3.2.5). Rings split
            // at tau leave such triangles outside every chain. A refused
            // cylinder leaves `band` as the chain for the cone below.
            Band cb;
            cb.tris = band.tris;
            cb.verts = band.verts;
            cb.pts = band.pts;
            cb.sheet = band.sheet;
            cb.seed = band.seed;
            const int nChain = (int)cb.tris.size();
            auto member = [&](int t, const gp_XYZ& o, const gp_XYZ& a, double R) {
                if (t < 0 || t >= nTri || !wm.freeTri[(std::size_t)t]) return false;
                double zLo = 1e300, zHi = -1e300;
                for (int k = 0; k < 3; ++k) {
                    const gp_XYZ d = triVert(mv, t, k) - o;
                    const double z = d.Dot(a);
                    if (!(std::fabs((d - a * z).Modulus() - R) <= tau)) return false;
                    zLo = std::min(zLo, z);
                    zHi = std::max(zHi, z);
                }
                return zHi - zLo > tau;
            };
            std::vector<gp_XYZ> bpts;
            bandPoints(cb, bpts);
            CoaxialCert cc =
                certifyCoaxialBand(fam.groups, bpts, fam.origin, fam.axis, tau, wm.maxIter);
            bool stable = false;
            int iters = 0;
            while (cc.ok && iters < wm.maxIter) {
                ++iters;
                const double R = cc.fit.radius.back();
                std::vector<char> inI((std::size_t)nTri, 0);
                std::vector<int> I, stack;
                for (int t : cb.tris)
                    if (member(t, cc.fit.origin, cc.fit.axis, R)) {
                        inI[(std::size_t)t] = 1;
                        stack.push_back(t);
                    }
                while (!stack.empty()) {
                    const int c = stack.back();
                    stack.pop_back();
                    I.push_back(c);
                    for (int n : wm.adj.nbr[(std::size_t)c]) {
                        if (n < 0 || n >= nTri || inI[(std::size_t)n]) continue;
                        if (!member(n, cc.fit.origin, cc.fit.axis, R)) continue;
                        inI[(std::size_t)n] = 1;
                        stack.push_back(n);
                    }
                }
                std::sort(I.begin(), I.end());
                if (I == cb.tris) {
                    stable = true;
                    break;
                }
                cb.tris = I;
                cb.seed = minTriOf(I);
                collectBand(mv, fr, wm.g2l, cb);
                bandPoints(cb, bpts);
                cc = certifyCoaxialBand(fam.groups, bpts, fam.origin, fam.axis, tau, wm.maxIter);
            }
            const double Rb = cc.fit.radius.empty() ? 0.0 : cc.fit.radius.back();
            const char* why = cc.why;
            int anchorProv = -1;
            double anchorDist = -1;
            const ClosedM ccl = cc.ok ? closedOf(mv, fr, cb.tris, cb.verts, wm.thetaBin) : ClosedM();
            if (cc.ok && !stable) why = "unstable";
            else if (cc.ok && ccl.closed) why = "closed-deferred";
            else if (cc.ok && coaxialSameRadius(cc.fit, fam.groups.size(), tau)) why = "same-surface";
            else if (cc.ok && lf.tol && lf.sp &&
                     (anchorProv = coaxialTangentPlane(wm, work, *lf.sp, *lf.tol, cb.tris,
                                                       cc.fit.origin, cc.fit.axis, Rb,
                                                       anchorDist)) < 0)
                why = "no-anchor";
            else if (cc.ok && sameSurfaceCyl(wm, work, cb, cc.fit.origin, cc.fit.axis, Rb))
                why = "same-surface";
            else if (cc.ok) why = "accept";
            const bool accept = std::strcmp(why, "accept") == 0;
            int straddle = 0;
            if (accept) {
                Region r;
                r.id = 0;
                placeCoaxial(mv, r, cb.tris, bpts, cc.fit.origin, cc.fit.axis, wm.thetaBin);
                straddle = countStraddle(work, cb.tris, nTri);
                std::vector<int> touched;
                trimProvisionals(mv, work, cb.tris, nTri, ProvClaim::ConsumedFillet, touched);
                if (cc.refit && !fam.bandRegion.empty()) fam.moved = true;
                fam.origin = cc.fit.origin;
                fam.axis = cc.fit.axis;
                fam.groups.push_back(bpts);
                fam.bandRegion.push_back((int)work.accepted.size());
                work.accepted.push_back(std::move(r));
                for (int t : cb.tris)
                    if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
                ++claimed;
            }
            if (diagOn()) {
                double famDev = 0;
                for (std::size_t g = 0; g + 1 < cc.fit.maxDev.size(); ++g)
                    famDev = std::max(famDev, cc.fit.maxDev[g]);
                const double bandDev = cc.fit.maxDev.empty() ? 0.0 : cc.fit.maxDev.back();
                DiagRec d = baseDiag("coaxial", donor, cb, tau);
                d.stations = cc.stations;
                d.levels = cc.levels;
                d.params = kParamsAnchored;
                d.fit = cc.refit ? "joint-ls" : "family";
                d.r0 = Rb;
                d.axis = cc.fit.axis;
                d.maxDev = bandDev;
                d.maxDevTau = bandDev / tau;
                d.lineResTau = cc.lineRes / tau;
                d.straddle = straddle;
                d.why = why;
                printBand(d);
                const gp_XYZ off = cc.fit.origin - donor.fr.origin;
                std::fprintf(stderr,
                             "DIAG_REL_COAX seed=%d donor=%d family=%d groups=%zu donors=%zu "
                             "refit=%d familyDevTau=%.4f bandDevTau=%.4f radius=%.9g "
                             "tilt=%.4g offset=%.4g nChain=%d nTri=%zu iters=%d anchorProv=%d "
                             "anchorDistTau=%.4f closed=%d why=%s\n",
                             cb.seed, donor.idx, famId, fam.groups.size() - (accept ? 1 : 0),
                             fam.nDonors, cc.refit ? 1 : 0, famDev / tau, bandDev / tau, Rb,
                             cc.fit.axis.Crossed(donor.fr.axis).Modulus(),
                             (off - donor.fr.axis * off.Dot(donor.fr.axis)).Modulus(),
                             nChain, cb.tris.size(), iters, anchorProv,
                             anchorDist >= 0 ? anchorDist / tau : -1.0, ccl.closed ? 1 : 0, why);
                if (accept) {
                    std::fprintf(stderr, "DIAG_REL_COAX_TRIS seed=%d tris=", cb.seed);
                    for (int t : cb.tris) std::fprintf(stderr, "%d,", t);
                    std::fprintf(stderr, "\n");
                }
            }
            if (accept) continue;
        }
        if (stations < 2 || !fit.hatOk) {
            refuse("underdetermined");
            continue;
        }
        if (!(fit.orthMax <= tau)) {
            refuse("vertex-off");
            continue;
        }
        double zMin = band.pts[0].z, zMax = zMin, rMin = band.pts[0].rho, rMax = rMin;
        for (const ProfPt& pp : band.pts) {
            zMin = std::min(zMin, pp.z);
            zMax = std::max(zMax, pp.z);
            rMin = std::min(rMin, pp.rho);
            rMax = std::max(rMax, pp.rho);
        }
        if (zMax - zMin <= tau && rMax - rMin <= tau) {
            refuse("underdetermined");
            continue;
        }
        if (zMax - zMin <= tau) {
            refuse("line");
            continue;
        }
        if (rMax - rMin <= tau) {
            refuse("same-surface");
            continue;
        }
        if (!(maxP <= tau)) {
            refuse("vertex-off");
            continue;
        }
        if (sameSurfaceFr(wm, work, fr, band, fit, false)) {
            refuse("same-surface");
            continue;
        }
        if (!cl.closed) {
            // D-fc-5.5. A band that does not close is a PARTIAL cone when it
            // certifies about its donor's certified coaxial family axis.
            if (famId < 0) {
                refuse("partial-cone");
                continue;
            }
            const CoaxFamily& fam = lf.fams[(std::size_t)famId];
            // Each edge-connected component of the band is certified and
            // claimed on its own (one surface piece, one outer loop), in
            // min-triangle order.
            std::vector<std::vector<int>> comps;
            {
                std::vector<char> seen((std::size_t)nTri, 0);
                for (int t0 : band.tris) {
                    if (t0 < 0 || t0 >= nTri || seen[(std::size_t)t0]) continue;
                    std::vector<int> comp, stack{t0};
                    seen[(std::size_t)t0] = 1;
                    while (!stack.empty()) {
                        const int c = stack.back();
                        stack.pop_back();
                        comp.push_back(c);
                        for (int n : wm.adj.nbr[(std::size_t)c]) {
                            if (n < 0 || n >= nTri || seen[(std::size_t)n] || !band.in[(std::size_t)n])
                                continue;
                            seen[(std::size_t)n] = 1;
                            stack.push_back(n);
                        }
                    }
                    std::sort(comp.begin(), comp.end());
                    comps.push_back(std::move(comp));
                }
                std::sort(comps.begin(), comps.end(),
                          [](const std::vector<int>& x, const std::vector<int>& y) {
                              return x.front() < y.front();
                          });
            }
            for (const std::vector<int>& comp : comps) {
                bool live = true;
                for (int t : comp)
                    if (t < 0 || t >= nTri || !wm.freeTri[(std::size_t)t]) live = false;
                if (!live) continue;  // an earlier component's preimage took it
                Band cb;
                cb.in.assign((std::size_t)nTri, 0);
                cb.tris = comp;
                cb.sheet = band.sheet;
                cb.seed = minTriOf(comp);
                for (int t : comp) cb.in[(std::size_t)t] = 1;
                collectBand(mv, fr, wm.g2l, cb);
                std::vector<gp_XYZ> bpts;
                bandPoints(cb, bpts);
                PartialConeCert pc =
                    certifyPartialCone(fam.groups, bpts, fam.origin, fam.axis, tau);
                // The claimed set is the complete tau-preimage of the certified
                // cone edge-connected to the chain (free triangles, 3 vertices
                // within tau of it, >= 2 axial levels), refit and re-certified to
                // a fixed point (D-fc-3.2.4/3.2.5, as D-fc-5.1): ring chains split
                // one surface into pieces.
                const Frame pf = makeFrame(fam.origin, fam.axis, gp_XYZ(1, 0, 0));
                auto onCone = [&](int t, const PartialConeCert& c) {
                    if (t < 0 || t >= nTri || !wm.freeTri[(std::size_t)t]) return false;
                    const double sc = std::sqrt(1.0 + c.k * c.k);
                    double zLo = 1e300, zHi = -1e300;
                    for (int k = 0; k < 3; ++k) {
                        const ProfPt pp = profOfPoint(triVert(mv, t, k), pf, 0);
                        // A vertex on the axis (the apex) has no azimuth.
                        if (!(pp.rho > tau)) return false;
                        if (!(std::fabs(pp.rho - c.rhoAt(pp.z)) / sc <= tau)) return false;
                        zLo = std::min(zLo, pp.z);
                        zHi = std::max(zHi, pp.z);
                    }
                    return zHi - zLo > tau;
                };
                bool pstable = !pc.ok;
                for (int it = 0; pc.ok && it < wm.maxIter; ++it) {
                    std::vector<char> inI((std::size_t)nTri, 0);
                    std::vector<int> I, stack;
                    for (int t : cb.tris)
                        if (onCone(t, pc)) {
                            inI[(std::size_t)t] = 1;
                            stack.push_back(t);
                        }
                    while (!stack.empty()) {
                        const int c = stack.back();
                        stack.pop_back();
                        I.push_back(c);
                        for (int n : wm.adj.nbr[(std::size_t)c]) {
                            if (n < 0 || n >= nTri || inI[(std::size_t)n] || !onCone(n, pc)) continue;
                            inI[(std::size_t)n] = 1;
                            stack.push_back(n);
                        }
                    }
                    std::sort(I.begin(), I.end());
                    if (I == cb.tris) {
                        pstable = true;
                        break;
                    }
                    if (I.empty()) break;
                    cb.tris = I;
                    cb.in.assign((std::size_t)nTri, 0);
                    for (int t : I) cb.in[(std::size_t)t] = 1;
                    cb.seed = minTriOf(I);
                    collectBand(mv, fr, wm.g2l, cb);
                    bandPoints(cb, bpts);
                    pc = certifyPartialCone(fam.groups, bpts, fam.origin, fam.axis, tau);
                }
                const char* why = pc.why;
                if (pc.ok && !pstable) why = "unstable";
                Frame ff = makeFrame(fam.origin, fam.axis, donor.fr.e1);
                if (pc.ok && pstable) {
                    Band fb = cb;
                    collectBand(mv, ff, wm.g2l, fb);
                    const SurfaceFit ffit = fitSurface(fb.pts, tau);
                    if (sameSurfaceFr(wm, work, ff, fb, ffit, false)) why = "same-surface";
                }
                double pspan = 0;
                int pstraddle = 0;
                const bool paccept = pc.ok && pstable && std::strcmp(why, "accept") == 0 &&
                                     emitPartialCone(wm, work, donor, cb, fam.origin, fam.axis,
                                                     pc, pspan, pstraddle);
                if (pc.ok && std::strcmp(why, "accept") == 0 && !paccept) why = "degenerate";
                if (diagOn()) {
                    DiagRec d = baseDiag("line", donor, cb, tau);
                    d.stations = pc.stations;
                    d.levels = pc.levels;
                    d.params = kParamsLine;
                    d.fit = "family";
                    d.axis = fam.axis;
                    d.maxDev = pc.lineRes;
                    d.maxDevTau = pc.lineRes / tau;
                    d.lineResTau = pc.planeRes / tau;
                    d.span = pspan;
                    d.halfAng = std::atan(std::fabs(pc.k)) * 180.0 / kPi;
                    d.dihedral = d.halfAng;
                    d.r0 = pc.rhoAt(pc.zLo);
                    d.r1 = pc.rhoAt(pc.zHi);
                    d.straddle = pstraddle;
                    d.side0 = sideBuf0;
                    d.side1 = sideBuf1;
                    d.why = paccept ? "accept"
                                    : (std::strcmp(why, "accept") == 0 ? "partial-cone" : why);
                    printBand(d);
                    std::fprintf(stderr,
                                 "DIAG_REL_PCONE seed=%d donor=%d family=%d nTri=%zu comps=%zu "
                                 "stations=%d levels=%d cols=%d famDevTau=%.4f lineResTau=%.4f "
                                 "planeResTau=%.4g r=%.9g/%.9g z=%.6f/%.6f span=%.6f why=%s\n",
                                 cb.seed, donor.idx, famId, cb.tris.size(), comps.size(),
                                 pc.stations, pc.levels, pc.cols, pc.famDev / tau,
                                 pc.lineRes / tau, pc.planeRes / tau, pc.rhoAt(pc.zLo),
                                 pc.rhoAt(pc.zHi), pc.zLo, pc.zHi, pspan,
                                 paccept ? "accept" : why);
                }
                if (!paccept) {
                    ++refused;
                    continue;
                }
                for (int t : cb.tris)
                    if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
                ++claimed;
            }
            continue;
        }
        {
            double best = 1e300;
            gp_XYZ bestR(0, 0, 0);
            bool have = false;
            for (const Vtx& v : band.verts) {
                gp_XYZ rad = v.p - fr.origin;
                rad -= fr.axis * rad.Dot(fr.axis);
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
            if (have) {
                // ρ is hypot(d·e1, d·e2). A later chain in this gather reads the
                // same donor frame, so e2 has to stay the complement of e1.
                donor.fr.e1 = bestR;
                gp_XYZ e2 = donor.fr.axis.Crossed(donor.fr.e1);
                const double em = e2.Modulus();
                if (em > gp::Resolution()) donor.fr.e2 = e2 / em;
            }
        }
        double ang0 = 0, ang1 = 0;
        int straddle = 0;
        if (!emitCone(wm, work, donor, band, fit, cl, ang0, ang1, straddle)) {
            refuse("degenerate");
            continue;
        }
        shareSeamFrame(wm, work, donor, band);
        updateDonor(work, donor);
        for (int t : band.tris)
            if (t >= 0 && t < nTri) wm.freeTri[(std::size_t)t] = 0;
        ++claimed;
        DiagRec d = baseDiag("line", donor, band, tau);
        d.stations = stations;
        d.levels = levels;
        d.params = kParamsLine;
        d.maxDev = maxP;
        d.maxDevTau = tau > 0 ? maxP / tau : 0;
        d.lineResTau = tau > 0 ? fit.orthMax / tau : 0;
        d.span = cl.span;
        d.spanTol = (cl.nBands > 0) ? 1.5 * kTwoPi / cl.nBands : 0;
        d.halfAng = std::atan(std::fabs(fit.b)) * 180.0 / kPi;
        d.dihedral = d.halfAng;
        d.r0 = work.accepted.back().radius;
        d.r1 = work.accepted.back().vMax;
        d.tan0 = ang0;
        d.tan1 = ang1;
        d.straddle = straddle;
        d.side0 = sideBuf0;
        d.side1 = sideBuf1;
        d.why = "accept";
        printBand(d);
    }

    // Triangles with some, but not all, vertices on rings are not a line band.
    // Score them (the elliptical cap) and do not claim them.
    std::vector<char> mixed((std::size_t)nTri, 0), seenT((std::size_t)nTri, 0);
    for (int t = 0; t < nTri; ++t) {
        if (!wm.freeTri[(std::size_t)t]) continue;
        const int on = vertOnRing(t);
        if (on > 0 && on < 3) mixed[(std::size_t)t] = 1;
    }
    for (int t = 0; t < nTri; ++t) {
        if (!mixed[(std::size_t)t] || seenT[(std::size_t)t]) continue;
        std::vector<int> comp, stack{t};
        seenT[(std::size_t)t] = 1;
        while (!stack.empty()) {
            const int c = stack.back();
            stack.pop_back();
            comp.push_back(c);
            for (int n : wm.adj.nbr[(std::size_t)c]) {
                if (n < 0 || n >= nTri || !mixed[(std::size_t)n] || seenT[(std::size_t)n]) continue;
                seenT[(std::size_t)n] = 1;
                stack.push_back(n);
            }
        }
        Band band;
        band.in.assign((std::size_t)nTri, 0);
        band.tris = std::move(comp);
        band.seed = minTriOf(band.tris);
        for (int c : band.tris)
            if (c >= 0 && c < nTri) band.in[(std::size_t)c] = 1;
        collectBand(mv, fr, wm.g2l, band);
        ++refused;
        DiagRec d = baseDiag("line", donor, band, tau);
        d.params = kParamsLine;
        d.why = "vertex-off";
        printBand(d);
    }

    scan.scanned = true;
    scans.push_back(scan);
}

// D-fc-4.6 F2: a free-band variant of sameSurfaceFr that does not count atoms.
// An atom is a candidate unit of the band, not a neighbour surface: inside the
// band it joins or leaves whole, beside it reclaimRelational(Torus) offers it
// the same test against the torus the band becomes.
bool sameSurfaceTorus(const WorkMesh& wm, const SegmentWork& work, const Frame& fr,
                      const Band& band, const SurfaceFit& fit, const std::vector<char>& atomDead) {
    const int nTri = (int)wm.mv->nTri;
    std::vector<char> bandTri((std::size_t)nTri, 0);
    for (int t : band.tris)
        if (t >= 0 && t < nTri) bandTri[(std::size_t)t] = 1;
    for (int i = 0; i < (int)work.accepted.size(); ++i) {
        const Region& r = work.accepted[(std::size_t)i];
        if (r.tris.empty() || torusAtom(r)) continue;
        if (i < (int)atomDead.size() && atomDead[(std::size_t)i]) continue;
        bool touch = false;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int n : wm.adj.nbr[(std::size_t)t])
                if (n >= 0 && n < nTri && bandTri[(std::size_t)n]) touch = true;
            if (touch) break;
        }
        if (!touch) continue;
        bool all = true;
        for (int t : r.tris) {
            if (t < 0 || t >= nTri) continue;
            for (int k = 0; k < 3 && all; ++k)
                if (!(circDist(fit, profOfPoint(triVert(*wm.mv, t, k), fr, 0)) <= wm.tau)) all = false;
            if (!all) break;
        }
        if (all) return true;
    }
    return false;
}

// D-fc-4.6 F1/F2. Torus mode gathers by axis like Line mode (D-fc-3.5.1).
// Rings are vertices at equal (rho, z) within tau about the donor axis, >= 3
// (kParamsLine + 1). Candidates are triangles with all three vertices on rings
// that are free or belong to an atom. A band is a chain of consecutive rings
// whose stations stay concyclic within tau (three stations fix a circle,
// kParamsCircle), split into edge-connected components; each component is
// certified by certifyTorus itself. An atom joins only whole and only with
// every vertex within tau of the band's torus, and over-determination is
// counted leave-one-region-out. Adjacency to the donor is not required.
void gatherTori(WorkMesh& wm, SegmentWork& work, Donor& donor, std::vector<AxisScan>& scans,
                const std::vector<int>& atomOf, std::vector<char>& atomDead, int& claimed,
                int& refused) {
    const MeshView& mv = *wm.mv;
    const int nTri = (int)mv.nTri;
    const double tau = wm.tau;
    if (!(tau > 0.0) || nTri <= 0) return;
    const double axisLen = std::max(donor.length, tau);
    const double reach = std::max(mv.diag, axisLen);
    for (const AxisScan& s : scans) {
        const double L = std::max(reach, std::max(s.length, tau));
        if (s.scanned && sameAxis(s.fr, donor.fr, tau, L)) return;
    }
    AxisScan scan;
    scan.fr = donor.fr;
    scan.length = donor.length;
    const Frame& fr = donor.fr;
    auto liveTri = [&](int t) {
        if (t < 0 || t >= nTri) return false;
        if (wm.freeTri[(std::size_t)t]) return true;
        const int a = atomOf[(std::size_t)t];
        if (a < 0 || a == donor.idx || atomDead[(std::size_t)a]) return false;
        return std::find(donor.members.begin(), donor.members.end(), a) == donor.members.end();
    };
    int maxV = 0;
    for (int t = 0; t < nTri; ++t)
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            maxV = std::max(maxV, mv.tris[gt][k]);
        }
    std::vector<double> rhoOf((std::size_t)maxV + 1, 0.0), zOf((std::size_t)maxV + 1, 0.0);
    std::vector<char> seenV((std::size_t)maxV + 1, 0);
    std::vector<int> allV;
    allV.reserve((std::size_t)nTri);
    for (int t = 0; t < nTri; ++t) {
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            if (v < 0 || v > maxV || seenV[(std::size_t)v]) continue;
            const ProfPt p = profOfPoint(mv.pts[v], fr, 0);
            if (!std::isfinite(p.rho) || !std::isfinite(p.z)) continue;
            seenV[(std::size_t)v] = 1;
            rhoOf[(std::size_t)v] = p.rho;
            zOf[(std::size_t)v] = p.z;
            allV.push_back(v);
        }
    }
    // Rings: the Line-mode cell rule (tau buckets, united while the union stays
    // within tau in rho and in z).
    struct Cell {
        int parent = 0;
        long long kr = 0, kz = 0;
        double rMin = 0, rMax = 0, zMin = 0, zMax = 0;
        std::vector<int> verts;
    };
    std::vector<Cell> cells;
    std::map<std::pair<long long, long long>, int> buckets;
    for (int v : allV) {
        const double rho = rhoOf[(std::size_t)v];
        const double z = zOf[(std::size_t)v];
        const auto key = std::make_pair((long long)std::llround(rho / tau), (long long)std::llround(z / tau));
        auto it = buckets.find(key);
        if (it == buckets.end()) {
            Cell c;
            c.parent = (int)cells.size();
            c.kr = key.first;
            c.kz = key.second;
            c.rMin = c.rMax = rho;
            c.zMin = c.zMax = z;
            c.verts.push_back(v);
            buckets.emplace(key, (int)cells.size());
            cells.push_back(std::move(c));
        } else {
            Cell& c = cells[(std::size_t)it->second];
            c.verts.push_back(v);
            c.rMin = std::min(c.rMin, rho);
            c.rMax = std::max(c.rMax, rho);
            c.zMin = std::min(c.zMin, z);
            c.zMax = std::max(c.zMax, z);
        }
    }
    auto findCell = [&](int x) {
        while (cells[(std::size_t)x].parent != x) {
            cells[(std::size_t)x].parent = cells[(std::size_t)cells[(std::size_t)x].parent].parent;
            x = cells[(std::size_t)x].parent;
        }
        return x;
    };
    auto uniteCell = [&](int a, int b) {
        a = findCell(a);
        b = findCell(b);
        if (a == b) return false;
        const double rMin = std::min(cells[(std::size_t)a].rMin, cells[(std::size_t)b].rMin);
        const double rMax = std::max(cells[(std::size_t)a].rMax, cells[(std::size_t)b].rMax);
        const double zMin = std::min(cells[(std::size_t)a].zMin, cells[(std::size_t)b].zMin);
        const double zMax = std::max(cells[(std::size_t)a].zMax, cells[(std::size_t)b].zMax);
        if (rMax - rMin > tau || zMax - zMin > tau) return false;
        cells[(std::size_t)b].parent = a;
        cells[(std::size_t)a].rMin = rMin;
        cells[(std::size_t)a].rMax = rMax;
        cells[(std::size_t)a].zMin = zMin;
        cells[(std::size_t)a].zMax = zMax;
        return true;
    };
    {
        bool changed = !cells.empty();
        for (int guard = 0; changed && guard < (int)cells.size() + 2; ++guard) {
            changed = false;
            for (int i = 0; i < (int)cells.size(); ++i)
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dz = -1; dz <= 1; ++dz) {
                        if (dr == 0 && dz == 0) continue;
                        const auto it = buckets.find(std::make_pair(cells[(std::size_t)i].kr + dr,
                                                                    cells[(std::size_t)i].kz + dz));
                        if (it != buckets.end() && uniteCell(i, it->second)) changed = true;
                    }
        }
    }
    struct Ring {
        double rho = 0, z = 0;
        int minV = 0;
    };
    std::vector<Ring> rings;
    std::vector<int> ringOf((std::size_t)maxV + 1, -1);
    {
        std::vector<std::vector<int>> groups(cells.size());
        for (int i = 0; i < (int)cells.size(); ++i) {
            auto& g = groups[(std::size_t)findCell(i)];
            g.insert(g.end(), cells[(std::size_t)i].verts.begin(), cells[(std::size_t)i].verts.end());
        }
        std::vector<std::pair<Ring, std::vector<int>>> made;
        for (auto& g : groups) {
            if ((int)g.size() < minStations(kParamsLine)) continue;
            Ring ring;
            ring.minV = g.front();
            double sr = 0, sz = 0;
            for (int v : g) {
                sr += rhoOf[(std::size_t)v];
                sz += zOf[(std::size_t)v];
                ring.minV = std::min(ring.minV, v);
            }
            ring.rho = sr / (double)g.size();
            ring.z = sz / (double)g.size();
            made.emplace_back(ring, std::move(g));
        }
        std::sort(made.begin(), made.end(), [](const auto& a, const auto& b) {
            if (a.first.z != b.first.z) return a.first.z < b.first.z;
            if (a.first.rho != b.first.rho) return a.first.rho < b.first.rho;
            return a.first.minV < b.first.minV;
        });
        for (int i = 0; i < (int)made.size(); ++i) {
            rings.push_back(made[(std::size_t)i].first);
            for (int v : made[(std::size_t)i].second) ringOf[(std::size_t)v] = i;
        }
    }
    const int nRings = (int)rings.size();
    struct Cand {
        int t = -1;
        int sheet = 0;
        int nR = 0;
        int ids[3] = {-1, -1, -1};
    };
    std::vector<Cand> cands;
    for (int t = 0; t < nTri; ++t) {
        if (!liveTri(t)) continue;
        Cand c;
        c.t = t;
        bool onAll = true;
        for (int k = 0; k < 3; ++k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            const int v = mv.tris[gt][k];
            const int r = (v >= 0 && v <= maxV) ? ringOf[(std::size_t)v] : -1;
            if (r < 0) {
                onAll = false;
                break;
            }
            bool dup = false;
            for (int i = 0; i < c.nR; ++i)
                if (c.ids[i] == r) dup = true;
            if (!dup) c.ids[c.nR++] = r;
        }
        if (!onAll || c.nR < 2) continue;
        c.sheet = sheetOf(mv, t, fr);
        if (c.sheet == 0) continue;
        cands.push_back(c);
    }
    if (cands.empty()) {
        scan.scanned = true;
        scans.push_back(scan);
        return;
    }
    struct Edge {
        int u = 0, v = 0, sheet = 0, minTri = 0;
        bool used = false;
    };
    std::vector<Edge> edges;
    std::map<std::tuple<int, int, int>, int> edgeKey;
    for (const Cand& c : cands)
        for (int a = 0; a < c.nR; ++a)
            for (int b = a + 1; b < c.nR; ++b) {
                const int u = std::min(c.ids[a], c.ids[b]), v = std::max(c.ids[a], c.ids[b]);
                const auto key = std::make_tuple(c.sheet, u, v);
                auto it = edgeKey.find(key);
                if (it == edgeKey.end()) {
                    Edge e;
                    e.u = u;
                    e.v = v;
                    e.sheet = c.sheet;
                    e.minTri = c.t;
                    edgeKey.emplace(key, (int)edges.size());
                    edges.push_back(e);
                } else if (c.t < edges[(std::size_t)it->second].minTri) {
                    edges[(std::size_t)it->second].minTri = c.t;
                }
            }
    std::vector<std::vector<int>> edgeAt((std::size_t)nRings);
    for (int i = 0; i < (int)edges.size(); ++i) {
        edgeAt[(std::size_t)edges[(std::size_t)i].u].push_back(i);
        edgeAt[(std::size_t)edges[(std::size_t)i].v].push_back(i);
    }
    auto concyclic = [&](const std::vector<int>& ids) {
        if ((int)ids.size() <= kParamsCircle) return true;
        std::vector<ProfPt> pts;
        pts.reserve(ids.size());
        for (int id : ids) {
            ProfPt p;
            p.rho = rings[(std::size_t)id].rho;
            p.z = rings[(std::size_t)id].z;
            pts.push_back(p);
        }
        double Rmaj = 0, z0 = 0, Rmin = 0;
        if (!fitProfileCircle(pts, tau, Rmaj, z0, Rmin) || !std::isfinite(Rmin)) return false;
        for (const ProfPt& p : pts)
            if (!(std::fabs(std::hypot(p.rho - Rmaj, p.z - z0) - Rmin) <= tau)) return false;
        return true;
    };
    std::vector<int> edgeOrder(edges.size());
    for (int i = 0; i < (int)edges.size(); ++i) edgeOrder[(std::size_t)i] = i;
    std::sort(edgeOrder.begin(), edgeOrder.end(), [&](int a, int b) {
        const Edge& A = edges[(std::size_t)a];
        const Edge& B = edges[(std::size_t)b];
        if (A.minTri != B.minTri) return A.minTri < B.minTri;
        if (A.u != B.u) return A.u < B.u;
        if (A.v != B.v) return A.v < B.v;
        return A.sheet < B.sheet;
    });
    std::set<std::vector<int>> tried;
    // D-fc-4.10 F1, refined by D-fc-4.11: the candidate sets of chains, and the
    // bands, that failed the GEOMETRIC certificate; no later band that overlaps
    // one may pass (a passing subset of a failing band is a refusal).
    std::vector<char> failedTri((std::size_t)nTri, 0);
    for (int ei0 : edgeOrder) {
        if (edges[(std::size_t)ei0].used) continue;
        const int sheet = edges[(std::size_t)ei0].sheet;
        std::vector<int> chain{edges[(std::size_t)ei0].u, edges[(std::size_t)ei0].v};
        std::vector<char> inChain((std::size_t)nRings, 0);
        inChain[(std::size_t)chain[0]] = inChain[(std::size_t)chain[1]] = 1;
        auto grow = [&](bool front) {
            while (true) {
                const int cur = front ? chain.front() : chain.back();
                int bestE = -1, bestN = -1, bestMin = std::numeric_limits<int>::max();
                for (int ei : edgeAt[(std::size_t)cur]) {
                    const Edge& e = edges[(std::size_t)ei];
                    if (e.used || e.sheet != sheet) continue;
                    const int nxt = (e.u == cur) ? e.v : e.u;
                    if (inChain[(std::size_t)nxt]) continue;
                    std::vector<int> trial = chain;
                    if (front) trial.insert(trial.begin(), nxt);
                    else trial.push_back(nxt);
                    if (!concyclic(trial)) continue;
                    if (e.minTri < bestMin || (e.minTri == bestMin && nxt < bestN)) {
                        bestMin = e.minTri;
                        bestE = ei;
                        bestN = nxt;
                    }
                }
                if (bestE < 0) break;
                inChain[(std::size_t)bestN] = 1;
                if (front) chain.insert(chain.begin(), bestN);
                else chain.push_back(bestN);
            }
        };
        grow(false);
        grow(true);
        {
            std::vector<int> key = chain;
            std::sort(key.begin(), key.end());
            key.push_back(sheet);
            if (!tried.insert(key).second) continue;
        }
        std::vector<int> bandTris;
        for (const Cand& c : cands) {
            if (c.sheet != sheet || !liveTri(c.t)) continue;
            bool inside = true;
            for (int i = 0; i < c.nR && inside; ++i)
                if (!inChain[(std::size_t)c.ids[i]]) inside = false;
            if (inside) bandTris.push_back(c.t);
        }
        if (bandTris.empty()) continue;
        // D-fc-4.8 (a): the chain's candidate set is the band's superset.
        const std::vector<int> chainCands = bandTris;
        // An atom joins only whole: one that is not wholly on the chain leaves it.
        {
            std::map<int, int> hit;
            for (int t : bandTris)
                if (atomOf[(std::size_t)t] >= 0) ++hit[atomOf[(std::size_t)t]];
            std::vector<char> drop(work.accepted.size(), 0);
            bool any = false;
            for (const auto& h : hit)
                if (h.second != (int)work.accepted[(std::size_t)h.first].tris.size()) {
                    drop[(std::size_t)h.first] = 1;
                    any = true;
                }
            if (any) {
                std::vector<int> keep;
                for (int t : bandTris) {
                    const int a = atomOf[(std::size_t)t];
                    if (a < 0 || !drop[(std::size_t)a]) keep.push_back(t);
                }
                bandTris.swap(keep);
            }
        }
        if (bandTris.empty()) continue;
        // D-fc-4.9 P2 (i): the chain is certified with its donor family. The
        // family axis is kept when it already certifies the chain's candidate set
        // (every donor vertex is within tau on it by construction); otherwise one
        // LS axis over the family's donor vertices plus the chain's vertices
        // decides, and the chain is accepted only when every donor vertex and
        // every chain vertex is within tau on it. The chain's bands are certified
        // and built on the axis the chain certified on. A band that is a strict
        // subset of a refused chain set is refused (D-fc-4.8 (a)).
        bool chainOk = false;
        Frame cfr = fr;
        Donor jointDonor;
        bool joint = false;
        bool chainGeomFail = false;  // D-fc-4.11: the candidate set failed geometrically
        TorusCert onFamilyCert;
        JointTorusCert jointCert;
        bool jointRan = false;
        {
            Band sup;
            sup.in.assign((std::size_t)nTri, 0);
            sup.tris = chainCands;
            sup.sheet = sheet;
            for (int t : chainCands) sup.in[(std::size_t)t] = 1;
            collectBand(mv, fr, wm.g2l, sup);
            const TorusCert onFamily = certifyTorus(sup.pts, tau);
            onFamilyCert = onFamily;
            chainOk = onFamily.ok;
            // certifyTorus needs kTorusProfPts distinct profile points on any axis:
            // a chain with fewer distinct vertices cannot certify jointly either.
            if (!chainOk && (int)sup.pts.size() >= kTorusProfPts) {
                std::vector<gp_XYZ> cpts;
                bandPoints(sup, cpts);
                const JointTorusCert jc = jointTorusCertify(donor.famPts, cpts, fr.origin, fr.axis,
                                                            fr.e1, tau, wm.maxIter);
                jointCert = jc;
                jointRan = true;
                if (diagOn()) {
                    double dMax = 0;
                    for (double m : jc.family.maxDev) dMax = std::max(dMax, m);
                    std::fprintf(stderr,
                                 "DIAG_REL_TORUS_JOINT seed=%d nCand=%zu nV=%zu donor=%d "
                                 "members=%zu familyWhy=%s familyDevTau=%.4f fitted=%d "
                                 "tilt=%.4g donorDevTau=%.4f chainDevTau=%.4f radii=%.8g/%.8g "
                                 "chainWhy=%s ok=%d\n",
                                 minTriOf(chainCands), chainCands.size(), sup.pts.size(),
                                 donor.idx, donor.famPts.size(), torusWhy(onFamily),
                                 onFamily.determinate ? onFamily.maxDev / tau : -1.0,
                                 jc.fitted ? 1 : 0, jc.tilt, dMax / tau,
                                 jc.chain.determinate ? jc.chain.maxDev / tau : -1.0,
                                 jc.chain.Rmaj, jc.chain.Rmin, torusWhy(jc.chain), jc.ok ? 1 : 0);
                }
                if (jc.ok) {
                    chainOk = true;
                    joint = true;
                    cfr = makeFrame(jc.family.origin, jc.family.axis, fr.e1);
                    jointDonor = donor;
                    jointDonor.fr = cfr;
                    jointDonor.radius = jc.family.radius.front();
                    jointDonor.maxDev = jc.family.maxDev.front();
                }
            }
            chainGeomFail = !chainOk && torusChainFailsGeometric(
                                            onFamily, jointRan ? &jointCert : nullptr, tau);
        }
        Donor& cdonor = joint ? jointDonor : donor;
        // D-fc-4.11 (F1 refined): the chain consumes its ring-graph edges after its
        // bands, when one of them was claimed (Line mode's partition) or when its
        // candidate set or one of its bands failed the geometric certificate. A
        // chain refused only for over-determination leaves its edges to later,
        // larger windows, which must themselves pass every check.
        bool claimedChain = false;
        bool bandGeomFail = false;
        int nBandClaimed = 0, nBandGeom = 0, nBandOverdet = 0, nBandOther = 0;
        // Each edge-connected component is one band.
        auto components = [&](const std::vector<int>& tris) {
            std::vector<char> in((std::size_t)nTri, 0), seen((std::size_t)nTri, 0);
            for (int t : tris) in[(std::size_t)t] = 1;
            std::vector<std::vector<int>> comps;
            for (int t : tris) {
                if (seen[(std::size_t)t]) continue;
                std::vector<int> comp, stack{t};
                seen[(std::size_t)t] = 1;
                while (!stack.empty()) {
                    const int c = stack.back();
                    stack.pop_back();
                    comp.push_back(c);
                    for (int n : wm.adj.nbr[(std::size_t)c])
                        if (n >= 0 && n < nTri && in[(std::size_t)n] && !seen[(std::size_t)n]) {
                            seen[(std::size_t)n] = 1;
                            stack.push_back(n);
                        }
                }
                std::sort(comp.begin(), comp.end());
                comps.push_back(std::move(comp));
            }
            return comps;
        };
        std::vector<std::vector<int>> queue = components(bandTris);
        std::sort(queue.begin(), queue.end(),
                  [](const std::vector<int>& a, const std::vector<int>& b) { return a.front() < b.front(); });
        for (std::size_t wi = 0; wi < queue.size(); ++wi) {
            std::vector<int> comp = queue[wi];
            Band band;
            band.in.assign((std::size_t)nTri, 0);
            band.tris = comp;
            band.sheet = sheet;
            band.seed = comp.front();
            for (int t : comp) band.in[(std::size_t)t] = 1;
            collectBand(mv, cfr, wm.g2l, band);
            if (band.pts.empty()) continue;
            const SurfaceFit fit = fitSurface(band.pts, tau);
            const ClosedM cl = closedOf(mv, cfr, band.tris, band.verts, wm.thetaBin);
            std::vector<double> zs, psis, rhos;
            for (const ProfPt& pp : band.pts) {
                zs.push_back(pp.z);
                psis.push_back(pp.psi);
                rhos.push_back(pp.rho);
            }
            const int levels = countLevels(zs, tau);
            const int stations = countProfPts(band.pts, tau);
            std::vector<int> atoms;
            for (int t : comp)
                if (atomOf[(std::size_t)t] >= 0) atoms.push_back(atomOf[(std::size_t)t]);
            std::sort(atoms.begin(), atoms.end());
            atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
            // geom: the band itself failed the geometric certificate (D-fc-4.11);
            // its triangles become a failing band. overdet: refused only for
            // over-determination.
            auto refuse = [&](const char* why, const TorusCert* c, bool geom = false,
                              bool overdet = false) {
                ++refused;
                if (geom) {
                    bandGeomFail = true;
                    ++nBandGeom;
                    for (int t : band.tris) failedTri[(std::size_t)t] = 1;
                } else if (overdet) {
                    ++nBandOverdet;
                } else {
                    ++nBandOther;
                }
                DiagRec d = baseDiag("torus", cdonor, band, tau);
                d.stations = stations;
                d.levels = levels;
                d.params = kParamsCircle;
                d.maxDev = (c && c->determinate) ? c->maxDev : fit.circMax;
                d.maxDevTau = tau > 0 ? d.maxDev / tau : 0;
                d.lineResTau = tau > 0 ? fit.orthMax / tau : 0;
                d.span = cl.span;
                d.spanTol = (cl.nBands > 0) ? 1.5 * kTwoPi / cl.nBands : 0;
                d.r0 = (c && c->Rmaj > 0) ? c->Rmaj : fit.Rmaj;
                d.r1 = (c && c->Rmin > 0) ? c->Rmin : fit.Rmin;
                d.why = why;
                printBand(d);
            };
            const TorusCert cert = certifyTorus(band.pts, tau);
            // Line-first is a geometric refusal only on determinate evidence: a
            // band below certifyTorus's own counts is under-determined first
            // (D-fc-1.6.2 certificate order).
            if (fit.line && fit.orthMax <= tau) {
                refuse("line", nullptr, cert.determinate, !cert.determinate);
                continue;
            }
            if (!cert.ok) {
                // An atom with a vertex off the band's circle leaves whole; what
                // remains is certified again as its own components.
                bool dropped = false;
                if (cert.determinate && cert.Rmin > 0.0 &&
                    std::strcmp(cert.why, "vertex-off-surface") == 0 && !atoms.empty()) {
                    std::vector<char> off(work.accepted.size(), 0);
                    for (int t : comp) {
                        const int a = atomOf[(std::size_t)t];
                        if (a < 0 || off[(std::size_t)a]) continue;
                        for (int k = 0; k < 3; ++k) {
                            const ProfPt p = profOfPoint(triVert(mv, t, k), cfr, 0);
                            if (std::fabs(std::hypot(p.rho - cert.Rmaj, p.z - cert.z0) - cert.Rmin) > tau) {
                                off[(std::size_t)a] = 1;
                                break;
                            }
                        }
                    }
                    std::vector<int> keep;
                    for (int t : comp) {
                        const int a = atomOf[(std::size_t)t];
                        if (a >= 0 && off[(std::size_t)a]) {
                            dropped = true;
                            continue;
                        }
                        keep.push_back(t);
                    }
                    if (dropped) {
                        std::vector<std::vector<int>> sub = components(keep);
                        std::sort(sub.begin(), sub.end(), [](const std::vector<int>& a,
                                                             const std::vector<int>& b) {
                            return a.front() < b.front();
                        });
                        for (auto& s : sub) queue.push_back(std::move(s));
                    }
                }
                if (!dropped)
                    refuse(torusWhy(cert), &cert, torusRefusalGeometric(cert),
                           !torusRefusalGeometric(cert));
                continue;
            }
            // Leave-one-region-out: the band without any one atom still meets
            // certifyTorus's own counts.
            {
                std::vector<std::array<ProfPt, 3>> triPts;
                std::vector<int> triAtom;
                triPts.reserve(comp.size());
                triAtom.reserve(comp.size());
                for (int t : comp) {
                    std::array<ProfPt, 3> tp;
                    for (int k = 0; k < 3; ++k) tp[(std::size_t)k] = profOfPoint(triVert(mv, t, k), cfr, localId(wm, t, k));
                    triPts.push_back(tp);
                    const int a = atomOf[(std::size_t)t];
                    int slot = -1;
                    if (a >= 0)
                        slot = (int)(std::lower_bound(atoms.begin(), atoms.end(), a) - atoms.begin());
                    triAtom.push_back(slot);
                }
                if (!leaveOneOutDeterminate(triPts, triAtom, (int)atoms.size(), tau)) {
                    refuse("underdetermined", &cert, false, true);
                    continue;
                }
            }
            if (band.tris.size() < chainCands.size() && !chainOk) {
                refuse("vertex-off", &cert, false, !chainGeomFail);
                continue;
            }
            if (torusWindowOverlapsRefused(band.tris, failedTri)) {
                refuse("vertex-off", &cert);
                continue;
            }
            if (sameSurfaceTorus(wm, work, cfr, band, fit, atomDead)) {
                refuse("same-surface", &cert);
                continue;
            }
            // D-fc-4.11 (a): a band that holds a triangle of a torus face the
            // builder could not ship is withheld. It is gathered as this claim
            // (its edges, triangles and atoms leave later windows, as they would
            // after the claim) but nothing is emitted, trimmed or updated: its
            // triangles and atoms keep their pre-T state.
            if (torusWindowOverlapsRefused(band.tris, wm.withheldTri)) {
                for (int a : atoms) {
                    atomDead[(std::size_t)a] = 1;
                    wm.atomHeld[(std::size_t)a] = 1;
                }
                for (int t : band.tris) wm.freeTri[(std::size_t)t] = 0;
                claimedChain = true;
                ++nBandOther;
                if (diagOn())
                    std::fprintf(stderr,
                                 "DIAG_REL_TORUS_HELD seed=%d nTri=%zu atoms=%zu chainRings=%zu "
                                 "donor=%d maxDevTau=%.4f radii=%.6f/%.6f\n",
                                 band.seed, band.tris.size(), atoms.size(), chain.size(),
                                 donor.idx, cert.maxDev / tau, cert.Rmaj, cert.Rmin);
                continue;
            }
            int straddle = countStraddle(work, band.tris, nTri);
            if (!emitTorus(wm, work, cdonor, band, cert, cl, straddle)) {
                refuse("degenerate", &cert);
                continue;
            }
            {
                Region& tor = work.accepted.back();
                if (tor.torusPending >= 0 && (std::size_t)tor.torusPending < work.torusPending.size()) {
                    TorusPendingReclaim& pack = work.torusPending[(std::size_t)tor.torusPending];
                    for (int a : atoms) {
                        pack.stashedCylinders.push_back(work.accepted[(std::size_t)a]);
                        pack.reclaimAccIdx.push_back(a);
                    }
                }
            }
            for (int a : atoms) atomDead[(std::size_t)a] = 1;
            for (int t : band.tris) wm.freeTri[(std::size_t)t] = 0;
            claimedChain = true;
            ++nBandClaimed;
            // F8: one seam azimuth per revolution axis. The coaxial cylinders the
            // torus touches take the donor frame it took.
            shareSeamFrame(wm, work, cdonor, band);
            updateDonor(work, cdonor);
            if (joint) donor.updated = donor.updated || jointDonor.updated;
            ++claimed;
            DiagRec d = baseDiag("torus", cdonor, band, tau);
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
            double ang0 = 0, ang1 = 0;
            endTangency(wm, cfr, band, fit, ang0, ang1);
            d.tan0 = ang0;
            d.tan1 = ang1;
            d.straddle = straddle;
            d.why = "accept";
            printBand(d);
            if (diagOn())
                std::fprintf(stderr,
                             "DIAG_REL_TORUS_BAND seed=%d nTri=%zu atoms=%zu chainRings=%zu donor=%d "
                             "donorR=%.6f donorClosed=%d\n",
                             band.seed, band.tris.size(), atoms.size(), chain.size(), donor.idx,
                             donor.radius, work.accepted[(std::size_t)donor.idx].closed360 ? 1 : 0);
        }
        const bool consumes = claimedChain || chainGeomFail || bandGeomFail;
        if (consumes)
            for (Edge& e : edges)
                if (e.sheet == sheet && inChain[(std::size_t)e.u] && inChain[(std::size_t)e.v])
                    e.used = true;
        if (chainGeomFail)
            for (int t : chainCands) failedTri[(std::size_t)t] = 1;
        if (diagOn()) {
            double donorDev = 0;
            for (double m : jointCert.family.maxDev) donorDev = std::max(donorDev, m);
            std::fprintf(stderr,
                         "DIAG_REL_TORUS_CHAIN seed=%d nCand=%zu rings=%zu donor=%d familyWhy=%s "
                         "familyDevTau=%.4f joint=%d jointOk=%d donorDevTau=%.4f chainDevTau=%.4f "
                         "jointWhy=%s chainOk=%d chainGeomFail=%d bands=%zu claimed=%d geom=%d "
                         "overdet=%d other=%d consumed=%d failingSuperset=%d\n",
                         minTriOf(chainCands), chainCands.size(), chain.size(), donor.idx,
                         torusWhy(onFamilyCert),
                         onFamilyCert.determinate ? onFamilyCert.maxDev / tau : -1.0,
                         jointRan ? 1 : 0, jointCert.ok ? 1 : 0, jointRan ? donorDev / tau : -1.0,
                         (jointRan && jointCert.chain.determinate) ? jointCert.chain.maxDev / tau
                                                                   : -1.0,
                         jointRan ? torusWhy(jointCert.chain) : "-", chainOk ? 1 : 0,
                         chainGeomFail ? 1 : 0, queue.size(), nBandClaimed, nBandGeom,
                         nBandOverdet, nBandOther, consumes ? 1 : 0, chainGeomFail ? 1 : 0);
        }
    }
    scan.scanned = true;
    scans.push_back(scan);
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
            // D-fc-4.6 F8: a cylinder donor whose stored frame already certifies its
            // own vertices within tau (the frame H0 refit and wrote, which its
            // coaxial cones carry) keeps it, so a torus on it is exactly coaxial
            // with them. Line mode is unchanged.
            if (mode == RelRevolveMode::Torus && r.type == SurfType::Cylinder) {
                const Frame sf = makeFrame(r.ax.Location().XYZ(), r.ax.Direction().XYZ(),
                                           r.ax.XDirection().XYZ());
                double sMax = 0;
                for (int g : gids) {
                    sMax = std::max(sMax, std::fabs(profOfPoint(mv.pts[g], sf, 0).rho - r.radius));
                    if (!(sMax <= tau)) break;
                }
                if (sMax <= tau && r.radius > tau) {
                    fr = sf;
                    radius = r.radius;
                }
            }
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
            d.gids = gids;
            donors.push_back(d);
        }
        std::sort(donors.begin(), donors.end(), [](const Donor& a, const Donor& b) {
            if (a.minTri != b.minTri) return a.minTri < b.minTri;
            return a.idx < b.idx;
        });

        int claimed = 0, refused = 0;
        std::vector<AxisScan> axisScans;

        if (mode != RelRevolveMode::Torus) {
            // D-fc-5.1 / D-fc-4.8 (b): the coaxial families of the Line donors,
            // formed once in scan order. They carry the joint certificate of a
            // coaxial band; the cone gather keeps each donor's own frame.
            LineFamilies lf;
            lf.famOfAcc.assign(work.accepted.size(), -1);
            lf.tol = &tol;
            lf.sp = &p;
            {
                std::vector<AxisDonor> ad(donors.size());
                std::vector<int> scanOrder;
                for (std::size_t i = 0; i < donors.size(); ++i) {
                    const Donor& d = donors[i];
                    ad[i].origin = d.fr.origin;
                    ad[i].axis = d.fr.axis;
                    ad[i].radius = d.radius;
                    ad[i].length = d.length;
                    for (int g : d.gids) ad[i].pts.push_back(mv.pts[g]);
                    scanOrder.push_back((int)i);
                }
                for (const CoaxialFamily& f : coaxialFamilies(ad, scanOrder, tau, wm.maxIter)) {
                    CoaxFamily cf;
                    cf.origin = f.fit.origin;
                    cf.axis = f.fit.axis;
                    for (int m : f.members) {
                        cf.groups.push_back(ad[(std::size_t)m].pts);
                        const int acc = donors[(std::size_t)m].idx;
                        if (acc >= 0 && acc < (int)lf.famOfAcc.size())
                            lf.famOfAcc[(std::size_t)acc] = (int)lf.fams.size();
                    }
                    cf.nDonors = cf.groups.size();
                    lf.fams.push_back(std::move(cf));
                }
            }
            for (Donor& donor : donors)
                gatherLines(wm, work, donor, axisScans, lf, claimed, refused);
            // One family, one axis: a band placed before a later band moved the
            // family axis takes the final one. The final joint certificate
            // already holds every band's vertices within tau of its radius on it.
            for (CoaxFamily& fam : lf.fams) {
                if (!fam.moved) continue;
                for (std::size_t k = 0; k < fam.bandRegion.size(); ++k) {
                    Region& r = work.accepted[(std::size_t)fam.bandRegion[k]];
                    const std::vector<int> tris = r.tris;
                    placeCoaxial(mv, r, tris, fam.groups[fam.nDonors + k], fam.origin, fam.axis,
                                 wm.thetaBin);
                }
            }
            // D-fc-5.2: section fillets need the circles certified above (the
            // donors and the coaxial bands), so they follow the Line gather.
            std::vector<SectionCircle> circles;
            for (const Donor& d : donors) {
                SectionCircle c;
                c.origin = d.fr.origin;
                c.axis = d.fr.axis;
                c.radius = d.radius;
                c.tol = anchorTol(tau, d.maxDev);
                c.length = d.length;
                c.acc = d.idx;
                circles.push_back(c);
            }
            for (const CoaxFamily& fam : lf.fams)
                for (int ri : fam.bandRegion) {
                    const Region& r = work.accepted[(std::size_t)ri];
                    SectionCircle c;
                    c.origin = fam.origin;
                    c.axis = fam.axis;
                    c.radius = r.radius;
                    c.tol = anchorTol(tau, r.maxVertexDev);
                    c.length = std::max(r.vMax - r.vMin, tau);
                    c.acc = ri;
                    circles.push_back(c);
                }
            claimSectionFillets(mv, p, tol, work, circles);
        } else {
            // D-fc-4.6 F2/F3. Atoms are candidate units of every scan but their
            // own; an atom a band absorbed leaves later scans and stops donating.
            // Certified closed donors (closed-360 cylinders, closed cones) scan
            // first, so their axes meet the atoms before the atoms' own axes do.
            std::vector<int> atomOf((std::size_t)nTri, -1);
            for (int i = 0; i < (int)work.accepted.size(); ++i) {
                const Region& r = work.accepted[(std::size_t)i];
                if (!torusAtom(r)) continue;
                for (int t : r.tris)
                    if (t >= 0 && t < nTri) atomOf[(std::size_t)t] = i;
            }
            std::vector<char> atomDead(work.accepted.size(), 0);
            wm.atomHeld.assign(work.accepted.size(), 0);
            wm.withheldTri.assign((std::size_t)nTri, 0);
            for (int t : p.torusWithheldTris)
                if (t >= 0 && t < nTri) wm.withheldTri[(std::size_t)t] = 1;
            std::vector<std::size_t> order(donors.size());
            for (std::size_t i = 0; i < donors.size(); ++i) order[i] = i;
            std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                const bool ca = work.accepted[(std::size_t)donors[a].idx].closed360;
                const bool cb = work.accepted[(std::size_t)donors[b].idx].closed360;
                return ca && !cb;
            });
            // D-fc-4.8 (b): coaxial donors whose certified vertices jointly certify
            // one axis at tau are one family, scanned once on that axis.
            std::vector<AxisDonor> ad(donors.size());
            std::vector<int> scanOrder;
            for (std::size_t i = 0; i < donors.size(); ++i) {
                const Donor& d = donors[i];
                ad[i].origin = d.fr.origin;
                ad[i].axis = d.fr.axis;
                ad[i].radius = d.radius;
                ad[i].length = d.length;
                for (int g : d.gids) ad[i].pts.push_back(mv.pts[g]);
            }
            for (std::size_t oi : order) scanOrder.push_back((int)oi);
            const std::vector<CoaxialFamily> fams =
                coaxialFamilies(ad, scanOrder, tau, wm.maxIter);
            for (const CoaxialFamily& fam : fams) {
                Donor donor = donors[(std::size_t)fam.members.front()];
                for (int m : fam.members) {
                    donor.members.push_back(donors[(std::size_t)m].idx);
                    donor.famPts.push_back(ad[(std::size_t)m].pts);
                }
                if (fam.members.size() > 1) {
                    donor.fr = makeFrame(fam.fit.origin, fam.fit.axis, donor.fr.e1);
                    donor.radius = fam.fit.radius.front();
                    donor.maxDev = fam.fit.maxDev.front();
                    double zMin = 1e300, zMax = -1e300, devMax = 0;
                    for (std::size_t k = 0; k < fam.members.size(); ++k) {
                        devMax = std::max(devMax, fam.fit.maxDev[k]);
                        for (const gp_XYZ& pt : ad[(std::size_t)fam.members[k]].pts) {
                            const double z = (pt - donor.fr.origin).Dot(donor.fr.axis);
                            zMin = std::min(zMin, z);
                            zMax = std::max(zMax, z);
                        }
                    }
                    donor.length = (zMax > zMin) ? (zMax - zMin) : tau;
                    if (diagOn()) {
                        std::fprintf(stderr,
                                     "DIAG_REL_TORUS_FAMILY lead=%d members=%zu refit=%d "
                                     "maxDevTau=%.4f axis=%.9g,%.9g,%.9g ids=",
                                     donor.idx, fam.members.size(), fam.refit ? 1 : 0,
                                     devMax / tau, donor.fr.axis.X(), donor.fr.axis.Y(),
                                     donor.fr.axis.Z());
                        for (std::size_t k = 0; k < fam.members.size(); ++k)
                            std::fprintf(stderr, "%s%d(R%.4f)", k ? "," : "",
                                         donors[(std::size_t)fam.members[k]].idx,
                                         fam.fit.radius[k]);
                        std::fprintf(stderr, "\n");
                    }
                }
                if (atomDead[(std::size_t)donor.idx]) continue;
                gatherTori(wm, work, donor, axisScans, atomOf, atomDead, claimed, refused);
            }
            // An atom a withheld band would have absorbed stays (D-fc-4.11 (a)).
            for (std::size_t i = 0; i < atomDead.size() && i < wm.atomHeld.size(); ++i)
                if (wm.atomHeld[i]) atomDead[i] = 0;
            int nDead = 0;
            for (char c : atomDead) nDead += c ? 1 : 0;
            if (nDead > 0) {
                const int nAcc = (int)work.accepted.size();
                std::vector<int> remap((std::size_t)nAcc, -1);
                std::vector<Region> kept;
                kept.reserve((std::size_t)(nAcc - nDead));
                for (int i = 0; i < nAcc; ++i) {
                    if (i < (int)atomDead.size() && atomDead[(std::size_t)i]) continue;
                    remap[(std::size_t)i] = (int)kept.size();
                    kept.push_back(std::move(work.accepted[(std::size_t)i]));
                }
                for (Region& r : kept)
                    if (r.origin == Origin::TorusBlend && r.filletNbrA >= 0 && r.filletNbrA < nAcc)
                        r.filletNbrA = remap[(std::size_t)r.filletNbrA];
                work.accepted.swap(kept);
            }
        }
        printSum(modeName, (int)donors.size(), claimed, refused + refusedDonors);
        (void)axisGap;
        return true;
    } catch (...) {
        return true;
    }
}

double sectionLineResid(const std::vector<gp_XYZ>& pts, const gp_XYZ& axisIn) {
    if (pts.size() < 2) return 0.0;
    gp_XYZ axis = axisIn;
    const double am = axis.Modulus();
    if (!(am > gp::Resolution())) return 0.0;
    axis /= am;
    const Frame f = makeFrame(gp_XYZ(0, 0, 0), axis, gp_XYZ(1, 0, 0));
    double mx = 0, my = 0;
    for (const gp_XYZ& p : pts) {
        mx += p.Dot(f.e1);
        my += p.Dot(f.e2);
    }
    mx /= (double)pts.size();
    my /= (double)pts.size();
    double sxx = 0, sxy = 0, syy = 0;
    for (const gp_XYZ& p : pts) {
        const double dx = p.Dot(f.e1) - mx, dy = p.Dot(f.e2) - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    // Principal direction by its angle (no eigenvalue cancellation); the line
    // normal is perpendicular to it.
    const double th = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    const double nx = -std::sin(th), ny = std::cos(th);
    double m = 0;
    for (const gp_XYZ& p : pts)
        m = std::max(m, std::fabs((p.Dot(f.e1) - mx) * nx + (p.Dot(f.e2) - my) * ny));
    return m;
}

CoaxialCert certifyCoaxialBand(const std::vector<std::vector<gp_XYZ>>& family,
                               const std::vector<gp_XYZ>& band, const gp_XYZ& origin,
                               const gp_XYZ& axisIn, double tau, int maxIters) {
    CoaxialCert c;
    if (!(tau > 0.0) || band.empty() || family.empty()) return c;
    gp_XYZ axis = axisIn;
    const double am = axis.Modulus();
    if (!(am > gp::Resolution())) return c;
    axis /= am;
    const Frame fr = makeFrame(origin, axis, gp_XYZ(1, 0, 0));
    std::vector<double> zs, psis, rhos;
    for (const gp_XYZ& p : band) {
        const ProfPt pp = profOfPoint(p, fr, 0);
        zs.push_back(pp.z);
        psis.push_back(pp.psi);
        rhos.push_back(pp.rho);
    }
    c.levels = countLevels(zs, tau);
    c.stations = countColumns(psis, rhos, tau);
    if (c.levels < 2 || c.stations < minStations(kParamsCircle)) {
        c.why = "underdetermined";
        return c;
    }
    c.lineRes = sectionLineResid(band, axis);
    if (c.lineRes <= tau) {
        c.why = "line";
        return c;
    }
    std::vector<std::vector<gp_XYZ>> groups = family;
    groups.push_back(band);
    c.fit = jointAxisAbout(groups, origin, axis);
    if (jointAxisCertifies(c.fit, tau)) {
        c.ok = true;
        c.why = "accept";
        return c;
    }
    JointAxisFit g;
    if (jointAxisFit(groups, origin, axis, maxIters, g)) {
        c.fit = g;
        c.refit = true;
        if (jointAxisCertifies(g, tau)) {
            c.ok = true;
            c.why = "accept";
            return c;
        }
    }
    c.why = "vertex-off";
    return c;
}

PartialConeCert certifyPartialCone(const std::vector<std::vector<gp_XYZ>>& family,
                                   const std::vector<gp_XYZ>& band, const gp_XYZ& origin,
                                   const gp_XYZ& axisIn, double tau) {
    PartialConeCert c;
    if (!(tau > 0.0) || band.empty() || family.empty()) return c;
    gp_XYZ axis = axisIn;
    const double am = axis.Modulus();
    if (!(am > gp::Resolution())) return c;
    axis /= am;
    const Frame fr = makeFrame(origin, axis, gp_XYZ(1, 0, 0));
    std::vector<ProfPt> pts;
    std::vector<double> zs, psis, rhos;
    for (const gp_XYZ& p : band) {
        const ProfPt pp = profOfPoint(p, fr, 0);
        pts.push_back(pp);
        zs.push_back(pp.z);
        psis.push_back(pp.psi);
        rhos.push_back(pp.rho);
    }
    // The profile line (total least squares in (rho, z)).
    double sr = 0, sz = 0;
    for (const ProfPt& pp : pts) {
        sr += pp.rho;
        sz += pp.z;
    }
    c.rm = sr / (double)pts.size();
    c.zm = sz / (double)pts.size();
    double srr = 0, srz = 0, szz = 0;
    for (const ProfPt& pp : pts) {
        const double dr = pp.rho - c.rm, dz = pp.z - c.zm;
        srr += dr * dr;
        srz += dr * dz;
        szz += dz * dz;
    }
    const double th = 0.5 * std::atan2(2.0 * srz, srr - szz);
    const double ur = std::cos(th), uz = std::sin(th);  // the line's direction
    // Stations of a line profile: the band's positions ALONG the line merged at
    // tau (countLevels on the line parameter); scatter across the line is the
    // residual, not a station.
    std::vector<double> ts;
    for (const ProfPt& pp : pts) ts.push_back((pp.rho - c.rm) * ur + (pp.z - c.zm) * uz);
    c.levels = countLevels(zs, tau);
    c.stations = countLevels(ts, tau);
    c.cols = countColumns(psis, rhos, tau);
    if (c.levels < 2 || c.stations < minStations(kParamsLine) ||
        c.cols < minStations(kParamsCircle)) {
        c.why = "underdetermined";
        return c;
    }
    // Line first: the vertices on one plane (least squares) within tau.
    {
        gp_XYZ m(0, 0, 0);
        for (const gp_XYZ& p : band) m += p;
        m /= (double)band.size();
        std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
        for (const gp_XYZ& p : band) {
            const double d[3] = {p.X() - m.X(), p.Y() - m.Y(), p.Z() - m.Z()};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) cov[(std::size_t)i][(std::size_t)j] += d[i] * d[j];
        }
        std::array<double, 3> ev{};
        std::array<std::array<double, 3>, 3> evec{};
        if (!jacobiEigenSymmetric3(cov, ev, evec)) return c;
        const gp_XYZ n(evec[0][0], evec[0][1], evec[0][2]);
        for (const gp_XYZ& p : band) c.planeRes = std::max(c.planeRes, std::fabs((p - m).Dot(n)));
        if (c.planeRes <= tau) {
            c.why = "line";
            return c;
        }
    }
    double rLo = 1e300, rHi = -1e300;
    c.zLo = 1e300;
    c.zHi = -1e300;
    for (const ProfPt& pp : pts) {
        c.zLo = std::min(c.zLo, pp.z);
        c.zHi = std::max(c.zHi, pp.z);
        rLo = std::min(rLo, pp.rho);
        rHi = std::max(rHi, pp.rho);
    }
    if (c.zHi - c.zLo <= tau) {
        c.why = "line";
        return c;
    }
    if (rHi - rLo <= tau || !(std::fabs(uz) > 0.0)) {
        c.why = "same-surface";
        return c;
    }
    c.k = ur / uz;
    for (const ProfPt& pp : pts)
        c.lineRes = std::max(c.lineRes, std::fabs(-(pp.rho - c.rm) * uz + (pp.z - c.zm) * ur));
    if (c.lineRes > tau) {
        c.why = "vertex-off";
        return c;
    }
    const JointAxisFit f = jointAxisAbout(family, origin, axis);
    for (double d : f.maxDev) c.famDev = std::max(c.famDev, d);
    if (!jointAxisCertifies(f, tau)) {
        c.why = "vertex-off";
        return c;
    }
    c.ok = true;
    c.why = "accept";
    return c;
}

bool coaxialSameRadius(const JointAxisFit& fit, std::size_t nFamily, double tau) {
    if (fit.radius.size() <= nFamily) return false;
    const double R = fit.radius.back();
    for (std::size_t g = 0; g < nFamily; ++g)
        if (std::fabs(fit.radius[g] - R) <= tau) return true;
    return false;
}

}  // namespace refit
}  // namespace stl2step
