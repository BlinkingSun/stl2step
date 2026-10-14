// Relational extrude (W1.1, D-fc-1.3 … D-fc-1.7, D-fc-1.15).
// Claims a band whose section perpendicular to one axis is one circle at τ:
// plane|plane fillets, U-turn hems, and closed pegs. Nothing else.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_profile_cert.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

bool diagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

gp_XYZ cornerPt(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.pts[mv.tris[gt][k]];
}

int cornerLocal(const std::vector<int>& g2l, const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    const int g = mv.tris[gt][k];
    if (g < 0 || g >= (int)g2l.size()) return -1;
    return g2l[g];
}

double triAreaNormal(const MeshView& mv, int t, gp_XYZ& nOut) {
    const gp_XYZ a = cornerPt(mv, t, 0);
    const gp_XYZ b = cornerPt(mv, t, 1);
    const gp_XYZ c = cornerPt(mv, t, 2);
    const gp_XYZ n = (b - a).Crossed(c - a);
    const double mag = n.Modulus();
    if (mag > 0.0) nOut = n / mag;
    else nOut = gp_XYZ(0, 0, 1);
    return 0.5 * mag;
}

double triArea(const MeshView& mv, int t) {
    gp_XYZ n;
    return triAreaNormal(mv, t, n);
}

bool frameOf(const gp_XYZ& aUnit, gp_XYZ& u, gp_XYZ& v) {
    const double nx = std::abs(aUnit.X());
    const double ny = std::abs(aUnit.Y());
    const double nz = std::abs(aUnit.Z());
    gp_XYZ w;
    if (nx <= ny && nx <= nz) w = gp_XYZ(1, 0, 0);
    else if (ny <= nz) w = gp_XYZ(0, 1, 0);
    else w = gp_XYZ(0, 0, 1);
    u = aUnit.Crossed(w);
    const double um = u.Modulus();
    if (um < 1e-15) return false;
    u.Divide(um);
    v = aUnit.Crossed(u);
    const double vm = v.Modulus();
    if (vm < 1e-15) return false;
    v.Divide(vm);
    return true;
}

double sectionDist(const gp_XYZ& p, const gp_XYZ& q, const gp_XYZ& axis) {
    return (p - q).Crossed(axis).Modulus();
}

double radialDist(const gp_XYZ& p, const gp_XYZ& center, const gp_XYZ& axis) {
    return (p - center).Crossed(axis).Modulus();
}

// Unsigned angle between two axes, in [0, π/2].
double axisAngle(const gp_XYZ& a, const gp_XYZ& b) {
    const double c = std::min(1.0, std::abs(a.Dot(b)));
    return std::acos(c);
}

bool solve3(const gp_XYZ& r0, const gp_XYZ& r1, const gp_XYZ& r2,
            double b0, double b1, double b2, gp_XYZ& x) {
    double A[3][4] = {
        {r0.X(), r0.Y(), r0.Z(), b0},
        {r1.X(), r1.Y(), r1.Z(), b1},
        {r2.X(), r2.Y(), r2.Z(), b2},
    };
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int row = col + 1; row < 3; ++row)
            if (std::abs(A[row][col]) > std::abs(A[piv][col])) piv = row;
        if (!(std::abs(A[piv][col]) > 1e-14)) return false;
        if (piv != col) {
            for (int k = col; k < 4; ++k) std::swap(A[piv][k], A[col][k]);
        }
        const double div = A[col][col];
        for (int k = col; k < 4; ++k) A[col][k] /= div;
        for (int row = 0; row < 3; ++row) {
            if (row == col) continue;
            const double f = A[row][col];
            for (int k = col; k < 4; ++k) A[row][k] -= f * A[col][k];
        }
    }
    x = gp_XYZ(A[0][3], A[1][3], A[2][3]);
    return std::isfinite(x.X()) && std::isfinite(x.Y()) && std::isfinite(x.Z());
}

bool invertSym(std::array<std::array<double, 4>, 4> A, std::array<std::array<double, 4>, 4>& inv) {
    double M[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M[i][j] = A[i][j];
        for (int j = 0; j < 4; ++j) M[i][4 + j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int row = col + 1; row < 4; ++row)
            if (std::abs(M[row][col]) > std::abs(M[piv][col])) piv = row;
        if (!(std::abs(M[piv][col]) > 1e-14)) return false;
        if (piv != col) {
            for (int k = 0; k < 8; ++k) std::swap(M[piv][k], M[col][k]);
        }
        const double div = M[col][col];
        for (int k = col; k < 8; ++k) M[col][k] /= div;
        for (int row = 0; row < 4; ++row) {
            if (row == col) continue;
            const double f = M[row][col];
            for (int k = col; k < 8; ++k) M[row][k] -= f * M[col][k];
        }
    }
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) inv[i][j] = M[i][4 + j];
    return true;
}

// D2.2 facet-normal clustering. Same rule the fillet stage uses for nSides.
int nBandsFromNormals(const std::vector<double>& psi, double thetaBin) {
    if (psi.empty()) return 0;
    std::vector<double> s = psi;
    std::sort(s.begin(), s.end());
    const int n = (int)s.size();
    if (n == 1) return 1;
    std::vector<double> gaps((std::size_t)n);
    for (int i = 0; i < n - 1; ++i)
        gaps[(std::size_t)i] = s[(std::size_t)i + 1] - s[(std::size_t)i];
    gaps[(std::size_t)n - 1] = (s[0] + kTwoPi) - s[(std::size_t)n - 1];
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
    int bands = 1;
    double prev = s[(std::size_t)start];
    for (int k = 1; k < n; ++k) {
        const int idx = (start + k) % n;
        const double cur = s[(std::size_t)idx];
        double gap = cur - prev;
        if (gap < 0.0) gap += kTwoPi;
        if (gap > th) ++bands;
        prev = cur;
    }
    return bands;
}

struct Vert {
    int id = -1;
    gp_XYZ p;
};

struct Flank {
    int prov = -1;
    double contact = 0;
    gp_XYZ n;
    gp_XYZ p0;
    double maxDev = 0;
    double tol = 0;
    bool certified = false;
};

struct Fit {
    bool ok = false;
    bool anchored = false;
    int params = 0;
    gp_XYZ axis;
    gp_XYZ center;
    double radius = 0;
    double maxDev = 0;
    double rms = 0;
};

struct BandDiag {
    const char* mode = "fillet";
    const char* why = "degenerate";
    const char* fit = "free";
    int seed = -1;
    int nTri = 0;
    int nV = 0;
    int stations = 0;
    int levels = 0;
    int params = 0;
    double radius = 0;
    gp_XYZ axis;
    double maxDev = 0;
    double lineRes = 0;
    double span = 0;
    double dihedral = 0;
    double spanTolV = 0;
    double tanDistA = 0, tanDistB = 0;
    double tanAngA = 0, tanAngB = 0;
    int straddle = 0;
    double tau = 0;
};

void emitDiag(const BandDiag& d) {
    if (!diagOn()) return;
    const double tau = d.tau > 0.0 ? d.tau : 1.0;
    std::fprintf(stderr,
                 "DIAG_REL_BAND stage=extrude mode=%s seed=%d nTri=%d nV=%d "
                 "stations=%d levels=%d params=%d fit=%s radius=%.9g halfAngle=0 "
                 "axis=%.9g,%.9g,%.9g maxDevMM=%.9g maxDevTau=%.9g lineResTau=%.9g "
                 "span=%.9g dihedral=%.9g spanTol=%.9g tanDistA=%.9g tanDistB=%.9g "
                 "tanAngA=%.9g tanAngB=%.9g straddle=%d why=%s\n",
                 d.mode, d.seed, d.nTri, d.nV, d.stations, d.levels, d.params, d.fit,
                 d.radius, d.axis.X(), d.axis.Y(), d.axis.Z(), d.maxDev, d.maxDev / tau,
                 d.lineRes / tau, d.span, d.dihedral, d.spanTolV, d.tanDistA, d.tanDistB,
                 d.tanAngA, d.tanAngB, d.straddle, d.why);
}

struct WhyCount {
    int bands = 0;
    int accept = 0;
    int line = 0;
    int vertexOff = 0;
    int notBlend = 0;
    int sameSurface = 0;
    int underdetermined = 0;
    int degenerate = 0;
    int noAnchor = 0;
};

void tallyWhy(WhyCount& c, const char* why) {
    ++c.bands;
    if (!why) return;
    if (std::strcmp(why, "accept") == 0) ++c.accept;
    else if (std::strcmp(why, "line") == 0) ++c.line;
    else if (std::strcmp(why, "vertex-off") == 0) ++c.vertexOff;
    else if (std::strcmp(why, "not-blend") == 0) ++c.notBlend;
    else if (std::strcmp(why, "same-surface") == 0) ++c.sameSurface;
    else if (std::strcmp(why, "underdetermined") == 0) ++c.underdetermined;
    else if (std::strcmp(why, "no-anchor") == 0) ++c.noAnchor;
    else ++c.degenerate;
}

bool edgeParallel(const gp_XYZ& dir, double len, const gp_XYZ& axis, double tau) {
    if (!(len > tau)) return false;
    const double s = std::min(1.0, dir.Crossed(axis).Modulus());
    const double ang = std::asin(s);
    // Two rulings of one strip. tolA = τ (certified floor): tanAngTol(τ, τ, L).
    return ang <= tanAngTol(tau, tau, len);
}

double sectionExtent(const std::vector<gp_XYZ>& pts, const gp_XYZ& axis) {
    if (pts.size() < 2) return 0.0;
    gp_XYZ u, v;
    if (!frameOf(axis, u, v)) return 0.0;
    double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
    for (const gp_XYZ& p : pts) {
        const double x = p.Dot(u);
        const double y = p.Dot(v);
        minx = std::min(minx, x);
        maxx = std::max(maxx, x);
        miny = std::min(miny, y);
        maxy = std::max(maxy, y);
    }
    return std::hypot(maxx - minx, maxy - miny);
}

double axialExtent(const std::vector<gp_XYZ>& pts, const gp_XYZ& axis) {
    if (pts.empty()) return 0.0;
    double mn = 1e300, mx = -1e300;
    for (const gp_XYZ& p : pts) {
        const double z = p.Dot(axis);
        mn = std::min(mn, z);
        mx = std::max(mx, z);
    }
    return mx - mn;
}

void uniqueVerts(const MeshView& mv, const std::vector<int>& g2l,
                 const std::vector<int>& tris, std::vector<Vert>& out) {
    std::vector<int> ids;
    ids.reserve(tris.size() * 3);
    for (int t : tris) {
        if (t < 0 || (std::size_t)t >= mv.nTri) continue;
        for (int k = 0; k < 3; ++k) {
            const int lv = cornerLocal(g2l, mv, t, k);
            if (lv >= 0) ids.push_back(lv);
        }
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    out.clear();
    out.reserve(ids.size());
    for (int lv : ids) {
        Vert v;
        v.id = lv;
        const int g = mv.compVtx ? mv.compVtx[lv] : lv;
        v.p = mv.pts[g];
        out.push_back(v);
    }
}

void provPoints(const MeshView& mv, const std::vector<int>& g2l, const Provisional& pr,
                std::vector<gp_XYZ>& pts) {
    std::vector<Vert> vs;
    uniqueVerts(mv, g2l, pr.tris, vs);
    pts.clear();
    pts.reserve(vs.size());
    for (const Vert& v : vs) pts.push_back(v.p);
}

bool gaussPole(const MeshView& mv, const std::vector<int>& tris, gp_XYZ& axis, bool& determined) {
    determined = false;
    std::array<std::array<double, 3>, 3> G{};
    for (auto& row : G) row.fill(0.0);
    double areaSum = 0.0;
    for (int t : tris) {
        gp_XYZ n;
        const double a = triAreaNormal(mv, t, n);
        if (!(a > 0.0)) continue;
        G[0][0] += a * n.X() * n.X();
        G[0][1] += a * n.X() * n.Y();
        G[0][2] += a * n.X() * n.Z();
        G[1][0] += a * n.Y() * n.X();
        G[1][1] += a * n.Y() * n.Y();
        G[1][2] += a * n.Y() * n.Z();
        G[2][0] += a * n.Z() * n.X();
        G[2][1] += a * n.Z() * n.Y();
        G[2][2] += a * n.Z() * n.Z();
        areaSum += a;
    }
    if (!(areaSum > 0.0)) return false;
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(G, eval, evec)) return false;
    axis = gp_XYZ(evec[0][0], evec[0][1], evec[0][2]);
    const double am = axis.Modulus();
    if (!(am > 1e-15)) return false;
    axis.Divide(am);
    // A point-like gauss map (one normal) does not name an axis.
    determined = eval[1] > 1e-12 * std::max(eval[2], 1e-300);
    return true;
}

void alignTo(gp_XYZ& axis, const gp_XYZ& hemi) {
    if (axis.Dot(hemi) < 0.0) axis.Reverse();
}

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
}

struct ChainGeom {
    std::vector<int> tris;
    std::vector<Vert> verts;
    std::vector<Flank> flanks;
    std::vector<Flank> caps;
    gp_XYZ ruling;
    double seedChord = 0;
};

// Stations: section points merged at τ (countProfPts). Levels: axial (countLevels).
void profileCounts(const std::vector<Vert>& verts, const gp_XYZ& axis, double tau,
                   int& stations, int& levels, double& lineRes) {
    stations = 0;
    levels = 0;
    lineRes = 0;
    gp_XYZ u, v;
    if (!frameOf(axis, u, v) || verts.empty()) return;
    std::vector<ProfPt> sec;
    std::vector<double> z;
    sec.reserve(verts.size());
    z.reserve(verts.size());
    const gp_XYZ o = verts.front().p;
    for (const Vert& vt : verts) {
        const gp_XYZ d = vt.p - o;
        ProfPt p;
        p.rho = d.Dot(u);
        p.z = d.Dot(v);
        sec.push_back(p);
        z.push_back(d.Dot(axis));
    }
    stations = countProfPts(sec, tau);
    levels = countLevels(z, tau);
    lineRes = profileLineMaxResid(sec);
}

double devStats(const std::vector<Vert>& verts, const gp_XYZ& center, const gp_XYZ& axis,
                double radius, double& rms) {
    double mx = 0, ss = 0;
    int n = 0;
    for (const Vert& v : verts) {
        const double d = std::abs(radialDist(v.p, center, axis) - radius);
        mx = std::max(mx, d);
        ss += d * d;
        ++n;
    }
    rms = (n > 0) ? std::sqrt(ss / (double)n) : 0.0;
    return mx;
}

bool fitAnchored(const gp_XYZ& nA, const gp_XYZ& pA, const gp_XYZ& nB, const gp_XYZ& pB,
                 const std::vector<Vert>& verts, int maxIter, gp_XYZ ruling,
                 gp_XYZ& axis, gp_XYZ& center, double& radius) {
    if (verts.empty()) return false;
    gp_XYZ a = nA.Crossed(nB);
    const double am = a.Modulus();
    if (!(am > 1e-12)) return false;
    a.Divide(am);
    alignTo(a, ruling);
    gp_XYZ C(0, 0, 0);
    for (const Vert& v : verts) C += v.p;
    C /= (double)verts.size();
    const double sA = ((C - pA).Dot(nA) >= 0.0) ? 1.0 : -1.0;
    const double sB = ((C - pB).Dot(nB) >= 0.0) ? 1.0 : -1.0;
    gp_XYZ c0, c1;
    if (!solve3(nA, nB, a, pA.Dot(nA), pB.Dot(nB), C.Dot(a), c0)) return false;
    if (!solve3(nA, nB, a, pA.Dot(nA) + sA, pB.Dot(nB) + sB, C.Dot(a), c1)) return false;
    const gp_XYZ cR = c1 - c0;
    const gp_XYZ dax = cR.Crossed(a);

    auto sse = [&](double r) {
        const gp_XYZ c = c0 + cR * r;
        double e = 0.0;
        for (const Vert& v : verts) {
            const double rho = radialDist(v.p, c, a);
            const double f = rho - r;
            e += f * f;
        }
        return e;
    };

    std::vector<double> roots;
    roots.reserve(verts.size());
    for (const Vert& v : verts) {
        const gp_XYZ w = (v.p - c0).Crossed(a);
        const double A = dax.Dot(dax) - 1.0;
        const double B = -2.0 * w.Dot(dax);
        const double Cq = w.Dot(w);
        if (std::abs(A) < 1e-12) {
            if (std::abs(B) > 1e-12) {
                const double r = -Cq / B;
                if (r > 0.0 && std::isfinite(r)) roots.push_back(r);
            }
            continue;
        }
        const double disc = B * B - 4.0 * A * Cq;
        if (disc < 0.0) continue;
        const double s = std::sqrt(disc);
        const double r1 = (-B + s) / (2.0 * A);
        const double r2 = (-B - s) / (2.0 * A);
        if (r1 > 0.0 && std::isfinite(r1)) roots.push_back(r1);
        if (r2 > 0.0 && std::isfinite(r2)) roots.push_back(r2);
    }
    if (roots.empty()) return false;
    std::sort(roots.begin(), roots.end());
    double r = roots[roots.size() / 2];
    double bestE = sse(r);
    // One least-squares start: the root of the lowest-id vertex with the smaller SSE.
    {
        const gp_XYZ w = (verts.front().p - c0).Crossed(a);
        const double A = dax.Dot(dax) - 1.0;
        const double B = -2.0 * w.Dot(dax);
        const double Cq = w.Dot(w);
        double cand[2] = {r, r};
        int nc = 1;
        if (std::abs(A) < 1e-12) {
            if (std::abs(B) > 1e-12) {
                const double rr = -Cq / B;
                if (rr > 0.0) cand[nc++] = rr;
            }
        } else {
            const double disc = B * B - 4.0 * A * Cq;
            if (disc >= 0.0) {
                const double s = std::sqrt(disc);
                const double r1 = (-B + s) / (2.0 * A);
                const double r2 = (-B - s) / (2.0 * A);
                nc = 0;
                if (r1 > 0.0) cand[nc++] = r1;
                if (r2 > 0.0 && nc < 2) cand[nc++] = r2;
            }
        }
        for (int i = 0; i < nc; ++i) {
            const double e = sse(cand[i]);
            if (e < bestE) {
                bestE = e;
                r = cand[i];
            }
        }
    }

    const int itCap = std::max(1, maxIter);
    for (int it = 0; it < itCap; ++it) {
        const gp_XYZ c = c0 + cR * r;
        double g = 0.0, h = 0.0;
        for (const Vert& v : verts) {
            const gp_XYZ w = (v.p - c0).Crossed(a);
            const gp_XYZ rad = w - dax * r;
            const double rho = rad.Modulus();
            if (!(rho > 1e-15)) continue;
            const double dr = -(rad.Dot(dax)) / rho;
            const double J = dr - 1.0;
            const double e = rho - r;
            g += J * e;
            h += J * J;
        }
        if (!(h > 1e-18)) break;
        const double step = -g / h;
        if (!std::isfinite(step)) break;
        double r2 = r + step;
        if (!(r2 > 0.0)) break;
        const double e2 = sse(r2);
        if (e2 <= bestE) {
            r = r2;
            bestE = e2;
        } else {
            break;
        }
    }
    if (!(r > 0.0) || !std::isfinite(r)) return false;
    axis = a;
    center = c0 + cR * r;
    center += a * ((C - center).Dot(a));
    radius = r;
    return true;
}

// O2: centre on the mid-line (the offset bisector of antiparallel planes), one
// parameter along that line. r is the mean radial distance and is checked, not imposed.
bool fitHem(const gp_XYZ& nA, const gp_XYZ& pA, const gp_XYZ& nB, const gp_XYZ& pB,
            const std::vector<Vert>& verts, const gp_XYZ& ruling, int maxIter,
            gp_XYZ& axis, gp_XYZ& center, double& radius) {
    if (verts.empty()) return false;
    gp_XYZ a = ruling;
    const double am = a.Modulus();
    if (!(am > 1e-15)) return false;
    a.Divide(am);
    gp_XYZ nGap = nA;
    if ((pB - pA).Dot(nA) < 0.0) nGap = nA * -1.0;
    const double gap = std::abs((pB - pA).Dot(nA));
    const gp_XYZ mid = pA + nGap * (0.5 * gap);
    gp_XYZ m = a.Crossed(nGap);
    const double mm = m.Modulus();
    if (!(mm > 1e-15)) return false;
    m.Divide(mm);
    gp_XYZ C(0, 0, 0);
    for (const Vert& v : verts) C += v.p;
    C /= (double)verts.size();
    double s = (C - mid).Dot(m);
    // Orient m so the bulge (centroid off the gap) is the positive side.
    if ((C - mid).Dot(m) < 0.0) {
        m.Reverse();
        s = -s;
    }
    const gp_XYZ dax = m.Crossed(a);
    const int itCap = std::max(1, maxIter);
    for (int it = 0; it < itCap; ++it) {
        std::vector<double> rho(verts.size(), 0.0);
        std::vector<double> dr(verts.size(), 0.0);
        double sum = 0.0, sumJ = 0.0;
        int nOk = 0;
        for (std::size_t i = 0; i < verts.size(); ++i) {
            const gp_XYZ w = (verts[i].p - mid).Crossed(a);
            const gp_XYZ rad = w - dax * s;
            const double mrad = rad.Modulus();
            if (!(mrad > 1e-15)) continue;
            rho[i] = mrad;
            dr[i] = -(rad.Dot(dax)) / mrad;
            sum += mrad;
            sumJ += dr[i];
            ++nOk;
        }
        if (nOk < 1) return false;
        const double rMean = sum / (double)nOk;
        const double jMean = sumJ / (double)nOk;
        double g = 0.0, h = 0.0;
        for (std::size_t i = 0; i < verts.size(); ++i) {
            if (!(rho[i] > 0.0)) continue;
            const double J = dr[i] - jMean;
            const double e = rho[i] - rMean;
            g += J * e;
            h += J * J;
        }
        if (!(h > 1e-18)) break;
        const double step = -g / h;
        if (!std::isfinite(step)) break;
        s += step;
    }
    const gp_XYZ c = mid + m * s;
    double sum = 0.0;
    int nOk = 0;
    for (const Vert& v : verts) {
        const double rho = radialDist(v.p, c, a);
        if (!(rho > 1e-15)) continue;
        sum += rho;
        ++nOk;
    }
    if (nOk < 1) return false;
    axis = a;
    center = c;
    radius = sum / (double)nOk;
    return radius > 0.0 && std::isfinite(radius);
}

// ℓ1 norm of the hat-matrix row of a new radial sample. Rank-deficient → 1
// (the one-parameter mean), so extrapBound falls back to q.
double hatRowL1(const std::vector<std::array<double, 4>>& J, int p,
                const std::array<double, 4>& jNew) {
    if (J.empty() || p < 1) return 1.0;
    std::array<std::array<double, 4>, 4> A{};
    for (auto& row : A) row.fill(0.0);
    for (const auto& row : J) {
        for (int a = 0; a < p; ++a)
            for (int b = 0; b < p; ++b) A[a][b] += row[a] * row[b];
    }
    std::array<std::array<double, 4>, 4> inv{};
    if (p == 4) {
        if (!invertSym(A, inv)) return 1.0;
    } else {
        if (!(std::abs(A[0][0]) > 1e-18)) return 1.0;
        inv[0][0] = 1.0 / A[0][0];
    }
    std::array<double, 4> t{};
    t.fill(0.0);
    for (int a = 0; a < p; ++a)
        for (int b = 0; b < p; ++b) t[a] += jNew[b] * inv[b][a];
    double l1 = 0.0;
    for (const auto& row : J) {
        double h = 0.0;
        for (int a = 0; a < p; ++a) h += t[a] * row[a];
        l1 += std::abs(h);
    }
    if (!std::isfinite(l1) || !(l1 > 0.0)) return 1.0;
    return l1;
}

void jacobianRows(const std::vector<Vert>& verts, const gp_XYZ& center, const gp_XYZ& axis,
                  double radius, bool anchored, const gp_XYZ& cR,
                  std::vector<std::array<double, 4>>& rows) {
    rows.clear();
    rows.reserve(verts.size());
    gp_XYZ u, v;
    if (!frameOf(axis, u, v)) return;
    if (anchored) {
        const gp_XYZ dax = cR.Crossed(axis);
        for (const Vert& vt : verts) {
            const double rho = radialDist(vt.p, center, axis);
            if (!(rho > 1e-15)) continue;
            const gp_XYZ w = (vt.p - (center - cR * radius)).Crossed(axis);
            const gp_XYZ rad = w - dax * radius;
            const double dr = -(rad.Dot(dax)) / rho;
            std::array<double, 4> row{{dr - 1.0, 0, 0, 0}};
            rows.push_back(row);
        }
        return;
    }
    // Free cylinder, r eliminated as the mean (cylinderFitLS): 4 columns.
    std::vector<std::array<double, 4>> raw;
    raw.reserve(verts.size());
    std::array<double, 4> mean{};
    mean.fill(0.0);
    for (const Vert& vt : verts) {
        const gp_XYZ d = vt.p - center;
        const double wd = d.Dot(axis);
        const gp_XYZ rad = d - axis * wd;
        const double m = rad.Modulus();
        if (!(m > 1e-15)) continue;
        const gp_XYZ rhat = rad / m;
        std::array<double, 4> row{{
            -wd * d.Dot(u) / m,
            -wd * d.Dot(v) / m,
            -rhat.Dot(u),
            -rhat.Dot(v),
        }};
        raw.push_back(row);
        for (int k = 0; k < 4; ++k) mean[k] += row[k];
    }
    if (raw.empty()) return;
    for (int k = 0; k < 4; ++k) mean[k] /= (double)raw.size();
    for (auto row : raw) {
        for (int k = 0; k < 4; ++k) row[k] -= mean[k];
        rows.push_back(row);
    }
    (void)radius;
}

double extrapAt(const gp_XYZ& p, const std::vector<std::array<double, 4>>& rows, int pcount,
                const gp_XYZ& center, const gp_XYZ& axis, double radius, bool anchored,
                const gp_XYZ& cR, double q) {
    std::array<double, 4> jn{{0, 0, 0, 0}};
    if (anchored) {
        const gp_XYZ dax = cR.Crossed(axis);
        const double rho = radialDist(p, center, axis);
        if (rho > 1e-15) {
            const gp_XYZ w = (p - (center - cR * radius)).Crossed(axis);
            const gp_XYZ rad = w - dax * radius;
            jn[0] = -(rad.Dot(dax)) / rho - 1.0;
        }
        return q * hatRowL1(rows, 1, jn);
    }
    gp_XYZ u, v;
    if (!frameOf(axis, u, v)) return q;
    const gp_XYZ d = p - center;
    const double wd = d.Dot(axis);
    const gp_XYZ rad = d - axis * wd;
    const double m = rad.Modulus();
    if (!(m > 1e-15)) return q;
    const gp_XYZ rhat = rad / m;
    jn[0] = -wd * d.Dot(u) / m;
    jn[1] = -wd * d.Dot(v) / m;
    jn[2] = -rhat.Dot(u);
    jn[3] = -rhat.Dot(v);
    return q * hatRowL1(rows, pcount, jn);
}

// Vertex azimuth span: 2π minus the largest gap. Closed-360 uses D-fc-1.6.13.
struct SpanInfo {
    double span = 0;
    double maxGap = 0;
    double uMin = 0;
    double uMax = 0;
    int nBands = 0;
    bool closed360 = false;
    double vMin = 0;
    double vMax = 0;
};

SpanInfo measureSpan(const MeshView& mv, const std::vector<int>& g2l,
                     const std::vector<int>& tris, const std::vector<Vert>& verts,
                     const gp_XYZ& center, const gp_XYZ& axis, const gp_XYZ& xdir,
                     double thetaBin) {
    SpanInfo s;
    const gp_XYZ vdir = axis.Crossed(xdir);
    std::vector<double> chi;
    chi.reserve(verts.size());
    double vMin = 1e300, vMax = -1e300;
    for (const Vert& vt : verts) {
        const gp_XYZ d = vt.p - center;
        vMin = std::min(vMin, d.Dot(axis));
        vMax = std::max(vMax, d.Dot(axis));
        const double ur = d.Dot(xdir);
        const double vr = d.Dot(vdir);
        if (ur * ur + vr * vr < 1e-24) continue;
        chi.push_back(std::atan2(vr, ur));
    }
    s.vMin = (vMin > vMax) ? 0.0 : vMin;
    s.vMax = (vMin > vMax) ? 0.0 : vMax;
    std::sort(chi.begin(), chi.end());
    double maxGap = 0.0;
    std::size_t gapAt = 0;
    if (chi.size() >= 2) {
        for (std::size_t i = 0; i + 1 < chi.size(); ++i) {
            const double g = chi[i + 1] - chi[i];
            if (g > maxGap) {
                maxGap = g;
                gapAt = i;
            }
        }
        const double wrap = (chi.front() + kTwoPi) - chi.back();
        if (wrap > maxGap) {
            maxGap = wrap;
            gapAt = chi.size() - 1;
        }
    } else {
        maxGap = kTwoPi;
    }
    s.maxGap = maxGap;
    s.span = kTwoPi - maxGap;
    if (s.span < 0.0) s.span = 0.0;
    if (chi.size() >= 2) {
        const std::size_t start = (gapAt + 1) % chi.size();
        s.uMin = chi[start];
        if (s.uMin > kPi) s.uMin -= kTwoPi;
        if (s.uMin <= -kPi) s.uMin += kTwoPi;
    }
    s.uMax = s.uMin + s.span;

    gp_XYZ u0, v0;
    if (!frameOf(axis, u0, v0)) return s;
    std::vector<double> psi;
    psi.reserve(tris.size());
    for (int t : tris) {
        gp_XYZ nt;
        triAreaNormal(mv, t, nt);
        psi.push_back(std::atan2(nt.Dot(v0), nt.Dot(u0)));
    }
    s.nBands = nBandsFromNormals(psi, thetaBin);
    if (s.nBands >= 3 && s.maxGap <= 1.5 * (kTwoPi / (double)s.nBands)) {
        s.closed360 = true;
        s.uMin = 0.0;
        s.uMax = kTwoPi;
        s.span = kTwoPi;
    }
    (void)g2l;
    return s;
}

gp_XYZ lowestXDir(const std::vector<Vert>& verts, const gp_XYZ& center, const gp_XYZ& axis) {
    for (const Vert& vt : verts) {
        gp_XYZ rad = vt.p - center;
        rad -= axis * rad.Dot(axis);
        const double m = rad.Modulus();
        if (m > 1e-12) return rad / m;
    }
    gp_XYZ w = (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    gp_XYZ x = axis.Crossed(w);
    const double xm = x.Modulus();
    if (xm > 1e-15) x.Divide(xm);
    else x = gp_XYZ(1, 0, 0);
    return x;
}

double areaOfTris(const MeshView& mv, const std::vector<int>& tris) {
    double a = 0.0;
    for (int t : tris) a += triArea(mv, t);
    return a;
}

void refitPlaneProv(const MeshView& mv, const std::vector<int>& g2l, Provisional& pr) {
    gp_Ax3 pl;
    if (!pr.tris.empty() && pcaPlane(mv, pr.tris, pl)) pr.plane = pl;
    pr.area = areaOfTris(mv, pr.tris);
    std::vector<Vert> vs;
    uniqueVerts(mv, g2l, pr.tris, vs);
    const gp_XYZ n = pr.plane.Direction().XYZ();
    const gp_XYZ o = pr.plane.Location().XYZ();
    double mx = 0.0, ss = 0.0;
    for (const Vert& v : vs) {
        const double d = std::abs((v.p - o).Dot(n));
        mx = std::max(mx, d);
        ss += d * d;
    }
    pr.maxVertexDev = mx;
    pr.rmsVertexDev = vs.empty() ? 0.0 : std::sqrt(ss / (double)vs.size());
}

}  // namespace

bool claimRelExtrude(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work) {
    try {
        const double tau = tauOf(mv);
        const double q = (mv.quantFloor > 0.0) ? mv.quantFloor : 0.0;
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0)
            return true;

        const int nTri = (int)mv.nTri;
        const int nVtx = (int)mv.nVtx;
        const int nEdge = (int)mv.nEdge;
        const int nProv = (int)work.provisionals.size();

        std::vector<int> g2l;
        int maxG = 0;
        for (int lv = 0; lv < nVtx; ++lv) {
            const int g = mv.compVtx ? mv.compVtx[lv] : lv;
            if (g > maxG) maxG = g;
        }
        for (int t = 0; t < nTri; ++t) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            for (int k = 0; k < 3; ++k) maxG = std::max(maxG, mv.tris[gt][k]);
        }
        g2l.assign((std::size_t)maxG + 1, -1);
        for (int lv = 0; lv < nVtx; ++lv) {
            const int g = mv.compVtx ? mv.compVtx[lv] : lv;
            if (g >= 0 && g < (int)g2l.size()) g2l[g] = lv;
        }

        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{{-1, -1}});
        std::vector<int> triProv(nTri, -1);
        for (int pi = 0; pi < nProv; ++pi) {
            for (int t : work.provisionals[pi].tris) {
                if (t >= 0 && t < nTri) triProv[t] = pi;
            }
        }
        for (int t = 0; t < nTri; ++t) {
            for (int s = 0; s < 3; ++s) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || e >= nEdge) continue;
                auto& et = edgeTris[(std::size_t)e];
                if (et[0] < 0) et[0] = t;
                else if (et[0] != t) et[1] = t;
            }
        }
        auto otherTri = [&](int e, int t) -> int {
            if (e < 0 || e >= nEdge) return -1;
            const int a = edgeTris[(std::size_t)e][0];
            const int b = edgeTris[(std::size_t)e][1];
            if (a == t) return b;
            if (b == t) return a;
            return -1;
        };

        std::vector<char> taken(nTri, 0);
        for (int t = 0; t < nTri; ++t) {
            const int pi = triProv[t];
            if (pi < 0 || work.provisionals[pi].claim != ProvClaim::Unclaimed) taken[t] = 1;
        }

        struct Seed {
            int prov;
            int minTri;
            gp_XYZ ruling;
            double chord;
        };
        std::vector<Seed> seeds;
        for (int pi = 0; pi < nProv; ++pi) {
            const Provisional& pr = work.provisionals[pi];
            if (pr.claim != ProvClaim::Unclaimed || pr.tris.empty()) continue;
            struct Fam {
                gp_XYZ dir;
                double len = 0;
                int n = 0;
                int minEdge = std::numeric_limits<int>::max();
            };
            std::vector<Fam> fams;
            bool anyNbr = false;
            for (int t : pr.tris) {
                if (t < 0 || t >= nTri || !mv.triEdges) continue;
                for (int s = 0; s < 3; ++s) {
                    const int e = mv.triEdges[t][s];
                    if (e < 0 || e >= nEdge || !mv.compEdges) continue;
                    const int ot = otherTri(e, t);
                    const int op = (ot >= 0 && ot < nTri) ? triProv[ot] : -1;
                    if (op == pi) continue;  // interior to the provisional
                    const int va = mv.compEdges[e].first;
                    const int vb = mv.compEdges[e].second;
                    const int ga = mv.compVtx ? mv.compVtx[va] : va;
                    const int gb = mv.compVtx ? mv.compVtx[vb] : vb;
                    gp_XYZ d = mv.pts[gb] - mv.pts[ga];
                    const double len = d.Modulus();
                    if (!(len > tau)) continue;
                    d.Divide(len);
                    if (op >= 0 && op != pi) anyNbr = true;
                    int slot = -1;
                    for (int f = 0; f < (int)fams.size(); ++f) {
                        gp_XYZ fd = fams[f].dir;
                        if (fd.Dot(d) < 0.0) fd.Reverse();
                        if (edgeParallel(d, len, fd, tau) || axisAngle(d, fams[f].dir) <=
                                                                  tanAngTol(tau, tau, len)) {
                            slot = f;
                            break;
                        }
                    }
                    if (slot < 0) {
                        Fam f;
                        f.dir = d;
                        f.len = len;
                        f.n = 1;
                        f.minEdge = e;
                        fams.push_back(f);
                    } else {
                        if (fams[slot].dir.Dot(d) < 0.0) d.Reverse();
                        const double nl = fams[slot].len + len;
                        fams[slot].dir = fams[slot].dir * fams[slot].len + d * len;
                        const double dm = fams[slot].dir.Modulus();
                        if (dm > 1e-15) fams[slot].dir.Divide(dm);
                        fams[slot].len = nl;
                        fams[slot].n += 1;
                        fams[slot].minEdge = std::min(fams[slot].minEdge, e);
                    }
                }
            }
            if (!anyNbr || fams.empty()) continue;
            int best = 0;
            for (int f = 1; f < (int)fams.size(); ++f) {
                if (fams[f].len > fams[best].len + 1e-15) best = f;
                else if (std::abs(fams[f].len - fams[best].len) <= 1e-15 &&
                         fams[f].minEdge < fams[best].minEdge)
                    best = f;
            }
            if (fams[best].n < 2) continue;
            std::vector<gp_XYZ> pts;
            provPoints(mv, g2l, pr, pts);
            const double chord = sectionExtent(pts, fams[best].dir);
            const double axial = axialExtent(pts, fams[best].dir);
            // A face whose section is wider than its ruling is not a strip seed.
            if (chord > axial + tau) continue;
            Seed sd;
            sd.prov = pi;
            sd.minTri = minTriOf(pr.tris);
            sd.ruling = fams[best].dir;
            sd.chord = chord;
            seeds.push_back(sd);
        }
        std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
            if (a.minTri != b.minTri) return a.minTri < b.minTri;
            return a.prov < b.prov;
        });

        WhyCount yc;
        std::vector<char> provInBand(nProv, 0);

        for (const Seed& sd : seeds) {
            if (sd.prov < 0 || sd.prov >= nProv) continue;
            if (work.provisionals[sd.prov].claim != ProvClaim::Unclaimed) continue;
            bool seedTaken = true;
            for (int t : work.provisionals[sd.prov].tris) {
                if (t >= 0 && t < nTri && !taken[t]) seedTaken = false;
            }
            if (seedTaken) continue;

            std::vector<int> band;
            std::vector<char> inBand(nTri, 0);
            provInBand.assign(nProv, 0);
            std::vector<int> queue;
            queue.push_back(sd.prov);
            provInBand[sd.prov] = 1;
            for (int t : work.provisionals[sd.prov].tris) {
                if (t >= 0 && t < nTri && !taken[t]) {
                    band.push_back(t);
                    inBand[t] = 1;
                }
            }
            if (band.empty()) continue;

            struct Touch {
                int prov;
                double contact;
            };
            std::vector<Touch> touches;
            auto addTouch = [&](int op, double len) {
                if (op < 0 || op >= nProv) return;
                for (Touch& tc : touches) {
                    if (tc.prov == op) {
                        tc.contact += len;
                        return;
                    }
                }
                touches.push_back({op, len});
            };

            std::vector<int> layer = queue;
            while (!layer.empty()) {
                std::vector<int> next;
                std::vector<std::pair<int, double>> cand;  // prov, contact
                for (int pi : layer) {
                    if (pi < 0 || pi >= nProv) continue;
                    for (int t : work.provisionals[pi].tris) {
                        if (t < 0 || t >= nTri || !inBand[t]) continue;
                        for (int s = 0; s < 3; ++s) {
                            const int e = mv.triEdges[t][s];
                            if (e < 0 || e >= nEdge || !mv.compEdges) continue;
                            const int va = mv.compEdges[e].first;
                            const int vb = mv.compEdges[e].second;
                            const int ga = mv.compVtx ? mv.compVtx[va] : va;
                            const int gb = mv.compVtx ? mv.compVtx[vb] : vb;
                            gp_XYZ ed = mv.pts[gb] - mv.pts[ga];
                            const double elen = ed.Modulus();
                            if (!edgeParallel(ed / std::max(elen, 1e-30), elen, sd.ruling, tau) &&
                                !(elen > tau && axisAngle(ed / elen, sd.ruling) <=
                                                    tanAngTol(tau, tau, elen)))
                                continue;
                            const int ot = otherTri(e, t);
                            if (ot < 0 || ot >= nTri || taken[ot]) continue;
                            const int op = triProv[ot];
                            if (op < 0 || op >= nProv || provInBand[op]) continue;
                            if (work.provisionals[op].claim != ProvClaim::Unclaimed) continue;
                            addTouch(op, elen);
                            bool seen = false;
                            for (auto& c : cand)
                                if (c.first == op) {
                                    c.second += elen;
                                    seen = true;
                                    break;
                                }
                            if (!seen) cand.push_back({op, elen});
                        }
                    }
                }
                std::sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) {
                    const int ma = a.first;
                    const int mb = b.first;
                    if (ma != mb) return ma < mb;
                    return a.second > b.second;
                });
                cand.erase(std::unique(cand.begin(), cand.end(),
                                       [](const auto& a, const auto& b) {
                                           return a.first == b.first;
                                       }),
                           cand.end());
                for (const auto& c : cand) {
                    const int op = c.first;
                    if (provInBand[op]) continue;
                    std::vector<gp_XYZ> pts;
                    provPoints(mv, g2l, work.provisionals[op], pts);
                    const double ext = sectionExtent(pts, sd.ruling);
                    // Anchor faces are wider in section than the seed strip.
                    // The slack is τ (one vertex quantum on a measured length).
                    if (ext > sd.chord + tau) continue;
                    bool near = false;
                    if (!pts.empty() && !band.empty()) {
                        near = true;
                        for (const gp_XYZ& pv : pts) {
                            double md = 1e300;
                            for (int bt : band) {
                                for (int k = 0; k < 3; ++k) {
                                    md = std::min(md, sectionDist(pv, cornerPt(mv, bt, k),
                                                                   sd.ruling));
                                }
                            }
                            if (md > sd.chord + tau) {
                                near = false;
                                break;
                            }
                        }
                    }
                    if (!near) continue;
                    provInBand[op] = 1;
                    next.push_back(op);
                    for (int t : work.provisionals[op].tris) {
                        if (t >= 0 && t < nTri && !taken[t] && !inBand[t]) {
                            band.push_back(t);
                            inBand[t] = 1;
                        }
                    }
                }
                layer.swap(next);
            }
            std::sort(band.begin(), band.end());
            band.erase(std::unique(band.begin(), band.end()), band.end());
            if (band.empty()) continue;

            std::vector<Vert> verts;
            uniqueVerts(mv, g2l, band, verts);

            std::vector<Flank> flanks;
            std::vector<Flank> caps;
            for (const Touch& tc : touches) {
                if (tc.prov < 0 || tc.prov >= nProv) continue;
                if (provInBand[tc.prov]) continue;
                const Provisional& ap = work.provisionals[tc.prov];
                if (ap.claim != ProvClaim::Unclaimed) continue;
                if (!a3WouldCommit(mv, tol, ap, false)) continue;
                Flank f;
                f.prov = tc.prov;
                f.contact = tc.contact;
                f.n = ap.plane.Direction().XYZ();
                const double nm = f.n.Modulus();
                if (!(nm > 1e-15)) continue;
                f.n.Divide(nm);
                f.p0 = ap.plane.Location().XYZ();
                f.maxDev = ap.maxVertexDev;
                f.tol = anchorTol(tau, ap.maxVertexDev);
                f.certified = anchorCertified(tau, ap.maxVertexDev);
                // Cap: plane normal parallel to the ruling (a floor the peg stands on).
                const double nang = axisAngle(f.n, sd.ruling);
                if (nang <= tanAngTol(f.tol, tau, std::max(axialExtent(
                    [&]() {
                        std::vector<gp_XYZ> ps;
                        provPoints(mv, g2l, ap, ps);
                        return ps;
                    }(),
                    sd.ruling), tau))) {
                    caps.push_back(f);
                } else {
                    flanks.push_back(f);
                }
            }
            std::sort(flanks.begin(), flanks.end(), [](const Flank& a, const Flank& b) {
                if (a.contact != b.contact) return a.contact > b.contact;
                return a.prov < b.prov;
            });
            if (flanks.size() > 2) flanks.resize(2);
            std::sort(flanks.begin(), flanks.end(),
                      [](const Flank& a, const Flank& b) { return a.prov < b.prov; });

            BandDiag dg;
            dg.tau = tau;
            dg.seed = sd.minTri;
            dg.nTri = (int)band.size();
            dg.nV = (int)verts.size();
            dg.axis = sd.ruling;
            dg.mode = "fillet";

            const bool two = flanks.size() == 2;
            bool anti = false;
            if (two) {
                const double cross = flanks[0].n.Crossed(flanks[1].n).Modulus();
                const double dot = flanks[0].n.Dot(flanks[1].n);
                const double Lref = std::max(axialExtent(
                    [&]() {
                        std::vector<gp_XYZ> ps;
                        ps.reserve(verts.size());
                        for (const Vert& v : verts) ps.push_back(v.p);
                        return ps;
                    }(),
                    sd.ruling), tau);
                const double atol = tanAngTol(std::max(flanks[0].tol, flanks[1].tol), tau, Lref);
                anti = dot < 0.0 && cross <= std::sin(std::min(1.0, atol));
                if (!anti && cross <= 1e-12 && dot < 0.0) anti = true;
            }
            const bool bothCert = two && flanks[0].certified && flanks[1].certified;
            const bool hem = two && anti;
            const bool filletPair = two && !anti;
            if (hem) dg.mode = "hem";
            else if (!filletPair && !caps.empty()) dg.mode = "peg";
            else if (!filletPair) dg.mode = "free";

            gp_XYZ hyp = sd.ruling;
            if (filletPair) {
                gp_XYZ cr = flanks[0].n.Crossed(flanks[1].n);
                const double cm = cr.Modulus();
                if (cm > 1e-12) {
                    cr.Divide(cm);
                    alignTo(cr, sd.ruling);
                    hyp = cr;
                }
            } else if (!caps.empty()) {
                gp_XYZ cn = caps[0].n;
                alignTo(cn, sd.ruling);
                hyp = cn;
            }

            const bool anchored = bothCert && (filletPair || hem);
            dg.fit = anchored ? "anchored" : "free";
            dg.params = anchored ? kParamsAnchored
                                 : (kParamsCircle + kParamsAxisDir);
            profileCounts(verts, hyp, tau, dg.stations, dg.levels, dg.lineRes);
            dg.axis = hyp;

            auto finish = [&](const char* why) {
                dg.why = why;
                emitDiag(dg);
                tallyWhy(yc, why);
            };

            if (dg.stations < minStations(dg.params) ||
                (!anchored && dg.levels < 2)) {
                finish("underdetermined");
                continue;
            }
            if (dg.lineRes <= tau) {
                finish("line");
                continue;
            }

            bool poleOk = false;
            gp_XYZ pole;
            const bool havePole = gaussPole(mv, band, pole, poleOk);
            if (havePole) alignTo(pole, hyp);
            std::vector<gp_XYZ> bandPts;
            bandPts.reserve(verts.size());
            for (const Vert& v : verts) bandPts.push_back(v.p);
            const double Lband = std::max(axialExtent(bandPts, hyp), tau);
            const double tolH = two ? std::max(flanks[0].tol, flanks[1].tol) : tau;
            if (poleOk && havePole && axisAngle(pole, hyp) > tanAngTol(tolH, tau, Lband)) {
                finish("degenerate");
                continue;
            }

            Fit fit;
            fit.anchored = anchored;
            fit.params = dg.params;
            gp_XYZ cR(0, 0, 0);
            if (anchored && filletPair) {
                gp_XYZ ax, c;
                double r = 0;
                fit.ok = fitAnchored(flanks[0].n, flanks[0].p0, flanks[1].n, flanks[1].p0,
                                     verts, p.maxRefineIters, hyp, ax, c, r);
                if (fit.ok) {
                    fit.axis = ax;
                    fit.center = c;
                    fit.radius = r;
                    // dc/dr for the hat row: rebuild cR from the same side choice.
                    gp_XYZ C(0, 0, 0);
                    for (const Vert& v : verts) C += v.p;
                    C /= (double)verts.size();
                    const double sA = ((C - flanks[0].p0).Dot(flanks[0].n) >= 0.0) ? 1.0 : -1.0;
                    const double sB = ((C - flanks[1].p0).Dot(flanks[1].n) >= 0.0) ? 1.0 : -1.0;
                    gp_XYZ c0, c1;
                    if (solve3(flanks[0].n, flanks[1].n, ax, flanks[0].p0.Dot(flanks[0].n),
                               flanks[1].p0.Dot(flanks[1].n), C.Dot(ax), c0) &&
                        solve3(flanks[0].n, flanks[1].n, ax,
                               flanks[0].p0.Dot(flanks[0].n) + sA,
                               flanks[1].p0.Dot(flanks[1].n) + sB, C.Dot(ax), c1))
                        cR = c1 - c0;
                }
            } else if (anchored && hem) {
                gp_XYZ ax, c;
                double r = 0;
                fit.ok = fitHem(flanks[0].n, flanks[0].p0, flanks[1].n, flanks[1].p0, verts,
                                hyp, p.maxRefineIters, ax, c, r);
                if (fit.ok) {
                    fit.axis = ax;
                    fit.center = c;
                    fit.radius = r;
                }
            } else {
                if ((int)verts.size() < 6) {
                    finish("underdetermined");
                    continue;
                }
                gp_Dir axOut;
                gp_Pnt cOut;
                double rOut = 0, mOut = 0;
                const gp_Dir seedDir(hyp.X(), hyp.Y(), hyp.Z());
                fit.ok = cylinderFitLS(mv, band, seedDir, axOut, cOut, rOut, mOut,
                                       p.maxRefineIters);
                if (fit.ok) {
                    fit.axis = axOut.XYZ();
                    alignTo(fit.axis, hyp);
                    fit.center = cOut.XYZ();
                    fit.radius = rOut;
                    if (poleOk && axisAngle(fit.axis, hyp) > tanAngTol(tolH, tau, Lband)) {
                        dg.axis = fit.axis;
                        dg.radius = fit.radius;
                        finish("degenerate");
                        continue;
                    }
                }
            }
            if (!fit.ok || !(fit.radius > 0.0)) {
                // Anchored failure is a refusal. Free failure likewise: no second model.
                finish(anchored ? "vertex-off" : "degenerate");
                continue;
            }
            fit.maxDev = devStats(verts, fit.center, fit.axis, fit.radius, fit.rms);
            dg.axis = fit.axis;
            dg.radius = fit.radius;
            dg.maxDev = fit.maxDev;
            if (fit.maxDev > tau) {
                finish("vertex-off");
                continue;
            }

            // Growth across shared edges. A round that leaves the τ certificate
            // is not admitted: the previous band is the maximal certified set.
            std::vector<char> isAnchor(nProv, 0);
            for (const Flank& f : flanks) isAnchor[f.prov] = 1;
            for (const Flank& f : caps) isAnchor[f.prov] = 1;
            const int growCap = std::max(1, p.maxRefineIters);
            for (int round = 0; round < growCap; ++round) {
                std::vector<std::array<double, 4>> rows;
                jacobianRows(verts, fit.center, fit.axis, fit.radius, anchored && filletPair,
                             cR, rows);
                const int pc = (anchored && filletPair) ? 1 : 4;
                std::vector<int> add;
                for (int t : band) {
                    for (int s = 0; s < 3; ++s) {
                        const int e = mv.triEdges[t][s];
                        const int ot = otherTri(e, t);
                        if (ot < 0 || ot >= nTri || inBand[ot] || taken[ot]) continue;
                        const int op = triProv[ot];
                        if (op >= 0 && op < nProv && isAnchor[op]) continue;
                        if (op >= 0 && op < nProv &&
                            work.provisionals[op].claim != ProvClaim::Unclaimed)
                            continue;
                        bool ok3 = true;
                        for (int k = 0; k < 3 && ok3; ++k) {
                            const gp_XYZ pt = cornerPt(mv, ot, k);
                            const double dev = std::abs(radialDist(pt, fit.center, fit.axis) -
                                                        fit.radius);
                            const double xb = extrapAt(pt, rows, pc, fit.center, fit.axis,
                                                       fit.radius, anchored && filletPair, cR, q);
                            if (dev > tau + xb) ok3 = false;
                        }
                        if (ok3) add.push_back(ot);
                    }
                }
                std::sort(add.begin(), add.end());
                add.erase(std::unique(add.begin(), add.end()), add.end());
                if (add.empty()) break;
                std::vector<int> trial = band;
                trial.insert(trial.end(), add.begin(), add.end());
                std::sort(trial.begin(), trial.end());
                std::vector<Vert> tv;
                uniqueVerts(mv, g2l, trial, tv);
                Fit tf = fit;
                bool tok = false;
                if (anchored && filletPair) {
                    gp_XYZ ax, c;
                    double r = 0;
                    tok = fitAnchored(flanks[0].n, flanks[0].p0, flanks[1].n, flanks[1].p0, tv,
                                      p.maxRefineIters, hyp, ax, c, r);
                    if (tok) {
                        tf.axis = ax;
                        tf.center = c;
                        tf.radius = r;
                    }
                } else if (anchored && hem) {
                    gp_XYZ ax, c;
                    double r = 0;
                    tok = fitHem(flanks[0].n, flanks[0].p0, flanks[1].n, flanks[1].p0, tv, hyp,
                                 p.maxRefineIters, ax, c, r);
                    if (tok) {
                        tf.axis = ax;
                        tf.center = c;
                        tf.radius = r;
                    }
                } else if ((int)tv.size() >= 6) {
                    gp_Dir axOut;
                    gp_Pnt cOut;
                    double rOut = 0, mOut = 0;
                    gp_Dir seedDir(fit.axis.X(), fit.axis.Y(), fit.axis.Z());
                    tok = cylinderFitLS(mv, trial, seedDir, axOut, cOut, rOut, mOut,
                                        p.maxRefineIters);
                    if (tok) {
                        tf.axis = axOut.XYZ();
                        alignTo(tf.axis, hyp);
                        tf.center = cOut.XYZ();
                        tf.radius = rOut;
                    }
                }
                if (!tok) break;
                double rms = 0;
                const double mx = devStats(tv, tf.center, tf.axis, tf.radius, rms);
                if (mx > tau) break;
                band.swap(trial);
                verts.swap(tv);
                fit = tf;
                fit.maxDev = mx;
                fit.rms = rms;
                inBand.assign(nTri, 0);
                for (int t : band) inBand[t] = 1;
                dg.nTri = (int)band.size();
                dg.nV = (int)verts.size();
                dg.radius = fit.radius;
                dg.maxDev = mx;
                dg.axis = fit.axis;
            }

            // Same-surface: an adjacent accepted cylinder whose vertices all lie on this fit.
            bool same = false;
            for (const Region& acc : work.accepted) {
                if (acc.type != SurfType::Cylinder || acc.tris.empty()) continue;
                bool adj = false;
                for (int t : acc.tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3 && !adj; ++s) {
                        const int ot = otherTri(mv.triEdges[t][s], t);
                        if (ot >= 0 && ot < nTri && inBand[ot]) adj = true;
                    }
                    if (adj) break;
                }
                if (!adj) continue;
                std::vector<Vert> av;
                uniqueVerts(mv, g2l, acc.tris, av);
                if (av.empty()) continue;
                bool all = true;
                for (const Vert& v : av) {
                    if (std::abs(radialDist(v.p, fit.center, fit.axis) - fit.radius) > tau) {
                        all = false;
                        break;
                    }
                }
                if (all) {
                    same = true;
                    break;
                }
            }
            if (same) {
                dg.radius = fit.radius;
                dg.maxDev = fit.maxDev;
                dg.axis = fit.axis;
                finish("same-surface");
                continue;
            }

            // Blend classification against the two generator flanks.
            bool blend = false;
            double dih = 0, span = 0, stol = 0;
            if (flanks.size() == 2 && fit.radius > tau) {
                dih = std::acos(std::max(-1.0, std::min(1.0, flanks[0].n.Dot(flanks[1].n))));
                const gp_XYZ xdir = lowestXDir(verts, fit.center, fit.axis);
                const SpanInfo sp = measureSpan(mv, g2l, band, verts, fit.center, fit.axis,
                                                xdir, tol.thetaBin);
                span = sp.span;
                stol = spanTol(flanks[0].tol, flanks[1].tol, fit.radius);
                const double dA = std::abs((fit.center - flanks[0].p0).Dot(flanks[0].n));
                const double dB = std::abs((fit.center - flanks[1].p0).Dot(flanks[1].n));
                dg.tanDistA = std::abs(dA - fit.radius);
                dg.tanDistB = std::abs(dB - fit.radius);
                dg.tanAngA = std::abs(fit.axis.Dot(flanks[0].n));
                dg.tanAngB = std::abs(fit.axis.Dot(flanks[1].n));
                dg.span = span;
                dg.dihedral = dih;
                dg.spanTolV = stol;
                const double angTolA = tanAngTol(flanks[0].tol, tau, Lband);
                const double angTolB = tanAngTol(flanks[1].tol, tau, Lband);
                blend = std::abs(span - dih) <= stol &&
                        dg.tanDistA <= tanDistTol(flanks[0].tol, tau) &&
                        dg.tanDistB <= tanDistTol(flanks[1].tol, tau) &&
                        dg.tanAngA <= angTolA && dg.tanAngB <= angTolB;
            }

            const gp_XYZ xdir = lowestXDir(verts, fit.center, fit.axis);
            SpanInfo sp = measureSpan(mv, g2l, band, verts, fit.center, fit.axis, xdir,
                                      tol.thetaBin);
            if (!blend && !sp.closed360) {
                dg.span = sp.span;
                finish(flanks.size() < 2 ? "no-anchor" : "not-blend");
                continue;
            }

            const bool asFillet = blend;
            // Area-weighted centroid, triangles in ascending id (already sorted).
            gp_XYZ cent(0, 0, 0);
            double areaSum = 0.0;
            for (int t : band) {
                const double at = triArea(mv, t);
                areaSum += at;
                gp_XYZ ctri = (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) /
                              3.0;
                cent += ctri * at;
            }
            if (areaSum > 0.0) cent /= areaSum;
            gp_XYZ loc = fit.center + fit.axis * ((cent - fit.center).Dot(fit.axis));

            double sigma = 0.0;
            for (int t : band) {
                gp_XYZ nt;
                const double at = triAreaNormal(mv, t, nt);
                gp_XYZ rho = (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) /
                                 3.0 -
                             loc;
                rho -= fit.axis * rho.Dot(fit.axis);
                const double rm = rho.Modulus();
                if (rm < 1e-15) continue;
                rho /= rm;
                sigma += at * nt.Dot(rho);
            }
            const bool outward = sigma > 0.0;
            const double useSpan = sp.closed360 ? kTwoPi : std::max(sp.span, 1e-12);
            int nBands = std::max(sp.nBands, 1);
            int nSides = (int)std::llround(kTwoPi * (double)nBands / useSpan);
            if (nSides < 1) nSides = 1;

            double rms = 0;
            const double mx = devStats(verts, loc, fit.axis, fit.radius, rms);

            int straddles = 0;
            for (int pi = 0; pi < nProv; ++pi) {
                int nin = 0, nout = 0;
                for (int t : work.provisionals[pi].tris) {
                    if (t >= 0 && t < nTri && inBand[t]) ++nin;
                    else ++nout;
                }
                if (nin > 0 && nout > 0) ++straddles;
            }

            Region r;
            r.id = 0;
            r.type = SurfType::Cylinder;
            r.origin = asFillet ? Origin::FilletStrip : Origin::CylGrow;
            r.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(fit.axis), gp_Dir(xdir));
            r.radius = fit.radius;
            r.uMin = sp.uMin;
            r.uMax = sp.uMax;
            r.vMin = sp.vMin;
            r.vMax = sp.vMax;
            r.closed360 = sp.closed360;
            r.outwardNormal = outward;
            r.tris = band;
            r.maxVertexDev = mx;
            r.rmsVertexDev = rms;
            r.nSides = nSides;
            r.chordSagitta = chordSagitta(fit.radius, nSides);
            r.dVolPredicted = dVolCylinderSector(areaSum, fit.radius, nSides, outward);
            r.maxVertexSnap = 0.0;
            r.reject = Reject::None;
            r.builtAs = BuiltAs::NotBuilt;
            r.filletNbrA = (asFillet && flanks.size() == 2) ? flanks[0].prov : -1;
            r.filletNbrB = (asFillet && flanks.size() == 2) ? flanks[1].prov : -1;

            const ProvClaim stamp = asFillet ? ProvClaim::ConsumedFillet
                                             : ProvClaim::ConsumedCylinder;
            // Triangle granularity: certified triangles leave their provisional.
            // A straddling remainder keeps A2's plane fit and splits into components.
            std::vector<int> fresh;
            for (int pi = 0; pi < nProv; ++pi) {
                Provisional& pr = work.provisionals[pi];
                std::vector<int> keep, drop;
                for (int t : pr.tris) {
                    if (t >= 0 && t < nTri && inBand[t]) drop.push_back(t);
                    else keep.push_back(t);
                }
                if (drop.empty()) continue;
                if (keep.empty()) {
                    pr.claim = stamp;
                    continue;
                }
                pr.tris.swap(keep);
                refitPlaneProv(mv, g2l, pr);
                // Connected components of the remainder, ascending min tri.
                std::vector<char> mark(pr.tris.size(), 0);
                std::vector<std::vector<int>> comps;
                for (std::size_t i = 0; i < pr.tris.size(); ++i) {
                    if (mark[i]) continue;
                    std::vector<int> comp;
                    std::vector<std::size_t> st = {i};
                    mark[i] = 1;
                    while (!st.empty()) {
                        const std::size_t u = st.back();
                        st.pop_back();
                        comp.push_back(pr.tris[u]);
                        const int tu = pr.tris[u];
                        for (int s = 0; s < 3; ++s) {
                            const int ot = otherTri(mv.triEdges[tu][s], tu);
                            if (ot < 0) continue;
                            for (std::size_t j = 0; j < pr.tris.size(); ++j) {
                                if (!mark[j] && pr.tris[j] == ot) {
                                    mark[j] = 1;
                                    st.push_back(j);
                                }
                            }
                        }
                    }
                    std::sort(comp.begin(), comp.end());
                    comps.push_back(std::move(comp));
                }
                std::sort(comps.begin(), comps.end(),
                          [](const std::vector<int>& a, const std::vector<int>& b) {
                              return a.front() < b.front();
                          });
                pr.tris = comps[0];
                pr.claim = ProvClaim::Unclaimed;
                refitPlaneProv(mv, g2l, pr);
                for (std::size_t c = 1; c < comps.size(); ++c) fresh.push_back((int)c);
                // Stash extra components on the provisional by appending after the loop
                // via a side list carried on `fresh` as encoded sizes — rebuilt below.
                // Re-store components in a local that outlives this iteration.
                // (handled by copying into side storage immediately)
                for (std::size_t c = 1; c < comps.size(); ++c) {
                    Provisional np = pr;
                    np.tris = comps[c];
                    np.claim = ProvClaim::Unclaimed;
                    refitPlaneProv(mv, g2l, np);
                    work.provisionals.push_back(std::move(np));
                }
            }
            (void)fresh;

            for (int t : band)
                if (t >= 0 && t < nTri) taken[t] = 1;

            dg.radius = r.radius;
            dg.maxDev = r.maxVertexDev;
            dg.axis = fit.axis;
            dg.span = sp.span;
            dg.straddle = straddles;
            dg.nTri = (int)band.size();
            dg.nV = (int)verts.size();
            profileCounts(verts, fit.axis, tau, dg.stations, dg.levels, dg.lineRes);
            finish("accept");
            work.accepted.push_back(std::move(r));
        }

        if (diagOn()) {
            std::fprintf(stderr,
                         "DIAG_REL_SUM stage=extrude bands=%d accept=%d line=%d "
                         "vertex-off=%d not-blend=%d same-surface=%d underdetermined=%d "
                         "degenerate=%d no-anchor=%d\n",
                         yc.bands, yc.accept, yc.line, yc.vertexOff, yc.notBlend,
                         yc.sameSurface, yc.underdetermined, yc.degenerate, yc.noAnchor);
        }
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
