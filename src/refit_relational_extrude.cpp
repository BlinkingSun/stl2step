// Relational extrude (W1.1, amended D-fc-3.2).
// Anchor-pair vote, slab τ-membership, chord-bounded coverage, then pegs.
// Hems (D-fc-5.11 F9) are built by the one-sided station walk, which finds
// the far plane. Growth and the ruled-chain seed are not used.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_profile_cert.hpp"
#include "refit_relational_arc.hpp"
#include "parallel.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
// Rounding bound of the forOnCircle superset test, in ulps of the coordinate
// magnitude. Not a geometric tolerance.
constexpr double kUlp = 256.0 * std::numeric_limits<double>::epsilon();

bool diagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

// A pair evaluated on a worker buffers its DIAG lines here. The pair walk
// writes each buffer in pair order, so stderr is the serial walk's.
thread_local std::string* tDiagSink = nullptr;

struct DiagSinkScope {
    std::string* prev;
    explicit DiagSinkScope(std::string* s) : prev(tDiagSink) { tDiagSink = s; }
    ~DiagSinkScope() { tDiagSink = prev; }
    DiagSinkScope(const DiagSinkScope&) = delete;
    DiagSinkScope& operator=(const DiagSinkScope&) = delete;
};

void diagPrintf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    if (tDiagSink) {
        va_list ap2;
        va_copy(ap2, ap);
        const int n = std::vsnprintf(nullptr, 0, fmt, ap2);
        va_end(ap2);
        if (n > 0) {
            const std::size_t old = tDiagSink->size();
            tDiagSink->resize(old + (std::size_t)n + 1);
            std::vsnprintf(&(*tDiagSink)[old], (std::size_t)n + 1, fmt, ap);
            tDiagSink->resize(old + (std::size_t)n);
        }
    } else {
        std::vfprintf(stderr, fmt, ap);
    }
    va_end(ap);
}

// REV ORCH-3. A wave-1 claim exists only when both plane anchors are
// anchorCertified (D-fc-1.4), so the fit has one free parameter. Free-mode
// fits and the peg walk stay compiled below and are not called. The predicate
// is the anchor state, not a per-part switch.
bool extrudeAnchoredPair(bool certA, bool certB) { return certA && certB; }

// ORCH-3. A chamfer has one facet normal, so its separation is 0 and it does
// not span the dihedral. A fillet's facets turn by more than spanTol.
// Reaching the dihedral itself is the vertex-span test (|span - dihedral| <=
// spanTol): discrete facets stop one step short of the walls, and that step
// is the mesh angle, which spanTol (a function of anchorTol/r) does not cover.
bool normalsSpanDihedral(double normSep, double spanTolV) { return normSep > spanTolV; }

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

// D-fc-4.1. The provisional's maximum-area facet, ties to the smallest local
// id; -1 when no facet has positive area.
int maxAreaFacet(const MeshView& mv, const std::vector<int>& tris, int nTri) {
    int best = -1;
    double bestA = 0.0;
    for (int t : tris) {
        if (t < 0 || t >= nTri) continue;
        const double a = triArea(mv, t);
        if (a > bestA || (a == bestA && best >= 0 && t < best)) {
            bestA = a;
            best = t;
        }
    }
    return best;
}

// D-fc-4.3. An accepted CylGrow region whose own maxVertexDev exceeds τ does
// not certify its surface (D-fc-1.4), so a certified relational band may rob it.
bool uncertifiedCylGrow(const Region& r, double tau) {
    return r.origin == Origin::CylGrow && r.maxVertexDev > tau;
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
    int nI = 0;
    int nCover = 0;
    int compsI = 0;
    int compsC = 0;
    double certDev = 0;
    double coverDev = 0;
    double sDev = 0;
    double anchorTolA = 0;
    double anchorTolB = 0;
};

void emitDiag(const BandDiag& d) {
    if (!diagOn()) return;
    const double tau = d.tau > 0.0 ? d.tau : 1.0;
    diagPrintf("DIAG_REL_BAND stage=extrude mode=%s seed=%d nTri=%d nV=%d "
               "stations=%d levels=%d params=%d fit=%s radius=%.9g halfAngle=0 "
               "axis=%.9g,%.9g,%.9g maxDevMM=%.9g maxDevTau=%.9g lineResTau=%.9g "
               "span=%.9g dihedral=%.9g spanTol=%.9g tanDistA=%.9g tanDistB=%.9g "
               "tanAngA=%.9g tanAngB=%.9g straddle=%d "
               "nI=%d nCover=%d compsI=%d compsC=%d certDevTau=%.9g coverDevTau=%.9g "
               "sTau=%.9g anchorTolATau=%.9g anchorTolBTau=%.9g why=%s\n",
               d.mode, d.seed, d.nTri, d.nV, d.stations, d.levels, d.params, d.fit,
               d.radius, d.axis.X(), d.axis.Y(), d.axis.Z(), d.maxDev, d.maxDev / tau,
               d.lineRes / tau, d.span, d.dihedral, d.spanTolV, d.tanDistA, d.tanDistB,
               d.tanAngA, d.tanAngB, d.straddle, d.nI, d.nCover, d.compsI, d.compsC,
               d.certDev / tau, d.coverDev / tau, d.sDev / tau, d.anchorTolA / tau,
               d.anchorTolB / tau, d.why);
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
    int unstable = 0;
    int hemDeferred = 0;
    int anchorInBand = 0;
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
    else if (std::strcmp(why, "unstable") == 0) ++c.unstable;
    else if (std::strcmp(why, "hem-deferred") == 0) ++c.hemDeferred;
    else if (std::strcmp(why, "anchor-in-band") == 0) ++c.anchorInBand;
    else ++c.degenerate;
}

void addWhy(WhyCount& a, const WhyCount& b) {
    a.bands += b.bands;
    a.accept += b.accept;
    a.line += b.line;
    a.vertexOff += b.vertexOff;
    a.notBlend += b.notBlend;
    a.sameSurface += b.sameSurface;
    a.underdetermined += b.underdetermined;
    a.degenerate += b.degenerate;
    a.noAnchor += b.noAnchor;
    a.unstable += b.unstable;
    a.hemDeferred += b.hemDeferred;
    a.anchorInBand += b.anchorInBand;
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
                 const std::vector<Vert>& verts, int maxIter, gp_XYZ ruling, double sA,
                 double sB, gp_XYZ& axis, gp_XYZ& center, double& radius) {
    if (verts.empty()) return false;
    gp_XYZ a = nA.Crossed(nB);
    const double am = a.Modulus();
    if (!(am > 1e-12)) return false;
    a.Divide(am);
    alignTo(a, ruling);
    gp_XYZ C(0, 0, 0);
    for (const Vert& v : verts) C += v.p;
    C /= (double)verts.size();
    if (!(std::abs(sA) > 0.0) || !(std::abs(sB) > 0.0)) return false;
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
        // A tangent generator sits on a double root. Cancellation leaves a
        // disc a few ulps below zero (case 10, 4° planes). That is rounding,
        // not a geometric tolerance.
        const double discRaw = B * B - 4.0 * A * Cq;
        const double discScale = std::max(B * B, std::abs(4.0 * A * Cq));
        const double discFloor =
            64.0 * std::numeric_limits<double>::epsilon() * std::max(discScale, 1.0);
        if (discRaw < -discFloor) continue;
        const double disc = std::max(0.0, discRaw);
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

// D-fc-1.4: fillet facet normals must span the anchor dihedral (chamfer ⇒ one normal).
double maxFacetNormalSep(const MeshView& mv, const std::vector<int>& tris) {
    std::vector<gp_XYZ> ns;
    ns.reserve(tris.size());
    for (int t : tris) {
        gp_XYZ nt;
        triAreaNormal(mv, t, nt);
        const double nm = nt.Modulus();
        if (!(nm > 1e-15)) continue;
        nt.Divide(nm);
        ns.push_back(nt);
    }
    double best = 0.0;
    for (std::size_t i = 0; i < ns.size(); ++i) {
        for (std::size_t j = i + 1; j < ns.size(); ++j) {
            const double d = std::max(-1.0, std::min(1.0, ns[i].Dot(ns[j])));
            best = std::max(best, std::acos(d));
        }
    }
    return best;
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

namespace {

const char* originName(Origin o) {
    switch (o) {
    case Origin::PlaneGrow: return "PlaneGrow";
    case Origin::CylGrow: return "CylGrow";
    case Origin::FilletStrip: return "FilletStrip";
    case Origin::NgonWall: return "NgonWall";
    case Origin::ChamferCone: return "ChamferCone";
    case Origin::TorusBlend: return "TorusBlend";
    }
    return "other";
}

void canonicalHemi(gp_XYZ& a) {
    if (a.X() < -1e-15) {
        a.Reverse();
        return;
    }
    if (std::abs(a.X()) <= 1e-15 && a.Y() < -1e-15) {
        a.Reverse();
        return;
    }
    if (std::abs(a.X()) <= 1e-15 && std::abs(a.Y()) <= 1e-15 && a.Z() < 0.0) a.Reverse();
}

double wrap2pi(double a) {
    const double t = std::fmod(a, kTwoPi);
    return (t < 0.0) ? t + kTwoPi : t;
}

// Ascending sort of positive finite doubles. Their IEEE-754 bit patterns order
// as the values do, so this is std::sort's sequence: LSD radix on the top four
// bytes, then an insertion pass that orders each run of equal top bytes by the
// full pattern.
void sortPositive(std::vector<double>& v, std::vector<std::uint64_t>& a,
                  std::vector<std::uint64_t>& b) {
    const std::size_t n = v.size();
    if (n < 2) return;
    a.resize(n);
    b.resize(n);
    std::memcpy(a.data(), v.data(), n * sizeof(double));
    std::size_t cnt[4][256];
    std::memset(cnt, 0, sizeof(cnt));
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint64_t k = a[i];
        for (int d = 0; d < 4; ++d) ++cnt[d][(k >> (32 + 8 * d)) & 0xff];
    }
    for (int d = 0; d < 4; ++d) {
        std::size_t* c = cnt[d];
        bool trivial = false;  // every key has this byte: the pass is the identity
        for (int x = 0; x < 256 && !trivial; ++x) trivial = c[x] == n;
        if (trivial) continue;
        std::size_t sum = 0;
        for (int x = 0; x < 256; ++x) {
            const std::size_t t = c[x];
            c[x] = sum;
            sum += t;
        }
        const int sh = 32 + 8 * d;
        for (std::size_t i = 0; i < n; ++i) b[c[(a[i] >> sh) & 0xff]++] = a[i];
        a.swap(b);
    }
    for (std::size_t i = 1; i < n; ++i) {
        const std::uint64_t k = a[i];
        std::size_t j = i;
        while (j > 0 && a[j - 1] > k) {
            a[j] = a[j - 1];
            --j;
        }
        a[j] = k;
    }
    std::memcpy(v.data(), a.data(), n * sizeof(double));
}

void pushCircleRoots(const gp_XYZ& p, const gp_XYZ& c0, const gp_XYZ& cR, const gp_XYZ& axis,
                     double tau, std::vector<double>& roots) {
    const gp_XYZ w = (p - c0).Crossed(axis);
    const gp_XYZ dax = cR.Crossed(axis);
    const double A = dax.Dot(dax) - 1.0;
    const double B = -2.0 * w.Dot(dax);
    const double Cq = w.Dot(w);
    auto consider = [&](double r) {
        if (!(r > 0.0) || !std::isfinite(r)) return;
        const gp_XYZ c = c0 + cR * r;
        if (std::abs(radialDist(p, c, axis) - r) <= tau) roots.push_back(r);
    };
    if (std::abs(A) < 1e-12) {
        if (std::abs(B) > 1e-12) consider(-Cq / B);
        return;
    }
    // Same ulp clamp as fitAnchored: a tangent double root cancels to a
    // disc a few ulps negative. Not a geometric tolerance.
    const double discRaw = B * B - 4.0 * A * Cq;
    const double discScale = std::max(B * B, std::abs(4.0 * A * Cq));
    const double discFloor =
        64.0 * std::numeric_limits<double>::epsilon() * std::max(discScale, 1.0);
    if (discRaw < -discFloor) return;
    const double disc = std::max(0.0, discRaw);
    const double s = std::sqrt(disc);
    consider((-B + s) / (2.0 * A));
    consider((-B - s) / (2.0 * A));
}

// D-fc-3.2.2. The radii r > 0 of the side's circle family (centre c0 + cR·r,
// radius r, about axis) whose circle, inflated by τ, contains p:
// |(p − c(r)) × axis| <= r + τ. Every r whose circle passes within τ of p is
// in it. Aq = |cR × axis|² − 1 > 0. The discriminant carries kUlp of its terms'
// magnitude and the ends kUlp of theirs, so the window is a superset of the
// exact one (a tangent double root included). False when it is empty.
bool radiusWindow(const gp_XYZ& p, const gp_XYZ& c0, const gp_XYZ& dax, double Aq,
                  const gp_XYZ& axis, double tau, double& lo, double& hi) {
    const gp_XYZ w = (p - c0).Crossed(axis);
    const double b = w.Dot(dax) + tau;
    const double C = w.Dot(w) - tau * tau;
    const double mag = std::max({b * b, std::abs(Aq * C), w.Dot(w) * (Aq + 1.0)});
    const double disc = b * b - Aq * C + kUlp * mag;
    if (!(disc >= 0.0)) return false;
    const double s = std::sqrt(disc);
    lo = (b - s) / Aq;
    hi = (b + s) / Aq;
    const double pad = kUlp * (std::abs(lo) + std::abs(hi));
    lo -= pad;
    hi += pad;
    return hi > 0.0;
}

// D-fc-4.9. Least squares of the one-sided blend's circle in section: centre
// b + g·t + nT·(s·r), tangent to the line through b along g (the anchor T),
// kParamsOneSided = 2 free (t, r). Gauss-Newton from the caller's (t, r),
// capped by maxIter; a step that raises the residual ends it.
bool fitOneSided(const std::vector<Vert>& verts, const gp_XYZ& axis, const gp_XYZ& b,
                 const gp_XYZ& g, const gp_XYZ& nT, double s, int maxIter, double& t,
                 double& r) {
    if (verts.empty() || !(r > 0.0)) return false;
    auto sse = [&](double tp, double rr) {
        const gp_XYZ c = b + g * tp + nT * (s * rr);
        double e = 0.0;
        for (const Vert& v : verts) {
            const double f = radialDist(v.p, c, axis) - rr;
            e += f * f;
        }
        return e;
    };
    double bestE = sse(t, r);
    const int itCap = std::max(1, maxIter);
    for (int it = 0; it < itCap; ++it) {
        double JtJ00 = 0, JtJ01 = 0, JtJ11 = 0, Jte0 = 0, Jte1 = 0;
        const gp_XYZ c = b + g * t + nT * (s * r);
        for (const Vert& v : verts) {
            gp_XYZ d = v.p - c;
            d -= axis * d.Dot(axis);
            const double rho = d.Modulus();
            if (!(rho > 0.0)) continue;
            const double e = rho - r;
            const double j0 = -d.Dot(g) / rho;
            const double j1 = -s * d.Dot(nT) / rho - 1.0;
            JtJ00 += j0 * j0;
            JtJ01 += j0 * j1;
            JtJ11 += j1 * j1;
            Jte0 += j0 * e;
            Jte1 += j1 * e;
        }
        const double det = JtJ00 * JtJ11 - JtJ01 * JtJ01;
        if (!(det > 0.0)) break;
        const double dt = -(JtJ11 * Jte0 - JtJ01 * Jte1) / det;
        const double dr = -(JtJ00 * Jte1 - JtJ01 * Jte0) / det;
        if (!std::isfinite(dt) || !std::isfinite(dr) || !(r + dr > 0.0)) break;
        const double e2 = sse(t + dt, r + dr);
        if (!(e2 <= bestE)) break;
        t += dt;
        r += dr;
        bestE = e2;
    }
    return r > 0.0 && std::isfinite(r) && std::isfinite(t);
}

}  // namespace

bool claimRelExtrude(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work) {
    try {
        const auto tStart = std::chrono::steady_clock::now();
        const double tau = tauOf(mv);
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;

        const int nTri = (int)mv.nTri;
        const int nVtx = (int)mv.nVtx;
        const int nEdge = (int)mv.nEdge;
        int nProv = (int)work.provisionals.size();

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
        for (const Region& acc : work.accepted) {
            for (int t : acc.tris) {
                if (t >= 0 && t < nTri) taken[t] = 1;
            }
        }
        // D-fc-4.3. A CylGrow region whose own maxVertexDev exceeds τ is not a
        // certified surface. Its triangles are candidates of the vote, and a
        // certified band takes those that are I or Cover of its own surface
        // (D-fc-3.2.7). The robbed region dissolves at commit. A region within
        // τ is never touched.
        std::vector<int> uncertOwner((std::size_t)nTri, -1);
        std::vector<int> dissolved;
        for (int ri = 0; ri < (int)work.accepted.size(); ++ri) {
            const Region& acc = work.accepted[(std::size_t)ri];
            if (!uncertifiedCylGrow(acc, tau)) continue;
            for (int t : acc.tris) {
                if (t < 0 || t >= nTri) continue;
                uncertOwner[(std::size_t)t] = ri;
                taken[t] = 0;
            }
        }

        auto candidate = [&](int t) -> bool {
            if (t < 0 || t >= nTri || taken[t]) return false;
            if (uncertOwner[(std::size_t)t] >= 0) return true;
            const int pi = triProv[t];
            if (pi < 0 || pi >= (int)work.provisionals.size()) return false;
            return work.provisionals[pi].claim == ProvClaim::Unclaimed;
        };

        auto triLevelCount = [&](int t, const gp_XYZ& axis) -> int {
            std::vector<double> z(3);
            for (int k = 0; k < 3; ++k) z[(std::size_t)k] = cornerPt(mv, t, k).Dot(axis);
            return countLevels(z, tau);
        };

        struct ARef {
            int prov;
            int minTri;
        };
        std::vector<ARef> arefs;
        for (int pi = 0; pi < nProv; ++pi) {
            const Provisional& pr = work.provisionals[pi];
            if (pr.claim != ProvClaim::Unclaimed || pr.tris.empty()) continue;
            if (!a3WouldCommit(mv, tol, pr, false)) continue;
            ARef ar;
            ar.prov = pi;
            ar.minTri = minTriOf(pr.tris);
            arefs.push_back(ar);
        }
        std::sort(arefs.begin(), arefs.end(), [](const ARef& a, const ARef& b) {
            if (a.minTri != b.minTri) return a.minTri < b.minTri;
            return a.prov < b.prov;
        });

        struct LiveAnchor {
            int prov = -1;
            int minTri = -1;
            gp_XYZ n;
            gp_XYZ p0;
            double maxDev = 0;
            double tolA = 0;
            double extent = 0;
            bool certified = false;
            bool peeled = false;
            int nAll = 0;
            std::vector<int> body;
            std::vector<gp_XYZ> boundary;
        };
        // Reloaded after every claim. A committed band rewrites provisionals.
        int anchorEpoch = 1;
        std::vector<int> anchorEpochSeen;
        std::vector<LiveAnchor> anchorStored;
        std::vector<char> anchorStoredOk;
        int nInlierAnchors = 0;
        std::vector<char> inlierCounted;

        auto triOnPlane = [&](int t, const gp_XYZ& n, const gp_XYZ& p0) -> bool {
            if (t < 0 || t >= nTri) return false;
            for (int k = 0; k < 3; ++k) {
                if (std::abs((cornerPt(mv, t, k) - p0).Dot(n)) > tau) return false;
            }
            return true;
        };
        auto fillBoundary = [&](LiveAnchor& A) {
            A.boundary.clear();
            if (A.body.empty()) return;
            std::vector<char> inB(nTri, 0);
            for (int t : A.body)
                if (t >= 0 && t < nTri) inB[t] = 1;
            std::vector<int> seen(nVtx, 0);
            for (int t : A.body) {
                if (t < 0 || t >= nTri) continue;
                for (int s = 0; s < 3; ++s) {
                    const int ot = otherTri(mv.triEdges[t][s], t);
                    if (ot >= 0 && ot < nTri && inB[ot]) continue;
                    for (int dk = 0; dk < 2; ++dk) {
                        const int k = (s + dk) % 3;
                        const int lv = cornerLocal(g2l, mv, t, k);
                        if (lv < 0 || lv >= nVtx || seen[lv]) continue;
                        seen[lv] = 1;
                        const int g = mv.compVtx ? mv.compVtx[lv] : lv;
                        A.boundary.push_back(mv.pts[g]);
                    }
                }
            }
        };
        // REV ORCH-4. A flank whose stored maxVertexDev exceeds τ often still
        // owns a wall: the extra is fillet facets. The anchor body is the
        // τ-preimage of the plane, refit until that preimage is stable.
        // anchorCertified is unchanged (body maxVertexDev ≤ τ). A provisional
        // that is already entirely within τ keeps its stored deviation, so a
        // caller-supplied uncertified plane is not promoted (case 9).
        // Pure in the current provisionals: no cache, tally or DIAG, so an
        // epoch's anchors can be computed on worker threads. storeAnchor
        // records a load in walk order.
        auto computeAnchor = [&](int prov, LiveAnchor& A) -> bool {
            if (prov < 0 || prov >= (int)work.provisionals.size()) return false;
            const Provisional& pr = work.provisionals[prov];
            // The anchor predicate itself (audit-fx-arch F-5: one factor).
            PlaneAnchor pa;
            if (!planeAnchorOf(mv, p, tol, pr, pa)) return false;
            A = LiveAnchor();
            A.prov = prov;
            A.n = pa.n;
            A.p0 = pa.p0;
            A.body.swap(pa.body);
            A.minTri = minTriOf(A.body);
            A.maxDev = pa.dev;
            A.tolA = pa.tol;
            A.certified = pa.certified;
            A.peeled = pa.peeled;
            A.nAll = pa.peeled ? (int)pr.tris.size() : 0;
            std::vector<gp_XYZ> pts;
            pts.reserve(A.body.size() * 3);
            for (int t : A.body) {
                if (t < 0 || t >= nTri) continue;
                for (int k = 0; k < 3; ++k) pts.push_back(cornerPt(mv, t, k));
            }
            A.extent = sectionExtent(pts, A.n);
            if (!(A.extent > 0.0)) A.extent = tau;
            fillBoundary(A);
            return true;
        };
        // The first load of prov in an epoch: cache it, count a peeled inlier
        // body once, print DIAG_REL_PEEL.
        auto storeAnchor = [&](int prov, bool ok, const LiveAnchor& A) -> bool {
            anchorEpochSeen[(std::size_t)prov] = anchorEpoch;
            anchorStoredOk[(std::size_t)prov] = ok ? 1 : 0;
            if (!ok) return false;
            anchorStored[(std::size_t)prov] = A;
            if (A.peeled) {
                if ((int)inlierCounted.size() <= prov)
                    inlierCounted.resize((std::size_t)prov + 1, 0);
                if (!inlierCounted[(std::size_t)prov]) {
                    inlierCounted[(std::size_t)prov] = 1;
                    ++nInlierAnchors;
                }
                if (diagOn()) {
                    std::fprintf(stderr,
                                 "DIAG_REL_PEEL stage=extrude prov=%d minTri=%d nBody=%d "
                                 "nAll=%d maxDevTau=%.9g n=%.6g,%.6g,%.6g\n",
                                 prov, A.minTri, (int)A.body.size(), A.nAll,
                                 (tau > 0.0 ? A.maxDev / tau : 0.0), A.n.X(), A.n.Y(), A.n.Z());
                }
            }
            return true;
        };
        // Cache slot for prov; false when prov is out of range.
        auto anchorSlot = [&](int prov) -> bool {
            if (prov < 0 || prov >= (int)work.provisionals.size()) return false;
            if ((int)anchorStored.size() <= prov) {
                anchorStored.resize((std::size_t)prov + 1);
                anchorStoredOk.resize((std::size_t)prov + 1, 0);
                anchorEpochSeen.resize((std::size_t)prov + 1, 0);
            }
            return true;
        };
        auto loadAnchor = [&](int prov, LiveAnchor& A) -> bool {
            if (!anchorSlot(prov)) return false;
            if (anchorEpochSeen[(std::size_t)prov] == anchorEpoch) {
                if (!anchorStoredOk[(std::size_t)prov]) return false;
                A = anchorStored[(std::size_t)prov];
                return true;
            }
            const bool ok = computeAnchor(prov, A);
            return storeAnchor(prov, ok, A);
        };

        WhyCount yc;
        int nPairs = 0;
        int nEvalSides = 0;
        int nPrunedSides = 0;

        auto edgeComps = [&](const std::vector<int>& tris) -> int {
            std::vector<char> in(nTri, 0);
            for (int t : tris)
                if (t >= 0 && t < nTri) in[t] = 1;
            std::vector<char> seen(nTri, 0);
            int comps = 0;
            std::vector<int> order = tris;
            std::sort(order.begin(), order.end());
            for (int t : order) {
                if (t < 0 || t >= nTri || !in[t] || seen[t]) continue;
                ++comps;
                std::vector<int> st = {t};
                seen[t] = 1;
                while (!st.empty()) {
                    const int u = st.back();
                    st.pop_back();
                    for (int s = 0; s < 3; ++s) {
                        const int ot = otherTri(mv.triEdges[u][s], u);
                        if (ot < 0 || ot >= nTri || !in[ot] || seen[ot]) continue;
                        seen[ot] = 1;
                        st.push_back(ot);
                    }
                }
            }
            return comps;
        };

        auto commitBand = [&](std::vector<int> band, const gp_XYZ& axis, const gp_XYZ& center,
                              double radius, bool asFillet, int flankA, int flankB, int seed,
                              BandDiag& dg) {
            std::sort(band.begin(), band.end());
            band.erase(std::unique(band.begin(), band.end()), band.end());
            std::vector<char> inBand(nTri, 0);
            for (int t : band)
                if (t >= 0 && t < nTri) inBand[t] = 1;
            std::vector<Vert> verts;
            uniqueVerts(mv, g2l, band, verts);
            const gp_XYZ xdir = lowestXDir(verts, center, axis);
            SpanInfo sp = measureSpan(mv, g2l, band, verts, center, axis, xdir, tol.thetaBin);
            gp_XYZ cent(0, 0, 0);
            double areaSum = 0.0;
            for (int t : band) {
                const double at = triArea(mv, t);
                areaSum += at;
                const gp_XYZ ctri =
                    (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) / 3.0;
                cent += ctri * at;
            }
            if (areaSum > 0.0) cent /= areaSum;
            gp_XYZ loc = center + axis * ((cent - center).Dot(axis));
            double sigma = 0.0;
            for (int t : band) {
                gp_XYZ nt;
                const double at = triAreaNormal(mv, t, nt);
                gp_XYZ rho =
                    (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) / 3.0 - loc;
                rho -= axis * rho.Dot(axis);
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
            double rms = 0.0;
            const double mx = devStats(verts, loc, axis, radius, rms);
            int straddles = 0;
            const int nProvNow = (int)work.provisionals.size();
            for (int pi = 0; pi < nProvNow; ++pi) {
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
            r.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(axis), gp_Dir(xdir));
            r.radius = radius;
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
            r.chordSagitta = chordSagitta(radius, nSides);
            r.dVolPredicted = dVolCylinderSector(areaSum, radius, nSides, outward);
            r.maxVertexSnap = 0.0;
            r.reject = Reject::None;
            r.builtAs = BuiltAs::NotBuilt;
            r.filletNbrA = asFillet ? flankA : -1;
            r.filletNbrB = asFillet ? flankB : -1;

            const ProvClaim stamp =
                asFillet ? ProvClaim::ConsumedFillet : ProvClaim::ConsumedCylinder;
            // D-fc-4.3. Every uncertified region the band robbed dissolves. Its
            // own provisionals return to Unclaimed (the rewrite below then
            // drops the band's triangles from them), and its other triangles
            // leave it.
            std::vector<int> robbed;
            for (int t : band)
                if (t >= 0 && t < nTri && uncertOwner[(std::size_t)t] >= 0)
                    robbed.push_back(uncertOwner[(std::size_t)t]);
            std::sort(robbed.begin(), robbed.end());
            robbed.erase(std::unique(robbed.begin(), robbed.end()), robbed.end());
            for (int ri : robbed) {
                Region& rr = work.accepted[(std::size_t)ri];
                std::vector<char> inR(nTri, 0);
                int nStolen = 0;
                for (int t : rr.tris) {
                    if (t < 0 || t >= nTri) continue;
                    inR[t] = 1;
                    if (inBand[t]) ++nStolen;
                }
                for (Provisional& pr : work.provisionals) {
                    if (pr.claim != ProvClaim::ConsumedCylinder || pr.tris.empty()) continue;
                    bool own = true;
                    for (int t : pr.tris) {
                        if (t < 0 || t >= nTri || !inR[t]) {
                            own = false;
                            break;
                        }
                    }
                    if (own) pr.claim = ProvClaim::Unclaimed;
                }
                for (int t : rr.tris) {
                    if (t < 0 || t >= nTri) continue;
                    uncertOwner[(std::size_t)t] = -1;
                    if (!inBand[t]) taken[t] = 0;
                }
                if (diagOn())
                    diagPrintf("DIAG_REL_ROB stage=extrude seed=%d region=%d radius=%.9g "
                               "devTau=%.9g nTri=%d stolen=%d released=%d\n",
                               seed, ri, rr.radius, rr.maxVertexDev / tau,
                               (int)rr.tris.size(), nStolen, (int)rr.tris.size() - nStolen);
                rr.tris.clear();
                dissolved.push_back(ri);
            }
            for (int pi = 0; pi < nProvNow; ++pi) {
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
                for (std::size_t c = 1; c < comps.size(); ++c) {
                    Provisional np = pr;
                    np.tris = comps[c];
                    np.claim = ProvClaim::Unclaimed;
                    refitPlaneProv(mv, g2l, np);
                    work.provisionals.push_back(std::move(np));
                }
            }
            for (int t : band)
                if (t >= 0 && t < nTri) taken[t] = 1;
            dg.radius = r.radius;
            dg.maxDev = r.maxVertexDev;
            dg.axis = axis;
            dg.span = sp.span;
            dg.straddle = straddles;
            dg.nTri = (int)band.size();
            dg.nV = (int)verts.size();
            dg.seed = seed;
            profileCounts(verts, axis, tau, dg.stations, dg.levels, dg.lineRes);
            dg.why = "accept";
            emitDiag(dg);
            tallyWhy(yc, "accept");
            work.accepted.push_back(std::move(r));
            // The claim rewrites provisionals. Later pairs reload anchors.
            ++anchorEpoch;
        };

        auto emitAdj = [&](int seed, const gp_XYZ& axis, const gp_XYZ& center, double radius,
                           const std::vector<char>& inC) {
            if (!diagOn()) return;
            const int nAcc = (int)work.accepted.size();
            for (int i = 0; i < nAcc; ++i) {
                const Region& acc = work.accepted[i];
                if (acc.type != SurfType::Cylinder || acc.tris.empty()) continue;
                int nTouch = 0;
                for (int t : acc.tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3; ++s) {
                        const int ot = otherTri(mv.triEdges[t][s], t);
                        if (ot >= 0 && ot < nTri && inC[ot]) ++nTouch;
                    }
                }
                if (nTouch < 1) continue;
                gp_Dir axOut = acc.ax.Direction();
                gp_Pnt cOut = acc.ax.Location();
                double rOut = acc.radius, mOut = acc.maxVertexDev;
                const bool fitOk =
                    cylinderFitLS(mv, acc.tris, acc.ax.Direction(), axOut, cOut, rOut, mOut,
                                  std::max(1, p.maxRefineIters));
                const double dev = fitOk ? mOut : acc.maxVertexDev;
                int within = 1;
                double maxOff = 0.0;
                std::vector<Vert> av;
                uniqueVerts(mv, g2l, acc.tris, av);
                for (const Vert& v : av) {
                    const double off = std::abs(radialDist(v.p, center, axis) - radius);
                    if (off > maxOff) maxOff = off;
                    if (off > tau) within = 0;
                }
                diagPrintf("DIAG_REL_ADJ stage=extrude seed=%d adjIndex=%d origin=%s "
                           "radius=%.9g refitDevTau=%.9g withinTau=%d maxOffTau=%.9g nTri=%d\n",
                           seed, i, originName(acc.origin), acc.radius, dev / tau, within,
                           (tau > 0.0 ? maxOff / tau : 0.0), (int)acc.tris.size());
            }
        };

        auto adjacentAcceptedSurf = [&](int t) -> bool {
            if (t < 0 || t >= nTri) return false;
            for (int s = 0; s < 3; ++s) {
                const int ot = otherTri(mv.triEdges[t][s], t);
                if (ot < 0) continue;
                for (const Region& acc : work.accepted) {
                    if (acc.tris.empty()) continue;
                    if (acc.type != SurfType::Torus && acc.type != SurfType::Cone &&
                        acc.type != SurfType::Cylinder)
                        continue;
                    for (int at : acc.tris) {
                        if (at == ot) return true;
                    }
                }
            }
            return false;
        };

        auto sameSurface = [&](const gp_XYZ& axis, const gp_XYZ& center, double radius,
                               const std::vector<char>& inC) -> bool {
            for (const Region& acc : work.accepted) {
                if (acc.type != SurfType::Cylinder || acc.tris.empty()) continue;
                bool adj = false;
                for (int t : acc.tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3 && !adj; ++s) {
                        const int ot = otherTri(mv.triEdges[t][s], t);
                        if (ot >= 0 && ot < nTri && inC[ot]) adj = true;
                    }
                    if (adj) break;
                }
                if (!adj) continue;
                std::vector<Vert> av;
                uniqueVerts(mv, g2l, acc.tris, av);
                if (av.empty()) continue;
                bool all = true;
                for (const Vert& v : av) {
                    if (std::abs(radialDist(v.p, center, axis) - radius) > tau) {
                        all = false;
                        break;
                    }
                }
                if (all) return true;
            }
            return false;
        };

        struct CandV {
            int lv;
            gp_XYZ p;
        };
        std::vector<CandV> candPts;
        double candMag = 0.0;
        int candEpoch = -1;
        // Rebuilt when a claim changes the candidate set. One vertex list per
        // epoch, not one triangle scan per pair-side.
        auto rebuildCand = [&]() {
            candPts.clear();
            candMag = 0.0;
            std::vector<int> seen(nVtx, 0);
            candPts.reserve((std::size_t)nTri);
            for (int t = 0; t < nTri; ++t) {
                if (!candidate(t)) continue;
                for (int k = 0; k < 3; ++k) {
                    const int lv = cornerLocal(g2l, mv, t, k);
                    if (lv < 0 || lv >= nVtx || seen[lv]) continue;
                    seen[lv] = 1;
                    const int g = mv.compVtx ? mv.compVtx[lv] : lv;
                    CandV cv;
                    cv.lv = lv;
                    cv.p = mv.pts[g];
                    candPts.push_back(cv);
                    candMag = std::max(candMag, cv.p.Modulus());
                }
            }
            candEpoch = anchorEpoch;
        };

        // Triangle incidence of each local vertex. A τ-member triangle is
        // incident to a vertex within τ, so an empty I is known without walking
        // every triangle. The slab scan still runs when that set is non-empty,
        // because coverage includes sector triangles that are not within τ.
        std::vector<std::vector<int>> vertTris((std::size_t)nVtx);
        for (int t = 0; t < nTri; ++t) {
            for (int k = 0; k < 3; ++k) {
                const int lv = cornerLocal(g2l, mv, t, k);
                if (lv < 0 || lv >= nVtx) continue;
                vertTris[(std::size_t)lv].push_back(t);
            }
        }

        const double sideSign[4][2] = {{1.0, 1.0}, {1.0, -1.0}, {-1.0, 1.0}, {-1.0, -1.0}};
        // FilletStrips accepted before this call. A strip this stage commits
        // is not one of them (plate E1/E2 meet the fillet just claimed).
        const int nPriorAcc = (int)work.accepted.size();
        std::vector<int> priorFilTri((std::size_t)nTri, -1);
        for (int ai = 0; ai < nPriorAcc; ++ai) {
            const Region& acc = work.accepted[(std::size_t)ai];
            if (acc.origin != Origin::FilletStrip) continue;
            for (int t : acc.tris)
                if (t >= 0 && t < nTri) priorFilTri[(std::size_t)t] = ai;
        }
        // 0 unknown, 1 own refit within τ, 2 own refit exceeds τ. Workers share
        // it; a racing pair computes the same value.
        std::unique_ptr<std::atomic<signed char>[]> priorFit(
            new std::atomic<signed char>[(std::size_t)std::max(nPriorAcc, 1)]);
        for (int ai = 0; ai < std::max(nPriorAcc, 1); ++ai)
            priorFit[(std::size_t)ai].store(0, std::memory_order_relaxed);
        // D-fc-1.6.10: corner triangles stay faceted. A boundary edge of C that
        // is not a generator (edgeParallel, tanAngTol(τ, τ, len)) and that
        // meets a prior FilletStrip is a corner. D-fc-1.6.5 keeps that corner
        // when the strip's own cylinder refit is within τ. A strip that misses
        // τ is not a blend partner, and the band is not-blend. Generator
        // contact stays. Same-surface is the separate all-within-τ refusal.
        auto priorCornerUncertified = [&](const std::vector<int>& C,
                                           const std::vector<char>& inC,
                                           const gp_XYZ& axis) -> bool {
            for (int t : C) {
                if (t < 0 || t >= nTri) continue;
                for (int s = 0; s < 3; ++s) {
                    const int ot = otherTri(mv.triEdges[t][s], t);
                    if (ot < 0 || ot >= nTri || inC[ot]) continue;
                    const int ai = priorFilTri[(std::size_t)ot];
                    if (ai < 0) continue;
                    const gp_XYZ pa = cornerPt(mv, t, s);
                    const gp_XYZ pb = cornerPt(mv, t, (s + 1) % 3);
                    gp_XYZ dir = pb - pa;
                    const double len = dir.Modulus();
                    if (!(len > tau)) continue;
                    dir /= len;
                    if (edgeParallel(dir, len, axis, tau)) continue;
                    signed char fit = priorFit[(std::size_t)ai].load(std::memory_order_relaxed);
                    if (fit == 0) {
                        const Region& acc = work.accepted[(std::size_t)ai];
                        gp_Dir axOut = acc.ax.Direction();
                        gp_Pnt cOut = acc.ax.Location();
                        double rOut = acc.radius;
                        double mOut = acc.maxVertexDev;
                        const bool fitOk = cylinderFitLS(
                            mv, acc.tris, acc.ax.Direction(), axOut, cOut, rOut, mOut,
                            std::max(1, p.maxRefineIters));
                        const double dev = fitOk ? mOut : acc.maxVertexDev;
                        fit = (dev <= tau) ? 1 : 2;
                        priorFit[(std::size_t)ai].store(fit, std::memory_order_relaxed);
                    }
                    if (fit == 2) return true;
                }
            }
            return false;
        };

        // The pair walk. A pair's evaluation reads the claim state and writes
        // only its PairOut, and it stops at a commit. Pairs are evaluated in
        // blocks on worker threads; the walk then replays them in pair order:
        // anchor loads, DIAG text, tallies, and the first commit. A commit
        // changes the state, so the walk resumes that pair from (si, guard)
        // with the anchors it started with, as the serial loop did, and
        // re-evaluates every later pair. Output does not depend on the worker
        // count; threads 1 runs it inline.
        struct PairScratch {
            std::vector<int> vertGen, triGen;
            int vertTick = 0, triTick = 0;
            std::vector<double> offA, offB, roots;
            std::vector<std::uint64_t> keyA, keyB;
            std::vector<std::size_t> hit;
            std::vector<std::pair<double, double>> winA, winB;
            std::vector<std::size_t> onA, onB;
        };
        struct PairOut {
            std::string diag;
            WhyCount yc;
            int nPairs = 0;
            int nEvalSides = 0;
            int nPrunedSides = 0;
            bool commit = false;
            int si = 0;
            int guard = 0;
            std::vector<int> C;
            gp_XYZ axis, center;
            double radius = 0;
            int seed = -1;
            BandDiag dg;
        };
        // D-fc-3.2.2 lossless prune. A committed band is tangent to A and B,
        // and its contact vertices lie on both anchor planes: boundary
        // vertices of the anchor bodies (watertight), or, when the band meets
        // another piece of the same wall, candidate vertices within
        // tanDistTol of the anchor's plane. Either is within τ of the band's
        // circle. Every radius of this side that brings such a vertex of A
        // within τ lies in the union of A's radiusWindow; likewise B. When the
        // two unions do not meet, no radius of this side can be a claim: the
        // side is skipped. Superset windows only.
        auto sideMeetsAnchors = [&](const LiveAnchor& A, const LiveAnchor& B, const gp_XYZ& c0,
                                    const gp_XYZ& cR, const gp_XYZ& axis,
                                    PairScratch& sc) -> bool {
            const gp_XYZ dax = cR.Crossed(axis);
            const double Aq = dax.Dot(dax) - 1.0;
            if (!(Aq > 0.0) || !std::isfinite(Aq)) return true;
            auto windows = [&](const LiveAnchor& X, const std::vector<std::size_t>& onX,
                               std::vector<std::pair<double, double>>& win) {
                win.clear();
                double lo = 0, hi = 0;
                for (const gp_XYZ& p : X.boundary)
                    if (radiusWindow(p, c0, dax, Aq, axis, tau, lo, hi)) win.emplace_back(lo, hi);
                for (std::size_t i : onX)
                    if (radiusWindow(candPts[i].p, c0, dax, Aq, axis, tau, lo, hi))
                        win.emplace_back(lo, hi);
                std::sort(win.begin(), win.end());
            };
            windows(A, sc.onA, sc.winA);
            if (sc.winA.empty()) return false;
            windows(B, sc.onB, sc.winB);
            std::size_t i = 0, j = 0;
            while (i < sc.winA.size() && j < sc.winB.size()) {
                const auto& a = sc.winA[i];
                const auto& b = sc.winB[j];
                if (a.second < b.first) {
                    ++i;
                } else if (b.second < a.first) {
                    ++j;
                } else {
                    return true;
                }
            }
            return false;
        };
        auto runPair = [&](std::size_t ia, std::size_t ib, const LiveAnchor& A,
                           const LiveAnchor& B, bool resume, int startSi, int startGuard,
                           PairScratch& sc, PairOut& out) {
            std::vector<int>& vertGen = sc.vertGen;
            std::vector<int>& triGen = sc.triGen;
            int& vertTick = sc.vertTick;
            int& triTick = sc.triTick;
            if (!resume) ++out.nPairs;
            const int seed = std::min(arefs[ia].minTri, arefs[ib].minTri);
            const double cross = A.n.Crossed(B.n).Modulus();
            const double rhs = 2.0 * (A.tolA + tau) / A.extent + 2.0 * (B.tolA + tau) / B.extent;
            if (cross <= rhs) {
                // A parallel pair has no n_A × n_B axis. The hem (D-fc-5.11 F9)
                // takes its axis from a third certified plane, in the one-sided
                // station walk below; this walk leaves the pair.
                BandDiag dg;
                dg.tau = tau;
                dg.mode = "hem";
                dg.fit = (A.certified && B.certified) ? "anchored" : "free";
                dg.seed = seed;
                dg.anchorTolA = A.tolA;
                dg.anchorTolB = B.tolA;
                dg.why = "hem-deferred";
                emitDiag(dg);
                tallyWhy(out.yc, "hem-deferred");
                return;
            }
            gp_XYZ axis0 = A.n.Crossed(B.n);
            const double am = axis0.Modulus();
            if (!(am > 1e-12)) return;
            axis0.Divide(am);
            canonicalHemi(axis0);
            const bool anchored = extrudeAnchoredPair(A.certified, B.certified);
            // Free mode stays in the fit branch below and is not entered.
            if (!anchored) return;

            // Signed offsets of the candidate vertices from each anchor plane.
            // A circle about axis0 lies within r of its centre along A.n and
            // B.n (both normal to axis0), so a vertex within τ of it is within
            // r + τ there. forOnCircle tests that first and the τ test still
            // decides, in candPts order: a superset filter, not a tolerance.
            // The extra τ covers rounding, which kUlp·mag/am bounds; when it
            // does not, every vertex takes the τ test.
            const double offA = A.p0.Dot(A.n);
            const double offB = B.p0.Dot(B.n);
            const std::size_t nCand = candPts.size();
            sc.offA.resize(nCand);
            sc.offB.resize(nCand);
            const double onTolA = tanDistTol(A.tolA, tau);
            const double onTolB = tanDistTol(B.tolA, tau);
            sc.onA.clear();
            sc.onB.clear();
            for (std::size_t i = 0; i < nCand; ++i) {
                sc.offA[i] = A.n.Dot(candPts[i].p) - offA;
                sc.offB[i] = B.n.Dot(candPts[i].p) - offB;
                if (std::abs(sc.offA[i]) <= onTolA) sc.onA.push_back(i);
                if (std::abs(sc.offB[i]) <= onTolB) sc.onB.push_back(i);
            }
            auto forOnCircle = [&](const gp_XYZ& c, const gp_XYZ& ax, double r, auto&& fn) {
                const double mag = std::max(
                    {candMag, c.Modulus(), std::abs(offA), std::abs(offB), std::abs(r)});
                const bool sup = ax.X() == axis0.X() && ax.Y() == axis0.Y() &&
                                 ax.Z() == axis0.Z() &&
                                 kUlp * mag <= tau * am;
                const double hA = A.n.Dot(c) - offA;
                const double hB = B.n.Dot(c) - offB;
                const double lim = r + 2.0 * tau;
                for (std::size_t i = 0; i < nCand; ++i) {
                    if (sup && (std::abs(sc.offA[i] - hA) > lim || std::abs(sc.offB[i] - hB) > lim))
                        continue;
                    const CandV& cv = candPts[i];
                    if (std::abs(radialDist(cv.p, c, ax) - r) <= tau) fn(i, cv);
                }
            };

            for (int si = resume ? startSi : 0; si < 4; ++si) {
                const double sA = sideSign[si][0];
                const double sB = sideSign[si][1];
                gp_XYZ c0, c1;
                if (!solve3(A.n, B.n, axis0, A.p0.Dot(A.n), B.p0.Dot(B.n), 0.0, c0)) continue;
                if (!solve3(A.n, B.n, axis0, A.p0.Dot(A.n) + sA, B.p0.Dot(B.n) + sB, 0.0, c1))
                    continue;
                const gp_XYZ cR = c1 - c0;
                const bool resumed = resume && si == startSi;
                if (!resumed && !sideMeetsAnchors(A, B, c0, cR, axis0, sc)) {
                    ++out.nPrunedSides;
                    continue;
                }
                if (!resumed) ++out.nEvalSides;

                // pushCircleRoots keeps a root r of p only when p is within τ of
                // the circle about c0 + cR·r, which is tangent to both anchors on
                // this side, so p is then within τ of the side's quadrant. Every
                // root is at most rMax: the centre runs off at |cR × axis0| > 1
                // per unit r. When the solve3 residuals at rMax and the rounding
                // (kUlp) stay within τ, a vertex more than 2τ outside the quadrant
                // has no root and is not tried. Otherwise every vertex is tried.
                bool quadSkip = false;
                {
                    const double D = cR.Crossed(axis0).Modulus();
                    if (D > 1.0) {
                        const double c0m = c0.Modulus();
                        const double rMax = (candMag + c0m + 2.0 * tau) / (D - 1.0);
                        const double mag = std::max({candMag, c0m + rMax * cR.Modulus(), rMax,
                                                     std::abs(offA), std::abs(offB)});
                        const double slack =
                            kUlp * mag + 2.0 * mag *
                                             std::max(std::abs(A.n.Dot(axis0)),
                                                      std::abs(B.n.Dot(axis0)));
                        const double errA = std::abs(A.n.Dot(c0) - offA) +
                                            rMax * std::abs(A.n.Dot(cR) - sA) + slack;
                        const double errB = std::abs(B.n.Dot(c0) - offB) +
                                            rMax * std::abs(B.n.Dot(cR) - sB) + slack;
                        quadSkip = std::isfinite(rMax) && std::isfinite(mag) && errA <= tau &&
                                   errB <= tau;
                    }
                }

                // candPts is current: the walk rebuilds it before every call.
                int guard = resumed ? startGuard : 0;
                while (guard++ < nTri + 1) {
                    std::vector<double>& roots = sc.roots;
                    roots.clear();
                    for (std::size_t i = 0; i < nCand; ++i) {
                        if (quadSkip &&
                            (sA * sc.offA[i] < -2.0 * tau || sB * sc.offB[i] < -2.0 * tau))
                            continue;
                        pushCircleRoots(candPts[i].p, c0, cR, axis0, tau, roots);
                    }
                    const bool noRoots = roots.empty();
                    double rHyp = 0.0;
                    BandDiag dg;
                    dg.tau = tau;
                    dg.mode = "fillet";
                    dg.fit = anchored ? "anchored" : "free";
                    dg.params = anchored ? kParamsAnchored : (kParamsCircle + kParamsAxisDir);
                    dg.seed = seed;
                    dg.axis = axis0;
                    dg.anchorTolA = A.tolA;
                    dg.anchorTolB = B.tolA;
                    auto finish = [&](const char* why) {
                        dg.why = why;
                        emitDiag(dg);
                        tallyWhy(out.yc, why);
                    };
                    // Hypothesis, one per pair-side (D-fc-3.2.3, prototype _vote_radius).
                    // Admissible roots are clustered into 2τ windows. A window is scored
                    // by how many candidate vertices lie within τ of its median circle.
                    // The densest root cluster can be a spurious quadratic mode, so the
                    // score is the geometric inlier count, not the root count. The
                    // maximum score is the hypothesis. A later certificate refusal ends
                    // the side: the window is not replaced.
                    {
                        std::vector<double>& pool = roots;
                        sortPositive(pool, sc.keyA, sc.keyB);
                        bool any = false;
                        int bestIn = -1;
                        double bestR = 0.0;
                        // Five residual modes, as _vote_radius. Not a geometric tolerance:
                        // hypothesis selection only. A certificate refusal still ends the side.
                        for (int mode = 0; mode < 5; ++mode) {
                            const int n = (int)pool.size();
                            if (n <= 0) break;
                            int j = 0, bn = 0, bj = 0, bi = 0;
                            for (int i = 0; i < n; ++i) {
                                while (pool[(std::size_t)i] - pool[(std::size_t)j] > 2.0 * tau)
                                    ++j;
                                const int nw = i - j + 1;
                                if (nw > bn) {
                                    bn = nw;
                                    bj = j;
                                    bi = i;
                                }
                            }
                            if (bn <= 0) break;
                            const int mid = bj + bn / 2;
                            const double rmed =
                                (bn % 2) != 0
                                    ? pool[(std::size_t)mid]
                                    : 0.5 * (pool[(std::size_t)mid - 1] + pool[(std::size_t)mid]);
                            const gp_XYZ hc = c0 + cR * rmed;
                            int nIn = 0;
                            forOnCircle(hc, axis0, rmed,
                                        [&](std::size_t, const CandV&) { ++nIn; });
                            if (!any || nIn > bestIn || (nIn == bestIn && rmed < bestR)) {
                                any = true;
                                bestIn = nIn;
                                bestR = rmed;
                            }
                            pool.erase(pool.begin() + bj, pool.begin() + bi + 1);
                        }
                        if (!any || bestIn <= 0) {
                            finish("degenerate");
                            break;
                        }
                        rHyp = bestR;
                    }
                    dg.radius = rHyp;
                    gp_XYZ axis = axis0;
                    gp_XYZ center = c0 + cR * rHyp;
                    double radius = rHyp;
                    // Line first (D-fc-1.6.3) on the hypothesis inliers. A sagitta
                    // within τ lies on two tangent radii; the smaller-r tie then
                    // misses the sector, so the slab is empty (case 6). The inliers
                    // are that section. A line residual within τ is a plane.
                    {
                        std::vector<Vert> inl;
                        forOnCircle(center, axis, radius, [&](std::size_t, const CandV& v) {
                            Vert vt;
                            vt.id = v.lv;
                            vt.p = v.p;
                            inl.push_back(vt);
                        });
                        int st = 0, lv = 0;
                        double lr = 0;
                        profileCounts(inl, axis, tau, st, lv, lr);
                        if (st >= minStations(kParamsLine) && lr <= tau) {
                            dg.stations = st;
                            dg.levels = lv;
                            dg.nV = (int)inl.size();
                            dg.lineRes = lr;
                            finish("line");
                            break;
                        }
                    }
                    const int itCap = std::max(1, p.maxRefineIters);
                    bool stable = false;
                    std::vector<int> Iset;
                    auto collectI = [&](const gp_XYZ& ax, const gp_XYZ& c, double r,
                                        std::vector<int>& slab, std::vector<int>& I) {
                        slab.clear();
                        I.clear();
                        gp_XYZ u, v;
                        if (!frameOf(ax, u, v) || !(r > 0.0)) return;
                        gp_XYZ dA = A.n * (-sA);
                        gp_XYZ dB = B.n * (-sB);
                        dA -= ax * dA.Dot(ax);
                        dB -= ax * dB.Dot(ax);
                        if (!(dA.Modulus() > 1e-15) || !(dB.Modulus() > 1e-15)) return;
                        const double azA = std::atan2(dA.Dot(v), dA.Dot(u));
                        const double azB = std::atan2(dB.Dot(v), dB.Dot(u));
                        const double wAB = wrap2pi(azB - azA);
                        // D-fc-3.3. The sector is fixed by the two contact generators.
                        // Of the two arcs, the fillet is the one nearer the anchor
                        // dihedral. That dihedral is <= pi, and the shorter arc is
                        // always the nearer one. A witness on a generator used to
                        // keep the reflex (rw2b audit, case 1).
                        double start = azA;
                        double width = wAB;
                        bool startB = false;
                        if (width > kPi) {
                            start = azB;
                            width = kTwoPi - width;
                            startB = true;
                        }
                        // D-fc-4.2. Each contact generator's slack is its
                        // propagated position bound: the contact on X moves
                        // along X by the other anchor's offset bound plus the
                        // radius bound, over sin(dihedral). cross > 0 here:
                        // hem-deferred took every pair with cross <= rhs.
                        const double slackA = contactSlack(B.tolA, tau, r, cross);
                        const double slackB = contactSlack(A.tolA, tau, r, cross);
                        const double slackStart = startB ? slackB : slackA;
                        start = wrap2pi(start - slackStart);
                        width += slackA + slackB;
                        if (width > kTwoPi) width = kTwoPi;
                        // Empty I needs no slab. Membership is the τ-preimage, so it is
                        // enough to look at triangles incident to a vertex within τ.
                        if (++vertTick > 2000000000) {
                            std::fill(vertGen.begin(), vertGen.end(), 0);
                            vertTick = 1;
                        }
                        if (++triTick > 2000000000) {
                            std::fill(triGen.begin(), triGen.end(), 0);
                            triTick = 1;
                        }
                        sc.hit.clear();
                        forOnCircle(c, ax, r, [&](std::size_t i, const CandV& cv) {
                            vertGen[(std::size_t)cv.lv] = vertTick;
                            sc.hit.push_back(i);
                        });
                        bool hinted = false;
                        for (std::size_t hi : sc.hit) {
                            const CandV& cv = candPts[hi];
                            for (int t : vertTris[(std::size_t)cv.lv]) {
                                if (t < 0 || t >= nTri || triGen[(std::size_t)t] == triTick ||
                                    !candidate(t))
                                    continue;
                                triGen[(std::size_t)t] = triTick;
                                bool inSecH = true, allTauH = true;
                                for (int k = 0; k < 3; ++k) {
                                    const int lv = cornerLocal(g2l, mv, t, k);
                                    if (lv < 0 || vertGen[(std::size_t)lv] != vertTick)
                                        allTauH = false;
                                    const gp_XYZ pnt = cornerPt(mv, t, k);
                                    const gp_XYZ d = pnt - c;
                                    const double az =
                                        std::atan2((d - ax * d.Dot(ax)).Dot(v),
                                                   (d - ax * d.Dot(ax)).Dot(u));
                                    if (!(width >= kTwoPi - 1e-9) &&
                                        wrap2pi(az - start) > width + 1e-9)
                                        inSecH = false;
                                }
                                if (inSecH && allTauH && triLevelCount(t, ax) >= 2) {
                                    hinted = true;
                                    break;
                                }
                            }
                            if (hinted) break;
                        }
                        if (!hinted) return;
                        for (int t = 0; t < nTri; ++t) {
                            if (!candidate(t)) continue;
                            bool inSec = true;
                            bool allTau = true;
                            for (int k = 0; k < 3; ++k) {
                                const gp_XYZ pnt = cornerPt(mv, t, k);
                                const gp_XYZ d = pnt - c;
                                const double az = std::atan2((d - ax * d.Dot(ax)).Dot(v),
                                                             (d - ax * d.Dot(ax)).Dot(u));
                                if (!(width >= kTwoPi - 1e-9) &&
                                    wrap2pi(az - start) > width + 1e-9)
                                    inSec = false;
                                if (std::abs(radialDist(pnt, c, ax) - r) > tau) allTau = false;
                            }
                            if (!inSec) continue;
                            slab.push_back(t);
                            if (allTau && triLevelCount(t, ax) >= 2) I.push_back(t);
                        }
                        std::sort(I.begin(), I.end());
                        std::sort(slab.begin(), slab.end());
                    };

                    std::vector<int> slab;
                    collectI(axis, center, radius, slab, Iset);
                    if (Iset.empty()) {
                        finish(noRoots ? "degenerate" : "vertex-off");
                        break;
                    }
                    for (int it = 0; it < itCap; ++it) {
                        std::vector<Vert> iv;
                        uniqueVerts(mv, g2l, Iset, iv);
                        bool ok = false;
                        if (anchored) {
                            gp_XYZ ax, c;
                            double r = 0;
                            ok = fitAnchored(A.n, A.p0, B.n, B.p0, iv, itCap, axis0, sA, sB,
                                             ax, c, r);
                            if (ok) {
                                axis = ax;
                                center = c;
                                radius = r;
                            }
                        } else {
                            if ((int)iv.size() < 6) {
                                finish("underdetermined");
                                ok = false;
                                Iset.clear();
                                break;
                            }
                            gp_Dir axOut;
                            gp_Pnt cOut;
                            double rOut = 0, mOut = 0;
                            const gp_Dir seedDir(axis.X(), axis.Y(), axis.Z());
                            ok = cylinderFitLS(mv, Iset, seedDir, axOut, cOut, rOut, mOut,
                                               itCap);
                            if (ok && rOut > 0.0) {
                                axis = axOut.XYZ();
                                alignTo(axis, axis0);
                                center = cOut.XYZ();
                                radius = rOut;
                            } else {
                                ok = false;
                            }
                        }
                        if (!ok) {
                            if (dg.why == nullptr || std::strcmp(dg.why, "degenerate") == 0)
                                finish(anchored ? "vertex-off" : "degenerate");
                            Iset.clear();
                            break;
                        }
                        std::vector<int> I2, slab2;
                        collectI(axis, center, radius, slab2, I2);
                        if (I2 == Iset) {
                            slab.swap(slab2);
                            stable = true;
                            break;
                        }
                        Iset.swap(I2);
                        slab.swap(slab2);
                        if (Iset.empty()) break;
                    }
                    if (!Iset.empty() && dg.why != nullptr &&
                        std::strcmp(dg.why, "degenerate") != 0 &&
                        std::strcmp(dg.why, "underdetermined") != 0 &&
                        std::strcmp(dg.why, "vertex-off") != 0) {
                        break;
                    }
                    if (Iset.empty()) {
                        if (dg.why == nullptr) finish(stable ? "unstable" : "vertex-off");
                        break;
                    }
                    if (!stable) {
                        finish("unstable");
                        break;
                    }
                    std::vector<Vert> iv;
                    uniqueVerts(mv, g2l, Iset, iv);
                    profileCounts(iv, axis, tau, dg.stations, dg.levels, dg.lineRes);
                    dg.nV = (int)iv.size();
                    dg.radius = radius;
                    dg.axis = axis;
                    dg.nI = (int)Iset.size();
                    dg.compsI = edgeComps(Iset);
                    double iRms = 0;
                    dg.certDev = devStats(iv, center, axis, radius, iRms);
                    dg.maxDev = dg.certDev;
                    if (dg.stations < minStations(dg.params) || (!anchored && dg.levels < 2)) {
                        finish("underdetermined");
                        break;
                    }
                    if (dg.lineRes <= tau) {
                        finish("line");
                        break;
                    }
                    if (dg.certDev > tau) {
                        finish("vertex-off");
                        break;
                    }
                    double zlo = 1e300, zhi = -1e300;
                    for (const Vert& v : iv) {
                        const double z = v.p.Dot(axis);
                        zlo = std::min(zlo, z);
                        zhi = std::max(zhi, z);
                    }
                    const double L = std::max(zhi - zlo, tau);
                    auto endDist = [&](const gp_XYZ& n, const gp_XYZ& p0, double z) {
                        const gp_XYZ pt = center + axis * (z - center.Dot(axis));
                        return std::abs((pt - p0).Dot(n));
                    };
                    const double dAlo = endDist(A.n, A.p0, zlo);
                    const double dAhi = endDist(A.n, A.p0, zhi);
                    const double dBlo = endDist(B.n, B.p0, zlo);
                    const double dBhi = endDist(B.n, B.p0, zhi);
                    const double dAmid = endDist(A.n, A.p0, 0.5 * (zlo + zhi));
                    const double dBmid = endDist(B.n, B.p0, 0.5 * (zlo + zhi));
                    dg.tanDistA = std::abs(dAmid - radius);
                    dg.tanDistB = std::abs(dBmid - radius);
                    dg.tanAngA = std::max(std::abs(dAlo - radius), std::abs(dAhi - radius));
                    dg.tanAngB = std::max(std::abs(dBlo - radius), std::abs(dBhi - radius));
                    dg.dihedral = std::acos(std::max(-1.0, std::min(1.0, A.n.Dot(B.n))));
                    dg.spanTolV = (radius > 0.0) ? spanTol(A.tolA, B.tolA, radius) : 0.0;
                    if (!anchored) {
                        const bool tang =
                            dg.tanDistA <= tanDistTol(A.tolA, tau) &&
                            dg.tanDistB <= tanDistTol(B.tolA, tau) &&
                            dg.tanAngA <= tanDistTol(A.tolA, tau) &&
                            dg.tanAngB <= tanDistTol(B.tolA, tau);
                        (void)L;
                        if (!tang) {
                            finish("not-blend");
                            break;
                        }
                    }
                    double sDev = 0.0;
                    for (int t : Iset) {
                        for (int e = 0; e < 3; ++e) {
                            const gp_XYZ mid =
                                0.5 * (cornerPt(mv, t, e) + cornerPt(mv, t, (e + 1) % 3));
                            sDev = std::max(sDev, std::abs(radialDist(mid, center, axis) - radius));
                        }
                    }
                    dg.sDev = sDev;
                    std::vector<char> inI(nTri, 0), inSlab(nTri, 0);
                    for (int t : Iset)
                        if (t >= 0 && t < nTri) inI[t] = 1;
                    for (int t : slab)
                        if (t >= 0 && t < nTri) inSlab[t] = 1;
                    std::vector<int> cover;
                    for (int t : slab) {
                        if (t < 0 || t >= nTri || inI[t] || !candidate(t)) continue;
                        bool one = false, allS = true, allZ = true;
                        for (int k = 0; k < 3; ++k) {
                            const gp_XYZ pnt = cornerPt(mv, t, k);
                            const double dev = std::abs(radialDist(pnt, center, axis) - radius);
                            const double z = pnt.Dot(axis);
                            if (dev <= tau) one = true;
                            if (dev > sDev) allS = false;
                            if (z < zlo - tau || z > zhi + tau) allZ = false;
                        }
                        // D-fc-3.2 item 1: fewer than 2 axial levels is never a member
                        // (section ear). Item 7's cover does not put it back.
                        if (one && allS && allZ && triLevelCount(t, axis) >= 2)
                            cover.push_back(t);
                    }
                    std::sort(cover.begin(), cover.end());
                    dg.nCover = (int)cover.size();
                    std::vector<Vert> cv;
                    uniqueVerts(mv, g2l, cover, cv);
                    double cRms = 0;
                    dg.coverDev = cv.empty() ? 0.0 : devStats(cv, center, axis, radius, cRms);
                    std::vector<int> C = Iset;
                    C.insert(C.end(), cover.begin(), cover.end());
                    std::sort(C.begin(), C.end());
                    dg.compsC = edgeComps(C);
                    dg.nTri = (int)C.size();
                    std::vector<char> inC(nTri, 0);
                    for (int t : C)
                        if (t >= 0 && t < nTri) inC[t] = 1;
                    auto outside = [&](int prov) -> bool {
                        if (prov < 0 || prov >= (int)work.provisionals.size()) return true;
                        const auto& tris = work.provisionals[prov].tris;
                        if (tris.empty()) return false;
                        for (int t : tris)
                            if (t < 0 || t >= nTri || !inC[t]) return true;
                        return false;
                    };
                    if (!outside(A.prov) || !outside(B.prov)) {
                        finish("anchor-in-band");
                        break;
                    }
                    emitAdj(seed, axis, center, radius, inC);
                    if (sameSurface(axis, center, radius, inC)) {
                        finish("same-surface");
                        break;
                    }
                    std::vector<Vert> cvAll;
                    uniqueVerts(mv, g2l, C, cvAll);
                    const gp_XYZ xdir = lowestXDir(cvAll, center, axis);
                    const SpanInfo sp =
                        measureSpan(mv, g2l, C, cvAll, center, axis, xdir, tol.thetaBin);
                    dg.span = sp.span;
                    const double normSep = maxFacetNormalSep(mv, C);
                    const bool blendVerts = std::abs(sp.span - dg.dihedral) <= dg.spanTolV;
                    const bool blendNorms = normalsSpanDihedral(normSep, dg.spanTolV);
                    const bool blend = blendVerts && blendNorms;
                    if (!blend || priorCornerUncertified(C, inC, axis)) {
                        finish("not-blend");
                        break;
                    }
                    // Certified claim. The walk commits it, then resumes here to
                    // re-vote the remaining vertices on this side.
                    out.commit = true;
                    out.si = si;
                    out.guard = guard;
                    out.C = std::move(C);
                    out.axis = axis;
                    out.center = center;
                    out.radius = radius;
                    out.seed = seed;
                    out.dg = dg;
                    return;
                }
            }
        };

        const unsigned nWork =
            std::max(1u, detail::resolveThreadCount(detail::currentRequestedThreads()));
        std::vector<PairScratch> scratch(nWork);
        for (PairScratch& sc : scratch) {
            sc.vertGen.assign((std::size_t)nVtx, 0);
            sc.triGen.assign((std::size_t)nTri, 0);
        }
        const std::size_t nA = arefs.size();
        // This epoch's anchors, by aref index.
        std::vector<LiveAnchor> anchorTab(nA);
        std::vector<char> anchorTabOk(nA, 0);
        int anchorTabEpoch = -1;
        auto parallelEach = [&](std::size_t n, const auto& fn) {
            const unsigned w = (unsigned)std::min<std::size_t>(nWork, std::max<std::size_t>(n, 1));
            std::atomic<std::size_t> next{0};
            detail::runPool(w, [&](unsigned wi) {
                for (;;) {
                    const std::size_t i = next.fetch_add(1);
                    if (i >= n) break;
                    if (!fn(i, wi)) break;
                }
            });
        };
        auto refreshAnchorTab = [&](std::size_t from) {
            if (anchorTabEpoch == anchorEpoch) return;
            parallelEach(nA - from, [&](std::size_t i, unsigned) {
                const std::size_t k = from + i;
                anchorTabOk[k] = computeAnchor(arefs[k].prov, anchorTab[k]) ? 1 : 0;
                return true;
            });
            anchorTabEpoch = anchorEpoch;
        };
        // loadAnchor's walk-order effects for aref k, from the table.
        auto touchAnchor = [&](std::size_t k) -> bool {
            const int prov = arefs[k].prov;
            if (!anchorSlot(prov)) return false;
            if (anchorEpochSeen[(std::size_t)prov] == anchorEpoch)
                return anchorStoredOk[(std::size_t)prov] != 0;
            return storeAnchor(prov, anchorTabOk[k] != 0, anchorTab[k]);
        };
        auto flushOut = [&](const PairOut& o) {
            if (!o.diag.empty()) std::fputs(o.diag.c_str(), stderr);
            addWhy(yc, o.yc);
            nPairs += o.nPairs;
            nEvalSides += o.nEvalSides;
            nPrunedSides += o.nPrunedSides;
        };
        // Scheduling grain only: pairs per block per worker.
        const std::size_t blockCap = (std::size_t)nWork * 1024;
        std::vector<std::pair<std::size_t, std::size_t>> blockPairs;
        std::vector<PairOut> outs;
        std::size_t curA = 0, curB = 1;
        auto advance = [&](std::size_t& a, std::size_t& b) {
            if (++b >= nA) {
                ++a;
                b = a + 1;
            }
        };
        while (curB < nA) {
            if (candEpoch != anchorEpoch) rebuildCand();
            refreshAnchorTab(curA);
            blockPairs.clear();
            {
                std::size_t a = curA, b = curB;
                while (b < nA && blockPairs.size() < blockCap) {
                    blockPairs.emplace_back(a, b);
                    advance(a, b);
                }
            }
            const std::size_t nb = blockPairs.size();
            outs.clear();
            outs.resize(nb);
            // Pairs past the first commit are stale; workers stop taking them.
            std::atomic<std::size_t> stopAt{nb};
            parallelEach(nb, [&](std::size_t k, unsigned wi) {
                if (k > stopAt.load()) return false;
                const std::size_t a = blockPairs[k].first, b = blockPairs[k].second;
                if (!anchorTabOk[a] || !anchorTabOk[b]) return true;
                DiagSinkScope sink(&outs[k].diag);
                runPair(a, b, anchorTab[a], anchorTab[b], false, 0, 0, scratch[wi], outs[k]);
                if (outs[k].commit) {
                    std::size_t cur = stopAt.load();
                    while (k < cur && !stopAt.compare_exchange_weak(cur, k)) {
                    }
                }
                return true;
            });
            std::size_t k = 0;
            for (; k < nb; ++k) {
                const std::size_t a = blockPairs[k].first, b = blockPairs[k].second;
                if (!touchAnchor(a) || !touchAnchor(b)) continue;
                flushOut(outs[k]);
                if (outs[k].commit) break;
            }
            if (k == nb) {
                curA = blockPairs.back().first;
                curB = blockPairs.back().second;
                advance(curA, curB);
                continue;
            }
            const std::size_t a = blockPairs[k].first, b = blockPairs[k].second;
            const LiveAnchor A = anchorTab[a];
            const LiveAnchor B = anchorTab[b];
            PairOut cur = std::move(outs[k]);
            // commitBand always appends the band, so the side always re-votes.
            while (cur.commit) {
                commitBand(std::move(cur.C), cur.axis, cur.center, cur.radius, true, A.prov,
                           B.prov, cur.seed, cur.dg);
                rebuildCand();
                PairOut nxt;
                {
                    DiagSinkScope sink(&nxt.diag);
                    runPair(a, b, A, B, true, cur.si, cur.guard, scratch[0], nxt);
                }
                flushOut(nxt);
                cur = std::move(nxt);
            }
            curA = a;
            curB = b;
            advance(curA, curB);
        }

        // D-fc-4.9. One-sided extrusion blend (spike-fc3-onesided §1). A
        // cylinder tangent to a certified plane T along a contact generator,
        // with axis ∥ n_T × n_X for a certified X, whose far end generator lies
        // on a certified plane X' parallel to the axis that it CROSSES. The
        // station is a boundary edge of T's body parallel to that axis whose
        // other triangle is a candidate off T: it supplies the position only.
        // Free after the anchors: kParamsOneSided (contact position, r). One
        // vote per station (D-fc-3.2.3 style); the side comes from the data.
        // Certificate: τ fixed point on I, stations >= minStations, levels >= 2,
        // line first, contact on T within tanDistTol, coverage D-fc-3.2.7,
        // anchor-in-band, trim X' (holds the far vertex and >= 2 levels within
        // tanDistTol, all of C on one side, parallel, crossed), vertex span ==
        // the fitted contact→X' arc within spanTol and facet normals span it,
        // same-surface refused. A far plane that is tangent and antiparallel
        // to T is a hem (D-fc-5.11 F9): its surface is fixed by the anchors
        // and certified on its own τ-preimage.
        //
        // T walks arefs order and its stations follow its edges, one serial
        // walk that restarts T after a commit. A station's evaluation reads the
        // claim state and writes only its OsbOut, so an epoch's stations are
        // evaluated on worker threads and replayed in walk order (DIAG text,
        // tallies, the first commit), as the pair walk does. Output does not
        // depend on the worker count; threads 1 runs it inline.
        {
            const auto tOsb = std::chrono::steady_clock::now();
            int osbEdges = 0, osbStations = 0, osbHyp = 0, osbAcc = 0, osbHem = 0;
            // The pair walk refreshes its table from curA on; load every row.
            anchorTabEpoch = -1;
            std::vector<int> osbInBody((std::size_t)nTri, 0);
            int osbBodyTick = 0;
            struct OsbEdge {
                gp_XYZ pa, pb;
                int ot;
                int far;
            };
            struct OsbSt {
                std::size_t it;
                int ot;
                gp_XYZ axis, b, nTs, g;
                double sg;
            };
            struct OsbT {
                std::size_t it;
                int nEdges;
                std::size_t first, last;
            };
            struct OsbOut {
                std::string diag;
                WhyCount yc;
                int hyp = 0;
                bool commit = false;
                std::vector<int> C;
                gp_XYZ axis, center;
                double radius = 0;
                int provT = -1, provX = -1;
                BandDiag dg;
            };
            // Stations of T = arefs[it] in the current epoch, appended to out.
            // False when T is not a certified anchor.
            auto osbStationsOf = [&](std::size_t it, int& nEdges, std::vector<OsbSt>& out) -> bool {
                nEdges = 0;
                if (!touchAnchor(it) || !anchorTab[it].certified) return false;
                const LiveAnchor& TA = anchorTab[it];
                // Contact edges: boundary edges of T's body whose other
                // triangle is a candidate with a corner off T's plane. far is
                // that triangle's corner off the edge: the band side.
                std::vector<OsbEdge> ce;
                if (++osbBodyTick > 2000000000) {
                    std::fill(osbInBody.begin(), osbInBody.end(), 0);
                    osbBodyTick = 1;
                }
                for (int t : TA.body)
                    if (t >= 0 && t < nTri) osbInBody[(std::size_t)t] = osbBodyTick;
                for (int t : TA.body) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3; ++s) {
                        const int ot = otherTri(mv.triEdges[t][s], t);
                        if (ot < 0 || ot >= nTri || osbInBody[(std::size_t)ot] == osbBodyTick ||
                            !candidate(ot))
                            continue;
                        if (triOnPlane(ot, TA.n, TA.p0)) continue;
                        const int la = cornerLocal(g2l, mv, t, s);
                        const int lb = cornerLocal(g2l, mv, t, (s + 1) % 3);
                        int far = -1;
                        for (int k = 0; k < 3; ++k) {
                            const int lk = cornerLocal(g2l, mv, ot, k);
                            if (lk != la && lk != lb) far = k;
                        }
                        if (far < 0) continue;
                        ce.push_back(OsbEdge{cornerPt(mv, t, s), cornerPt(mv, t, (s + 1) % 3), ot, far});
                    }
                }
                std::stable_sort(ce.begin(), ce.end(),
                                 [](const OsbEdge& x, const OsbEdge& y) { return x.ot < y.ot; });
                nEdges = (int)ce.size();
                if (ce.empty()) return true;
                // Axis sources: certified X in arefs order, not parallel to T
                // by the hem test.
                std::vector<gp_XYZ> axes;
                for (std::size_t ix = 0; ix < nA; ++ix) {
                    if (ix == it || !touchAnchor(ix) || !anchorTab[ix].certified) continue;
                    const LiveAnchor& XA = anchorTab[ix];
                    const double cross = TA.n.Crossed(XA.n).Modulus();
                    const double rhs =
                        2.0 * (TA.tolA + tau) / TA.extent + 2.0 * (XA.tolA + tau) / XA.extent;
                    if (cross <= rhs) continue;
                    gp_XYZ a = TA.n.Crossed(XA.n);
                    a.Divide(cross);
                    canonicalHemi(a);
                    axes.push_back(a);
                }
                if (axes.empty()) return true;
                // (axis, section point, band side) already listed: collinear
                // contact edges are one station.
                const std::size_t first = out.size();
                for (const OsbEdge& e : ce) {
                    gp_XYZ dir = e.pb - e.pa;
                    const double len = dir.Modulus();
                    if (!(len > tau)) continue;
                    dir /= len;
                    std::size_t axk = axes.size();
                    for (std::size_t k = 0; k < axes.size(); ++k) {
                        if (edgeParallel(dir, len, axes[k], tau)) {
                            axk = k;
                            break;
                        }
                    }
                    if (axk == axes.size()) continue;
                    OsbSt st;
                    st.it = it;
                    st.ot = e.ot;
                    st.axis = axes[axk];
                    st.nTs = TA.n - st.axis * TA.n.Dot(st.axis);
                    const double nTm = st.nTs.Modulus();
                    if (!(nTm > 0.0)) continue;
                    st.nTs /= nTm;
                    st.g = st.axis.Crossed(st.nTs);
                    // Section point of the station on T's plane (axial coordinate 0).
                    st.b = e.pa - st.axis * e.pa.Dot(st.axis);
                    st.b -= st.nTs * (st.b - TA.p0).Dot(st.nTs);
                    // The band leaves the contact on the side of its own
                    // triangle. A triangle normal to T there is no tangent strip.
                    gp_XYZ w = cornerPt(mv, e.ot, e.far) - st.b;
                    w -= st.axis * w.Dot(st.axis);
                    const double gw = st.g.Dot(w);
                    if (std::abs(gw) <= tau) continue;
                    st.sg = gw > 0.0 ? 1.0 : -1.0;
                    bool dup = false;
                    for (std::size_t k = first; k < out.size(); ++k) {
                        const OsbSt& d = out[k];
                        if (d.axis.X() == st.axis.X() && d.axis.Y() == st.axis.Y() &&
                            d.axis.Z() == st.axis.Z() && d.sg == st.sg &&
                            (d.b - st.b).Modulus() <= tau) {
                            dup = true;
                            break;
                        }
                    }
                    if (!dup) out.push_back(st);
                }
                return true;
            };
            // One station: the vote, the fixed point and the certificate.
            // Reads the claim state only; a certified band is returned in out.
            auto evalStation = [&](const OsbSt& st, PairScratch& sc, OsbOut& out) {
                const LiveAnchor& TA = anchorTab[st.it];
                const gp_XYZ& axis = st.axis;
                const gp_XYZ& b = st.b;
                const gp_XYZ& nTs = st.nTs;
                const gp_XYZ& g = st.g;
                const double sg = st.sg;
                // Vote r per side s = sign(n_T · w): the circle tangent to T
                // at b through each candidate vertex past the contact.
                std::vector<double> pool[2];
                for (const CandV& cv : candPts) {
                    gp_XYZ w = cv.p - b;
                    w -= axis * w.Dot(axis);
                    const double h = nTs.Dot(w);
                    if (std::abs(h) <= tau) continue;
                    if (sg * g.Dot(w) < -tau) continue;
                    const double r = w.Dot(w) / (2.0 * std::abs(h));
                    if (!(r > 0.0) || !std::isfinite(r)) continue;
                    pool[h > 0.0 ? 0 : 1].push_back(r);
                }
                int bestIn = -1, bestS = 0;
                double bestR = 0.0;
                for (int si = 0; si < 2; ++si) {
                    std::vector<double>& pl = pool[si];
                    if (pl.empty()) continue;
                    std::sort(pl.begin(), pl.end());
                    const int n = (int)pl.size();
                    int j = 0, bn = 0, bj = 0;
                    for (int i = 0; i < n; ++i) {
                        while (pl[(std::size_t)i] - pl[(std::size_t)j] > 2.0 * tau) ++j;
                        if (i - j + 1 > bn) {
                            bn = i - j + 1;
                            bj = j;
                        }
                    }
                    const int mid = bj + bn / 2;
                    const double rmed =
                        (bn % 2) != 0
                            ? pl[(std::size_t)mid]
                            : 0.5 * (pl[(std::size_t)mid - 1] + pl[(std::size_t)mid]);
                    const double s = si == 0 ? 1.0 : -1.0;
                    const gp_XYZ c = b + nTs * (s * rmed);
                    int nIn = 0;
                    for (const CandV& cv : candPts) {
                        gp_XYZ w = cv.p - b;
                        w -= axis * w.Dot(axis);
                        if (sg * g.Dot(w) < -tau || s * nTs.Dot(w) < -tau) continue;
                        if (std::abs(radialDist(cv.p, c, axis) - rmed) <= tau) ++nIn;
                    }
                    if (nIn > bestIn || (nIn == bestIn && rmed < bestR)) {
                        bestIn = nIn;
                        bestR = rmed;
                        bestS = si;
                    }
                }
                // Fewer inlier vertices than minStations cannot certify (lossless).
                if (bestIn < minStations(kParamsOneSided)) return;
                out.hyp = 1;
                const double s = bestS == 0 ? 1.0 : -1.0;
                double tpos = 0.0, radius = bestR;
                auto centerOf = [&](double tp, double r) {
                    return b + g * tp + nTs * (s * r);
                };
                // Slab: past the contact and on the centre side of T, τ
                // slack each. Every I or Cover triangle has a corner in the
                // slab within τ of the circle, so only triangles incident to
                // such a vertex are examined.
                auto collect = [&](double tp, double r, std::vector<int>& slab,
                                   std::vector<int>& I) {
                    slab.clear();
                    I.clear();
                    const gp_XYZ bc = b + g * tp;
                    const gp_XYZ c = centerOf(tp, r);
                    // Past the contact: the fitted contact OR the station
                    // generator, where T's certified body ends (D-fc-5.5 F-1).
                    auto inSlab = [&](const gp_XYZ& pnt) {
                        gp_XYZ w = pnt - bc;
                        w -= axis * w.Dot(axis);
                        gp_XYZ w0 = pnt - b;
                        w0 -= axis * w0.Dot(axis);
                        return osbPastContact(g.Dot(w), g.Dot(w0), sg, tau) &&
                               !(s * nTs.Dot(w) < -tau);
                    };
                    if (++sc.triTick > 2000000000) {
                        std::fill(sc.triGen.begin(), sc.triGen.end(), 0);
                        sc.triTick = 1;
                    }
                    for (const CandV& cv : candPts) {
                        if (std::abs(radialDist(cv.p, c, axis) - r) > tau || !inSlab(cv.p))
                            continue;
                        for (int t : vertTris[(std::size_t)cv.lv]) {
                            if (t < 0 || t >= nTri || sc.triGen[(std::size_t)t] == sc.triTick ||
                                !candidate(t))
                                continue;
                            sc.triGen[(std::size_t)t] = sc.triTick;
                            bool in = true, allTau = true;
                            for (int k = 0; k < 3 && in; ++k) {
                                const gp_XYZ pnt = cornerPt(mv, t, k);
                                if (!inSlab(pnt)) in = false;
                                if (std::abs(radialDist(pnt, c, axis) - r) > tau) allTau = false;
                            }
                            if (!in) continue;
                            slab.push_back(t);
                            if (allTau && triLevelCount(t, axis) >= 2) I.push_back(t);
                        }
                    }
                    std::sort(slab.begin(), slab.end());
                    std::sort(I.begin(), I.end());
                };
                BandDiag dg;
                dg.tau = tau;
                dg.mode = "onesided";
                dg.fit = "anchored";
                dg.params = kParamsOneSided;
                dg.seed = st.ot;
                dg.axis = axis;
                dg.radius = radius;
                dg.anchorTolA = TA.tolA;
                auto finish = [&](const char* why) {
                    dg.why = why;
                    emitDiag(dg);
                    tallyWhy(out.yc, why);
                };
                std::vector<int> slab, Iset;
                collect(tpos, radius, slab, Iset);
                bool stable = false;
                const int itCap = std::max(1, p.maxRefineIters);
                for (int iter = 0; iter < itCap && !Iset.empty(); ++iter) {
                    std::vector<Vert> iv;
                    uniqueVerts(mv, g2l, Iset, iv);
                    if (!fitOneSided(iv, axis, b, g, nTs, s, itCap, tpos, radius)) {
                        Iset.clear();
                        break;
                    }
                    std::vector<int> I2, slab2;
                    collect(tpos, radius, slab2, I2);
                    if (I2 == Iset) {
                        slab.swap(slab2);
                        stable = true;
                        break;
                    }
                    Iset.swap(I2);
                    slab.swap(slab2);
                }
                dg.radius = radius;
                if (Iset.empty()) {
                    finish("vertex-off");
                    return;
                }
                if (!stable) {
                    finish("unstable");
                    return;
                }
                const gp_XYZ center = centerOf(tpos, radius);
                std::vector<Vert> iv;
                uniqueVerts(mv, g2l, Iset, iv);
                profileCounts(iv, axis, tau, dg.stations, dg.levels, dg.lineRes);
                dg.nV = (int)iv.size();
                dg.nI = (int)Iset.size();
                dg.compsI = edgeComps(Iset);
                double iRms = 0;
                dg.certDev = devStats(iv, center, axis, radius, iRms);
                dg.maxDev = dg.certDev;
                if (dg.stations < minStations(kParamsOneSided) || dg.levels < 2) {
                    finish("underdetermined");
                    return;
                }
                // D-fc-5.5 F-2: over-determination counted leave-one-provisional-out
                // on I -- without any one source (an Unclaimed provisional, or the
                // uncertified region robbed, D-fc-4.3) I keeps minStations.
                {
                    std::vector<std::array<gp_XYZ, 3>> tp;
                    std::vector<int> ts;
                    tp.reserve(Iset.size());
                    ts.reserve(Iset.size());
                    for (int t : Iset) {
                        tp.push_back({cornerPt(mv, t, 0), cornerPt(mv, t, 1), cornerPt(mv, t, 2)});
                        const int ow = uncertOwner[(std::size_t)t];
                        ts.push_back(ow >= 0 ? nProv + ow : triProv[t]);
                    }
                    if (leaveOneSourceOutStations(tp, ts, axis, tau) < minStations(kParamsOneSided)) {
                        finish("underdetermined");
                        return;
                    }
                }
                if (dg.lineRes <= tau) {
                    finish("line");
                    return;
                }
                if (dg.certDev > tau) {
                    finish("vertex-off");
                    return;
                }
                // Contact on T.
                double minOnT = std::numeric_limits<double>::infinity();
                for (const Vert& vt : iv)
                    minOnT = std::min(minOnT, std::abs((vt.p - TA.p0).Dot(TA.n)));
                dg.tanDistA = minOnT;
                if (minOnT > tanDistTol(TA.tolA, tau)) {
                    finish("no-anchor");
                    return;
                }
                // Coverage (D-fc-3.2.7).
                double zlo = std::numeric_limits<double>::infinity();
                double zhi = -std::numeric_limits<double>::infinity();
                for (const Vert& vt : iv) {
                    zlo = std::min(zlo, vt.p.Dot(axis));
                    zhi = std::max(zhi, vt.p.Dot(axis));
                }
                double sDev = 0.0;
                for (int t : Iset) {
                    for (int k = 0; k < 3; ++k) {
                        const gp_XYZ mid =
                            0.5 * (cornerPt(mv, t, k) + cornerPt(mv, t, (k + 1) % 3));
                        sDev = std::max(sDev, std::abs(radialDist(mid, center, axis) - radius));
                    }
                }
                dg.sDev = sDev;
                std::vector<char> inI((std::size_t)nTri, 0);
                for (int t : Iset) inI[(std::size_t)t] = 1;
                std::vector<int> cover;
                for (int t : slab) {
                    if (inI[(std::size_t)t] || !candidate(t)) continue;
                    bool one = false, allS = true, allZ = true;
                    for (int k = 0; k < 3; ++k) {
                        const gp_XYZ pnt = cornerPt(mv, t, k);
                        const double dv = std::abs(radialDist(pnt, center, axis) - radius);
                        const double z = pnt.Dot(axis);
                        if (dv <= tau) one = true;
                        if (dv > sDev) allS = false;
                        if (z < zlo - tau || z > zhi + tau) allZ = false;
                    }
                    // D-fc-3.2 item 1: fewer than 2 axial levels is never a member.
                    if (one && allS && allZ && triLevelCount(t, axis) >= 2) cover.push_back(t);
                }
                dg.nCover = (int)cover.size();
                {
                    std::vector<Vert> cvv;
                    uniqueVerts(mv, g2l, cover, cvv);
                    double cRms = 0;
                    dg.coverDev = cvv.empty() ? 0.0 : devStats(cvv, center, axis, radius, cRms);
                }
                std::vector<int> C = Iset;
                C.insert(C.end(), cover.begin(), cover.end());
                std::sort(C.begin(), C.end());
                dg.compsC = edgeComps(C);
                dg.nTri = (int)C.size();
                std::vector<char> inC((std::size_t)nTri, 0);
                for (int t : C) inC[(std::size_t)t] = 1;
                {
                    bool tOut = false;
                    for (int t : work.provisionals[(std::size_t)TA.prov].tris) {
                        if (t < 0 || t >= nTri || !inC[(std::size_t)t]) {
                            tOut = true;
                            break;
                        }
                    }
                    // D-fc-5.5 (F-1's premise): T's certified body ends at the
                    // station generator. A band holding any of T's certified
                    // triangles runs back over T (one surface, one face).
                    if (!tOut || bandHoldsAny(TA.body, inC)) {
                        finish("anchor-in-band");
                        return;
                    }
                }
                std::vector<Vert> cvAll;
                uniqueVerts(mv, g2l, C, cvAll);
                gp_XYZ u0, v0;
                frameOf(axis, u0, v0);
                auto azOf = [&](const gp_XYZ& pnt) {
                    gp_XYZ d = pnt - center;
                    d -= axis * d.Dot(axis);
                    return std::atan2(d.Dot(v0), d.Dot(u0));
                };
                // Far end: the C vertex at the largest arc from the contact.
                const gp_XYZ contact = center - nTs * (s * radius);
                const double azC = azOf(contact);
                gp_XYZ farP = center;
                double farArc = -1.0, azFar = azC;
                for (const Vert& vt : cvAll) {
                    const double az = azOf(vt.p);
                    const double da = std::abs(std::remainder(az - azC, kTwoPi));
                    if (da > farArc) {
                        farArc = da;
                        azFar = az;
                        farP = vt.p;
                    }
                }
                // Trim X': the first certified anchor (arefs order) that owns
                // triangles outside C, holds the far vertex and >= 2 levels of C
                // within tanDistTol, has all of C on one side, is parallel to the
                // axis, and is crossed.
                std::size_t trimIdx = nA, hemIdx = nA;
                bool farTangent = false;
                for (std::size_t ix = 0; ix < nA && !farTangent; ++ix) {
                    if (ix == st.it || !anchorTabOk[ix] || !anchorTab[ix].certified) continue;
                    const LiveAnchor& XA = anchorTab[ix];
                    {
                        bool outC = false;  // D-fc-3.2.9
                        for (int t : XA.body) {
                            if (t >= 0 && t < nTri && !inC[(std::size_t)t]) {
                                outC = true;
                                break;
                            }
                        }
                        if (!outC) continue;
                    }
                    const double tolX = tanDistTol(XA.tolA, tau);
                    if (std::abs((farP - XA.p0).Dot(XA.n)) > tolX) continue;
                    std::vector<double> zOn;
                    int nPos = 0, nNeg = 0;
                    for (const Vert& vt : cvAll) {
                        const double h = (vt.p - XA.p0).Dot(XA.n);
                        if (std::abs(h) <= tolX) zOn.push_back(vt.p.Dot(axis));
                        else if (h > 0.0) ++nPos;
                        else ++nNeg;
                    }
                    if (nPos > 0 && nNeg > 0) continue;
                    if (countLevels(zOn, tau) < 2) continue;
                    const gp_XYZ plo = center + axis * (zlo - center.Dot(axis));
                    const gp_XYZ phi = center + axis * (zhi - center.Dot(axis));
                    const double dlo = std::abs((plo - XA.p0).Dot(XA.n));
                    const double dhi = std::abs((phi - XA.p0).Dot(XA.n));
                    if (std::abs(dlo - dhi) > tolX) continue;
                    if (std::abs(0.5 * (dlo + dhi) - radius) <= tolX) {
                        // Tangent at the far end: two-sided, or a hem (O2, D-fc-5.11 F9).
                        farTangent = true;
                        hemIdx = ix;
                        break;
                    }
                    if (trimIdx < nA) continue;
                    trimIdx = ix;
                    dg.tanDistB = std::abs(0.5 * (dlo + dhi) - radius);
                    dg.anchorTolB = XA.tolA;
                }
                if (farTangent) {
                    const LiveAnchor& XH = anchorTab[hemIdx];
                    const double crossTX = TA.n.Crossed(XH.n).Modulus();
                    const double rhsTX =
                        2.0 * (TA.tolA + tau) / TA.extent + 2.0 * (XH.tolA + tau) / XH.extent;
                    if (!hemAntiparallel(crossTX, rhsTX, TA.n.Dot(XH.n))) {
                        // Tangent, not antiparallel: the pair walk's two-sided blend.
                        finish("not-blend");
                        return;
                    }
                    // D-fc-5.11 F9 (un-defers D-fc-3.1). The far plane closes a
                    // half-turn: a hem. Its surface comes from the anchors alone
                    // (kParamsHem): the axis from (T, X), r = half the T|X' gap,
                    // the centre midway, the contact at the station generator
                    // where T's certified body ends. The band is that surface's
                    // τ-preimage in the slab; nothing is fitted on it. The
                    // one-sided fit above only found the band and its far plane.
                    BandDiag hg;
                    hg.tau = tau;
                    hg.mode = "hem";
                    hg.fit = "anchored";
                    hg.params = kParamsHem;
                    hg.seed = st.ot;
                    hg.axis = axis;
                    hg.anchorTolA = TA.tolA;
                    hg.anchorTolB = XH.tolA;
                    hg.dihedral = kPi;
                    auto finishH = [&](const char* why) {
                        hg.why = why;
                        emitDiag(hg);
                        tallyWhy(out.yc, why);
                    };
                    // Gap from the contact on T, along the centre side, to X'.
                    const gp_XYZ toC = nTs * s;
                    const double den = XH.n.Dot(toC);
                    const double gap = std::abs(den) > 0.0 ? (XH.p0 - b).Dot(XH.n) / den : 0.0;
                    double rH = 0.0;
                    if (!hemRadiusOf(gap, tau, rH)) {
                        finishH("degenerate");
                        return;
                    }
                    hg.radius = rH;
                    const gp_XYZ cH = centerOf(0.0, rH);
                    std::vector<int> slabH, IH;
                    collect(0.0, rH, slabH, IH);
                    {
                        // The hem lies between its two planes: the slab is
                        // also on the centre side of X' (τ slack), as it is of
                        // T, so either plane's station gives one band.
                        const double sX = (cH - XH.p0).Dot(XH.n) < 0.0 ? -1.0 : 1.0;
                        auto insideX = [&](int t) {
                            for (int k = 0; k < 3; ++k)
                                if (sX * (cornerPt(mv, t, k) - XH.p0).Dot(XH.n) < -tau) return false;
                            return true;
                        };
                        slabH.erase(std::remove_if(slabH.begin(), slabH.end(),
                                                   [&](int t) { return !insideX(t); }),
                                    slabH.end());
                        IH.erase(std::remove_if(IH.begin(), IH.end(),
                                                [&](int t) { return !insideX(t); }),
                                 IH.end());
                    }
                    if (IH.empty()) {
                        finishH("vertex-off");
                        return;
                    }
                    std::vector<Vert> ivH;
                    uniqueVerts(mv, g2l, IH, ivH);
                    profileCounts(ivH, axis, tau, hg.stations, hg.levels, hg.lineRes);
                    hg.nV = (int)ivH.size();
                    hg.nI = (int)IH.size();
                    hg.compsI = edgeComps(IH);
                    double hRms = 0;
                    hg.certDev = devStats(ivH, cH, axis, rH, hRms);
                    hg.maxDev = hg.certDev;
                    if (hg.stations < minStations(kParamsHem) || hg.levels < 2) {
                        finishH("underdetermined");
                        return;
                    }
                    {
                        // D-fc-5.5 F-2, counted on the hem's own I.
                        std::vector<std::array<gp_XYZ, 3>> tp;
                        std::vector<int> ts;
                        tp.reserve(IH.size());
                        ts.reserve(IH.size());
                        for (int t : IH) {
                            tp.push_back({cornerPt(mv, t, 0), cornerPt(mv, t, 1), cornerPt(mv, t, 2)});
                            const int ow = uncertOwner[(std::size_t)t];
                            ts.push_back(ow >= 0 ? nProv + ow : triProv[t]);
                        }
                        if (leaveOneSourceOutStations(tp, ts, axis, tau) < minStations(kParamsHem)) {
                            finishH("underdetermined");
                            return;
                        }
                    }
                    if (hg.lineRes <= tau) {
                        finishH("line");
                        return;
                    }
                    if (hg.certDev > tau) {
                        finishH("vertex-off");
                        return;
                    }
                    // Both contacts: I meets T and holds >= 2 levels on X'.
                    const double tolXH = tanDistTol(XH.tolA, tau);
                    double minOnTH = std::numeric_limits<double>::infinity();
                    std::vector<double> zOnX;
                    double zloH = std::numeric_limits<double>::infinity();
                    double zhiH = -std::numeric_limits<double>::infinity();
                    for (const Vert& vt : ivH) {
                        minOnTH = std::min(minOnTH, std::abs((vt.p - TA.p0).Dot(TA.n)));
                        if (std::abs((vt.p - XH.p0).Dot(XH.n)) <= tolXH)
                            zOnX.push_back(vt.p.Dot(axis));
                        zloH = std::min(zloH, vt.p.Dot(axis));
                        zhiH = std::max(zhiH, vt.p.Dot(axis));
                    }
                    hg.tanDistA = minOnTH;
                    hg.tanDistB = std::abs(std::abs((cH - XH.p0).Dot(XH.n)) - rH);
                    if (minOnTH > tanDistTol(TA.tolA, tau) || countLevels(zOnX, tau) < 2) {
                        finishH("no-anchor");
                        return;
                    }
                    // Coverage (D-fc-3.2.7) of the hem's own I.
                    double sDevH = 0.0;
                    for (int t : IH) {
                        for (int k = 0; k < 3; ++k) {
                            const gp_XYZ mid =
                                0.5 * (cornerPt(mv, t, k) + cornerPt(mv, t, (k + 1) % 3));
                            sDevH = std::max(sDevH, std::abs(radialDist(mid, cH, axis) - rH));
                        }
                    }
                    hg.sDev = sDevH;
                    std::vector<char> inIH((std::size_t)nTri, 0);
                    for (int t : IH) inIH[(std::size_t)t] = 1;
                    std::vector<int> CH = IH;
                    for (int t : slabH) {
                        if (inIH[(std::size_t)t] || !candidate(t)) continue;
                        bool one = false, allS = true, allZ = true;
                        for (int k = 0; k < 3; ++k) {
                            const gp_XYZ pnt = cornerPt(mv, t, k);
                            const double dv = std::abs(radialDist(pnt, cH, axis) - rH);
                            const double z = pnt.Dot(axis);
                            if (dv <= tau) one = true;
                            if (dv > sDevH) allS = false;
                            if (z < zloH - tau || z > zhiH + tau) allZ = false;
                        }
                        if (one && allS && allZ && triLevelCount(t, axis) >= 2) {
                            CH.push_back(t);
                            ++hg.nCover;
                        }
                    }
                    std::sort(CH.begin(), CH.end());
                    if (hg.nCover > 0) {
                        std::vector<int> cov;
                        for (int t : CH)
                            if (!inIH[(std::size_t)t]) cov.push_back(t);
                        std::vector<Vert> cvv;
                        uniqueVerts(mv, g2l, cov, cvv);
                        double cRms = 0;
                        hg.coverDev = devStats(cvv, cH, axis, rH, cRms);
                    }
                    hg.compsC = edgeComps(CH);
                    hg.nTri = (int)CH.size();
                    std::vector<char> inCH((std::size_t)nTri, 0);
                    for (int t : CH) inCH[(std::size_t)t] = 1;
                    // D-fc-3.2.9 / D-fc-5.8 AF-1 on both anchors.
                    auto ownsOutside = [&](int prov) {
                        for (int t : work.provisionals[(std::size_t)prov].tris)
                            if (t < 0 || t >= nTri || !inCH[(std::size_t)t]) return true;
                        return false;
                    };
                    if (!ownsOutside(TA.prov) || bandHoldsAny(TA.body, inCH) ||
                        !ownsOutside(XH.prov) || bandHoldsAny(XH.body, inCH)) {
                        finishH("anchor-in-band");
                        return;
                    }
                    std::vector<Vert> cvH;
                    uniqueVerts(mv, g2l, CH, cvH);
                    const gp_XYZ xdirH = lowestXDir(cvH, cH, axis);
                    const SpanInfo spH =
                        measureSpan(mv, g2l, CH, cvH, cH, axis, xdirH, tol.thetaBin);
                    hg.span = spH.span;
                    hg.spanTolV = spanTol(TA.tolA, XH.tolA, rH);
                    if (!hemSpansHalfTurn(spH.span, hg.spanTolV) ||
                        !normalsSpanDihedral(maxFacetNormalSep(mv, CH), hg.spanTolV)) {
                        finishH("not-blend");
                        return;
                    }
                    emitAdj(st.ot, axis, cH, rH, inCH);
                    if (sameSurface(axis, cH, rH, inCH)) {
                        finishH("same-surface");
                        return;
                    }
                    if (diagOn())
                        diagPrintf("DIAG_REL_HEM_COMMIT stage=extrude seed=%d provT=%d minTriT=%d "
                                   "nT=%.6g,%.6g,%.6g dT=%.9g provX=%d minTriX=%d "
                                   "nX=%.6g,%.6g,%.6g dX=%.9g gap=%.9g\n",
                                   st.ot, TA.prov, TA.minTri, TA.n.X(), TA.n.Y(), TA.n.Z(),
                                   TA.p0.Dot(TA.n), XH.prov, XH.minTri, XH.n.X(), XH.n.Y(),
                                   XH.n.Z(), XH.p0.Dot(XH.n), gap);
                    out.commit = true;
                    out.C = std::move(CH);
                    out.axis = axis;
                    out.center = cH + axis * (0.5 * (zloH + zhiH) - cH.Dot(axis));
                    out.radius = rH;
                    out.provT = TA.prov;
                    out.provX = XH.prov;
                    out.dg = hg;
                    return;
                }
                if (trimIdx == nA) {
                    finish("no-anchor");
                    return;
                }
                const LiveAnchor& XT = anchorTab[trimIdx];
                // Expected arc: the contact azimuth to the circle's crossing
                // with X' nearer the far vertex.
                double expArc = -1.0;
                {
                    gp_XYZ nX = XT.n - axis * XT.n.Dot(axis);
                    const double nXm = nX.Modulus();
                    if (nXm > 0.0) {
                        nX /= nXm;
                        const double h = (XT.p0 - center).Dot(nX);
                        const double half2 = radius * radius - h * h;
                        if (half2 >= 0.0) {
                            const gp_XYZ gx = axis.Crossed(nX);
                            const gp_XYZ base = center + nX * h;
                            const double half = std::sqrt(half2);
                            const double a1 = azOf(base + gx * half);
                            const double a2 = azOf(base - gx * half);
                            const double d1 = std::abs(std::remainder(a1 - azFar, kTwoPi));
                            const double d2 = std::abs(std::remainder(a2 - azFar, kTwoPi));
                            const double aX = d1 <= d2 ? a1 : a2;
                            expArc = std::abs(std::remainder(aX - azC, kTwoPi));
                        }
                    }
                }
                dg.dihedral = expArc;
                const gp_XYZ xdir = lowestXDir(cvAll, center, axis);
                const SpanInfo sp =
                    measureSpan(mv, g2l, C, cvAll, center, axis, xdir, tol.thetaBin);
                dg.span = sp.span;
                dg.spanTolV = spanTol(TA.tolA, XT.tolA, radius);
                const double normSep = maxFacetNormalSep(mv, C);
                if (!(expArc >= 0.0) || std::abs(sp.span - expArc) > dg.spanTolV ||
                    !normalsSpanDihedral(normSep, dg.spanTolV)) {
                    finish("not-blend");
                    return;
                }
                emitAdj(st.ot, axis, center, radius, inC);
                if (sameSurface(axis, center, radius, inC)) {
                    finish("same-surface");
                    return;
                }
                const int provT = TA.prov, provX = XT.prov;
                if (diagOn()) {
                    // The anchors and every provisional the band splits
                    // (D-fc-1.6.12), so a downstream change can be traced.
                    std::string sp;
                    for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
                        int nin = 0, nout = 0;
                        for (int t : work.provisionals[(std::size_t)pi].tris) {
                            if (t >= 0 && t < nTri && inC[(std::size_t)t]) ++nin;
                            else ++nout;
                        }
                        if (nin == 0) continue;
                        const gp_XYZ pn = work.provisionals[(std::size_t)pi].plane.Direction().XYZ();
                        char buf[160];
                        std::snprintf(buf, sizeof buf, " prov%d(minTri=%d,in=%d,out=%d,n=%.6g,%.6g,%.6g)",
                                      pi, minTriOf(work.provisionals[(std::size_t)pi].tris), nin, nout,
                                      pn.X(), pn.Y(), pn.Z());
                        sp += buf;
                    }
                    diagPrintf("DIAG_REL_OSB_COMMIT stage=extrude seed=%d provT=%d minTriT=%d "
                               "nT=%.6g,%.6g,%.6g dT=%.9g provX=%d minTriX=%d "
                               "nX=%.6g,%.6g,%.6g dX=%.9g bandProvs=%s\n",
                               st.ot, provT, TA.minTri, TA.n.X(), TA.n.Y(), TA.n.Z(),
                               TA.p0.Dot(TA.n), provX, XT.minTri, XT.n.X(), XT.n.Y(),
                               XT.n.Z(), XT.p0.Dot(XT.n), sp.c_str());
                }
                out.commit = true;
                out.C = std::move(C);
                out.axis = axis;
                out.center = center + axis * (0.5 * (zlo + zhi) - center.Dot(axis));
                out.radius = radius;
                out.provT = provT;
                out.provX = provX;
                out.dg = dg;
            };
            std::vector<OsbT> osbTs;
            std::vector<OsbSt> osbSts;
            std::vector<OsbOut> osbOuts;
            std::size_t curIt = 0;
            int guard = 0;
            while (curIt < nA && guard++ < nTri + 1) {
                if (candEpoch != anchorEpoch) rebuildCand();
                refreshAnchorTab(0);
                osbTs.clear();
                osbSts.clear();
                for (std::size_t it = curIt; it < nA; ++it) {
                    OsbT tr;
                    tr.it = it;
                    tr.first = osbSts.size();
                    if (!osbStationsOf(it, tr.nEdges, osbSts)) continue;
                    tr.last = osbSts.size();
                    osbTs.push_back(tr);
                }
                const std::size_t ns = osbSts.size();
                osbOuts.clear();
                osbOuts.resize(ns);
                // Stations past the first commit are stale; workers stop taking them.
                std::atomic<std::size_t> stopAt{ns};
                parallelEach(ns, [&](std::size_t k, unsigned wi) {
                    if (k > stopAt.load()) return false;
                    DiagSinkScope sink(&osbOuts[k].diag);
                    evalStation(osbSts[k], scratch[wi], osbOuts[k]);
                    if (osbOuts[k].commit) {
                        std::size_t cur = stopAt.load();
                        while (k < cur && !stopAt.compare_exchange_weak(cur, k)) {
                        }
                    }
                    return true;
                });
                bool committed = false;
                for (const OsbT& tr : osbTs) {
                    osbEdges += tr.nEdges;
                    for (std::size_t k = tr.first; k < tr.last; ++k) {
                        OsbOut& o = osbOuts[k];
                        ++osbStations;
                        osbHyp += o.hyp;
                        if (!o.diag.empty()) std::fputs(o.diag.c_str(), stderr);
                        addWhy(yc, o.yc);
                        if (!o.commit) continue;
                        ++osbAcc;
                        if (std::strcmp(o.dg.mode, "hem") == 0) ++osbHem;
                        commitBand(std::move(o.C), o.axis, o.center, o.radius, true, o.provT,
                                   o.provX, osbSts[k].ot, o.dg);
                        committed = true;
                        // The commit changed the claim state: T restarts.
                        curIt = tr.it;
                        break;
                    }
                    if (committed) break;
                }
                if (!committed) break;
            }
            if (diagOn()) {
                const double seconds =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - tOsb).count();
                std::fprintf(stderr,
                             "DIAG_REL_OSB stage=extrude contactEdges=%d stations=%d "
                             "hypotheses=%d accept=%d seconds=%.6f hems=%d\n",
                             osbEdges, osbStations, osbHyp, osbAcc, seconds, osbHem);
            }
        }

        // Pegs. Axis is the plane anchor normal. Certified anchors keep that axis;
        // an uncertified anchor is a free cylinder fit.
        struct SeenCyl {
            gp_XYZ c;
            gp_XYZ a;
            double r;
        };
        std::vector<SeenCyl> seenCyl;
        auto onSeen = [&](const std::vector<Vert>& vs) -> bool {
            if (vs.empty()) return false;
            for (const SeenCyl& s : seenCyl) {
                bool all = true;
                for (const Vert& v : vs) {
                    if (std::abs(radialDist(v.p, s.c, s.a) - s.r) > tau) {
                        all = false;
                        break;
                    }
                }
                if (all) return true;
            }
            return false;
        };
        for (const ARef& ar : arefs) {
            LiveAnchor A;
            if (!loadAnchor(ar.prov, A)) continue;
            // A peg has one plane anchor. extrudeAnchoredPair requires two, so
            // this walk is not called. The body stays compiled (REV ORCH-3).
            if (!extrudeAnchoredPair(A.certified, false)) continue;
            gp_XYZ axis = A.n;
            canonicalHemi(axis);
            std::vector<char> pegTried(nTri, 0);
            for (int seed = 0; seed < nTri; ++seed) {
                if (pegTried[seed] || !candidate(seed) || triLevelCount(seed, axis) < 2) continue;
                std::vector<Vert> seedVs;
                uniqueVerts(mv, g2l, std::vector<int>{seed}, seedVs);
                if (onSeen(seedVs)) continue;
                std::vector<int> grown = {seed};
                std::vector<char> inG(nTri, 0);
                inG[seed] = 1;
                auto stationCount = [&]() -> int {
                    std::vector<Vert> vs;
                    uniqueVerts(mv, g2l, grown, vs);
                    int st = 0, lv = 0;
                    double lr = 0;
                    profileCounts(vs, axis, tau, st, lv, lr);
                    return st;
                };
                int growGuard = 0;
                while (stationCount() < 3 && growGuard++ < nTri) {
                    int best = -1;
                    for (int u : grown) {
                        for (int s = 0; s < 3; ++s) {
                            const int ot = otherTri(mv.triEdges[u][s], u);
                            if (!candidate(ot) || inG[ot]) continue;
                            if (best < 0 || ot < best) best = ot;
                        }
                    }
                    if (best < 0) break;
                    inG[best] = 1;
                    grown.push_back(best);
                }
                for (int gt : grown)
                    if (gt >= 0 && gt < nTri) pegTried[gt] = 1;
                if (stationCount() < 3) continue;
                gp_Pnt cOut;
                double r0 = 0;
                if (!eberlyCenterRadius(mv, grown, gp_Dir(axis), cOut, r0) || !(r0 > 0.0))
                    continue;
                gp_XYZ center = cOut.XYZ();
                std::vector<Vert> gvs;
                uniqueVerts(mv, g2l, grown, gvs);
                if (onSeen(gvs)) continue;
                seenCyl.push_back(SeenCyl{center, axis, r0});

                auto gather = [&](const gp_XYZ& ax, const gp_XYZ& c, double r) {
                    std::vector<int> I;
                    for (int t = 0; t < nTri; ++t) {
                        if (!candidate(t) || triLevelCount(t, ax) < 2) continue;
                        bool all = true;
                        for (int k = 0; k < 3; ++k) {
                            if (std::abs(radialDist(cornerPt(mv, t, k), c, ax) - r) > tau)
                                all = false;
                        }
                        if (all) I.push_back(t);
                    }
                    std::sort(I.begin(), I.end());
                    return I;
                };
                std::vector<int> Iset = gather(axis, center, r0);
                const bool pegAnchored = A.certified;
                const int itCap = std::max(1, p.maxRefineIters);
                bool stable = false;
                gp_XYZ pAxis = axis;
                for (int it = 0; it < itCap && !Iset.empty(); ++it) {
                    bool ok = false;
                    if (pegAnchored) {
                        gp_Pnt c2;
                        double r2 = 0;
                        ok = eberlyCenterRadius(mv, Iset, gp_Dir(axis), c2, r2) && r2 > 0.0;
                        if (ok) {
                            center = c2.XYZ();
                            r0 = r2;
                            pAxis = axis;
                        }
                    } else {
                        std::vector<Vert> iv;
                        uniqueVerts(mv, g2l, Iset, iv);
                        if ((int)iv.size() < 6) break;
                        gp_Dir axOut;
                        gp_Pnt c2;
                        double r2 = 0, m2 = 0;
                        ok = cylinderFitLS(mv, Iset, gp_Dir(pAxis), axOut, c2, r2, m2, itCap) &&
                             r2 > 0.0;
                        if (ok) {
                            pAxis = axOut.XYZ();
                            alignTo(pAxis, axis);
                            center = c2.XYZ();
                            r0 = r2;
                        }
                    }
                    if (!ok) break;
                    std::vector<int> I2 = gather(pAxis, center, r0);
                    if (I2 == Iset) {
                        stable = true;
                        break;
                    }
                    Iset.swap(I2);
                }
                BandDiag dg;
                dg.tau = tau;
                dg.mode = "peg";
                dg.fit = pegAnchored ? "anchored" : "free";
                dg.params = pegAnchored ? kParamsCircle : (kParamsCircle + kParamsAxisDir);
                dg.seed = seed;
                dg.axis = pAxis;
                dg.radius = r0;
                dg.anchorTolA = A.tolA;
                dg.nI = (int)Iset.size();
                dg.compsI = edgeComps(Iset);
                dg.compsC = dg.compsI;
                auto finish = [&](const char* why) {
                    dg.why = why;
                    emitDiag(dg);
                    tallyWhy(yc, why);
                };
                if (!stable || Iset.empty()) {
                    finish(Iset.empty() ? "vertex-off" : "unstable");
                    continue;
                }
                std::vector<Vert> iv;
                uniqueVerts(mv, g2l, Iset, iv);
                profileCounts(iv, pAxis, tau, dg.stations, dg.levels, dg.lineRes);
                dg.nV = (int)iv.size();
                double rms = 0;
                dg.certDev = devStats(iv, center, pAxis, r0, rms);
                dg.maxDev = dg.certDev;
                if (dg.stations < minStations(dg.params) || dg.levels < 2) {
                    finish("underdetermined");
                    continue;
                }
                if (dg.certDev > tau) {
                    finish("vertex-off");
                    continue;
                }
                std::vector<char> inC(nTri, 0);
                for (int t : Iset)
                    if (t >= 0 && t < nTri) inC[t] = 1;
                const gp_XYZ xdir = lowestXDir(iv, center, pAxis);
                const SpanInfo sp =
                    measureSpan(mv, g2l, Iset, iv, center, pAxis, xdir, tol.thetaBin);
                dg.span = sp.span;
                if (!sp.closed360) {
                    finish("not-blend");
                    continue;
                }
                const int pegMinTri = 2 * minStations(dg.params);
                if ((int)Iset.size() < pegMinTri) {
                    finish("underdetermined");
                    continue;
                }
                bool pegAdj = false;
                for (int t : Iset) {
                    if (adjacentAcceptedSurf(t)) {
                        pegAdj = true;
                        break;
                    }
                }
                if (pegAdj) {
                    finish("not-blend");
                    continue;
                }
                emitAdj(seed, pAxis, center, r0, inC);
                if (sameSurface(pAxis, center, r0, inC)) {
                    finish("same-surface");
                    continue;
                }
                seenCyl.push_back(SeenCyl{center, pAxis, r0});
                commitBand(Iset, pAxis, center, r0, false, A.prov, -1, seed, dg);
            }
        }

        // D-fc-1.6.11 observation. An accepted cylinder whose own vertices lie
        // within τ, with an adjacent unclaimed triangle also within τ, is the
        // same surface. Nothing is claimed. The stored axis is the certificate
        // when it already meets τ; a free refit is only a fallback. A second
        // fit of the unclaimed side must put every accepted vertex back within τ.
        {
            const int nAcc = (int)work.accepted.size();
            for (int i = 0; i < nAcc; ++i) {
                const Region& acc = work.accepted[i];
                if (acc.type != SurfType::Cylinder || acc.tris.empty()) continue;
                gp_XYZ ax = acc.ax.Direction().XYZ();
                gp_XYZ c = acc.ax.Location().XYZ();
                double radius = acc.radius;
                std::vector<Vert> av;
                uniqueVerts(mv, g2l, acc.tris, av);
                if (av.empty() || !(radius > 0.0)) continue;
                bool stored = true;
                for (const Vert& v : av) {
                    if (std::abs(radialDist(v.p, c, ax) - radius) > tau) stored = false;
                }
                if (!stored) {
                    gp_Dir axOut;
                    gp_Pnt cOut;
                    double rOut = 0, mOut = 0;
                    if (!cylinderFitLS(mv, acc.tris, acc.ax.Direction(), axOut, cOut, rOut, mOut,
                                       std::max(1, p.maxRefineIters)) ||
                        !(rOut > 0.0) || mOut > tau)
                        continue;
                    ax = axOut.XYZ();
                    c = cOut.XYZ();
                    radius = rOut;
                }
                std::vector<int> U;
                std::vector<char> seenU(nTri, 0);
                for (int t : acc.tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3; ++s) {
                        const int ot = otherTri(mv.triEdges[t][s], t);
                        if (!candidate(ot) || seenU[ot]) continue;
                        if (triLevelCount(ot, ax) < 2) continue;
                        bool all = true;
                        for (int k = 0; k < 3; ++k) {
                            if (std::abs(radialDist(cornerPt(mv, ot, k), c, ax) - radius) > tau)
                                all = false;
                        }
                        if (!all) continue;
                        seenU[ot] = 1;
                        U.push_back(ot);
                    }
                }
                if (U.empty()) continue;
                // Every vertex of U is already within τ of this cylinder, and
                // every accepted vertex is too. That is the same surface.
                // A free refit of a semicircle is not required: it can fail
                // while the vertices still lie on the accepted cylinder.
                BandDiag dg;
                dg.tau = tau;
                dg.mode = "fillet";
                dg.fit = "free";
                dg.seed = minTriOf(U);
                dg.axis = ax;
                dg.radius = radius;
                dg.nTri = (int)U.size();
                dg.nI = (int)U.size();
                dg.maxDev = 0.0;
                dg.why = "same-surface";
                emitDiag(dg);
                tallyWhy(yc, "same-surface");
            }
        }

        if (diagOn()) {
            const double seconds = std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - tStart)
                                       .count();
            std::fprintf(stderr,
                         "DIAG_REL_SUM stage=extrude bands=%d accept=%d line=%d "
                         "vertex-off=%d not-blend=%d same-surface=%d underdetermined=%d "
                         "degenerate=%d no-anchor=%d unstable=%d hem-deferred=%d "
                         "anchor-in-band=%d pairs=%d evaluatedSides=%d inlierAnchors=%d "
                         "seconds=%.6f prunedSides=%d\n",
                         yc.bands, yc.accept, yc.line, yc.vertexOff, yc.notBlend, yc.sameSurface,
                         yc.underdetermined, yc.degenerate, yc.noAnchor, yc.unstable,
                         yc.hemDeferred, yc.anchorInBand, nPairs, nEvalSides, nInlierAnchors,
                         seconds, nPrunedSides);
        }

        // D-fc-4.3. A dissolved region leaves the accepted list. A chamfer
        // cone names its cylinder by accepted index (pre-D), so that index
        // follows the region, or becomes -1 when its cylinder dissolved.
        if (!dissolved.empty()) {
            const int nAcc = (int)work.accepted.size();
            std::vector<char> gone((std::size_t)nAcc, 0);
            for (int ri : dissolved)
                if (ri >= 0 && ri < nAcc) gone[(std::size_t)ri] = 1;
            std::vector<int> newIdx((std::size_t)nAcc, -1);
            std::vector<Region> keep;
            keep.reserve((std::size_t)nAcc);
            for (int i = 0; i < nAcc; ++i) {
                if (gone[(std::size_t)i]) continue;
                newIdx[(std::size_t)i] = (int)keep.size();
                keep.push_back(std::move(work.accepted[(std::size_t)i]));
            }
            for (Region& r : keep) {
                if (r.origin != Origin::ChamferCone || r.filletNbrA < 0 || r.filletNbrA >= nAcc)
                    continue;
                r.filletNbrA = newIdx[(std::size_t)r.filletNbrA];
            }
            work.accepted.swap(keep);
        }
        return true;
    } catch (...) {
        return false;
    }
}

int leaveOneSourceOutStations(const std::vector<std::array<gp_XYZ, 3>>& triPts,
                              const std::vector<int>& triSrc, const gp_XYZ& axis, double tau) {
    if (triPts.size() != triSrc.size()) return 0;
    std::vector<int> srcs(triSrc);
    std::sort(srcs.begin(), srcs.end());
    srcs.erase(std::unique(srcs.begin(), srcs.end()), srcs.end());
    if (srcs.empty()) return 0;
    int fewest = std::numeric_limits<int>::max();
    std::vector<Vert> vs;
    for (int k : srcs) {
        vs.clear();
        for (std::size_t i = 0; i < triPts.size(); ++i) {
            if (triSrc[i] == k) continue;
            for (const gp_XYZ& c : triPts[i]) vs.push_back(Vert{-1, c});
        }
        int st = 0, lv = 0;
        double lr = 0;
        profileCounts(vs, axis, tau, st, lv, lr);
        fewest = std::min(fewest, st);
    }
    return fewest;
}

bool planeAnchorOf(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                   const Provisional& pr, PlaneAnchor& A) {
    A = PlaneAnchor();
    const double tau = tauOf(mv);
    const int nTri = (int)mv.nTri;
    if (!(tau > 0.0) || pr.claim != ProvClaim::Unclaimed || pr.tris.empty()) return false;
    if (!a3WouldCommit(mv, tol, pr, false)) return false;
    A.n = pr.plane.Direction().XYZ();
    const double nm = A.n.Modulus();
    if (!(nm > 1e-15)) return false;
    A.n.Divide(nm);
    A.p0 = pr.plane.Location().XYZ();
    A.body = pr.tris;
    A.dev = pr.maxVertexDev;
    A.tol = anchorTol(tau, pr.maxVertexDev);
    A.certified = anchorCertified(tau, pr.maxVertexDev);
    if (A.certified) return true;
    auto triOnPlane = [&](int t, const gp_XYZ& n, const gp_XYZ& p0) -> bool {
        if (t < 0 || t >= nTri) return false;
        for (int k = 0; k < 3; ++k)
            if (std::abs((cornerPt(mv, t, k) - p0).Dot(n)) > tau) return false;
        return true;
    };
    std::vector<int> cur = pr.tris;
    gp_XYZ n = A.n;
    gp_XYZ p0 = A.p0;
    const int seedT = maxAreaFacet(mv, pr.tris, nTri);
    if (seedT >= 0) {
        gp_XYZ sn;
        triAreaNormal(mv, seedT, sn);
        if (sn.Dot(A.n) < 0.0) sn.Reverse();
        n = sn;
        p0 = cornerPt(mv, seedT, 0);
    }
    bool peeled = false;
    const int itCap = std::max(1, p.maxRefineIters);
    for (int it = 0; it < itCap; ++it) {
        std::vector<int> in;
        for (int t : cur)
            if (triOnPlane(t, n, p0)) in.push_back(t);
        if (in.empty()) {
            peeled = false;
            break;
        }
        if (in.size() == cur.size()) break;
        peeled = true;
        gp_Ax3 np;
        if (!pcaPlane(mv, in, np)) break;
        gp_XYZ nn = np.Direction().XYZ();
        const double nnm = nn.Modulus();
        if (!(nnm > 1e-15)) break;
        nn.Divide(nnm);
        if (nn.Dot(A.n) < 0.0) nn.Reverse();
        n = nn;
        p0 = np.Location().XYZ();
        cur.swap(in);
    }
    if (!peeled) return true;
    double dev = 0.0;
    for (int t : cur)
        for (int k = 0; k < 3; ++k) dev = std::max(dev, std::abs((cornerPt(mv, t, k) - p0).Dot(n)));
    Provisional body = pr;
    body.tris = cur;
    body.area = areaOfTris(mv, cur);
    body.maxVertexDev = dev;
    if (anchorCertified(tau, dev) && a3WouldCommit(mv, tol, body, false)) {
        A.n = n;
        A.p0 = p0;
        A.body.swap(cur);
        A.dev = dev;
        A.tol = anchorTol(tau, dev);
        A.certified = true;
        A.peeled = true;
    }
    return true;
}

// D-fc-5.2. Section fillets between two certified extruded anchors.

bool sectionFilletCentre(const SectionAnchor& A, const SectionAnchor& B, double r, double hx,
                         double hy, double& px, double& py) {
    if (!(r > 0.0)) return false;
    if (!A.circle && !B.circle) return false;  // plane|plane: the extrude pair walk
    double qx[2] = {0, 0}, qy[2] = {0, 0};
    if (!A.circle || !B.circle) {
        const SectionAnchor& L = A.circle ? B : A;
        const SectionAnchor& C = A.circle ? A : B;
        const double rho = C.R + C.side * r;
        if (!(rho > 0.0)) return false;
        const double tx = -L.ny, ty = L.nx;
        const double h = L.d + L.side * r - (L.nx * C.cx + L.ny * C.cy);
        const double disc = rho * rho - h * h;
        if (disc < 0.0) return false;
        const double s = std::sqrt(disc);
        for (int k = 0; k < 2; ++k) {
            const double sg = k ? 1.0 : -1.0;
            qx[k] = C.cx + h * L.nx + sg * s * tx;
            qy[k] = C.cy + h * L.ny + sg * s * ty;
        }
    } else {
        const double r1 = A.R + A.side * r, r2 = B.R + B.side * r;
        if (!(r1 > 0.0) || !(r2 > 0.0)) return false;
        const double dx = B.cx - A.cx, dy = B.cy - A.cy;
        const double d = std::hypot(dx, dy);
        if (!(d > 0.0)) return false;
        const double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        const double h2 = r1 * r1 - a * a;
        if (h2 < 0.0) return false;
        const double h = std::sqrt(h2);
        const double ux = dx / d, uy = dy / d;
        for (int k = 0; k < 2; ++k) {
            const double sg = k ? 1.0 : -1.0;
            qx[k] = A.cx + a * ux - sg * h * uy;
            qy[k] = A.cy + a * uy + sg * h * ux;
        }
    }
    const int k = std::hypot(qx[1] - hx, qy[1] - hy) < std::hypot(qx[0] - hx, qy[0] - hy) ? 1 : 0;
    px = qx[k];
    py = qy[k];
    return true;
}

bool sectionConcentric(const SectionAnchor& A, const SectionAnchor& B) {
    if (!A.circle || !B.circle) return false;
    return std::hypot(A.cx - B.cx, A.cy - B.cy) <= A.tol + B.tol;
}

namespace {

struct SecFrame {
    gp_XYZ o, a, e1, e2;
};

SecFrame secFrame(const gp_XYZ& o, const gp_XYZ& axisIn) {
    SecFrame f;
    f.o = o;
    f.a = axisIn / axisIn.Modulus();
    gp_XYZ w = (std::abs(f.a.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    f.e1 = w - f.a * w.Dot(f.a);
    f.e1 /= f.e1.Modulus();
    f.e2 = f.a.Crossed(f.e1);
    return f;
}

struct Sec {
    double x = 0, y = 0, z = 0;
};

Sec toSec(const SecFrame& f, const gp_XYZ& p) {
    const gp_XYZ d = p - f.o;
    Sec s;
    s.x = d.Dot(f.e1);
    s.y = d.Dot(f.e2);
    s.z = d.Dot(f.a);
    return s;
}

// The fillet model: two anchors in section with their sides. Concentric
// circle anchors fix r at half their gap; the parameter is then the centre's
// azimuth on the mid circle. Otherwise the parameter is r.
struct FilletModel {
    SectionAnchor A, B;
    bool conc = false;
    double rConc = 0, rhoMid = 0, cx = 0, cy = 0;
};

bool modelAt(const FilletModel& m, double param, double hx, double hy, double& r, double& px,
             double& py) {
    if (m.conc) {
        r = m.rConc;
        px = m.cx + m.rhoMid * std::cos(param);
        py = m.cy + m.rhoMid * std::sin(param);
        return true;
    }
    r = param;
    return sectionFilletCentre(m.A, m.B, param, hx, hy, px, py);
}

double modelCost(const FilletModel& m, double param, double hx, double hy,
                 const std::vector<Sec>& q, double& r, double& px, double& py, bool& ok) {
    ok = modelAt(m, param, hx, hy, r, px, py);
    if (!ok) return 0.0;
    double e = 0;
    for (const Sec& v : q) {
        const double f = std::hypot(v.x - px, v.y - py) - r;
        e += f * f;
    }
    return e;
}

// Least squares on the one free parameter (Gauss-Newton, derivative by a
// central difference at the quantum: tau for r, tau / rho for an azimuth),
// keeping only steps that lower the cost, capped at maxIter.
bool fitModel(const FilletModel& m, const std::vector<Sec>& q, double tau, int maxIter,
              double& param, double& hx, double& hy) {
    double r = 0, px = 0, py = 0;
    bool ok = false;
    double e = modelCost(m, param, hx, hy, q, r, px, py, ok);
    if (!ok) return false;
    hx = px;
    hy = py;
    const double h = m.conc ? tau / m.rhoMid : tau;
    for (int it = 0; it < maxIter; ++it) {
        double r1, p1x, p1y, r2, p2x, p2y;
        if (!modelAt(m, param - h, hx, hy, r1, p1x, p1y) ||
            !modelAt(m, param + h, hx, hy, r2, p2x, p2y))
            break;
        double jj = 0, jf = 0;
        for (const Sec& v : q) {
            const double f = std::hypot(v.x - px, v.y - py) - r;
            const double f1 = std::hypot(v.x - p1x, v.y - p1y) - r1;
            const double f2 = std::hypot(v.x - p2x, v.y - p2y) - r2;
            const double j = (f2 - f1) / (2.0 * h);
            jj += j * j;
            jf += j * f;
        }
        if (!(jj > 0.0)) break;
        const double next = param - jf / jj;
        double rn, pnx, pny;
        bool okn = false;
        const double en = modelCost(m, next, hx, hy, q, rn, pnx, pny, okn);
        if (!okn || !(en < e)) break;
        param = next;
        e = en;
        r = rn;
        px = pnx;
        py = pny;
        hx = px;
        hy = py;
    }
    return true;
}

// The anchors' corner in section: their intersection nearest (px, py). A
// blend rounds that corner, so its arc runs between the two contacts on the
// side facing it (D-fc-1.6.7: span = the anchors' dihedral). False when the
// anchors do not meet (concentric circles: the hem, whose two arcs are both pi).
bool sectionCorner(const SectionAnchor& A, const SectionAnchor& B, double px, double py,
                   double& kx, double& ky) {
    std::vector<std::pair<double, double>> K;
    if (!A.circle && !B.circle) return false;
    if (!A.circle || !B.circle) {
        const SectionAnchor& L = A.circle ? B : A;
        const SectionAnchor& C = A.circle ? A : B;
        const double h = L.d - (L.nx * C.cx + L.ny * C.cy);
        const double disc = C.R * C.R - h * h;
        if (disc < 0.0) return false;
        const double s = std::sqrt(disc);
        for (int k = 0; k < 2; ++k) {
            const double sg = k ? 1.0 : -1.0;
            K.emplace_back(C.cx + h * L.nx - sg * s * L.ny, C.cy + h * L.ny + sg * s * L.nx);
        }
    } else {
        const double dx = B.cx - A.cx, dy = B.cy - A.cy;
        const double d = std::hypot(dx, dy);
        if (!(d > 0.0)) return false;
        const double a = (A.R * A.R - B.R * B.R + d * d) / (2.0 * d);
        const double h2 = A.R * A.R - a * a;
        if (h2 < 0.0) return false;
        const double h = std::sqrt(h2), ux = dx / d, uy = dy / d;
        for (int k = 0; k < 2; ++k) {
            const double sg = k ? 1.0 : -1.0;
            K.emplace_back(A.cx + a * ux - sg * h * uy, A.cy + a * uy + sg * h * ux);
        }
    }
    const int k = std::hypot(K[1].first - px, K[1].second - py) <
                          std::hypot(K[0].first - px, K[0].second - py)
                      ? 1
                      : 0;
    kx = K[(std::size_t)k].first;
    ky = K[(std::size_t)k].second;
    return true;
}

double wrapTwoPi(double a) {
    while (a < 0.0) a += kTwoPi;
    while (a >= kTwoPi) a -= kTwoPi;
    return a;
}

struct LineAnc {
    gp_XYZ n;
    double d = 0, tol = 0;
    int prov = -1;
};

bool diagSectionOn() { return diagOn(); }

}  // namespace

bool claimSectionFillets(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                         SegmentWork& work, const std::vector<SectionCircle>& circles) {
    try {
        const double tau = tauOf(mv);
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;
        const int nTri = (int)mv.nTri;
        const int nVtx = (int)mv.nVtx;
        const int maxIter = p.maxRefineIters > 0 ? p.maxRefineIters : 1;
        std::vector<int> g2l;
        {
            int maxG = 0;
            for (int t = 0; t < nTri; ++t) {
                const int gt = mv.compTris ? mv.compTris[t] : t;
                for (int k = 0; k < 3; ++k) maxG = std::max(maxG, mv.tris[gt][k]);
            }
            for (int lv = 0; lv < nVtx; ++lv)
                maxG = std::max(maxG, mv.compVtx ? mv.compVtx[lv] : lv);
            g2l.assign((std::size_t)maxG + 1, -1);
            for (int lv = 0; lv < nVtx; ++lv) {
                const int g = mv.compVtx ? mv.compVtx[lv] : lv;
                if (g >= 0) g2l[(std::size_t)g] = lv;
            }
        }
        std::vector<std::vector<int>> nbr((std::size_t)nTri);
        {
            std::vector<std::array<int, 2>> et((std::size_t)std::max((int)mv.nEdge, 0),
                                               std::array<int, 2>{{-1, -1}});
            for (int t = 0; t < nTri; ++t)
                for (int s = 0; s < 3; ++s) {
                    const int e = mv.triEdges[t][s];
                    if (e < 0 || e >= (int)mv.nEdge) continue;
                    if (et[(std::size_t)e][0] < 0) et[(std::size_t)e][0] = t;
                    else if (et[(std::size_t)e][0] != t) et[(std::size_t)e][1] = t;
                }
            for (const auto& e : et)
                if (e[0] >= 0 && e[1] >= 0) {
                    nbr[(std::size_t)e[0]].push_back(e[1]);
                    nbr[(std::size_t)e[1]].push_back(e[0]);
                }
        }
        // Candidates: Unclaimed provisional triangles held by no region.
        std::vector<char> freeTri((std::size_t)nTri, 0);
        for (const Provisional& pr : work.provisionals) {
            if (pr.claim != ProvClaim::Unclaimed) continue;
            for (int t : pr.tris)
                if (t >= 0 && t < nTri) freeTri[(std::size_t)t] = 1;
        }
        std::vector<int> atomOf((std::size_t)nTri, -1);
        for (int i = 0; i < (int)work.accepted.size(); ++i) {
            const Region& r = work.accepted[(std::size_t)i];
            for (int t : r.tris)
                if (t >= 0 && t < nTri) {
                    freeTri[(std::size_t)t] = 0;
                    if (torusAtom(r)) atomOf[(std::size_t)t] = i;
                }
        }
        // Line anchors (D-fc-1.4): certified planes A3 would commit.
        std::vector<LineAnc> lines;
        std::vector<std::vector<int>> lineBody;
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
            PlaneAnchor pa;
            if (!planeAnchorOf(mv, p, tol, work.provisionals[(std::size_t)pi], pa) || !pa.certified)
                continue;
            LineAnc L;
            L.n = pa.n;
            L.d = pa.n.Dot(pa.p0);
            L.tol = pa.tol;
            L.prov = pi;
            lines.push_back(L);
            lineBody.push_back(pa.body);
        }
        auto vertsOf = [&](const std::vector<int>& tris) {
            std::vector<int> g;
            for (int t : tris) {
                if (t < 0 || t >= nTri) continue;
                const int gt = mv.compTris ? mv.compTris[t] : t;
                for (int k = 0; k < 3; ++k) g.push_back(mv.tris[gt][k]);
            }
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            return g;
        };
        std::vector<char> dead(work.accepted.size(), 0);
        int nAccept = 0, nEval = 0;
        const int nAcc0 = (int)work.accepted.size();
        for (int ai = 0; ai < nAcc0; ++ai) {
            if (dead[(std::size_t)ai]) continue;
            const Region& atom = work.accepted[(std::size_t)ai];
            if (!torusAtom(atom)) continue;
            std::vector<int> atris = atom.tris;
            std::sort(atris.begin(), atris.end());
            gp_Dir axOut = atom.ax.Direction();
            gp_Pnt cen;
            double rAtom = 0, devAtom = 0;
            if (!cylinderFitLS(mv, atris, atom.ax.Direction(), axOut, cen, rAtom, devAtom, maxIter))
                continue;
            if (!(rAtom > tau)) continue;
            const gp_XYZ aAtom = axOut.XYZ();
            const double tolR = anchorTol(tau, devAtom);
            const std::vector<int> agids = vertsOf(atris);
            double zLo = 1e300, zHi = -1e300;
            for (int g : agids) {
                const double z = (mv.pts[g] - cen.XYZ()).Dot(aAtom);
                zLo = std::min(zLo, z);
                zHi = std::max(zHi, z);
            }
            const double Lat = std::max(zHi - zLo, tau);
            // Anchors tangent to the atom's circle, in a fixed order: circles,
            // then lines. Each records its side.
            struct Tan {
                bool circle = false;
                int idx = -1;
                int side = 1;
            };
            std::vector<Tan> tans;
            for (int ci = 0; ci < (int)circles.size(); ++ci) {
                const SectionCircle& C = circles[(std::size_t)ci];
                if (C.acc == ai) continue;
                // Hypothesis: tangency within both fits' own tolerances.
                if (C.axis.Crossed(aAtom).Modulus() > (C.tol + tolR) / Lat) continue;
                const SecFrame f = secFrame(C.origin, C.axis);
                const Sec sc = toSec(f, cen.XYZ());
                const double dc = std::hypot(sc.x, sc.y);
                if (dc <= C.tol + tolR) continue;  // coaxial with the atom
                Tan t;
                t.circle = true;
                t.idx = ci;
                if (std::fabs(dc - (C.radius + rAtom)) <= C.tol + tolR) t.side = 1;
                else if (std::fabs(dc - (C.radius - rAtom)) <= C.tol + tolR) t.side = -1;
                else continue;
                tans.push_back(t);
            }
            if (tans.empty()) continue;  // plane|plane belongs to the extrude pair walk
            for (int li = 0; li < (int)lines.size(); ++li) {
                const LineAnc& L = lines[(std::size_t)li];
                if (std::fabs(L.n.Dot(aAtom)) > (L.tol + tolR) / Lat) continue;
                const double dl = L.n.Dot(cen.XYZ()) - L.d;
                if (std::fabs(std::fabs(dl) - rAtom) > L.tol + tolR) continue;
                Tan t;
                t.idx = li;
                t.side = dl >= 0.0 ? 1 : -1;
                tans.push_back(t);
            }
            bool done = false;
            for (std::size_t ia = 0; ia < tans.size() && !done; ++ia) {
                if (!tans[ia].circle) continue;
                for (std::size_t ib = ia + 1; ib < tans.size() && !done; ++ib) {
                    ++nEval;
                    const SectionCircle& C0 = circles[(std::size_t)tans[ia].idx];
                    const SecFrame f = secFrame(C0.origin, C0.axis);
                    auto toAnchor = [&](const Tan& t) {
                        SectionAnchor s;
                        s.side = t.side;
                        if (t.circle) {
                            const SectionCircle& C = circles[(std::size_t)t.idx];
                            const Sec sc = toSec(f, C.origin);
                            s.circle = true;
                            s.cx = sc.x;
                            s.cy = sc.y;
                            s.R = C.radius;
                            s.tol = C.tol;
                        } else {
                            const LineAnc& L = lines[(std::size_t)t.idx];
                            double nx = L.n.Dot(f.e1), ny = L.n.Dot(f.e2);
                            const double nm = std::hypot(nx, ny);
                            nx /= nm;
                            ny /= nm;
                            s.nx = nx;
                            s.ny = ny;
                            s.d = (L.d - L.n.Dot(f.o)) / nm;
                            s.tol = L.tol;
                        }
                        return s;
                    };
                    // Both anchors parallel to the certified extrusion axis over
                    // the atom's length (tanAngTol, D-fc-1.6.7).
                    const Tan& tb = tans[ib];
                    if (tb.circle) {
                        const SectionCircle& C1 = circles[(std::size_t)tb.idx];
                        if (C1.axis.Crossed(C0.axis).Modulus() > tanAngTol(C1.tol, tau, Lat)) continue;
                    } else if (std::fabs(lines[(std::size_t)tb.idx].n.Dot(C0.axis)) >
                               tanAngTol(lines[(std::size_t)tb.idx].tol, tau, Lat)) {
                        continue;
                    }
                    FilletModel m;
                    m.A = toAnchor(tans[ia]);
                    m.B = toAnchor(tans[ib]);
                    const Sec hc = toSec(f, cen.XYZ());
                    double param = rAtom;
                    if (sectionConcentric(m.A, m.B)) {
                        // The fillet lies outside the inner circle and inside the
                        // outer one; any other pair of sides has no tangent circle.
                        const SectionAnchor& in = (m.A.R < m.B.R) ? m.A : m.B;
                        const SectionAnchor& out = (m.A.R < m.B.R) ? m.B : m.A;
                        if (in.side != 1 || out.side != -1) continue;
                        m.conc = true;
                        m.cx = in.cx;
                        m.cy = in.cy;
                        m.rConc = 0.5 * (out.R - in.R);
                        m.rhoMid = in.R + m.rConc;
                        if (!(m.rConc > tau)) continue;
                        param = std::atan2(hc.y - m.cy, hc.x - m.cx);
                    }
                    const char* why = "accept";
                    double hx = hc.x, hy = hc.y;
                    std::vector<int> band, joined;
                    double r = 0, px = 0, py = 0, phiS = 0, span = 0;
                    double maxDev = 0, rms = 0, lineRes = 0, vSpan = 0, spanTolV = 0;
                    int stations = 0, levels = 0;
                    bool stable = false;
                    {
                        std::vector<Sec> q;
                        for (int g : agids) q.push_back(toSec(f, mv.pts[g]));
                        if (!fitModel(m, q, tau, maxIter, param, hx, hy)) why = "degenerate";
                    }
                    auto contactAngle = [&](const SectionAnchor& s, double cpx, double cpy) {
                        double qx, qy;
                        if (s.circle) {
                            const double dx = cpx - s.cx, dy = cpy - s.cy;
                            const double dm = std::hypot(dx, dy);
                            qx = s.cx + dx / dm * s.R;
                            qy = s.cy + dy / dm * s.R;
                        } else {
                            qx = cpx - s.side * r * s.nx;
                            qy = cpy - s.side * r * s.ny;
                        }
                        return std::atan2(qy - cpy, qx - cpx);
                    };
                    for (int it = 0; it < maxIter && std::strcmp(why, "accept") == 0; ++it) {
                        if (!modelAt(m, param, hx, hy, r, px, py)) {
                            why = "degenerate";
                            break;
                        }
                        // The fillet arc runs between its two contact points, on
                        // the side of the atom's vertices.
                        const double phA = contactAngle(m.A, px, py);
                        const double phB = contactAngle(m.B, px, py);
                        // The blend arc faces the anchors' corner. Concentric
                        // anchors (a hem) have no corner and two arcs of pi:
                        // the side is the atom's (from the data).
                        double psiM = 0;
                        double kx = 0, ky = 0;
                        if (sectionCorner(m.A, m.B, px, py, kx, ky)) {
                            psiM = std::atan2(ky - py, kx - px);
                        } else if (m.conc) {
                            double sx = 0, sy = 0;
                            for (int g : agids) {
                                const Sec s = toSec(f, mv.pts[g]);
                                const double dm = std::hypot(s.x - px, s.y - py);
                                if (dm > 0.0) {
                                    sx += (s.x - px) / dm;
                                    sy += (s.y - py) / dm;
                                }
                            }
                            psiM = std::atan2(sy, sx);
                        } else {
                            why = "not-blend";
                            break;
                        }
                        const double sAB = wrapTwoPi(phB - phA);
                        const bool fromA = wrapTwoPi(psiM - phA) <= sAB;
                        phiS = fromA ? phA : phB;
                        span = fromA ? sAB : kTwoPi - sAB;
                        const double slackS = sectionContactSlack(fromA ? m.A.tol : m.B.tol, tau, r);
                        const double slackE = sectionContactSlack(fromA ? m.B.tol : m.A.tol, tau, r);
                        auto inSlab = [&](const gp_XYZ& pt, double& dev, double& z) {
                            const Sec s = toSec(f, pt);
                            z = s.z;
                            dev = std::fabs(std::hypot(s.x - px, s.y - py) - r);
                            double u = wrapTwoPi(std::atan2(s.y - py, s.x - px) - phiS + slackS);
                            return u <= span + slackS + slackE;
                        };
                        auto member = [&](int t) {
                            double zl = 1e300, zh = -1e300;
                            for (int k = 0; k < 3; ++k) {
                                double dev, z;
                                if (!inSlab(cornerPt(mv, t, k), dev, z) || !(dev <= tau)) return false;
                                zl = std::min(zl, z);
                                zh = std::max(zh, z);
                            }
                            return zh - zl > tau;
                        };
                        std::vector<int> I;
                        for (int t = 0; t < nTri; ++t)
                            if (freeTri[(std::size_t)t] && member(t)) I.push_back(t);
                        // Atoms join only whole, every vertex in the slab within tau.
                        std::vector<int> J;
                        std::vector<char> tested(work.accepted.size(), 0);
                        for (int t = 0; t < nTri; ++t) {
                            const int aj = atomOf[(std::size_t)t];
                            if (aj < 0 || dead[(std::size_t)aj] || tested[(std::size_t)aj]) continue;
                            tested[(std::size_t)aj] = 1;
                            double dev, z;
                            if (!inSlab(cornerPt(mv, t, 0), dev, z) || !(dev <= tau)) continue;
                            bool all = true;
                            for (int g : vertsOf(work.accepted[(std::size_t)aj].tris)) {
                                if (!inSlab(mv.pts[g], dev, z) || !(dev <= tau)) {
                                    all = false;
                                    break;
                                }
                            }
                            if (all) J.push_back(aj);
                        }
                        std::sort(J.begin(), J.end());
                        if (std::find(J.begin(), J.end(), ai) == J.end()) {
                            why = "vertex-off";
                            break;
                        }
                        std::vector<int> B = I;
                        for (int aj : J)
                            for (int t : work.accepted[(std::size_t)aj].tris) B.push_back(t);
                        std::sort(B.begin(), B.end());
                        B.erase(std::unique(B.begin(), B.end()), B.end());
                        if (B == band && J == joined) {
                            stable = true;
                            break;
                        }
                        band = B;
                        joined = J;
                        std::vector<Sec> q;
                        for (int g : vertsOf(band)) q.push_back(toSec(f, mv.pts[g]));
                        if (!fitModel(m, q, tau, maxIter, param, hx, hy)) why = "degenerate";
                    }
                    if (std::strcmp(why, "accept") == 0 && !stable) why = "unstable";
                    // A band of joined atoms alone claims nothing new: merging
                    // same-surface regions is reclaim's (D-fc-1.6.11, D-fc-4.4).
                    if (std::strcmp(why, "accept") == 0) {
                        bool anyFree = false;
                        for (int t : band)
                            if (freeTri[(std::size_t)t]) anyFree = true;
                        if (!anyFree) why = "same-surface";
                    }
                    std::vector<int> bgids;
                    std::vector<double> psis, rhos, zs;
                    if (std::strcmp(why, "accept") == 0) {
                        bgids = vertsOf(band);
                        double s2 = 0;
                        for (int g : bgids) {
                            const Sec s = toSec(f, mv.pts[g]);
                            const double dm = std::hypot(s.x - px, s.y - py);
                            const double d = std::fabs(dm - r);
                            maxDev = std::max(maxDev, d);
                            s2 += d * d;
                            psis.push_back(wrapTwoPi(std::atan2(s.y - py, s.x - px)));
                            rhos.push_back(dm);
                            zs.push_back(s.z);
                        }
                        rms = bgids.empty() ? 0.0 : std::sqrt(s2 / (double)bgids.size());
                        stations = countColumns(psis, rhos, tau);
                        levels = countLevels(zs, tau);
                        // D-fc-4.6 F2: over-determination counted leave-one-atom-out.
                        bool det = stations >= minStations(kParamsAnchored);
                        for (int aj : joined) {
                            if (!det) break;
                            std::vector<char> inA((std::size_t)nTri, 0);
                            for (int t : work.accepted[(std::size_t)aj].tris)
                                if (t >= 0 && t < nTri) inA[(std::size_t)t] = 1;
                            std::vector<int> rest;
                            for (int t : band)
                                if (!inA[(std::size_t)t]) rest.push_back(t);
                            std::vector<double> ps2, rs2;
                            for (int g : vertsOf(rest)) {
                                const Sec s = toSec(f, mv.pts[g]);
                                ps2.push_back(wrapTwoPi(std::atan2(s.y - py, s.x - px)));
                                rs2.push_back(std::hypot(s.x - px, s.y - py));
                            }
                            if (countColumns(ps2, rs2, tau) < minStations(kParamsAnchored)) det = false;
                        }
                        std::vector<gp_XYZ> bp;
                        for (int g : bgids) bp.push_back(mv.pts[g]);
                        lineRes = sectionLineResid(bp, f.a);
                        if (!det) why = "underdetermined";
                        else if (lineRes <= tau) why = "line";
                        else if (!(maxDev <= tau)) why = "vertex-off";
                    }
                    if (std::strcmp(why, "accept") == 0) {
                        // Contact: the band shares a vertex with each anchor on it.
                        for (int side = 0; side < 2 && std::strcmp(why, "accept") == 0; ++side) {
                            const Tan& t = side ? tans[ib] : tans[ia];
                            const SectionAnchor& s = side ? m.B : m.A;
                            const std::vector<int> atrisX =
                                t.circle ? work.accepted[(std::size_t)circles[(std::size_t)t.idx].acc].tris
                                         : lineBody[(std::size_t)t.idx];
                            const std::vector<int> xg = vertsOf(atrisX);
                            bool touch = false;
                            for (int g : bgids) {
                                if (!std::binary_search(xg.begin(), xg.end(), g)) continue;
                                const Sec sp = toSec(f, mv.pts[g]);
                                const double dx = s.circle ? std::fabs(std::hypot(sp.x - s.cx, sp.y - s.cy) - s.R)
                                                           : std::fabs(s.nx * sp.x + s.ny * sp.y - s.d);
                                if (dx <= tanDistTol(s.tol, tau)) {
                                    touch = true;
                                    break;
                                }
                            }
                            if (!touch) why = "no-anchor";
                        }
                    }
                    if (std::strcmp(why, "accept") == 0) {
                        // Normals span the tangent-to-tangent arc (D-fc-1.6.7).
                        std::vector<double> ps = psis;
                        std::sort(ps.begin(), ps.end());
                        double gap = ps.empty() ? kTwoPi : (ps.front() + kTwoPi) - ps.back();
                        for (std::size_t i = 0; i + 1 < ps.size(); ++i) gap = std::max(gap, ps[i + 1] - ps[i]);
                        vSpan = kTwoPi - gap;
                        spanTolV = spanTol(m.A.tol, m.B.tol, r);
                        if (std::fabs(vSpan - span) > spanTolV ||
                            !normalsSpanDihedral(maxFacetNormalSep(mv, band), spanTolV))
                            why = "not-blend";
                    }
                    if (std::strcmp(why, "accept") == 0) {
                        // D-fc-1.6.11 on the claimed set.
                        std::vector<char> inB((std::size_t)nTri, 0);
                        for (int t : band) inB[(std::size_t)t] = 1;
                        for (int ri = 0; ri < (int)work.accepted.size() && std::strcmp(why, "accept") == 0; ++ri) {
                            if (dead[(std::size_t)ri] || std::binary_search(joined.begin(), joined.end(), ri)) continue;
                            const Region& rr = work.accepted[(std::size_t)ri];
                            bool touch = false;
                            for (int t : rr.tris) {
                                if (t < 0 || t >= nTri) continue;
                                for (int n : nbr[(std::size_t)t])
                                    if (inB[(std::size_t)n]) touch = true;
                                if (touch) break;
                            }
                            if (!touch || rr.tris.empty()) continue;
                            bool all = true;
                            for (int g : vertsOf(rr.tris)) {
                                const Sec s = toSec(f, mv.pts[g]);
                                if (!(std::fabs(std::hypot(s.x - px, s.y - py) - r) <= tau)) {
                                    all = false;
                                    break;
                                }
                            }
                            if (all) why = "same-surface";
                        }
                    }
                    const bool accept = std::strcmp(why, "accept") == 0;
                    int straddle = 0;
                    if (accept) {
                        const gp_XYZ P = f.o + f.e1 * px + f.e2 * py;
                        std::vector<char> inB((std::size_t)nTri, 0);
                        for (int t : band) inB[(std::size_t)t] = 1;
                        gp_XYZ cent(0, 0, 0);
                        double area = 0, sigma = 0;
                        for (int t : band) {
                            gp_XYZ nt;
                            const double at = triAreaNormal(mv, t, nt);
                            const gp_XYZ c3 = (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) / 3.0;
                            area += at;
                            cent += c3 * at;
                        }
                        if (area > 0.0) cent /= area;
                        const gp_XYZ loc = P + f.a * ((cent - P).Dot(f.a));
                        for (int t : band) {
                            gp_XYZ nt;
                            const double at = triAreaNormal(mv, t, nt);
                            gp_XYZ rad = (cornerPt(mv, t, 0) + cornerPt(mv, t, 1) + cornerPt(mv, t, 2)) / 3.0 - loc;
                            rad -= f.a * rad.Dot(f.a);
                            const double rm = rad.Modulus();
                            if (rm > 0.0) sigma += at * nt.Dot(rad / rm);
                        }
                        std::vector<Vert> verts;
                        uniqueVerts(mv, g2l, band, verts);
                        const gp_XYZ xdir = lowestXDir(verts, loc, f.a);
                        const SpanInfo sp = measureSpan(mv, g2l, band, verts, loc, f.a, xdir, tol.thetaBin);
                        const bool outward = sigma > 0.0;
                        const double useSpan = std::max(sp.span, 1e-12);
                        int nSides = (int)std::llround(kTwoPi * (double)std::max(sp.nBands, 1) / useSpan);
                        if (nSides < 1) nSides = 1;
                        for (const Provisional& pr : work.provisionals) {
                            if (pr.claim != ProvClaim::Unclaimed) continue;
                            int nin = 0, nout = 0;
                            for (int t : pr.tris) (t >= 0 && t < nTri && inB[(std::size_t)t]) ? ++nin : ++nout;
                            if (nin > 0 && nout > 0) ++straddle;
                        }
                        const int host = joined.front();
                        Region& R = work.accepted[(std::size_t)host];
                        R.type = SurfType::Cylinder;
                        R.origin = Origin::FilletStrip;
                        R.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(f.a), gp_Dir(xdir));
                        R.radius = r;
                        R.uMin = sp.uMin;
                        R.uMax = sp.uMax;
                        R.vMin = sp.vMin;
                        R.vMax = sp.vMax;
                        R.closed360 = false;
                        R.outwardNormal = outward;
                        R.tris = band;
                        R.maxVertexDev = maxDev;
                        R.rmsVertexDev = rms;
                        R.nSides = nSides;
                        R.chordSagitta = chordSagitta(r, nSides);
                        R.dVolPredicted = dVolCylinderSector(area, r, nSides, outward);
                        R.maxVertexSnap = 0.0;
                        R.reject = Reject::None;
                        R.builtAs = BuiltAs::NotBuilt;
                        R.filletNbrA = tans[ia].circle ? -1 : lines[(std::size_t)tans[ia].idx].prov;
                        R.filletNbrB = tans[ib].circle ? -1 : lines[(std::size_t)tans[ib].idx].prov;
                        for (std::size_t k = 1; k < joined.size(); ++k) {
                            dead[(std::size_t)joined[k]] = 1;
                            work.accepted[(std::size_t)joined[k]].tris.clear();
                        }
                        for (int t : band) {
                            freeTri[(std::size_t)t] = 0;
                            atomOf[(std::size_t)t] = -1;
                        }
                        // The unclaimed triangles leave their provisionals; a
                        // remainder is refit and split into its components.
                        const int nProvNow = (int)work.provisionals.size();
                        for (int pi = 0; pi < nProvNow; ++pi) {
                            Provisional& pr = work.provisionals[(std::size_t)pi];
                            if (pr.claim != ProvClaim::Unclaimed) continue;
                            std::vector<int> keep;
                            bool hit = false;
                            for (int t : pr.tris) {
                                if (t >= 0 && t < nTri && inB[(std::size_t)t]) hit = true;
                                else keep.push_back(t);
                            }
                            if (!hit) continue;
                            if (keep.empty()) {
                                pr.claim = ProvClaim::ConsumedFillet;
                                continue;
                            }
                            std::vector<char> inK((std::size_t)nTri, 0);
                            for (int t : keep) inK[(std::size_t)t] = 1;
                            std::vector<char> seen((std::size_t)nTri, 0);
                            std::vector<std::vector<int>> comps;
                            for (int t0 : keep) {
                                if (seen[(std::size_t)t0]) continue;
                                std::vector<int> comp, st{t0};
                                seen[(std::size_t)t0] = 1;
                                while (!st.empty()) {
                                    const int c = st.back();
                                    st.pop_back();
                                    comp.push_back(c);
                                    for (int n : nbr[(std::size_t)c])
                                        if (inK[(std::size_t)n] && !seen[(std::size_t)n]) {
                                            seen[(std::size_t)n] = 1;
                                            st.push_back(n);
                                        }
                                }
                                std::sort(comp.begin(), comp.end());
                                comps.push_back(std::move(comp));
                            }
                            std::sort(comps.begin(), comps.end());
                            pr.tris = comps[0];
                            refitPlaneProv(mv, g2l, pr);
                            for (std::size_t c = 1; c < comps.size(); ++c) {
                                Provisional np = pr;
                                np.tris = comps[c];
                                refitPlaneProv(mv, g2l, np);
                                work.provisionals.push_back(std::move(np));
                            }
                        }
                        ++nAccept;
                        done = true;
                    }
                    if (diagSectionOn()) {
                        diagPrintf("DIAG_REL_BAND stage=extrude mode=section seed=%d nTri=%zu nV=%zu "
                                   "stations=%d levels=%d params=%d fit=anchored radius=%.9g "
                                   "axis=%.6g,%.6g,%.6g maxDevMM=%.6g maxDevTau=%.6g lineResTau=%.6g "
                                   "span=%.9g dihedral=%.9g spanTol=%.6g straddle=%d atom=%d "
                                   "atoms=%zu kindA=%s kindB=%s sideA=%d sideB=%d hem=%d why=%s\n",
                                   minTriOf(band.empty() ? atris : band), band.size(), bgids.size(),
                                   stations, levels, kParamsAnchored, r, f.a.X(), f.a.Y(), f.a.Z(),
                                   maxDev, maxDev / tau, lineRes / tau, vSpan, span, spanTolV,
                                   straddle, ai, joined.size(), tans[ia].circle ? "circle" : "line",
                                   tans[ib].circle ? "circle" : "line", tans[ia].side, tans[ib].side,
                                   m.conc ? 1 : 0, why);
                    }
                    // audit-fx-arch F-3: a GEOMETRIC refusal of the atom's band ends
                    // the atom; another anchor pair is not tried on it.
                    if (lineChainGeometric(why) || std::strcmp(why, "not-blend") == 0) done = true;
                }
            }
        }
        if (diagSectionOn())
            diagPrintf("DIAG_REL_SUM stage=section pairs=%d accept=%d circles=%zu lines=%zu\n", nEval,
                       nAccept, circles.size(), lines.size());
        // A joined atom other than the host leaves the accepted list (D-fc-4.6
        // F2). A chamfer cone names its cylinder by accepted index.
        bool anyDead = false;
        for (char c : dead) anyDead = anyDead || c;
        if (anyDead) {
            const int nAcc = (int)work.accepted.size();
            std::vector<int> newIdx((std::size_t)nAcc, -1);
            std::vector<Region> keep;
            for (int i = 0; i < nAcc; ++i) {
                if (i < (int)dead.size() && dead[(std::size_t)i]) continue;
                newIdx[(std::size_t)i] = (int)keep.size();
                keep.push_back(std::move(work.accepted[(std::size_t)i]));
            }
            for (Region& r : keep)
                if (r.origin == Origin::ChamferCone && r.filletNbrA >= 0 && r.filletNbrA < nAcc)
                    r.filletNbrA = newIdx[(std::size_t)r.filletNbrA];
            work.accepted.swap(keep);
        }
        return true;
    } catch (...) {
        return true;
    }
}

}  // namespace refit
}  // namespace stl2step
