// Sphere corners by relation (D-fc-5.11 F8). A corner where three certified
// rounds of one radius meet is a sphere: its centre is the point the rounds'
// spines share and its radius is theirs, so nothing is fitted to the corner's
// own triangles. The claimed set is every free triangle inside the corner
// sector whose three vertices lie within tau of that sphere; it is claimed only
// when it is one disk bounded by the three rounds alone. Runs after the tori
// (H1), before A3, so a torus round is a round like a cylinder.
// SPDX-License-Identifier: MIT

#include "refit_relational_sphere.hpp"

#include "refit_internal.hpp"
#include "refit_profile_cert.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {

namespace {
constexpr double kTwoPi = 6.28318530717958647692528676655900577;
}  // namespace

double spineDistance(const RoundSpine& s, const gp_XYZ& x) {
    const gp_XYZ d = x - s.o;
    const double h = d.Dot(s.a);
    const gp_XYZ rho = d - s.a * h;
    if (!s.circle) return rho.Modulus();
    return std::hypot(rho.Modulus() - s.R, h);
}

bool spineTangentAt(const RoundSpine& s, const gp_XYZ& x, gp_XYZ& t) {
    if (!s.circle) {
        t = s.a;
        return true;
    }
    const gp_XYZ d = x - s.o;
    const gp_XYZ rho = d - s.a * d.Dot(s.a);
    const double m = rho.Modulus();
    if (!(m > 0.0)) return false;
    t = s.a.Crossed(rho / m);
    return true;
}

bool lineSpinesPoint(const std::vector<RoundSpine>& spines, gp_XYZ& c) {
    double M[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double b[3] = {0, 0, 0};
    int n = 0;
    for (const RoundSpine& s : spines) {
        if (s.circle) continue;
        const double a[3] = {s.a.X(), s.a.Y(), s.a.Z()};
        const double o[3] = {s.o.X(), s.o.Y(), s.o.Z()};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const double pij = (i == j ? 1.0 : 0.0) - a[i] * a[j];
                M[i][j] += pij;
                b[i] += pij * o[j];
            }
        ++n;
    }
    if (n < 2) return false;
    const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
                       M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                       M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (!(std::fabs(det) > 0.0) || !std::isfinite(det)) return false;
    double x[3];
    for (int k = 0; k < 3; ++k) {
        double A[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) A[i][j] = (j == k) ? b[i] : M[i][j];
        x[k] = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) /
               det;
    }
    if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2])) return false;
    c = gp_XYZ(x[0], x[1], x[2]);
    return true;
}

SphereCornerFit sphereOfRounds(const std::array<RoundSpine, 3>& rounds, double tau) {
    SphereCornerFit f;
    if (!(tau > 0.0)) return f;
    std::vector<RoundSpine> v(rounds.begin(), rounds.end());
    if (!lineSpinesPoint(v, f.c)) return f;
    double lo = rounds[0].r, hi = rounds[0].r, sum = 0.0;
    for (const RoundSpine& s : rounds) {
        lo = std::min(lo, s.r);
        hi = std::max(hi, s.r);
        sum += s.r;
        f.spineDev = std::max(f.spineDev, spineDistance(s, f.c));
    }
    f.r = sum / 3.0;
    f.radiusSpread = hi - lo;
    if (f.radiusSpread > tau || f.spineDev > tau || !(f.r > tau)) {
        f.why = "not-blend";
        return f;
    }
    f.ok = true;
    f.why = "accept";
    return f;
}

double planeResidual(const std::vector<gp_XYZ>& pts) {
    if (pts.size() < 4) return 0.0;
    gp_XYZ m(0, 0, 0);
    for (const gp_XYZ& p : pts) m += p;
    m /= (double)pts.size();
    std::array<std::array<double, 3>, 3> C{};
    for (const gp_XYZ& p : pts) {
        const double d[3] = {p.X() - m.X(), p.Y() - m.Y(), p.Z() - m.Z()};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) C[i][j] += d[i] * d[j];
    }
    std::array<double, 3> ev{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(C, ev, evec)) return HUGE_VAL;
    int k = 0;
    for (int i = 1; i < 3; ++i)
        if (ev[i] < ev[k]) k = i;
    gp_XYZ n(evec[k][0], evec[k][1], evec[k][2]);  // row k: eigenvector of eval[k]
    const double nm = n.Modulus();
    if (!(nm > 0.0)) return HUGE_VAL;
    n /= nm;
    double mx = 0.0;
    for (const gp_XYZ& p : pts) mx = std::max(mx, std::fabs((p - m).Dot(n)));
    return mx;
}

bool pointsPlanar(const std::vector<gp_XYZ>& pts, double tau) {
    return planeResidual(pts) <= tau;
}

RoundCert certifyCylinderRound(const RoundSpine& spine,
                               const std::vector<std::array<gp_XYZ, 3>>& tris, double tau) {
    RoundCert rc;
    if (spine.circle || !(tau > 0.0) || tris.empty()) return rc;
    gp_XYZ e1 = std::fabs(spine.a.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    e1 = e1 - spine.a * e1.Dot(spine.a);
    e1 /= e1.Modulus();
    const gp_XYZ e2 = spine.a.Crossed(e1);
    auto dev = [&](const gp_XYZ& x) { return std::fabs(spineDistance(spine, x) - spine.r); };
    std::vector<char> inI(tris.size(), 0);
    std::vector<double> psi, rho, zs;
    double zlo = HUGE_VAL, zhi = -HUGE_VAL;
    for (std::size_t i = 0; i < tris.size(); ++i) {
        bool in = true;
        for (const gp_XYZ& x : tris[i]) in = in && dev(x) <= tau;
        if (!in) continue;
        inI[i] = 1;
        ++rc.nI;
        for (int k = 0; k < 3; ++k) {
            const gp_XYZ& x = tris[i][(std::size_t)k];
            const gp_XYZ d = x - spine.o;
            const double z = d.Dot(spine.a);
            const gp_XYZ q = d - spine.a * z;
            psi.push_back(std::atan2(q.Dot(e2), q.Dot(e1)));
            rho.push_back(q.Modulus());
            zs.push_back(z);
            zlo = std::min(zlo, z);
            zhi = std::max(zhi, z);
            rc.certDev = std::max(rc.certDev, dev(x));
            const gp_XYZ mid = (x + tris[i][(std::size_t)((k + 1) % 3)]) * 0.5;
            rc.s = std::max(rc.s, dev(mid));
        }
    }
    if (rc.nI == 0) return rc;
    rc.stations = countColumns(psi, rho, tau);
    rc.levels = countLevels(zs, tau);
    for (std::size_t i = 0; i < tris.size(); ++i) {
        if (inI[i]) continue;
        bool one = false, all = true;
        for (const gp_XYZ& x : tris[i]) {
            const double d = dev(x);
            one = one || d <= tau;
            all = all && d <= rc.s;
            const double z = (x - spine.o).Dot(spine.a);
            all = all && z >= zlo - tau && z <= zhi + tau;
            rc.coverDev = std::max(rc.coverDev, d);
        }
        if (one && all) ++rc.nCover;
        else ++rc.nOther;
    }
    rc.ok = rc.nOther == 0 && rc.stations >= minStations(kParamsCircle) && rc.levels >= 2;
    return rc;
}

bool sphereCylinderContact(const gp_XYZ& c, double r, const gp_XYZ& o, const gp_XYZ& a,
                           double R, double tol, double& v) {
    RoundSpine s;
    s.o = o;
    s.a = a;
    s.r = R;
    if (!(tol >= 0.0) || spineDistance(s, c) > tol || std::fabs(R - r) > tol) return false;
    v = (c - o).Dot(a);
    return true;
}

double sphereIsoBound(const gp_XYZ& centre, double r, const gp_XYZ& sx, const gp_XYZ& sy,
                      const gp_XYZ& sz, const gp_XYZ& c0, double rho, const gp_XYZ& xc,
                      const gp_XYZ& yc, bool meridian, double u0, double v0, double s) {
    gp_XYZ s0, sc, sn;
    if (!meridian) {
        const gp_XYZ m = sx * std::cos(u0) + sy * std::sin(u0);
        const gp_XYZ w = sx * (-std::sin(u0)) + sy * std::cos(u0);
        s0 = centre + sz * (r * std::sin(v0));
        sc = m * (r * std::cos(v0));
        sn = w * (r * std::cos(v0) * s);
    } else {
        const gp_XYZ m = sx * std::cos(u0) + sy * std::sin(u0);
        s0 = centre;
        sc = (m * std::cos(v0) + sz * std::sin(v0)) * r;
        sn = (m * (-std::sin(v0)) + sz * std::cos(v0)) * (r * s);
    }
    const gp_XYZ P = c0 - s0;
    const gp_XYZ Q = xc * rho - sc;
    const gp_XYZ R = yc * rho - sn;
    const double a = Q.Dot(Q), b = Q.Dot(R), c = R.Dot(R);
    const double lam = 0.5 * (a + c) + std::sqrt(0.25 * (a - c) * (a - c) + b * b);
    return P.Modulus() + std::sqrt(std::max(lam, 0.0));
}

SphereIsoLine sphereIsoLineOf(const gp_XYZ& centre, double r, const gp_XYZ& sx, const gp_XYZ& sy,
                              const gp_XYZ& sz, const gp_XYZ& c0, double rho, const gp_XYZ& xc,
                              const gp_XYZ& yc, double tMid) {
    SphereIsoLine L;
    const gp_XYZ n = xc.Crossed(yc);
    auto wrap = [](double u) {
        while (u < 0.0) u += kTwoPi;
        while (u >= kTwoPi) u -= kTwoPi;
        return u;
    };
    if (std::fabs(n.Dot(sz)) >= std::hypot(n.Dot(sx), n.Dot(sy))) {
        // Parallel: u(t) = azimuth of Xc + s t.
        L.meridian = false;
        L.s = n.Dot(sz) >= 0.0 ? 1.0 : -1.0;
        L.u0 = wrap(std::atan2(xc.Dot(sy), xc.Dot(sx)));
        L.v0 = std::asin(std::max(-1.0, std::min(1.0, (c0 - centre).Dot(sz) / r)));
    } else {
        // Meridian through the circle's mid point: v(t) = beta + s t.
        const gp_XYZ pm = c0 + (xc * std::cos(tMid) + yc * std::sin(tMid)) * rho - centre;
        gp_XYZ m = pm - sz * pm.Dot(sz);
        const double mm = m.Modulus();
        if (!(mm > 0.0)) return L;
        m /= mm;
        L.meridian = true;
        L.u0 = wrap(std::atan2(m.Dot(sy), m.Dot(sx)));
        const double beta = std::atan2(xc.Dot(sz), xc.Dot(m));
        const gp_XYZ yPlus = m * (-std::sin(beta)) + sz * std::cos(beta);
        L.s = yc.Dot(yPlus) >= 0.0 ? 1.0 : -1.0;
        L.v0 = beta;
    }
    L.bound = sphereIsoBound(centre, r, sx, sy, sz, c0, rho, xc, yc, L.meridian, L.u0, L.v0, L.s);
    // The line's own branch: the value at tMid in [0, 2 pi) (u of a parallel)
    // or in [-pi, pi) (v of a meridian). A whole-period shift is exact.
    if (!L.meridian) {
        const double um = L.u0 + L.s * tMid;
        L.u0 += wrap(um) - um;
    } else {
        const double vm = L.v0 + L.s * tMid;
        L.v0 += (wrap(vm + 0.5 * kTwoPi) - 0.5 * kTwoPi) - vm;
    }
    return L;
}

double sphereCapVolume(const gp_XYZ& centre, double r, const gp_XYZ& p, const gp_XYZ& q,
                       const gp_XYZ& s) {
    const gp_XYZ a = p - centre, b = q - centre, d = s - centre;
    const double la = a.Modulus(), lb = b.Modulus(), ld = d.Modulus();
    const double triple = std::fabs(a.Dot(b.Crossed(d)));
    const double den = la * lb * ld + a.Dot(b) * ld + a.Dot(d) * lb + b.Dot(d) * la;
    const double omega = 2.0 * std::atan2(triple, den);
    return r * r * r * omega / 3.0 - triple / 6.0;
}

namespace {

bool diagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

void diagPrintf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
}

gp_XYZ ptOf(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.pts[mv.tris[gt][k]];
}

double triAreaNormal(const MeshView& mv, int t, gp_XYZ& n) {
    const gp_XYZ a = ptOf(mv, t, 0), b = ptOf(mv, t, 1), c = ptOf(mv, t, 2);
    gp_XYZ x = (b - a).Crossed(c - a);
    const double m = x.Modulus();
    n = m > 0.0 ? x / m : gp_XYZ(0, 0, 0);
    return 0.5 * m;
}

// D-fc-4.3. An accepted CylGrow region whose own maxVertexDev exceeds tau does
// not certify its surface (D-fc-1.4); a certified band may rob it.
bool uncertifiedCylGrow(const Region& r, double tau) {
    return r.origin == Origin::CylGrow && r.maxVertexDev > tau;
}

void refitPlaneProv(const MeshView& mv, Provisional& pr) {
    gp_Ax3 pl;
    if (!pr.tris.empty() && pcaPlane(mv, pr.tris, pl)) pr.plane = pl;
    double area = 0.0;
    std::vector<int> vs;
    for (int t : pr.tris) {
        gp_XYZ n;
        area += triAreaNormal(mv, t, n);
        const int gt = mv.compTris ? mv.compTris[t] : t;
        for (int k = 0; k < 3; ++k) vs.push_back(mv.tris[gt][k]);
    }
    std::sort(vs.begin(), vs.end());
    vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
    pr.area = area;
    const gp_XYZ n = pr.plane.Direction().XYZ();
    const gp_XYZ o = pr.plane.Location().XYZ();
    double mx = 0.0, ss = 0.0;
    for (int g : vs) {
        const double d = std::fabs((mv.pts[g] - o).Dot(n));
        mx = std::max(mx, d);
        ss += d * d;
    }
    pr.maxVertexDev = mx;
    pr.rmsVertexDev = vs.empty() ? 0.0 : std::sqrt(ss / (double)vs.size());
}

struct Round {
    int acc = -1;
    RoundSpine s;
    bool certified = false;
    RoundCert cert;
};

}  // namespace

bool claimSphereCorners(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                        SegmentWork& work) {
    (void)p;
    (void)tol;
    try {
        const double tau = tauOf(mv);
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;
        const int nTri = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const bool dg = diagOn();

        // --- rounds: cylinders and torus blends; certified per D-fc-1.4 ------
        std::vector<Round> rounds;
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai) {
            const Region& R = work.accepted[(std::size_t)ai];
            if (R.tris.empty()) continue;
            Round rd;
            rd.acc = ai;
            rd.s.o = R.ax.Location().XYZ();
            rd.s.a = R.ax.Direction().XYZ();
            if (R.type == SurfType::Cylinder) {
                rd.s.circle = false;
                rd.s.r = R.radius;
            } else if (R.type == SurfType::Torus && R.origin == Origin::TorusBlend) {
                rd.s.circle = true;
                rd.s.R = R.radius;
                rd.s.r = R.radius2;
            } else {
                continue;
            }
            if (!(rd.s.r > tau)) continue;
            rd.certified = anchorCertified(tau, R.maxVertexDev);
            // A relational band's surface is certified on its tau-preimage I;
            // its chord-bounded cover does not un-certify it (D-fc-3.2.7). A
            // B1 region is certified by its own residual alone (D-fc-4.3).
            if (!rd.certified && !rd.s.circle && R.origin != Origin::CylGrow) {
                std::vector<std::array<gp_XYZ, 3>> pts;
                pts.reserve(R.tris.size());
                for (int t : R.tris) {
                    if (t < 0 || t >= nTri) continue;
                    pts.push_back({ptOf(mv, t, 0), ptOf(mv, t, 1), ptOf(mv, t, 2)});
                }
                rd.cert = certifyCylinderRound(rd.s, pts, tau);
                rd.certified = rd.cert.ok;
            }
            rounds.push_back(rd);
        }
        if (rounds.size() < 3) return true;

        // --- corners: two line spines that meet, and every round through it --
        std::vector<std::array<int, 3>> corners;
        for (std::size_t i = 0; i < rounds.size(); ++i) {
            if (rounds[i].s.circle) continue;
            for (std::size_t j = i + 1; j < rounds.size(); ++j) {
                if (rounds[j].s.circle) continue;
                if (std::fabs(rounds[i].s.r - rounds[j].s.r) > tau) continue;
                gp_XYZ c;
                if (!lineSpinesPoint({rounds[i].s, rounds[j].s}, c)) continue;
                if (spineDistance(rounds[i].s, c) > tau || spineDistance(rounds[j].s, c) > tau)
                    continue;
                for (std::size_t k = 0; k < rounds.size(); ++k) {
                    if (k == i || k == j) continue;
                    if (std::fabs(rounds[k].s.r - rounds[i].s.r) > tau) continue;
                    if (spineDistance(rounds[k].s, c) > tau) continue;
                    std::array<int, 3> key = {(int)i, (int)j, (int)k};
                    std::sort(key.begin(), key.end());
                    corners.push_back(key);
                }
            }
        }
        std::sort(corners.begin(), corners.end());
        corners.erase(std::unique(corners.begin(), corners.end()), corners.end());

        int nAccept = 0, nRefused = 0;
        if (corners.empty()) {
            if (dg)
                diagPrintf("DIAG_REL_SUM stage=sphere rounds=%zu corners=0 accept=0 refused=0\n",
                           rounds.size());
            return true;
        }

        // --- triangle state -----------------------------------------------
        std::vector<int> triProv((std::size_t)nTri, -1);
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi)
            for (int t : work.provisionals[(std::size_t)pi].tris)
                if (t >= 0 && t < nTri) triProv[(std::size_t)t] = pi;
        std::vector<int> triAcc((std::size_t)nTri, -1);
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai)
            for (int t : work.accepted[(std::size_t)ai].tris)
                if (t >= 0 && t < nTri) triAcc[(std::size_t)t] = ai;
        std::vector<char> robbable(work.accepted.size(), 0);
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai)
            robbable[(std::size_t)ai] = uncertifiedCylGrow(work.accepted[(std::size_t)ai], tau);
        auto candidate = [&](int t) -> bool {
            const int ai = triAcc[(std::size_t)t];
            if (ai >= 0) return robbable[(std::size_t)ai] != 0;
            const int pi = triProv[(std::size_t)t];
            return pi >= 0 && work.provisionals[(std::size_t)pi].claim == ProvClaim::Unclaimed;
        };

        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{{-1, -1}});
        for (int t = 0; t < nTri; ++t)
            for (int s = 0; s < 3; ++s) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || e >= nEdge) continue;
                auto& et = edgeTris[(std::size_t)e];
                if (et[0] < 0) et[0] = t;
                else if (et[0] != t) et[1] = t;
            }
        auto otherTri = [&](int e, int t) -> int {
            if (e < 0 || e >= nEdge) return -1;
            const auto& et = edgeTris[(std::size_t)e];
            if (et[0] == t) return et[1];
            if (et[1] == t) return et[0];
            return -1;
        };

        std::vector<int> dissolved;
        std::vector<char> claimed((std::size_t)nTri, 0);

        for (const std::array<int, 3>& key : corners) {
            const Round* rd[3] = {&rounds[(std::size_t)key[0]], &rounds[(std::size_t)key[1]],
                                  &rounds[(std::size_t)key[2]]};
            const SphereCornerFit fit = sphereOfRounds({rd[0]->s, rd[1]->s, rd[2]->s}, tau);
            const gp_XYZ c = fit.c;
            const double r = fit.r;
            int accIdx[3];
            for (int q = 0; q < 3; ++q) accIdx[q] = rd[q]->acc;
            const char* why = fit.why;
            std::vector<int> band;
            int nV = 0, comps = 0, chi = 0, foreign = 0, nRob = 0;
            int touch[3] = {0, 0, 0};
            double maxDev = 0.0, planeRes = 0.0;
            auto emit = [&](const char* w) {
                if (!dg) return;
                int seed = band.empty() ? -1 : band.front();
                diagPrintf("DIAG_REL_BAND stage=sphere mode=corner seed=%d nTri=%zu nV=%d "
                           "stations=%d levels=0 params=%d fit=anchored radius=%.9g "
                           "centre=%.6f,%.6f,%.6f maxDev=%.6g maxDevTau=%.6g lineRes=%.6g "
                           "spineDevTau=%.6g radiusSpreadTau=%.6g rounds=%d,%d,%d "
                           "roundOrigin=%d,%d,%d roundDevTau=%.4g,%.4g,%.4g touch=%d,%d,%d "
                           "roundCertI=%d,%d,%d roundCover=%d,%d,%d comps=%d chi=%d foreign=%d rob=%d why=%s\n",
                           seed, band.size(), nV, nV, kParamsSphere, r, c.X(), c.Y(), c.Z(),
                           maxDev, maxDev / tau, planeRes / tau, fit.spineDev / tau,
                           fit.radiusSpread / tau, accIdx[0], accIdx[1], accIdx[2],
                           (int)work.accepted[(std::size_t)accIdx[0]].origin,
                           (int)work.accepted[(std::size_t)accIdx[1]].origin,
                           (int)work.accepted[(std::size_t)accIdx[2]].origin,
                           work.accepted[(std::size_t)accIdx[0]].maxVertexDev / tau,
                           work.accepted[(std::size_t)accIdx[1]].maxVertexDev / tau,
                           work.accepted[(std::size_t)accIdx[2]].maxVertexDev / tau, touch[0],
                           touch[1], touch[2], rd[0]->cert.nI, rd[1]->cert.nI, rd[2]->cert.nI,
                           rd[0]->cert.nCover, rd[1]->cert.nCover, rd[2]->cert.nCover, comps, chi,
                           foreign, nRob, w);
            };
            auto refuse = [&](const char* w) {
                ++nRefused;
                emit(w);
            };
            if (!fit.ok) {
                refuse(why);
                continue;
            }
            // D-fc-1.4: every round is a certified surface.
            if (!rd[0]->certified || !rd[1]->certified || !rd[2]->certified) {
                refuse("no-anchor");
                continue;
            }
            // Each round ends at the corner: its contact plane (through c, normal
            // to its spine there) has the round's whole body on one side.
            gp_XYZ tan[3];
            double side[3] = {0, 0, 0};
            bool ends = true;
            for (int q = 0; q < 3 && ends; ++q) {
                if (!spineTangentAt(rd[q]->s, c, tan[q])) {
                    ends = false;
                    break;
                }
                const Region& R = work.accepted[(std::size_t)rd[q]->acc];
                double sum = 0.0;
                for (int t : R.tris)
                    for (int k = 0; k < 3; ++k) sum += (ptOf(mv, t, k) - c).Dot(tan[q]);
                side[q] = sum >= 0.0 ? 1.0 : -1.0;
                for (int t : R.tris)
                    for (int k = 0; k < 3; ++k)
                        if (!onRoundSide(ptOf(mv, t, k), c, tan[q], side[q], tau)) ends = false;
            }
            if (!ends) {
                refuse("not-blend");
                continue;
            }
            // I: every free triangle in the sector with 3 vertices within tau.
            for (int t = 0; t < nTri; ++t) {
                if (claimed[(std::size_t)t] || !candidate(t)) continue;
                bool in = true;
                for (int k = 0; k < 3 && in; ++k) {
                    const gp_XYZ x = ptOf(mv, t, k);
                    if (std::fabs((x - c).Modulus() - r) > tau) in = false;
                    for (int q = 0; q < 3 && in; ++q)
                        if (!onPatchSide(x, c, tan[q], side[q], tau)) in = false;
                }
                if (in) band.push_back(t);
            }
            if (band.empty()) {
                refuse("underdetermined");
                continue;
            }
            std::vector<char> inBand((std::size_t)nTri, 0);
            for (int t : band) inBand[(std::size_t)t] = 1;
            std::vector<int> gv;
            for (int t : band) {
                const int gt = mv.compTris ? mv.compTris[t] : t;
                for (int k = 0; k < 3; ++k) gv.push_back(mv.tris[gt][k]);
            }
            std::sort(gv.begin(), gv.end());
            gv.erase(std::unique(gv.begin(), gv.end()), gv.end());
            nV = (int)gv.size();
            std::vector<gp_XYZ> vp;
            vp.reserve(gv.size());
            double ss = 0.0;
            for (int g : gv) {
                vp.push_back(mv.pts[g]);
                const double d = std::fabs((mv.pts[g] - c).Modulus() - r);
                maxDev = std::max(maxDev, d);
                ss += d * d;
            }
            if (nV < minStations(kParamsSphere)) {
                refuse("underdetermined");
                continue;
            }
            // Plane first (fewest parameters).
            planeRes = planeResidual(vp);
            if (pointsPlanar(vp, tau)) {
                refuse("line");
                continue;
            }
            // One disk bounded by the three rounds alone.
            int nE = 0, nBoundary = 0, contactOff = 0;
            {
                std::vector<int> es;
                for (int t : band)
                    for (int s = 0; s < 3; ++s) es.push_back(mv.triEdges[t][s]);
                std::sort(es.begin(), es.end());
                es.erase(std::unique(es.begin(), es.end()), es.end());
                nE = (int)es.size();
                chi = nV - nE + (int)band.size();
                std::vector<char> seen((std::size_t)nTri, 0);
                for (int t0 : band) {
                    if (seen[(std::size_t)t0]) continue;
                    ++comps;
                    std::vector<int> st = {t0};
                    seen[(std::size_t)t0] = 1;
                    while (!st.empty()) {
                        const int u = st.back();
                        st.pop_back();
                        for (int s = 0; s < 3; ++s) {
                            const int o = otherTri(mv.triEdges[u][s], u);
                            if (o >= 0 && inBand[(std::size_t)o] && !seen[(std::size_t)o]) {
                                seen[(std::size_t)o] = 1;
                                st.push_back(o);
                            }
                        }
                    }
                }
                bool sameSurf = false;
                for (int t : band)
                    for (int s = 0; s < 3; ++s) {
                        const int o = otherTri(mv.triEdges[t][s], t);
                        if (o >= 0 && inBand[(std::size_t)o]) continue;
                        ++nBoundary;
                        const int ai = o >= 0 ? triAcc[(std::size_t)o] : -1;
                        int which = -1;
                        for (int q = 0; q < 3; ++q)
                            if (ai >= 0 && ai == rd[q]->acc) which = q;
                        if (which >= 0) {
                            ++touch[which];
                            // The round meets the sphere on its certified
                            // surface: its contact triangle is within tau.
                            for (int k = 0; k < 3; ++k)
                                if (std::fabs(spineDistance(rd[which]->s, ptOf(mv, o, k)) -
                                              rd[which]->s.r) > tau)
                                    ++contactOff;
                            continue;
                        }
                        ++foreign;
                        if (ai >= 0 && !robbable[(std::size_t)ai]) {
                            // A certified neighbour lying on this sphere is a second
                            // face of one surface (D-fc-1.6(11)).
                            bool on = true;
                            for (int t2 : work.accepted[(std::size_t)ai].tris)
                                for (int k = 0; k < 3 && on; ++k)
                                    if (std::fabs((ptOf(mv, t2, k) - c).Modulus() - r) > tau)
                                        on = false;
                            if (on) {
                                sameSurf = true;
                                if (dg) {
                                    const Region& sr = work.accepted[(std::size_t)ai];
                                    diagPrintf("DIAG_REL_SPHERE_SAMESURF seed=%d region=%d "
                                               "origin=%d type=%d nTri=%zu radius=%.9g "
                                               "devTau=%.6g\n",
                                               band.front(), ai, (int)sr.origin, (int)sr.type,
                                               sr.tris.size(), sr.radius,
                                               sr.maxVertexDev / tau);
                                }
                            }
                        }
                    }
                if (sameSurf) {
                    refuse("same-surface");
                    continue;
                }
            }
            if (comps != 1 || chi != 1 || foreign != 0 || touch[0] == 0 || touch[1] == 0 ||
                touch[2] == 0) {
                refuse("not-blend");
                continue;
            }
            if (contactOff != 0) {
                refuse("vertex-off");
                continue;
            }

            // --- accept -----------------------------------------------------
            gp_XYZ msum(0, 0, 0);
            double sigma = 0.0, areaSum = 0.0, dvol = 0.0;
            for (int t : band) {
                gp_XYZ nt;
                const double at = triAreaNormal(mv, t, nt);
                const gp_XYZ a = ptOf(mv, t, 0), b = ptOf(mv, t, 1), d = ptOf(mv, t, 2);
                gp_XYZ dir = (a + b + d) / 3.0 - c;
                const double dm = dir.Modulus();
                if (dm > 0.0) dir /= dm;
                msum += dir * at;
                sigma += at * nt.Dot(dir);
                areaSum += at;
                dvol += sphereCapVolume(c, r, a, b, d);
            }
            const double mm = msum.Modulus();
            if (!(mm > 0.0)) {
                refuse("degenerate");
                continue;
            }
            const gp_XYZ mdir = msum / mm;
            // The sphere's frame. When one round's spine tangent is
            // perpendicular to the other two (within tau / r: the tilt that
            // moves a contact circle by tau), the axis is that tangent pointing
            // into the patch: its contact is the equator, the other two are
            // meridians and meet at the pole -- every contact is an iso-line,
            // so every pcurve is an exact line. Otherwise the patch's middle
            // sits on the equator. Either way the patch is centred on u = pi,
            // away from the seam.
            int polar = -1;
            for (int q = 0; q < 3 && polar < 0; ++q)
                if (std::fabs(tan[q].Dot(tan[(q + 1) % 3])) <= tau / r &&
                    std::fabs(tan[q].Dot(tan[(q + 2) % 3])) <= tau / r)
                    polar = q;
            gp_XYZ zd, xd;
            if (polar >= 0) {
                zd = tan[polar] * (-side[polar]);
                zd /= zd.Modulus();
                gp_XYZ mp = mdir - zd * mdir.Dot(zd);
                const double mpm = mp.Modulus();
                if (!(mpm > 0.0)) {
                    refuse("degenerate");
                    continue;
                }
                xd = (mp / mpm).Reversed();
            } else {
                const double ax3[3] = {std::fabs(mdir.X()), std::fabs(mdir.Y()),
                                       std::fabs(mdir.Z())};
                int ke = 0;
                for (int i = 1; i < 3; ++i)
                    if (ax3[i] < ax3[ke]) ke = i;
                gp_XYZ e(ke == 0 ? 1.0 : 0.0, ke == 1 ? 1.0 : 0.0, ke == 2 ? 1.0 : 0.0);
                zd = mdir.Crossed(e);
                zd /= zd.Modulus();
                xd = mdir.Reversed();
            }
            const bool outward = sigma > 0.0;

            Region reg;
            reg.id = 0;
            reg.type = SurfType::Sphere;
            reg.origin = Origin::SphereCorner;
            reg.ax = gp_Ax3(gp_Pnt(c), gp_Dir(zd), gp_Dir(xd));
            reg.radius = r;
            reg.radius2 = 0.0;
            {
                double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300;
                const gp_XYZ yd = zd.Crossed(xd);
                for (const gp_XYZ& x : vp) {
                    const gp_XYZ d = x - c;
                    const double u = std::atan2(d.Dot(yd), d.Dot(xd));
                    const double uu = u < 0.0 ? u + kTwoPi : u;
                    const double v = std::atan2(d.Dot(zd), std::hypot(d.Dot(xd), d.Dot(yd)));
                    u0 = std::min(u0, uu);
                    u1 = std::max(u1, uu);
                    v0 = std::min(v0, v);
                    v1 = std::max(v1, v);
                }
                reg.uMin = u0;
                reg.uMax = u1;
                reg.vMin = v0;
                reg.vMax = v1;
            }
            reg.closed360 = false;
            reg.outwardNormal = outward;
            reg.tris = band;
            reg.maxVertexDev = maxDev;
            reg.rmsVertexDev = std::sqrt(ss / (double)gv.size());
            reg.nSides = 0;
            reg.chordSagitta = 0.0;
            reg.dVolPredicted = (outward ? 1.0 : -1.0) * dvol;
            reg.maxVertexSnap = 0.0;
            reg.reject = Reject::None;
            reg.builtAs = BuiltAs::NotBuilt;
            reg.filletNbrA = -1;
            reg.filletNbrB = -1;

            // D-fc-4.3: every uncertified region the patch robbed dissolves.
            std::vector<int> robbed;
            for (int t : band)
                if (triAcc[(std::size_t)t] >= 0) robbed.push_back(triAcc[(std::size_t)t]);
            std::sort(robbed.begin(), robbed.end());
            robbed.erase(std::unique(robbed.begin(), robbed.end()), robbed.end());
            nRob = (int)robbed.size();
            for (int ai : robbed) {
                Region& rr = work.accepted[(std::size_t)ai];
                std::vector<char> inR((std::size_t)nTri, 0);
                int nStolen = 0;
                for (int t : rr.tris) {
                    if (t < 0 || t >= nTri) continue;
                    inR[(std::size_t)t] = 1;
                    if (inBand[(std::size_t)t]) ++nStolen;
                }
                for (Provisional& pr : work.provisionals) {
                    if (pr.claim != ProvClaim::ConsumedCylinder || pr.tris.empty()) continue;
                    bool own = true;
                    for (int t : pr.tris)
                        if (t < 0 || t >= nTri || !inR[(std::size_t)t]) {
                            own = false;
                            break;
                        }
                    if (own) pr.claim = ProvClaim::Unclaimed;
                }
                for (int t : rr.tris)
                    if (t >= 0 && t < nTri) triAcc[(std::size_t)t] = -1;
                if (dg)
                    diagPrintf("DIAG_REL_ROB stage=sphere seed=%d region=%d radius=%.9g "
                               "devTau=%.9g nTri=%d stolen=%d released=%d\n",
                               band.front(), ai, rr.radius, rr.maxVertexDev / tau,
                               (int)rr.tris.size(), nStolen, (int)rr.tris.size() - nStolen);
                rr.tris.clear();
                robbable[(std::size_t)ai] = 0;
                dissolved.push_back(ai);
            }
            // Triangle granularity (D-fc-1.6(12)): the patch leaves every
            // provisional; a remainder is refit and split into components.
            const int nProvNow = (int)work.provisionals.size();
            for (int pi = 0; pi < nProvNow; ++pi) {
                Provisional& pr = work.provisionals[(std::size_t)pi];
                std::vector<int> keep;
                bool any = false;
                for (int t : pr.tris) {
                    if (t >= 0 && t < nTri && inBand[(std::size_t)t]) any = true;
                    else keep.push_back(t);
                }
                if (!any) continue;
                if (keep.empty()) {
                    pr.claim = ProvClaim::ConsumedCylinder;
                    continue;
                }
                std::vector<char> inKeep((std::size_t)nTri, 0);
                for (int t : keep) inKeep[(std::size_t)t] = 1;
                std::vector<std::vector<int>> comp;
                std::vector<char> mark((std::size_t)nTri, 0);
                for (int t0 : keep) {
                    if (mark[(std::size_t)t0]) continue;
                    std::vector<int> cc, st = {t0};
                    mark[(std::size_t)t0] = 1;
                    while (!st.empty()) {
                        const int u = st.back();
                        st.pop_back();
                        cc.push_back(u);
                        for (int s = 0; s < 3; ++s) {
                            const int o = otherTri(mv.triEdges[u][s], u);
                            if (o >= 0 && inKeep[(std::size_t)o] && !mark[(std::size_t)o]) {
                                mark[(std::size_t)o] = 1;
                                st.push_back(o);
                            }
                        }
                    }
                    std::sort(cc.begin(), cc.end());
                    comp.push_back(std::move(cc));
                }
                std::sort(comp.begin(), comp.end(),
                          [](const std::vector<int>& a, const std::vector<int>& b) {
                              return a.front() < b.front();
                          });
                const ProvClaim keepClaim =
                    pr.claim == ProvClaim::ConsumedCylinder ? ProvClaim::Unclaimed : pr.claim;
                pr.tris = comp[0];
                pr.claim = keepClaim;
                refitPlaneProv(mv, pr);
                for (std::size_t q = 1; q < comp.size(); ++q) {
                    Provisional np = work.provisionals[(std::size_t)pi];
                    np.tris = comp[q];
                    np.claim = keepClaim;
                    refitPlaneProv(mv, np);
                    work.provisionals.push_back(std::move(np));
                }
            }
            for (int t : band) {
                claimed[(std::size_t)t] = 1;
                triAcc[(std::size_t)t] = (int)work.accepted.size();
            }
            triProv.assign((std::size_t)nTri, -1);
            for (int pi = 0; pi < (int)work.provisionals.size(); ++pi)
                for (int t : work.provisionals[(std::size_t)pi].tris)
                    if (t >= 0 && t < nTri) triProv[(std::size_t)t] = pi;
            ++nAccept;
            emit("accept");
            work.accepted.push_back(std::move(reg));
            robbable.push_back(0);
        }

        if (dg)
            diagPrintf("DIAG_REL_SUM stage=sphere rounds=%zu corners=%zu accept=%d refused=%d "
                       "dissolved=%zu\n",
                       rounds.size(), corners.size(), nAccept, nRefused, dissolved.size());

        // D-fc-4.3. A dissolved region leaves the accepted list. Indices into
        // it (a chamfer cone's cylinder, a pending torus's absorbed regions)
        // follow the regions that stay.
        if (!dissolved.empty()) {
            const int nAcc = (int)work.accepted.size();
            std::vector<char> gone((std::size_t)nAcc, 0);
            for (int ai : dissolved)
                if (ai >= 0 && ai < nAcc) gone[(std::size_t)ai] = 1;
            std::vector<int> newIdx((std::size_t)nAcc, -1);
            std::vector<Region> keep;
            keep.reserve((std::size_t)nAcc);
            for (int i = 0; i < nAcc; ++i) {
                if (gone[(std::size_t)i]) continue;
                newIdx[(std::size_t)i] = (int)keep.size();
                keep.push_back(std::move(work.accepted[(std::size_t)i]));
            }
            for (Region& rg : keep) {
                if (rg.origin != Origin::ChamferCone || rg.filletNbrA < 0 || rg.filletNbrA >= nAcc)
                    continue;
                rg.filletNbrA = newIdx[(std::size_t)rg.filletNbrA];
            }
            for (TorusPendingReclaim& pack : work.torusPending)
                for (int& ai : pack.reclaimAccIdx)
                    if (ai >= 0 && ai < nAcc) ai = newIdx[(std::size_t)ai];
            work.accepted.swap(keep);
        }
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
