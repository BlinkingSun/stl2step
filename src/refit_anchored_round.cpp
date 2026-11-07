// Anchored quarter-round torus (D-gid-2.6 / D-gid-3.5).
//
// certifyTorus needs kTorusProfPts = 5 distinct profile points because a FREE
// profile circle has three parameters. A round whose circle is tangent to
// certified anchors has fewer: tangency to a certified coaxial cylinder or a
// certified plane normal to the axis fixes one parameter each. A three-row
// quarter-round between such a wall and such a plane has four stations and
// one free parameter (r): it is over-determined in the anchored mode although
// certifyTorus calls it under-determined. certifyTorus is not edited.
//
// Every threshold is tau = 2q or a count of parameters (D-fc-1.6).
// SPDX-License-Identifier: MIT

#include "refit_anchored_round.hpp"

#include "refit_internal.hpp"
#include "refit_profile_cert.hpp"
#include "refit_relational.hpp"
#include "refit_relational_sphere.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

#include <gp.hxx>

namespace stl2step {
namespace refit {

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

// Iteration cap of the anchored least-squares fits: fitProfileCircle's own
// (I5 determinism cap, not an acceptance threshold).
constexpr int kFitIters = 16;

// I5 cap on the dVolPredicted subdivision of one triangle (per side).
constexpr int kMaxVolSub = 64;

double circDevOf(const ProfPt& p, double cr, double cz, double r) {
    return std::fabs(std::hypot(p.rho - cr, p.z - cz) - r);
}

// Unsigned angle between two 2-D vectors.
double angleBetween(double ax, double az, double bx, double bz) {
    return std::atan2(std::fabs(ax * bz - az * bx), ax * bx + az * bz);
}

bool diagOn() {
    static const bool on = [] {
        const char* e = std::getenv("STL2STEP_DIAG_REL");
        return e && e[0] && e[0] != '0';
    }();
    return on;
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

int gidOf(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.tris[gt][k];
}

double triAreaNormal(const MeshView& mv, int t, gp_XYZ& n) {
    const gp_XYZ a = ptOf(mv, t, 0), b = ptOf(mv, t, 1), c = ptOf(mv, t, 2);
    gp_XYZ x = (b - a).Crossed(c - a);
    const double m = x.Modulus();
    if (m > 0.0) n = x / m;
    else n = gp_XYZ(0, 0, 0);
    return 0.5 * m;
}

// The round's frame: the donor cylinder's axis and XDirection (D-fc-4.6 F8: one
// seam azimuth per revolution axis).
struct RoundFrame {
    gp_XYZ o, a, e1, e2;
};

ProfPt profOf(const RoundFrame& f, const gp_XYZ& x, int lv) {
    const gp_XYZ d = x - f.o;
    ProfPt p;
    p.z = d.Dot(f.a);
    const gp_XYZ r = d - f.a * p.z;
    p.rho = r.Modulus();
    p.psi = std::atan2(r.Dot(f.e2), r.Dot(f.e1));
    if (p.psi < 0.0) p.psi += kTwoPi;
    p.lv = lv;
    return p;
}

// Closest point of the torus (frame f, tube centre (Rmaj, z0), radius r) to x.
gp_XYZ torusFoot(const RoundFrame& f, double Rmaj, double z0, double r, const gp_XYZ& x) {
    const gp_XYZ d = x - f.o;
    const double z = d.Dot(f.a);
    gp_XYZ rad = d - f.a * z;
    const double rho = rad.Modulus();
    if (rho > gp::Resolution()) rad /= rho;
    else rad = f.e1;
    const double dr = rho - Rmaj, dz = z - z0;
    const double m = std::hypot(dr, dz);
    if (!(m > gp::Resolution())) return x;
    const double fr = Rmaj + r * dr / m, fz = z0 + r * dz / m;
    return f.o + rad * fr + f.a * fz;
}

double tetVol(const gp_XYZ& a, const gp_XYZ& b, const gp_XYZ& c, const gp_XYZ& d) {
    return (b - a).Crossed(c - a).Dot(d - a) / 6.0;
}

// I9 for the round: the signed volume between one mesh triangle and the torus
// over it (positive where the torus lies on the triangle's outward side). The
// triangle is subdivided m x m and each sub-triangle's corners are projected
// onto the torus; a flat sub-triangle on the torus misses the surface by its
// sagitta L^2 / (8 r), so m is the smallest count that brings that sagitta to
// tau (the certificate's own resolution): the error is <= area * tau.
double torusTriVolume(const RoundFrame& f, double Rmaj, double z0, double r, double tau,
                      const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2) {
    const double L = std::max({(p1 - p0).Modulus(), (p2 - p1).Modulus(), (p0 - p2).Modulus()});
    int m = 1;
    if (tau > 0.0 && r > 0.0) {
        const double need = L / std::sqrt(8.0 * r * tau);
        if (std::isfinite(need)) m = std::max(1, std::min(kMaxVolSub, (int)std::ceil(need)));
    }
    const gp_XYZ e1 = (p1 - p0) / (double)m, e2 = (p2 - p0) / (double)m;
    std::vector<gp_XYZ> X, Y;
    X.reserve((std::size_t)((m + 1) * (m + 2) / 2));
    auto idx = [m](int i, int j) {  // row i (along e1), column j (along e2), i + j <= m
        return i * (m + 1) - i * (i - 1) / 2 + j;
    };
    for (int i = 0; i <= m; ++i)
        for (int j = 0; i + j <= m; ++j) {
            const gp_XYZ x = p0 + e1 * (double)i + e2 * (double)j;
            X.push_back(x);
            Y.push_back(torusFoot(f, Rmaj, z0, r, x));
        }
    auto prism = [&](int a, int b, int c) {
        const gp_XYZ &Xa = X[(std::size_t)a], &Xb = X[(std::size_t)b], &Xc = X[(std::size_t)c];
        const gp_XYZ &Ya = Y[(std::size_t)a], &Yb = Y[(std::size_t)b], &Yc = Y[(std::size_t)c];
        return tetVol(Xa, Xb, Xc, Ya) + tetVol(Xb, Xc, Ya, Yb) + tetVol(Xc, Ya, Yb, Yc);
    };
    double v = 0.0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; i + j < m; ++j) {
            v += prism(idx(i, j), idx(i + 1, j), idx(i, j + 1));
            if (i + j + 1 < m) v += prism(idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1));
        }
    return std::isfinite(v) ? v : 0.0;
}

}  // namespace

AnchoredRoundCert certifyAnchoredRound(const std::vector<ProfPt>& P,
                                       const std::vector<ProfileAnchor>& anchors, double tau) {
    AnchoredRoundCert c;
    c.anchors = (int)anchors.size();
    if (!(tau > 0.0) || P.empty()) {
        c.why = "underdetermined";
        return c;
    }
    if (anchors.empty()) {
        c.why = "free-round";
        return c;
    }
    if (anchors.size() > 2) {
        c.why = "over-anchored";
        return c;
    }
    c.params = kParamsCircle - (int)anchors.size();
    {
        std::vector<double> ps, rs;
        ps.reserve(P.size());
        rs.reserve(P.size());
        for (const ProfPt& p : P) {
            ps.push_back(p.psi);
            rs.push_back(p.rho);
        }
        c.stations = countProfPts(P, tau);
        c.cols = countColumns(ps, rs, tau);
    }
    if (c.stations < minStations(c.params) || c.cols < 2) {
        c.why = "underdetermined";
        return c;
    }
    c.determinate = true;
    // Fewest parameters first: a profile line is a cylinder, cone or plane.
    if (profileLineMaxResid(P) <= tau) {
        c.why = "line";
        return c;
    }
    // Convexity from the fit: the side of each anchor the free circle's centre
    // lies on.
    double fR = 0, fZ = 0, fr = 0;
    if (!fitProfileCircle(P, tau, fR, fZ, fr)) {
        c.why = "fit-singular";
        return c;
    }
    double side[2] = {0, 0};
    for (std::size_t i = 0; i < anchors.size(); ++i) {
        const ProfileAnchor& A = anchors[i];
        const double s = A.nr * fR + A.nz * fZ - A.d;
        if (!std::isfinite(s) || s == 0.0) {
            c.why = "degenerate";
            return c;
        }
        side[i] = s > 0.0 ? 1.0 : -1.0;
    }
    double cr = 0, cz = 0, r = 0;
    if (anchors.size() == 2) {
        // n_i . c = d_i + side_i r  ->  c(r) = c0 + r w.
        const ProfileAnchor &A = anchors[0], &B = anchors[1];
        const double det = A.nr * B.nz - A.nz * B.nr;
        if (!(std::fabs(det) > gp::Resolution())) {
            c.why = "parallel-anchors";
            return c;
        }
        const double c0r = (A.d * B.nz - A.nz * B.d) / det;
        const double c0z = (A.nr * B.d - A.d * B.nr) / det;
        const double wr = (side[0] * B.nz - A.nz * side[1]) / det;
        const double wz = (A.nr * side[1] - side[0] * B.nr) / det;
        r = fr;
        for (int it = 0; it < kFitIters; ++it) {
            double sjr = 0, sjj = 0;
            for (const ProfPt& p : P) {
                const double qx = p.rho - (c0r + r * wr), qz = p.z - (c0z + r * wz);
                const double dd = std::hypot(qx, qz);
                if (!(dd > gp::Resolution())) continue;
                const double res = dd - r;
                const double j = -(qx * wr + qz * wz) / dd - 1.0;
                sjr += j * res;
                sjj += j * j;
            }
            if (!(sjj > 0.0)) break;
            const double dr = -sjr / sjj;
            if (!std::isfinite(dr)) break;
            r += dr;
            if (std::fabs(dr) <= tau) break;
        }
        cr = c0r + r * wr;
        cz = c0z + r * wz;
    } else {
        // centre = foot + t T + side r N, N = (nr, nz), T = (-nz, nr).
        const ProfileAnchor& A = anchors[0];
        const double Tr = -A.nz, Tz = A.nr;
        const double footR = A.d * A.nr, footZ = A.d * A.nz;
        double t = Tr * fR + Tz * fZ;
        r = side[0] * (A.nr * fR + A.nz * fZ - A.d);
        for (int it = 0; it < kFitIters; ++it) {
            double h11 = 0, h12 = 0, h22 = 0, g1 = 0, g2 = 0;
            const double ccr = footR + t * Tr + side[0] * r * A.nr;
            const double ccz = footZ + t * Tz + side[0] * r * A.nz;
            for (const ProfPt& p : P) {
                const double qx = p.rho - ccr, qz = p.z - ccz;
                const double dd = std::hypot(qx, qz);
                if (!(dd > gp::Resolution())) continue;
                const double res = dd - r;
                const double jt = -(qx * Tr + qz * Tz) / dd;
                const double jr = -side[0] * (qx * A.nr + qz * A.nz) / dd - 1.0;
                h11 += jt * jt;
                h12 += jt * jr;
                h22 += jr * jr;
                g1 += jt * res;
                g2 += jr * res;
            }
            const double hd = h11 * h22 - h12 * h12;
            if (!(std::fabs(hd) > gp::Resolution())) break;
            const double dt = -(g1 * h22 - g2 * h12) / hd;
            const double dr = -(h11 * g2 - h12 * g1) / hd;
            if (!std::isfinite(dt) || !std::isfinite(dr)) break;
            t += dt;
            r += dr;
            if (std::hypot(dt, dr) <= tau) break;
        }
        cr = footR + t * Tr + side[0] * r * A.nr;
        cz = footZ + t * Tz + side[0] * r * A.nz;
    }
    if (!std::isfinite(cr) || !std::isfinite(cz) || !std::isfinite(r)) {
        c.why = "fit-singular";
        return c;
    }
    c.Rmaj = cr;
    c.z0 = cz;
    c.Rmin = r;
    if (!(r > tau)) {
        c.why = "degenerate";
        return c;
    }
    if (!(cr > r)) {
        c.why = "spindle-torus";
        return c;
    }
    double sum2 = 0.0;
    for (const ProfPt& p : P) {
        const double d = circDevOf(p, cr, cz, r);
        c.maxDev = std::max(c.maxDev, d);
        sum2 += d * d;
    }
    c.rmsDev = std::sqrt(sum2 / (double)P.size());
    if (c.maxDev > tau) {
        c.why = "vertex-off-surface";
        return c;
    }
    // The normals span the anchor dihedral. The circle's normal at its tangent
    // point on anchor i is that anchor's normal; a station there is the round
    // meeting the anchor.
    double tR[2] = {0, 0}, tZ[2] = {0, 0};
    for (std::size_t i = 0; i < anchors.size(); ++i) {
        tR[i] = cr - side[i] * r * anchors[i].nr;
        tZ[i] = cz - side[i] * r * anchors[i].nz;
        bool meets = false;
        for (const ProfPt& p : P)
            if (std::hypot(p.rho - tR[i], p.z - tZ[i]) <= tau) {
                meets = true;
                break;
            }
        if (!meets) {
            c.why = "not-blend";
            return c;
        }
    }
    if (anchors.size() == 2) {
        // The fillet arc runs between the two tangent points on the side of the
        // anchors' corner c0: the bisector of the arc points at c0.
        const ProfileAnchor &A = anchors[0], &B = anchors[1];
        const double det = A.nr * B.nz - A.nz * B.nr;
        const double c0r = (A.d * B.nz - A.nz * B.d) / det;
        const double c0z = (A.nr * B.d - A.d * B.nr) / det;
        const double br = c0r - cr, bz = c0z - cz;
        if (!(std::hypot(br, bz) > gp::Resolution())) {
            c.why = "degenerate";
            return c;
        }
        const double half = angleBetween(tR[0] - cr, tZ[0] - cz, br, bz);
        const double slack = tau / r;
        for (const ProfPt& p : P)
            if (angleBetween(p.rho - cr, p.z - cz, br, bz) > half + slack) {
                c.why = "not-blend";
                return c;
            }
    }
    c.ok = true;
    c.why = "accept";
    return c;
}

// D-gid-2.6. Runs in runTorusStages after H1 (claimRelRevolve / reclaimRelational
// in Torus mode), before the sphere corners and A3.
//
// Donor: an accepted cylinder that certifies its surface (anchorCertified).
// Band: an edge-connected component of free triangles (no accepted region; an
// Unclaimed provisional that is not itself a plane anchor) touching the donor.
// Anchors: the donor (rho = R) and every certified plane normal to the donor's
// axis that the band touches (an accepted plane, or a provisional A3 would
// commit whose own residual is within tau): z = h. A provisional is an atom: it
// joins only whole, and only when every one of its vertices is within tau of
// the anchored torus (D-gid-3.5); an atom with a vertex off leaves whole and
// what remains is certified again as its own components (the H1 rule).
// The band ships as one closed TorusBlend when it is one annulus about the axis
// (two boundary loops, each winding once) and certifyAnchoredRound accepts it
// with every vertex of the band within tau.
bool claimAnchoredRounds(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                         SegmentWork& work) {
    try {
        const double tau = tauOf(mv);
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;
        const int nTri = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const bool dg = diagOn();

        std::vector<int> triAcc((std::size_t)nTri, -1);
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai)
            for (int t : work.accepted[(std::size_t)ai].tris)
                if (t >= 0 && t < nTri) triAcc[(std::size_t)t] = ai;
        std::vector<int> triProv((std::size_t)nTri, -1);
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi)
            for (int t : work.provisionals[(std::size_t)pi].tris)
                if (t >= 0 && t < nTri) triProv[(std::size_t)t] = pi;
        std::vector<char> provPlane(work.provisionals.size(), 0);
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
            const Provisional& pr = work.provisionals[(std::size_t)pi];
            provPlane[(std::size_t)pi] = !pr.tris.empty() && a3WouldCommit(mv, tol, pr, false) &&
                                         anchorCertified(tau, pr.maxVertexDev);
        }
        std::vector<char> withheld((std::size_t)nTri, 0);
        for (int t : p.torusWithheldTris)
            if (t >= 0 && t < nTri) withheld[(std::size_t)t] = 1;
        std::vector<char> taken((std::size_t)nTri, 0);
        auto candidate = [&](int t) -> bool {
            if (t < 0 || t >= nTri || taken[(std::size_t)t] || triAcc[(std::size_t)t] >= 0)
                return false;
            const int pi = triProv[(std::size_t)t];
            return pi >= 0 && work.provisionals[(std::size_t)pi].claim == ProvClaim::Unclaimed &&
                   !provPlane[(std::size_t)pi];
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
        std::vector<char> mark((std::size_t)nTri, 0);
        auto components = [&](const std::vector<int>& tris) {
            for (int t : tris) mark[(std::size_t)t] = 1;
            std::vector<std::vector<int>> comps;
            std::vector<char> seen((std::size_t)nTri, 0);
            for (int t0 : tris) {
                if (seen[(std::size_t)t0]) continue;
                std::vector<int> cc, st = {t0};
                seen[(std::size_t)t0] = 1;
                while (!st.empty()) {
                    const int u = st.back();
                    st.pop_back();
                    cc.push_back(u);
                    for (int s = 0; s < 3; ++s) {
                        const int o = otherTri(mv.triEdges[u][s], u);
                        if (o >= 0 && mark[(std::size_t)o] && !seen[(std::size_t)o]) {
                            seen[(std::size_t)o] = 1;
                            st.push_back(o);
                        }
                    }
                }
                std::sort(cc.begin(), cc.end());
                comps.push_back(std::move(cc));
            }
            for (int t : tris) mark[(std::size_t)t] = 0;
            std::sort(comps.begin(), comps.end(),
                      [](const std::vector<int>& a, const std::vector<int>& b) {
                          return a.front() < b.front();
                      });
            return comps;
        };

        int nDonors = 0, nBands = 0, nAccept = 0, nRefused = 0;
        const int nAcc0 = (int)work.accepted.size();
        for (int di = 0; di < nAcc0; ++di) {
            // A copy: the claims below append to work.accepted.
            const Region D = work.accepted[(std::size_t)di];
            if (D.type != SurfType::Cylinder || D.tris.empty()) continue;
            if (!(D.radius > tau)) continue;
            // D-fc-1.4: the donor certifies its surface. A relational band's
            // surface is certified on its tau-preimage I; its chord-bounded cover
            // does not un-certify it (D-fc-3.2.7). A B1 region is certified by
            // its own residual alone (D-fc-4.3). The sphere corners' rule.
            bool donorCert = anchorCertified(tau, D.maxVertexDev);
            RoundCert rc;
            if (!donorCert && D.origin != Origin::CylGrow) {
                RoundSpine sp;
                sp.circle = false;
                sp.o = D.ax.Location().XYZ();
                sp.a = D.ax.Direction().XYZ();
                sp.r = D.radius;
                std::vector<std::array<gp_XYZ, 3>> tp;
                tp.reserve(D.tris.size());
                for (int t : D.tris)
                    if (t >= 0 && t < nTri) tp.push_back({ptOf(mv, t, 0), ptOf(mv, t, 1), ptOf(mv, t, 2)});
                rc = certifyCylinderRound(sp, tp, tau);
                donorCert = rc.ok;
            }
            if (dg)
                diagPrintf("DIAG_REL_ANCHORED_DONOR acc=%d origin=%d R=%.9g nTri=%zu devTau=%.4f "
                           "certI=%d nI=%d nCover=%d nOther=%d stations=%d levels=%d "
                           "certDevTau=%.4f coverDevTau=%.4f sTau=%.4f certified=%d\n",
                           di, (int)D.origin, D.radius, D.tris.size(), D.maxVertexDev / tau,
                           rc.ok ? 1 : 0, rc.nI, rc.nCover, rc.nOther, rc.stations, rc.levels,
                           rc.certDev / tau, rc.coverDev / tau, rc.s / tau, donorCert ? 1 : 0);
            if (!donorCert) continue;
            RoundFrame fr;
            fr.o = D.ax.Location().XYZ();
            fr.a = D.ax.Direction().XYZ();
            fr.e1 = D.ax.XDirection().XYZ();
            fr.e1 -= fr.a * fr.e1.Dot(fr.a);
            if (!(fr.e1.Modulus() > gp::Resolution())) continue;
            fr.e1 /= fr.e1.Modulus();
            fr.e2 = fr.a.Crossed(fr.e1);
            // Free triangles across the donor's boundary.
            std::vector<int> free;
            {
                std::vector<int> seeds;
                for (int t : D.tris)
                    for (int s = 0; s < 3; ++s) {
                        const int o = otherTri(mv.triEdges[t][s], t);
                        if (candidate(o)) seeds.push_back(o);
                    }
                if (seeds.empty()) continue;
                std::sort(seeds.begin(), seeds.end());
                seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
                std::vector<char> seen((std::size_t)nTri, 0);
                std::vector<int> st = seeds;
                for (int t : seeds) seen[(std::size_t)t] = 1;
                while (!st.empty()) {
                    const int u = st.back();
                    st.pop_back();
                    free.push_back(u);
                    for (int s = 0; s < 3; ++s) {
                        const int o = otherTri(mv.triEdges[u][s], u);
                        if (candidate(o) && !seen[(std::size_t)o]) {
                            seen[(std::size_t)o] = 1;
                            st.push_back(o);
                        }
                    }
                }
                std::sort(free.begin(), free.end());
            }
            ++nDonors;
            std::vector<std::vector<int>> queue = components(free);
            // Per queue entry: the band an atom-off came from (its triangles are
            // the candidates of the tau-preimage below), and whether the entry IS
            // that preimage (it splits atoms by construction, D-gid-5.15).
            std::vector<std::vector<int>> parentOf(queue.size());
            std::vector<char> preimageOf(queue.size(), 0);
            auto enqueue = [&](std::vector<int> tris, const std::vector<int>& parent, bool pre) {
                queue.push_back(std::move(tris));
                parentOf.push_back(parent);
                preimageOf.push_back(pre ? 1 : 0);
            };
            for (std::size_t qi = 0; qi < queue.size(); ++qi) {
                std::vector<int> band = queue[qi];
                const std::vector<int> parent = parentOf[qi];
                const bool preimage = preimageOf[qi] != 0;
                ++nBands;
                std::vector<char> inBand((std::size_t)nTri, 0);
                for (int t : band) inBand[(std::size_t)t] = 1;
                std::vector<char> inParent((std::size_t)nTri, 0);
                for (int t : parent) inParent[(std::size_t)t] = 1;
                // Atoms join only whole.
                if (!preimage) {
                    std::vector<int> keep;
                    bool drop = false;
                    for (int t : band) {
                        const Provisional& pr =
                            work.provisionals[(std::size_t)triProv[(std::size_t)t]];
                        bool whole = true;
                        for (int u : pr.tris)
                            if (u < 0 || u >= nTri || !inBand[(std::size_t)u]) whole = false;
                        if (whole) keep.push_back(t);
                        else drop = true;
                    }
                    if (drop) {
                        for (auto& s : components(keep)) enqueue(std::move(s), parent, false);
                        continue;
                    }
                }
                // Profile of every band vertex about the donor axis.
                std::vector<int> gv;
                for (int t : band)
                    for (int k = 0; k < 3; ++k) gv.push_back(gidOf(mv, t, k));
                std::sort(gv.begin(), gv.end());
                gv.erase(std::unique(gv.begin(), gv.end()), gv.end());
                std::vector<ProfPt> P;
                P.reserve(gv.size());
                double rhoMax = 0.0;
                for (int g : gv) {
                    P.push_back(profOf(fr, mv.pts[g], g));
                    rhoMax = std::max(rhoMax, P.back().rho);
                }
                // Anchors: the donor and the certified planes normal to its axis
                // that the band touches.
                std::vector<ProfileAnchor> anchors;
                bool touchDonor = false;
                std::vector<int> rimOwner;  // per boundary edge: acc idx, -(prov + 2), or -1
                std::vector<int> rimTri;    // per boundary edge: the triangle across it, or -1
                std::vector<int> anchorOwner = {di};
                {
                    std::vector<int> accN, provN;
                    for (int t : band)
                        for (int s = 0; s < 3; ++s) {
                            const int o = otherTri(mv.triEdges[t][s], t);
                            if (o >= 0 && inBand[(std::size_t)o]) continue;
                            rimTri.push_back(o);
                            if (o < 0) {
                                rimOwner.push_back(-1);
                            } else if (triAcc[(std::size_t)o] >= 0) {
                                accN.push_back(triAcc[(std::size_t)o]);
                                rimOwner.push_back(triAcc[(std::size_t)o]);
                            } else if (triProv[(std::size_t)o] >= 0 &&
                                       provPlane[(std::size_t)triProv[(std::size_t)o]]) {
                                provN.push_back(triProv[(std::size_t)o]);
                                rimOwner.push_back(-(triProv[(std::size_t)o] + 2));
                            } else {
                                rimOwner.push_back(-1);
                            }
                        }
                    std::sort(accN.begin(), accN.end());
                    accN.erase(std::unique(accN.begin(), accN.end()), accN.end());
                    std::sort(provN.begin(), provN.end());
                    provN.erase(std::unique(provN.begin(), provN.end()), provN.end());
                    // L: the ring's diameter, over which a tilted plane moves.
                    const double L = std::max(2.0 * rhoMax, tau);
                    auto addPlane = [&](const gp_XYZ& n, const gp_XYZ& o, double devA) {
                        const double tolA = anchorTol(tau, devA);
                        if (!(n.Crossed(fr.a).Modulus() <= tanAngTol(tolA, tau, L))) return false;
                        const ProfileAnchor pa{0.0, 1.0, (o - fr.o).Dot(fr.a)};
                        for (const ProfileAnchor& q : anchors)
                            if (q.nr == pa.nr && q.nz == pa.nz &&
                                std::fabs(q.d - pa.d) <= tanDistTol(tolA, tau))
                                return true;
                        anchors.push_back(pa);
                        return true;
                    };
                    for (int ai : accN) {
                        if (ai == di) {
                            touchDonor = true;
                            continue;
                        }
                        const Region& A = work.accepted[(std::size_t)ai];
                        if (A.type != SurfType::Plane || A.tris.empty() ||
                            !anchorCertified(tau, A.maxVertexDev))
                            continue;
                        if (addPlane(A.ax.Direction().XYZ(), A.ax.Location().XYZ(), A.maxVertexDev))
                            anchorOwner.push_back(ai);
                    }
                    for (int pi : provN) {
                        const Provisional& pr = work.provisionals[(std::size_t)pi];
                        if (addPlane(pr.plane.Direction().XYZ(), pr.plane.Location().XYZ(),
                                     pr.maxVertexDev))
                            anchorOwner.push_back(-(pi + 2));
                    }
                    anchors.insert(anchors.begin(), ProfileAnchor{1.0, 0.0, D.radius});
                }
                AnchoredRoundCert cert;
                const char* why = "accept";
                int loops = 0, wind0 = 0, wind1 = 0;
                double dVol = 0.0;
                auto emitDiag = [&](const char* w) {
                    if (!dg) return;
                    diagPrintf("DIAG_REL_BAND stage=anchored mode=round seed=%d nTri=%zu nV=%zu "
                               "stations=%d cols=%d params=%d anchors=%d donor=%d donorR=%.9g "
                               "radii=%.9g/%.9g z0=%.9g maxDev=%.6g maxDevTau=%.4f "
                               "loops=%d winding=%d,%d dVol=%.9g why=%s\n",
                               band.front(), band.size(), gv.size(), cert.stations, cert.cols,
                               cert.params, (int)anchors.size(), di, D.radius, cert.Rmaj,
                               cert.Rmin, cert.z0, cert.maxDev, cert.maxDev / tau, loops, wind0,
                               wind1, dVol, w);
                };
                if (!touchDonor) {
                    ++nRefused;
                    emitDiag("no-donor-contact");
                    continue;
                }
                cert = certifyAnchoredRound(P, anchors, tau);
                if (!preimage && !cert.ok && cert.determinate && cert.Rmin > 0.0 &&
                    std::strcmp(cert.why, "vertex-off-surface") == 0) {
                    // An atom with a vertex off the anchored circle leaves whole;
                    // what remains is certified again as its own components.
                    std::vector<int> keep;
                    bool dropped = false;
                    std::map<int, bool> off;
                    for (int t : band) {
                        const int pi = triProv[(std::size_t)t];
                        auto it = off.find(pi);
                        if (it == off.end()) {
                            bool o = false;
                            for (int u : work.provisionals[(std::size_t)pi].tris)
                                for (int k = 0; k < 3 && !o; ++k)
                                    if (circDevOf(profOf(fr, ptOf(mv, u, k), -1), cert.Rmaj,
                                                  cert.z0, cert.Rmin) > tau)
                                        o = true;
                            it = off.emplace(pi, o).first;
                        }
                        if (it->second) dropped = true;
                        else keep.push_back(t);
                    }
                    if (dropped && !keep.empty()) {
                        emitDiag("atom-off");
                        const std::vector<int>& from = parent.empty() ? band : parent;
                        for (auto& s : components(keep)) enqueue(std::move(s), from, false);
                        continue;
                    }
                }
                if (!cert.ok) {
                    ++nRefused;
                    emitDiag(cert.why);
                    continue;
                }
                // D-gid-5.15: the face claims the tau-preimage of the certified
                // circle -- every triangle of the band the atom-off came from whose
                // three vertices lie within tau of it, edge-connected to this band.
                // An atom whose other triangles have a vertex off is split: its
                // preimage triangles join, the rest stay facets (its atom is taken
                // as pending below, so A3 does not commit them) and the face is
                // bounded there by their mesh polyline. No vertex past tau joins.
                if (!preimage && !parent.empty()) {
                    auto onCircle = [&](int t) {
                        for (int k = 0; k < 3; ++k)
                            if (circDevOf(profOf(fr, ptOf(mv, t, k), -1), cert.Rmaj, cert.z0,
                                          cert.Rmin) > tau)
                                return false;
                        return true;
                    };
                    std::vector<char> inPre((std::size_t)nTri, 0);
                    for (int t : parent)
                        if (inBand[(std::size_t)t] || (candidate(t) && onCircle(t)))
                            inPre[(std::size_t)t] = 1;
                    std::vector<int> grown, st = band;
                    std::vector<char> seen((std::size_t)nTri, 0);
                    for (int t : band) seen[(std::size_t)t] = 1;
                    while (!st.empty()) {
                        const int u = st.back();
                        st.pop_back();
                        grown.push_back(u);
                        for (int s = 0; s < 3; ++s) {
                            const int o = otherTri(mv.triEdges[u][s], u);
                            if (o >= 0 && inPre[(std::size_t)o] && !seen[(std::size_t)o]) {
                                seen[(std::size_t)o] = 1;
                                st.push_back(o);
                            }
                        }
                    }
                    std::sort(grown.begin(), grown.end());
                    if (dg)
                        diagPrintf("DIAG_REL_ANCHORED_PREIMAGE seed=%d nParent=%zu nBand=%zu "
                                   "nPreimage=%zu\n",
                                   band.front(), parent.size(), band.size(), grown.size());
                    if (grown.size() > band.size()) {
                        enqueue(std::move(grown), parent, true);
                        continue;
                    }
                }
                // D-fc-1.11 (ii)(a): the band without any one atom still meets the
                // anchored station count.
                {
                    std::vector<int> atoms;
                    for (int t : band) atoms.push_back(triProv[(std::size_t)t]);
                    std::sort(atoms.begin(), atoms.end());
                    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
                    bool det = true;
                    if (atoms.size() > 1)
                        for (int a : atoms) {
                            std::vector<int> g2;
                            for (int t : band)
                                if (triProv[(std::size_t)t] != a)
                                    for (int k = 0; k < 3; ++k) g2.push_back(gidOf(mv, t, k));
                            std::sort(g2.begin(), g2.end());
                            g2.erase(std::unique(g2.begin(), g2.end()), g2.end());
                            std::vector<ProfPt> P2;
                            P2.reserve(g2.size());
                            for (int g : g2) P2.push_back(profOf(fr, mv.pts[g], g));
                            if (countProfPts(P2, tau) < minStations(cert.params)) {
                                det = false;
                                break;
                            }
                        }
                    else
                        det = false;
                    if (!det) {
                        ++nRefused;
                        emitDiag("underdetermined");
                        continue;
                    }
                }
                // One annulus about the axis: two boundary loops, each winding once
                // (D-fc-5.7 O-1).
                {
                    std::map<int, std::vector<int>> next;  // boundary edge a -> b
                    for (int t : band)
                        for (int s = 0; s < 3; ++s) {
                            const int o = otherTri(mv.triEdges[t][s], t);
                            if (o >= 0 && inBand[(std::size_t)o]) continue;
                            // Edge s of triangle t: find its two corners in the
                            // triangle's own order.
                            const int e = mv.triEdges[t][s];
                            if (e < 0 || e >= nEdge || !mv.compEdges) continue;
                            const std::pair<int, int>& ce = mv.compEdges[e];
                            const int ga = mv.compVtx ? mv.compVtx[ce.first] : ce.first;
                            const int gb = mv.compVtx ? mv.compVtx[ce.second] : ce.second;
                            int ka = -1, kb = -1;
                            for (int k = 0; k < 3; ++k) {
                                if (gidOf(mv, t, k) == ga) ka = k;
                                if (gidOf(mv, t, k) == gb) kb = k;
                            }
                            if (ka < 0 || kb < 0) continue;
                            if ((ka + 1) % 3 == kb) next[ga].push_back(gb);
                            else next[gb].push_back(ga);
                        }
                    bool pinch = false;
                    for (const auto& kv : next)
                        if (kv.second.size() != 1) pinch = true;
                    std::vector<int> winds;
                    if (!pinch) {
                        std::map<int, char> done;
                        for (const auto& kv : next) {
                            if (done[kv.first]) continue;
                            std::vector<double> az;
                            int v = kv.first;
                            bool closedLoop = false;
                            for (std::size_t guard = 0; guard <= next.size(); ++guard) {
                                done[v] = 1;
                                const gp_XYZ d = mv.pts[v] - fr.o;
                                az.push_back(std::atan2(d.Dot(fr.e2), d.Dot(fr.e1)));
                                auto it = next.find(v);
                                if (it == next.end()) break;
                                v = it->second.front();
                                if (v == kv.first) {
                                    closedLoop = true;
                                    break;
                                }
                            }
                            winds.push_back(closedLoop ? windingNumber(az) : 0);
                        }
                    }
                    loops = (int)winds.size();
                    wind0 = loops > 0 ? winds[0] : 0;
                    wind1 = loops > 1 ? winds[1] : 0;
                    if (pinch || !annulusAboutAxis(1, winds)) {
                        ++nRefused;
                        emitDiag("not-annulus");
                        continue;
                    }
                }
                // D-140-2 / D-fc-4.6 F9: a rim edge whose partner is not one of the
                // band's anchors (facets, an uncertified region, an atom that left)
                // cannot carry the anchor's circle, so the ring cannot ship; its
                // triangles keep their state.
                // D-gid-5.15: the one exception is a preimage's own remainder -- a
                // triangle of the parent band whose atom the preimage split. It
                // stays a facet and the face is bounded there by its mesh polyline.
                std::vector<char> splitAtom(work.provisionals.size(), 0);
                if (preimage)
                    for (int t : band) splitAtom[(std::size_t)triProv[(std::size_t)t]] = 1;
                int foreign = 0, facetRim = 0;
                for (std::size_t ei = 0; ei < rimOwner.size(); ++ei) {
                    const int ow = rimOwner[ei];
                    if (std::find(anchorOwner.begin(), anchorOwner.end(), ow) != anchorOwner.end())
                        continue;
                    const int o = ei < rimTri.size() ? rimTri[ei] : -1;
                    if (preimage && o >= 0 && inParent[(std::size_t)o] && triAcc[(std::size_t)o] < 0 &&
                        triProv[(std::size_t)o] >= 0 && splitAtom[(std::size_t)triProv[(std::size_t)o]])
                        ++facetRim;
                    else
                        ++foreign;
                }
                if (dg && preimage)
                    diagPrintf("DIAG_REL_ANCHORED_RIM seed=%d rimEdges=%zu foreign=%d facetRim=%d "
                               "anchors=%zu\n",
                               band.front(), rimOwner.size(), foreign, facetRim, anchorOwner.size());
                if (foreign > 0) {
                    ++nRefused;
                    if (dg)
                        diagPrintf("DIAG_REL_ANCHORED_RIM seed=%d rimEdges=%zu foreign=%d "
                                   "anchors=%zu\n",
                                   band.front(), rimOwner.size(), foreign, anchorOwner.size());
                    emitDiag("foreign-rim");
                    continue;
                }
                // Orientation: every triangle faces the same way about the tube.
                double sigma = 0.0;
                int nPos = 0, nNeg = 0;
                for (int t : band) {
                    gp_XYZ n;
                    const double at = triAreaNormal(mv, t, n);
                    const gp_XYZ cen = (ptOf(mv, t, 0) + ptOf(mv, t, 1) + ptOf(mv, t, 2)) / 3.0;
                    const gp_XYZ foot = torusFoot(fr, cert.Rmaj, cert.z0, cert.Rmin, cen);
                    // Tube-centre direction at the centroid's azimuth.
                    const gp_XYZ d = cen - fr.o;
                    gp_XYZ rad = d - fr.a * d.Dot(fr.a);
                    const double rm = rad.Modulus();
                    if (!(rm > gp::Resolution())) continue;
                    rad /= rm;
                    gp_XYZ mhat = foot - (fr.o + rad * cert.Rmaj + fr.a * cert.z0);
                    const double mm = mhat.Modulus();
                    if (!(mm > gp::Resolution())) continue;
                    const double s = n.Dot(mhat / mm);
                    sigma += at * s;
                    if (s > 0.0) ++nPos;
                    else if (s < 0.0) ++nNeg;
                }
                if (nPos > 0 && nNeg > 0) {
                    ++nRefused;
                    emitDiag("not-blend");
                    continue;
                }
                bool held = false;
                for (int t : band)
                    if (withheld[(std::size_t)t]) held = true;
                if (held) {
                    // U-3 per torus (D-fc-4.11 (a)): the builder could not ship this
                    // band; its triangles keep their pre-T state.
                    ++nRefused;
                    emitDiag("withheld");
                    continue;
                }
                for (int t : band)
                    dVol += torusTriVolume(fr, cert.Rmaj, cert.z0, cert.Rmin, tau, ptOf(mv, t, 0),
                                           ptOf(mv, t, 1), ptOf(mv, t, 2));
                std::vector<double> vs;
                vs.reserve(P.size());
                for (const ProfPt& pp : P) vs.push_back(std::atan2(pp.z - cert.z0, pp.rho - cert.Rmaj));
                double vLo = 0, vHi = 0;
                {
                    // The profile arc's extent: the complement of its largest gap.
                    std::vector<double> a = vs;
                    for (double& x : a)
                        if (x < 0.0) x += kTwoPi;
                    std::sort(a.begin(), a.end());
                    std::size_t ig = a.size() - 1;
                    double maxGap = (a.front() + kTwoPi) - a.back();
                    for (std::size_t i = 0; i + 1 < a.size(); ++i)
                        if (a[i + 1] - a[i] > maxGap) {
                            maxGap = a[i + 1] - a[i];
                            ig = i;
                        }
                    vLo = a[(ig + 1) % a.size()];
                    vHi = a[ig];
                    while (vLo > kPi) vLo -= kTwoPi;
                    while (vLo <= -kPi) vLo += kTwoPi;
                    while (vHi > kPi) vHi -= kTwoPi;
                    while (vHi <= -kPi) vHi += kTwoPi;
                    if (vHi < vLo) vHi += kTwoPi;
                }
                const int nSides = torusSidesOf(cert.cols, kTwoPi, true);
                const gp_XYZ centre = fr.o + fr.a * cert.z0;
                Region r;
                r.id = 0;
                r.type = SurfType::Torus;
                r.origin = Origin::TorusBlend;
                r.ax = gp_Ax3(gp_Pnt(centre.X(), centre.Y(), centre.Z()),
                              gp_Dir(fr.a.X(), fr.a.Y(), fr.a.Z()),
                              gp_Dir(fr.e1.X(), fr.e1.Y(), fr.e1.Z()));
                r.radius = cert.Rmaj;
                r.radius2 = cert.Rmin;
                r.uMin = 0.0;
                r.uMax = kTwoPi;
                r.vMin = vLo;
                r.vMax = vHi;
                r.closed360 = true;
                r.outwardNormal = sigma >= 0.0;
                r.tris = band;
                r.maxVertexDev = cert.maxDev;
                r.rmsVertexDev = cert.rmsDev;
                r.chordSagitta = chordSagitta(cert.Rmaj, nSides);
                r.nSides = nSides;
                r.dVolPredicted = dVol;
                r.reject = Reject::None;
                r.builtAs = BuiltAs::NotBuilt;
                r.filletNbrA = di;
                r.filletNbrB = -1;
                r.publishRimFirst = true;
                r.meshPolylineCap = facetRim > 0;
                // Every atom is whole inside the band: each provisional is taken
                // as a pending torus (A3 skips it; U-3 commits or reverts).
                TorusPendingReclaim pend;
                {
                    std::vector<int> atoms;
                    for (int t : band) atoms.push_back(triProv[(std::size_t)t]);
                    std::sort(atoms.begin(), atoms.end());
                    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
                    for (int a : atoms) {
                        work.provisionals[(std::size_t)a].claim = ProvClaim::PendingTorus;
                        pend.provIdx.push_back(a);
                    }
                }
                r.torusPending = (int)work.torusPending.size();
                work.torusPending.push_back(std::move(pend));
                work.torusAdmitted = true;
                for (int t : band) {
                    taken[(std::size_t)t] = 1;
                    triAcc[(std::size_t)t] = (int)work.accepted.size();
                }
                ++nAccept;
                emitDiag("accept");
                work.accepted.push_back(std::move(r));
            }
        }
        if (dg)
            diagPrintf("DIAG_REL_SUM stage=anchored donors=%d bands=%d accept=%d refused=%d\n",
                       nDonors, nBands, nAccept, nRefused);
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
