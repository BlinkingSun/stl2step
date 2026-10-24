// fc/fx-t2-walls (D-fc-6.8 RQ-1 and RQ-3, D-fc-6.11). Predicates of the coaxial Line
// stage's closed walls and tangent-anchored partial walls. Inline; no tolerance of
// its own: every bound is tau = 2q, the chord bound s of the certified set
// (D-fc-3.2.7), tanDistTol / tanAngTol (D-fc-1.6.7) or B1's seed floor thetaCylLo
// (read only).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_COAX_CLOSED_HPP
#define STL2STEP_REFIT_COAX_CLOSED_HPP

#include "refit_axis_seam.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {

// |rho(p) - R| about the axis through o along the unit vector a.
inline double coaxRadialDev(const gp_XYZ& p, const gp_XYZ& o, const gp_XYZ& a, double R) {
    const gp_XYZ d = p - o;
    return std::fabs((d - a * d.Dot(a)).Modulus() - R);
}

// D-fc-3.2.7. The chord bound s of the certified set I: the largest radial
// deviation of the midpoint of any edge of an I triangle. It is a deviation I
// itself exhibits (its rim chords), so a cover triangle within s deviates no more
// than the certified triangles' own chords already do.
inline double coaxChordBound(const std::vector<std::array<gp_XYZ, 3>>& I, const gp_XYZ& o,
                             const gp_XYZ& a, double R) {
    double s = 0;
    for (const std::array<gp_XYZ, 3>& t : I)
        for (int k = 0; k < 3; ++k)
            s = std::max(s, coaxRadialDev(0.5 * (t[(std::size_t)k] + t[(std::size_t)((k + 1) % 3)]),
                                          o, a, R));
    return s;
}

// D-fc-3.2.7 Cover inside the coaxial gather (D-fc-6.8 RQ-1): a triangle with at
// least one vertex within tau of the certified cylinder, every vertex within s,
// two axial levels (a section ear is never a member, and cover does not put it
// back) and every vertex inside I's axial extent [zLo - tau, zHi + tau] (cover
// fills, never extends). tau <= 0 admits nothing.
inline bool coaxCoverAdmits(const std::array<gp_XYZ, 3>& v, const gp_XYZ& o, const gp_XYZ& a,
                            double R, double tau, double s, double zLo, double zHi) {
    if (!(tau > 0.0)) return false;
    bool one = false;
    double lo = 1e300, hi = -1e300;
    for (const gp_XYZ& p : v) {
        const double d = coaxRadialDev(p, o, a, R);
        if (!(d <= s)) return false;
        if (d <= tau) one = true;
        const double z = (p - o).Dot(a);
        if (z < zLo - tau || z > zHi + tau) return false;
        lo = std::min(lo, z);
        hi = std::max(hi, z);
    }
    return one && hi - lo > tau;
}

// D-fc-6.11 (a). A fill of a closed coaxial wall (D-fc-6.8 RQ-1): a triangle ENCLOSED
// by I + Cover -- each of its edges e[k] shared with an I + Cover triangle
// (shared[k]) or lying on the band's own rims (both ends within tau of one rim
// level, zLo or zHi) -- whose every vertex is within s. It spans two axial levels
// inside I's extent (a section ear is never a member; cover does not put it back).
// A fill closes a hole in the closure; it never extends the band's outline.
inline bool coaxFillAdmits(const std::array<gp_XYZ, 3>& v,
                           const std::array<std::array<gp_XYZ, 2>, 3>& e,
                           const std::array<bool, 3>& shared, const gp_XYZ& o, const gp_XYZ& a,
                           double R, double tau, double s, double zLo, double zHi) {
    if (!(tau > 0.0)) return false;
    double lo = 1e300, hi = -1e300;
    for (const gp_XYZ& p : v) {
        if (!(coaxRadialDev(p, o, a, R) <= s)) return false;
        const double z = (p - o).Dot(a);
        if (z < zLo - tau || z > zHi + tau) return false;
        lo = std::min(lo, z);
        hi = std::max(hi, z);
    }
    if (!(hi - lo > tau)) return false;
    for (int k = 0; k < 3; ++k) {
        if (shared[(std::size_t)k]) continue;
        const double z0 = (e[(std::size_t)k][0] - o).Dot(a), z1 = (e[(std::size_t)k][1] - o).Dot(a);
        const bool onLo = std::fabs(z0 - zLo) <= tau && std::fabs(z1 - zLo) <= tau;
        const bool onHi = std::fabs(z0 - zHi) <= tau && std::fabs(z1 - zHi) <= tau;
        if (!onLo && !onHi) return false;
    }
    return true;
}

// Partial-wall cover stays inside the certified set's own angular extent: the
// complement of the largest gap [gapLo, gapHi) of I's vertex azimuths (windowGap,
// the convention of Region::uMin/uMax), with tau/R slack.
inline bool coaxInArc(double u, double gapLo, double gapHi, double slack) {
    const double twoPi = 2.0 * std::acos(-1.0);
    double w = u;
    while (w < gapLo) w += twoPi;
    while (w >= gapLo + twoPi) w -= twoPi;
    return !(w > gapLo + slack && w < gapHi - slack);
}

// The dihedral fold between two triangle normals (radians).
inline double coaxFold(const gp_XYZ& n1, const gp_XYZ& n2) {
    const double m = n1.Modulus() * n2.Modulus();
    if (!(m > 0.0)) return 0.0;
    return std::acos(std::max(-1.0, std::min(1.0, n1.Dot(n2) / m)));
}

// D-fc-5.11 F1 scope (the ring path, refit_relational_reclaim.cpp claimRingCylinders):
// the lift names walls B1 cannot seed, every fold of the certified set under B1's
// seed floor thetaCylLo. A wall that folds at or above it is B1 / NgonWall /
// tess-law's to decide: its vertices alone cannot tell a prism from a cylinder.
inline bool coaxFoldInScope(double maxFold, double thetaCylLo) { return maxFold < thetaCylLo; }

// D-fc-7.3 RQ-HUB-1. A closed wall that folds at or above B1's seed floor is still
// the lift's when B1 / NgonWall / tess-law left every one of its triangles free
// (everyTriFree) AND it carries the D-fc-6.8 RQ-3 (i) certified tangent parallel
// cylinder anchor on a shared contact vertex (anchored): the stages the fold scope
// hands such a wall to did not claim it, and the anchor is the relation its
// vertices alone cannot give.
inline bool coaxClosedLiftInScope(double maxFold, double thetaCylLo, bool everyTriFree,
                                  bool anchored) {
    return coaxFoldInScope(maxFold, thetaCylLo) || (everyTriFree && anchored);
}

// D-fc-7.3 RQ-HUB-1 corollary. A band whose closure closes (D2) is a closed wall or
// nothing, whichever attempt of the stage measured that closure: a set inside a
// D2-closed closure (each sorted) is a passing subset of it, never a partial wall.
// `set` is sorted.
inline bool coaxInsideClosure(const std::vector<int>& set,
                              const std::vector<std::vector<int>>& closures) {
    if (set.empty()) return false;
    for (const std::vector<int>& c : closures)
        if (std::includes(c.begin(), c.end(), set.begin(), set.end())) return true;
    return false;
}

// D-fc-6.8 RQ-1. The lift is the F1 route's: a closed band whose tau-preimage I is
// already ONE edge-connected annulus is the F1 ring path's to claim, unchanged. The
// coaxial stage commits only the closure I alone cannot make (chord triangles split
// it), with the D-fc-3.2.7 cover.
inline bool coaxLiftRoutes(bool preimageIsAnnulus) { return !preimageIsAnnulus; }

// D-fc-6.8 RQ-1. A closed wall: closed by the stage's own D2 rule (closedOf), ONE
// edge component, its boundary the two rims, each winding once about the family
// axis (annulusAboutAxis, D-fc-5.7 O-1). A band with windows (more boundary loops)
// is coaxClosedWallWithWindows' to decide.
inline bool coaxClosedWall(bool d2Closed, int edgeComponents, const std::vector<int>& windings) {
    return d2Closed && annulusAboutAxis(edgeComponents, windings);
}

// D-fc-7.3 RQ-HUB-2 (the D-fc-6.8 RQ-1 window clause). A closed wall may carry
// windows: D2-closed, ONE edge component, its loops the two rims winding once in
// opposite senses plus loops winding zero times (closedWallWithWindows) that one
// seam generator clears (seamClearOfWindows over each window's loop azimuths
// windowU) -- the D-fc-4.4 ring closure's test. The windows are the face's inner
// wires. A band with no window is coaxClosedWall's (seamU = 0).
inline bool coaxClosedWallWithWindows(bool d2Closed, int edgeComponents,
                                      const std::vector<int>& windings,
                                      const std::vector<std::vector<double>>& windowU,
                                      double& seamU) {
    seamU = 0.0;
    if (coaxClosedWall(d2Closed, edgeComponents, windings)) return true;
    return d2Closed && edgeComponents == 1 && closedWallWithWindows(windings) &&
           seamClearOfWindows(windowU, seamU);
}

// D-fc-6.8 RQ-3 (i). A partial coaxial wall (family axis oF/aF, radius R) is
// anchored by a certified cylinder (oA/aA, r) tangent to it with a PARALLEL axis:
// |aF x aA| <= tanAngTol(tolA, tau, L) and |dist(axes) - (R + r)| or
// |dist(axes) - |R - r|| <= tanDistTol(tolA, tau). The anchor's axis is distinct
// from the family's (a coaxial cylinder is a family member, never a tangent
// anchor). Contact on a shared vertex is the caller's test.
struct CoaxTangent {
    bool ok = false;
    bool internal = false;  // |R - r|: the anchor lies inside the wall's circle
    double cross = 0;       // |aF x aA|
    double dist = 0;        // distance between the parallel axes
    double tanDist = 0;     // the smaller of |dist - (R + r)|, |dist - |R - r||
};
inline CoaxTangent coaxTangentCylinder(const gp_XYZ& oF, const gp_XYZ& aFIn, double R,
                                       const gp_XYZ& oA, const gp_XYZ& aAIn, double r,
                                       double tolA, double tau, double L) {
    CoaxTangent c;
    if (!(tau > 0.0) || !(R > tau) || !(r > tau) || !(L > 0.0)) return c;
    const double mf = aFIn.Modulus(), ma = aAIn.Modulus();
    if (!(mf > 0.0) || !(ma > 0.0)) return c;
    const gp_XYZ aF = aFIn / mf, aA = aAIn / ma;
    c.cross = aF.Crossed(aA).Modulus();
    const gp_XYZ d = oA - oF;
    c.dist = (d - aF * d.Dot(aF)).Modulus();
    const double ext = std::fabs(c.dist - (R + r));
    const double inn = std::fabs(c.dist - std::fabs(R - r));
    c.internal = inn < ext;
    c.tanDist = std::min(ext, inn);
    c.ok = c.cross <= tanAngTol(tolA, tau, L) && c.tanDist <= tanDistTol(tolA, tau) &&
           c.dist > tanDistTol(tolA, tau);
    return c;
}

// D-fc-6.8 RQ-3. The contact: a vertex the wall and the anchor share lies within
// tanDistTol of both surfaces.
inline bool coaxContact(const gp_XYZ& p, const gp_XYZ& oF, const gp_XYZ& aF, double R,
                        const gp_XYZ& oA, const gp_XYZ& aA, double r, double tolA, double tau) {
    const double m = aA.Modulus();
    if (!(m > 0.0)) return false;
    return coaxRadialDev(p, oF, aF, R) <= tanDistTol(tolA, tau) &&
           coaxRadialDev(p, oA, aA / m, r) <= tanDistTol(tolA, tau);
}

// D-fc-7.6. `round` was fitted exactly tangent (Region::tangentExact) to the very
// surface `wall` carries: the same axis point, direction and radius as recorded at
// that fit (an identity, not a tolerance: any later change to the wall breaks it).
inline bool roundTangentToWall(const Region& round, const Region& wall) {
    if (!round.tangentExact || wall.type != SurfType::Cylinder) return false;
    const gp_XYZ o = wall.ax.Location().XYZ(), a = wall.ax.Direction().XYZ();
    return o.X() == round.tangentWallOrigin.X() && o.Y() == round.tangentWallOrigin.Y() &&
           o.Z() == round.tangentWallOrigin.Z() && a.X() == round.tangentWallAxis.X() &&
           a.Y() == round.tangentWallAxis.Y() && a.Z() == round.tangentWallAxis.Z() &&
           wall.radius == round.tangentWallRadius;
}

// D-fc-7.6 (T). A round tangent to a certified coaxial wall (axis through o along
// the unit a, radius R) fitted in the D-fc-1.4 anchored mode: its axis is held
// parallel to the wall's and the tangency exact, |o_wall o_round| = R + r
// (external) or |R - r| (internal), so its centre in the wall's cross-section is
// (R + s r)(cos phi e1 + sin phi e2), s = +1 / -1. phi is free; r is free when
// fitRadius, else held at the free fit's r0 (fewest parameters first, D-fc-1.6:
// the caller tries phi alone before phi and r). Gauss-Newton from the free fit
// (centre c0 on the round's axis, radius r0) over the round's vertices pts, at
// most maxIter steps. maxDev is the largest radial deviation of pts from the
// anchored cylinder; the caller certifies it at tau.
struct AnchoredTangentFit {
    bool ok = false;
    gp_XYZ centre{0, 0, 0};  // on the wall's cross-section through o (axis point of the round)
    double radius = 0;
    double maxDev = 0, rmsDev = 0;
};
inline AnchoredTangentFit anchoredTangentFit(const std::vector<gp_XYZ>& pts, const gp_XYZ& o,
                                             const gp_XYZ& aIn, double R, bool internal,
                                             const gp_XYZ& c0, double r0, int maxIter,
                                             bool fitRadius = true) {
    AnchoredTangentFit f;
    const double am = aIn.Modulus();
    if (pts.size() < 3 || !(am > 0.0) || !(R > 0.0) || !(r0 > 0.0) || maxIter < 1) return f;
    const gp_XYZ a = aIn / am;
    gp_XYZ e1 = std::fabs(a.X()) < 0.9 ? gp_XYZ(1, 0, 0).Crossed(a) : gp_XYZ(0, 1, 0).Crossed(a);
    e1 /= e1.Modulus();
    const gp_XYZ e2 = a.Crossed(e1);
    const double s = internal ? -1.0 : 1.0;
    gp_XYZ d0 = c0 - o;
    d0 -= a * d0.Dot(a);
    double phi = std::atan2(d0.Dot(e2), d0.Dot(e1));
    double r = r0;
    if (R + s * r < 0.0) phi += std::acos(-1.0);  // the centre lies past the axis
    std::vector<gp_XYZ> q;
    q.reserve(pts.size());
    for (const gp_XYZ& p : pts) {
        gp_XYZ v = p - o;
        q.push_back(v - a * v.Dot(a));
    }
    for (int it = 0; it < maxIter; ++it) {
        const gp_XYZ u(e1 * std::cos(phi) + e2 * std::sin(phi));
        const gp_XYZ t(e2 * std::cos(phi) - e1 * std::sin(phi));
        const gp_XYZ c = u * (R + s * r);
        double h00 = 0, h01 = 0, h11 = 0, g0 = 0, g1 = 0;
        for (const gp_XYZ& p : q) {
            const gp_XYZ v = p - c;
            const double n = v.Modulus();
            if (!(n > 0.0)) return f;
            const gp_XYZ w = v / n;
            const double res = n - r;
            const double j0 = -w.Dot(t * (R + s * r));
            const double j1 = -w.Dot(u * s) - 1.0;
            h00 += j0 * j0;
            h01 += j0 * j1;
            h11 += j1 * j1;
            g0 += j0 * res;
            g1 += j1 * res;
        }
        double dphi = 0, dr = 0;
        if (fitRadius) {
            const double det = h00 * h11 - h01 * h01;
            if (!(std::fabs(det) > 0.0)) break;
            dphi = -(h11 * g0 - h01 * g1) / det;
            dr = -(h00 * g1 - h01 * g0) / det;
        } else {
            if (!(h00 > 0.0)) break;
            dphi = -g0 / h00;
        }
        if (!std::isfinite(dphi) || !std::isfinite(dr)) return f;
        phi += dphi;
        r += dr;
        if (dphi == 0.0 && dr == 0.0) break;
    }
    if (!(r > 0.0)) return f;
    const gp_XYZ c = (e1 * std::cos(phi) + e2 * std::sin(phi)) * (R + s * r);
    double mx = 0, s2 = 0;
    for (const gp_XYZ& p : q) {
        const double d = std::fabs((p - c).Modulus() - r);
        mx = std::max(mx, d);
        s2 += d * d;
    }
    f.ok = true;
    f.centre = o + c;
    f.radius = r;
    f.maxDev = mx;
    f.rmsDev = std::sqrt(s2 / (double)q.size());
    return f;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_COAX_CLOSED_HPP
