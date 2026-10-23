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

// D-fc-6.8 RQ-1. The lift is the F1 route's: a closed band whose tau-preimage I is
// already ONE edge-connected annulus is the F1 ring path's to claim, unchanged. The
// coaxial stage commits only the closure I alone cannot make (chord triangles split
// it), with the D-fc-3.2.7 cover.
inline bool coaxLiftRoutes(bool preimageIsAnnulus) { return !preimageIsAnnulus; }

// D-fc-6.8 RQ-1. A closed wall: closed by the stage's own D2 rule (closedOf), ONE
// edge component, its boundary the two rims, each winding once about the family
// axis (annulusAboutAxis, D-fc-5.7 O-1). A band with windows (more boundary loops)
// is not committed here (closedWallWithWindows needs oriented loop windings).
inline bool coaxClosedWall(bool d2Closed, int edgeComponents, const std::vector<int>& windings) {
    return d2Closed && annulusAboutAxis(edgeComponents, windings);
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

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_COAX_CLOSED_HPP
