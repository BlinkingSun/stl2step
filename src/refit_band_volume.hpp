// Chord volume of a mesh band on a surface of revolution (D-fc-5.7 O-3).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_BAND_VOLUME_HPP
#define STL2STEP_REFIT_BAND_VOLUME_HPP

#include <gp_XYZ.hxx>

#include <cmath>

namespace stl2step {
namespace refit {

// The volume between a band of planar mesh triangles and the surface of
// revolution their vertices lie on, for a section radius that is affine in the
// axial coordinate: rho_S(u) = rho0 + slope * u, u = (p - origin) . axis.
// slope = 0 is a cylinder, any other slope a cone.
//
// THE DERIVATION. In cylindrical coordinates (rho, theta, u) about the axis,
// the field F = rho e_rho / 2 has divergence 1 and is tangent to every
// half-plane theta = const. For a triangle T, its axis wedge W(T) (the
// points between the axis and T along rays perpendicular to the axis) is
// bounded by T, the axis and surfaces made of such rays, so
//   vol W(T) = |flux of F through T| = |c_perp . S_perp| / 2,
// with c the centroid of T and S the area vector (p1-p0)x(p2-p0)/2. F is
// affine on T, so the centroid rule is exact. The surface over the same
// footprint holds  integral rho_S(u)^2 / 2  d theta du.  On the (theta, u)
// triangle of T's vertices that integrand is a quadratic in u, so the
// edge-midpoint rule is exact:
//   patch(T) = |A_thetau| * (g(m01) + g(m12) + g(m20)) / 3,  g = rho_S^2 / 2.
// A shared edge has one footprint in both triangles, so over a band the
// patches tile the band's footprint. Each rim that lies in a plane
// perpendicular to the axis is a straight line u = const in that footprint.
// For a closed band between two such rims, the sum
//   sum_T patch(T) - vol W(T)
// is therefore exactly (frustum volume) - (volume inside the band's facets).
// No side count enters: rims and rows may carry different numbers of
// vertices. A regular N-gon prism or frustum reproduces dVolCylinderSector /
// coneFrustumChordVolume. The open side of a partial band is a
// chord in (theta, u), which is exact to third order in its azimuth step.
//
// The term of ONE triangle is not its own chord volume (the patch uses the
// straight footprint). Only the sum over a band is meaningful. e1, e2 and
// axis must be orthonormal. The result is >= 0 for vertices on the surface.
inline double revolvedChordTerm(const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2,
                                const gp_XYZ& origin, const gp_XYZ& axis, const gp_XYZ& e1,
                                const gp_XYZ& e2, double rho0, double slope) {
    constexpr double kPiL = 3.14159265358979323846;
    const gp_XYZ d0 = p0 - origin, d1 = p1 - origin, d2 = p2 - origin;
    const gp_XYZ c = (d0 + d1 + d2) / 3.0;
    const double thc = std::atan2(c.Dot(e2), c.Dot(e1));
    auto dth = [&](const gp_XYZ& d) {
        double a = std::atan2(d.Dot(e2), d.Dot(e1)) - thc;
        while (a > kPiL) a -= 2.0 * kPiL;
        while (a <= -kPiL) a += 2.0 * kPiL;
        return a;
    };
    const double t0 = dth(d0), t1 = dth(d1), t2 = dth(d2);
    const double u0 = d0.Dot(axis), u1 = d1.Dot(axis), u2 = d2.Dot(axis);
    const double aTU = 0.5 * std::fabs((t1 - t0) * (u2 - u0) - (t2 - t0) * (u1 - u0));
    auto g = [&](double u) {
        const double r = rho0 + slope * u;
        return 0.5 * r * r;
    };
    const double patch =
        aTU * (g(0.5 * (u0 + u1)) + g(0.5 * (u1 + u2)) + g(0.5 * (u2 + u0))) / 3.0;
    const gp_XYZ S = (d1 - d0).Crossed(d2 - d0) * 0.5;
    const gp_XYZ cPerp = c - axis * c.Dot(axis);
    const gp_XYZ sPerp = S - axis * S.Dot(axis);
    const double wedge = 0.5 * std::fabs(cPerp.Dot(sPerp));
    const double v = patch - wedge;
    return std::isfinite(v) ? v : 0.0;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_BAND_VOLUME_HPP
