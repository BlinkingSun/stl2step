// Chord volume of a mesh band on a surface of revolution (D-fc-5.7 O-3).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_BAND_VOLUME_HPP
#define STL2STEP_REFIT_BAND_VOLUME_HPP

#include <gp_XYZ.hxx>

#include <algorithm>
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

// D-gid-8.8. The same volume for a cone band that holds its own apex (a closed
// drill point: a fan about a vertex on the axis). revolvedChordTerm cannot take
// it: its footprint is the straight (theta, u) triangle of the three vertices,
// and the apex vertex has no azimuth.
//
// THE BOUNDARY FORM. With H(u) = integral from the apex to u of rho_S^2 / 2,
// Green's theorem turns the patch over a footprint into  - (loop of H du-theta),
//   patch(T) = - sum over T's edges p -> q (winding order) of dtheta(p, q) * Hbar(p, q),
// dtheta = the edge's own wrapped azimuth step, Hbar = H's mean along the edge.
// u is linear along a straight footprint edge and H is a cubic in u, so
// Simpson's rule gives Hbar exactly. The wedge is the signed flux of
// F = rho e_rho / 2 through T, c_perp . S_perp / 2 (F is linear on T).
//
// Over a band, an edge shared by two triangles enters with opposite sense and
// the same dtheta, so it cancels whatever its endpoints' azimuths are. Every
// edge at the apex is such an edge, so only the band's rim survives:
//   sum_T term(T) = - (rim loop of H dtheta) - (flux of F through the facets).
// H(apex) = 0, so the footprint runs from the apex to the rim. That is the
// volume between the cone and the facets: F is tangent to the rays through the
// rim and vanishes on the axis. A rim in a plane perpendicular to the axis is
// u = const, a straight footprint edge, so the sum is exact; any other rim is
// exact to third order in its azimuth step, as for revolvedChordTerm.
//
// rho_S(u) = slope * u, u = (p - apex) . axis. The term of ONE triangle is not
// its own volume; only the sum over a closed fan about the axis is meaningful.
// The sum is positive where the cone lies on the side the facets' winding
// normals point to. e1, e2 and axis must be orthonormal.
inline double apexConeChordTerm(const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2,
                                const gp_XYZ& apex, const gp_XYZ& axis, const gp_XYZ& e1,
                                const gp_XYZ& e2, double slope) {
    constexpr double kPiL = 3.14159265358979323846;
    const gp_XYZ d[3] = {p0 - apex, p1 - apex, p2 - apex};
    const double k = slope * slope / 6.0;
    auto H = [&](double u) { return k * u * u * u; };
    double th[3], u[3];
    for (int i = 0; i < 3; ++i) {
        th[i] = std::atan2(d[i].Dot(e2), d[i].Dot(e1));
        u[i] = d[i].Dot(axis);
    }
    double patch = 0.0;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        double dt = th[j] - th[i];
        while (dt > kPiL) dt -= 2.0 * kPiL;
        while (dt <= -kPiL) dt += 2.0 * kPiL;
        const double hBar = (H(u[i]) + 4.0 * H(0.5 * (u[i] + u[j])) + H(u[j])) / 6.0;
        patch -= dt * hBar;
    }
    const gp_XYZ c = (d[0] + d[1] + d[2]) / 3.0;
    const gp_XYZ S = (d[1] - d[0]).Crossed(d[2] - d[0]) * 0.5;
    const gp_XYZ cPerp = c - axis * c.Dot(axis);
    const gp_XYZ sPerp = S - axis * S.Dot(axis);
    const double v = patch - 0.5 * cPerp.Dot(sPerp);
    return std::isfinite(v) ? v : 0.0;
}

// D-gid-8.18. The rim term of revolvedChordTerm for one band rim edge pi -> pj
// (in the sense of the triangle that owns it). revolvedChordTerm takes each
// edge's image on the surface as straight in (theta, u); the field's rays carry
// a chord point to the surface at its own azimuth and height, so the true image
// of the rim chord is theta(s) = atan2 of the chord's projection, u(s) linear.
// The two differ when the edge is oblique (it changes both azimuth and height):
// a generator-crossing end of a partial band, or a rim on a plane that is not
// perpendicular to the axis. Interior edges cancel in a band. With
// G(u) = integral rho_S(u)^2 / 2 du, Green on the (theta, u) footprint gives
//   rim = -integral_0^1 G(u(s)) (theta'(s) - dtheta) ds,
//   theta'(s) = (x_i dy - y_i dx) / (x(s)^2 + y(s)^2),
// for the footprint oriented counter-clockwise; revolvedChordTerm takes every
// footprint with a positive area, so the caller multiplies by the sign of the
// owning triangle's (theta, u) area (revolvedParamSign). The integrand is
// smooth on a chord that does not pass the axis; a fixed 16-point
// Gauss-Legendre rule evaluates it, no step or tolerance is chosen. Zero for a
// meridian (dtheta = 0 and theta' = 0) or a parallel (du = 0) edge.
inline double coneRimTerm(const gp_XYZ& pi, const gp_XYZ& pj, const gp_XYZ& origin,
                          const gp_XYZ& axis, const gp_XYZ& e1, const gp_XYZ& e2, double rho0,
                          double slope) {
    constexpr double kPiL = 3.14159265358979323846;
    const gp_XYZ di = pi - origin, dj = pj - origin;
    const double xi = di.Dot(e1), yi = di.Dot(e2), ui = di.Dot(axis);
    const double xj = dj.Dot(e1), yj = dj.Dot(e2), uj = dj.Dot(axis);
    double dth = std::atan2(yj, xj) - std::atan2(yi, xi);
    while (dth > kPiL) dth -= 2.0 * kPiL;
    while (dth <= -kPiL) dth += 2.0 * kPiL;
    const double dx = xj - xi, dy = yj - yi, du = uj - ui;
    const double cross = xi * dy - yi * dx;
    auto G = [&](double u) {
        return 0.5 * rho0 * rho0 * u + 0.5 * rho0 * slope * u * u + slope * slope * u * u * u / 6.0;
    };
    static const double kX[8] = {0.095012509837637441, 0.28160355077925892, 0.45801677765722737,
                                 0.61787624440264377,  0.75540440835500311, 0.86563120238783164,
                                 0.9445750230732326,   0.98940093499164994};
    static const double kW[8] = {0.18945061045506859, 0.18260341504492361, 0.16915651939500262,
                                 0.14959598881657665, 0.12462897125553395, 0.095158511682492994,
                                 0.062253523938647609, 0.027152459411754055};
    double acc = 0.0;
    for (int k = 0; k < 8; ++k)
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            const double s = 0.5 * (1.0 + sgn * kX[k]);
            const double w = 0.5 * kW[k];
            const double x = xi + s * dx, y = yi + s * dy, u = ui + s * du;
            const double q = x * x + y * y;
            if (!(q > 0.0)) return 0.0;
            acc += w * G(u) * (cross / q - dth);
        }
    const double out = -acc;
    return std::isfinite(out) ? out : 0.0;
}

// The sign of the (theta, u) area revolvedChordTerm takes in absolute value
// (+1 counter-clockwise about axis = e1 x e2), 0 for a degenerate footprint.
inline int revolvedParamSign(const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2,
                             const gp_XYZ& origin, const gp_XYZ& axis, const gp_XYZ& e1,
                             const gp_XYZ& e2) {
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
    const double a = (t1 - t0) * (u2 - u0) - (t2 - t0) * (u1 - u0);
    return a > 0.0 ? 1 : (a < 0.0 ? -1 : 0);
}

// The volume between one planar mesh triangle and a torus (D-gid-7.6), by the
// same field F = rho e_rho / 2 as revolvedChordTerm (divergence 1, tangent to
// every ray perpendicular to the axis).
//
// THE TORUS. x = o + (Rmaj + r cos psi)(cos phi e1 + sin phi e2) + r sin psi a,
// with (e1, e2, a) right-handed. x_phi x x_psi = r (Rmaj + r cos psi) m, m the
// outward tube normal, so a triangle counter-clockwise in (phi, psi) faces out
// of the tube. The flux of F through the surface is
//   h(psi) dphi dpsi,  h = (r / 2) (Rmaj + r cos psi)^2 cos psi
//                        = a0 + a1 cos psi + a2 cos 2psi + a3 cos 3psi,
//   a0 = a2 = Rmaj r^2 / 2,  a1 = (r / 2)(Rmaj^2 + 3 r^2 / 4),  a3 = r^3 / 8.
// By Green, the flux over a straight-edged (phi, psi) triangle is the edge sum
// of  integral phi h(psi) dpsi,  which is closed form on a straight edge:
//   dpsi [a0 phim + sum_k a_k (phim cos(k psim) S0(k dpsi)
//                               - dphi sin(k psim) S1(k dpsi))],
//   S0(b) = int_{-1/2}^{1/2} cos(t b) dt,  S1(b) = int_{-1/2}^{1/2} t sin(t b) dt.
// The flux through the triangle itself is c_perp . S_perp / 2 (F is linear).
// The difference P - W is the volume between the triangle and the torus over
// the triangle's footprint, signed by the triangle's winding (positive where
// the torus lies on the side its winding normal points to).
//
// A vertex maps to the torus point at its own height and azimuth (the ray the
// field is tangent to): psi = atan2(v, +-sqrt(r^2 - v^2)), the branch of its
// side of the tube, v clamped to [-r, r]. Over a band, every interior edge
// enters twice with opposite sense, so only the band's rim vertices and edges
// matter: the sum is exact for a rim made of meridians (constant phi) and
// parallels (constant height); an oblique rim edge adds torusRimTerm. No side
// count, row count or tessellation regularity enters.
// side: +1 / -1 takes every vertex to the outer / inner half of the tube (the
// band's own half, when all of it lies on one); 0 takes each vertex's own.
// paramArea receives the triangle's signed (phi, psi) area; its sign over a
// band is the band's winding relative to the outward tube normal.
struct TorusChordTerm {
    double vol = 0.0;
    double paramArea = 0.0;
};

inline TorusChordTerm torusChordTerm(const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2,
                                     const gp_XYZ& origin, const gp_XYZ& axis,
                                     const gp_XYZ& e1, const gp_XYZ& e2, double Rmaj,
                                     double rMin, int side = 0) {
    constexpr double kPiL = 3.14159265358979323846;
    TorusChordTerm out;
    if (!(rMin > 0.0) || !std::isfinite(Rmaj) || !std::isfinite(rMin)) return out;
    const gp_XYZ* P[3] = {&p0, &p1, &p2};
    double phi[3], psi[3];
    for (int i = 0; i < 3; ++i) {
        const gp_XYZ d = *P[i] - origin;
        const double v = d.Dot(axis);
        const double x = d.Dot(e1), y = d.Dot(e2);
        phi[i] = std::atan2(y, x);
        const double vc = std::max(-rMin, std::min(rMin, v));
        // (r - v)(r + v): exact zero on the rims |v| = r, whatever the FP contraction.
        const double c = std::sqrt(std::max(0.0, (rMin - vc) * (rMin + vc)));
        const bool outerV = side != 0 ? side > 0 : std::hypot(x, y) >= Rmaj;
        psi[i] = std::atan2(vc, outerV ? c : -c);
    }
    auto wrap = [&](double a) {
        while (a > kPiL) a -= 2.0 * kPiL;
        while (a <= -kPiL) a += 2.0 * kPiL;
        return a;
    };
    double fr[3], sr[3];
    for (int i = 0; i < 3; ++i) {
        fr[i] = wrap(phi[i] - phi[0]);
        sr[i] = wrap(psi[i] - psi[0]);
    }
    // S0, S1 in closed form; below |b| = 0.1 by their Taylor series (the same
    // functions; the closed form of S1 cancels there). Truncation < 1e-16.
    auto S0 = [](double b) {
        if (std::fabs(b) < 0.1) {
            const double b2 = b * b;
            return 1.0 - b2 / 24.0 + b2 * b2 / 1920.0 - b2 * b2 * b2 / 322560.0;
        }
        return 2.0 * std::sin(0.5 * b) / b;
    };
    auto S1 = [](double b) {
        if (std::fabs(b) < 0.1) {
            const double b2 = b * b;
            return b * (1.0 / 12.0 - b2 / 480.0 + b2 * b2 / 53760.0 - b2 * b2 * b2 / 11612160.0);
        }
        return 2.0 * std::sin(0.5 * b) / (b * b) - std::cos(0.5 * b) / b;
    };
    const double r = rMin;
    const double a0 = 0.5 * Rmaj * r * r;
    const double ak[4] = {a0, 0.5 * r * (Rmaj * Rmaj + 0.75 * r * r), a0, 0.125 * r * r * r};
    double patch = 0.0;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const double dphi = fr[j] - fr[i], dpsi = sr[j] - sr[i];
        const double phim = 0.5 * (fr[i] + fr[j]);
        const double psim = psi[0] + 0.5 * (sr[i] + sr[j]);
        double e = ak[0] * phim;
        for (int k = 1; k <= 3; ++k) {
            const double kk = (double)k;
            e += ak[k] * (phim * std::cos(kk * psim) * S0(kk * dpsi) -
                          dphi * std::sin(kk * psim) * S1(kk * dpsi));
        }
        patch += dpsi * e;
    }
    const gp_XYZ d0 = p0 - origin, d1 = p1 - origin, d2 = p2 - origin;
    const gp_XYZ c = (d0 + d1 + d2) / 3.0;
    const gp_XYZ S = (d1 - d0).Crossed(d2 - d0) * 0.5;
    const gp_XYZ cPerp = c - axis * c.Dot(axis);
    const gp_XYZ sPerp = S - axis * S.Dot(axis);
    const double wedge = 0.5 * cPerp.Dot(sPerp);
    out.vol = patch - wedge;
    out.paramArea = 0.5 * (fr[1] * sr[2] - fr[2] * sr[1]);
    if (!std::isfinite(out.vol) || !std::isfinite(out.paramArea)) out = TorusChordTerm{};
    return out;
}

// The rim term of torusChordTerm for one band rim edge pi -> pj (in the sense
// of the triangle that owns it). torusChordTerm takes each edge's image on the
// torus as straight in (phi, psi). The field's rays carry a chord point to the
// torus at its own azimuth and height, so the true image of the rim chord is
//   phi(s) = atan2 of the chord's projection,  v(s) linear,  psi = Psi(v(s)),
// which differs from the straight edge when the edge crosses rows obliquely
// (phi(s) is not linear because the chord's radius varies along a torus
// profile, and psi is not linear in v). Interior edges cancel in a band, so
// only rim edges carry this term. With H(psi) = integral h dpsi,
//   rim = -integral_0^1 [H(Psi(v(s))) phi'(s) - dphi H(psi_lin(s))] ds,
//   phi'(s) = (x_i dy - y_i dx) / (x(s)^2 + y(s)^2),
// the difference of the Green terms over the two curves (same endpoints).
// A chord point above the profile's top (or below its bottom), |v| > r, has no
// ray image; it maps to the rim circle at its azimuth, Q(s), and the ribbon
// X(s, t) = p(s) + t (Q(s) - p(s)) is then not made of the field's rays, so its
// flux, integral F . (X_s x X_t) ds dt (oriented p_i, p_j, Q_j, Q_i), is added.
// [0, 1] is split where |v(s)| = r; on each piece the integrands are smooth
// (H(Psi(v)) has derivative rho_S(v)^2 / 2 in v) and a fixed 16-point
// Gauss-Legendre rule evaluates them, the ribbon's t-integral (quadratic) by
// 2-point Gauss; no step or tolerance is chosen. side as in torusChordTerm;
// with side 0 the term is zero when the two ends lie on different halves.
inline double torusRimTerm(const gp_XYZ& pi, const gp_XYZ& pj, const gp_XYZ& origin,
                           const gp_XYZ& axis, const gp_XYZ& e1, const gp_XYZ& e2, double Rmaj,
                           double rMin, int side = 0) {
    constexpr double kPiL = 3.14159265358979323846;
    if (!(rMin > 0.0) || !std::isfinite(Rmaj) || !std::isfinite(rMin)) return 0.0;
    const gp_XYZ di = pi - origin, dj = pj - origin;
    const double xi = di.Dot(e1), yi = di.Dot(e2), vi = di.Dot(axis);
    const double xj = dj.Dot(e1), yj = dj.Dot(e2), vj = dj.Dot(axis);
    bool outer = side > 0;
    if (side == 0) {
        outer = std::hypot(xi, yi) >= Rmaj;
        if (outer != (std::hypot(xj, yj) >= Rmaj)) return 0.0;
    }
    const double r = rMin;
    auto Psi = [&](double v) {
        const double s = std::asin(std::max(-1.0, std::min(1.0, v / r)));
        return outer ? s : kPiL - s;
    };
    const double a0 = 0.5 * Rmaj * r * r;
    const double a1 = 0.5 * r * (Rmaj * Rmaj + 0.75 * r * r), a3 = 0.125 * r * r * r;
    auto H = [&](double p) {
        return a0 * p + a1 * std::sin(p) + 0.5 * a0 * std::sin(2.0 * p) +
               (a3 / 3.0) * std::sin(3.0 * p);
    };
    double dphi = std::atan2(yj, xj) - std::atan2(yi, xi);
    while (dphi > kPiL) dphi -= 2.0 * kPiL;
    while (dphi <= -kPiL) dphi += 2.0 * kPiL;
    const double dx = xj - xi, dy = yj - yi, dv = vj - vi;
    const double cross = xi * dy - yi * dx;
    const double Pi = Psi(vi), Pj = Psi(vj);
    static const double kX[8] = {0.095012509837637441, 0.28160355077925892, 0.45801677765722737,
                                 0.61787624440264377,  0.75540440835500311, 0.86563120238783164,
                                 0.9445750230732326,   0.98940093499164994};
    static const double kW[8] = {0.18945061045506859, 0.18260341504492361, 0.16915651939500262,
                                 0.14959598881657665, 0.12462897125553395, 0.095158511682492994,
                                 0.062253523938647609, 0.027152459411754055};
    // Pieces of [0, 1] between the crossings of v = +r and v = -r.
    double cuts[4] = {0.0, 1.0, 1.0, 1.0};
    int nCut = 1;
    if (dv != 0.0)
        for (double lim : {r, -r}) {
            const double sc = (lim - vi) / dv;
            if (sc > 0.0 && sc < 1.0) cuts[nCut++] = sc;
        }
    cuts[nCut++] = 1.0;
    std::sort(cuts, cuts + nCut);
    const double t2[2] = {0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0)};
    double acc = 0.0, rib = 0.0;
    for (int piece = 0; piece + 1 < nCut; ++piece) {
        const double s0 = cuts[piece], L = cuts[piece + 1] - s0;
        if (!(L > 0.0)) continue;
        const double vMid = vi + (s0 + 0.5 * L) * dv;
        const bool clamped = vMid > r || vMid < -r;
        const double vc = vMid > r ? r : -r;
        for (int k = 0; k < 8; ++k) {
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double s = s0 + L * 0.5 * (1.0 + sgn * kX[k]);
                const double w = 0.5 * kW[k] * L;
                const double x = xi + s * dx, y = yi + s * dy, v = vi + s * dv;
                const double q = x * x + y * y;
                if (!(q > 0.0)) return 0.0;
                acc += w * (H(Psi(v)) * cross / q - dphi * H(Pi + s * (Pj - Pi)));
                if (!clamped) continue;
                // Ribbon to the rim circle Q(s) = Rmaj rho_hat(phi(s)) + vc a.
                const double rho = std::sqrt(q);
                const double ux = x / rho, uy = y / rho;
                const double dphids = cross / q;
                const gp_XYZ P(x, y, v), dP(dx, dy, dv);
                const gp_XYZ Q(Rmaj * ux, Rmaj * uy, vc);
                const gp_XYZ dQ(-Rmaj * uy * dphids, Rmaj * ux * dphids, 0.0);
                const gp_XYZ D = Q - P;
                for (double t : t2) {
                    const gp_XYZ X = P + D * t;
                    const gp_XYZ Xs = dP + (dQ - dP) * t;
                    const gp_XYZ Fv(0.5 * X.X(), 0.5 * X.Y(), 0.0);
                    rib += w * 0.5 * Fv.Dot(Xs.Crossed(D));
                }
            }
        }
    }
    const double out = rib - acc;
    return std::isfinite(out) ? out : 0.0;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_BAND_VOLUME_HPP
