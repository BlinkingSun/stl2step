// stl2step torus math — the closed-form primitive behind the 1.4.0 toroidal
// round (lane torus-math, SPEC-torus PART B, D-140-6 §3(4)).  Signatures are
// published in src/refit_torus_math.hpp; the derivation is stated there and in
// _team/reports/torus-math.md.
//
// This is deliberately NOT part of src/refit_math.cpp, for the reason
// refit_cone_math.cpp is not (0d67eba): that file is one of the five P1 sources
// the D5.3 include allowlist covers, and the allowlist — a gate instrument, not
// a style rule — admits no project header beyond refit.hpp / refit_internal.hpp.
// Nothing here touches MeshView, Region or any engine state: it is pure gp_
// geometry, which is also what lets tests/unit/torus_math_test.cpp compile this
// one file and certify the math without linking the engine.
//
// Everything below is closed form.  Nothing samples the surface, iterates, or
// approximates, because D-130-2 only lets an edge ship analytic when the program
// computes the supremum deviation in closed form.
//
// The one identity the whole file rests on: the surface is a surface of
// revolution, so the nearest surface point to p lies in p's own meridian
// half-plane, and there the RING torus (R > r) is exactly the profile circle
// (rho - R)^2 + z^2 = r^2 — the whole circle, because R > r keeps it inside
// rho >= 0.  So dist(p, S) = | hypot(rho_p - R, z_p) - r |, with no clamping and
// no "nearest point is on the axis" case.  For a spindle (R < r) that identity
// fails — part of the profile circle has rho < 0, and the untrimmed surface
// there is its mirror image, so the formula is only an upper bound — and the
// spindle is refused rather than certified with a value that is not the
// supremum.
//
// SPDX-License-Identifier: MIT

#include "refit_torus_math.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <gp.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Torus.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

// Resolves the torus into (R, r).  False when the torus is not a ring: non-finite
// radii, no tube (r <= gp::Resolution()), or R - r <= gp::Resolution() — the
// latter is gp_Torus::SetMajorRadius / SetMinorRadius's own construction guard,
// mirrored so the profile circle is guaranteed to lie in rho >= 0 (the header
// identity) and nothing here can see a spindle, a horn or a sphere.
bool torusParams(const gp_Torus& torus, double& R, double& r) {
    R = torus.MajorRadius();
    r = torus.MinorRadius();
    if (!std::isfinite(R) || !std::isfinite(r)) return false;
    if (!(r > gp::Resolution())) return false;
    if (!(R - r > gp::Resolution())) return false;
    return true;
}

bool tolsOk(double angTol, double linTol) {
    return std::isfinite(angTol) && angTol >= 0.0 && std::isfinite(linTol) && linTol >= 0.0;
}

// Distance from x to the closed interval [lo, hi]; 0 inside.
double distToInterval(double x, double lo, double hi) {
    return std::max(0.0, std::max(lo - x, x - hi));
}

constexpr double kTorusNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

const char* torusDevClassName(TorusDevClass c) {
    switch (c) {
        case TorusDevClass::Coaxial:   return "circle-on-torus";
        case TorusDevClass::Unhandled: break;
    }
    return "unhandled-circle-on-torus";
}

bool torusDevClassIsExact(TorusDevClass c) {
    return c == TorusDevClass::Coaxial;
}

bool torusIsRing(const gp_Torus& torus) {
    double R, r;
    return torusParams(torus, R, r);
}

double torusAxialCoord(const gp_Torus& torus, const gp_Pnt& p) {
    return (p.XYZ() - torus.Location().XYZ()).Dot(torus.Axis().Direction().XYZ());
}

double torusRadialCoord(const gp_Torus& torus, const gp_Pnt& p) {
    const gp_XYZ a = torus.Axis().Direction().XYZ();
    gp_XYZ d = p.XYZ() - torus.Location().XYZ();
    d -= a * d.Dot(a);
    return d.Modulus();
}

double torusVOfProfilePoint(const gp_Torus& torus, double rho, double z) {
    double R, r;
    if (!torusParams(torus, R, r)) return kTorusNaN;
    if (!std::isfinite(rho) || !std::isfinite(z)) return kTorusNaN;
    return std::atan2(z, rho - R);
}

bool torusVIsoCircle(const gp_Torus& torus, double v, gp_Circ& out) {
    double R, r;
    if (!torusParams(torus, R, r)) return false;
    if (!std::isfinite(v)) return false;
    const double radius = R + r * std::cos(v);  // > 0: R > r >= r*|cos v|
    const double offset = r * std::sin(v);
    if (!std::isfinite(radius) || !(radius > 0.0) || !std::isfinite(offset)) return false;
    const gp_Dir za = torus.Axis().Direction();
    const gp_Pnt c = torus.Location().Translated(gp_Vec(za) * offset);
    out = gp_Circ(gp_Ax2(c, za, torus.Position().XDirection()), radius);
    return true;
}

double pointTorusDist(const gp_Torus& torus, const gp_Pnt& p) {
    double R, r;
    if (!torusParams(torus, R, r)) return -1.0;
    const double z = torusAxialCoord(torus, p);
    const double rho = torusRadialCoord(torus, p);
    if (!std::isfinite(z) || !std::isfinite(rho)) return -1.0;
    return std::fabs(std::hypot(rho - R, z) - r);
}

// Supremum over the circle of the distance to the untrimmed ring torus.
//
// Exact formula per point (the identity at the head of this file):
//     d = | hypot(rho - R, z) - r |
// A COAXIAL circle of radius r_c at axial offset h has (rho, z) = (r_c, h) at
// every point, so d is constant along it and the supremum is that constant.
//
// The circle is admitted as coaxial when its normal is parallel to the axis
// within angTol and its centre is on the axis within linTol.  Inside those
// tolerances the (rho, z) trajectory is not literally a point, so the value
// returned is the exact supremum over its bounding box —
//   z   in [h - r_c*sin(beta), h + r_c*sin(beta)]        (tight)
//   rho in [max(0, r_c*cos(beta) - d0, d0 - r_c), d0 + r_c]
// (the same two intervals circleOnConeMax uses, for the same reasons) —
// g = hypot(rho - R, z) is convex in (rho, z), so over the box its maximum is
// at a corner and its minimum at the projection of (R, 0) onto the box, and
//     sup |g - r| = max(g_max - r, r - g_min).
// With beta = 0 and d0 = 0 the box is the point (r_c, h), g_max = g_min =
// hypot(r_c - R, h), and the expression IS the constant above — not a bound.
// Otherwise it is >= the truth, which is the only direction D-130-2 permits.
double circleOnTorusMax(const gp_Torus& torus, const gp_Circ& circ, TorusDevClass* clsOut,
                        double angTol, double linTol) {
    if (clsOut) *clsOut = TorusDevClass::Unhandled;
    double R, r;
    if (!torusParams(torus, R, r)) return -1.0;
    if (!tolsOk(angTol, linTol)) return -1.0;
    const double rc = circ.Radius();
    if (!std::isfinite(rc) || rc < 0.0) return -1.0;

    const double h = torusAxialCoord(torus, circ.Location());
    const double d0 = torusRadialCoord(torus, circ.Location());
    if (!std::isfinite(h) || !std::isfinite(d0)) return -1.0;

    // sin(beta) from the CROSS product (exactly zero for a normal that is the
    // axis direction or its negation; only the axis LINE matters), never from
    // sqrt(1 - dot*dot) — see refit_cone_math.cpp for the measured reason.
    // cos(beta) is clamped, never raised: understating it only widens rhoLo.
    gp_XYZ cn = circ.Axis().Direction().XYZ();
    double dot = cn.Dot(torus.Axis().Direction().XYZ());
    if (dot < 0.0) {
        cn = cn.Multiplied(-1.0);
        dot = -dot;
    }
    const double cosB = std::min(1.0, dot);
    const double sinB = cn.Crossed(torus.Axis().Direction().XYZ()).Modulus();
    if (!(sinB <= angTol)) return -1.0;  // tilted: not coaxial (tier 2, counted)
    if (!(d0 <= linTol)) return -1.0;    // off-axis: not coaxial (tier 2, counted)

    const double zLo = h - rc * sinB;
    const double zHi = h + rc * sinB;
    const double rhoLo = std::max(0.0, std::max(rc * cosB - d0, d0 - rc));
    const double rhoHi = d0 + rc;

    const double gMax = std::hypot(std::max(std::fabs(rhoLo - R), std::fabs(rhoHi - R)),
                                   std::max(std::fabs(zLo), std::fabs(zHi)));
    const double gMin = std::hypot(distToInterval(R, rhoLo, rhoHi), distToInterval(0.0, zLo, zHi));

    if (clsOut) *clsOut = TorusDevClass::Coaxial;
    return std::max(gMax - r, r - gMin);
}

}  // namespace refit
}  // namespace stl2step
