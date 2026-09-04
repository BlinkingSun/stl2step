// stl2step cone math — the closed-form primitives behind the 1.3.0 chamfer
// frustum (lane 130-CONE-MATH, SPEC-130-cone-math).  Signatures are published in
// src/refit_cone_math.hpp; the derivations are in
// docs/1.3.0-linkage-detectors.md.
//
// This is deliberately NOT part of src/refit_math.cpp.  That file is one of the
// five P1 sources the D5.3 include allowlist covers, and the allowlist — a gate
// instrument, not a style rule — admits no bare gp.hxx and no project header
// beyond refit.hpp / refit_internal.hpp.  Both are needed here, so the cone math
// lives in its own TU and P1 keeps exactly the surface it was measured with.
// Nothing here touches MeshView, Region or any engine state: it is pure gp_
// geometry, which is also what lets tests/unit/cone_math_test.cpp compile this
// one file and certify the math without linking the engine.
//
// Everything below is closed form.  Nothing samples the surface, iterates, or
// approximates, because D-130-2 only lets an edge ship analytic when the program
// computes the supremum deviation in closed form.
//
// The one identity the whole file rests on: the meridian plane through a point
// meets the untrimmed (double-nappe) cone in TWO full lines, and both lines lie
// entirely on the OCCT surface — S(u, v) for v past the apex has a negative
// radius factor, which is the point at meridian u + PI.  So the perpendicular
// foot is always a genuine surface point, there is no "nearest point is the
// apex" case to clamp, and the meridian trace of the surface is simply
// rho = |rho_cone(z)|.
//
// SPDX-License-Identifier: MIT

#include "refit_cone_math.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <gp.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;

// Resolves the cone into (R, tan Ang, cos Ang, sin Ang). False when the cone
// could not have come from a legal gp_Cone: non-finite parameters, a negative
// reference radius, or a semi-angle outside ]-PI/2, PI/2[ \ {0} (gp_Cone's own
// construction guard, mirrored so nothing here can divide by zero).
bool coneParams(const gp_Cone& cone, double& R, double& ta, double& ca, double& sa) {
    R = cone.RefRadius();
    const double ang = cone.SemiAngle();
    if (!std::isfinite(R) || R < 0.0 || !std::isfinite(ang)) return false;
    const double aa = std::fabs(ang);
    if (aa <= gp::Resolution() || aa >= kPi * 0.5 - gp::Resolution()) return false;
    ca = std::cos(ang);
    sa = std::sin(ang);
    ta = sa / ca;
    return std::isfinite(ca) && std::isfinite(sa) && std::isfinite(ta) && ca > 0.0 &&
           std::fabs(sa) > 0.0;
}

bool tolsOk(double angTol, double linTol) {
    return std::isfinite(angTol) && angTol >= 0.0 && std::isfinite(linTol) && linTol >= 0.0;
}

// The circle centred on the cone axis at axial height z, oriented off the
// cone's own frame so repeated calls are byte-identical (I5).
gp_Circ coaxialCircle(const gp_Cone& cone, double z, double radius) {
    const gp_Dir za = cone.Axis().Direction();
    const gp_Pnt c = cone.Location().Translated(gp_Vec(za) * z);
    return gp_Circ(gp_Ax2(c, za, cone.Position().XDirection()), radius);
}

constexpr double kConeNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

const char* coneDevClassName(ConeDevClass c) {
    switch (c) {
        case ConeDevClass::CoaxialRim:   return "circle-on-cone";
        case ConeDevClass::PerpPlane:    return "circle-on-cone-perp";
        case ConeDevClass::GeneralBound: return "circle-on-cone-bound";
        case ConeDevClass::Unhandled:    break;
    }
    return "unhandled-circle-on-cone";
}

bool coneDevClassIsExact(ConeDevClass c) {
    return c == ConeDevClass::CoaxialRim || c == ConeDevClass::PerpPlane;
}

double coneAxialCoord(const gp_Cone& cone, const gp_Pnt& p) {
    return (p.XYZ() - cone.Location().XYZ()).Dot(cone.Axis().Direction().XYZ());
}

double coneRadialCoord(const gp_Cone& cone, const gp_Pnt& p) {
    const gp_XYZ a = cone.Axis().Direction().XYZ();
    gp_XYZ d = p.XYZ() - cone.Location().XYZ();
    d -= a * d.Dot(a);
    return d.Modulus();
}

// Every scalar helper below returns NaN on a cone that coneParams rejects: a
// silent 0 or inf would be indistinguishable from a real answer downstream.
double coneRadiusAt(const gp_Cone& cone, double z) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return kConeNaN;
    return R + z * ta;
}

double coneApexZ(const gp_Cone& cone) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return kConeNaN;
    return -R / ta;
}

double coneApexV(const gp_Cone& cone) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return kConeNaN;
    return -R / sa;
}

double coneVAtZ(const gp_Cone& cone, double z) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return kConeNaN;
    return z / ca;
}

double coneVAtRadius(const gp_Cone& cone, double rho) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return kConeNaN;
    return (rho - R) / sa;
}

bool coneFrustumVRange(const gp_Cone& cone, double zA, double zB, double& vLo, double& vHi) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return false;
    if (!std::isfinite(zA) || !std::isfinite(zB)) return false;
    if (!(std::fabs(zA - zB) > 0.0)) return false;
    // rho_cone(z) = tan(Ang) * (z - z_apex): two strictly positive rim radii put
    // both rims on the same side of the apex, which is exactly apex exclusion.
    if (!(R + zA * ta > 0.0) || !(R + zB * ta > 0.0)) return false;
    const double a = zA / ca;
    const double b = zB / ca;
    vLo = std::min(a, b);
    vHi = std::max(a, b);
    const double vApex = -R / sa;
    if (vApex >= vLo && vApex <= vHi) return false;  // defensive: unreachable above
    return true;
}

bool coneFrustumVRangeFromRadii(const gp_Cone& cone, double rhoA, double rhoB, double& vLo,
                                double& vHi) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return false;
    if (!std::isfinite(rhoA) || !std::isfinite(rhoB)) return false;
    if (!(rhoA > 0.0) || !(rhoB > 0.0)) return false;
    if (!(std::fabs(rhoA - rhoB) > 0.0)) return false;
    const double a = (rhoA - R) / sa;
    const double b = (rhoB - R) / sa;
    vLo = std::min(a, b);
    vHi = std::max(a, b);
    const double vApex = -R / sa;
    if (vApex >= vLo && vApex <= vHi) return false;  // defensive: unreachable above
    return true;
}

bool coneFromRims(const gp_Ax3& rimA, double rhoA, double rhoB, double height, gp_Cone& out) {
    if (!std::isfinite(rhoA) || !std::isfinite(rhoB) || !std::isfinite(height)) return false;
    if (!(rhoA > 0.0) || !(rhoB > 0.0)) return false;
    if (!(std::fabs(height) > 0.0)) return false;
    const double dR = rhoB - rhoA;
    if (!(std::fabs(dR) > 0.0)) return false;
    const double ang = std::atan(dR / height);  // in ]-PI/2, PI/2[ for any finite ratio
    if (!std::isfinite(ang)) return false;
    // gp_Cone's own construction guard, checked before we can trip it.
    const double aa = std::fabs(ang);
    if (aa <= gp::Resolution() || aa >= kPi * 0.5 - gp::Resolution()) return false;
    out = gp_Cone(rimA, ang, rhoA);
    return true;
}

double pointConeDist(const gp_Cone& cone, const gp_Pnt& p) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return -1.0;
    const double z = coneAxialCoord(cone, p);
    const double rho = coneRadialCoord(cone, p);
    if (!std::isfinite(z) || !std::isfinite(rho)) return -1.0;
    return ca * std::fabs(rho - std::fabs(R + z * ta));
}

// Supremum over the circle of the distance to the untrimmed cone.
//
// Exact formula per point (see the identity at the head of this file):
//     d = cos(Ang) * | rho - |rho_cone(z)| |
// so the supremum is the maximum of |rho - xi| over the circle's (rho, z)
// trajectory, with xi = |rho_cone(z)|. Both coordinates are bounded EXACTLY in
// closed form:
//   z   in [z0 - r*sin(beta), z0 + r*sin(beta)]   (tight; the circle's axial
//       spread is r*sin(beta) by construction, beta = normal-to-axis angle)
//   rho in [max(0, r*cos(beta) - d0, d0 - r), d0 + r]
//       The circle projects onto the plane perpendicular to the axis as an
//       ellipse of semi-axes r and r*cos(beta) centred at distance d0 from the
//       axis, so the projected offset has modulus in [r*cos(beta), r] and the
//       triangle inequality gives both ends.
// The maximum of |rho - xi| over that box is attained at a corner, hence
//     sup <= cos(Ang) * max(rho_hi - xi_min, xi_max - rho_lo).
//
// When the circle's plane is perpendicular to the axis, sin(beta) = 0 collapses
// the z interval to a point and cos(beta) = 1 collapses the rho interval to the
// EXACT span [|d0 - r|, d0 + r] — the inequality above becomes an equality and
// the result is the supremum itself, not a bound. That is the only difference
// between the exact classes and the general one: one evaluator, two regimes.
double circleOnConeMax(const gp_Cone& cone, const gp_Circ& circ, ConeDevClass* clsOut,
                       double angTol, double linTol) {
    if (clsOut) *clsOut = ConeDevClass::Unhandled;
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa)) return -1.0;
    if (!tolsOk(angTol, linTol)) return -1.0;
    const double r = circ.Radius();
    if (!std::isfinite(r) || r < 0.0) return -1.0;

    const double z0 = coneAxialCoord(cone, circ.Location());
    const double d0 = coneRadialCoord(cone, circ.Location());
    if (!std::isfinite(z0) || !std::isfinite(d0)) return -1.0;

    // sin(beta) comes from the CROSS product, never from sqrt(1 - dot*dot).
    // Two unit directions that are bitwise (anti)parallel still dot to 1 -+ 1ulp,
    // and sqrt of that ulp is ~1e-8 — eight orders above Precision::Angular(), so
    // the sqrt form reports a spurious tilt for a third of all cone orientations
    // and demotes an exact rim to GeneralBound (measured: 669 of 2000 random
    // placements). The cross product of identical coordinates is exactly zero,
    // and near-parallel it carries absolute, not square-rooted, error.
    // Only the axis LINE matters, so fold the normal onto the axis first; the
    // negation is exact and leaves the cross modulus bit-identical.
    gp_XYZ cn = circ.Axis().Direction().XYZ();
    double dot = cn.Dot(cone.Axis().Direction().XYZ());
    if (dot < 0.0) {
        cn = cn.Multiplied(-1.0);
        dot = -dot;
    }
    // cosB may fall 1ulp short of 1; that only shrinks rhoLo, which enlarges the
    // bound. Over-stating it would shrink the bound, so it is clamped, not raised.
    const double cosB = std::min(1.0, dot);
    const double sinB = cn.Crossed(cone.Axis().Direction().XYZ()).Modulus();

    const double zLo = z0 - r * sinB;
    const double zHi = z0 + r * sinB;
    const double xiA = std::fabs(R + zLo * ta);
    const double xiB = std::fabs(R + zHi * ta);
    const double zApex = -R / ta;
    const double xiMax = std::max(xiA, xiB);
    const double xiMin = (zApex >= zLo && zApex <= zHi) ? 0.0 : std::min(xiA, xiB);

    const double rhoLo = std::max(0.0, std::max(r * cosB - d0, d0 - r));
    const double rhoHi = d0 + r;

    if (clsOut) {
        *clsOut = (sinB > angTol)     ? ConeDevClass::GeneralBound
                  : (d0 > linTol)     ? ConeDevClass::PerpPlane
                                      : ConeDevClass::CoaxialRim;
    }
    return ca * std::max(rhoHi - xiMin, xiMax - rhoLo);
}

bool coneIntPlaneCircle(const gp_Cone& cone, const gp_Pln& pln, gp_Circ& out, double angTol,
                        double linTol) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa) || !tolsOk(angTol, linTol)) return false;
    const gp_Dir za = cone.Axis().Direction();
    // Perpendicular plane <=> its normal is parallel (or antiparallel) to the axis.
    if (!pln.Axis().Direction().IsParallel(za, angTol)) return false;
    const double zp = coneAxialCoord(cone, pln.Location());
    if (!std::isfinite(zp)) return false;
    const double rho = R + zp * ta;
    if (!(std::fabs(rho) > linTol)) return false;  // apex plane: a point, not a circle
    out = coaxialCircle(cone, zp, std::fabs(rho));
    return true;
}

bool coneIntCylCircle(const gp_Cone& cone, const gp_Cylinder& cyl, gp_Circ& out, double angTol,
                      double linTol) {
    double R, ta, ca, sa;
    if (!coneParams(cone, R, ta, ca, sa) || !tolsOk(angTol, linTol)) return false;
    const gp_Dir za = cone.Axis().Direction();
    if (!cyl.Axis().Direction().IsParallel(za, angTol)) return false;
    if (!(coneRadialCoord(cone, cyl.Location()) <= linTol)) return false;  // coaxial only
    const double Rc = cyl.Radius();
    if (!std::isfinite(Rc) || !(Rc > linTol)) return false;
    // rho_cone(z) = R + z*tan(Ang) = Rc. The far-nappe root -Rc lies past the
    // apex and so outside every frustum's v-range; it is not returned.
    const double z = (Rc - R) / ta;
    if (!std::isfinite(z)) return false;
    out = coaxialCircle(cone, z, Rc);
    return true;
}

bool coneIntConeCircle(const gp_Cone& a, const gp_Cone& b, gp_Circ& out, double angTol,
                       double linTol) {
    double Ra, taA, caA, saA;
    double Rb, taB, caB, saB;
    if (!coneParams(a, Ra, taA, caA, saA) || !coneParams(b, Rb, taB, caB, saB)) return false;
    if (!tolsOk(angTol, linTol)) return false;
    const gp_Dir za = a.Axis().Direction();
    if (!b.Axis().Direction().IsParallel(za, angTol)) return false;
    if (!(coneRadialCoord(a, b.Location()) <= linTol)) return false;  // coaxial only

    // Re-express b in a's frame: a point at a-height z is at b-height
    // sgn*(z - zb0), so rho_b(z) = (Rb - zb0*tb) + z*tb with tb = sgn*tan(AngB).
    const double sgn = (b.Axis().Direction().Dot(za) >= 0.0) ? 1.0 : -1.0;
    const double zb0 = coneAxialCoord(a, b.Location());
    if (!std::isfinite(zb0)) return false;
    const double tb = sgn * taB;
    const double RbA = Rb - zb0 * tb;
    const double dt = taA - tb;
    // |d tan| >= |d angle| for |Ang| < PI/2, so this refuses generators that are
    // parallel within angTol (and with them the "same cone" case).
    if (!(std::fabs(dt) > angTol)) return false;
    const double z = (RbA - Ra) / dt;
    if (!std::isfinite(z)) return false;
    const double rho = Ra + z * taA;
    if (!(rho > linTol)) return false;  // crossing at or behind an apex: not shippable
    out = coaxialCircle(a, z, rho);
    return true;
}


// ---------------------------------------------------------------------------
// 7. facet-to-frustum volume  (D-130-11(2))
// ---------------------------------------------------------------------------

namespace {

// Shared guard and the two factors both public entry points need:
//   segFac = (gamma - sin gamma), the circular-segment shape factor,
//   slant  = hypot(h, R2 - R1),   the length of one N-gon generator.
// False on anything that could not describe a frustum, so neither entry point
// can return a plausible-looking number for input that has no volume.
bool frustumVolParams(double R1, double R2, double height, int nSides, double& gamma,
                      double& segFac, double& slant) {
    if (nSides < 3) return false;
    if (!std::isfinite(R1) || !std::isfinite(R2) || !std::isfinite(height)) return false;
    if (!(R1 > 0.0) || !(R2 > 0.0)) return false;
    const double h = std::fabs(height);
    if (!(h > 0.0)) return false;
    gamma = 2.0 * kPi / static_cast<double>(nSides);
    if (!(gamma < kPi)) return false;
    segFac = gamma - std::sin(gamma);
    // The facet's own height, NOT the generator's length. A facet is the
    // trapezoid between the two parallel chords; its chords sit at axis
    // distance rho*cos(gamma/2), so the perpendicular offset between them is
    // hypot(h, dR*cos(gamma/2)). The generator hypot(h, dR) joins two rim
    // VERTICES and is longer -- measured on the plate's frustum the difference
    // is 3.1e-4 relative, which the unit test's triangulated skin sees.
    slant = std::hypot(h, (R2 - R1) * std::cos(0.5 * gamma));
    return std::isfinite(segFac) && std::isfinite(slant) && slant > 0.0;
}

}  // namespace

double coneFrustumChordVolume(double R1, double R2, double height, int nSides) {
    double gamma = 0.0, segFac = 0.0, slant = 0.0;
    if (!frustumVolParams(R1, R2, height, nSides, gamma, segFac, slant)) return 0.0;
    // integral of rho^2 dz over the frustum = |h| * (R1^2 + R1 R2 + R2^2)/3.
    const double iRho2 = std::fabs(height) * (R1 * R1 + R1 * R2 + R2 * R2) / 3.0;
    const double v = static_cast<double>(nSides) * 0.5 * segFac * iRho2;
    return std::isfinite(v) ? v : 0.0;
}

double coneFrustumDVolDensity(double rho, double R1, double R2, double height, int nSides) {
    double gamma = 0.0, segFac = 0.0, slant = 0.0;
    if (!frustumVolParams(R1, R2, height, nSides, gamma, segFac, slant)) return 0.0;
    if (!std::isfinite(rho)) return 0.0;
    const double s = std::sin(0.5 * gamma);
    if (!(std::fabs(s) > 1e-15)) return 0.0;
    // The facet skin carries area w(z) * slant per unit of the axial parameter
    // while the gap carries A_seg(z) per unit of z, so the density that makes
    // the two integrals agree is A_seg / (w * slant / |h|)
    //   = rho * segFac * |h| / (4 * sin(gamma/2) * slant).
    const double k = rho * segFac * std::fabs(height) / (4.0 * s * slant);
    return std::isfinite(k) ? k : 0.0;
}

}  // namespace refit
}  // namespace stl2step
