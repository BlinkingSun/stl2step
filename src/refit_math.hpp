// stl2step cone math — closed-form primitives for the 1.3.0 chamfer frustum.
//
// SPEC-130-cone-math, lane 130-CONE-MATH.  Decisions in force:
//   * D-130-2  an edge ships analytic only when the program computes, IN CLOSED
//              FORM, the supremum deviation between the shipped curve and the
//              exact intersection of the two shipped surfaces.  Everything in
//              this header is that closed form for the cone; nothing here
//              samples, iterates, or approximates.
//   * D-130-3  chamfer frustums first.  No cone face ships until
//              `circle-on-cone` is a certified `exactMaxAtBind` class — this
//              header publishes that class and its supremum.
//
// Why a second header for `refit_math.cpp` (whose cross-TU signatures are frozen
// in `refit_internal.hpp`): every declaration below is pure `gp_` geometry with
// no `MeshView`, no `Region`, no engine state, so `tests/unit/cone_math_test.cpp`
// can compile `src/refit_math.cpp` on its own and certify this math without
// linking — or waiting for — the engine.  Nothing in `refit_internal.hpp` moves.
//
// FRAME AND PARAMETRISATION (OCCT, `Geom_ConicalSurface` / `gp_Cone`):
//
//   S(u, v) = O + (R + v*sin(Ang)) * (cos(u)*XDir + sin(u)*YDir) + v*cos(Ang)*ZDir
//
// with `R = cone.RefRadius() >= 0`, `Ang = cone.SemiAngle()` in ]-PI/2, PI/2[
// and NON-ZERO (`Ang` may be negative: the apex is then on the positive side of
// ZDir).  Throughout this header:
//
//   z    axial coordinate of a point, (p - O) . ZDir
//   rho  radial coordinate of a point, its distance to the cone axis (>= 0)
//   v    the surface parameter above; v = z / cos(Ang), so v and z are affine
//        images of one another and the apex is excluded from both or neither.
//
// The SIGNED cone radius at height z is `rho_cone(z) = R + z*tan(Ang)`, which is
// what `coneRadiusAt` returns.  It is negative past the apex, where the surface
// point sits on the opposite meridian (u + PI) — the untrimmed OCCT surface is
// the complete DOUBLE cone, and every formula here is stated for that surface.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_MATH_HPP
#define STL2STEP_REFIT_MATH_HPP

#include <cstdint>

#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

namespace stl2step {
namespace refit {

// ---------------------------------------------------------------------------
// 1. circle-on-cone deviation classes
// ---------------------------------------------------------------------------
//
// `circleOnConeMax` reports which of these it computed.  The distinction is the
// whole point of D-130-2: an EXACT class may be recorded as the bind tolerance;
// a BOUND may only be compared against one (and, being >= the truth, is safe to
// compare), and it is never presented as the deviation itself.
enum class ConeDevClass : std::uint8_t {
    Unhandled = 0,  // degenerate input; the returned value is -1.0
    CoaxialRim,     // circle plane PERP axis, centre ON axis  — EXACT, constant in u
    PerpPlane,      // circle plane PERP axis, centre OFF axis — EXACT supremum
    GeneralBound    // tilted and/or off-axis                  — PROVABLE UPPER BOUND
};

// `exactMaxAtBind`-style class name, for the edge census D-130-2 requires:
//   CoaxialRim   -> "circle-on-cone"        (the class D-130-3 gates cone faces on)
//   PerpPlane    -> "circle-on-cone-perp"
//   GeneralBound -> "circle-on-cone-bound"
//   Unhandled    -> "unhandled-circle-on-cone"
const char* coneDevClassName(ConeDevClass c);

// True only for CoaxialRim and PerpPlane, i.e. only where the returned value IS
// the supremum rather than an upper bound on it.
bool coneDevClassIsExact(ConeDevClass c);

// ---------------------------------------------------------------------------
// 2. cone frame helpers
// ---------------------------------------------------------------------------

// z = (p - cone.Location()) . cone.Axis().Direction().
double coneAxialCoord(const gp_Cone& cone, const gp_Pnt& p);

// rho = distance from p to the cone axis (>= 0).
double coneRadialCoord(const gp_Cone& cone, const gp_Pnt& p);

// ---------------------------------------------------------------------------
// 3. frustum helpers  (SPEC item 3 — CONE-BUILD excludes the apex with these)
// ---------------------------------------------------------------------------

// The five scalar helpers below return NaN — never a plausible-looking 0 or
// inf — for a cone that could not have come from a legal gp_Cone (non-finite
// parameters, negative reference radius, semi-angle outside ]-PI/2, PI/2[).

// SIGNED cone radius at axial height z:  rho_cone(z) = R + z*tan(Ang).
// Zero exactly at the apex, negative beyond it.
double coneRadiusAt(const gp_Cone& cone, double z);

// Axial coordinate of the apex: z_apex = -R / tan(Ang).  Agrees with
// gp_Cone::Apex() by construction.
double coneApexZ(const gp_Cone& cone);

// Surface parameter of the apex: v_apex = -R / sin(Ang) = coneVAtZ(coneApexZ()).
double coneApexV(const gp_Cone& cone);

// v <-> z and v <-> rho, both exact inverses of the parametrisation above.
double coneVAtZ(const gp_Cone& cone, double z);
double coneVAtRadius(const gp_Cone& cone, double rho);

// v-range of the frustum whose two rims sit at axial heights zA and zB.
// Returns false — and writes nothing — unless the frustum is well posed:
// zA != zB and the SIGNED cone radius is strictly positive at both heights.
// That condition is exactly apex exclusion: rho_cone(z) = tan(Ang)*(z - z_apex),
// so two positive radii put both rims on the same side of the apex, and
// [vLo, vHi] therefore never contains v_apex.  vLo < vHi on success.
bool coneFrustumVRange(const gp_Cone& cone, double zA, double zB, double& vLo, double& vHi);

// Same, keyed on the two rim RADII (what detector C measures) instead of their
// heights.  Requires rhoA > 0, rhoB > 0, rhoA != rhoB.
bool coneFrustumVRangeFromRadii(const gp_Cone& cone, double rhoA, double rhoB, double& vLo,
                                double& vHi);

// Build the cone carrying a frustum given as: rim A of radius rhoA lying in the
// XY plane of `rimA` (whose Location() is on the axis and whose Direction() is
// the axis), and rim B of radius rhoB at axial height `height` from it.
// Ang = atan((rhoB - rhoA) / height), RefRadius = rhoA, Position = rimA.
// Returns false when the frustum is degenerate (non-finite input, height ~ 0,
// rhoA ~ rhoB, or a non-positive rim radius) — i.e. exactly when gp_Cone would
// raise Standard_ConstructionError.
bool coneFromRims(const gp_Ax3& rimA, double rhoA, double rhoB, double height, gp_Cone& out);

// ---------------------------------------------------------------------------
// 4. exact point-to-surface distance
// ---------------------------------------------------------------------------

// Exact distance from p to the untrimmed (double-nappe) conical surface:
//
//     d(p) = cos(Ang) * | rho_p - |rho_cone(z_p)| |
//
// Closed form, no clamping: the meridian plane through p meets the double cone
// in two full lines, both entirely on the surface, so the perpendicular foot is
// always a real surface point and the apex is never the nearest point except in
// the limit.  Returns -1.0 on degenerate input.
double pointConeDist(const gp_Cone& cone, const gp_Pnt& p);

// ---------------------------------------------------------------------------
// 5. THE SUPREMUM — circle on cone  (SPEC item 1)
// ---------------------------------------------------------------------------

// Supremum over the circle's whole parameter range of the distance from the
// circle to the untrimmed conical surface.
//
//   (a) circle plane PERP to the axis (CoaxialRim / PerpPlane): EXACT.
//       Every point of the circle has the same z, so with
//       xi = |rho_cone(z)| and the circle's radial span [rho_lo, rho_hi]
//       (exactly [|d0 - r|, d0 + r] for a centre at axis distance d0):
//
//           sup = cos(Ang) * max(rho_hi - xi, xi - rho_lo)
//
//       which for a centre ON the axis collapses to the constant
//       cos(Ang) * |r - |rho_cone(z)||.
//
//   (b) any other circle (tilted and/or off-axis): a PROVABLE UPPER BOUND,
//       never an approximation.  Interval evaluation of the same exact formula:
//       z is confined to [z0 -+ r*sin(beta)] (tight, beta = axis/normal angle)
//       and rho to [max(0, r*cos(beta) - d0, d0 - r), d0 + r] (the circle
//       projects onto the plane perpendicular to the axis as an ellipse of
//       semi-axes r and r*cos(beta), so the projected offset has modulus in
//       [r*cos(beta), r] and the triangle inequality gives both ends), and the
//       maximum of |rho - xi| over that box is taken at a corner.  The
//       perpendicular case is the same expression with sin(beta) = 0 and
//       cos(beta) = 1, which is why one evaluator serves both regimes.
//
// Returns -1.0 with clsOut = Unhandled only for degenerate input (non-finite
// values, negative radii, a semi-angle outside ]-PI/2, PI/2[ \ {0}).
// `angTol` decides "plane perpendicular to the axis"; `linTol` decides
// "centre on the axis".  clsOut may be null.
//
// The class is INVARIANT under rigid placement: a circle whose normal is the
// cone's own axis direction classes CoaxialRim/PerpPlane for every orientation
// of that cone in world space, never GeneralBound.  (That is not automatic —
// sin(beta) must come from the cross product; see the note at the definition.)
double circleOnConeMax(const gp_Cone& cone, const gp_Circ& circ, ConeDevClass* clsOut = nullptr,
                       double angTol = Precision::Angular(),
                       double linTol = Precision::Confusion());

// ---------------------------------------------------------------------------
// 6. tier-1 IntAna cone cases  (SPEC item 2, D-130-2 tier 1)
// ---------------------------------------------------------------------------
//
// Each returns the EXACT gp_Circ of the intersection, or false ("none").  The
// supremum deviation of these curves from the true intersection is 0 by
// construction, which is what makes them tier 1.  Non-coaxial and oblique
// configurations return false and fall to tier 2 (mesh polyline, counted) —
// never to an approximated conic.
//
// All three produce a circle whose gp_Ax2 takes its Z from the cone axis and
// its X from the cone's own XDirection, so repeated calls are byte-identical
// (I5 determinism).

// Cone with a plane PERPENDICULAR to the cone axis.  false when the plane is
// not perpendicular within `angTol`, or passes through the apex (a point, not a
// circle) within `linTol`.
bool coneIntPlaneCircle(const gp_Cone& cone, const gp_Pln& pln, gp_Circ& out,
                        double angTol = Precision::Angular(),
                        double linTol = Precision::Confusion());

// Cone with a COAXIAL cylinder: the circle at the height where the cone radius
// equals the cylinder radius, z = (R_cyl - R) / tan(Ang).  false unless the
// axes are parallel within `angTol` AND coincident within `linTol`, and
// R_cyl > linTol.
//
// The double cone also meets the cylinder on its far nappe, at
// z' = (-R_cyl - R) / tan(Ang); that root is below the apex and therefore
// outside the v-range of any frustum (whose rims both have positive radius), so
// it is never the edge a shipped chamfer face needs.  This function returns the
// near-nappe root only.
bool coneIntCylCircle(const gp_Cone& cone, const gp_Cylinder& cyl, gp_Circ& out,
                      double angTol = Precision::Angular(),
                      double linTol = Precision::Confusion());

// Two COAXIAL cones: the circle at the height where their signed radii agree,
// after cone `b` is re-expressed in cone `a`'s frame.  false unless the axes
// are parallel within `angTol` and coincident within `linTol`, the generators
// are not parallel, or the crossing radius is not strictly positive.  The
// branch solved is rho_a(z) == rho_b(z) with both SIGNED radii positive; the
// two double cones also meet where the signed radii are equal and NEGATIVE (a
// real circle, but on the far nappe of both) and where rho_a == -rho_b (one
// nappe each).  Neither can lie in the v-range of a frustum, whose rims are
// both positive, so neither is returned.
bool coneIntConeCircle(const gp_Cone& a, const gp_Cone& b, gp_Circ& out,
                       double angTol = Precision::Angular(),
                       double linTol = Precision::Confusion());

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_MATH_HPP
