// stl2step torus math — the closed-form primitive behind the 1.4.0 toroidal
// round (lane torus-math, SPEC-torus PART B).  Decisions in force:
//   * D-130-2  an edge ships analytic only when the program computes, IN CLOSED
//              FORM, the supremum deviation between the shipped curve and the
//              exact intersection of the two shipped surfaces.  Everything in
//              this header is that closed form for the torus; nothing here
//              samples, iterates, or approximates.
//   * D-140-6  §3(4), the D-130-3 precedent applied verbatim: no torus face
//              ships until `circle-on-torus` is a certified `exactMaxAtBind`
//              class.  This header publishes that class and its supremum.
//
// Why its own translation unit (`src/refit_torus_math.cpp`) rather than the P1
// math TU: exactly the reason `refit_cone_math` has one (0d67eba).
// `src/refit_math.cpp` is one of the five P1 sources the D5.3 include allowlist
// covers, and that allowlist — a gate instrument, not a style rule — admits
// `refit.hpp` and `refit_internal.hpp` as the only project headers.  This math
// needs `gp_Torus.hxx` and nothing of the engine, so it lives beside the cone
// math, outside the allowlist by construction; the gate stays untouched.
// Nothing here touches `MeshView`, `Region` or any engine state, which is what
// lets `tests/unit/torus_math_test.cpp` compile this one TU and certify the
// math without linking — or waiting for — the engine.
//
// FRAME AND PARAMETRISATION (OCCT, `Geom_ToroidalSurface` / `gp_Torus`):
//
//   S(u, v) = O + (R + r*cos(v)) * (cos(u)*XDir + sin(u)*YDir) + r*sin(v)*ZDir
//
// with `R = torus.MajorRadius()` and `r = torus.MinorRadius()`.  Throughout:
//
//   z    axial coordinate of a point, (p - O) . ZDir        (h in SPEC-torus)
//   rho  radial coordinate of a point, its distance to the axis (>= 0)
//   u    azimuth about the axis; v tube angle, measured from the outer equator.
//
// The meridian half-plane {rho >= 0} meets the surface in the PROFILE CIRCLE
//
//   (rho - R)^2 + z^2 = r^2,
//
// and that circle lies entirely in rho >= 0 exactly when R >= r.  Everything in
// this header is stated for the RING torus, R > r, which is also the only torus
// OCCT's own `gp_Torus::SetMajorRadius` / `SetMinorRadius` guards accept
// (`R - r > gp::Resolution()`).  A spindle torus (R < r) is refused: its profile
// circle crosses the axis, the untrimmed surface there is the circle's MIRROR
// image, and the closed form below is then only an upper bound, not the
// supremum — it reads `2` on the spindle's own mirrored equator, a circle that
// lies ON the surface (the unit test exhibits that witness) — and D-130-2 lets
// only an exact class be recorded as a bind tolerance.  `R -> 0` is a sphere and
// a different class (D-140-6 §1(4): `Rmaj <= tau` routes to SphereNYI), never a
// torus.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_TORUS_MATH_HPP
#define STL2STEP_REFIT_TORUS_MATH_HPP

#include <cstdint>

#include <Precision.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>
#include <gp_Torus.hxx>

namespace stl2step {
namespace refit {

// ---------------------------------------------------------------------------
// 1. circle-on-torus deviation classes
// ---------------------------------------------------------------------------
//
// `circleOnTorusMax` reports which of these it computed.  There is one exact
// class and one refusal — D-140-6 §3(4) certifies COAXIAL circles only; a circle
// that is not coaxial (a general plane section of the torus) is `Unhandled`
// here, ships tier 2 (mesh polyline, counted) under D-130-2, and this lane does
// not widen that.
enum class TorusDevClass : std::uint8_t {
    Unhandled = 0,  // refused input; the returned value is -1.0
    Coaxial         // circle plane PERP axis, centre ON axis — EXACT, constant along the circle
};

// `exactMaxAtBind`-style class name, for the edge census D-130-2 requires:
//   Coaxial   -> "circle-on-torus"        (the class D-140-6 §3(4) gates torus faces on)
//   Unhandled -> "unhandled-circle-on-torus"
const char* torusDevClassName(TorusDevClass c);

// True only for Coaxial, i.e. only where the returned value IS the supremum.
bool torusDevClassIsExact(TorusDevClass c);

// ---------------------------------------------------------------------------
// 2. torus frame helpers
// ---------------------------------------------------------------------------

// The admission predicate every entry point below applies: finite radii,
// `r > gp::Resolution()` (a tube exists) and `R - r > gp::Resolution()` (a ring,
// OCCT's own guard).  False for a horn (R == r), a spindle (R < r), a sphere
// (R == 0) and a degenerate circle (r == 0).
bool torusIsRing(const gp_Torus& torus);

// z = (p - torus.Location()) . torus.Axis().Direction().
double torusAxialCoord(const gp_Torus& torus, const gp_Pnt& p);

// rho = distance from p to the torus axis (>= 0).
double torusRadialCoord(const gp_Torus& torus, const gp_Pnt& p);

// The tube angle of the profile point nearest to (rho, z): atan2(z, rho - R),
// in ]-PI, PI].  Inverse of `torusVIsoCircle` on the circle's own (radius,
// offset).  For the two tangency circles of a mouth round this is 0 or PI on
// the cylinder side (rho = R +- r, z = 0) and +-PI/2 on the plane side
// (rho = R, z = +-r).  NaN on a refused torus.
double torusVOfProfilePoint(const gp_Torus& torus, double rho, double z);

// The v-iso circle of the surface: radius `R + r*cos(v)` at axial offset
// `r*sin(v)`, oriented off the torus's own frame (Z = the axis, X = the torus
// XDirection) so repeated calls are byte-identical (I5).  These are the tier-1
// bounds a blend patch carries (D-140-6 §3(3)): for a mouth round, v = 0 or PI
// is the cylinder-side tangency circle and v = +-PI/2 the plane-side one.
// False — and writes nothing — on a refused torus.
bool torusVIsoCircle(const gp_Torus& torus, double v, gp_Circ& out);

// ---------------------------------------------------------------------------
// 3. exact point-to-surface distance
// ---------------------------------------------------------------------------

// Exact distance from p to the untrimmed ring torus:
//
//     d(p) = | hypot(rho_p - R, z_p) - r |
//
// PROOF.  The surface is invariant under rotation about the axis, so for any
// surface point Q with cylindrical coordinates (rho_Q, phi_Q, z_Q),
//   |p - Q|^2 = rho_p^2 + rho_Q^2 - 2 rho_p rho_Q cos(phi_Q - phi_p) + (z_p - z_Q)^2
//             >= (rho_p - rho_Q)^2 + (z_p - z_Q)^2,
// with equality at phi_Q = phi_p; so the nearest surface point lies in p's own
// meridian half-plane, and there the surface is the profile circle of centre
// (R, 0) and radius r, whose distance from the point (rho_p, z_p) is the
// modulus above.  It needs R >= r so that the whole profile circle is surface
// (see the header note); `torusIsRing` guarantees it.  Returns -1.0 on a
// refused torus or a non-finite point.
double pointTorusDist(const gp_Torus& torus, const gp_Pnt& p);

// ---------------------------------------------------------------------------
// 4. THE SUPREMUM — circle on torus  (SPEC-torus B.1, D-140-6 §3(4))
// ---------------------------------------------------------------------------

// Supremum over the circle's whole parameter range of the distance from the
// circle to the untrimmed toroidal surface, for a circle COAXIAL with the torus:
// plane perpendicular to the axis within `angTol`, centre on the axis within
// `linTol`.
//
//   For a coaxial circle of radius r_c at axial offset h (z of its centre)
//   every point has the same (rho, z) = (r_c, h), so by `pointTorusDist`
//   every point is at the SAME distance from the surface and the supremum is
//   that constant:
//
//       sup = | hypot(r_c - R, h) - r |                          (EXACT)
//
//   It is exactly 0 on every v-iso circle of the surface — in particular on
//   the two tangency circles a mouth round is bounded by: r_c = R +- r at
//   h = 0 (the cylinder side) and r_c = R at h = +-r (the plane side).  The
//   sign of (hypot - r) says whether the circle runs inside (-) or outside (+)
//   the tube; the modulus drops it, and the result is a distance, >= 0.
//
//   Degenerate circle inputs are NOT refused, because the identity holds for
//   them: r_c = 0 is a point on the axis at distance |hypot(R, h) - r|, and
//   h = 0 is the equatorial plane, where the value is ||r_c - R| - r| and the
//   zeros are the two equatorial circles r_c = R +- r.
//
//   Within the admission tolerances the value returned is the exact supremum
//   over the circle's true (rho, z) trajectory bounded as an interval box —
//   z in [h -+ r_c*sin(beta)], rho in [max(0, r_c*cos(beta) - d0, d0 - r_c),
//   d0 + r_c] (beta = axis/normal angle, d0 = centre-to-axis distance), and
//   hypot(rho - R, z) is convex so its extremes over the box are closed form —
//   which collapses to the constant above when beta = 0 and d0 = 0 and is
//   otherwise >= the truth (D-130-2: never a value below the supremum).
//
// Returns -1.0 with clsOut = Unhandled for: a torus `torusIsRing` refuses; a
// non-finite or negative circle radius; a circle whose plane is not
// perpendicular to the axis within `angTol` or whose centre is off the axis by
// more than `linTol` (non-coaxial — tier 2, counted, not this lane's to widen);
// a negative or non-finite tolerance.  clsOut may be null.
//
// The class is INVARIANT under rigid placement: a circle whose normal is the
// torus's own axis direction classes Coaxial for every orientation of that
// torus in world space.  sin(beta) is taken from the cross product, never from
// sqrt(1 - dot*dot) — see refit_cone_math for the measured reason.
double circleOnTorusMax(const gp_Torus& torus, const gp_Circ& circ,
                        TorusDevClass* clsOut = nullptr,
                        double angTol = Precision::Angular(),
                        double linTol = Precision::Confusion());

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_TORUS_MATH_HPP
