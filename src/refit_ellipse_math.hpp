// stl2step ellipse math — the closed forms behind the ELLIPSE bind class
// (lane ellipse-math, SPEC-ellipse PART A, D-140-9 §2; D-130-2, D-S3-111).
//
// The third instance of the 130-CONE-MATH / 140-torus-math shape: a pure gp_
// translation unit beside refit_cone_math (outside the P1 D5.3 include
// allowlist by construction — it needs both refit.hpp-class headers and bare
// gp_, exactly as refit_cone_math.hpp documents), no MeshView, no Region, no
// engine state, compiled and certified by tests/unit/ellipse_math_test.cpp
// WITHOUT linking the engine.  Nothing in src/refit_build.cpp calls it until
// the bind lane (PART B) does; `exactMaxAtBind` still returns `unhandled-other`
// for an ellipse when this TU lands, so B0 holds by construction.
//
// WHAT WAS MEASURED (SPEC-ellipse PART 0, do not re-derive): every one of the
// 55 `unhandled` edges on the 9-fixture edge-class red line is a plane∩cylinder
// ELLIPSE.  Its plane-side pcurve is the exact 2-D ELLIPSE `GeomAPI::To2d`
// writes (dev 1.3e-15); its cylinder-side pcurve is one degree-8 Bézier whose
// u-poles are in exact arithmetic progression and whose v-poles approximate a
// sinusoid (dev 4.05e-10, all of it axial).  Hence TWO classes:
//
//   ellipse-on-plane           EXACT supremum          (A.1)
//   ellipse-on-cylinder(-bound) certified UPPER BOUND  (A.2)
//
// and no exact cylinder class in general — no NURBS reproduces cos t in t, so
// the axial residual is a theorem, not an engine limitation.
//
// NOTATION.  E(t) = C + a cos t Û + b sin t V̂ (gp_Elips: Û = XAxis, V̂ = YAxis,
// a = MajorRadius, b = MinorRadius), t in [f, l].  A plane is (O; X̂, Ŷ, N̂) with
// S(x, y) = O + x X̂ + y Ŷ.  A cylinder is (O; X̂, Ŷ, d̂) of radius R with
// S(u, v) = O + R (cos u X̂ + sin u Ŷ) + v d̂; for a point, rho = distance to the
// axis, theta = atan2(·Ŷ, ·X̂), z = ·d̂.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_ELLIPSE_MATH_HPP
#define STL2STEP_REFIT_ELLIPSE_MATH_HPP

#include <cstdint>

#include <Precision.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace stl2step {
namespace refit {

// ---------------------------------------------------------------------------
// 1. classes
// ---------------------------------------------------------------------------
//
// The distinction is D-130-2's: an EXACT class may be recorded as the bind
// tolerance; a BOUND is >= the truth and may be compared against one, but is
// never presented as the deviation itself.
enum class EllipseDevClass : std::uint8_t {
    Unhandled = 0,   // refused; the returned value is -1.0
    OnPlane,         // 3-D ellipse vs its To2d image on a plane — EXACT supremum
    OnCylinder,      // exact: the pcurve reproduces the section identically
    OnCylinderBound  // PROVABLE UPPER BOUND on the projected pcurve
};

//   OnPlane -> "ellipse-on-plane" | OnCylinder -> "ellipse-on-cylinder"
//   OnCylinderBound -> "ellipse-on-cylinder-bound" | Unhandled -> "unhandled-ellipse"
const char* ellipseDevClassName(EllipseDevClass c);

// true for OnPlane and OnCylinder only
bool ellipseDevClassIsExact(EllipseDevClass c);

// Every refusal is NAMED (D-140-9 §3: "by named clause, never by silence").
// Both class functions take an optional trailing `clauseOut`; on refusal it
// receives one of the strings below, on admission the class name.  The strings
// are the census `clsA`/`clsB` values the bind TU forwards verbatim.
//   plane:    "unhandled-ellipse-plane-param"      lifted 2-D frame not in `pl`
//             "unhandled-ellipse-plane-degenerate" non-finite, a|b <= Resolution,
//                                                  l - f <= PConfusion
//   cylinder: "unhandled-ellipse-cyl-oncyl"        Step 1 fails (not a section)
//             "unhandled-ellipse-cyl-uaffine"      Step 2: eps_u exceeds the margin
//             "unhandled-ellipse-cyl-span"         nPoles/degree/knots/range clause
//             "unhandled-ellipse-cyl-degenerate"   non-finite, R|a|b <= Resolution
//   grazing:  "unhandled-ellipse-grazing"          the bind site's refusal when
//                                                  ellipseVertexBound > meshTolCap
extern const char* const kEllipseClausePlaneParam;
extern const char* const kEllipseClausePlaneDegenerate;
extern const char* const kEllipseClauseCylOnCyl;
extern const char* const kEllipseClauseCylUAffine;
extern const char* const kEllipseClauseCylSpan;
extern const char* const kEllipseClauseCylDegenerate;
extern const char* const kEllipseClauseGrazing;

// ---------------------------------------------------------------------------
// 2. the plane class — EXACT supremum  (SPEC-ellipse A.1)
// ---------------------------------------------------------------------------
//
//   sup_{t in [f, l]} |E(t) - S(pc(t))|,   S(pc(t)) = P0 + a2 cos t Û2 + b2 sin t V̂2
//
// where (P0; Û2, V̂2; a2, b2) is the 2-D ellipse `GeomAPI::To2d` wrote, LIFTED
// into 3-space through the plane's own frame: P0 = O + x0 X̂ + y0 Ŷ,
// Û2 = xu X̂ + yu Ŷ, V̂2 = xv X̂ + yv Ŷ.  (The unit test asserts against OCCT that
// To2d is parameter-preserving — the 2-D point at parameter t IS the image of
// the 3-D point at parameter t — on 13 placements; that is a certificate
// clause, not an assumption.  Lifting is the caller's, and the "plane-param"
// clause refuses a lifted frame that is not in `pl`: P0 off the plane by more
// than Precision::Confusion(), or Û2/V̂2 out of the plane or not orthogonal by
// more than Precision::Angular().)
//
// THE DERIVATION.  E(t) - S(pc(t)) = W + A cos t + B sin t with
//   W = C - P0,  A = a Û - a2 Û2,  B = b V̂ - b2 V̂2.
// g(t) = |W + A cos t + B sin t|^2 has
//   g'(t)/2 = (|B|^2 - |A|^2) sin t cos t + (A·B)(cos^2 t - sin^2 t)
//             - (W·A) sin t + (W·B) cos t  =: h(t),
// and under tau = tan(t/2), (1 + tau^2)^2 h(t) is the QUARTIC
//   (Q - Rb) tau^4 - 2 (P + Ra) tau^3 - 6 Q tau^2 + 2 (P - Ra) tau + (Q + Rb),
//   P = |B|^2 - |A|^2,  Q = A·B,  Ra = W·A,  Rb = W·B.
// A continuous function on a compact interval attains its supremum at an
// endpoint or at a critical point, so the supremum is EXACTLY
//   max( |·|(f), |·|(l), max over the quartic's real roots t* in (f, l) ).
// The quartic is solved by Ferrari (resolvent cubic by Cardano/trigonometric
// form); the substitution is taken about a shift t0 chosen — threshold-free —
// as the point of largest |h| among eight equispaced samples, so the leading
// coefficient h(t0 + pi) is never the one that vanishes, and every closed-form
// root is polished by three fixed Newton steps on h (deterministic; the
// supremum is second-order insensitive to the root's position, the polish only
// removes Ferrari's cancellation).  The lost point t0 + pi and the eight
// samples are added as candidates: extra candidates never change the maximum
// of a set that already contains every critical point and both endpoints.
//
// Branches (each unit-tested, SPEC A.1):
//   W = 0 and l - f >= pi  ->  ellipseTrigSupZeroOffset(A, B): |A cos t + B sin t|
//                              has period pi, so on any interval of length >= pi
//                              its supremum is sigma_max([A B]).  Bitwise the
//                              same value as calling that function directly.
//   A = B = W = 0          ->  exactly 0.0 (the SHIPPED case: the ellipse IS the
//                              plane section and To2d is parameter-preserving).
//   a = b (a circle)       ->  the value is <= circleOnPlaneMax(pl, circ) (the
//                              D-S3-111 triangle bound) and >= the sampled truth
//                              — the sanity relation between the two classes.
double ellipseOnPlaneMax(const gp_Pln& pl, const gp_Elips& el, double f, double l,
                         const gp_Pnt& p0, const gp_Dir& u2, const gp_Dir& v2,
                         double a2, double b2,
                         EllipseDevClass* clsOut = nullptr,
                         const char** clauseOut = nullptr);

// Degenerate branch, published separately so the bind site and the unit can
// both call it: W = 0  =>  sup_t |A cos t + B sin t| = the largest singular
// value of the 3x2 matrix [A B],
//   sigma_max = sqrt( (|A|^2 + |B|^2 + sqrt((|A|^2 - |B|^2)^2 + 4 (A·B)^2)) / 2 ).
// Exactly 0.0 for A = B = 0.  Also the residual norm the cylinder class uses
// for its section identity (Step 1 below), with 2-D vectors lifted to z = 0.
double ellipseTrigSupZeroOffset(const gp_Vec& A, const gp_Vec& B);

// ---------------------------------------------------------------------------
// 3. the cylinder section identity  (SPEC-ellipse A.2 Step 1)
// ---------------------------------------------------------------------------
//
// True iff E lies ON cyl identically: centre on the axis (<= linTol), and with
// Û⊥ = Û - (Û·d̂) d̂, V̂⊥ = V̂ - (V̂·d̂) d̂:  a|Û⊥| = b|V̂⊥| = R (<= linTol) and
// Û⊥ ⟂ V̂⊥ (cosine <= Precision::Angular(); ellipseOnCylMax passes its own
// angTol).  Under it, IDENTICALLY:
//   rho(t) == R,   theta(t) = theta0 + sigma t (sigma = +-1),
//   z(t) = z0 + a (Û·d̂) cos t + b (V̂·d̂) sin t = z0 + Lambda cos(t - psi),
//   Lambda = sqrt((a Û·d̂)^2 + (b V̂·d̂)^2).
// theta0 is the azimuth of Û⊥ in the cylinder's (X̂, Ŷ); sigma the sense of
// (Û⊥ x V̂⊥)·d̂; z0 = (C - O)·d̂; psi = atan2(b V̂·d̂, a Û·d̂).  No case split on
// which axis is major.
//
// Inside the tolerances the identity holds up to a residual the cylinder class
// does NOT ignore: with M = [a Û⊥ | b V̂⊥] in the cylinder's (X̂, Ŷ) and
// M0 = R [cos theta0, -sigma sin theta0; sin theta0, sigma cos theta0], the
// perturbation Delta = M - M0 satisfies sup_t |Delta (cos t, sin t)| =
// sigma_max(Delta) (ellipseTrigSupZeroOffset), so |rho - R| <= sigma_max and
// |theta - (theta0 + sigma t)| <= asin(sigma_max / R).  Both enter the bound
// below; both are exactly 0 on an exact section.
bool ellipseIsCylinderSection(const gp_Cylinder& cyl, const gp_Elips& el,
                              double& theta0, double& sigma, double& z0,
                              double& lambda, double& psi,
                              double linTol = Precision::Confusion());

// ---------------------------------------------------------------------------
// 4. the grazing predicate  (SPEC-ellipse A.3, D-S3-111's open predicate)
// ---------------------------------------------------------------------------
//
// gamma = the angle between the plane normal N̂ and the cylinder normal rhô,
// sin gamma = |N̂ x rhô|.  Along the section it is smallest at the ellipse's
// major vertices, where rhô is the in-plane component of N̂'s direction, and
// there sin gamma = |N̂·d̂| = cos alpha = b/a (alpha = the cut's tilt from the
// axis).  A plane perpendicular to the axis (a circle) returns 1; a plane
// tangent to the cylinder returns 0.  Returns -1.0 on non-finite input.
double ellipseCylGrazingSin(const gp_Cylinder& cyl, const gp_Pln& pl);

// A mesh vertex P within deltaPlane of the shipped plane and deltaCyl of the
// shipped cylinder (both <= tau = 2q by the region certificate).  Writing
// P = Q + s m̂ + h N̂ with Q on the ellipse and m̂ the in-plane normal of the
// ellipse, rho(P) - R = s sin gamma + h cos gamma to first order, so
//   |h| <= deltaPlane,
//   |s| <= (deltaCyl + deltaPlane |cos gamma|) / sin gamma,
//   d(P, E) <= sqrt(s^2 + h^2) + s^2 / (2 rhoMinCurv),
// rhoMinCurv = b^2/a = R cos alpha the ellipse's tightest radius of curvature.
// Returns exactly that expression; +inf when sin gamma <= 0 (a tangential cut
// is not an ellipse the program can certify); -1.0 on non-finite or negative
// input.  THE PREDICATE: the class refuses (`unhandled-ellipse-grazing`) when
// this bound exceeds the edge's meshTolCap — a margin derived from q, R, alpha
// and meshTolCap, never a hand angle.
double ellipseVertexBound(double deltaPlane, double deltaCyl, double sinGamma,
                          double cosGamma, double rhoMinCurv);

// ---------------------------------------------------------------------------
// 5. the cylinder class — certified UPPER BOUND  (SPEC-ellipse A.2)
// ---------------------------------------------------------------------------
//
// Certified bound on sup_t |E(t) - S(pc(t))| for a 2-D B-spline / Bézier
// pcurve pc(t) = (u(t), v(t)) on the cylinder, parametrised by the SAME t as E
// (the BRep SameParameter contract).  Poles / knots are passed as plain arrays
// so this TU stays free of Geom2d_: `degree`, `nPoles` poles (poleU[i],
// poleV[i]), `nKnots` distinct knots with multiplicities `mults` (a Bézier is
// nKnots = 2, mults = {degree+1, degree+1}).
//
// Step 1  the section identity (above)                       -> "-oncyl"
// Step 2  the azimuthal residual, EXACT from the poles.  The affine function
//         theta0 + sigma t has B-spline coefficients theta0 + sigma xi_i at the
//         Greville abscissae xi_i (for a Bézier, f + (i/n)(l - f)), so
//           eps_u = max_i |u_i - (theta0 + sigma xi_i)|
//         and by the convex-hull property (basis >= 0, sums to 1)
//           |u(t) - (theta0 + sigma t)| <= eps_u   for all t.
//         theta0 is taken modulo 2 pi to the branch the poles use (S is
//         2 pi-periodic in u).  eps_u is carried EXACTLY into Step 4, never
//         dropped and never gated by an angular tolerance (measured on
//         GeomProjLib::Curve2d output it ranges 1e-15 .. 1e-9 with the arc);
//         refused ("-uaffine") only at the derived margin eps_u + eta >= pi,
//         beyond which the azimuthal term stops being monotone in eps_u.
// Step 3  the axial residual, per knot span, CLOSED FORM.  On a span of width
//         H, v(t) is one polynomial of degree n and z(t) = z0 + Lambda cos(t -
//         psi) has every derivative bounded by Lambda.  With x_j the n+1
//         Chebyshev nodes of the span and p_n the interpolant of z there,
//           z - v = (z - p_n) + (p_n - v),
//           |z - p_n|   <= Lambda H^(n+1) / (2^(2n+1) (n+1)!)      (classical:
//                          max|omega| = 2 (H/4)^(n+1) for Chebyshev nodes),
//           |p_n - v|   <= L_n max_j |v(x_j) - z(x_j)|              (p_n - v is a
//                          degree-n polynomial interpolating the residuals,
//                          L_n <= (2/pi) ln(n+1) + 1 its Lebesgue constant).
//         The node set is finite and fixed by the pcurve's own degree — the
//         "certified-sample" door D-130-2(b) names.  A midpoint-Taylor
//         remainder is REFUSED (2^n looser; measured 256x on S11-b).
// Step 4  two points on one cylinder, rho1 = R + d:  |P1 - P2|^2 =
//         4 R (R + d) sin^2(dTheta/2) + d^2 + dz^2 — exact, no approximation —
//         so with the Step-1 residual sigma_max (0 on an exact section)
//           dev <= sqrt( 4 R (R + sigma_max) sin^2((eps_u + asin(sigma_max/R))/2)
//                        + sigma_max^2 + (max over spans |dz|)^2 ).
//
// Classes: OnCylinder (exactly 0.0) when every residual term is exactly zero —
// the pcurve reproduces the section identically; OnCylinderBound otherwise.
// Refusals, each named: "-oncyl", "-uaffine", "-span" (nPoles != degree + 1 +
// (sum of interior multiplicities), knots not increasing, mults out of range,
// or the knot range != [f, l] within Precision::PConfusion() — the range
// clause is D-140-9 §6's and is load-bearing), "-degenerate".
double ellipseOnCylMax(const gp_Cylinder& cyl, const gp_Elips& el, double f, double l,
                       int degree, int nPoles, const double* poleU, const double* poleV,
                       int nKnots, const double* knots, const int* mults,
                       EllipseDevClass* clsOut = nullptr,
                       double angTol = Precision::Angular(),
                       double linTol = Precision::Confusion(),
                       const char** clauseOut = nullptr);

// The Step-3 remainder alone, published so the unit can certify its
// monotonicity in n and H (halving H divides it by 2^(n+1)):
//   Lambda H^(n+1) / (2^(2n+1) (n+1)!),  0 <= n <= 25.  -1.0 outside.
double ellipseChebyshevRemainder(double lambda, double H, int degree);

// The Lebesgue-constant bound L_n <= (2/pi) ln(n+1) + 1 used in Step 3.
double ellipseChebyshevLebesgue(int degree);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_ELLIPSE_MATH_HPP
