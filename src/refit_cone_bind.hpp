// stl2step cone BIND site — the certified `circle-on-cone` edge class and the
// cone pcurves it certifies.  Lane 130-CONE-BUILD, D-130-3 lane 2 of 2.
//
// 130-CONE-MATH (`refit_cone_math.hpp`) publishes `circleOnConeMax`, the
// supremum distance from a circle to the untrimmed conical SURFACE.  That is
// not the quantity a bind site records.  `exactMaxAtBind` records
//
//     sup_t | C3d(t) - S(pcurve(t)) |
//
// which also carries the PARAMETRISATION: a rim displaced radially by delta
// sits `delta*cos(Ang)` from the surface but `delta` from the point its own
// pcurve names.  Omitting that factor would record a tolerance smaller than
// the gap the edge actually has.  This header is that composition, and it is
// closed form throughout — nothing here samples, iterates or approximates.
//
// D-130-2 in force: an edge ships analytic only where the program computes the
// supremum in closed form.  So each entry point either CERTIFIES (returns the
// exact supremum and an exact class name) or REFUSES (-1.0 and an `unhandled-`
// class), and refusal means tier 2 — mesh polyline, counted — never a widened
// tolerance.  The certificates are stated at each declaration; they are
// decidable from the curve, the pcurve and the surface alone.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_CONE_BIND_HPP
#define STL2STEP_REFIT_CONE_BIND_HPP

#include <Geom2d_Curve.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Standard_Handle.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Lin.hxx>

namespace stl2step {
namespace refit {

// The conical surface under `srf`, seeing through a Geom_RectangularTrimmedSurface.
// Null when `srf` is not a cone.
Handle(Geom_ConicalSurface) coneBasisOf(const Handle(Geom_Surface)& srf);

// ---------------------------------------------------------------------------
// pcurves on a cone
// ---------------------------------------------------------------------------
//
// Both are the cone analogue of the cylinder branch of `makePCurveOnSurf`, and
// both exploit the fact that OCCT's cone parametrisation is UNIT SPEED in v
// (dS/dv = sin(Ang)*rho_hat + cos(Ang)*Z is a unit vector) and ANGULAR in u, so
// a straight 2d line reproduces a rim circle and a generator line exactly.

// Rim circle -> Geom2d_Line (phi, v0) + t*(sgn, 0):  u runs with the circle's
// own parameter, v is pinned to the circle's axial height via `coneVAtZ`, and
// `sgn` matches the circle's rotation sense to the cone frame's.  Null unless
// the circle's axis is parallel to the cone's.
Handle(Geom2d_Curve) conePCurveForCircle(const gp_Cone& cone, const gp_Circ& circ);

// Generator line -> Geom2d_Line (u0, v0) + t*(0, sgn).  Null unless the line's
// direction is a generator direction of the cone (within `angTol`).
Handle(Geom2d_Curve) conePCurveForLine(const gp_Cone& cone, const gp_Lin& lin, double tMid,
                                       double angTol = 1e-9);

// ---------------------------------------------------------------------------
// the certified bind-site supremum
// ---------------------------------------------------------------------------

// CIRCLE on cone.  Certificate (all four decidable in closed form):
//   (1) the pcurve is a 2d line whose direction has |dv| <= Precision::PConfusion()
//       — v is constant, so the pcurve's image is a circle of the cone,
//   (2) the circle's axis is parallel to the cone's within Precision::Angular(),
//   (3) the circle's centre is on the cone axis within Precision::Confusion(),
//   (4) the two traversals share an angular velocity VECTOR: sign(du) * Zcone
//       equals the circle's own axis direction.
// Under (1)-(4) the two curves are concentric circles in parallel planes swept
// at one angular rate, so |C(t) - S(pc(t))| is CONSTANT in t and the value at
// either end IS the supremum — no sampling, and no assumption that the
// perpendicular distance and the parametric gap agree (they differ by exactly
// cos(Ang) for a radial displacement, which is the whole reason this function
// exists).  Class "circle-on-cone".
//
// Anything else — a tilted rim, an off-axis rim, a projected or reparametrised
// pcurve — returns -1.0 with class "unhandled-circle-on-cone".  That is tier 2
// by D-130-2 and must be COUNTED, never absorbed.
double coneCircleBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                         const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                         const TopLoc_Location& loc, const char** clsOut);

// LINE on cone (the seam generator).  Certificate: the pcurve is a 2d line whose
// direction has |du| <= Precision::PConfusion(), so u is constant and
// S(pc(t)) = O + (R + v(t)sinA)*rho_hat(u0) + v(t)cosA*Z is AFFINE in t.  C3d is
// affine too, so their difference is affine and |difference| is convex: the
// supremum is at an endpoint, exactly.  Class "line-on-cone"; otherwise -1.0
// and "unhandled-line-on-cone".
double coneLineBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                       const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                       const TopLoc_Location& loc, const char** clsOut);

// The single entry `exactMaxAtBind` calls once it has established that `srf` is
// a cone: dispatches on the 3d curve type and returns -1.0 with
// "unhandled-curve-on-cone" for anything that is neither a circle nor a line.
double coneBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                   const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                   const TopLoc_Location& loc, const char** clsOut);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_CONE_BIND_HPP
