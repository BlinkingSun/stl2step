// stl2step ELLIPSE BIND site — the certified `ellipse-on-plane` and
// `ellipse-on-cylinder[-bound]` edge classes.  Lane `ellipse-bind`, D-140-9 §3,
// SPEC-ellipse PART B.2.  The `refit_cone_bind.hpp` pattern, verbatim.
//
// `refit_ellipse_math.hpp` (PART A) publishes the SURFACE-side supremum and
// bound.  That is not the quantity a bind site records.  `exactMaxAtBind`
// records
//
//     sup_t | C3d(t) - S(pcurve(t)) |
//
// which carries the PARAMETRISATION as well as the geometry.  This TU is that
// composition: it reads the parametrisation the bind site actually stored --
// the 2-D ellipse on the plane side, the 2-D poles/knots on the cylinder side --
// hands the math TU plain arrays (so the math TU stays free of `Geom2d_`), and
// either CERTIFIES (the exact value or the certified bound, with an exact class
// name) or REFUSES (-1.0 with a named `unhandled-ellipse-*` class).  There is
// no third outcome, and a refusal is tier 2 by D-130-2 -- counted, never a
// widened tolerance.
//
// D-140-9 §3: `makePCurveOnSurf` is NOT edited by this lane.  The projector's
// output is already certifiable; rewriting the pcurve to make a bound easier is
// widening a tolerance to pass, and it is forbidden.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_ELLIPSE_BIND_HPP
#define STL2STEP_REFIT_ELLIPSE_BIND_HPP

#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Standard_Handle.hxx>
#include <TopLoc_Location.hxx>

namespace stl2step {
namespace refit {

// ELLIPSE on a PLANE.  Certificate, all decidable from the curve, the pcurve and
// the surface alone:
//   (P1) the 3-D basis curve is a `Geom_Ellipse` with a > Resolution,
//        b > Resolution and l - f > PConfusion,
//   (P2) `srf` is a `Geom_Plane`, so S(u, v) = O + u X + v Y is AFFINE and the
//        image of any 2-D point lifts exactly,
//   (P3) the stored pcurve is a `Geom2d_Ellipse`.
// Under (P1)-(P3) the difference is, IDENTICALLY IN t and with no assumption
// about `GeomAPI::To2d`'s convention,
//     E(t) - S(pc(t)) = W + A cos t + B sin t
// with W = C - P0, A = a U - a2 U2, B = b V - b2 V2 -- because BOTH sides are
// evaluated at the SAME t and both parametrisations are the gp_ cosine/sine
// form.  Parameter preservation is therefore MEASURED by the resulting
// supremum, not assumed: a To2d that reparametrised would make W, A, B large
// and the class would report that, honestly, as a large exact deviation.
// `ellipseOnPlaneMax` returns that supremum in closed form.  Class
// "ellipse-on-plane"; otherwise -1.0 and "unhandled-ellipse-plane-param"
// (P3 fails: the pcurve is not the 2-D ellipse the plane branch writes) or
// "unhandled-ellipse-plane-degenerate" (P1/P2 fail).
double ellipseBindSupOnPlane(const Handle(Geom_Curve)& c3, double f, double l,
                             const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                             const TopLoc_Location& loc, const char** clsOut);

// ELLIPSE on a CYLINDER.  Certificate:
//   (C1) the 3-D basis curve is a `Geom_Ellipse` and `srf` a cylinder,
//   (C2) the ellipse IS a section of that cylinder -- centre on the axis, and
//        a|Uperp| = b|Vperp| = R with Uperp.Vperp = 0 (`ellipseIsCylinderSection`,
//        refusal "unhandled-ellipse-cyl-oncyl"),
//   (C3) the stored pcurve is a NON-RATIONAL 2-D Bezier / B-spline in
//        piecewise-Bezier form (interior knot multiplicity = degree,
//        nPoles = degree*nSpans + 1) whose parameter range equals [f, l] within
//        `Precision::PConfusion()` -- THE RANGE CLAUSE, refusal
//        "unhandled-ellipse-cyl-span",
//   (C4) the pole azimuths track the section's affine theta within the derived
//        margin (refusal "unhandled-ellipse-cyl-uaffine"),
//   (C5) the section is not grazing (refusal "unhandled-ellipse-grazing").
// Under (C1)-(C5) `ellipseOnCylMax` returns the certified bound.  Class
// "ellipse-on-cylinder" when the residual is identically zero, otherwise
// "ellipse-on-cylinder-bound" -- a bound may be compared against, never
// presented as an exact deviation (`ellipseDevClassIsExact`).
double ellipseBindSupOnCyl(const Handle(Geom_Curve)& c3, double f, double l,
                           const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                           const TopLoc_Location& loc, const char** clsOut);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_ELLIPSE_BIND_HPP
