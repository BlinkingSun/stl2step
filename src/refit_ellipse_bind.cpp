// stl2step ELLIPSE BIND site — see refit_ellipse_bind.hpp for the certificates.
//
// SPDX-License-Identifier: MIT

#include "refit_ellipse_bind.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <gp.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax22d.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include "refit_ellipse_math.hpp"

namespace stl2step {
namespace refit {
namespace {

const double kTwoPi = 8.0 * std::atan(1.0);

// The SEAM sub-clause of the range clause (D-130-16, restated by SPEC-ellipse
// B.4).  A full-turn ellipse is closed at ONE vertex by `makeFullEllipse`
// (refit_build.cpp:11194), and D-130-16 requires a 360 degree rim to collapse
// into TWO arcs of the same ellipse meeting at the chain's own mesh vertex, so
// that both faces share the edges and no vertex is invented.  B.4: "state the
// clause, do not build for it."  This lane therefore does not certify a single
// full-turn ellipse edge -- it refuses it, by the same named clause the range
// carries, and the census counts it.
//
// B.4 recorded "measured: none of the 55 is a full ellipse", measured on S11-b.
// SPEC-ellipse B.5's own census (this lane's first act) measured otherwise on
// S16-R1-explode-success: `#82 = EDGE_CURVE('',#60,#60,#83,.T.)` closes on one
// vertex over a 2 pi range, and its cylinder-side pcurve is a degree-1 segment
// at CONSTANT v -- a circle at one height, not the section.  B.5 rules that
// case: "the certificate's clauses -- not the bound's constants -- are what
// must widen, and that is a report, not a tune."  Nothing numeric moves here.
bool isFullTurn(double f, double l) { return (l - f) >= kTwoPi - Precision::PConfusion(); }

Handle(Geom_Curve) basis3d(const Handle(Geom_Curve)& c) {
    Handle(Geom_TrimmedCurve) t = Handle(Geom_TrimmedCurve)::DownCast(c);
    if (!t.IsNull() && !t->BasisCurve().IsNull()) return t->BasisCurve();
    return c;
}

Handle(Geom2d_Curve) basis2d(const Handle(Geom2d_Curve)& c) {
    Handle(Geom2d_TrimmedCurve) t = Handle(Geom2d_TrimmedCurve)::DownCast(c);
    if (!t.IsNull() && !t->BasisCurve().IsNull()) return t->BasisCurve();
    return c;
}

Handle(Geom_Surface) basisSurf(const Handle(Geom_Surface)& s) {
    Handle(Geom_RectangularTrimmedSurface) t =
        Handle(Geom_RectangularTrimmedSurface)::DownCast(s);
    if (!t.IsNull() && !t->BasisSurface().IsNull()) return t->BasisSurface();
    return s;
}

// The world ellipse the bind site actually records. `c3` is already in world
// space at this site (`curveSurfDevAtBind` evaluates it without `loc`), so no
// transform is applied here -- the surface is the side that carries `loc`.
bool worldEllipse(const Handle(Geom_Curve)& c3, gp_Elips& out) {
    Handle(Geom_Ellipse) ge = Handle(Geom_Ellipse)::DownCast(basis3d(c3));
    if (ge.IsNull()) return false;
    out = ge->Elips();
    return true;
}

}  // namespace

double ellipseBindSupOnPlane(const Handle(Geom_Curve)& c3, double f, double l,
                             const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                             const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled-ellipse-plane-degenerate";
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return -1.0;
    Handle(Geom_Plane) gp_ = Handle(Geom_Plane)::DownCast(basisSurf(srf));
    if (gp_.IsNull()) return -1.0;
    gp_Elips el;
    if (!worldEllipse(c3, el)) return -1.0;
    if (!(l - f > Precision::PConfusion())) return -1.0;
    if (isFullTurn(f, l)) {
        if (clsOut) *clsOut = "unhandled-ellipse-plane-param";
        return -1.0;
    }

    // (P3) the stored pcurve is the 2-D ellipse the plane branch of
    // `makePCurveOnSurf` writes through `GeomAPI::To2d`.
    Handle(Geom2d_Ellipse) e2 = Handle(Geom2d_Ellipse)::DownCast(basis2d(c2d));
    if (e2.IsNull()) {
        if (clsOut) *clsOut = "unhandled-ellipse-plane-param";
        return -1.0;
    }

    // (P2) lift the 2-D ellipse's frame through the plane's own AFFINE map,
    // S(u, v) = O + u X + v Y, with the face's location carried on the surface
    // side exactly as `curveSurfDevAtBind` does.
    gp_Pln pln = gp_->Pln();
    if (!loc.IsIdentity()) pln.Transform(loc.Transformation());
    const gp_Ax3 pos = pln.Position();
    const gp_Vec X(pos.XDirection()), Y(pos.YDirection());
    const gp_Ax22d f2 = e2->Position();
    const gp_Pnt2d c2 = f2.Location();
    const gp_Dir2d xd = f2.XDirection(), yd = f2.YDirection();
    const gp_Pnt p0 = pos.Location().Translated(X * c2.X() + Y * c2.Y());
    const gp_Vec u2v = X * xd.X() + Y * xd.Y();
    const gp_Vec v2v = X * yd.X() + Y * yd.Y();
    if (!(u2v.Magnitude() > gp::Resolution()) || !(v2v.Magnitude() > gp::Resolution()))
        return -1.0;

    EllipseDevClass cls = EllipseDevClass::Unhandled;
    const double m = ellipseOnPlaneMax(pln, el, f, l, p0, gp_Dir(u2v), gp_Dir(v2v),
                                       e2->MajorRadius(), e2->MinorRadius(), &cls);
    if (!(m >= 0.0) || cls != EllipseDevClass::OnPlane) return -1.0;
    if (clsOut) *clsOut = ellipseDevClassName(cls);
    return m;
}

double ellipseBindSupOnCyl(const Handle(Geom_Curve)& c3, double f, double l,
                           const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                           const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled-ellipse-cyl-degenerate";
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return -1.0;
    Handle(Geom_CylindricalSurface) gcy =
        Handle(Geom_CylindricalSurface)::DownCast(basisSurf(srf));
    if (gcy.IsNull()) return -1.0;
    gp_Elips el;
    if (!worldEllipse(c3, el)) return -1.0;
    if (!(l - f > Precision::PConfusion())) return -1.0;
    gp_Cylinder cyl = gcy->Cylinder();
    if (!loc.IsIdentity()) cyl.Transform(loc.Transformation());

    // (C5) the grazing predicate (D-S3-111's open predicate, SPEC-ellipse A.3).
    // The section plane is the ellipse's own plane, so the published
    // `ellipseCylGrazingSin` is called on the pair the edge actually carries.
    // A tangential cut is not an ellipse this program can certify.
    const gp_Pln secPln(el.Location(), el.Position().Direction());
    const double sinG = ellipseCylGrazingSin(cyl, secPln);
    if (!(sinG > gp::Resolution())) {
        if (clsOut) *clsOut = "unhandled-ellipse-grazing";
        return -1.0;
    }

    // (C3) the stored pcurve, as plain arrays: non-rational, piecewise-Bezier,
    // and the RANGE CLAUSE. `ellipseOnCylMax` re-checks the structure from the
    // arrays; this side is only the extraction.
    int degree = 0, nPoles = 0, nKnots = 0;
    std::vector<double> pu, pv, kn;
    std::vector<int> mu;
    const Handle(Geom2d_Curve) b2 = basis2d(c2d);
    Handle(Geom2d_BSplineCurve) bs = Handle(Geom2d_BSplineCurve)::DownCast(b2);
    Handle(Geom2d_BezierCurve) bz = Handle(Geom2d_BezierCurve)::DownCast(b2);
    if (!bs.IsNull()) {
        if (bs->IsRational()) {
            if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
            return -1.0;
        }
        degree = bs->Degree();
        nPoles = bs->NbPoles();
        nKnots = bs->NbKnots();
        if (degree < 1 || nPoles < 2 || nKnots < 2) {
            if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
            return -1.0;
        }
        pu.resize(nPoles);
        pv.resize(nPoles);
        for (int i = 1; i <= nPoles; i++) {
            const gp_Pnt2d P = bs->Pole(i);
            pu[i - 1] = P.X();
            pv[i - 1] = P.Y();
        }
        kn.resize(nKnots);
        mu.resize(nKnots);
        for (int i = 1; i <= nKnots; i++) {
            kn[i - 1] = bs->Knot(i);
            mu[i - 1] = bs->Multiplicity(i);
        }
    } else if (!bz.IsNull()) {
        if (bz->IsRational()) {
            if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
            return -1.0;
        }
        degree = bz->Degree();
        nPoles = bz->NbPoles();
        if (degree < 1 || nPoles != degree + 1) {
            if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
            return -1.0;
        }
        pu.resize(nPoles);
        pv.resize(nPoles);
        for (int i = 1; i <= nPoles; i++) {
            const gp_Pnt2d P = bz->Pole(i);
            pu[i - 1] = P.X();
            pv[i - 1] = P.Y();
        }
        // A Bezier's own range is [0, 1]; the range clause below decides.
        nKnots = 2;
        kn.assign({bz->FirstParameter(), bz->LastParameter()});
        mu.assign({degree + 1, degree + 1});
    } else {
        if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
        return -1.0;
    }

    // The range clause, stated at the bind site so the refusal it produces is
    // named `-span` rather than hidden inside the math TU's return. Its seam
    // sub-clause (`isFullTurn`, above) is part of the same clause and carries
    // the same name: no new census class string enters the vocabulary.
    if (isFullTurn(f, l) || std::fabs(kn.front() - f) > Precision::PConfusion() ||
        std::fabs(kn.back() - l) > Precision::PConfusion()) {
        if (clsOut) *clsOut = "unhandled-ellipse-cyl-span";
        return -1.0;
    }

    // (C2) the section identity, so the refusal it produces is named `-oncyl`.
    {
        double th0 = 0.0, sg = 1.0, z0 = 0.0, lam = 0.0, psi = 0.0;
        if (!ellipseIsCylinderSection(cyl, el, th0, sg, z0, lam, psi)) {
            if (clsOut) *clsOut = "unhandled-ellipse-cyl-oncyl";
            return -1.0;
        }
    }

    EllipseDevClass cls = EllipseDevClass::Unhandled;
    const double m = ellipseOnCylMax(cyl, el, f, l, degree, nPoles, pu.data(), pv.data(),
                                     nKnots, kn.data(), mu.data(), &cls);
    if (!(m >= 0.0) || cls == EllipseDevClass::Unhandled) {
        // (C4) everything the math TU still refuses after (C2)/(C3) have passed
        // is the azimuthal clause: the poles do not track the section's affine
        // theta within the derived margin.
        if (clsOut) *clsOut = "unhandled-ellipse-cyl-uaffine";
        return -1.0;
    }
    if (clsOut) *clsOut = ellipseDevClassName(cls);
    return m;
}

}  // namespace refit
}  // namespace stl2step
