// stl2step cone BIND site — see refit_cone_bind.hpp for the certificates.
//
// SPDX-License-Identifier: MIT

#include "refit_cone_bind.hpp"

#include <algorithm>
#include <cmath>

#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Lin2d.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include "refit_cone_math.hpp"

namespace stl2step {
namespace refit {
namespace {

Handle(Geom_Curve) basisOf3d(const Handle(Geom_Curve)& c) {
    Handle(Geom_TrimmedCurve) t = Handle(Geom_TrimmedCurve)::DownCast(c);
    if (!t.IsNull() && !t->BasisCurve().IsNull()) return t->BasisCurve();
    return c;
}

Handle(Geom2d_Curve) basisOf2d(const Handle(Geom2d_Curve)& c) {
    Handle(Geom2d_TrimmedCurve) t = Handle(Geom2d_TrimmedCurve)::DownCast(c);
    if (!t.IsNull() && !t->BasisCurve().IsNull()) return t->BasisCurve();
    return c;
}

// |C3d(t) - S(pcurve(t))| at one parameter. Both operands are in world space:
// the surface point carries `loc`, exactly as the bind site's own
// curveSurfDevAtBind does.
double gapAt(const Handle(Geom_Curve)& c3, double t, const Handle(Geom_Surface)& srf,
             const Handle(Geom2d_Curve)& c2d, const TopLoc_Location& loc) {
    const gp_Pnt P = c3->Value(t);
    const gp_Pnt2d uv = c2d->Value(t);
    gp_Pnt Q = srf->Value(uv.X(), uv.Y());
    if (!loc.IsIdentity()) Q.Transform(loc.Transformation());
    return P.Distance(Q);
}

// The world-space cone under `srf`, or false.
bool worldCone(const Handle(Geom_Surface)& srf, const TopLoc_Location& loc, gp_Cone& out) {
    Handle(Geom_ConicalSurface) cs = coneBasisOf(srf);
    if (cs.IsNull()) return false;
    out = cs->Cone();
    if (!loc.IsIdentity()) out.Transform(loc.Transformation());
    return true;
}

}  // namespace

Handle(Geom_ConicalSurface) coneBasisOf(const Handle(Geom_Surface)& srf) {
    if (srf.IsNull()) return {};
    Handle(Geom_ConicalSurface) c = Handle(Geom_ConicalSurface)::DownCast(srf);
    if (!c.IsNull()) return c;
    Handle(Geom_RectangularTrimmedSurface) t =
        Handle(Geom_RectangularTrimmedSurface)::DownCast(srf);
    if (!t.IsNull()) return Handle(Geom_ConicalSurface)::DownCast(t->BasisSurface());
    return {};
}

Handle(Geom2d_Curve) conePCurveForCircle(const gp_Cone& cone, const gp_Circ& circ) {
    const gp_Ax3 pos = cone.Position();
    const double cosA = std::cos(cone.SemiAngle());
    if (!(std::fabs(cosA) > Precision::Confusion())) return {};
    const double dotAx = circ.Axis().Direction().Dot(pos.Direction());
    if (std::fabs(dotAx) < 1.0 - Precision::Angular()) return {};
    const gp_Dir cx = circ.XAxis().Direction();
    const double phi = std::atan2(cx.Dot(pos.YDirection()), cx.Dot(pos.XDirection()));
    const double z = gp_Vec(pos.Location(), circ.Location()).Dot(pos.Direction());
    const double v0 = coneVAtZ(cone, z);
    if (!std::isfinite(v0)) return {};
    // u must sweep with the circle's own rotation sense: du * Zcone == circle axis.
    const double sgn = (dotAx >= 0.0) ? 1.0 : -1.0;
    return new Geom2d_Line(gp_Pnt2d(phi, v0), gp_Dir2d(sgn, 0.0));
}

Handle(Geom2d_Curve) conePCurveForLine(const gp_Cone& cone, const gp_Lin& lin, double tMid,
                                       double angTol) {
    const gp_Ax3 pos = cone.Position();
    const double sinA = std::sin(cone.SemiAngle());
    const double cosA = std::cos(cone.SemiAngle());
    if (!(std::fabs(cosA) > Precision::Confusion())) return {};
    // u0 from where the line actually sits, taken at its own mid parameter so a
    // seam that starts on the axis (it cannot, on a frustum) is not consulted at
    // a degenerate point.
    const gp_Pnt pmid = lin.Location().Translated(gp_Vec(lin.Direction()) * tMid);
    gp_Vec rho(pos.Location(), pmid);
    const gp_Vec ax(pos.Direction());
    gp_Vec rad = rho - ax * rho.Dot(ax);
    if (!(rad.Magnitude() > Precision::Confusion())) return {};
    const double u0 = std::atan2(rad.Dot(gp_Vec(pos.YDirection())),
                                 rad.Dot(gp_Vec(pos.XDirection())));
    // The generator through u0, unit speed in v.
    const gp_Vec rhoHat(std::cos(u0) * pos.XDirection().XYZ() +
                        std::sin(u0) * pos.YDirection().XYZ());
    const gp_Vec gen = rhoHat * sinA + ax * cosA;
    const double d = gp_Vec(lin.Direction()).Dot(gen);
    if (std::fabs(std::fabs(d) - 1.0) > angTol) return {};
    const double v0 = gp_Vec(pos.Location(), lin.Location()).Dot(ax) / cosA;
    if (!std::isfinite(v0)) return {};
    const double sgn = (d >= 0.0) ? 1.0 : -1.0;
    return new Geom2d_Line(gp_Pnt2d(u0, v0), gp_Dir2d(0.0, sgn));
}

double coneCircleBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                         const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                         const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled-circle-on-cone";
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return -1.0;
    gp_Cone cone;
    if (!worldCone(srf, loc, cone)) return -1.0;
    Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(basisOf3d(c3));
    if (gc.IsNull()) return -1.0;
    Handle(Geom2d_Line) pl2 = Handle(Geom2d_Line)::DownCast(basisOf2d(c2d));
    if (pl2.IsNull()) return -1.0;

    // (1) v constant along the pcurve.
    const gp_Dir2d d2 = pl2->Lin2d().Direction();
    if (std::fabs(d2.Y()) > Precision::PConfusion()) return -1.0;
    const double du = d2.X();
    if (!(std::fabs(du) > 0.0)) return -1.0;

    const gp_Circ circ = gc->Circ();
    const gp_Ax3 pos = cone.Position();
    // (2) axes parallel and (4) matching angular velocity vector.
    const double dotAx = circ.Axis().Direction().Dot(pos.Direction());
    if (std::fabs(dotAx) < 1.0 - Precision::Angular()) return -1.0;
    if (du * dotAx <= 0.0) return -1.0;
    // (3) centre on the cone axis.
    {
        const gp_Vec w(pos.Location(), circ.Location());
        const gp_Vec ax(pos.Direction());
        if ((w - ax * w.Dot(ax)).Magnitude() > Precision::Confusion()) return -1.0;
    }

    // Certified: the gap is constant in t. Evaluate at both ends (they agree by
    // the certificate; taking the max costs nothing and is the safe direction)
    // and carry the surface-distance supremum published by 130-CONE-MATH, which
    // this quantity dominates -- a perpendicular distance is never larger than
    // the parametric gap that contains it.
    double m = std::max(gapAt(c3, f, srf, c2d, loc), gapAt(c3, l, srf, c2d, loc));
    ConeDevClass cls = ConeDevClass::Unhandled;
    const double surfSup = circleOnConeMax(cone, circ, &cls);
    if (!(surfSup >= 0.0) || !coneDevClassIsExact(cls)) return -1.0;
    m = std::max(m, surfSup);
    if (clsOut) *clsOut = "circle-on-cone";
    return m;
}

double coneLineBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                       const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                       const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled-line-on-cone";
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return -1.0;
    gp_Cone cone;
    if (!worldCone(srf, loc, cone)) return -1.0;
    if (Handle(Geom_Line)::DownCast(basisOf3d(c3)).IsNull()) return -1.0;
    Handle(Geom2d_Line) pl2 = Handle(Geom2d_Line)::DownCast(basisOf2d(c2d));
    if (pl2.IsNull()) return -1.0;
    // u constant => the surface image is the generator, affine in the parameter.
    if (std::fabs(pl2->Lin2d().Direction().X()) > Precision::PConfusion()) return -1.0;
    const double m = std::max(gapAt(c3, f, srf, c2d, loc), gapAt(c3, l, srf, c2d, loc));
    if (clsOut) *clsOut = "line-on-cone";
    return m;
}

double coneBindSup(const Handle(Geom_Curve)& c3, double f, double l,
                   const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                   const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled-curve-on-cone";
    if (c3.IsNull()) return -1.0;
    const Handle(Geom_Curve) src = basisOf3d(c3);
    if (!Handle(Geom_Circle)::DownCast(src).IsNull())
        return coneCircleBindSup(c3, f, l, srf, c2d, loc, clsOut);
    if (!Handle(Geom_Line)::DownCast(src).IsNull())
        return coneLineBindSup(c3, f, l, srf, c2d, loc, clsOut);
    return -1.0;
}

}  // namespace refit
}  // namespace stl2step
