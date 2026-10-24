// D-fc-6.1(d) / RQ-8. A mesh-polyline segment of an analytic plane face that
// enters a certified circle already on the same wire, by more than tau, and
// shares that circle's end vertex, is replaced by the circle through the same
// two vertices that is tangent to the certified circle there. The certified
// circle and both vertices stay where they are. No new threshold: tau = tauOf
// (2q), and the replacement's sagitta may not exceed the larger of the edge's
// existing tolerance and the D-fc-6.7 G5 bound 0.1*R of the certified circle.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_WIRE_CLIP_HPP
#define STL2STEP_REFIT_WIRE_CLIP_HPP

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include <cmath>

namespace stl2step {
namespace refit {

// Deepest entry of the open segment AB into the disk (O, R). Positive is
// inside. Zero when the closest point of the segment is an endpoint, so a
// chord that only touches the circle at its shared vertex does not count.
inline double segmentDiskPenetration(const gp_Pnt& O, double R, const gp_Pnt& A,
                                     const gp_Pnt& B) {
    const gp_Vec AB(A, B);
    const double L2 = AB.SquareMagnitude();
    if (!(L2 > 0.0) || !(R > 0.0)) return 0.0;
    const double t = gp_Vec(A, O).Dot(AB) / L2;
    if (!(t > 1e-9 && t < 1.0 - 1e-9)) return 0.0;
    const gp_Pnt C = A.Translated(AB * t);
    const double pen = R - C.Distance(O);
    return pen > 0.0 ? pen : 0.0;
}

// D-fc-7.3 RQ-HUB-3 (b), D-fc-7.6. Two circles of one plane wire (centres o0, o1,
// radii r0, r1) tangent within `tol` -- externally (|o0 o1| = r0 + r1) or
// internally (|o0 o1| = |r0 - r1|) -- meet at their tangency point on the line of
// centres: the midpoint of the two circles' own points there. Their intersection
// roots lie +-sqrt(2 r_eff delta) along the common tangent for an overlap delta
// (ill-conditioned); the tangency point is not. False when the centres coincide
// within `tol` or the circles are not tangent within it.
inline bool circleTangencyPoint(const gp_Pnt& o0, double r0, const gp_Pnt& o1, double r1,
                                double tol, gp_Pnt& T) {
    if (!(tol > 0.0) || !(r0 > 0.0) || !(r1 > 0.0)) return false;
    const double d = o0.Distance(o1);
    if (!(d > tol)) return false;
    const double ext = std::fabs(d - (r0 + r1));
    const double inn = std::fabs(d - std::fabs(r0 - r1));
    if (!(std::min(ext, inn) <= tol)) return false;
    const gp_Vec u = gp_Vec(o0, o1) / d;
    // External: o0's point faces o1 and o1's faces o0. Internal: both points lie
    // on the side of the larger circle's centre away from the smaller one.
    double s0 = 1.0, s1 = -1.0;
    if (inn < ext) {
        s0 = r0 >= r1 ? 1.0 : -1.0;
        s1 = s0;
    }
    const gp_Pnt t0 = o0.Translated(u * (s0 * r0));
    const gp_Pnt t1 = o1.Translated(u * (s1 * r1));
    T = gp_Pnt(0.5 * (t0.XYZ() + t1.XYZ()));
    return true;
}

// D-fc-7.6 (the D-fc-6.11 (c) mint on a rim published as several arcs). Split the
// circular arc edge `e` at the vertex V standing on its circle: two edges on the same
// circle, e's first vertex to V and V to e's last vertex (parameter order of e's own
// curve), each FORWARD. False when e is not an arc of a circle or V's parameter is
// not strictly inside the arc. No geometry moves; e's tolerance is kept.
inline bool splitArcAt(const TopoDS_Edge& e, const TopoDS_Vertex& V, TopoDS_Edge& e1,
                       TopoDS_Edge& e2) {
    if (e.IsNull() || V.IsNull()) return false;
    try {
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull() || !(l > f)) return false;
        Handle(Geom_Curve) b = c;
        if (b->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
            b = Handle(Geom_TrimmedCurve)::DownCast(b)->BasisCurve();
        Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(b);
        if (gc.IsNull()) return false;
        constexpr double kTwoPi = 6.28318530717958647692;
        double t = ElCLib::Parameter(gc->Circ(), BRep_Tool::Pnt(V));
        while (t < f) t += kTwoPi;
        while (t >= f + kTwoPi) t -= kTwoPi;
        if (!(t > f + Precision::PConfusion() && t < l - Precision::PConfusion())) return false;
        const TopoDS_Edge ef = TopoDS::Edge(e.Oriented(TopAbs_FORWARD));
        TopoDS_Vertex vf, vl;
        TopExp::Vertices(ef, vf, vl);
        if (vf.IsNull() || vl.IsNull()) return false;
        BRepBuilderAPI_MakeEdge m1(gc, vf, V, f, t), m2(gc, V, vl, t, l);
        if (!m1.IsDone() || !m2.IsDone()) return false;
        e1 = m1.Edge();
        e2 = m2.Edge();
        BRep_Builder B;
        B.UpdateEdge(e1, BRep_Tool::Tolerance(e));
        B.UpdateEdge(e2, BRep_Tool::Tolerance(e));
        return true;
    } catch (const Standard_Failure&) {
        return false;
    }
}

// D-fc-6.7. An edge tolerance may not exceed a tenth of the certified radius.
inline double g5EdgeTolCap(double radius) {
    return (radius > 0.0 && std::isfinite(radius)) ? 0.1 * radius : 0.0;
}

struct TangentClipArc {
    bool ok = false;
    gp_Circ circ;
    double penetration = 0;
    double sagitta = 0;
    double sweep = 0;  // signed parameter of P, V at 0; |sweep| <= pi
};

// The unique circle through V and P whose centre lies on the certified
// circle's radial through V. Refuses unless V is on the certified circle
// within tau, P lies outside it by more than tau, the segment enters by more
// than tau, and the short arc stays outside the certified disk within tau.
// Building that circle does not move V or P: a point farther than `snap` from
// the circle refuses the clip.
inline TangentClipArc tangentClipArc(const gp_Circ& cert, const gp_Pnt& V, const gp_Pnt& P,
                                     double tau, double snap) {
    TangentClipArc out;
    if (!(tau > 0.0) || !(cert.Radius() > 0.0) || !(snap >= 0.0)) return out;
    const gp_Pnt O = cert.Location();
    const double R = cert.Radius();
    const gp_Dir axis = cert.Axis().Direction();
    gp_Vec radial(O, V);
    radial -= gp_Vec(axis) * radial.Dot(axis);
    const double rV = radial.Magnitude();
    if (!(rV > Precision::Confusion())) return out;
    if (std::fabs(rV - R) > tau) return out;
    const gp_Vec n = radial / rV;
    if (P.Distance(O) < R + tau) return out;
    const double pen = segmentDiskPenetration(O, R, V, P);
    out.penetration = pen;
    if (!(pen > tau)) return out;
    const gp_Vec OP(O, P);
    const double denom = 2.0 * (n.Dot(OP) - rV);
    if (!(std::fabs(denom) > tau)) return out;
    const double t = (OP.SquareMagnitude() - rV * rV) / denom;
    if (!std::isfinite(t)) return out;
    const gp_Pnt Q = O.Translated(n * t);
    const double rho = Q.Distance(V);
    if (!(rho > Precision::Confusion()) || !std::isfinite(rho)) return out;
    gp_Vec xv(Q, V);
    if (xv.Magnitude() <= Precision::Confusion()) return out;
    xv.Normalize();
    const gp_Circ arc(gp_Ax2(Q, axis, gp_Dir(xv)), rho);
    auto offCircle = [&](const gp_Pnt& X) {
        const gp_Vec d(Q, X);
        const double axial = std::fabs(d.Dot(axis));
        const double radial = (d - gp_Vec(axis) * d.Dot(axis)).Magnitude();
        return std::max(axial, std::fabs(radial - rho));
    };
    if (offCircle(V) > snap || offCircle(P) > snap) return out;
    double du = ElCLib::Parameter(arc, P);
    constexpr double kPi = 3.14159265358979323846;
    while (du > kPi) du -= 2.0 * kPi;
    while (du <= -kPi) du += 2.0 * kPi;
    if (!(std::fabs(du) > 1e-12) || !std::isfinite(du)) return out;
    const gp_Pnt chordMid(0.5 * (V.X() + P.X()), 0.5 * (V.Y() + P.Y()), 0.5 * (V.Z() + P.Z()));
    for (int i = 1; i <= 4; ++i) {
        const gp_Pnt s = ElCLib::Value(du * (double)i / 4.0, arc);
        if (R - s.Distance(O) > tau) return out;
        if (i == 2) out.sagitta = s.Distance(chordMid);
    }
    if (!std::isfinite(out.sagitta)) return out;
    out.circ = arc;
    out.sweep = du;
    out.ok = true;
    return out;
}

// Degree-1 pcurve of the chord on `pln`, parameterized on the 3d curve's own
// [f, l]. Geom2d_Line cannot carry this: its direction is unit, so a parameter
// that is an angle walks one unit of length per radian and the gap becomes
// the chord. The spline's gap on a plane that contains the chord is the
// sagitta of the arc, which the caller compares with g5EdgeTolCap.
inline Handle(Geom2d_Curve) planeSameParameterChord(const gp_Pln& pln, const gp_Pnt& pF,
                                                   const gp_Pnt& pL, double f, double l) {
    if (!(l > f) || !std::isfinite(f) || !std::isfinite(l)) return {};
    try {
        Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0;
        ElSLib::Parameters(pln, pF, u0, v0);
        ElSLib::Parameters(pln, pL, u1, v1);
        const gp_Vec2d d(u1 - u0, v1 - v0);
        if (!(d.SquareMagnitude() > Precision::PConfusion() * Precision::PConfusion())) return {};
        TColgp_Array1OfPnt2d poles(1, 2);
        poles.SetValue(1, gp_Pnt2d(u0, v0));
        poles.SetValue(2, gp_Pnt2d(u1, v1));
        TColStd_Array1OfReal knots(1, 2);
        knots.SetValue(1, f);
        knots.SetValue(2, l);
        TColStd_Array1OfInteger mults(1, 2);
        mults.SetValue(1, 2);
        mults.SetValue(2, 2);
        return new Geom2d_BSplineCurve(poles, knots, mults, 1);
    } catch (const Standard_Failure&) {
        return {};
    }
}

// Max |C3d(t) - S(C2d(t))| at 9 uniform parameters. A null curve is unusable.
inline double curveSurfaceDeviation(const Handle(Geom_Curve)& c3, double f, double l,
                                    const Handle(Geom_Surface)& srf,
                                    const Handle(Geom2d_Curve)& c2d) {
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull() || !(l > f)) return 1.0e300;
    double m = 0.0;
    try {
        for (int i = 0; i <= 8; ++i) {
            const double t = f + (l - f) * (double)i / 8.0;
            const gp_Pnt P = c3->Value(t);
            const gp_Pnt2d uv = c2d->Value(t);
            m = std::max(m, P.Distance(srf->Value(uv.X(), uv.Y())));
        }
    } catch (const Standard_Failure&) {
        return 1.0e300;
    }
    return m;
}

// An edge on the clip arc between the two existing vertices. Null when the
// clip does not apply, the sagitta exceeds tolCap, or the edge would mint a
// new vertex. Tolerance becomes the larger of `tol` and the sagitta.
inline TopoDS_Edge tangentClipEdge(const gp_Circ& cert, const TopoDS_Vertex& vV,
                                   const TopoDS_Vertex& vP, double tau, double tol,
                                   double tolCap) {
    if (vV.IsNull() || vP.IsNull()) return TopoDS_Edge();
    const double snap = std::max(tau, std::max(BRep_Tool::Tolerance(vV), BRep_Tool::Tolerance(vP)));
    const TangentClipArc clip = tangentClipArc(cert, BRep_Tool::Pnt(vV), BRep_Tool::Pnt(vP), tau, snap);
    if (!clip.ok) return TopoDS_Edge();
    const double need = std::max(tol, clip.sagitta);
    if (!(tolCap > 0.0) || need > tolCap) return TopoDS_Edge();
    Handle(Geom_Circle) gc = new Geom_Circle(clip.circ);
    const double u0 = std::min(0.0, clip.sweep);
    const double u1 = std::max(0.0, clip.sweep);
    Handle(Geom_TrimmedCurve) tc = new Geom_TrimmedCurve(gc, u0, u1, Standard_True);
    const bool forward = clip.sweep >= 0.0;
    BRepBuilderAPI_MakeEdge me(tc, forward ? vV : vP, forward ? vP : vV);
    if (!me.IsDone()) return TopoDS_Edge();
    TopoDS_Edge e = me.Edge();
    TopoDS_Vertex a, b;
    TopExp::Vertices(e, a, b);
    const bool same = (a.IsSame(vV) && b.IsSame(vP)) || (a.IsSame(vP) && b.IsSame(vV));
    if (!same) return TopoDS_Edge();
    BRep_Builder B;
    B.UpdateEdge(e, std::max(need, Precision::Confusion()));
    return e;
}

}  // namespace refit
}  // namespace stl2step

#endif
