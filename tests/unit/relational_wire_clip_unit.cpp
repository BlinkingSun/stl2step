// D-fc-6.1(d) / RQ-8. A planar face whose polyline chord enters a tangent
// certified circle by more than tau is clipped at the shared vertex and the
// face is BRepCheck-valid. A chord that stays outside within tau is left alone.
// Synthetic geometry only. SPDX-License-Identifier: MIT

#include "refit_wire_clip.hpp"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ElSLib.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>

#include <cmath>
#include <cstdio>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846;
static constexpr double kTau = 1e-5;

static TopoDS_Vertex vtx(const gp_Pnt& p) {
    TopoDS_Vertex v;
    BRep_Builder B;
    B.MakeVertex(v, p, 1e-7);
    return v;
}

static bool selfIntersects(const TopoDS_Face& face, const TopoDS_Wire& wire) {
    if (face.IsNull() || wire.IsNull()) return false;
    BRepCheck_Wire chk(wire);
    chk.GeometricControls(Standard_True);
    TopoDS_Edge e1, e2;
    return chk.SelfIntersect(face, e1, e2, Standard_False) == BRepCheck_SelfIntersectingWire;
}

// Quarter circle V->A, then A->Q->P, then the chord P->V. All in z = 0.
struct Bowl {
    gp_Circ cert;
    TopoDS_Vertex vV, vA, vP, vQ;
    TopoDS_Edge arc, aq, qp, chord;
};

static Bowl bowl(const gp_Pnt& P) {
    Bowl b;
    b.cert = gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 1.0);
    const gp_Pnt V(1, 0, 0), A(0, 1, 0), Q(2, 2, 0);
    b.vV = vtx(V);
    b.vA = vtx(A);
    b.vP = vtx(P);
    b.vQ = vtx(Q);
    Handle(Geom_TrimmedCurve) tc =
        new Geom_TrimmedCurve(new Geom_Circle(b.cert), 0.0, 0.5 * kPi, Standard_True);
    b.arc = BRepBuilderAPI_MakeEdge(tc, b.vV, b.vA).Edge();
    b.aq = BRepBuilderAPI_MakeEdge(b.vA, b.vQ).Edge();
    b.qp = BRepBuilderAPI_MakeEdge(b.vQ, b.vP).Edge();
    b.chord = BRepBuilderAPI_MakeEdge(b.vP, b.vV).Edge();
    return b;
}

static TopoDS_Face faceOf(const TopoDS_Edge& arc, const TopoDS_Edge& aq, const TopoDS_Edge& qp,
                          const TopoDS_Edge& closing, TopoDS_Wire& wire) {
    BRepBuilderAPI_MakeWire mw;
    mw.Add(arc);
    mw.Add(aq);
    mw.Add(qp);
    mw.Add(closing);
    if (!mw.IsDone()) return TopoDS_Face();
    wire = mw.Wire();
    BRepBuilderAPI_MakeFace mf(wire, Standard_True);
    return mf.IsDone() ? mf.Face() : TopoDS_Face();
}

static bool sameEnds(const TopoDS_Edge& e, const TopoDS_Vertex& a, const TopoDS_Vertex& b) {
    TopoDS_Vertex u, v;
    TopExp::Vertices(e, u, v);
    return (u.IsSame(a) && v.IsSame(b)) || (u.IsSame(b) && v.IsSame(a));
}

int main() {
    expect(std::fabs(g5EdgeTolCap(3.0) - 0.3) < 1e-15, "g5: tenth of the certified radius");
    expect(g5EdgeTolCap(0.0) == 0.0 && g5EdgeTolCap(-1.0) == 0.0, "g5: non-positive radius is 0");

    const gp_Pnt V(1, 0, 0);
    const gp_Pnt P(0.90, 0.50, 0);
    const gp_Circ cert(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 1.0);
    const double pen = segmentDiskPenetration(cert.Location(), cert.Radius(), V, P);
    expect(pen > kTau, "crossing: the open chord enters by more than tau");

    const gp_Pnt Pout(1.02, 0.50, 0);
    expect(segmentDiskPenetration(cert.Location(), cert.Radius(), V, Pout) == 0.0,
           "outside: a chord that misses the disk has no penetration");
    const TangentClipArc miss =
        tangentClipArc(cert, V, Pout, kTau, kTau);
    expect(!miss.ok, "outside: a chord outside the circle within tau is not clipped");
    expect(tangentClipEdge(cert, vtx(V), vtx(Pout), kTau, 1e-7, g5EdgeTolCap(1.0)).IsNull(),
           "outside: the edge builder leaves that chord untouched");

    const TangentClipArc shallow = tangentClipArc(cert, V, P, pen, kTau);
    expect(!shallow.ok, "within tau: a penetration that is not greater than tau is not clipped");

    const TangentClipArc hit = tangentClipArc(cert, V, P, kTau, kTau);
    expect(hit.ok && hit.penetration > kTau && std::fabs(hit.sweep) <= kPi,
           "crossing: the tangent arc exists and is the short one");
    expect(hit.sagitta > 0.0 && hit.sagitta < g5EdgeTolCap(cert.Radius()),
           "crossing: the sagitta stays under the G5 cap");

    Bowl deep = bowl(P);
    expect(!deep.arc.IsNull() && !deep.chord.IsNull(), "face: the certified arc and the chord build");
    TopoDS_Wire rawWire;
    const TopoDS_Face raw = faceOf(deep.arc, deep.aq, deep.qp, deep.chord, rawWire);
    expect(!raw.IsNull() && selfIntersects(raw, rawWire),
           "face: the unclipped chord makes a self-intersecting wire");
    expect(!BRepCheck_Analyzer(raw, Standard_True).IsValid(),
           "face: the unclipped face is not BRepCheck-valid");

    const double cap = g5EdgeTolCap(cert.Radius());
    TopoDS_Edge clipped = tangentClipEdge(cert, deep.vV, deep.vP, kTau, 1e-7, cap);
    expect(!clipped.IsNull(), "face: the crossing chord is replaced");
    expect(sameEnds(clipped, deep.vV, deep.vP),
           "face: the clip reuses the shared vertex and the far vertex");
    Standard_Real f = 0, l = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(clipped, f, l);
    expect(!c3.IsNull() && c3->DynamicType() != STANDARD_TYPE(Geom_Line),
           "face: the replacement is not the polyline line");
    expect(BRep_Tool::Tolerance(clipped) <= cap + 1e-12, "face: the edge tolerance stays within G5");
    // The certified arc's curve is the circle we built; the clip does not replace it.
    Standard_Real af = 0, al = 0;
    Handle(Geom_Curve) ac = BRep_Tool::Curve(deep.arc, af, al);
    expect(!ac.IsNull() && ac->DynamicType() == STANDARD_TYPE(Geom_Circle) &&
               std::fabs(Handle(Geom_Circle)::DownCast(ac)->Radius() - 1.0) < 1e-12,
           "face: the certified circle is not moved");

    TopoDS_Wire fixedWire;
    const TopoDS_Face fixed = faceOf(deep.arc, deep.aq, deep.qp, clipped, fixedWire);
    expect(!fixed.IsNull() && !selfIntersects(fixed, fixedWire) &&
               BRepCheck_Analyzer(fixed, Standard_True).IsValid(),
           "face: clipped at the shared tangent vertex, the face is BRepCheck-valid");

    expect(tangentClipEdge(cert, deep.vV, deep.vP, kTau, 1e-7, hit.sagitta * 0.5).IsNull(),
           "cap: a sagitta above the tolerance cap refuses the clip");

    // An angle-parameterized arc whose chord is much longer than its parameter
    // span. A unit-speed Geom2d_Line misses by about the chord; the degree-1
    // spline on the curve's own knots misses by the sagitta.
    {
        constexpr double kRho = 20.0;
        constexpr double kSweep = 0.05;
        const gp_Circ big(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), kRho);
        const gp_Pnt BV = ElCLib::Value(0.0, big);
        const gp_Pnt BP = ElCLib::Value(kSweep, big);
        TopoDS_Vertex vB = vtx(BV), vP2 = vtx(BP);
        Handle(Geom_TrimmedCurve) tc =
            new Geom_TrimmedCurve(new Geom_Circle(big), 0.0, kSweep, Standard_True);
        TopoDS_Edge arcE = BRepBuilderAPI_MakeEdge(tc, vB, vP2).Edge();
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(arcE, f, l);
        const gp_Pln host(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)));
        const gp_Vec chord(BV, BP);
        const gp_Dir nPart(chord.Crossed(gp_Vec(0, 0, 1)));
        const gp_Pln partner(BV, nPart);
        Handle(Geom_Plane) hostS = new Geom_Plane(host);
        Handle(Geom_Plane) partS = new Geom_Plane(partner);
        Handle(Geom2d_Curve) hostC = planeSameParameterChord(host, c3->Value(f), c3->Value(l), f, l);
        Handle(Geom2d_Curve) partC =
            planeSameParameterChord(partner, c3->Value(f), c3->Value(l), f, l);
        const double sag = ElCLib::Value(0.5 * kSweep, big).Distance(gp_Pnt(0.5 * (BV.X() + BP.X()),
                                                                             0.5 * (BV.Y() + BP.Y()),
                                                                             0.5 * (BV.Z() + BP.Z())));
        const double dHost = curveSurfaceDeviation(c3, f, l, hostS, hostC);
        const double dPart = curveSurfaceDeviation(c3, f, l, partS, partC);
        std::fprintf(stderr, "chord-spline f=%.6f l=%.6f sag=%.6e dHost=%.6e dPart=%.6e\n", f, l,
                     sag, dHost, dPart);
        expect(dHost <= sag + 1e-9 && dPart <= sag + 1e-9,
               "pcurve: the same-parameter chord stays within the sagitta on both planes");

        Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0;
        ElSLib::Parameters(host, c3->Value(f), u0, v0);
        ElSLib::Parameters(host, c3->Value(l), u1, v1);
        gp_Vec2d d2(u1 - u0, v1 - v0);
        d2 /= (l - f);
        Handle(Geom2d_Curve) unitLine =
            new Geom2d_Line(gp_Pnt2d(u0 - f * d2.X(), v0 - f * d2.Y()), gp_Dir2d(d2));
        const double dUnit = curveSurfaceDeviation(c3, f, l, hostS, unitLine);
        std::fprintf(stderr, "unit-speed-line dHost=%.6e cap=%.3f\n", dUnit, g5EdgeTolCap(3.0));
        expect(dUnit > g5EdgeTolCap(3.0) && dUnit > dHost * 10.0,
               "pcurve: a unit-speed line on an angular arc exceeds the G5 cap");

        const double tolE = std::max(sag, 1e-7);
        BRep_Builder B;
        const TopLoc_Location loc;
        B.UpdateEdge(arcE, hostC, hostS, loc, tolE);
        B.Range(arcE, hostS, loc, f, l);
        B.UpdateEdge(arcE, partC, partS, loc, tolE);
        B.Range(arcE, partS, loc, f, l);
        B.UpdateEdge(arcE, tolE);
        B.SameParameter(arcE, Standard_True);
        BRepLib::SameParameter(arcE, 0.0013, Standard_True);
        Standard_Real f2 = 0, l2 = 0;
        Handle(Geom_Curve) c3b = BRep_Tool::Curve(arcE, f2, l2);
        Standard_Boolean storedH = Standard_False, storedP = Standard_False;
        Handle(Geom2d_Curve) pcH = BRep_Tool::CurveOnSurface(arcE, hostS, loc, f2, l2, &storedH);
        Handle(Geom2d_Curve) pcP = BRep_Tool::CurveOnSurface(arcE, partS, loc, f2, l2, &storedP);
        const double afterH = curveSurfaceDeviation(c3b, f, l, hostS, pcH);
        const double afterP = curveSurfaceDeviation(c3b, f, l, partS, pcP);
        const double tolAfter = BRep_Tool::Tolerance(arcE);
        std::fprintf(stderr,
                     "after-SameParameter storedH=%d storedP=%d dHost=%.6e dPart=%.6e tol=%.6e\n",
                     storedH ? 1 : 0, storedP ? 1 : 0, afterH, afterP, tolAfter);
        expect(storedH && storedP && afterH <= g5EdgeTolCap(3.0) && afterP <= g5EdgeTolCap(3.0) &&
                   tolAfter <= g5EdgeTolCap(3.0) && afterH <= tolAfter + 1e-7 &&
                   afterP <= tolAfter + 1e-7,
               "pcurve: SameParameter at sew scale keeps both planes inside the G5 cap");
    }

    // D-fc-7.3 RQ-HUB-3 (b), D-fc-7.6: the tangency point of two circles of one plane.
    {
        const gp_Pnt o0(0, 0, 0);
        const double tol = 2 * kTau;
        gp_Pnt T;
        const gp_Pnt o1(12.0 - 0.5 * kTau, 0, 0);
        expect(circleTangencyPoint(o0, 9.0, o1, 3.0, tol, T) &&
                   std::fabs(T.X() - (9.0 - 0.25 * kTau)) < 1e-12 && std::fabs(T.Y()) < 1e-15,
               "tangency: external pair meets on the line of centres");
        const gp_Pnt o2(0, 6.0 + 0.3 * kTau, 0);
        expect(circleTangencyPoint(o0, 10.0, o2, 4.0, tol, T) &&
                   std::fabs(T.Y() - (10.0 + 0.15 * kTau)) < 1e-12,
               "tangency: internal pair meets beyond the small centre");
        expect(circleTangencyPoint(o2, 4.0, o0, 10.0, tol, T) &&
                   std::fabs(T.Y() - (10.0 + 0.15 * kTau)) < 1e-12,
               "tangency: the order of the pair does not matter");
        expect(!circleTangencyPoint(o0, 9.0, gp_Pnt(12.0 + 3 * tol, 0, 0), 3.0, tol, T),
               "tangency: circles apart by more than tanDistTol are not tangent");
        expect(!circleTangencyPoint(o0, 9.0, gp_Pnt(0.5 * tol, 0, 0), 9.0, tol, T),
               "tangency: concentric circles have no line of centres");
    }
    // D-fc-7.6: an arc of one exact rim split at a vertex on its circle.
    {
        const gp_Circ c(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 2.0);
        const TopoDS_Vertex a = vtx(ElCLib::Value(0.0, c)), b = vtx(ElCLib::Value(0.5 * kPi, c));
        Handle(Geom_TrimmedCurve) tc =
            new Geom_TrimmedCurve(new Geom_Circle(c), 0.0, 0.5 * kPi, Standard_True);
        const TopoDS_Edge arc = BRepBuilderAPI_MakeEdge(tc, a, b).Edge();
        const TopoDS_Vertex m = vtx(ElCLib::Value(0.25 * kPi, c));
        TopoDS_Edge e1, e2;
        expect(splitArcAt(arc, m, e1, e2) && sameEnds(e1, a, m) && sameEnds(e2, m, b),
               "split arc: two arcs through the shared vertices");
        Standard_Real f1 = 0, l1 = 0, f2 = 0, l2 = 0;
        Handle(Geom_Curve) c1 = BRep_Tool::Curve(e1, f1, l1);
        Handle(Geom_Curve) c2 = BRep_Tool::Curve(e2, f2, l2);
        expect(!c1.IsNull() && !c2.IsNull() && std::fabs(f1) < 1e-12 &&
                   std::fabs(l1 - 0.25 * kPi) < 1e-12 && std::fabs(f2 - 0.25 * kPi) < 1e-12 &&
                   std::fabs(l2 - 0.5 * kPi) < 1e-12,
               "split arc: the same circle, the parameter range divided at the vertex");
        expect(!splitArcAt(arc, vtx(ElCLib::Value(kPi, c)), e1, e2),
               "split arc: a vertex off the arc's range does not split it");
        expect(!splitArcAt(BRepBuilderAPI_MakeEdge(a, b).Edge(), m, e1, e2),
               "split arc: a line is not an arc");
    }

    if (gFail) std::fprintf(stderr, "%d FAIL\n", gFail);
    else std::fprintf(stderr, "ALL PASS\n");
    return gFail ? 1 : 0;
}
