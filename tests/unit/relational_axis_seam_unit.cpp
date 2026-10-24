// Synthetic cases for the torus-family seam predicates (D-fc-4.8), the closed-wall /
// closed-torus build predicates of fx-closed360, the crossing-cap predicates of
// D-fc-5.7 O-1, the U-3 naming predicates and the split-cap / single-wire
// orientation rules of D-fc-5.7 O-4. No mesh.
// D-fc-5.7 O-1, the U-3 naming predicates and the closed-ellipse rim construction of
// D-fc-5.13 and the D-fc-5.20 refusal of a non-positive quantum. No mesh file.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_axis_seam.hpp"
#include "refit_edge_bound.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <Geom2d_Line.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static const char* brepCheckStatusName(int st) {
    switch (st) {
        case BRepCheck_NoError: return "NoError";
        case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
        case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
        case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
        case BRepCheck_No3DCurve: return "No3DCurve";
        case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
        case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
        case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
        case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
        case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
        case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
        case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
        case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
        case BRepCheck_FreeEdge: return "FreeEdge";
        case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
        case BRepCheck_InvalidRange: return "InvalidRange";
        case BRepCheck_EmptyWire: return "EmptyWire";
        case BRepCheck_RedundantEdge: return "RedundantEdge";
        case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
        case BRepCheck_NoSurface: return "NoSurface";
        case BRepCheck_InvalidWire: return "InvalidWire";
        case BRepCheck_RedundantWire: return "RedundantWire";
        case BRepCheck_IntersectingWires: return "IntersectingWires";
        case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
        case BRepCheck_EmptyShell: return "EmptyShell";
        case BRepCheck_RedundantFace: return "RedundantFace";
        case BRepCheck_InvalidImbricationOfShells: return "InvalidImbricationOfShells";
        case BRepCheck_UnorientableShape: return "UnorientableShape";
        case BRepCheck_NotClosed: return "NotClosed";
        case BRepCheck_NotConnected: return "NotConnected";
        case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
        case BRepCheck_BadOrientation: return "BadOrientation";
        case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
        case BRepCheck_InvalidPolygonOnTriangulation: return "InvalidPolygonOnTriangulation";
        case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
        case BRepCheck_EnclosedRegion: return "EnclosedRegion";
        case BRepCheck_CheckFail: return "CheckFail";
        default: return "?";
    }
}

static void fprintfBrepResultStatuses(const char* label, const Handle(BRepCheck_Result)& res) {
    if (res.IsNull()) {
        std::fprintf(stderr, "  %s: (null BRepCheck_Result)\n", label);
        return;
    }
    std::fprintf(stderr, "  %s:", label);
    bool any = false;
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next()) {
        const int st = static_cast<int>(it.Value());
        std::fprintf(stderr, "%s%s", any ? "," : " ", brepCheckStatusName(st));
        any = true;
    }
    if (!any) std::fprintf(stderr, " NoError");
    std::fprintf(stderr, "\n");
}

static void dumpFaceBrepCheck(const TopoDS_Face& face) {
    BRepCheck_Analyzer an(face, Standard_True);
    std::fprintf(stderr, "BRepCheck dump (IsValid=%d):\n", an.IsValid() ? 1 : 0);
    fprintfBrepResultStatuses("face", an.Result(face));
    for (TopExp_Explorer wx(face, TopAbs_WIRE); wx.More(); wx.Next())
        fprintfBrepResultStatuses("wire", an.Result(wx.Current()));
    int ei = 0;
    for (TopExp_Explorer ex(face, TopAbs_EDGE); ex.More(); ex.Next(), ++ei) {
        const TopoDS_Edge e = TopoDS::Edge(ex.Current());
        std::fprintf(stderr, "  edge[%d] tol=%.9g:", ei, BRep_Tool::Tolerance(e));
        const Handle(BRepCheck_Result) er = an.Result(e);
        if (er.IsNull()) {
            std::fprintf(stderr, " (null)\n");
            continue;
        }
        bool any = false;
        for (BRepCheck_ListOfStatus::Iterator it(er->Status()); it.More(); it.Next()) {
            const int st = static_cast<int>(it.Value());
            std::fprintf(stderr, "%s%s", any ? "," : " ", brepCheckStatusName(st));
            any = true;
        }
        if (!any) std::fprintf(stderr, " NoError");
        std::fprintf(stderr, "\n");
    }
}

static void expectFaceValidBeforeHeal(const TopoDS_Face& face, Standard_Boolean wantValid,
                                     const char* name) {
    BRepCheck_Analyzer before(face, Standard_True);
    const bool got = before.IsValid() == Standard_True;
    const bool want = wantValid == Standard_True;
    if (got != want) dumpFaceBrepCheck(face);
    expect(got == want, name);
}

// Same dense sampling and margin as fitEdgeTolerancesToFace (refit_build.cpp:11551-11560).
static double maxCurveOnSurfaceDeviation(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                                         const Handle(Geom_Surface)& srf,
                                         const Handle(Geom2d_Curve)& pc,
                                         const TopLoc_Location& loc) {
    if (c3.IsNull() || pc.IsNull() || srf.IsNull() || l <= f) return 0.0;
    Handle(Geom_Surface) sb =
        Handle(Geom_Surface)::DownCast(srf->Transformed(loc.Transformation()));
    if (sb.IsNull()) sb = srf;
    const gp_Trsf tr = loc.Transformation();
    double maxD = 0.0;
    auto atT = [&](double t) {
        const gp_Pnt2d uv = pc->Value(t);
        gp_Pnt ps = sb->Value(uv.X(), uv.Y());
        if (!loc.IsIdentity()) ps.Transform(tr);
        maxD = std::max(maxD, c3->Value(t).Distance(ps));
    };
    for (int k = 0; k <= 64; ++k) atT(f + (l - f) * (static_cast<double>(k) / 64.0));
    atT(f);
    atT(l);
    atT(0.5 * (f + l));
    return maxD;
}

static double edgeTolFromMeasuredDeviation(double d) { return d * 1.001 + Precision::Confusion(); }

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kTau = 1e-4;

static gp_Pnt onCircle(double R, double u, double z) {
    return gp_Pnt(R * std::cos(u), R * std::sin(u), z);
}

// fx-hubtori: a mesh chord on a torus with a straight UV pcurve (FixAddPCurve's
// shape for a short chord), parametrised as the 3D line on [0, L].
struct ChordOnTorus {
    TopoDS_Face face;
    TopoDS_Edge edge;
    double L = 0;
};

static ChordOnTorus chordOnTorus(const Handle(Geom_ToroidalSurface)& S, double u0, double v0,
                                 double u1, double v1, double tol) {
    ChordOnTorus c;
    BRep_Builder B;
    B.MakeFace(c.face, S, Precision::Confusion());
    const gp_Pnt P0 = S->Value(u0, v0), P1 = S->Value(u1, v1);
    c.L = P0.Distance(P1);
    Handle(Geom_Line) ln = new Geom_Line(P0, gp_Dir(gp_Vec(P0, P1)));
    TColgp_Array1OfPnt2d poles(1, 2);
    poles(1) = gp_Pnt2d(u0, v0);
    poles(2) = gp_Pnt2d(u1, v1);
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0.0;
    knots(2) = c.L;
    TColStd_Array1OfInteger mults(1, 2);
    mults(1) = mults(2) = 2;
    Handle(Geom2d_BSplineCurve) pc = new Geom2d_BSplineCurve(poles, knots, mults, 1);
    B.MakeEdge(c.edge, ln, tol);
    TopoDS_Vertex a, b;
    B.MakeVertex(a, P0, tol);
    B.MakeVertex(b, P1, tol);
    B.Add(c.edge, a.Oriented(TopAbs_FORWARD));
    B.Add(c.edge, b.Oriented(TopAbs_REVERSED));
    B.Range(c.edge, 0.0, c.L);
    B.UpdateEdge(c.edge, pc, c.face, tol);
    return c;
}

static double sampledDev(const Handle(Geom_ToroidalSurface)& S, const ChordOnTorus& c, int n) {
    Standard_Real f = 0, l = 0, a = 0, b = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(c.edge, f, l);
    Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(c.edge, c.face, a, b);
    double m = 0;
    for (int i = 0; i < n; ++i) {
        const double t = f + (l - f) * (double)i / (double)(n - 1);
        const gp_Pnt2d q = pc->Value(t);
        m = std::max(m, S->Value(q.X(), q.Y()).Distance(c3->Value(t)));
    }
    return m;
}

static void caseEdgeBound() {
    Handle(Geom_ToroidalSurface) S =
        new Geom_ToroidalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 12.0, 3.0);
    // A diagonal chord between two rows of a fillet tube (u and v both move).
    const ChordOnTorus d = chordOnTorus(S, 0.30, 2.62, 0.43, 2.88, Precision::Confusion());
    const double ex = exactEdgeDeviationOnFace(d.edge, d.face);
    const double grid = sampledDev(S, d, 23);
    const double dense = sampledDev(S, d, 4001);
    const double mid = sampledDev(S, d, 3);
    expect(ex > 0.0, "edge bound: exact deviation measured");
    expect(ex >= grid && ex >= dense && ex >= mid, "edge bound: exact >= every sampled reading");
    expect(ex <= dense + 2.0 * Precision::Confusion(), "edge bound: exact is the maximum, not a bound");
    expect(edgeBoundOf(ex, mid, 1.0) == EdgeBound::Raise, "edge bound: midpoint tolerance is raised");
    expect(edgeBoundOf(ex, ex, 1.0) == EdgeBound::Keep, "edge bound: an exact tolerance is kept");
    expect(edgeBoundOf(ex, 10.0 * ex, 1.0) == EdgeBound::Keep, "edge bound: never lowered");
    expect(edgeBoundOf(ex, mid, 0.5 * ex) == EdgeBound::OverCap, "edge bound: never past cap");
    expect(edgeBoundOf(-1.0, mid, 1.0) == EdgeBound::Keep, "edge bound: unmeasured edge kept");
    // A chord along one row (u only): its deviation is the sagitta at the middle.
    const ChordOnTorus r = chordOnTorus(S, 0.30, kPi, 0.45, kPi, Precision::Confusion());
    const double exR = exactEdgeDeviationOnFace(r.edge, r.face);
    const double sag = 9.0 * (1.0 - std::cos(0.075));
    expect(std::fabs(exR - Precision::Confusion() - sag) <= 1e-9, "edge bound: row chord = its sagitta");
}

// A single wire on surface S through the (u, v) corners in the given order, each side an
// iso line built from its pcurve; the face is made bare (no MakeFace orientation pass).
static TopoDS_Face faceThroughUV(const Handle(Geom_Surface)& S, const std::vector<gp_Pnt2d>& c) {
    BRepBuilderAPI_MakeWire mw;
    for (size_t i = 0; i < c.size(); ++i) {
        const gp_Pnt2d a = c[i], b = c[(i + 1) % c.size()];
        Handle(Geom2d_Line) l = new Geom2d_Line(a, gp_Dir2d(gp_Vec2d(a, b)));
        BRepBuilderAPI_MakeEdge me(l, S, 0.0, a.Distance(b));
        TopoDS_Edge e = me.Edge();
        BRepLib::BuildCurves3d(e);
        mw.Add(e);
    }
    TopoDS_Face f;
    BRep_Builder B;
    B.MakeFace(f, S, Precision::Confusion());
    B.Add(f, mw.Wire());
    return f;
}

static int edgeCount(const TopoDS_Shape& s) {
    int n = 0;
    for (TopExp_Explorer ex(s, TopAbs_EDGE); ex.More(); ex.Next()) ++n;
    return n;
}

static void caseInsideOut() {
    // A patch of a convex torus tube (u and v ranges of a partial fillet).
    Handle(Geom_ToroidalSurface) S =
        new Geom_ToroidalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 10.0, 2.5);
    const std::vector<gp_Pnt2d> ccw = {gp_Pnt2d(0.5, 0.2), gp_Pnt2d(2.0, 0.2), gp_Pnt2d(2.0, 1.3),
                                       gp_Pnt2d(0.5, 1.3)};
    const std::vector<gp_Pnt2d> cw = {ccw[0], ccw[3], ccw[2], ccw[1]};
    TopoDS_Face good = faceThroughUV(S, ccw), bad = faceThroughUV(S, cw);
    expect(!singleWireInsideOut(good), "inside out: a counter-clockwise wire is not");
    expect(singleWireInsideOut(bad), "inside out: a clockwise single wire is");
    TopoDS_Face badR = TopoDS::Face(bad.Reversed());
    expect(singleWireInsideOut(badR), "inside out: read on the FORWARD face whatever the flag");
    // The reversal: same loop, same edges, same orientation flag; now valid by BRepCheck's rule.
    TopoDS_Face fixedF = bad;
    expect(reverseInsideOutWire(fixedF), "inside out: the clockwise wire is reversed");
    expect(!singleWireInsideOut(fixedF), "inside out: the reversed wire leaves infinity OUT");
    expect(fixedF.Orientation() == bad.Orientation(), "inside out: the face keeps its flag");
    expect(edgeCount(fixedF) == edgeCount(bad), "inside out: the same edges");
    {
        TopTools_IndexedMapOfShape a, b;
        TopExp::MapShapes(bad, TopAbs_EDGE, a);
        TopExp::MapShapes(fixedF, TopAbs_EDGE, b);
        bool same = a.Extent() == b.Extent();
        for (int i = 1; same && i <= a.Extent(); ++i) same = b.Contains(a(i));
        expect(same, "inside out: every edge keeps its TShape");
    }
    expect(BRepCheck_Analyzer(fixedF).IsValid() == Standard_True, "inside out: BRepCheck accepts it");
    TopoDS_Face goodCopy = good;
    expect(!reverseInsideOutWire(goodCopy) && goodCopy.IsEqual(good),
           "inside out: a healthy face is left untouched");
    TopoDS_Face badRCopy = badR;
    expect(reverseInsideOutWire(badRCopy) && badRCopy.Orientation() == TopAbs_REVERSED &&
               !singleWireInsideOut(badRCopy),
           "inside out: a REVERSED face is fixed on its FORWARD reading");
    // Two wires (a hole) are not this rule's: classification needs both.
    {
        const std::vector<gp_Pnt2d> in = {gp_Pnt2d(1.0, 0.5), gp_Pnt2d(1.0, 0.9), gp_Pnt2d(1.4, 0.9),
                                          gp_Pnt2d(1.4, 0.5)};
        TopoDS_Face two = faceThroughUV(S, cw);
        TopoDS_Face hole = faceThroughUV(S, in);
        BRep_Builder B;
        for (TopExp_Explorer wx(hole, TopAbs_WIRE); wx.More(); wx.Next()) B.Add(two, wx.Current());
        expect(!singleWireInsideOut(two), "inside out: a two-wire face is not read");
    }
    expect(!singleWireInsideOut(TopoDS_Face()), "inside out: a null face is not");
}

int main() {
    caseInsideOut();
    caseEdgeBound();
    const gp_Ax1 fam(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));

    // axisOnFamilyAxis: line distance and direction error times lever, each at tau.
    expect(axisOnFamilyAxis(gp_Ax1(gp_Pnt(0, 0, 7), gp_Dir(0, 0, 1)), fam, 50.0, kTau),
           "axis: same line elsewhere along it");
    expect(axisOnFamilyAxis(gp_Ax1(gp_Pnt(0.5 * kTau, 0, 0), gp_Dir(0, 0, 1)), fam, 50.0, kTau),
           "axis: offset 0.5 tau is on");
    expect(!axisOnFamilyAxis(gp_Ax1(gp_Pnt(2.0 * kTau, 0, 0), gp_Dir(0, 0, 1)), fam, 50.0, kTau),
           "axis: offset 2 tau is off");
    expect(axisOnFamilyAxis(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1)), fam, 50.0, kTau),
           "axis: antiparallel is the same axis");
    {
        const double lever = 40.0;
        const double aIn = 0.5 * kTau / lever, aOut = 2.0 * kTau / lever;
        expect(axisOnFamilyAxis(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(std::sin(aIn), 0, std::cos(aIn))), fam,
                                lever, kTau),
               "axis: tilt x lever 0.5 tau is on");
        expect(!axisOnFamilyAxis(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(std::sin(aOut), 0, std::cos(aOut))),
                                 fam, lever, kTau),
               "axis: tilt x lever 2 tau is off");
    }
    expect(!axisOnFamilyAxis(fam, fam, 1.0, 0.0), "axis: tau 0 admits nothing");

    // circleOnFamilyAxis: the lever is the circle's own radius.
    {
        const double R = 10.0;
        const double aIn = 0.5 * kTau / R, aOut = 2.0 * kTau / R;
        const gp_Circ cIn(gp_Ax2(gp_Pnt(0, 0, 3), gp_Dir(std::sin(aIn), 0, std::cos(aIn))), R);
        const gp_Circ cOut(gp_Ax2(gp_Pnt(0, 0, 3), gp_Dir(std::sin(aOut), 0, std::cos(aOut))), R);
        expect(circleOnFamilyAxis(cIn, fam, kTau), "circle: tilt within tau at R");
        expect(!circleOnFamilyAxis(cOut, fam, kTau), "circle: tilt beyond tau at R");
    }

    // frameOnRefX: X follows refX projected normal to the axis, handedness kept.
    {
        const gp_Ax3 direct(gp_Pnt(1, 2, 3), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        gp_Ax3 indirect = direct;
        indirect.YReverse();
        const gp_Dir refX(std::cos(0.7), std::sin(0.7), 0.3);
        gp_Ax3 t1, t2;
        const bool ok1 = frameOnRefX(direct, refX, t1);
        const bool ok2 = frameOnRefX(indirect, refX, t2);
        const gp_Dir want(std::cos(0.7), std::sin(0.7), 0.0);
        expect(ok1 && t1.XDirection().IsEqual(want, 1e-12), "frame: X is refX projected");
        expect(ok1 && t1.Direction().IsEqual(direct.Direction(), 1e-15) &&
                   t1.Location().IsEqual(direct.Location(), 0.0),
               "frame: axis line unchanged");
        expect(ok1 && ok2 && t1.Direct() && !t2.Direct(), "frame: handedness kept");
        gp_Ax3 t3;
        expect(!frameOnRefX(direct, gp_Dir(0, 0, 1), t3), "frame: refX along the axis refused");
    }

    // frameNeedsTurn: the angle acts on the arm.
    {
        const gp_Ax3 a(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const double arm = 25.0;
        const double dIn = 0.5 * kTau / arm, dOut = 2.0 * kTau / arm;
        const gp_Ax3 bIn(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(std::cos(dIn), std::sin(dIn), 0));
        const gp_Ax3 bOut(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(std::cos(dOut), std::sin(dOut), 0));
        expect(!frameNeedsTurn(a, bIn, arm, kTau), "turn: within tau at the arm stays");
        expect(frameNeedsTurn(a, bOut, arm, kTau), "turn: beyond tau at the arm turns");
        expect(!frameNeedsTurn(a, bOut, arm, 0.0), "turn: tau 0 never turns");
    }

    // seamHalfPlaneOn: the exact circle point on the half-plane.
    {
        const gp_Circ c(gp_Ax2(gp_Pnt(0, 0, 5), gp_Dir(0, 0, 1), gp_Dir(0, 1, 0)), 12.0);
        SeamHalfPlane h;
        const bool ok = seamHalfPlaneOn(c, gp_Dir(1, 0, 0.4), h);
        expect(ok && h.p.IsEqual(gp_Pnt(12.0, 0, 5), 1e-12), "half-plane: P = centre + R x");
        expect(ok && std::fabs(h.x.Dot(gp_Vec(0, 0, 1))) < 1e-15 &&
                   h.y.IsEqual(gp_Vec(0, 1, 0), 1e-12, 1e-12),
               "half-plane: x in the circle plane, y = n x x");
        SeamHalfPlane h2;
        expect(!seamHalfPlaneOn(c, gp_Dir(0, 0, 1), h2), "half-plane: refX along the normal refused");
    }

    // nearestToSeamHalfPlane: smallest |azimuth|, first on ties.
    {
        const gp_Circ c(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 10.0);
        SeamHalfPlane h;
        seamHalfPlaneOn(c, gp_Dir(1, 0, 0), h);
        std::vector<gp_Pnt> pts = {onCircle(10, 0.3, 0), onCircle(10, -0.02, 0), onCircle(10, 0.05, 0),
                                   onCircle(10, kPi, 0)};
        expect(nearestToSeamHalfPlane(c, h, pts) == 1, "nearest: -0.02 beats 0.05 and 0.3");
        std::vector<gp_Pnt> tie = {onCircle(10, 0.04, 0), onCircle(10, -0.04, 0)};
        expect(nearestToSeamHalfPlane(c, h, tie) == 0, "nearest: tie keeps the first");
        expect(nearestToSeamHalfPlane(c, h, {}) == -1, "nearest: empty rim");
    }

    // rimSeamOffset: the chord to the half-plane at the nearest point's own radius.
    {
        const gp_Dir refX(1, 0, 0);
        std::vector<gp_Pnt> on = {onCircle(20, 0.0, 4), onCircle(20, 1.0, 4)};
        expect(rimSeamOffset(fam, refX, on) < 1e-12, "offset: a point on the half-plane is 0");
        const double d = 0.01;
        std::vector<gp_Pnt> off = {onCircle(20, 2.0, 4), onCircle(20, -d, 4)};
        expect(std::fabs(rimSeamOffset(fam, refX, off) - 2.0 * 20.0 * std::sin(0.5 * d)) < 1e-12,
               "offset: 2 rho sin(du/2) of the nearest point");
        const double s = 0.4 * kTau / 20.0;
        std::vector<gp_Pnt> within = {onCircle(20, s, 4)};
        expect(rimSeamOffset(fam, refX, within) <= kTau, "offset: an aligned rim is within tau");
    }

    // rimBuiltBeforeTorus: closed-360 cylinders build first, then region id order.
    expect(rimBuiltBeforeTorus(true, 40, 7), "built-before: a closed cylinder builds first");
    expect(rimBuiltBeforeTorus(false, 3, 7), "built-before: a lower id builds first");
    expect(!rimBuiltBeforeTorus(false, 9, 7), "built-before: a higher id builds after");

    // rimRunsWithSurfaceU: a rim circle coaxial with the torus, axis along or against it.
    {
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        Handle(Geom_Surface) tor = new Geom_ToroidalSurface(ax, 14.0, 3.0);
        const double v = kPi;  // inner equator, rho = 11
        Handle(Geom_Curve) with = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 11.0);
        Handle(Geom_Curve) against = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0)), 11.0);
        expect(rimRunsWithSurfaceU(with, 0.0, tor, 0.0, v), "rim sense: coaxial circle runs with u");
        expect(!rimRunsWithSurfaceU(against, 0.0, tor, 0.0, v), "rim sense: antiparallel circle runs -u");
        expect(rimRunsWithSurfaceU(Handle(Geom_Curve)(), 0.0, tor, 0.0, v), "rim sense: no curve keeps +u");
        expect(rimRunsWithSurfaceU(against, 0.0, Handle(Geom_Surface)(), 0.0, v), "rim sense: no surface keeps +u");
    }

    // splitCapHealApplies: certified at tau and one exact rim, both required.
    expect(splitCapHealApplies(0.5 * kTau, kTau, true), "split-cap: certified wall with an exact rim");
    expect(!splitCapHealApplies(2.0 * kTau, kTau, true), "split-cap: uncertified wall keeps its path");
    expect(!splitCapHealApplies(0.5 * kTau, kTau, false), "split-cap: all-polyline rims keep their path");
    expect(!splitCapHealApplies(0.0, 0.0, true), "split-cap: tau 0 admits nothing");
    expect(!splitCapHealApplies(0.0, -kTau, true), "split-cap: tau < 0 admits nothing");
    expect(splitCapHealApplies(kTau, kTau, true), "split-cap: every vertex within tau, bound inclusive");
    expect(!splitCapHealApplies(std::nextafter(kTau, 1.0), kTau, true), "split-cap: just past tau refused");

    // crossingCapApplies (D-fc-5.7 O-1): interior certified, a rim vertex off tau, one
    // exact rim and at least one crossing chain -- all four required.
    expect(crossingCapApplies(0.5 * kTau, 300.0 * kTau, kTau, true, 15),
           "crossing cap: rim-off wall with an exact rim and a crossing rim");
    expect(!crossingCapApplies(2.0 * kTau, 300.0 * kTau, kTau, true, 15),
           "crossing cap: an interior vertex off tau keeps the path");
    expect(!crossingCapApplies(0.5 * kTau, 0.5 * kTau, kTau, true, 15),
           "crossing cap: a wall certified on every vertex is the split-cap heal's");
    expect(!crossingCapApplies(0.5 * kTau, 300.0 * kTau, kTau, false, 15),
           "crossing cap: no exact rim keeps the path");
    expect(!crossingCapApplies(0.5 * kTau, 300.0 * kTau, kTau, true, 0),
           "crossing cap: no crossing chain keeps the path");
    expect(!crossingCapApplies(0.0, 1.0, 0.0, true, 15), "crossing cap: tau 0 admits nothing");

    // capPointOnGenerator: where a generator meets a cap's plane -- an oblique cut
    // (normal (1, 0, 2)) through a cylinder of radius 4 about z, generator at 60 degrees.
    {
        const gp_Ax2 cap(gp_Pnt(0, 0, 3), gp_Dir(1.0 / std::sqrt(5.0), 0, 2.0 / std::sqrt(5.0)));
        const double R = 4.0, u = kPi / 3.0;
        const gp_Lin gen(onCircle(R, u, -40.0), gp_Dir(0, 0, 1));
        gp_Pnt x;
        const bool ok = capPointOnGenerator(cap, gen, x);
        const double zWant = 3.0 - 0.5 * R * std::cos(u);  // n.(x - o) = 0, n ~ (1, 0, 2)
        expect(ok && std::fabs(x.Z() - zWant) < 1e-12 && gen.Distance(x) < 1e-12,
               "cap point: on the generator and on the cap plane");
        const gp_Ax2 side(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0));
        expect(!capPointOnGenerator(side, gp_Lin(gp_Pnt(0, 2, 0), gp_Dir(0, 0, 1)), x),
               "cap point: a generator parallel to the plane meets it nowhere");
    }

    // capLoopUnshareable: every chain of the cap loop must border a face that can carry
    // the rim circle; only the first chain can hand over a built rim edge.
    {
        const std::vector<int> buildable = {3, 5};  // chain ids whose partner can carry it
        auto ok = [&](int ci) {
            for (int b : buildable)
                if (b == ci) return true;
            return false;
        };
        expect(!capLoopUnshareable(std::vector<int>{3}, false, ok), "cap loop: one buildable chain");
        expect(capLoopUnshareable(std::vector<int>{4}, false, ok), "cap loop: one facet partner");
        expect(!capLoopUnshareable(std::vector<int>{4}, true, ok),
               "cap loop: a reused front edge is shared");
        expect(capLoopUnshareable(std::vector<int>{3, 4}, false, ok),
               "cap loop: a facet partner after the first chain");
        expect(capLoopUnshareable(std::vector<int>{4, 3}, true, ok) == false,
               "cap loop: reused front, buildable rest");
        expect(capLoopUnshareable(std::vector<int>{4, 4}, true, ok),
               "cap loop: a reused front does not excuse a later facet partner");
        expect(capLoopUnshareable(std::vector<int>{3, 5, 7}, false, ok),
               "cap loop: the last chain decides too");
        expect(!capLoopUnshareable(std::vector<int>{3, -1, 5}, false, ok),
               "cap loop: a chain with no index is skipped");
        expect(!capLoopUnshareable(std::vector<int>{}, false, ok), "cap loop: no chain, nothing to share");
    }

    // torusNamingEnd / torusRoundDecided: a failed-torus round stops after the last
    // region (in build order) that can name a torus.
    expect(torusNamingEnd(std::vector<char>{}) == 0, "naming end: empty order");
    expect(torusNamingEnd(std::vector<char>{0, 0, 0}) == 0, "naming end: no naming region");
    expect(torusNamingEnd(std::vector<char>{0, 1, 0, 1, 0}) == 4, "naming end: one past the last");
    expect(torusNamingEnd(std::vector<char>{1}) == 1, "naming end: a single torus");
    expect(!torusRoundDecided(false, 9, 4), "round decided: no failed torus builds on");
    expect(!torusRoundDecided(true, 3, 4), "round decided: a naming region is still unbuilt");
    expect(torusRoundDecided(true, 4, 4), "round decided: every naming region built");
    expect(torusRoundDecided(true, 7, 4), "round decided: later regions cannot name a torus");

    // closedTorusPartnerFailed: a failed region beside a closed torus built whole
    // decides the round like a failed torus.
    {
        auto besideBuilt = [](int id) { return id == 5 || id == 9; };
        expect(!closedTorusPartnerFailed(std::vector<int>{}, besideBuilt),
               "partner failed: no failed region");
        expect(!closedTorusPartnerFailed(std::vector<int>{1, 2, 7}, besideBuilt),
               "partner failed: no failed region beside a built closed torus");
        expect(closedTorusPartnerFailed(std::vector<int>{1, 9}, besideBuilt),
               "partner failed: the last failed region is a cap partner");
        expect(closedTorusPartnerFailed(std::vector<int>{5}, besideBuilt),
               "partner failed: a single cap partner");
        const bool partner = closedTorusPartnerFailed(std::vector<int>{3, 5}, besideBuilt);
        expect(!torusRoundDecided(false || partner, 3, 4),
               "partner failed: a naming region is still unbuilt");
        expect(torusRoundDecided(false || partner, 4, 4),
               "partner failed: every naming region built");
        expect(!torusRoundDecided(false || closedTorusPartnerFailed(std::vector<int>{3}, besideBuilt),
                                  9, 4),
               "partner failed: no partner, no failed torus, the round builds on");
    }

    // closedWallWithWindows: two rims winding once in opposite senses, every other
    // loop a window that does not wind.
    expect(closedWallWithWindows(std::vector<int>{1, -1}), "wall: two rims, no window");
    expect(closedWallWithWindows(std::vector<int>{0, -1, 0, 1}), "wall: two rims and two windows");
    expect(!closedWallWithWindows(std::vector<int>{0}), "wall: a disk is no wall");
    expect(!closedWallWithWindows(std::vector<int>{0, 0}), "wall: two windows and no rim");
    expect(!closedWallWithWindows(std::vector<int>{1, 1}), "wall: two rims in one sense");
    expect(!closedWallWithWindows(std::vector<int>{1, -1, 1}), "wall: a third winding loop");
    expect(!closedWallWithWindows(std::vector<int>{2, -2}), "wall: a loop winding twice");
    expect(!closedWallWithWindows(std::vector<int>{}), "wall: no loop");

    // windowGap: the complement of the largest gap of a window's own azimuths.
    {
        double lo = 0, hi = 0;
        expect(windowGap(std::vector<double>{0.4, 0.5, 0.6}, lo, hi) &&
                   std::fabs(lo - 0.6) < 1e-12 && std::fabs(hi - (0.4 + 2 * kPi)) < 1e-12,
               "window gap: a window at 0.4..0.6 leaves 0.6..0.4+2pi");
        expect(windowGap(std::vector<double>{kPi - 0.1, -kPi + 0.1}, lo, hi) &&
                   std::fabs(lo - (-kPi + 0.1)) < 1e-12 && std::fabs(hi - (kPi - 0.1)) < 1e-12,
               "window gap: a window across +-pi leaves the side opposite");
        expect(!windowGap(std::vector<double>{1.0}, lo, hi), "window gap: one azimuth is no window");
    }

    // seamClearOfWindows: the seam stands in the widest gap of the pooled windows and
    // must lie inside every window's own gap.
    {
        double u = -1;
        expect(seamClearOfWindows(std::vector<std::vector<double>>{}, u) && u == 0.0,
               "seam: no window, u = 0");
        // Two slot windows facing each other (35 deg wide about +-pi/2): the seam stands
        // between them, at 0 or pi.
        const double h = 0.3;
        const std::vector<std::vector<double>> slots = {
            {0.5 * kPi - h, 0.5 * kPi, 0.5 * kPi + h}, {-0.5 * kPi - h, -0.5 * kPi, -0.5 * kPi + h}};
        expect(seamClearOfWindows(slots, u) &&
                   (std::fabs(u) < 1e-9 || std::fabs(std::fabs(u) - kPi) < 1e-9),
               "seam: two opposed windows, seam between them");
        // The builder stands the seam in the widest POOLED gap; when that gap is the
        // inside of a sparsely sampled window the seam would cross it: refused.
        const std::vector<std::vector<double>> blocked = {
            {-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0}, {1.5, 2.0, 2.5, 3.0}};
        expect(!seamClearOfWindows(blocked, u), "seam: widest pooled gap inside a window");
        // Windows that cover every generator between them.
        const std::vector<std::vector<double>> full = {{-3.0, -1.0, 0.5}, {0.4, 2.0, 3.3}};
        expect(!seamClearOfWindows(full, u), "seam: windows that cover every generator");
        // Windows side by side: the pooled widest gap is clear of both.
        const std::vector<std::vector<double>> side = {{0.1, 0.2}, {0.3, 0.4}};
        expect(seamClearOfWindows(side, u) && std::fabs(u - std::remainder(0.25 + kPi, 2 * kPi)) < 1e-9,
               "seam: adjacent windows, seam opposite them");
        expect(!seamClearOfWindows(std::vector<std::vector<double>>{{0.3}}, u),
               "seam: a one-azimuth window is refused");
    }

    // closedLoopFromTerminal: the loop starts at its terminal and closes on it.
    expect(closedLoopFromTerminal(std::vector<int>{7, 8, 9, 4, 5}, 4) == std::vector<int>({4, 5, 7, 8, 9, 4}),
           "rim loop: rotated to its terminal");
    expect(closedLoopFromTerminal(std::vector<int>{7, 8, 9}, 7) == std::vector<int>({7, 8, 9, 7}),
           "rim loop: terminal first is unchanged");
    expect(closedLoopFromTerminal(std::vector<int>{7, 8, 9}, 3) == std::vector<int>({7, 8, 9, 7}),
           "rim loop: a terminal off the loop keeps the order");
    expect(closedLoopFromTerminal(std::vector<int>{}, 3).empty(), "rim loop: empty stays empty");

    // capPathsRunOpposite (D-fc-5.7 O-4): the plain split-cap wire closes in (u, v) only
    // when both cap paths sweep the same way round the axis.
    expect(!capPathsRunOpposite(-2.0 * kPi, -2.0 * kPi), "cap paths: both -u close the plain wire");
    expect(!capPathsRunOpposite(2.0 * kPi, 2.0 * kPi), "cap paths: both +u close the plain wire");
    expect(capPathsRunOpposite(-2.0 * kPi, 2.0 * kPi), "cap paths: -u high, +u low jump a period");
    expect(capPathsRunOpposite(2.0 * kPi, -2.0 * kPi), "cap paths: +u high, -u low jump a period");
    expect(!capPathsRunOpposite(0.0, 2.0 * kPi), "cap paths: a path with no sweep decides nothing");
    expect(!capPathsRunOpposite(-2.0 * kPi, 0.0), "cap paths: no sweep on the low path either");

    // D-fc-5.13: edgeSampledSignedAreaUV reads every curve under its orientation flag;
    // closedEllipseHalves / closedAnalyticRim ship and recognise a closed ellipse rim as
    // two half-turn arcs. An ellipse 5 x 3 about z on the z = 0 plane.
    {
        const double a = 5.0, b = 3.0;
        const gp_Elips el(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), a, b);
        const gp_Pln up(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)));
        const gp_Pln down(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0)));
        const TopoDS_Edge full = BRepBuilderAPI_MakeEdge(el).Edge();
        const TopoDS_Edge fullR = TopoDS::Edge(full.Reversed());
        // 32 samples over a full turn: the inscribed 32-gon, (32 / 2) a b sin(2 pi / 32).
        const double gon = 16.0 * a * b * std::sin(2.0 * kPi / 32.0);
        const double aF = edgeSampledSignedAreaUV(full, up);
        expect(std::fabs(aF - gon) < 1e-9, "sampled area: closed ellipse FORWARD, counter-clockwise");
        expect(std::fabs(edgeSampledSignedAreaUV(fullR, up) + aF) < 1e-12,
               "sampled area: the same ellipse REVERSED is clockwise");
        expect(std::fabs(edgeSampledSignedAreaUV(full, down) + aF) < 1e-12,
               "sampled area: read about the opposite normal");
        const TopoDS_Edge line =
            BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(4, 1, 0)).Edge();
        expect(std::fabs(edgeSampledSignedAreaUV(line, up)) < 1e-12 &&
                   std::fabs(edgeSampledSignedAreaUV(TopoDS::Edge(line.Reversed()), up)) < 1e-12,
               "sampled area: a line encloses nothing either way");

        const std::vector<TopoDS_Edge> h = closedEllipseHalves(full);
        bool shape = h.size() == 2;
        TopoDS_Vertex V0, V1, A0, A1;
        if (shape) {
            TopExp::Vertices(full, V0, V1, Standard_True);
            TopoDS_Vertex h0f, h0l, h1f, h1l;
            TopExp::Vertices(h[0], h0f, h0l, Standard_True);
            TopExp::Vertices(h[1], h1f, h1l, Standard_True);
            Standard_Real f0 = 0, l0 = 0, f1 = 0, l1 = 0;
            const Handle(Geom_Curve) c0 = BRep_Tool::Curve(h[0], f0, l0);
            const Handle(Geom_Curve) c1 = BRep_Tool::Curve(h[1], f1, l1);
            shape = h0f.IsSame(V0) && h1l.IsSame(V0) && h0l.IsSame(h1f) && !h0l.IsSame(V0) &&
                    std::fabs((l0 - f0) - kPi) < 1e-12 && std::fabs((l1 - f1) - kPi) < 1e-12 &&
                    !Handle(Geom_Ellipse)::DownCast(c0).IsNull() && c0 == c1;
            A0 = h0l;
        }
        expect(shape, "halves: vertex to antipode and back, half a turn each, one ellipse");
        expect(!A0.IsNull() && BRep_Tool::Pnt(A0).Distance(gp_Pnt(-a, 0, 0)) < 1e-12,
               "halves: the antipode is the exact curve point half a turn on");
        expect(h.size() == 2 && edgeSampledSignedAreaUV(h[0], up) > 0.0 &&
                   std::fabs(edgeSampledSignedAreaUV(h[0], up) -
                             edgeSampledSignedAreaUV(h[1], up)) < 1e-9,
               "halves: each half runs with the whole, and they are equal");
        const std::vector<TopoDS_Edge> hR = closedEllipseHalves(fullR);
        bool rev = hR.size() == 2;
        if (rev) {
            TopoDS_Vertex r0f, r0l, r1f, r1l;
            TopExp::Vertices(hR[0], r0f, r0l, Standard_True);
            TopExp::Vertices(hR[1], r1f, r1l, Standard_True);
            rev = r0f.IsSame(V0) && r1l.IsSame(V0) && r0l.IsSame(r1f) &&
                  edgeSampledSignedAreaUV(hR[0], up) < 0.0 && edgeSampledSignedAreaUV(hR[1], up) < 0.0;
        }
        expect(rev, "halves: a REVERSED ellipse keeps its traversal (clockwise halves)");
        // The antipode a caller already has a vertex for (the chain's own mesh vertex):
        // the halves meet AT that vertex, on its own parameter, and mint nothing.
        gp_Pnt pA;
        expect(closedEllipseAntipode(full, pA) && pA.Distance(gp_Pnt(-a, 0, 0)) < 1e-12,
               "antipode: half a turn from the edge's vertex");
        {
            TopoDS_Vertex mesh;
            BRep_Builder().MakeVertex(mesh, gp_Pnt(-a, 2e-7, 0), 1e-6);
            const std::vector<TopoDS_Edge> hm = closedEllipseHalves(full, mesh);
            bool meet = hm.size() == 2;
            if (meet) {
                TopoDS_Vertex m0f, m0l, m1f, m1l;
                TopExp::Vertices(hm[0], m0f, m0l, Standard_True);
                TopExp::Vertices(hm[1], m1f, m1l, Standard_True);
                Standard_Real f0 = 0, l0 = 0;
                BRep_Tool::Range(hm[0], f0, l0);
                meet = m0l.IsSame(mesh) && m1f.IsSame(mesh) && m0f.IsSame(V0) && m1l.IsSame(V0) &&
                       std::fabs((l0 - f0) - kPi) < 1e-6 && closedAnalyticRim(hm);
            }
            expect(meet, "halves: meet at the caller's vertex standing at the antipode");
            TopoDS_Vertex v0f = V0;
            expect(closedEllipseHalves(full, v0f).empty(),
                   "halves: the edge's own vertex cannot be its antipode");
        }
        const TopoDS_Edge circ =
            BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 4.0)).Edge();
        expect(closedEllipseHalves(circ).empty(), "halves: a circle ships whole");
        expect(!closedEllipseAntipode(circ, pA), "antipode: a circle has none here");
        const TopoDS_Edge arc = BRepBuilderAPI_MakeEdge(el, 0.0, kPi).Edge();
        expect(closedEllipseHalves(arc).empty(), "halves: an ellipse arc is not a closed rim");
        expect(closedEllipseHalves(TopoDS_Edge()).empty(), "halves: a null edge");

        expect(closedAnalyticRim({full}), "rim: one closed ellipse");
        expect(closedAnalyticRim({circ}), "rim: one closed circle");
        expect(closedAnalyticRim(h) && closedAnalyticRim(hR), "rim: the two halves, either sense");
        expect(h.size() == 2 && closedAnalyticRim({h[1], h[0]}), "rim: the halves from the antipode");
        expect(!closedAnalyticRim({arc}), "rim: an open arc");
        expect(h.size() == 2 && !closedAnalyticRim({h[0]}), "rim: one half alone");
        const std::vector<TopoDS_Edge> hOther = closedEllipseHalves(BRepBuilderAPI_MakeEdge(
            gp_Elips(gp_Ax2(gp_Pnt(0, 0, 1), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), a, b)).Edge());
        expect(h.size() == 2 && hOther.size() == 2 && !closedAnalyticRim({h[0], hOther[1]}),
               "rim: halves of two ellipses");
        const TopoDS_Edge c1 = BRepBuilderAPI_MakeEdge(
            gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 4.0), 0.0, kPi).Edge();
        TopoDS_Vertex cv0, cv1;
        TopExp::Vertices(c1, cv0, cv1, Standard_True);
        const TopoDS_Edge c2 = BRepBuilderAPI_MakeEdge(
            new Geom_Circle(gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 4.0)), cv1, cv0, kPi,
            2.0 * kPi).Edge();
        expect(!closedAnalyticRim({c1, c2}), "rim: two circle arcs are not this construction");
        expect(h.size() == 2 && !closedAnalyticRim({h[0], h[1], full}), "rim: three edges");
        expect(!closedAnalyticRim({}), "rim: no edge");
    }

    // D-fc-5.20 B-1. The review's closed elliptical rim: radii 5 and 3, eight
    // vertices on the ellipse, chain start at the major-axis vertex, the opposite
    // vertex exactly at the antipode. q ≤ 0 and a non-finite q refuse the split.
    {
        const double a = 5.0, b = 3.0;
        const gp_Elips el(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), a, b);
        const TopoDS_Edge full = BRepBuilderAPI_MakeEdge(el).Edge();
        TopoDS_Vertex edgeV0, edgeV1;
        TopExp::Vertices(full, edgeV0, edgeV1, Standard_True);
        expect(!edgeV0.IsNull() && BRep_Tool::Pnt(edgeV0).Distance(gp_Pnt(a, 0, 0)) < 1e-12,
               "b1 mesh: the edge starts at the major-axis vertex");
        gp_Pnt pAnti;
        expect(closedEllipseAntipode(full, pAnti), "b1 mesh: the antipode exists");
        std::vector<TopoDS_Vertex> verts(8);
        std::vector<int> ids(8);
        BRep_Builder bb;
        for (int k = 0; k < 8; ++k) {
            const double t = k * (2.0 * kPi / 8.0);
            const gp_Pnt p = (k == 4) ? pAnti : ElCLib::Value(t, el);
            bb.MakeVertex(verts[(size_t)k], p, Precision::Confusion());
            ids[(size_t)k] = k;
        }
        const int iTerm = 0;
        const int iAnti = 4;
        expect(BRep_Tool::Pnt(verts[(size_t)iAnti]).Distance(pAnti) == 0.0,
               "b1 mesh: the opposite vertex is the antipode");
        expect(BRep_Tool::Pnt(verts[(size_t)iTerm]).Distance(gp_Pnt(a, 0, 0)) < 1e-9,
               "b1 mesh: the chain starts at the major-axis vertex");

        const double refused[] = {0.0, -0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::infinity(),
                                  -std::numeric_limits<double>::infinity()};
        for (double q : refused) {
            char name[80];
            std::snprintf(name, sizeof name, "b1: q=%g admits nothing", q);
            expect(!ellipseRimSplitAdmitted(q), name);
            std::snprintf(name, sizeof name, "b1: q=%g mints nothing", q);
            expect(closedEllipseChainSplit(full, q, TopoDS_Vertex()).empty(), name);
            std::snprintf(name, sizeof name, "b1: q=%g refuses the split with the antipode in hand", q);
            expect(closedEllipseChainSplit(full, q, verts[(size_t)iAnti]).empty(), name);
        }
        // Asking the constructor directly still mints. That mint is a second
        // vertex at the mesh vertex's place: the twin the gate refuses.
        {
            const std::vector<TopoDS_Edge> minted = closedEllipseHalves(full);
            bool twin = minted.size() == 2;
            if (twin) {
                TopoDS_Vertex m0f, m0l;
                TopExp::Vertices(minted[0], m0f, m0l, Standard_True);
                twin = !m0l.IsSame(verts[(size_t)iAnti]) &&
                       BRep_Tool::Pnt(m0l).Distance(BRep_Tool::Pnt(verts[(size_t)iAnti])) < 1e-12;
            }
            expect(twin, "b1 witness: an unguarded mint twins the antipode vertex");
        }

        const double q = 1.0e-6;
        const double tau = 2.0 * q;
        expect(ellipseRimSplitAdmitted(q), "b1: a positive finite q admits");
        const int lv = closedEllipseAntipodePick(full, verts, ids, iTerm, tau);
        expect(lv == iAnti, "b1: positive q picks the mesh antipode");
        expect(closedEllipseAntipodePick(full, verts, ids, iTerm, 0.0) < 0, "b1: tau 0 picks nothing");
        expect(closedEllipseAntipodePick(full, verts, ids, iTerm,
                                         std::numeric_limits<double>::quiet_NaN()) < 0,
               "b1: NaN tau picks nothing");
        const std::vector<TopoDS_Edge> h =
            (lv >= 0) ? closedEllipseChainSplit(full, q, verts[(size_t)lv]) : std::vector<TopoDS_Edge>{};
        bool meet = h.size() == 2;
        if (meet) {
            TopoDS_Vertex h0f, h0l, h1f, h1l;
            TopExp::Vertices(h[0], h0f, h0l, Standard_True);
            TopExp::Vertices(h[1], h1f, h1l, Standard_True);
            Standard_Real f0 = 0, l0 = 0, f1 = 0, l1 = 0;
            BRep_Tool::Range(h[0], f0, l0);
            BRep_Tool::Range(h[1], f1, l1);
            meet = h0l.IsSame(verts[(size_t)iAnti]) && h1f.IsSame(verts[(size_t)iAnti]) &&
                   std::fabs((l0 - f0) - kPi) < 1e-8 && std::fabs((l1 - f1) - kPi) < 1e-8;
        }
        expect(meet, "b1: positive q meets at the mesh vertex, half a turn each");

        std::vector<TopoDS_Vertex> far(1);
        bb.MakeVertex(far[0], gp_Pnt(0, 0, 40), Precision::Confusion());
        expect(closedEllipseAntipodePick(full, far, std::vector<int>{0}, -1, tau) < 0,
               "b1: a vertex outside tau is not the antipode");
        expect(closedEllipseChainSplit(full, q, TopoDS_Vertex()).size() == 2,
               "b1: positive q with no in-range vertex still mints the exact point");
    }

    // c1 (D-130-16, crossing-cap rim): an exact rim point is the chain's own mesh
    // vertex when one stands within tau of it. The rim's mesh vertices sit off the
    // fitted ellipse by a fraction of tau (the fit's residual); a vertex minted at
    // the exact point would be a second vertex at that mesh vertex's place.
    {
        const double a = 6.0, b = 4.0;
        const gp_Elips el(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), a, b);
        const double q = 1.0e-5;
        const double tau = 2.0 * q;
        const int n = 12;
        std::vector<TopoDS_Vertex> verts((size_t)n);
        std::vector<int> ids((size_t)n);
        BRep_Builder bb;
        for (int k = 0; k < n; ++k) {
            const double t = k * (2.0 * kPi / n);
            const gp_Pnt p = ElCLib::Value(t, el).Translated(gp_Vec(0, 0, 0.4 * tau));
            bb.MakeVertex(verts[(size_t)k], p, Precision::Confusion());
            ids[(size_t)k] = k;
        }
        const int kE = 3, kA = 9;  // half a turn apart
        const gp_Pnt pE = ElCLib::Value(kE * (2.0 * kPi / n), el);
        const gp_Pnt pA = ElCLib::Value(kA * (2.0 * kPi / n), el);
        expect(chainVertexWithinTau(pE, verts, ids, TopoDS_Vertex(), tau) == kE,
               "c1: the mesh vertex within tau of the seam point is reused");
        const int lvE = chainVertexWithinTau(pE, verts, ids, TopoDS_Vertex(), tau);
        expect(lvE >= 0 && chainVertexWithinTau(pA, verts, ids, verts[(size_t)lvE], tau) == kA,
               "c1: the mesh vertex within tau of the antipode is reused");
        expect(chainVertexWithinTau(pE, verts, ids, verts[(size_t)kE], tau) < 0,
               "c1: the skipped vertex is never its own antipode");
        const gp_Pnt pMid = ElCLib::Value((kE + 0.5) * (2.0 * kPi / n), el);
        expect(chainVertexWithinTau(pMid, verts, ids, TopoDS_Vertex(), tau) < 0,
               "c1: no mesh vertex within tau -> mint (-1)");

        const double refused[] = {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::infinity()};
        for (double t : refused) {
            char name[80];
            std::snprintf(name, sizeof name, "c1: tau=%g reuses nothing", t);
            expect(chainVertexWithinTau(pE, verts, ids, TopoDS_Vertex(), t) < 0, name);
        }

        // Exactly tau away is within tau; just beyond is not.
        std::vector<TopoDS_Vertex> edge(2);
        bb.MakeVertex(edge[0], gp_Pnt(0, 0, tau), Precision::Confusion());
        bb.MakeVertex(edge[1], gp_Pnt(0, 0, -tau * (1.0 + 1e-6)), Precision::Confusion());
        expect(chainVertexWithinTau(gp_Pnt(0, 0, 0), edge, std::vector<int>{1, 0},
                                    TopoDS_Vertex(), tau) == 0,
               "c1: a vertex exactly tau away is within tau");
        expect(chainVertexWithinTau(gp_Pnt(0, 0, 0), edge, std::vector<int>{1},
                                    TopoDS_Vertex(), tau) < 0,
               "c1: a vertex beyond tau is not reused");

        // Nearest wins; an exact tie goes to the later meshVerts entry, so the pick
        // is fixed by the chain's walk, never by container order elsewhere.
        std::vector<TopoDS_Vertex> two(3);
        bb.MakeVertex(two[0], gp_Pnt(0.5 * tau, 0, 0), Precision::Confusion());
        bb.MakeVertex(two[1], gp_Pnt(-0.5 * tau, 0, 0), Precision::Confusion());
        bb.MakeVertex(two[2], gp_Pnt(0, 0.25 * tau, 0), Precision::Confusion());
        expect(chainVertexWithinTau(gp_Pnt(0, 0, 0), two, std::vector<int>{0, 1, 2},
                                    TopoDS_Vertex(), tau) == 2,
               "c1: the nearest vertex within tau wins");
        expect(chainVertexWithinTau(gp_Pnt(0, 0, 0), two, std::vector<int>{0, 1},
                                    TopoDS_Vertex(), tau) == 1 &&
                   chainVertexWithinTau(gp_Pnt(0, 0, 0), two, std::vector<int>{1, 0},
                                        TopoDS_Vertex(), tau) == 0,
               "c1: a tie goes to the later meshVerts entry");

        // Slots outside verts[] or null are not candidates.
        std::vector<TopoDS_Vertex> holes(2);
        bb.MakeVertex(holes[1], gp_Pnt(0, 0, 0), Precision::Confusion());
        expect(chainVertexWithinTau(gp_Pnt(0, 0, 0), holes, std::vector<int>{-1, 0, 7, 1},
                                    TopoDS_Vertex(), tau) == 1,
               "c1: null and out-of-range slots are skipped");

        // closedEllipseAntipodePick is the same rule at the antipode.
        std::vector<TopoDS_Vertex> onEl = verts;
        bb.MakeVertex(onEl[0], ElCLib::Value(0.0, el), Precision::Confusion());
        BRepBuilderAPI_MakeEdge mFull(el, onEl[0], onEl[0]);
        gp_Pnt pAnti;
        if (mFull.IsDone() && closedEllipseAntipode(mFull.Edge(), pAnti)) {
            const int pick = closedEllipseAntipodePick(mFull.Edge(), onEl, ids, 0, tau);
            expect(pick == n / 2 && pick == chainVertexWithinTau(pAnti, onEl, ids, onEl[0], tau),
                   "c1: the antipode pick is chainVertexWithinTau at the antipode");
        } else {
            expect(false, "c1: the full ellipse has an antipode");
        }
    }

    // D-fc-7.2 (4). A partial cylinder whose u-span straddles 0, bounded by the
    // rot-trim / rev-hole edges (circle, two generators, one chord), is one valid
    // face when every edge takes the trim's own u-chart.
    {
        const double R = 5.0;
        const double uLo = -0.45;
        const double uHi = 0.70;
        const double vLo = 0.0;
        const double vHi = 4.0;
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_CylindricalSurface) basis =
            new Geom_CylindricalSurface(gp_Cylinder(ax, R));
        const Handle(Geom_RectangularTrimmedSurface) trim =
            new Geom_RectangularTrimmedSurface(basis, uLo, uHi, vLo, vHi);
        auto at = [&](double u, double v) { return gp_Pnt(R * std::cos(u), R * std::sin(u), v); };
        BRep_Builder bb;
        TopoDS_Vertex v00, v10, v11, v01;
        bb.MakeVertex(v00, at(uLo, vLo), Precision::Confusion());
        bb.MakeVertex(v10, at(uHi, vLo), Precision::Confusion());
        bb.MakeVertex(v11, at(uHi, vHi), Precision::Confusion());
        bb.MakeVertex(v01, at(uLo, vHi), Precision::Confusion());
        const gp_Ax2 cax(gp_Pnt(0, 0, vLo), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_TrimmedCurve) arc = new Geom_TrimmedCurve(new Geom_Circle(cax, R), uLo, uHi);
        const Handle(Geom_Line) genHi = new Geom_Line(at(uHi, vLo), gp_Dir(0, 0, 1));
        const gp_Vec chordDir(at(uHi, vHi), at(uLo, vHi));
        const Handle(Geom_Line) chord = new Geom_Line(at(uHi, vHi), gp_Dir(chordDir));
        const Handle(Geom_Line) genLo = new Geom_Line(at(uLo, vHi), gp_Dir(0, 0, -1));
        BRepBuilderAPI_MakeEdge eArc(arc, v00, v10);
        BRepBuilderAPI_MakeEdge eGenHi(genHi, v10, v11);
        BRepBuilderAPI_MakeEdge eChord(chord, v11, v01);
        BRepBuilderAPI_MakeEdge eGenLo(genLo, v01, v00);
        expect(eArc.IsDone() && eGenHi.IsDone() && eChord.IsDone() && eGenLo.IsDone(),
               "partial chart: rot-trim / rev-hole edges built");
        TopoDS_Edge edges[4] = {eArc.Edge(), eGenHi.Edge(), eChord.Edge(), eGenLo.Edge()};
        bool charts = true;
        for (const TopoDS_Edge& e : edges) {
            Standard_Real f = 0, l = 0;
            Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
            Handle(Geom2d_Curve) c2d = partialTrimChartPCurve(e, c3, f, l, trim);
            if (c2d.IsNull()) {
                charts = false;
                continue;
            }
            Standard_Real buLo = 0, buHi = 0, bvLo = 0, bvHi = 0;
            trim->Bounds(buLo, buHi, bvLo, bvHi);
            double dev = 0.0;
            for (int i = 0; i <= 32; ++i) {
                const double t = f + (l - f) * (static_cast<double>(i) / 32.0);
                const gp_Pnt2d uv = c2d->Value(t);
                if (uv.X() < buLo - Precision::Confusion() || uv.X() > buHi + Precision::Confusion())
                    charts = false;
                dev = std::max(dev, c3->Value(t).Distance(trim->Value(uv.X(), uv.Y())));
            }
            const double tol = std::max(Precision::Confusion(), dev);
            const TopLoc_Location loc;
            bb.UpdateEdge(e, c2d, trim, loc, tol);
            bb.Range(e, trim, loc, f, l);
        }
        expect(charts, "partial chart: every edge lies on the straddling trim");
        // A trim of the same span that does not contain the vertices is not the
        // chart (the rotated-basis trim whose parameter window is the mirror).
        {
            const Handle(Geom_RectangularTrimmedSurface) mirror =
                new Geom_RectangularTrimmedSurface(basis, 0.0, 0.40, vLo, vHi);
            TopoDS_Vertex na, nb;
            bb.MakeVertex(na, at(-0.40, vLo), Precision::Confusion());
            bb.MakeVertex(nb, at(-0.40, vHi), Precision::Confusion());
            const Handle(Geom_Line) genNeg = new Geom_Line(at(-0.40, vLo), gp_Dir(0, 0, 1));
            BRepBuilderAPI_MakeEdge eNeg(genNeg, na, nb);
            expect(eNeg.IsDone(), "partial chart: generator off the mirror trim");
            if (eNeg.IsDone()) {
                Standard_Real f = 0, l = 0;
                Handle(Geom_Curve) c3 = BRep_Tool::Curve(eNeg.Edge(), f, l);
                expect(partialTrimChartPCurve(eNeg.Edge(), c3, f, l, mirror).IsNull(),
                       "partial chart: mirrored trim is refused");
                const Handle(Geom_RectangularTrimmedSurface) own =
                    new Geom_RectangularTrimmedSurface(basis, -0.40, 0.0, vLo, vHi);
                expect(!partialTrimChartPCurve(eNeg.Edge(), c3, f, l, own).IsNull(),
                       "partial chart: the region's own trim is accepted");
            }
        }
        BRepBuilderAPI_MakeWire mkW;
        for (const TopoDS_Edge& e : edges) mkW.Add(e);
        expect(mkW.IsDone(), "partial chart: wire closed");
        if (mkW.IsDone()) {
            BRepBuilderAPI_MakeFace mkF(trim, mkW.Wire());
            expect(mkF.IsDone(), "partial chart: face built");
            if (mkF.IsDone()) {
                BRepCheck_Analyzer an(mkF.Face(), Standard_True);
                expect(an.IsValid() == Standard_True, "partial chart: one valid face");
                const TopoDS_Edge& arcE = edges[0];
                const TopoDS_Edge& genE = edges[1];
                Standard_Real fa = 0, la = 0, fg = 0, lg = 0;
                const TopLoc_Location loc;
                Handle(Geom2d_Curve) pa = BRep_Tool::CurveOnSurface(arcE, trim, loc, fa, la);
                Handle(Geom2d_Curve) pg = BRep_Tool::CurveOnSurface(genE, trim, loc, fg, lg);
                expect(!pa.IsNull() && !pg.IsNull() &&
                           pa->Value(la).Distance(pg->Value(fg)) <= Precision::Confusion(),
                       "partial chart: circle and generator meet at the shared vertex");
            }
        }
    }

    // D-fc-7.10. A partial cylinder that is already valid keeps its own pcurves
    // when the heal runs. The chart object is not written onto the face.
    {
        const double R = 5.0;
        const double uLo = 0.40;
        const double uHi = 1.10;
        const double vLo = 0.0;
        const double vHi = 3.0;
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_CylindricalSurface) basis =
            new Geom_CylindricalSurface(gp_Cylinder(ax, R));
        const Handle(Geom_RectangularTrimmedSurface) trim =
            new Geom_RectangularTrimmedSurface(basis, uLo, uHi, vLo, vHi);
        auto at = [&](double u, double v) { return gp_Pnt(R * std::cos(u), R * std::sin(u), v); };
        BRep_Builder bb;
        const double vt = Precision::Confusion();
        TopoDS_Vertex v00, v10, v11, v01;
        bb.MakeVertex(v00, at(uLo, vLo), vt);
        bb.MakeVertex(v10, at(uHi, vLo), vt);
        bb.MakeVertex(v11, at(uHi, vHi), vt);
        bb.MakeVertex(v01, at(uLo, vHi), vt);
        const gp_Ax2 cax(gp_Pnt(0, 0, vLo), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_TrimmedCurve) arc =
            new Geom_TrimmedCurve(new Geom_Circle(cax, R), uLo, uHi);
        const Handle(Geom_Line) genHi = new Geom_Line(at(uHi, vLo), gp_Dir(0, 0, 1));
        const gp_Vec chordDir(at(uHi, vHi), at(uLo, vHi));
        const Handle(Geom_Line) chord = new Geom_Line(at(uHi, vHi), gp_Dir(chordDir));
        const Handle(Geom_Line) genLo = new Geom_Line(at(uLo, vHi), gp_Dir(0, 0, -1));
        BRepBuilderAPI_MakeEdge eArc(arc, v00, v10);
        BRepBuilderAPI_MakeEdge eGenHi(genHi, v10, v11);
        BRepBuilderAPI_MakeEdge eChord(chord, v11, v01);
        BRepBuilderAPI_MakeEdge eGenLo(genLo, v01, v00);
        expect(eArc.IsDone() && eGenHi.IsDone() && eChord.IsDone() && eGenLo.IsDone(),
               "valid partial: edges built");
        TopoDS_Edge edges[4] = {eArc.Edge(), eGenHi.Edge(), eChord.Edge(), eGenLo.Edge()};
        const TopLoc_Location loc;
        bool placed = true;
        for (const TopoDS_Edge& e : edges) {
            Standard_Real f = 0, l = 0;
            Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
            Handle(Geom2d_Curve) src = partialTrimChartPCurve(e, c3, f, l, trim);
            if (src.IsNull()) {
                placed = false;
                continue;
            }
            Handle(Geom2d_Curve) own = Handle(Geom2d_Curve)::DownCast(src->Copy());
            if (own.IsNull() || own.get() == src.get()) {
                placed = false;
                continue;
            }
            const double dev = maxCurveOnSurfaceDeviation(c3, f, l, trim, own, loc);
            const double tol = std::max(vt, edgeTolFromMeasuredDeviation(dev));
            bb.UpdateEdge(e, own, trim, loc, tol);
            bb.Range(e, trim, loc, f, l);
        }
        expect(placed, "valid partial: own pcurves, not the chart object");
        BRepBuilderAPI_MakeWire mkW;
        for (const TopoDS_Edge& e : edges) mkW.Add(e);
        expect(mkW.IsDone(), "valid partial: wire closed");
        if (mkW.IsDone() && placed) {
            TopoDS_Face face;
            bb.MakeFace(face, trim, vt);
            bb.Add(face, mkW.Wire());
            expectFaceValidBeforeHeal(face, Standard_True, "valid partial: face valid before the heal");
            const Geom2d_Curve* raw[4] = {};
            double tol0[4] = {};
            gp_Pnt2d end0[4][2] = {};
            for (int i = 0; i < 4; ++i) {
                Standard_Real a = 0, b = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(edges[i], face, a, b);
                raw[i] = pc.get();
                tol0[i] = BRep_Tool::Tolerance(edges[i]);
                if (!pc.IsNull()) {
                    end0[i][0] = pc->Value(a);
                    end0[i][1] = pc->Value(b);
                }
            }
            const bool healed = healInvalidPartialTrimChart(face, kTau);
            expect(!healed, "valid partial: the heal does not rebuild a valid face");
            const double vtBefore = BRep_Tool::Tolerance(v00);
            const bool raised = raiseMeasuredCornerTol(face, kTau);
            expect(!raised && std::fabs(BRep_Tool::Tolerance(v00) - vtBefore) <= Precision::Confusion(),
                   "valid partial: a valid face does not take the corner tolerance");
            bool same = true;
            for (int i = 0; i < 4; ++i) {
                Standard_Real a = 0, b = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(edges[i], face, a, b);
                if (pc.get() != raw[i]) same = false;
                if (std::fabs(BRep_Tool::Tolerance(edges[i]) - tol0[i]) > Precision::Confusion())
                    same = false;
                if (pc.IsNull() || pc->Value(a).Distance(end0[i][0]) > Precision::Confusion() ||
                    pc->Value(b).Distance(end0[i][1]) > Precision::Confusion())
                    same = false;
                Standard_Real f = 0, l = 0;
                Handle(Geom_Curve) c3 = BRep_Tool::Curve(edges[i], f, l);
                Handle(Geom2d_Curve) chart = partialTrimChartPCurve(edges[i], c3, f, l, trim);
                if (!chart.IsNull() && chart.get() == pc.get()) same = false;
            }
            BRepCheck_Analyzer after(face, Standard_True);
            expect(same && after.IsValid() == Standard_True,
                   "valid partial: pcurves byte-identical, chart not installed");
        }
    }

    // D-fc-7.10. A partial-cylinder wire that still self-intersects after the
    // chart (the two diagonals cross) keeps today's pcurves.
    {
        const double R = 5.0;
        const double uLo = 0.20;
        const double uHi = 1.00;
        const double vLo = 0.0;
        const double vHi = 2.0;
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_CylindricalSurface) basis =
            new Geom_CylindricalSurface(gp_Cylinder(ax, R));
        const Handle(Geom_RectangularTrimmedSurface) trim =
            new Geom_RectangularTrimmedSurface(basis, uLo, uHi, vLo, vHi);
        auto at = [&](double u, double v) { return gp_Pnt(R * std::cos(u), R * std::sin(u), v); };
        BRep_Builder bb;
        const double vt = Precision::Confusion();
        TopoDS_Vertex a, b, c, d;
        bb.MakeVertex(a, at(uLo, vLo), vt);
        bb.MakeVertex(b, at(uHi, vHi), vt);
        bb.MakeVertex(c, at(uHi, vLo), vt);
        bb.MakeVertex(d, at(uLo, vHi), vt);
        auto lineEdge = [&](const TopoDS_Vertex& p, const TopoDS_Vertex& q) {
            const gp_Pnt P = BRep_Tool::Pnt(p);
            const gp_Pnt Q = BRep_Tool::Pnt(q);
            return BRepBuilderAPI_MakeEdge(new Geom_Line(P, gp_Dir(gp_Vec(P, Q))), p, q);
        };
        BRepBuilderAPI_MakeEdge e0 = lineEdge(a, b);
        BRepBuilderAPI_MakeEdge e1 = lineEdge(b, c);
        BRepBuilderAPI_MakeEdge e2 = lineEdge(c, d);
        BRepBuilderAPI_MakeEdge e3 = lineEdge(d, a);
        expect(e0.IsDone() && e1.IsDone() && e2.IsDone() && e3.IsDone(),
               "unrepaired partial: bowtie edges built");
        TopoDS_Edge edges[4] = {e0.Edge(), e1.Edge(), e2.Edge(), e3.Edge()};
        const gp_Pnt2d uvCorner[4][2] = {{gp_Pnt2d(uLo, vLo), gp_Pnt2d(uHi, vHi)},
                                         {gp_Pnt2d(uHi, vHi), gp_Pnt2d(uHi, vLo)},
                                         {gp_Pnt2d(uHi, vLo), gp_Pnt2d(uLo, vHi)},
                                         {gp_Pnt2d(uLo, vHi), gp_Pnt2d(uLo, vLo)}};
        const TopLoc_Location loc;
        for (int i = 0; i < 4; ++i) {
            Standard_Real f = 0, l = 0;
            Handle(Geom_Curve) c3 = BRep_Tool::Curve(edges[i], f, l);
            const gp_Vec2d step(uvCorner[i][1].X() - uvCorner[i][0].X(),
                                uvCorner[i][1].Y() - uvCorner[i][0].Y());
            const double df = l - f;
            Handle(Geom2d_Line) own = new Geom2d_Line(
                gp_Pnt2d(uvCorner[i][0].X() - f * step.X() / df,
                         uvCorner[i][0].Y() - f * step.Y() / df),
                gp_Dir2d(step));
            bb.UpdateEdge(edges[i], own, trim, loc, vt);
            bb.Range(edges[i], trim, loc, f, l);
        }
        BRepBuilderAPI_MakeWire mkW;
        for (const TopoDS_Edge& e : edges) mkW.Add(e);
        expect(mkW.IsDone(), "unrepaired partial: wire built");
        if (mkW.IsDone()) {
            TopoDS_Face face;
            bb.MakeFace(face, trim, vt);
            bb.Add(face, mkW.Wire());
            expectFaceValidBeforeHeal(face, Standard_False,
                                      "unrepaired partial: face invalid before the heal");
            const Geom2d_Curve* raw[4] = {};
            bool lines = true;
            for (int i = 0; i < 4; ++i) {
                Standard_Real a0 = 0, b0 = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(edges[i], face, a0, b0);
                raw[i] = pc.get();
                if (Handle(Geom2d_Line)::DownCast(pc).IsNull()) lines = false;
            }
            expect(lines, "unrepaired partial: today's pcurves are lines");
            const bool healed = healInvalidPartialTrimChart(face, kTau);
            expect(!healed, "unrepaired partial: the chart does not validate the bowtie");
            bool kept = true;
            for (int i = 0; i < 4; ++i) {
                Standard_Real a0 = 0, b0 = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(edges[i], face, a0, b0);
                if (pc.get() != raw[i]) kept = false;
                if (Handle(Geom2d_Line)::DownCast(pc).IsNull()) kept = false;
                if (std::fabs(BRep_Tool::Tolerance(edges[i]) - vt) > Precision::Confusion())
                    kept = false;
            }
            BRepCheck_Analyzer after(face, Standard_True);
            expect(kept && after.IsValid() == Standard_False,
                   "unrepaired partial: today's pcurves kept, face still invalid");
        }
    }

    // D-fc-7.10 fallback. A shared corner whose pcurves meet within tau of
    // the vertex is eligible for the measured tolerance. The tolerance is
    // written only when BRepCheck then passes; this face stays invalid, so
    // the vertex tolerance is restored and today's path stands. The same
    // corner beyond tau is refused. A face that is already valid is not
    // touched (covered above).
    {
        const double R = 5.0;
        const double vHi = 2.0;
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_CylindricalSurface) basis =
            new Geom_CylindricalSurface(gp_Cylinder(ax, R));
        auto at = [&](double u, double v) { return gp_Pnt(R * std::cos(u), R * std::sin(u), v); };
        auto linePc = [&](double u0, double v0, double u1, double v1, double f, double l) {
            const double df = l - f;
            const gp_Vec2d d((u1 - u0) / df, (v1 - v0) / df);
            return Handle(Geom2d_Line)(
                new Geom2d_Line(gp_Pnt2d(u0 - f * d.X(), v0 - f * d.Y()), gp_Dir2d(d)));
        };
        // Build the face twice: once with the crossing inside tau, once past it.
        // The arc runs nearly a full turn so the crossing parameter sits far
        // from the shared vertex's parameter (the checker's yawn then fails)
        // while the 3D point stays within the requested gap of that vertex.
        auto build = [&](double gap, TopoDS_Face& face, TopoDS_Vertex& corner, TopoDS_Edge& shifted,
                         double& measured) -> bool {
            const double eps = gap / R;
            const double uEnd = 2.0 * kPi - 0.5 * eps;
            const double uCross = 2.0 * kPi - eps;
            const Handle(Geom_RectangularTrimmedSurface) trim =
                new Geom_RectangularTrimmedSurface(basis, 0.0, uEnd, 0.0, vHi);
            BRep_Builder bb;
            const double vt = Precision::Confusion();
            TopoDS_Vertex A, B, C, D;
            bb.MakeVertex(A, at(0.0, 0.0), vt);
            bb.MakeVertex(B, at(uEnd, 0.0), vt);
            bb.MakeVertex(C, at(uEnd, vHi), vt);
            bb.MakeVertex(D, at(0.0, vHi), vt);
            const gp_Ax2 cax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
            const Handle(Geom_TrimmedCurve) arc =
                new Geom_TrimmedCurve(new Geom_Circle(cax, R), 0.0, uEnd);
            BRepBuilderAPI_MakeEdge eArc(arc, A, B);
            BRepBuilderAPI_MakeEdge eGen(new Geom_Line(at(uEnd, 0.0), gp_Dir(0, 0, 1)), B, C);
            const gp_Pnt top0 = at(uEnd, vHi), top1 = at(0.0, vHi);
            BRepBuilderAPI_MakeEdge eTop(new Geom_Line(top0, gp_Dir(gp_Vec(top0, top1))), C, D);
            BRepBuilderAPI_MakeEdge eShift(new Geom_Line(at(0.0, vHi), gp_Dir(0, 0, -1)), D, A);
            if (!eArc.IsDone() || !eGen.IsDone() || !eTop.IsDone() || !eShift.IsDone()) return false;
            const TopLoc_Location loc;
            auto place = [&](const TopoDS_Edge& e, const Handle(Geom2d_Curve)& pc, double tol) {
                Standard_Real f = 0, l = 0;
                Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
                (void)c3;
                bb.UpdateEdge(e, pc, trim, loc, tol);
                bb.Range(e, trim, loc, f, l);
            };
            Standard_Real fa = 0, la = 0, fg = 0, lg = 0, ft = 0, lt = 0, fs = 0, ls = 0;
            BRep_Tool::Curve(eArc.Edge(), fa, la);
            BRep_Tool::Curve(eGen.Edge(), fg, lg);
            BRep_Tool::Curve(eTop.Edge(), ft, lt);
            BRep_Tool::Curve(eShift.Edge(), fs, ls);
            place(eArc.Edge(), linePc(0.0, 0.0, uEnd, 0.0, fa, la), vt);
            place(eGen.Edge(), linePc(uEnd, 0.0, uEnd, vHi, fg, lg), vt);
            place(eTop.Edge(), linePc(uEnd, vHi, 2.0 * kPi, vHi, ft, lt), vt);
            // The shifted generator's pcurve meets the arc inside the arc's
            // parameter, gap millimetres from the shared vertex. Its edge
            // tolerance covers that curve-on-surface gap; the vertex does not.
            place(eShift.Edge(), linePc(uCross, vHi, uCross, 0.0, fs, ls), gap);
            bb.MakeFace(face, trim, vt);
            BRepBuilderAPI_MakeWire mkW;
            mkW.Add(eArc.Edge());
            mkW.Add(eGen.Edge());
            mkW.Add(eTop.Edge());
            mkW.Add(eShift.Edge());
            if (!mkW.IsDone()) return false;
            bb.Add(face, mkW.Wire());
            corner = A;
            shifted = eShift.Edge();
            measured = at(0.0, 0.0).Distance(at(uCross, 0.0));
            return true;
        };
        const double gapIn = 0.5 * kTau;
        TopoDS_Face faceIn;
        TopoDS_Vertex cornerIn;
        TopoDS_Edge shiftIn;
        double measuredIn = 0;
        expect(build(gapIn, faceIn, cornerIn, shiftIn, measuredIn) && measuredIn <= kTau &&
                   measuredIn > Precision::Confusion(),
               "corner tol: a within-tau gap is built");
        if (!faceIn.IsNull()) {
            const double vt0 = BRep_Tool::Tolerance(cornerIn);
            const double et0 = BRep_Tool::Tolerance(shiftIn);
            BRepCheck_Analyzer before(faceIn, Standard_True);
            expect(before.IsValid() == Standard_False, "corner tol: within-tau face is invalid");
            const bool healed = raiseMeasuredCornerTol(faceIn, kTau);
            BRepCheck_Analyzer after(faceIn, Standard_True);
            const double vt1 = BRep_Tool::Tolerance(cornerIn);
            expect(!healed && after.IsValid() == Standard_False,
                   "corner tol: within tau, still invalid, today's path kept");
            expect(std::fabs(vt1 - vt0) <= Precision::Confusion(),
                   "corner tol: an unvalidated raise restores the vertex tolerance");
            expect(std::fabs(BRep_Tool::Tolerance(shiftIn) - et0) <= Precision::Confusion(),
                   "corner tol: the edge tolerance is not rewritten");
        }
        const double gapOut = 2.0 * kTau;
        TopoDS_Face faceOut;
        TopoDS_Vertex cornerOut;
        TopoDS_Edge shiftOut;
        double measuredOut = 0;
        expect(build(gapOut, faceOut, cornerOut, shiftOut, measuredOut) && measuredOut > kTau,
               "corner tol: a beyond-tau gap is built");
        if (!faceOut.IsNull()) {
            const double vt0 = BRep_Tool::Tolerance(cornerOut);
            const bool healed = raiseMeasuredCornerTol(faceOut, kTau);
            BRepCheck_Analyzer after(faceOut, Standard_True);
            expect(!healed && after.IsValid() == Standard_False &&
                       std::fabs(BRep_Tool::Tolerance(cornerOut) - vt0) <= Precision::Confusion(),
                   "corner tol: beyond tau today's vertex tolerance stands");
        }
    }

    // D-fc-7.16. The rotated frame is a measured function of the trim, and it
    // applies only to a partial cylinder that touches the seam or already
    // carries a seam edge.
    {
        const double uHi = 1.2;
        const double twoPi = 2.0 * kPi;
        expect(partialTrimTouchesPeriodicSeam(0.0, uHi), "seam frame: u = 0 touches the seam");
        expect(partialTrimTouchesPeriodicSeam(0.2, twoPi), "seam frame: u = 2pi touches the seam");
        expect(!partialTrimTouchesPeriodicSeam(0.4, 1.1), "seam frame: an interior trim does not");
        expect(std::fabs(partialCylSeamFrameAngle(0.0, uHi) - (0.5 * uHi - kPi)) <=
                   Precision::PConfusion(),
               "seam frame: the rotation centres the trim at pi");
        expect(partialCylSeamFrameApplies(0.0, uHi, 0), "seam frame: a seam-touching partial applies");
        expect(partialCylSeamFrameApplies(0.4, 1.1, 3),
               "seam frame: a partial with a seam edge applies");
        expect(!partialCylSeamFrameApplies(0.4, 1.1, 0),
               "seam frame: an interior partial with no seam edge does not");
        expect(!partialCylSeamFrameApplies(0.0, twoPi, 3),
               "seam frame: a full turn is not this heal");
    }

    // D-fc-7.16. A partial cylinder whose trim starts at u = 0 and whose
    // closing edges are parameterised on the periodic image at 2π is invalid.
    // The same edges, charted on the frame rotated by (uLo + uHi)/2 − π, are
    // one valid face. A trim that starts at u = 0 and is already valid is not
    // rewritten.
    {
        const double R = 3.0;
        const double uLo = 0.0;
        const double uHi = 1.2;
        const double vLo = 0.0;
        const double vHi = 3.0;
        const double twoPi = 2.0 * kPi;
        const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_CylindricalSurface) basis =
            new Geom_CylindricalSurface(gp_Cylinder(ax, R));
        auto at = [&](double u, double v) { return gp_Pnt(R * std::cos(u), R * std::sin(u), v); };
        auto place = [&](BRep_Builder& bb, const TopoDS_Edge& e,
                         const Handle(Geom_Surface)& s, double ua, double va, double ub, double vb,
                         double tol) {
            Standard_Real f = 0, l = 0;
            BRep_Tool::Curve(e, f, l);
            const double df = l - f;
            const gp_Vec2d d((ub - ua) / df, (vb - va) / df);
            Handle(Geom2d_Curve) pc;
            if (std::abs(d.Magnitude() - 1.0) <= Precision::PConfusion())
                pc = new Geom2d_Line(gp_Pnt2d(ua - f * d.X(), va - f * d.Y()), gp_Dir2d(d));
            else {
                TColgp_Array1OfPnt2d poles(1, 2);
                poles.SetValue(1, gp_Pnt2d(ua, va));
                poles.SetValue(2, gp_Pnt2d(ub, vb));
                TColStd_Array1OfReal knots(1, 2);
                knots.SetValue(1, f);
                knots.SetValue(2, l);
                TColStd_Array1OfInteger mults(1, 2);
                mults.SetValue(1, 2);
                mults.SetValue(2, 2);
                pc = new Geom2d_BSplineCurve(poles, knots, mults, 1);
            }
            const TopLoc_Location loc;
            bb.UpdateEdge(e, pc, s, loc, tol);
            bb.Range(e, s, loc, f, l);
        };
        BRep_Builder bb;
        const double vt = Precision::Confusion();
        TopoDS_Vertex v00, v10, v11, v01;
        bb.MakeVertex(v00, at(uLo, vLo), vt);
        bb.MakeVertex(v10, at(uHi, vLo), vt);
        bb.MakeVertex(v11, at(uHi, vHi), vt);
        bb.MakeVertex(v01, at(uLo, vHi), vt);
        const gp_Ax2 cLo(gp_Pnt(0, 0, vLo), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const gp_Ax2 cHi(gp_Pnt(0, 0, vHi), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        const Handle(Geom_TrimmedCurve) arcLo =
            new Geom_TrimmedCurve(new Geom_Circle(cLo, R), uLo, uHi);
        // Sense false walks the short arc from uHi back to uLo.
        const Handle(Geom_TrimmedCurve) arcHi =
            new Geom_TrimmedCurve(new Geom_Circle(cHi, R), uLo, uHi, Standard_False);
        BRepBuilderAPI_MakeEdge eArc(arcLo, v00, v10);
        BRepBuilderAPI_MakeEdge eGenHi(new Geom_Line(at(uHi, vLo), gp_Dir(0, 0, 1)), v10, v11);
        BRepBuilderAPI_MakeEdge eArcHi(arcHi, v11, v01);
        BRepBuilderAPI_MakeEdge eGenLo(new Geom_Line(at(uLo, vHi), gp_Dir(0, 0, -1)), v01, v00);
        expect(eArc.IsDone() && eGenHi.IsDone() && eArcHi.IsDone() && eGenLo.IsDone(),
               "seam frame: boundary edges built");
        if (eArc.IsDone() && eGenHi.IsDone() && eArcHi.IsDone() && eGenLo.IsDone()) {
            const Handle(Geom_RectangularTrimmedSurface) trim =
                new Geom_RectangularTrimmedSurface(basis, uLo, uHi, vLo, vHi);
            // The closing pair sits on the periodic image at 2π, not on u = 0.
            place(bb, eArc.Edge(), trim, uLo, vLo, uHi, vLo, vt);
            place(bb, eGenHi.Edge(), trim, uHi, vLo, uHi, vHi, vt);
            place(bb, eArcHi.Edge(), trim, uHi, vHi, twoPi, vHi, vt);
            place(bb, eGenLo.Edge(), trim, twoPi, vHi, twoPi, vLo, vt);
            BRepBuilderAPI_MakeWire mkW;
            mkW.Add(eArc.Edge());
            mkW.Add(eGenHi.Edge());
            mkW.Add(eArcHi.Edge());
            mkW.Add(eGenLo.Edge());
            expect(mkW.IsDone(), "seam frame: wire closed");
            if (mkW.IsDone()) {
                TopoDS_Face face;
                bb.MakeFace(face, trim, vt);
                bb.Add(face, mkW.Wire());
                BRepCheck_Analyzer before(face, Standard_True);
                expect(before.IsValid() == Standard_False,
                       "seam frame: trim at u = 0 with the periodic image is invalid");
                const bool healed = healPartialCylSeamFrame(face, kTau);
                BRepCheck_Analyzer after(face, Standard_True);
                expect(healed && after.IsValid() == Standard_True,
                       "seam frame: the rotated frame is a valid face");
                Handle(Geom_RectangularTrimmedSurface) got =
                    Handle(Geom_RectangularTrimmedSurface)::DownCast(BRep_Tool::Surface(face));
                expect(!got.IsNull(), "seam frame: the shipped surface is the rotated trim");
                if (!got.IsNull()) {
                    Standard_Real gu0 = 0, gu1 = 0, gv0 = 0, gv1 = 0;
                    got->Bounds(gu0, gu1, gv0, gv1);
                    expect(std::fabs(0.5 * (gu0 + gu1) - kPi) <= Precision::PConfusion(),
                           "seam frame: shipped trim is centred at pi");
                    expect(gu0 > Precision::PConfusion() &&
                               gu1 < twoPi - Precision::PConfusion(),
                           "seam frame: shipped trim is clear of the seam");
                    Handle(Geom_CylindricalSurface) gc =
                        Handle(Geom_CylindricalSurface)::DownCast(got->BasisSurface());
                    expect(!gc.IsNull() && std::fabs(gc->Radius() - R) <= Precision::Confusion() &&
                               gc->Position().Location().Distance(ax.Location()) <=
                                   Precision::Confusion() &&
                               gc->Position().Direction().IsEqual(ax.Direction(),
                                                                  Precision::Angular()),
                           "seam frame: axis, location and radius are unchanged");
                    bool onePc = !gc.IsNull();
                    for (TopExp_Explorer ex(face, TopAbs_EDGE); ex.More() && onePc; ex.Next()) {
                        Standard_Real a = 0, b = 0;
                        Handle(Geom2d_Curve) pc =
                            BRep_Tool::CurveOnSurface(TopoDS::Edge(ex.Current()), face, a, b);
                        if (pc.IsNull()) {
                            onePc = false;
                            continue;
                        }
                        const double uA = pc->Value(a).X();
                        const double uB = pc->Value(b).X();
                        if (uA < gu0 - Precision::PConfusion() || uA > gu1 + Precision::PConfusion() ||
                            uB < gu0 - Precision::PConfusion() || uB > gu1 + Precision::PConfusion())
                            onePc = false;
                    }
                    expect(onePc, "seam frame: one pcurve per edge, inside the rotated trim");
                    TopoDS_Shell sh;
                    bb.MakeShell(sh);
                    bb.Add(sh, face);
                    BRepCheck_Analyzer anSh(sh, Standard_True);
                    expect(anSh.IsValid() == Standard_True,
                           "seam frame: the shell that contains the face is valid");
                }
            }
            // The same boundary with pcurves on the trim itself is valid, and
            // the heal must not replace that surface or those pcurves.
            TopoDS_Vertex w00, w10, w11, w01;
            bb.MakeVertex(w00, at(uLo, vLo), vt);
            bb.MakeVertex(w10, at(uHi, vLo), vt);
            bb.MakeVertex(w11, at(uHi, vHi), vt);
            bb.MakeVertex(w01, at(uLo, vHi), vt);
            BRepBuilderAPI_MakeEdge fArc(arcLo, w00, w10);
            BRepBuilderAPI_MakeEdge fHi(new Geom_Line(at(uHi, vLo), gp_Dir(0, 0, 1)), w10, w11);
            BRepBuilderAPI_MakeEdge fTop(arcHi, w11, w01);
            BRepBuilderAPI_MakeEdge fLo(new Geom_Line(at(uLo, vHi), gp_Dir(0, 0, -1)), w01, w00);
            expect(fArc.IsDone() && fHi.IsDone() && fTop.IsDone() && fLo.IsDone(),
                   "seam frame: valid-at-first edges built");
            if (fArc.IsDone() && fHi.IsDone() && fTop.IsDone() && fLo.IsDone()) {
                const Handle(Geom_RectangularTrimmedSurface) okTrim =
                    new Geom_RectangularTrimmedSurface(basis, uLo, uHi, vLo, vHi);
                place(bb, fArc.Edge(), okTrim, uLo, vLo, uHi, vLo, vt);
                place(bb, fHi.Edge(), okTrim, uHi, vLo, uHi, vHi, vt);
                place(bb, fTop.Edge(), okTrim, uHi, vHi, uLo, vHi, vt);
                place(bb, fLo.Edge(), okTrim, uLo, vHi, uLo, vLo, vt);
                BRepBuilderAPI_MakeWire okW;
                okW.Add(fArc.Edge());
                okW.Add(fHi.Edge());
                okW.Add(fTop.Edge());
                okW.Add(fLo.Edge());
                expect(okW.IsDone(), "seam frame: valid-at-first wire closed");
                if (okW.IsDone()) {
                    TopoDS_Face okFace;
                    bb.MakeFace(okFace, okTrim, vt);
                    bb.Add(okFace, okW.Wire());
                    BRepCheck_Analyzer okBefore(okFace, Standard_True);
                    expect(okBefore.IsValid() == Standard_True,
                           "seam frame: a trim at u = 0 with its own pcurves is valid");
                    const Geom_Surface* surf0 = BRep_Tool::Surface(okFace).get();
                    const Geom2d_Curve* pcs[4] = {};
                    TopoDS_Edge eds[4] = {fArc.Edge(), fHi.Edge(), fTop.Edge(), fLo.Edge()};
                    for (int i = 0; i < 4; ++i) {
                        Standard_Real a = 0, b = 0;
                        pcs[i] = BRep_Tool::CurveOnSurface(eds[i], okFace, a, b).get();
                    }
                    const bool healed = healPartialCylSeamFrame(okFace, kTau);
                    const Geom_Surface* surf1 = BRep_Tool::Surface(okFace).get();
                    bool same = surf0 == surf1;
                    for (int i = 0; i < 4 && same; ++i) {
                        Standard_Real a = 0, b = 0;
                        if (BRep_Tool::CurveOnSurface(eds[i], okFace, a, b).get() != pcs[i])
                            same = false;
                    }
                    BRepCheck_Analyzer okAfter(okFace, Standard_True);
                    expect(!healed && same && okAfter.IsValid() == Standard_True,
                           "seam frame: a face valid at first build is untouched");
                }
            }
        }
    }

    std::fprintf(stderr, "relational_axis_seam_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
