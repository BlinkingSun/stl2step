// Exact edge-on-face deviation for edges whose pcurve FixAddPCurve projected
// (fx-hubtori, fx-tori-r1 Open 2). No threshold of its own: the measure is
// OCCT's exact curve-on-surface maximum, the bound is the caller's mesh cap.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_EDGE_BOUND_HPP
#define STL2STEP_REFIT_EDGE_BOUND_HPP

#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>

#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRep_Builder.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Tool.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>

namespace stl2step {
namespace refit {

// The exact maximum of |C(t) - S(pc(t))| over the edge's range on face f, plus
// Precision::Confusion() so a sampled check (BRepCheck's control points) can
// never read more than the edge carries. -1 when it cannot be measured.
inline double exactEdgeDeviationOnFace(const TopoDS_Edge& e, const TopoDS_Face& f) {
    if (e.IsNull() || f.IsNull() || BRep_Tool::Degenerated(e)) return -1.0;
    try {
        BRepLib_CheckCurveOnSurface chk(e, f);
        chk.Perform();
        if (!chk.IsDone()) return -1.0;
        return chk.MaxDistance() + Precision::Confusion();
    } catch (const Standard_Failure&) {
        return -1.0;
    }
}

enum class EdgeBound { Keep, Raise, OverCap };

// A tolerance is only ever raised, and never past cap.
inline EdgeBound edgeBoundOf(double exactDev, double tol, double cap) {
    if (!(exactDev > tol)) return EdgeBound::Keep;
    if (exactDev > cap) return EdgeBound::OverCap;
    return EdgeBound::Raise;
}

// BRepCheck's single-wire rule: on the FORWARD face the one wire must leave the
// infinite point of the (u, v) plane OUT, i.e. run counter-clockwise about the
// surface normal. BRepLib_MakeFace tests that itself (CheckInside) only for
// edges that already carry pcurves; a mesh-rim wire gets its pcurves later
// (addPcurvesOnFace). True when f has exactly one wire and that wire is
// inside out; its reversal is the same closed loop.
inline bool singleWireInsideOut(const TopoDS_Face& f) {
    if (f.IsNull()) return false;
    const TopoDS_Face fF = TopoDS::Face(f.Oriented(TopAbs_FORWARD));
    int nW = 0;
    for (TopoDS_Iterator it(fF); it.More(); it.Next())
        if (it.Value().ShapeType() == TopAbs_WIRE) ++nW;
    if (nW != 1) return false;
    try {
        BRepTopAdaptor_FClass2d cl(fF, Precision::PConfusion());
        return cl.PerformInfinitePoint() != TopAbs_OUT;
    } catch (const Standard_Failure&) {
        return false;
    }
}

// A single-wire face whose wire is inside out is rebuilt with that wire reversed, on
// the same surface and with the same orientation. Every edge keeps its curves and
// pcurves; nothing moves. False (f untouched) when the wire is not inside out.
inline bool reverseInsideOutWire(TopoDS_Face& f) {
    if (!singleWireInsideOut(f)) return false;
    try {
        const TopoDS_Face fF = TopoDS::Face(f.Oriented(TopAbs_FORWARD));
        TopoDS_Wire w;
        for (TopoDS_Iterator it(fF); it.More(); it.Next())
            if (it.Value().ShapeType() == TopAbs_WIRE) w = TopoDS::Wire(it.Value());
        TopoDS_Face nf = TopoDS::Face(fF.EmptyCopied());
        BRep_Builder B;
        B.Add(nf, w.Reversed());
        nf.Orientation(f.Orientation());
        f = nf;
        return true;
    } catch (const Standard_Failure&) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_EDGE_BOUND_HPP
