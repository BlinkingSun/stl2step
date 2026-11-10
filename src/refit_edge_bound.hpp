// Exact edge-on-face deviation for edges whose pcurve FixAddPCurve projected
// (fx-hubtori, fx-tori-r1 Open 2). No threshold of its own: the measure is
// OCCT's exact curve-on-surface maximum, the bound is the caller's mesh cap.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_EDGE_BOUND_HPP
#define STL2STEP_REFIT_EDGE_BOUND_HPP

#include <algorithm>
#include <vector>

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
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

// D-gid-8.28. The tolerance an edge carries into the shell's forced
// SameParameter is the one the builder wrote for it (the binder's exact
// curve-on-surface bound, a measured seam gap, ... but also a placeholder such
// as the sewing tolerance a pcurve is bound with). SameParameter re-measures the
// edge on its own samples and can write less than the edge's exact deviation,
// and a valid face can then fail its wire check (D-fc-6.9: a tolerance below
// its measured value is forbidden). One row per edge TShape: `before` is the
// tolerance it carried in, `after` the one SameParameter left, `kept` the one
// the edge leaves with (filled by keepSnapshotEdgeTolerances).
struct EdgeTolRow {
    TopoDS_Edge e;
    double before = 0.0;
    double after = 0.0;
    double kept = 0.0;
    double measured = -1.0;
};

inline std::vector<EdgeTolRow> snapshotEdgeTolerances(const TopoDS_Shape& s) {
    std::vector<EdgeTolRow> rows;
    if (s.IsNull()) return rows;
    TopTools_IndexedMapOfShape emap;
    TopExp::MapShapes(s, TopAbs_EDGE, emap);
    rows.reserve((size_t)emap.Extent());
    for (int i = 1; i <= emap.Extent(); i++) {
        const TopoDS_Edge& e = TopoDS::Edge(emap(i));
        const double t = BRep_Tool::Tolerance(e);
        rows.push_back({e, t, t, t, -1.0});
    }
    return rows;
}

// After the shell's SameParameter, every edge SameParameter lowered keeps
//   max(its tolerance now, min(the tolerance it carried in, measure(e))),
// where measure(e) is the edge's measured deviation (< 0: not measurable, the
// edge is left as SameParameter wrote it). A measure that always returns +inf
// keeps the carried-in tolerance itself. Tolerances only rise, never past the
// value the edge already had; no curve, pcurve, range or vertex is touched.
// Returns every edge SameParameter lowered, raised or not (kept == after when
// it was not raised).
template <class Measure>
inline std::vector<EdgeTolRow> keepSnapshotEdgeTolerances(const std::vector<EdgeTolRow>& snap,
                                                          Measure&& measure) {
    std::vector<EdgeTolRow> lowered;
    BRep_Builder B;
    for (const EdgeTolRow& r : snap) {
        if (r.e.IsNull()) continue;
        const double now = BRep_Tool::Tolerance(r.e);
        if (!(r.before > now)) continue;
        const double m = measure(r.e);
        const double keep = (m >= 0.0) ? std::min(r.before, m) : now;
        if (keep > now) B.UpdateEdge(r.e, keep);
        lowered.push_back({r.e, r.before, now, BRep_Tool::Tolerance(r.e), m});
    }
    return lowered;
}

// The measure the builder bounds an edge with (boundEdgesExactOnFace): the
// exact curve-on-surface maximum on every face of `faces` that the edge bounds,
// plus Precision::Confusion(). -1 when no face measures it.
inline double exactEdgeDeviationOnFaces(const TopoDS_Edge& e, const TopTools_ListOfShape& faces) {
    double m = -1.0;
    for (TopTools_ListOfShape::Iterator it(faces); it.More(); it.Next()) {
        if (it.Value().ShapeType() != TopAbs_FACE) continue;
        const double d = exactEdgeDeviationOnFace(e, TopoDS::Face(it.Value()));
        if (d > m) m = d;
    }
    return m;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_EDGE_BOUND_HPP
