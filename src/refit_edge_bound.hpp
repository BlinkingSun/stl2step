// Exact edge-on-face deviation for edges whose pcurve FixAddPCurve projected
// (fx-hubtori, fx-tori-r1 Open 2). No threshold of its own: the measure is
// OCCT's exact curve-on-surface maximum, the bound is the caller's mesh cap.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_EDGE_BOUND_HPP
#define STL2STEP_REFIT_EDGE_BOUND_HPP

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <BRepLib_CheckCurveOnSurface.hxx>
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

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_EDGE_BOUND_HPP
