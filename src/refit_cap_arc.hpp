// Builder predicates for the cap arcs and corner faces of relational blends.
// Pure OCCT geometry and BRepCheck statuses, no thresholds.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_CAP_ARC_HPP
#define STL2STEP_REFIT_CAP_ARC_HPP

#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <ElCLib.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>

#include <cmath>

namespace stl2step {
namespace refit {

constexpr double kCapArcPi = 3.14159265358979323846264338327950288;

// The cap arc of a two-terminal chain (one mesh edge) on a partial cylinder
// (D-fc-3.2.10). pa and pb are the chain's terminals on the cap circle, hint is
// the region's fitted u-span midpoint on that circle, [u0, u1] (u1 >= u0,
// unwrapped) is the region's fitted u-span and azimuth(p) is the region's own
// azimuth of a point on any 2π branch. A chord's two arcs are its short arc and
// the complement. True, with shortMid the short arc's middle, when the hint is
// off the short arc and the short arc lies inside the fitted span: one edge of
// the region's triangles then spans the short arc, never its complement. False
// (the hint stands) otherwise, and for coincident terminals.
template <class Azimuth>
inline bool shortArcInsideSpanOffMid(const gp_Circ& circ, const gp_Pnt& pa, const gp_Pnt& pb,
                                     const gp_Pnt& hint, double u0, double u1, Azimuth&& azimuth,
                                     gp_Pnt& shortMid) {
    const double kTwoPi = 2.0 * kCapArcPi;
    auto wrap0 = [&](double t) {
        t = std::fmod(t, kTwoPi);
        if (t < 0.0) t += kTwoPi;
        return t;
    };
    const double p1 = ElCLib::Parameter(circ, pa);
    const double p2 = ElCLib::Parameter(circ, pb);
    const double pm = ElCLib::Parameter(circ, hint);
    double df = wrap0(p2 - p1);
    if (df > kCapArcPi) df -= kTwoPi;
    if (df == 0.0) return false;
    const double dm = (df > 0.0) ? wrap0(pm - p1) : wrap0(p1 - pm);
    if (dm <= std::fabs(df)) return false;
    const gp_Pnt m = ElCLib::Value(p1 + 0.5 * df, circ);
    const double t = u0 + wrap0(azimuth(m) - u0);
    if (t > u1) return false;
    shortMid = m;
    return true;
}

// The face's BRepCheck statuses, split: orientOnFace counts the face's own
// orientation statuses (UnorientableShape, BadOrientation,
// BadOrientationOfSubshape), otherBad every other non-NoError status on the
// face or on any wire, edge or vertex of it, in the face's context.
inline void faceCheckSplit(const TopoDS_Face& f, int& orientOnFace, int& otherBad) {
    orientOnFace = 0;
    otherBad = 0;
    if (f.IsNull()) {
        ++otherBad;
        return;
    }
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        auto tally = [&](const Handle(BRepCheck_Result)& res, bool faceLevel) {
            if (res.IsNull()) return;
            for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next()) {
                const BRepCheck_Status st = it.Value();
                if (st == BRepCheck_NoError) continue;
                const bool orient = st == BRepCheck_UnorientableShape ||
                                    st == BRepCheck_BadOrientation ||
                                    st == BRepCheck_BadOrientationOfSubshape;
                if (faceLevel && orient) ++orientOnFace;
                else ++otherBad;
            }
            for (res->InitContextIterator(); res->MoreShapeInContext(); res->NextShapeInContext()) {
                if (!res->ContextualShape().IsSame(f)) continue;
                for (BRepCheck_ListOfStatus::Iterator it(res->StatusOnShape()); it.More();
                     it.Next())
                    if (it.Value() != BRepCheck_NoError) ++otherBad;
            }
        };
        tally(an.Result(f), true);
        for (TopAbs_ShapeEnum k : {TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX})
            for (TopExp_Explorer ex(f, k); ex.More(); ex.Next()) tally(an.Result(ex.Current()), false);
    } catch (const Standard_Failure&) {
        ++otherBad;
    }
}

// A partial torus face (D-fc-4.6 build) whose only BRepCheck failures are the
// face's own orientation statuses has a wire the walk closed in the opposite
// sense to the surface's (u, v); every edge and pcurve checks clean. The face
// is retried once on the reversed wire. A face with any other failure is not.
inline bool partialTorusWireRetry(int orientStatusesOnFace, int otherStatuses) {
    return orientStatusesOnFace > 0 && otherStatuses == 0;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_CAP_ARC_HPP
