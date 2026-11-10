// D-gid-8.7. The periodic branch of a projected pcurve on a partial cylinder.
//
// GeomProjLib frames a projection onto a periodic surface into
// [UFirst, UFirst + 2π). An edge a hair below the face chart's lower end -- a
// generator that IS the chart boundary but lies off an uncertified cylinder,
// so the closed-form generator pcurve is refused -- comes back a whole turn
// away. The wire then jumps 2π, BRepCheck still passes on the periodic
// surface, and the face closes on its complement.
//
// The projection is kept as it is; only its period is chosen, by a whole turn
// (exact), as the image nearest the face chart (chartBranchShift). The chart is
// the surface's own when it is trimmed in u, otherwise the region's [uMin, uMax]
// carried into the surface's frame (cylChartInFrame). No tolerance is involved.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_PCURVE_BRANCH_HPP
#define STL2STEP_REFIT_PCURVE_BRANCH_HPP

#include <Geom2d_Curve.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Curve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <gp_Ax3.hxx>
#include <gp_Vec2d.hxx>

#include <cmath>

#include "refit_partial_cyl.hpp"

namespace stl2step {
namespace refit {

// The region chart [uMin, uMax] (uMax < uMin means it wraps once) in the frame
// `pos` of one of the region's cylinder surfaces. `pos` must be a turn of the
// region frame `ax` about the same axis, same sense and handedness; then
// u_pos = u_ax - rot with rot = the angle of pos's XDirection in ax. Returns
// false otherwise, or for a full turn (no branch to choose).
inline bool cylChartInFrame(const gp_Ax3& ax, double uMin, double uMax, const gp_Ax3& pos,
                            double chart[2]) {
    constexpr double kPiLocal = 3.14159265358979323846264338327950288;
    if (pos.Direct() != ax.Direct()) return false;
    if (pos.Direction().Dot(ax.Direction()) <= 0.0) return false;
    double u0 = uMin, u1 = uMax;
    if (u1 < u0) u1 += 2.0 * kPiLocal;
    if (!(u1 - u0 < 2.0 * kPiLocal)) return false;
    const gp_Dir xv = pos.XDirection();
    const double rot = std::atan2(xv.Dot(ax.YDirection()), xv.Dot(ax.XDirection()));
    chart[0] = u0 - rot;
    chart[1] = u1 - rot;
    return true;
}

// The face chart of a pcurve on `srf`: a sheet trimmed in u (span < 2π) is its
// own chart, its bounds; else `uChart` when the caller has one; else none.
inline bool faceChartU(const Handle(Geom_Surface)& srf, const double* uChart, double& c0,
                       double& c1, const char** from = nullptr) {
    constexpr double kPiLocal = 3.14159265358979323846264338327950288;
    Handle(Geom_RectangularTrimmedSurface) rt =
        Handle(Geom_RectangularTrimmedSurface)::DownCast(srf);
    if (!rt.IsNull()) {
        Standard_Real a = 0, b = 0, va = 0, vb = 0;
        rt->Bounds(a, b, va, vb);
        if (std::isfinite(a) && std::isfinite(b) && b >= a && b - a < 2.0 * kPiLocal) {
            c0 = a;
            c1 = b;
            if (from) *from = "trim";
            return true;
        }
    }
    if (uChart) {
        c0 = uChart[0];
        c1 = uChart[1];
        if (from) *from = "region";
        return true;
    }
    return false;
}

struct ChartBranchInfo {
    const char* from = "none";  // trim | region | none
    double c0 = 0.0, c1 = 0.0;  // the chart used
    double uMid = 0.0;          // projected u at the curve's mid parameter
    double shift = 0.0;         // whole-period shift applied (0 = unchanged)
};

// GeomProjLib::Curve2d, then the whole-period shift onto the face chart's
// branch. Without a chart the projection is returned unchanged.
inline Handle(Geom2d_Curve) projectOnChartBranch(const Handle(Geom_Curve)& c3, Standard_Real f,
                                                 Standard_Real l, const Handle(Geom_Surface)& srf,
                                                 const double* uChart,
                                                 ChartBranchInfo* info = nullptr) {
    Handle(Geom2d_Curve) pj = GeomProjLib::Curve2d(c3, f, l, srf);
    if (pj.IsNull()) return pj;
    ChartBranchInfo in;
    if (!faceChartU(srf, uChart, in.c0, in.c1, &in.from)) {
        if (info) *info = in;
        return pj;
    }
    in.uMid = pj->Value(0.5 * (f + l)).X();
    in.shift = chartBranchShift(in.uMid, in.c0, in.c1);
    if (info) *info = in;
    if (in.shift == 0.0) return pj;
    Handle(Geom2d_Curve) moved =
        Handle(Geom2d_Curve)::DownCast(pj->Translated(gp_Vec2d(in.shift, 0.0)));
    return moved.IsNull() ? pj : moved;
}

}  // namespace refit
}  // namespace stl2step

#endif
