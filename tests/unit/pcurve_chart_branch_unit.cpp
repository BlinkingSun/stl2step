// D-gid-8.7: the fallback projection of an edge that lies off a partial
// cylinder keeps the face chart's periodic branch. Synthetic geometry only: a
// cylinder of radius 50 about +Z and generators parallel to its axis, set a
// small radial offset off the surface (so the closed-form generator pcurve
// would be refused and the projector is used), at the chart's lower end -- the
// cylinder's u = 0 seam lies just inside the chart.
// SPDX-License-Identifier: MIT

#include "refit_partial_cyl.hpp"
#include "refit_pcurve_branch.hpp"

#include <ElSLib.hxx>
#include <GeomProjLib.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <cstdio>
#include <cstring>

using namespace stl2step::refit;

namespace {

constexpr double kPiT = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPiT;
constexpr double kEq = 1e-9;  // test comparison only (exact arithmetic up to rounding)

int gFail = 0;
int gPass = 0;

void expect(bool ok, const char* name) {
    if (ok) {
        ++gPass;
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

// A line parallel to the axis of `ax` at azimuth u and distance rho from it.
Handle(Geom_Line) generatorAt(const gp_Ax3& ax, double u, double rho) {
    const gp_Vec rad =
        gp_Vec(ax.XDirection()) * std::cos(u) + gp_Vec(ax.YDirection()) * std::sin(u);
    const gp_Pnt p = ax.Location().Translated(rad * rho);
    return new Geom_Line(gp_Lin(p, ax.Direction()));
}

double uAt(const Handle(Geom2d_Curve)& pc, double t) { return pc->Value(t).X(); }
double vAt(const Handle(Geom2d_Curve)& pc, double t) { return pc->Value(t).Y(); }

// u of a 3D point in a cylinder frame, normalised to [0, 2π) by ElSLib.
double uInFrame(const gp_Ax3& pos, double R, const gp_Pnt& p) {
    Standard_Real u = 0, v = 0;
    ElSLib::CylinderParameters(pos, R, p, u, v);
    return u;
}

double modTwoPiDist(double a, double b) {
    double d = std::fmod(std::fabs(a - b), kTwoPi);
    return std::min(d, kTwoPi - d);
}

void testShiftPure() {
    const double u0 = -0.01, u1 = 0.6, mid = 0.5 * (u0 + u1);
    expect(chartBranchShift(mid, u0, u1) == 0.0, "shift_on_chart_stays");
    expect(chartBranchShift(u0, u0, u1) == 0.0, "shift_chart_lower_end_stays");
    expect(chartBranchShift(u1, u0, u1) == 0.0, "shift_chart_upper_end_stays");
    expect(chartBranchShift(u0 - 1e-6, u0, u1) == 0.0, "shift_hair_below_stays");
    // The projector's image of a curve a hair below u0: u0 + 2π - ε.
    expect(std::fabs(chartBranchShift(u0 + kTwoPi - 1e-6, u0, u1) + kTwoPi) < kEq,
           "shift_turn_above_pulled_back");
    expect(std::fabs(chartBranchShift(mid + 2.0 * kTwoPi, u0, u1) + 2.0 * kTwoPi) < kEq,
           "shift_two_turns_above_pulled_back");
    expect(std::fabs(chartBranchShift(mid - kTwoPi, u0, u1) - kTwoPi) < kEq,
           "shift_turn_below_pushed_up");
    // partialCylPeriodShift only looks at images in [u0, u0 + 2π): it leaves
    // the turn-above image where it is. That is the gap this lane closes.
    expect(partialCylPeriodShift(u0 + kTwoPi - 1e-6, u0, u1) == 0.0,
           "partial_shift_blind_to_turn_above");
    // The fold is the antipode of the chart centre: an image off the chart
    // moves only to an image that is nearer to it.
    expect(chartBranchShift(3.0, 0.0, 1.0) == 0.0, "shift_nearer_branch_kept");
    expect(std::fabs(chartBranchShift(4.5, 0.0, 1.0) + kTwoPi) < kEq, "shift_farther_branch_moved");
    // No branch to choose on a full turn, nor on a malformed chart.
    expect(chartBranchShift(mid + kTwoPi, 0.0, kTwoPi) == 0.0, "shift_full_turn_none");
    expect(chartBranchShift(mid + kTwoPi, u1, u0) == 0.0, "shift_reversed_chart_none");
    expect(chartBranchShift(std::nan(""), u0, u1) == 0.0, "shift_nan_none");
}

// The S29-class defect: an off-surface generator on the chart's lower end, on
// the untrimmed cylinder whose chart is [u0, u1] with u0 < 0 < u1.
void testUntrimmedSeam() {
    const double R = 50.0, off = 0.04, u0 = -0.01, u1 = 0.6, f = 0.0, l = 5.0;
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    Handle(Geom_CylindricalSurface) cyl = new Geom_CylindricalSurface(gp_Cylinder(ax, R));
    Handle(Geom_Line) gen = generatorAt(ax, u0, R + off);

    // Precondition (OCCT 7.9.3): the projector frames u into [0, 2π).
    Handle(Geom2d_Curve) raw = GeomProjLib::Curve2d(gen, f, l, cyl);
    expect(!raw.IsNull(), "raw_projection_exists");
    if (raw.IsNull()) return;
    expect(std::fabs(uAt(raw, 0.5 * (f + l)) - (u0 + kTwoPi)) < kEq,
           "raw_projection_is_a_turn_off_the_chart");

    double chart[2] = {0.0, 0.0};
    expect(cylChartInFrame(ax, u0, u1, cyl->Position(), chart), "chart_in_own_frame");
    expect(std::fabs(chart[0] - u0) < kEq && std::fabs(chart[1] - u1) < kEq,
           "chart_in_own_frame_is_the_chart");

    ChartBranchInfo info;
    Handle(Geom2d_Curve) pc = projectOnChartBranch(gen, f, l, cyl, chart, &info);
    expect(!pc.IsNull(), "branch_projection_exists");
    if (pc.IsNull()) return;
    expect(std::strcmp(info.from, "region") == 0, "branch_uses_region_chart");
    expect(std::fabs(info.shift + kTwoPi) < kEq, "branch_shift_is_one_turn");
    expect(std::fabs(uAt(pc, f) - u0) < kEq && std::fabs(uAt(pc, l) - u0) < kEq,
           "branch_u_on_chart_lower_end");
    expect(std::fabs(vAt(pc, f) - vAt(raw, f)) < kEq && std::fabs(vAt(pc, l) - vAt(raw, l)) < kEq,
           "branch_v_unchanged");
    // Exact: the same surface points, so the edge's deviation is unchanged.
    double maxMove = 0.0;
    for (int i = 0; i <= 8; ++i) {
        const double t = f + (l - f) * i / 8.0;
        const gp_Pnt2d a = raw->Value(t), b = pc->Value(t);
        maxMove = std::max(maxMove, cyl->Value(a.X(), a.Y()).Distance(cyl->Value(b.X(), b.Y())));
    }
    expect(maxMove < kEq, "branch_shift_is_exact_on_surface");

    // A generator inside the chart is untouched.
    Handle(Geom_Line) inner = generatorAt(ax, 0.3, R + off);
    Handle(Geom2d_Curve) pcIn = projectOnChartBranch(inner, f, l, cyl, chart, &info);
    expect(!pcIn.IsNull() && info.shift == 0.0 && std::fabs(uAt(pcIn, f) - 0.3) < kEq,
           "interior_generator_unchanged");

    // No chart: the projection is returned as it was.
    Handle(Geom2d_Curve) pcNone = projectOnChartBranch(gen, f, l, cyl, nullptr, &info);
    expect(!pcNone.IsNull() && std::strcmp(info.from, "none") == 0 &&
               std::fabs(uAt(pcNone, f) - uAt(raw, f)) < kEq,
           "no_chart_unchanged");
}

// A sheet trimmed in u is its own chart. OCCT re-frames a trim on a periodic
// basis, so the bounds can sit a turn away from the region's numbers; the
// pcurve must land inside them.
void testTrimmedSheets() {
    const double R = 50.0, off = 0.04, u0 = -0.01, u1 = 0.6, f = 0.0, l = 5.0;
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    Handle(Geom_CylindricalSurface) cyl = new Geom_CylindricalSurface(gp_Cylinder(ax, R));
    Handle(Geom_RectangularTrimmedSurface) rect =
        new Geom_RectangularTrimmedSurface(cyl, u0, u1, f, l);
    Standard_Real b0 = 0, b1 = 0, v0 = 0, v1 = 0;
    rect->Bounds(b0, b1, v0, v1);
    Handle(Geom_Line) gen = generatorAt(ax, u0, R + off);
    ChartBranchInfo info;
    double regionChart[2] = {u0, u1};
    Handle(Geom2d_Curve) pc = projectOnChartBranch(gen, f, l, rect, regionChart, &info);
    expect(!pc.IsNull() && std::strcmp(info.from, "trim") == 0, "rect_trim_uses_own_bounds");
    if (!pc.IsNull()) {
        const double u = uAt(pc, 0.5 * (f + l));
        expect(u >= b0 - kEq && u <= b1 + kEq, "rect_trim_pcurve_inside_bounds");
        expect(modTwoPiDist(u, u0) < kEq, "rect_trim_pcurve_same_azimuth");
    }

    // The rotated sheet [0, span]: a generator a hair below its u = 0.
    gp_Ax3 rotAx = ax;
    rotAx.Rotate(gp_Ax1(ax.Location(), ax.Direction()), u0);
    Handle(Geom_CylindricalSurface) rotCyl = new Geom_CylindricalSurface(gp_Cylinder(rotAx, R));
    Handle(Geom_RectangularTrimmedSurface) rot =
        new Geom_RectangularTrimmedSurface(rotCyl, 0.0, u1 - u0, f, l);
    const double hair = -1e-7;
    Handle(Geom_Line) genRot = generatorAt(rotAx, hair, R + off);
    Handle(Geom2d_Curve) pcRot = projectOnChartBranch(genRot, f, l, rot, nullptr, &info);
    expect(!pcRot.IsNull() && std::strcmp(info.from, "trim") == 0, "rot_trim_uses_own_bounds");
    if (!pcRot.IsNull())
        expect(std::fabs(uAt(pcRot, f) - hair) < kEq, "rot_trim_hair_below_stays_on_chart");
}

// The region chart carried into a turned frame, direct and indirect.
void testFrames() {
    const double u0 = -0.01, u1 = 0.6, R = 50.0;
    for (int indirect = 0; indirect <= 1; ++indirect) {
        gp_Ax3 ax(gp_Pnt(1, 2, 3), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        if (indirect) ax.YReverse();
        for (double theta : {0.0, 0.25, -1.3, 2.9}) {
            gp_Ax3 pos = ax;
            pos.Rotate(gp_Ax1(ax.Location(), ax.Direction()), theta);
            double chart[2] = {0.0, 0.0};
            const bool ok = cylChartInFrame(ax, u0, u1, pos, chart);
            expect(ok, indirect ? "frame_indirect_ok" : "frame_direct_ok");
            if (!ok) continue;
            expect(std::fabs((chart[1] - chart[0]) - (u1 - u0)) < kEq, "frame_span_kept");
            // The chart ends name the same 3D points in either frame.
            gp_Cylinder cb(ax, R);
            const gp_Pnt p0 = ElSLib::Value(u0, 0.0, cb), p1 = ElSLib::Value(u1, 0.0, cb);
            expect(modTwoPiDist(uInFrame(pos, R, p0), chart[0]) < kEq &&
                       modTwoPiDist(uInFrame(pos, R, p1), chart[1]) < kEq,
                   indirect ? "frame_indirect_maps_chart" : "frame_direct_maps_chart");
        }
    }
    // Not a turn of the region's frame: no chart.
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Ax3 flipped(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
    double chart[2] = {0.0, 0.0};
    expect(!cylChartInFrame(ax, u0, u1, flipped, chart), "frame_reversed_axis_none");
    gp_Ax3 other = ax;
    other.YReverse();
    expect(!cylChartInFrame(ax, u0, u1, other, chart), "frame_other_handedness_none");
    expect(!cylChartInFrame(ax, 0.0, kTwoPi, ax, chart), "frame_full_turn_none");
}

// A chart across u = π on the region frame, and the same arc written below
// -π: the branch follows the chart, not the projector's [0, 2π).
void testStraddle() {
    const double R = 50.0, off = 0.04, f = 0.0, l = 5.0;
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    Handle(Geom_CylindricalSurface) cyl = new Geom_CylindricalSurface(gp_Cylinder(ax, R));
    ChartBranchInfo info;
    {
        double chart[2] = {2.9, 3.5};
        Handle(Geom2d_Curve) pc =
            projectOnChartBranch(generatorAt(ax, 3.5, R + off), f, l, cyl, chart, &info);
        expect(!pc.IsNull() && info.shift == 0.0 && std::fabs(uAt(pc, f) - 3.5) < kEq,
               "straddle_upper_end_kept");
    }
    {
        double chart[2] = {-3.5, -2.9};
        Handle(Geom2d_Curve) pc =
            projectOnChartBranch(generatorAt(ax, -2.9, R + off), f, l, cyl, chart, &info);
        expect(!pc.IsNull() && std::fabs(uAt(pc, f) + 2.9) < kEq, "straddle_negative_chart_followed");
    }
}

}  // namespace

int main() {
    testShiftPure();
    testUntrimmedSeam();
    testTrimmedSheets();
    testFrames();
    testStraddle();
    std::printf("pcurve_chart_branch_unit: %d passed, %d failed\n", gPass, gFail);
    return gFail == 0 ? 0 : 1;
}
