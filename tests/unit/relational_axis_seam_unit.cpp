// Synthetic cases for the torus-family seam predicates (D-fc-4.8), the closed-wall /
// closed-torus build predicates of fx-closed360 and the U-3 naming predicates. No mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_axis_seam.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

#include <Geom_Circle.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kTau = 1e-4;

static gp_Pnt onCircle(double R, double u, double z) {
    return gp_Pnt(R * std::cos(u), R * std::sin(u), z);
}

int main() {
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
    }

    // splitCapHealApplies: certified at tau and one exact rim, both required.
    expect(splitCapHealApplies(0.5 * kTau, kTau, true), "split-cap: certified wall with an exact rim");
    expect(!splitCapHealApplies(2.0 * kTau, kTau, true), "split-cap: uncertified wall keeps its path");
    expect(!splitCapHealApplies(0.5 * kTau, kTau, false), "split-cap: all-polyline rims keep their path");
    expect(!splitCapHealApplies(0.0, 0.0, true), "split-cap: tau 0 admits nothing");

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

    std::fprintf(stderr, "relational_axis_seam_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
