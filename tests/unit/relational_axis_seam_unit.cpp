// Synthetic cases for the torus-family seam predicates (D-fc-4.8), the closed-wall /
// closed-torus build predicates of fx-closed360, the crossing-cap predicates of
// D-fc-5.7 O-1, the U-3 naming predicates and the closed-ellipse rim construction of
// D-fc-5.13 and the D-fc-5.20 refusal of a non-positive quantum. No mesh file.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_axis_seam.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
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

    std::fprintf(stderr, "relational_axis_seam_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
