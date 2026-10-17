// Synthetic cases for the relational cap-arc / corner-face build predicates
// (refit_cap_arc.hpp) and the D-fc-4.3 robbing predicates (refit_relational.hpp).
// No mesh, no STL paths.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_cap_arc.hpp"
#include "refit_relational.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

static gp_Pnt onCircle(double R, double u) { return gp_Pnt(R * std::cos(u), R * std::sin(u), 0.0); }

int main() {
    const double R = 3.0;
    const gp_Circ circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R);
    auto azimuth = [](const gp_Pnt& p) { return std::atan2(p.Y(), p.X()); };

    // shortArcInsideSpanOffMid. A strip of two edge-components fitted as one
    // region over u in [-1.5, 0.06]: the second component's one-edge cap spans
    // [0, 0.06] and the span midpoint (-0.72) is off that chord.
    {
        const double u0 = -1.5, u1 = 0.06;
        const gp_Pnt pa = onCircle(R, 0.0), pb = onCircle(R, 0.06);
        const gp_Pnt hint = onCircle(R, 0.5 * (u0 + u1));
        gp_Pnt m;
        const bool got = shortArcInsideSpanOffMid(circ, pa, pb, hint, u0, u1, azimuth, m);
        expect(got && std::fabs(azimuth(m) - 0.03) < 1e-9,
               "capArc: span midpoint in the gap between components -> the chord's short arc");
        gp_Pnt m2;
        const bool got2 = shortArcInsideSpanOffMid(circ, pb, pa, hint, u0, u1, azimuth, m2);
        expect(got2 && m2.Distance(m) < 1e-12, "capArc: terminal order does not matter");
    }
    {
        // The whole cap of a one-component region is the chord: the span
        // midpoint is on its short arc and stands.
        const double u0 = 0.2, u1 = 0.4;
        const gp_Pnt hint = onCircle(R, 0.3);
        gp_Pnt m;
        expect(!shortArcInsideSpanOffMid(circ, onCircle(R, u0), onCircle(R, u1), hint, u0, u1,
                                         azimuth, m),
               "capArc: midpoint on the short arc -> hint stands");
    }
    {
        // A 300 deg region whose cap chord joins its two span ends: the short
        // arc is the 60 deg gap outside the span, so the hint (on the long arc)
        // stands.
        const double u0 = 0.0, u1 = 300.0 * kPi / 180.0;
        const gp_Pnt hint = onCircle(R, 0.5 * (u0 + u1));
        gp_Pnt m;
        expect(!shortArcInsideSpanOffMid(circ, onCircle(R, u0), onCircle(R, u1), hint, u0, u1,
                                         azimuth, m),
               "capArc: short arc outside the fitted span -> hint stands");
    }
    {
        gp_Pnt m;
        const gp_Pnt p = onCircle(R, 1.0);
        expect(!shortArcInsideSpanOffMid(circ, p, p, onCircle(R, 2.0), 0.0, 3.0, azimuth, m),
               "capArc: coincident terminals -> hint stands");
    }
    {
        // A span that crosses the azimuth branch cut at pi: [2.9, 3.5] unwrapped.
        const double u0 = 2.9, u1 = 3.5;
        const gp_Pnt pa = onCircle(R, 3.4), pb = onCircle(R, 3.45);
        const gp_Pnt hint = onCircle(R, 0.5 * (u0 + u1));
        gp_Pnt m;
        const bool got = shortArcInsideSpanOffMid(circ, pa, pb, hint, u0, u1, azimuth, m);
        expect(got && m.Distance(onCircle(R, 3.425)) < 1e-9,
               "capArc: span across the branch cut is read unwrapped");
    }

    // partialTorusWireRetry: only a face whose sole failures are its own
    // orientation statuses is retried.
    expect(partialTorusWireRetry(1, 0), "torusRetry: orientation only -> retry");
    expect(!partialTorusWireRetry(0, 0), "torusRetry: clean face -> no retry");
    expect(!partialTorusWireRetry(1, 1), "torusRetry: orientation plus an edge status -> no retry");
    expect(!partialTorusWireRetry(0, 2), "torusRetry: edge statuses only -> no retry");

    // faceCheckSplit on a horn-torus corner patch (major = minor = r: the blend
    // where a round of radius 2r meets an end round of radius r), u in
    // [0, 70 deg], v in [-pi/2, 0]: the face on its own outer wire checks clean,
    // the face on the reversed wire fails on the face's orientation only, and
    // the reversed wire's reverse checks clean again.
    {
        Handle(Geom_ToroidalSurface) tor =
            new Geom_ToroidalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), 3.0, 3.0 - 1e-6);
        const double uSpan = 70.0 * kPi / 180.0;
        Handle(Geom_RectangularTrimmedSurface) trim =
            new Geom_RectangularTrimmedSurface(tor, 0.0, uSpan, -0.5 * kPi, 0.0);
        BRepBuilderAPI_MakeFace mf0(trim, 1e-7);
        int o0 = -1, b0 = -1, o1 = -1, b1 = -1, o2 = -1, b2 = -1;
        if (mf0.IsDone()) {
            const TopoDS_Face f0 = mf0.Face();
            faceCheckSplit(f0, o0, b0);
            // The wire's edges keep their pcurves on `trim`; a face that takes the
            // wire as given (no MakeFace re-orientation) is the builder's case.
            const TopoDS_Wire w = BRepTools::OuterWire(f0);
            auto faceOn = [&](const TopoDS_Wire& ww) {
                BRep_Builder B;
                TopoDS_Face f;
                B.MakeFace(f, trim, 1e-7);
                B.Add(f, ww);
                return f;
            };
            faceCheckSplit(faceOn(TopoDS::Wire(w.Reversed())), o1, b1);
            faceCheckSplit(faceOn(w), o2, b2);
        }
        std::fprintf(stderr, "faceCheckSplit own=(%d,%d) reversed=(%d,%d) again=(%d,%d)\n", o0, b0,
                     o1, b1, o2, b2);
        expect(o0 == 0 && b0 == 0, "faceCheckSplit: corner patch on its own wire is clean");
        expect(partialTorusWireRetry(o1, b1),
               "faceCheckSplit: reversed wire fails on orientation only -> retried");
        expect(o2 == 0 && b2 == 0, "faceCheckSplit: the retry's wire checks clean");
    }
    {
        int o = -1, b = -1;
        faceCheckSplit(TopoDS_Face(), o, b);
        expect(o == 0 && b > 0, "faceCheckSplit: a null face is not an orientation-only failure");
    }

    // D-fc-4.3 robbableUncertified: CylGrow and FilletStrip cylinders above tau.
    {
        const double tau = 2.0e-4;
        auto reg = [&](Origin o, SurfType t, double dev) {
            Region r;
            r.origin = o;
            r.type = t;
            r.maxVertexDev = dev;
            return r;
        };
        expect(robbableUncertified(reg(Origin::CylGrow, SurfType::Cylinder, 2.0 * tau), tau),
               "robbable: uncertified CylGrow");
        expect(robbableUncertified(reg(Origin::FilletStrip, SurfType::Cylinder, 50.0 * tau), tau),
               "robbable: uncertified FilletStrip");
        expect(!robbableUncertified(reg(Origin::FilletStrip, SurfType::Cylinder, tau), tau),
               "robbable: FilletStrip at tau is certified");
        expect(!robbableUncertified(reg(Origin::CylGrow, SurfType::Cylinder, 0.5 * tau), tau),
               "robbable: certified CylGrow is never touched");
        expect(!robbableUncertified(reg(Origin::ChamferCone, SurfType::Cone, 9.0 * tau), tau),
               "robbable: a cone is not a robbable cylinder");
    }

    // blendPartnerGone: robbed by this band, or dissolved by an earlier one.
    {
        Region live;
        live.tris = {4, 5};
        Region gone;
        const std::vector<int> robs = {2, 7};
        expect(blendPartnerGone(robs, 7, live), "partnerGone: robbed by this band");
        expect(blendPartnerGone(robs, 3, gone), "partnerGone: dissolved earlier");
        expect(!blendPartnerGone(robs, 3, live), "partnerGone: a live prior strip stays a partner");
        expect(!blendPartnerGone({}, 3, live), "partnerGone: no robs, live strip");
    }

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_capface_unit ok\n");
    return 0;
}
