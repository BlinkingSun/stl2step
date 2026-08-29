// refit::segment — parameter derivation, stage sequencing, threading policy, failure
// handling. Algorithms live in refit_grow / refit_fillet / refit_chains / refit_math.
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <string>

#include "refit.hpp"
#include "refit_internal.hpp"

namespace stl2step {
namespace refit {
namespace {

constexpr double kDegToRad = M_PI / 180.0;

void adaptCoarseSegmentParams(const MeshView& mv, SegmentParams& p) {
    if (!coarseFusionBand(mv)) return;
    p.thetaPlaneDeg = std::max(p.thetaPlaneDeg, 15.0);
    p.thetaCylHiDeg = std::max(p.thetaCylHiDeg, 70.0);
}

DerivedTols deriveTols(const MeshView& mv, const SegmentParams& p) {
    DerivedTols tol;
    tol.epsMesh = (p.epsMesh > 0.0)
                      ? p.epsMesh
                      : std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
    tol.epsPlane = (p.epsPlane > 0.0) ? p.epsPlane
                                      : std::max({tol.epsMesh, mv.sewTol, 0.02});
    tol.thetaPlane = p.thetaPlaneDeg * kDegToRad;
    tol.thetaSharp = p.thetaSharpDeg * kDegToRad;
    tol.thetaCylLo = p.thetaCylLoDeg * kDegToRad;
    tol.thetaCylHi = p.thetaCylHiDeg * kDegToRad;
    tol.thetaBin = p.thetaBinDeg * kDegToRad;
    return tol;
}

bool runStages(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
               RegionSet& out) {
    SegmentWork work;

    if (!chartsA1(mv, p, tol, work)) return false;

    if (!growProvisionalA2(mv, p, tol, work)) return false;
    if (!claimLawBandsL(mv, p, tol, work)) return false;
    if (!claimCylindersB1(mv, p, tol, work)) return false;

    if (p.doFillets) {
        if (!claimFilletsC1(mv, p, tol, work)) return false;
    }

    if (!commitPlanesA3(mv, p, tol, work)) return false;
    if (!buildTopologyD(mv, p, tol, work, out)) return false;

    return true;
}

}  // namespace

bool segment(const MeshView& mv, const SegmentParams& p, RegionSet& out, WarnFn warn) {
    try {
        out = RegionSet{};

        SegmentParams params = p;
        adaptCoarseSegmentParams(mv, params);
        const DerivedTols tol = deriveTols(mv, params);

        if (!runStages(mv, params, tol, out)) {
            out = RegionSet{};
            if (warn) warn("refit::segment: stage failed");
            return false;
        }

        return true;
    } catch (...) {
        out = RegionSet{};
        if (warn) warn("refit::segment: unknown error");
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
