// refit::segment — parameter derivation, stage sequencing, threading policy, failure
// handling. Algorithms live in refit_grow / refit_fillet / refit_chains / refit_math.
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "refit.hpp"
#include "refit_internal.hpp"

namespace stl2step {
namespace refit {

// Landed by sibling agents; declared here so runStages can wire them without
// touching refit_internal.hpp / grow / fillet.
bool claimNgonWallsA(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work);
bool claimChamferConesC(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                        SegmentWork& work);

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
    if (!claimNgonWallsA(mv, p, tol, work)) return false;
    if (!claimCylindersB1(mv, p, tol, work)) return false;
    if (!claimChamferConesC(mv, p, tol, work)) return false;

    if (p.doFillets) {
        if (!claimFilletsC1(mv, p, tol, work)) return false;
    }

    if (!commitPlanesA3(mv, p, tol, work)) return false;
    if (!buildTopologyD(mv, p, tol, work, out)) return false;

    if (const char* d130 = std::getenv("STL2STEP_DIAG_130"); d130 && d130[0] && d130[0] != '0') {
        int nNgon = 0, nCone = 0, nCyl = 0, nPl = 0, nFil = 0;
        for (const Region& r : out.regions) {
            if (r.origin == Origin::NgonWall) nNgon++;
            if (r.origin == Origin::ChamferCone) nCone++;
            if (r.type == SurfType::Cylinder) nCyl++;
            if (r.type == SurfType::Plane) nPl++;
            if (r.origin == Origin::FilletStrip) nFil++;
        }
        std::fprintf(stderr,
                     "DIAG_130_CENSUS nTri=%zu nReg=%zu planes=%d cyl=%d ngon=%d "
                     "cone=%d fillet=%d chains=%zu\n",
                     mv.nTri, out.regions.size(), nPl, nCyl, nNgon, nCone, nFil,
                     out.chains.size());
        for (const Region& r : out.regions) {
            if (r.type == SurfType::Cylinder || r.type == SurfType::Cone)
                std::fprintf(stderr,
                             "  DIAG_130_REG id=%d type=%d origin=%d R=%.4f nTri=%zu "
                             "nSides=%d closed360=%d v=[%.4f,%.4f]\n",
                             r.id, (int)r.type, (int)r.origin, r.radius, r.tris.size(),
                             r.nSides, r.closed360 ? 1 : 0, r.vMin, r.vMax);
        }
    }

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
