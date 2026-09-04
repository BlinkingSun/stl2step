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

    // D-130-11(1) cross-check: the same floor measured on the plane regions that
    // actually reached the RegionSet, so the stage-L number (taken from the A2
    // provisionals, which is where the absorb needs it) can be compared with the
    // set the decision names rather than assumed equal to it.
    if (const char* lb = std::getenv("STL2STEP_LAWBAND_DIAG"); lb && lb[0] && lb[0] != '0') {
        std::vector<double> devs;
        for (const Region& r : out.regions)
            if (r.type == SurfType::Plane && std::isfinite(r.maxVertexDev))
                devs.push_back(r.maxVertexDev);
        std::sort(devs.begin(), devs.end());
        std::fprintf(stderr,
                     "DIAG_NOISEFLOOR_ACCEPTED nTri=%zu planes=%zu floor=%.6g median=%.6g\n",
                     mv.nTri, devs.size(), devs.empty() ? 0.0 : devs.back(),
                     devs.empty() ? 0.0 : devs[devs.size() / 2]);
    }

    if (const char* d130 = std::getenv("STL2STEP_DIAG_130"); d130 && d130[0] && d130[0] != '0') {
        int nNgon = 0, nCone = 0, nCyl = 0, nPl = 0, nFil = 0, nTorus = 0;
        for (const Region& r : out.regions) {
            if (r.origin == Origin::NgonWall) nNgon++;
            if (r.origin == Origin::ChamferCone) nCone++;
            if (r.type == SurfType::Cylinder) nCyl++;
            if (r.type == SurfType::Plane) nPl++;
            if (r.origin == Origin::FilletStrip) nFil++;
            // D-130-9(a): a strip whose neighbour is a RECOGNISED cylinder is a
            // sliver lying on a torus, not a fillet. C1 already refuses it
            // (nbrIsCylinder -> Reject::TorusNYI); the count is what 1.4
            // inherits, so it is reported rather than left to be inferred from
            // the difference between two fillet counts.
            if (r.reject == Reject::TorusNYI) nTorus++;
        }
        for (const Region& r : work.rejected)
            if (r.reject == Reject::TorusNYI) nTorus++;
        std::fprintf(stderr,
                     "DIAG_130_CENSUS nTri=%zu nReg=%zu planes=%d cyl=%d ngon=%d "
                     "cone=%d fillet=%d torusNYI=%d chains=%zu\n",
                     mv.nTri, out.regions.size(), nPl, nCyl, nNgon, nCone, nFil, nTorus,
                     out.chains.size());
        for (const Region& r : work.rejected) {
            if (r.reject != Reject::TorusNYI) continue;
            std::fprintf(stderr, "  DIAG_130_TORUSNYI nTri=%zu firstTri=%d nbrA=%d nbrB=%d\n",
                         r.tris.size(), r.tris.empty() ? -1 : r.tris.front(), r.filletNbrA,
                         r.filletNbrB);
        }
        {
            static const char* kRejectName[] = {
                "None",           "GaussPlanarity", "VertexResidual", "ChordConsistency",
                "RadiusSanity",   "Span",           "FilletConsensus", "NeighborNotAnalytic",
                "StripWidth",     "TorusNYI",       "ConeNYI",        "SphereNYI",
                "DirtyComponent", "FaceBuildFailed", "ChainUnstable"};
            int hist[15] = {0};
            for (const Region& r : work.rejected) {
                const int i = (int)r.reject;
                if (i >= 0 && i < 15) hist[i]++;
            }
            std::fprintf(stderr, "DIAG_130_REJECT nRejected=%zu", work.rejected.size());
            for (int i = 0; i < 15; i++)
                if (hist[i]) std::fprintf(stderr, " %s=%d", kRejectName[i], hist[i]);
            std::fprintf(stderr, "\n");
            // Where every triangle of the mesh ended up. A region count alone
            // cannot say whether a strip was REFUSED or simply never offered:
            // both read as "fewer fillets".
            size_t tPl = 0, tCyl = 0, tCone = 0, tFil = 0, tRej = 0;
            std::vector<char> seen(mv.nTri, 0);
            for (const Region& r : out.regions) {
                for (int t : r.tris)
                    if (t >= 0 && (size_t)t < mv.nTri) seen[(size_t)t] = 1;
                if (r.origin == Origin::FilletStrip) tFil += r.tris.size();
                else if (r.origin == Origin::ChamferCone) tCone += r.tris.size();
                else if (r.type == SurfType::Cylinder) tCyl += r.tris.size();
                else if (r.type == SurfType::Plane) tPl += r.tris.size();
            }
            for (const Region& r : work.rejected) tRej += r.tris.size();
            size_t tFree = 0;
            for (size_t t = 0; t < mv.nTri; t++)
                if (!seen[t]) tFree++;
            std::fprintf(stderr,
                         "DIAG_130_TRICLASS nTri=%zu inPlane=%zu inCyl=%zu inCone=%zu "
                         "inFillet=%zu inNoRegion=%zu rejectedTris=%zu\n",
                         mv.nTri, tPl, tCyl, tCone, tFil, tFree, tRej);
            // One character per triangle, so two runs can be diffed triangle by
            // triangle instead of by aggregate: P plane, C cylinder, K cone,
            // F fillet strip, '.' no region.
            std::string map(mv.nTri, '.');
            for (const Region& r : out.regions) {
                const char ch = (r.origin == Origin::FilletStrip)   ? 'F'
                                : (r.origin == Origin::ChamferCone) ? 'K'
                                : (r.type == SurfType::Cylinder)    ? 'C'
                                : (r.type == SurfType::Plane)       ? 'P'
                                                                    : '?';
                for (int t : r.tris)
                    if (t >= 0 && (size_t)t < mv.nTri) map[(size_t)t] = ch;
            }
            std::fprintf(stderr, "DIAG_130_TRIMAP %s\n", map.c_str());
            for (const Region& r : out.regions)
                std::fprintf(stderr, "  DIAG_130_ALLREG id=%d origin=%d type=%d R=%.4f nTri=%zu "
                                     "t0=%d\n",
                             r.id, (int)r.origin, (int)r.type, r.radius, r.tris.size(),
                             r.tris.empty() ? -1 : r.tris.front());
            for (const Region& r : out.regions) {
                if (r.origin != Origin::FilletStrip) continue;
                std::fprintf(stderr,
                             "  DIAG_130_FIL id=%d R=%.4f nTri=%zu t0=%d nbrA=%d nbrB=%d\n",
                             r.id, r.radius, r.tris.size(), r.tris.empty() ? -1 : r.tris.front(),
                             r.filletNbrA, r.filletNbrB);
            }
        }
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
