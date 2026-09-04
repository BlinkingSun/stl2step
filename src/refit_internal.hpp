// stl2step P1 internal cross-TU contract (five translation units).
// Declarations only — FROZEN on landing. Not installed; not public API.
//
// Pipeline (DECISION-p1-math D1 §1.2):
//   A1 charts        refit_grow.cpp   — union-find over edges with |phi| <= theta_sharp
//   A2 provisional   refit_grow.cpp   — running-PCA plane growth per chart; TOTAL partition
//   L  law bands     refit_grow.cpp   — Stage L after A2, before A(ngon)
//   A  ngon walls    (new TU)         — claimNgonWallsA; Origin::NgonWall
//   B1 cyl claim     refit_grow.cpp   — merge whole provisional regions (never triangles)
//   C  chamfer cone  (new TU)         — claimChamferConesC; Origin::ChamferCone
//   C1 fillet claim  refit_fillet.cpp — strips of <= 3 unclaimed provisionals, D7 radius
//   A3 plane commit  refit_grow.cpp   — every still-unclaimed provisional -> PlaneGrow plane
//   D  topology      refit_chains.cpp — chains, I8, loops/LoopRole, ids, RefitStats
//   glue             refit_segment.cpp — resolve SegmentParams/DerivedTols, sequencing, threads
//   math             refit_math.cpp   — Jacobi/PCA/Eberly/Pratt/Gauss-map/dVol (no Eigen)
//
// B1 claims whole provisional regions, not individual triangles, because the
// discriminating vertex-residual signal is maximal on a designed tangent plane
// (large extent -> immediate blow-up) and ~0 on a true cylinder band (D1 §1.2).
// RegionSet::regions is not materialised until stage D, so no published region id
// is ever rewritten.
//
// DerivedTols: segment() resolves epsMesh/epsPlane/angles once (D5.2) into a
// DerivedTols passed to every stage entry point. R-dependent tolerances
// (epsCylGrow, epsCylAccept, epsFillet) are inline helpers on that struct, not
// stored fields — they stay out of SegmentParams per D5.2.
//
// Determinism (I5) — obligation on every stage before a TU boundary:
//   * sort every produced collection by the D5.6 keys (each key ends in an id);
//   * accumulate every floating-point reduction in ascending local triangle id;
//   * never sort on a key that can tie without a trailing id field.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_INTERNAL_HPP
#define STL2STEP_REFIT_INTERNAL_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include "refit.hpp"

namespace stl2step { namespace refit {

// --- A2 currency; internal only — no Region id, never dumped -----------------

enum class ProvClaim : uint8_t {
    Unclaimed,          // eligible for B1 seed/grow, C1, or A3
    InCylinderClaim,    // transient B1 member; rollback restores Unclaimed
    ConsumedCylinder,   // B1 committed; permanently ineligible
    InFilletClaim,      // transient C1 member; rollback restores Unclaimed
    ConsumedFillet,     // C1 committed; permanently ineligible
    CommittedPlane      // A3 committed as PlaneGrow
};

struct Provisional {
    int chartId = -1;
    std::vector<int> tris;          // LOCAL triangle ids, ascending (I5)
    gp_Ax3 plane;                   // A2 running-PCA plane (point + normal)
    double area = 0;
    double maxVertexDev = 0;
    double rmsVertexDev = 0;
    ProvClaim claim = ProvClaim::Unclaimed;
    bool seedTried = false;         // B1: after a failed claim involving this region,
                                    // may be absorbed but never seed again (D1 §1.3)
};

// D5.2 derived-not-stored constants — computed once in refit_segment.cpp, read-only
// everywhere else. R-dependent helpers are functions, not cached fields.
struct DerivedTols {
    double epsMesh = 0;
    double epsPlane = 0;
    double thetaPlane = 0;    // radians
    double thetaSharp = 0;
    double thetaCylLo = 0;
    double thetaCylHi = 0;
    double thetaBin = 0;

    static constexpr double kGaussPlanarity = 0.05;
    static constexpr double kG3Lo = 0.35;
    static constexpr double kG3Hi = 2.00;
    static constexpr double kG5SpanClosedDeg = 300.0;
    static constexpr double kG5SpanPartialDeg = 40.0;
    static constexpr int    kG5NSidesMin = 6;
    static constexpr int    kG5NBandsMin = 4;

    static inline double gaussAxisTiltSin() {
        return std::sin(3.0 * M_PI / 180.0);
    }
    inline double epsCylGrow(double R) const {
        return std::max(epsMesh, 0.05 * R);
    }
    inline double epsCylAccept(double R) const {
        return std::max(epsMesh, 0.01 * R);
    }
    inline double epsFilletTol(double R) const {
        return std::max(epsMesh, 0.05 * R);
    }
};

// Coarse Fusion/STLB export band (handle-lock @ 908 tris). Lane B coarse gates
// and adaptCoarseSegmentParams must stay in sync — unscoped B (nTri<=1200 or
// global evaluateCommit/A3 changes) breaks p2real S13/S14 fillet fixtures.
// DO NOT WIDEN: S13/S14 (16/20 tris, fillet R=2/2.5) sit outside this band.
inline bool coarseFusionBand(const MeshView& mv) {
    return mv.nTri >= 500 && mv.nTri <= 1200;
}

// Arch-chain detector applicability (chaingen-math). Overlaps coarseFusionBand
// on [500, 1200]; upper bound ~8000 excludes Body11 file (15300 tris) and its
// 12060-tri body. S13/S14 stay outside both predicates. Grow commit remains
// gated by coarseFusionBand until the grow lane wires this predicate.
inline bool archChainBand(const MeshView& mv) {
    return mv.nTri >= 500 && mv.nTri <= 8000;
}

// Working state threaded A1 -> D. Region ids and loops are filled in buildTopologyD.
struct SegmentWork {
    std::vector<int> triChart;          // size mv.nTri
    int nCharts = 0;
    std::vector<Provisional> provisionals;
    std::vector<Region> accepted;       // geometry from B1/C1/A3; topology from D
    std::vector<Region> rejected;       // B1 G1-G4 rollbacks only (G5 silent)
};

// --- Stage entry points (one per D1 §1.2 stage; bool = stage ok, never throws) -

bool chartsA1(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
              SegmentWork& work);

bool growProvisionalA2(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                       SegmentWork& work);

bool claimCylindersB1(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                      SegmentWork& work);

bool claimFilletsC1(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                    SegmentWork& work);

bool commitPlanesA3(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                    SegmentWork& work);

bool buildTopologyD(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                    SegmentWork& work, RegionSet& out);

// --- refit_math.cpp primitives (written in-TU; no Eigen) ---------------------

bool jacobiEigenSymmetric3(const std::array<std::array<double, 3>, 3>& m,
                           std::array<double, 3>& eval,
                           std::array<std::array<double, 3>, 3>& evec);

bool pcaPlane(const MeshView& mv, const std::vector<int>& tris, gp_Ax3& plane);

bool eberlyCenterRadius(const MeshView& mv, const std::vector<int>& tris,
                        const gp_Dir& axis, gp_Pnt& center, double& radius);

// Geometric least-squares cylinder over the unique vertices of `tris`, with the
// axis DIRECTION among the unknowns (D-130-12). `axisSeed` is the starting
// direction and fixes the returned hemisphere; `centerOut` is the axis point
// nearest the vertex mean; `maxResidOut` is the max |radial distance - radius|
// over those vertices. Levenberg-Marquardt, hard-capped at `maxIters` accepted
// steps (I5). Returns false rather than a plausible-looking answer.
bool cylinderFitLS(const MeshView& mv, const std::vector<int>& tris, const gp_Dir& axisSeed,
                   gp_Dir& axisOut, gp_Pnt& centerOut, double& radiusOut, double& maxResidOut,
                   int maxIters = 16);

bool prattCircleFit(const gp_Pnt* points, std::size_t n, double spanRad,
                    gp_Pnt& center, double& radius);

bool gaussMapAxis(const MeshView& mv, const std::vector<int>& tris, gp_Dir& axis);

double chordSagitta(double radius, int nSides);

// Chord-sagitta / inscribed-polygon radius recovery (refit_math.cpp).
double radiusFromChordLength(double chordLen, int nSides);
double circumradiusFromInscribed(double rInscribed, int nSides);
double radiusFromChordSagitta(double halfChord, double sagitta);
int estimateFullCircleSides(const MeshView& mv, const std::vector<int>& tris);
bool refineCylinderRadius(const MeshView& mv, const std::vector<int>& tris,
                          const gp_Dir& axis, gp_Pnt& center, double& radius,
                          int nSides, double spanRad, double rHint = 0.0);

double dVolCylinderSector(double areaReg, double radius, int nSides, bool outwardNormal);

double dVolPlaneRegion(const MeshView& mv, const std::vector<int>& tris, const gp_Ax3& ax);

// Large-R arc strip (lane F): monotonic normal rotation about a common axis, or a
// static-normal vertex ring on coarse tessellation. Used by peelLargeArcStripsA2b.
struct ArcStripDetect {
    bool ok = false;
    gp_Dir axis;
    gp_Pnt center;
    double radius = 0;
    double spanRad = 0;
    bool staticNormals = false;
    double chainScore = 0.0;  // arch-chain signal strength (lane ARCHCHAINS)
    bool fromArchChain = false;
    double areaCV = 0.0;      // relative area CV of the ordered chain
    double angCV = 0.0;       // relative dihedral CV of arc links
    int chainN = 0;           // triangle count of the ordered chain
};

bool detectLargeArcStrip(const MeshView& mv, const std::vector<int>& tris,
                         const DerivedTols& tol, ArcStripDetect& out);

// Point-to-point arch chain (lane ARCHCHAINS): equal-area strip chain with uniform
// dihedral steps; R = w/(2 sin(θ/2)). Returns chainScore in out.chainScore.
bool detectArchChain(const MeshView& mv, const std::vector<int>& tris,
                     const DerivedTols& tol, ArcStripDetect& out);

// Core chain-chord radius from an ordered tri path; axis filters circumferential edges.
bool radiusFromArchChain(const MeshView& mv, const std::vector<int>& chain,
                         const gp_Dir& axis, double& radiusOut, double& chainScoreOut,
                         double rHint = 0.0);

// Build path chain from patch tris and compute arch-chain radius (evaluateCommit hook).
bool archChainRadiusFromPatch(const MeshView& mv, const std::vector<int>& tris,
                              const gp_Dir& axis, double& radiusOut,
                              double& chainScoreOut, double rHint = 0.0);

bool peelLargeArcStripsA2b(const MeshView& mv, const DerivedTols& tol, SegmentWork& work);

// src/refit_internal.hpp  — Stage L (law-band recognition), D4 §1-2.
// d and alpha are PER-EXPORT unknowns: intervals, never points (RULE 4.2b).
struct TessLawInterval {
    double dLo = 0.0, dHi = 0.0;        // mm, chordal deviation; empty when dHi <= dLo
    double alphaLo = 0.0, alphaHi = 0.0;// rad, adjacent-normal cap; alphaHi<=alphaLo => unconstrained
    int    nDLimited = 0, nAlphaLimited = 0;
    bool   empty = true;                // RULE 4.2c: true => decline wholesale
};

struct LawBand {
    std::vector<int>    tris;           // owning triangles
    std::vector<double> theta;          // per-strip generator angle, rad
    std::vector<double> w;              // per-strip circumferential chord, mm
    double  R = 0.0;                    // median w_i/(2 sin(theta_i/2))  -- RULE 4.1d
    gp_Ax1  axis;                       // from Delta-theta ~ 0 generator edges
    double  phi = 0.0;                  // sum theta_i, or 2*pi when closed
    int     N = 0;                      // strips == theta.size()
    bool    closed360 = false;
    double  cvTheta = 0.0, cvR = 0.0, maxVertResid = 0.0;  // Tier-1 residuals
    bool    lowConfidence = false;      // RULE 4.2d: N == 2
};

// Tier 1 -- parameter-free. No d, no alpha, no degree threshold.
bool lawChainAccept(const MeshView& mv, const std::vector<int>& tris,
                    const DerivedTols& tol, LawBand& out);

// Tier 2 -- constraint intersection over accepted bands. Never gates Tier 1.
TessLawInterval lawCalibrate(const std::vector<LawBand>& bands);

// Advisory only (RULE 4.2b). Returns the SET of N consistent with the interval.
bool lawNConsistent(const LawBand& b, const TessLawInterval& li);

// RULE 4.1b/4.1c. Merge iff coaxial + R-consistent + shared generator + still equal-theta.
bool lawBandsMergeable(const LawBand& a, const LawBand& b, const DerivedTols& tol);

// Stage L entry (owned by L2). Runs AFTER growProvisionalA2, BEFORE claimCylindersB1
// and commitPlanesA3. Gated to archChainBand(mv). RULE 4.1a.
bool claimLawBandsL(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                    SegmentWork& work);

bool claimNgonWallsA(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work);
bool claimChamferConesC(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                        SegmentWork& work);

bool tryPlaneLoopCircles(RegionSet& rs, const MeshView& mv, double sewTol);

// Pratt / LS circle in the owning plane. N≥6 unique loop verts; accept iff
// max vertex residual ≤ max(sewTol, chordSagitta(R,N)). Stadium/slot mixed
// loops fail residual. Builder (agent 06) calls this during collapse to emit
// Geom_Circle; do not force-polyline when this returns false if IntAna already
// produced a curve.
bool loopIsCircle(const RegionSet& rs, const Loop& loop, const MeshView& mv,
                  gp_Ax2& ax, double& R, double sewTol);

// Cache filled by tryPlaneLoopCircles (thread_local). loopIndex indexes
// Region::loops. Returns false if that loop was not tagged.
bool planeLoopCircleOf(int regionId, int loopIndex, gp_Ax2& ax, double& R);

// Same cache, keyed by RegionSet::chains index. A multi-chain circular loop
// shares one (ax,R); each chain is an arc of that circle, not a full circle.
bool planeLoopCircleForChain(int chainIndex, gp_Ax2& ax, double& R);

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_INTERNAL_HPP
