// stl2step P1 internal cross-TU contract (five translation units).
// Declarations only — FROZEN on landing. Not installed; not public API.
//
// Pipeline (DECISION-p1-math D1 §1.2):
//   A1 charts        refit_grow.cpp   — union-find over edges with |phi| <= theta_sharp
//   A2 provisional   refit_grow.cpp   — running-PCA plane growth per chart; TOTAL partition
//   B1 cyl claim     refit_grow.cpp   — merge whole provisional regions (never triangles)
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
                          int nSides, double spanRad);

double dVolCylinderSector(double areaReg, double radius, int nSides, bool outwardNormal);

double dVolPlaneRegion(const MeshView& mv, const std::vector<int>& tris, const gp_Ax3& ax);

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_INTERNAL_HPP
