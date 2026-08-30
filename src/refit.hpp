// stl2step internal surface-refit contract (P1 <-> P2).
// Not part of the public API. Never installed, never included from
// include/stl2step/stl2step.hpp.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_HPP
#define STL2STEP_REFIT_HPP

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_XYZ.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

// P1 guarantees (invariants P2 may rely on without re-checking).
// I1-I6, I8, I9 from DECISION-insertion §3.1. I7 amended and I7b added
// from DECISION-p1-math §3.4.
//
// I1. triRegion[t] >= 0 XOR triIsland[t] >= 0 for every local triangle — a
//     total, non-overlapping partition.
// I2. A mesh edge whose two triangles are in different regions/islands
//     appears in exactly one chain; a mesh edge interior to a region
//     appears in none.
// I3. Chains are ordered and oriented so that walking meshVerts keeps
//     regA on the left.
// I4. Every vertex of every triangle of an accepted region satisfies
//     dev <= maxVertexDev <= eps for that region's fitted surface.
// I5. Region ids, chain ids and island ids are byte-for-byte independent
//     of --threads.
// I6. No accepted region contains a triangle from a dirty component
//     (P1 is never called on one).
// I7 (amended). For a region with closed360 == false, Region::loops is complete
// and closed: exactly one loop with role == LoopRole::Outer, plus zero or more with
// role == LoopRole::Inner. Every chain touching that region appears in exactly one of its
// loops.
// I7b (new). For a region with closed360 == true, Region::loops contains exactly
// one LoopRole::CapLow (at vMin) and exactly one LoopRole::CapHigh (at vMax), plus
// zero or more LoopRole::Inner, and no LoopRole::Outer. Every chain touching that
// region appears in exactly one of its loops. No loop corresponds to the seam.
// I8. Chains are split at every vertex incident to >=3 distinct
//     regions/islands, and at every vertex touching a non-manifold or
//     open mesh edge. (§2.2 — the collapse-safety linchpin.)
// I9. dVolPredicted is signed and physically meaningful: positive where
//     the analytic surface bulges outward relative to the chords (bosses),
//     negative for holes.

namespace stl2step { namespace refit {

// Progress/diagnostic sink. Binds to the engine's Conv::warn(const std::string&).
// MAY BE NULL - always null-check. P1 calls this ONLY when segment() is about to
// return false (R3). Every reject/NYI/faceted outcome is silent and counted, because
// a warning becomes Result::warnings and therefore exit code 2.
using WarnFn = std::function<void(const std::string& msg)>;

struct SegmentParams {
    // --- tolerances, mm; 0 => auto-derive from MeshView -------------------
    double epsMesh       = 0.0;   // 0 -> max(mv.weldTol, 1e-4*mv.diag, 1e-3)
    double epsPlane      = 0.0;   // 0 -> max(epsMesh, mv.sewTol, 0.02)   [Options::smoothTolMM]
    // --- angles, DEGREES --------------------------------------------------
    double thetaPlaneDeg = 2.0;   // Phase A normal gate               [Options::smoothAngleDeg]
    double thetaSharpDeg = 30.0;  // chart / growth hard stop
    double thetaCylLoDeg = 5.0;   // Phase B seed band, INCLUSIVE
    double thetaCylHiDeg = 60.0;  // Phase B seed band, INCLUSIVE       (D5.4)
    double thetaBinDeg   = 0.25;  // nSides band-clustering floor       (D2.2)
    // --- switches ---------------------------------------------------------
    bool   doFillets     = true;  // Phase C                           [Options::smoothFillets]
    // --- determinism (I5): every iterative step is hard-capped ------------
    int    maxRefineIters = 8;    // geometric (c,R) refinement cap
};

enum class SurfType : uint8_t { Plane, Cylinder, /* v2 */ Cone, Sphere, Torus };
enum class Origin   : uint8_t { PlaneGrow, CylGrow, FilletStrip };
enum class BuiltAs  : uint8_t { NotBuilt, Single, Seamed360, TwoHalves, ExplodedToFacets };
enum class Reject   : uint8_t { None, GaussPlanarity, VertexResidual, ChordConsistency,
                                RadiusSanity, Span, FilletConsensus, NeighborNotAnalytic,
                                StripWidth, TorusNYI, ConeNYI, SphereNYI, DirtyComponent,
                                FaceBuildFailed, ChainUnstable };

enum class LoopRole : uint8_t {
    Outer,    // the single outer boundary of a region with closed360 == false
    Inner,    // an inner (hole) boundary; legal on closed360 regions too (cross-hole)
    CapLow,   // closed360 only: the cap loop at v == vMin
    CapHigh   // closed360 only: the cap loop at v == vMax
};

// --- the engine fills this from its existing arrays; refit owns nothing. -------
struct MeshView {
    const gp_XYZ*                        pts;        // global welded points (post-scale, mm)
    const int (*tris)[3];                            // global triangles (a,b,c)
    const int*                           compTris;   // local tri -> global tri id
    const int*                           compVtx;    // local vtx -> global point id
    const std::pair<int,int>*            compEdges;  // local edge -> (local vLo, vHi)
    const std::array<int,3>*             triEdges;   // local tri -> 3 local edge ids
    const uint8_t*                       triDirs;    // local tri -> direction bits
    size_t nTri, nVtx, nEdge;
    double diag, weldTol, sewTol;                    // derived tolerance inputs
};

struct Loop {                       // a face boundary loop, in order
    std::vector<int>     chainIdx;  // indices into RegionSet::chains
    std::vector<uint8_t> reversed;  // 1 = traverse that chain backwards
    LoopRole             role;      // see I7 / I7b; NOT a bool - closed360 has two caps
};

struct Region {
    int      id;                    // dense; deterministic (area desc, then min local tri id)
    SurfType type;
    Origin   origin;                // which stage committed the region:
                                    //   PlaneGrow   - Phase A3 (provisional region unclaimed)
                                    //   CylGrow     - Phase B1 (provisional-region merge)
                                    //   FilletStrip - Phase C1 (<=3 bands, two neighbours)
                                    // Phase A produces a PROVISIONAL oversegmentation; B and
                                    // C claim whole provisional regions; A3 commits the rest.
                                    // G3 runs on CylGrow ONLY.
    // geometry, mm, same frame as pts[]
    gp_Ax3   ax;                    // Plane: point+normal. Cylinder: axis placement, Loc on axis
                                    // Cylinder: Location() = axis point nearest the
                                    // area-weighted vertex centroid; XDirection() is set to
                                    // a REAL MESH VERTEX azimuth (the region vertex with the
                                    // lowest local id), so a 360° seam at u=0 lands on a
                                    // facet generator instead of bisecting a facet.
    double   radius = 0;            // Cylinder only
    double   uMin = 0, uMax = 0;    // Cylinder: angular start + end from ax.XDirection().
                                    // closed360 => exactly 0 and 2*pi. Partial => the
                                    // complement of the largest vertex-azimuth gap.
    double   vMin = 0, vMax = 0;    // axial extent: min/max of (p_i - ax.Location()).a
    bool     closed360 = false;     // nBands >= 3 AND largest vertex-azimuth gap
                                    // <= 1.5*(2*pi/nBands). Drives the §2.4 seam ladder.
    bool     outwardNormal = true;  // false => face must be Reversed() (hole/concave fillet).
                                    // Definition: the material lies on the AXIS side.
                                    // Computed as sign of sum_t area(t)*(n_t . rho_hat(cent_t)),
                                    // rho_hat = unit radial direction from the axis.
                                    // This is exactly the sigma of dVolPredicted (I9).
    bool     lawBand = false;
    // membership
    std::vector<int> tris;          // LOCAL triangle indices
    std::vector<Loop> loops;        // P1 owns loop extraction (mesh topology work).
                                    // closed360 == false: exactly one Outer + N Inner (I7).
                                    // closed360 == true : exactly one CapLow (at vMin) and
                                    //   one CapHigh (at vMax) + N Inner, and NO Outer (I7b).
                                    // No loop is ever the seam - P2 constructs that.
    // quality / budget
    double maxVertexDev = 0, rmsVertexDev = 0;   // region vertices -> fitted surface
    double chordSagitta = 0;        // delta = R * (1 - cos(pi / nSides)); 0 for planes.
                                    // G3's reference length (Origin::CylGrow only).
    int    nSides = 0;              // tessellation count around a FULL circle:
                                    //   nSides = round(2*pi * nBands / span)
                                    // nBands from facet-normal azimuth clustering,
                                    // span from VERTEX azimuths (see DECISION-p1-math D2).
                                    // NEVER round(2*pi / median_dihedral) - that estimator
                                    // is withdrawn (it doubles on quad-split cylinders).
                                    // 0 for planes. For FilletStrip it is the equivalent
                                    // full-circle count, so chordSagitta/dVolPredicted are
                                    // one formula for Phase B and Phase C alike.
    double dVolPredicted = 0;       // SIGNED mm^3 the analytic surface adds vs the chords.
                                    // Cylinder/fillet (exact circular-segment prism):
                                    //   gamma = 2*pi/nSides
                                    //   sigma * A_reg * R * (gamma - sin gamma)
                                    //                     / (4 * sin(gamma/2))
                                    //   sigma = outwardNormal ? +1 : -1
                                    // Plane: -sum_t A_t * (h_a + h_b + h_c)/3, h along
                                    //   ax.Direction(). NOT pi*R^2*H*(delta/R) - that
                                    //   approximation is 0.75x low and is withdrawn.
    double maxVertexSnap = 0;       // max distance a shared vertex was moved (0 if none)
    // diagnostics
    Reject  reject  = Reject::None; // != None => never built; reported only
    BuiltAs builtAs = BuiltAs::NotBuilt;   // written by P2
    int     filletNbrA = -1, filletNbrB = -1;  // Origin::FilletStrip ONLY; -1 otherwise.
                                    // Provenance, not a build instruction: P2 keys the
                                    // constructed-generator rule off BoundaryChain::tangent.
};

struct BoundaryChain {
    int  regA = -1, regB = -1;      // region ids; -1 => faceted island
    int  islandA = -1, islandB = -1;// island ids when the corresponding reg* == -1
    bool tangent    = false;        // G1 => construct the generator, never intersect.
                                    // Geometric test, independent of Origin: both sides are
                                    // accepted analytic regions AND
                                    //   plane|cyl : |dist(axis,plane) - R| <= epsPlane
                                    //               and |a . n_plane| <= sin(3 deg)
                                    //   cyl|cyl   : axes parallel within 3 deg and
                                    //               |dist(axes) - |R1 +/- R2|| <= epsPlane
                                    //   plane|pln : n_L . n_R >= cos(3 deg)
                                    //   any|island: false
    bool closedLoop = false;        // e.g. a hole cap ring
    std::vector<int> meshEdges;     // ordered LOCAL edge ids
    std::vector<int> meshVerts;     // ordered LOCAL vertex ids; size == meshEdges+1 (open) or == meshEdges (closed)
};

struct RefitStats {
    int planes = 0, cylinders = 0, fillets = 0, rejected = 0;
    int facetIslands = 0, facetTriangles = 0;
    int distinctRadii = 0;
    double maxVertexDev = 0, maxEdgeTol = 0, dVolPredicted = 0;
};

struct RegionSet {                  // exactly one per CLEAN component
    int compRoot = -1;
    std::vector<Region>        regions;   // accepted only (reject == None)
    std::vector<Region>        rejected;  // diagnostics; never built
    std::vector<BoundaryChain> chains;
    std::vector<int>           triRegion; // local tri -> region id, else -1
    std::vector<int>           triIsland; // local tri -> island id, else -1
    int                        nIslands = 0;
    RefitStats                 stats;
};

// P1 entry point. Pure math. Links no OCCT topology: gp_ value types and math_ solvers
// ONLY. Banned in refit_segment.cpp: TopoDS_*, BRep*, ShapeFix_*, ShapeUpgrade_*,
// ShapeAnalysis_* (incl. ShapeAnalysis_Geom), IntAna_*, Geom_*/Geom2d_* handles,
// GeomAPI_*, GeomConvert_*, GProp_*, Poly_*, STEPControl_*. PCA/Jacobi and both circle
// fits are written in-TU; no Eigen. P0 greps the include list against this rule.
bool segment(const MeshView& mv, const SegmentParams& p, RegionSet& out, WarnFn warn);

// P2 entry point. Consumes RegionSet, produces the SAME vector<TopoDS_Face> the
// per-triangle loop produces today. Returns false => caller falls back to R2 (§5).
bool buildFaces(const MeshView& mv, RegionSet& rs,
                const std::vector<TopoDS_Vertex>& verts,   // engine's existing slot array
                std::vector<TopoDS_Face>& out, WarnFn warn);
}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_HPP
