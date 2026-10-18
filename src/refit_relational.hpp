// Shared relational rulings (D-fc-1.4, D-fc-1.6). Inline; no algorithm.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_HPP
#define STL2STEP_REFIT_RELATIONAL_HPP

#include "refit_internal.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace stl2step {
namespace refit {

// D-fc-1.6(1). τ = 2q. q <= 0 admits nothing.
inline double tauOf(const MeshView& mv) {
    return mv.quantFloor > 0 ? 2 * mv.quantFloor : 0;
}

// D-fc-1.4. anchorTol(A) = max(τ, A.maxVertexDev).
inline double anchorTol(double tau, double maxVertexDev) {
    return std::max(tau, maxVertexDev);
}

// D-fc-1.4. anchorCertified(A) iff A.maxVertexDev <= τ.
inline bool anchorCertified(double tau, double maxVertexDev) {
    return maxVertexDev <= tau;
}

// D-fc-1.6(2). Model parameter counts. Under-determined when stations < params + 1.
constexpr int kParamsAnchored = 1, kParamsLine = 2, kParamsCircle = 3, kParamsAxisDir = 2;

// D-fc-4.9. One-sided extrusion blend: axis ∥ n_T × n_X and tangency to T
// imposed; free are the contact position along T and r.
constexpr int kParamsOneSided = 2;

// D-fc-5.11 F9. Hem: the half-cylinder tangent to two certified antiparallel
// planes. Axis ∥ n_T × n_X (certified), r = half the gap, centre midway, the
// contact where T's certified body ends: zero free parameters.
constexpr int kParamsHem = 0;

// D-fc-5.17. A horn torus fixes Rmaj = r (the profile circle passes through
// the axis). Free parameters are z0 and r — two, one fewer than the profile
// circle. When r is the certified edge-round radius over 2 (Rmaj = R_edge - r),
// only z0 is free and the count is kParamsAnchored.
constexpr int kParamsHorn = kParamsLine;

// D-fc-1.6(2).
inline int minStations(int params) { return params + 1; }

// D-fc-5.11 F9. T and X' are antiparallel: parallel by the pair walk's test
// (|n_T × n_X'| <= rhs, D-fc-3.2.2) with opposed normals. Co-directed
// parallel planes are a step, never one half-turn.
inline bool hemAntiparallel(double crossTX, double rhs, double dotTX) {
    return crossTX <= rhs && dotTX < 0.0;
}

// D-fc-5.11 F9. Tangent to both planes: r = half the gap, measured from the
// contact on T along the centre side to X'. A gap within τ is one sheet.
inline bool hemRadiusOf(double gap, double tau, double& r) {
    r = 0.5 * gap;
    return gap > tau && std::isfinite(gap);
}

// D-fc-5.11 F9 / D-fc-1.5 O2. The hem turns through π: the vertex span
// matches π within spanTol (D-fc-1.6.7 with dihedral π).
inline bool hemSpansHalfTurn(double span, double spanTolV) {
    constexpr double kHalfTurn = 3.14159265358979323846264338327950288;
    return std::abs(span - kHalfTurn) <= spanTolV;
}

// D-fc-5.5 F-1. One-sided blend slab, "past the contact": a section offset is
// on the band side of the fitted contact generator (gw) OR of the station
// generator, where T's certified body ends (gw0), each with τ slack. Offsets
// are along the contact direction g; sg is the band side (±1).
inline bool osbPastContact(double gw, double gw0, double sg, double tau) {
    return !(sg * gw < -tau) || !(sg * gw0 < -tau);
}

// D-fc-5.5. True when the band (inBand[t] != 0) holds any triangle of `tris`
// -- e.g. of an anchor's certified body, which is that anchor's surface.
inline bool bandHoldsAny(const std::vector<int>& tris, const std::vector<char>& inBand) {
    for (int t : tris)
        if (t >= 0 && (std::size_t)t < inBand.size() && inBand[(std::size_t)t]) return true;
    return false;
}

// D-fc-1.4. Exact predicate of commitPlanesA3 (refit_grow.cpp:2293-2297).
// areaMin = tol.epsPlane * mv.diag (refit_grow.cpp:2198). Shatter-class
// (n <= 2 && area <= areaMin) is not committed unless filletFlank. A
// provisional that is not Unclaimed is not committed.
inline bool a3WouldCommit(const MeshView& mv, const DerivedTols& tol,
                          const Provisional& prov, bool filletFlank) {
    if (prov.claim != ProvClaim::Unclaimed) return false;
    const double areaMin = tol.epsPlane * mv.diag;
    const int n = static_cast<int>(prov.tris.size());
    // refit_grow.cpp:2295
    if (n <= 2 && prov.area <= areaMin && !filletFlank) return false;
    return true;
}

// D-fc-1.6(7). spanTol = sqrt(2·anchorTol(A)/r) + sqrt(2·anchorTol(B)/r).
inline double spanTol(double tolA, double tolB, double r) {
    return std::sqrt(2 * tolA / r) + std::sqrt(2 * tolB / r);
}

// D-fc-4.2. Slab slack at the contact generator on anchor X, as an angle:
// (anchorTol(other) + τ) / (r · |n_A × n_B|).
inline double contactSlack(double tolOther, double tau, double r, double sinDihedral) {
    return (tolOther + tau) / (r * sinDihedral);
}

// D-fc-1.6(7). |dist(axis, A) − r| <= anchorTol(A) + τ.
inline double tanDistTol(double tolA, double tau) { return tolA + tau; }

// D-fc-1.6(7). |a·n_A| <= (anchorTol(A) + τ) / L.
inline double tanAngTol(double tolA, double tau, double L) { return (tolA + tau) / L; }

// A closed wall between two rims is an annulus: two boundary loops (topology).
constexpr int kAnnulusLoops = 2;

// D-fc-5.7 O-1. Winding number about an axis of a closed vertex cycle, from its
// vertex azimuths in traversal order: the wrapped steps summed over 2π. The sum is
// an integer multiple of 2π, read off by rounding; no threshold.
inline int windingNumber(const std::vector<double>& azimuths) {
    const double kTwoPi = 6.283185307179586476925286766559;
    if (azimuths.size() < 3) return 0;
    double sum = 0;
    for (std::size_t i = 0; i < azimuths.size(); ++i) {
        double d = azimuths[(i + 1) % azimuths.size()] - azimuths[i];
        d = std::remainder(d, kTwoPi);
        sum += d;
    }
    return static_cast<int>(std::lround(sum / kTwoPi));
}

// D-fc-5.7 O-1. A band is a closed wall by construction when it is ONE edge-component
// bounded by kAnnulusLoops cycles, each winding once about the axis: every generator
// crosses it, whatever its tessellation (D2's band count cannot see that on a wall
// tessellated coarsely at one rim and finely at the other).
inline bool annulusAboutAxis(int edgeComponents, const std::vector<int>& cycleWindings) {
    if (edgeComponents != 1 || (int)cycleWindings.size() != kAnnulusLoops) return false;
    for (int w : cycleWindings)
        if (w != 1 && w != -1) return false;
    return true;
}

// D-fc-4.3. An accepted cylinder region the extrude stage finds on entry whose
// own maxVertexDev exceeds τ does not certify its surface (D-fc-1.4: B1
// phantoms, and plane|plane strips of the segment's fillet stage that miss τ),
// so a certified relational band may rob it; the robbed region dissolves. A
// region within τ is never touched.
inline bool robbableUncertified(const Region& r, double tau) {
    return (r.origin == Origin::CylGrow || r.origin == Origin::FilletStrip) &&
           r.type == SurfType::Cylinder && r.maxVertexDev > tau;
}

// D-fc-4.3 with D-fc-1.6.10. A prior strip the band robs dissolves at its commit,
// and a strip an earlier band dissolved holds no triangle: neither is a blend
// partner at a corner of the band. robs is sorted.
inline bool blendPartnerGone(const std::vector<int>& robs, int ai, const Region& acc) {
    return acc.tris.empty() || std::binary_search(robs.begin(), robs.end(), ai);
}

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_RELATIONAL_HPP
