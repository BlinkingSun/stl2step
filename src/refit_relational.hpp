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

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_RELATIONAL_HPP
