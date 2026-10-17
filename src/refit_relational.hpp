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

// D-fc-1.6(2).
inline int minStations(int params) { return params + 1; }

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
