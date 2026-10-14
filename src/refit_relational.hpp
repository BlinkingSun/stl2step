// Shared relational rulings (D-fc-1.4, D-fc-1.6). Inline; no algorithm.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_HPP
#define STL2STEP_REFIT_RELATIONAL_HPP

#include "refit_internal.hpp"

#include <algorithm>
#include <cmath>

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

// D-fc-1.6(2).
inline int minStations(int params) { return params + 1; }

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

// D-fc-1.6(7). |dist(axis, A) − r| <= anchorTol(A) + τ.
inline double tanDistTol(double tolA, double tau) { return tolA + tau; }

// D-fc-1.6(7). |a·n_A| <= (anchorTol(A) + τ) / L.
inline double tanAngTol(double tolA, double tau, double L) { return (tolA + tau) / L; }

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_RELATIONAL_HPP
