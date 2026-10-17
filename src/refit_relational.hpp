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

// D-fc-2.1 / D-fc-5.17 (a). Each azimuth is known to τ/r, a step to 2τ/r, so
// two true-equal steps differ by at most 4τ/r. No other band.
inline double equalStepTol(double tau, double r) {
    return (r > 0.0 && tau > 0.0) ? (4.0 * tau / r) : 0.0;
}

// A rail of k chords is this section's equal-step law when k over-determines
// a free circle and the steps agree within equalStepTol.
inline bool equalStepRail(int k, double minStep, double maxStep, double tau, double r) {
    return k >= minStations(kParamsCircle) && maxStep >= minStep &&
           (maxStep - minStep) <= equalStepTol(tau, r);
}

// D-fc-5.17 (a) / D-fc-1.6.7. A section circle is tangent to one certified
// line anchor when the extrusion is parallel to the plane and the centre
// stands one radius off it.
inline bool sectionTangentLine(double axisDotN, double centreDist, double r, double tolA,
                               double tau, double L) {
    if (!(r > 0.0) || !(L > 0.0)) return false;
    return std::fabs(axisDotN) <= tanAngTol(tolA, tau, L) &&
           std::fabs(std::fabs(centreDist) - r) <= tanDistTol(tolA, tau);
}

// Tangent to one certified circle. side is +1 when the centres are R + r
// apart (external) and -1 when they are |R − r| apart (internal).
inline bool sectionTangentCircle(double axisCross, double centreDist, double R, double r,
                                 double tolA, double tau, double L, int& side) {
    side = 0;
    if (!(r > 0.0) || !(R > 0.0) || !(L > 0.0)) return false;
    if (axisCross > tanAngTol(tolA, tau, L)) return false;
    const double tol = tanDistTol(tolA, tau);
    const double ext = std::fabs(centreDist - (R + r));
    const double inn = std::fabs(centreDist - std::fabs(R - r));
    if (ext <= tol && ext <= inn) {
        side = 1;
        return true;
    }
    if (inn <= tol) {
        side = -1;
        return true;
    }
    return false;
}

// Tangent plane of a section circle at azimuth phi: n · x = d, |n| = 1.
// Zero free parameters (the circle and the end azimuth).
inline void sectionTangentPlane(double cx, double cy, double r, double phi, double& nx,
                                double& ny, double& d) {
    nx = std::cos(phi);
    ny = std::sin(phi);
    d = nx * (cx + r * nx) + ny * (cy + r * ny);
}

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_RELATIONAL_HPP
