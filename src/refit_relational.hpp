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

// D-fc-1.6.7. Branch-cut on a circle: curve parameter u (radians) sits on k·2π
// within arc-length τ when |u − k·2π| ≤ τ/r (parametric confusion on the natural
// circle chart, not a hand epsilon).
inline double circlePeriodParamTol(double tau, double radius) {
    const double r = std::max(radius, 0.0);
    return (tau > 0.0 && r > 0.0) ? (tau / r) : 0.0;
}

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
// region within τ is never touched. D-gid-7.6: nor is a wall that passed the
// ring certificate; its maxVertexDev is the chord sag that certificate excuses.
inline bool robbableUncertified(const Region& r, double tau) {
    return (r.origin == Origin::CylGrow || r.origin == Origin::FilletStrip) &&
           r.type == SurfType::Cylinder && !r.ringCertified && r.maxVertexDev > tau;
}

// F-3 (D-fc-4.12) with D-gid-8.19. An OPEN chain between torus T and a partner
// that holds T's tangent v-iso (torusPlaneContactCircle /
// torusCylinderContactCircle) takes that circle as its edge. A partial torus
// takes it on the partner's hold alone. A closed torus takes it only when T
// is certified, one cap of T is the D-gid-5.15 mesh polyline (the facets its
// τ-preimage left interrupt that rim, so the chains there are open), and every
// chain vertex is within τ of the circle: chainDev is their largest distance
// to it, < 0 when not measured.
inline bool openChainTakesTorusContact(const Region& T, double tau, double chainDev) {
    if (!T.closed360) return true;
    return T.meshPolylineCap && anchorCertified(tau, T.maxVertexDev) && chainDev >= 0.0 &&
           chainDev <= tau;
}

// D-fc-4.3. A dissolved region leaves the accepted list (the sphere-corner
// stage's compaction, refit_relational_sphere.cpp). Indices into it -- a
// chamfer cone's or a torus blend's neighbour, a pending torus's absorbed
// regions -- follow the regions that stay; a dissolved slot clears the link.
inline void compactDissolved(SegmentWork& work, const std::vector<int>& dissolved) {
    if (dissolved.empty()) return;
    const int nAcc = (int)work.accepted.size();
    std::vector<char> gone((std::size_t)nAcc, 0);
    for (int ai : dissolved)
        if (ai >= 0 && ai < nAcc) gone[(std::size_t)ai] = 1;
    std::vector<int> newIdx((std::size_t)nAcc, -1);
    std::vector<Region> keep;
    keep.reserve((std::size_t)nAcc);
    for (int i = 0; i < nAcc; ++i) {
        if (gone[(std::size_t)i]) continue;
        newIdx[(std::size_t)i] = (int)keep.size();
        keep.push_back(std::move(work.accepted[(std::size_t)i]));
    }
    for (Region& rg : keep) {
        const bool nbr = rg.origin == Origin::ChamferCone || rg.origin == Origin::TorusBlend;
        if (!nbr || rg.filletNbrA < 0 || rg.filletNbrA >= nAcc) continue;
        rg.filletNbrA = newIdx[(std::size_t)rg.filletNbrA];
    }
    for (TorusPendingReclaim& pack : work.torusPending)
        for (int& ai : pack.reclaimAccIdx)
            if (ai >= 0 && ai < nAcc) ai = newIdx[(std::size_t)ai];
    work.accepted.swap(keep);
}

// D-fc-4.3 with D-fc-1.6.10. A prior strip the band robs dissolves at its commit,
// and a strip an earlier band dissolved holds no triangle: neither is a blend
// partner at a corner of the band. robs is sorted.
inline bool blendPartnerGone(const std::vector<int>& robs, int ai, const Region& acc) {
    return acc.tris.empty() || std::binary_search(robs.begin(), robs.end(), ai);
}

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
