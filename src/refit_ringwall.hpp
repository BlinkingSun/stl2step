// D-fc-6.10. Ring-wall predicates (RQ-H1, RQ-H2, RQ-H4).
// A chord vertex is a T-junction on a facet edge: it carries no circle
// information, and it certifies on the chord rather than on the circle.
// No new threshold. tau is the caller's tauOf (2q). The corner count is
// minStations(kParamsCircle), passed in. SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RINGWALL_HPP
#define STL2STEP_REFIT_RINGWALL_HPP

#include <gp_Circ.hxx>
#include <gp_XYZ.hxx>

#include <cmath>
#include <vector>

namespace stl2step {
namespace refit {

// Perpendicular distance from P to the infinite line through A and B.
// Negative when the line is degenerate (A == B): that vertex is not a chord.
inline double chordLineDistance(const gp_XYZ& P, const gp_XYZ& A, const gp_XYZ& B) {
    const gp_XYZ u = B - A;
    const double L = u.Modulus();
    if (!(L > 0.0)) return -1.0;
    return (P - A).Crossed(u).Modulus() / L;
}

// RQ-H1, C-1 (D-fc-6.12 F1). A ring vertex is a chord vertex when its
// distance to the segment joining its two cycle neighbours is <= tau.
// The foot parameter t on that chord must lie in [0, 1]; a foot past the
// segment is a fold-back, not a chord. A degenerate neighbour pair is not
// a chord. tau <= 0 admits nothing. t in [0, 1] is the segment parameter,
// not a new tolerance.
inline bool isChordVertex(const gp_XYZ& P, const gp_XYZ& A, const gp_XYZ& B, double tau) {
    if (!(tau > 0.0)) return false;
    const gp_XYZ u = B - A;
    const double L = u.Modulus();
    if (!(L > 0.0)) return false;
    const double t = (P - A).Dot(u) / (L * L);
    if (t < 0.0 || t > 1.0) return false;
    const double d = chordLineDistance(P, A, B);
    return d >= 0.0 && d <= tau;
}

// How a chart vertex certifies against a ring wall (D-fc-6.10).
// Corner: on a certified ring, not a chord vertex.
// Chord: within tau of its neighbours' chord (RQ-H1); its cylinder residual
//         is the facet sagitta and is not the accept test.
// OnCylinder: an interior vertex within tau of the joint fit.
// Off: none of those.
enum class RingWallVert : unsigned char { Corner, Chord, OnCylinder, Off };

inline RingWallVert classifyRingWallVertex(bool onRingCorner, bool onRingChord, double cylDev,
                                           double tau) {
    if (onRingCorner) return RingWallVert::Corner;
    if (onRingChord) return RingWallVert::Chord;
    if (tau > 0.0 && cylDev <= tau) return RingWallVert::OnCylinder;
    return RingWallVert::Off;
}

inline bool ringWallVertexCertifies(RingWallVert v) { return v != RingWallVert::Off; }

// A triangle of a certified ring wall is covered when every vertex certifies.
inline bool ringWallTriangleCovered(RingWallVert a, RingWallVert b, RingWallVert c) {
    return ringWallVertexCertifies(a) && ringWallVertexCertifies(b) && ringWallVertexCertifies(c);
}

// RQ-H2. One vertex on the circle, the other two certified chord vertices.
// Distinct from the narrow-heal cover (kCoverOnVertices = 2), which stands.
inline bool oneVertexChordCover(RingWallVert a, RingWallVert b, RingWallVert c) {
    const RingWallVert v[3] = {a, b, c};
    int corners = 0, chords = 0;
    for (int i = 0; i < 3; ++i) {
        if (v[i] == RingWallVert::Corner) ++corners;
        else if (v[i] == RingWallVert::Chord) ++chords;
        else return false;
    }
    return corners == 1 && chords == 2;
}

// 3-D distance from a point to a circle: radial and axial offset, hypot.
inline double circlePointResidual(const gp_Circ& c, const gp_XYZ& p) {
    const gp_XYZ d = p - c.Location().XYZ();
    const gp_XYZ N = c.Axis().Direction().XYZ();
    const double h = d.Dot(N);
    const double rho = (d - N * h).Modulus();
    return std::hypot(rho - c.Radius(), h);
}

struct RimChordCert {
    bool ok = false;
    int corners = 0;
    int chords = 0;
};

// RQ-H4, C-2. A plane|cylinder chain (the caller has already required a
// closed360 cylinder and an IntAna circle) is a certified contact when every
// vertex is a corner (within tau of the circle) or a chord vertex (within tau
// of the segment between its nearest corners along the chain, foot on that
// segment), and the corners number at least minCorners
// (minStations(kParamsCircle)). An open chain must end on corners. The sewTol
// vertex test does not apply to a certified chain.
inline RimChordCert rimChordCertified(const gp_Circ& c, const std::vector<gp_XYZ>& P,
                                      bool closedLoop, double tau, int minCorners) {
    RimChordCert out;
    const int n = (int)P.size();
    if (!(tau > 0.0) || minCorners < 1 || n < minCorners) return out;
    std::vector<char> on((std::size_t)n, 0);
    for (int i = 0; i < n; ++i) {
        if (circlePointResidual(c, P[(std::size_t)i]) <= tau) {
            on[(std::size_t)i] = 1;
            ++out.corners;
        }
    }
    if (out.corners < minCorners) return out;
    if (!closedLoop && (!on.front() || !on.back())) return out;
    for (int i = 0; i < n; ++i) {
        if (on[(std::size_t)i]) continue;
        int a = i, b = i;
        for (int k = 1; k < n; ++k) {
            const int j = closedLoop ? (i - k + n) % n : i - k;
            if (j < 0) return out;
            if (on[(std::size_t)j]) {
                a = j;
                break;
            }
        }
        for (int k = 1; k < n; ++k) {
            const int j = closedLoop ? (i + k) % n : i + k;
            if (j >= n) return out;
            if (on[(std::size_t)j]) {
                b = j;
                break;
            }
        }
        if (a == i || b == i || a == b) return out;
        const gp_XYZ u = P[(std::size_t)b] - P[(std::size_t)a];
        const double L = u.Modulus();
        if (!(L > 0.0)) return out;
        const double t = (P[(std::size_t)i] - P[(std::size_t)a]).Dot(u) / (L * L);
        if (t < 0.0 || t > 1.0) return out;
        if ((P[(std::size_t)i] - P[(std::size_t)a]).Crossed(u).Modulus() / L > tau) return out;
        ++out.chords;
    }
    out.ok = true;
    return out;
}

// C-3. The seam vertex is the lowest id among ring-0 corners that has a
// ring-1 corner on the same generator (|delta azimuth| * R <= tau). A chord
// vertex is never the seam. -1 when no pair exists; the caller keeps the
// axis frame it already had.
inline int cornerPairSeamVertex(const std::vector<int>& id0, const std::vector<char>& corner0,
                                const std::vector<double>& az0, const std::vector<int>& id1,
                                const std::vector<char>& corner1, const std::vector<double>& az1,
                                double tau, double radius) {
    if (!(tau > 0.0) || !(radius > 0.0)) return -1;
    if (id0.size() != corner0.size() || id0.size() != az0.size()) return -1;
    if (id1.size() != corner1.size() || id1.size() != az1.size()) return -1;
    auto wrap = [](double t) {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kTwoPi = 2.0 * kPi;
        while (t <= -kPi) t += kTwoPi;
        while (t > kPi) t -= kTwoPi;
        return t;
    };
    int pick = -1;
    for (std::size_t i = 0; i < id0.size(); ++i) {
        if (!corner0[i]) continue;
        const int lv = id0[i];
        if (pick >= 0 && lv >= pick) continue;
        for (std::size_t j = 0; j < id1.size(); ++j) {
            if (!corner1[j]) continue;
            if (std::abs(wrap(az1[j] - az0[i])) * radius <= tau) {
                pick = lv;
                break;
            }
        }
    }
    return pick;
}

}  // namespace refit
}  // namespace stl2step

#endif
