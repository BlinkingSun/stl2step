// D-fc-6.10. Chord vertex (RQ-H1), one-vertex cover (RQ-H2), plane|cylinder
// rim certificate (RQ-H4), corner-pair seam (C-3). Synthetic geometry only.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_ringwall.hpp"

#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static gp_XYZ xyz(double x, double y, double z) { return gp_XYZ(x, y, z); }

static void testChordVertex() {
    const double tau = 1e-4;
    const gp_XYZ A = xyz(0, 0, 0), B = xyz(1, 0, 0);
    const gp_XYZ mid = xyz(0.5, 0, 0);
    expect(isChordVertex(mid, A, B, tau), "h1: midpoint of the neighbour chord");
    expect(isChordVertex(xyz(0.5, tau, 0), A, B, tau), "h1: offset of exactly tau is a chord");
    expect(!isChordVertex(xyz(0.5, 2 * tau, 0), A, B, tau), "h1: offset of 2 tau is not a chord");
    expect(!isChordVertex(mid, A, B, 0.0), "h1: tau <= 0 admits nothing");
    expect(!isChordVertex(mid, A, A, tau), "h1: a degenerate neighbour pair is not a chord");
    // D-fc-6.12 F1. Within tau of the infinite line, foot past the segment.
    expect(!isChordVertex(xyz(2.0, 0.5 * tau, 0), A, B, tau),
           "h1: a foot past the segment is not a chord vertex");

    // A true corner of a regular octagon sits at the sagitta of its neighbours'
    // chord, far above tau, so it stays a corner.
    const double R = 1.0;
    const double a = 3.14159265358979323846 / 4.0;
    const gp_XYZ left = xyz(R * std::cos(a), R * std::sin(a), 0);
    const gp_XYZ right = xyz(R * std::cos(-a), R * std::sin(-a), 0);
    const gp_XYZ corner = xyz(R, 0, 0);
    const double sag = chordLineDistance(corner, left, right);
    expect(sag > 1000 * tau, "h1: octagon corner sagitta is many tau");
    expect(!isChordVertex(corner, left, right, tau), "h1: a true corner is not a chord vertex");

    // A T-junction on that same chord, lifted by half tau, is a chord vertex.
    const gp_XYZ junction = xyz(R * std::cos(a), 0, 0.5 * tau);
    expect(isChordVertex(junction, left, right, tau), "h1: T-junction on the corner chord");
}

static void testOneVertexCover() {
    const double tau = 1e-4;
    expect(classifyRingWallVertex(true, false, 1000 * tau, tau) == RingWallVert::Corner,
           "h2: a ring corner stays a corner when the cylinder residual is large");
    expect(classifyRingWallVertex(false, true, 1000 * tau, tau) == RingWallVert::Chord,
           "h2: a chord vertex certifies off the cylinder");
    expect(classifyRingWallVertex(false, false, tau, tau) == RingWallVert::OnCylinder,
           "h2: an interior vertex within tau is on the cylinder");
    expect(classifyRingWallVertex(false, false, tau * 1.01, tau) == RingWallVert::Off,
           "h2: an interior vertex past tau is off");
    expect(classifyRingWallVertex(false, false, 0.0, 0.0) == RingWallVert::Off,
           "h2: tau <= 0 does not certify an interior vertex");

    const RingWallVert one[3] = {RingWallVert::Corner, RingWallVert::Chord, RingWallVert::Chord};
    expect(oneVertexChordCover(one[0], one[1], one[2]), "h2: one corner and two chords");
    expect(oneVertexChordCover(one[1], one[0], one[2]), "h2: the corner may be any vertex");
    expect(oneVertexChordCover(one[1], one[2], one[0]), "h2: the corner may be the last vertex");
    expect(ringWallTriangleCovered(one[0], one[1], one[2]), "h2: that triangle is covered");
    expect(!oneVertexChordCover(RingWallVert::Corner, RingWallVert::Chord, RingWallVert::OnCylinder),
           "h2: a non-chord partner is not the one-vertex rule");
    expect(!oneVertexChordCover(RingWallVert::Corner, RingWallVert::Corner, RingWallVert::Chord),
           "h2: two corners is not the one-vertex rule");
    expect(ringWallTriangleCovered(RingWallVert::Corner, RingWallVert::Corner, RingWallVert::Corner),
           "h2: three corners are covered");
    expect(!ringWallTriangleCovered(RingWallVert::Corner, RingWallVert::Chord, RingWallVert::Off),
           "h2: an off vertex is not covered");
}

static gp_Circ unitCircle() {
    return gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 1.0);
}

static gp_XYZ onCircle(double ang) { return xyz(std::cos(ang), std::sin(ang), 0); }

static void testRimCertificate() {
    const double tau = 1e-4;
    const int minC = minStations(kParamsCircle);
    expect(minC == 4, "h4: a circle needs minStations(kParamsCircle) corners");
    const gp_Circ c = unitCircle();

    std::vector<gp_XYZ> corners;
    for (int i = 0; i < 8; ++i) corners.push_back(onCircle(i * 3.14159265358979323846 / 4.0));
    const RimChordCert clean = rimChordCertified(c, corners, true, tau, minC);
    expect(clean.ok && clean.corners == 8 && clean.chords == 0, "h4: eight corners, no chords");

    std::vector<gp_XYZ> withChord = corners;
    const gp_XYZ chord = xyz(0.5 * (corners[0].X() + corners[1].X()),
                             0.5 * (corners[0].Y() + corners[1].Y()), 0.5 * tau);
    withChord.insert(withChord.begin() + 1, chord);
    const RimChordCert mixed = rimChordCertified(c, withChord, true, tau, minC);
    expect(mixed.ok && mixed.corners == 8 && mixed.chords == 1, "h4: a chord vertex on the rim certifies");

    std::vector<gp_XYZ> offRim = corners;
    offRim.insert(offRim.begin() + 1, xyz(0.5 * (corners[0].X() + corners[1].X()),
                                          0.5 * (corners[0].Y() + corners[1].Y()), 10 * tau));
    expect(!rimChordCertified(c, offRim, true, tau, minC).ok, "h4: a vertex off the chord is refused");

    std::vector<gp_XYZ> few(corners.begin(), corners.begin() + (minC - 1));
    expect(!rimChordCertified(c, few, true, tau, minC).ok, "h4: fewer than minStations corners is refused");

    expect(!rimChordCertified(c, corners, true, 0.0, minC).ok, "h4: tau <= 0 admits nothing");

    // Open chain: the chord vertex sits between corners, and the chain ends on corners.
    std::vector<gp_XYZ> open = {corners[0], chord, corners[1], corners[2], corners[3], corners[4]};
    const RimChordCert op = rimChordCertified(c, open, false, tau, minC);
    expect(op.ok && op.chords == 1, "h4: an open chain ending on corners certifies");

    std::vector<gp_XYZ> badEnds = {chord, corners[0], corners[1], corners[2], corners[3]};
    expect(!rimChordCertified(c, badEnds, false, tau, minC).ok,
           "h4: an open chain that does not end on a corner is refused");

    // Foot past the segment: close to the infinite line, not on the corner segment.
    std::vector<gp_XYZ> past = corners;
    past.insert(past.begin() + 1, xyz(corners[0].X() + 0.5, corners[0].Y(), 0.0));
    expect(!rimChordCertified(c, past, true, tau, minC).ok,
           "h4: a foot outside the corner segment is refused");
}

static void testCornerPairSeam() {
    const double tau = 1e-4;
    const double R = 2.0;
    // Two ring-0 corners share a generator with the ring-1 corner. The lower id
    // wins. A still-lower chord vertex on that same generator is not eligible.
    const std::vector<int> id0 = {10, 4, 2};
    const std::vector<char> c0 = {1, 1, 0};
    const std::vector<double> az0 = {0.0, 0.0, 0.0};
    const std::vector<int> id1 = {3, 8};
    const std::vector<char> c1 = {1, 0};
    const std::vector<double> az1 = {0.5 * tau / R, 1.0};
    expect(cornerPairSeamVertex(id0, c0, az0, id1, c1, az1, tau, R) == 4,
           "c3: lowest corner id with a mate on the other ring");
    expect(cornerPairSeamVertex(id0, c0, az0, id1, c1, az1, tau, R) != 2,
           "c3: a chord vertex is not the seam");
    const std::vector<double> azFar = {1.0, 1.0};
    expect(cornerPairSeamVertex(id0, c0, az0, id1, c1, azFar, tau, R) == -1,
           "c3: no generator pair leaves the frame unchanged");
    expect(cornerPairSeamVertex(id0, c0, az0, id1, c1, az1, 0.0, R) == -1,
           "c3: tau <= 0 pairs nothing");
}

int main() {
    testChordVertex();
    testOneVertexCover();
    testRimCertificate();
    testCornerPairSeam();
    if (gFail) std::fprintf(stderr, "relational_ringwall_unit: %d FAIL\n", gFail);
    return gFail ? 1 : 0;
}
