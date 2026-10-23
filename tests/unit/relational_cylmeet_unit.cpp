// Partial counterbore wall predicates (D-fc-5.17). Public quantities only.
// SPDX-License-Identifier: MIT

#include "refit_cyl_meet.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static void testBudget() {
    expect(kCascadeResidFactor == 3.0, "cascade factor is the existing kR");
    expect(kTightRelFloor == 1e-4, "tight floor is the existing relative term");
    const double b = partialSheetVolBudget(1000.0, 4, 2.0);
    expect(std::fabs(b - 6.0) <= 1e-12, "budget takes the larger of floor and kR*|dVol|");
    const double floorOnly = partialSheetVolBudget(1000.0, 4, 0.0);
    expect(std::fabs(floorOnly - 0.025) <= 1e-12, "zero predictor leaves the relative floor");
    expect(partialSheetVolumeAgrees(10.0, 8.0, 1.5, 6.0), "residual inside the budget agrees");
    // Public-fixture sheet: face volume and chord volume disagree by ~104 mm^3
    // against a predictor of -1.630 and a budget of 4.89.
    const double vFace = -1296.30848, vChord = -1191.85768, dVol = -1.630, budget = 4.89;
    expect(!partialSheetVolumeAgrees(vFace, vChord, dVol, budget),
           "rect-trim sheet outside the budget is refused");
    expect(partialSheetVolumeAgrees(vChord + dVol, vChord, dVol, budget),
           "a face that matches the predictor agrees");
}

static void testSector() {
    const double budget = partialSheetVolBudget(66570.18737, 20, 1.630);
    const double inflated = cylinderSectorGap(10.6442, 10.0, 2.611, 4.5);
    expect(inflated > budget, "stored-versus-certified sector exceeds the sheet budget");
    const double tau = 2e-6;
    const double tiny = cylinderSectorGap(10.0 + tau, 10.0, 2.611, 4.5);
    expect(tiny < budget, "a tau-level radius error stays inside the sheet budget");
    expect(cylinderSectorGap(10.0, 10.0, 2.611, 4.5) == 0.0, "equal radii have no sector gap");
}

static void testAxesAndPlane() {
    expect(axesCross(1.0, 1e-6, 10.0), "perpendicular axes cross");
    expect(!axesCross(0.0, 1e-6, 10.0), "parallel axes do not cross");
    expect(!axesCross(1.0, 1e-6, 1e-9), "a length inside tau is not a crossing");
    expect(planeNormalParallelAxis(1.0, 10.0, 1e-6), "a normal along the axis is a cap plane");
    expect(!planeNormalParallelAxis(0.0, 10.0, 1e-6), "a side plane is not a cap");
}

static void testIntervene() {
    CylMeetQuantities q;
    q.cylinder = true;
    q.closed360 = false;
    q.tau = 1e-6;
    q.storedDev = 0.644;
    q.rStored = 10.6442;
    q.rFit = 10.0;
    q.fitOk = true;
    q.crossingMate = true;
    q.gap = 78.0;
    q.budget = 4.89;
    expect(cylMeetIntervene(q), "over-budget open wall with a crossing mate intervenes");
    q.closed360 = true;
    expect(!cylMeetIntervene(q), "a closed wall does not intervene");
    q.closed360 = false;
    q.storedDev = q.tau;
    expect(!cylMeetIntervene(q), "a surface already within tau does not intervene");
    q.storedDev = 0.644;
    q.gap = q.budget;
    expect(!cylMeetIntervene(q), "a gap that meets the budget does not intervene");
    q.gap = 78.0;
    q.crossingMate = false;
    expect(!cylMeetIntervene(q), "no crossing mate does not intervene");
    q.crossingMate = true;
    q.fitOk = false;
    expect(!cylMeetIntervene(q), "an uncertified section does not intervene");
    q.fitOk = true;
    q.rStored = 9.0;
    q.rFit = 10.0;
    expect(!cylMeetIntervene(q), "a stored radius inside the certificate does not intervene");
}

static void testOuterRho() {
    const double tau = 1e-6;
    std::vector<double> rho, z;
    for (int level = 0; level < 2; ++level) {
        const double zz = level == 0 ? 0.0 : 4.5;
        for (int k = 0; k < 8; ++k) {
            rho.push_back(10.0);
            z.push_back(zz);
        }
    }
    // Bore-cut samples sit inside the outer radius.
    rho.push_back(9.4);
    z.push_back(1.0);
    rho.push_back(8.2);
    z.push_back(2.0);
    RhoCertificate rc;
    expect(certifyOuterRho(rho.data(), z.data(), (int)rho.size(), tau, rc),
           "outer rho of a bore-cut wall certifies");
    expect(rc.ok && std::fabs(rc.R - 10.0) <= tau && rc.inliers == 16 && rc.levels >= 2,
           "the certificate is the maximum rho, not the cut");
    rho[0] = 10.0 + 10.0 * tau;
    RhoCertificate shifted;
    expect(!certifyOuterRho(rho.data(), z.data(), (int)rho.size(), tau, shifted),
           "one sample outside tau does not certify the wall");
    expect(std::fabs(shifted.R - (10.0 + 10.0 * tau)) <= 1e-12 && shifted.inliers == 1 && !shifted.ok,
           "that sample is the outer radius and stands alone");
}

static void testSectorCentre() {
    const double open[] = {-1.0, 0.0, 1.0};
    SectorCentre sc;
    expect(cylMeetSectorCentre(open, 3, sc), "three samples on an open arc have a centre");
    expect(sc.ok && std::fabs(sc.centre) <= 1e-12 && std::fabs(sc.span - 2.0) <= 1e-12,
           "the centre is the midpoint of that arc");
    const double side[] = {1.5, 2.0, 2.5};
    expect(cylMeetSectorCentre(side, 3, sc) && std::fabs(sc.centre - 2.0) <= 1e-12 &&
               std::fabs(sc.span - 1.0) <= 1e-12,
           "a sector away from the angle origin keeps its midpoint");
    const double wrap[] = {-3.0, -2.7, 2.2, 2.6};
    expect(cylMeetSectorCentre(wrap, 4, sc), "a sector across the branch has a centre");
    const double span = 1.3831853071795862;
    const double centre = 2.891592653589793;
    expect(std::fabs(sc.span - span) <= 1e-12 && std::fabs(sc.centre - centre) <= 1e-12,
           "the centre is the midpoint of the wrapped arc");
    expect(!cylMeetSectorCentre(open, 1, sc), "one sample has no sector");
}

// D-fc-6.4(a): the RULE 2.2 shell budget keeps the literals it replaced.
static void testShellBudget() {
    const double meshAbs = 66570.187372, shippedAbs = 1.25;
    expect(shellVolBudget(meshAbs, shippedAbs) == std::max(1e-4 * meshAbs, 3.0 * shippedAbs),
           "shell budget is max(1e-4 * |mesh|, 3 * shipped |dVol|), bit for bit");
    expect(shellVolBudget(meshAbs, 0.0) == 1e-4 * meshAbs, "no predictor leaves the floor");
    expect(shellVolBudget(100.0, 10.0) == 30.0, "a large predictor sets the budget");
    expect(partialSheetVolBudget(meshAbs, 1, shippedAbs) == shellVolBudget(meshAbs, shippedAbs),
           "one analytic region gets the shell budget");
}

// D-fc-6.4(b): the chord census measures the chord midpoint, as segment() does.
static void testChordCensus() {
    const gp_Ax3 zAxis(gp_Pnt(30, 30, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const double R = 10.0, half = 3.14159265358979323846 / 48.0;
    const gp_XYZ p(30 + R * std::cos(-half), 30 + R * std::sin(-half), 15.0);
    const gp_XYZ q(30 + R * std::cos(half), 30 + R * std::sin(half), 15.0);
    const double sag = R * (1.0 - std::cos(half));
    expect(std::fabs(chordMidDevCylinder(zAxis, R, p, q) - sag) <= 1e-12,
           "a wall chord deviates from its cylinder by the sagitta");
    expect(std::fabs(chordMidDevCylinder(zAxis, R + 0.644, p, q) - (0.644 + sag)) <= 1e-12,
           "the same chord against an inflated radius carries the radius error");
    const gp_Ax3 floorPlane(gp_Pnt(0, 0, 14), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    expect(std::fabs(chordMidDevPlane(floorPlane, p, q) - 1.0) <= 1e-12,
           "a chord one unit above a plane deviates by one");
    const gp_XYZ a(0, 0, 14), b(5, 3, 14);
    expect(chordMidDevPlane(floorPlane, a, b) == 0.0, "a chord in the plane does not deviate");
}

int main() {
    testBudget();
    testSector();
    testAxesAndPlane();
    testIntervene();
    testOuterRho();
    testSectorCentre();
    testShellBudget();
    testChordCensus();
    if (gFail) {
        std::fprintf(stderr, "relational_cylmeet_unit FAILED %d\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_cylmeet_unit OK\n");
    return 0;
}
