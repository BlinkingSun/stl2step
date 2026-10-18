// cyl|cyl contact predicates (D-fc-5.17). No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_cylcyl_contact.hpp"

#include <cmath>
#include <cstdio>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static gp_Ax3 axZ() { return gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)); }
static gp_Ax3 axX() { return gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0), gp_Dir(0, 0, 1)); }

static gp_Pnt onElips(const gp_Elips& e, double t) {
    const gp_Vec x(e.XAxis().Direction());
    const gp_Vec y(e.YAxis().Direction());
    return e.Location().Translated(x * (e.MajorRadius() * std::cos(t)) +
                                   y * (e.MinorRadius() * std::sin(t)));
}

static void caseSeparation() {
    const AxisSep meet = cylinderAxisSeparation(axZ(), axX());
    expect(!meet.parallel && meet.distance <= Precision::Confusion(),
           "perpendicular axes through one point meet");

    gp_Ax3 shifted = axX();
    shifted.SetLocation(gp_Pnt(0, 0.01, 0));
    const AxisSep skew = cylinderAxisSeparation(axZ(), shifted);
    expect(!skew.parallel && std::fabs(skew.distance - 0.01) <= 1e-12,
           "perpendicular axes 0.01 apart");

    gp_Ax3 parallel = axZ();
    parallel.SetLocation(gp_Pnt(1, 0, 0));
    const AxisSep par = cylinderAxisSeparation(axZ(), parallel);
    expect(par.parallel && std::fabs(par.distance - 1.0) <= 1e-12, "translated axis is parallel");
}

static void caseAdmissible() {
    const double tau = 1e-4;
    const AxisSep meet = cylinderAxisSeparation(axZ(), axX());
    expect(steinmetzAdmissible(1.0, 1.0, meet, tau), "equal meeting axes within tau");
    expect(!steinmetzAdmissible(1.0, 1.0, meet, 0.0), "tau <= 0 refuses");
    expect(steinmetzAdmissible(1.0, 1.0 + tau, meet, tau), "radius gap equal to tau admits");
    expect(!steinmetzAdmissible(1.0, 1.0 + tau + tau, meet, tau), "radius gap of 2 tau refuses");
    expect(steinmetzAdmissible(1.0, 1.0 + 0.5 * tau, meet, tau), "half radius gap admits");

    gp_Ax3 shifted = axX();
    shifted.SetLocation(gp_Pnt(0, 0.25 * tau, 0));
    const AxisSep sep = cylinderAxisSeparation(axZ(), shifted);
    expect(steinmetzAdmissible(1.0, 1.0 + 0.25 * tau, sep, tau),
           "radius gap 0.25 tau plus separation 0.25 tau admits");
    expect(!steinmetzAdmissible(1.0, 1.0, sep, 0.2 * tau), "separation above tau refuses");

    gp_Ax3 parallel = axZ();
    parallel.SetLocation(gp_Pnt(0.1 * tau, 0, 0));
    const AxisSep par = cylinderAxisSeparation(axZ(), parallel);
    expect(!steinmetzAdmissible(1.0, 1.0, par, tau), "parallel axes are not Steinmetz");
    expect(steinmetzSections(CylAxis{axZ(), 1.0}, CylAxis{parallel, 1.0}, tau).n == 0,
           "parallel pair builds no section");
}

static void caseExactMeeting() {
    const double tau = 1e-6;
    const double R = 1.0;
    const SteinmetzPair sp = steinmetzSections(CylAxis{axZ(), R}, CylAxis{axX(), R}, tau);
    expect(sp.n == 2, "meeting equal cylinders give two sections");
    if (sp.n != 2) return;
    const gp_Elips els[2] = {sp.e0.elips, sp.e1.elips};
    const CylAxis a{axZ(), R};
    const CylAxis b{axX(), R};
    bool on = true;
    bool shape = true;
    for (const gp_Elips& e : els) {
        if (!(std::fabs(e.MinorRadius() - R) <= 1e-12)) shape = false;
        if (!(std::fabs(e.MajorRadius() - R * std::sqrt(2.0)) <= 1e-9)) shape = false;
        for (int i = 0; i < 32; i++) {
            const gp_Pnt p = onElips(e, (2.0 * cylcylPi() * i) / 32.0);
            if (distToCylinder(a.ax, a.radius, p) > 1e-9) on = false;
            if (distToCylinder(b.ax, b.radius, p) > 1e-9) on = false;
            if (!withinTauOfBoth(p, a, b, tau)) on = false;
        }
    }
    expect(shape, "both sections have minor R and major R*sqrt(2)");
    expect(on, "section samples lie on both cylinders");
    const gp_Pnt far(10, 10, 10);
    expect(!withinTauOfBoth(far, a, b, tau), "a far point is outside tau of both");
}

static void caseChordBound() {
    const double R = 5.0;
    const double tau = 1e-4;
    const gp_Ax3 ax = axZ();
    const double h = polylineHalfAzimuth(R, tau);
    expect(h > 0.0 && h < 0.5 * cylcylPi(), "half-azimuth exists inside the budget");
    const gp_Pnt p0(R, 0, 0);
    // Full azimuth h means the binder's half-angle is h/2, inside the solved step.
    const double a = h;
    const gp_Pnt p1(R * std::cos(a), R * std::sin(a), 0);
    expect(polylineOnCylSup(p0, p1, ax, R) <= tau, "chord inside the solved step is within tau");
    const double coarse = std::min(0.5 * cylcylPi() - Precision::Angular(), h * 16.0);
    const gp_Pnt p2(R * std::cos(coarse), R * std::sin(coarse), 0);
    expect(coarse > a && polylineOnCylSup(p0, p2, ax, R) > tau,
           "a coarse chord exceeds tau");
    expect(polylineHalfAzimuth(R, 0.0) == 0.0, "no budget, no step");
    const gp_Pnt pushed(R + 2.0 * tau, 0, 0);
    expect(distToCylinder(ax, R, pushed) > tau, "a point pushed off the cylinder exceeds tau");
}

static void caseContactChordBand() {
    const gp_Ax3 a = axZ();
    const gp_Ax3 b = axX();
    const double R = 5.0;
    const gp_Pnt p0(R, 0, 0);
    const double h = 0.05;
    const gp_Pnt p1(R * std::cos(h), R * std::sin(h), 0);
    const gp_Pnt fine[2] = {p0, p1};
    const double band = contactChordBand(fine, 2, false, a, R, b, R);
    const double sum = polylineOnCylSup(p0, p1, a, R) + polylineOnCylSup(p0, p1, b, R);
    expect(std::fabs(band - sum) <= 1e-12, "band is the sum of the two chord certificates");
    expect(band > 0.0, "a real chord has a positive band");
    const double h2 = 0.4;
    const gp_Pnt p2(R * std::cos(h2), R * std::sin(h2), 0);
    const gp_Pnt coarse[2] = {p0, p2};
    expect(contactChordBand(coarse, 2, false, a, R, b, R) > band, "a coarser chord widens the band");
    expect(contactChordBand(fine, 1, true, a, R, b, R) == 0.0, "fewer than two points has no band");
    const gp_Pnt tri[3] = {p0, p1, p2};
    const double openB = contactChordBand(tri, 3, false, a, R, b, R);
    const double closedB = contactChordBand(tri, 3, true, a, R, b, R);
    const double closing = polylineOnCylSup(p2, p0, a, R) + polylineOnCylSup(p2, p0, b, R);
    expect(closedB + 1e-15 >= openB, "closing the loop cannot shrink the band");
    expect(closedB + 1e-12 >= closing, "the closing edge is in the closed band");
}

static void caseIntssGate() {
    const double tau = 1e-4;
    const double q = intssTolOf(tau);
    expect(q == tau * 0.5, "IntSS tolerance is tau/2");
    expect(intssTolOf(0.0) == 0.0, "non-positive tau has no IntSS tolerance");
    expect(intssTolOf(-tau) == 0.0, "negative tau has no IntSS tolerance");

    const AxisSep meet = cylinderAxisSeparation(axZ(), axX());
    expect(axesMeetForContact(8.0, 5.0, meet, tau), "unequal perpendicular axes meet");
    expect(!axesMeetForContact(8.0, 5.0, meet, 0.0), "tau <= 0 axes do not meet");

    gp_Ax3 parallel = axZ();
    parallel.SetLocation(gp_Pnt(1.0, 0.0, 0.0));
    const AxisSep par = cylinderAxisSeparation(axZ(), parallel);
    expect(!axesMeetForContact(5.0, 5.0, par, tau), "parallel axes do not meet");

    const double sum = 8.0 + 5.0 + tau;
    gp_Ax3 atSum = axX();
    atSum.SetLocation(gp_Pnt(0.0, sum, 0.0));
    const AxisSep sumSep = cylinderAxisSeparation(axZ(), atSum);
    expect(axesMeetForContact(8.0, 5.0, sumSep, tau), "gap equal to rA+rB+tau meets");
    gp_Ax3 past = axX();
    past.SetLocation(gp_Pnt(0.0, sum + tau, 0.0));
    const AxisSep pastSep = cylinderAxisSeparation(axZ(), past);
    expect(!axesMeetForContact(8.0, 5.0, pastSep, tau), "one more tau of gap does not meet");

    const double band = 0.034;
    const double snap = 2.3;
    const double clearance = 0.017;
    expect(intssAdmissible(8.0, 5.0, meet, tau, clearance, band, snap, -1.0, q),
           "unequal meeting axes inside the band admit IntSS");
    expect(intssAdmissible(8.0, 5.0, meet, tau, band, band, snap, -1.0, q),
           "clearance equal to the band admits IntSS");
    expect(!intssAdmissible(8.0, 5.0, meet, tau, band + tau, band, snap, -1.0, q),
           "clearance above the band refuses IntSS");
    expect(!intssAdmissible(8.0, 5.0, meet, tau, snap + tau, band + snap, snap, -1.0, q),
           "clearance above the snap ceiling refuses IntSS");
    expect(!intssAdmissible(8.0, 5.0, par, tau, clearance, band, snap, -1.0, q),
           "parallel axes do not admit IntSS");
    expect(!intssAdmissible(8.0, 5.0, pastSep, tau, clearance, band, snap, -1.0, q),
           "axes that do not meet do not admit IntSS");

    expect(steinmetzAdmissible(5.0, 5.0, meet, tau), "equal meeting radii are Steinmetz");
    expect(!intssAdmissible(5.0, 5.0, meet, tau, clearance, band, snap, -1.0, q),
           "equal radii with no section do not call IntSS");
    const double gate = band;
    expect(!intssAdmissible(5.0, 5.0, meet, tau, clearance, band, snap, gate + q + q, q),
           "a section that misses the gate by more than q refuses IntSS");
    expect(intssAdmissible(5.0, 5.0, meet, tau, clearance, band, snap, gate + q, q),
           "a section within q of the gate still admits IntSS");

    const double tightSnap = 0.02;
    const double wideBand = 5.0;
    expect(intssAdmissible(8.0, 5.0, meet, tau, 0.01, wideBand, tightSnap, -1.0, q),
           "clearance inside the snap ceiling admits IntSS when the band is wider");
    expect(!intssAdmissible(8.0, 5.0, meet, tau, tightSnap + tau, wideBand, tightSnap, -1.0, q),
           "clearance above a tighter snap ceiling refuses IntSS");

    const SteinmetzPair sp = steinmetzSections(CylAxis{axZ(), 1.0}, CylAxis{axX(), 1.0}, tau);
    expect(sp.n == 2, "section for the clearance sample");
    if (sp.n == 2) {
        const gp_Pnt on = onElips(sp.e0.elips, 0.3);
        const gp_Pnt one[1] = {on};
        const double c = chainSurfaceClearance(one, 1, axZ(), 1.0, axX(), 1.0);
        expect(c >= 0.0 && c <= 1e-9, "a section point has no cylinder clearance");
        const gp_Pnt far[1] = {gp_Pnt(10.0, 10.0, 10.0)};
        expect(chainSurfaceClearance(far, 1, axZ(), 1.0, axX(), 1.0) > tau,
               "a far point's clearance exceeds tau");
    }
    expect(chainSurfaceClearance(nullptr, 1, axZ(), 1.0, axX(), 1.0) == 0.0,
           "no points has no clearance");
}

static void caseSegmentGap() {
    const double tau = 1e-4;
    const gp_Pnt a(0, 0, 0);
    const gp_Pnt b(2, 0, 0);
    expect(segmentGap(a, b, gp_Pnt(1, 0, 0)) <= Precision::Confusion(),
           "midpoint of a segment has no gap");
    expect(segmentGap(a, b, gp_Pnt(1, tau, 0)) <= tau, "offset equal to tau admits");
    expect(segmentGap(a, b, gp_Pnt(1, tau + tau, 0)) > tau, "offset of 2 tau refuses");
    expect(segmentGap(a, a, gp_Pnt(0, 0, 1)) > tau,
           "coincident ends do not hide a point off the chord");
}

int main() {
    caseSeparation();
    caseAdmissible();
    caseExactMeeting();
    caseChordBound();
    caseContactChordBand();
    caseIntssGate();
    caseSegmentGap();
    std::fprintf(stderr, "relational_cylcyl_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
