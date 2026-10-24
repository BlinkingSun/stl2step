// Same-surface re-host (D-fc-7.3 RQ-HUB-5). A short arc within tau of a
// certified coaxial member takes that member's surface; one vertex outside
// tau keeps the region's own fit. A cylinder already on that surface is not
// written a second time. No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_coax_rehost.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kTau = 2e-5;
static constexpr double kRMember = 25.0;
static constexpr double kROwn = 25.151;

struct OwnFit {
    double radius = kROwn;
    gp_XYZ axis = gp_XYZ(0, 0, 1);
};

static std::vector<gp_XYZ> shortArc() {
    std::vector<gp_XYZ> v;
    for (double z : {0.0, 1.0}) {
        for (double a : {-0.02, 0.0, 0.02})
            v.emplace_back(kRMember * std::cos(a), kRMember * std::sin(a), z);
    }
    return v;
}

static CoaxRehostMember memberAt(int acc, double radius, bool certified) {
    CoaxRehostMember m;
    m.acc = acc;
    m.origin = gp_XYZ(0, 0, 0);
    m.axis = gp_XYZ(0, 0, 1);
    m.radius = radius;
    m.certified = certified;
    return m;
}

static void apply(OwnFit& own, const CoaxRehostDecision& d) {
    if (!d.take) return;
    own.radius = d.radius;
    const double am = d.axis.Modulus();
    own.axis = am > 0.0 ? d.axis / am : d.axis;
}

static void caseTakesMember() {
    OwnFit own;
    const std::vector<gp_XYZ> verts = shortArc();
    const CoaxRehostDecision d =
        coaxRehostChoose(verts, {memberAt(3, kRMember, true)}, kTau);
    apply(own, d);
    expect(d.take, "short arc within tau certifies on the member");
    expect(d.acc == 3 && std::fabs(own.radius - kRMember) <= 0.0,
           "short arc takes the member radius");
    expect(std::fabs(own.axis.Z() - 1.0) <= 1e-15 && std::fabs(own.axis.X()) <= 1e-15,
           "short arc takes the member axis");
}

static void caseOneVertexOutside() {
    OwnFit own;
    std::vector<gp_XYZ> verts = shortArc();
    verts[0] += gp_XYZ(3.0 * kTau, 0, 0);
    const CoaxRehostDecision d =
        coaxRehostChoose(verts, {memberAt(3, kRMember, true)}, kTau);
    apply(own, d);
    expect(!d.take, "one vertex outside tau does not certify");
    expect(std::fabs(own.radius - kROwn) <= 0.0, "one vertex outside tau keeps its own fit");
}

static void caseNotCertified() {
    OwnFit own;
    const CoaxRehostDecision d =
        coaxRehostChoose(shortArc(), {memberAt(1, kRMember, false)}, kTau);
    apply(own, d);
    expect(!d.take && std::fabs(own.radius - kROwn) <= 0.0,
           "an uncertified member is not a re-host target");
}

static void caseLowestIndex() {
    const CoaxRehostDecision d = coaxRehostChoose(
        shortArc(), {memberAt(5, kRMember, true), memberAt(2, kRMember, true)}, kTau);
    expect(d.take && d.acc == 2, "two certifying members: the lower accepted index");
}

static void caseMutualKeepsLower() {
    const double tau = 1e-4;
    const double rLo = 10.0;
    const double rHi = rLo + 0.5 * tau;
    std::vector<gp_XYZ> onLo, onHi;
    for (double z : {0.0, 1.0}) {
        for (double a : {0.0, 0.4, 0.8, 1.2}) {
            onLo.emplace_back(rLo * std::cos(a), rLo * std::sin(a), z);
            onHi.emplace_back(rHi * std::cos(a), rHi * std::sin(a), z);
        }
    }
    CoaxRehostMember lo = memberAt(1, rLo, true);
    CoaxRehostMember hi = memberAt(4, rHi, true);
    lo.radius = rLo;
    hi.radius = rHi;
    const std::vector<CoaxRehostMember> members{lo, hi};
    expect(coaxRehostStable(onLo, lo.acc, members, tau),
           "the lower member of a mutual pair keeps its own fit");
    expect(!coaxRehostStable(onHi, hi.acc, members, tau),
           "the higher member is not a second surface");
    OwnFit own;
    own.radius = rHi;
    const CoaxRehostDecision d = coaxRehostChoose(onHi, {lo}, tau);
    apply(own, d);
    expect(d.take && std::fabs(own.radius - rLo) <= 0.0,
           "the higher member takes the lower member's surface");
}

static void caseRefusesEmpty() {
    expect(!coaxSameSurfaceRehost({}, gp_XYZ(0, 0, 0), gp_XYZ(0, 0, 1), kRMember, kTau),
           "an empty vertex set does not re-host");
    expect(!coaxRehostChoose(shortArc(), {memberAt(1, kRMember, true)}, 0.0).take,
           "tau <= 0 admits nothing");
}

// Caller shape: if (coaxRehostAlready(...)) continue; else write the member.
// The stored origin may sit anywhere on the member axis (placement is the
// region's centroid). That slide is still the same surface.
static void caseAlreadySkipsSecondWrite() {
    const double length = 1.0;
    const CoaxRehostDecision d =
        coaxRehostChoose(shortArc(), {memberAt(3, kRMember, true)}, kTau);
    expect(d.take, "already-guard: the short arc certifies on the member");

    gp_XYZ origin(0, 0, 0.4);
    gp_XYZ storedAxis(0, 0, 1);
    double storedR = kROwn;
    const bool first = coaxRehostAlready(origin, storedAxis, storedR, length, d.origin, d.axis,
                                         d.radius, kTau);
    expect(!first, "own radius is not already the member, so the first write runs");
    if (!first) {
        storedR = d.radius;
        const double am = d.axis.Modulus();
        storedAxis = am > 0.0 ? d.axis / am : d.axis;
    }
    expect(std::fabs(storedR - kRMember) <= 0.0, "the first write stores the member radius");

    const bool second = coaxRehostAlready(origin, storedAxis, storedR, length, d.origin, d.axis,
                                          d.radius, kTau);
    expect(second, "an axial slide of the written cylinder is already the member");
    const double rFreeze = storedR;
    const double zFreeze = storedAxis.Z();
    const double oFreeze = origin.Z();
    if (!second) {
        storedR = d.radius;
        storedAxis = d.axis;
        origin = d.origin;
    }
    expect(std::fabs(storedR - rFreeze) <= 0.0 && std::fabs(storedAxis.Z() - zFreeze) <= 0.0 &&
               std::fabs(origin.Z() - oFreeze) <= 0.0,
           "a second pass does not write the region again");
}

static void caseAlreadyBothBranches() {
    const double length = 1.0;
    const gp_XYZ oM(0, 0, 0);
    const gp_XYZ aM(0, 0, 1);
    const gp_XYZ oSlide(0, 0, 0.4);
    const double lim = tanAngTol(kTau, kTau, length);
    const auto tilt = [](double cross) {
        return gp_XYZ(cross, 0, std::sqrt(1.0 - cross * cross));
    };

    expect(coaxRehostAlready(oSlide, aM, kRMember, length, oM, aM, kRMember, kTau),
           "parallel, same radius, origin slid along the axis: already");
    expect(coaxRehostAlready(oSlide, gp_XYZ(0, 0, 2), kRMember, length, oM, gp_XYZ(0, 0, 3),
                             kRMember, kTau),
           "non-unit axes on the same line still count as already");
    expect(coaxRehostAlready(gp_XYZ(0.5 * kTau, 0, 0.4), aM, kRMember + 0.5 * kTau, length, oM,
                             aM, kRMember, kTau),
           "perpendicular gap and radius inside tau: already");
    expect(coaxRehostAlready(gp_XYZ(kTau, 0, 0), aM, kRMember + kTau, length, oM, aM, kRMember,
                             kTau),
           "perpendicular gap and radius equal to tau: already");
    expect(coaxRehostAlready(oSlide, tilt(0.5 * lim), kRMember, length, oM, aM, kRMember, kTau),
           "tilt inside tanAngTol: already");
    expect(coaxRehostAlready(oM, aM, kRMember, 0.0, oM, aM, kRMember, kTau),
           "a zero length floors at tau and a matching cylinder is already");

    expect(!coaxRehostAlready(oSlide, aM, kROwn, length, oM, aM, kRMember, kTau),
           "the region's own radius is not already the member");
    expect(!coaxRehostAlready(oSlide, aM, kRMember + 2.0 * kTau, length, oM, aM, kRMember, kTau),
           "radius past tau is not already");
    expect(!coaxRehostAlready(gp_XYZ(2.0 * kTau, 0, 0), aM, kRMember, length, oM, aM, kRMember,
                              kTau),
           "perpendicular gap past tau is not already");
    expect(!coaxRehostAlready(oSlide, tilt(2.0 * lim), kRMember, length, oM, aM, kRMember, kTau),
           "tilt past tanAngTol is not already");
    expect(!coaxRehostAlready(oSlide, aM, kRMember, length, oM, aM, kRMember, 0.0),
           "tau <= 0 is not already");
    expect(!coaxRehostAlready(oSlide, gp_XYZ(0, 0, 0), kRMember, length, oM, aM, kRMember, kTau),
           "a zero region axis is not already");
    expect(!coaxRehostAlready(oSlide, aM, kRMember, length, oM, gp_XYZ(0, 0, 0), kRMember, kTau),
           "a zero member axis is not already");
}

int main() {
    caseTakesMember();
    caseOneVertexOutside();
    caseNotCertified();
    caseLowestIndex();
    caseMutualKeepsLower();
    caseRefusesEmpty();
    caseAlreadySkipsSecondWrite();
    caseAlreadyBothBranches();
    if (gFail) {
        std::fprintf(stderr, "%d FAILED\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "all passed\n");
    return 0;
}
