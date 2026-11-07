// Class-B heal predicates (D-fc-6.7).
// SPDX-License-Identifier: MIT

#include "refit_classb.hpp"

#include <cmath>
#include <cstdio>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static void testDissolveConjuncts() {
    const double tau = 2e-5;
    const double dVol = 8.83;
    const double area = 100.0;
    const double swept = 2000.0;
    const double maxDev = 0.04;
    expect(classBDissolveQualifies(true, maxDev, tau, swept, dVol, area),
           "all conjuncts: qualifies");
    expect(!classBDissolveQualifies(false, maxDev, tau, swept, dVol, area),
           "not shipped: refused");
    expect(!classBDissolveQualifies(true, tau * 0.5, tau, swept, dVol, area),
           "certified surface: refused");
    expect(!classBDissolveQualifies(true, maxDev, tau, 1.0, dVol, area),
           "swept below bound: refused");
    const double bound = classBVolumeBound(dVol, area, tau);
    expect(swept > bound, "S29-scale swept exceeds bound");
}

static void testStopRule() {
    const int ids[] = {2, 3, 4};
    expect(classBHealPassHasCandidates(ids, 3), "non-empty pass has candidates");
    expect(!classBHealPassHasCandidates(ids, 0), "empty pass stops heal");
}

static void testCertFailPass() {
    const double tau = 2e-5;
    expect(classBNoQualCertFailPass(0), "nQual 0 selects the cert-fail pass");
    expect(!classBNoQualCertFailPass(1), "positive nQual stays on the swept loop");
    expect(!classBNoQualCertFailPass(-1), "a negative count is not an empty set");
    expect(classBCertFailPassMember(true, tau * 2.0, tau, 0),
           "shipped over-tau region is a member when nQual is 0");
    expect(!classBCertFailPassMember(true, tau * 2.0, tau, 3),
           "swept candidates suppress the cert-fail pass");
    expect(!classBCertFailPassMember(false, tau * 2.0, tau, 0),
           "an unshipped region is not a member");
    expect(!classBCertFailPassMember(true, tau * 0.5, tau, 0),
           "a certified region stays");
    expect(!classBCertFailPassMember(true, tau, tau, 0),
           "a residual equal to tau stays");
}

int main() {
    testDissolveConjuncts();
    testStopRule();
    testCertFailPass();
    if (gFail) {
        std::fprintf(stderr, "relational_classb_unit: %d FAIL\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_classb_unit: all PASS\n");
    return 0;
}
