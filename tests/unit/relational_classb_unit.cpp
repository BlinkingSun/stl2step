// Class-B heal predicates (D-fc-6.7); ring-certified walls spared (D-gid-7.13).
// SPDX-License-Identifier: MIT

#include "refit.hpp"
#include "refit_classb.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

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
    expect(classBDissolveQualifies(true, maxDev, false, tau, swept, dVol, area),
           "all conjuncts: qualifies");
    expect(!classBDissolveQualifies(false, maxDev, false, tau, swept, dVol, area),
           "not shipped: refused");
    expect(!classBDissolveQualifies(true, tau * 0.5, false, tau, swept, dVol, area),
           "certified surface: refused");
    expect(!classBDissolveQualifies(true, maxDev, false, tau, 1.0, dVol, area),
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
    expect(classBCertFailPassMember(true, tau * 2.0, false, tau, 0),
           "shipped over-tau region is a member when nQual is 0");
    expect(!classBCertFailPassMember(true, tau * 2.0, false, tau, 3),
           "swept candidates suppress the cert-fail pass");
    expect(!classBCertFailPassMember(false, tau * 2.0, false, tau, 0),
           "an unshipped region is not a member");
    expect(!classBCertFailPassMember(true, tau * 0.5, false, tau, 0),
           "a certified region stays");
    expect(!classBCertFailPassMember(true, tau, false, tau, 0),
           "a residual equal to tau stays");
}

// D-gid-7.13. A closed hole wall faceted as a regular n-gon (corners on the
// circle) carries the chord sagitta R(1 - cos(pi/n)) as its cylinder residual.
// When closeRingWalls certified it (ringCertified), neither class-B pass counts
// it as failing; the same row without the certificate fails as before.
static void testRingCertifiedSpared() {
    const double kPi = 3.14159265358979323846;
    const double tau = 2e-5;
    const int n = 12;
    const double R = 1.0, H = 2.0;
    const double sag = R * (1.0 - std::cos(kPi / n));
    expect(sag > 100.0 * tau, "d-gid-7.13: the chord sag of the wall is many tau");
    const double area = 2.0 * n * R * std::sin(kPi / n) * H;  // the n facets
    // Swept volume between facets and cylinder: n circular segments times H.
    const double swept = 0.5 * R * R * (2.0 * kPi - n * std::sin(2.0 * kPi / n)) * H;
    const double dVol = 0.0;
    expect(swept > classBVolumeBound(dVol, area, tau),
           "d-gid-7.13: its swept volume exceeds the class-B bound");

    expect(classBOwnCertFails(sag, false, tau), "d-gid-7.13: uncertified wall fails own cert");
    expect(!classBOwnCertFails(sag, true, tau), "d-gid-7.13: ring-certified wall does not");
    expect(classBDissolveQualifies(true, sag, false, tau, swept, dVol, area),
           "d-gid-7.13: swept pass dissolves the uncertified wall");
    expect(!classBDissolveQualifies(true, sag, true, tau, swept, dVol, area),
           "d-gid-7.13: swept pass spares the ring-certified wall");
    expect(classBCertFailPassMember(true, sag, false, tau, 0),
           "d-gid-7.13: cert-fail pass takes the uncertified wall");
    expect(!classBCertFailPassMember(true, sag, true, tau, 0),
           "d-gid-7.13: cert-fail pass spares the ring-certified wall");

    // The two passes as the driver runs them, over two rows copied from their
    // regions the way measureClassBShipMetrics copies them: the certified wall
    // and the same wall without the certificate.
    Region wall;
    wall.id = 7;
    wall.type = SurfType::Cylinder;
    wall.origin = Origin::CylGrow;
    wall.closed360 = true;
    wall.radius = R;
    wall.maxVertexDev = sag;
    wall.ringCertified = true;
    Region bare = wall;
    bare.id = 8;
    bare.ringCertified = false;
    ClassBRegionMetric rows[2];
    const Region* regs[2] = {&wall, &bare};
    for (int i = 0; i < 2; ++i) {
        rows[i].rid = regs[i]->id;
        rows[i].maxVertexDev = regs[i]->maxVertexDev;
        rows[i].ringCertified = regs[i]->ringCertified;
        rows[i].dVolPredicted = dVol;
        rows[i].shippedAnalytic = true;
        rows[i].area = area;
        rows[i].swept = swept;
    }
    auto qualifying = [&](std::vector<int>& ids) {
        ids.clear();
        for (const auto& w : rows)
            if (classBDissolveQualifies(w.shippedAnalytic, w.maxVertexDev, w.ringCertified, tau,
                                        w.swept, w.dVolPredicted, w.area))
                ids.push_back(w.rid);
        return static_cast<int>(ids.size());
    };
    std::vector<int> ids;
    expect(qualifying(ids) == 1 && ids[0] == 8,
           "d-gid-7.13: swept pass qualifies only the uncertified row");
    // Swept volumes inside the bound: nobody qualifies, the cert-fail pass runs.
    for (auto& w : rows) w.swept = 0.0;
    const int nQual = qualifying(ids);
    expect(nQual == 0 && classBNoQualCertFailPass(nQual),
           "d-gid-7.13: with swept inside the bound the cert-fail pass is selected");
    ids.clear();
    for (const auto& w : rows)
        if (classBCertFailPassMember(w.shippedAnalytic, w.maxVertexDev, w.ringCertified, tau,
                                     nQual))
            ids.push_back(w.rid);
    expect(ids.size() == 1 && ids[0] == 8,
           "d-gid-7.13: cert-fail pass dissolves only the uncertified row");
    // The certificate excuses the residual only; a ring-certified row that did
    // not ship is not a member either way.
    expect(!classBCertFailPassMember(false, sag, true, tau, 0),
           "d-gid-7.13: an unshipped ring-certified row is not a member");
}

int main() {
    testDissolveConjuncts();
    testStopRule();
    testCertFailPass();
    testRingCertifiedSpared();
    if (gFail) {
        std::fprintf(stderr, "relational_classb_unit: %d FAIL\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_classb_unit: all PASS\n");
    return 0;
}
