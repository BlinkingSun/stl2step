// Synthetic anchored quarter-round cases (D-gid-2.6): certifyAnchoredRound on a
// three-row quarter-round between a coaxial wall and a plane normal to the
// axis. Four profile stations certify with two anchors (one free parameter)
// and with one anchor (two free), and refuse as a free round; certifyTorus
// refuses the same stations (profPts < 5). No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_anchored_round.hpp"
#include "refit_profile_cert.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kTau = 2e-6;

// A concave quarter-round: wall rho = Rw below, plane z = zp beside it, tube
// centre (Rw + r, zp - r). Stations at profile angles 180 - 90 * k / rows
// degrees, k = 0..rows (k = 0 on the wall, k = rows on the plane), each
// repeated at `cols` azimuths. `bump` moves station `bumpK` off the circle.
static std::vector<ProfPt> round(double Rw, double zp, double r, int rows, int cols,
                                 double bump = 0.0, int bumpK = -1) {
    std::vector<ProfPt> P;
    const double cr = Rw + r, cz = zp - r;
    for (int k = 0; k <= rows; ++k) {
        const double a = kPi - 0.5 * kPi * (double)k / (double)rows;
        const double rr = r + (k == bumpK ? bump : 0.0);
        for (int j = 0; j < cols; ++j) {
            ProfPt p;
            p.rho = cr + rr * std::cos(a);
            p.z = cz + rr * std::sin(a);
            p.psi = 2.0 * kPi * (double)j / (double)cols;
            P.push_back(p);
        }
    }
    return P;
}

int main() {
    // Generic synthetic dimensions (no part's radii or heights).
    const double Rw = 10.0, zp = 3.0, r = 1.0;
    const ProfileAnchor wall{1.0, 0.0, Rw};
    const ProfileAnchor plane{0.0, 1.0, zp};

    // Three rows: four stations.
    const std::vector<ProfPt> P4 = round(Rw, zp, r, 3, 8);
    {
        const TorusCert free = certifyTorus(P4, kTau);
        expect(!free.ok && std::strcmp(free.why, "profile-underdetermined") == 0,
               "certifyTorus refuses 4 stations (profPts < 5)");
    }
    {
        const AnchoredRoundCert c = certifyAnchoredRound(P4, {wall, plane}, kTau);
        expect(c.ok, "two anchors: 4 stations certify");
        expect(c.params == 1 && c.stations == 4, "two anchors: 1 free parameter, 4 stations");
        expect(std::fabs(c.Rmin - r) <= kTau && std::fabs(c.Rmaj - (Rw + r)) <= kTau &&
                   std::fabs(c.z0 - (zp - r)) <= kTau,
               "two anchors: circle recovered within tau");
        expect(c.maxDev <= kTau, "two anchors: every station within tau");
    }
    {
        const AnchoredRoundCert c = certifyAnchoredRound(P4, {plane}, kTau);
        expect(c.ok && c.params == 2, "one anchor: 4 stations certify (2 free)");
    }
    {
        const AnchoredRoundCert c = certifyAnchoredRound(P4, {}, kTau);
        expect(!c.ok && std::strcmp(c.why, "free-round") == 0,
               "free 4-station round refuses (certifyTorus governs free fits)");
    }
    // One row: two stations. Two anchors leave one parameter: 2 >= minStations(1).
    // One anchor leaves two: 2 < 3.
    {
        const std::vector<ProfPt> P2 = round(Rw, zp, r, 1, 8);
        const AnchoredRoundCert c1 = certifyAnchoredRound(P2, {plane}, kTau);
        expect(!c1.ok && std::strcmp(c1.why, "underdetermined") == 0,
               "one anchor: 2 stations underdetermined");
    }
    // A station 3 tau off the circle refuses.
    {
        const std::vector<ProfPt> Pb = round(Rw, zp, r, 3, 8, 3.0 * kTau, 2);
        const AnchoredRoundCert c = certifyAnchoredRound(Pb, {wall, plane}, kTau);
        expect(!c.ok && c.determinate && std::strcmp(c.why, "vertex-off-surface") == 0,
               "two anchors: a station 3 tau off refuses vertex-off-surface");
    }
    // The wrong corner: the same circle offered a wall on the far side of the
    // tube centre. The round does not meet that anchor.
    {
        const ProfileAnchor farWall{1.0, 0.0, Rw + 2.0 * r};
        const AnchoredRoundCert c = certifyAnchoredRound(P4, {farWall, plane}, kTau);
        expect(!c.ok, "two anchors: a wall the round does not meet refuses");
    }
    // A three-quarter round past the plane's tangent point: stations off the
    // arc between the two tangent points.
    {
        std::vector<ProfPt> P = round(Rw, zp, r, 3, 8);
        const double cr = Rw + r, cz = zp - r;
        for (int j = 0; j < 8; ++j) {
            ProfPt p;
            const double a = -0.25 * kPi;
            p.rho = cr + r * std::cos(a);
            p.z = cz + r * std::sin(a);
            p.psi = 2.0 * kPi * (double)j / 8.0;
            P.push_back(p);
        }
        const AnchoredRoundCert c = certifyAnchoredRound(P, {wall, plane}, kTau);
        expect(!c.ok && std::strcmp(c.why, "not-blend") == 0,
               "two anchors: a station past a tangent point refuses not-blend");
    }
    // A profile line (a 45 deg chamfer through the same corner) is a cone, not
    // a round: line first.
    {
        std::vector<ProfPt> P;
        for (int k = 0; k <= 3; ++k)
            for (int j = 0; j < 8; ++j) {
                ProfPt p;
                p.rho = Rw + r * (double)k / 3.0;
                p.z = zp - r + r * (double)k / 3.0;
                p.psi = 2.0 * kPi * (double)j / 8.0;
                P.push_back(p);
            }
        const AnchoredRoundCert c = certifyAnchoredRound(P, {wall, plane}, kTau);
        expect(!c.ok && std::strcmp(c.why, "line") == 0, "a profile line refuses as line");
    }
    // One azimuth only: the revolution is not evidenced.
    {
        const std::vector<ProfPt> P1 = round(Rw, zp, r, 3, 1);
        const AnchoredRoundCert c = certifyAnchoredRound(P1, {wall, plane}, kTau);
        expect(!c.ok && std::strcmp(c.why, "underdetermined") == 0,
               "one column: underdetermined");
    }
    std::fprintf(stderr, "%s (%d failure%s)\n", gFail ? "FAIL" : "PASS", gFail,
                 gFail == 1 ? "" : "s");
    return gFail ? 1 : 0;
}
