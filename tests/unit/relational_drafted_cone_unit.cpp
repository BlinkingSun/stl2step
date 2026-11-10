// Synthetic drafted-wall cases (D-gid-8.18). certifyDraftedCone on points of a
// cone whose axis is a given plane normal: a wall between one level row and a
// row on a tilted plane certifies and recovers its parameters; two flat rows,
// a row spread within 2 tau (one level), two stations per level, an
// unresolvable draft, a cylinder, a plane, a torus zone and two rows of a
// sphere (the D-gid-8.29 ambiguity, equal parameter count) refuse.
// The chord volume of a partial band (revolvedChordTerm + coneRimTerm) matches
// a ray-integrated reference (and misses it without the rim term); certifyAnchoredRound takes a cone's profile line as
// its anchor. Generic dimensions only, no training mesh.
// SPDX-License-Identifier: MIT

#include "refit_anchored_round.hpp"
#include "refit_band_volume.hpp"
#include "refit_drafted_cone.hpp"
#include "refit_profile_cert.hpp"

#include <gp_XYZ.hxx>

#include <array>
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

// Generic synthetic cone: axis +z through (cx, cy), radius rho0 at z = 0,
// draft half-angle a (radians): rho(z) = rho0 + tan(a) z.
struct Cone {
    double cx = 3.0, cy = -2.0, rho0 = 5.0, a = 2.0 * kPi / 180.0;
    double rho(double z) const { return rho0 + std::tan(a) * z; }
    gp_XYZ at(double phi, double z) const {
        return gp_XYZ(cx + rho(z) * std::cos(phi), cy + rho(z) * std::sin(phi), z);
    }
};

// A wall patch: a bottom row on the level z = z0 at nB azimuths over [p0, p1],
// a top row at nT other azimuths whose heights follow a plane tilted by `tilt`
// about the x axis (z = z1 + tilt * y).
static std::vector<gp_XYZ> wall(const Cone& c, double z0, double z1, double tilt, int nB, int nT,
                                double p0, double p1) {
    std::vector<gp_XYZ> P;
    for (int i = 0; i < nB; ++i) P.push_back(c.at(p0 + (p1 - p0) * i / (nB - 1), z0));
    for (int i = 0; i < nT; ++i) {
        const double phi = p0 + (p1 - p0) * (i + 0.5) / nT;
        // Solve z = z1 + tilt * y(phi, z) by fixed point (converges: |tilt| small).
        double z = z1;
        for (int it = 0; it < 60; ++it) z = z1 + tilt * c.at(phi, z).Y();
        P.push_back(c.at(phi, z));
    }
    return P;
}

int main() {
    const Cone c;
    const gp_XYZ n(0, 0, 1), p0(0, 0, 0);

    // 1. One level row + a row on a tilted plane: more than two levels.
    {
        const std::vector<gp_XYZ> P = wall(c, 1.0, 9.0, 0.01, 5, 5, 0.2, 1.4);
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(f.ok && std::strcmp(f.why, "accept") == 0, "tilted-top wall certifies");
        expect(f.levels >= 3 && f.levelStations >= 3 && f.stations >= 5,
               "tilted-top wall: > 2 levels, a determined level circle, > 4 stations");
        expect(std::fabs(f.rhoP - c.rho0) <= kTau && std::fabs(f.k - std::tan(c.a)) * 9.0 <= kTau,
               "tilted-top wall: radius at the plane and draft recovered within tau");
        expect(std::fabs(f.c0.X() - c.cx) <= kTau && std::fabs(f.c0.Y() - c.cy) <= kTau,
               "tilted-top wall: centre recovered within tau");
        expect(f.maxDev <= kTau, "tilted-top wall: every vertex within tau");
        // The same wall with its axis taken from a plane normal off by 1e-3 rad.
        const gp_XYZ nOff(std::sin(1e-3), 0, std::cos(1e-3));
        const DraftedConeFit g = certifyDraftedCone(P, nOff, p0, kTau);
        expect(!g.ok, "a normal that is not the axis does not certify the wall");
        // A vertex moved off the cone by 3 tau.
        std::vector<gp_XYZ> Q = P;
        Q[2] = Q[2] + gp_XYZ(std::cos(0.5), std::sin(0.5), 0) * (3.0 * kTau);
        const DraftedConeFit h = certifyDraftedCone(Q, n, p0, kTau);
        expect(!h.ok && std::strcmp(h.why, "vertex-off-surface") == 0,
               "a vertex 3 tau off refuses vertex-off-surface");
    }
    // 2. Two flat rows: two levels only, the draft is exactly determined.
    {
        const std::vector<gp_XYZ> P = wall(c, 1.0, 9.0, 0.0, 5, 5, 0.2, 1.4);
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(!f.ok && std::strcmp(f.why, "underdetermined") == 0 && f.levels == 2,
               "two flat rows refuse underdetermined (2 levels)");
    }
    // 2b. D-gid-8.25 (a): a top row whose heights spread over 1.5 tau is one
    // level (heights within the certificate's diameter 2 tau are one level).
    {
        const std::vector<gp_XYZ> P = wall(c, 1.0, 9.0, 1.5 * kTau / 3.9, 5, 5, 0.2, 1.4);
        double zLo = 1e300, zHi = -1e300;
        for (std::size_t i = 5; i < P.size(); ++i) {
            zLo = std::min(zLo, P[i].Z());
            zHi = std::max(zHi, P[i].Z());
        }
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(zHi - zLo > kTau && zHi - zLo <= 2.0 * kTau,
               "top row spread is between tau and 2 tau");
        expect(!f.ok && std::strcmp(f.why, "underdetermined") == 0 && f.levels == 2,
               "a top row within the certificate's diameter is one level: underdetermined");
    }
    // 2c. D-gid-8.25 (b): a draft whose radius change over the band is within
    // 2 tau is not resolvable -- a cylinder to the certificate -- even when one
    // vertex 0.6 tau off keeps the cylinder fit over tau.
    {
        Cone cd = c;
        cd.a = std::atan(1.9 * kTau / 8.0);
        std::vector<gp_XYZ> P = wall(cd, 1.0, 9.0, 0.01, 5, 5, 0.2, 1.4);
        P[0] = P[0] + gp_XYZ(std::cos(0.2), std::sin(0.2), 0) * (0.6 * kTau);
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        std::fprintf(stderr, "  unresolvable draft: cylinder fit %.3f tau\n", f.cylDev / kTau);
        expect(f.cylDev > kTau, "the cylinder about the normal misses one vertex by more than tau");
        expect(!f.ok && std::strcmp(f.why, "cylinder") == 0,
               "a draft within 2 tau over the band is a cylinder, not a cone");
    }
    // 3. Two stations per level on three levels: the centre slides on the bisector.
    {
        std::vector<gp_XYZ> P;
        for (double z : {1.0, 2.0, 3.0})
            for (double phi : {0.3, 0.6}) P.push_back(c.at(phi, z));
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(!f.ok && std::strcmp(f.why, "underdetermined") == 0 && f.levelStations == 2,
               "two stations per level refuse underdetermined");
    }
    // 4. A cylinder about the normal (no draft) is a cylinder, not a cone.
    {
        Cone cy = c;
        cy.a = 0.0;
        const std::vector<gp_XYZ> P = wall(cy, 1.0, 9.0, 0.01, 5, 5, 0.2, 1.4);
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(!f.ok && std::strcmp(f.why, "cylinder") == 0, "an undrafted wall is a cylinder");
    }
    // 5. Coplanar points are a plane.
    {
        std::vector<gp_XYZ> P;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 3; ++j) P.push_back(gp_XYZ(i, 0.5 * j, 0.25 * i));
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(!f.ok && std::strcmp(f.why, "plane") == 0, "coplanar points are a plane");
    }
    // 6. A torus zone about the normal (a curved profile that is not a sphere)
    // on four levels refuses: no family certifies it.
    {
        std::vector<gp_XYZ> P;
        const double Rmaj = 4.0, r = 2.0;
        for (double z : {0.5, 1.0, 1.5, 1.9}) {
            const double rho = Rmaj + std::sqrt(r * r - z * z);
            for (int i = 0; i < 4; ++i) {
                const double phi = 0.2 + 0.3 * i;
                P.push_back(gp_XYZ(c.cx + rho * std::cos(phi), c.cy + rho * std::sin(phi), z));
            }
        }
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(!f.ok && std::strcmp(f.why, "vertex-off-surface") == 0,
               "a torus zone on four levels refuses vertex-off-surface");
    }
    // 6b. D-gid-8.29 ambiguity, the sphere tie: two rows of a sphere (radius 3,
    // centre on the axis) beside its equator -- a level row and a tilted row
    // 0.04 above it, heights spread 0.0004 -- are also a 0.4 degree drafted
    // cone within tau. A sphere has the cone's four parameters, so the band is
    // ambiguous and refuses "sphere".
    {
        const double Rs = 3.0, H = 0.04, spread = 0.0004;
        std::vector<gp_XYZ> P;
        std::vector<double> rho, h;
        for (int i = 0; i < 7; ++i) {
            const double phi = 0.1 * i;
            P.push_back(gp_XYZ(c.cx + Rs * std::cos(phi), c.cy + Rs * std::sin(phi), 0.0));
            rho.push_back(Rs);
            h.push_back(0.0);
        }
        for (int i = 0; i < 7; ++i) {
            const double phi = 0.02 + 0.56 * i / 6.0;
            const double z = H + spread * i / 6.0;
            const double rr = std::sqrt(Rs * Rs - z * z);
            P.push_back(gp_XYZ(c.cx + rr * std::cos(phi), c.cy + rr * std::sin(phi), z));
            rho.push_back(rr);
            h.push_back(z);
        }
        // Independently: the cone about the axis through the sphere's centre,
        // rho = a + b h by least squares, orthogonal residual.
        double sh = 0, shh = 0, sr = 0, shr = 0;
        const double nn = (double)h.size();
        for (std::size_t i = 0; i < h.size(); ++i) {
            sh += h[i];
            shh += h[i] * h[i];
            sr += rho[i];
            shr += h[i] * rho[i];
        }
        const double b = (nn * shr - sh * sr) / (nn * shh - sh * sh);
        const double a = (sr - b * sh) / nn;
        double coneMax = 0.0;
        for (std::size_t i = 0; i < h.size(); ++i)
            coneMax = std::max(coneMax, std::fabs(rho[i] - a - b * h[i]) / std::sqrt(1.0 + b * b));
        std::fprintf(stderr, "  sphere tie: cone %.3f tau, draft %.4f deg\n", coneMax / kTau,
                     std::atan(b) * 180.0 / kPi);
        expect(coneMax <= kTau, "the two sphere rows are also a drafted cone within tau");
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(f.determinate && f.levels >= 3 && f.cylDev > kTau && f.planeDev > kTau,
               "sphere tie: determinate, not a plane, not a cylinder");
        expect(!f.ok && std::strcmp(f.why, "sphere") == 0 && f.sphereDev <= kTau,
               "a band that is also a sphere at tau is ambiguous: refuses sphere");
    }
    // 6d. D-gid-8.30: the ambiguity refusal does not depend on the frame. Two
    // rows of a sphere (radius 3, centre on the axis), the level row at 0.12
    // azimuth steps, the tilted row 0.04 above it over a 0.7 azimuth span,
    // heights spread 0.0004; the frame origin 8 m off the rows' plane and far
    // off the axis, then on the plane.
    {
        const double Rs = 3.0, H = 0.04, spread = 0.0004;
        std::vector<gp_XYZ> P;
        for (int i = 0; i < 7; ++i) {
            const double phi = 0.1 + 0.12 * i;
            P.push_back(gp_XYZ(c.cx + Rs * std::cos(phi), c.cy + Rs * std::sin(phi), 0.0));
        }
        for (int i = 0; i < 7; ++i) {
            const double phi = 0.13 + 0.7 * i / 6.0;
            const double z = H + spread * i / 6.0;
            const double rr = std::sqrt(Rs * Rs - z * z);
            P.push_back(gp_XYZ(c.cx + rr * std::cos(phi), c.cy + rr * std::sin(phi), z));
        }
        for (const gp_XYZ& o : {gp_XYZ(25000.0, -17000.0, 8000.0), gp_XYZ(25000.0, -17000.0, 0.0),
                                gp_XYZ(0.0, 0.0, 0.0)}) {
            const DraftedConeFit g = certifyDraftedCone(P, n, o, kTau);
            std::fprintf(stderr, "  frame origin (%g, %g, %g): why=%s sphere %.3g tau\n", o.X(),
                         o.Y(), o.Z(), g.why, g.sphereDev / kTau);
            expect(!g.ok && std::strcmp(g.why, "sphere") == 0 && g.sphereDev <= kTau,
                   "sphere-tie refusal is independent of the frame origin");
        }
    }
    // 6c. The real drafted wall of case 1 is not a sphere: it still certifies.
    {
        const std::vector<gp_XYZ> P = wall(c, 1.0, 9.0, 0.01, 5, 5, 0.2, 1.4);
        const DraftedConeFit f = certifyDraftedCone(P, n, p0, kTau);
        expect(f.ok && f.sphereDev > kTau, "a drafted wall misses every sphere by more than tau");
    }
    // 7. Chord volume of a partial band with oblique ends and a tilted top rim:
    // revolvedChordTerm + coneRimTerm on the rim against a ray-integrated reference.
    {
        const int nB = 6, nT = 5;
        const double z0 = 1.0, z1 = 9.0, tilt = 0.01, pA = 0.2, pB = 1.4;
        std::vector<gp_XYZ> B, T;
        for (int i = 0; i < nB; ++i) B.push_back(c.at(pA + (pB - pA) * i / (nB - 1), z0));
        for (int i = 0; i < nT; ++i) {
            const double phi = pA + 0.03 + (pB - pA - 0.06) * i / (nT - 1);
            double z = z1;
            for (int it = 0; it < 60; ++it) z = z1 + tilt * c.at(phi, z).Y();
            T.push_back(c.at(phi, z));
        }
        // Triangulate the strip between the rows, both rows ascending in phi,
        // every triangle counter-clockwise seen from outside.
        std::vector<std::array<gp_XYZ, 3>> tris;
        int i = 0, j = 0;
        while (i < nB - 1 || j < nT - 1) {
            const bool advB = (j == nT - 1) || (i < nB - 1 && (pA + (pB - pA) * (i + 1) / (nB - 1)) <=
                                                                  (pA + 0.03 + (pB - pA - 0.06) * (j + 1) / (nT - 1)));
            if (advB) {
                tris.push_back({B[i], B[i + 1], T[j]});
                ++i;
            } else {
                tris.push_back({B[i], T[j + 1], T[j]});
                ++j;
            }
        }
        const gp_XYZ o(c.cx, c.cy, 0.0), ax(0, 0, 1), e1(1, 0, 0), e2(0, 1, 0);
        const double slope = std::tan(c.a);
        double vol = 0.0, volNoRim = 0.0;
        // Rim edges: every edge used once.
        std::vector<std::array<gp_XYZ, 2>> rim;
        auto same = [](const gp_XYZ& a, const gp_XYZ& b) { return (a - b).Modulus() < 1e-12; };
        for (std::size_t q = 0; q < tris.size(); ++q) {
            const auto& t = tris[q];
            vol += revolvedChordTerm(t[0], t[1], t[2], o, ax, e1, e2, c.rho0, slope);
            volNoRim += revolvedChordTerm(t[0], t[1], t[2], o, ax, e1, e2, c.rho0, slope);
            const int sg = revolvedParamSign(t[0], t[1], t[2], o, ax, e1, e2);
            for (int k = 0; k < 3; ++k) {
                const gp_XYZ& a = t[k];
                const gp_XYZ& b = t[(k + 1) % 3];
                int uses = 0;
                for (const auto& u : tris)
                    for (int m = 0; m < 3; ++m) {
                        const gp_XYZ& x = u[m];
                        const gp_XYZ& y = u[(m + 1) % 3];
                        if ((same(x, a) && same(y, b)) || (same(x, b) && same(y, a))) ++uses;
                    }
                if (uses == 1) vol += sg * coneRimTerm(a, b, o, ax, e1, e2, c.rho0, slope);
            }
        }
        // Reference, independent of the flux identity: the volume between the
        // facets and the cone along the horizontal rays from the axis,
        //   integral (rho_cone(z)^2 - rho_facet(theta, z)^2) / 2  dtheta dz,
        // by the midpoint rule on a fine (theta, z) grid; a ray that meets no
        // facet adds nothing.
        auto rayVolume = [&](int nTh, int nZ) {
            const double th0 = pA - 0.01, th1 = pB + 0.01, zz0 = z0 - 0.01, zz1 = z1 + 0.2;
            const double dth = (th1 - th0) / nTh, dz = (zz1 - zz0) / nZ;
            double v = 0.0;
            for (int a = 0; a < nTh; ++a) {
                const double th = th0 + (a + 0.5) * dth;
                const gp_XYZ dir(std::cos(th), std::sin(th), 0.0);
                for (int b = 0; b < nZ; ++b) {
                    const double z = zz0 + (b + 0.5) * dz;
                    const gp_XYZ base(c.cx, c.cy, z);
                    for (const auto& t : tris) {
                        const gp_XYZ nn = (t[1] - t[0]).Crossed(t[2] - t[0]);
                        const double den = nn.Dot(dir);
                        if (std::fabs(den) < 1e-15) continue;
                        const double rf = nn.Dot(t[0] - base) / den;
                        if (!(rf > 0.0)) continue;
                        const gp_XYZ x = base + dir * rf;
                        // inside: same side of all three edges as the normal
                        bool in = true;
                        for (int k = 0; k < 3 && in; ++k)
                            if ((t[(k + 1) % 3] - t[k]).Crossed(x - t[k]).Dot(nn) < 0.0) in = false;
                        if (!in) continue;
                        const double rc = c.rho(z);
                        v += 0.5 * (rc * rc - rf * rf) * dth * dz;
                        break;
                    }
                }
            }
            return v;
        };
        const double ref = rayVolume(1600, 1600);
        double area = 0.0;
        for (const auto& t : tris) area += 0.5 * (t[1] - t[0]).Crossed(t[2] - t[0]).Modulus();
        std::fprintf(stderr, "  chord volume %.9g reference %.9g area %.4g\n", vol, ref, area);
        expect(std::fabs(std::fabs(vol) - std::fabs(ref)) <= area * kTau,
               "partial band chord volume = reference within area * tau");
        expect(std::fabs(std::fabs(volNoRim) - std::fabs(ref)) > area * kTau,
               "without the rim term the oblique rim misses the reference");
    }
    // 8. certifyAnchoredRound with a cone's profile line as its anchor: a round
    // tangent to the cone wall and to a plane normal to the axis.
    {
        const double slope = std::tan(c.a), zp = -1.0, r = 0.8;
        const double s = std::sqrt(1.0 + slope * slope);
        const ProfileAnchor coneLine{1.0 / s, -slope / s, c.rho0 / s};
        const ProfileAnchor plane{0.0, 1.0, zp};
        // Centre inside the cone (side -1) and above the plane (side +1).
        // n_c . C = d_c - r, n_p . C = zp + r.
        const double Cz = zp + r;
        const double Cr = (c.rho0 / s - r - coneLine.nz * Cz) / coneLine.nr;
        // Tangent points: on the plane (Cr, zp); on the cone C + r n_c.
        const double a0 = std::atan2(coneLine.nz, coneLine.nr);
        std::vector<ProfPt> P;
        for (int k = 0; k <= 3; ++k) {
            const double t = -0.5 * kPi + (a0 + 0.5 * kPi) * (double)k / 3.0;
            for (int j = 0; j < 4; ++j) {
                ProfPt q;
                q.rho = Cr + r * std::cos(t);
                q.z = Cz + r * std::sin(t);
                q.psi = 0.3 * j;
                P.push_back(q);
            }
        }
        const AnchoredRoundCert cert = certifyAnchoredRound(P, {coneLine, plane}, kTau);
        expect(cert.ok && cert.params == 1, "round on a cone anchor and a plane certifies (1 free)");
        expect(std::fabs(cert.Rmin - r) <= kTau && std::fabs(cert.Rmaj - Cr) <= kTau &&
                   std::fabs(cert.z0 - Cz) <= kTau,
               "round on a cone anchor: circle recovered within tau");
    }
    std::fprintf(stderr, "%s\n", gFail ? "FAILED" : "ALL PASS");
    return gFail ? 1 : 0;
}
