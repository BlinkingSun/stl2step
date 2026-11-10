// Synthetic cases for sphereOutranksCylinder (D-gid-8.3, R-k31-2): a cylinder
// region's triangles belong to a sphere corner only when every vertex is within
// tau of BOTH the region's own cylinder and the sphere. No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_relational_sphere.hpp"

#include <gp_XYZ.hxx>

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kQ = 1e-6;
static constexpr double kTau = 2 * kQ;

// Sphere (centre C, radius R) and a cylinder through C along z with the same
// radius: the two surfaces share the equator, and a point on the sphere at
// height h is R - sqrt(R^2 - h^2) ~ h^2 / (2R) off the cylinder.
static const gp_XYZ C(1.0, 2.0, 3.0);
static constexpr double R = 3.0;

static gp_XYZ onSphere(double az, double h) {
    const double rho = std::sqrt(R * R - h * h);
    return C + gp_XYZ(rho * std::cos(az), rho * std::sin(az), h);
}
static gp_XYZ onCylinder(double az, double h) {
    return C + gp_XYZ(R * std::cos(az), R * std::sin(az), h);
}

// A strip of triangles between heights h0 and h1 over azimuths [a0, a1].
static std::vector<std::array<gp_XYZ, 3>> strip(double a0, double a1, int n, double h0, double h1,
                                                gp_XYZ (*at)(double, double)) {
    std::vector<std::array<gp_XYZ, 3>> t;
    for (int i = 0; i < n; ++i) {
        const double u0 = a0 + (a1 - a0) * i / n, u1 = a0 + (a1 - a0) * (i + 1) / n;
        t.push_back({at(u0, h0), at(u1, h0), at(u1, h1)});
        t.push_back({at(u0, h0), at(u1, h1), at(u0, h1)});
    }
    return t;
}

int main() {
    RoundSpine cyl;
    cyl.circle = false;
    cyl.o = C;
    cyl.a = gp_XYZ(0, 0, 1);
    cyl.r = R;

    // Heights: 1e-3 puts a sphere point 1.7e-7 (0.08 tau) off the cylinder;
    // 1e-2 puts it 1.7e-5 (8 tau) off.
    const double hIn = 1e-3, hOut = 1e-2;

    const auto band = strip(0.0, kPi / 12, 6, -hIn, hIn, onSphere);
    expect(sphereOutranksCylinder(C, R, cyl, band, kTau),
           "equator band within tau of both surfaces -> owned by the sphere");

    auto tallS = band;
    tallS.back()[2] = onSphere(kPi / 24, hOut);
    expect(!sphereOutranksCylinder(C, R, cyl, tallS, kTau),
           "one vertex on the sphere but 8 tau off its own cylinder -> not owned");

    auto tallC = band;
    tallC.front()[1] = onCylinder(kPi / 24, hOut);
    expect(!sphereOutranksCylinder(C, R, cyl, tallC, kTau),
           "one vertex on its own cylinder but 8 tau off the sphere -> not owned");

    const auto cylOnly = strip(0.0, kPi / 12, 6, -hOut, hOut, onCylinder);
    expect(!sphereOutranksCylinder(C, R, cyl, cylOnly, kTau),
           "a cylinder band taller than the shared equator -> not owned");

    expect(!sphereOutranksCylinder(C, R + 10 * kTau, cyl, band, kTau),
           "sphere radius 10 tau off -> not owned");

    expect(!sphereOutranksCylinder(C, R, cyl, {}, kTau), "no triangles -> false");
    expect(!sphereOutranksCylinder(C, R, cyl, band, 0.0), "tau 0 admits nothing");
    expect(!sphereOutranksCylinder(C, R, cyl, band, -kTau), "negative tau admits nothing");
    expect(!sphereOutranksCylinder(C, R, cyl, band, std::numeric_limits<double>::quiet_NaN()),
           "NaN tau admits nothing");

    RoundSpine circle = cyl;
    circle.circle = true;
    circle.R = 5.0;
    expect(!sphereOutranksCylinder(C, R, circle, band, kTau), "a circle spine is not a cylinder");

    RoundSpine noAxis = cyl;
    noAxis.a = gp_XYZ(0, 0, 0);
    expect(!sphereOutranksCylinder(C, R, noAxis, band, kTau), "a zero axis is refused");

    if (gFail) {
        std::fprintf(stderr, "relational_sphere_owns_unit: %d FAILED\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_sphere_owns_unit: all passed\n");
    return 0;
}
