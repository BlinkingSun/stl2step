// Pins the shared closed forms (D-fc-1.4, D-fc-1.6). No extrude algorithm.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"

#include <cmath>
#include <cstdio>

static int gFail = 0;

static void expect(bool ok, const char* name) {
    if (ok) {
        std::fprintf(stderr, "PASS %s\n", name);
        return;
    }
    std::fprintf(stderr, "FAIL %s\n", name);
    ++gFail;
}

int main() {
    using namespace stl2step::refit;

    MeshView qHalf{};
    qHalf.quantFloor = 0.5;
    MeshView qZero{};
    qZero.quantFloor = 0.0;
    const double tauHalf = tauOf(qHalf);
    const double tauZero = tauOf(qZero);
    std::fprintf(stderr, "tauOf(0.5)=%.17g tauOf(0)=%.17g\n", tauHalf, tauZero);
    expect(tauHalf == 1.0, "tauOf(0.5) = 1");
    expect(tauZero == 0.0, "tauOf(0) = 0");

    const double a1 = anchorTol(1.0, 3.0);
    const double a2 = anchorTol(4.0, 1.0);
    std::fprintf(stderr, "anchorTol(1,3)=%.17g anchorTol(4,1)=%.17g\n", a1, a2);
    expect(a1 == 3.0, "anchorTol(1, 3) = 3");
    expect(a2 == 4.0, "anchorTol(4, 1) = 4");

    const double s1 = spanTol(2.0, 8.0, 1.0);
    const double s2 = spanTol(0.5, 0.5, 1.0);
    std::fprintf(stderr, "spanTol(2,8,1)=%.17g spanTol(0.5,0.5,1)=%.17g\n", s1, s2);
    expect(s1 == 6.0, "spanTol(2, 8, 1) = 6");
    expect(s2 == 2.0, "spanTol(0.5, 0.5, 1) = 2");

    // areaMin = epsPlane * diag = 0.25 * 8 = 2. A 2-triangle provisional at
    // area 1 is shatter-class. A 3-triangle provisional at the same area is not.
    MeshView mv{};
    mv.diag = 8.0;
    DerivedTols tol;
    tol.epsPlane = 0.25;
    Provisional two;
    two.claim = ProvClaim::Unclaimed;
    two.tris = {0, 1};
    two.area = 1.0;
    Provisional three = two;
    three.tris = {0, 1, 2};
    const bool c2 = a3WouldCommit(mv, tol, two, false);
    const bool c3 = a3WouldCommit(mv, tol, three, false);
    std::fprintf(stderr, "a3WouldCommit(n=2,area=1)=%d a3WouldCommit(n=3,area=1)=%d\n",
                 c2 ? 1 : 0, c3 ? 1 : 0);
    expect(!c2, "a3WouldCommit 2-triangle tiny = false");
    expect(c3, "a3WouldCommit 3-triangle = true");

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_extrude_unit ok\n");
    return 0;
}
