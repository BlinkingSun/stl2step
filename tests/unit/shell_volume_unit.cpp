// D-fc-6.2: the caller shell-volume gate compares magnitude with magnitude.
// A signed-negative shell of equal magnitude passes. A magnitude miss beyond
// the budget still fails, whether the shell mass is negative or positive.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "shell_volume.hpp"

#include <cmath>
#include <cstdio>

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

int main() {
    // Equal magnitude, inward shell. The budget is an argument, not a threshold
    // this test invents: 0.25 is exact in binary so the boundary is exact.
    const double mesh = 8.0;
    const double budget = 0.25;
    expect(stl2step::shellMagnitudesAgree(-mesh, mesh, budget), "neg-equal-passes");
    expect(stl2step::shellMagnitudesAgree(mesh, mesh, budget), "pos-equal-passes");
    expect(!stl2step::shellMagnitudesAgree(-7.0, mesh, budget), "neg-miss-fails");
    expect(!stl2step::shellMagnitudesAgree(7.0, mesh, budget), "pos-miss-fails");
    // The caller fails only when the miss is strictly greater than the budget.
    expect(stl2step::shellMagnitudesAgree(-(mesh + budget), mesh, budget), "boundary-passes");
    expect(!stl2step::shellMagnitudesAgree(-(mesh + budget + 0.01), mesh, budget), "beyond-fails");
    expect(!stl2step::shellMagnitudesAgree(mesh + budget + 0.01, mesh, budget), "pos-beyond-fails");

    // Measured Class A shell (public fixture volume) against its |mesh|.
    expect(stl2step::shellMagnitudesAgree(-384.0, 384.0, 0.0384), "class-a-neg-passes");
    // Both masses positive and apart by more than the budget: still a miss.
    expect(!stl2step::shellMagnitudesAgree(22571.732569, 19706.499583, 73.560550),
           "magnitude-miss-still-fails");

    if (gFail) {
        std::fprintf(stderr, "FAILED %d\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "OK\n");
    return 0;
}
