// D-fc-7.4 Q2. Period shift for a partial cylinder's pcurves.
// A circle edge and a generator of the same chart can be written a whole turn
// apart. The shift that pulls the edge onto the chart is exact (one period).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_PARTIAL_CYL_HPP
#define STL2STEP_REFIT_PARTIAL_CYL_HPP

#include <algorithm>
#include <cmath>

namespace stl2step {
namespace refit {

// Whole-period shift (added to u) that pulls a pcurve midpoint onto the
// partial chart [u0, u1], u1 >= u0, span < 2π. Zero when the midpoint is
// already on that period, the chart is not a partial turn, or the image would
// not land closer to the chart by more than half a turn (the cone chart's
// branch fold — a full period is the only move, and half a turn is the fold
// of that period, not an acceptance band).
inline double partialCylPeriodShift(double uMid, double u0, double u1) {
    constexpr double kPiLocal = 3.14159265358979323846264338327950288;
    const double twoPi = 2.0 * kPiLocal;
    const double span = u1 - u0;
    if (!(span > 0.0) || !(span < twoPi)) return 0.0;
    const double k = std::floor((uMid - u0) / twoPi);
    if (k == 0.0) return 0.0;
    const double uShifted = uMid - k * twoPi;
    auto distChart = [&](double u) {
        if (u >= u0 && u <= u1) return 0.0;
        return std::min(std::fabs(u - u0), std::fabs(u - u1));
    };
    if (!(distChart(uShifted) + kPiLocal < distChart(uMid))) return 0.0;
    return -k * twoPi;
}

}  // namespace refit
}  // namespace stl2step

#endif
