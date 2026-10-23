// D-fc-6.2 (Class A). The caller shell gate compares magnitude with magnitude.
// shellVol is the signed shell mass; meshAbs is the mesh volume (the caller
// passes |mesh|); budget is the caller's existing
// max(1e-4·|mesh|, 3·Σ|dVol|). This predicate adds no threshold.
// A signed-negative shell whose magnitude matches |mesh| agrees; a magnitude
// miss strictly beyond the budget does not. The shipped solid is then oriented
// by the same makeSolids step as the off-path (Verbatim) solid.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_SHELL_VOLUME_HPP
#define STL2STEP_SHELL_VOLUME_HPP

#include <cmath>

namespace stl2step {

inline bool shellMagnitudesAgree(double shellVol, double meshAbs, double budget) {
    return !(std::fabs(std::fabs(shellVol) - std::fabs(meshAbs)) > budget);
}

}  // namespace stl2step

#endif
