// Relational reclaim stub (W1.0 contract). Algorithm belongs to lane b1-reclaim.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

namespace stl2step {
namespace refit {

// Lane b1-reclaim (W1.3). No-op stub; admits nothing (D-fc-1.11).
bool reclaimRelational(const MeshView&, const SegmentParams&, const DerivedTols&,
                       SegmentWork&, RelRevolveMode) {
    return true;
}

}}  // namespace stl2step::refit
