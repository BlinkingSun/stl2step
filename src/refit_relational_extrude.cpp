// Relational extrude stub (W1.0 contract). Algorithm belongs to lane rel-extrude.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

namespace stl2step {
namespace refit {

// Lane rel-extrude (W1.1). No-op stub; admits nothing (D-fc-1.3).
bool claimRelExtrude(const MeshView&, const SegmentParams&, const DerivedTols&,
                     SegmentWork&) {
    return true;
}

}}  // namespace stl2step::refit
