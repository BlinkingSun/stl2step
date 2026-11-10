// Class-B volume heal predicates (D-fc-6.7). Header-only quantities; no OCCT.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_CLASSB_HPP
#define STL2STEP_REFIT_CLASSB_HPP

#include "refit_cyl_meet.hpp"

#include <cmath>

namespace stl2step {
namespace refit {

// Per-region share of the caller's RULE 2.2 budget (kR = kCascadeResidFactor).
inline double classBVolumeBound(double dVolPredicted, double area, double tau) {
    if (!(tau > 0.0) || !std::isfinite(tau)) return 0.0;
    const double pred =
        std::isfinite(dVolPredicted) ? kCascadeResidFactor * std::fabs(dVolPredicted) : 0.0;
    const double floor = std::isfinite(area) && area > 0.0 ? area * tau : 0.0;
    return pred + floor;
}

// Own certificate fails: vertex residual over τ. D-gid-7.13: a ring-certified
// wall (Region::ringCertified, D-gid-7.6) passed the ring certificate in
// closeRingWalls; its maxVertexDev is the chord sagitta that certificate
// excuses, so it does not fail here.
inline bool classBOwnCertFails(double maxVertexDev, bool ringCertified, double tau) {
    return !ringCertified && tau > 0.0 && std::isfinite(maxVertexDev) && maxVertexDev > tau;
}

// dissolve(r) ⇔ shipped ∧ cert fails ∧ |swept| > bound.
inline bool classBDissolveQualifies(bool shippedAnalytic, double maxVertexDev,
                                    bool ringCertified, double tau, double swept,
                                    double dVolPredicted, double area) {
    if (!shippedAnalytic || !(tau > 0.0)) return false;
    if (!classBOwnCertFails(maxVertexDev, ringCertified, tau)) return false;
    if (!std::isfinite(swept)) return false;
    const double bound = classBVolumeBound(dVolPredicted, area, tau);
    return std::fabs(swept) > bound;
}

// Stop rule: a pass with an empty qualifying set ends the heal (revert with named reason).
inline bool classBHealPassHasCandidates(const int* qualifyingRids, int nQualifying) {
    return nQualifying > 0;
}

// D-gid-2.2 / D-gid-3.8. The swept test qualified nobody (nQual = 0), so the
// caller runs one classBOwnCertFails dissolve from the pristine plan. A
// positive count stays on the swept-test loop. A negative count is not an
// empty set.
inline bool classBNoQualCertFailPass(int nQualifying) {
    return nQualifying == 0;
}

// Member of that one pass: shipped analytic and own certificate fails
// (classBOwnCertFails). nQualifying is the finished swept-test count. The pass
// cannot be selected while that count is still positive, and it does not rank
// or subset the ids.
inline bool classBCertFailPassMember(bool shippedAnalytic, double maxVertexDev,
                                     bool ringCertified, double tau, int nQualifying) {
    return classBNoQualCertFailPass(nQualifying) && shippedAnalytic &&
           classBOwnCertFails(maxVertexDev, ringCertified, tau);
}

}  // namespace refit
}  // namespace stl2step

#endif
