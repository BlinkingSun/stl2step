// Anchored quarter-round torus (D-gid-2.6 / D-gid-3.5). Body in
// refit_anchored_round.cpp. Every threshold is tau = 2q or a count of
// parameters (D-fc-1.6). certifyTorus is not edited: it governs FREE torus
// fits; this certificate is the anchored mode beside it.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_ANCHORED_ROUND_HPP
#define STL2STEP_REFIT_ANCHORED_ROUND_HPP

#include "refit_profile_cert.hpp"

#include <vector>

namespace stl2step {
namespace refit {

// A certified anchor seen in the (rho, z) profile half-plane of the round's
// revolution axis: the line nr * rho + nz * z = d, (nr, nz) unit. A coaxial
// cylinder of radius R is (1, 0, R); a plane normal to the axis at height h is
// (0, 1, h).
struct ProfileAnchor {
    double nr = 0, nz = 0, d = 0;
};

// D-gid-2.6. The profile circle of a round tangent to its certified anchors.
// Tangency to an anchor fixes one of the circle's three parameters, so one
// anchor leaves 2 free (contact position along it and r) and two anchors leave
// 1 (r). Zero anchors is a FREE fit: certifyTorus governs it (why "free-round").
// Order (D-fc-1.6.2):
//  - stations (countProfPts) >= minStations(free parameters), columns >= 2
//    (the profile is seen at two azimuths), else "underdetermined";
//  - the profile line first: every station within tau of one line -> "line";
//  - convexity from the fit: each anchor's side is the side the free circle's
//    centre (fitProfileCircle) lies on; no side is assumed;
//  - the anchored least-squares circle; r <= tau "degenerate", centre on or
//    inside the axis (Rmaj <= r) "spindle-torus";
//  - every station within tau of that circle (certifyTorus's circle residual),
//    else "vertex-off-surface" (determinate);
//  - the normals span the anchor dihedral: every anchor's tangent point is a
//    station within tau (the round meets each anchor), and with two anchors
//    every station lies on the arc between the two tangent points that faces
//    the anchors' corner, within tau of it, else "not-blend".
struct AnchoredRoundCert {
    bool ok = false;
    bool determinate = false;
    const char* why = "none";
    double Rmaj = 0, z0 = 0, Rmin = 0;
    double maxDev = 0, rmsDev = 0;
    int stations = 0, cols = 0, params = 0, anchors = 0;
};
AnchoredRoundCert certifyAnchoredRound(const std::vector<ProfPt>& P,
                                       const std::vector<ProfileAnchor>& anchors, double tau);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_ANCHORED_ROUND_HPP
