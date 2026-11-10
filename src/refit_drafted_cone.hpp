// Drafted wall about a certified plane normal (D-gid-8.18). Body in
// refit_drafted_cone.cpp. Every threshold is tau = 2q or a count of
// parameters (D-fc-1.6).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_DRAFTED_CONE_HPP
#define STL2STEP_REFIT_DRAFTED_CONE_HPP

#include "refit_profile_cert.hpp"

#include <gp_XYZ.hxx>

#include <cstddef>
#include <vector>

namespace stl2step {
namespace refit {

// D-gid-8.18. A drafted wall whose axis direction is a certified plane's
// normal (the pull direction). The direction is the anchor's, so the free
// parameters are the plane-section circle (centre 2, radius 1) and the draft
// (1).
constexpr int kParamsDraftCone = 4;

// D-gid-8.18. The certificate on a point set (every vertex of a band), the
// axis direction `n` taken from a certified plane through `p0`. Order
// (D-fc-1.6.2, fewest parameters first):
//  - stations (section points merged at tau, countProfPts) >=
//    minStations(kParamsDraftCone); the stations lie on more than two axial
//    levels (countLevelClusters >= minStations(kParamsLine): the profile line,
//    radius and draft, is over-determined; D-gid-8.25 (a) two heights are one
//    level unless more than levelSep(tau) = 2 tau apart); and one level holds
//    kParamsCircle stations, so that level's circle fixes the centre (with two
//    stations per level the centre slides along their bisector and every
//    further station is absorbed by it), else "underdetermined";
//  - D-gid-8.29 ambiguity: every point within tau of one plane -> "plane";
//    of one cylinder about n -> "cylinder"; of one sphere (four parameters,
//    the cone's count) -> "sphere" (ambiguousWithCheaperFamily);
//  - the least-squares cone about n (orthogonal distance, iteration cap
//    kDraftFitIters); D-gid-8.25 (b) its radius changes over the band by more
//    than levelSep(tau) (the draft is resolvable), else "cylinder"; the cone's
//    radius at both axial ends of the points exceeds tau (no apex inside the
//    band), else "apex";
//  - every point within tau of that cone, else "vertex-off-surface";
//  - else "accept".
// `rhoP` and `k` describe the profile rho(h) = rhoP + k * h, h the height
// above `p0` along `axis`; `c0` is the axis point at h = 0.
// D-gid-8.25 (a). The certificate's diameter: a point within tau of the
// surface lies anywhere in a slab of width 2 tau, so two heights (or two radii)
// closer than this are one to the certificate.
inline double levelSep(double tau) { return 2.0 * tau; }

struct DraftedConeFit {
    bool ok = false;
    bool determinate = false;
    const char* why = "underdetermined";
    gp_XYZ axis{0, 0, 1};
    gp_XYZ c0{0, 0, 0};
    double rhoP = 0, k = 0;
    double hLo = 0, hHi = 0;
    double maxDev = 0, rmsDev = 0;
    double planeDev = 0, cylDev = 0, sphereDev = 0;
    int stations = 0, levels = 0, levelStations = 0;
};

// D-gid-8.25 (a). The axial levels of heights h: sorted, a new level at a gap
// wider than levelSep(tau). Point indices per level, ascending height.
std::vector<std::vector<std::size_t>> levelClusters(const std::vector<double>& h, double tau);
int countLevelClusters(const std::vector<double>& h, double tau);
// The most section stations (countProfPts) on one level of levelClusters.
int maxLevelStations(const std::vector<double>& h, const std::vector<ProfPt>& sec, double tau);

DraftedConeFit certifyDraftedCone(const std::vector<gp_XYZ>& pts, const gp_XYZ& n,
                                  const gp_XYZ& p0, double tau);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_DRAFTED_CONE_HPP
