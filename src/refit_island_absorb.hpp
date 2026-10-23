// D-fc-6.1(d). Absorb an unclaimed facet island into the one closed360 cylinder
// that already owns every one of its vertices. hostSagitta is the cylinder's own
// mesh sagitta; the lateral gate is sin(thetaSharp), the B1 member quantity.
// No refit. Bodies live in refit_island_absorb.cpp: refit_chains.cpp is on the
// D5.3 include allowlist and cannot include this header.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_ISLAND_ABSORB_HPP
#define STL2STEP_REFIT_ISLAND_ABSORB_HPP

#include "refit.hpp"

#include <vector>

namespace stl2step {
namespace refit {

// max over t in S of |dist(centroid(t), axis_S) - R_S|. A function of the mesh.
double hostSagitta(const MeshView& mv, const Region& S);

// One host. False when any vertex is new, the island sagitta exceeds hostSagitta,
// a triangle is not lateral, or the island shares no edge with S.
bool f1WouldAbsorb(const MeshView& mv, const Region& S, const std::vector<int>& island,
                   double lateralSin, double* sagOut);

// Index of the single closed360 cylinder that would absorb `island`, or false.
bool uniqueF1Host(const MeshView& mv, const std::vector<Region>& regions,
                  const std::vector<int>& island, double lateralSin, int& host,
                  double& sag, double& hostSag);

}  // namespace refit
}  // namespace stl2step

#endif
