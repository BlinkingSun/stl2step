// Torus relation predicates shared by H1 (D-fc-4.6). Bodies in
// refit_relational_reclaim.cpp. Every threshold is tau = 2q or a count of
// certifyTorus itself.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_TORUS_HPP
#define STL2STEP_REFIT_RELATIONAL_TORUS_HPP

#include "refit_internal.hpp"
#include "refit_profile_cert.hpp"

#include <array>
#include <vector>

namespace stl2step {
namespace refit {

// F2. An atom is an accepted, non-closed CylGrow cylinder (a B1 commit or a
// stage-L law band). It joins a certified torus only whole, and only when every
// one of its vertices is within tau of that torus.
bool torusAtom(const Region& r);

// Distance of p from the torus of a TorusBlend region (ax.Location is the centre).
double torusSurfaceDev(const Region& torus, const gp_XYZ& p);

// Largest distance of any vertex of `region` from `torus`.
double regionTorusDev(const MeshView& mv, const Region& torus, const Region& region);

// F2, D-fc-1.11(ii)(a) amended: over-determination counted leave-one-region-out.
// triPts[i] holds the three profile points of band triangle i, triAtom[i] the
// atom it belongs to (-1: unclaimed). True when, for every atom k < nAtoms, the
// band without atom k's triangles still meets certifyTorus's counts.
bool leaveOneOutDeterminate(const std::vector<std::array<ProfPt, 3>>& triPts,
                            const std::vector<int>& triAtom, int nAtoms, double tau);

// F4, the rail rule. A torus is a full ring when the largest azimuth gap of its
// own vertices is no larger than the largest azimuth step of one of its own
// boundary edges: the uncovered azimuth is then the inside of one edge.
struct TorusRails {
    double gap = 0;      // largest vertex-azimuth gap (rad)
    double maxStep = 0;  // largest boundary-edge azimuth step (rad)
    bool full = false;
};
TorusRails torusRails(const MeshView& mv, const Region& torus);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_RELATIONAL_TORUS_HPP
