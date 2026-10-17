// Torus relation predicates shared by H1 (D-fc-4.6). Bodies in
// refit_relational_reclaim.cpp. Every threshold is tau = 2q or a count of
// certifyTorus itself.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_TORUS_HPP
#define STL2STEP_REFIT_RELATIONAL_TORUS_HPP

#include "refit_internal.hpp"
#include "refit_profile_cert.hpp"

#include <array>
#include <gp_XYZ.hxx>
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

// D-fc-4.8 (b): one axis for several donors. Least squares over the donors' own
// certified vertex sets: one axis (4 DOF) and one radius per donor (the mean rho
// of its own vertices about the axis). maxDev[g] is the largest |rho - radius[g]|
// over donor g's vertices. False when the fit is not finite.
struct JointAxisFit {
    gp_XYZ origin;  // the axis point nearest the vertex mean
    gp_XYZ axis;    // unit, in the seed's hemisphere
    std::vector<double> radius, maxDev;
};
bool jointAxisFit(const std::vector<std::vector<gp_XYZ>>& groups, const gp_XYZ& originSeed,
                  const gp_XYZ& axisSeed, int maxIters, JointAxisFit& out);

// The radii and deviations of every group about a GIVEN axis (no refit).
JointAxisFit jointAxisAbout(const std::vector<std::vector<gp_XYZ>>& groups, const gp_XYZ& origin,
                            const gp_XYZ& axis);

// The joint certificate: every group's vertices within tau of its own radius.
bool jointAxisCertifies(const JointAxisFit& f, double tau);

// D-fc-4.8 (b): coaxial donor families, formed once in scan order. A donor not
// yet in a family leads one on its own certified frame. A later donor joins when
// the family's axis already certifies every member (it included) within tau, or
// else when the joint least-squares axis over all of them does -- which then
// becomes the family axis. A donor whose own axis differs from the family's by
// more than tanAngTol at the shorter length, or whose axis line stands farther
// from the family axis than its radius plus the family's largest radius, is not
// tested (it cannot enclose that axis); it is never admitted untested.
struct AxisDonor {
    std::vector<gp_XYZ> pts;  // its own certified vertices
    gp_XYZ origin, axis;      // its own certified frame
    double radius = 0, length = 0;
};
struct CoaxialFamily {
    std::vector<int> members;  // donor indices in scan order; front() leads
    JointAxisFit fit;          // the family axis; radius / maxDev per member
    bool refit = false;        // the axis is the joint LS (not the leader's own frame)
};
std::vector<CoaxialFamily> coaxialFamilies(const std::vector<AxisDonor>& donors,
                                           const std::vector<int>& order, double tau,
                                           int maxIters);

// D-fc-4.8 (a): in Torus mode the D-fc-1.6.9 superset of a band is its ring
// chain's candidate set. A band that is a strict subset of a candidate set that
// certifyTorus refuses is a passing subset of a failing band.
bool torusChainSupersetFails(const std::vector<ProfPt>& candidatePts, std::size_t nCandidates,
                             std::size_t nBand, double tau);

// D-fc-4.10 F1: a Torus-mode band that shares a triangle with the candidate set
// of an earlier refused chain (failedTri[t] != 0) is a passing window of a
// failing band, and is refused.
bool torusWindowOverlapsRefused(const std::vector<int>& bandTris,
                                const std::vector<char>& failedTri);

// D-fc-4.11 (F1 refined): a certifyTorus refusal is GEOMETRIC when its evidence
// was determinate (certifyTorus's own counts met) and a geometric test failed:
// line-first, sphere, degenerate, spindle, or a vertex off tau. A refusal for
// over-determination (profile-underdetermined: levels / cols / profPts) or a
// singular circle fit is not.
bool torusRefusalGeometric(const TorusCert& c);

// The full-turn side count of a torus band's azimuthal tessellation, for its
// chord sagitta and dVolPredicted: a closed band's cols columns cover the whole
// turn; a partial band's cols columns are cols - 1 steps over its own azimuth
// span (the rule partial cylinders use, 2 pi nBands / span).
int torusSidesOf(int cols, double span, bool closed);

// The profile (rho, z, psi) of each point about the frame (origin, axis, xHint):
// the Torus scan's frame, psi measured from xHint projected off the axis.
std::vector<ProfPt> profileAbout(const std::vector<gp_XYZ>& pts, const gp_XYZ& origin,
                                 const gp_XYZ& axis, const gp_XYZ& xHint);

// D-fc-4.9 P2 option (i): the over-determined family certificate of a Torus-mode
// ring chain. One least-squares axis (4 DOF) over the coaxial donor family's
// certified vertices (one radius per donor, projected out) PLUS the chain's
// vertices (one profile circle Rmaj, z0, Rmin about that axis), from the seed
// frame. The chain is accepted only when every donor vertex stays within tau of
// its donor's radius about the joint axis AND certifyTorus itself passes on the
// chain's profile about it (every chain vertex within tau of the torus). Never a
// band-only refit: the donors are always in the fit and in the certificate.
struct JointTorusCert {
    JointAxisFit family;  // the joint axis; radius / maxDev per donor
    TorusCert chain;      // certifyTorus on the chain's profile about the joint axis
    double tilt = 0;      // angle between the seed axis and the joint axis (rad)
    bool fitted = false;  // the joint least squares produced a finite axis
    bool ok = false;      // jointAxisCertifies(family, tau) && chain.ok
};
JointTorusCert jointTorusCertify(const std::vector<std::vector<gp_XYZ>>& donors,
                                 const std::vector<gp_XYZ>& chain, const gp_XYZ& originSeed,
                                 const gp_XYZ& axisSeed, const gp_XYZ& xHint, double tau,
                                 int maxIters);

// D-fc-4.11 (F1 refined): a Torus-mode ring chain's candidate set fails the
// GEOMETRIC certificate. When the joint family certificate ran on a finite axis
// it decides: on a chain that certifyTorus finds determinate about the joint
// axis, a donor vertex off tau of its donor's radius or a geometric refusal of
// certifyTorus. Otherwise the family-axis certificate decides. A chain refused
// only for over-determination (or that certifies) does not fail it.
bool torusChainFailsGeometric(const TorusCert& onFamily, const JointTorusCert* joint, double tau);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_RELATIONAL_TORUS_HPP
