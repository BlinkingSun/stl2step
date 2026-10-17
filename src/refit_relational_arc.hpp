// Arcs of an extrusion by relation (D-fc-5.1, D-fc-5.2). Bodies in
// refit_relational_revolve.cpp (coaxial family) and refit_relational_extrude.cpp
// (section fillets). Every threshold is tau = 2q or a count of parameters.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_ARC_HPP
#define STL2STEP_REFIT_RELATIONAL_ARC_HPP

#include "refit_relational_torus.hpp"

#include <gp_XYZ.hxx>
#include <vector>

namespace stl2step {
namespace refit {

// D-fc-5.1. The certificate of a band of rings at one rho about a coaxial
// family's axis. In order: under-determined (fewer than 2 axial levels, or
// fewer than minStations(kParamsCircle) azimuth stations: the band alone must
// over-determine a free section circle); line (the band's section points lie
// on one line within tau: a planar strip); then the joint family certificate
// of D-fc-4.9 option (i) in Line mode -- the family axis is kept when every
// family vertex and every band vertex is within tau of its own radius about
// it, else one least-squares axis over the family's certified vertices plus
// the band's vertices decides (vertex-off when it does not certify).
// fit.radius / fit.maxDev hold one entry per family group, the band last.
struct CoaxialCert {
    JointAxisFit fit;
    bool refit = false;  // the band certified on a new least-squares axis
    bool ok = false;
    const char* why = "underdetermined";
    int stations = 0, levels = 0;
    double lineRes = 0;  // section line residual of the band (mm)
};
CoaxialCert certifyCoaxialBand(const std::vector<std::vector<gp_XYZ>>& family,
                               const std::vector<gp_XYZ>& band, const gp_XYZ& origin,
                               const gp_XYZ& axis, double tau, int maxIters);

// D-fc-5.1 / D-fc-1.6.11. A certified coaxial band whose radius is within tau
// of a family member's radius about the same axis is another piece of that
// member's surface: one surface is one face (D-130-14), so it is refused.
bool coaxialSameRadius(const JointAxisFit& fit, std::size_t nFamily, double tau);

// Largest distance of a 2-D point set from its own least-squares line: the
// section line test (line first, fewest parameters).
double sectionLineResid(const std::vector<gp_XYZ>& pts, const gp_XYZ& axis);

// D-fc-5.2. In the cross-section of an extrusion an anchor is a line (a
// certified plane parallel to the extrusion) or a circle (a certified cylinder
// parallel to it). A section fillet is a circle of radius r tangent to both on
// the sides the data give: line n.p - d = side * r; circle |p - c| = R + side * r.
struct SectionAnchor {
    bool circle = false;
    double nx = 0, ny = 0, d = 0;  // line: n.x = d, |n| = 1
    double cx = 0, cy = 0, R = 0;  // circle
    int side = 1;
    double tol = 0;  // anchorTol (D-fc-1.4)
};

// The centre of the fillet of radius r tangent to A and B: the intersection
// of their r-offsets nearest (hx, hy). False when the offsets do not meet.
// Concentric circles have no point intersection (sectionConcentric).
bool sectionFilletCentre(const SectionAnchor& A, const SectionAnchor& B, double r, double hx,
                         double hy, double& px, double& py);

// Two circle anchors whose centres are within the sum of their tolerances are
// concentric: their r-offsets meet only when r is half their gap, and the free
// parameter is then the centre's azimuth on the mid circle (a hem; D-fc-1.5 O2).
bool sectionConcentric(const SectionAnchor& A, const SectionAnchor& B);

// D-fc-5.2 contact slack: the angular uncertainty of a contact point on anchor
// X seen from the fillet centre, tanDistTol(anchorTol(X), tau) / r.
inline double sectionContactSlack(double tolX, double tau, double r) {
    return r > 0.0 ? (tolX + tau) / r : 0.0;
}

// D-fc-1.4 + D-fc-4.1: the plane anchor of an Unclaimed provisional A3 would
// commit -- certified on its own maxVertexDev, or else on the tau-preimage of
// its maximum-area facet's plane refit to a fixed point (the REV ORCH-4 peel,
// whose body is then the anchor). The predicate of claimRelExtrude's anchors.
struct PlaneAnchor {
    gp_XYZ n, p0;
    double dev = 0, tol = 0;
    bool certified = false;
    std::vector<int> body;
};
bool planeAnchorOf(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                   const Provisional& pr, PlaneAnchor& out);

// A certified cylinder that may anchor a section fillet (a Line-mode donor
// or a coaxial band, D-fc-1.4).
struct SectionCircle {
    gp_XYZ origin, axis;
    double radius = 0, tol = 0, length = 0;
    int acc = -1;  // work.accepted index
};

// D-fc-5.2. Runs after the coaxial family claims (it needs their circles).
// Hypotheses come from accepted non-closed CylGrow regions (D-fc-4.6 F2 atoms)
// tangent to two anchors, at least one a circle (plane|plane is the extrude
// pair walk). The fillet is fit with one free parameter, its claimed set is the
// complete tau-preimage of the unclaimed candidates and atoms in the slab
// between its contact generators (atoms join whole, over-determination counted
// leave-one-atom-out), refit to a fixed point and certified in the D-fc-3.2.6
// order. Never throws.
bool claimSectionFillets(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                         SegmentWork& work, const std::vector<SectionCircle>& circles);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_RELATIONAL_ARC_HPP
