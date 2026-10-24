// D-fc-7.3 RQ-HUB-5 (F-E). Same-surface re-host onto a certified coaxial
// family member. Inline; no tolerance of its own. tau = tauOf = 2q.
// Not a merge: the caller writes the member's (axis, R) onto the region
// and leaves its triangles where they are (one surface per family, D-130-14,
// faces stay separate). No subset: one vertex outside tau keeps the region's
// own fit.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_COAX_REHOST_HPP
#define STL2STEP_REFIT_COAX_REHOST_HPP

#include "refit_coax_closed.hpp"

#include <vector>

#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {

// An accepted region that coaxialFamilies grouped, certified on the surface
// it ships (anchorCertified: maxVertexDev <= tau). acc is its accepted index
// so the choice among several certifying members is deterministic.
struct CoaxRehostMember {
    int acc = -1;
    gp_XYZ origin;
    gp_XYZ axis;
    double radius = 0;
    bool certified = false;
};

struct CoaxRehostDecision {
    bool take = false;
    int acc = -1;
    gp_XYZ origin;
    gp_XYZ axis;
    double radius = 0;
};

// Every vertex lies within tau of the cylinder (axis through o along axisIn,
// radius R). tau <= 0, a non-positive radius, or an empty set admits nothing.
inline bool coaxSameSurfaceRehost(const std::vector<gp_XYZ>& verts, const gp_XYZ& o,
                                  const gp_XYZ& axisIn, double R, double tau) {
    if (!(tau > 0.0) || !(R > tau) || verts.empty()) return false;
    const double am = axisIn.Modulus();
    if (!(am > 0.0)) return false;
    const gp_XYZ a = axisIn / am;
    for (const gp_XYZ& p : verts)
        if (!(coaxRadialDev(p, o, a, R) <= tau)) return false;
    return true;
}

// A certified member keeps its own fit when no lower-index certified member
// already carries every one of its vertices. The higher member is then not a
// surface to copy: both would take each other and the faces would swap.
inline bool coaxRehostStable(const std::vector<gp_XYZ>& ownVerts, int acc,
                             const std::vector<CoaxRehostMember>& members, double tau) {
    if (!(tau > 0.0)) return false;
    for (const CoaxRehostMember& o : members) {
        if (!o.certified || o.acc >= acc) continue;
        if (coaxSameSurfaceRehost(ownVerts, o.origin, o.axis, o.radius, tau)) return false;
    }
    return true;
}

// Lowest accepted index among certified members whose cylinder every vertex
// certifies on. None => take stays false and the caller keeps the own fit.
// Pass only stable members (coaxRehostStable): a member that itself lies on a
// lower member is not a surface of its own.
inline CoaxRehostDecision coaxRehostChoose(const std::vector<gp_XYZ>& verts,
                                           const std::vector<CoaxRehostMember>& members,
                                           double tau) {
    CoaxRehostDecision d;
    if (!(tau > 0.0)) return d;
    for (const CoaxRehostMember& m : members) {
        if (!m.certified || !(m.radius > tau)) continue;
        if (d.take && m.acc >= d.acc) continue;
        if (!coaxSameSurfaceRehost(verts, m.origin, m.axis, m.radius, tau)) continue;
        d.take = true;
        d.acc = m.acc;
        d.origin = m.origin;
        d.axis = m.axis;
        d.radius = m.radius;
    }
    return d;
}

// The stored cylinder is already that surface: axis line within tau, direction
// within tanAngTol over the region's own axial length, radius within tau.
// Rewriting it would only slide the placement origin along the same line.
inline bool coaxRehostAlready(const gp_XYZ& oR, const gp_XYZ& aR, double radiusR, double length,
                              const gp_XYZ& oM, const gp_XYZ& aM, double radiusM, double tau) {
    if (!(tau > 0.0)) return false;
    const double mr = aR.Modulus(), mm = aM.Modulus();
    if (!(mr > 0.0) || !(mm > 0.0)) return false;
    const gp_XYZ ur = aR / mr, um = aM / mm;
    const double L = std::max(length, tau);
    if (ur.Crossed(um).Modulus() > tanAngTol(tau, tau, L)) return false;
    const gp_XYZ dlt = oR - oM;
    if ((dlt - um * dlt.Dot(um)).Modulus() > tau) return false;
    return std::fabs(radiusR - radiusM) <= tau;
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_COAX_REHOST_HPP
