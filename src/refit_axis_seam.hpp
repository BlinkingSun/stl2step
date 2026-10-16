// One seam azimuth per revolution axis (D-fc-4.8). A closed torus's coaxial
// family -- closed-360 cylinders, chamfer cones, closed tori -- shares one frame
// XDirection, and every closed rim of the family starts at the exact circle
// point on that half-plane, so every seam of the family is one generator /
// meridian through exact rim points. Every threshold is tau = 2q.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_AXIS_SEAM_HPP
#define STL2STEP_REFIT_AXIS_SEAM_HPP

#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <Precision.hxx>

#include <cmath>
#include <vector>

namespace stl2step {
namespace refit {

// An axis is on the family axis within tau when its line passes within tau of
// the family line and its direction error, acting on `lever` (the farthest of
// its own rim points from the family's axis point), moves no point by more than
// tau. Antiparallel directions are the same axis.
inline bool axisOnFamilyAxis(const gp_Ax1& a, const gp_Ax1& fam, double lever, double tau) {
    if (!(tau > 0.0)) return false;
    if (gp_Lin(fam).Distance(a.Location()) > tau) return false;
    const double sinA = gp_Vec(a.Direction()).Crossed(gp_Vec(fam.Direction())).Magnitude();
    return sinA * lever <= tau;
}

// A rim circle lies on the family axis when its centre is within tau of the
// family line and its plane's tilt moves no rim point by more than tau.
inline bool circleOnFamilyAxis(const gp_Circ& c, const gp_Ax1& fam, double tau) {
    return axisOnFamilyAxis(c.Axis(), fam, c.Radius(), tau);
}

// `ax` with its XDirection turned onto refX (refX projected normal to ax's own
// axis). Handedness is kept. False when refX is parallel to the axis.
inline bool frameOnRefX(const gp_Ax3& ax, const gp_Dir& refX, gp_Ax3& out) {
    const gp_Vec n(ax.Direction());
    const gp_Vec x = gp_Vec(refX) - n * gp_Vec(refX).Dot(n);
    if (x.Magnitude() <= Precision::Confusion()) return false;
    out = ax;
    out.SetXDirection(gp_Dir(x));
    return true;
}

// A frame is turned when its own u = 0 stands more than tau off the family
// half-plane at `arm` (the largest radius its rims reach).
inline bool frameNeedsTurn(const gp_Ax3& ax, const gp_Ax3& turned, double arm, double tau) {
    return tau > 0.0 && ax.XDirection().Angle(turned.XDirection()) * arm > tau;
}

// The family half-plane on a rim circle: the unit radial x (refX projected into
// the circle's plane), y = n x x, and the exact circle point P = centre + R x.
struct SeamHalfPlane {
    gp_Vec x, y;
    gp_Pnt p;
};
inline bool seamHalfPlaneOn(const gp_Circ& c, const gp_Dir& refX, SeamHalfPlane& out) {
    const gp_Vec n(c.Axis().Direction());
    gp_Vec x = gp_Vec(refX) - n * gp_Vec(refX).Dot(n);
    if (x.Magnitude() <= Precision::Confusion()) return false;
    x.Normalize();
    out.x = x;
    out.y = n.Crossed(x);
    out.p = c.Location().Translated(x * c.Radius());
    return true;
}

// The rim point nearest the half-plane in azimuth (first on ties); -1 if none.
inline int nearestToSeamHalfPlane(const gp_Circ& c, const SeamHalfPlane& h,
                                  const std::vector<gp_Pnt>& pts) {
    int best = -1;
    double bestA = 1e300;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const gp_Vec d(c.Location(), pts[i]);
        const double a = std::fabs(std::atan2(d.Dot(h.y), d.Dot(h.x)));
        if (a < bestA) {
            bestA = a;
            best = (int)i;
        }
    }
    return best;
}

// How far a rim's seam vertex stands off the family half-plane before any
// construction: the chord, at its own radius, from the azimuth-nearest rim
// point to the half-plane through the family axis.
inline double rimSeamOffset(const gp_Ax1& axis, const gp_Dir& refX, const std::vector<gp_Pnt>& pts) {
    const gp_Circ unit(gp_Ax2(axis.Location(), axis.Direction()), 1.0);
    SeamHalfPlane h;
    if (!seamHalfPlaneOn(unit, refX, h)) return 0.0;
    const int k = nearestToSeamHalfPlane(unit, h, pts);
    if (k < 0) return 0.0;
    const gp_Vec d(axis.Location(), pts[(std::size_t)k]);
    const double x = d.Dot(h.x), y = d.Dot(h.y);
    return 2.0 * std::hypot(x, y) * std::sin(0.5 * std::fabs(std::atan2(y, x)));
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_AXIS_SEAM_HPP
