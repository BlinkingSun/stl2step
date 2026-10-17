// Sphere corners by relation (D-fc-5.11 F8). Body in refit_relational_sphere.cpp.
// Every threshold is tau = 2q or a count of parameters (D-fc-1.6).
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_SPHERE_HPP
#define STL2STEP_REFIT_RELATIONAL_SPHERE_HPP

#include <gp_XYZ.hxx>

#include <array>
#include <cmath>
#include <vector>

namespace stl2step {
namespace refit {

// D-fc-5.11 F8. A sphere has four parameters (centre and radius). At a corner
// the three rounds fix all four, so the sphere has none left to fit.
constexpr int kParamsSphere = 4;

// The spine of a certified round: the locus of its section-circle centres. A
// cylinder's spine is its axis line; a torus blend's is its centre circle
// (radius R_major in the plane through `o` normal to `a`). `r` is the round's
// own radius (cylinder radius, torus R_minor).
struct RoundSpine {
    bool circle = false;
    gp_XYZ o;
    gp_XYZ a;       // unit
    double R = 0;   // circle only
    double r = 0;
};

// Distance from x to the spine.
double spineDistance(const RoundSpine& s, const gp_XYZ& x);

// Unit tangent of the spine at the spine point nearest x. False when it is
// undefined (x on a circle spine's axis).
bool spineTangentAt(const RoundSpine& s, const gp_XYZ& x, gp_XYZ& t);

// Least-squares point of the LINE spines in `spines`: the point minimising the
// sum of squared distances to the lines. False with fewer than two lines or
// when the lines are parallel (the normal matrix is singular).
bool lineSpinesPoint(const std::vector<RoundSpine>& spines, gp_XYZ& c);

// D-fc-5.11 F8. The sphere of three rounds meeting at a corner. The centre is
// the point their spines share -- the least-squares point of the line spines,
// which needs at least two of them -- and the radius is the rounds' mean
// radius. It certifies only when the three radii agree within tau and the
// centre lies within tau of every spine: zero free parameters, nothing fitted
// to the corner's own triangles. why: "accept" | "underdetermined" (fewer than
// two line spines) | "not-blend" (radii or spines disagree).
struct SphereCornerFit {
    bool ok = false;
    const char* why = "underdetermined";
    gp_XYZ c;
    double r = 0;
    double spineDev = 0;      // max distance of c from the three spines
    double radiusSpread = 0;  // max - min of the three radii
};
SphereCornerFit sphereOfRounds(const std::array<RoundSpine, 3>& rounds, double tau);

// Largest distance of the points from their own least-squares plane (0 for
// fewer than four points).
double planeResidual(const std::vector<gp_XYZ>& pts);

// Plane first (fewest parameters, D-fc-1.6(3)): true when every point is within
// tau of the points' own least-squares plane.
bool pointsPlanar(const std::vector<gp_XYZ>& pts, double tau);

// I9 for a sphere patch: the volume between the spherical triangle over
// (p, q, s) -- three points on the sphere (centre, r) -- and the flat triangle:
// the spherical sector r^3 * Omega / 3 minus the tetrahedron (centre, p, q, s).
// Unsigned.
double sphereCapVolume(const gp_XYZ& centre, double r, const gp_XYZ& p, const gp_XYZ& q,
                       const gp_XYZ& s);

// D-fc-1.4 / D-fc-3.2.6-7 read on a cylinder round's own triangles: the round's
// surface is certified on its tau-preimage I (every triangle whose three
// vertices are within tau of it) when I over-determines the section circle
// (>= minStations(kParamsCircle) azimuth stations, >= 2 axial levels), and every
// other triangle of the round is a chord-bounded fill of I (>= 1 vertex within
// tau, every vertex within s = I's own largest edge-midpoint deviation, inside
// I's axial extent widened by tau). A round with no off-tau triangle is I.
struct RoundCert {
    bool ok = false;
    int nI = 0, nCover = 0, nOther = 0, stations = 0, levels = 0;
    double certDev = 0, coverDev = 0, s = 0;
};
RoundCert certifyCylinderRound(const RoundSpine& spine,
                               const std::vector<std::array<gp_XYZ, 3>>& tris, double tau);

// The corner sector. A round meeting the corner ends on the plane through the
// centre normal to its spine tangent t there; its body lies on the side
// side * (x - c).t >= 0. The sphere patch lies on the other side of all three
// planes. Each test carries tau slack (a vertex on a contact circle lies on
// its plane).
inline bool onRoundSide(const gp_XYZ& x, const gp_XYZ& c, const gp_XYZ& t, double side,
                        double tau) {
    return side * (x - c).Dot(t) >= -tau;
}
inline bool onPatchSide(const gp_XYZ& x, const gp_XYZ& c, const gp_XYZ& t, double side,
                        double tau) {
    return side * (x - c).Dot(t) <= tau;
}

// D-fc-5.11 F8 build. A sphere corner (centre c) meets a cylinder round (axis
// point o, unit axis a, radius R) along the cylinder's v-iso through c -- a
// great circle of the sphere when the radii agree. It exists only when c lies
// within tol of the axis and |R - r| <= tol; v is the iso's axial parameter.
bool sphereCylinderContact(const gp_XYZ& c, double r, const gp_XYZ& o, const gp_XYZ& a,
                           double R, double tol, double& v);

// D-fc-5.11 F8 build. A circle C(t) = C0 + rho (cos t Xc + sin t Yc) against a
// sphere iso-curve traced by a 2-D line pcurve of unit speed in the circle's
// own parameter: a parallel (u = u0 + s t, v = v0) or a meridian (u = u0,
// v = v0 + s t). Their difference is P + Q cos t + R sin t, so
// sup_t |C(t) - S(t)| <= |P| + sqrt(lambda_max([Q.Q, Q.R; Q.R, R.R])) -- a
// closed-form bound over the whole circle. `meridian` selects the iso kind;
// (sx, sy, sz) is the sphere's frame (sz its axis).
double sphereIsoBound(const gp_XYZ& centre, double r, const gp_XYZ& sx, const gp_XYZ& sy,
                      const gp_XYZ& sz, const gp_XYZ& c0, double rho, const gp_XYZ& xc,
                      const gp_XYZ& yc, bool meridian, double u0, double v0, double s);

// The iso line of a sphere nearest a circle: the parallel when the circle's
// plane normal is nearer the sphere axis than the plane, else the meridian
// through the circle's point at tMid. Returns its (u0, v0, s, meridian) and
// sphereIsoBound for it; u0 is put in [0, 2 pi).
struct SphereIsoLine {
    bool meridian = false;
    double u0 = 0, v0 = 0, s = 1;
    double bound = HUGE_VAL;
};
SphereIsoLine sphereIsoLineOf(const gp_XYZ& centre, double r, const gp_XYZ& sx, const gp_XYZ& sy,
                              const gp_XYZ& sz, const gp_XYZ& c0, double rho, const gp_XYZ& xc,
                              const gp_XYZ& yc, double tMid);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_RELATIONAL_SPHERE_HPP
