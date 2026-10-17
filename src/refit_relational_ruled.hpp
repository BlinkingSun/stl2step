// RULED QUAD between certified rails (D-fc-5.9). Stage body in
// refit_relational_ruled.cpp. A chart whose boundary is four certified
// straight lines -- two rails (wall|wall creases) and two caps lying in
// anchorCertified planes -- is the degree-(1,1) doubly ruled patch they span.
// Its corners are rail ∩ cap plane: zero free parameters. Every threshold is
// tau = 2q or a count of parameters.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_RELATIONAL_RULED_HPP
#define STL2STEP_REFIT_RELATIONAL_RULED_HPP

#include "refit_relational.hpp"

#include <gp_XYZ.hxx>
#include <array>
#include <vector>

namespace stl2step {
namespace refit {

// D-fc-5.9. The patch is fixed by its rails and cap planes; nothing is fitted
// on the chart itself. Its evidence is the chart's interior vertices, so the
// chart is under-determined with fewer than minStations(kParamsRuled) of them.
constexpr int kParamsRuled = 0;

// Corner poles {P00, P10, P01, P11} of S(u,v), u, v in [0, 1].
using RuledPoles = std::array<gp_XYZ, 4>;

// Line-first certificate of one boundary chain (D-fc-5.9 step 2, D-fc-1.6.2):
// the least-squares line of the chain's points; stations are the points
// projected on it, merged at tau. ok iff stations >= minStations(kParamsLine)
// and every point is within tau of the line. Points are accumulated in the
// order given (callers pass them ascending by local vertex id, so the two
// walls of one rail fit the identical line).
struct ChainLineCert {
    gp_XYZ p, d;          // a point on the line (the centroid) and its unit direction
    double maxDev = 0;    // largest point-to-line distance (mm)
    int stations = 0;
    bool ok = false;
    const char* why = "underdetermined";
};
ChainLineCert certifyChainLine(const std::vector<gp_XYZ>& pts, double tau);

// A cap chain lies in a plane anchor (D-fc-1.4): an Unclaimed provisional A3
// would commit whose own maxVertexDev is within tau.
bool ruledCapAnchor(bool unclaimed, bool a3Commits, double maxVertexDev, double tau);

// Corner = rail line ∩ cap plane (closed form). False when the rail is parallel
// to the plane (no point).
bool railCapCorner(const gp_XYZ& railP, const gp_XYZ& railD, const gp_XYZ& planeP,
                   const gp_XYZ& planeN, gp_XYZ& out);

// Line first on the surface: the largest distance of the four corners from
// their least-squares plane. Within tau the chart is a plane (A3 owns it).
double ruledCornerPlanarity(const RuledPoles& P);

gp_XYZ ruledPoint(const RuledPoles& P, double u, double v);
// S_u x S_v at (u, v); not normalised.
gp_XYZ ruledNormal(const RuledPoles& P, double u, double v);

// Distance from x to S. Gauss-Newton on |S - x|^2 from the closed-form inverse
// of the bilinear map in the patch's centre tangent plane, at most maxIters
// steps. With clampDomain the
// distance is measured to the face (u, v clamped to [0, 1]); otherwise to the
// unbounded surface. The value is the distance to an actual surface point, so
// it never understates the true distance (an unconverged step refuses, it
// never admits).
double ruledDistance(const RuledPoles& P, const gp_XYZ& x, int maxIters, bool clampDomain,
                     double* uOut = nullptr, double* vOut = nullptr);

// Vertex certificate (D-fc-1.6.5): every point within tau of the face. Max,
// not RMS.
struct RuledVertexCert {
    double maxDev = 0, rmsDev = 0;
    bool ok = false;
};
RuledVertexCert certifyRuledVertices(const RuledPoles& P, const std::vector<gp_XYZ>& pts,
                                     double tau, int maxIters);

// Same surface (D-130-14): every point of an adjacent chart within tau of the
// unbounded S. One surface is one face, so such a chart joins the quad.
bool ruledSameSurface(const RuledPoles& P, const std::vector<gp_XYZ>& pts, double tau,
                      int maxIters);

// dVolPredicted of one mesh triangle (a, b, c in its outward winding): the
// surface-minus-chord integral (S(proj p) - p) . n_t over the triangle by the
// 7-point (degree-3) rule, signed per I9. sagOut receives the largest |S - p|
// at the rule's points.
double dVolRuledTriangle(const RuledPoles& P, const gp_XYZ& a, const gp_XYZ& b,
                         const gp_XYZ& c, int maxIters, double* sagOut = nullptr);

bool claimRuledQuads(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_RELATIONAL_RULED_HPP
