// One seam azimuth per revolution axis (D-fc-4.8). A closed torus's coaxial
// family -- closed-360 cylinders, chamfer cones, closed tori -- shares one frame
// XDirection, and every closed rim of the family starts at the exact circle
// point on that half-plane, so every seam of the family is one generator /
// meridian through exact rim points. Every threshold is tau = 2q.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_AXIS_SEAM_HPP
#define STL2STEP_REFIT_AXIS_SEAM_HPP

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec2d.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <Precision.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
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

// The builder builds closed-360 cylinders first, then regions in id order. A rim
// shared with a face built before the torus is that face's own circle, minted at
// that face's vertex; the closed-torus path can re-mint only a rim built after it.
inline bool rimBuiltBeforeTorus(bool partnerClosedCylinder, int partnerId, int torusId) {
    return partnerClosedCylinder || partnerId < torusId;
}

// A rim curve runs with the surface's u when its tangent at parameter f has no
// negative component along the surface's du at (u0, v). A coaxial wall fitted with
// an antiparallel axis gives its rims the opposite sense.
inline bool rimRunsWithSurfaceU(const Handle(Geom_Curve)& c, double f,
                                const Handle(Geom_Surface)& s, double u0, double v) {
    if (c.IsNull() || s.IsNull()) return true;
    gp_Pnt P, S;
    gp_Vec T, Du, Dv;
    c->D1(f, P, T);
    s->D1(u0, v, S, Du, Dv);
    return !(T.Dot(Du) < 0.0);
}

// Heal-not-revert for a closed wall whose split cap meets an unbuilt partner: only a
// wall certified at tau (every vertex within tau, D-fc-1.4) that still shares one
// exact rim with an analytic neighbour ships with its mesh-polyline rim.
inline bool splitCapHealApplies(double maxVertexDev, double tau, bool oneExactRim) {
    return tau > 0.0 && maxVertexDev <= tau && oneExactRim;
}

// D-fc-5.7 O-1: a closed wall whose rim crosses its neighbours along mesh edges -- a
// bore leaving through a round, its rim a non-planar mesh polyline -- ships that rim
// as the polyline beside one exact rim. Scope: the wall's own vertices leave tau only
// on its rims (every interior vertex within tau, D-fc-1.4; a rim vertex beyond it is
// where the mesh leaves the surface, the chord-bounded cover of D-fc-4.5). A wall
// certified on every vertex keeps the split-cap heal above.
inline bool crossingCapApplies(double interiorDev, double maxVertexDev, double tau,
                               bool oneExactRim, int crossingChains) {
    return tau > 0.0 && interiorDev <= tau && maxVertexDev > tau && oneExactRim &&
           crossingChains > 0;
}

// Where a generator meets the plane of a cap circle or ellipse; the caller takes the
// curve point at that point's parameter (exact on the curve). False when the line
// runs parallel to the plane.
inline bool capPointOnGenerator(const gp_Ax2& capPlane, const gp_Lin& gen, gp_Pnt& x) {
    const gp_XYZ n = capPlane.Direction().XYZ();
    const double den = n.Dot(gen.Direction().XYZ());
    if (std::fabs(den) <= gp::Resolution()) return false;
    const double t = n.Dot(capPlane.Location().XYZ() - gen.Location().XYZ()) / den;
    x = gp_Pnt(gen.Location().XYZ() + t * gen.Direction().XYZ());
    return true;
}

// D-fc-4.6 F9 over the whole rim. A closed torus's rim is one 2π circle, so every
// chain of its cap loop must border a face that can carry that circle. A chain whose
// partner cannot (facets, an island, an exploded or E'-filled region) leaves the
// circle beside a mesh polyline that no later round can share. Only the loop's first
// chain can hand the torus a rim edge already built (frontReused).
template <class PartnerBuildable>
inline bool capLoopUnshareable(const std::vector<int>& chainIdx, bool frontReused,
                               PartnerBuildable&& partnerBuildable) {
    for (std::size_t k = 0; k < chainIdx.size(); ++k) {
        const int ci = chainIdx[k];
        if (ci < 0 || (k == 0 && frontReused)) continue;
        if (!partnerBuildable(ci)) return true;
    }
    return false;
}

// U-3 per torus (D-fc-4.11 (a)). A build round in which a torus failed ships nothing:
// the builder returns the tori it names and the driver withholds them. That set is
// decided by each torus's own build and by the cap-partner rule over the round's
// failures beside a closed torus; namesTorus marks those regions in build order.
// torusNamingEnd is one past the last of them (0 when there is none).
inline std::size_t torusNamingEnd(const std::vector<char>& namesTorus) {
    for (std::size_t k = namesTorus.size(); k > 0; --k)
        if (namesTorus[k - 1]) return k;
    return 0;
}

// Once a torus failed, the round is decided when every region that can name a torus
// has been built: the regions after them cannot change the named set.
inline bool torusRoundDecided(bool torusFailed, std::size_t built, std::size_t namingEnd) {
    return torusFailed && built >= namingEnd;
}

// The cap-partner rule ends a round with no failed torus the same way: a failed
// region beside a closed torus built whole (builtClosedTorusBeside) leaves that
// torus's 2π rim beside a mesh polyline, so the round ships nothing and names the
// torus. Every region beside a closed torus names one, so once the naming regions
// are built the regions after them cannot change the outcome (torusRoundDecided).
template <class BuiltClosedTorusBeside>
inline bool closedTorusPartnerFailed(const std::vector<int>& failedIds,
                                     BuiltClosedTorusBeside&& builtClosedTorusBeside) {
    for (int id : failedIds)
        if (builtClosedTorusBeside(id)) return true;
    return false;
}

// A closed rim chain read from its terminal (D-130-16 two-arc rim collapse). The
// chain's terminal is the seam column of the wall it bounds, which need not be
// its first mesh vertex; the loop is rotated to start there and closed by
// repeating it, so the split vertex's index divides the loop into the two arcs
// the collapse builds. A terminal not on the loop leaves the order as it is.
inline std::vector<int> closedLoopFromTerminal(std::vector<int> loop, int terminal) {
    if (loop.empty()) return loop;
    const auto at = std::find(loop.begin(), loop.end(), terminal);
    if (at != loop.end()) std::rotate(loop.begin(), at, loop.end());
    loop.push_back(loop.front());
    return loop;
}

// A closed wall with windows (D-fc-4.4 ring closure). Each boundary loop of a
// claimed triangle set winds round the axis a whole number of times. The set is
// a closed wall when exactly two loops wind once, in opposite senses (its two
// rims), and every other loop winds zero times (a window: an inner wire of the
// face). Topology only; no tolerance.
inline bool closedWallWithWindows(const std::vector<int>& windings) {
    int up = 0, down = 0;
    for (int w : windings) {
        if (w == 1) ++up;
        else if (w == -1) ++down;
        else if (w != 0) return false;
    }
    return up == 1 && down == 1;
}

// The arc a window covers is the complement of the largest gap of its own
// azimuths; returns that gap as [lo, hi) with hi - lo in (0, 2pi]. False for
// fewer than two azimuths.
inline bool windowGap(std::vector<double> us, double& lo, double& hi) {
    if (us.size() < 2) return false;
    const double twoPi = 2.0 * std::acos(-1.0);
    for (double& u : us) u = std::remainder(u, twoPi);
    std::sort(us.begin(), us.end());
    lo = us.back();
    hi = us.front() + twoPi;
    for (std::size_t i = 1; i < us.size(); ++i)
        if (us[i] - us[i - 1] > hi - lo) {
            lo = us[i - 1];
            hi = us[i];
        }
    return hi > lo;
}

// The seam of a closed wall is ONE generator of the face and may not enter a
// window (OCCT: every inner wire inside the chart the seam names). It stands in
// the middle of the widest gap of the windows' pooled azimuths -- the builder's
// own placement (trySeamed360) -- and it is clear when it lies strictly inside
// every window's own gap. With no window every generator is clear (seamU = 0).
inline bool seamClearOfWindows(const std::vector<std::vector<double>>& windowU, double& seamU) {
    const double twoPi = 2.0 * std::acos(-1.0);
    std::vector<double> pooled;
    for (const std::vector<double>& w : windowU) pooled.insert(pooled.end(), w.begin(), w.end());
    seamU = 0.0;
    if (pooled.empty()) return true;
    double lo = 0, hi = 0;
    if (!windowGap(pooled, lo, hi)) return false;
    seamU = std::remainder(0.5 * (lo + hi), twoPi);
    for (const std::vector<double>& w : windowU) {
        double wl = 0, wh = 0;
        if (!windowGap(w, wl, wh)) return false;
        double u = seamU;
        while (u < wl) u += twoPi;
        while (u >= wl + twoPi) u -= twoPi;
        if (!(u > wl && u < wh)) return false;
    }
    return true;
}

// D-fc-5.7 O-4. A split-cap closed wall's plain wire runs seam up, the high cap path,
// seam down, the low cap path reversed. It closes in (u, v) only when the two paths
// sweep the same way round the axis. Each cap loop is oriented by its own chains, so
// the two can arrive sweeping opposite ways. That wire then jumps a whole period at
// each seam copy, and the wall takes the composite UV-sense and seam-pair binding.
// A path that reads no sweep decides nothing.
inline bool capPathsRunOpposite(double sweepH, double sweepL) {
    return sweepH * sweepL < 0.0;
}

// D-S3-171 case (c), the sampled branch: the signed area (plane UV, about the plane's
// normal) of an edge's curve sampled densely along its own parameter and closed back
// to its start. The edge's traversal is that run read under its orientation flag
// (D-fc-5.13: a plane whose inner loop is one closed ellipse handed REVERSED).
inline double edgeSampledSignedAreaUV(const TopoDS_Edge& e, const gp_Pln& pln) {
    if (e.IsNull()) return 0.0;
    try {
        BRepAdaptor_Curve ac(e);
        const Standard_Real f = ac.FirstParameter();
        const Standard_Real l = ac.LastParameter();
        if (l - f <= Precision::PConfusion()) return 0.0;
        std::vector<gp_Pnt2d> uv;
        uv.reserve(32);
        for (int i = 0; i < 32; i++) {
            const Standard_Real t = f + (l - f) * (Standard_Real)i / 32.0;
            Standard_Real u = 0, vv = 0;
            ElSLib::Parameters(pln, ac.Value(t), u, vv);
            const gp_Pnt2d q(u, vv);
            if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
            uv.push_back(q);
        }
        if (uv.size() < 3) return 0.0;
        double a = 0.0;
        for (std::size_t i = 0; i < uv.size(); i++) {
            const gp_Pnt2d& p0 = uv[i];
            const gp_Pnt2d& p1 = uv[(i + 1) % uv.size()];
            a += p0.X() * p1.Y() - p1.X() * p0.Y();
        }
        if (e.Orientation() == TopAbs_REVERSED) a = -a;
        return 0.5 * a;
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

inline Handle(Geom_Curve) rimBasisCurve(const Handle(Geom_Curve)& c) {
    const Handle(Geom_TrimmedCurve) t = Handle(Geom_TrimmedCurve)::DownCast(c);
    return (!t.IsNull() && !t->BasisCurve().IsNull()) ? t->BasisCurve() : c;
}

// The point half a turn on from a closed ellipse edge's vertex (its antipode); false
// when the edge is not one full turn of an ellipse.
inline bool closedEllipseAntipode(const TopoDS_Edge& e, gp_Pnt& p) {
    const double kHalfTurn = 3.14159265358979323846;
    if (e.IsNull()) return false;
    try {
        Standard_Real f = 0, l = 0;
        const Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull()) return false;
        const Handle(Geom_Ellipse) ge = Handle(Geom_Ellipse)::DownCast(rimBasisCurve(c));
        if (ge.IsNull() || std::fabs((l - f) - 2.0 * kHalfTurn) > Precision::PConfusion())
            return false;
        p = ge->Value(f + kHalfTurn);
        return true;
    } catch (const Standard_Failure&) {
    }
    return false;
}

// SPEC-common-fc: if q ≤ 0, admit nothing. A non-finite quantum is the same
// refusal. D-fc-5.20: the closed-ellipse split (reuse a mesh vertex, or mint
// the antipode) runs only for a positive finite q. Otherwise the full-turn
// edge stays and no twin is minted.
inline bool ellipseRimSplitAdmitted(double q) {
    return std::isfinite(q) && q > 0.0;
}

// D-130-16: index of the chain vertex nearest the exact rim point `p` and within
// tau, skipping `skip` (a vertex the rim already carries). Tie: the later
// meshVerts entry. -1 when tau is not a positive finite length or no vertex
// lies within tau; the caller then mints `p`. A vertex minted where a mesh
// vertex stands is that vertex's twin.
inline int chainVertexWithinTau(const gp_Pnt& p, const std::vector<TopoDS_Vertex>& verts,
                                const std::vector<int>& meshVerts, const TopoDS_Vertex& skip,
                                double tau) {
    if (!(std::isfinite(tau) && tau > 0.0)) return -1;
    int best = -1;
    double dBest = tau;
    for (int lv : meshVerts) {
        if (lv < 0 || (size_t)lv >= verts.size() || verts[(size_t)lv].IsNull()) continue;
        if (!skip.IsNull() && verts[(size_t)lv].IsSame(skip)) continue;
        const double d = BRep_Tool::Pnt(verts[(size_t)lv]).Distance(p);
        if (d <= dBest) {
            dBest = d;
            best = lv;
        }
    }
    return best;
}

// Index of the chain vertex nearest the antipode and within tau, skipping the
// chain terminal `skip` (the edge's own vertex). Tie: the later meshVerts
// entry. -1 when tau is not a positive finite length, the edge is not one
// full turn of an ellipse, or no other vertex lies within tau.
inline int closedEllipseAntipodePick(const TopoDS_Edge& e,
                                     const std::vector<TopoDS_Vertex>& verts,
                                     const std::vector<int>& meshVerts, int skip, double tau) {
    if (!(std::isfinite(tau) && tau > 0.0)) return -1;
    gp_Pnt pA;
    if (!closedEllipseAntipode(e, pA)) return -1;
    const bool haveSkip = skip >= 0 && (size_t)skip < verts.size();
    return chainVertexWithinTau(pA, verts, meshVerts,
                                haveSkip ? verts[(size_t)skip] : TopoDS_Vertex(), tau);
}

// D-fc-5.13 (the D-fc-5.7 O-1 construction): a closed ellipse edge as two half-turn
// arcs of its own curve, from its vertex to the antipode and back, in the edge's own
// traversal. The antipode is `atAntipode` when the caller has a vertex standing there
// (the chain's own mesh vertex, D-130-16: a second vertex at a mesh vertex's place is
// a twin), else the exact curve point half a turn on, minted. A full-turn ellipse has
// no certified bind class (refit_ellipse_bind's range clause); a half-turn arc has.
// Empty when the edge is not one full turn of an ellipse closed at one vertex.
inline std::vector<TopoDS_Edge> closedEllipseHalves(const TopoDS_Edge& e,
                                                    const TopoDS_Vertex& atAntipode = TopoDS_Vertex()) {
    const double kHalfTurn = 3.14159265358979323846;
    if (e.IsNull()) return {};
    try {
        Standard_Real f = 0, l = 0;
        const Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull()) return {};
        const Handle(Geom_Ellipse) ge = Handle(Geom_Ellipse)::DownCast(rimBasisCurve(c));
        if (ge.IsNull() || std::fabs((l - f) - 2.0 * kHalfTurn) > Precision::PConfusion())
            return {};
        TopoDS_Vertex V, V2;
        TopExp::Vertices(e, V, V2, Standard_False);
        if (V.IsNull() || !V.IsSame(V2) || atAntipode.IsSame(V)) return {};
        V.Orientation(TopAbs_FORWARD);
        TopoDS_Vertex VA;
        double tA = f + kHalfTurn;
        if (!atAntipode.IsNull()) {
            // The vertex's own parameter on the curve, on the edge's turn [f, l).
            VA = TopoDS::Vertex(atAntipode.Oriented(TopAbs_FORWARD));
            tA = ElCLib::Parameter(ge->Elips(), BRep_Tool::Pnt(VA));
            while (tA < f) tA += 2.0 * kHalfTurn;
            while (tA >= l) tA -= 2.0 * kHalfTurn;
            if (!(tA > f + Precision::PConfusion() && tA < l - Precision::PConfusion())) return {};
        } else {
            BRep_Builder().MakeVertex(VA, ge->Value(tA), Precision::Confusion());
        }
        BRepBuilderAPI_MakeEdge m1(ge, V, VA, f, tA);
        BRepBuilderAPI_MakeEdge m2(ge, VA, V, tA, l);
        if (!m1.IsDone() || !m2.IsDone()) return {};
        const TopoDS_Edge a1 = m1.Edge(), a2 = m2.Edge();
        if (e.Orientation() == TopAbs_REVERSED)
            return {TopoDS::Edge(a2.Reversed()), TopoDS::Edge(a1.Reversed())};
        return {a1, a2};
    } catch (const Standard_Failure&) {
    }
    return {};
}

// The two half-turn arcs, or empty. Empty when the quantum does not admit the
// split (D-fc-5.20): closedEllipseHalves is not called, so a null atAntipode
// does not mint. When the quantum admits, a null atAntipode mints the exact
// antipode and a mesh vertex standing there is reused.
inline std::vector<TopoDS_Edge> closedEllipseChainSplit(const TopoDS_Edge& e, double q,
                                                        const TopoDS_Vertex& atAntipode) {
    if (!ellipseRimSplitAdmitted(q)) return {};
    return closedEllipseHalves(e, atAntipode);
}

// One closed analytic rim: a single closed circle or ellipse edge, or the two half-turn
// arcs of one ellipse that closedEllipseHalves makes, joined end to end.
inline bool closedAnalyticRim(const std::vector<TopoDS_Edge>& edges) {
    if (edges.empty() || edges.size() > 2) return false;
    Handle(Geom_Curve) basis;
    for (const TopoDS_Edge& e : edges) {
        if (e.IsNull()) return false;
        Standard_Real f = 0, l = 0;
        const Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
        if (c3.IsNull()) return false;
        const Handle(Geom_Curve) c = rimBasisCurve(c3);
        const bool ellipse = !Handle(Geom_Ellipse)::DownCast(c).IsNull();
        if (!ellipse && Handle(Geom_Circle)::DownCast(c).IsNull()) return false;
        if (edges.size() == 2 && !ellipse) return false;
        if (!basis.IsNull() && basis != c) return false;
        basis = c;
    }
    TopoDS_Vertex a, b;
    TopExp::Vertices(edges.front(), a, b, Standard_True);
    if (a.IsNull() || b.IsNull()) return false;
    if (edges.size() == 1) return a.IsSame(b);
    TopoDS_Vertex c2, d2;
    TopExp::Vertices(edges.back(), c2, d2, Standard_True);
    return !a.IsSame(b) && b.IsSame(c2) && d2.IsSame(a);
}

// D-fc-7.10. The u-chart of one edge on a partial-cylinder rectangular trim:
// the straight UV segment through the edge's own vertices on that trim's
// [uLo, uHi]. Two edges that share a vertex then meet at that vertex's UV.
// A full-turn image is used only when both ends land in the trim. The slack
// is the edge's own tolerance over the cylinder radius (the angle that
// tolerance subtends). Null when the edge is not on this trim — the caller
// keeps its other construction. Installed only by healInvalidPartialTrimChart,
// and only on a face whose BRepCheck already failed.
inline Handle(Geom2d_Curve) partialTrimChartPCurve(const TopoDS_Edge& e, const Handle(Geom_Curve)& c3,
                                                   Standard_Real f, Standard_Real l,
                                                   const Handle(Geom_Surface)& srf) {
    if (e.IsNull() || c3.IsNull() || srf.IsNull()) return {};
    Handle(Geom_RectangularTrimmedSurface) trim =
        Handle(Geom_RectangularTrimmedSurface)::DownCast(srf);
    if (trim.IsNull()) return {};
    Handle(Geom_CylindricalSurface) cyl =
        Handle(Geom_CylindricalSurface)::DownCast(trim->BasisSurface());
    if (cyl.IsNull()) return {};
    Standard_Real uLo = 0, uHi = 0, vLo = 0, vHi = 0;
    trim->Bounds(uLo, uHi, vLo, vHi);
    const double halfTurn = std::acos(-1.0);
    const double twoPi = 2.0 * halfTurn;
    const double span = uHi - uLo;
    if (!(span > Precision::Confusion() && span < twoPi - Precision::Confusion())) return {};
    if (l - f <= Precision::PConfusion()) return {};
    TopoDS_Vertex vf, vl;
    TopExp::Vertices(e, vf, vl, Standard_False);
    if (vf.IsNull() || vl.IsNull() || vf.IsSame(vl)) return {};
    const gp_Ax3 ax = cyl->Position();
    const double radius = cyl->Radius();
    if (!(radius > Precision::Confusion())) return {};
    auto uvOf = [&](const gp_Pnt& p) {
        gp_Vec rho(ax.Location(), p);
        const double v = rho.Dot(ax.Direction());
        gp_Vec rad = rho - gp_Vec(ax.Direction()) * v;
        double u = 0.0;
        if (rad.Magnitude() > Precision::Confusion())
            u = std::atan2(rad.Dot(ax.YDirection()), rad.Dot(ax.XDirection()));
        return gp_Pnt2d(u, v);
    };
    const gp_Pnt2d a0 = uvOf(BRep_Tool::Pnt(vf));
    gp_Pnt2d b0 = uvOf(BRep_Tool::Pnt(vl));
    double du = b0.X() - a0.X();
    Handle(Geom_TrimmedCurve) tr = Handle(Geom_TrimmedCurve)::DownCast(c3);
    Handle(Geom_Curve) src = (!tr.IsNull() && !tr->BasisCurve().IsNull()) ? tr->BasisCurve() : c3;
    Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(src);
    if (!gc.IsNull()) {
        const double sgn = gc->Circ().Axis().Direction().Dot(ax.Direction()) >= 0.0 ? 1.0 : -1.0;
        const double want = sgn * (l - f);
        du += twoPi * std::round((want - du) / twoPi);
    } else {
        while (du > halfTurn) du -= twoPi;
        while (du < -halfTurn) du += twoPi;
    }
    b0.SetX(a0.X() + du);
    const double uSlack =
        std::max(BRep_Tool::Tolerance(e), std::max(BRep_Tool::Tolerance(vf), BRep_Tool::Tolerance(vl))) /
        radius;
    auto inside = [&](double u) { return u >= uLo - uSlack && u <= uHi + uSlack; };
    int useK = 0;
    bool fit = false;
    for (int kk : {0, 1, -1}) {
        if (inside(a0.X() + twoPi * kk) && inside(b0.X() + twoPi * kk)) {
            useK = kk;
            fit = true;
            break;
        }
    }
    if (!fit) return {};
    const gp_Pnt2d a(a0.X() + twoPi * useK, a0.Y());
    const gp_Pnt2d b(b0.X() + twoPi * useK, b0.Y());
    const double df = l - f;
    const gp_Vec2d d2((b.X() - a.X()) / df, (b.Y() - a.Y()) / df);
    if (d2.SquareMagnitude() <= Precision::PConfusion() * Precision::PConfusion()) return {};
    try {
        // Geom2d_Line is unit-speed. Iso-curves of the cylinder (circle, generator)
        // have UV speed 1 on the 3D parameter, so the line hits both vertices.
        // A chord's UV step is not that parameter length; a degree-1 spline carries
        // the same segment on [f, l].
        if (std::abs(d2.Magnitude() - 1.0) <= Precision::PConfusion())
            return new Geom2d_Line(gp_Pnt2d(a.X() - f * d2.X(), a.Y() - f * d2.Y()), gp_Dir2d(d2));
        TColgp_Array1OfPnt2d poles(1, 2);
        poles.SetValue(1, a);
        poles.SetValue(2, b);
        TColStd_Array1OfReal knots(1, 2);
        knots.SetValue(1, f);
        knots.SetValue(2, l);
        TColStd_Array1OfInteger mults(1, 2);
        mults.SetValue(1, 2);
        mults.SetValue(2, 2);
        return new Geom2d_BSplineCurve(poles, knots, mults, 1);
    } catch (const Standard_Failure&) {
        return {};
    }
}

// D-fc-7.10. `cap` is the caller's existing partial-face budget. A face that
// is already BRepCheck-valid is returned untouched — the chart is not
// installed. An invalid face whose surface is a partial-cylinder trim is
// rebuilt by installing partialTrimChartPCurve on each edge that lies on that
// trim within `cap`. The edge tolerance is the measured deviation when that
// is larger, and it is restored with the pcurves if the rebuild fails.
// The rebuild is kept
// only when BRepCheck then passes; otherwise every pcurve, range, tolerance
// and parameter flag is restored and today's path stands. True only when an
// invalid face became valid.
// `tau` is tauOf (2q). When the chart is installed and the face is still
// invalid, the D-fc-7.10 corner fallback runs on that face before the
// pcurves are restored. tau <= 0 skips the fallback.
bool healInvalidPartialTrimChart(TopoDS_Face& face, double cap, double tau = 0);

// D-fc-7.10 fallback, used only when the chart above cannot validate the
// face. `tau` is tauOf (2q). A face that is already BRepCheck-valid is
// returned untouched. On an invalid face, a shared corner whose pcurves
// cross away from the vertex takes that vertex's tolerance equal to the
// measured 3D distance from the vertex to the crossing, and only when that
// distance is <= tau. A larger distance is refused, not clamped. The saved
// tolerance is restored unless BRepCheck then passes. True only when an
// invalid face became valid.
bool raiseMeasuredCornerTol(TopoDS_Face& face, double tau);

// D-fc-7.16. Either trim bound lies within Precision::PConfusion() of the
// periodic seam (0 and 2π are the same angle). The tolerance is OCCT's
// parametric confusion, not a new length.
inline bool partialTrimTouchesPeriodicSeam(double uLo, double uHi) {
    const double twoPi = 2.0 * std::acos(-1.0);
    const double tol = Precision::PConfusion();
    auto nearSeam = [&](double u) { return std::fabs(std::remainder(u, twoPi)) <= tol; };
    return nearSeam(uLo) || nearSeam(uHi);
}

// D-fc-7.16. XDirection rotation about the cylinder axis that centres this
// trim at u = π. Measured from the trim's own bounds. No new threshold.
inline double partialCylSeamFrameAngle(double uLo, double uHi) {
    return 0.5 * (uLo + uHi) - std::acos(-1.0);
}

// D-fc-7.16. The rotated frame applies only to a partial (not full-turn)
// cylinder whose trim touches the seam, or whose face already carries a seam
// edge (nSeamE > 0). A full turn is a closed cylinder and is not this heal.
inline bool partialCylSeamFrameApplies(double uLo, double uHi, int nSeamEdges) {
    const double twoPi = 2.0 * std::acos(-1.0);
    const double span = uHi - uLo;
    if (!(span > Precision::Confusion() && span < twoPi - Precision::Confusion())) return false;
    return partialTrimTouchesPeriodicSeam(uLo, uHi) || nSeamEdges > 0;
}

// D-fc-7.16. A face that is already BRepCheck-valid is returned untouched.
// An invalid partial cylinder whose trim touches the periodic seam (or whose
// wire carries a seam edge) is rebuilt on the same cylinder — same axis,
// radius and location — with XDirection rotated by partialCylSeamFrameAngle
// so the trim is centred at π. Every boundary edge gets one pcurve from
// partialTrimChartPCurve on that rotated trim. The rebuild ships only when
// the face and a shell that contains it both pass BRepCheck; otherwise the
// added pcurves and any raised edge tolerance are restored. True only when
// an invalid face became valid.
bool healPartialCylSeamFrame(TopoDS_Face& face, double cap);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_AXIS_SEAM_HPP
