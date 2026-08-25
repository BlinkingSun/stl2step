// stl2step P2 — refit::buildFaces (analytic face construction, shared topology).
// Consumes a RegionSet + the engine's vertex slot array and produces the same
// vector<TopoDS_Face> the per-triangle loop produces today, except accepted
// regions become analytic (Geom_Plane / Geom_CylindricalSurface).
//
// Binding: SPEC-P2, DECISION-insertion §2 / §2.3 / §2.4 / §2.5 / §3.2 / §5,
// DECISION-p1-math D3. Shared topology, no sewing pass (J5).
//
// SPDX-License-Identifier: MIT

#include "refit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <ElCLib.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_ResultType.hxx>
#include <Precision.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Vec2d.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

const double kPi = 3.14159265358979323846;

void emit(WarnFn warn, const std::string& msg) {
    if (warn) warn(msg);
}

gp_Pnt pntOf(const MeshView& mv, int localVtx) {
    if (localVtx < 0 || (size_t)localVtx >= mv.nVtx || !mv.compVtx || !mv.pts)
        return gp_Pnt(0.0, 0.0, 0.0);
    return gp_Pnt(mv.pts[mv.compVtx[localVtx]]);
}

// R2: a RegionSet whose indices do not address this MeshView/verts[] must
// not crash. Return false (abandon) on any out-of-range id.
bool regionSetConsistent(const MeshView& mv, const RegionSet& rs,
                         const std::vector<TopoDS_Vertex>& verts) {
    if (!mv.pts || !mv.tris || !mv.compTris || !mv.compVtx || !mv.compEdges || !mv.triEdges ||
        !mv.triDirs)
        return false;
    if (mv.nTri == 0 || verts.size() < mv.nVtx) return false;
    for (size_t i = 0; i < mv.nTri; i++) {
        if (mv.compTris[i] < 0) return false;
        for (int s = 0; s < 3; s++) {
            int e = mv.triEdges[i][s];
            if (e < 0 || (size_t)e >= mv.nEdge) return false;
        }
    }
    for (size_t i = 0; i < mv.nEdge; i++) {
        int a = mv.compEdges[i].first, b = mv.compEdges[i].second;
        if (a < 0 || b < 0 || (size_t)a >= mv.nVtx || (size_t)b >= mv.nVtx) return false;
    }
    auto knownRegion = [&](int id) -> bool {
        if (id < 0) return true;
        for (const Region& r : rs.regions)
            if (r.id == id) return true;
        return false;
    };
    for (const Region& r : rs.regions) {
        if (r.id < 0) return false;
        for (int t : r.tris)
            if (t < 0 || (size_t)t >= mv.nTri) return false;
        for (const Loop& lp : r.loops) {
            if (lp.chainIdx.size() != lp.reversed.size()) return false;
            for (int ci : lp.chainIdx)
                if (ci < 0 || (size_t)ci >= rs.chains.size()) return false;
        }
    }
    for (const BoundaryChain& ch : rs.chains) {
        if (!knownRegion(ch.regA) || !knownRegion(ch.regB)) return false;
        for (int e : ch.meshEdges)
            if (e < 0 || (size_t)e >= mv.nEdge) return false;
        for (int v : ch.meshVerts)
            if (v < 0 || (size_t)v >= mv.nVtx) return false;
    }
    size_t nMap = std::min(mv.nTri, rs.triRegion.size());
    for (size_t k = 0; k < nMap; k++)
        if (!knownRegion(rs.triRegion[k])) return false;
    return true;
}

bool edgeSpansFullCircle(const TopoDS_Edge& e) {
    if (e.IsNull()) return false;
    Standard_Real f = 0, l = 0;
    Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
    if (c.IsNull()) return false;
    return std::fabs(std::fabs(l - f) - 2.0 * kPi) <= 0.05;
}

gp_Pln asPlane(const Region& r) { return gp_Pln(r.ax); }
gp_Cylinder asCyl(const Region& r) { return gp_Cylinder(r.ax, r.radius); }

bool isAnalytic(const Region* r) {
    return r && (r->type == SurfType::Plane || r->type == SurfType::Cylinder);
}

const Region* regionById(const RegionSet& rs, int id) {
    if (id < 0) return nullptr;
    for (const Region& r : rs.regions)
        if (r.id == id) return &r;
    return nullptr;
}

Region* regionByIdMut(RegionSet& rs, int id) {
    if (id < 0) return nullptr;
    for (Region& r : rs.regions)
        if (r.id == id) return &r;
    return nullptr;
}

int otherReg(const BoundaryChain& ch, int id) {
    if (ch.regA == id) return ch.regB;
    if (ch.regB == id) return ch.regA;
    return -1;
}

bool chainTouches(const BoundaryChain& ch, int id) {
    return ch.regA == id || ch.regB == id;
}

bool bothAnalytic(const BoundaryChain& ch, const RegionSet& rs) {
    return isAnalytic(regionById(rs, ch.regA)) && isAnalytic(regionById(rs, ch.regB));
}

bool meshComponentClosed(const MeshView& mv) {
    std::vector<int> use(mv.nEdge, 0);
    for (size_t k = 0; k < mv.nTri; k++) {
        for (int s = 0; s < 3; s++) {
            int e = mv.triEdges[k][s];
            if (e >= 0 && (size_t)e < mv.nEdge) use[(size_t)e]++;
        }
    }
    for (size_t i = 0; i < mv.nEdge; i++)
        if (use[i] == 1) return false;
    return mv.nTri > 0;
}

double azimuthOf(const Region& r, const gp_Pnt& p) {
    gp_Vec rho(r.ax.Location(), p);
    double v = rho.Dot(r.ax.Direction());
    (void)v;
    gp_Vec rad = rho - gp_Vec(r.ax.Direction()) * rho.Dot(r.ax.Direction());
    if (rad.Magnitude() < Precision::Confusion()) return 0.0;
    double x = rad.Dot(r.ax.XDirection());
    double y = rad.Dot(r.ax.YDirection());
    double u = std::atan2(y, x);
    if (u < 0.0) u += 2.0 * kPi;
    return u;
}

int seamVertexOf(const MeshView& mv, const Region& cyl, const BoundaryChain& ch) {
    // Vertex whose azimuth is closest to 0 (XDirection is already a mesh-vertex azimuth).
    int best = ch.meshVerts.empty() ? -1 : ch.meshVerts.front();
    double bestU = 1e300;
    for (int lv : ch.meshVerts) {
        double u = azimuthOf(cyl, pntOf(mv, lv));
        double d = std::min(u, 2.0 * kPi - u);
        if (d < bestU) {
            bestU = d;
            best = lv;
        }
    }
    return best;
}

int vertexClosestToU(const MeshView& mv, const Region& cyl, const BoundaryChain& ch, double target) {
    int best = ch.meshVerts.empty() ? -1 : ch.meshVerts.front();
    double bestD = 1e300;
    for (int lv : ch.meshVerts) {
        if (lv < 0 || (size_t)lv >= mv.nVtx) continue;
        double u = azimuthOf(cyl, pntOf(mv, lv));
        double d = std::fabs(u - target);
        d = std::min(d, 2.0 * kPi - d);
        if (d < bestD) {
            bestD = d;
            best = lv;
        }
    }
    return best;
}

int vertexClosestToUOnLoop(const MeshView& mv, const Region& cyl, const Loop& lp,
                           const RegionSet& rs, double target) {
    int best = -1;
    double bestD = 1e300;
    for (int ci : lp.chainIdx) {
        if (ci < 0 || (size_t)ci >= rs.chains.size()) continue;
        int v = vertexClosestToU(mv, cyl, rs.chains[(size_t)ci], target);
        if (v < 0) continue;
        double u = azimuthOf(cyl, pntOf(mv, v));
        double d = std::fabs(u - target);
        d = std::min(d, 2.0 * kPi - d);
        if (d < bestD) {
            bestD = d;
            best = v;
        }
    }
    return best;
}

bool faceIsValid(const TopoDS_Face& f) {
    if (f.IsNull()) return false;
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        return an.IsValid() == Standard_True;
    } catch (const Standard_Failure&) {
        return false;
    }
}

// F3: DECISION §5 post-fit residual gate. Chord-midpoint radial deviation on
// a cylinder is the sagitta; sagitta/R ≤ 1−cos(π/4) admits nSides≥4 and
// rejects the 3-sided closed360 fixtures (r1_success / r8 / r1_round2).
bool cylinderPostFitOk(const Region& r, const MeshView& mv, const RegionSet& rs) {
    if (r.type != SurfType::Cylinder || r.radius <= 0.0) return true;
    const double lim = r.radius * (1.0 - std::cos(kPi / 4.0)) + Precision::Confusion();
    gp_Cylinder cyl = asCyl(r);
    auto radDev = [&](const gp_Pnt& p) {
        gp_Vec v(cyl.Location(), p);
        return std::fabs(gp_Vec(cyl.Axis().Direction()).Crossed(v).Magnitude() - cyl.Radius());
    };
    double d = 0.0;
    for (int lt : r.tris) {
        if (lt < 0 || (size_t)lt >= mv.nTri) continue;
        const int* T = mv.tris[mv.compTris[lt]];
        for (int i = 0; i < 3; i++) d = std::max(d, radDev(gp_Pnt(mv.pts[T[i]])));
        if (d > lim) return false;
    }
    for (const BoundaryChain& ch : rs.chains) {
        if (ch.regA != r.id && ch.regB != r.id) continue;
        for (int eid : ch.meshEdges) {
            if (eid < 0 || (size_t)eid >= mv.nEdge) continue;
            const auto& ev = mv.compEdges[eid];
            gp_Pnt a = pntOf(mv, ev.first), b = pntOf(mv, ev.second);
            d = std::max(d, radDev(gp_Pnt(0.5 * (a.XYZ() + b.XYZ()))));
            if (d > lim) return false;
        }
    }
    return true;
}

// F4: a face whose boundary fell back to mesh chords must absorb sagitta in
// edge/vertex tolerance or fail (never emit BRepCheck-invalid).
bool ensureFaceValid(TopoDS_Face& f) {
    if (f.IsNull()) return false;
    if (faceIsValid(f)) return true;
    try {
        BRep_Builder B;
        BRepAdaptor_Surface s(f, Standard_False);
        auto devPnt = [&](const gp_Pnt& p) {
            if (s.GetType() == GeomAbs_Plane) return s.Plane().Distance(p);
            if (s.GetType() == GeomAbs_Cylinder) {
                gp_Cylinder c = s.Cylinder();
                gp_Vec v(c.Location(), p);
                return std::fabs(gp_Vec(c.Axis().Direction()).Crossed(v).Magnitude() - c.Radius());
            }
            return 0.0;
        };
        for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
            const TopoDS_Vertex& vv = TopoDS::Vertex(vx.Current());
            double d = devPnt(BRep_Tool::Pnt(vv));
            if (d > BRep_Tool::Tolerance(vv))
                B.UpdateVertex(vv, d * 1.001 + Precision::Confusion());
        }
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
            double d = 0.0;
            for (TopExp_Explorer vx(e, TopAbs_VERTEX); vx.More(); vx.Next())
                d = std::max(d, devPnt(BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()))));
            Standard_Real fp = 0, lp = 0;
            Handle(Geom_Curve) c = BRep_Tool::Curve(e, fp, lp);
            if (!c.IsNull()) d = std::max(d, devPnt(c->Value(0.5 * (fp + lp))));
            if (d > BRep_Tool::Tolerance(e))
                B.UpdateEdge(e, d * 1.001 + Precision::Confusion());
        }
    } catch (const Standard_Failure&) {
        return false;
    }
    return faceIsValid(f);
}

void addPcurvesOnFace(TopoDS_Face& f, double sewTol, bool /*closed360*/) {
    ShapeFix_Edge sfe;
    TopTools_IndexedMapOfShape emap;
    TopExp::MapShapes(f, TopAbs_EDGE, emap);
    // Count how many times each TShape appears on the face (seam = 2).
    std::vector<int> count(static_cast<size_t>(emap.Extent()) + 1, 0);
    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
        int i = emap.FindIndex(ex.Current());
        if (i > 0) count[static_cast<size_t>(i)]++;
    }
    for (int i = 1; i <= emap.Extent(); i++) {
        const TopoDS_Edge& e = TopoDS::Edge(emap(i));
        Standard_Boolean isSeam = count[static_cast<size_t>(i)] >= 2 ? Standard_True
                                                                    : Standard_False;
        sfe.FixAddPCurve(e, f, isSeam, sewTol);
    }
}

void addPcurvesOnSurface(const Handle(Geom_Surface)& surf, TopoDS_Edge& e, bool isSeam,
                         double sewTol) {
    ShapeFix_Edge sfe;
    sfe.FixAddPCurve(e, surf, TopLoc_Location(), isSeam ? Standard_True : Standard_False, sewTol);
}

TopoDS_Edge orientEdgeFromTo(const TopoDS_Edge& e, const TopoDS_Vertex& from) {
    TopoDS_Vertex v1, v2;
    TopExp::Vertices(e, v1, v2, Standard_True);  // with orientation
    if (v1.IsNull()) return e;
    if (v1.IsSame(from)) return e;
    return TopoDS::Edge(e.Reversed());
}

// Walk a chain's local vertex ids in loop order.
std::vector<int> walkVerts(const BoundaryChain& ch, bool reversed) {
    std::vector<int> v = ch.meshVerts;
    if (reversed) std::reverse(v.begin(), v.end());
    return v;
}

int edgeConnecting(const MeshView& mv, const BoundaryChain& ch, int va, int vb) {
    for (int eid : ch.meshEdges) {
        if (eid < 0 || (size_t)eid >= mv.nEdge) continue;
        const auto& e = mv.compEdges[eid];
        if ((e.first == va && e.second == vb) || (e.first == vb && e.second == va)) return eid;
    }
    return -1;
}

struct ChainGeom {
    bool collapsed = false;           // analytic (possibly 1 or 2 edges)
    std::vector<TopoDS_Edge> edges;   // collapsed edges, in walk-forward (regA-left) order
};

struct AnalyticCurve {
    enum Kind { None, Lin, Circ, Elips } kind = None;
    gp_Lin lin;
    gp_Circ circ;
    gp_Elips elips;
};

double curveResidual(const AnalyticCurve& c, const gp_Pnt& p) {
    try {
        if (c.kind == AnalyticCurve::Lin)
            return c.lin.Distance(p);
        if (c.kind == AnalyticCurve::Circ) {
            double t = ElCLib::Parameter(c.circ, p);
            return ElCLib::Value(t, c.circ).Distance(p);
        }
        if (c.kind == AnalyticCurve::Elips) {
            double t = ElCLib::Parameter(c.elips, p);
            return ElCLib::Value(t, c.elips).Distance(p);
        }
    } catch (const Standard_Failure&) {
        return 1e300;
    }
    return 1e300;
}

double chainResidual(const AnalyticCurve& c, const MeshView& mv, const BoundaryChain& ch) {
    if (ch.meshVerts.empty()) return 1e300;
    double s = 0;
    for (int lv : ch.meshVerts) s += curveResidual(c, pntOf(mv, lv));
    return s / (double)ch.meshVerts.size();
}

AnalyticCurve pickIntAna(const IntAna_QuadQuadGeo& iq, const MeshView& mv, const BoundaryChain& ch,
                         double sewTol) {
    AnalyticCurve best;
    if (!iq.IsDone()) return best;
    IntAna_ResultType ty = iq.TypeInter();
    if (ty == IntAna_Empty || ty == IntAna_Same || ty == IntAna_NoGeometricSolution) return best;

    double bestR = 1e300;
    auto consider = [&](const AnalyticCurve& cand) {
        double r = chainResidual(cand, mv, ch);
        if (r < bestR) {
            bestR = r;
            best = cand;
        }
    };

    int n = iq.NbSolutions();
    if (ty == IntAna_Line || ty == IntAna_PointAndCircle) {
        for (int i = 1; i <= n; i++) {
            AnalyticCurve c;
            try {
                c.kind = AnalyticCurve::Lin;
                c.lin = iq.Line(i);
                consider(c);
            } catch (const Standard_Failure&) {
            }
        }
    }
    if (ty == IntAna_Circle || ty == IntAna_PointAndCircle) {
        for (int i = 1; i <= n; i++) {
            AnalyticCurve c;
            try {
                c.kind = AnalyticCurve::Circ;
                c.circ = iq.Circle(i);
                consider(c);
            } catch (const Standard_Failure&) {
            }
        }
    }
    if (ty == IntAna_Ellipse) {
        for (int i = 1; i <= n; i++) {
            AnalyticCurve c;
            try {
                c.kind = AnalyticCurve::Elips;
                c.elips = iq.Ellipse(i);
                consider(c);
            } catch (const Standard_Failure&) {
            }
        }
    }
    if (best.kind != AnalyticCurve::None && bestR > std::max(sewTol * 50.0, 1.0)) {
        // Intersection exists but is the wrong branch relative to this chain.
        best.kind = AnalyticCurve::None;
    }
    return best;
}

AnalyticCurve constructedGenerator(const Region& cylR, const Region& plnR) {
    AnalyticCurve out;
    gp_Cylinder cyl = asCyl(cylR);
    gp_Pln pln = asPlane(plnR);
    gp_Pnt axp = cyl.Location();
    gp_Dir a = cyl.Axis().Direction();
    gp_Dir n = pln.Axis().Direction();
    double sd = gp_Vec(pln.Location(), axp).Dot(gp_Vec(n));
    gp_XYZ toward = (sd >= 0.0) ? n.XYZ().Reversed() : n.XYZ();
    double mag = toward.Modulus();
    if (mag < Precision::Confusion()) return out;
    toward /= mag;
    gp_Pnt origin(axp.XYZ() + toward * cyl.Radius());
    out.kind = AnalyticCurve::Lin;
    out.lin = gp_Lin(origin, a);
    return out;
}

AnalyticCurve constructedCylCylGenerator(const Region& a, const Region& b) {
    AnalyticCurve out;
    gp_Cylinder c1 = asCyl(a), c2 = asCyl(b);
    gp_Pnt p1 = c1.Location(), p2 = c2.Location();
    gp_Dir d1 = c1.Axis().Direction();
    gp_Vec delta(p1, p2);
    delta -= gp_Vec(d1) * delta.Dot(d1);
    double dist = delta.Magnitude();
    if (dist < Precision::Confusion()) return out;
    gp_Dir n(delta);
    double R1 = c1.Radius(), R2 = c2.Radius();
    gp_XYZ origin;
    if (std::fabs(dist - (R1 + R2)) <= std::max(1e-6, 1e-3 * dist)) {
        origin = p1.XYZ() + n.XYZ() * R1;
    } else if (std::fabs(dist - std::fabs(R1 - R2)) <= std::max(1e-6, 1e-3 * dist)) {
        double s = (R1 >= R2) ? R1 : -R1;
        origin = p1.XYZ() + n.XYZ() * s;
    } else {
        origin = p1.XYZ() + n.XYZ() * R1;
    }
    out.kind = AnalyticCurve::Lin;
    out.lin = gp_Lin(gp_Pnt(origin), d1);
    return out;
}

AnalyticCurve intersectSurfaces(const Region& A, const Region& B, const MeshView& mv,
                                const BoundaryChain& ch, double sewTol, WarnFn warn) {
    const double tolAng = Precision::Angular();
    const double tol = std::max(sewTol, Precision::Confusion());
    try {
        if (A.type == SurfType::Plane && B.type == SurfType::Plane) {
            IntAna_QuadQuadGeo iq(asPlane(A), asPlane(B), tolAng, tol);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: IntAna plane|plane empty/same — keeping mesh polyline");
            return c;
        }
        if (A.type == SurfType::Plane && B.type == SurfType::Cylinder) {
            double H = std::fabs(B.vMax - B.vMin);
            IntAna_QuadQuadGeo iq(asPlane(A), asCyl(B), tolAng, tol, H);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: IntAna plane|cyl empty/same — keeping mesh polyline");
            return c;
        }
        if (A.type == SurfType::Cylinder && B.type == SurfType::Plane) {
            return intersectSurfaces(B, A, mv, ch, sewTol, warn);
        }
        if (A.type == SurfType::Cylinder && B.type == SurfType::Cylinder) {
            IntAna_QuadQuadGeo iq(asCyl(A), asCyl(B), tol);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: IntAna cyl|cyl empty/same — keeping mesh polyline");
            return c;
        }
    } catch (const Standard_Failure&) {
        emit(warn, "smooth: IntAna threw — keeping mesh polyline");
    }
    return {};
}

TopoDS_Edge makeEdgeFromCurve(const AnalyticCurve& c, const TopoDS_Vertex& v1,
                              const TopoDS_Vertex& v2, bool closedFull) {
    try {
        if (c.kind == AnalyticCurve::Lin) {
            BRepBuilderAPI_MakeEdge me(c.lin, v1, v2);
            if (me.IsDone()) return me.Edge();
        } else if (c.kind == AnalyticCurve::Circ) {
            if (closedFull || v1.IsSame(v2)) {
                BRepBuilderAPI_MakeEdge me(c.circ, v1, v1);
                if (me.IsDone()) return me.Edge();
                BRepBuilderAPI_MakeEdge me2(c.circ);
                if (me2.IsDone()) return me2.Edge();
            } else {
                // Mesh chain is the selector: pick the arc whose midpoint matches the chain.
                BRepBuilderAPI_MakeEdge me(c.circ, v1, v2);
                if (me.IsDone()) return me.Edge();
            }
        } else if (c.kind == AnalyticCurve::Elips) {
            if (closedFull || v1.IsSame(v2)) {
                BRepBuilderAPI_MakeEdge me(c.elips, v1, v1);
                if (me.IsDone()) return me.Edge();
            } else {
                BRepBuilderAPI_MakeEdge me(c.elips, v1, v2);
                if (me.IsDone()) return me.Edge();
            }
        }
    } catch (const Standard_Failure&) {
    }
    return TopoDS_Edge();
}

TopoDS_Edge makeArc(const gp_Circ& circ, const TopoDS_Vertex& vA, const TopoDS_Vertex& vB,
                    const gp_Pnt& midHint) {
    // Bound the circle with the shared verts[] slots; the mesh mid-vertex picks
    // which of the two arcs (critical at 180°).
    const double twopi = 2.0 * kPi;
    auto wrap = [&](double t) {
        while (t < 0.0) t += twopi;
        while (t >= twopi) t -= twopi;
        return t;
    };
    try {
        double p1 = wrap(ElCLib::Parameter(circ, BRep_Tool::Pnt(vA)));
        double p2 = wrap(ElCLib::Parameter(circ, BRep_Tool::Pnt(vB)));
        double pm = wrap(ElCLib::Parameter(circ, midHint));
        double df = p2 - p1;
        while (df <= 0.0) df += twopi;
        double dm = pm - p1;
        while (dm < 0.0) dm += twopi;
        Handle(Geom_Circle) gc = new Geom_Circle(circ);
        BRep_Builder B;
        TopoDS_Edge e;
        B.MakeEdge(e, gc, Precision::Confusion());
        if (dm <= df + 1e-9) {
            B.Add(e, vA.Oriented(TopAbs_FORWARD));
            B.Add(e, vB.Oriented(TopAbs_REVERSED));
            B.Range(e, p1, p1 + df);
        } else {
            // Complementary arc: first param is V2, last is V1, then Reverse
            // so the edge still walks vA → vB along the midHint side.
            double db = twopi - df;
            B.Add(e, vB.Oriented(TopAbs_FORWARD));
            B.Add(e, vA.Oriented(TopAbs_REVERSED));
            B.Range(e, p1 - db, p1);
            e.Reverse();
        }
        return e;
    } catch (const Standard_Failure&) {
    }
    try {
        BRepBuilderAPI_MakeEdge me(circ, vA, vB);
        if (me.IsDone()) return me.Edge();
    } catch (const Standard_Failure&) {
    }
    return TopoDS_Edge();
}

bool triangleOnExplodedOrIsland(const RegionSet& rs, const std::vector<char>& exploded, size_t k) {
    if (k >= rs.triIsland.size()) return false;
    if (rs.triIsland[k] >= 0) return true;
    int rid = (k < rs.triRegion.size()) ? rs.triRegion[k] : -1;
    if (rid < 0) return true;
    if (rid >= 0 && (size_t)rid < exploded.size() && exploded[(size_t)rid]) return true;
    // exploded is indexed by region id, which is dense
    if ((size_t)rid < exploded.size() && exploded[(size_t)rid]) return true;
    return false;
}

bool regionExploded(const std::vector<char>& exploded, int id) {
    return id >= 0 && (size_t)id < exploded.size() && exploded[(size_t)id] != 0;
}

TopoDS_Face makeFacet(const MeshView& mv, const std::vector<TopoDS_Vertex>& verts,
                      const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                      size_t k) {
    (void)verts;
    if (k >= mv.nTri || !mv.compTris || !mv.tris || !mv.pts || !mv.triEdges || !mv.triDirs)
        return TopoDS_Face();
    const int gt = mv.compTris[k];
    if (gt < 0) return TopoDS_Face();
    const int* T = mv.tris[gt];
    const gp_XYZ &A = mv.pts[T[0]], &Bp = mv.pts[T[1]], &Cp = mv.pts[T[2]];
    gp_XYZ n = (Bp - A).Crossed(Cp - A);
    double mag = n.Modulus();
    double l2 = std::max((Bp - A).SquareModulus(), (Cp - A).SquareModulus());
    if (mag < 1e-12 || mag * mag < l2 * l2 * 1e-20) return TopoDS_Face();
    BRep_Builder bb;
    TopoDS_Wire w;
    bb.MakeWire(w);
    for (int s = 0; s < 3; s++) {
        int id = mv.triEdges[k][s];
        if (id < 0 || (size_t)id >= meshE.size() || !edgeOk[(size_t)id]) return TopoDS_Face();
        bool fwd = (mv.triDirs[k] >> s) & 1;
        bb.Add(w, fwd ? meshE[(size_t)id] : TopoDS::Edge(meshE[(size_t)id].Reversed()));
    }
    w.Closed(Standard_True);
    // BRep_Builder::MakeFace + Add keeps the wire's verts[] slots (J2).
    // MakeFace(gp_Pln, wire) copies vertices and is the explode-path twin source.
    try {
        Handle(Geom_Plane) pln = new Geom_Plane(gp_Pln(gp_Pnt(A), gp_Dir(n)));
        TopoDS_Face f;
        bb.MakeFace(f, pln, Precision::Confusion());
        bb.Add(f, w);
        return f;
    } catch (const Standard_Failure&) {
        return TopoDS_Face();
    }
}

bool appendPolyline(BRep_Builder& B, TopoDS_Wire& w, const MeshView& mv, const BoundaryChain& ch,
                    bool reversed, const std::vector<TopoDS_Edge>& meshE,
                    const std::vector<char>& edgeOk) {
    std::vector<int> vs = walkVerts(ch, reversed);
    if (vs.size() < 2 && !(ch.closedLoop && vs.size() == 1)) {
        // closed with verts==edges; need >= 2 for a polyline
    }
    size_t nSeg = ch.meshEdges.size();
    if (vs.empty() || nSeg == 0) return false;
    for (size_t i = 0; i < nSeg; i++) {
        int va = vs[i % vs.size()];
        int vb = vs[(i + 1) % vs.size()];
        if (!ch.closedLoop && i + 1 >= vs.size()) return false;
        int eid = edgeConnecting(mv, ch, va, vb);
        if (eid < 0) {
            // fall back to stored order
            size_t idx = reversed ? (nSeg - 1 - i) : i;
            if (idx >= ch.meshEdges.size()) return false;
            eid = ch.meshEdges[idx];
        }
        if (eid < 0 || (size_t)eid >= meshE.size() || !edgeOk[(size_t)eid]) return false;
        const auto& ev = mv.compEdges[eid];
        bool fwd = (ev.first == va);
        B.Add(w, fwd ? meshE[(size_t)eid] : TopoDS::Edge(meshE[(size_t)eid].Reversed()));
    }
    return true;
}

bool appendCollapsed(BRep_Builder& B, TopoDS_Wire& w, const ChainGeom& g, bool reversed) {
    if (g.edges.empty()) return false;
    if (!reversed) {
        for (const auto& e : g.edges) B.Add(w, e);
    } else {
        for (int i = (int)g.edges.size() - 1; i >= 0; i--)
            B.Add(w, TopoDS::Edge(g.edges[(size_t)i].Reversed()));
    }
    return true;
}

bool buildLoopWire(TopoDS_Wire& w, const Loop& loop, const RegionSet& rs, const MeshView& mv,
                   const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                   const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk) {
    BRep_Builder B;
    B.MakeWire(w);
    if (loop.chainIdx.size() != loop.reversed.size()) return false;
    for (size_t i = 0; i < loop.chainIdx.size(); i++) {
        int ci = loop.chainIdx[i];
        if (ci < 0 || (size_t)ci >= rs.chains.size()) return false;
        bool rev = loop.reversed[i] != 0;
        bool ok;
        if ((size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
            (size_t)ci < geom.size() && geom[(size_t)ci].collapsed && !geom[(size_t)ci].edges.empty())
            ok = appendCollapsed(B, w, geom[(size_t)ci], rev);
        else
            ok = appendPolyline(B, w, mv, rs.chains[(size_t)ci], rev, meshE, edgeOk);
        if (!ok) return false;
    }
    w.Closed(Standard_True);
    return true;
}

// Fix the call: appendPolyline signature uses BoundaryChain not RegionSet. Oops — I
// inlined the wrong name above. The definition takes BoundaryChain.

void setFaceOutward(TopoDS_Face& f, bool outwardNormal) {
    // Same discipline as trySeamed360: SET orientation from outwardNormal,
    // never toggle on a failed BRepCheck (that undoes a hole's Reverse).
    if (f.IsNull()) return;
    if (!outwardNormal && f.Orientation() == TopAbs_FORWARD) f.Reverse();
    else if (outwardNormal && f.Orientation() == TopAbs_REVERSED) f.Reverse();
}

double regionU(const Region& r, const gp_Pnt& p) {
    double u = azimuthOf(r, p);  // [0, 2π)
    // Map into (−π, π] so a half-hole with uMin=-π, uMax=0 is a simple interval
    // and does not cross the cylinder's periodic seam at 0/2π in 2d.
    if (u > kPi) u -= 2.0 * kPi;
    return u;
}

double regionV(const Region& r, const gp_Pnt& p) {
    return gp_Vec(r.ax.Location(), p).Dot(r.ax.Direction());
}

void bindCylPCurves(TopoDS_Wire& w, const Handle(Geom_Surface)& surf, const Region& r,
                    double sewTol) {
    BRep_Builder B;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        TopoDS_Edge e = TopoDS::Edge(it.Value());
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        gp_Pnt pa, pb;
        if (!c3.IsNull()) {
            pa = c3->Value(f3);
            pb = c3->Value(l3);
        } else {
            TopoDS_Vertex va, vb;
            TopExp::Vertices(e, va, vb, Standard_False);
            if (va.IsNull() || vb.IsNull()) continue;
            pa = BRep_Tool::Pnt(va);
            pb = BRep_Tool::Pnt(vb);
        }
        double u1 = regionU(r, pa), v1 = regionV(r, pa);
        double u2 = regionU(r, pb), v2 = regionV(r, pb);
        while (u2 - u1 > kPi) u2 -= 2.0 * kPi;
        while (u1 - u2 > kPi) u2 += 2.0 * kPi;
        gp_Vec2d duv(u2 - u1, v2 - v1);
        double mag = duv.Magnitude();
        if (mag < Precision::PConfusion()) continue;
        Handle(Geom2d_Line) ln = new Geom2d_Line(gp_Pnt2d(u1, v1), gp_Dir2d(duv));
        B.UpdateEdge(e, ln, surf, TopLoc_Location(), sewTol);
        B.Range(e, surf, TopLoc_Location(), 0.0, mag);
    }
}

bool addPcurvesThenFace(const Handle(Geom_Surface)& surf, TopoDS_Wire& outer,
                        const std::vector<TopoDS_Wire>& inners, double sewTol, bool outward,
                        TopoDS_Face& outF) {
    // Pcurves on every edge before MakeFace on a non-plane (spike-occt §3.2).
    auto pcurveWire = [&](TopoDS_Wire& ww) {
        for (TopoDS_Iterator it(ww); it.More(); it.Next()) {
            TopoDS_Edge e = TopoDS::Edge(it.Value());
            addPcurvesOnSurface(surf, e, false, sewTol);
        }
    };
    pcurveWire(outer);
    for (auto& iw : inners) pcurveWire(const_cast<TopoDS_Wire&>(iw));

    try {
        BRepBuilderAPI_MakeFace mf(surf, outer, Standard_True);
        if (!mf.IsDone()) return false;
        for (const auto& iw : inners) mf.Add(iw);
        if (!mf.IsDone()) return false;
        outF = mf.Face();
        setFaceOutward(outF, outward);
        addPcurvesOnFace(outF, sewTol, false);
        return !outF.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool buildPlanarFace(const Region& r, const RegionSet& rs, const MeshView& mv,
                     const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                     const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                     TopoDS_Face& outF) {
    const Loop* outer = nullptr;
    std::vector<const Loop*> inners;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::Outer) outer = &lp;
        else if (lp.role == LoopRole::Inner) inners.push_back(&lp);
    }
    if (!outer) return false;
    TopoDS_Wire ow;
    if (!buildLoopWire(ow, *outer, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
    auto addInners = [&](BRepBuilderAPI_MakeFace& mf) -> bool {
        for (const Loop* ip : inners) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, *ip, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
            mf.Add(iw);
        }
        return mf.IsDone() == Standard_True;
    };
    try {
        BRepBuilderAPI_MakeFace mf(asPlane(r), ow, Standard_True);
        if (!mf.IsDone()) return false;
        if (!addInners(mf)) return false;
        outF = mf.Face();
        setFaceOutward(outF, r.outwardNormal);
        if (!faceIsValid(outF)) {
            // Mixed loops can be two endpoint-sharing paths rather than one
            // circulating walk. Reconnect by vertex identity (MakeWire) without
            // copying TShapes — s09; skip when the original wire already works
            // (slotted_stadium planes).
            try {
                BRepBuilderAPI_MakeWire mw;
                for (TopoDS_Iterator it(ow); it.More(); it.Next())
                    mw.Add(TopoDS::Edge(it.Value()));
                if (mw.IsDone()) {
                    TopoDS_Wire ow2 = mw.Wire();
                    ow2.Closed(Standard_True);
                    BRepBuilderAPI_MakeFace mf2(asPlane(r), ow2, Standard_True);
                    if (mf2.IsDone() && addInners(mf2)) {
                        TopoDS_Face f2 = mf2.Face();
                        setFaceOutward(f2, r.outwardNormal);
                        if (faceIsValid(f2)) {
                            outF = f2;
                            ow = ow2;
                        }
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF)) {
            // Outer-wire winding vs plane normal: retry reversed wire.
            try {
                BRepBuilderAPI_MakeFace mfR(asPlane(r), TopoDS::Wire(ow.Reversed()),
                                            Standard_True);
                if (mfR.IsDone()) {
                    for (const Loop* ip : inners) {
                        TopoDS_Wire iw;
                        if (!buildLoopWire(iw, *ip, rs, mv, geom, collapsed, meshE, edgeOk))
                            break;
                        mfR.Add(iw);
                    }
                    if (mfR.IsDone()) {
                        TopoDS_Face fR = mfR.Face();
                        setFaceOutward(fR, r.outwardNormal);
                        if (faceIsValid(fR)) outF = fR;
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF) && !inners.empty()) {
            // Retry with inner wires flipped — imbrication/orientation of holes.
            try {
                BRepBuilderAPI_MakeFace mf2(asPlane(r), ow, Standard_True);
                if (mf2.IsDone()) {
                    for (const Loop* ip : inners) {
                        TopoDS_Wire iw;
                        if (!buildLoopWire(iw, *ip, rs, mv, geom, collapsed, meshE, edgeOk))
                            break;
                        mf2.Add(TopoDS::Wire(iw.Reversed()));
                    }
                    if (mf2.IsDone()) {
                        TopoDS_Face f2 = mf2.Face();
                        if (!r.outwardNormal) f2.Reverse();
                        if (faceIsValid(f2)) outF = f2;
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF)) {
            try {
                ShapeFix_Face sff(outF);
                sff.FixOrientationMode() = 1;
                sff.FixAddNaturalBoundMode() = 0;
                sff.FixMissingSeamMode() = 0;
                sff.Perform();
                TopoDS_Shape res = sff.Result();
                if (res.IsNull()) res = sff.Face();
                for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                    outF = TopoDS::Face(fx.Current());
                    break;
                }
            } catch (const Standard_Failure&) {
            }
        }
        // Prefer Region::outwardNormal, but keep the BRepCheck-valid orientation
        // if forcing the flag would leave Unorientable/BadOrientation (s09).
        {
            TopoDS_Face want = outF, other = outF;
            setFaceOutward(want, r.outwardNormal);
            setFaceOutward(other, !r.outwardNormal);
            if (faceIsValid(want)) outF = want;
            else if (faceIsValid(other)) outF = other;
            else outF = want;
        }
        ensureFaceValid(outF);
        return !outF.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool buildPartialCylinder(const Region& r, const RegionSet& rs, const MeshView& mv,
                          const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                          const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                          double sewTol, TopoDS_Face& outF) {
    const Loop* outer = nullptr;
    std::vector<TopoDS_Wire> inners;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::Outer) outer = &lp;
        else if (lp.role == LoopRole::Inner) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, lp, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
            inners.push_back(iw);
        }
    }
    if (!outer) return false;
    TopoDS_Wire ow;
    if (!buildLoopWire(ow, *outer, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
    Handle(Geom_CylindricalSurface) surf = new Geom_CylindricalSurface(asCyl(r));
    bindCylPCurves(ow, surf, r, sewTol);
    for (auto& iw : inners) bindCylPCurves(iw, surf, r, sewTol);
    try {
        BRepBuilderAPI_MakeFace mf(surf, ow, Standard_True);
        if (!mf.IsDone()) return false;
        for (const auto& iw : inners) mf.Add(iw);
        if (!mf.IsDone()) return false;
        outF = mf.Face();
    } catch (const Standard_Failure&) {
        return false;
    }
    setFaceOutward(outF, r.outwardNormal);
    addPcurvesOnFace(outF, sewTol, false);
    if (!faceIsValid(outF)) {
        try {
            ShapeFix_Face sff(outF);
            sff.FixOrientationMode() = 1;
            sff.FixAddNaturalBoundMode() = 0;
            sff.FixMissingSeamMode() = 0;
            sff.Perform();
            TopoDS_Shape res = sff.Result();
            if (res.IsNull()) res = sff.Face();
            for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                outF = TopoDS::Face(fx.Current());
                break;
            }
        } catch (const Standard_Failure&) {
        }
        setFaceOutward(outF, r.outwardNormal);
        addPcurvesOnFace(outF, sewTol, false);
    }
    // F4: never emit a BRepCheck-invalid face — absorb sagitta or fail (R1).
    if (!ensureFaceValid(outF)) return false;
    return faceIsValid(outF);
}

TopoDS_Edge makeFullCircle(const gp_Circ& circ, const TopoDS_Vertex& V) {
    try {
        BRepBuilderAPI_MakeEdge me(circ, V, V);
        if (me.IsDone()) {
            TopoDS_Edge e = me.Edge();
            e.Closed(Standard_True);
            return e;
        }
    } catch (const Standard_Failure&) {
    }
    try {
        Handle(Geom_Circle) gc = new Geom_Circle(circ);
        BRep_Builder B;
        TopoDS_Edge e;
        B.MakeEdge(e, gc, Precision::Confusion());
        B.Add(e, V.Oriented(TopAbs_FORWARD));
        B.Add(e, V.Oriented(TopAbs_REVERSED));
        B.Range(e, 0.0, 2.0 * kPi);
        e.Closed(Standard_True);
        return e;
    } catch (const Standard_Failure&) {
    }
    return TopoDS_Edge();
}

bool trySeamed360(const Region& r, const RegionSet& rs, const MeshView& mv,
                  const std::vector<TopoDS_Vertex>& verts, std::vector<ChainGeom>& geom,
                  std::vector<char>& collapsed, const std::vector<TopoDS_Edge>& meshE,
                  const std::vector<char>& edgeOk, double sewTol, WarnFn warn, TopoDS_Face& outF) {
    const Loop *capL = nullptr, *capH = nullptr;
    std::vector<TopoDS_Wire> inners;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::CapLow) capL = &lp;
        else if (lp.role == LoopRole::CapHigh) capH = &lp;
        else if (lp.role == LoopRole::Inner) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, lp, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
            inners.push_back(std::move(iw));
        }
    }
    if (!capL || !capH || capL->chainIdx.empty() || capH->chainIdx.empty()) return false;
    int cL = capL->chainIdx.front();
    int cH = capH->chainIdx.front();
    if (cL < 0 || cH < 0 || (size_t)cL >= rs.chains.size() || (size_t)cH >= rs.chains.size())
        return false;

    int vL = seamVertexOf(mv, r, rs.chains[(size_t)cL]);
    int vH = seamVertexOf(mv, r, rs.chains[(size_t)cH]);
    if (vL < 0 || vH < 0 || (size_t)vL >= verts.size() || (size_t)vH >= verts.size()) return false;

    Handle(Geom_CylindricalSurface) surf = new Geom_CylindricalSurface(asCyl(r));
    gp_Ax2 axLow(r.ax.Location().Translated(gp_Vec(r.ax.Direction()) * r.vMin), r.ax.Direction(),
                 r.ax.XDirection());
    gp_Ax2 axHigh(r.ax.Location().Translated(gp_Vec(r.ax.Direction()) * r.vMax), r.ax.Direction(),
                  r.ax.XDirection());
    gp_Circ circL(axLow, r.radius);
    gp_Circ circH(axHigh, r.radius);

    TopoDS_Edge eL, eH;
    // F2: the walk-direction heuristic may Reverse() a collapsed cap. The
    // pcurve written below is u: 0→2π on a FORWARD edge — normalise first.
    // A collapsed cap that is not a single 2π circle (split / 180° arc) is
    // rung 2 (TwoHalves), not a seamed face with a desynced Range(0,2π).
    auto takeFullCap = [&](int ci, const gp_Circ& circ, const TopoDS_Vertex& V,
                           TopoDS_Edge& e) -> bool {
        if ((size_t)ci < geom.size() && geom[(size_t)ci].collapsed) {
            if (geom[(size_t)ci].edges.size() != 1) return false;
            e = TopoDS::Edge(geom[(size_t)ci].edges[0].Oriented(TopAbs_FORWARD));
            return edgeSpansFullCircle(e);
        }
        e = makeFullCircle(circ, V);
        return !e.IsNull();
    };
    if (!takeFullCap(cL, circL, verts[(size_t)vL], eL) ||
        !takeFullCap(cH, circH, verts[(size_t)vH], eH)) {
        emit(warn, "seamed360: cap is not a full 2π circle — try TwoHalves");
        return false;
    }

    TopoDS_Edge eSeam;
    try {
        gp_Lin lin(pntOf(mv, vL), r.ax.Direction());
        BRepBuilderAPI_MakeEdge ms(lin, verts[(size_t)vL], verts[(size_t)vH]);
        if (!ms.IsDone()) {
            emit(warn, "seamed360: seam MakeEdge failed");
            return false;
        }
        eSeam = ms.Edge();
    } catch (const Standard_Failure&) {
        emit(warn, "seamed360: seam threw");
        return false;
    }

    // Pcurves: caps are V-isos; seam is a U-iso at 0 and 2π.
    try {
        BRep_Builder B;
        Handle(Geom2d_Line) pcL = new Geom2d_Line(gp_Pnt2d(0.0, r.vMin), gp_Dir2d(1.0, 0.0));
        Handle(Geom2d_Line) pcH = new Geom2d_Line(gp_Pnt2d(0.0, r.vMax), gp_Dir2d(1.0, 0.0));
        Handle(Geom2d_Line) pcS0 = new Geom2d_Line(gp_Pnt2d(0.0, r.vMin), gp_Dir2d(0.0, 1.0));
        Handle(Geom2d_Line) pcS1 = new Geom2d_Line(gp_Pnt2d(2.0 * kPi, r.vMin), gp_Dir2d(0.0, 1.0));
        B.UpdateEdge(eL, pcL, surf, TopLoc_Location(), sewTol);
        B.Range(eL, 0.0, 2.0 * kPi);
        B.UpdateEdge(eH, pcH, surf, TopLoc_Location(), sewTol);
        B.Range(eH, 0.0, 2.0 * kPi);
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(eSeam, f3, l3);
        B.UpdateEdge(eSeam, pcS0, pcS1, surf, TopLoc_Location(), sewTol);
        if (!c3.IsNull()) B.Range(eSeam, f3, l3);
        eSeam.Closed(Standard_False);

        TopoDS_Wire w;
        B.MakeWire(w);
        B.Add(w, eSeam);
        B.Add(w, eH);
        B.Add(w, TopoDS::Edge(eSeam.Reversed()));
        B.Add(w, TopoDS::Edge(eL.Reversed()));
        w.Closed(Standard_True);

        BRepBuilderAPI_MakeFace mf(surf, w, Standard_True);
        TopoDS_Face got;
        if (mf.IsDone()) got = mf.Face();

        if (got.IsNull()) {
            // Spec ladder rung 1: UV-box + ShapeFix_Face(FixMissingSeamMode=1).
            BRepBuilderAPI_MakeFace box(asCyl(r), 0.0, 2.0 * kPi, r.vMin, r.vMax);
            if (!box.IsDone()) {
                emit(warn, "seamed360: UV-box MakeFace not done");
                return false;
            }
            ShapeFix_Face sff(box.Face());
            sff.FixMissingSeamMode() = 1;
            sff.FixAddNaturalBoundMode() = 0;
            sff.Add(w);
            for (const auto& iw : inners) sff.Add(iw);
            sff.Perform();
            TopoDS_Shape res = sff.Result();
            if (res.IsNull()) res = sff.Face();
            int nF = 0;
            for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                nF++;
                got = TopoDS::Face(fx.Current());
            }
            if (nF != 1) {
                emit(warn, "seamed360: ShapeFix produced " + std::to_string(nF) + " faces");
                return false;
            }
        } else {
            ShapeFix_Face sff(got);
            sff.FixMissingSeamMode() = 1;
            sff.FixAddNaturalBoundMode() = 0;
            for (const auto& iw : inners) sff.Add(iw);
            sff.Perform();
            TopoDS_Shape res = sff.Result();
            if (!res.IsNull()) {
                int nF = 0;
                TopoDS_Face g2;
                for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                    nF++;
                    g2 = TopoDS::Face(fx.Current());
                }
                if (nF == 1) got = g2;
            }
        }

        if (got.IsNull()) return false;
        if (!r.outwardNormal && got.Orientation() == TopAbs_FORWARD) got.Reverse();
        else if (r.outwardNormal && got.Orientation() == TopAbs_REVERSED) got.Reverse();
        addPcurvesOnFace(got, sewTol, true);
        if (!faceIsValid(got)) {
            emit(warn, "seamed360: BRepCheck invalid on seamed face");
            return false;
        }
        // Publish the cap circles so neighbour planes share the same TShape.
        geom[(size_t)cL].collapsed = true;
        geom[(size_t)cL].edges = {eL};
        collapsed[(size_t)cL] = 1;
        geom[(size_t)cH].collapsed = true;
        geom[(size_t)cH].edges = {eH};
        collapsed[(size_t)cH] = 1;
        outF = got;
        return true;
    } catch (const Standard_Failure& e) {
        emit(warn, std::string("seamed360: ") + (e.GetMessageString() ? e.GetMessageString() : "throw"));
        return false;
    }
}

// F9: rung 2 fires when rung 1 cannot take a single 2π cap (split / partial
// collapsed circle) and F3 still passes. nSides=3 still fail F3 on both rungs
// and explode (r8 / r1_success).
bool tryTwoHalves(const Region& r, const MeshView& mv, const std::vector<TopoDS_Vertex>& verts,
                  const Loop* capL, const Loop* capH, const RegionSet& rs,
                  std::vector<ChainGeom>& geom, std::vector<char>& collapsed, double sewTol,
                  std::vector<TopoDS_Face>& faces) {
    if (!capL || !capH) return false;
    if (capL->chainIdx.empty() || capH->chainIdx.empty()) return false;
    int cL = capL->chainIdx.front();
    int cH = capH->chainIdx.front();
    if (cL < 0 || cH < 0 || (size_t)cL >= rs.chains.size() || (size_t)cH >= rs.chains.size())
        return false;

    int vL0 = vertexClosestToUOnLoop(mv, r, *capL, rs, 0.0);
    int vLpi = vertexClosestToUOnLoop(mv, r, *capL, rs, kPi);
    int vH0 = vertexClosestToUOnLoop(mv, r, *capH, rs, 0.0);
    int vHpi = vertexClosestToUOnLoop(mv, r, *capH, rs, kPi);
    if (vL0 < 0 || vLpi < 0 || vH0 < 0 || vHpi < 0) return false;
    if (vL0 == vLpi || vH0 == vHpi) return false;
    if (std::fabs(r.vMax - r.vMin) < Precision::Confusion()) return false;
    if ((size_t)vL0 >= verts.size() || (size_t)vLpi >= verts.size() ||
        (size_t)vH0 >= verts.size() || (size_t)vHpi >= verts.size())
        return false;

    gp_Ax2 ax2(r.ax.Location(), r.ax.Direction(), r.ax.XDirection());
    gp_Circ circL(gp_Ax2(r.ax.Location().Translated(gp_Vec(r.ax.Direction()) * r.vMin),
                         r.ax.Direction(), r.ax.XDirection()),
                  r.radius);
    gp_Circ circH(gp_Ax2(r.ax.Location().Translated(gp_Vec(r.ax.Direction()) * r.vMax),
                         r.ax.Direction(), r.ax.XDirection()),
                  r.radius);
    (void)ax2;

    gp_Pnt midL0(ElCLib::Value(kPi * 0.5, circL));
    gp_Pnt midL1(ElCLib::Value(kPi * 1.5, circL));
    gp_Pnt midH0(ElCLib::Value(kPi * 0.5, circH));
    gp_Pnt midH1(ElCLib::Value(kPi * 1.5, circH));

    TopoDS_Edge aL0 = makeArc(circL, verts[(size_t)vL0], verts[(size_t)vLpi], midL0);
    TopoDS_Edge aL1 = makeArc(circL, verts[(size_t)vLpi], verts[(size_t)vL0], midL1);
    TopoDS_Edge aH0 = makeArc(circH, verts[(size_t)vH0], verts[(size_t)vHpi], midH0);
    TopoDS_Edge aH1 = makeArc(circH, verts[(size_t)vHpi], verts[(size_t)vH0], midH1);
    TopoDS_Edge e0, ePi;
    try {
        BRepBuilderAPI_MakeEdge m0(verts[(size_t)vL0], verts[(size_t)vH0]);
        BRepBuilderAPI_MakeEdge mPi(verts[(size_t)vLpi], verts[(size_t)vHpi]);
        if (!m0.IsDone() || !mPi.IsDone()) return false;
        e0 = m0.Edge();
        ePi = mPi.Edge();
    } catch (const Standard_Failure&) {
        return false;
    }
    if (aL0.IsNull() || aL1.IsNull() || aH0.IsNull() || aH1.IsNull()) return false;

    auto oneHalf = [&](const TopoDS_Edge& genA, const TopoDS_Edge& arcH, const TopoDS_Edge& genB,
                       const TopoDS_Edge& arcL, double u0, double u1, TopoDS_Face& f) -> bool {
        (void)u0;
        (void)u1;
        try {
            BRepBuilderAPI_MakeWire mw;
            mw.Add(genA);
            mw.Add(arcH);
            mw.Add(TopoDS::Edge(genB.Reversed()));
            mw.Add(TopoDS::Edge(arcL.Reversed()));
            if (!mw.IsDone()) return false;
            TopoDS_Wire ow = mw.Wire();
            Handle(Geom_CylindricalSurface) surf = new Geom_CylindricalSurface(asCyl(r));
            bindCylPCurves(ow, surf, r, sewTol);
            BRepBuilderAPI_MakeFace mf(surf, ow, Standard_True);
            if (!mf.IsDone()) return false;
            f = mf.Face();
            setFaceOutward(f, r.outwardNormal);
            addPcurvesOnFace(f, sewTol, false);
            return !f.IsNull();
        } catch (const Standard_Failure&) {
            return false;
        }
    };

    TopoDS_Face f0, f1;
    if (!oneHalf(e0, aH0, ePi, aL0, 0.0, kPi, f0)) return false;
    if (!oneHalf(ePi, aH1, e0, aL1, kPi, 2.0 * kPi, f1)) return false;
    auto fixHalf = [&](TopoDS_Face& f) {
        if (faceIsValid(f)) return;
        try {
            ShapeFix_Face sff(f);
            sff.FixOrientationMode() = 1;
            sff.FixAddNaturalBoundMode() = 0;
            sff.FixMissingSeamMode() = 0;
            sff.Perform();
            TopoDS_Shape res = sff.Result();
            if (res.IsNull()) res = sff.Face();
            for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                f = TopoDS::Face(fx.Current());
                break;
            }
        } catch (const Standard_Failure&) {
        }
        setFaceOutward(f, r.outwardNormal);
        ensureFaceValid(f);
    };
    fixHalf(f0);
    fixHalf(f1);
    if (!faceIsValid(f0) || !faceIsValid(f1)) return false;

    // Split the cap chains to the two arcs so neighbour planar caps share them.
    geom[(size_t)cL].collapsed = true;
    geom[(size_t)cL].edges = {aL0, aL1};
    geom[(size_t)cH].collapsed = true;
    geom[(size_t)cH].edges = {aH0, aH1};
    collapsed[(size_t)cL] = 1;
    collapsed[(size_t)cH] = 1;

    faces.push_back(f0);
    faces.push_back(f1);
    return true;
}

void uncollapseRegionChains(int rid, RegionSet& rs, std::vector<char>& collapsed,
                            std::vector<ChainGeom>& geom) {
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        if (!chainTouches(rs.chains[ci], rid)) continue;
        collapsed[ci] = 0;
        geom[ci].collapsed = false;
        geom[ci].edges.clear();
    }
}

void orientFaceWalk(std::vector<TopoDS_Face>& faces) {
    // Adjacent faces must use opposite orientations of each shared TShape.
    // Reverse() only, seeded from face 0 (J5: no shell-fix / sew / boolean).
    if (faces.size() < 2) return;
    struct Pair {
        int a, b;
        bool aFwd, bFwd;
    };
    auto pairs = [&]() {
        struct Rec {
            int fi;
            bool fwd;
        };
        std::vector<std::pair<TopoDS_Shape, std::vector<Rec>>> acc;
        auto find = [&](const TopoDS_Shape& e) -> size_t {
            for (size_t i = 0; i < acc.size(); i++)
                if (acc[i].first.IsSame(e)) return i;
            acc.push_back({e, {}});
            return acc.size() - 1;
        };
        for (int fi = 0; fi < (int)faces.size(); fi++) {
            for (TopExp_Explorer ex(faces[(size_t)fi], TopAbs_EDGE); ex.More(); ex.Next()) {
                const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
                acc[find(e)].second.push_back({fi, e.Orientation() == TopAbs_FORWARD});
            }
        }
        std::vector<Pair> out;
        for (auto& kv : acc) {
            if (kv.second.size() != 2) continue;
            if (kv.second[0].fi == kv.second[1].fi) continue;  // seam (same face twice)
            out.push_back({kv.second[0].fi, kv.second[1].fi, kv.second[0].fwd, kv.second[1].fwd});
        }
        return out;
    };

    std::vector<char> seen(faces.size(), 0);
    std::vector<int> q{0};
    seen[0] = 1;
    auto ps = pairs();
    auto flipFace = [&](int fj) {
        faces[(size_t)fj].Reverse();
        for (Pair& p : ps) {
            if (p.a == fj) p.aFwd = !p.aFwd;
            if (p.b == fj) p.bFwd = !p.bFwd;
        }
    };
    for (size_t qi = 0; qi < q.size(); qi++) {
        int fi = q[qi];
        for (const Pair& p : ps) {
            int fj = -1;
            bool agree = false;
            if (p.a == fi) {
                fj = p.b;
                agree = p.aFwd == p.bFwd;
            } else if (p.b == fi) {
                fj = p.a;
                agree = p.aFwd == p.bFwd;
            }
            if (fj < 0 || seen[(size_t)fj]) continue;
            seen[(size_t)fj] = 1;
            q.push_back(fj);
            if (agree) flipFace(fj);
        }
        if (qi + 1 == q.size()) {
            for (int k = 0; k < (int)faces.size(); k++)
                if (!seen[(size_t)k]) {
                    seen[(size_t)k] = 1;
                    q.push_back(k);
                    break;
                }
        }
    }
}

double mixedEdgeDeviation(const MeshView& mv, const TopoDS_Vertex& v1, const TopoDS_Vertex& v2,
                          const Region* analytic) {
    if (!analytic) return 0.0;
    gp_Pnt a = BRep_Tool::Pnt(v1), b = BRep_Tool::Pnt(v2);
    gp_Pnt m(0.5 * (a.XYZ() + b.XYZ()));
    auto dev = [&](const gp_Pnt& p) {
        if (analytic->type == SurfType::Plane) return asPlane(*analytic).Distance(p);
        if (analytic->type == SurfType::Cylinder) {
            gp_Cylinder cyl = asCyl(*analytic);
            gp_Vec v(cyl.Location(), p);
            gp_Vec cr = gp_Vec(cyl.Axis().Direction()).Crossed(v);
            return std::fabs(cr.Magnitude() - cyl.Radius());
        }
        return 0.0;
    };
    return std::max(dev(m), std::max(dev(a), dev(b)));
}

}  // namespace

void fitAnalyticTolerances(const TopoDS_Shape& shape) {
    // Sibling of fitPlanarTolerances (src/stl2step.cpp). Cylinder radial
    // deviation d = |‖axis × (p − loc)‖ − R|, vertices AND edge midpoints
    // (a chord's max deviation from a cylinder is at its midpoint). Same
    // ×1.001 + Precision::Confusion() policy. P3 adds the guarded call.
    BRep_Builder B;
    for (TopExp_Explorer fx(shape, TopAbs_FACE); fx.More(); fx.Next()) {
        const TopoDS_Face& f = TopoDS::Face(fx.Current());
        BRepAdaptor_Surface s(f, Standard_False);
        if (s.GetType() != GeomAbs_Cylinder) continue;
        gp_Cylinder cyl = s.Cylinder();
        gp_Pnt loc = cyl.Location();
        gp_Dir ad = cyl.Axis().Direction();
        double R = cyl.Radius();
        auto radDev = [&](const gp_Pnt& p) {
            gp_Vec v(loc, p);
            gp_Vec cross = gp_Vec(ad).Crossed(v);
            return std::fabs(cross.Magnitude() - R);
        };
        for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
            const TopoDS_Vertex& v = TopoDS::Vertex(vx.Current());
            double d = radDev(BRep_Tool::Pnt(v));
            if (d > BRep_Tool::Tolerance(v))
                B.UpdateVertex(v, d * 1.001 + Precision::Confusion());
        }
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
            double d = 0;
            for (TopExp_Explorer vx(e, TopAbs_VERTEX); vx.More(); vx.Next())
                d = std::max(d, radDev(BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()))));
            Standard_Real fpar = 0, lpar = 0;
            Handle(Geom_Curve) c = BRep_Tool::Curve(e, fpar, lpar);
            if (!c.IsNull()) {
                gp_Pnt mid = c->Value(0.5 * (fpar + lpar));
                d = std::max(d, radDev(mid));
            }
            if (d > BRep_Tool::Tolerance(e))
                B.UpdateEdge(e, d * 1.001 + Precision::Confusion());
        }
    }
}

bool buildFaces(const MeshView& mv, RegionSet& rs, const std::vector<TopoDS_Vertex>& verts,
                std::vector<TopoDS_Face>& out, WarnFn warn) {
    out.clear();
    try {
        if (mv.nTri == 0) return false;
        if (verts.size() < mv.nVtx) return false;
        if (!regionSetConsistent(mv, rs, verts)) {
            for (Region& r : rs.regions) {
                r.reject = Reject::ChainUnstable;
                r.builtAs = BuiltAs::NotBuilt;
            }
            return false;
        }

        const double sewTol = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
        const bool wasClosed = meshComponentClosed(mv);

        int maxId = -1;
        for (const Region& r : rs.regions) maxId = std::max(maxId, r.id);
        std::vector<char> exploded((size_t)std::max(0, maxId) + 1, 0);

        std::vector<TopoDS_Edge> meshE(mv.nEdge);
        std::vector<char> edgeOk(mv.nEdge, 0);
        for (size_t i = 0; i < mv.nEdge; i++) {
            int a = mv.compEdges[i].first, b = mv.compEdges[i].second;
            if (a < 0 || b < 0 || (size_t)a >= verts.size() || (size_t)b >= verts.size()) continue;
            try {
                BRepBuilderAPI_MakeEdge me(verts[(size_t)a], verts[(size_t)b]);
                if (me.IsDone()) {
                    meshE[i] = me.Edge();
                    edgeOk[i] = 1;
                }
            } catch (const Standard_Failure&) {
            }
        }

        std::vector<char> collapsed(rs.chains.size(), 0);
        std::vector<ChainGeom> geom(rs.chains.size());

        auto rebuildCollapsed = [&]() {
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                collapsed[ci] = 0;
                geom[ci] = ChainGeom{};
                const BoundaryChain& ch = rs.chains[ci];
                const Region* A = regionById(rs, ch.regA);
                const Region* B = regionById(rs, ch.regB);
                if (A && regionExploded(exploded, A->id)) A = nullptr;
                if (B && regionExploded(exploded, B->id)) B = nullptr;
                if (!isAnalytic(A) || !isAnalytic(B)) {
                    // mixed analytic|faceted, or island|island: polyline verbatim
                    if ((isAnalytic(A) && !isAnalytic(B)) || (!isAnalytic(A) && isAnalytic(B))) {
                        const Region* an = isAnalytic(A) ? A : B;
                        for (int eid : ch.meshEdges) {
                            if (eid < 0 || (size_t)eid >= mv.nEdge || !edgeOk[(size_t)eid]) continue;
                            const auto& ev = mv.compEdges[eid];
                            double d = mixedEdgeDeviation(mv, verts[(size_t)ev.first],
                                                          verts[(size_t)ev.second], an);
                            rs.stats.maxEdgeTol = std::max(rs.stats.maxEdgeTol, d);
                        }
                    }
                    continue;
                }
                // analytic | analytic
                AnalyticCurve curve;
                if (ch.tangent) {
                    if (A->type == SurfType::Cylinder && B->type == SurfType::Plane)
                        curve = constructedGenerator(*A, *B);
                    else if (B->type == SurfType::Cylinder && A->type == SurfType::Plane)
                        curve = constructedGenerator(*B, *A);
                    else if (A->type == SurfType::Cylinder && B->type == SurfType::Cylinder)
                        curve = constructedCylCylGenerator(*A, *B);
                    else
                        curve = intersectSurfaces(*A, *B, mv, ch, sewTol, warn);
                } else {
                    curve = intersectSurfaces(*A, *B, mv, ch, sewTol, warn);
                }
                if (curve.kind == AnalyticCurve::None) continue;

                // Terminals: open chain = first/last meshVert; closed = seam vertex (or first).
                int ia, ib;
                bool full = ch.closedLoop;
                if (ch.closedLoop) {
                    const Region* cyl = (A->type == SurfType::Cylinder) ? A
                                       : (B->type == SurfType::Cylinder) ? B
                                                                         : nullptr;
                    int sv = cyl ? seamVertexOf(mv, *cyl, ch) : ch.meshVerts.front();
                    ia = ib = sv;
                } else {
                    if (ch.meshVerts.size() < 2) continue;
                    ia = ch.meshVerts.front();
                    ib = ch.meshVerts.back();
                }
                if (ia < 0 || ib < 0 || (size_t)ia >= verts.size() || (size_t)ib >= verts.size())
                    continue;
                // F1: closed360 cap circles must be parameterised on the region's
                // u=0 frame (XDirection is a mesh-vertex azimuth). IntAna picks
                // its own XDirection; leaving it makes Range(0,2π) desync the
                // 3d curve from the u:0→2π pcurve trySeamed360 then writes.
                if (curve.kind == AnalyticCurve::Circ) {
                    const Region* cylR = (A->type == SurfType::Cylinder) ? A
                                       : (B->type == SurfType::Cylinder) ? B
                                                                         : nullptr;
                    if (cylR) {
                        curve.circ = gp_Circ(
                            gp_Ax2(curve.circ.Location(), cylR->ax.Direction(),
                                   cylR->ax.XDirection()),
                            curve.circ.Radius());
                    }
                }
                TopoDS_Edge e;
                if (!full && curve.kind == AnalyticCurve::Circ && ch.meshVerts.size() >= 3) {
                    // Open circular chain (partial-cyl cap): 180° MakeEdge(V,V) is
                    // ambiguous — the mesh mid-vertex selects the notch, not the
                    // complementary half (UnorientableShape / wrong volume sign).
                    gp_Pnt mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
                    e = makeArc(curve.circ, verts[(size_t)ia], verts[(size_t)ib], mid);
                } else {
                    e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
                }
                if (e.IsNull()) {
                    emit(warn, "smooth: analytic MakeEdge failed — keeping mesh polyline");
                    continue;
                }
                if (full && curve.kind == AnalyticCurve::Circ && ch.meshVerts.size() >= 2) {
                    double t0 = ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[0]));
                    double t1 = ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[1]));
                    double dt = t1 - t0;
                    while (dt < 0.0) dt += 2.0 * kPi;
                    if (dt > kPi) e.Reverse();
                }
                geom[ci].collapsed = true;
                geom[ci].edges = {e};
                collapsed[ci] = 1;
            }
        };

        auto buildOneRegion = [&](Region& r, std::vector<TopoDS_Face>& acc) -> bool {
            if (regionExploded(exploded, r.id)) return true;
            r.builtAs = BuiltAs::NotBuilt;
            if (r.type == SurfType::Cone || r.type == SurfType::Sphere || r.type == SurfType::Torus) {
                r.reject = (r.type == SurfType::Cone)     ? Reject::ConeNYI
                           : (r.type == SurfType::Sphere) ? Reject::SphereNYI
                                                          : Reject::TorusNYI;
                return false;
            }
            if (r.closed360 && r.type == SurfType::Cylinder) {
                TopoDS_Face f360;
                if (trySeamed360(r, rs, mv, verts, geom, collapsed, meshE, edgeOk, sewTol, warn,
                                 f360) &&
                    cylinderPostFitOk(r, mv, rs)) {
                    r.builtAs = BuiltAs::Seamed360;
                    acc.push_back(f360);
                    return true;
                }
                const Loop *capL = nullptr, *capH = nullptr;
                for (const Loop& lp : r.loops) {
                    if (lp.role == LoopRole::CapLow) capL = &lp;
                    if (lp.role == LoopRole::CapHigh) capH = &lp;
                }
                std::vector<TopoDS_Face> halves;
                if (tryTwoHalves(r, mv, verts, capL, capH, rs, geom, collapsed, sewTol, halves) &&
                    cylinderPostFitOk(r, mv, rs)) {
                    bool allOk = true;
                    for (auto& hf : halves)
                        if (!ensureFaceValid(hf)) allOk = false;
                    if (allOk) {
                        r.builtAs = BuiltAs::TwoHalves;
                        for (auto& f : halves) acc.push_back(f);
                        return true;
                    }
                }
                return false;
            }
            if (r.type == SurfType::Cylinder) {
                TopoDS_Face f;
                if (!buildPartialCylinder(r, rs, mv, geom, collapsed, meshE, edgeOk, sewTol, f))
                    return false;
                r.builtAs = BuiltAs::Single;
                acc.push_back(f);
                return true;
            }
            if (r.type == SurfType::Plane) {
                TopoDS_Face f;
                if (!buildPlanarFace(r, rs, mv, geom, collapsed, meshE, edgeOk, f)) return false;
                r.builtAs = BuiltAs::Single;
                acc.push_back(f);
                return true;
            }
            return false;
        };

        int rounds = 0;
        const int kMaxRounds = 2;
        bool unstable = false;
        std::vector<TopoDS_Face> built;
        while (true) {
            rebuildCollapsed();
            built.clear();
            bool anyFail = false;
            std::vector<int> failedIds;
            for (Region& r : rs.regions) {
                if (regionExploded(exploded, r.id)) continue;
                std::vector<TopoDS_Face> acc;
                if (!buildOneRegion(r, acc)) {
                    anyFail = true;
                    failedIds.push_back(r.id);
                    r.reject = Reject::FaceBuildFailed;
                    r.builtAs = BuiltAs::NotBuilt;
                } else {
                    for (auto& f : acc) built.push_back(f);
                }
            }
            if (!anyFail) break;
            if (rounds >= kMaxRounds) {
                unstable = true;
                break;
            }
            // R1: explode failed regions, un-collapse their chains, rebuild neighbours.
            for (int id : failedIds) {
                if (id >= 0 && (size_t)id < exploded.size()) exploded[(size_t)id] = 1;
                Region* rr = regionByIdMut(rs, id);
                if (rr) {
                    rr->builtAs = BuiltAs::ExplodedToFacets;
                    rr->reject = Reject::FaceBuildFailed;
                }
                uncollapseRegionChains(id, rs, collapsed, geom);
            }
            rounds++;
        }

        if (unstable) {
            for (Region& r : rs.regions) {
                if (r.reject == Reject::FaceBuildFailed) r.reject = Reject::ChainUnstable;
            }
            out.clear();
            return false;
        }

        // Facet islands + exploded-region triangles.
        for (size_t k = 0; k < mv.nTri; k++) {
            int rid = (k < rs.triRegion.size()) ? rs.triRegion[k] : -1;
            int iid = (k < rs.triIsland.size()) ? rs.triIsland[k] : -1;
            bool exp = (rid >= 0 && regionExploded(exploded, rid));
            if (iid < 0 && !exp) continue;
            TopoDS_Face f = makeFacet(mv, verts, meshE, edgeOk, k);
            if (!f.IsNull()) built.push_back(f);
        }

        if (built.empty()) return false;

        for (auto& f : built) {
            try {
                // Only repair faces that are currently BRepCheck-invalid.
                // Skipping this on open meshes left s09's mixed plane
                // UnorientableShape (wire order); running it on already-valid
                // faces (esp. cylinders / island facets) copies edges and
                // undoes hole Reverse / breaks TShape sharing.
                {
                    BRepCheck_Analyzer an(f, Standard_False);
                    if (an.IsValid()) continue;
                    BRepAdaptor_Surface sa(f, Standard_False);
                    int nE = 0;
                    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) nE++;
                    // Triangle facets already share meshE / verts[]. ShapeFix_Face
                    // copies them and is the R1 explode J2-twin source.
                    if (nE == 3 && sa.GetType() == GeomAbs_Plane) continue;
                }
                ShapeFix_Face sff(f);
                sff.FixOrientationMode() = 1;
                sff.FixAddNaturalBoundMode() = 0;
                sff.FixMissingSeamMode() = 0;
                if (sff.Perform()) {
                    TopoDS_Shape res = sff.Result();
                    if (!res.IsNull()) {
                        for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                            f = TopoDS::Face(fx.Current());
                            break;
                        }
                    } else if (!sff.Face().IsNull()) {
                        f = sff.Face();
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }

        orientFaceWalk(built);

        // Re-assert outwardNormal on *partial* (not closed360) cylinders after
        // ShapeFix/walk. Seamed 360 faces already match AC4; touching them
        // regresses full_360_hole.
        {
            bool partialHole = false, partialBoss = false;
            for (const Region& rg : rs.regions) {
                if (rg.type != SurfType::Cylinder || rg.closed360) continue;
                if (regionExploded(exploded, rg.id)) continue;
                if (rg.outwardNormal) partialBoss = true;
                else partialHole = true;
            }
            if (partialHole != partialBoss) {
                bool wantOut = partialBoss;
                for (auto& f : built) {
                    BRepAdaptor_Surface sa(f, Standard_False);
                    if (sa.GetType() != GeomAbs_Cylinder) continue;
                    // Skip faces that wrap a full 2π (seamed 360).
                    double umin = 0, umax = 0, vmin = 0, vmax = 0;
                    BRepTools::UVBounds(f, umin, umax, vmin, vmax);
                    if (umax - umin > 1.5 * kPi) continue;
                    setFaceOutward(f, wantOut);
                }
            }
        }

        // J4: assemble a shell and SameParameter(forced).
        BRep_Builder B;
        TopoDS_Shell sh;
        B.MakeShell(sh);
        for (auto& f : built) B.Add(sh, f);
        try {
            BRepLib::SameParameter(sh, sewTol, /*forced=*/Standard_True);
            bool needRelax = false;
            for (const Region& rg : rs.regions) {
                if (rg.type == SurfType::Cylinder && !rg.closed360 && !rg.outwardNormal &&
                    !regionExploded(exploded, rg.id))
                    needRelax = true;
            }
            if (needRelax) {
                BRep_Builder bb;
                bool any = false;
                for (TopExp_Explorer ex(sh, TopAbs_EDGE); ex.More(); ex.Next()) {
                    const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
                    if (BRep_Tool::SameParameter(e)) continue;
                    any = true;
                    bb.UpdateEdge(e, std::max(BRep_Tool::Tolerance(e), sewTol) * 10.0
                                         + Precision::Confusion());
                }
                if (any) BRepLib::SameParameter(sh, sewTol * 10.0, /*forced=*/Standard_True);
            }
        } catch (const Standard_Failure&) {
            return false;
        }

        // J6
        bool shClosed = BRep_Tool::IsClosed(sh);
        sh.Closed(shClosed ? Standard_True : Standard_False);
        if (wasClosed && !shClosed) return false;

        fitAnalyticTolerances(sh);

        out = std::move(built);
        return true;
    } catch (const Standard_Failure&) {
        out.clear();
        return false;
    } catch (const std::exception&) {
        out.clear();
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
