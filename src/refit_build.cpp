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
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
// BRepCheck_ListIteratorOfListOfStatus.hxx and
// TopTools_ListIteratorOfListOfShape.hxx are typedef-only headers removed in
// OCCT 8.0; iterate via ListOfStatus::Iterator / ListOfShape::Iterator.
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <ElCLib.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_CurveType.hxx>
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
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax1.hxx>
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
#include <GeomAPI.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomProjLib.hxx>

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

Handle(Geom_CylindricalSurface) cylSurfaceForRegion(const Region& r);

gp_Cylinder cylForIntersect(const Region& r) {
    if (r.type != SurfType::Cylinder) return asCyl(r);
    return cylSurfaceForRegion(r)->Cylinder();
}

// H1b / D3S-3: build-scope E′ demotion side table. Consulted by isAnalytic so
// all five call sites see one truth. Never mutates Region::type. thread_local
// because convert() fans out components.
thread_local std::vector<char> gEPrimeDemoted;

bool isAnalytic(const Region* r) {
    if (!r) return false;
    if (r->type != SurfType::Plane && r->type != SurfType::Cylinder) return false;
    if (r->id >= 0 && (size_t)r->id < gEPrimeDemoted.size() && gEPrimeDemoted[(size_t)r->id])
        return false;
    return true;
}

void installEPrimeDemotion(const MeshView& mv, const RegionSet& rs) {
    gEPrimeDemoted.clear();
    int maxId = -1;
    std::vector<int> cylIdx;
    cylIdx.reserve(rs.regions.size());
    for (size_t i = 0; i < rs.regions.size(); ++i) {
        const Region& r = rs.regions[i];
        maxId = std::max(maxId, r.id);
        if (r.type == SurfType::Cylinder) cylIdx.push_back((int)i);
    }
    if (maxId >= 0) gEPrimeDemoted.assign((size_t)maxId + 1, 0);
    if (cylIdx.empty()) return;

    // RULE 4.2a, recomputed locally (Stage-P TUs stay frozen).
    const double weld = mv.weldTol > 0.0 ? mv.weldTol : 0.0;
    const double diag = std::isfinite(mv.diag) && mv.diag > 0.0 ? mv.diag : 0.0;
    const double tauSurf = std::max(5e-5, std::max(4.0 * weld, 1e-6 * diag));
    double hMin = 0.0;
    for (int i : cylIdx) {
        const Region& r = rs.regions[(size_t)i];
        const double h = std::abs(r.vMax - r.vMin);
        if (!(h > 0.0) || !std::isfinite(h)) continue;
        if (hMin <= 0.0 || h < hMin) hMin = h;
    }
    double tauAx = 1e-6;
    if (hMin > 0.0) tauAx = std::max(1e-6, 2.0 * tauSurf / hMin);

    auto axisSin = [](const gp_Dir& a, const gp_Dir& b) {
        return gp_Vec(a).Crossed(gp_Vec(b)).Magnitude();
    };

    const int n = (int)cylIdx.size();
    std::vector<int> parent((size_t)n);
    for (int i = 0; i < n; ++i) parent[(size_t)i] = i;
    auto find = [&](int x) {
        int r = x;
        while (parent[(size_t)r] != r) r = parent[(size_t)r];
        while (parent[(size_t)x] != r) {
            const int nxt = parent[(size_t)x];
            parent[(size_t)x] = r;
            x = nxt;
        }
        return r;
    };
    auto unite = [&](int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (a < b) parent[(size_t)b] = a;
        else parent[(size_t)a] = b;
    };
    for (int i = 0; i < n; ++i) {
        const gp_Dir di = rs.regions[(size_t)cylIdx[i]].ax.Direction();
        for (int j = i + 1; j < n; ++j) {
            const gp_Dir dj = rs.regions[(size_t)cylIdx[j]].ax.Direction();
            if (axisSin(di, dj) <= tauAx) unite(i, j);
        }
    }

    struct Cl {
        int tris = 0;
        int minRid = INT_MAX;
        int minIdx = -1;
    };
    std::vector<Cl> stats((size_t)n);
    for (int i = 0; i < n; ++i) {
        const Region& r = rs.regions[(size_t)cylIdx[i]];
        const int p = find(i);
        stats[(size_t)p].tris += (int)r.tris.size();
        if (r.id < stats[(size_t)p].minRid) {
            stats[(size_t)p].minRid = r.id;
            stats[(size_t)p].minIdx = i;
        }
    }
    int best = -1;
    for (int i = 0; i < n; ++i) {
        if (find(i) != i) continue;
        if (best < 0) {
            best = i;
            continue;
        }
        if (stats[(size_t)i].tris > stats[(size_t)best].tris ||
            (stats[(size_t)i].tris == stats[(size_t)best].tris &&
             stats[(size_t)i].minRid < stats[(size_t)best].minRid))
            best = i;
    }
    if (best < 0 || stats[(size_t)best].minIdx < 0) return;
    const gp_Dir ahat =
        rs.regions[(size_t)cylIdx[stats[(size_t)best].minIdx]].ax.Direction();

    const char* ev = std::getenv("STL2STEP_EPRIME_DIAG");
    const bool diagOn = ev && ev[0] && ev[0] != '0';
    for (int i = 0; i < n; ++i) {
        const Region& r = rs.regions[(size_t)cylIdx[i]];
        const double sinVs = axisSin(r.ax.Direction(), ahat);
        const bool demote = r.radius <= 3.5 && !r.closed360 && sinVs > tauAx;
        if (demote && r.id >= 0 && (size_t)r.id < gEPrimeDemoted.size())
            gEPrimeDemoted[(size_t)r.id] = 1;
        if (diagOn)
            std::fprintf(stderr,
                         "DIAG_EPRIME rid=%d R=%.6f nTri=%zu closed360=%d "
                         "sinVsDom=%.6e demoted=%d\n",
                         r.id, r.radius, r.tris.size(), r.closed360 ? 1 : 0, sinVs,
                         demote ? 1 : 0);
    }
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
// Live nSides=4 fillets sit exactly on the limit; 0.1% slack keeps them in
// and still rejects nSides=3 (sagitta/R = 0.5).
bool cylinderPostFitOk(const Region& r, const MeshView& mv, const RegionSet& rs) {
    if (r.type != SurfType::Cylinder || r.radius <= 0.0) return true;
    const double lim = r.radius * (1.0 - std::cos(kPi / 4.0)) * 1.001 + Precision::Confusion();
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

double meshTolCap(const MeshView& mv, const Region* r) {
    double c = std::max(mv.sewTol, Precision::Confusion());
    if (r) {
        c = std::max(c, r->chordSagitta);
        if (r->radius > 0.0)
            c = std::max(c, r->radius * (1.0 - std::cos(kPi / 4.0)));
    }
    return c * 1.001 + Precision::Confusion();
}

// Same derivation as refit_segment.cpp DerivedTols::epsPlane (frozen TU).
double derivedEpsPlane(const MeshView& mv) {
    const double diag = mv.diag > 0.0 ? mv.diag : 1.0;
    const double epsMesh = std::max(std::max(mv.weldTol, 1e-4 * diag), 1e-3);
    return std::max(std::max(epsMesh, mv.sewTol), 0.02);
}

// Partial-cylinder face validity: prefer measured chord sag over the nSides≥4
// π/4 floor when nSides was counted (handle-lock rid=5: 2.14 mm vs 4.7 mm).
double partialFaceTolCap(const MeshView& mv, const Region& r) {
    double c = std::max(mv.sewTol, Precision::Confusion());
    c = std::max(c, r.maxVertexDev);
    c = std::max(c, derivedEpsPlane(mv));
    if (r.nSides > 0 && r.chordSagitta > 0.0)
        c = std::max(c, r.chordSagitta);
    else if (r.radius > 0.0)
        c = std::max(c, r.radius * (1.0 - std::cos(kPi / 4.0)));
    return c * 1.001 + Precision::Confusion();
}

// Snap budget for an accepted IntAna / constructed curve. Floor is SPEC-F2
// max(sewTol, region.maxVertexDev, epsPlane). Ceiling is pickIntAna's
// wrong-branch gate (max(sewTol*50, 1 mm)): a curve we kept is allowed
// that residual. Does not move shared verts; only bumps tolerance.
double analyticSnapCap(const MeshView& mv, const Region* A, const Region* B) {
    double c = std::max(meshTolCap(mv, A), meshTolCap(mv, B));
    if (A) c = std::max(c, A->maxVertexDev);
    if (B) c = std::max(c, B->maxVertexDev);
    c = std::max(c, derivedEpsPlane(mv));
    return c;
}

// F4: absorb sagitta in edge/vertex tolerance, but NEVER inflate a shared
// TShape past the mesh budget (adjudication F5). Beyond that the fit is
// wrong — reject the region rather than poison every adjacent face.
bool ensureFaceValid(TopoDS_Face& f, double cap) {
    if (f.IsNull()) return false;
    if (faceIsValid(f)) return true;
    if (!(cap > 0.0)) cap = Precision::Confusion();
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
            if (!(d > 0.0)) continue;
            if (d > cap) return false;
            if (d > BRep_Tool::Tolerance(vv))
                B.UpdateVertex(vv, std::min(d * 1.001 + Precision::Confusion(), cap));
        }
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
            double d = 0.0;
            for (TopExp_Explorer vx(e, TopAbs_VERTEX); vx.More(); vx.Next())
                d = std::max(d, devPnt(BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()))));
            Standard_Real fp = 0, lp = 0;
            Handle(Geom_Curve) c = BRep_Tool::Curve(e, fp, lp);
            if (!c.IsNull()) d = std::max(d, devPnt(c->Value(0.5 * (fp + lp))));
            if (d > cap) return false;
            if (d > BRep_Tool::Tolerance(e))
                B.UpdateEdge(e, std::min(d * 1.001 + Precision::Confusion(), cap));
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

void bumpVertexTol(const TopoDS_Vertex& v, double d) {
    if (v.IsNull() || !(d > 0.0)) return;
    if (d > BRep_Tool::Tolerance(v)) {
        BRep_Builder B;
        B.UpdateVertex(v, d * 1.001 + Precision::Confusion());
    }
}

void snapVertexToCurve(const TopoDS_Vertex& v, const AnalyticCurve& c, double cap) {
    if (v.IsNull()) return;
    double d = curveResidual(c, BRep_Tool::Pnt(v));
    if (!std::isfinite(d) || d > cap) return;
    bumpVertexTol(v, d);
}

gp_Circ cylinderIsoCircle(const Region& cyl, double v) {
    gp_Pnt loc = cyl.ax.Location().Translated(gp_Vec(cyl.ax.Direction()) * v);
    return gp_Circ(gp_Ax2(loc, cyl.ax.Direction(), cyl.ax.XDirection()), cyl.radius);
}

bool planePerpCylinder(const Region& pln, const Region& cyl) {
    double c = std::fabs(pln.ax.Direction().Dot(cyl.ax.Direction()));
    return c >= std::cos(3.0 * kPi / 180.0);
}

double planeVOnCylinder(const Region& pln, const Region& cyl) {
    gp_Dir n = pln.ax.Direction();
    gp_Dir a = cyl.ax.Direction();
    double na = n.Dot(a);
    if (std::fabs(na) < 1e-12) return 0.0;
    return gp_Vec(cyl.ax.Location(), pln.ax.Location()).Dot(n) / na;
}

// Ordered edges of a closed wire, rotated so the first edge starts at V.
bool wireEdgesFromVertex(const TopoDS_Wire& w, const TopoDS_Vertex& V,
                         std::vector<TopoDS_Edge>& edges) {
    edges.clear();
    if (w.IsNull() || V.IsNull()) return false;
    for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next())
        edges.push_back(TopoDS::Edge(ex.Current()));
    if (edges.empty()) return false;
    if (edges.size() == 1) {
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(edges[0], v1, v2, Standard_False);
        return v1.IsSame(V) || v2.IsSame(V);
    }
    for (size_t i = 0; i < edges.size(); i++) {
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(edges[i], v1, v2, Standard_True);
        if (v1.IsSame(V)) {
            std::rotate(edges.begin(), edges.begin() + (std::ptrdiff_t)i, edges.end());
            return true;
        }
    }
    return false;
}

double intAnaAcceptResidual(const MeshView& mv, const BoundaryChain& ch, double sewTol,
                            const Region* A, const Region* B) {
    double acceptR = std::max(sewTol * 50.0, 1.0);
    if (A && B) {
        const double fit = std::max(A->maxVertexDev, B->maxVertexDev);
        if (ch.meshVerts.size() <= 3)
            acceptR = std::max(acceptR, fit * 3.0 + derivedEpsPlane(mv));
        else if (ch.meshVerts.size() <= 20)
            acceptR = std::max(acceptR, fit * 2.0 + derivedEpsPlane(mv));
        // Skew/near-tangent cyl|cyl: IntAna conics land off the mesh band on
        // coarse exports — widen acceptance vs the fitted residual.
        if (A->type == SurfType::Cylinder && B->type == SurfType::Cylinder) {
            acceptR = std::max(acceptR, fit * 4.0 + derivedEpsPlane(mv));
            if (ch.meshVerts.size() <= 3)
                acceptR = std::max(acceptR, fit * 5.0 + derivedEpsPlane(mv) * 2.0);
        }
        if ((A->type == SurfType::Plane && B->type == SurfType::Cylinder) ||
            (A->type == SurfType::Cylinder && B->type == SurfType::Plane)) {
            if (ch.meshVerts.size() <= 3)
                acceptR = std::max(acceptR, fit * 5.0 + derivedEpsPlane(mv) * 2.0);
        }
    }
    return acceptR;
}

AnalyticCurve pickIntAna(const IntAna_QuadQuadGeo& iq, const MeshView& mv, const BoundaryChain& ch,
                         double sewTol, double acceptResidual = -1.0) {
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
    if (best.kind != AnalyticCurve::None) {
        const double acceptR =
            acceptResidual > 0.0 ? acceptResidual : std::max(sewTol * 50.0, 1.0);
        if (bestR > acceptR) best.kind = AnalyticCurve::None;
    }
    return best;
}

bool diagP2Enabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_P2_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

bool diagJ6Enabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_J6_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

bool diagCoverEnabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_COVER_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// Default-off B0/B0.1 provenance. Thread-local so concurrent convert() calls
// do not share stamps. Written only when STL2STEP_COVER_DIAG is set.
thread_local int gDiagMeshEGen = 0;
thread_local std::unordered_map<const void*, const char*> gDiagEdgeSite;
thread_local std::unordered_map<const void*, int> gDiagEdgeGen;

const void* diagTShapePtr(const TopoDS_Shape& s) {
    if (s.IsNull()) return nullptr;
    return (const void*)s.TShape().get();
}

void diagNoteEdge(const TopoDS_Edge& e, const char* site) {
    if (!diagCoverEnabled() || e.IsNull() || !site) return;
    const void* p = diagTShapePtr(e);
    if (!p) return;
    gDiagEdgeSite[p] = site;
    auto it = gDiagEdgeGen.find(p);
    if (it == gDiagEdgeGen.end()) gDiagEdgeGen[p] = gDiagMeshEGen;
}

void diagMaybeDegenEmit(const TopoDS_Edge& e, const char* site, int ci, int ridA, int ridB,
                        int extra) {
    if (!diagCoverEnabled() || e.IsNull()) return;
    diagNoteEdge(e, site);
    TopoDS_Vertex va, vb;
    TopExp::Vertices(e, va, vb, Standard_True);
    gp_Pnt pa = BRep_Tool::Pnt(va), pb = BRep_Tool::Pnt(vb);
    const double len = pa.Distance(pb);
    const int flagged = BRep_Tool::Degenerated(e) ? 1 : 0;
    if (len > 1e-4 && !flagged) return;
    std::fprintf(stderr,
                 "DIAG_DEGEN_EMIT site=%s ci=%d ridA=%d ridB=%d extra=%d "
                 "len=%.6f degen=%d sameV=%d tshape=%p "
                 "pa=(%.6f,%.6f,%.6f) pb=(%.6f,%.6f,%.6f)\n",
                 site, ci, ridA, ridB, extra, len, flagged, va.IsSame(vb) ? 1 : 0,
                 diagTShapePtr(e), pa.X(), pa.Y(), pa.Z(), pb.X(), pb.Y(), pb.Z());
}

void diagJ6FreeEdges(const MeshView& mv, const TopoDS_Shell& sh,
                     const std::vector<TopoDS_Face>& built, const std::vector<int>& builtRid,
                     const RegionSet& rs, const std::vector<char>& collapsed) {
    try {
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            int idx = faceMap.FindIndex(built[i]);
            if (idx > 0 && i < builtRid.size()) faceRid[(size_t)idx] = builtRid[i];
        }
        int nFree = 0;
        for (int i = 1; i <= anc.Extent(); i++) {
            if (anc(i).Extent() >= 2) continue;
            const TopoDS_Edge& e = TopoDS::Edge(anc.FindKey(i));
            TopoDS_Vertex va, vb;
            TopExp::Vertices(e, va, vb, Standard_True);
            gp_Pnt pa = BRep_Tool::Pnt(va), pb = BRep_Tool::Pnt(vb);
            int frid = -1;
            if (anc(i).Extent() == 1) {
                const TopoDS_Face& f = TopoDS::Face(anc(i).First());
                int fi = faceMap.FindIndex(f);
                if (fi > 0) frid = faceRid[(size_t)fi];
            }
            int ciMatch = -1;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                if (ci >= collapsed.size() || !collapsed[ci]) continue;
                const BoundaryChain& ch = rs.chains[ci];
                if (ch.meshVerts.size() < 2) continue;
                gp_Pnt p0 = pntOf(mv, ch.meshVerts.front());
                gp_Pnt p1 = pntOf(mv, ch.meshVerts.back());
                if (pa.Distance(p0) < 0.05 && pb.Distance(p1) < 0.05 ||
                    pa.Distance(p1) < 0.05 && pb.Distance(p0) < 0.05)
                    ciMatch = (int)ci;
            }
            std::fprintf(stderr,
                         "DIAG_J6 freeE#%d faceRid=%d ci=%d pa=(%.3f,%.3f,%.3f) "
                         "pb=(%.3f,%.3f,%.3f) len=%.4f\n",
                         nFree, frid, ciMatch, pa.X(), pa.Y(), pa.Z(), pb.X(), pb.Y(), pb.Z(),
                         pa.Distance(pb));
            nFree++;
        }
        int analyticPoly = 0, analyticCol = 0;
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            const BoundaryChain& ch = rs.chains[ci];
            const Region* A = regionById(rs, ch.regA);
            const Region* B = regionById(rs, ch.regB);
            if (!isAnalytic(A) || !isAnalytic(B)) continue;
            if (ci < collapsed.size() && collapsed[ci]) analyticCol++;
            else analyticPoly++;
        }
        std::fprintf(stderr, "DIAG_J6 summary freeEdges=%d analyticCollapsed=%d "
                             "analyticPolyline=%d\n",
                     nFree, analyticCol, analyticPoly);
    } catch (const Standard_Failure&) {
    }
}

// B0/B0.1: default-off coverage / free-edge attribution. STL2STEP_COVER_DIAG=1.
// H1: twin predicate requires len > tol (no zero-length self-match); fourth
//     class degenerate-zero-length.
// H2: coverage keyed on builtRid membership, not Region::builtAs.
// H3: emit TShape identity, meshE generation, construction site per row.
void diagCoverDump(const MeshView& mv, const RegionSet& rs, const std::vector<char>& exploded,
                   const std::vector<int>& builtRid, const std::vector<TopoDS_Face>& built,
                   const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& collapsed,
                   const std::vector<ChainGeom>& geom, const TopoDS_Shell& sh, int j6Pass) {
    if (!diagCoverEnabled()) return;
    try {
        const size_t nTri = mv.nTri;
        int maxRid = -1;
        for (int id : builtRid) maxRid = std::max(maxRid, id);
        for (const Region& r : rs.regions) maxRid = std::max(maxRid, r.id);
        std::vector<char> shipped((size_t)std::max(0, maxRid) + 1, 0);
        int nShippedRids = 0;
        for (int id : builtRid) {
            if (id >= 0 && (size_t)id < shipped.size() && !shipped[(size_t)id]) {
                shipped[(size_t)id] = 1;
                nShippedRids++;
            }
        }

        auto triClass = [&](size_t k) -> int {
            // 0 analytic, 1 island, 2 exploded, 3 uncovered
            int rid = (k < rs.triRegion.size()) ? rs.triRegion[k] : -1;
            int iid = (k < rs.triIsland.size()) ? rs.triIsland[k] : -1;
            bool exp = rid >= 0 && (size_t)rid < exploded.size() && exploded[(size_t)rid];
            if (exp) return 2;
            if (rid >= 0 && (size_t)rid < shipped.size() && shipped[(size_t)rid]) return 0;
            if (iid >= 0) return 1;
            return 3;
        };

        int nA = 0, nI = 0, nE = 0, nU = 0, nFillSkip = 0, nCvWouldAdd = 0, nI1Hole = 0;
        std::vector<int> uByRid, eByRid;
        const int ridCap = 4096;
        uByRid.assign(ridCap, 0);
        eByRid.assign(ridCap, 0);
        int uNoRid = 0;
        for (size_t k = 0; k < nTri; k++) {
            int rid = (k < rs.triRegion.size()) ? rs.triRegion[k] : -1;
            int iid = (k < rs.triIsland.size()) ? rs.triIsland[k] : -1;
            bool exp = rid >= 0 && (size_t)rid < exploded.size() && exploded[(size_t)rid];
            int cls = triClass(k);
            if (cls == 0) nA++;
            else if (cls == 1) nI++;
            else if (cls == 2) {
                nE++;
                if (rid >= 0 && rid < ridCap) eByRid[(size_t)rid]++;
            } else {
                nU++;
                if (rid >= 0 && rid < ridCap) uByRid[(size_t)rid]++;
                else uNoRid++;
            }
            if (iid < 0 && !exp) nFillSkip++;
            const bool cvFacet = !(rid >= 0 && (size_t)rid < shipped.size() && shipped[(size_t)rid]);
            if (cvFacet && iid < 0 && !exp) nCvWouldAdd++;
            if (rid < 0 && iid < 0) nI1Hole++;
        }
        int nRejTris = 0, nRejIsland = 0, nRejRegion = 0;
        for (const Region& rr : rs.rejected) {
            for (int t : rr.tris) {
                nRejTris++;
                if (t >= 0 && (size_t)t < rs.triIsland.size() && rs.triIsland[(size_t)t] >= 0)
                    nRejIsland++;
                if (t >= 0 && (size_t)t < rs.triRegion.size() && rs.triRegion[(size_t)t] >= 0)
                    nRejRegion++;
            }
        }
        std::fprintf(stderr,
                     "DIAG_COVER_SUM pass=%d analytic=%d island=%d exploded=%d uncovered=%d "
                     "total=%zu fillSkip=%d cvWouldAdd=%d i1Holes=%d rejectedList=%d "
                     "rejectedAlsoIsland=%d rejectedAlsoRegion=%d builtRidKeys=%d meshEGen=%d\n",
                     j6Pass, nA, nI, nE, nU, nTri, nFillSkip, nCvWouldAdd, nI1Hole, nRejTris,
                     nRejIsland, nRejRegion, nShippedRids, gDiagMeshEGen);
        if (nA) std::fprintf(stderr, "DIAG_COVER class=analytic n=%d\n", nA);
        if (nI) std::fprintf(stderr, "DIAG_COVER class=island n=%d\n", nI);
        for (int rid = 0; rid < ridCap; rid++) {
            if (eByRid[(size_t)rid] > 0)
                std::fprintf(stderr, "DIAG_COVER class=exploded rid=%d n=%d\n", rid,
                             eByRid[(size_t)rid]);
        }
        if (nE) std::fprintf(stderr, "DIAG_COVER class=exploded n=%d\n", nE);
        for (int rid = 0; rid < ridCap; rid++) {
            if (uByRid[(size_t)rid] > 0)
                std::fprintf(stderr, "DIAG_COVER class=uncovered rid=%d n=%d\n", rid,
                             uByRid[(size_t)rid]);
        }
        if (uNoRid)
            std::fprintf(stderr, "DIAG_COVER class=uncovered rid=-1 n=%d\n", uNoRid);
        if (nU) std::fprintf(stderr, "DIAG_COVER class=uncovered n=%d\n", nU);

        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            int idx = faceMap.FindIndex(built[i]);
            if (idx > 0 && i < builtRid.size()) faceRid[(size_t)idx] = builtRid[i];
        }

        std::vector<int> edgeT0(mv.nEdge, -1), edgeT1(mv.nEdge, -1);
        if (mv.triEdges) {
            for (size_t k = 0; k < nTri; k++) {
                for (int s = 0; s < 3; s++) {
                    int eid = mv.triEdges[k][s];
                    if (eid < 0 || (size_t)eid >= mv.nEdge) continue;
                    if (edgeT0[(size_t)eid] < 0) edgeT0[(size_t)eid] = (int)k;
                    else edgeT1[(size_t)eid] = (int)k;
                }
            }
        }

        auto matchMeshEid = [&](const TopoDS_Edge& e, const gp_Pnt& pa, const gp_Pnt& pb) -> int {
            for (size_t eid = 0; eid < meshE.size(); eid++) {
                if (!meshE[eid].IsNull() && e.IsSame(meshE[eid])) return (int)eid;
            }
            int best = -1;
            double bestD = 1e300;
            const double tol = 0.05;
            if (!mv.compEdges) return -1;
            for (size_t eid = 0; eid < mv.nEdge; eid++) {
                const auto& ev = mv.compEdges[eid];
                gp_Pnt p0 = pntOf(mv, ev.first), p1 = pntOf(mv, ev.second);
                const double d =
                    std::min(pa.Distance(p0) + pb.Distance(p1), pa.Distance(p1) + pb.Distance(p0));
                if (d < bestD) {
                    bestD = d;
                    best = (int)eid;
                }
            }
            return (best >= 0 && bestD <= 2.0 * tol) ? best : -1;
        };

        auto matchGeomCi = [&](const TopoDS_Edge& e) -> int {
            for (size_t ci = 0; ci < geom.size(); ci++) {
                for (const TopoDS_Edge& ge : geom[ci].edges) {
                    if (!ge.IsNull() && e.IsSame(ge)) return (int)ci;
                }
            }
            return -1;
        };

        int nUnc = 0, nTwin = 0, nNoCol = 0, nDegen = 0, nFree = 0;
        const double twinTol = 0.05;
        // %.4f printed 0.0000 in contain-r1; Precision::Confusion is 1e-7 and
        // would miss slightly-longer collapsed chords that still round to 0.
        const double degenTol = 5e-5;
        for (int i = 1; i <= anc.Extent(); i++) {
            if (anc(i).Extent() >= 2) continue;
            const TopoDS_Edge& e = TopoDS::Edge(anc.FindKey(i));
            TopoDS_Vertex va, vb;
            TopExp::Vertices(e, va, vb, Standard_True);
            gp_Pnt pa = BRep_Tool::Pnt(va), pb = BRep_Tool::Pnt(vb);
            const double len = pa.Distance(pb);
            int frid = -1;
            if (anc(i).Extent() == 1) {
                const TopoDS_Face& f = TopoDS::Face(anc(i).First());
                int fi = faceMap.FindIndex(f);
                if (fi > 0) frid = faceRid[(size_t)fi];
            }
            int meshEid = matchMeshEid(e, pa, pb);
            int otherTri = -1, otherCls = -1;
            if (meshEid >= 0 && (size_t)meshEid < edgeT0.size()) {
                int t0 = edgeT0[(size_t)meshEid], t1 = edgeT1[(size_t)meshEid];
                auto hostish = [&](int t) {
                    if (t < 0) return false;
                    int rid = ((size_t)t < rs.triRegion.size()) ? rs.triRegion[(size_t)t] : -1;
                    if (frid >= 0 && rid == frid) return true;
                    int cls = triClass((size_t)t);
                    if (frid < 0 && (cls == 1 || cls == 2)) return true;
                    return false;
                };
                if (hostish(t0) && !hostish(t1)) otherTri = t1;
                else if (hostish(t1) && !hostish(t0)) otherTri = t0;
                else if (t0 >= 0 && t1 >= 0) {
                    int c0 = triClass((size_t)t0), c1 = triClass((size_t)t1);
                    if (c0 == 3 && c1 != 3) otherTri = t0;
                    else if (c1 == 3 && c0 != 3) otherTri = t1;
                    else otherTri = t1;
                } else {
                    otherTri = (t0 >= 0 && !hostish(t0)) ? t0 : t1;
                }
                if (otherTri >= 0) otherCls = triClass((size_t)otherTri);
            }
            int ciMatch = -1;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                if (ci >= collapsed.size() || !collapsed[ci]) continue;
                const BoundaryChain& ch = rs.chains[ci];
                if (ch.meshVerts.size() < 2) continue;
                gp_Pnt p0 = pntOf(mv, ch.meshVerts.front());
                gp_Pnt p1 = pntOf(mv, ch.meshVerts.back());
                if ((pa.Distance(p0) < 0.05 && pb.Distance(p1) < 0.05) ||
                    (pa.Distance(p1) < 0.05 && pb.Distance(p0) < 0.05))
                    ciMatch = (int)ci;
            }
            const int geomCi = matchGeomCi(e);
            bool twin = false;
            if (len > twinTol) {
                for (int j = 1; j <= anc.Extent(); j++) {
                    if (j == i) continue;
                    const TopoDS_Edge& e2 = TopoDS::Edge(anc.FindKey(j));
                    if (e.IsSame(e2)) continue;
                    TopoDS_Vertex va2, vb2;
                    TopExp::Vertices(e2, va2, vb2, Standard_True);
                    gp_Pnt pa2 = BRep_Tool::Pnt(va2), pb2 = BRep_Tool::Pnt(vb2);
                    const double d2 = pa2.Distance(pb2);
                    if (d2 <= twinTol) continue;
                    const double d = std::min(pa.Distance(pa2) + pb.Distance(pb2),
                                              pa.Distance(pb2) + pb.Distance(pa2));
                    if (d <= 2.0 * twinTol) {
                        twin = true;
                        break;
                    }
                }
            }
            const int flagged = BRep_Tool::Degenerated(e) ? 1 : 0;
            const char* miss = "no-collapsed-partner";
            if (len <= degenTol || flagged) {
                miss = "degenerate-zero-length";
                nDegen++;
            } else if (otherCls == 3) {
                miss = "uncovered-tri";
                nUnc++;
            } else if (twin || ciMatch >= 0) {
                miss = "tshape-mismatch";
                nTwin++;
            } else {
                nNoCol++;
            }
            const void* tsp = diagTShapePtr(e);
            const void* meshTp = nullptr;
            int sameMeshE = 0;
            if (meshEid >= 0 && (size_t)meshEid < meshE.size() && !meshE[(size_t)meshEid].IsNull()) {
                meshTp = diagTShapePtr(meshE[(size_t)meshEid]);
                sameMeshE = e.IsSame(meshE[(size_t)meshEid]) ? 1 : 0;
            }
            const char* site = "?";
            int gen = -1;
            if (tsp) {
                auto sit = gDiagEdgeSite.find(tsp);
                if (sit != gDiagEdgeSite.end()) site = sit->second;
                auto git = gDiagEdgeGen.find(tsp);
                if (git != gDiagEdgeGen.end()) gen = git->second;
            }
            if (geomCi >= 0 && site[0] == '?') site = "appendCollapsed";
            else if (sameMeshE && site[0] == '?') site = "meshE-slot";
            const char* cty = "unk";
            double cFirst = 0, cLast = 0;
            try {
                BRepAdaptor_Curve ac(e);
                cFirst = ac.FirstParameter();
                cLast = ac.LastParameter();
                switch (ac.GetType()) {
                    case GeomAbs_Line:
                        cty = "line";
                        break;
                    case GeomAbs_Circle:
                        cty = "circle";
                        break;
                    case GeomAbs_Ellipse:
                        cty = "ellipse";
                        break;
                    case GeomAbs_BSplineCurve:
                        cty = "bspline";
                        break;
                    case GeomAbs_BezierCurve:
                        cty = "bezier";
                        break;
                    default:
                        cty = "other";
                        break;
                }
            } catch (const Standard_Failure&) {
            }
            if (std::strcmp(miss, "no-collapsed-partner") == 0)
                std::fprintf(stderr,
                             "DIAG_B02 #%d hostRid=%d edgeSite=%s nbrRid=%d otherCls=%d "
                             "copyKind=%s meshEid=%d sameMeshE=%d geomCi=%d\n",
                             nFree, frid, site, (otherTri >= 0 && (size_t)otherTri < rs.triRegion.size())
                                                      ? rs.triRegion[(size_t)otherTri]
                                                      : -1,
                             otherCls,
                             (sameMeshE ? "live-meshE"
                                        : (geomCi >= 0 ? "independent-geom" : "meshE-copy")),
                             meshEid, sameMeshE, geomCi);
            std::fprintf(stderr,
                         "DIAG_FREEEDGE #%d faceRid=%d miss=%s meshEid=%d otherTri=%d "
                         "otherCls=%d ci=%d geomCi=%d twin=%d pa=(%.4f,%.4f,%.4f) "
                         "pb=(%.4f,%.4f,%.4f) len=%.4f degen=%d sameV=%d tshape=%p "
                         "meshETshape=%p sameMeshE=%d gen=%d site=%s cty=%s "
                         "range=%.6f..%.6f\n",
                         nFree, frid, miss, meshEid, otherTri, otherCls, ciMatch, geomCi,
                         twin ? 1 : 0, pa.X(), pa.Y(), pa.Z(), pb.X(), pb.Y(), pb.Z(), len,
                         flagged, va.IsSame(vb) ? 1 : 0, tsp, meshTp, sameMeshE, gen, site, cty,
                         cFirst, cLast);
            nFree++;
        }
        std::fprintf(stderr,
                     "DIAG_FREEEDGE_SUM pass=%d freeEdges=%d uncovered-tri=%d "
                     "tshape-mismatch=%d no-collapsed-partner=%d degenerate-zero-length=%d "
                     "sum=%d\n",
                     j6Pass, nFree, nUnc, nTwin, nNoCol, nDegen,
                     nUnc + nTwin + nNoCol + nDegen);
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_COVER exception\n");
    }
}

double distAxesCyl(const gp_Ax3& a, const gp_Ax3& b) {
    const gp_XYZ u = a.Direction().XYZ();
    const gp_XYZ v = b.Direction().XYZ();
    const gp_XYZ w = b.Location().XYZ() - a.Location().XYZ();
    const gp_XYZ cr = u.Crossed(v);
    const double crn = cr.Modulus();
    if (crn <= std::sin(8.0 * kPi / 180.0) + 1e-15) return w.Crossed(u).Modulus();
    return std::fabs(w.Dot(cr)) / crn;
}

bool axesNearParallel8(const gp_Dir& a, const gp_Dir& b) {
    return a.XYZ().Crossed(b.XYZ()).Modulus() <= std::sin(8.0 * kPi / 180.0) + 1e-15;
}

bool cylCylTangentContact(const Region& a, const Region& b, double epsPlane) {
    if (!axesNearParallel8(a.ax.Direction(), b.ax.Direction())) return false;
    const double tol = std::max(epsPlane, std::max(a.maxVertexDev, b.maxVertexDev));
    const double d = distAxesCyl(a.ax, b.ax);
    const double ext = a.radius + b.radius;
    const double inn = std::fabs(a.radius - b.radius);
    return std::fabs(d - ext) <= tol || std::fabs(d - inn) <= tol;
}

bool cylCylParallelOffset(const Region& a, const Region& b, double epsPlane) {
    if (!axesNearParallel8(a.ax.Direction(), b.ax.Direction())) return false;
    const double tol = std::max(epsPlane, std::max(a.maxVertexDev, b.maxVertexDev));
    return distAxesCyl(a.ax, b.ax) > tol;
}

bool planeCylSideContact(const Region& plnR, const Region& cylR, double epsPlane) {
    const double adn = std::fabs(cylR.ax.Direction().Dot(plnR.ax.Direction()));
    if (adn > std::sin(3.0 * kPi / 180.0) + 1e-15) return false;
    const double d = std::fabs(plnR.ax.Direction().XYZ().Dot(
        cylR.ax.Location().XYZ() - plnR.ax.Location().XYZ()));
    const double tol =
        std::max(epsPlane, std::max(plnR.maxVertexDev, cylR.maxVertexDev));
    return std::fabs(d - cylR.radius) <= tol;
}

AnalyticCurve constructedPlaneCylCap(const Region& plnR, const Region& cylR) {
    AnalyticCurve out;
    if (!planePerpCylinder(plnR, cylR)) return out;
    out.kind = AnalyticCurve::Circ;
    out.circ = cylinderIsoCircle(cylR, planeVOnCylinder(plnR, cylR));
    return out;
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

AnalyticCurve meshAnchoredCylGenerator(const Region& cyl, const MeshView& mv,
                                       const BoundaryChain& ch) {
    AnalyticCurve out;
    if (ch.meshVerts.empty()) return out;
    gp_Dir axis = cyl.ax.Direction();
    gp_Pnt meshP = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
    gp_Vec toP(cyl.ax.Location(), meshP);
    toP -= gp_Vec(axis) * toP.Dot(axis);
    if (toP.Magnitude() < Precision::Confusion()) return out;
    gp_Dir rad(toP);
    gp_Pnt origin = cyl.ax.Location().Translated(gp_Vec(rad) * cyl.radius);
    out.kind = AnalyticCurve::Lin;
    out.lin = gp_Lin(origin, axis);
    return out;
}

AnalyticCurve bestCylCylConstructed(const Region& a, const Region& b, const MeshView& mv,
                                    const BoundaryChain& ch, double acceptR) {
    AnalyticCurve best;
    double bestR = 1e300;
    auto consider = [&](const AnalyticCurve& c) {
        if (c.kind == AnalyticCurve::None) return;
        const double r = chainResidual(c, mv, ch);
        if (r <= acceptR && r < bestR) {
            bestR = r;
            best = c;
        }
    };
    consider(constructedCylCylGenerator(a, b));
    consider(constructedCylCylGenerator(b, a));
    if (axesNearParallel8(a.ax.Direction(), b.ax.Direction())) {
        consider(meshAnchoredCylGenerator(a, mv, ch));
        consider(meshAnchoredCylGenerator(b, mv, ch));
    }
    return best;
}

const char* cylCylClassName(const Region& a, const Region& b, double epsPlane) {
    if (!axesNearParallel8(a.ax.Direction(), b.ax.Direction())) return "skew";
    const double d = distAxesCyl(a.ax, b.ax);
    if (d <= Precision::Confusion()) return "coax";
    if (cylCylTangentContact(a, b, epsPlane)) return "nearTan";
    if (cylCylParallelOffset(a, b, epsPlane)) return "parOff";
    return "other";
}

// A constructed generator is only shareable if it is a generator of BOTH
// cylinders. Off-axis (8° "parallel") lines sit on one surface and orphan
// on the other — Body9's MakeEdge-succeeds / J6-opens class.
bool constructedLinOnBothCylinders(const AnalyticCurve& c, const Region& a, const Region& b,
                                   double tol) {
    if (c.kind != AnalyticCurve::Lin) return false;
    if (!(tol > 0.0) || !std::isfinite(tol)) return false;
    if (!axesNearParallel8(c.lin.Direction(), a.ax.Direction())) return false;
    if (!axesNearParallel8(c.lin.Direction(), b.ax.Direction())) return false;
    gp_Ax3 linAx(c.lin.Location(), c.lin.Direction());
    const double dA = distAxesCyl(linAx, a.ax);
    const double dB = distAxesCyl(linAx, b.ax);
    return std::fabs(dA - a.radius) <= tol && std::fabs(dB - b.radius) <= tol;
}

double cylCylOnSurfaceTol(const Region& a, const Region& b, const MeshView& mv, double sewTol) {
    const double eps = derivedEpsPlane(mv);
    const double fit = std::max(a.maxVertexDev, b.maxVertexDev);
    const double sew = (sewTol > 0.0) ? sewTol : Precision::Confusion();
    return std::max(std::max(eps, fit) * 2.0, sew);
}

// Grow stamps filletNbrA = -2 on in-band arch-chain commits (refit_grow.cpp).
// Those boundaries are new edge classes the r2 fallback-only guard never saw.
constexpr int kArchChainSewMark = -2;

bool regionIsArchChainCommit(const Region* r) {
    // filletNbrA is remapped to -1 in refit_chains; maxVertexSnap < 0 is the
    // grow handshake that survives topology (in-band closed arch-chain only).
    return r && r->origin == Origin::CylGrow
           && (r->filletNbrA == kArchChainSewMark || r->maxVertexSnap < 0.0);
}

bool chainTouchesArchCommit(const Region* A, const Region* B) {
    return regionIsArchChainCommit(A) || regionIsArchChainCommit(B);
}

bool regionHasChain(const Region& r, int ci) {
    for (const Loop& lp : r.loops) {
        for (int idx : lp.chainIdx)
            if (idx == ci) return true;
    }
    return false;
}

bool bothRegionsReferenceChain(const RegionSet& rs, const BoundaryChain& ch, int ci) {
    const Region* A = regionById(rs, ch.regA);
    const Region* B = regionById(rs, ch.regB);
    if (!A || !B) return false;
    return regionHasChain(*A, ci) && regionHasChain(*B, ci);
}

bool faceHasEdgeTShape(const TopoDS_Face& f, const TopoDS_Edge& e) {
    if (f.IsNull() || e.IsNull()) return false;
    try {
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            if (ex.Current().IsSame(e)) return true;
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

int matchCollapsedChainToSegment(const MeshView& mv, const RegionSet& rs,
                                 const std::vector<char>& collapsed, const gp_Pnt& pa,
                                 const gp_Pnt& pb, double tol) {
    int best = -1;
    double bestD = 1e300;
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        if (ci >= collapsed.size() || !collapsed[ci]) continue;
        const BoundaryChain& ch = rs.chains[ci];
        if (ch.meshVerts.size() < 2) continue;
        gp_Pnt p0 = pntOf(mv, ch.meshVerts.front());
        gp_Pnt p1 = pntOf(mv, ch.meshVerts.back());
        const double d00 = pa.Distance(p0), d01 = pa.Distance(p1);
        const double d10 = pb.Distance(p0), d11 = pb.Distance(p1);
        const double d = std::min(std::min(d00 + d11, d01 + d10), std::min(d00 + d10, d01 + d11));
        if (d < bestD) {
            bestD = d;
            best = (int)ci;
        }
    }
    return (best >= 0 && bestD <= 2.0 * tol) ? best : -1;
}

void emitDiagFallback(int ci, const BoundaryChain& ch, const char* cls, const char* outcome,
                      int nA, int nB, double residual, size_t nV) {
    if (!diagJ6Enabled()) return;
    std::fprintf(stderr,
                 "DIAG_FALLBACK ci=%d class=%s outcome=%s nFaceA=%d nFaceB=%d "
                 "regA=%d regB=%d residual=%.4f nV=%zu\n",
                 ci, cls, outcome, nA, nB, ch.regA, ch.regB, residual, nV);
}

AnalyticCurve intersectSurfaces(const Region& A, const Region& B, const MeshView& mv,
                                const BoundaryChain& ch, double sewTol, WarnFn warn,
                                bool* usedCylCylFallback = nullptr) {
    if (usedCylCylFallback) *usedCylCylFallback = false;
    const double tolAng = Precision::Angular();
    const double tol = std::max(sewTol, Precision::Confusion());
    try {
        if (A.type == SurfType::Plane && B.type == SurfType::Plane) {
            IntAna_QuadQuadGeo iq(asPlane(A), asPlane(B), tolAng, tol);
            const double acceptR = intAnaAcceptResidual(mv, ch, sewTol, &A, &B);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol, acceptR);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: IntAna plane|plane empty/same — keeping mesh polyline");
            return c;
        }
        if (A.type == SurfType::Plane && B.type == SurfType::Cylinder) {
            double H = std::fabs(B.vMax - B.vMin);
            IntAna_QuadQuadGeo iq(asPlane(A), cylForIntersect(B), tolAng, tol, H);
            const double acceptR = intAnaAcceptResidual(mv, ch, sewTol, &A, &B);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol, acceptR);
            // Oblique plane|cyl: IntAna ellipse on 2-vertex coarse chains often
            // exceeds the first-pass residual gate but still matches the mesh band.
            if (c.kind == AnalyticCurve::None && iq.IsDone() &&
                iq.TypeInter() == IntAna_Ellipse && ch.meshVerts.size() <= 3) {
                const double fit = std::max(A.maxVertexDev, B.maxVertexDev);
                const double loose = std::max(acceptR, fit * 8.0 + derivedEpsPlane(mv) * 3.0);
                c = pickIntAna(iq, mv, ch, sewTol, loose);
            }
            const double epsPl = derivedEpsPlane(mv);
            // G4: side-grazing plane|cyl => generator line (§2.5). Re-evaluate
            // geometry here — ch.tangent can miss when |dist-R| is within fit
            // residual but above the epsPlane floor.
            if (c.kind == AnalyticCurve::None && planeCylSideContact(A, B, epsPl))
                c = constructedGenerator(B, A);
            // Cap circle from the fitted cylinder when IntAna misses or picks
            // the wrong branch (coarse meshes: pickIntAna residual gate).
            if (c.kind == AnalyticCurve::None && planePerpCylinder(A, B))
                c = constructedPlaneCylCap(A, B);
            // Ellipse branch on a cap plane is often a noisy IntAna circle; prefer
            // the fitted iso-circle when the plane is perpendicular to the axis.
            if (c.kind == AnalyticCurve::Elips && planePerpCylinder(A, B) && !B.closed360) {
                AnalyticCurve cap = constructedPlaneCylCap(A, B);
                if (cap.kind == AnalyticCurve::Circ) c = cap;
            }
            if (c.kind == AnalyticCurve::None) {
                if (diagP2Enabled()) {
                    IntAna_ResultType ty = iq.IsDone() ? iq.TypeInter() : IntAna_Empty;
                    std::fprintf(stderr,
                                 "DIAG_P2 plane|cyl regA=%d regB=%d R=%.4f ty=%d "
                                 "tangent=%d side=%d perp=%d nV=%zu\n",
                                 ch.regA, ch.regB, B.radius, (int)ty, (int)ch.tangent,
                                 (int)planeCylSideContact(A, B, epsPl),
                                 (int)planePerpCylinder(A, B), ch.meshVerts.size());
                }
                emit(warn, "smooth: IntAna plane|cyl empty/same — keeping mesh polyline");
            }
            return c;
        }
        if (A.type == SurfType::Cylinder && B.type == SurfType::Plane) {
            return intersectSurfaces(B, A, mv, ch, sewTol, warn);
        }
        if (A.type == SurfType::Cylinder && B.type == SurfType::Cylinder) {
            IntAna_QuadQuadGeo iq(cylForIntersect(A), cylForIntersect(B), tol);
            const double acceptR = intAnaAcceptResidual(mv, ch, sewTol, &A, &B);
            AnalyticCurve c = pickIntAna(iq, mv, ch, sewTol, acceptR);
            const double epsPl = derivedEpsPlane(mv);
            bool constructed = false;
            // Near-tangent fillet|bore and arch|arch: IntAna empty on coarse
            // meshes — construct the shared generator from fitted geometry.
            if (c.kind == AnalyticCurve::None && cylCylTangentContact(A, B, epsPl)) {
                c = constructedCylCylGenerator(A, B);
                constructed = (c.kind != AnalyticCurve::None);
            }
            if (c.kind == AnalyticCurve::None && cylCylParallelOffset(A, B, epsPl)) {
                c = bestCylCylConstructed(A, B, mv, ch, acceptR);
                constructed = (c.kind != AnalyticCurve::None);
            }
            // Log off-surface constructions; do not discard here. A line can
            // still be dual-referenced (Body18). Orphans are rejected after
            // both faces exist (dual-face / J6 free-edge ledger).
            if (constructed) {
                if (usedCylCylFallback) *usedCylCylFallback = true;
                const double onTol = cylCylOnSurfaceTol(A, B, mv, sewTol);
                if (diagJ6Enabled() && !constructedLinOnBothCylinders(c, A, B, onTol))
                    std::fprintf(stderr,
                                 "DIAG_FALLBACK class=%s outcome=off-surface "
                                 "regA=%d regB=%d onTol=%.4f nV=%zu\n",
                                 cylCylClassName(A, B, epsPl), ch.regA, ch.regB, onTol,
                                 ch.meshVerts.size());
            }
            if (c.kind == AnalyticCurve::None) {
                if (diagP2Enabled()) {
                    IntAna_ResultType ty = iq.IsDone() ? iq.TypeInter() : IntAna_Empty;
                    std::fprintf(stderr,
                                 "DIAG_P2 cyl|cyl regA=%d regB=%d R=%.4f/%.4f ty=%d "
                                 "tangent=%d tanContact=%d parOff=%d nV=%zu\n",
                                 ch.regA, ch.regB, A.radius, B.radius, (int)ty, (int)ch.tangent,
                                 (int)cylCylTangentContact(A, B, epsPl),
                                 (int)cylCylParallelOffset(A, B, epsPl), ch.meshVerts.size());
                }
                emit(warn, "smooth: IntAna cyl|cyl empty/same — keeping mesh polyline");
            }
            return c;
        }
    } catch (const Standard_Failure&) {
        emit(warn, "smooth: IntAna threw — keeping mesh polyline");
    }
    return {};
}

TopoDS_Edge makeFullCircle(const gp_Circ& circ, const TopoDS_Vertex& V);
TopoDS_Edge makeFullEllipse(const gp_Elips& el, const TopoDS_Vertex& V);

// SPEC-F2: MakeEdge from projected parameters on the IntAna curve (chain is
// selector). Bypasses BRepBuilderAPI_MakeEdge's PointProjectionFailed when
// the vertex sits outside its current tolerance but on an accepted curve.
TopoDS_Edge bindEdgeByParam(const Handle(Geom_Curve)& gc, const TopoDS_Vertex& v1,
                            const TopoDS_Vertex& v2, double p1, double p2) {
    if (gc.IsNull() || v1.IsNull() || v2.IsNull()) return {};
    try {
        BRepBuilderAPI_MakeEdge me(gc, v1, v2, p1, p2);
        if (me.IsDone()) return me.Edge();
    } catch (const Standard_Failure&) {
    }
    return {};
}

TopoDS_Edge makeEdgeFromCurve(const AnalyticCurve& c, const TopoDS_Vertex& v1,
                              const TopoDS_Vertex& v2, bool closedFull) {
    try {
        if (c.kind == AnalyticCurve::Lin) {
            BRepBuilderAPI_MakeEdge me(c.lin, v1, v2);
            if (me.IsDone()) return me.Edge();
            const double p1 = ElCLib::Parameter(c.lin, BRep_Tool::Pnt(v1));
            const double p2 = ElCLib::Parameter(c.lin, BRep_Tool::Pnt(v2));
            if (std::fabs(p1 - p2) <= Precision::Confusion()) return {};
            Handle(Geom_Line) gl = new Geom_Line(c.lin);
            TopoDS_Edge pe = bindEdgeByParam(gl, v1, v2, p1, p2);
            if (!pe.IsNull()) return pe;
        } else if (c.kind == AnalyticCurve::Circ) {
            if (closedFull || v1.IsSame(v2)) {
                BRepBuilderAPI_MakeEdge me(c.circ, v1, v1);
                if (me.IsDone()) return me.Edge();
                BRepBuilderAPI_MakeEdge me2(c.circ);
                if (me2.IsDone()) return me2.Edge();
            } else {
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
            // Complementary arc. Keep the 3d range non-negative: a negative
            // Range on Geom_Circle (r9 180° seam-straddle) is
            // BRepCheck_BadOrientationOfSubshape and ShapeFix cannot recover.
            const double db = twopi - df;
            const double t0 = wrap(p1 - db);  // parameter of vB, in [0, 2π)
            B.Add(e, vB.Oriented(TopAbs_FORWARD));
            B.Add(e, vA.Oriented(TopAbs_REVERSED));
            B.Range(e, t0, t0 + db);
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
    try {
        GC_MakeArcOfCircle mk(BRep_Tool::Pnt(vA), midHint, BRep_Tool::Pnt(vB));
        if (mk.IsDone()) {
            BRepBuilderAPI_MakeEdge me(mk.Value(), vA, vB);
            if (me.IsDone()) return me.Edge();
        }
    } catch (const Standard_Failure&) {
    }
    return TopoDS_Edge();
}

// ME_CYLPLN_ELLIPSE_PROJ: same mid-vertex selector as makeArc, on the IntAna
// ellipse. Open arcs keep the chain as the selector.
TopoDS_Edge makeEllipseArc(const gp_Elips& el, const TopoDS_Vertex& vA, const TopoDS_Vertex& vB,
                           const gp_Pnt& midHint) {
    const double twopi = 2.0 * kPi;
    auto wrap = [&](double t) {
        while (t < 0.0) t += twopi;
        while (t >= twopi) t -= twopi;
        return t;
    };
    try {
        double p1 = wrap(ElCLib::Parameter(el, BRep_Tool::Pnt(vA)));
        double p2 = wrap(ElCLib::Parameter(el, BRep_Tool::Pnt(vB)));
        double pm = wrap(ElCLib::Parameter(el, midHint));
        double df = p2 - p1;
        while (df <= 0.0) df += twopi;
        double dm = pm - p1;
        while (dm < 0.0) dm += twopi;
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(el);
        BRep_Builder B;
        TopoDS_Edge e;
        B.MakeEdge(e, ge, Precision::Confusion());
        if (dm <= df + 1e-9) {
            B.Add(e, vA.Oriented(TopAbs_FORWARD));
            B.Add(e, vB.Oriented(TopAbs_REVERSED));
            B.Range(e, p1, p1 + df);
            B.UpdateVertex(vA, p1, e, std::max(BRep_Tool::Tolerance(vA), Precision::Confusion()));
            B.UpdateVertex(vB, p1 + df, e, std::max(BRep_Tool::Tolerance(vB), Precision::Confusion()));
        } else {
            const double db = twopi - df;
            B.Add(e, vA.Oriented(TopAbs_FORWARD));
            B.Add(e, vB.Oriented(TopAbs_REVERSED));
            B.Range(e, p1, p1 + db);
            B.UpdateVertex(vA, p1, e, std::max(BRep_Tool::Tolerance(vA), Precision::Confusion()));
            B.UpdateVertex(vB, p1 + db, e, std::max(BRep_Tool::Tolerance(vB), Precision::Confusion()));
        }
        return e;
    } catch (const Standard_Failure&) {
    }
    try {
        BRepBuilderAPI_MakeEdge me(el, vA, vB);
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
        diagNoteEdge(meshE[(size_t)id], "makeFacet");
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
        diagNoteEdge(meshE[(size_t)eid], "appendPolyline");
        B.Add(w, fwd ? meshE[(size_t)eid] : TopoDS::Edge(meshE[(size_t)eid].Reversed()));
    }
    return true;
}

bool appendCollapsed(BRep_Builder& B, TopoDS_Wire& w, const ChainGeom& g, bool reversed) {
    if (g.edges.empty()) return false;
    if (!reversed) {
        for (const auto& e : g.edges) {
            diagNoteEdge(e, "appendCollapsed");
            B.Add(w, e);
        }
    } else {
        for (int i = (int)g.edges.size() - 1; i >= 0; i--) {
            TopoDS_Edge rev = TopoDS::Edge(g.edges[(size_t)i].Reversed());
            diagNoteEdge(rev, "appendCollapsed");
            B.Add(w, rev);
        }
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

bool collectWireEdges(const TopoDS_Wire& w, std::vector<TopoDS_Edge>& edges) {
    edges.clear();
    if (w.IsNull()) return false;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next())
            edges.push_back(TopoDS::Edge(ex.Current()));
    } catch (const Standard_Failure&) {
        edges.clear();
    }
    if (edges.empty()) {
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            edges.push_back(TopoDS::Edge(it.Value()));
        }
    }
    return !edges.empty();
}

bool rotateEdgesToVertex(std::vector<TopoDS_Edge>& edges, const TopoDS_Vertex& V) {
    if (edges.empty() || V.IsNull()) return false;
    auto startOf = [&](const TopoDS_Edge& e) -> TopoDS_Vertex {
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(e, v1, v2, Standard_True);
        return v1;
    };
    for (size_t i = 0; i < edges.size(); i++) {
        if (startOf(edges[i]).IsSame(V)) {
            std::rotate(edges.begin(), edges.begin() + (std::ptrdiff_t)i, edges.end());
            return true;
        }
    }
    for (size_t i = 0; i < edges.size(); i++) {
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(edges[i], v1, v2, Standard_False);
        if (v1.IsSame(V) || v2.IsSame(V)) {
            if (!startOf(edges[i]).IsSame(V)) edges[i] = TopoDS::Edge(edges[i].Reversed());
            std::rotate(edges.begin(), edges.begin() + (std::ptrdiff_t)i, edges.end());
            return true;
        }
    }
    return false;
}

void setFaceOutward(TopoDS_Face& f, bool outwardNormal);

bool diagPlatesEnabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_DIAG_PLATES");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

bool collapseDiagEnabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_COLLAPSE_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// STL2STEP_FAIL_RID — parsed once (G0.1). rid 0 is legal; do not copy the
// DIAG `v[0] != '0'` idiom. Unset/empty ⇒ disabled. Malformed ⇒ disabled.
struct FailRidKnob {
    bool ready = false;
    bool enabled = false;
    bool malformed = false;
    bool warned = false;
    std::vector<int> ids;
};

FailRidKnob& failRidKnob() {
    static FailRidKnob k;
    if (k.ready) return k;
    k.ready = true;
    const char* v = std::getenv("STL2STEP_FAIL_RID");
    if (!v || !v[0]) return k;
    const char* p = v;
    std::vector<int> ids;
    while (*p) {
        while (*p == ' ' || *p == '\t') ++p;
        if (!*p) break;
        if (*p == ',') {
            k.malformed = true;
            break;
        }
        errno = 0;
        char* end = nullptr;
        const long x = std::strtol(p, &end, 10);
        if (end == p || errno == ERANGE || x < static_cast<long>(INT_MIN) ||
            x > static_cast<long>(INT_MAX)) {
            k.malformed = true;
            break;
        }
        ids.push_back(static_cast<int>(x));
        p = end;
        while (*p == ' ' || *p == '\t') ++p;
        if (*p == ',') {
            ++p;
            continue;
        }
        if (*p) {
            k.malformed = true;
            break;
        }
    }
    if (k.malformed || ids.empty()) {
        k.enabled = false;
        k.ids.clear();
        return k;
    }
    k.enabled = true;
    k.ids = std::move(ids);
    return k;
}

bool failRidHits(int id) {
    const FailRidKnob& k = failRidKnob();
    if (!k.enabled) return false;
    for (int x : k.ids)
        if (x == id) return true;
    return false;
}

void emitFailRidWarningOnce(WarnFn warn) {
    FailRidKnob& k = failRidKnob();
    if (!k.malformed || k.warned) return;
    k.warned = true;
    emit(warn, "STL2STEP_FAIL_RID: malformed value, ignored");
}

void diagCascadeInject() {
    if (!collapseDiagEnabled()) return;
    const FailRidKnob& k = failRidKnob();
    if (!k.enabled) return;
    std::fprintf(stderr, "DIAG_CASCADE inject=");
    for (size_t i = 0; i < k.ids.size(); i++) {
        if (i) std::fputc(',', stderr);
        std::fprintf(stderr, "%d", k.ids[i]);
    }
    std::fputc('\n', stderr);
}

const char* brepCheckName(int st) {
    switch (st) {
        case BRepCheck_NoError: return "NoError";
        case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
        case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
        case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
        case BRepCheck_No3DCurve: return "No3DCurve";
        case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
        case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
        case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
        case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
        case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
        case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
        case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
        case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
        case BRepCheck_FreeEdge: return "FreeEdge";
        case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
        case BRepCheck_InvalidRange: return "InvalidRange";
        case BRepCheck_EmptyWire: return "EmptyWire";
        case BRepCheck_RedundantEdge: return "RedundantEdge";
        case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
        case BRepCheck_NoSurface: return "NoSurface";
        case BRepCheck_InvalidWire: return "InvalidWire";
        case BRepCheck_RedundantWire: return "RedundantWire";
        case BRepCheck_IntersectingWires: return "IntersectingWires";
        case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
        case BRepCheck_EmptyShell: return "EmptyShell";
        case BRepCheck_RedundantFace: return "RedundantFace";
        case BRepCheck_InvalidImbricationOfShells: return "InvalidImbricationOfShells";
        case BRepCheck_UnorientableShape: return "UnorientableShape";
        case BRepCheck_NotClosed: return "NotClosed";
        case BRepCheck_NotConnected: return "NotConnected";
        case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
        case BRepCheck_BadOrientation: return "BadOrientation";
        case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
        case BRepCheck_InvalidPolygonOnTriangulation: return "InvalidPolygonOnTriangulation";
        case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
        case BRepCheck_EnclosedRegion: return "EnclosedRegion";
        case BRepCheck_CheckFail: return "CheckFail";
        default: return "?";
    }
}

void dumpBRepStatuses(const Handle(BRepCheck_Result)& res) {
    if (res.IsNull()) return;
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next()) {
        const int st = (int)it.Value();
        if (st == (int)BRepCheck_NoError) continue;
        std::fprintf(stderr, " st=%d(%s)", st, brepCheckName(st));
    }
}

double pcurveSignedArea(const TopoDS_Face& f, const TopoDS_Wire& w) {
    double A = 0.0;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
            const TopoDS_Edge e = ex.Current();
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
            if (pc.IsNull()) continue;
            const int N = 12;
            gp_Pnt2d p0 = pc->Value(a);
            for (int i = 1; i <= N; i++) {
                const double t = a + (b - a) * (double)i / (double)N;
                gp_Pnt2d p1 = pc->Value(t);
                A += p0.X() * p1.Y() - p1.X() * p0.Y();
                p0 = p1;
            }
        }
    } catch (const Standard_Failure&) {
    }
    return 0.5 * A;
}

void dumpFaceTopo(int rid, const TopoDS_Face& f, const char* tag, const char* phase) {
    if (!collapseDiagEnabled() || f.IsNull()) return;
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        int nW = 0, nE = 0, nPc = 0, nClosedW = 0, nWalkStall = 0;
        double uvA = 0.0;
        for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next()) {
            nW++;
            const TopoDS_Wire w = TopoDS::Wire(wx.Current());
            if (w.Closed() || BRep_Tool::IsClosed(w)) nClosedW++;
            int nWe = 0, nIt = 0;
            for (TopoDS_Iterator it(w); it.More(); it.Next())
                if (it.Value().ShapeType() == TopAbs_EDGE) nIt++;
            int nWalk = 0;
            try {
                for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) nWalk++;
            } catch (const Standard_Failure&) {
            }
            if (nWalk < nIt) nWalkStall++;
            nWe = nIt;
            nE += nWe;
            uvA += pcurveSignedArea(f, w);
            for (TopExp_Explorer ex(w, TopAbs_EDGE); ex.More(); ex.Next()) {
                const TopoDS_Edge e = TopoDS::Edge(ex.Current());
                Standard_Real a = 0, b = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
                if (!pc.IsNull()) nPc++;
            }
        }
        std::fprintf(stderr,
                     "DIAG_PARTIAL_TOPO rid=%d phase=%s tag=%s valid=%d nW=%d nE=%d closedW=%d "
                     "walkStall=%d pc=%d/%d uvA=%.4f faceOri=%s",
                     rid, phase ? phase : "-", tag ? tag : "-", an.IsValid() ? 1 : 0, nW, nE,
                     nClosedW, nWalkStall, nPc, nE, uvA,
                     f.Orientation() == TopAbs_FORWARD ? "F" : "R");
        dumpBRepStatuses(an.Result(f));
        std::fprintf(stderr, "\n");
        int wi = 0;
        for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), wi++) {
            const TopoDS_Wire w = TopoDS::Wire(wx.Current());
            Handle(BRepCheck_Result) wr = an.Result(w);
            std::fprintf(stderr, "DIAG_PARTIAL_WIRE rid=%d wi=%d closed=%d uvA=%.4f", rid, wi,
                         (w.Closed() || BRep_Tool::IsClosed(w)) ? 1 : 0, pcurveSignedArea(f, w));
            dumpBRepStatuses(wr);
            std::fprintf(stderr, "\n");
            int ei = 0;
            try {
                for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next(), ei++) {
                    const TopoDS_Edge e = ex.Current();
                    Handle(BRepCheck_Result) er = an.Result(e);
                    Standard_Real a = 0, b = 0;
                    Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
                    Standard_Real f3 = 0, l3 = 0;
                    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
                    const char* cty = "none";
                    if (!c3.IsNull()) {
                        if (c3->DynamicType() == STANDARD_TYPE(Geom_Line)) cty = "line";
                        else if (c3->DynamicType() == STANDARD_TYPE(Geom_Circle)) cty = "circ";
                        else if (c3->DynamicType() == STANDARD_TYPE(Geom_Ellipse)) cty = "elips";
                        else cty = "other";
                    }
                    gp_Pnt2d uv0, uv1;
                    if (!pc.IsNull()) {
                        uv0 = pc->Value(a);
                        uv1 = pc->Value(b);
                    }
                    TopoDS_Vertex va = ex.CurrentVertex();
                    gp_Pnt pa = va.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(va);
                    std::fprintf(stderr,
                                 "DIAG_PARTIAL_EDGE rid=%d wi=%d ei=%d hasPC=%d cty=%s "
                                 "pcUV=(%.4f,%.4f)->(%.4f,%.4f) c3Range=[%.4f,%.4f] eOri=%s "
                                 "v=(%.3f,%.3f,%.3f)",
                                 rid, wi, ei, pc.IsNull() ? 0 : 1, cty,
                                 pc.IsNull() ? 0.0 : uv0.X(), pc.IsNull() ? 0.0 : uv0.Y(),
                                 pc.IsNull() ? 0.0 : uv1.X(), pc.IsNull() ? 0.0 : uv1.Y(), f3, l3,
                                 e.Orientation() == TopAbs_FORWARD ? "F" : "R", pa.X(), pa.Y(),
                                 pa.Z());
                    dumpBRepStatuses(er);
                    std::fprintf(stderr, "\n");
                }
            } catch (const Standard_Failure&) {
            }
        }
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_PARTIAL_TOPO rid=%d phase=%s tag=%s threw\n", rid,
                     phase ? phase : "-", tag ? tag : "-");
    }
}

void dumpShellCheck(const TopoDS_Shell& sh, const std::vector<TopoDS_Face>& built,
                    const std::vector<int>& builtRid, const RegionSet& rs) {
    if (!collapseDiagEnabled() || sh.IsNull()) return;
    try {
        BRepCheck_Analyzer an(sh, Standard_True);
        std::fprintf(stderr, "DIAG_SHELL valid=%d closed=%d nF=%zu", an.IsValid() ? 1 : 0,
                     BRep_Tool::IsClosed(sh) ? 1 : 0, built.size());
        dumpBRepStatuses(an.Result(sh));
        std::fprintf(stderr, "\n");
        for (size_t i = 0; i < built.size(); i++) {
            const TopoDS_Face& f = built[i];
            if (f.IsNull()) continue;
            const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
            Handle(BRepCheck_Result) fr = an.Result(f);
            bool faceBad = false;
            BRepCheck_Analyzer fan(f, Standard_True);
            if (!fan.IsValid()) faceBad = true;
            if (!fr.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(fr->Status()); it.More(); it.Next())
                    if (it.Value() != BRepCheck_NoError) faceBad = true;
            }
            if (!faceBad) continue;
            const Region* rr = regionById(rs, rid);
            const char* ty = "facet";
            double R = 0.0;
            int nTri = 0;
            if (rr) {
                ty = (rr->type == SurfType::Cylinder) ? "cyl"
                     : (rr->type == SurfType::Plane)  ? "pln"
                                                      : "oth";
                R = rr->radius;
                nTri = (int)rr->tris.size();
            }
            std::fprintf(stderr,
                         "DIAG_SHELL_FACE rid=%d type=%s R=%.4f nTri=%d faceValid=%d ori=%s", rid,
                         ty, R, nTri, fan.IsValid() ? 1 : 0,
                         f.Orientation() == TopAbs_FORWARD ? "F" : "R");
            dumpBRepStatuses(fr);
            dumpBRepStatuses(fan.Result(f));
            std::fprintf(stderr, "\n");
        }
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            int idx = faceMap.FindIndex(built[i]);
            if (idx > 0 && i < builtRid.size()) faceRid[(size_t)idx] = builtRid[i];
        }
        int nBadE = 0;
        for (int i = 1; i <= anc.Extent(); i++) {
            const TopoDS_Edge& e = TopoDS::Edge(anc.FindKey(i));
            Handle(BRepCheck_Result) er = an.Result(e);
            bool bad = false;
            if (!er.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(er->Status()); it.More(); it.Next())
                    if (it.Value() != BRepCheck_NoError) bad = true;
            }
            if (!bad) continue;
            if (nBadE >= 40) {
                std::fprintf(stderr, "DIAG_SHELL_EDGE ... truncated\n");
                break;
            }
            int r0 = -1, r1 = -1;
            bool aFwd = false, bFwd = false;
            const int nAnc = anc(i).Extent();
            if (nAnc >= 1) {
                const TopoDS_Face& fa = TopoDS::Face(anc(i).First());
                int ia = faceMap.FindIndex(fa);
                if (ia > 0) r0 = faceRid[(size_t)ia];
                for (TopExp_Explorer ex(fa, TopAbs_EDGE); ex.More(); ex.Next()) {
                    if (ex.Current().IsSame(e)) {
                        aFwd = ex.Current().Orientation() == TopAbs_FORWARD;
                        break;
                    }
                }
            }
            if (nAnc >= 2) {
                const TopoDS_Face& fb = TopoDS::Face(anc(i).Last());
                int ib = faceMap.FindIndex(fb);
                if (ib > 0) r1 = faceRid[(size_t)ib];
                for (TopExp_Explorer ex(fb, TopAbs_EDGE); ex.More(); ex.Next()) {
                    if (ex.Current().IsSame(e)) {
                        bFwd = ex.Current().Orientation() == TopAbs_FORWARD;
                        break;
                    }
                }
            }
            std::fprintf(stderr,
                         "DIAG_SHELL_EDGE rids=%d,%d nAnc=%d sameOri=%d aFwd=%d bFwd=%d", r0, r1,
                         nAnc, (aFwd == bFwd) ? 1 : 0, aFwd ? 1 : 0, bFwd ? 1 : 0);
            dumpBRepStatuses(er);
            std::fprintf(stderr, "\n");
            nBadE++;
        }
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_SHELL threw\n");
    }
}

void diagPlateMakeFaceFail(int rid, const TopoDS_Face& f) {
    if (!diagPlatesEnabled()) return;
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        std::fprintf(stderr, "DIAG_PLATE rid=%d valid=%d", rid, an.IsValid() ? 1 : 0);
        if (!an.IsValid() && !f.IsNull()) {
            Handle(BRepCheck_Result) res = an.Result(f);
            if (!res.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More();
                     it.Next())
                    std::fprintf(stderr, " st=%d", (int)it.Value());
            }
        }
        std::fprintf(stderr, "\n");
    } catch (const Standard_Failure&) {
    }
}

void emitFaceSurfDiag(int ridHint, const TopoDS_Face& F);
Handle(Geom2d_Curve) makePCurveOnSurf(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                                      const Handle(Geom_Surface)& srf, const char* kind);
thread_local std::unordered_set<const void*> gSeamTShapes;
thread_local int gCompNTri = 0;

// BRep_Builder::MakeFace + Add keeps the input wire's TShapes (J2).
// BRepBuilderAPI_MakeFace(plane/surf, wire) may project/copy edges and is the
// mixed-analytic free-edge source on live multi-inner plates.
bool makeFaceKeep(const Handle(Geom_Surface)& surf, const TopoDS_Wire& outer,
                  const std::vector<TopoDS_Wire>& inners, bool outward, TopoDS_Face& outF) {
    if (surf.IsNull() || outer.IsNull()) return false;
    try {
        BRep_Builder B;
        TopoDS_Face f;
        B.MakeFace(f, surf, Precision::Confusion());
        B.Add(f, outer);
        for (const auto& iw : inners) {
            if (iw.IsNull()) return false;
            B.Add(f, iw);
        }
        setFaceOutward(f, outward);
        if (gCompNTri < 10000) {
            TopLoc_Location loc;
            Handle(Geom_Surface) sF = BRep_Tool::Surface(f, loc);
            if (!sF.IsNull()) {
                BRep_Builder Bb;
                for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
                    TopoDS_Edge eW = TopoDS::Edge(ex.Current());
                    const void* ts = diagTShapePtr(eW);
                    if (!ts || !gSeamTShapes.count(ts)) continue;
                    Standard_Boolean hasPc = Standard_False;
                    Standard_Real pf = 0, pl = 0;
                    (void)BRep_Tool::CurveOnSurface(eW, sF, loc, pf, pl, &hasPc);
                    if (hasPc) continue;
                    Standard_Real a = 0, b = 0;
                    Handle(Geom_Curve) c3 = BRep_Tool::Curve(eW, a, b);
                    if (c3.IsNull() || b - a <= Precision::PConfusion()) continue;
                    Handle(Geom2d_Curve) c2d = makePCurveOnSurf(c3, a, b, sF, nullptr);
                    if (c2d.IsNull()) continue;
                    double tol = Precision::Confusion();
                    const double existing = BRep_Tool::Tolerance(eW);
                    if (existing > tol) tol = existing;
                    Bb.UpdateEdge(eW, c2d, sF, loc, tol);
                    Bb.Range(eW, sF, loc, a, b);
                }
            }
        }
        outF = f;
        emitFaceSurfDiag(-1, outF);
        return !outF.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

Handle(Geom_CylindricalSurface) cylSurfaceForRegion(const Region& r) {
    gp_Ax3 ax = r.ax;
    if (!r.closed360) {
        double u0 = r.uMin;
        while (u0 > kPi) u0 -= 2.0 * kPi;
        while (u0 <= -kPi) u0 += 2.0 * kPi;
        if (std::fabs(u0) > 1e-14) {
            ax.Rotate(gp_Ax1(ax.Location(), ax.Direction()), u0);
        }
    }
    return new Geom_CylindricalSurface(gp_Cylinder(ax, r.radius));
}

bool singleChain360Cap(const RegionSet& rs, int ci, const Region*& cylOut, double& vOut) {
    cylOut = nullptr;
    if (ci < 0) return false;
    for (const Region& r : rs.regions) {
        if (!r.closed360 || r.type != SurfType::Cylinder) continue;
        for (const Loop& lp : r.loops) {
            if (lp.role != LoopRole::CapLow && lp.role != LoopRole::CapHigh) continue;
            if (lp.chainIdx.size() != 1 || lp.chainIdx[0] != ci) continue;
            cylOut = &r;
            vOut = (lp.role == LoopRole::CapLow) ? r.vMin : r.vMax;
            return true;
        }
    }
    return false;
}

bool shellIsValid(const TopoDS_Shape& sh) {
    if (sh.IsNull()) return false;
    try {
        BRepCheck_Analyzer an(sh, Standard_True);
        return an.IsValid() == Standard_True;
    } catch (const Standard_Failure&) {
        return false;
    }
}

void setFaceOutward(TopoDS_Face& f, bool outwardNormal) {
    // Same discipline as trySeamed360: SET orientation from outwardNormal,
    // never toggle on a failed BRepCheck (that undoes a hole's Reverse).
    if (f.IsNull()) return;
    if (!outwardNormal && f.Orientation() == TopAbs_FORWARD) f.Reverse();
    else if (outwardNormal && f.Orientation() == TopAbs_REVERSED) f.Reverse();
}

double wrapToPi(double t) {
    while (t <= -kPi) t += 2.0 * kPi;
    while (t > kPi) t -= 2.0 * kPi;
    return t;
}

// DECISION-p1-growx RULING-b: shift chi into [uMin, uMin+2π). Do not wrap uMax
// into (−π, π] — MakeFace uses the unwrapped interval.
double shiftIntoUSpan(double chi, double uMin, double uMax) {
    const double twoPi = 2.0 * kPi;
    // RULING-b: t = chi + 2π·ceil((uMin-chi)/2π) into [uMin, uMin+2π).
    double k = std::ceil((uMin - chi) / twoPi);
    double t = chi + twoPi * k;
    if (t < uMin) t += twoPi;
    if (t >= uMin + twoPi) t -= twoPi;
    // Prefer the 2π-image inside [uMin, uMax] so a noisy vertex just
    // outside uMin does not jump to uMin+2π and turn a 180° arc into 2π.
    if (t > uMax) {
        const double t2 = t - twoPi;
        const double dHi = t - uMax;
        const double dLo = uMin - t2;
        if (dLo <= dHi) t = t2;
    }
    return t;
}

bool seamStraddleU(const Region& r) {
    if (r.type != SurfType::Cylinder || r.closed360) return false;
    double u0 = r.uMin, u1 = r.uMax;
    if (u1 < u0) u1 += 2.0 * kPi;
    // Crosses the U=0 seam in the interior (S05 r9). A 180° slot centred on
    // U=0 with uMin exactly -π/2 (hand slotted_stadium) is NOT this case.
    return (u0 < -0.5 * kPi && u1 > 0.0) || (u0 < 2.0 * kPi && u1 > 2.0 * kPi);
}

// Mid-vertex selector for partial-cylinder cap arcs when the chain has only
// two terminals (no mesh interior point). Uses the fitted u-span midpoint on
// the iso-circle at the cap V — not OCCT's shortest arc between endpoints.
gp_Pnt partialCylCapArcMid(const Region& cyl, const Region* plnR, const gp_Circ& circ) {
    double v = gp_Vec(cyl.ax.Location(), circ.Location()).Dot(cyl.ax.Direction());
    if (plnR && planePerpCylinder(*plnR, cyl)) v = planeVOnCylinder(*plnR, cyl);
    double u0 = cyl.uMin, u1 = cyl.uMax;
    if (u1 < u0) u1 += 2.0 * kPi;
    double um = 0.5 * (u0 + u1);
    double ang = um;
    if (!seamStraddleU(cyl)) {
        while (ang > kPi) ang -= 2.0 * kPi;
        while (ang <= -kPi) ang += 2.0 * kPi;
    }
    gp_Vec rad(gp_Vec(cyl.ax.XDirection()) * std::cos(ang) +
               gp_Vec(cyl.ax.YDirection()) * std::sin(ang));
    return cyl.ax.Location().Translated(gp_Vec(cyl.ax.Direction()) * v + rad * cyl.radius);
}

double regionU(const Region& r, const gp_Pnt& p) {
    double chi = wrapToPi(azimuthOf(r, p));
    // RULING-b unwrapped sheet only when the fitted interval straddles the
    // cylinder seam; otherwise (−π, π] (r8 / quarter blends / half-holes).
    if (seamStraddleU(r)) return shiftIntoUSpan(chi, r.uMin, r.uMax);
    return chi;
}

double regionV(const Region& r, const gp_Pnt& p) {
    return gp_Vec(r.ax.Location(), p).Dot(r.ax.Direction());
}

double azimuthOnAx(const gp_Ax3& ax, const gp_Pnt& p) {
    gp_Vec rho(ax.Location(), p);
    gp_Vec rad = rho - gp_Vec(ax.Direction()) * rho.Dot(ax.Direction());
    if (rad.Magnitude() < Precision::Confusion()) return 0.0;
    double u = std::atan2(rad.Dot(ax.YDirection()), rad.Dot(ax.XDirection()));
    if (u < 0.0) u += 2.0 * kPi;
    if (u > kPi) u -= 2.0 * kPi;
    return u;
}

double vOnAx(const gp_Ax3& ax, const gp_Pnt& p) {
    return gp_Vec(ax.Location(), p).Dot(ax.Direction());
}

double cylRadialDev(const Region& r, const gp_Pnt& p) {
    gp_Vec rho(r.ax.Location(), p);
    gp_Vec ax(r.ax.Direction());
    gp_Vec rad = rho - ax * rho.Dot(ax);
    return std::fabs(rad.Magnitude() - r.radius);
}

gp_Pnt projectPntOnCylinder(const Region& r, const gp_Pnt& p) {
    gp_Vec rho(r.ax.Location(), p);
    gp_Vec ax(r.ax.Direction());
    double v = rho.Dot(ax);
    gp_Vec rad = rho - ax * v;
    const double mag = rad.Magnitude();
    if (mag < Precision::Confusion())
        rad = gp_Vec(r.ax.XDirection()) * r.radius;
    else
        rad.Scale(r.radius / mag);
    return r.ax.Location().Translated(ax * v + rad);
}

bool edgeUsesLinearCylPCurve(const TopoDS_Edge& e) {
    Standard_Real f = 0, l = 0;
    Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
    if (c.IsNull()) return true;
    if (c->DynamicType() == STANDARD_TYPE(Geom_Ellipse)) return false;
    if (c->DynamicType() == STANDARD_TYPE(Geom_TrimmedCurve)) {
        Handle(Geom_TrimmedCurve) tc = Handle(Geom_TrimmedCurve)::DownCast(c);
        if (!tc.IsNull() && tc->BasisCurve()->DynamicType() == STANDARD_TYPE(Geom_Ellipse))
            return false;
    }
    return true;
}

void rebindCirclePCurvesOnWire(TopoDS_Wire& w, const Handle(Geom_Surface)& surf, double sewTol) {
    for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
        TopoDS_Edge eW = ex.Current();
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c = BRep_Tool::Curve(eW, f, l);
        if (!c.IsNull() && c->DynamicType() == STANDARD_TYPE(Geom_Circle))
            addPcurvesOnSurface(surf, eW, false, sewTol);
    }
}

// D3H-1: construct-once per-(Region, variant) surface table. Faces fetch, never construct.
enum class SurfVar : int {
    Plane = 0,
    CylBase = 1,
    CylRotAx = 2,
    CylRotU1 = 3,
    CylRectTrim = 4,
    CylRotTrim = 5
};

const char* surfVarName(SurfVar v) {
    switch (v) {
    case SurfVar::Plane: return "Plane";
    case SurfVar::CylBase: return "CylBase";
    case SurfVar::CylRotAx: return "CylRotAx";
    case SurfVar::CylRotU1: return "CylRotU1";
    case SurfVar::CylRectTrim: return "CylRectTrim";
    case SurfVar::CylRotTrim: return "CylRotTrim";
    }
    return "?";
}

thread_local std::unordered_map<long long, Handle(Geom_Surface)> gRegionSurf;
thread_local std::unordered_map<const void*, std::pair<int, SurfVar>> gSurfIdent;
thread_local std::unordered_set<const void*> gMeshTShapes;

long long regionSurfKey(int id, SurfVar v) { return (static_cast<long long>(id) << 8) | static_cast<int>(v); }

void registerRegionSurf(int id, SurfVar v, const Handle(Geom_Surface)& s) {
    if (s.IsNull()) return;
    gRegionSurf[regionSurfKey(id, v)] = s;
    gSurfIdent[static_cast<const void*>(s.get())] = {id, v};
}

Handle(Geom_Surface) regionSurf(const Region& r, SurfVar v) {
    const long long k = regionSurfKey(r.id, v);
    auto it = gRegionSurf.find(k);
    if (it != gRegionSurf.end()) return it->second;
    Handle(Geom_Surface) s;
    try {
        if (v == SurfVar::Plane) {
            if (r.type == SurfType::Plane) s = new Geom_Plane(asPlane(r));
        } else if (r.type == SurfType::Cylinder) {
            double u0 = r.uMin, u1 = r.uMax;
            if (u1 < u0) u1 += 2.0 * kPi;
            const double span = u1 - u0;
            if (v == SurfVar::CylBase) {
                s = new Geom_CylindricalSurface(asCyl(r));
            } else if (v == SurfVar::CylRotAx) {
                s = cylSurfaceForRegion(r);
            } else if (v == SurfVar::CylRotU1) {
                gp_Ax3 ax = r.ax;
                ax.Rotate(gp_Ax1(ax.Location(), ax.Direction()), u1);
                s = new Geom_CylindricalSurface(gp_Cylinder(ax, r.radius));
            } else if (v == SurfVar::CylRectTrim) {
                Handle(Geom_CylindricalSurface) base =
                    Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylBase));
                if (!base.IsNull() && span > Precision::Confusion())
                    s = new Geom_RectangularTrimmedSurface(base, u0, u1, r.vMin, r.vMax);
            } else if (v == SurfVar::CylRotTrim) {
                Handle(Geom_CylindricalSurface) basis;
                if (!r.outwardNormal && r.uMin >= -1e-10 && u1 <= 0.5 * kPi + 0.02)
                    basis = Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylRotU1));
                else
                    basis = Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylRotAx));
                if (!basis.IsNull() && span > Precision::Confusion())
                    s = new Geom_RectangularTrimmedSurface(basis, 0.0, span, r.vMin, r.vMax);
            }
        }
    } catch (const Standard_Failure&) {
        s.Nullify();
    }
    registerRegionSurf(r.id, v, s);
    return s;
}

void variantsForRegion(const Region& r, SurfVar* out, int& n) {
    n = 0;
    if (r.type == SurfType::Plane) {
        out[n++] = SurfVar::Plane;
        return;
    }
    if (r.type != SurfType::Cylinder) return;
    out[n++] = SurfVar::CylBase;
    out[n++] = SurfVar::CylRotAx;
    out[n++] = SurfVar::CylRotU1;
    if (!r.closed360) {
        out[n++] = SurfVar::CylRectTrim;
        out[n++] = SurfVar::CylRotTrim;
    }
}

Handle(Geom_Surface) untrimmedBasisOf(const Handle(Geom_Surface)& srf) {
    if (srf.IsNull()) return {};
    Handle(Geom_RectangularTrimmedSurface) t = Handle(Geom_RectangularTrimmedSurface)::DownCast(srf);
    if (!t.IsNull() && !t->BasisSurface().IsNull()) return t->BasisSurface();
    return srf;
}

Handle(Geom_CylindricalSurface) basisCylOf(const Handle(Geom_Surface)& srf) {
    Handle(Geom_CylindricalSurface) c = Handle(Geom_CylindricalSurface)::DownCast(srf);
    if (!c.IsNull()) return c;
    return Handle(Geom_CylindricalSurface)::DownCast(untrimmedBasisOf(srf));
}

Handle(Geom_Plane) basisPlaneOf(const Handle(Geom_Surface)& srf) {
    Handle(Geom_Plane) pl = Handle(Geom_Plane)::DownCast(srf);
    if (!pl.IsNull()) return pl;
    return Handle(Geom_Plane)::DownCast(untrimmedBasisOf(srf));
}

Handle(Geom_Curve) basisCurveOf(const Handle(Geom_Curve)& c) {
    Handle(Geom_TrimmedCurve) tr = Handle(Geom_TrimmedCurve)::DownCast(c);
    if (!tr.IsNull() && !tr->BasisCurve().IsNull()) return tr->BasisCurve();
    return c;
}

double pcurveDev(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                 const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d);

Handle(Geom2d_Curve) makePCurveOnSurf(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                                      const Handle(Geom_Surface)& srf, const char* kind) {
    if (c3.IsNull() || srf.IsNull()) return {};
    Handle(Geom_Plane) pl = Handle(Geom_Plane)::DownCast(srf);
    if (!pl.IsNull()) {
        Handle(Geom_Curve) src = basisCurveOf(c3);
        try {
            return GeomAPI::To2d(src, pl->Pln());
        } catch (const Standard_Failure&) {
            return {};
        }
    }
    Handle(Geom_Curve) src = basisCurveOf(c3);
    const bool isCirc = (kind && std::strcmp(kind, "circ") == 0) ||
                        (!src.IsNull() && src->DynamicType() == STANDARD_TYPE(Geom_Circle));
    const bool isLin = (kind && std::strcmp(kind, "lin") == 0) ||
                       (!src.IsNull() && src->DynamicType() == STANDARD_TYPE(Geom_Line));
    const bool isElips = (kind && std::strcmp(kind, "elips") == 0) ||
                         (!src.IsNull() && src->DynamicType() == STANDARD_TYPE(Geom_Ellipse));
    Handle(Geom_CylindricalSurface) cyl = basisCylOf(srf);
    if (cyl.IsNull()) return {};
    const gp_Ax3 pos = cyl->Position();
    try {
        if (isCirc && !isElips) {
            Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(src);
            if (gc.IsNull()) return {};
            const gp_Circ circ = gc->Circ();
            const gp_Dir cx = circ.XAxis().Direction();
            const double phi = std::atan2(cx.Dot(pos.YDirection()), cx.Dot(pos.XDirection()));
            const double v0 = gp_Vec(pos.Location(), circ.Location()).Dot(pos.Direction());
            const double sgn = circ.Axis().Direction().Dot(pos.Direction()) >= 0.0 ? 1.0 : -1.0;
            return new Geom2d_Line(gp_Pnt2d(phi, v0), gp_Dir2d(sgn, 0.0));
        }
        if (isLin) {
            Handle(Geom_Line) gl = Handle(Geom_Line)::DownCast(src);
            Handle(Geom2d_Curve) c2dLin;
            if (!gl.IsNull()) {
                const gp_Lin lin = gl->Lin();
                const double zdot = lin.Direction().Dot(pos.Direction());
                if (std::fabs(std::fabs(zdot) - 1.0) <= 1e-9) {
                    const gp_Pnt pmid = c3->Value(0.5 * (f + l));
                    gp_Vec rho(pos.Location(), pmid);
                    gp_Vec rad = rho - gp_Vec(pos.Direction()) * rho.Dot(pos.Direction());
                    double u0 = 0.0;
                    if (rad.Magnitude() > Precision::Confusion())
                        u0 = std::atan2(rad.Dot(pos.YDirection()), rad.Dot(pos.XDirection()));
                    gp_Vec rho0(pos.Location(), lin.Location());
                    const double v0 = rho0.Dot(pos.Direction());
                    const double sgn = zdot >= 0.0 ? 1.0 : -1.0;
                    c2dLin = new Geom2d_Line(gp_Pnt2d(u0, v0), gp_Dir2d(0.0, sgn));
                }
            }
            if (!c2dLin.IsNull() && pcurveDev(c3, f, l, srf, c2dLin) <= Precision::Confusion())
                return c2dLin;
            return GeomProjLib::Curve2d(c3, f, l, srf);
        }
        if (isElips || !src.IsNull()) {
            return GeomProjLib::Curve2d(c3, f, l, srf);
        }
    } catch (const Standard_Failure&) {
        return {};
    }
    return {};
}

void pcurveDevTN(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                 const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d, double& devT,
                 double& devN, double* devOut = nullptr, const char** projOut = nullptr,
                 bool* unprojOut = nullptr) {
    devT = 0.0;
    devN = 0.0;
    double dev = 0.0;
    const char* proj = "closedform";
    bool unproj = false;
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) {
        devT = 1e300;
        devN = 1e300;
        dev = 1e300;
        proj = "iterative";
        unproj = true;
        if (devOut) *devOut = dev;
        if (projOut) *projOut = proj;
        if (unprojOut) *unprojOut = unproj;
        return;
    }
    const int N = 16;
    const double df = l - f;
    if (!(df > Precision::PConfusion())) {
        if (devOut) *devOut = 0.0;
        if (projOut) *projOut = proj;
        if (unprojOut) *unprojOut = false;
        return;
    }
    try {
        // P* on the BASIS (untrimmed) surface. Evaluating Q on the table
        // handle S (possibly trimmed) is required; projecting onto a
        // RectangularTrimmedSurface clamps to the trim and inflates devT.
        // pproj = Q is forbidden.
        Handle(Geom_Surface) basis = untrimmedBasisOf(srf);
        Handle(Geom_Plane) pl = basisPlaneOf(srf);
        Handle(Geom_CylindricalSurface) cyl = basisCylOf(srf);
        if (pl.IsNull() && cyl.IsNull()) proj = "iterative";
        for (int i = 0; i < N; ++i) {
            const double t = f + df * static_cast<double>(i) / static_cast<double>(N - 1);
            const gp_Pnt P = c3->Value(t);  // stored C3d; never BasisCurve / mesh
            const gp_Pnt2d uv = c2d->Value(t);
            const gp_Pnt Q = srf->Value(uv.X(), uv.Y());
            gp_Pnt Pstar = P;
            double n = 0.0;
            bool got = false;
            if (!pl.IsNull()) {
                const gp_Pln g = pl->Pln();
                n = std::fabs(g.Distance(P));
                const gp_Vec nn(g.Axis().Direction());
                const double signedD = gp_Vec(g.Location(), P).Dot(nn);
                Pstar = P.Translated(nn * (-signedD));
                got = true;
            } else if (!cyl.IsNull()) {
                const gp_Ax3 pos = cyl->Position();
                const double R = cyl->Radius();
                gp_Vec rho(pos.Location(), P);
                const double v = rho.Dot(pos.Direction());
                gp_Vec rad = rho - gp_Vec(pos.Direction()) * v;
                const double rmag = rad.Magnitude();
                n = std::fabs(rmag - R);
                if (rmag > Precision::Confusion())
                    Pstar = pos.Location().Translated(gp_Vec(pos.Direction()) * v +
                                                      rad * (R / rmag));
                else
                    Pstar = pos.Location().Translated(gp_Vec(pos.Direction()) * v +
                                                      gp_Vec(pos.XDirection()) * R);
                got = true;
            } else if (!basis.IsNull()) {
                GeomAPI_ProjectPointOnSurf projector(P, basis);
                if (projector.IsDone() && projector.NbPoints() >= 1) {
                    Pstar = projector.NearestPoint();
                    n = P.Distance(Pstar);
                    got = true;
                }
            }
            if (!got) {
                unproj = true;
                continue;
            }
            devN = std::max(devN, n);
            devT = std::max(devT, Q.Distance(Pstar));
            dev = std::max(dev, P.Distance(Q));
        }
        if (unproj && proj[0] == 'c') proj = "iterative";
    } catch (const Standard_Failure&) {
        unproj = true;
        proj = "iterative";
    }
    if (devOut) *devOut = dev;
    if (projOut) *projOut = proj;
    if (unprojOut) *unprojOut = unproj;
}

double pcurveDev(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                 const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d) {
    double t = 0, n = 0;
    pcurveDevTN(c3, f, l, srf, c2d, t, n);
    return std::max(t, n);
}

double bindTolFromDev(double dev, double cap) {
    double t = Precision::Confusion();
    if (std::isfinite(dev) && dev > t) t = dev;
    if (t > cap) t = cap;
    return t;
}

void emitFaceSurfDiag(int ridHint, const TopoDS_Face& F) {
    if (!diagP2Enabled() || F.IsNull()) return;
    TopLoc_Location loc;
    Handle(Geom_Surface) s = BRep_Tool::Surface(F, loc);
    int rid = ridHint;
    int fromTable = 0;
    const char* varName = "unknown";
    if (!s.IsNull()) {
        auto it = gSurfIdent.find(static_cast<const void*>(s.get()));
        if (it != gSurfIdent.end()) {
            fromTable = 1;
            rid = it->second.first;
            varName = surfVarName(it->second.second);
        }
    }
    const bool facePlane = (s && Handle(Geom_Plane)::DownCast(s));
    int nEdges = 0, nHasPc = 0, nOrphan = 0, nSynth = 0, nI = 0, nII = 0, nIII = 0;
    for (TopExp_Explorer wx(F, TopAbs_WIRE); wx.More(); wx.Next()) {
        const TopoDS_Wire w = TopoDS::Wire(wx.Current());
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            nEdges++;
            const TopoDS_Edge e = TopoDS::Edge(it.Value());
            Standard_Real f = 0, l = 0;
            Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
            Standard_Boolean stored = Standard_False;
            Standard_Real pf = 0, pl = 0;
            if (!s.IsNull()) (void)BRep_Tool::CurveOnSurface(e, s, loc, pf, pl, &stored);
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) synth = BRep_Tool::CurveOnSurface(e, F, a, b);
            if (!stored && !synth.IsNull()) nSynth++;
            if (stored) nHasPc++;
            else if (!c3.IsNull() && (l - f) > Precision::PConfusion()) {
                nOrphan++;
                const void* ts = diagTShapePtr(e);
                if (gSeamTShapes.count(ts)) nI++;
                else if (facePlane && gMeshTShapes.count(ts)) nII++;
                else nIII++;
            }
        }
    }
    std::fprintf(stderr,
                 "DIAG_FACESURF rid=%d variant=%s fromTable=%d nEdges=%d nHasPc=%d nOrphan=%d "
                 "nSynth=%d nI=%d nII=%d nIII=%d\n",
                 rid, varName, fromTable, nEdges, nHasPc, nOrphan, nSynth, nI, nII, nIII);
}

struct SeamBindCounts {
    int nWrite = 0;
    int nGuardHit = 0;
    int nDoubleWrite = 0;
};

void bindAllVariants(const TopoDS_Edge& e, const Region& R, const char* kind, const MeshView& mv,
                     int ci, SeamBindCounts& acc) {
    Standard_Real f = 0, l = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
    if (c3.IsNull() || l - f <= Precision::PConfusion()) return;
    if (const void* ts = diagTShapePtr(e)) gSeamTShapes.insert(ts);
    const TopLoc_Location loc;
    const double cap = meshTolCap(mv, &R);
    SurfVar vars[6];
    int nv = 0;
    variantsForRegion(R, vars, nv);
    for (int i = 0; i < nv; ++i) {
        Handle(Geom_Surface) srf = regionSurf(R, vars[i]);
        if (srf.IsNull()) continue;
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        (void)BRep_Tool::CurveOnSurface(e, srf, loc, pf, pl, &hasPc);
        if (hasPc) {
            acc.nGuardHit++;
            continue;
        }
        Handle(Geom2d_Curve) c2d = makePCurveOnSurf(c3, f, l, srf, kind);
        if (c2d.IsNull()) continue;
        double devT = 0, devN = 0, dev = 0;
        const char* proj = "closedform";
        bool unproj = false;
        pcurveDevTN(c3, f, l, srf, c2d, devT, devN, &dev, &proj, &unproj);
        if (diagP2Enabled() &&
            (R.type == SurfType::Cylinder || R.type == SurfType::Plane))
            std::fprintf(stderr,
                         "DIAG_PCDEV ci=%d rid=%d variant=%s kind=%s dev=%.6g devT=%.6g "
                         "devN=%.6g proj=%s\n",
                         ci, R.id, surfVarName(vars[i]),
                         unproj ? "unproj" : (kind ? kind : "poly"), dev, devT, devN,
                         proj ? proj : "iterative");
        double tol = bindTolFromDev(std::max(devT, devN), cap);
        const double existing = BRep_Tool::Tolerance(e);
        if (existing > tol) tol = existing;
        BRep_Builder B;
        B.UpdateEdge(e, c2d, srf, loc, tol);
        B.Range(e, srf, loc, f, l);
        acc.nWrite++;
    }
}

// U3: consumer/assertion. Never creates. nSkippedNoCurve = no 3D curve only (D3H-3).
void bindCylPCurves(TopoDS_Wire& w, const Handle(Geom_Surface)& surf, const Region& r,
                    double /*sewTol*/) {
    if (w.IsNull() || surf.IsNull()) return;
    const TopLoc_Location loc;
    int nEdgesIter = 0, nEdgesExp = 0;
    int nBound = 0, nSkippedHasPc = 0, nSkippedNoCurve = 0, nOrphan = 0;
    const int nCreated = 0;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() == TopAbs_EDGE) nEdgesIter++;
    }
    for (TopExp_Explorer cex(w, TopAbs_EDGE); cex.More(); cex.Next()) nEdgesExp++;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() != TopAbs_EDGE) continue;
        TopoDS_Edge eW = TopoDS::Edge(it.Value());
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        (void)BRep_Tool::CurveOnSurface(eW, surf, loc, pf, pl, &hasPc);
        if (hasPc) {
            nSkippedHasPc++;
            continue;
        }
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(eW, f, l);
        if (c3.IsNull() || l - f <= Precision::PConfusion()) {
            nSkippedNoCurve++;
            continue;
        }
        nOrphan++;
    }
    if (diagP2Enabled()) {
        const int covered = nBound + nSkippedHasPc + nSkippedNoCurve + nOrphan;
        const int g7 = (covered == nEdgesIter) ? 1 : 0;
        std::fprintf(stderr,
                     "DIAG_CYLG7 rid=%d nEdgesIter=%d nEdgesExp=%d nBound=%d nSkippedHasPc=%d "
                     "nSkippedNoCurve=%d covered=%d g7=%d nCreated=%d nOrphan=%d\n",
                     r.id, nEdgesIter, nEdgesExp, nBound, nSkippedHasPc, nSkippedNoCurve, covered,
                     g7, nCreated, nOrphan);
    }
}

// U4: plane-side consumer/assertion (D3H-3/D3H-6). Never creates.
void bindPlanePCurves(TopoDS_Wire& w, const Handle(Geom_Plane)& thePlane, const Region& r,
                      double /*tol*/, int& nBound, int& nSkippedHasPc, int& nSkippedNoCurve) {
    if (w.IsNull() || thePlane.IsNull()) return;
    const TopLoc_Location loc;
    int nEdgesIter = 0, nEdgesExp = 0;
    int wBound = 0, wSkipPc = 0, wSkipNc = 0, wOrphan = 0;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() == TopAbs_EDGE) nEdgesIter++;
    }
    for (TopExp_Explorer cex(w, TopAbs_EDGE); cex.More(); cex.Next()) nEdgesExp++;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() != TopAbs_EDGE) continue;
        TopoDS_Edge eW = TopoDS::Edge(it.Value());
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        (void)BRep_Tool::CurveOnSurface(eW, thePlane, loc, pf, pl, &hasPc);
        if (hasPc) {
            nSkippedHasPc++;
            wSkipPc++;
            continue;
        }
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(eW, f, l);
        if (c3.IsNull() || l - f <= Precision::PConfusion()) {
            nSkippedNoCurve++;
            wSkipNc++;
            continue;
        }
        wOrphan++;
    }
    (void)wBound;
    (void)nBound;
    if (diagP2Enabled()) {
        const int covered = wBound + wSkipPc + wSkipNc + wOrphan;
        const int g7 = (covered == nEdgesIter) ? 1 : 0;
        std::fprintf(stderr,
                     "DIAG_G7 rid=%d nEdgesIter=%d nEdgesExp=%d nBound=%d nSkippedHasPc=%d "
                     "nSkippedNoCurve=%d covered=%d g7=%d nOrphan=%d\n",
                     r.id, nEdgesIter, nEdgesExp, wBound, wSkipPc, wSkipNc, covered, g7, wOrphan);
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
            if (ip->chainIdx.size() == 1) {
                const int ci = ip->chainIdx[0];
                if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                    (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                    geom[(size_t)ci].edges.size() == 1 &&
                    edgeSpansFullCircle(geom[(size_t)ci].edges[0]))
                    iw = TopoDS::Wire(iw.Reversed());
            }
            if (ip->chainIdx.size() > 6) {
                try {
                    BRepBuilderAPI_MakeWire mw;
                    for (TopoDS_Iterator it(iw); it.More(); it.Next())
                        mw.Add(TopoDS::Edge(it.Value()));
                    if (mw.IsDone()) {
                        TopoDS_Wire iw2 = mw.Wire();
                        iw2.Closed(Standard_True);
                        iw = iw2;
                    }
                } catch (const Standard_Failure&) {
                }
            }
            mf.Add(iw);
        }
        return mf.IsDone() == Standard_True;
    };
    try {
        Handle(Geom_Plane) gpl = Handle(Geom_Plane)::DownCast(regionSurf(r, SurfVar::Plane));
        std::vector<TopoDS_Wire> innerW;
        bool innersOk = true;
        for (const Loop* ip : inners) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, *ip, rs, mv, geom, collapsed, meshE, edgeOk)) {
                innersOk = false;
                break;
            }
            // Hole loops must bind shared Seamed360 cap circles with reversed
            // wire orientation (same TShape, flipped in the wire).
            if (ip->chainIdx.size() == 1) {
                const int ci = ip->chainIdx[0];
                if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                    (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                    geom[(size_t)ci].edges.size() == 1 &&
                    edgeSpansFullCircle(geom[(size_t)ci].edges[0]))
                    iw = TopoDS::Wire(iw.Reversed());
            }
            innerW.push_back(iw);
        }
        if (innersOk) {
            int nBound = 0, nSkippedHasPc = 0, nSkippedNoCurve = 0;
            bindPlanePCurves(ow, gpl, r, 0.0, nBound, nSkippedHasPc, nSkippedNoCurve);
            for (auto& iw : innerW)
                bindPlanePCurves(iw, gpl, r, 0.0, nBound, nSkippedHasPc, nSkippedNoCurve);
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_PLANEPCBIND rid=%d nBound=%d nSkippedHasPc=%d nSkippedNoCurve=%d\n",
                             r.id, nBound, nSkippedHasPc, nSkippedNoCurve);
        }
        if (innersOk && makeFaceKeep(gpl, ow, innerW, r.outwardNormal, outF)) {
            if (faceIsValid(outF) || ensureFaceValid(outF, meshTolCap(mv, &r)))
                return !outF.IsNull();
            diagPlateMakeFaceFail(r.id, outF);
        } else if (innersOk) {
            diagPlateMakeFaceFail(r.id, outF);
        }
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
        if (!faceIsValid(outF) && inners.empty()) {
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
        ensureFaceValid(outF, meshTolCap(mv, &r));
        return !outF.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool buildPartialCylinder(const Region& r, const RegionSet& rs, const MeshView& mv,
                          const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                          const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                          double sewTol, TopoDS_Face& outF) {
    auto diagPartial = [&](const char* why) {
        if (collapseDiagEnabled())
            std::fprintf(stderr, "DIAG_PARTIAL rid=%d %s\n", r.id, why);
    };
    const char* surfTag = "unset";
    const char* wireTag = "fwd";
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
    if (!outer) {
        diagPartial("no-outer");
        return false;
    }
    TopoDS_Wire ow;
    auto refreshOuterWire = [&]() -> bool {
        if (!buildLoopWire(ow, *outer, rs, mv, geom, collapsed, meshE, edgeOk)) return false;
        if (seamStraddleU(r) || outer->chainIdx.size() > 6) {
            try {
                BRepBuilderAPI_MakeWire mw;
                for (TopoDS_Iterator it(ow); it.More(); it.Next())
                    mw.Add(TopoDS::Edge(it.Value()));
                if (mw.IsDone()) {
                    TopoDS_Wire ow2 = mw.Wire();
                    ow2.Closed(Standard_True);
                    ow = ow2;
                }
            } catch (const Standard_Failure&) {
            }
        }
        return true;
    };
    if (!refreshOuterWire()) {
        diagPartial("loop-wire");
        return false;
    }
    if (collapseDiagEnabled()) {
        double u0 = r.uMin, u1 = r.uMax;
        if (u1 < u0) u1 += 2.0 * kPi;
        int nOuterCh = outer ? (int)outer->chainIdx.size() : 0;
        std::fprintf(stderr,
                     "DIAG_PARTIAL_RID rid=%d R=%.4f nTri=%d outward=%d closed360=%d "
                     "u=[%.4f,%.4f] span=%.4f v=[%.4f,%.4f] nOuterCh=%d nInner=%zu "
                     "dVol=%.3f nSides=%d\n",
                     r.id, r.radius, (int)r.tris.size(), r.outwardNormal ? 1 : 0,
                     r.closed360 ? 1 : 0, r.uMin, r.uMax, u1 - u0, r.vMin, r.vMax, nOuterCh,
                     inners.size(), r.dVolPredicted, r.nSides);
    }
    const double cap = std::max(partialFaceTolCap(mv, r), analyticSnapCap(mv, &r, nullptr));
    auto finishPartial = [&](TopoDS_Face cand) -> bool {
        if (cand.IsNull()) return false;
        char tagbuf[96];
        std::snprintf(tagbuf, sizeof(tagbuf), "%s/%s", surfTag, wireTag);
        dumpFaceTopo(r.id, cand, tagbuf, "raw");
        {
            TopoDS_Face want = cand, other = cand;
            setFaceOutward(want, r.outwardNormal);
            setFaceOutward(other, !r.outwardNormal);
            const bool wantOk = faceIsValid(want);
            const bool otherOk = faceIsValid(other);
            if (wantOk) cand = want;
            else if (otherOk) cand = other;
            else cand = want;
            if (collapseDiagEnabled() && otherOk && !wantOk)
                std::fprintf(stderr, "DIAG_PARTIAL rid=%d flipped-outward tag=%s\n", r.id, tagbuf);
        }
        if (!faceIsValid(cand)) {
            dumpFaceTopo(r.id, cand, tagbuf, "pre-fix");
            addPcurvesOnFace(cand, sewTol, false);
        }
        if (!faceIsValid(cand)) {
            try {
                ShapeFix_Face sff(cand);
                sff.FixOrientationMode() = 1;
                sff.FixAddNaturalBoundMode() = 0;
                sff.FixMissingSeamMode() = seamStraddleU(r) ? 1 : 0;
                sff.FixWireMode() = 1;
                sff.Perform();
                TopoDS_Shape res = sff.Result();
                if (res.IsNull()) res = sff.Face();
                for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                    cand = TopoDS::Face(fx.Current());
                    break;
                }
            } catch (const Standard_Failure&) {
            }
            setFaceOutward(cand, r.outwardNormal);
            addPcurvesOnFace(cand, sewTol, false);
            dumpFaceTopo(r.id, cand, tagbuf, "post-fix");
        }
        if (!ensureFaceValid(cand, cap)) {
            dumpFaceTopo(r.id, cand, tagbuf, "ensure");
            if (collapseDiagEnabled()) {
                double maxD = 0.0;
                for (TopExp_Explorer vx(cand, TopAbs_VERTEX); vx.More(); vx.Next())
                    maxD = std::max(maxD, cylRadialDev(r, BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()))));
                std::fprintf(stderr,
                             "DIAG_PARTIAL rid=%d maxCylDev=%.4f cap=%.4f valid=%d tag=%s", r.id,
                             maxD, cap, faceIsValid(cand) ? 1 : 0, tagbuf);
                try {
                    BRepCheck_Analyzer an(cand, Standard_True);
                    if (!an.IsValid()) dumpBRepStatuses(an.Result(cand));
                } catch (const Standard_Failure&) {
                }
                std::fprintf(stderr, "\n");
            }
            diagPartial(faceIsValid(cand) ? "ensure-only" : "ensure-invalid");
            return false;
        }
        if (!faceIsValid(cand)) {
            dumpFaceTopo(r.id, cand, tagbuf, "face-invalid");
            diagPartial("face-invalid");
            return false;
        }
        if (collapseDiagEnabled())
            std::fprintf(stderr, "DIAG_PARTIAL rid=%d valid=1 tag=%s ori=%s\n", r.id, tagbuf,
                         cand.Orientation() == TopAbs_FORWARD ? "F" : "R");
        outF = cand;
        return true;
    };
    auto makeFaceBound = [&](const Handle(Geom_Surface)& surf, TopoDS_Face& cand) -> bool {
        if (mv.nTri >= 500 && !refreshOuterWire()) return false;
        bindCylPCurves(ow, surf, r, sewTol);
        for (auto& iw : inners) bindCylPCurves(iw, surf, r, sewTol);
        auto fillMissingPc = [&](TopoDS_Wire& ww) {
            if (gCompNTri >= 10000) return;
            if (ww.IsNull() || surf.IsNull()) return;
            const TopLoc_Location loc;
            BRep_Builder Bb;
            for (TopoDS_Iterator it(ww); it.More(); it.Next()) {
                if (it.Value().ShapeType() != TopAbs_EDGE) continue;
                TopoDS_Edge eW = TopoDS::Edge(it.Value());
                Standard_Boolean hasPc = Standard_False;
                Standard_Real pf = 0, pl = 0;
                (void)BRep_Tool::CurveOnSurface(eW, surf, loc, pf, pl, &hasPc);
                if (hasPc) continue;
                Standard_Real f = 0, l = 0;
                Handle(Geom_Curve) c3 = BRep_Tool::Curve(eW, f, l);
                if (c3.IsNull() || l - f <= Precision::PConfusion()) continue;
                Handle(Geom2d_Curve) c2d = makePCurveOnSurf(c3, f, l, surf, nullptr);
                if (c2d.IsNull()) continue;
                double tol = Precision::Confusion();
                const double existing = BRep_Tool::Tolerance(eW);
                if (existing > tol) tol = existing;
                Bb.UpdateEdge(eW, c2d, surf, loc, tol);
                Bb.Range(eW, surf, loc, f, l);
            }
        };
        fillMissingPc(ow);
        for (auto& iw : inners) fillMissingPc(iw);
        auto attempt = [&]() -> bool {
            try {
                BRepBuilderAPI_MakeFace mf(surf, ow, Standard_True);
                if (!mf.IsDone()) return false;
                for (const auto& iw : inners) mf.Add(iw);
                if (!mf.IsDone()) return false;
                cand = mf.Face();
                emitFaceSurfDiag(r.id, cand);
                return !cand.IsNull();
            } catch (const Standard_Failure&) {
                return false;
            }
        };
        if (attempt()) return true;
        rebindCirclePCurvesOnWire(ow, surf, sewTol);
        for (auto& iw : inners) rebindCirclePCurvesOnWire(iw, surf, sewTol);
        return attempt();
    };
    auto tryPartial = [&](const Handle(Geom_CylindricalSurface)& surf) -> bool {
        TopoDS_Face cand;
        bool got = makeFaceBound(surf, cand);
        if (!got) {
            if (mv.nTri >= 500 && !refreshOuterWire()) {
                diagPartial("makeface-null");
                return false;
            }
            got = makeFaceKeep(surf, ow, inners, r.outwardNormal, cand);
        }
        if (!got) {
            diagPartial("makeface-null");
            return false;
        }
        return finishPartial(cand);
    };
    auto wantRotTrim = [&]() -> bool {
        if (mv.nTri < 500) return false;
        double u0 = r.uMin, u1 = r.uMax;
        if (u1 < u0) u1 += 2.0 * kPi;
        const double span = u1 - u0;
        if (span <= Precision::Confusion()) return false;
        return seamStraddleU(r) ||
               (r.uMin < -1e-10 && r.uMax > 1e-10 && span < kPi) ||
               (!r.outwardNormal && span <= 0.5 * kPi + 0.02) ||
               (r.radius > 0.0 && r.radius < 2.0 && span <= 0.5 * kPi + 0.02);
    };
    auto tryRotTrimmedSheet = [&]() -> bool {
        surfTag = "rot-trim";
        if (r.closed360 || !wantRotTrim()) return false;
        double u0 = r.uMin, u1 = r.uMax;
        if (u1 < u0) u1 += 2.0 * kPi;
        const double span = u1 - u0;
        if (!refreshOuterWire()) return false;
        try {
            Handle(Geom_CylindricalSurface) surfR;
            // Positive-u quarter holes (rid=20) mirror negative-u (rid=19) via uMax
            // rotation onto the same [0, span] trimmed sheet.
            if (!r.outwardNormal && r.uMin >= -1e-10 && u1 <= 0.5 * kPi + 0.02) {
                surfR = Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylRotU1));
            } else {
                surfR = Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylRotAx));
            }
            Handle(Geom_RectangularTrimmedSurface) trim =
                Handle(Geom_RectangularTrimmedSurface)::DownCast(regionSurf(r, SurfVar::CylRotTrim));
            (void)surfR;
            (void)span;
            if (trim.IsNull()) return false;
            TopoDS_Face cand;
            bool got = makeFaceBound(trim, cand);
            if (!got && (mv.nTri < 500 || refreshOuterWire()))
                got = makeFaceKeep(trim, ow, inners, r.outwardNormal, cand);
            if (!got) return false;
            return finishPartial(cand);
        } catch (const Standard_Failure&) {
            return false;
        }
    };
    Handle(Geom_CylindricalSurface) surf0 =
        Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylBase));
    auto tryAllSurfs = [&](TopoDS_Wire& wire) -> bool {
        TopoDS_Wire owSaved = ow;
        ow = wire;
        if (wantRotTrim() && tryRotTrimmedSheet()) {
            ow = owSaved;
            return true;
        }
        double u0 = r.uMin, u1 = r.uMax;
        if (u1 < u0) u1 += 2.0 * kPi;
        try {
            surfTag = "rect-trim";
            Handle(Geom_RectangularTrimmedSurface) trim =
                Handle(Geom_RectangularTrimmedSurface)::DownCast(regionSurf(r, SurfVar::CylRectTrim));
            TopoDS_Face cand;
            bool got = makeFaceBound(trim, cand);
            if (!got && (mv.nTri < 500 || refreshOuterWire()))
                got = makeFaceKeep(trim, ow, inners, r.outwardNormal, cand);
            if (got && finishPartial(cand)) {
                ow = owSaved;
                return true;
            }
        } catch (const Standard_Failure&) {
        }
        if (seamStraddleU(r)) {
            try {
                surfTag = "seam-box";
                BRepBuilderAPI_MakeFace box(asCyl(r), u0, u1, r.vMin, r.vMax);
                if (box.IsDone()) {
                    Handle(Geom_Surface) sbox = BRep_Tool::Surface(box.Face());
                    TopoDS_Face cand;
                    BRep_Builder Bb;
                    Bb.MakeFace(cand, sbox, Precision::Confusion());
                    Bb.Add(cand, ow);
                    for (const auto& iw : inners) Bb.Add(cand, iw);
                    if (finishPartial(cand)) {
                        ow = owSaved;
                        return true;
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }
        surfTag = "untrim";
        if (tryPartial(surf0)) {
            ow = owSaved;
            return true;
        }
        if (!r.closed360 && std::fabs(r.uMin) > 1e-14) {
            surfTag = "rot-ax";
            Handle(Geom_CylindricalSurface) surf1 =
                Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylRotAx));
            if (tryPartial(surf1)) {
                ow = owSaved;
                return true;
            }
        }
        ow = owSaved;
        return false;
    };
    if (mv.nTri < 500) {
        wireTag = "fwd";
        if (tryAllSurfs(ow)) return true;
    } else {
        if (!r.outwardNormal) {
            if (refreshOuterWire()) {
                TopoDS_Wire rev = ow;
                rev.Reverse();
                wireTag = "rev-hole";
                if (tryAllSurfs(rev)) return true;
            }
        }
        wireTag = "fwd";
        if (refreshOuterWire() && tryAllSurfs(ow)) return true;
    }
    {
        TopoDS_Wire rev = ow;
        rev.Reverse();
        wireTag = "rev";
        if (tryAllSurfs(rev)) return true;
    }
    diagPartial("all-failed");
    return false;
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
        BRepBuilderAPI_MakeEdge me2(circ);
        if (me2.IsDone()) {
            TopoDS_Edge e = me2.Edge();
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

// ME_CYLPLN_DIFF_PTS_CLOSED on an ellipse: MakeEdge(elips, V, V) after
// snapping one vertex, not two distinct projections.
TopoDS_Edge makeFullEllipse(const gp_Elips& el, const TopoDS_Vertex& V) {
    try {
        BRepBuilderAPI_MakeEdge me(el, V, V);
        if (me.IsDone()) {
            TopoDS_Edge e = me.Edge();
            e.Closed(Standard_True);
            return e;
        }
    } catch (const Standard_Failure&) {
    }
    try {
        BRepBuilderAPI_MakeEdge me2(el);
        if (me2.IsDone()) {
            TopoDS_Edge e = me2.Edge();
            e.Closed(Standard_True);
            return e;
        }
    } catch (const Standard_Failure&) {
    }
    try {
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(el);
        BRep_Builder B;
        TopoDS_Edge e;
        B.MakeEdge(e, ge, Precision::Confusion());
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
                  const std::vector<char>& edgeOk, double sewTol, WarnFn warn, TopoDS_Face& outF,
                  const std::vector<char>* exploded = nullptr) {
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

    int vL = vertexClosestToUOnLoop(mv, r, *capL, rs, 0.0);
    int vH = vertexClosestToUOnLoop(mv, r, *capH, rs, 0.0);
    if (vL < 0 || vH < 0 || (size_t)vL >= verts.size() || (size_t)vH >= verts.size()) return false;

    Handle(Geom_CylindricalSurface) surf =
        Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylBase));
    gp_Circ circL = cylinderIsoCircle(r, r.vMin);
    gp_Circ circH = cylinderIsoCircle(r, r.vMax);
    AnalyticCurve acL, acH;
    acL.kind = AnalyticCurve::Circ;
    acL.circ = circL;
    acH.kind = AnalyticCurve::Circ;
    acH.circ = circH;
    const double snapCap = meshTolCap(mv, &r);
    snapVertexToCurve(verts[(size_t)vL], acL, snapCap);
    snapVertexToCurve(verts[(size_t)vH], acH, snapCap);

    auto takeFullCap = [&](int ci, const gp_Circ& circ, const TopoDS_Vertex& V,
                           TopoDS_Edge& e) -> bool {
        if (ci >= 0 && (size_t)ci < geom.size() && geom[(size_t)ci].collapsed) {
            if (geom[(size_t)ci].edges.size() != 1) return false;
            e = TopoDS::Edge(geom[(size_t)ci].edges[0].Oriented(TopAbs_FORWARD));
            return edgeSpansFullCircle(e);
        }
        e = makeFullCircle(circ, V);
        return !e.IsNull();
    };

    auto bindIsoPCurves = [&](TopoDS_Edge& eCap, double vIso, double u0, TopoDS_Edge& eSeam,
                             bool writeSeam) {
        BRep_Builder B;
        // Do not B.Range() a cap circle that is already shared with a neighbour
        // plane (J2). Range mutates the 3d bounds and opens the shell (J6)
        // when the plate faces were built first. Pcurves follow the existing
        // 3d range; u0 is the seam azimuth used only for the seam pcurves.
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(eCap, f, l);
        if (c3.IsNull()) {
            f = u0;
            l = u0 + 2.0 * kPi;
        }
        Handle(Geom2d_Line) pc = new Geom2d_Line(gp_Pnt2d(f, vIso), gp_Dir2d(1.0, 0.0));
        B.UpdateEdge(eCap, pc, surf, TopLoc_Location(), sewTol);
        if (!writeSeam) return;
        Standard_Real fs = 0, ls = 0;
        Handle(Geom_Curve) cs = BRep_Tool::Curve(eSeam, fs, ls);
        Handle(Geom2d_Line) pcS0 = new Geom2d_Line(gp_Pnt2d(f, r.vMin), gp_Dir2d(0.0, 1.0));
        Handle(Geom2d_Line) pcS1 =
            new Geom2d_Line(gp_Pnt2d(f + 2.0 * kPi, r.vMin), gp_Dir2d(0.0, 1.0));
        B.UpdateEdge(eSeam, pcS0, pcS1, surf, TopLoc_Location(), sewTol);
        if (!cs.IsNull()) B.Range(eSeam, fs, ls);
        eSeam.Closed(Standard_False);
    };

    auto finishFace = [&](TopoDS_Face got, bool publishSimple, int ciL, int ciH,
                          const TopoDS_Edge& eL, const TopoDS_Edge& eH) -> bool {
        if (got.IsNull()) return false;
        setFaceOutward(got, r.outwardNormal);
        addPcurvesOnFace(got, sewTol, true);
        if (!faceIsValid(got)) {
            try {
                ShapeFix_Face sff(got);
                sff.FixMissingSeamMode() = 1;
                sff.FixAddNaturalBoundMode() = 0;
                sff.FixOrientationMode() = 1;
                sff.Perform();
                TopoDS_Shape res = sff.Result();
                if (res.IsNull()) res = sff.Face();
                int nF = 0;
                TopoDS_Face g2;
                for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                    nF++;
                    g2 = TopoDS::Face(fx.Current());
                }
                if (nF == 1) got = g2;
            } catch (const Standard_Failure&) {
            }
            setFaceOutward(got, r.outwardNormal);
            addPcurvesOnFace(got, sewTol, true);
        }
        if (!faceIsValid(got) && !ensureFaceValid(got, meshTolCap(mv, &r))) {
            emit(warn, "seamed360: BRepCheck invalid on seamed face");
            return false;
        }
        if (publishSimple && ciL >= 0 && ciH >= 0 && (size_t)ciL < geom.size() &&
            (size_t)ciH < geom.size()) {
            geom[(size_t)ciL].collapsed = true;
            geom[(size_t)ciL].edges = {eL};
            collapsed[(size_t)ciL] = 1;
            geom[(size_t)ciH].collapsed = true;
            geom[(size_t)ciH].edges = {eH};
            collapsed[(size_t)ciH] = 1;
        }
        outF = got;
        return true;
    };

    TopoDS_Edge eSeam;
    try {
        gp_Lin lin(pntOf(mv, vL), r.ax.Direction());
        AnalyticCurve acs;
        acs.kind = AnalyticCurve::Lin;
        acs.lin = lin;
        snapVertexToCurve(verts[(size_t)vL], acs, snapCap);
        snapVertexToCurve(verts[(size_t)vH], acs, snapCap);
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

    // Simple path: each cap is one collapsed (or constructed) full 2π circle.
    TopoDS_Edge eL, eH;
    const int ciL0 = capL->chainIdx.front();
    const int ciH0 = capH->chainIdx.front();
    const bool simple =
        capL->chainIdx.size() == 1 && capH->chainIdx.size() == 1 &&
        takeFullCap(ciL0, circL, verts[(size_t)vL], eL) &&
        takeFullCap(ciH0, circH, verts[(size_t)vH], eH);

    try {
        if (simple) {
            const double u0 = azimuthOf(r, BRep_Tool::Pnt(verts[(size_t)vL]));
            bindIsoPCurves(eL, r.vMin, u0, eSeam, true);
            bindIsoPCurves(eH, r.vMax, u0, eSeam, false);
            BRep_Builder B;
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
            }
            return finishFace(got, true, ciL0, ciH0, eL, eH);
        }

        // Live I8 splits a cap into N chains (S04 boss-top torus junction). Build
        // the actual cap wires and insert the seam at u≈0; do not invent a
        // single circle that would steal neighbour TShapes.
        TopoDS_Wire wL, wH;
        if (!buildLoopWire(wL, *capL, rs, mv, geom, collapsed, meshE, edgeOk) ||
            !buildLoopWire(wH, *capH, rs, mv, geom, collapsed, meshE, edgeOk)) {
            emit(warn, "seamed360: composite cap wire failed — try TwoHalves");
            return false;
        }
        std::vector<TopoDS_Edge> pathL, pathH;
        if (!collectWireEdges(wH, pathH) || !collectWireEdges(wL, pathL) ||
            !rotateEdgesToVertex(pathH, verts[(size_t)vH]) ||
            !rotateEdgesToVertex(pathL, verts[(size_t)vL])) {
            emit(warn, "seamed360: cap wire does not pass seam vertex — try TwoHalves");
            return false;
        }
        BRep_Builder Bw;
        TopoDS_Wire ow;
        Bw.MakeWire(ow);
        Bw.Add(ow, eSeam);
        for (const auto& e : pathH) Bw.Add(ow, e);
        Bw.Add(ow, TopoDS::Edge(eSeam.Reversed()));
        for (int i = (int)pathL.size() - 1; i >= 0; --i)
            Bw.Add(ow, TopoDS::Edge(pathL[(size_t)i].Reversed()));
        ow.Closed(Standard_True);
        bindCylPCurves(ow, surf, r, sewTol);
        TopoDS_Face got;
        if (!makeFaceKeep(surf, ow, inners, r.outwardNormal, got)) {
            BRepBuilderAPI_MakeFace mf(surf, ow, Standard_True);
            if (!mf.IsDone()) {
                emit(warn, "seamed360: composite MakeFace not done — try TwoHalves");
                return false;
            }
            for (const auto& iw : inners) mf.Add(iw);
            if (!mf.IsDone()) return false;
            got = mf.Face();
        }
        return finishFace(got, false, -1, -1, eL, eH);
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
            Handle(Geom_CylindricalSurface) surf =
                Handle(Geom_CylindricalSurface)::DownCast(regionSurf(r, SurfVar::CylBase));
            bindCylPCurves(ow, surf, r, sewTol);
            BRepBuilderAPI_MakeFace mf(surf, ow, Standard_True);
            if (!mf.IsDone()) return false;
            f = mf.Face();
            setFaceOutward(f, r.outwardNormal);
            addPcurvesOnFace(f, sewTol, false);
            emitFaceSurfDiag(r.id, f);
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
        ensureFaceValid(f, meshTolCap(mv, &r));
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

// --- cascade ladder (D1 §1). Shared helper callable from A/B/C; wave-2 uses B. ---

enum class CascadeRung { None, U0, U1, U2 };

struct CascadeState {
    int u0Rounds = 0;
    bool u1Done = false;
    bool u2Done = false;
};

struct CascadeHit {
    int rid = -1;
    const char* src = "face";  // face | shell | resid
    bool faceValid = true;
    double residErr = 0;
};

struct CascadePlan {
    CascadeRung rung = CascadeRung::None;
    std::vector<int> explode;
    bool hostR2 = false;
};

bool brepStatusBad(const Handle(BRepCheck_Result)& res) {
    if (res.IsNull()) return false;
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next())
        if (it.Value() != BRepCheck_NoError) return true;
    return false;
}

void addCascadeHit(std::vector<CascadeHit>& hits, int rid, const char* src, bool faceValid,
                   double residErr = 0) {
    if (rid < 0) return;
    for (CascadeHit& h : hits) {
        if (h.rid != rid) continue;
        if (residErr > h.residErr) h.residErr = residErr;
        if (h.faceValid && !faceValid) h.faceValid = false;
        if (std::strcmp(h.src, "face") != 0 && std::strcmp(src, "face") == 0) h.src = src;
        return;
    }
    hits.push_back({rid, src, faceValid, residErr});
}

int countAnalyticCylinders(const RegionSet& rs) {
    int n = 0;
    for (const Region& r : rs.regions)
        if (r.type == SurfType::Cylinder) ++n;
    return n;
}

int countExplodedCylinders(const RegionSet& rs, const std::vector<char>& exploded) {
    int n = 0;
    for (const Region& r : rs.regions)
        if (r.type == SurfType::Cylinder && regionExploded(exploded, r.id)) ++n;
    return n;
}

int countIncidentChains(const RegionSet& rs, int rid) {
    int n = 0;
    for (const BoundaryChain& ch : rs.chains)
        if (chainTouches(ch, rid)) ++n;
    return n;
}

std::vector<int> regionNeighbours(const RegionSet& rs, int rid) {
    std::vector<int> out;
    for (const BoundaryChain& ch : rs.chains) {
        if (!chainTouches(ch, rid)) continue;
        const int o = otherReg(ch, rid);
        if (o < 0) continue;
        if (std::find(out.begin(), out.end(), o) == out.end()) out.push_back(o);
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool seamedClosed360(const Region* r) {
    return r && r->type == SurfType::Cylinder && r->closed360 &&
           (r->builtAs == BuiltAs::Seamed360 || r->builtAs == BuiltAs::TwoHalves);
}

std::vector<int> cylinderOneHop(const RegionSet& rs, int rid, const std::vector<char>& exploded) {
    std::vector<int> out;
    for (const BoundaryChain& ch : rs.chains) {
        if (!chainTouches(ch, rid)) continue;
        const int o = otherReg(ch, rid);
        const Region* rr = regionById(rs, o);
        if (!rr || rr->type != SurfType::Cylinder) continue;
        if (regionExploded(exploded, o)) continue;
        if (rr->closed360) continue;  // RULE 3.2 / U2: closed360 never U1 collateral
        if (std::find(out.begin(), out.end(), o) == out.end()) out.push_back(o);
    }
    std::sort(out.begin(), out.end());
    return out;
}

void collectFaceCulprits(const std::vector<TopoDS_Face>& built, const std::vector<int>& builtRid,
                         const std::vector<char>& exploded, std::vector<CascadeHit>& hits) {
    for (size_t i = 0; i < built.size(); i++) {
        const int id = (i < builtRid.size()) ? builtRid[i] : -1;
        if (id < 0 || regionExploded(exploded, id)) continue;
        const bool ok = faceIsValid(built[i]);
        if (!ok) addCascadeHit(hits, id, "face", false);
    }
}

void collectShellCulprits(const TopoDS_Shape& sh, const std::vector<TopoDS_Face>& built,
                          const std::vector<int>& builtRid, const RegionSet& rs,
                          const std::vector<char>& exploded, std::vector<CascadeHit>& hits) {
    if (sh.IsNull()) return;
    try {
        BRepCheck_Analyzer an(sh, Standard_True);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        std::vector<char> faceOk((size_t)faceMap.Extent() + 1, 1);
        for (size_t i = 0; i < built.size(); i++) {
            const int idx = faceMap.FindIndex(built[i]);
            if (idx <= 0 || i >= builtRid.size()) continue;
            faceRid[(size_t)idx] = builtRid[i];
            faceOk[(size_t)idx] = faceIsValid(built[i]) ? 1 : 0;
        }
        auto markFace = [&](const TopoDS_Shape& fsh) {
            const int idx = faceMap.FindIndex(fsh);
            if (idx <= 0) return;
            const int rid = faceRid[(size_t)idx];
            if (rid < 0 || regionExploded(exploded, rid)) return;
            addCascadeHit(hits, rid, "shell", faceOk[(size_t)idx] != 0);
        };
        for (int i = 1; i <= faceMap.Extent(); i++) {
            if (brepStatusBad(an.Result(faceMap(i)))) markFace(faceMap(i));
        }
        TopTools_IndexedDataMapOfShapeListOfShape wanc, eanc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_WIRE, TopAbs_FACE, wanc);
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, eanc);
        for (int i = 1; i <= wanc.Extent(); i++) {
            if (!brepStatusBad(an.Result(wanc.FindKey(i)))) continue;
            for (TopTools_ListOfShape::Iterator it(wanc(i)); it.More(); it.Next())
                markFace(it.Value());
        }
        for (int i = 1; i <= eanc.Extent(); i++) {
            if (!brepStatusBad(an.Result(eanc.FindKey(i)))) continue;
            for (TopTools_ListOfShape::Iterator it(eanc(i)); it.More(); it.Next())
                markFace(it.Value());
        }
        (void)rs;
    } catch (const Standard_Failure&) {
    }
}

template <typename Fn>
void cascadeParallelFor(size_t n, Fn fn) {
    if (n == 0) return;
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    const size_t workers = std::min((size_t)hw, n);
    if (workers <= 1) {
        for (size_t i = 0; i < n; i++) fn(i);
        return;
    }
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (size_t w = 0; w < workers; w++) {
        pool.emplace_back([=]() {
            for (size_t i = w; i < n; i += workers) fn(i);
        });
    }
    for (auto& t : pool) t.join();
}

double meshTriChordVol(const MeshView& mv, int localTri) {
    if (localTri < 0 || (size_t)localTri >= mv.nTri || !mv.compTris || !mv.tris || !mv.pts)
        return 0.0;
    const int gt = mv.compTris[localTri];
    if (gt < 0) return 0.0;
    const int* t = mv.tris[gt];
    const gp_XYZ a = mv.pts[t[0]], b = mv.pts[t[1]], c = mv.pts[t[2]];
    return a.Dot(b.Crossed(c)) / 6.0;
}

double meshViewVolume(const MeshView& mv) {
    double v = 0;
    for (size_t k = 0; k < mv.nTri; k++) v += meshTriChordVol(mv, (int)k);
    return v;
}

double regionChordVol(const MeshView& mv, const Region& r) {
    double v = 0;
    for (int t : r.tris) v += meshTriChordVol(mv, t);
    return v;
}

double faceVolumeContribution(const TopoDS_Face& f) {
    if (f.IsNull()) return 0.0;
    try {
        BRepAdaptor_Surface sa(f, Standard_False);
        if (sa.GetType() == GeomAbs_Plane) {
            // VolumeProperties(OnlyClosed=false) returns ~0 on a lone plane.
            // Pyramid: (1/3) (c · n) Area, n follows face orientation.
            GProp_GProps sp;
            BRepGProp::SurfaceProperties(f, sp);
            const double area = sp.Mass();
            gp_Dir n = sa.Plane().Axis().Direction();
            if (f.Orientation() == TopAbs_REVERSED) n.Reverse();
            return sp.CentreOfMass().XYZ().Dot(gp_XYZ(n.X(), n.Y(), n.Z())) * area / 3.0;
        }
        GProp_GProps gp;
        BRepGProp::VolumeProperties(f, gp, Standard_False);
        return gp.Mass();
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

double shapeVolSafe(const TopoDS_Shape& s) {
    if (s.IsNull()) return 0.0;
    try {
        GProp_GProps gp;
        BRepGProp::VolumeProperties(s, gp);
        return gp.Mass();
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

bool regionShippedAnalytic(const Region& r, const std::vector<char>& exploded) {
    if (regionExploded(exploded, r.id)) return false;
    if (!isAnalytic(&r)) return false;
    return r.builtAs == BuiltAs::Single || r.builtAs == BuiltAs::Seamed360 ||
           r.builtAs == BuiltAs::TwoHalves;
}

void collectResidualCulprits(const MeshView& mv, const RegionSet& rs, const TopoDS_Shape& sh,
                             const std::vector<TopoDS_Face>& built, const std::vector<int>& builtRid,
                             const std::vector<char>& exploded, double meshVol,
                             std::vector<CascadeHit>& hits, double& allAbs, double& shippedAbs) {
    allAbs = 0;
    shippedAbs = 0;
    int nAnalytic = 0;
    for (const Region& r : rs.regions) {
        allAbs += std::fabs(r.dVolPredicted);
        if (regionShippedAnalytic(r, exploded)) {
            shippedAbs += std::fabs(r.dVolPredicted);
            ++nAnalytic;
        }
    }
            if (collapseDiagEnabled())
        std::fprintf(stderr, "DIAG_CASCADE budget all=%.6f shipped=%.6f mesh=%.6f\n", allAbs,
                     shippedAbs, meshVol);
    if (nAnalytic <= 0) return;
    const double meshAbs = std::fabs(meshVol);
    const double kR = 3.0;      // unchanged
    const double kArc = 1.05;   // D4 §3: max plate ratio 235.37/224.2 = 1.050, +5% margin
    struct ResidRow {
        int rid = -1;
        double vface = 0, vchord = 0, dvol = 0, sub = 0, err = 0;
        bool shipped = false, pass = true, planeSkip = false;
    };
    std::vector<ResidRow> rows(rs.regions.size());
    cascadeParallelFor(rs.regions.size(), [&](size_t i) {
        const Region& r = rs.regions[i];
        ResidRow& row = rows[i];
        row.rid = r.id;
        if (!regionShippedAnalytic(r, exploded)) return;
        row.shipped = true;
        double vface = 0;
        for (size_t fi = 0; fi < built.size(); fi++) {
            if ((fi < builtRid.size() ? builtRid[fi] : -1) != r.id) continue;
            vface += faceVolumeContribution(built[fi]);
        }
        const double vchord = regionChordVol(mv, r);
        const double resid = vface - vchord;
        const bool predictorBlind = (r.type == SurfType::Plane
                                     && std::fabs(r.dVolPredicted) < 1e-6);
        const double arcTerm = predictorBlind ? kArc * std::fabs(resid) : 0.0;
        const double sub = std::max({1e-4 * meshAbs / (double)nAnalytic,
                                     kR * std::fabs(r.dVolPredicted),
                                     arcTerm});
        row.vface = vface;
        row.vchord = vchord;
        row.dvol = r.dVolPredicted;
        row.sub = sub;
        row.err = std::fabs(resid - r.dVolPredicted);
        row.pass = row.err <= sub;
        row.planeSkip = (r.type == SurfType::Plane && std::fabs(r.dVolPredicted) < 1e-6);
        // BRepGProp returned ~0 while the mesh patch has real volume — not a
        // landmine, a measurement miss. Do not explode.
        if (std::fabs(vface) < 1e-3 && std::fabs(vchord) > 1.0) {
            row.pass = true;
            row.planeSkip = true;
        }
    });
    for (const ResidRow& row : rows) {
        if (!row.shipped) continue;
        if (collapseDiagEnabled())
            std::fprintf(stderr,
                         "DIAG_CASCADE resid rid=%d vface=%.6f vchord=%.6f dVolPred=%.6f "
                         "sub=%.6f pass=%d\n",
                         row.rid, row.vface, row.vchord, row.dvol, row.sub, row.pass ? 1 : 0);
        // Actuate only the +121% landmine class (|dVol| large). Fixture-scale
        // cylinders must not die on a Vface/Vchord convention miss.
        const Region* rr = regionById(rs, row.rid);
        if (rr && rr->closed360) continue;
        if (!row.pass && !row.planeSkip && std::fabs(row.dvol) > 100.0)
            addCascadeHit(hits, row.rid, "resid", true, row.err);
    }
    (void)sh;
}

void sortHitsWorstFirst(std::vector<CascadeHit>& hits) {
    std::sort(hits.begin(), hits.end(), [](const CascadeHit& a, const CascadeHit& b) {
        if (a.residErr != b.residErr) return a.residErr > b.residErr;
        return a.rid < b.rid;
    });
}

void uniqueSorted(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

bool wouldSaturate(const RegionSet& rs, const std::vector<char>& exploded,
                   const std::vector<int>& extra) {
    const int nCyl = countAnalyticCylinders(rs);
    if (nCyl <= 0) return false;
    int would = countExplodedCylinders(rs, exploded);
    for (int id : extra) {
        const Region* r = regionById(rs, id);
        if (!r || r->type != SurfType::Cylinder) continue;
        if (regionExploded(exploded, id)) continue;
        ++would;
    }
    // RULE 1.5: escalate to U2 when exploded set would exceed 50% (ties round up).
    const int cap = (nCyl + 1) / 2;
    return would > cap;
}

bool touchesClosed360Cyl(const RegionSet& rs, int rid) {
    for (const BoundaryChain& ch : rs.chains) {
        if (!chainTouches(ch, rid)) continue;
        const Region* o = regionById(rs, otherReg(ch, rid));
        if (!o || o->type != SurfType::Cylinder || !o->closed360) continue;
        if (o->builtAs == BuiltAs::Seamed360 || o->builtAs == BuiltAs::TwoHalves ||
            o->builtAs == BuiltAs::Single)
            return true;
    }
    return false;
}

bool isHubPlane(const RegionSet& rs, int rid) {
    const Region* r = regionById(rs, rid);
    if (!r || r->type == SurfType::Cylinder) return false;
    const int nCyl = countAnalyticCylinders(rs);
    if (nCyl <= 0) return false;
    int hop = 0;
    for (const BoundaryChain& ch : rs.chains) {
        if (!chainTouches(ch, rid)) continue;
        const Region* o = regionById(rs, otherReg(ch, rid));
        if (o && o->type == SurfType::Cylinder) ++hop;
    }
    // Plane 0 shares a chain with 15/16 cylinders — exploding it is the blanket.
    return hop * 2 > nCyl;
}

// U0 set: face-invalid cylinders first (the true unit), then non-hub planes,
// then one hub plane (smallest chain count) only if nothing else remains.
// Exploding every hub plate in one pass opens the shell (census 0, RULE 1.4).
std::vector<int> selectU0Explode(const std::vector<CascadeHit>& culprits, const RegionSet& rs,
                                const std::vector<char>& exploded) {
    std::vector<int> cyls, smallPln, hubs;
    std::vector<int> badPlanes;
    for (const CascadeHit& h : culprits) {
        if (h.rid < 0 || regionExploded(exploded, h.rid)) continue;
        const Region* r = regionById(rs, h.rid);
        if (!r) continue;
        if (r->type == SurfType::Cylinder) {
            if (r->closed360) continue;
            cyls.push_back(h.rid);
        }
        else if (touchesClosed360Cyl(rs, h.rid))
            continue;  // RULE 3.2/S03: exploding plates against seamed holes steals TShapes
        else if (isHubPlane(rs, h.rid))
            hubs.push_back(h.rid);
        else
            smallPln.push_back(h.rid);
    }
    // Join-suspect: cylinder 1-hop from >=2 *non-hub* invalid planes.
    // Hub plates (0/1) touch almost every cylinder; they must not recruit the
    // blanket. rid=11 touches invalid plates 2 and 9 (spike M4).
    if (!smallPln.empty()) {
        for (const Region& r : rs.regions) {
            if (r.type != SurfType::Cylinder || regionExploded(exploded, r.id)) continue;
            if (r.closed360) continue;
            int nBad = 0;
            for (const BoundaryChain& ch : rs.chains) {
                if (!chainTouches(ch, r.id)) continue;
                const int o = otherReg(ch, r.id);
                if (std::find(smallPln.begin(), smallPln.end(), o) != smallPln.end())
                    ++nBad;
            }
            if (nBad >= 2) cyls.push_back(r.id);
        }
    }
    uniqueSorted(cyls);
    uniqueSorted(smallPln);
    uniqueSorted(hubs);
    if (!cyls.empty()) {
        uniqueSorted(cyls);
        if (wouldSaturate(rs, exploded, cyls)) {
            std::sort(cyls.begin(), cyls.end(), [&](int a, int b) {
                const int ca = countIncidentChains(rs, a), cb = countIncidentChains(rs, b);
                if (ca != cb) return ca < cb;
                return a < b;
            });
            std::vector<int> keep;
            for (int id : cyls) {
                keep.push_back(id);
                if (wouldSaturate(rs, exploded, keep)) {
                    keep.pop_back();
                    break;
                }
            }
            return keep.empty() ? std::vector<int>{cyls.front()} : keep;
        }
        return cyls;
    }
    if (!smallPln.empty()) return smallPln;
    if (hubs.empty()) return {};
    std::sort(hubs.begin(), hubs.end(), [&](int a, int b) {
        const int ca = countIncidentChains(rs, a), cb = countIncidentChains(rs, b);
        if (ca != cb) return ca < cb;
        return a < b;
    });
    return {hubs.front()};
}

std::vector<int> u2BlanketRids(const RegionSet& rs, const std::vector<char>& exploded) {
    std::vector<int> out;
    for (const Region& rg : rs.regions) {
        if (regionExploded(exploded, rg.id)) continue;
        if (rg.type != SurfType::Cylinder || rg.closed360) continue;
        out.push_back(rg.id);
    }
    return out;
}

// Shared helper — A/B/C may call; wave-2 site B is the only caller.
CascadePlan cascadeLadderPlan(CascadeState& st, const std::vector<CascadeHit>& culprits,
                              const RegionSet& rs, const std::vector<char>& exploded,
                              bool residActuator) {
    CascadePlan plan;
    if (st.u2Done) {
        plan.hostR2 = true;
        return plan;
    }
    std::vector<int> u0;
    if (residActuator) {
        for (const CascadeHit& h : culprits) {
            if (h.rid < 0 || regionExploded(exploded, h.rid)) continue;
            u0.push_back(h.rid);
            break;  // RULE 2.3: explode-worst, one per pass
        }
    } else {
        u0 = selectU0Explode(culprits, rs, exploded);
    }

    const bool satU0 = wouldSaturate(rs, exploded, u0);
    if (!satU0 && st.u0Rounds < 2 && !u0.empty()) {
        plan.rung = CascadeRung::U0;
        plan.explode = std::move(u0);
        st.u0Rounds++;
        return plan;
    }

    std::vector<int> u1;
    if (!st.u1Done) {
        for (const CascadeHit& h : culprits) {
            if (h.rid < 0 || regionExploded(exploded, h.rid)) continue;
            const auto nbrs = cylinderOneHop(rs, h.rid, exploded);
            u1.insert(u1.end(), nbrs.begin(), nbrs.end());
        }
        uniqueSorted(u1);
    }
    const bool satU1 = st.u1Done || wouldSaturate(rs, exploded, u1);
    if (!satU1 && (st.u0Rounds >= 2 || u0.empty()) && !st.u1Done) {
        st.u1Done = true;
        if (!u1.empty()) {
            plan.rung = CascadeRung::U1;
            plan.explode = std::move(u1);
            return plan;
        }
        // empty U1 → fall through to U2 in the same decide
    } else if (!st.u1Done && satU1) {
        st.u1Done = true;
    }

    st.u2Done = true;
    plan.rung = CascadeRung::U2;
    plan.explode = u2BlanketRids(rs, exploded);
    if (plan.explode.empty()) plan.hostR2 = true;
    return plan;
}

const char* cascadeRungName(CascadeRung r) {
    switch (r) {
        case CascadeRung::U0: return "U0";
        case CascadeRung::U1: return "U1";
        case CascadeRung::U2: return "U2";
        default: return "-";
    }
}

void diagCascadeCulprits(const std::vector<CascadeHit>& hits) {
    if (!collapseDiagEnabled()) return;
    for (const CascadeHit& h : hits)
        std::fprintf(stderr, "DIAG_CASCADE culprit rid=%d src=%s faceValid=%d\n", h.rid, h.src,
                     h.faceValid ? 1 : 0);
}

void diagCascadeExplode(int rid, CascadeRung rung, const RegionSet& rs) {
    if (!collapseDiagEnabled()) return;
    const auto nbrs = regionNeighbours(rs, rid);
    std::fprintf(stderr, "DIAG_CASCADE explode rid=%d rung=%s chains=%d nbrs=[", rid,
                 cascadeRungName(rung), countIncidentChains(rs, rid));
    for (size_t i = 0; i < nbrs.size(); i++) {
        if (i) std::fputc(',', stderr);
        std::fprintf(stderr, "%d", nbrs[i]);
    }
    std::fprintf(stderr, "]\n");
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
        double cap = Precision::Confusion();
        for (TopExp_Explorer vx0(f, TopAbs_VERTEX); vx0.More(); vx0.Next()) {
            gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(vx0.Current()));
            cap = std::max(cap, p.Distance(loc) * 0.01);
        }
        cap = std::max(cap, 1e-3);
        for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
            const TopoDS_Vertex& v = TopoDS::Vertex(vx.Current());
            double d = radDev(BRep_Tool::Pnt(v));
            if (d > cap) continue;
            if (d > BRep_Tool::Tolerance(v))
                B.UpdateVertex(v, std::min(d * 1.001 + Precision::Confusion(), cap));
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
            if (d > cap) continue;
            if (d > BRep_Tool::Tolerance(e))
                B.UpdateEdge(e, std::min(d * 1.001 + Precision::Confusion(), cap));
        }
    }
}

bool tryStageP(const MeshView& mv, RegionSet& rs, std::vector<TopoDS_Face>& out);

bool buildFaces(const MeshView& mv, RegionSet& rs, const std::vector<TopoDS_Vertex>& verts,
                std::vector<TopoDS_Face>& out, WarnFn warn) {
    out.clear();
    try {
        gCompNTri = (int)mv.nTri;
        if (mv.nTri == 0) return false;
        if (verts.size() < mv.nVtx) return false;
        if (!regionSetConsistent(mv, rs, verts)) {
            for (Region& r : rs.regions) {
                r.reject = Reject::ChainUnstable;
                r.builtAs = BuiltAs::NotBuilt;
            }
            return false;
        }

        installEPrimeDemotion(mv, rs);

        // Stage P routing (RULE 5.1a): prismatic -> profiles/prisms/union;
        // decline or build failure falls through byte-identically.
        if (tryStageP(mv, rs, out)) return true;

        const double sewTol = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
        const bool wasClosed = meshComponentClosed(mv);

        struct VSnap {
            gp_Pnt p;
            double t = 0;
        };
        std::vector<VSnap> vsnap(verts.size());
        for (size_t i = 0; i < verts.size(); i++) {
            if (verts[i].IsNull()) continue;
            vsnap[i].p = BRep_Tool::Pnt(verts[i]);
            vsnap[i].t = BRep_Tool::Tolerance(verts[i]);
        }

        int maxId = -1;
        for (const Region& r : rs.regions) maxId = std::max(maxId, r.id);
        std::vector<char> exploded((size_t)std::max(0, maxId) + 1, 0);

        std::vector<TopoDS_Edge> meshE(mv.nEdge);
        std::vector<char> edgeOk(mv.nEdge, 0);
        auto rebuildMeshEdges = [&]() {
            meshE.assign(mv.nEdge, TopoDS_Edge());
            edgeOk.assign(mv.nEdge, 0);
            if (diagCoverEnabled()) gDiagMeshEGen++;
            for (size_t i = 0; i < mv.nEdge; i++) {
                int a = mv.compEdges[i].first, b = mv.compEdges[i].second;
                if (a < 0 || b < 0 || (size_t)a >= verts.size() || (size_t)b >= verts.size())
                    continue;
                try {
                    BRepBuilderAPI_MakeEdge me(verts[(size_t)a], verts[(size_t)b]);
                    if (me.IsDone()) {
                        meshE[i] = me.Edge();
                        edgeOk[i] = 1;
                        if (diagCoverEnabled()) {
                            diagNoteEdge(meshE[i], "rebuildMeshEdges");
                            gDiagEdgeGen[diagTShapePtr(meshE[i])] = gDiagMeshEGen;
                            gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)a]);
                            gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)b]);
                            if (pa.Distance(pb) <= 1e-4 || a == b)
                                diagMaybeDegenEmit(meshE[i], "rebuildMeshEdges", -1, a, b,
                                                   (int)i);
                        }
                    }
                } catch (const Standard_Failure&) {
                }
            }
        };
        rebuildMeshEdges();
        auto restoreShared = [&]() {
            // OCCT UpdateVertex will not decrease tolerance. Recreate the
            // slot TShapes so R1/R2 do not inherit fat tols from a discarded
            // analytic attempt (adjudication F6). Assign a fresh handle —
            // in-place MakeVertex leaves snapped TShapes alive when leftover
            // faces still reference them, which dropped R2 triangles (Body9
            // volΔ 0.030937 → 0.030954).
            if (diagCoverEnabled())
                std::fprintf(stderr, "DIAG_RESTORE meshEGen-before=%d\n", gDiagMeshEGen);
            BRep_Builder Bb;
            for (size_t i = 0; i < verts.size(); i++) {
                if (verts[i].IsNull() && vsnap[i].t <= 0) continue;
                try {
                    TopoDS_Vertex nv;
                    Bb.MakeVertex(nv, vsnap[i].p,
                                  vsnap[i].t > 0 ? vsnap[i].t : Precision::Confusion());
                    const_cast<TopoDS_Vertex&>(verts[i]) = nv;
                } catch (const Standard_Failure&) {
                }
            }
            rebuildMeshEdges();
        };

        std::vector<char> collapsed(rs.chains.size(), 0);
        std::vector<ChainGeom> geom(rs.chains.size());
        // Per-chain edge-failure ledger (DECISION §5). Producer: chain build.
        // Consumer: J6 recoverPass 0 only — first-pass R1 enrollment on these
        // ids cascades to zero cylinders on Body11. IA_CYLCYL_NOGEOM (cyl|cyl
        // IntAna none) is legal polyline and is never enrolled.
        std::vector<char> chainEdgeFail(rs.chains.size(), 0);
        std::vector<char> fallbackBanned(rs.chains.size(), 0);
        std::vector<char> collapseBanned(rs.chains.size(), 0);
        std::vector<char> fallbackUsed(rs.chains.size(), 0);
        int recoverPass = 0;
        int rounds = 0;
        int fallbackGuardPass = 0;
        int j6UncollapsePass = 0;
        CascadeState cascadeSt;
        emitFailRidWarningOnce(warn);
        diagCascadeInject();

        auto rebuildCollapsed = [&]() {
            chainEdgeFail.assign(rs.chains.size(), 0);
            fallbackUsed.assign(rs.chains.size(), 0);
            int diagCollapse = 0;
            if (const char* v = std::getenv("STL2STEP_COLLAPSE_DIAG"); v && v[0] && v[0] != '0')
                diagCollapse = 1;
            int nMix = 0, nNone = 0, nFail = 0, nOk = 0;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                collapsed[ci] = 0;
                geom[ci] = ChainGeom{};
                const BoundaryChain& ch = rs.chains[ci];
                const Region* A = regionById(rs, ch.regA);
                const Region* B = regionById(rs, ch.regB);
                if (A && regionExploded(exploded, A->id)) A = nullptr;
                if (B && regionExploded(exploded, B->id)) B = nullptr;
                if (!isAnalytic(A) || !isAnalytic(B)) {
                    nMix++;
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
                static int lastCanonClearPass = -1;
                if (recoverPass != lastCanonClearPass) {
                    gRegionSurf.clear();
                    gSurfIdent.clear();
                    gSeamTShapes.clear();
                    gMeshTShapes.clear();
                    for (const auto& me : meshE) {
                        if (!me.IsNull()) gMeshTShapes.insert(diagTShapePtr(me));
                    }
                    lastCanonClearPass = recoverPass;
                }
                const bool bothCyl =
                    A->type == SurfType::Cylinder && B->type == SurfType::Cylinder;
                // Banned collapse: keep mesh polyline. No IntAna warning (AC #2)
                // and no chainEdgeFail enrollment (Body11).
                if (collapseBanned[ci] || (bothCyl && fallbackBanned[ci])) {
                    nNone++;
                    continue;
                }
                AnalyticCurve curve;
                bool usedFallback = false;
                WarnFn chainWarn =
                    (recoverPass == 0 && rounds == 0 && fallbackGuardPass == 0 &&
                     j6UncollapsePass == 0)
                        ? warn
                        : nullptr;
                if (ch.tangent) {
                    if (A->type == SurfType::Cylinder && B->type == SurfType::Plane)
                        curve = constructedGenerator(*A, *B);
                    else if (B->type == SurfType::Cylinder && A->type == SurfType::Plane)
                        curve = constructedGenerator(*B, *A);
                    else if (bothCyl) {
                        curve = bestCylCylConstructed(
                            *A, *B, mv, ch, intAnaAcceptResidual(mv, ch, sewTol, A, B));
                        usedFallback = (curve.kind != AnalyticCurve::None);
                    } else
                        curve = intersectSurfaces(*A, *B, mv, ch, sewTol, chainWarn);
                } else {
                    curve = intersectSurfaces(*A, *B, mv, ch, sewTol, chainWarn, &usedFallback);
                }
                if (usedFallback && curve.kind != AnalyticCurve::None) {
                    const double epsPl = derivedEpsPlane(mv);
                    const char* cls = cylCylClassName(*A, *B, epsPl);
                    if (!bothRegionsReferenceChain(rs, ch, (int)ci)) {
                        emitDiagFallback((int)ci, ch, cls, "reject-loop", 0, 0, 0.0,
                                         ch.meshVerts.size());
                        curve = {};
                        usedFallback = false;
                    } else if (curve.kind == AnalyticCurve::Lin &&
                               !constructedLinOnBothCylinders(
                                   curve, *A, *B, cylCylOnSurfaceTol(*A, *B, mv, sewTol))) {
                        emitDiagFallback((int)ci, ch, cls, "off-surface", 0, 0,
                                         chainResidual(curve, mv, ch), ch.meshVerts.size());
                    }
                }
                fallbackUsed[ci] = usedFallback ? 1 : 0;
                // IntAna none: cyl|cyl is legal polyline (IA_CYLCYL_NOGEOM).
                // plane|plane and plane|cyl misses are post-fit edge failures.
                if (curve.kind == AnalyticCurve::None) {
                    nNone++;
                    if ((A->type == SurfType::Plane && B->type == SurfType::Plane) ||
                        (A->type == SurfType::Plane && B->type == SurfType::Cylinder) ||
                        (A->type == SurfType::Cylinder && B->type == SurfType::Plane))
                        chainEdgeFail[ci] = 1;
                    continue;
                }

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
                const Region* cylR = (A->type == SurfType::Cylinder) ? A
                                   : (B->type == SurfType::Cylinder) ? B
                                                                     : nullptr;
                const Region* plnR = (A->type == SurfType::Plane) ? A
                                   : (B->type == SurfType::Plane) ? B
                                                                 : nullptr;
                // F1: closed360 cap circles must be V-isos of the fitted cylinder
                // (axis + R + XDirection). IntAna's own radius/location is ~1e-7
                // off the fitted R on live N>=12 prisms; that desyncs the 3d
                // circle from the cylindrical surface and BRepCheck-fails the
                // seamed face. Project onto the axis and use Region::radius.
                if (curve.kind == AnalyticCurve::Circ) {
                    if (cylR) {
                        double v = gp_Vec(cylR->ax.Location(), curve.circ.Location())
                                       .Dot(cylR->ax.Direction());
                        if (plnR && planePerpCylinder(*plnR, *cylR))
                            v = planeVOnCylinder(*plnR, *cylR);
                        curve.circ = cylinderIsoCircle(*cylR, v);
                    }
                }
                double snapCap = analyticSnapCap(mv, A, B);
                // ME_PLNPLN_*_PROJ: terminal-to-line residual (median 0.08 mm,
                // p90 0.45) exceeds sewTol. Bump to residual up to the
                // pickIntAna acceptance (1 mm) on plane|plane only — cyl
                // snaps stay on the I4 floor so extra MakeFace fails do not
                // eat the 127-cylinder small component.
                if (A->type == SurfType::Plane && B->type == SurfType::Plane) {
                    const double dA = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ia]));
                    const double dB = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ib]));
                    const double sew = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
                    const double acceptCap = std::max(sew * 50.0, 1.0);
                    snapCap = std::max(snapCap,
                                       std::min(std::max(std::isfinite(dA) ? dA : 0.0,
                                                         std::isfinite(dB) ? dB : 0.0),
                                                acceptCap));
                }
                if (A->type == SurfType::Cylinder && B->type == SurfType::Cylinder &&
                    curve.kind == AnalyticCurve::Lin) {
                    const double dA = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ia]));
                    const double dB = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ib]));
                    const double sew = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
                    const double acceptCap = intAnaAcceptResidual(mv, ch, sewTol, A, B);
                    snapCap = std::max(snapCap,
                                       std::min(std::max(std::isfinite(dA) ? dA : 0.0,
                                                         std::isfinite(dB) ? dB : 0.0),
                                                acceptCap));
                }
                snapVertexToCurve(verts[(size_t)ia], curve, snapCap);
                snapVertexToCurve(verts[(size_t)ib], curve, snapCap);
                TopoDS_Edge e;
                const Region* cylHint = (A && A->type == SurfType::Cylinder) ? A
                                       : (B && B->type == SurfType::Cylinder) ? B
                                                                             : nullptr;
                const Region* plnHint = (A && A->type == SurfType::Plane) ? A
                                       : (B && B->type == SurfType::Plane) ? B
                                                                         : nullptr;
                if (full && curve.kind == AnalyticCurve::Circ) {
                    // Plane inner wires around Seamed360 holes: MakeEdge(circ,V,V)
                    // fails on coarse loops — use the F1 iso-circle + makeFullCircle
                    // ladder (same as cap rings).
                    e = makeFullCircle(curve.circ, verts[(size_t)ia]);
                    if (e.IsNull())
                        e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
                } else {
                    const bool partialCylArc =
                        cylHint && !cylHint->closed360 && curve.kind == AnalyticCurve::Circ;
                    if (partialCylArc || (!full && curve.kind == AnalyticCurve::Circ &&
                                          ch.meshVerts.size() >= 3)) {
                        gp_Pnt mid;
                        if (partialCylArc && ch.meshVerts.size() < 3)
                            mid = partialCylCapArcMid(*cylHint, plnHint, curve.circ);
                        else {
                            // Open circular chain (partial-cyl cap): 180° MakeEdge(V,V) is
                            // ambiguous — pick the mesh vertex whose azimuth is closest
                            // to the region's u-mid (S05 nSides=26 is a true 180°).
                            mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
                            if (cylHint) {
                                const double um = 0.5 * (cylHint->uMin + cylHint->uMax);
                                double best = 1e300;
                                for (int lv : ch.meshVerts) {
                                    double du =
                                        std::fabs(regionU(*cylHint, pntOf(mv, lv)) - um);
                                    if (du < best) {
                                        best = du;
                                        mid = pntOf(mv, lv);
                                    }
                                }
                            }
                        }
                        e = makeArc(curve.circ, verts[(size_t)ia], verts[(size_t)ib], mid);
                        if (e.IsNull())
                            e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib],
                                                  full);
                    } else if (!full && curve.kind == AnalyticCurve::Elips &&
                               ch.meshVerts.size() >= 2) {
                        gp_Pnt mid;
                        if (ch.meshVerts.size() >= 3)
                            mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
                        else {
                            gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)ia]);
                            gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)ib]);
                            mid = gp_Pnt(0.5 * (pa.X() + pb.X()), 0.5 * (pa.Y() + pb.Y()),
                                         0.5 * (pa.Z() + pb.Z()));
                        }
                        e = makeEllipseArc(curve.elips, verts[(size_t)ia], verts[(size_t)ib], mid);
                        if (e.IsNull())
                            e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib],
                                                  full);
                    } else {
                        e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
                    }
                }
                if (e.IsNull()) {
                    const bool identicLine =
                        curve.kind == AnalyticCurve::Lin &&
                        (ia == ib || BRep_Tool::Pnt(verts[(size_t)ia])
                                         .Distance(BRep_Tool::Pnt(verts[(size_t)ib])) <=
                                         Precision::Confusion());
                    // IDENTIC_POINTS: keep polyline, do not invent an edge, do
                    // not enroll edge-failure R1 (the chain's legal fallback).
                    if (!identicLine) chainEdgeFail[ci] = 1;
                    if (!identicLine) nFail++;
                    emit(chainWarn, "smooth: analytic MakeEdge failed — keeping mesh polyline");
                    continue;
                }
                if (full && curve.kind == AnalyticCurve::Circ && ch.meshVerts.size() >= 2) {
                    double t0 = ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[0]));
                    double t1 = ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[1]));
                    double dt = t1 - t0;
                    while (dt < 0.0) dt += 2.0 * kPi;
                    if (dt > kPi) e.Reverse();
                }
                if (curve.kind == AnalyticCurve::Lin) {
                    auto partialCyl = [](const Region* R) {
                        return R && R->type == SurfType::Cylinder && !R->closed360;
                    };
                    if (partialCyl(A) || partialCyl(B))
                        e = orientEdgeFromTo(e, verts[(size_t)ia]);
                }
                {
                    const Region* ownerR = nullptr;
                    const Region* consumerR = nullptr;
                    if (A->type == SurfType::Plane && B->type == SurfType::Cylinder) {
                        ownerR = A;
                        consumerR = B;
                    } else if (B->type == SurfType::Plane && A->type == SurfType::Cylinder) {
                        ownerR = B;
                        consumerR = A;
                    } else if (A->id <= B->id) {
                        ownerR = A;
                        consumerR = B;
                    } else {
                        ownerR = B;
                        consumerR = A;
                    }
                    const char* kindStr = "poly";
                    if (curve.kind == AnalyticCurve::Circ) kindStr = "circ";
                    else if (curve.kind == AnalyticCurve::Lin) kindStr = "lin";
                    else if (curve.kind == AnalyticCurve::Elips) kindStr = "elips";
                    SeamBindCounts acc;
                    bindAllVariants(e, *ownerR, kindStr, mv, (int)ci, acc);
                    bindAllVariants(e, *consumerR, kindStr, mv, (int)ci, acc);
                    bindAllVariants(e, *ownerR, kindStr, mv, (int)ci, acc);
                    if (diagP2Enabled()) {
                        std::fprintf(stderr,
                                     "DIAG_SEAMBIND ci=%d nWrite=%d nGuardHit=%d nDoubleWrite=%d\n",
                                     (int)ci, acc.nWrite, acc.nGuardHit, acc.nDoubleWrite);
                        auto surfTypeStr = [](const Region* R) -> const char* {
                            if (!R) return "facet";
                            if (R->type == SurfType::Plane) return "plane";
                            if (R->type == SurfType::Cylinder) return "cyl";
                            return "facet";
                        };
                        Standard_Real f = 0, l = 0;
                        (void)BRep_Tool::Curve(e, f, l);
                        const TopLoc_Location loc;
                        Handle(Geom_Surface) surfO = ownerR->type == SurfType::Plane
                            ? regionSurf(*ownerR, SurfVar::Plane)
                            : regionSurf(*ownerR, SurfVar::CylBase);
                        Handle(Geom_Surface) surfC = consumerR->type == SurfType::Plane
                            ? regionSurf(*consumerR, SurfVar::Plane)
                            : regionSurf(*consumerR, SurfVar::CylBase);
                        Standard_Boolean hasO = Standard_False, hasC = Standard_False;
                        Standard_Real pf = 0, pl = 0;
                        (void)BRep_Tool::CurveOnSurface(e, surfO, loc, pf, pl, &hasO);
                        (void)BRep_Tool::CurveOnSurface(e, surfC, loc, pf, pl, &hasC);
                        std::fprintf(stderr,
                                     "DIAG_SEAM ci=%d owner=%d consumer=%d ownerType=%s "
                                     "consumerType=%s kind=%s f=%.9g l=%.9g storedOwner=%d "
                                     "storedConsumer=%d\n",
                                     (int)ci, ownerR->id, consumerR->id, surfTypeStr(ownerR),
                                     surfTypeStr(consumerR), kindStr, f, l, hasO ? 1 : 0,
                                     hasC ? 1 : 0);
                    }
                }
                geom[ci].collapsed = true;
                geom[ci].edges = {e};
                collapsed[ci] = 1;
                nOk++;
                if (diagCoverEnabled()) {
                    const char* esite = "rebuildCollapsed";
                    if (curve.kind == AnalyticCurve::Circ) esite = "rebuildCollapsed:circ";
                    else if (curve.kind == AnalyticCurve::Lin) esite = "rebuildCollapsed:lin";
                    else if (curve.kind == AnalyticCurve::Elips) esite = "rebuildCollapsed:elips";
                    diagMaybeDegenEmit(e, esite, (int)ci, A ? A->id : -1, B ? B->id : -1, ia);
                }
            }
            if (diagCollapse)
                std::fprintf(stderr,
                             "DIAG_COLLAPSE mix=%d none=%d fail=%d ok=%d total=%zu recover=%d "
                             "rounds=%d\n",
                             nMix, nNone, nFail, nOk, rs.chains.size(), recoverPass, rounds);
        };

        auto buildOneRegion = [&](Region& r, std::vector<TopoDS_Face>& acc) -> bool {
            if (failRidHits(r.id)) {
                r.reject = Reject::FaceBuildFailed;
                r.builtAs = BuiltAs::NotBuilt;
                return false;
            }
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
                                 f360, &exploded) &&
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
                        if (!ensureFaceValid(hf, meshTolCap(mv, &r))) allOk = false;
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

        auto explodeRegion = [&](int id) {
            if (id < 0 || (size_t)id >= exploded.size() || exploded[(size_t)id]) return;
            exploded[(size_t)id] = 1;
            Region* rr = regionByIdMut(rs, id);
            if (rr) {
                rr->builtAs = BuiltAs::ExplodedToFacets;
                rr->reject = Reject::FaceBuildFailed;
            }
            uncollapseRegionChains(id, rs, collapsed, geom);
            restoreShared();
        };
        auto explodeNonSeamed = [&]() -> bool {
            bool any = false;
            for (const Region& r : rs.regions) {
                if (regionExploded(exploded, r.id)) continue;
                if (r.builtAs == BuiltAs::Seamed360) continue;
                explodeRegion(r.id);
                any = true;
            }
            return any;
        };
        auto explodeAll = [&]() -> bool {
            bool any = false;
            for (const Region& r : rs.regions) {
                if (regionExploded(exploded, r.id)) continue;
                explodeRegion(r.id);
                any = true;
            }
            return any;
        };

        bool unstable = false;
        std::vector<TopoDS_Face> built;
        std::vector<int> builtRid;
        TopoDS_Shell sh;
    try_rebuild:
        rounds = 0;
        unstable = false;
        built.clear();
        builtRid.clear();
        const int kMaxRounds = 2;
        while (true) {
            rebuildCollapsed();
            built.clear();
            builtRid.clear();
            bool anyFail = false;
            std::vector<int> failedIds;
            std::vector<size_t> buildOrder(rs.regions.size());
            for (size_t i = 0; i < buildOrder.size(); i++) buildOrder[i] = i;
            std::stable_sort(buildOrder.begin(), buildOrder.end(), [&](size_t a, size_t b) {
                const Region& ra = rs.regions[a], &rb = rs.regions[b];
                const bool a360 = ra.closed360 && ra.type == SurfType::Cylinder;
                const bool b360 = rb.closed360 && rb.type == SurfType::Cylinder;
                if (a360 != b360) return a360;
                return ra.id < rb.id;
            });
            for (size_t oi : buildOrder) {
                Region& r = rs.regions[oi];
                if (regionExploded(exploded, r.id)) continue;
                std::vector<TopoDS_Face> acc;
                if (!buildOneRegion(r, acc)) {
                    anyFail = true;
                    failedIds.push_back(r.id);
                    r.reject = Reject::FaceBuildFailed;
                    r.builtAs = BuiltAs::NotBuilt;
                } else {
                    for (auto& f : acc) {
                        built.push_back(f);
                        builtRid.push_back(r.id);
                    }
                }
            }
            // Edge-failure R1 (DECISION §5). A region whose analytic reconstruction
            // failed — MakeFace not done, or the face is BRepCheck-invalid after
            // the projected-param / snap-cap MakeEdge path — is a post-fit
            // failure. Explode it, uncollapse its chains, rebuild neighbours.
            // MakeEdge-fail that still yields a BRepCheck-valid mixed face is
            // the legal polyline fallback (IA_* empty/same, IDENTIC_POINTS)
            // and is not exploded. Cap 2 rounds, then R2.
            if (!anyFail) break;
            std::sort(failedIds.begin(), failedIds.end());
            failedIds.erase(std::unique(failedIds.begin(), failedIds.end()), failedIds.end());
            for (int id : failedIds) explodeRegion(id);
            rounds++;
            if (rounds > kMaxRounds) {
                // 2-round cap: leftover post-fit fails were just exploded.
                // Rebuild once more. Remaining fails ⇒ ChainUnstable / R2.
                rebuildCollapsed();
                built.clear();
                builtRid.clear();
                anyFail = false;
                failedIds.clear();
                for (size_t oi : buildOrder) {
                    Region& r = rs.regions[oi];
                    if (regionExploded(exploded, r.id)) continue;
                    std::vector<TopoDS_Face> acc;
                    if (!buildOneRegion(r, acc)) {
                        anyFail = true;
                        failedIds.push_back(r.id);
                        r.reject = Reject::FaceBuildFailed;
                        r.builtAs = BuiltAs::NotBuilt;
                    } else {
                        for (auto& f : acc) {
                            built.push_back(f);
                            builtRid.push_back(r.id);
                        }
                    }
                }
                if (anyFail) unstable = true;
                break;
            }
        }

        // Dual-face incidence guard (SPEC): a constructed cyl|cyl fallback
        // lands only when both adjacent analytic faces reference the same
        // TShape. Asymmetric incidence is an orphan — revert that chain to
        // the pre-353f0fe polyline and rebuild. Both-missing (MakeFace copy)
        // is logged, not rejected, so handle-lock TShape sharing is not
        // second-guessed. Local sew probe is not used.
        // AC-SEAMFIX: the same ledger covers collapsed chains that touch an
        // in-band arch-chain commit (grow filletNbrA mark). Those plane|cyl
        // / cyl|cyl edges are new classes r2 never guarded.
        if (!unstable && fallbackGuardPass < 2) {
            bool rejected = false;
            int nSewed = 0, nOrphan = 0, nCopied = 0;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                if (!collapsed[ci] || geom[ci].edges.empty()) continue;
                const BoundaryChain& ch = rs.chains[ci];
                const Region* A = regionById(rs, ch.regA);
                const Region* B = regionById(rs, ch.regB);
                // Commit band (Body18/20 = 1948): guard every collapsed
                // analytic edge, not only constructed cyl|cyl fallbacks.
                // Grow's in-band applies add IntAna classes r2 never saw;
                // those 24 extra collapses are the 264→288 / 18→24 spike.
                const bool commitBand = mv.nTri > 1200 && mv.nTri <= 2500;
                if (!fallbackUsed[ci] && !chainTouchesArchCommit(A, B) && !commitBand)
                    continue;
                const double epsPl = derivedEpsPlane(mv);
                const char* cls = "other";
                if (A && B && A->type == SurfType::Cylinder && B->type == SurfType::Cylinder)
                    cls = cylCylClassName(*A, *B, epsPl);
                else if (chainTouchesArchCommit(A, B))
                    cls = "arch-chain";
                bool aBuilt = false, bBuilt = false, aHas = false, bHas = false;
                for (size_t i = 0; i < built.size(); i++) {
                    const bool hit = faceHasEdgeTShape(built[i], geom[ci].edges[0]);
                    if (builtRid[i] == ch.regA) {
                        aBuilt = true;
                        if (hit) aHas = true;
                    }
                    if (builtRid[i] == ch.regB) {
                        bBuilt = true;
                        if (hit) bHas = true;
                    }
                }
                const double resid = 0.0;
                if (aBuilt && bBuilt && aHas && bHas) {
                    nSewed++;
                    emitDiagFallback((int)ci, ch, cls, "sewed", 1, 1, resid,
                                     ch.meshVerts.size());
                } else if (aBuilt && bBuilt && aHas != bHas) {
                    nOrphan++;
                    emitDiagFallback((int)ci, ch, cls, "orphan", aHas ? 1 : 0, bHas ? 1 : 0,
                                     resid, ch.meshVerts.size());
                    fallbackBanned[ci] = 1;
                    collapseBanned[ci] = 1;  // plane|cyl too — fallbackBanned alone is cyl|cyl
                    collapsed[ci] = 0;
                    geom[ci] = {};
                    rejected = true;
                } else if (aBuilt && bBuilt && !aHas && !bHas) {
                    nCopied++;
                    emitDiagFallback((int)ci, ch, cls, "copied", 0, 0, resid,
                                     ch.meshVerts.size());
                } else {
                    nOrphan++;
                    emitDiagFallback((int)ci, ch, cls, "orphan-unbuilt", aHas ? 1 : 0,
                                     bHas ? 1 : 0, resid, ch.meshVerts.size());
                    fallbackBanned[ci] = 1;
                    collapseBanned[ci] = 1;
                    collapsed[ci] = 0;
                    geom[ci] = {};
                    rejected = true;
                }
            }
            if (diagJ6Enabled())
                std::fprintf(stderr,
                             "DIAG_FALLBACK_SUM sewed=%d orphan=%d copied=%d "
                             "guardPass=%d rejected=%d\n",
                             nSewed, nOrphan, nCopied, fallbackGuardPass, (int)rejected);
            if (rejected) {
                restoreShared();
                fallbackGuardPass++;
                goto try_rebuild;
            }
        }

        if (unstable) {
            for (Region& r : rs.regions) {
                if (r.reject == Reject::FaceBuildFailed) r.reject = Reject::ChainUnstable;
            }
            restoreShared();
            explodeAll();
            built.clear();
            builtRid.clear();
            for (size_t k = 0; k < mv.nTri; k++) {
                TopoDS_Face f = makeFacet(mv, verts, meshE, edgeOk, k);
                if (!f.IsNull()) {
                    built.push_back(f);
                    builtRid.push_back(-1);
                }
            }
            if (built.empty()) {
                out.clear();
                return false;
            }
            // Fall through to shell assembly of the R2 faceted baseline.
        }

        // Facet islands + exploded-region triangles. R2 already emitted every
        // triangle above.
        if (!unstable) {
            for (size_t k = 0; k < mv.nTri; k++) {
                int rid = (k < rs.triRegion.size()) ? rs.triRegion[k] : -1;
                int iid = (k < rs.triIsland.size()) ? rs.triIsland[k] : -1;
                bool exp = (rid >= 0 && regionExploded(exploded, rid));
                if (iid < 0 && !exp) continue;
                TopoDS_Face f = makeFacet(mv, verts, meshE, edgeOk, k);
                if (!f.IsNull()) {
                    ensureFaceValid(f, meshTolCap(mv, nullptr));
                    built.push_back(f);
                    builtRid.push_back(exp ? rid : -1);
                }
            }
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
                    // copies them and is the R1 explode J2-twin source. Cylinders
                    // and inner-wired plates share collapsed cap TShapes; copying
                    // those opens the shell (live S03/S10).
                    if (nE == 3 && sa.GetType() == GeomAbs_Plane) continue;
                    if (sa.GetType() == GeomAbs_Cylinder) continue;
                    if (nE > 6) continue;
                }
                // F7: do not ShapeFix_Face here — it displaces shared verts[]
                // (measured 18.8 mm on S02). Invalid faces explode via R1.
            } catch (const Standard_Failure&) {
            }
        }

        {
            orientFaceWalk(built);
        }

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
        B.MakeShell(sh);
        for (auto& f : built) B.Add(sh, f);
        const bool allTriFacets = [&]() {
            if (built.empty()) return false;
            for (const auto& f : built) {
                int nE = 0;
                for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) nE++;
                if (nE != 3) return false;
            }
            return true;
        }();
        double spCap = meshTolCap(mv, nullptr);
        for (const Region& rg : rs.regions) spCap = std::max(spCap, meshTolCap(mv, &rg));
        spCap = std::max(spCap, 0.05 * (mv.diag > 0.0 ? mv.diag : 1.0));
        spCap = std::max(spCap, 25.0);  // same floor as live tolOk; blocks 10^3 mm poison
        if (mv.nTri < 10000) {
            BRep_Builder Bb;
            for (auto& F : built) {
                if (F.IsNull()) continue;
                TopLoc_Location loc;
                Handle(Geom_Surface) sF = BRep_Tool::Surface(F, loc);
                if (sF.IsNull()) continue;
                for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next()) {
                    TopoDS_Edge eW = TopoDS::Edge(ex.Current());
                    const void* ts = diagTShapePtr(eW);
                    if (!ts || !gSeamTShapes.count(ts)) continue;
                    Standard_Boolean hasPc = Standard_False;
                    Standard_Real pf = 0, pl = 0;
                    (void)BRep_Tool::CurveOnSurface(eW, sF, loc, pf, pl, &hasPc);
                    if (hasPc) continue;
                    Standard_Real a = 0, b = 0;
                    Handle(Geom_Curve) c3 = BRep_Tool::Curve(eW, a, b);
                    if (c3.IsNull() || b - a <= Precision::PConfusion()) continue;
                    Handle(Geom2d_Curve) c2d = makePCurveOnSurf(c3, a, b, sF, nullptr);
                    if (c2d.IsNull()) continue;
                    double tol = Precision::Confusion();
                    const double existing = BRep_Tool::Tolerance(eW);
                    if (existing > tol) tol = existing;
                    Bb.UpdateEdge(eW, c2d, sF, loc, tol);
                    Bb.Range(eW, sF, loc, a, b);
                }
            }
        }
        if (diagP2Enabled()) {
            const bool preClosed = BRep_Tool::IsClosed(sh) == Standard_True;
            int freeE = 0, nI = 0, nII = 0, nIII = 0;
            try {
                TopTools_IndexedDataMapOfShapeListOfShape anc;
                TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
                for (int i = 1; i <= anc.Extent(); i++)
                    if (anc(i).Extent() < 2) freeE++;
                for (const auto& F : built) {
                    if (F.IsNull()) continue;
                    TopLoc_Location loc;
                    Handle(Geom_Surface) s = BRep_Tool::Surface(F, loc);
                    const bool facePlane = (s && Handle(Geom_Plane)::DownCast(s));
                    const bool faceCyl = (s && basisCylOf(s));
                    if (!facePlane && !faceCyl) continue;
                    std::unordered_set<const void*> seen;
                    for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next()) {
                        const TopoDS_Edge e = TopoDS::Edge(ex.Current());
                        const void* ts = diagTShapePtr(e);
                        if (!ts || !seen.insert(ts).second) continue;
                        Standard_Real f = 0, l = 0;
                        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
                        if (c3.IsNull() || l - f <= Precision::PConfusion()) continue;
                        Standard_Boolean stored = Standard_False;
                        Standard_Real pf = 0, pl = 0;
                        if (!s.IsNull())
                            (void)BRep_Tool::CurveOnSurface(e, s, loc, pf, pl, &stored);
                        if (stored) continue;
                        if (gSeamTShapes.count(ts)) nI++;
                        else if (facePlane && gMeshTShapes.count(ts)) nII++;
                        else nIII++;
                    }
                }
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr,
                         "DIAG_PREJ4 closed=%d freeEdges=%d nOrphanI=%d nOrphanII=%d nOrphanIII=%d\n",
                         preClosed ? 1 : 0, freeE, nI, nII, nIII);
        }
        try {
            BRepLib::SameParameter(sh, std::min(sewTol, spCap), /*forced=*/Standard_True);
            bool needRelax = false;
            for (const Region& rg : rs.regions) {
                if (rg.type == SurfType::Cylinder && !regionExploded(exploded, rg.id))
                    needRelax = true;
            }
            if (needRelax && !allTriFacets) {
                BRep_Builder bb;
                bool any = false;
                for (TopExp_Explorer ex(sh, TopAbs_EDGE); ex.More(); ex.Next()) {
                    const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
                    if (BRep_Tool::SameParameter(e)) continue;
                    any = true;
                    const double fat =
                        std::min(std::max(BRep_Tool::Tolerance(e), sewTol) * 10.0 +
                                     Precision::Confusion(),
                                 spCap);
                    bb.UpdateEdge(e, fat);
                }
                if (any)
                    BRepLib::SameParameter(sh, std::min(sewTol * 10.0, spCap),
                                           /*forced=*/Standard_True);
            }
        } catch (const Standard_Failure&) {
            return false;
        }

        // J6
        bool shClosed = BRep_Tool::IsClosed(sh);
        sh.Closed(shClosed ? Standard_True : Standard_False);
        if (wasClosed && !shClosed) {
            int freeE = 0;
            TopTools_IndexedDataMapOfShapeListOfShape anc;
            TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
            for (int i = 1; i <= anc.Extent(); i++)
                if (anc(i).Extent() < 2) freeE++;
            if (diagCoverEnabled())
                diagCoverDump(mv, rs, exploded, builtRid, built, meshE, collapsed, geom, sh,
                              j6UncollapsePass);
            if (diagJ6Enabled())
                diagJ6FreeEdges(mv, sh, built, builtRid, rs, collapsed);
            emit(warn, "J6: shell not closed freeEdges=" + std::to_string(freeE) +
                           " faces=" + std::to_string(built.size()) +
                           " recover=" + std::to_string(recoverPass));
            // Targeted J6 heal: uncollapse analytic chains whose terminals
            // match a free edge. Those TShapes did not pair. Mesh polyline
            // shares meshE and is the pre-353f0fe path for that chain only.
            // Skip large multi-hole parts (Body11, 15300 tris): mixed
            // uncollapse + abort-recover zeroes the 127-cyl floor. The
            // 8-file J6 set is all < 8000 tris.
            if (j6UncollapsePass < 1 && freeE > 0 && mv.nTri < 10000) {
                const double matchTol =
                    std::max(std::max(0.5, sewTol * 20.0), derivedEpsPlane(mv) * 4.0);
                std::vector<int> heal;
                auto addHeal = [&](int ci) {
                    if (ci < 0 || (size_t)ci >= collapseBanned.size()) return;
                    if (std::find(heal.begin(), heal.end(), ci) != heal.end()) return;
                    heal.push_back(ci);
                };
                for (int i = 1; i <= anc.Extent(); i++) {
                    if (anc(i).Extent() >= 2) continue;
                    const TopoDS_Edge& e = TopoDS::Edge(anc.FindKey(i));
                    TopoDS_Vertex va, vb;
                    TopExp::Vertices(e, va, vb, Standard_True);
                    gp_Pnt pa = BRep_Tool::Pnt(va), pb = BRep_Tool::Pnt(vb);
                    addHeal(matchCollapsedChainToSegment(mv, rs, collapsed, pa, pb, matchTol));
                }
                if (!heal.empty()) {
                    for (int ci : heal) {
                        collapseBanned[(size_t)ci] = 1;
                        if ((size_t)ci < fallbackBanned.size()) fallbackBanned[(size_t)ci] = 1;
                        collapsed[(size_t)ci] = 0;
                        geom[(size_t)ci] = {};
                        if (diagJ6Enabled()) {
                            const BoundaryChain& ch = rs.chains[(size_t)ci];
                            const Region* A = regionById(rs, ch.regA);
                            const Region* B = regionById(rs, ch.regB);
                            const char* cls = (A && B && A->type == SurfType::Cylinder &&
                                               B->type == SurfType::Cylinder)
                                                  ? cylCylClassName(*A, *B, derivedEpsPlane(mv))
                                                  : "mixed";
                            emitDiagFallback(ci, ch, cls, "j6-uncollapse", 0, 0, 0.0,
                                             ch.meshVerts.size());
                        }
                    }
                    if (diagJ6Enabled())
                        std::fprintf(stderr, "DIAG_FALLBACK_SUM j6-uncollapse=%zu freeWas=%d\n",
                                     heal.size(), freeE);
                    restoreShared();
                    j6UncollapsePass++;
                    goto try_rebuild;
                }
            }
            // Revert-class (Body9/12/18/20): after the heal, do not explode
            // recover — that dirties vertex TShapes and worsens R2 volΔ.
            // Body11 never enters this branch (nTri >= 10000 skips the heal).
            if (j6UncollapsePass > 0) {
                restoreShared();
                out.clear();
                return false;
            }
            // Clean R1 to all-facets only when no Seamed360 survived. Rebuilding
            // seamed holes against exploded plates invents cap circles and
            // steals neighbour TShapes (S03/S16). That is AC4 escalation.
            bool keepCyl = false;
            for (const Region& rg : rs.regions) {
                if (rg.type == SurfType::Cylinder &&
                    (rg.builtAs == BuiltAs::Single || rg.builtAs == BuiltAs::Seamed360 ||
                     rg.builtAs == BuiltAs::TwoHalves) &&
                    !regionExploded(exploded, rg.id))
                    keepCyl = true;
            }
            if (recoverPass < 1) {
                bool did = false;
                // J6 deferred edge-failure R1: explode regions whose analytic
                // chains failed (ledger from rebuildCollapsed). Not enrolled on
                // the MakeFace R1 pass — that blanket cascade zeroed cylinders.
                {
                    std::vector<int> edgeFailIds;
                    auto addId = [&](int id) {
                        if (id < 0 || regionExploded(exploded, id)) return;
                        if (std::find(edgeFailIds.begin(), edgeFailIds.end(), id) ==
                            edgeFailIds.end())
                            edgeFailIds.push_back(id);
                    };
                    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                        if (!chainEdgeFail[ci]) continue;
                        addId(rs.chains[ci].regA);
                        addId(rs.chains[ci].regB);
                    }
                    std::sort(edgeFailIds.begin(), edgeFailIds.end());
                    for (int id : edgeFailIds) {
                        explodeRegion(id);
                        did = true;
                    }
                }
                if (!did) {
                    if (!keepCyl)
                        did = explodeAll();
                    else {
                        // Keep Seamed360 holes; drop junk partials (S03 cone
                        // fragments) that left the shell open.
                        bool any = false;
                        for (const Region& rg : rs.regions) {
                            if (regionExploded(exploded, rg.id)) continue;
                            if (rg.type != SurfType::Cylinder || rg.closed360) continue;
                            explodeRegion(rg.id);
                            any = true;
                        }
                        did = any;
                    }
                }
                if (did) {
                    recoverPass++;
                    goto try_rebuild;
                }
            } else if (recoverPass < 2 && explodeAll()) {
                recoverPass++;
                goto try_rebuild;
            }
            restoreShared();
            out.clear();
            return false;
        }

        auto applyCascadePlan = [&](const CascadePlan& plan) -> bool {
            bool any = false;
            for (int id : plan.explode) {
                if (id < 0 || regionExploded(exploded, id)) continue;
                diagCascadeExplode(id, plan.rung, rs);
                explodeRegion(id);
                any = true;
            }
            return any;
        };

        const bool shValid = shellIsValid(sh);
        // Site B is closed-but-invalid only. Open shells are site C (untouched).
        if (BRep_Tool::IsClosed(sh) && !shValid) {
            dumpShellCheck(sh, built, builtRid, rs);
            std::vector<CascadeHit> culprits;
            collectFaceCulprits(built, builtRid, exploded, culprits);
            collectShellCulprits(sh, built, builtRid, rs, exploded, culprits);
            diagCascadeCulprits(culprits);
            // RULE 1.4: empty culprit set + invalid shell ⇒ escalate, never ship.
            CascadePlan plan = cascadeLadderPlan(cascadeSt, culprits, rs, exploded, false);
            if (applyCascadePlan(plan)) {
                recoverPass++;
                goto try_rebuild;
            }
            if (plan.hostR2 || cascadeSt.u2Done || !shValid) {
                restoreShared();
                out.clear();
                return false;
            }
        } else if (BRep_Tool::IsClosed(sh) && shValid) {
            // RULE 2.2 / 2.3 / 2.4 — detector + explode-worst on a closed valid shell.
            const double meshVol = meshViewVolume(mv);
            std::vector<CascadeHit> residHits;
            double allAbs = 0, shippedAbs = 0;
            collectResidualCulprits(mv, rs, sh, built, builtRid, exploded, meshVol, residHits,
                                    allAbs, shippedAbs);
            const double shellVol = shapeVolSafe(sh);
            const double meshAbs = std::fabs(meshVol);
            const double tightBudget = std::max(1e-4 * meshAbs, 3.0 * shippedAbs);
            const bool tightFail = std::fabs(shellVol - meshVol) > tightBudget;
            if (tightFail && residHits.empty()) {
                // Aggregate-only miss: explode the landmine-class worst |dVol|.
                // Never pick a closed360 hole (S03/S16 / full_360_hole).
                int worst = -1;
                double worstAbs = -1;
                for (const Region& r : rs.regions) {
                    if (!regionShippedAnalytic(r, exploded)) continue;
                    if (r.closed360) continue;
                    const double a = std::fabs(r.dVolPredicted);
                    if (a <= 100.0) continue;
                    if (a > worstAbs) {
                        worstAbs = a;
                        worst = r.id;
                    }
                }
                if (worst >= 0) addCascadeHit(residHits, worst, "resid", true, worstAbs);
            }
            if (!residHits.empty()) {
                sortHitsWorstFirst(residHits);
                diagCascadeCulprits(residHits);
                CascadePlan plan = cascadeLadderPlan(cascadeSt, residHits, rs, exploded, true);
                if (applyCascadePlan(plan)) {
                    recoverPass++;
                    goto try_rebuild;
                }
                if (plan.hostR2 && tightFail && !cascadeSt.u2Done) {
                    restoreShared();
                    out.clear();
                    return false;
                }
                // After U2 the host D4.5 probe (RULE 2.0) is the volume authority.
                // Do not R2 a closed+valid U2 shell here — that is census 0.
            }
        }

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
