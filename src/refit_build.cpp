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
#include "refit_cone_bind.hpp"
#include "refit_cone_math.hpp"
#include "refit_internal.hpp"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_CurveRepresentation.hxx>
#include <BRep_ListOfCurveRepresentation.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepLib_ValidateEdge.hxx>
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
#include <BRepCheck_Face.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <BRep_CurveOnSurface.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
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
    const bool planeOrCyl = (r->type == SurfType::Plane || r->type == SurfType::Cylinder);
    const bool chamferCone =
        (r->type == SurfType::Cone && r->origin == Origin::ChamferCone);
    if (!planeOrCyl && !chamferCone) return false;
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
    if (maxId >= 0) {
        gEPrimeDemoted.assign((size_t)maxId + 1, 0);
    }
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
    // E′′ v2 measurement (default-off, STL2STEP_EPRIME_DIAG only; no product change).
    const double diagLen = mv.diag > 0.0 ? mv.diag : 1.0;
    const double epsMesh = std::max(std::max(mv.weldTol, 1e-4 * diagLen), 1e-3);
    const double epsPlane = std::max(std::max(epsMesh, mv.sewTol), 0.02);
    auto meshTolCapInline = [&](const Region& rr) {
        double c = std::max(mv.sewTol, Precision::Confusion());
        c = std::max(c, rr.chordSagitta);
        if (rr.radius > 0.0)
            c = std::max(c, rr.radius * (1.0 - std::cos(kPi / 4.0)));
        return c * 1.001 + Precision::Confusion();
    };
    auto axisSeparation = [](const Region& a, const Region& b) {
        const gp_Vec dloc(a.ax.Location(), b.ax.Location());
        return gp_Vec(a.ax.Direction()).Crossed(dloc).Magnitude();
    };
    auto fitResiduals = [&](const Region& rr, double& maxV, double& chordMid) {
        maxV = 0.0;
        chordMid = 0.0;
        if (!(rr.radius > 0.0)) return;
        gp_Cylinder cyl = asCyl(rr);
        auto radDev = [&](const gp_Pnt& p) {
            gp_Vec v(cyl.Location(), p);
            return std::fabs(gp_Vec(cyl.Axis().Direction()).Crossed(v).Magnitude() -
                             cyl.Radius());
        };
        for (int lt : rr.tris) {
            if (lt < 0 || (size_t)lt >= mv.nTri) continue;
            const int* T = mv.tris[mv.compTris[lt]];
            for (int k = 0; k < 3; k++) maxV = std::max(maxV, radDev(gp_Pnt(mv.pts[T[k]])));
        }
        for (const BoundaryChain& ch : rs.chains) {
            if (ch.regA != rr.id && ch.regB != rr.id) continue;
            for (int eid : ch.meshEdges) {
                if (eid < 0 || (size_t)eid >= mv.nEdge) continue;
                const auto& evv = mv.compEdges[eid];
                gp_Pnt a = pntOf(mv, evv.first), b = pntOf(mv, evv.second);
                chordMid = std::max(chordMid, radDev(gp_Pnt(0.5 * (a.XYZ() + b.XYZ()))));
            }
        }
    };
    auto circleFitSigmaR = [&](const Region& rr, double& sigmaR, double& Rfit) {
        sigmaR = 1e300;
        Rfit = rr.radius;
        if (!(rr.radius > 0.0)) return;
        gp_Cylinder cyl = asCyl(rr);
        const gp_Pnt origin = cyl.Location();
        const gp_Dir axD = cyl.Axis().Direction();
        std::vector<double> res;
        res.reserve(rr.tris.size() * 3);
        for (int lt : rr.tris) {
            if (lt < 0 || (size_t)lt >= mv.nTri) continue;
            const int* T = mv.tris[mv.compTris[lt]];
            for (int k = 0; k < 3; k++) {
                gp_Vec v(origin, gp_Pnt(mv.pts[T[k]]));
                res.push_back(gp_Vec(axD).Crossed(v).Magnitude() - rr.radius);
            }
        }
        const std::size_t np = res.size();
        if (np < 2) return;
        double mean = 0.0;
        for (double e : res) mean += e;
        mean /= static_cast<double>(np);
        double ss = 0.0;
        for (double e : res) {
            const double d = e - mean;
            ss += d * d;
        }
        const double sigmaE =
            std::sqrt(ss / std::max(1.0, static_cast<double>(np) - 1.0));
        const double spanU =
            rr.closed360 ? 2.0 * kPi : std::max(0.0, rr.uMax - rr.uMin);
        const double halfAng = std::min(0.5 * spanU, kPi);
        double arcDen = 1.0 - std::cos(halfAng);
        if (arcDen < 1e-9) arcDen = 1e-9;
        sigmaR = sigmaE / (arcDen * std::sqrt(static_cast<double>(np)));
    };
    auto findRegById = [&](int id) -> const Region* {
        if (id < 0) return nullptr;
        for (const Region& rr : rs.regions)
            if (rr.id == id) return &rr;
        return nullptr;
    };
    struct S4Row {
        double spanU = 0.0;
        double sigmaR = 0.0;
        double sigmaROverR = 0.0;
        double clusterScatter = 0.0;
        int clusterId = -1;
        int chainN = 0;
        int chainBestRid = -1;
        double chainBestR = -1.0;
        int chainBestNTri = -1;
        double chainBestAbsDR = -1.0;
        int confirmRid = -1;
        int confirmNTri = -1;
        double confirmR = -1.0;
        int sibRelRid = -1;
        double sibRelR = -1.0;
        int sibRelNTri = -1;
        double sibRelSep = -1.0;
        double sibRelAbsDR = -1.0;
        double sibRelAbsDROverCap = -1.0;
        int absorbParent = -1;
        const char* absorbFail = "related";
        double absorbCarryDev = -1.0;
        double absorbCarryLim = -1.0;
    };
    std::vector<S4Row> s4((size_t)n);
    {
        std::vector<int> clParent((size_t)n);
        for (int i = 0; i < n; ++i) clParent[(size_t)i] = i;
        auto clFind = [&](int x) {
            int r = x;
            while (clParent[(size_t)r] != r) r = clParent[(size_t)r];
            while (clParent[(size_t)x] != r) {
                const int nxt = clParent[(size_t)x];
                clParent[(size_t)x] = r;
                x = nxt;
            }
            return r;
        };
        auto clUnite = [&](int a, int b) {
            a = clFind(a);
            b = clFind(b);
            if (a == b) return;
            if (a < b) clParent[(size_t)b] = a;
            else clParent[(size_t)a] = b;
        };
        for (int i = 0; i < n; ++i) {
            const Region& ri = rs.regions[(size_t)cylIdx[i]];
            for (int j = i + 1; j < n; ++j) {
                const Region& rj = rs.regions[(size_t)cylIdx[j]];
                if (axisSin(ri.ax.Direction(), rj.ax.Direction()) > tauAx) continue;
                const double sep = axisSeparation(ri, rj);
                const double sepTol =
                    std::max(meshTolCapInline(ri), meshTolCapInline(rj));
                if (sep > sepTol) continue;
                clUnite(i, j);
            }
        }
        std::unordered_map<int, std::vector<int>> clusters;
        for (int i = 0; i < n; ++i) clusters[clFind(i)].push_back(i);
        for (int i = 0; i < n; ++i) {
            const Region& r = rs.regions[(size_t)cylIdx[i]];
            S4Row& row = s4[(size_t)i];
            row.spanU = r.closed360 ? 2.0 * kPi : std::max(0.0, r.uMax - r.uMin);
            double Rfit = r.radius;
            circleFitSigmaR(r, row.sigmaR, Rfit);
            row.sigmaROverR =
                r.radius > 0.0 ? row.sigmaR / r.radius : row.sigmaR;
            row.clusterId = clFind(i);
            const auto& memb = clusters[row.clusterId];
            double rMin = 1e300, rMax = -1e300, rSum = 0.0;
            for (int mi : memb) {
                const double Rm = rs.regions[(size_t)cylIdx[mi]].radius;
                rMin = std::min(rMin, Rm);
                rMax = std::max(rMax, Rm);
                rSum += Rm;
            }
            const double rMean = memb.empty() ? r.radius : rSum / memb.size();
            row.clusterScatter =
                rMean > 0.0 ? (rMax - rMin) / rMean : 0.0;
            for (const BoundaryChain& ch : rs.chains) {
                int other = -1;
                if (ch.regA == r.id) other = ch.regB;
                else if (ch.regB == r.id) other = ch.regA;
                if (other < 0) continue;
                const Region* orr = findRegById(other);
                if (!orr || orr->type != SurfType::Cylinder) continue;
                row.chainN++;
                const double dR = std::fabs(r.radius - orr->radius);
                if (row.chainBestRid < 0 || (int)orr->tris.size() > row.chainBestNTri ||
                    ((int)orr->tris.size() == row.chainBestNTri && dR < row.chainBestAbsDR)) {
                    row.chainBestRid = orr->id;
                    row.chainBestR = orr->radius;
                    row.chainBestNTri = (int)orr->tris.size();
                    row.chainBestAbsDR = dR;
                }
            }
            const double agreeTol = row.sigmaR;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                const Region& s = rs.regions[(size_t)cylIdx[j]];
                if ((int)s.tris.size() <= (int)r.tris.size()) continue;
                bool related = false;
                for (const BoundaryChain& ch : rs.chains) {
                    if ((ch.regA == r.id && ch.regB == s.id) ||
                        (ch.regB == r.id && ch.regA == s.id)) {
                        related = true;
                        break;
                    }
                }
                if (!related && axisSin(r.ax.Direction(), s.ax.Direction()) > tauAx)
                    continue;
                double sSig = 1e300, sRfit = s.radius;
                circleFitSigmaR(s, sSig, sRfit);
                const double tol = std::max(row.sigmaR, sSig);
                const double dR = std::fabs(r.radius - s.radius);
                if (dR > tol) continue;
                row.confirmRid = s.id;
                row.confirmNTri = (int)s.tris.size();
                row.confirmR = s.radius;
                break;
            }
            const double sew = mv.sewTol > 0.0 ? mv.sewTol : Precision::Confusion();
            // D-S3-27: [vMin,vMax] overlap/abut within sewTol, in the sibling's
            // axis frame. Local-frame compare is meaningless (each region's
            // Location is its own origin, so stacked holes both look like [0,h]).
            auto vOnS = [](const Region& rr, const Region& ss, double vLocal) {
                const gp_XYZ oR = rr.ax.Location().XYZ();
                const gp_XYZ oS = ss.ax.Location().XYZ();
                const gp_XYZ aS = ss.ax.Direction().XYZ();
                const gp_XYZ aR = rr.ax.Direction().XYZ();
                const double sign = aR.Dot(aS) >= 0.0 ? 1.0 : -1.0;
                return (oR - oS).Dot(aS) + sign * vLocal;
            };
            auto vAbut = [&](const Region& a, const Region& b) {
                const double a0 = vOnS(a, b, a.vMin);
                const double a1 = vOnS(a, b, a.vMax);
                const double aLo = std::min(a0, a1), aHi = std::max(a0, a1);
                const double bLo = std::min(b.vMin, b.vMax), bHi = std::max(b.vMin, b.vMax);
                if (aHi < bLo - sew) return false;
                if (bHi < aLo - sew) return false;
                return true;
            };
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                const Region& s = rs.regions[(size_t)cylIdx[j]];
                if ((int)s.tris.size() <= (int)r.tris.size()) continue;
                if (axisSin(r.ax.Direction(), s.ax.Direction()) > tauAx) continue;
                const double sep = axisSeparation(r, s);
                const double sepTol = meshTolCapInline(s);
                if (sep > sepTol) continue;
                if (!vAbut(r, s)) continue;
                const double dR = std::fabs(r.radius - s.radius);
                if (!(dR > 0.0)) continue;
                const double sCap = meshTolCapInline(s);
                const double oC = sCap > 0.0 ? dR / sCap : dR;
                if (dR > sCap) continue;
                if (row.sibRelRid < 0 || oC < row.sibRelAbsDROverCap) {
                    row.sibRelRid = s.id;
                    row.sibRelR = s.radius;
                    row.sibRelNTri = (int)s.tris.size();
                    row.sibRelSep = sep;
                    row.sibRelAbsDR = dR;
                    row.sibRelAbsDROverCap = oC;
                }
            }
        }
        auto chainAdj = [&](int idA, int idB) {
            for (const BoundaryChain& ch : rs.chains) {
                if ((ch.regA == idA && ch.regB == idB) || (ch.regB == idA && ch.regA == idB))
                    return true;
            }
            return false;
        };
        auto relatedTo = [&](const Region& rr, const Region& ss) {
            if (chainAdj(rr.id, ss.id)) return true;
            if (axisSin(rr.ax.Direction(), ss.ax.Direction()) > tauAx) return false;
            return axisSeparation(rr, ss) <= meshTolCapInline(ss);
        };
        auto carryDev = [&](const Region& rr, const Region& ss) {
            gp_Cylinder cyl = asCyl(ss);
            double maxV = 0.0;
            auto radDev = [&](const gp_Pnt& p) {
                gp_Vec v(cyl.Location(), p);
                return std::fabs(gp_Vec(cyl.Axis().Direction()).Crossed(v).Magnitude() -
                                 cyl.Radius());
            };
            for (int lt : rr.tris) {
                if (lt < 0 || (size_t)lt >= mv.nTri) continue;
                const int* T = mv.tris[mv.compTris[lt]];
                for (int k = 0; k < 3; k++) maxV = std::max(maxV, radDev(gp_Pnt(mv.pts[T[k]])));
            }
            return maxV;
        };
        const double sewAbs = mv.sewTol > 0.0 ? mv.sewTol : Precision::Confusion();
        for (int i = 0; i < n; ++i) {
            const Region& r = rs.regions[(size_t)cylIdx[i]];
            S4Row& row = s4[(size_t)i];
            int bestJ = -1;
            bool anyRelated = false, anyWitness = false;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                const Region& s = rs.regions[(size_t)cylIdx[j]];
                if (!relatedTo(r, s)) continue;
                anyRelated = true;
                if (!(s4[(size_t)j].sigmaR < s4[(size_t)i].sigmaR &&
                      (int)s.tris.size() > (int)r.tris.size()))
                    continue;
                anyWitness = true;
                const double dev = carryDev(r, s);
                const double lim = std::max(s.maxVertexDev, sewAbs);
                if (row.absorbCarryDev < 0.0 || dev < row.absorbCarryDev) {
                    row.absorbCarryDev = dev;
                    row.absorbCarryLim = lim;
                }
                if (dev > lim) continue;
                auto better = [&](int ja, int jb) {
                    const Region& a = rs.regions[(size_t)cylIdx[ja]];
                    const Region& b = rs.regions[(size_t)cylIdx[jb]];
                    const int na = (int)a.tris.size(), nb = (int)b.tris.size();
                    if (na != nb) return na > nb;
                    if (s4[(size_t)ja].sigmaR != s4[(size_t)jb].sigmaR)
                        return s4[(size_t)ja].sigmaR < s4[(size_t)jb].sigmaR;
                    return a.id < b.id;
                };
                if (bestJ < 0 || better(j, bestJ)) bestJ = j;
            }
            if (bestJ >= 0) {
                const Region& s = rs.regions[(size_t)cylIdx[bestJ]];
                row.absorbParent = s.id;
                row.absorbFail = "none";
                row.absorbCarryDev = carryDev(r, s);
                row.absorbCarryLim = std::max(s.maxVertexDev, sewAbs);
            } else {
                row.absorbParent = -1;
                row.absorbFail = !anyRelated ? "related" : (!anyWitness ? "witness" : "carry");
            }
        }
        if (diagOn)
            std::fprintf(stderr, "DIAG_S4FIT_META nCyl=%d tauAx=%.6e epsPlane=%.6e\n", n, tauAx,
                         epsPlane);
    }
    // D-S3-26 absorb map is computed above (DIAG_ABSORB). Topology-merge
    // opened HP PREJ4 (freeEdges=50). Leftover bind of rid 193 alone
    // opened PREJ4 (freeEdges=22). Residue fill (below) is the specified
    // E″ v3 admission; it does not close HP (eprime-leak-r1-d1 mixed
    // edges). Struck E′ (R≤3.5 ∧ ¬closed360 ∧ sinVs>tauAx) is removed.
    for (int i = 0; i < n; ++i) {
        const Region& r = rs.regions[(size_t)cylIdx[i]];
        const double sinVs = axisSin(r.ax.Direction(), ahat);
        const S4Row& row = s4[(size_t)i];
        const bool absorbed = row.absorbParent >= 0;
        const bool unident = row.sigmaROverR > 0.003;
        const bool confirmed = row.confirmRid >= 0;
        // E″ v3 residue: unidentifiable ∧ ¬confirmed ∧ ¬∃absorb.
        // Absorb construction (merge / leftover bind) opened PREJ4; residue
        // fill is the remaining admission. siblingInconsistent without a
        // parent is not leftover-bound (handle-lock rid 12).
        const bool residue = !absorbed && unident && !confirmed;
        const bool demote = residue;
        if (demote && r.id >= 0 && (size_t)r.id < gEPrimeDemoted.size())
            gEPrimeDemoted[(size_t)r.id] = 1;
        if (diagOn)
            std::fprintf(stderr,
                         "DIAG_EPRIME rid=%d R=%.6f nTri=%zu closed360=%d "
                         "sinVsDom=%.6e demoted=%d absorbed=%d leftoverSib=%d residue=%d\n",
                         r.id, r.radius, r.tris.size(), r.closed360 ? 1 : 0, sinVs,
                         demote ? 1 : 0, absorbed ? 1 : 0, row.absorbParent, residue ? 1 : 0);
        if (diagOn) {
            double maxV = 0.0, chordMid = 0.0;
            fitResiduals(r, maxV, chordMid);
            const double cap = meshTolCapInline(r);
            const double overR = r.radius > 0.0 ? maxV / r.radius : 0.0;
            const double overCap = cap > 0.0 ? maxV / cap : 0.0;
            const double chordOverCap = cap > 0.0 ? chordMid / cap : 0.0;
            double minSibAbsDROverR = -1.0;
            double minSibAbsDROverCap = -1.0;
            int sibRid = -1, sibNTri = -1;
            double sibR = -1.0;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                const Region& s = rs.regions[(size_t)cylIdx[j]];
                if ((int)s.tris.size() <= (int)r.tris.size()) continue;
                if (axisSin(r.ax.Direction(), s.ax.Direction()) > tauAx) continue;
                const double sep = axisSeparation(r, s);
                if (sep > epsPlane) continue;
                const double dR = std::fabs(r.radius - s.radius);
                if (!(dR > 0.0)) continue;
                const double oR = r.radius > 0.0 ? dR / r.radius : dR;
                const double sCap = meshTolCapInline(s);
                const double oC = sCap > 0.0 ? dR / sCap : dR;
                if (minSibAbsDROverR < 0.0 || oR < minSibAbsDROverR) {
                    minSibAbsDROverR = oR;
                    minSibAbsDROverCap = oC;
                    sibRid = s.id;
                    sibNTri = (int)s.tris.size();
                    sibR = s.radius;
                }
            }
            const S4Row& row = s4[(size_t)i];
            const int unident = row.sigmaROverR > 0.003 ? 1 : 0;
            const int confirmed = row.confirmRid >= 0 ? 1 : 0;
            const int sibIncon = row.sibRelRid >= 0 ? 1 : 0;
            std::fprintf(stderr,
                         "DIAG_S4FIT rid=%d R=%.9f nTri=%zu nSides=%d closed360=%d "
                         "sinVs=%.6e spanU=%.6f sigmaR=%.9e sigmaROverR=%.9e "
                         "clusterId=%d clusterScatter=%.9e chainN=%d chainBestRid=%d "
                         "chainBestR=%.9f chainBestNTri=%d chainBestAbsDR=%.9e "
                         "confirmRid=%d confirmR=%.9f confirmNTri=%d "
                         "sibRelRid=%d sibRelR=%.9f sibRelNTri=%d sibRelSep=%.9e "
                         "sibRelAbsDR=%.9e sibRelAbsDROverCap=%.9e "
                         "unidentifiable=%d confirmed=%d siblingInconsistent=%d "
                         "maxVertexDev=%.9f maxVertexDevOverR=%.9e "
                         "chordSagitta=%.9f meshTolCap=%.9f minSibAbsDROverR=%.9e "
                         "demoted=%d maxVertexDevOverCap=%.9e chordOverCap=%.9e "
                         "minSibAbsDROverCap=%.9e sibRid=%d sibR=%.9f sibNTri=%d "
                         "sagStored=%.9f maxDevStored=%.9f\n",
                         r.id, r.radius, r.tris.size(), r.nSides, r.closed360 ? 1 : 0, sinVs,
                         row.spanU, row.sigmaR, row.sigmaROverR, row.clusterId,
                         row.clusterScatter, row.chainN, row.chainBestRid, row.chainBestR,
                         row.chainBestNTri, row.chainBestAbsDR, row.confirmRid, row.confirmR,
                         row.confirmNTri, row.sibRelRid, row.sibRelR, row.sibRelNTri,
                         row.sibRelSep, row.sibRelAbsDR, row.sibRelAbsDROverCap, unident,
                         confirmed, sibIncon, maxV, overR, chordMid, cap,
                         minSibAbsDROverR, demote ? 1 : 0, overCap, chordOverCap,
                         minSibAbsDROverCap, sibRid, sibR, sibNTri, r.chordSagitta,
                         r.maxVertexDev);
            std::fprintf(stderr,
                         "DIAG_ABSORB rid=%d parent=%d fail=%s nTri=%zu sigmaR=%.9e "
                         "carryDev=%.9e carryLim=%.9e unidentifiable=%d confirmed=%d "
                         "siblingInconsistent=%d\n",
                         r.id, row.absorbParent, row.absorbFail, r.tris.size(), row.sigmaR,
                         row.absorbCarryDev, row.absorbCarryLim, unident, confirmed, sibIncon);
        }
    }
    if (diagOn) {
        for (const Region& r : rs.rejected) {
            if (r.type != SurfType::Cylinder) continue;
            std::fprintf(stderr,
                         "DIAG_S4REJ rid=%d R=%.9f nTri=%zu nSides=%d closed360=%d "
                         "reject=%d maxDevStored=%.9f\n",
                         r.id, r.radius, r.tris.size(), r.nSides, r.closed360 ? 1 : 0,
                         (int)r.reject, r.maxVertexDev);
        }
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

// D-130: plate-face validity budget. The plane analogue of partialFaceTolCap.
// A plate's boundary edges can never sit closer to the fitted plane than the
// region's own vertices do, so the face-validity cap must carry the region's
// measured planar residual (Region::maxVertexDev, the same quantity the partial
// cylinder cap already uses). meshTolCap alone is the sew budget and silently
// refuses a cap plane whose rim rides an off-plane mesh vertex.
double plateFaceTolCap(const MeshView& mv, const Region& r) {
    double c = meshTolCap(mv, &r);
    if (r.maxVertexDev > 0.0) c = std::max(c, r.maxVertexDev * 1.001 + Precision::Confusion());
    return c;
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

extern thread_local int gTolRewriteFire;
void stampTolWriter(const TopoDS_Edge& e, const char* site);
void fireTolRewriteEdge(BRep_Builder& B, const TopoDS_Edge& e, double after, const char* site,
                        double meshCap = -1.0);

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
            if (d > BRep_Tool::Tolerance(e)) {
                fireTolRewriteEdge(B, e,
                                   std::min(d * 1.001 + Precision::Confusion(), cap), "other",
                                   cap);
            }
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

struct AnalyticCurve {
    enum Kind { None, Lin, Circ, Elips } kind = None;
    gp_Lin lin;
    gp_Circ circ;
    gp_Elips elips;
};

struct ChainGeom {
    bool collapsed = false;           // analytic (possibly 1 or 2 edges)
    std::vector<TopoDS_Edge> edges;   // collapsed edges, in walk-forward (regA-left) order
    AnalyticCurve curve;              // birth curve (rebuildCollapsed)
    int ia = -1, ib = -1;             // terminal verts[] slots
    gp_Pnt midHint;
    int regA = -1, regB = -1;
    double chainSnapCap = 0;
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
    // F4/F5: never inflate a shared TShape past the mesh budget (1 mm,
    // same ceiling as plane|plane acceptCap). 8–13 mm tols collapse distinct
    // hole vertices (dist=0) and MakeEdge of Circ arcs fails.
    d = std::min(d, 1.0);
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

int intersectLinLin(const gp_Lin& a, const gp_Lin& b, gp_Pnt& out) {
    const gp_Vec d1(a.Direction()), d2(b.Direction());
    const gp_Vec n = d1.Crossed(d2);
    const double n2 = n.SquareMagnitude();
    if (n2 < Precision::Confusion() * Precision::Confusion()) return 0;
    const gp_Vec w(a.Location(), b.Location());
    const double t = w.Crossed(d2).Dot(n) / n2;
    out = a.Location().Translated(d1 * t);
    return 1;
}

int intersectLinCirc(const gp_Lin& lin, const gp_Circ& circ, gp_Pnt& pA, gp_Pnt& pB) {
    const gp_Pnt c = circ.Location();
    const gp_Dir ax = circ.Axis().Direction();
    const gp_Dir D = lin.Direction();
    const gp_Pnt o = lin.Location();
    const double tFoot = gp_Vec(o, c).Dot(D);
    const gp_Pnt foot = o.Translated(gp_Vec(D) * tFoot);
    const double off = gp_Vec(c, foot).Dot(ax);
    const gp_Pnt footP = foot.Translated(gp_Vec(ax) * (-off));
    const double dist = footP.Distance(c);
    const double R = circ.Radius();
    if (dist > R + Precision::Confusion()) return 0;
    gp_Vec along = gp_Vec(D) - gp_Vec(ax) * D.Dot(ax);
    if (along.Magnitude() < Precision::Confusion()) return 0;
    along.Normalize();
    try {
        if (dist >= R - Precision::Confusion()) {
            gp_Vec rad(c, footP);
            if (rad.Magnitude() < Precision::Confusion()) return 0;
            pA = ElCLib::Value(ElCLib::Parameter(circ, c.Translated(rad.Normalized() * R)), circ);
            return 1;
        }
        const double h = std::sqrt(std::max(0.0, R * R - dist * dist));
        pA = ElCLib::Value(ElCLib::Parameter(circ, footP.Translated(along * (-h))), circ);
        pB = ElCLib::Value(ElCLib::Parameter(circ, footP.Translated(along * h)), circ);
        return 2;
    } catch (const Standard_Failure&) {
        return 0;
    }
}

int intersectCircCirc(const gp_Circ& c1, const gp_Circ& c2, gp_Pnt& pA, gp_Pnt& pB) {
    gp_Pnt o1 = c1.Location(), o2 = c2.Location();
    const double r1 = c1.Radius(), r2 = c2.Radius();
    gp_Vec d(o1, o2);
    const double dist = d.Magnitude();
    if (dist < Precision::Confusion()) return 0;
    if (dist > r1 + r2 + Precision::Confusion()) return 0;
    if (dist < std::fabs(r1 - r2) - Precision::Confusion()) return 0;
    const double a = (r1 * r1 - r2 * r2 + dist * dist) / (2.0 * dist);
    const double h2 = r1 * r1 - a * a;
    if (h2 < -Precision::Confusion() * Precision::Confusion()) return 0;
    const double h = h2 <= 0.0 ? 0.0 : std::sqrt(h2);
    gp_Dir dir(d);
    gp_Pnt mid = o1.Translated(gp_Vec(dir) * a);
    gp_Vec perp = gp_Vec(c1.Axis().Direction()) ^ gp_Vec(dir);
    if (perp.Magnitude() < Precision::Confusion()) return 0;
    perp.Normalize();
    pA = mid.Translated(perp * h);
    pB = mid.Translated(perp * (-h));
    return h <= Precision::Confusion() ? 1 : 2;
}

int intersectAnalyticCurves(const AnalyticCurve& a, const AnalyticCurve& b, gp_Pnt& pA, gp_Pnt& pB) {
    if (a.kind == AnalyticCurve::Lin && b.kind == AnalyticCurve::Lin) {
        const int n = intersectLinLin(a.lin, b.lin, pA);
        pB = pA;
        return n;
    }
    if (a.kind == AnalyticCurve::Lin && b.kind == AnalyticCurve::Circ)
        return intersectLinCirc(a.lin, b.circ, pA, pB);
    if (a.kind == AnalyticCurve::Circ && b.kind == AnalyticCurve::Lin)
        return intersectLinCirc(b.lin, a.circ, pA, pB);
    if (a.kind == AnalyticCurve::Circ && b.kind == AnalyticCurve::Circ)
        return intersectCircCirc(a.circ, b.circ, pA, pB);
    return 0;
}

// D-S3-72: per-surface deviation of a point vs a fitted Region.
double regionPointDev(const Region* R, const gp_Pnt& p) {
    if (!R) return 1e300;
    try {
        if (R->type == SurfType::Plane)
            return std::fabs(gp_Vec(R->ax.Location(), p).Dot(R->ax.Direction()));
        if (R->type == SurfType::Cylinder) {
            const gp_XYZ a(R->ax.Direction().X(), R->ax.Direction().Y(), R->ax.Direction().Z());
            const gp_XYZ c(R->ax.Location().X(), R->ax.Location().Y(), R->ax.Location().Z());
            const gp_XYZ px(p.X(), p.Y(), p.Z());
            const double radial = a.Crossed(px - c).Modulus();
            return std::fabs(radial - R->radius);
        }
    } catch (const Standard_Failure&) {
    }
    return 1e300;
}

double regionDevBound(const MeshView& mv, const Region* R) {
    double b = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
    if (R && R->maxVertexDev > b) b = R->maxVertexDev;
    return b;
}

double dihedralBetweenRegions(const Region* A, const Region* B, const gp_Pnt& p) {
    if (!A || !B) return -1.0;
    gp_Dir nA, nB;
    try {
        if (A->type == SurfType::Plane)
            nA = A->ax.Direction();
        else if (A->type == SurfType::Cylinder) {
            gp_Vec rho(A->ax.Location(), p);
            rho -= gp_Vec(A->ax.Direction()) * rho.Dot(A->ax.Direction());
            if (rho.Magnitude() < Precision::Confusion()) return 0.0;
            nA = gp_Dir(rho);
        } else
            return -1.0;
        if (B->type == SurfType::Plane)
            nB = B->ax.Direction();
        else if (B->type == SurfType::Cylinder) {
            gp_Vec rho(B->ax.Location(), p);
            rho -= gp_Vec(B->ax.Direction()) * rho.Dot(B->ax.Direction());
            if (rho.Magnitude() < Precision::Confusion()) return 0.0;
            nB = gp_Dir(rho);
        } else
            return -1.0;
        const double c = std::clamp(nA.Dot(nB), -1.0, 1.0);
        return std::acos(c);
    } catch (const Standard_Failure&) {
        return -1.0;
    }
}

bool pointOnCurveSpan(const AnalyticCurve& c, const gp_Pnt& q, const gp_Pnt& pa, const gp_Pnt& pb) {
    if (curveResidual(c, q) > Precision::Confusion() * 10.0) return false;
    if (c.kind == AnalyticCurve::Lin) {
        const double t = ElCLib::Parameter(c.lin, q);
        double t0 = ElCLib::Parameter(c.lin, pa);
        double t1 = ElCLib::Parameter(c.lin, pb);
        if (t0 > t1) std::swap(t0, t1);
        return t >= t0 - Precision::Confusion() && t <= t1 + Precision::Confusion();
    }
    if (c.kind == AnalyticCurve::Circ) {
        const double chord = pa.Distance(pb);
        const double viaQ = pa.Distance(q) + q.Distance(pb);
        return viaQ <= chord + c.circ.Radius() * 0.01 + Precision::Confusion();
    }
    return false;
}

bool projectOntoCurveSpan(const AnalyticCurve& c, const gp_Pnt& p, const gp_Pnt& pa, const gp_Pnt& pb,
                          gp_Pnt& out) {
    try {
        if (c.kind == AnalyticCurve::Lin) {
            double t = ElCLib::Parameter(c.lin, p);
            double t0 = ElCLib::Parameter(c.lin, pa);
            double t1 = ElCLib::Parameter(c.lin, pb);
            if (t0 > t1) std::swap(t0, t1);
            if (t1 - t0 <= Precision::Confusion()) return false;
            t = std::max(t0, std::min(t1, t));
            out = ElCLib::Value(t, c.lin);
            return true;
        }
        if (c.kind == AnalyticCurve::Circ) {
            const gp_Pnt ctr = c.circ.Location();
            const gp_Vec ax(c.circ.Axis().Direction());
            const gp_Vec rho(ctr, p);
            const gp_Vec radial = rho - ax * rho.Dot(ax);
            if (radial.Magnitude() <= Precision::Confusion()) return false;
            const double t = ElCLib::Parameter(c.circ, p);
            out = ElCLib::Value(t, c.circ);
            if (pa.Distance(pb) <= Precision::Confusion()) return true;
            return pointOnCurveSpan(c, out, pa, pb);
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

gp_Circ cylinderIsoCircle(const Region& cyl, double v) {
    gp_Pnt loc = cyl.ax.Location().Translated(gp_Vec(cyl.ax.Direction()) * v);
    return gp_Circ(gp_Ax2(loc, cyl.ax.Direction(), cyl.ax.XDirection()), cyl.radius);
}

bool isChamferConeR(const Region& r) {
    return r.type == SurfType::Cone && r.origin == Origin::ChamferCone;
}

// SIGNED axial offset from the R_lo rim (the region's Location) to the R_hi rim,
// measured along +Direction. Negative when the taper runs against the canonical
// axis; detector C never flips the axis, it carries the sign here, so a cone and
// the cylinder it meets always share one frame (refit_chamfer_cone.cpp).
double chamferHeightOf(const Region& r) {
    if (std::fabs(r.maxVertexSnap) > Precision::Confusion()) return r.maxVertexSnap;
    return std::fabs(r.vMax - r.vMin);
}

double chamferRhiOf(const Region& r) {
    if (r.vMax > r.radius + Precision::Confusion()) return r.vMax;
    return r.radius + std::fabs(chamferHeightOf(r));
}

gp_Circ coneIsoCircle(const Region& cone, double vAxial) {
    const double h = chamferHeightOf(cone);
    const double Rhi = chamferRhiOf(cone);
    // t is a fraction of the SIGNED height, so t = 1 is the R_hi rim for either
    // taper sign and rho(v) = R_lo + t*(R_hi - R_lo) needs no case analysis.
    const double t = (std::fabs(h) > Precision::Confusion()) ? (vAxial / h) : 0.0;
    double R = cone.radius + t * (Rhi - cone.radius);
    if (!(R > Precision::Confusion())) R = cone.radius;
    gp_Pnt loc = cone.ax.Location().Translated(gp_Vec(cone.ax.Direction()) * vAxial);
    return gp_Circ(gp_Ax2(loc, cone.ax.Direction(), cone.ax.XDirection()), R);
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

AnalyticCurve constructedPlaneConeCap(const Region& plnR, const Region& coneR) {
    AnalyticCurve out;
    if (!planePerpCylinder(plnR, coneR)) return out;
    out.kind = AnalyticCurve::Circ;
    out.circ = coneIsoCircle(coneR, planeVOnCylinder(plnR, coneR));
    return out;
}

AnalyticCurve constructedCylConeRim(const Region& cylR, const Region& coneR) {
    AnalyticCurve out;
    const double d0 = std::fabs(cylR.radius - coneR.radius);
    const double d1 = std::fabs(cylR.radius - chamferRhiOf(coneR));
    out.kind = AnalyticCurve::Circ;
    if (d0 <= d1)
        out.circ = coneIsoCircle(coneR, 0.0);
    else
        out.circ = coneIsoCircle(coneR, chamferHeightOf(coneR));
    return out;
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

// D-130-14 / D-130-16 -- the union's landing stage (see refit_grow.cpp).
// With it off every face this file builds is exactly what the branch tip built.
bool unionBuildOn() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_UNION");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

bool diagP2Enabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_P2_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// D-130-16 -- THE TWO-ARC RIM COLLAPSE.
//
// A closed360 cylinder that carries an enclosed interruption ships as ONE face
// with that interruption as an inner wire (D-130-14). The face's seam is a
// GENERATOR of the same face, so it must stand at an azimuth the inner wire
// does not cover, on a vertex both cap wires actually have. On the mesh that
// produced this ruling the cap's 288 deg rim chain collapses into a SINGLE
// analytic arc, so the cap wire keeps only six vertices -- and all six are the
// tooth-arc junctions, every one of them inside the inner wire's span. The seam
// has nowhere to stand and the face is UnorientableShape.
//
// The rim's own mesh has twenty vertices outside that span; the collapse is
// what discards them. So the rim collapses into TWO arcs of the same circle
// instead of one, meeting at the chain's own mesh vertex nearest the seam
// target. Nothing is invented and nothing is widened: the meeting point is a
// real mesh vertex of that chain, the rim stays analytic (two arcs of one
// circle, not a polyline), both faces still share the same edges, and the split
// vertex must lie on the chain's own analytic curve within the chain's own
// snap cap -- the admission the terminals already get. The neighbouring plane's
// hole gains exactly that one real vertex.
//
// Returns the local mesh vertex to split at, or -1 when the chain is not such a
// rim (then the collapse is the single arc it has always been).
int rimSplitVertex(const RegionSet& rs, const MeshView& mv, const BoundaryChain& ch, int ci,
                   const Region* A, const Region* B, double snapCap, const AnalyticCurve& curve) {
    if (ch.closedLoop || ch.meshVerts.size() < 3) return -1;
    if (curve.kind != AnalyticCurve::Circ) return -1;
    for (const Region* R : {A, B}) {
        if (!R || R->type != SurfType::Cylinder || !R->closed360) continue;
        // ... carrying an inner wire, and this chain on one of its rims.
        bool onCap = false;
        std::vector<double> innerU;
        for (const Loop& lp : R->loops) {
            if (lp.role == LoopRole::CapLow || lp.role == LoopRole::CapHigh) {
                for (int k : lp.chainIdx)
                    if (k == ci) onCap = true;
            } else if (lp.role == LoopRole::Inner) {
                for (int k : lp.chainIdx) {
                    if (k < 0 || (size_t)k >= rs.chains.size()) continue;
                    for (int lv : rs.chains[(size_t)k].meshVerts)
                        innerU.push_back(azimuthOf(*R, pntOf(mv, lv)));
                }
            }
        }
        if (!onCap || innerU.empty()) continue;
        // The seam target is the middle of the widest azimuth gap between the
        // interruptions -- the same measurement trySeamed360 makes, of this
        // region and nothing else. "Outside the inner wire's span" is exactly
        // "strictly inside that gap".
        std::sort(innerU.begin(), innerU.end());
        double gap = innerU.front() + 2.0 * kPi - innerU.back();
        double gLo = innerU.back(), gHi = innerU.front() + 2.0 * kPi;
        for (size_t i = 1; i < innerU.size(); i++) {
            const double g = innerU[i] - innerU[i - 1];
            if (g > gap) {
                gap = g;
                gLo = innerU[i - 1];
                gHi = innerU[i];
            }
        }
        const double target = 0.5 * (gLo + gHi);
        int best = -1;
        double bestD = 1e300;
        for (size_t i = 1; i + 1 < ch.meshVerts.size(); i++) {  // interior vertices only
            const int lv = ch.meshVerts[i];
            const gp_Pnt p = pntOf(mv, lv);
            if (curveResidual(curve, p) > snapCap) continue;  // not on this rim's own circle
            double u = azimuthOf(*R, p);
            while (u < gLo) u += 2.0 * kPi;
            if (!(u > gLo && u < gHi)) continue;              // inside the interruption's span
            const double d = std::fabs(u - target);
            if (d < bestD || (d == bestD && lv < best)) {
                bestD = d;
                best = lv;
            }
        }
        if (best >= 0) {
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_RIMSPLIT ci=%d rid=%d nV=%zu gap=[%.5f,%.5f] target=%.5f "
                             "splitLv=%d u=%.5f\n",
                             ci, R->id, ch.meshVerts.size(), gLo, gHi, target, best,
                             azimuthOf(*R, pntOf(mv, best)));
            return best;
        }
    }
    return -1;
}

bool diag130Enabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_DIAG_130");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// Isolated TShape copy for print-only try-B / BRepCheck. Never fall back to live.
TopoDS_Shape diagIsolatedCopy(const TopoDS_Shape& s) {
    if (s.IsNull()) return {};
    try {
        BRepBuilderAPI_Copy cop(s, Standard_True);
        if (cop.IsDone() && !cop.Shape().IsNull()) return cop.Shape();
    } catch (const Standard_Failure&) {
    }
    return {};
}

void dumpDiagPlateLine(int rid, const TopoDS_Face& f) {
    const TopoDS_Shape iso = diagIsolatedCopy(f);
    if (iso.IsNull()) return;
    try {
        BRepCheck_Analyzer an(iso, Standard_True);
        std::fprintf(stderr, "DIAG_PLATE rid=%d valid=%d", rid, an.IsValid() ? 1 : 0);
        if (!an.IsValid()) {
            Handle(BRepCheck_Result) res = an.Result(iso);
            if (!res.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next())
                    std::fprintf(stderr, " st=%d", (int)it.Value());
            }
        }
        std::fprintf(stderr, "\n");
        try {
            GProp_GProps sp;
            BRepGProp::SurfaceProperties(f, sp);
            std::fprintf(stderr, "DIAG_GPROP_FACE rid=%d area=%.6f\n", rid, sp.Mass());
        } catch (const Standard_Failure&) {
        }
    } catch (const Standard_Failure&) {
    }
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
thread_local bool gHubWireRowB = false;
thread_local std::unordered_set<int> gChainSewFallbackCi;
thread_local std::unordered_map<const void*, int> gWirePopCi;
thread_local std::unordered_map<const void*, const char*> gWirePopSite;
thread_local int gReconnectFires = 0;
thread_local int gOriNotDetermined = 0;
struct WireOriTally {
    int reversedByRule = 0;
    int notDetermined = 0;
    int seedReversed = 0;
    int flipped = 0;
    int ci0 = -1;
    const char* seedOrientedBy = "none";
    const char* loopLabel = "?";
};
thread_local std::vector<WireOriTally> gPlateWireOri;

void setHubWireRowB() { gHubWireRowB = true; }
const char* hubWireRowLabel() { return gHubWireRowB ? "B" : "A"; }
void noteChainSewFallback(int ci) {
    if (ci >= 0) gChainSewFallbackCi.insert(ci);
}

const void* diagTShapePtr(const TopoDS_Shape& s) {
    if (s.IsNull()) return nullptr;
    return (const void*)s.TShape().get();
}

void noteWirePop(const TopoDS_Edge& e, const char* site, int ci) {
    if (e.IsNull() || !site) return;
    const void* ts = diagTShapePtr(e);
    if (!ts) return;
    gWirePopSite[ts] = site;
    gWirePopCi[ts] = ci;
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
        if (!collapsed.empty() && (ci >= collapsed.size() || !collapsed[ci])) continue;
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
        if ((A.type == SurfType::Plane && isChamferConeR(B)) ||
            (B.type == SurfType::Plane && isChamferConeR(A))) {
            const Region& pln = (A.type == SurfType::Plane) ? A : B;
            const Region& cn = isChamferConeR(A) ? A : B;
            AnalyticCurve c = constructedPlaneConeCap(pln, cn);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: plane|cone empty — keeping mesh polyline");
            return c;
        }
        if ((A.type == SurfType::Cylinder && isChamferConeR(B)) ||
            (B.type == SurfType::Cylinder && isChamferConeR(A))) {
            const Region& cyl = (A.type == SurfType::Cylinder) ? A : B;
            const Region& cn = isChamferConeR(A) ? A : B;
            AnalyticCurve c = constructedCylConeRim(cyl, cn);
            if (c.kind == AnalyticCurve::None)
                emit(warn, "smooth: cyl|cone empty — keeping mesh polyline");
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
                double p1 = ElCLib::Parameter(c.circ, BRep_Tool::Pnt(v1));
                double p2 = ElCLib::Parameter(c.circ, BRep_Tool::Pnt(v2));
                if (p2 < p1) p2 += 2.0 * kPi;
                if (std::fabs(p2 - p1) <= Precision::Angular()) return {};
                Handle(Geom_Circle) gc = new Geom_Circle(c.circ);
                TopoDS_Edge pe = bindEdgeByParam(gc, v1, v2, p1, p2);
                if (!pe.IsNull()) return pe;
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
        bumpVertexTol(vA, Precision::Confusion() * 10.0);
        bumpVertexTol(vB, Precision::Confusion() * 10.0);
        if (dm <= df + 1e-9) {
            TopoDS_Edge pe = bindEdgeByParam(gc, vA, vB, p1, p1 + df);
            if (!pe.IsNull()) return pe;
        } else {
            const double db = twopi - df;
            const double t0 = wrap(p1 - db);
            TopoDS_Edge pe = bindEdgeByParam(gc, vB, vA, t0, t0 + db);
            if (!pe.IsNull()) {
                pe.Reverse();
                return pe;
            }
        }
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

// True when the increasing-parameter arc vA -> vB is the one the mesh mid
// vertex sits on. Same selector as makeEllipseArc's own branch, factored out so
// the fallback cannot answer with the complement.
bool midHintOnForwardArc(const gp_Elips& el, const TopoDS_Vertex& vA, const TopoDS_Vertex& vB,
                         const gp_Pnt& midHint) {
    const double twopi = 2.0 * kPi;
    auto wrap = [&](double t) {
        while (t < 0.0) t += twopi;
        while (t >= twopi) t -= twopi;
        return t;
    };
    try {
        const double p1 = wrap(ElCLib::Parameter(el, BRep_Tool::Pnt(vA)));
        const double p2 = wrap(ElCLib::Parameter(el, BRep_Tool::Pnt(vB)));
        const double pm = wrap(ElCLib::Parameter(el, midHint));
        double df = p2 - p1;
        while (df <= 0.0) df += twopi;
        double dm = pm - p1;
        while (dm < 0.0) dm += twopi;
        return dm <= df + 1e-9;
    } catch (const Standard_Failure&) {
    }
    return true;
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
        if (dm <= df + 1e-9) {
            TopoDS_Edge pe = bindEdgeByParam(ge, vA, vB, p1, p1 + df);
            if (!pe.IsNull()) return pe;
        } else {
            // Complementary arc: the mid vertex says the sweep runs the other
            // way round, so vB owns p1 - db, not p1 + db. Naming p1 + db binds
            // vB to a point it does not occupy (90 deg away on S11-b), which
            // BRepBuilderAPI_MakeEdge refuses -- and the arc-blind fallback
            // below then hands back the MAJOR arc. Bind vB -> vA over its own
            // range and reverse, exactly as makeArc does for a circle.
            const double db = twopi - df;
            const double t0 = wrap(p1 - db);
            TopoDS_Edge pe = bindEdgeByParam(ge, vB, vA, t0, t0 + db);
            if (!pe.IsNull()) {
                pe.Reverse();
                return pe;
            }
        }
    } catch (const Standard_Failure&) {
    }
    // ME_ELLIPSE_FALLBACK: this edge is the increasing-parameter arc from vA to
    // vB, so it is only the arc the mesh asked for when the mid vertex lies on
    // it. Handing back the complement instead is a silent geometry swap; the
    // chain must demote to its mesh polyline, which every neighbour can share.
    if (!midHintOnForwardArc(el, vA, vB, midHint)) return TopoDS_Edge();
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

bool collapseDiagEnabled();

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
    // The three points this wire actually joins. They are NOT mv.pts[]: a shared
    // vertex on the boundary of an analytic region has been snapped onto that
    // region's chain curve (snapVertexToCurve), and the mesh array does not move
    // with it.
    gp_XYZ wp[3];
    bool haveWp = true;
    for (int s = 0; s < 3; s++) {
        int id = mv.triEdges[k][s];
        if (id < 0 || (size_t)id >= meshE.size() || !edgeOk[(size_t)id]) return TopoDS_Face();
        bool fwd = (mv.triDirs[k] >> s) & 1;
        diagNoteEdge(meshE[(size_t)id], "makeFacet");
        const TopoDS_Edge oe =
            fwd ? meshE[(size_t)id] : TopoDS::Edge(meshE[(size_t)id].Reversed());
        const TopoDS_Vertex v0 = TopExp::FirstVertex(oe, Standard_True);
        if (v0.IsNull()) haveWp = false;
        else wp[s] = BRep_Tool::Pnt(v0).XYZ();
        bb.Add(w, oe);
    }
    w.Closed(Standard_True);
    // A facet is filled on the plane of its OWN wire. Taking the plane from
    // mv.pts[] instead leaves the face not containing the edges it is built
    // from: measured 3.21725e-07 mm on the plate's cross-bore island facet
    // (vertex snapped 4.4e-06 mm onto the R10 generator), which survives
    // construction only because the shared edge still carries the mesh
    // tolerance. BRepLib::SameParameter(forced) then recomputes that tolerance
    // from the STORED representations -- and a facet plane stores no pcurve
    // (ban N3) -- so the edge drops to Precision::Confusion(), the face reads
    // InvalidCurveOnSurface in its own context and BRepCheck_Face cannot orient
    // its single wire: st=27 UnorientableShape. Refill locally, and only when a
    // wire vertex has actually left the mesh plane by more than the tolerance
    // the face will carry; where nothing was snapped the two planes are the same
    // numbers and the face is unchanged. The mesh normal keeps the casting vote:
    // a refill may never flip the facet's side or degenerate it.
    gp_XYZ pA = A;
    gp_XYZ pn = n;
    if (haveWp) {
        const gp_XYZ nHat = n / mag;
        double off = 0.0;
        for (int s = 0; s < 3; s++) off = std::max(off, std::abs((wp[s] - A).Dot(nHat)));
        if (off > Precision::Confusion()) {
            const gp_XYZ n2 = (wp[1] - wp[0]).Crossed(wp[2] - wp[0]);
            const double mag2 = n2.Modulus();
            if (mag2 > 0.5 * mag && n2.Dot(n) > 0.0) {
                pA = wp[0];
                pn = n2;
                if (collapseDiagEnabled())
                    std::fprintf(stderr,
                                 "DIAG_FACETREFILL k=%zu gt=%d off=%.6g p=(%.7f,%.7f,%.7f)\n", k,
                                 gt, off, wp[0].X(), wp[0].Y(), wp[0].Z());
            }
        }
    }
    // BRep_Builder::MakeFace + Add keeps the wire's verts[] slots (J2).
    // MakeFace(gp_Pln, wire) copies vertices and is the explode-path twin source.
    try {
        Handle(Geom_Plane) pln = new Geom_Plane(gp_Pln(gp_Pnt(pA), gp_Dir(pn)));
        TopoDS_Face f;
        bb.MakeFace(f, pln, Precision::Confusion());
        bb.Add(f, w);
        return f;
    } catch (const Standard_Failure&) {
        return TopoDS_Face();
    }
}

void bindEdgePcurveOnInternedPlane(const TopoDS_Edge& e, const Region& plnR, const MeshView& mv,
                                   const char* kind, int ci);
extern thread_local std::unordered_set<const void*> gSeamTShapes;

struct WireAppendOri {
    bool havePrev = false;
    TopoDS_Edge prev;
    bool seedPreset = false;
    TopoDS_Edge seedPresetEdge;
    bool wholeLoopAsBuilt = false;
    bool useSlotOri = false;
    std::vector<TopoDS_Edge> slotOri;
    size_t slotConsume = 0;
    const char* seedOrientedBy = "none";
    int reversedByRule = 0;
    int notDetermined = 0;
    int seedReversed = 0;
};

bool sameVtxTShape(const TopoDS_Vertex& a, const TopoDS_Vertex& b) {
    if (a.IsNull() || b.IsNull()) return false;
    return diagTShapePtr(a) == diagTShapePtr(b);
}

// OP2b hypothesis print (D-S3-129): bind/pin after Reversed() local handle?
// bindAfter=1 only if UpdateEdge/pin ran on the handed (possibly reversed) handle.
void dumpOp2bBind(const TopoDS_Edge& eHanded, const TopoDS_Edge& eNative, int ci, const char* site,
                  int reversedThis, int bindAfter, const char* bindKind, const Region* plateR,
                  const TopoDS_Vertex& contVtx) {
    if (!diagP2Enabled() || !plateR || plateR->id > 5) return;
    const TopoDS_Vertex nF = TopExp::FirstVertex(eNative, Standard_False);
    const TopoDS_Vertex nL = TopExp::LastVertex(eNative, Standard_False);
    const TopoDS_Vertex oF = TopExp::FirstVertex(eHanded, Standard_True);
    const TopoDS_Vertex oL = TopExp::LastVertex(eHanded, Standard_True);
    const int nativeMatch = (sameVtxTShape(contVtx, nF) || sameVtxTShape(contVtx, nL)) ? 1 : 0;
    const int oriMatch = (sameVtxTShape(contVtx, oF) || sameVtxTShape(contVtx, oL)) ? 1 : 0;
    int pcAgree = -1;
    try {
        const gp_Pln pln(plateR->ax);
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(eHanded, f, l);
        if (!c3.IsNull() && l - f > Precision::PConfusion() && !oF.IsNull() && !oL.IsNull()) {
            Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0, uf = 0, vf = 0, ul = 0, vl = 0;
            ElSLib::Parameters(pln, BRep_Tool::Pnt(oF), u0, v0);
            ElSLib::Parameters(pln, BRep_Tool::Pnt(oL), u1, v1);
            ElSLib::Parameters(pln, c3->Value(f), uf, vf);
            ElSLib::Parameters(pln, c3->Value(l), ul, vl);
            const double dot = (u1 - u0) * (ul - uf) + (v1 - v0) * (vl - vf);
            pcAgree = (dot >= 0.0) ? 1 : 0;
        }
    } catch (const Standard_Failure&) {
        pcAgree = -1;
    }
    const char* handedOri = (eHanded.Orientation() == TopAbs_REVERSED) ? "R" : "F";
    std::fprintf(stderr,
                 "DIAG_OP2B_BIND rid=%d ci=%d site=%s handedOri=%s reversedThis=%d bindAfter=%d "
                 "bindKind=%s pcAgree=%d nativeFirstT=%p nativeLastT=%p oriFirstT=%p oriLastT=%p "
                 "contVtxT=%p nativeMatch=%d oriMatch=%d\n",
                 plateR->id, ci, site ? site : "?", handedOri, reversedThis, bindAfter,
                 bindKind ? bindKind : "none", pcAgree, diagTShapePtr(nF), diagTShapePtr(nL),
                 diagTShapePtr(oF), diagTShapePtr(oL), diagTShapePtr(contVtx), nativeMatch,
                 oriMatch);
}

int loopStartMeshVert(const Loop& loop, const RegionSet& rs) {
    if (loop.chainIdx.empty() || loop.reversed.empty()) return -1;
    const int ci = loop.chainIdx[0];
    if (ci < 0 || (size_t)ci >= rs.chains.size()) return -1;
    const BoundaryChain& ch = rs.chains[(size_t)ci];
    if (ch.meshVerts.empty()) return -1;
    const bool rev = loop.reversed[0] != 0;
    if (ch.closedLoop) return ch.meshVerts.front();
    return rev ? ch.meshVerts.back() : ch.meshVerts.front();
}

// D-S3-130 (b'): signed area in the plate UV frame (ElSLib on r.ax). No pcurves.
double loopSignedAreaUV(const Loop& loop, const RegionSet& rs, const MeshView& mv,
                        const Region& plate) {
    const gp_Pln pln(plate.ax);
    std::vector<gp_Pnt2d> uv;
    uv.reserve(16);
    for (size_t i = 0; i < loop.chainIdx.size() && i < loop.reversed.size(); i++) {
        const int ci = loop.chainIdx[i];
        if (ci < 0 || (size_t)ci >= rs.chains.size()) continue;
        std::vector<int> vs = walkVerts(rs.chains[(size_t)ci], loop.reversed[i] != 0);
        for (int lv : vs) {
            if (lv < 0 || (size_t)lv >= mv.nVtx) continue;
            const int gv = mv.compVtx[lv];
            const gp_XYZ p = mv.pts[gv];
            Standard_Real u = 0, v = 0;
            ElSLib::Parameters(pln, gp_Pnt(p), u, v);
            const gp_Pnt2d q(u, v);
            if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
            uv.push_back(q);
        }
    }
    if (uv.size() >= 2 && uv.front().Distance(uv.back()) <= Precision::PConfusion()) uv.pop_back();
    if (uv.size() < 3) return 0.0;
    double a = 0.0;
    for (size_t i = 0; i < uv.size(); i++) {
        const gp_Pnt2d& p0 = uv[i];
        const gp_Pnt2d& p1 = uv[(i + 1) % uv.size()];
        a += p0.X() * p1.Y() - p1.X() * p0.Y();
    }
    return 0.5 * a;
}

double wireSignedAreaUV(const TopoDS_Wire& w, const Region& plate) {
    if (w.IsNull()) return 0.0;
    const gp_Pln pln(plate.ax);
    std::vector<gp_Pnt2d> uv;
    uv.reserve(16);
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
            const TopoDS_Edge e = TopoDS::Edge(ex.Current());
            const TopoDS_Vertex v = TopExp::FirstVertex(e, Standard_True);
            if (v.IsNull()) continue;
            Standard_Real u = 0, vv = 0;
            ElSLib::Parameters(pln, BRep_Tool::Pnt(v), u, vv);
            const gp_Pnt2d q(u, vv);
            if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
            uv.push_back(q);
        }
    } catch (const Standard_Failure&) {
        return 0.0;
    }
    if (uv.size() >= 2 && uv.front().Distance(uv.back()) <= Precision::PConfusion()) uv.pop_back();
    if (uv.size() < 3) return 0.0;
    double a = 0.0;
    for (size_t i = 0; i < uv.size(); i++) {
        const gp_Pnt2d& p0 = uv[i];
        const gp_Pnt2d& p1 = uv[(i + 1) % uv.size()];
        a += p0.X() * p1.Y() - p1.X() * p0.Y();
    }
    return 0.5 * a;
}

// Minimum UV edge length along the loop polygon (geometric resolution of its vertices).
double loopUvMinEdgeLen(const Loop& loop, const RegionSet& rs, const MeshView& mv,
                        const Region& plate) {
    const gp_Pln pln(plate.ax);
    std::vector<gp_Pnt2d> uv;
    uv.reserve(16);
    for (size_t i = 0; i < loop.chainIdx.size() && i < loop.reversed.size(); i++) {
        const int ci = loop.chainIdx[i];
        if (ci < 0 || (size_t)ci >= rs.chains.size()) continue;
        std::vector<int> vs = walkVerts(rs.chains[(size_t)ci], loop.reversed[i] != 0);
        for (int lv : vs) {
            if (lv < 0 || (size_t)lv >= mv.nVtx) continue;
            const int gv = mv.compVtx[lv];
            const gp_XYZ p = mv.pts[gv];
            Standard_Real u = 0, v = 0;
            ElSLib::Parameters(pln, gp_Pnt(p), u, v);
            const gp_Pnt2d q(u, v);
            if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
            uv.push_back(q);
        }
    }
    if (uv.size() >= 2 && uv.front().Distance(uv.back()) <= Precision::PConfusion()) uv.pop_back();
    if (uv.size() < 2) return 0.0;
    double minE = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < uv.size(); i++) {
        const double d = uv[i].Distance(uv[(i + 1) % uv.size()]);
        if (d < minE) minE = d;
    }
    return minE;
}

void reverseLoopWalk(Loop& loop) {
    std::reverse(loop.chainIdx.begin(), loop.chainIdx.end());
    std::reverse(loop.reversed.begin(), loop.reversed.end());
    for (uint8_t& r : loop.reversed) r = r ? 0 : 1;
}

size_t chainWireEdgeCount(int ci, const RegionSet& rs, const std::vector<ChainGeom>& geom,
                          const std::vector<char>& collapsed) {
    if (ci < 0 || (size_t)ci >= rs.chains.size()) return 0;
    if ((size_t)ci < collapsed.size() && collapsed[(size_t)ci] && (size_t)ci < geom.size() &&
        geom[(size_t)ci].collapsed)
        return geom[(size_t)ci].edges.size();
    return rs.chains[(size_t)ci].meshEdges.size();
}

size_t loopWireEdgeCount(const Loop& loop, const RegionSet& rs, const std::vector<ChainGeom>& geom,
                         const std::vector<char>& collapsed) {
    size_t n = 0;
    for (size_t i = 0; i < loop.chainIdx.size(); i++) {
        int ci = loop.chainIdx[i];
        n += chainWireEdgeCount(ci, rs, geom, collapsed);
    }
    return n;
}

bool loopPolylineSegmentEdge(const BoundaryChain& ch, bool reversed, size_t segIx,
                             const MeshView& mv, const std::vector<TopoDS_Edge>& meshE,
                             const std::vector<char>& edgeOk, TopoDS_Edge& out) {
    std::vector<int> vs = walkVerts(ch, reversed);
    if (vs.empty() || ch.meshEdges.empty()) return false;
    const size_t nSeg = ch.meshEdges.size();
    if (segIx >= nSeg) return false;
    const int va = vs[segIx % vs.size()];
    const int vb = vs[(segIx + 1) % vs.size()];
    if (!ch.closedLoop && segIx + 1 >= vs.size()) return false;
    int eid = edgeConnecting(mv, ch, va, vb);
    if (eid < 0) {
        const size_t idx = reversed ? (nSeg - 1 - segIx) : segIx;
        if (idx >= ch.meshEdges.size()) return false;
        eid = ch.meshEdges[idx];
    }
    if (eid < 0 || (size_t)eid >= meshE.size() || !edgeOk[(size_t)eid]) return false;
    const auto& ev = mv.compEdges[eid];
    const bool fwd = (ev.first == va);
    out = fwd ? meshE[(size_t)eid] : TopoDS::Edge(meshE[(size_t)eid].Reversed());
    return true;
}

bool loopChainHandedEdge(size_t chainPos, size_t edgeIx, const Loop& loop, const RegionSet& rs,
                         const MeshView& mv, const std::vector<ChainGeom>& geom,
                         const std::vector<char>& collapsed, const std::vector<TopoDS_Edge>& meshE,
                         const std::vector<char>& edgeOk, TopoDS_Edge& out) {
    if (chainPos >= loop.chainIdx.size()) return false;
    const int ci = loop.chainIdx[chainPos];
    const bool rev = loop.reversed[chainPos] != 0;
    if (ci < 0 || (size_t)ci >= rs.chains.size()) return false;
    if ((size_t)ci < collapsed.size() && collapsed[(size_t)ci] && (size_t)ci < geom.size() &&
        geom[(size_t)ci].collapsed && !geom[(size_t)ci].edges.empty()) {
        const auto& g = geom[(size_t)ci];
        if (edgeIx >= g.edges.size()) return false;
        if (!rev)
            out = g.edges[edgeIx];
        else
            out = TopoDS::Edge(g.edges[g.edges.size() - 1 - edgeIx].Reversed());
        return !out.IsNull();
    }
    return loopPolylineSegmentEdge(rs.chains[(size_t)ci], rev, edgeIx, mv, meshE, edgeOk, out);
}

bool loopWireEdgeAt(size_t wireIx, const Loop& loop, const RegionSet& rs, const MeshView& mv,
                    const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                    const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                    TopoDS_Edge& out, int* outCi) {
    size_t pos = 0;
    for (size_t chainPos = 0; chainPos < loop.chainIdx.size(); chainPos++) {
        const int ci = loop.chainIdx[chainPos];
        const size_t nE = chainWireEdgeCount(ci, rs, geom, collapsed);
        if (wireIx < pos + nE) {
            if (outCi) *outCi = ci;
            return loopChainHandedEdge(chainPos, wireIx - pos, loop, rs, mv, geom, collapsed, meshE,
                                       edgeOk, out);
        }
        pos += nE;
    }
    return false;
}

bool loopWireSuccessorHandedEdge(const Loop& loop, const RegionSet& rs, const MeshView& mv,
                                 const std::vector<ChainGeom>& geom,
                                 const std::vector<char>& collapsed,
                                 const std::vector<TopoDS_Edge>& meshE,
                                 const std::vector<char>& edgeOk, TopoDS_Edge& out, int* outCi) {
    return loopWireEdgeAt(1, loop, rs, mv, geom, collapsed, meshE, edgeOk, out, outCi);
}

bool bothEndsSharedTShape(const TopoDS_Edge& a, const TopoDS_Edge& b) {
    const TopoDS_Vertex aF = TopExp::FirstVertex(a, Standard_True);
    const TopoDS_Vertex aL = TopExp::LastVertex(a, Standard_True);
    const TopoDS_Vertex bF = TopExp::FirstVertex(b, Standard_True);
    const TopoDS_Vertex bL = TopExp::LastVertex(b, Standard_True);
    return (sameVtxTShape(aF, bF) && sameVtxTShape(aL, bL)) ||
           (sameVtxTShape(aF, bL) && sameVtxTShape(aL, bF));
}

struct LoopWireSlot {
    size_t wireIx = 0;
    int ci = -1;
    TopoDS_Edge edge;
};

struct LoopConnWalk {
    bool ok = false;
    bool notDetermined = false;
    int failWireIx = -1;
    std::vector<size_t> order;
    std::vector<char> forward;
};

double vtxMinDist3d(const TopoDS_Vertex& a, const TopoDS_Vertex& b) {
    if (a.IsNull() || b.IsNull()) return 1e300;
    return BRep_Tool::Pnt(a).Distance(BRep_Tool::Pnt(b));
}

double edgeEndpointMinDist(const TopoDS_Vertex& vtx, const TopoDS_Edge& e) {
    const TopoDS_Vertex f = TopExp::FirstVertex(e, Standard_True);
    const TopoDS_Vertex l = TopExp::LastVertex(e, Standard_True);
    return std::min(vtxMinDist3d(vtx, f), vtxMinDist3d(vtx, l));
}

void collectJunctionCandidates(const TopoDS_Vertex& vtx, const std::vector<LoopWireSlot>& slots,
                               const std::vector<char>& visited, size_t fromIx,
                               std::vector<std::pair<size_t, double>>& tsHits,
                               std::vector<std::pair<size_t, double>>& geoNear) {
    tsHits.clear();
    geoNear.clear();
    for (size_t j = 0; j < slots.size(); j++) {
        if (visited[j] || j == fromIx) continue;
        const TopoDS_Edge& e = slots[j].edge;
        const TopoDS_Vertex f = TopExp::FirstVertex(e, Standard_True);
        const TopoDS_Vertex l = TopExp::LastVertex(e, Standard_True);
        if (sameVtxTShape(f, vtx) || sameVtxTShape(l, vtx)) {
            tsHits.push_back({j, 0.0});
            continue;
        }
        const double d = edgeEndpointMinDist(vtx, e);
        if (d < 1e299) geoNear.push_back({j, d});
    }
}

// D-S3-164: walk starts at chainIdx[0] in AS-BUILT orientation and steps from
// that edge's HEAD (Last vertex). Direction is never derived from a start vertex.
bool loopConnectivityWalk(const std::vector<LoopWireSlot>& slots, int rid,
                          const char* loopLabel, LoopConnWalk& out) {
    out = {};
    const size_t n = slots.size();
    if (n == 0) return false;
    const TopoDS_Edge& seed = slots[0].edge;
    if (n == 1) {
        const TopoDS_Vertex sF = TopExp::FirstVertex(seed, Standard_True);
        const TopoDS_Vertex sL = TopExp::LastVertex(seed, Standard_True);
        if (sameVtxTShape(sF, sL)) {
            out.ok = true;
            out.order = {0};
            out.forward = {1};
            return true;
        }
        out.notDetermined = true;
        std::fprintf(stderr,
                     "DIAG_OP2G_JUNCTION rid=%d loop=%s ci=%d candidates=0 dists=single-open\n",
                     rid, loopLabel ? loopLabel : "?", slots[0].ci);
        return false;
    }
    if (n == 2 && bothEndsSharedTShape(slots[0].edge, slots[1].edge)) {
        out.ok = true;
        out.order = {0, 1};
        out.forward = {1, 1};
        return true;
    }
    std::vector<char> visited(n, 0);
    std::vector<size_t> ord;
    std::vector<char> fwd;
    visited[0] = 1;
    ord.push_back(0);
    fwd.push_back(1);
    TopoDS_Vertex freeVtx = TopExp::LastVertex(seed, Standard_True);
    while (ord.size() < n) {
        std::vector<std::pair<size_t, double>> tsHits, geoNear;
        collectJunctionCandidates(freeVtx, slots, visited, ord.back(), tsHits, geoNear);
        if (tsHits.size() != 1) {
            out.notDetermined = true;
            out.failWireIx = (int)ord.back();
            std::fprintf(stderr,
                         "DIAG_OP2G_JUNCTION rid=%d loop=%s ci=%d candidates=%zu dists=",
                         rid, loopLabel ? loopLabel : "?", slots[ord.back()].ci, tsHits.size());
            if (tsHits.empty() && !geoNear.empty()) {
                for (size_t k = 0; k < geoNear.size(); k++)
                    std::fprintf(stderr, "%s%.9f", k ? "," : "", geoNear[k].second);
            } else {
                for (size_t k = 0; k < tsHits.size(); k++)
                    std::fprintf(stderr, "%s0", k ? "," : "");
            }
            std::fprintf(stderr, "\n");
            return false;
        }
        const size_t nextIx = tsHits[0].first;
        visited[nextIx] = 1;
        ord.push_back(nextIx);
        const TopoDS_Edge& ne = slots[nextIx].edge;
        const TopoDS_Vertex nf = TopExp::FirstVertex(ne, Standard_True);
        const bool nextFwd = sameVtxTShape(nf, freeVtx);
        fwd.push_back(nextFwd ? 1 : 0);
        TopoDS_Edge neDir = nextFwd ? ne : TopoDS::Edge(ne.Reversed());
        freeVtx = TopExp::LastVertex(neDir, Standard_True);
    }
    out.ok = true;
    out.order = std::move(ord);
    out.forward = std::move(fwd);
    return true;
}

bool buildLoopWireSlots(const Loop& loop, const RegionSet& rs, const MeshView& mv,
                        const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                        const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                        std::vector<LoopWireSlot>& slots) {
    slots.clear();
    size_t wireIx = 0;
    for (size_t chainPos = 0; chainPos < loop.chainIdx.size(); chainPos++) {
        const int ci = loop.chainIdx[chainPos];
        const size_t nE = chainWireEdgeCount(ci, rs, geom, collapsed);
        for (size_t edgeIx = 0; edgeIx < nE; edgeIx++) {
            LoopWireSlot s;
            s.wireIx = wireIx++;
            s.ci = ci;
            if (!loopChainHandedEdge(chainPos, edgeIx, loop, rs, mv, geom, collapsed, meshE, edgeOk,
                                     s.edge))
                return false;
            slots.push_back(std::move(s));
        }
    }
    return !slots.empty();
}

TopoDS_Edge directedSlotEdge(const LoopWireSlot& s, char forward) {
    return forward ? s.edge : TopoDS::Edge(s.edge.Reversed());
}

// D-S3-168 step 0: signed area of the walk's cyclic vertex sequence in plate UV.
double walkSignedAreaUV(const std::vector<LoopWireSlot>& slots, const LoopConnWalk& conn,
                        const Region& plate, size_t* nUvOut) {
    if (nUvOut) *nUvOut = 0;
    if (!conn.ok || conn.order.size() != conn.forward.size() || conn.order.empty()) return 0.0;
    const gp_Pln pln(plate.ax);
    std::vector<gp_Pnt2d> uv;
    uv.reserve(conn.order.size() + 1);
    for (size_t k = 0; k < conn.order.size(); k++) {
        const size_t ix = conn.order[k];
        if (ix >= slots.size()) continue;
        const TopoDS_Edge e = directedSlotEdge(slots[ix], conn.forward[k]);
        const TopoDS_Vertex v = TopExp::FirstVertex(e, Standard_True);
        if (v.IsNull()) continue;
        Standard_Real u = 0, vv = 0;
        ElSLib::Parameters(pln, BRep_Tool::Pnt(v), u, vv);
        const gp_Pnt2d q(u, vv);
        if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
        uv.push_back(q);
    }
    if (uv.size() >= 2 && uv.front().Distance(uv.back()) <= Precision::PConfusion()) uv.pop_back();
    if (nUvOut) *nUvOut = uv.size();
    if (uv.size() < 3) return 0.0;
    double a = 0.0;
    for (size_t i = 0; i < uv.size(); i++) {
        const gp_Pnt2d& p0 = uv[i];
        const gp_Pnt2d& p1 = uv[(i + 1) % uv.size()];
        a += p0.X() * p1.Y() - p1.X() * p0.Y();
    }
    return 0.5 * a;
}

// D-S3-171 case (c): signed area of one edge from the curve, not a 1-vertex polygon.
// Circle: sign(axis · plate normal) × traversal orientation × π r².
// Else: dense samples along the adaptor's own parameter range.
double edgeCurveSignedAreaUV(const TopoDS_Edge& eIn, const Region& plate) {
    if (eIn.IsNull()) return 0.0;
    const gp_Pln pln(plate.ax);
    const gp_Dir n = pln.Axis().Direction();
    try {
        BRepAdaptor_Curve ac(eIn);
        if (ac.GetType() == GeomAbs_Circle) {
            const gp_Circ circ = ac.Circle();
            double s = circ.Axis().Direction().Dot(n);
            if (eIn.Orientation() == TopAbs_REVERSED) s = -s;
            const double r = circ.Radius();
            return s * kPi * r * r;
        }
        const Standard_Real f = ac.FirstParameter();
        const Standard_Real l = ac.LastParameter();
        if (l - f <= Precision::PConfusion()) return 0.0;
        std::vector<gp_Pnt2d> uv;
        uv.reserve(32);
        for (int i = 0; i < 32; i++) {
            const Standard_Real t =
                f + (l - f) * (Standard_Real)i / 32.0;
            Standard_Real u = 0, vv = 0;
            ElSLib::Parameters(pln, ac.Value(t), u, vv);
            const gp_Pnt2d q(u, vv);
            if (!uv.empty() && uv.back().Distance(q) <= Precision::PConfusion()) continue;
            uv.push_back(q);
        }
        if (uv.size() < 3) return 0.0;
        double a = 0.0;
        for (size_t i = 0; i < uv.size(); i++) {
            const gp_Pnt2d& p0 = uv[i];
            const gp_Pnt2d& p1 = uv[(i + 1) % uv.size()];
            a += p0.X() * p1.Y() - p1.X() * p0.Y();
        }
        return 0.5 * a;
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

double walkCurveSignedAreaUV(const std::vector<LoopWireSlot>& slots, const LoopConnWalk& conn,
                             const Region& plate) {
    if (conn.ok && conn.order.size() == conn.forward.size() && !conn.order.empty()) {
        double a = 0.0;
        for (size_t k = 0; k < conn.order.size(); k++) {
            const size_t ix = conn.order[k];
            if (ix >= slots.size()) continue;
            a += edgeCurveSignedAreaUV(directedSlotEdge(slots[ix], conn.forward[k]), plate);
        }
        return a;
    }
    if (slots.size() == 1) return edgeCurveSignedAreaUV(slots[0].edge, plate);
    return 0.0;
}

double wireCurveSignedAreaUV(const TopoDS_Wire& w, const Region& plate) {
    if (w.IsNull()) return 0.0;
    double a = 0.0;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next())
            a += edgeCurveSignedAreaUV(TopoDS::Edge(ex.Current()), plate);
    } catch (const Standard_Failure&) {
        return 0.0;
    }
    return a;
}

// D-S3-168 step 0.5: reverse the cyclic walk; keep chainIdx[0] as start.
void reverseConnTraversal(LoopConnWalk& conn, const std::vector<LoopWireSlot>& slots) {
    const size_t n = conn.order.size();
    if (n < 1 || n != conn.forward.size() || conn.order[0] >= slots.size()) return;
    if (n == 1) {
        conn.forward[0] = conn.forward[0] ? 0 : 1;
        return;
    }
    if (n == 2 && bothEndsSharedTShape(slots[conn.order[0]].edge, slots[conn.order[1]].edge)) {
        for (size_t k = 0; k < n; k++)
            conn.forward[k] = conn.forward[k] ? 0 : 1;
        return;
    }
    std::reverse(conn.order.begin() + 1, conn.order.end());
    conn.forward.assign(n, 1);
    const TopoDS_Edge& seed = slots[conn.order[0]].edge;
    const TopoDS_Edge& succ = slots[conn.order[1]].edge;
    const TopoDS_Vertex sF = TopExp::FirstVertex(seed, Standard_True);
    const TopoDS_Vertex sL = TopExp::LastVertex(seed, Standard_True);
    const TopoDS_Vertex nF = TopExp::FirstVertex(succ, Standard_True);
    const TopoDS_Vertex nL = TopExp::LastVertex(succ, Standard_True);
    const bool headAtSucc = sameVtxTShape(sL, nF) || sameVtxTShape(sL, nL);
    const bool tailAtSucc = sameVtxTShape(sF, nF) || sameVtxTShape(sF, nL);
    conn.forward[0] = (headAtSucc && !tailAtSucc) ? 1 : 0;
    TopoDS_Edge cur = conn.forward[0] ? seed : TopoDS::Edge(seed.Reversed());
    TopoDS_Vertex freeVtx = TopExp::LastVertex(cur, Standard_True);
    for (size_t k = 1; k < n; k++) {
        const TopoDS_Edge& e = slots[conn.order[k]].edge;
        const TopoDS_Vertex f = TopExp::FirstVertex(e, Standard_True);
        conn.forward[k] = sameVtxTShape(f, freeVtx) ? 1 : 0;
        TopoDS_Edge ed = conn.forward[k] ? e : TopoDS::Edge(e.Reversed());
        freeVtx = TopExp::LastVertex(ed, Standard_True);
    }
}

bool walkOrderMatchesAppend(const LoopConnWalk& walk) {
    if (walk.order.size() != walk.forward.size()) return false;
    for (size_t i = 0; i < walk.order.size(); i++)
        if (walk.order[i] != i) return false;
    return true;
}

bool collectWireEdges(const TopoDS_Wire& w, std::vector<TopoDS_Edge>& edges);

void registerShippedLoopOrder(int rid, const char* loopLabel, const TopoDS_Wire& w,
                              const LoopConnWalk& walk, const std::vector<LoopWireSlot>& slots) {
    if (!walk.ok || walk.order.empty()) return;
    std::vector<TopoDS_Edge> shipped;
    if (!collectWireEdges(w, shipped)) return;
    const size_t n = walk.order.size();
    if (shipped.size() != n) {
        std::fprintf(stderr,
                     "REG_OP2G_LOOPORDER rid=%d loop=%s shipped=%zu walk=%zu mismatch=count\n", rid,
                     loopLabel ? loopLabel : "?", shipped.size(), n);
        return;
    }
    bool same = true;
    for (size_t i = 0; i < n; i++) {
        const TopoDS_Edge want = directedSlotEdge(slots[walk.order[i]], walk.forward[i]);
        const void* a = diagTShapePtr(shipped[i]);
        const void* b = diagTShapePtr(want);
        if (a != b) {
            same = false;
            break;
        }
    }
    std::fprintf(stderr, "REG_OP2G_LOOPORDER rid=%d loop=%s shippedVsWalk=%s\n", rid,
                 loopLabel ? loopLabel : "?", same ? "same" : "reordered");
}

const char* loopRoleLabel(LoopRole role) {
    switch (role) {
    case LoopRole::Outer:
        return "outer";
    case LoopRole::Inner:
        return "inner";
    case LoopRole::CapLow:
        return "capLow";
    case LoopRole::CapHigh:
        return "capHigh";
    }
    return "?";
}

double minEndpointDist3d(const TopoDS_Vertex& a0, const TopoDS_Vertex& a1,
                          const TopoDS_Vertex& b0, const TopoDS_Vertex& b1) {
    auto d = [](const TopoDS_Vertex& a, const TopoDS_Vertex& b) -> double {
        if (a.IsNull() || b.IsNull()) return 1e300;
        return BRep_Tool::Pnt(a).Distance(BRep_Tool::Pnt(b));
    };
    return std::min(std::min(d(a0, b0), d(a0, b1)), std::min(d(a1, b0), d(a1, b1)));
}

// D-S3-147 / D-S3-159: orient seed so its head (Last) is the single shared junction
// TShape with the connectivity successor; both/two shared -> area-only (b).
bool orientSeedBySuccessor(TopoDS_Edge& seed, const TopoDS_Edge& succ, int seedCi, int succCi,
                           WireAppendOri& st, int rid, const char* loopLabel) {
    const TopoDS_Vertex sF = TopExp::FirstVertex(seed, Standard_True);
    const TopoDS_Vertex sL = TopExp::LastVertex(seed, Standard_True);
    const TopoDS_Vertex nF = TopExp::FirstVertex(succ, Standard_True);
    const TopoDS_Vertex nL = TopExp::LastVertex(succ, Standard_True);
    if (bothEndsSharedTShape(seed, succ)) {
        st.seedOrientedBy = "area-only";
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2E_SEED seedCi=%d succCi=%d shared=both seedOrientedBy=area-only "
                         "seedTailT=%p seedHeadT=%p succTailT=%p succHeadT=%p\n",
                         seedCi, succCi, diagTShapePtr(sF), diagTShapePtr(sL), diagTShapePtr(nF),
                         diagTShapePtr(nL));
        return true;
    }
    const void* juncTs = nullptr;
    TopoDS_Vertex junc;
    int nSharedTs = 0;
    auto noteShare = [&](const TopoDS_Vertex& a, const TopoDS_Vertex& b) {
        if (!sameVtxTShape(a, b)) return;
        const void* ts = diagTShapePtr(a);
        if (!ts) return;
        if (juncTs && ts != juncTs)
            nSharedTs = 2;
        else {
            juncTs = ts;
            junc = a;
            nSharedTs = 1;
        }
    };
    noteShare(sF, nF);
    noteShare(sF, nL);
    noteShare(sL, nF);
    noteShare(sL, nL);
    if (nSharedTs == 2) {
        st.seedOrientedBy = "area-only";
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2E_SEED seedCi=%d succCi=%d shared=two seedOrientedBy=area-only "
                         "seedTailT=%p seedHeadT=%p succTailT=%p succHeadT=%p\n",
                         seedCi, succCi, diagTShapePtr(sF), diagTShapePtr(sL), diagTShapePtr(nF),
                         diagTShapePtr(nL));
        return true;
    }
    if (nSharedTs == 0) {
        const double dist = minEndpointDist3d(sF, sL, nF, nL);
        std::fprintf(stderr,
                     "DIAG_OP2G_JUNCTION rid=%d loop=%s ci=%d candidates=0 dists=%.9f\n", rid,
                     loopLabel ? loopLabel : "?", seedCi, dist);
        st.notDetermined++;
        gOriNotDetermined++;
        st.seedOrientedBy = "notDetermined";
        return false;
    }
    const bool headAtJ = sameVtxTShape(sL, junc);
    const bool tailAtJ = sameVtxTShape(sF, junc);
    if (headAtJ && !tailAtJ) {
        st.seedOrientedBy = "successor";
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2E_SEED seedCi=%d succCi=%d shared=one junctionT=%p "
                         "seedOrientedBy=successor seedReversed=0\n",
                         seedCi, succCi, juncTs);
        return true;
    }
    if (tailAtJ && !headAtJ) {
        seed = TopoDS::Edge(seed.Reversed());
        st.reversedByRule++;
        st.seedReversed = 1;
        st.seedOrientedBy = "successor";
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2E_SEED seedCi=%d succCi=%d shared=one junctionT=%p "
                         "seedOrientedBy=successor seedReversed=1\n",
                         seedCi, succCi, juncTs);
        return true;
    }
    st.notDetermined++;
    gOriNotDetermined++;
    st.seedOrientedBy = "notDetermined";
    std::fprintf(stderr,
                 "DIAG_WIREORI2_NOTDET site=seed-successor seedCi=%d succCi=%d shared=one "
                 "ambiguous headAt=%d tailAt=%d junctionT=%p\n",
                 seedCi, succCi, headAtJ ? 1 : 0, tailAtJ ? 1 : 0, juncTs);
    return false;
}

// D-S3-117/118/119: one rule at append. Local handle only (N4). Seed = loop start

TopoDS_Vertex loopStartVertex(const Loop& loop, const RegionSet& rs, const MeshView& mv,
                              const std::vector<TopoDS_Edge>& meshE) {
    const int startLv = loopStartMeshVert(loop, rs);
    if (startLv < 0 || loop.chainIdx.empty()) return {};
    const int ci = loop.chainIdx[0];
    if (ci < 0 || (size_t)ci >= rs.chains.size()) return {};
    const BoundaryChain& ch = rs.chains[(size_t)ci];
    for (int eid : ch.meshEdges) {
        if (eid < 0 || (size_t)eid >= mv.nEdge || (size_t)eid >= meshE.size()) continue;
        if (meshE[(size_t)eid].IsNull()) continue;
        const auto& ev = mv.compEdges[eid];
        if (ev.first != startLv && ev.second != startLv) continue;
        const bool fwd = (ev.first == startLv);
        const TopoDS_Edge me =
            fwd ? meshE[(size_t)eid] : TopoDS::Edge(meshE[(size_t)eid].Reversed());
        return TopExp::FirstVertex(me, Standard_True);
    }
    return {};
}

// D-S3-147/117/159: seed preset or as-built fallback; later edges by TShape continuity.
TopoDS_Edge orientAppendedEdge(const TopoDS_Edge& eIn, WireAppendOri& st, int ci,
                               const char* site) {
    TopoDS_Edge e = eIn;
    if (e.IsNull()) return e;
    if (st.useSlotOri && st.slotConsume < st.slotOri.size()) {
        e = st.slotOri[st.slotConsume++];
        st.havePrev = true;
        st.prev = e;
        return e;
    }
    if (st.wholeLoopAsBuilt) {
        st.havePrev = true;
        st.prev = e;
        return e;
    }
    if (!st.havePrev) {
        if (st.seedPreset) {
            e = st.seedPresetEdge;
            st.seedPreset = false;
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_OP2D_SEED_EDGE ci=%d seedOrientedBy=%s seedReversed=%d\n", ci,
                             st.seedOrientedBy ? st.seedOrientedBy : "?", st.seedReversed);
            st.havePrev = true;
            st.prev = e;
            return e;
        }
        st.notDetermined++;
        gOriNotDetermined++;
        st.seedOrientedBy = "notDetermined";
        std::fprintf(stderr, "DIAG_WIREORI2_NOTDET site=seed-missing ci=%d\n", ci);
        st.havePrev = true;
        st.prev = e;
        return e;
    }
    const TopoDS_Vertex prevHead = TopExp::LastVertex(st.prev, Standard_True);
    const TopoDS_Vertex first = TopExp::FirstVertex(e, Standard_True);
    const TopoDS_Vertex last = TopExp::LastVertex(e, Standard_True);
    if (sameVtxTShape(first, prevHead)) {
        // continuous
    } else if (sameVtxTShape(last, prevHead)) {
        e = TopoDS::Edge(e.Reversed());
        st.reversedByRule++;
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2G_REV ci=%d site=%s tailT=%p headT=%p prevHeadT=%p\n", ci,
                         site ? site : "continuity", diagTShapePtr(first), diagTShapePtr(last),
                         diagTShapePtr(prevHead));
    } else {
        st.notDetermined++;
        gOriNotDetermined++;
        std::fprintf(stderr,
                     "DIAG_WIREORI2_NOTDET site=continuity ci=%d neither-vertex-matches\n", ci);
    }
    st.prev = e;
    return e;
}

bool appendPolyline(BRep_Builder& B, TopoDS_Wire& w, const MeshView& mv, const BoundaryChain& ch,
                    bool reversed, const std::vector<TopoDS_Edge>& meshE,
                    const std::vector<char>& edgeOk, const Region* plateR, int ci,
                    WireAppendOri& ori) {
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
        noteWirePop(meshE[(size_t)eid], "appendPolyline", ci);
        bool boundBefore = false;
        if (plateR && plateR->type == SurfType::Plane && isAnalytic(plateR)) {
            const void* ts = diagTShapePtr(meshE[(size_t)eid]);
            if (!ts || !gSeamTShapes.count(ts)) {
                bindEdgePcurveOnInternedPlane(meshE[(size_t)eid], *plateR, mv, "poly", ci);
                boundBefore = true;
            }
        }
        TopoDS_Edge native = meshE[(size_t)eid];
        TopoDS_Edge addE = fwd ? native : TopoDS::Edge(native.Reversed());
        const TopoDS_Vertex contVtx =
            ori.havePrev ? TopExp::LastVertex(ori.prev, Standard_True) : TopoDS_Vertex();
        const int revBefore = ori.reversedByRule;
        addE = orientAppendedEdge(addE, ori, ci, "appendPolyline");
        dumpOp2bBind(addE, native, ci, "appendPolyline", ori.reversedByRule - revBefore, 0,
                     boundBefore ? "interned-before" : "none", plateR, contVtx);
        B.Add(w, addE);
    }
    return true;
}

bool appendCollapsed(BRep_Builder& B, TopoDS_Wire& w, const ChainGeom& g, bool reversed, int ci,
                     WireAppendOri& ori, const Region* plateR) {
    if (g.edges.empty()) return false;
    auto addOne = [&](const TopoDS_Edge& native, const TopoDS_Edge& handed) {
        diagNoteEdge(handed, "appendCollapsed");
        noteWirePop(handed, "appendCollapsed", ci);
        const TopoDS_Vertex contVtx =
            ori.havePrev ? TopExp::LastVertex(ori.prev, Standard_True) : TopoDS_Vertex();
        const int revBefore = ori.reversedByRule;
        TopoDS_Edge addE = orientAppendedEdge(handed, ori, ci, "appendCollapsed");
        dumpOp2bBind(addE, native, ci, "appendCollapsed", ori.reversedByRule - revBefore, 0,
                     "bindAllVariants-before", plateR, contVtx);
        B.Add(w, addE);
    };
    if (!reversed) {
        for (const auto& e : g.edges) addOne(e, e);
    } else {
        for (int i = (int)g.edges.size() - 1; i >= 0; i--) {
            TopoDS_Edge native = g.edges[(size_t)i];
            addOne(native, TopoDS::Edge(native.Reversed()));
        }
    }
    return true;
}

bool buildLoopWire(TopoDS_Wire& w, const Loop& loop, const RegionSet& rs, const MeshView& mv,
                   const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                   const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                   const Region* plateR) {
    BRep_Builder B;
    B.MakeWire(w);
    if (loop.chainIdx.size() != loop.reversed.size()) return false;
    Loop walk = loop;
    int flipped = 0;
    double area = 0.0;
    double vtxRes = 0.0;
    double areaFloor = 0.0;
    int wanted = 0;
    if (plateR && plateR->type == SurfType::Plane) {
        area = loopSignedAreaUV(walk, rs, mv, *plateR);
        vtxRes = loopUvMinEdgeLen(walk, rs, mv, *plateR);
        areaFloor = Precision::SquareConfusion();
        const bool degen = std::fabs(area) <= areaFloor;
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2C_SENSE rid=%d role=%d area=%.6g absArea=%.6g vtxRes=%.6g "
                         "areaFloor=%.6g degen=%d\n",
                         plateR->id, (int)walk.role, area, std::fabs(area), vtxRes, areaFloor,
                         degen ? 1 : 0);
        if (degen) {
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_OP2B_DEGEN rid=%d role=%d area=%.6g vtxRes=%.6g -- A4 not built\n",
                             plateR->id, (int)walk.role, area, vtxRes);
            return false;
        }
        const gp_Dir z = plateR->ax.Direction();
        const gp_Dir outward = plateR->outwardNormal ? z : gp_Dir(z.Reversed());
        const int signZ = (z.Dot(outward) >= 0.0) ? 1 : -1;
        // Outer CCW about outward => area*signZ > 0; inner opposite.
        const bool inner = (walk.role == LoopRole::Inner);
        wanted = inner ? -signZ : signZ;
        if (diagP2Enabled())
            std::fprintf(stderr,
                         "DIAG_OP2B_SENSE rid=%d role=%d area=%.6g wanted=%d flipped=%d\n",
                         plateR->id, (int)walk.role, area, wanted, flipped);
    }
    WireAppendOri ori;
    const int ci0 = walk.chainIdx.empty() ? -1 : walk.chainIdx[0];
    const int rev0 = walk.chainIdx.empty() ? 0 : (int)walk.reversed[0];
    const int startLv = loopStartMeshVert(walk, rs);
    const size_t nWireE = loopWireEdgeCount(walk, rs, geom, collapsed);
    const char* loopLbl = loopRoleLabel(walk.role);
    const int rid = plateR ? plateR->id : -1;
    if (nWireE == 0) return false;
    std::vector<LoopWireSlot> slots;
    if (!buildLoopWireSlots(walk, rs, mv, geom, collapsed, meshE, edgeOk, slots)) return false;
    LoopConnWalk conn;
    bool useAsBuilt = false;
    int succCi = -1;
    int loopNotDet = 0;
    const char* walkCmp = "same";
    int travReversed = 0;
    bool senseFromCurve = false;
    size_t senseNVtx = 0;
    if (nWireE == 1) {
        // D-S3-171 case (c): single closed edge → area-only from the EDGE curve.
        ori.seedOrientedBy = "area-only";
        TopoDS_Edge seedE = slots[0].edge;
        if (plateR && plateR->type == SurfType::Plane && wanted != 0) {
            const double curveArea = edgeCurveSignedAreaUV(seedE, *plateR);
            area = curveArea;
            senseFromCurve = true;
            senseNVtx = 1;
            std::fprintf(stderr,
                         "DIAG_OP2K_SENSE rid=%d loop=%s nVtx=1 areaSource=curve area=%.6f "
                         "wanted=%+d\n",
                         rid, loopLbl, curveArea, wanted);
            if ((curveArea > 0.0) != (wanted > 0)) {
                seedE = TopoDS::Edge(seedE.Reversed());
                ori.seedReversed = 1;
                ori.reversedByRule++;
                travReversed = 1;
            }
        }
        ori.seedPresetEdge = seedE;
        ori.seedPreset = true;
    } else {
        if (!loopConnectivityWalk(slots, rid, loopLbl, conn)) {
            useAsBuilt = true;
            loopNotDet = 1;
            ori.wholeLoopAsBuilt = true;
            ori.seedOrientedBy = "as-built";
            gOriNotDetermined++;
        } else {
            walkCmp = walkOrderMatchesAppend(conn) ? "same" : "reordered";
            // D-S3-168 step 0.5: reverse the traversal if the walk polygon's
            // sign disagrees with the loop role. |area|<=SquareConfusion → A4
            // only when the sequence is a real polygon (>=3 vertices).
            if (plateR && plateR->type == SurfType::Plane && wanted != 0) {
                size_t nUv = 0;
                const double seqArea = walkSignedAreaUV(slots, conn, *plateR, &nUv);
                if (nUv >= 3) {
                    const double areaFloorWalk = Precision::SquareConfusion();
                    if (std::fabs(seqArea) <= areaFloorWalk) {
                        if (diagP2Enabled())
                            std::fprintf(stderr,
                                         "DIAG_OP2B_DEGEN rid=%d role=%d area=%.6g vtxRes=%.6g "
                                         "-- A4 not built (walk 0.5)\n",
                                         plateR->id, (int)walk.role, seqArea, vtxRes);
                        return false;
                    }
                    area = seqArea;
                    senseNVtx = nUv;
                    std::fprintf(stderr,
                                 "DIAG_OP2K_SENSE rid=%d loop=%s nVtx=%zu areaSource=polygon "
                                 "area=%.6f wanted=%+d\n",
                                 rid, loopLbl, nUv, seqArea, wanted);
                    if ((seqArea > 0.0) != (wanted > 0)) {
                        reverseConnTraversal(conn, slots);
                        travReversed = 1;
                    }
                } else {
                    const double curveArea = walkCurveSignedAreaUV(slots, conn, *plateR);
                    area = curveArea;
                    senseFromCurve = true;
                    senseNVtx = nUv;
                    std::fprintf(stderr,
                                 "DIAG_OP2K_SENSE rid=%d loop=%s nVtx=%zu areaSource=curve "
                                 "area=%.6f wanted=%+d\n",
                                 rid, loopLbl, nUv, curveArea, wanted);
                    if ((curveArea > 0.0) != (wanted > 0)) {
                        reverseConnTraversal(conn, slots);
                        travReversed = 1;
                    }
                }
            }
            walkCmp = walkOrderMatchesAppend(conn) ? "same" : "reordered";
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_OP2J_TRAV rid=%d loop=%s reversed=%d area=%.6g wanted=%d\n",
                             rid, loopLbl, travReversed, area, wanted);
            const size_t succWireIx = conn.order.size() > 1 ? conn.order[1] : 0;
            TopoDS_Edge seedE = slots[0].edge;
            TopoDS_Edge succE = slots[succWireIx].edge;
            succCi = slots[succWireIx].ci;
            if (!orientSeedBySuccessor(seedE, succE, ci0, succCi, ori, rid, loopLbl)) {
                useAsBuilt = true;
                loopNotDet = 1;
                ori = {};
                ori.wholeLoopAsBuilt = true;
                ori.seedOrientedBy = "as-built";
            } else {
                ori.seedPresetEdge = seedE;
                ori.seedPreset = true;
                const bool areaOnly =
                    strcmp(ori.seedOrientedBy ? ori.seedOrientedBy : "", "area-only") == 0;
                if (areaOnly && conn.ok && conn.order.size() == slots.size()) {
                    // Two-edge both-shared: apply the (possibly reversed) walk
                    // forwards — same outcome as area-only, via the sequence.
                    ori.slotOri.assign(slots.size(), TopoDS_Edge());
                    for (size_t k = 0; k < conn.order.size(); k++) {
                        const size_t ix = conn.order[k];
                        TopoDS_Edge e = directedSlotEdge(slots[ix], conn.forward[k]);
                        if (!conn.forward[k]) {
                            ori.reversedByRule++;
                            if (k == 0) ori.seedReversed = 1;
                        }
                        ori.slotOri[ix] = e;
                        if (k == 0) {
                            ori.seedPresetEdge = e;
                            seedE = e;
                        }
                    }
                    ori.useSlotOri = true;
                } else if (conn.ok && conn.order.size() == slots.size()) {
                    ori.slotOri.assign(slots.size(), TopoDS_Edge());
                    ori.slotOri[conn.order[0]] = seedE;
                    TopoDS_Edge prevE = seedE;
                    for (size_t k = 1; k < conn.order.size(); k++) {
                        const size_t ix = conn.order[k];
                        TopoDS_Edge e = slots[ix].edge;
                        const TopoDS_Vertex prevHead = TopExp::LastVertex(prevE, Standard_True);
                        const TopoDS_Vertex first = TopExp::FirstVertex(e, Standard_True);
                        const TopoDS_Vertex last = TopExp::LastVertex(e, Standard_True);
                        if (sameVtxTShape(first, prevHead)) {
                            // already continuous
                        } else if (sameVtxTShape(last, prevHead)) {
                            e = TopoDS::Edge(e.Reversed());
                            ori.reversedByRule++;
                            if (diagP2Enabled())
                                std::fprintf(stderr,
                                             "DIAG_OP2G_REV ci=%d site=walk-pred tailT=%p "
                                             "headT=%p prevHeadT=%p\n",
                                             slots[ix].ci, diagTShapePtr(first),
                                             diagTShapePtr(last), diagTShapePtr(prevHead));
                        } else {
                            ori.notDetermined++;
                            gOriNotDetermined++;
                            std::fprintf(stderr,
                                         "DIAG_WIREORI2_NOTDET site=walk-pred ci=%d "
                                         "neither-vertex-matches\n",
                                         slots[ix].ci);
                        }
                        ori.slotOri[ix] = e;
                        prevE = e;
                    }
                    ori.useSlotOri = true;
                }
            }
        }
    }
    if (diagP2Enabled())
        std::fprintf(stderr,
                     "DIAG_OP2D_SEED_SETUP rid=%d ci0=%d nWireE=%zu seedOrientedBy=%s "
                     "seedReversed=%d asBuilt=%d walk=%s\n",
                     rid, ci0, nWireE, ori.seedOrientedBy ? ori.seedOrientedBy : "?",
                     ori.seedReversed, useAsBuilt ? 1 : 0, walkCmp);
    for (size_t i = 0; i < walk.chainIdx.size(); i++) {
        int ci = walk.chainIdx[i];
        if (ci < 0 || (size_t)ci >= rs.chains.size()) return false;
        bool rev = walk.reversed[i] != 0;
        bool ok;
        if ((size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
            (size_t)ci < geom.size() && geom[(size_t)ci].collapsed && !geom[(size_t)ci].edges.empty())
            ok = appendCollapsed(B, w, geom[(size_t)ci], rev, (int)ci, ori, plateR);
        else
            ok = appendPolyline(B, w, mv, rs.chains[(size_t)ci], rev, meshE, edgeOk, plateR,
                                (int)ci, ori);
        if (!ok) return false;
    }
    // D-S3-168 step 3: wire-level sense is an assertion. On a well-built loop
    // it is a no-op; print (do not reverse) when it is not.
    if (!useAsBuilt && plateR && plateR->type == SurfType::Plane && wanted != 0) {
        const double warea =
            senseFromCurve ? wireCurveSignedAreaUV(w, *plateR) : wireSignedAreaUV(w, *plateR);
        area = warea;
        if ((warea > 0.0) != (wanted > 0)) {
            flipped = 1;
            std::fprintf(stderr,
                         "DIAG_OP2J_ASSERT_SENSE rid=%d loop=%s wireArea=%.6g wanted=%d "
                         "travReversed=%d -- not a no-op\n",
                         rid, loopRoleLabel(walk.role), warea, wanted, travReversed);
        }
    }
    const char* loopLblDone = loopRoleLabel(walk.role);
    const int continuityRev = ori.reversedByRule - ori.seedReversed;
    const int loopNotDetTotal = loopNotDet + ori.notDetermined;
    std::fprintf(stderr,
                 "DIAG_OP2G_LOOP rid=%d loop=%s nE=%zu walkOrderFromAppend=%s seedCi=%d "
                 "succCi(connectivity)=%d seedOrientedBy=%s seedReversed=%d continuityRev=%d "
                 "flipped=%d notDetermined=%d area=%.6g\n",
                 rid, loopLblDone, nWireE, walkCmp, ci0, succCi,
                 ori.seedOrientedBy ? ori.seedOrientedBy : "?", ori.seedReversed, continuityRev,
                 flipped, loopNotDetTotal, area);
    if (conn.ok && plateR)
        registerShippedLoopOrder(rid, loopLblDone, w, conn, slots);
    if (rid == 2 && walk.role == LoopRole::Outer && diagP2Enabled()) {
        std::fprintf(stderr,
                     "DIAG_OP2G_HUB2_OUTER walkOrderFromAppend=%s seedCi=%d seedReversed=%d "
                     "reversedByRule=%d continuityRev=%d (seed odd-one-out if continuityRev>>seed)\n",
                     walkCmp, ci0, ori.seedReversed, ori.reversedByRule, continuityRev);
    }
    if (diagP2Enabled() && plateR) {
        std::fprintf(stderr,
                     "DIAG_OP2D_SEED rid=%d wireIx=%zu ci0=%d rev0=%d startLv=%d "
                     "seedOrientedBy=%s seedReversed=%d reversedByRule=%d notDetermined=%d "
                     "flipped=%d area=%.6g\n",
                     plateR->id, gPlateWireOri.size(), ci0, rev0, startLv,
                     ori.seedOrientedBy ? ori.seedOrientedBy : "?", ori.seedReversed,
                     ori.reversedByRule, loopNotDetTotal, flipped, area);
    }
    gPlateWireOri.push_back({ori.reversedByRule, loopNotDetTotal, ori.seedReversed, flipped, ci0,
                             ori.seedOrientedBy ? ori.seedOrientedBy : "none", loopLblDone});
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

double pcurveSignedArea(const TopoDS_Face& f, const TopoDS_Wire& w);

void dumpBRepStatuses(const Handle(BRepCheck_Result)& res) {
    if (res.IsNull()) return;
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next()) {
        const int st = (int)it.Value();
        if (st == (int)BRepCheck_NoError) continue;
        std::fprintf(stderr, " st=%d(%s)", st, brepCheckName(st));
    }
}

void dumpBRepStatusList(const BRepCheck_ListOfStatus& lst) {
    for (BRepCheck_ListOfStatus::Iterator it(lst); it.More(); it.Next()) {
        const int st = (int)it.Value();
        if (st == (int)BRepCheck_NoError) continue;
        std::fprintf(stderr, " st=%d(%s)", st, brepCheckName(st));
    }
}

bool diagWireEnabled() { return diagP2Enabled(); }

const char* topAbsStateName(TopAbs_State s) {
    switch (s) {
        case TopAbs_IN: return "IN";
        case TopAbs_OUT: return "OUT";
        case TopAbs_ON: return "ON";
        default: return "?";
    }
}

bool segIntersect2d(const gp_Pnt2d& a0, const gp_Pnt2d& a1, const gp_Pnt2d& b0, const gp_Pnt2d& b1,
                    double tol) {
    const double dax = a1.X() - a0.X(), day = a1.Y() - a0.Y();
    const double dbx = b1.X() - b0.X(), dby = b1.Y() - b0.Y();
    const double denom = dax * dby - day * dbx;
    if (std::fabs(denom) < tol * tol) return false;
    const double t = ((b0.X() - a0.X()) * dby - (b0.Y() - a0.Y()) * dbx) / denom;
    const double u = ((b0.X() - a0.X()) * day - (b0.Y() - a0.Y()) * dax) / denom;
    if (t < -tol || t > 1.0 + tol || u < -tol || u > 1.0 + tol) return false;
    if (a0.Distance(b0) <= tol || a0.Distance(b1) <= tol || a1.Distance(b0) <= tol ||
        a1.Distance(b1) <= tol)
        return false;
    return true;
}

void wireDiagMetrics(const TopoDS_Face& f, const TopoDS_Wire& w, double sewTol, double& areaUV,
                     int& nSelfX, double& maxUVgap, int& nDup, double* maxUVgapIterOut = nullptr,
                     double* maxGap3dOut = nullptr) {
    areaUV = pcurveSignedArea(f, w);
    nSelfX = 0;
    maxUVgap = 0.0;
    nDup = 0;
    double maxUVgapIter = 0.0;
    double maxGap3d = 0.0;
    const double tol = std::max(sewTol, Precision::PConfusion());
    struct Seg2d {
        gp_Pnt2d a, b;
        int ei;
    };
    std::vector<Seg2d> segs;
    std::vector<gp_Pnt> mids;
    const int nSample = 12;
    try {
        struct EdgeEnds {
            gp_Pnt2d uv0, uv1;
            gp_Pnt p0, p1;
            bool havePc = false;
        };
        std::vector<EdgeEnds> expEnds;
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
            const TopoDS_Edge e = ex.Current();
            EdgeEnds ee;
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
            if (!pc.IsNull()) {
                ee.havePc = true;
                if (e.Orientation() == TopAbs_REVERSED) std::swap(a, b);
                ee.uv0 = pc->Value(a);
                ee.uv1 = pc->Value(b);
            }
            TopoDS_Vertex v1, v2;
            TopExp::Vertices(e, v1, v2, Standard_True);
            ee.p0 = BRep_Tool::Pnt(v1);
            ee.p1 = BRep_Tool::Pnt(v2);
            expEnds.push_back(ee);
        }
        for (size_t i = 0; i + 1 < expEnds.size(); i++) {
            const double g3 = expEnds[i].p1.Distance(expEnds[i + 1].p0);
            if (g3 > maxGap3d) maxGap3d = g3;
            if (expEnds[i].havePc && expEnds[i + 1].havePc) {
                const double guv = expEnds[i].uv1.Distance(expEnds[i + 1].uv0);
                if (guv > maxUVgap) maxUVgap = guv;
            }
        }
        if (!expEnds.empty() && (w.Closed() || BRep_Tool::IsClosed(w))) {
            const size_t last = expEnds.size() - 1;
            const double g3 = expEnds[last].p1.Distance(expEnds[0].p0);
            if (g3 > maxGap3d) maxGap3d = g3;
            if (expEnds[last].havePc && expEnds[0].havePc) {
                const double guv = expEnds[last].uv1.Distance(expEnds[0].uv0);
                if (guv > maxUVgap) maxUVgap = guv;
            }
        }
        int ei = 0;
        gp_Pnt2d prevEnd;
        bool havePrev = false;
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next(), ei++) {
            const TopoDS_Edge e = ex.Current();
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
            gp_Pnt2d uv0, uv1;
            if (!pc.IsNull()) {
                uv0 = pc->Value(a);
                uv1 = pc->Value(b);
                gp_Pnt2d p0 = uv0;
                for (int i = 1; i <= nSample; i++) {
                    const double t = a + (b - a) * (double)i / (double)nSample;
                    gp_Pnt2d p1 = pc->Value(t);
                    segs.push_back({p0, p1, ei});
                    p0 = p1;
                }
            }
            if (havePrev && !pc.IsNull()) {
                const double gap = prevEnd.Distance(uv0);
                if (gap > maxUVgapIter) maxUVgapIter = gap;
            }
            if (!pc.IsNull()) {
                prevEnd = uv1;
                havePrev = true;
            }
            Standard_Real f3 = 0, l3 = 0;
            Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
            if (!c3.IsNull() && l3 > f3) {
                const double tm = 0.5 * (f3 + l3);
                mids.push_back(c3->Value(tm));
            }
        }
        if (havePrev && (w.Closed() || BRep_Tool::IsClosed(w)) && !segs.empty()) {
            gp_Pnt2d uvStart = segs.front().a;
            const double gap = prevEnd.Distance(uvStart);
            if (gap > maxUVgapIter) maxUVgapIter = gap;
        }
        for (size_t i = 0; i < segs.size(); i++) {
            for (size_t j = i + 1; j < segs.size(); j++) {
                if (segs[i].ei == segs[j].ei) continue;
                if (std::abs(segs[i].ei - segs[j].ei) == 1) continue;
                if (segIntersect2d(segs[i].a, segs[i].b, segs[j].a, segs[j].b, tol)) nSelfX++;
            }
        }
        for (size_t i = 0; i < mids.size(); i++) {
            for (size_t j = i + 1; j < mids.size(); j++) {
                if (mids[i].Distance(mids[j]) <= tol) nDup++;
            }
        }
    } catch (const Standard_Failure&) {
    }
    if (maxGap3dOut) *maxGap3dOut = maxGap3d;
    if (maxUVgapIterOut) *maxUVgapIterOut = maxUVgapIter;
}

void dumpDiagWire(int rid, int iw, const TopoDS_Face& f, const TopoDS_Wire& w,
                  const BRepCheck_Analyzer& an, double sewTol, const char* stage,
                  const TopTools_IndexedDataMapOfShapeListOfShape* anc,
                  const TopTools_IndexedMapOfShape* faceMap, const std::vector<int>* faceRid,
                  const std::vector<char>* eprimeFill, const RegionSet* rs) {
    int nE = 0, pcurveMissing = 0, nCylE = 0, nFillE = 0, nPlnE = 0, nSeamFree = 0;
    TopLoc_Location loc;
    Handle(Geom_Surface) sF;
    try {
        sF = BRep_Tool::Surface(f, loc);
    } catch (const Standard_Failure&) {
    }
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() != TopAbs_EDGE) continue;
        nE++;
        const TopoDS_Edge e = TopoDS::Edge(it.Value());
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        if (!sF.IsNull()) {
            try {
                (void)BRep_Tool::CurveOnSurface(e, sF, loc, pf, pl, &hasPc);
            } catch (const Standard_Failure&) {
                hasPc = Standard_False;
            }
        }
        if (!hasPc) pcurveMissing++;
        if (!anc || !faceMap || !faceRid || !rs) continue;
        int eidx = anc->FindIndex(e);
        if (eidx <= 0) {
            nSeamFree++;
            continue;
        }
        const int nAnc = anc->FindFromIndex(eidx).Extent();
        if (nAnc < 2) {
            nSeamFree++;
            continue;
        }
        for (TopTools_ListOfShape::Iterator ait(anc->FindFromIndex(eidx)); ait.More(); ait.Next()) {
            const TopoDS_Face fa = TopoDS::Face(ait.Value());
            if (fa.IsSame(f)) continue;
            int ia = faceMap->FindIndex(fa);
            const int orid = (ia > 0 && (size_t)ia < faceRid->size()) ? (*faceRid)[(size_t)ia] : -1;
            const bool fill = eprimeFill && orid >= 0 && (size_t)orid < eprimeFill->size() &&
                              (*eprimeFill)[(size_t)orid];
            const Region* rr = regionById(*rs, orid);
            if (fill) nFillE++;
            else if (rr && rr->type == SurfType::Cylinder) nCylE++;
            else if (rr && rr->type == SurfType::Plane) nPlnE++;
            else nSeamFree++;
        }
    }
    const char* infPt = "?";
    try {
        BRepTopAdaptor_FClass2d cl(f, Precision::PConfusion());
        infPt = topAbsStateName(cl.PerformInfinitePoint());
    } catch (const Standard_Failure&) {
    }
    double areaUV = 0.0, maxUVgap = 0.0, maxUVgapIter = 0.0, maxGap3d = 0.0;
    int nSelfX = 0, nDup = 0;
    wireDiagMetrics(f, w, sewTol, areaUV, nSelfX, maxUVgap, nDup, &maxUVgapIter, &maxGap3d);
    std::fprintf(stderr,
                 "DIAG_WIRE rid=%d iw=%d nE=%d ori=%s faceOri=%s infPt=%s areaUV=%.6e nSelfX=%d "
                 "maxUVgap=%.6e maxUVgapIter=%.6e maxGap3d=%.6e nDup=%d pcurveMissing=%d nCylE=%d "
                 "nFillE=%d nPlnE=%d nSeamFree=%d sewTol=%.6e stage=%s",
                 rid, iw, nE, w.Orientation() == TopAbs_FORWARD ? "F" : "R",
                 f.Orientation() == TopAbs_FORWARD ? "F" : "R", infPt, areaUV, nSelfX, maxUVgap,
                 maxUVgapIter, maxGap3d, nDup, pcurveMissing, nCylE, nFillE, nPlnE, nSeamFree,
                 sewTol, stage ? stage : "-");
    try {
        Handle(BRepCheck_Result) res = an.Result(w);
        if (!res.IsNull()) {
            if (res->IsStatusOnShape(f)) dumpBRepStatusList(res->StatusOnShape(f));
            dumpBRepStatuses(res);
        }
        Handle(BRepCheck_Result) fr = an.Result(f);
        if (!fr.IsNull()) dumpBRepStatuses(fr);
    } catch (const Standard_Failure&) {
    }
    std::fprintf(stderr, "\n");
}

void dumpDiagWiresOfFace(int rid, const TopoDS_Face& f, double sewTol, const char* stage = nullptr,
                         const TopTools_IndexedDataMapOfShapeListOfShape* anc = nullptr,
                         const TopTools_IndexedMapOfShape* faceMap = nullptr,
                         const std::vector<int>* faceRid = nullptr,
                         const std::vector<char>* eprimeFill = nullptr,
                         const RegionSet* rs = nullptr) {
    if (f.IsNull()) return;
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        int iw = 0;
        for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), iw++)
            dumpDiagWire(rid, iw, f, TopoDS::Wire(wx.Current()), an, sewTol, stage, anc, faceMap,
                         faceRid, eprimeFill, rs);
    } catch (const Standard_Failure&) {
    }
}

const char* geomCurveKind(const Handle(Geom_Curve)& c) {
    if (c.IsNull()) return "none";
    Handle(Geom_Curve) b = c;
    Handle(Geom_TrimmedCurve) tr = Handle(Geom_TrimmedCurve)::DownCast(c);
    if (!tr.IsNull() && !tr->BasisCurve().IsNull()) b = tr->BasisCurve();
    if (b->DynamicType() == STANDARD_TYPE(Geom_Circle)) return "circ";
    if (b->DynamicType() == STANDARD_TYPE(Geom_Line)) return "lin";
    if (b->DynamicType() == STANDARD_TYPE(Geom_Ellipse)) return "elips";
    return "other";
}

int findSharedCylRid(const TopoDS_Edge& e, const std::vector<TopoDS_Face>& built,
                     const std::vector<int>& builtRid, const RegionSet& rs, int selfRid) {
    const void* ts = e.TShape().get();
    if (!ts) return -1;
    for (size_t i = 0; i < built.size(); i++) {
        const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
        if (rid == selfRid) continue;
        const Region* rr = regionById(rs, rid);
        if (!rr || rr->type != SurfType::Cylinder) continue;
        for (TopExp_Explorer ex(built[i], TopAbs_EDGE); ex.More(); ex.Next()) {
            if (ex.Current().TShape().get() == ts) return rid;
        }
    }
    return -1;
}

void dumpDiagM5(int rid, const TopoDS_Edge& e, const TopoDS_Face& f, int sharedCylRid, int st) {
    Standard_Real f3 = 0, l3 = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
    const char* cty = geomCurveKind(c3);
    Standard_Real f2 = 0, l2 = 0;
    Standard_Boolean stored = Standard_False;
    Handle(Geom2d_Curve) c2;
    if (!f.IsNull()) c2 = BRep_Tool::CurveOnSurface(e, f, f2, l2, &stored);
    const char* cls = stored ? "bound" : (c2.IsNull() ? "unbound" : "synth");
    const double span3 = (double)(l3 - f3);
    const double span2 = (double)(l2 - f2);
    std::fprintf(stderr,
                 "DIAG_M5 rid=%d class=%s cty=%s span3=%.6f span2=%.6f stored=%d "
                 "sharedCylRid=%d st=%d\n",
                 rid, cls, cty, span3, span2, stored ? 1 : 0, sharedCylRid, st);
}

void dumpM5OnShell(const BRepCheck_Analyzer& an, const std::vector<TopoDS_Face>& built,
                   const std::vector<int>& builtRid, const RegionSet& rs) {
    if (!diagP2Enabled()) return;
    try {
        for (size_t i = 0; i < built.size(); i++) {
            const TopoDS_Face& f = built[i];
            if (f.IsNull()) continue;
            const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
            for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
                const TopoDS_Edge e = TopoDS::Edge(ex.Current());
                Handle(BRepCheck_Result) er = an.Result(e);
                if (er.IsNull()) continue;
                for (BRepCheck_ListOfStatus::Iterator it(er->Status()); it.More(); it.Next()) {
                    if ((int)it.Value() != 11) continue;
                    dumpDiagM5(rid, e, f, findSharedCylRid(e, built, builtRid, rs, rid), 11);
                }
            }
        }
    } catch (const Standard_Failure&) {
    }
}

const char* builtAsName(BuiltAs a) {
    switch (a) {
        case BuiltAs::Single: return "Single";
        case BuiltAs::Seamed360: return "Seamed360";
        case BuiltAs::TwoHalves: return "TwoHalves";
        case BuiltAs::ExplodedToFacets: return "ExplodedToFacets";
        default: return "NotBuilt";
    }
}

bool sliverCensusEnabled() { return diagP2Enabled() || diagPlatesEnabled(); }

const char* sliverCensusPartName(int nTri) {
    switch (nTri) {
        case 3338: return "HP";
        case 908: return "handle-lock";
        case 15300: return "Body11";
        case 412: return "S02";
        case 44: return "S11-b";
        default: {
            static char buf[32];
            std::snprintf(buf, sizeof(buf), "nTri%d", nTri);
            return buf;
        }
    }
}

double regionMeshArea(const MeshView& mv, const Region& r) {
    double a = 0;
    for (int lt : r.tris) {
        if (lt < 0 || (size_t)lt >= mv.nTri || !mv.compTris || !mv.tris) continue;
        const int gt = mv.compTris[lt];
        const int* T = mv.tris[gt];
        const gp_Pnt p0 = pntOf(mv, T[0]);
        const gp_Pnt p1 = pntOf(mv, T[1]);
        const gp_Pnt p2 = pntOf(mv, T[2]);
        a += 0.5 * gp_Vec(p0, p1).Crossed(gp_Vec(p0, p2)).Magnitude();
    }
    return a;
}

int regionMeshVertCount(const MeshView& mv, const Region& r) {
    std::unordered_set<int> v;
    for (int lt : r.tris) {
        if (lt < 0 || (size_t)lt >= mv.nTri || !mv.compTris || !mv.tris) continue;
        const int gt = mv.compTris[lt];
        const int* T = mv.tris[gt];
        v.insert(T[0]);
        v.insert(T[1]);
        v.insert(T[2]);
    }
    return (int)v.size();
}

thread_local int gRid100PrevValid = -1;

void dumpSliverRegionCensus(const MeshView& mv, const RegionSet& rs, double sewTol) {
    if (!sliverCensusEnabled()) return;
    const char* part = sliverCensusPartName((int)mv.nTri);
    int nExam = 0;
    std::fprintf(stderr, "DIAG_SLIVERCENSUS part=%s sewTol=%.6f nRegions=%zu\n", part, sewTol,
                 rs.regions.size());
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane && r.type != SurfType::Cylinder) continue;
        ++nExam;
        const double ratioSew =
            sewTol > 0.0 ? r.maxVertexDev / sewTol : r.maxVertexDev;
        const char* typ = (r.type == SurfType::Plane) ? "plane" : "cyl";
        std::fprintf(stderr,
                     "DIAG_SLIVERCENSUS part=%s rid=%d type=%s nTri=%zu nVert=%d "
                     "maxVertexDev=%.6f ratioSew=%.3f rmsDev=%.6f analyticNow=%d builtAs=%s "
                     "area=%.4f\n",
                     part, r.id, typ, r.tris.size(), regionMeshVertCount(mv, r),
                     r.maxVertexDev, ratioSew, r.rmsVertexDev, isAnalytic(&r) ? 1 : 0,
                     builtAsName(r.builtAs), regionMeshArea(mv, r));
    }
    std::fprintf(stderr, "DIAG_SLIVERCENSUS_SUM part=%s examined=%d present=%zu\n", part, nExam,
                 rs.regions.size());
}

void dumpChainFitCensus(const MeshView& mv, const RegionSet& rs,
                        const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                        double sewTol) {
    if (!sliverCensusEnabled() || mv.nTri != 3338) return;
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        if (geom[ci].curve.kind == AnalyticCurve::None) continue;
        const BoundaryChain& ch = rs.chains[ci];
        double maxDist = 0;
        for (int lv : ch.meshVerts)
            maxDist = std::max(maxDist, curveResidual(geom[ci].curve, pntOf(mv, lv)));
        const double snapCap = geom[ci].chainSnapCap;
        const double ratioSew = sewTol > 0.0 ? maxDist / sewTol : maxDist;
        const double ratioCap = snapCap > 0.0 ? maxDist / snapCap : maxDist;
        const char* kind = (geom[ci].curve.kind == AnalyticCurve::Circ)   ? "circ"
                           : (geom[ci].curve.kind == AnalyticCurve::Lin) ? "lin"
                                                                         : "elips";
        std::fprintf(stderr,
                     "DIAG_CHAINFIT ci=%d kind=%s regA=%d regB=%d nMeshV=%zu "
                     "maxDistToOwnVerts=%.6f ratioSew=%.3f snapCap=%.6f ratioCap=%.3f collapsed=%d\n",
                     (int)ci, kind, geom[ci].regA, geom[ci].regB, ch.meshVerts.size(), maxDist,
                     ratioSew, snapCap, ratioCap, collapsed[ci] ? 1 : 0);
    }
}

void dumpDiagTheta1136(const MeshView& mv, const RegionSet& rs) {
    if (!sliverCensusEnabled() || mv.nTri != 3338) return;
    const Region* r1 = regionById(rs, 1);
    const Region* r136 = regionById(rs, 136);
    if (!r1 || !r136) return;
    const int vid = 227;
    if (vid < 0 || (size_t)vid >= mv.nVtx) return;
    const gp_Pnt pV = pntOf(mv, vid);
    const double theta = dihedralBetweenRegions(r1, r136, pV);
    const double d136 = regionPointDev(r136, pV);
    std::fprintf(stderr,
                 "DIAG_THETA ridA=1 ridB=136 vid=%d theta_rad=%.8f theta_deg=%.6f dPlane136=%.8f "
                 "rid136_nTri=%zu rid136_maxVertexDev=%.6f\n",
                 vid, theta, theta * 180.0 / kPi, d136, r136->tris.size(), r136->maxVertexDev);
}

double chainMaxDistToOwnVerts(const AnalyticCurve& c, const MeshView& mv, const BoundaryChain& ch) {
    double maxDist = 0;
    for (int lv : ch.meshVerts) maxDist = std::max(maxDist, curveResidual(c, pntOf(mv, lv)));
    return maxDist;
}

bool fittedPlaneIntersection(const Region* A, const Region* B, gp_Lin& out) {
    if (!A || !B || A->type != SurfType::Plane || B->type != SurfType::Plane) return false;
    try {
        IntAna_QuadQuadGeo iq(asPlane(*A), asPlane(*B), Precision::Angular(),
                              Precision::Confusion());
        if (!iq.IsDone() || iq.NbSolutions() < 1) return false;
        if (iq.TypeInter() != IntAna_Line && iq.TypeInter() != IntAna_PointAndCircle)
            return false;
        out = iq.Line(1);
        return true;
    } catch (const Standard_Failure&) {
        return false;
    }
}

void dumpBoundCredCensus(const MeshView& mv, const RegionSet& rs,
                         const std::vector<ChainGeom>& geom, const std::vector<char>& /*collapsed*/,
                         double sewTol) {
    if (!sliverCensusEnabled()) return;
    const char* part = sliverCensusPartName((int)mv.nTri);
    const size_t nCh = rs.chains.size();
    int nExam = 0, nIsland = 0;
    std::fprintf(stderr, "DIAG_BOUNDCRED part=%s sewTol=%.6f nRegions=%zu nChains=%zu\n", part,
                 sewTol, rs.regions.size(), nCh);
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane && r.type != SurfType::Cylinder) continue;
        ++nExam;
        int nChains = 0, nChainsSew = 0, nChainsCap = 0;
        for (size_t ci = 0; ci < nCh; ci++) {
            const BoundaryChain& ch = rs.chains[ci];
            if (ch.regA != r.id && ch.regB != r.id) continue;
            ++nChains;
            if (ci >= geom.size() || geom[ci].curve.kind == AnalyticCurve::None) continue;
            const double maxDist = chainMaxDistToOwnVerts(geom[ci].curve, mv, ch);
            const double snapCap = geom[ci].chainSnapCap;
            if (maxDist <= sewTol) ++nChainsSew;
            if (snapCap > 0.0 && maxDist <= snapCap) ++nChainsCap;
        }
        const int allFailSew = (nChainsSew == 0) ? 1 : 0;
        const int allFailCap = (nChainsCap == 0) ? 1 : 0;
        if (nChains == 0) ++nIsland;
        const double ratioSew = sewTol > 0.0 ? r.maxVertexDev / sewTol : r.maxVertexDev;
        const char* typ = (r.type == SurfType::Plane) ? "plane" : "cyl";
        std::fprintf(stderr,
                     "DIAG_BOUNDCRED part=%s rid=%d type=%s nTri=%zu nVert=%d nChains=%d "
                     "nChainsSew=%d nChainsCap=%d allFailSew=%d allFailCap=%d "
                     "maxVertexDev=%.6f ratioSew=%.3f area=%.4f analyticNow=%d\n",
                     part, r.id, typ, r.tris.size(), regionMeshVertCount(mv, r), nChains,
                     nChainsSew, nChainsCap, allFailSew, allFailCap, r.maxVertexDev, ratioSew,
                     regionMeshArea(mv, r), isAnalytic(&r) ? 1 : 0);
    }
    std::fprintf(stderr, "DIAG_BOUNDCRED_SUM part=%s examined=%d present=%zu nIsland=%d\n", part,
                 nExam, rs.regions.size(), nIsland);
}

void dumpChainLineVsIntersection(const MeshView& mv, const RegionSet& rs,
                                 const std::vector<ChainGeom>& geom, int ci, int plA, int plB,
                                 int vidHint) {
    if (ci < 0 || (size_t)ci >= geom.size() || (size_t)ci >= rs.chains.size()) return;
    const AnalyticCurve& c = geom[(size_t)ci].curve;
    const char* kind = (c.kind == AnalyticCurve::Circ)   ? "circ"
                       : (c.kind == AnalyticCurve::Lin) ? "lin"
                                                         : (c.kind == AnalyticCurve::Elips) ? "elips"
                                                                                            : "none";
    const Region* A = regionById(rs, plA);
    const Region* B = regionById(rs, plB);
    gp_Lin ix;
    const int ixOk = fittedPlaneIntersection(A, B, ix) ? 1 : 0;
    gp_Pnt orig(0, 0, 0);
    gp_Dir dir(0, 0, 1);
    if (c.kind == AnalyticCurve::Lin) {
        orig = c.lin.Location();
        dir = c.lin.Direction();
    }
    double angleDeg = -1.0;
    if (c.kind == AnalyticCurve::Lin && ixOk) {
        double ang = dir.Angle(ix.Direction());
        if (ang > kPi * 0.5) ang = kPi - ang;
        angleDeg = ang * 180.0 / kPi;
    }
    gp_Pnt pV(0, 0, 0);
    int vid = vidHint;
    if (vid >= 0 && (size_t)vid < mv.nVtx)
        pV = pntOf(mv, vid);
    else if (!rs.chains[(size_t)ci].meshVerts.empty()) {
        vid = rs.chains[(size_t)ci].meshVerts.front();
        pV = pntOf(mv, vid);
    }
    const double dChain = (c.kind == AnalyticCurve::Lin) ? c.lin.Distance(pV) : -1.0;
    const double dIx = ixOk ? ix.Distance(pV) : -1.0;
    const char* cls = (c.kind == AnalyticCurve::Lin && ixOk && angleDeg >= 0.0 && angleDeg < 1.0)
                          ? "planeIntersection"
                          : "fit";
    std::fprintf(stderr,
                 "DIAG_CHAINLINE ci=%d kind=%s regA=%d regB=%d plA=%d plB=%d "
                 "orig=(%.6f,%.6f,%.6f) dir=(%.8f,%.8f,%.8f) "
                 "ixOk=%d ixOrig=(%.6f,%.6f,%.6f) ixDir=(%.8f,%.8f,%.8f) "
                 "angleDeg=%.6f vid=%d dChain=%.8f dIx=%.8f class=%s\n",
                 ci, kind, geom[(size_t)ci].regA, geom[(size_t)ci].regB, plA, plB, orig.X(),
                 orig.Y(), orig.Z(), dir.X(), dir.Y(), dir.Z(), ixOk, ix.Location().X(),
                 ix.Location().Y(), ix.Location().Z(), ix.Direction().X(), ix.Direction().Y(),
                 ix.Direction().Z(), angleDeg, vid, dChain, dIx, cls);
}

void dumpChainLineCensus(const MeshView& mv, const RegionSet& rs,
                         const std::vector<ChainGeom>& geom) {
    if (!sliverCensusEnabled() || mv.nTri != 3338) return;
    dumpChainLineVsIntersection(mv, rs, geom, 39, 1, 136, 227);
    dumpChainLineVsIntersection(mv, rs, geom, 122, 2, 174, -1);
}

void formatStatusList(const Handle(BRepCheck_Result)& res, char* buf, size_t n) {
    buf[0] = '\0';
    if (res.IsNull() || n < 4) return;
    size_t used = 0;
    buf[used++] = '[';
    buf[used] = '\0';
    bool any = false;
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next()) {
        const int st = (int)it.Value();
        if (st == (int)BRepCheck_NoError) continue;
        int w = std::snprintf(buf + used, n - used, "%s%d", any ? "," : "", st);
        if (w < 0 || (size_t)w >= n - used) break;
        used += (size_t)w;
        any = true;
    }
    if (used + 2 < n) {
        buf[used++] = ']';
        buf[used] = '\0';
    }
}

void diagRid100Status(const char* pass, int rid, const TopoDS_Face& f) {
    if (rid != 100 || !sliverCensusEnabled() || !pass) return;
    int valid = 0;
    char stbuf[96];
    stbuf[0] = '\0';
    if (f.IsNull()) {
        std::snprintf(stbuf, sizeof(stbuf), "null");
    } else {
        try {
            BRepCheck_Analyzer an(f, Standard_True);
            valid = an.IsValid() ? 1 : 0;
            Handle(BRepCheck_Result) res = an.Result(f);
            formatStatusList(res, stbuf, sizeof(stbuf));
            if (!stbuf[0]) std::snprintf(stbuf, sizeof(stbuf), valid ? "valid" : "invalid");
        } catch (const Standard_Failure&) {
            valid = 0;
            std::snprintf(stbuf, sizeof(stbuf), "exception");
        }
    }
    std::fprintf(stderr, "DIAG_RID100 pass=%s status=%s valid=%d\n", pass, stbuf, valid);
    if (gRid100PrevValid == 1 && valid == 0)
        std::fprintf(stderr, "DIAG_RID100 first-invalid pass=%s\n", pass);
    if (!f.IsNull()) gRid100PrevValid = valid;
}

int countWireEdges(const TopoDS_Wire& w) {
    int n = 0;
    if (w.IsNull()) return 0;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) n++;
    } catch (const Standard_Failure&) {
        for (TopoDS_Iterator it(w); it.More(); it.Next())
            if (it.Value().ShapeType() == TopAbs_EDGE) n++;
    }
    return n;
}

int pcurveMissingOnWireSurf(const TopoDS_Wire& w, const Handle(Geom_Surface)& surf) {
    int n = 0;
    if (w.IsNull() || surf.IsNull()) return 0;
    const TopLoc_Location loc;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
            const TopoDS_Edge e = TopoDS::Edge(ex.Current());
            Standard_Boolean hasPc = Standard_False;
            Standard_Real pf = 0, pl = 0;
            (void)BRep_Tool::CurveOnSurface(e, surf, loc, pf, pl, &hasPc);
            if (!hasPc) n++;
        }
    } catch (const Standard_Failure&) {
    }
    return n;
}

int pcurveMissingOnFace(const TopoDS_Face& f) {
    int n = 0;
    if (f.IsNull()) return 0;
    TopLoc_Location loc;
    Handle(Geom_Surface) sF = BRep_Tool::Surface(f, loc);
    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
        const TopoDS_Edge e = TopoDS::Edge(ex.Current());
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        if (!sF.IsNull()) {
            try {
                (void)BRep_Tool::CurveOnSurface(e, sF, loc, pf, pl, &hasPc);
            } catch (const Standard_Failure&) {
            }
        }
        if (!hasPc) n++;
    }
    return n;
}

int countFaceEdges(const TopoDS_Face& f) {
    int n = 0;
    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) n++;
    return n;
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
                    const std::vector<int>& builtRid, const RegionSet& rs, double sewTol) {
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
            // Name the face and the sub-shape that made it invalid, the same way
            // DIAG_R2SUB does at the R2 probe (ee834db): a face reported with an
            // empty own-status list is invalid only in the CONTEXT of a wire or
            // an edge, and the context list is the only place that says so.
            // A rid=-1 face carries no region, so its shape is the only name it
            // has: wire count, edge count, the wire's closure and its signed
            // area in UV -- a zero/negative UV area on the single outer wire is
            // exactly what BRepCheck_Face reports as UnorientableShape.
            {
                int nW = 0, nE = 0;
                for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next()) nW++;
                for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) nE++;
                std::fprintf(stderr, "DIAG_SHELL_FACESHAPE rid=%d nW=%d nE=%d", rid, nW, nE);
                int wi = 0;
                for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), wi++) {
                    const TopoDS_Wire w = TopoDS::Wire(wx.Current());
                    int nWE = 0;
                    for (TopExp_Explorer ex(w, TopAbs_EDGE); ex.More(); ex.Next()) nWE++;
                    std::fprintf(stderr, " w%d=[nE=%d closed=%d ori=%s uvA=%.6g]", wi, nWE,
                                 (w.Closed() || BRep_Tool::IsClosed(w)) ? 1 : 0,
                                 w.Orientation() == TopAbs_FORWARD ? "F" : "R",
                                 pcurveSignedArea(f, w));
                }
                gp_Pnt p0;
                bool haveP = false;
                for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
                    p0 = BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()));
                    haveP = true;
                    break;
                }
                if (haveP)
                    std::fprintf(stderr, " p=(%.5f,%.5f,%.5f)", p0.X(), p0.Y(), p0.Z());
                {
                    TopLoc_Location fl3;
                    Handle(Geom_Surface) fs3 = BRep_Tool::Surface(f, fl3);
                    Handle(Geom_Plane) pl3 = Handle(Geom_Plane)::DownCast(fs3);
                    if (!pl3.IsNull()) {
                        gp_Pln pp = pl3->Pln();
                        if (!fl3.IsIdentity()) pp.Transform(fl3.Transformation());
                        std::fprintf(stderr, " plnP=(%.7f,%.7f,%.7f) plnN=(%.7f,%.7f,%.7f)",
                                     pp.Location().X(), pp.Location().Y(), pp.Location().Z(),
                                     pp.Axis().Direction().X(), pp.Axis().Direction().Y(),
                                     pp.Axis().Direction().Z());
                    }
                    std::fprintf(stderr, " ftol=%.6g", BRep_Tool::Tolerance(f));
                }
                std::fprintf(stderr, "\n");
                TopTools_IndexedMapOfShape subMap;
                TopExp::MapShapes(f, subMap);
                for (int si = 1; si <= subMap.Extent(); si++) {
                    const TopoDS_Shape& sh = subMap(si);
                    Handle(BRepCheck_Result) sr = fan.Result(sh);
                    if (sr.IsNull()) continue;
                    char ctxSt[128];
                    ctxSt[0] = '\0';
                    size_t cu = 0;
                    for (sr->InitContextIterator(); sr->MoreShapeInContext();
                         sr->NextShapeInContext()) {
                        for (BRepCheck_ListOfStatus::Iterator it(sr->StatusOnShape()); it.More();
                             it.Next()) {
                            if (it.Value() == BRepCheck_NoError) continue;
                            const int wn =
                                std::snprintf(ctxSt + cu, sizeof(ctxSt) - cu, "%s%s",
                                              cu ? "," : "", brepCheckName((int)it.Value()));
                            if (wn < 0 || (size_t)wn >= sizeof(ctxSt) - cu) break;
                            cu += (size_t)wn;
                        }
                    }
                    char ownSt[128];
                    formatStatusList(sr, ownSt, sizeof(ownSt));
                    const bool ownClean =
                        (ownSt[0] == '\0' || (ownSt[0] == '[' && ownSt[1] == ']'));
                    if (ownClean && ctxSt[0] == '\0') continue;
                    gp_Pnt sp;
                    double dev = -1.0;
                    double etol = 0.0;
                    int sameP = -1;
                    if (sh.ShapeType() == TopAbs_VERTEX) {
                        sp = BRep_Tool::Pnt(TopoDS::Vertex(sh));
                        etol = BRep_Tool::Tolerance(TopoDS::Vertex(sh));
                    } else if (sh.ShapeType() == TopAbs_EDGE) {
                        const TopoDS_Edge ee = TopoDS::Edge(sh);
                        sameP = BRep_Tool::SameParameter(ee) ? 1 : 0;
                        etol = BRep_Tool::Tolerance(ee);
                        Standard_Real ef = 0, el = 0;
                        Handle(Geom_Curve) ec = BRep_Tool::Curve(ee, ef, el);
                        if (!ec.IsNull()) sp = ec->Value(ef);
                        TopLoc_Location fl;
                        Handle(Geom_Surface) fs2 = BRep_Tool::Surface(f, fl);
                        // No stored pcurve on a plane: measure the 3d curve
                        // against the SURFACE itself (the distance a synthesised
                        // pcurve would have to absorb), which is what
                        // InvalidCurveOnSurface reports.
                        if (!ec.IsNull() && !fs2.IsNull()) {
                            dev = 0.0;
                            GeomAPI_ProjectPointOnSurf pr;
                            for (int kk = 0; kk <= 64; kk++) {
                                const double t = ef + (el - ef) * (double)kk / 64.0;
                                gp_Pnt q = ec->Value(t);
                                if (!fl.IsIdentity())
                                    q.Transform(fl.Transformation().Inverted());
                                pr.Init(q, fs2);
                                if (pr.NbPoints() > 0) dev = std::max(dev, pr.LowerDistance());
                            }
                        }
                    }
                    if (sh.ShapeType() == TopAbs_EDGE) {
                        const TopoDS_Edge ee = TopoDS::Edge(sh);
                        Standard_Real ef = 0, el = 0;
                        Handle(Geom_Curve) ec = BRep_Tool::Curve(ee, ef, el);
                        const char* cty = "none";
                        if (!ec.IsNull()) {
                            if (ec->DynamicType() == STANDARD_TYPE(Geom_Line)) cty = "line";
                            else if (ec->DynamicType() == STANDARD_TYPE(Geom_Circle)) cty = "circ";
                            else if (ec->DynamicType() == STANDARD_TYPE(Geom_Ellipse)) cty = "elips";
                            else if (ec->DynamicType() == STANDARD_TYPE(Geom_BSplineCurve))
                                cty = "bspl";
                            else cty = "other";
                        }
                        Standard_Real pf = 0, pl = 0;
                        Standard_Boolean stored = Standard_False;
                        TopLoc_Location fl2;
                        Handle(Geom_Surface) fsur = BRep_Tool::Surface(f, fl2);
                        if (!fsur.IsNull())
                            (void)BRep_Tool::CurveOnSurface(ee, fsur, fl2, pf, pl, &stored);
                        gp_Pnt pa = ec.IsNull() ? gp_Pnt() : ec->Value(ef);
                        gp_Pnt pb = ec.IsNull() ? gp_Pnt() : ec->Value(el);
                        std::fprintf(stderr,
                                     "DIAG_SHELL_SUBE rid=%d i=%d cty=%s stored=%d "
                                     "range=[%.9g,%.9g] a=(%.7f,%.7f,%.7f) b=(%.7f,%.7f,%.7f) "
                                     "len=%.9g\n",
                                     rid, si, cty, stored ? 1 : 0, ef, el, pa.X(), pa.Y(), pa.Z(),
                                     pb.X(), pb.Y(), pb.Z(), pa.Distance(pb));
                    }
                    std::fprintf(stderr,
                                 "DIAG_SHELL_SUB rid=%d i=%d type=%s own=%s ctx=[%s] "
                                 "dev=%.6g sameP=%d tol=%.6g p=(%.5f,%.5f,%.5f)\n",
                                 rid, si,
                                 sh.ShapeType() == TopAbs_EDGE     ? "EDGE"
                                 : sh.ShapeType() == TopAbs_VERTEX ? "VERTEX"
                                 : sh.ShapeType() == TopAbs_WIRE   ? "WIRE"
                                 : sh.ShapeType() == TopAbs_FACE   ? "FACE"
                                                                   : "OTHER",
                                 ownSt, ctxSt, dev, sameP, etol, sp.X(), sp.Y(), sp.Z());
                }
            }
            if (diagP2Enabled()) dumpDiagWiresOfFace(rid, f, sewTol, "shell-invalid");
        }
        if (diagP2Enabled()) dumpM5OnShell(an, built, builtRid, rs);
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

extern thread_local std::unordered_set<const void*> gSeamTShapes;

bool wireStatusOnFace(const Handle(BRepCheck_Result)& res, const TopoDS_Face& f,
                      BRepCheck_Status want) {
    if (res.IsNull()) return false;
    if (res->IsStatusOnShape(f)) {
        for (BRepCheck_ListOfStatus::Iterator it(res->StatusOnShape(f)); it.More(); it.Next())
            if (it.Value() == want) return true;
    }
    for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More(); it.Next())
        if (it.Value() == want) return true;
    return false;
}

int wireEdgeIterIndex(const TopoDS_Wire& w, const TopoDS_Edge& target) {
    int i = 0;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() != TopAbs_EDGE) continue;
        const TopoDS_Edge e = TopoDS::Edge(it.Value());
        if (e.IsSame(target) || diagTShapePtr(e) == diagTShapePtr(target)) return i;
        i++;
    }
    return -1;
}

bool edgeShareVtxTShape(const TopoDS_Edge& e1, const TopoDS_Edge& e2) {
    TopoDS_Vertex a0, a1, b0, b1;
    TopExp::Vertices(e1, a0, a1, Standard_True);
    TopExp::Vertices(e2, b0, b1, Standard_True);
    const void* tsA[] = {diagTShapePtr(a0), diagTShapePtr(a1)};
    const void* tsB[] = {diagTShapePtr(b0), diagTShapePtr(b1)};
    for (const void* x : tsA)
        for (const void* y : tsB)
            if (x && x == y) return true;
    return false;
}

bool edgeStretchOverlap(const TopoDS_Edge& ea, const TopoDS_Edge& eb, double sewTol,
                        gp_Pnt& at) {
    try {
        BRepAdaptor_Curve ca(ea), cb(eb);
        if (ca.LastParameter() <= ca.FirstParameter() || cb.LastParameter() <= cb.FirstParameter())
            return false;
        const bool swap = (ca.LastParameter() - ca.FirstParameter()) >
                          (cb.LastParameter() - cb.FirstParameter());
        const BRepAdaptor_Curve& sh = swap ? cb : ca;
        const BRepAdaptor_Curve& lg = swap ? ca : cb;
        const double sf = sh.FirstParameter(), sl = sh.LastParameter();
        const double lf = lg.FirstParameter(), ll = lg.LastParameter();
        int nHit = 0;
        for (int k = 0; k < 5; k++) {
            const double t = sf + (sl - sf) * (double)k / 4.0;
            const gp_Pnt p = sh.Value(t);
            double bestD = 1e300;
            gp_Pnt bestQ;
            for (int j = 0; j <= 20; j++) {
                const double u = lf + (ll - lf) * (double)j / 20.0;
                const gp_Pnt q = lg.Value(u);
                const double d = p.Distance(q);
                if (d < bestD) {
                    bestD = d;
                    bestQ = q;
                }
            }
            if (bestD <= sewTol) {
                nHit++;
                if (k == 2) at = bestQ;
            }
        }
        return nHit >= 3;
    } catch (const Standard_Failure&) {
        return false;
    }
}

double edgeLen3d(const TopoDS_Edge& e) {
    try {
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull() || l <= f) return 0.0;
        return c->Value(f).Distance(c->Value(l));
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

int countStretchPairs3d(const TopoDS_Wire& w, double sewTol) {
    std::vector<TopoDS_Edge> edges;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() == TopAbs_EDGE) edges.push_back(TopoDS::Edge(it.Value()));
    }
    int n = 0;
    const double tol = std::max(sewTol, Precision::Confusion());
    gp_Pnt at;
    for (size_t i = 0; i < edges.size(); i++) {
        for (size_t j = i + 1; j < edges.size(); j++) {
            if (edgeShareVtxTShape(edges[i], edges[j])) continue;
            if (edgeStretchOverlap(edges[i], edges[j], tol, at) ||
                edgeStretchOverlap(edges[j], edges[i], tol, at))
                n++;
        }
    }
    return n;
}

int countOrientedPcSelfX(const TopoDS_Face& f, const TopoDS_Wire& w, double tol) {
    struct Seg {
        gp_Pnt2d a, b;
        int ei;
    };
    std::vector<Seg> segs;
    int ei = 0;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next(), ei++) {
            const TopoDS_Edge e = ex.Current();
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
            if (pc.IsNull()) continue;
            if (e.Orientation() == TopAbs_REVERSED) std::swap(a, b);
            gp_Pnt2d p0 = pc->Value(a);
            for (int s = 1; s <= 8; s++) {
                const double t = a + (b - a) * (double)s / 8.0;
                gp_Pnt2d p1 = pc->Value(t);
                segs.push_back({p0, p1, ei});
                p0 = p1;
            }
        }
    } catch (const Standard_Failure&) {
        return 0;
    }
    int n = 0;
    for (size_t i = 0; i < segs.size(); i++) {
        for (size_t j = i + 1; j < segs.size(); j++) {
            if (segs[i].ei == segs[j].ei) continue;
            if (std::abs(segs[i].ei - segs[j].ei) == 1) continue;
            if (segIntersect2d(segs[i].a, segs[i].b, segs[j].a, segs[j].b, tol)) n++;
        }
    }
    return n;
}

int sharedVtxTShapeCount(const TopoDS_Edge& e1, const TopoDS_Edge& e2) {
    TopoDS_Vertex a0, a1, b0, b1;
    TopExp::Vertices(e1, a0, a1, Standard_True);
    TopExp::Vertices(e2, b0, b1, Standard_True);
    const void* tsA[] = {diagTShapePtr(a0), diagTShapePtr(a1)};
    const void* tsB[] = {diagTShapePtr(b0), diagTShapePtr(b1)};
    int n = 0;
    for (const void* x : tsA) {
        if (!x) continue;
        for (const void* y : tsB)
            if (x == y) {
                n++;
                break;
            }
    }
    return n;
}

gp_Pnt sharedVtxPnt(const TopoDS_Edge& e1, const TopoDS_Edge& e2) {
    TopoDS_Vertex a0, a1, b0, b1;
    TopExp::Vertices(e1, a0, a1, Standard_True);
    TopExp::Vertices(e2, b0, b1, Standard_True);
    const void* tsB[] = {diagTShapePtr(b0), diagTShapePtr(b1)};
    TopoDS_Vertex vs[] = {a0, a1};
    for (int i = 0; i < 2; i++) {
        const void* x = diagTShapePtr(vs[i]);
        if (!x) continue;
        for (const void* y : tsB)
            if (x == y) return BRep_Tool::Pnt(vs[i]);
    }
    return gp_Pnt();
}

gp_Vec edgeDir3d(const TopoDS_Edge& e) {
    try {
        BRepAdaptor_Curve c(e);
        const double f = c.FirstParameter(), l = c.LastParameter();
        gp_Pnt p0 = c.Value(f), p1 = c.Value(l);
        if (e.Orientation() == TopAbs_REVERSED) std::swap(p0, p1);
        gp_Vec v(p0, p1);
        if (v.Magnitude() <= Precision::Confusion()) return gp_Vec(0, 0, 0);
        v.Normalize();
        return v;
    } catch (const Standard_Failure&) {
        return gp_Vec(0, 0, 0);
    }
}

double overlapLen3dTrue(const TopoDS_Edge& ea, const TopoDS_Edge& eb, double sewTol) {
    try {
        BRepAdaptor_Curve ca(ea), cb(eb);
        const bool swap = (ca.LastParameter() - ca.FirstParameter()) >
                          (cb.LastParameter() - cb.FirstParameter());
        const BRepAdaptor_Curve& sh = swap ? cb : ca;
        const BRepAdaptor_Curve& lg = swap ? ca : cb;
        const double sf = sh.FirstParameter(), sl = sh.LastParameter();
        const double lf = lg.FirstParameter(), ll = lg.LastParameter();
        const double len = sh.Value(sf).Distance(sh.Value(sl));
        int nHit = 0;
        const int nS = 21;
        for (int k = 0; k < nS; k++) {
            const double t = sf + (sl - sf) * (double)k / (double)(nS - 1);
            const gp_Pnt p = sh.Value(t);
            double bestD = 1e300;
            for (int j = 0; j <= 40; j++) {
                const double u = lf + (ll - lf) * (double)j / 40.0;
                const double d = p.Distance(lg.Value(u));
                if (d < bestD) bestD = d;
            }
            if (bestD <= sewTol) nHit++;
        }
        return len * (double)nHit / (double)nS;
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

bool edgeUvOnPlane(const TopoDS_Face& f, const TopoDS_Edge& e, gp_Pnt2d& uv0, gp_Pnt2d& uv1,
                   int& stored) {
    stored = 0;
    try {
        Standard_Real a = 0, b = 0;
        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
        if (!pc.IsNull()) {
            stored = 1;
            if (e.Orientation() == TopAbs_REVERSED) std::swap(a, b);
            uv0 = pc->Value(a);
            uv1 = pc->Value(b);
            return true;
        }
        TopLoc_Location loc;
        Handle(Geom_Plane) pl = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(f, loc));
        if (pl.IsNull()) return false;
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        if (c3.IsNull() || l3 <= f3) return false;
        Handle(Geom_Curve) src = c3;
        Handle(Geom_TrimmedCurve) tr = Handle(Geom_TrimmedCurve)::DownCast(c3);
        if (!tr.IsNull() && !tr->BasisCurve().IsNull()) src = tr->BasisCurve();
        Handle(Geom2d_Curve) c2 = GeomAPI::To2d(src, pl->Pln());
        if (c2.IsNull()) return false;
        if (e.Orientation() == TopAbs_REVERSED) std::swap(f3, l3);
        uv0 = c2->Value(f3);
        uv1 = c2->Value(l3);
        stored = 0;
        return true;
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool uvSegIntersect(const gp_Pnt2d& a0, const gp_Pnt2d& a1, const gp_Pnt2d& b0, const gp_Pnt2d& b1,
                    gp_Pnt2d& out) {
    const double dax = a1.X() - a0.X(), day = a1.Y() - a0.Y();
    const double dbx = b1.X() - b0.X(), dby = b1.Y() - b0.Y();
    const double denom = dax * dby - day * dbx;
    if (std::fabs(denom) < 1e-18) return false;
    const double t = ((b0.X() - a0.X()) * dby - (b0.Y() - a0.Y()) * dbx) / denom;
    const double u = ((b0.X() - a0.X()) * day - (b0.Y() - a0.Y()) * dax) / denom;
    if (t < -1e-9 || t > 1.0 + 1e-9 || u < -1e-9 || u > 1.0 + 1e-9) return false;
    out = gp_Pnt2d(a0.X() + t * dax, a0.Y() + t * day);
    return true;
}

int seamCylinderRid(const TopoDS_Edge& e, const RegionSet& rs, const std::vector<ChainGeom>& geom) {
    const void* ts = diagTShapePtr(e);
    gp_Pnt p0, p1;
    TopoDS_Vertex v0, v1;
    TopExp::Vertices(e, v0, v1, Standard_True);
    p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
    p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
    gp_Vec dir(p0, p1);
    const double len = dir.Magnitude();
    if (len > Precision::Confusion()) dir.Normalize();
    gp_Pnt mid(0.5 * (p0.X() + p1.X()), 0.5 * (p0.Y() + p1.Y()), 0.5 * (p0.Z() + p1.Z()));
    for (size_t ci = 0; ci < geom.size(); ci++) {
        for (const TopoDS_Edge& ge : geom[ci].edges) {
            if (ge.IsNull()) continue;
            if (!(e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts))) continue;
            if (ci >= rs.chains.size()) continue;
            const BoundaryChain& ch = rs.chains[ci];
            const Region* a = regionById(rs, ch.regA);
            const Region* b = regionById(rs, ch.regB);
            if (a && a->type == SurfType::Cylinder) return a->id;
            if (b && b->type == SurfType::Cylinder) return b->id;
        }
    }
    int best = -1;
    double bestD = 1e300;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder) continue;
        const gp_Dir ax = r.ax.Direction();
        if (len > Precision::Confusion() && std::fabs(dir.Dot(ax)) < 0.99) continue;
        gp_Vec rho(r.ax.Location(), mid);
        rho -= gp_Vec(ax) * rho.Dot(ax);
        const double dR = std::fabs(rho.Magnitude() - r.radius);
        if (dR < bestD) {
            bestD = dR;
            best = r.id;
        }
    }
    return (best >= 0 && bestD <= 0.5) ? best : -1;
}

void dumpDiagPlateSelfx(int rid, const TopoDS_Face& f, const TopoDS_Wire& ow, const RegionSet& rs,
                        const MeshView& mv, const std::vector<ChainGeom>& geom,
                        const std::vector<char>& collapsed, const std::vector<TopoDS_Edge>& meshE,
                        double sewTol) {
    if (!diagPlatesEnabled() && !diagP2Enabled()) return;
    if (f.IsNull() || ow.IsNull()) return;
    const double tol = std::max(sewTol, Precision::Confusion());
    try {
        auto matchChain = [&](const TopoDS_Edge& e, const gp_Pnt& p0, const gp_Pnt& p1) -> int {
            const void* ts = diagTShapePtr(e);
            for (size_t ci = 0; ci < geom.size(); ci++) {
                for (const TopoDS_Edge& ge : geom[ci].edges) {
                    if (!ge.IsNull() && (e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts)))
                        return (int)ci;
                }
            }
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                for (int eid : rs.chains[ci].meshEdges) {
                    if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull())
                        continue;
                    if (e.IsSame(meshE[(size_t)eid]) ||
                        (ts && diagTShapePtr(meshE[(size_t)eid]) == ts))
                        return (int)ci;
                }
            }
            int best = -1;
            double bestD = 1e300;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                const BoundaryChain& ch = rs.chains[ci];
                if (ch.meshVerts.size() < 2) continue;
                gp_Pnt a = pntOf(mv, ch.meshVerts.front());
                gp_Pnt b = pntOf(mv, ch.meshVerts.back());
                const double d = std::min(p0.Distance(a) + p1.Distance(b),
                                          p0.Distance(b) + p1.Distance(a));
                if (d < bestD) {
                    bestD = d;
                    best = (int)ci;
                }
            }
            return (best >= 0 && bestD <= 2.0 * tol) ? best : -1;
        };

        auto edgeSelfxCls = [&](const TopoDS_Edge& e, int ownerRid) -> std::pair<const char*, int> {
            const void* ts = diagTShapePtr(e);
            if (ts && gSeamTShapes.count(ts)) {
                // D-130: name the chain the analytic edge came from (identity on the
                // collapsed edge TShape) instead of reporting ci=-1 for every seam.
                for (size_t ci = 0; ci < geom.size(); ci++) {
                    for (const auto& ge : geom[ci].edges) {
                        if (ge.IsNull() || diagTShapePtr(ge) != ts) continue;
                        return {"seam", (int)ci};
                    }
                }
                TopoDS_Vertex s0, s1;
                TopExp::Vertices(e, s0, s1, Standard_True);
                const int ciEnd = matchChain(e, s0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(s0),
                                             s1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(s1));
                if (ciEnd < 0) {
                    static thread_local int nUnattrib = 0;
                    if (nUnattrib < 6) {
                        nUnattrib++;
                        int nCollapsed = 0, nGeomE = 0;
                        for (size_t ci = 0; ci < geom.size(); ci++) {
                            if (geom[ci].collapsed) nCollapsed++;
                            nGeomE += (int)geom[ci].edges.size();
                        }
                        TopoDS_Vertex q0, q1;
                        TopExp::Vertices(e, q0, q1, Standard_True);
                        const gp_Pnt a = q0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(q0);
                        const gp_Pnt b = q1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(q1);
                        std::fprintf(stderr,
                                     "DIAG_130_UNATTRIB rid=%d eT=%p inSeamTShapes=%d "
                                     "nChains=%zu nCollapsed=%d nGeomEdges=%d "
                                     "p0=(%.4f,%.4f,%.4f) p1=(%.4f,%.4f,%.4f)\n",
                                     ownerRid, ts, gSeamTShapes.count(ts) ? 1 : 0,
                                     rs.chains.size(), nCollapsed, nGeomE, a.X(), a.Y(), a.Z(),
                                     b.X(), b.Y(), b.Z());
                    }
                }
                return {"seam", ciEnd};
            }
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(e, v0, v1, Standard_True);
            gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
            gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
            const int ci = matchChain(e, p0, p1);
            if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                !geom[(size_t)ci].edges.empty()) {
                const BoundaryChain& ch = rs.chains[(size_t)ci];
                int other = (ch.regA == ownerRid) ? ch.regB : ch.regA;
                if (other < 0) other = (ch.regA == ownerRid) ? ch.regB : ch.regA;
                const Region* rr = regionById(rs, other);
                if (rr && rr->type == SurfType::Cylinder) return {"cyl", other};
                if (rr && rr->type == SurfType::Plane) return {"pln", other};
                return {"span", other};
            }
            if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                const BoundaryChain& ch = rs.chains[(size_t)ci];
                int other = (ch.regA == ownerRid) ? ch.regB : ch.regA;
                if (other < 0) other = ch.regA >= 0 ? ch.regA : ch.regB;
                const Region* rr = regionById(rs, other);
                if (rr && rr->type == SurfType::Cylinder) return {"cyl", other};
                if (rr && rr->type == SurfType::Plane) return {"pln", other};
                return {"fill", other};
            }
            return {"fill", -1};
        };

        auto storedPcOnFace = [&](const TopoDS_Edge& e) -> int {
            TopLoc_Location loc;
            Handle(Geom_Surface) sF = BRep_Tool::Surface(f, loc);
            if (sF.IsNull()) return 0;
            Standard_Boolean hasPc = Standard_False;
            Standard_Real pf = 0, pl = 0;
            (void)BRep_Tool::CurveOnSurface(e, sF, loc, pf, pl, &hasPc);
            return hasPc ? 1 : 0;
        };

        BRepCheck_Analyzer an(f, Standard_True);
        int iw = 0;
        for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), iw++) {
            const TopoDS_Wire w = TopoDS::Wire(wx.Current());
            const TopoDS_Wire& wUse = w.IsSame(ow) ? ow : w;
            Handle(BRepCheck_Result) res = an.Result(w);
            if (!wireStatusOnFace(res, f, BRepCheck_SelfIntersectingWire)) continue;
            BRepCheck_Wire chk(w);
            chk.GeometricControls(Standard_True);
            TopoDS_Edge e1, e2;
            const BRepCheck_Status si = chk.SelfIntersect(f, e1, e2, Standard_False);
            const char* status = brepCheckName((int)si);
            if (!e1.IsNull()) {
                auto c1 = edgeSelfxCls(e1, rid);
                auto c2 = edgeSelfxCls(e2, rid);
                const int i1 = wireEdgeIterIndex(wUse, e1);
                const int i2 = e2.IsNull() ? -1 : wireEdgeIterIndex(wUse, e2);
                gp_Pnt at;
                const char* overlap =
                    (!e2.IsNull() && edgeStretchOverlap(e1, e2, tol, at)) ? "stretch" : "point";
                if (e2.IsNull()) overlap = "point";
                const double len1 = edgeLen3d(e1);
                const double len2 = e2.IsNull() ? 0.0 : edgeLen3d(e2);
                std::fprintf(stderr,
                             "DIAG_SELFX rid=%d iw=%d status=%s e1=%d e2=%d cls1=%s cls2=%s "
                             "other1=%d other2=%d storedPc1=%d storedPc2=%d overlap=%s "
                             "len1=%.4f len2=%.4f p=(%.4f,%.4f,%.4f) shareVtx=%d\n",
                             rid, iw, status, i1, i2, c1.first, e2.IsNull() ? "-" : c2.first,
                             c1.second, e2.IsNull() ? -1 : c2.second, storedPcOnFace(e1),
                             e2.IsNull() ? 0 : storedPcOnFace(e2), overlap, len1, len2, at.X(),
                             at.Y(), at.Z(), e2.IsNull() ? 0 : (edgeShareVtxTShape(e1, e2) ? 1 : 0));
                if (!e2.IsNull()) {
                    auto classifyCi = [&](const TopoDS_Edge& e, int& ci, int& coll, int& partner,
                                          const char*& why) {
                        ci = -1;
                        coll = 0;
                        partner = -1;
                        why = "free";
                        TopoDS_Vertex va, vb;
                        TopExp::Vertices(e, va, vb, Standard_True);
                        gp_Pnt pa = va.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(va);
                        gp_Pnt pb = vb.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(vb);
                        ci = matchChain(e, pa, pb);
                        if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                            const BoundaryChain& ch = rs.chains[(size_t)ci];
                            coll = ((size_t)ci < collapsed.size() && collapsed[(size_t)ci]) ? 1 : 0;
                            partner = (ch.regA == rid) ? ch.regB : ch.regA;
                            if (partner < 0 && (ch.islandA >= 0 || ch.islandB >= 0))
                                why = "island";
                            else
                                why = "chain";
                            return;
                        }
                        const void* ts = diagTShapePtr(e);
                        if (ts && gSeamTShapes.count(ts)) {
                            const int cr = seamCylinderRid(e, rs, geom);
                            why = "seam";
                            partner = cr;
                            return;
                        }
                        why = "free";
                    };
                    int ci1 = -1, ci2 = -1, col1 = 0, col2 = 0, p1 = -1, p2 = -1;
                    const char *w1 = "-", *w2 = "-";
                    classifyCi(e1, ci1, col1, p1, w1);
                    classifyCi(e2, ci2, col2, p2, w2);
                    gp_Pnt2d u10, u11, u20, u21;
                    int st1 = 0, st2 = 0;
                    const bool have1 = edgeUvOnPlane(f, e1, u10, u11, st1);
                    const bool have2 = edgeUvOnPlane(f, e2, u20, u21, st2);
                    const double uvL1 = have1 ? u10.Distance(u11) : 0.0;
                    const double uvL2 = have2 ? u20.Distance(u21) : 0.0;
                    gp_Dir nPl(0, 0, 1);
                    try {
                        TopLoc_Location locN;
                        Handle(Geom_Plane) gpl =
                            Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(f, locN));
                        if (!gpl.IsNull()) nPl = gpl->Pln().Axis().Direction();
                    } catch (const Standard_Failure&) {
                    }
                    const gp_Vec d1 = edgeDir3d(e1), d2 = edgeDir3d(e2);
                    const double dot12 = d1.Dot(d2);
                    const int nShare = sharedVtxTShapeCount(e1, e2);
                    double ang = -1.0;
                    if (nShare > 0 && d1.Magnitude() > 0 && d2.Magnitude() > 0) {
                        double c = std::max(-1.0, std::min(1.0, dot12));
                        ang = std::acos(c) * 180.0 / 3.14159265358979323846;
                    }
                    const int colinear = (std::fabs(dot12) >= 0.99) ? 1 : 0;
                    const int antipar = (dot12 <= -0.99) ? 1 : 0;
                    const double dn1 = d1.Dot(nPl), dn2 = d2.Dot(nPl);
                    const double ovl = overlapLen3dTrue(e1, e2, tol);
                    gp_Pnt ip = (nShare > 0) ? sharedVtxPnt(e1, e2) : gp_Pnt();
                    if (nShare == 0 && have1 && have2) {
                        gp_Pnt2d uxi;
                        if (uvSegIntersect(u10, u11, u20, u21, uxi)) {
                            try {
                                TopLoc_Location locN;
                                Handle(Geom_Plane) gpl =
                                    Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(f, locN));
                                if (!gpl.IsNull()) ip = gpl->Value(uxi.X(), uxi.Y());
                            } catch (const Standard_Failure&) {
                            }
                        }
                    }
                    std::fprintf(stderr,
                                 "DIAG_SELFX2 rid=%d e1=%d e2=%d ci1=%d ci2=%d collapsed1=%d "
                                 "collapsed2=%d partner1=%d partner2=%d "
                                 "uv1=(%.4f,%.4f)->(%.4f,%.4f) uv2=(%.4f,%.4f)->(%.4f,%.4f) "
                                 "uvLen1=%.5f uvLen2=%.5f angleAtSharedVtx=%.2f colinear=%d "
                                 "antiparallel=%d dir1DotN=%.3f dir2DotN=%.3f overlapLen3d=%.5f "
                                 "sharedVtxCount=%d p=(%.4f,%.4f,%.4f) pcSrc1=%s pcSrc2=%s "
                                 "ciWhy1=%s ciWhy2=%s\n",
                                 rid, i1, i2, ci1, ci2, col1, col2, p1, p2, u10.X(), u10.Y(),
                                 u11.X(), u11.Y(), u20.X(), u20.Y(), u21.X(), u21.Y(), uvL1, uvL2,
                                 ang, colinear, antipar, dn1, dn2, ovl, nShare, ip.X(), ip.Y(),
                                 ip.Z(), st1 ? "stored" : "proj", st2 ? "stored" : "proj", w1, w2);
                }
            }
            const int nPcSelfX = countOrientedPcSelfX(f, wUse, tol);
            int nLegacySelfX = 0;
            double areaUV = 0, maxUVgap = 0;
            int nDup = 0;
            wireDiagMetrics(f, wUse, sewTol, areaUV, nLegacySelfX, maxUVgap, nDup);
            std::fprintf(stderr,
                         "DIAG_SELFX_SUM rid=%d iw=%d brepFirst=%d orientedPcSelfX=%d "
                         "legacyIterSelfX=%d note=legacy-skips-adjacent-unoriented\n",
                         rid, iw, e1.IsNull() ? 0 : 1, nPcSelfX, nLegacySelfX);
        }

        if (rid >= 0 && rid <= 3) {
            int nE = 0, pcurveMissing = 0, storedPc = 0, projectedPc = 0;
            TopLoc_Location loc;
            Handle(Geom_Surface) sF = BRep_Tool::Surface(f, loc);
            for (TopoDS_Iterator it(ow); it.More(); it.Next()) {
                if (it.Value().ShapeType() != TopAbs_EDGE) continue;
                nE++;
                const TopoDS_Edge e = TopoDS::Edge(it.Value());
                Standard_Boolean hasPc = Standard_False;
                Standard_Real pf = 0, pl = 0;
                if (!sF.IsNull())
                    (void)BRep_Tool::CurveOnSurface(e, sF, loc, pf, pl, &hasPc);
                if (hasPc)
                    storedPc++;
                else
                    pcurveMissing++;
                if (!hasPc && !sF.IsNull()) {
                    Standard_Real f3 = 0, l3 = 0;
                    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
                    if (!c3.IsNull() && l3 > f3) projectedPc++;
                }
            }
            const int nStretch3d = countStretchPairs3d(ow, sewTol);
            int nLegacySelfX = 0;
            double areaUV = 0, maxUVgap = 0;
            int nDup = 0;
            wireDiagMetrics(f, ow, sewTol, areaUV, nLegacySelfX, maxUVgap, nDup);
            char faceSt[64] = "[]";
            formatStatusList(an.Result(f), faceSt, sizeof(faceSt));
            std::fprintf(stderr,
                         "DIAG_HUBSELFX_SUM rid=%d nE=%d pcurveMissing=%d storedPc=%d "
                         "projectedPc=%d nStretch3d=%d legacySelfX=%d faceSt=%s valid=%d\n",
                         rid, nE, pcurveMissing, storedPc, projectedPc, nStretch3d, nLegacySelfX,
                         faceSt, an.IsValid() ? 1 : 0);
            gp_Dir nPl(0, 0, 1);
            Handle(Geom_Plane) gpl;
            try {
                TopLoc_Location locN;
                gpl = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(f, locN));
                if (!gpl.IsNull()) nPl = gpl->Pln().Axis().Direction();
            } catch (const Standard_Failure&) {
            }
            const double faceTol = f.IsNull() ? sewTol : std::max(BRep_Tool::Tolerance(f), sewTol);
            const char* infPt = "?";
            try {
                BRepTopAdaptor_FClass2d cl(f, Precision::PConfusion());
                infPt = topAbsStateName(cl.PerformInfinitePoint());
            } catch (const Standard_Failure&) {
            }
            int nWires = 0, outerIw = 0, maxNE = -1;
            int nCylKept = 0, nCylFill = 0, nPln = 0, nIsland = 0, nSeamCyl = 0, nFree = 0;
            int iwS = 0;
            for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), iwS++) {
                const TopoDS_Wire w = TopoDS::Wire(wx.Current());
                nWires++;
                int nEw = 0, nSpur = 0, nVertE = 0;
                double uvClosureGap = 0.0;
                gp_Pnt2d firstUv, prevUv;
                bool havePrev = false, haveFirst = false;
                int ei = 0, prevEi = -1, firstEi = -1, gapI = -1, gapJ = -1;
                int lastEi = -1;
                gp_Pnt2d lastUv;
                for (TopoDS_Iterator it(w); it.More(); it.Next(), ei++) {
                    if (it.Value().ShapeType() != TopAbs_EDGE) continue;
                    const TopoDS_Edge e = TopoDS::Edge(it.Value());
                    nEw++;
                    gp_Pnt2d uv0, uv1;
                    int stored = 0;
                    if (edgeUvOnPlane(f, e, uv0, uv1, stored)) {
                        const double uvL = uv0.Distance(uv1);
                        if (uvL < faceTol) nSpur++;
                        if (!haveFirst) {
                            firstUv = uv0;
                            firstEi = ei;
                            haveFirst = true;
                        }
                        if (havePrev) {
                            const double g = prevUv.Distance(uv0);
                            if (g > uvClosureGap) {
                                uvClosureGap = g;
                                gapI = prevEi;
                                gapJ = ei;
                            }
                        }
                        prevUv = uv1;
                        prevEi = ei;
                        lastUv = uv1;
                        lastEi = ei;
                        havePrev = true;
                    }
                    const gp_Vec d = edgeDir3d(e);
                    const double dn = d.Dot(nPl);
                    if (std::fabs(dn) > 0.99) nVertE++;
                    TopoDS_Vertex va, vb;
                    TopExp::Vertices(e, va, vb, Standard_True);
                    gp_Pnt pa = va.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(va);
                    gp_Pnt pb = vb.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(vb);
                    const int ci = matchChain(e, pa, pb);
                    const void* ts = diagTShapePtr(e);
                    const bool isSeam = ts && gSeamTShapes.count(ts);
                    int partner = -1;
                    int island = 0;
                    if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                        const BoundaryChain& ch = rs.chains[(size_t)ci];
                        partner = (ch.regA == rid) ? ch.regB : ch.regA;
                        if (partner < 0 && (ch.islandA >= 0 || ch.islandB >= 0)) island = 1;
                    }
                    const Region* pr = regionById(rs, partner);
                    const char* cls = "free";
                    int cylRid = -1;
                    if (isSeam || std::fabs(dn) > 0.99) {
                        cylRid = seamCylinderRid(e, rs, geom);
                        cls = "seam";
                        nSeamCyl++;
                    } else if (island) {
                        cls = "island";
                        nIsland++;
                    } else if (pr && pr->type == SurfType::Cylinder && isAnalytic(pr)) {
                        cls = "cyl-kept";
                        nCylKept++;
                    } else if (pr && pr->type == SurfType::Cylinder) {
                        cls = "cyl-filled";
                        nCylFill++;
                    } else if (pr && pr->type == SurfType::Plane) {
                        cls = "plane";
                        nPln++;
                    } else if (partner < 0) {
                        cls = "island";
                        nIsland++;
                    } else {
                        nFree++;
                    }
                    if (std::strcmp(cls, "seam") == 0 || std::strcmp(cls, "free") == 0) {
                        std::fprintf(stderr,
                                     "DIAG_HUBEDGE_SEAMFREE rid=%d iw=%d i=%d cls=%s cylRid=%d "
                                     "len3d=%.4f dirDotN=%.3f partner=%d ci=%d\n",
                                     rid, iwS, ei, cls, cylRid, edgeLen3d(e), dn, partner, ci);
                    }
                }
                if (havePrev && haveFirst) {
                    const double g = prevUv.Distance(firstUv);
                    if (g > uvClosureGap) {
                        uvClosureGap = g;
                        gapI = lastEi;
                        gapJ = firstEi;
                    }
                }
                if (nEw > maxNE) {
                    maxNE = nEw;
                    outerIw = iwS;
                }
                const char *c2s = "?", *oris = "?", *sis = "?";
                try {
                    BRepCheck_Wire chk(w);
                    chk.GeometricControls(Standard_True);
                    c2s = brepCheckName((int)chk.Closed2d(f, Standard_False));
                    oris = brepCheckName((int)chk.Orientation(f, Standard_False));
                    TopoDS_Edge se1, se2;
                    sis = brepCheckName((int)chk.SelfIntersect(f, se1, se2, Standard_False));
                } catch (const Standard_Failure&) {
                }
                std::fprintf(stderr,
                             "DIAG_WIRESTAT rid=%d iw=%d closed2d=%s orientation=%s "
                             "selfIntersect=%s nE=%d uvClosureGap=%.5f nSpur=%d nVerticalE=%d "
                             "infPt=%s\n",
                             rid, iwS, c2s, oris, sis, nEw, uvClosureGap, nSpur, nVertE, infPt);
                std::fprintf(stderr,
                             "DIAG_CLOSUREGAP rid=%d iw=%d i=%d j=%d gapUv=%.5f sewTol=%.5f "
                             "nE=%d examined=%d\n",
                             rid, iwS, gapI, gapJ, uvClosureGap, sewTol, nEw, nEw);
            }
            const char* cause = "other";
            try {
                BRepCheck_Face fc(f);
                fc.GeometricControls(Standard_True);
                const BRepCheck_Status iwires = fc.IntersectWires(Standard_False);
                const BRepCheck_Status imb = fc.ClassifyWires(Standard_False);
                const BRepCheck_Status oriW = fc.OrientationOfWires(Standard_False);
                if (iwires != BRepCheck_NoError) cause = "intersecting-wires";
                else if (imb != BRepCheck_NoError) cause = "imbrication";
                else if (oriW != BRepCheck_NoError) cause = "wire-orientation";
                else {
                    int iwC = 0;
                    for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), iwC++) {
                        BRepCheck_Wire chk(TopoDS::Wire(wx.Current()));
                        if (chk.Closed2d(f, Standard_False) != BRepCheck_NoError) {
                            cause = "wire-open-2d";
                            break;
                        }
                        if (chk.Orientation(f, Standard_False) != BRepCheck_NoError) {
                            cause = "wire-orientation";
                            break;
                        }
                    }
                }
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr,
                         "DIAG_FACESTAT rid=%d nWires=%d outerIw=%d st=%s unorientableCause=%s\n",
                         rid, nWires, outerIw, faceSt, cause);
            std::fprintf(stderr,
                         "DIAG_HUBEDGEHIST rid=%d nE=%d cylKept=%d cylFilled=%d plane=%d "
                         "island=%d seamCyl=%d free=%d\n",
                         rid, nE, nCylKept, nCylFill, nPln, nIsland, nSeamCyl, nFree);
            if (rid == 2) {
                int iwB = 0;
                bool named = false;
                for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More() && !named; wx.Next(), iwB++) {
                    const TopoDS_Wire w = TopoDS::Wire(wx.Current());
                    Handle(BRepCheck_Result) wr = an.Result(w);
                    if (wr.IsNull()) continue;
                    if (wireStatusOnFace(wr, f, BRepCheck_BadOrientationOfSubshape) ||
                        wireStatusOnFace(wr, f, BRepCheck_BadOrientation)) {
                        int ei = 0;
                        for (TopoDS_Iterator it(w); it.More(); it.Next(), ei++) {
                            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
                            const TopoDS_Edge e = TopoDS::Edge(it.Value());
                            Handle(BRepCheck_Result) er = an.Result(e);
                            if (er.IsNull()) continue;
                            bool bad = false;
                            for (BRepCheck_ListOfStatus::Iterator sit(er->Status()); sit.More();
                                 sit.Next()) {
                                if (sit.Value() == BRepCheck_BadOrientationOfSubshape ||
                                    sit.Value() == BRepCheck_BadOrientation)
                                    bad = true;
                            }
                            if (bad) {
                                std::fprintf(stderr,
                                             "DIAG_BADORI rid=2 kind=edge iw=%d ei=%d\n", iwB,
                                             ei);
                                named = true;
                                break;
                            }
                        }
                        if (!named) {
                            std::fprintf(stderr, "DIAG_BADORI rid=2 kind=wire iw=%d ei=-1\n",
                                         iwB);
                            named = true;
                        }
                    }
                }
                if (!named)
                    std::fprintf(stderr, "DIAG_BADORI rid=2 kind=face iw=-1 ei=-1\n");
            }
        }
    } catch (const Standard_Failure&) {
    }
}


const char* topAbsOriName(TopAbs_Orientation o) {
    switch (o) {
        case TopAbs_FORWARD: return "FORWARD";
        case TopAbs_REVERSED: return "REVERSED";
        case TopAbs_INTERNAL: return "INTERNAL";
        case TopAbs_EXTERNAL: return "EXTERNAL";
        default: return "?";
    }
}

std::string formatBRepStatusListAll(const BRepCheck_ListOfStatus& lst) {
    std::string s = "[";
    bool any = false;
    for (BRepCheck_ListOfStatus::Iterator it(lst); it.More(); it.Next()) {
        if (any) s += ',';
        any = true;
        const int st = (int)it.Value();
        char buf[96];
        std::snprintf(buf, sizeof(buf), "%d:%s", st, brepCheckName(st));
        s += buf;
    }
    s += ']';
    return s;
}

Handle(BRepCheck_Result) analyzerResultOf(const BRepCheck_Analyzer& an, const TopoDS_Shape& s) {
    try {
        return an.Result(s);
    } catch (const Standard_Failure&) {
        return Handle(BRepCheck_Result)();
    }
}

std::string analyzerStatusOwn(const Handle(BRepCheck_Result)& res) {
    if (res.IsNull()) return "[]";
    try {
        return formatBRepStatusListAll(res->Status());
    } catch (const Standard_Failure&) {
        return "[]";
    }
}

std::string analyzerStatusOnFace(const Handle(BRepCheck_Result)& res, const TopoDS_Face& f) {
    if (res.IsNull()) return "[]";
    try {
        if (res->IsStatusOnShape(f)) return formatBRepStatusListAll(res->StatusOnShape(f));
        res->InitContextIterator();
        for (; res->MoreShapeInContext(); res->NextShapeInContext()) {
            if (res->ContextualShape().IsSame(f) || res->ContextualShape().IsEqual(f))
                return formatBRepStatusListAll(res->StatusOnShape());
        }
    } catch (const Standard_Failure&) {
    }
    return "[]";
}

bool firstNonNoErrorInList(const BRepCheck_ListOfStatus& lst, BRepCheck_Status& st) {
    for (BRepCheck_ListOfStatus::Iterator it(lst); it.More(); it.Next()) {
        if (it.Value() != BRepCheck_NoError) {
            st = it.Value();
            return true;
        }
    }
    return false;
}

bool firstInContextError(const Handle(BRepCheck_Result)& res, const TopoDS_Face& f,
                         BRepCheck_Status& st) {
    if (res.IsNull()) return false;
    try {
        if (res->IsStatusOnShape(f)) return firstNonNoErrorInList(res->StatusOnShape(f), st);
        res->InitContextIterator();
        for (; res->MoreShapeInContext(); res->NextShapeInContext()) {
            if (res->ContextualShape().IsSame(f) || res->ContextualShape().IsEqual(f))
                return firstNonNoErrorInList(res->StatusOnShape(), st);
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

bool firstOwnError(const Handle(BRepCheck_Result)& res, BRepCheck_Status& st) {
    if (res.IsNull()) return false;
    try {
        return firstNonNoErrorInList(res->Status(), st);
    } catch (const Standard_Failure&) {
        return false;
    }
}

// D-S3-99 print-only: walk BRepCheck_Analyzer's result map. No geometry change.
void dumpDiagAnalyzerMap(int rid, char tryAB, const TopoDS_Face& f, const TopoDS_Wire& ow) {
    if (f.IsNull()) return;
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        int nFaceExam = 0, nFacePres = 0, nWireExam = 0, nWirePres = 0;
        int nEdgeExam = 0, nEdgePres = 0, nVertExam = 0, nVertPres = 0;
        char firstSub[64] = "none";
        BRepCheck_Status firstSt = BRepCheck_NoError;
        bool firstMissed = false;
        bool haveFirst = false;

        auto considerFirst = [&](const char* kind, int idx, const Handle(BRepCheck_Result)& res) {
            if (haveFirst) return;
            BRepCheck_Status st = BRepCheck_NoError;
            if (!firstInContextError(res, f, st)) return;
            haveFirst = true;
            firstSt = st;
            if (idx < 0)
                std::snprintf(firstSub, sizeof(firstSub), "%s", kind);
            else
                std::snprintf(firstSub, sizeof(firstSub), "%s %d", kind, idx);
            BRepCheck_Status own = BRepCheck_NoError;
            firstMissed = !firstOwnError(res, own);
        };

        {
            nFaceExam++;
            Handle(BRepCheck_Result) res = analyzerResultOf(an, f);
            if (!res.IsNull()) nFacePres++;
            const std::string own = analyzerStatusOwn(res);
            const std::string ctx = analyzerStatusOnFace(res, f);
            std::fprintf(stderr,
                         "DIAG_ANALYZER rid=%d try=%c sub=face status=%s inContext=%s\n", rid,
                         tryAB, own.c_str(), ctx.c_str());
        }

        TopoDS_Wire outerW;
        int outerIw = -1, maxNE = -1;
        int iw = 0;
        for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next(), iw++) {
            nWireExam++;
            const TopoDS_Wire w = TopoDS::Wire(wx.Current());
            Handle(BRepCheck_Result) res = analyzerResultOf(an, w);
            if (!res.IsNull()) nWirePres++;
            const std::string own = analyzerStatusOwn(res);
            const std::string ctx = analyzerStatusOnFace(res, f);
            std::fprintf(stderr,
                         "DIAG_ANALYZER rid=%d try=%c sub=wire %d status=%s inContext=%s\n", rid,
                         tryAB, iw, own.c_str(), ctx.c_str());
            considerFirst("wire", iw, res);
            int nE = 0;
            for (TopoDS_Iterator it(w); it.More(); it.Next())
                if (it.Value().ShapeType() == TopAbs_EDGE) nE++;
            if (nE > maxNE) {
                maxNE = nE;
                outerW = w;
                outerIw = iw;
            }
        }

        TopTools_IndexedMapOfShape emap, vmap;
        TopExp::MapShapes(f, TopAbs_EDGE, emap);
        TopExp::MapShapes(f, TopAbs_VERTEX, vmap);
        for (int ie = 1; ie <= emap.Extent(); ie++) {
            nEdgeExam++;
            const TopoDS_Shape& e = emap(ie);
            Handle(BRepCheck_Result) res = analyzerResultOf(an, e);
            if (!res.IsNull()) nEdgePres++;
            const std::string own = analyzerStatusOwn(res);
            const std::string ctx = analyzerStatusOnFace(res, f);
            std::fprintf(stderr,
                         "DIAG_ANALYZER rid=%d try=%c sub=edge %d status=%s inContext=%s\n", rid,
                         tryAB, ie - 1, own.c_str(), ctx.c_str());
            considerFirst("edge", ie - 1, res);
        }
        for (int iv = 1; iv <= vmap.Extent(); iv++) {
            nVertExam++;
            const TopoDS_Shape& v = vmap(iv);
            Handle(BRepCheck_Result) res = analyzerResultOf(an, v);
            if (!res.IsNull()) nVertPres++;
            const std::string own = analyzerStatusOwn(res);
            const std::string ctx = analyzerStatusOnFace(res, f);
            std::fprintf(stderr,
                         "DIAG_ANALYZER rid=%d try=%c sub=vertex %d status=%s inContext=%s\n", rid,
                         tryAB, iv - 1, own.c_str(), ctx.c_str());
            considerFirst("vertex", iv - 1, res);
        }

        if (!haveFirst) {
            Handle(BRepCheck_Result) fres = analyzerResultOf(an, f);
            BRepCheck_Status st = BRepCheck_NoError;
            if (firstInContextError(fres, f, st) || firstOwnError(fres, st)) {
                haveFirst = true;
                firstSt = st;
                std::snprintf(firstSub, sizeof(firstSub), "face");
                BRepCheck_Status own = BRepCheck_NoError;
                firstMissed = firstInContextError(fres, f, st) && !firstOwnError(fres, own);
            }
        }

        std::fprintf(stderr,
                     "DIAG_ANALYZER_COV rid=%d try=%c face=%d/%d wire=%d/%d edge=%d/%d "
                     "vertex=%d/%d\n",
                     rid, tryAB, nFacePres, nFaceExam, nWirePres, nWireExam, nEdgePres, nEdgeExam,
                     nVertPres, nVertExam);
        std::fprintf(stderr,
                     "DIAG_ANALYZER_VERDICT rid=%d try=%c first=%s:%s replicaMissed=%d "
                     "missBecause=%s\n",
                     rid, tryAB, haveFirst ? firstSub : "none",
                     haveFirst ? brepCheckName((int)firstSt) : "NoError", firstMissed ? 1 : 0,
                     firstMissed ? "Status() not StatusOnShape(face)" : "Status() already had it");

        if (tryAB == 'A' && !outerW.IsNull()) {
            const char *c2s = "?", *oris = "?";
            try {
                BRepCheck_Wire chk(outerW);
                chk.GeometricControls(Standard_True);
                c2s = brepCheckName((int)chk.Closed2d(f, Standard_False));
                oris = brepCheckName((int)chk.Orientation(f, Standard_False));
            } catch (const Standard_Failure&) {
            }
            const TopAbs_Orientation wOri = outerW.Orientation();
            const TopAbs_Orientation fOri = f.Orientation();
            const TopAbs_Orientation rel = TopAbs::Compose(fOri, wOri);
            (void)ow;
            std::fprintf(stderr,
                         "DIAG_ANALYZER_OUTER rid=%d try=A iw=%d closed2d=%s orientation=%s "
                         "wireOri=%s faceOri=%s relative=%s\n",
                         rid, outerIw, c2s, oris, topAbsOriName(wOri), topAbsOriName(fOri),
                         topAbsOriName(rel));
        }
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_ANALYZER rid=%d try=%c sub=face status=[] inContext=[]\n", rid,
                     tryAB);
    }
}

void dumpDiagAnalyzerBothTries(int rid, const TopoDS_Face& f, const TopoDS_Wire& ow) {
    static thread_local std::unordered_set<int> dumped;
    if (dumped.count(rid)) return;
    dumped.insert(rid);
    dumpDiagAnalyzerMap(rid, 'A', f, ow);
    try {
        const TopoDS_Face fB = TopoDS::Face(f.Reversed());
        dumpDiagAnalyzerMap(rid, 'B', fB, ow.IsNull() ? ow : TopoDS::Wire(ow.Reversed()));
    } catch (const Standard_Failure&) {
    }
}


void diagPlateMakeFaceFail(int rid, const TopoDS_Face& f, const TopoDS_Wire& ow,
                           const RegionSet& rs, const MeshView& mv,
                           const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                           const std::vector<TopoDS_Edge>& meshE, double sewTol) {
    if (!diagPlatesEnabled() && !diagP2Enabled()) return;
    const TopoDS_Shape isoF = diagIsolatedCopy(f);
    const TopoDS_Shape isoW = diagIsolatedCopy(ow);
    if (isoF.IsNull()) return;
    const TopoDS_Face fChk = TopoDS::Face(isoF);
    const TopoDS_Wire wChk = isoW.IsNull() ? TopoDS_Wire() : TopoDS::Wire(isoW);
    try {
        BRepCheck_Analyzer an(fChk, Standard_True);
        std::fprintf(stderr, "DIAG_PLATE rid=%d valid=%d", rid, an.IsValid() ? 1 : 0);
        if (!an.IsValid()) {
            Handle(BRepCheck_Result) res = an.Result(fChk);
            if (!res.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More();
                     it.Next())
                    std::fprintf(stderr, " st=%d", (int)it.Value());
            }
        }
        std::fprintf(stderr, "\n");
        if (!an.IsValid() && !wChk.IsNull())
            dumpDiagPlateSelfx(rid, fChk, wChk, rs, mv, geom, collapsed, meshE, sewTol);
    } catch (const Standard_Failure&) {
    }
}

void emitFaceSurfDiag(int ridHint, const TopoDS_Face& F);
Handle(Geom2d_Curve) makePCurveOnSurf(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                                      const Handle(Geom_Surface)& srf, const char* kind);
thread_local std::unordered_set<const void*> gSeamTShapes;
thread_local int lastCanonClearPass = -1;
thread_local int gCompNTri = 0;

int countSeamEdgesOnWire(const TopoDS_Wire& w) {
    int n = 0;
    if (w.IsNull()) return 0;
    try {
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) {
            const void* ts = diagTShapePtr(TopoDS::Edge(ex.Current()));
            if (ts && gSeamTShapes.count(ts)) n++;
        }
    } catch (const Standard_Failure&) {
    }
    return n;
}

void dumpDiagCylFail(const Region& r, const Loop* outer, const char* stage, const char* reason,
                     const TopoDS_Wire& ow, const Handle(Geom_Surface)& surf,
                     const TopoDS_Face& cand) {
    if (!diagP2Enabled()) return;
    const int nE = countWireEdges(ow);
    const int nChains = outer ? (int)outer->chainIdx.size() : 0;
    const int pcurveMissing = pcurveMissingOnWireSurf(ow, surf);
    const int nSeamE = countSeamEdgesOnWire(ow);
    char faceSt[64] = "[]";
    char wireSt[64] = "[]";
    if (!cand.IsNull()) {
        try {
            const TopoDS_Shape iso = diagIsolatedCopy(cand);
            if (!iso.IsNull()) {
                BRepCheck_Analyzer an(iso, Standard_True);
                formatStatusList(an.Result(iso), faceSt, sizeof(faceSt));
                for (TopExp_Explorer wx(iso, TopAbs_WIRE); wx.More(); wx.Next()) {
                    formatStatusList(an.Result(TopoDS::Wire(wx.Current())), wireSt, sizeof(wireSt));
                    break;
                }
            }
        } catch (const Standard_Failure&) {
        }
    }
    std::fprintf(stderr,
                 "DIAG_CYLFAIL rid=%d R=%.4f nTri=%d closed360=%d nChains=%d nE=%d stage=%s "
                 "reason=%s faceSt=%s wireSt=%s pcurveMissing=%d nSeamE=%d\n",
                 r.id, r.radius, (int)r.tris.size(), r.closed360 ? 1 : 0, nChains, nE, stage,
                 reason, faceSt, wireSt, pcurveMissing, nSeamE);
}

// BRep_Builder::MakeFace + Add keeps the input wire's TShapes (J2).
// BRepBuilderAPI_MakeFace(plane/surf, wire) copies edges — never use for Regions.
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
                    stampTolWriter(eW, "other");
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

// Cylinder Region faces: BRepBuilderAPI_MakeFace on the interned regionSurf handle
// orients the wire (D-S3-43a); planes stay on makeFaceKeep.
bool makeFaceCopy(const Handle(Geom_Surface)& surf, const TopoDS_Wire& outer,
                  const std::vector<TopoDS_Wire>& inners, bool outward, TopoDS_Face& outF) {
    if (surf.IsNull() || outer.IsNull()) return false;
    try {
        BRepBuilderAPI_MakeFace mf(surf, outer, Standard_True);
        if (!mf.IsDone()) return false;
        for (const auto& iw : inners) mf.Add(iw);
        if (!mf.IsDone()) return false;
        outF = mf.Face();
        setFaceOutward(outF, outward);
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
    CylRotTrim = 5,
    ConeBase = 6
};

const char* surfVarName(SurfVar v) {
    switch (v) {
    case SurfVar::Plane: return "Plane";
    case SurfVar::CylBase: return "CylBase";
    case SurfVar::CylRotAx: return "CylRotAx";
    case SurfVar::CylRotU1: return "CylRotU1";
    case SurfVar::CylRectTrim: return "CylRectTrim";
    case SurfVar::CylRotTrim: return "CylRotTrim";
    case SurfVar::ConeBase: return "ConeBase";
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
        } else if (isChamferConeR(r) && v == SurfVar::ConeBase) {
            // The frustum's surface comes from the REGION's own fitted numbers
            // (Location on the R_lo rim, Direction toward R_hi, RefRadius R_lo,
            // semi-angle from the two rim radii and the height) so that the
            // shared rim circles -- which coneIsoCircle builds from the same
            // fields -- lie on it exactly. One handle serves the face and the
            // bind site, exactly as SurfVar::CylBase does for a cylinder.
            gp_Cone cone;
            if (coneFromRims(r.ax, r.radius, chamferRhiOf(r), chamferHeightOf(r), cone))
                s = new Geom_ConicalSurface(cone);
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
    if (isChamferConeR(r)) {
        out[n++] = SurfVar::ConeBase;
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

bool regionSurfOnFace(const Region& r, const Handle(Geom_Surface)& onFace) {
    if (onFace.IsNull()) return false;
    const void* p = onFace.get();
    auto it = gSurfIdent.find(p);
    if (it != gSurfIdent.end() && it->second.first == r.id) return true;
    SurfVar vars[8];
    int nv = 0;
    variantsForRegion(r, vars, nv);
    for (int i = 0; i < nv; i++) {
        Handle(Geom_Surface) reg = regionSurf(r, vars[i]);
        if (!reg.IsNull() && reg.get() == p) return true;
    }
    if (r.type == SurfType::Cylinder) {
        Handle(Geom_Surface) rt = regionSurf(r, SurfVar::CylRectTrim);
        if (!rt.IsNull() && rt.get() == p) return true;
        Handle(Geom_Surface) rot = regionSurf(r, SurfVar::CylRotTrim);
        if (!rot.IsNull() && rot.get() == p) return true;
    }
    return false;
}

void dumpFillIdentity(const TopoDS_Shell& sh, const std::vector<TopoDS_Face>& built,
                      const std::vector<int>& builtRid, const std::vector<char>& eprimeFill,
                      const RegionSet& rs) {
    if (!diagP2Enabled() || sh.IsNull()) return;
    int surfMismatch = 0, nRegionFaces = 0;
    try {
        for (size_t i = 0; i < built.size(); i++) {
            const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
            if (rid < 0) continue;
            if ((size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid]) continue;
            const Region* rr = regionById(rs, rid);
            if (!rr || !isAnalytic(rr)) continue;
            if (rr->type != SurfType::Plane && rr->type != SurfType::Cylinder) continue;
            if (rr->builtAs != BuiltAs::Single && rr->builtAs != BuiltAs::Seamed360 &&
                rr->builtAs != BuiltAs::TwoHalves)
                continue;
            nRegionFaces++;
            Handle(Geom_Surface) s = BRep_Tool::Surface(built[i]);
            if (!regionSurfOnFace(*rr, s)) surfMismatch++;
        }
        int seamOrphan = 0, seamChecked = 0;
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            int idx = faceMap.FindIndex(built[i]);
            if (idx > 0 && i < builtRid.size()) faceRid[(size_t)idx] = builtRid[i];
        }
        auto fillRid = [&](int rid) {
            return rid >= 0 && (size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid];
        };
        for (int ei = 1; ei <= anc.Extent(); ei++) {
            const int nAnc = anc(ei).Extent();
            bool touchFill = false, touchAnalytic = false;
            for (TopTools_ListOfShape::Iterator it(anc(ei)); it.More(); it.Next()) {
                const TopoDS_Face& fa = TopoDS::Face(it.Value());
                int ia = faceMap.FindIndex(fa);
                const int rid = (ia > 0) ? faceRid[(size_t)ia] : -1;
                if (fillRid(rid)) touchFill = true;
                const Region* rr = regionById(rs, rid);
                if (rr && isAnalytic(rr) &&
                    (rr->type == SurfType::Plane || rr->type == SurfType::Cylinder))
                    touchAnalytic = true;
            }
            if (!touchFill || !touchAnalytic) continue;
            seamChecked++;
            if (nAnc < 2) seamOrphan++;
        }
        std::fprintf(stderr,
                     "DIAG_FILLID surfMismatch=%d nRegionFaces=%d seamOrphan=%d seamChecked=%d\n",
                     surfMismatch, nRegionFaces, seamOrphan, seamChecked);
    } catch (const Standard_Failure&) {
    }
}

int countShellFreeEdges(const TopoDS_Shell& sh) {
    if (sh.IsNull()) return 0;
    int freeE = 0;
    try {
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, anc);
        for (int i = 1; i <= anc.Extent(); i++)
            if (anc(i).Extent() < 2) freeE++;
    } catch (const Standard_Failure&) {
    }
    return freeE;
}

int countAdmittedRegions(const RegionSet& rs, const std::vector<char>& eprimeFill) {
    int n = 0;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane && r.type != SurfType::Cylinder) continue;
        if (r.id < 0) continue;
        if ((size_t)r.id < eprimeFill.size() && eprimeFill[(size_t)r.id]) continue;
        if (!isAnalytic(&r)) continue;
        n++;
    }
    return n;
}

bool regionBuiltAnalytic(int rid, const RegionSet& rs, const std::vector<char>& eprimeFill,
                         const std::vector<char>& exploded) {
    if (rid < 0) return false;
    if ((size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid]) return false;
    if (regionExploded(exploded, rid)) return false;
    const Region* r = regionById(rs, rid);
    if (!r) return false;
    return r->builtAs == BuiltAs::Single || r->builtAs == BuiltAs::Seamed360 ||
           r->builtAs == BuiltAs::TwoHalves;
}

bool regionClosureHealEligible(int rid, const RegionSet& rs, const std::vector<char>& eprimeFill,
                               const std::vector<char>& exploded) {
    if (!regionBuiltAnalytic(rid, rs, eprimeFill, exploded)) return false;
    const Region* r = regionById(rs, rid);
    return r && r->builtAs == BuiltAs::Single;
}

int matchFreeEdgeChain(const TopoDS_Edge& e, const gp_Pnt& p0, const gp_Pnt& p1,
                       const MeshView& mv, const RegionSet& rs,
                       const std::vector<TopoDS_Edge>& meshE, const std::vector<ChainGeom>& geom,
                       double tol) {
    const void* ts = diagTShapePtr(e);
    for (size_t ci = 0; ci < geom.size(); ci++) {
        for (const TopoDS_Edge& ge : geom[ci].edges) {
            if (!ge.IsNull() && (e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts))) return (int)ci;
        }
    }
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        for (int eid : rs.chains[ci].meshEdges) {
            if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
            if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts))
                return (int)ci;
        }
    }
    static const std::vector<char> kMatchAllChains;
    return matchCollapsedChainToSegment(mv, rs, kMatchAllChains, p0, p1, tol);
}

struct ClosureHealAttrib {
    int rid = -1;
    const char* from = nullptr;
};

std::vector<ClosureHealAttrib> collectClosureHealDemotions(
    const TopoDS_Shell& sh, const std::vector<TopoDS_Face>& built,
    const std::vector<int>& builtRid, const std::vector<char>& eprimeFill, const RegionSet& rs,
    const MeshView& mv, const std::vector<TopoDS_Edge>& meshE, const std::vector<ChainGeom>& geom,
    const std::vector<char>& collapsed, const std::vector<char>& exploded, double sewTol) {
    std::vector<ClosureHealAttrib> out;
    if (sh.IsNull()) return out;
    const double tol = std::max(sewTol, Precision::Confusion());
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
        std::unordered_set<int> seen;
        for (int ei = 1; ei <= anc.Extent(); ei++) {
            if (anc(ei).Extent() >= 2) continue;
            const TopoDS_Edge e = TopoDS::Edge(anc.FindKey(ei));
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(e, v0, v1, Standard_True);
            gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
            gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
            int ownerRid = -1;
            if (anc(ei).Extent() == 1) {
                const TopoDS_Face f = TopoDS::Face(anc(ei).First());
                int fi = faceMap.FindIndex(f);
                if (fi > 0) ownerRid = faceRid[(size_t)fi];
            }
            const int ci = matchFreeEdgeChain(e, p0, p1, mv, rs, meshE, geom, tol);
            int otherRid = -1;
            if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                const BoundaryChain& ch = rs.chains[(size_t)ci];
                if (ch.regA == ownerRid)
                    otherRid = ch.regB;
                else if (ch.regB == ownerRid)
                    otherRid = ch.regA;
                else if (ownerRid < 0)
                    otherRid = ch.regA >= 0 ? ch.regA : ch.regB;
            }
            int healRid = -1;
            const char* from = nullptr;
            if (regionClosureHealEligible(ownerRid, rs, eprimeFill, exploded)) {
                healRid = ownerRid;
                from = "owner";
            } else if (regionClosureHealEligible(otherRid, rs, eprimeFill, exploded)) {
                healRid = otherRid;
                from = "neighbour";
            }
            if (healRid < 0 || seen.count(healRid)) continue;
            seen.insert(healRid);
            out.push_back({healRid, from});
        }
    } catch (const Standard_Failure&) {
    }
    return out;
}

const char* freeEdgeRegionKind(int rid, const RegionSet& rs, const std::vector<char>& eprimeFill) {
    if (rid < 0) return "none";
    if ((size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid]) return "fill";
    const Region* r = regionById(rs, rid);
    if (!r) return "none";
    if (r->builtAs == BuiltAs::ExplodedToFacets) return "facet";
    if (r->type == SurfType::Plane) return "plane";
    if (r->type == SurfType::Cylinder) return "cyl";
    return "facet";
}

bool edgeEndsMatch(const gp_Pnt& a0, const gp_Pnt& a1, const gp_Pnt& b0, const gp_Pnt& b1,
                   double tol) {
    return (a0.Distance(b0) <= tol && a1.Distance(b1) <= tol) ||
           (a0.Distance(b1) <= tol && a1.Distance(b0) <= tol);
}

void dumpDiagFreeEdges(const TopoDS_Shell& sh, const std::vector<TopoDS_Face>& built,
                       const std::vector<int>& builtRid, const std::vector<char>& eprimeFill,
                       const RegionSet& rs, const MeshView& mv,
                       const std::vector<TopoDS_Edge>& meshE, const std::vector<ChainGeom>& geom,
                       const std::vector<char>& collapsed, double sewTol) {
    if (!diagP2Enabled() || sh.IsNull()) return;
    const double tol = std::max(sewTol, Precision::Confusion());
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

        auto kindOf = [&](int rid) { return freeEdgeRegionKind(rid, rs, eprimeFill); };
        auto isFill = [&](int rid) {
            return rid >= 0 && (size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid];
        };
        auto isDemoted = [&](int rid) {
            return rid >= 0 && (size_t)rid < gEPrimeDemoted.size() && gEPrimeDemoted[(size_t)rid];
        };

        auto faceHasEnds = [&](int rid, const gp_Pnt& p0, const gp_Pnt& p1) -> int {
            if (rid < 0) return 0;
            for (size_t i = 0; i < built.size(); i++) {
                if ((i < builtRid.size() ? builtRid[i] : -1) != rid) continue;
                for (TopExp_Explorer ex(built[i], TopAbs_EDGE); ex.More(); ex.Next()) {
                    const TopoDS_Edge e = TopoDS::Edge(ex.Current());
                    TopoDS_Vertex va, vb;
                    TopExp::Vertices(e, va, vb, Standard_True);
                    if (va.IsNull() || vb.IsNull()) continue;
                    if (edgeEndsMatch(p0, p1, BRep_Tool::Pnt(va), BRep_Tool::Pnt(vb), tol))
                        return 1;
                }
            }
            return 0;
        };

        auto nTriAlongChain = [&](int ci, int otherRid) -> int {
            if (ci < 0 || (size_t)ci >= rs.chains.size() || otherRid < 0) return 0;
            const BoundaryChain& ch = rs.chains[(size_t)ci];
            if (!mv.triEdges) {
                const Region* o = regionById(rs, otherRid);
                return o ? (int)o->tris.size() : 0;
            }
            std::unordered_set<int> seen;
            for (int eid : ch.meshEdges) {
                if (eid < 0) continue;
                for (size_t k = 0; k < mv.nTri; k++) {
                    bool hit = false;
                    for (int s = 0; s < 3; s++)
                        if (mv.triEdges[k][s] == eid) hit = true;
                    if (!hit) continue;
                    const int tr =
                        ((size_t)k < rs.triRegion.size()) ? rs.triRegion[(size_t)k] : -1;
                    if (tr == otherRid) seen.insert((int)k);
                }
            }
            return (int)seen.size();
        };

        auto matchChain = [&](const TopoDS_Edge& e, const gp_Pnt& p0, const gp_Pnt& p1) -> int {
            return matchFreeEdgeChain(e, p0, p1, mv, rs, meshE, geom, tol);
        };

        int iFree = 0;
        int nSpanFill = 0, nMeshFill = 0, nOther = 0;
        for (int ei = 1; ei <= anc.Extent(); ei++) {
            if (anc(ei).Extent() >= 2) continue;
            const TopoDS_Edge e = TopoDS::Edge(anc.FindKey(ei));
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(e, v0, v1, Standard_True);
            gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
            gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
            const double len = p0.Distance(p1);
            int ownerRid = -1;
            if (anc(ei).Extent() == 1) {
                const TopoDS_Face f = TopoDS::Face(anc(ei).First());
                int fi = faceMap.FindIndex(f);
                if (fi > 0) ownerRid = faceRid[(size_t)fi];
            }
            const Region* owner = regionById(rs, ownerRid);
            const char* ownerKind = kindOf(ownerRid);
            const char* ownerBuilt = owner ? builtAsName(owner->builtAs) : "NotBuilt";
            const void* ts = diagTShapePtr(e);
            const char* edgeKind = "meshE";
            const int ci = matchChain(e, p0, p1);
            if (ts && gSeamTShapes.count(ts))
                edgeKind = "seam";
            else if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                     (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                     !geom[(size_t)ci].edges.empty())
                edgeKind = "span";
            int otherRid = -1;
            if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                const BoundaryChain& ch = rs.chains[(size_t)ci];
                if (ch.regA == ownerRid)
                    otherRid = ch.regB;
                else if (ch.regB == ownerRid)
                    otherRid = ch.regA;
                else if (ownerRid < 0)
                    otherRid = ch.regA >= 0 ? ch.regA : ch.regB;
                else
                    otherRid = (ch.regA == ownerRid) ? ch.regB : ch.regA;
            }
            const char* otherKind = kindOf(otherRid);
            const int otherHas = faceHasEnds(otherRid, p0, p1);
            const int nAlong = nTriAlongChain(ci, otherRid >= 0 ? otherRid : ownerRid);
            if (std::strcmp(edgeKind, "span") == 0 && std::strcmp(otherKind, "fill") == 0)
                nSpanFill++;
            else if (std::strcmp(edgeKind, "meshE") == 0 && std::strcmp(otherKind, "fill") == 0)
                nMeshFill++;
            else
                nOther++;
            std::fprintf(stderr,
                         "DIAG_FREEEDGE i=%d ownerRid=%d ownerKind=%s ownerBuiltAs=%s "
                         "edgeKind=%s otherRid=%d otherKind=%s otherHasEdge=%d nTriAlong=%d "
                         "len=%.4f p0=(%.4f,%.4f,%.4f) p1=(%.4f,%.4f,%.4f)\n",
                         iFree, ownerRid, ownerKind, ownerBuilt, edgeKind, otherRid, otherKind,
                         otherHas, nAlong, len, p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z());
            iFree++;
        }
        std::fprintf(stderr,
                     "DIAG_FREEEDGE_CLS n=%d span-vs-fill=%d meshE-vs-fill=%d other=%d\n",
                     iFree, nSpanFill, nMeshFill, nOther);

        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            const BoundaryChain& ch = rs.chains[ci];
            int hub = -1, other = -1;
            if (ch.regA == 1) {
                hub = 1;
                other = ch.regB;
            } else if (ch.regB == 1) {
                hub = 1;
                other = ch.regA;
            }
            if (hub != 1) continue;
            const char* ek = "meshE";
            if ((size_t)ci < collapsed.size() && collapsed[ci] && (size_t)ci < geom.size() &&
                geom[ci].collapsed && !geom[ci].edges.empty())
                ek = "span";
            const void* seamTs = nullptr;
            if ((size_t)ci < geom.size()) {
                for (const auto& ge : geom[ci].edges) {
                    const void* t = diagTShapePtr(ge);
                    if (t && gSeamTShapes.count(t)) {
                        ek = "seam";
                        seamTs = t;
                        break;
                    }
                }
            }
            (void)seamTs;
            std::fprintf(stderr,
                         "DIAG_HUB1NBR otherRid=%d otherKind=%s demoted=%d fill=%d "
                         "edgeKind=%s collapsed=%d nMeshE=%zu nTriOther=%d\n",
                         other, kindOf(other), isDemoted(other) ? 1 : 0, isFill(other) ? 1 : 0,
                         ek, ((size_t)ci < collapsed.size() && collapsed[ci]) ? 1 : 0,
                         ch.meshEdges.size(),
                         regionById(rs, other) ? (int)regionById(rs, other)->tris.size() : 0);
        }
    } catch (const Standard_Failure&) {
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

Handle(Geom_ConicalSurface) basisConeOf(const Handle(Geom_Surface)& srf) {
    Handle(Geom_ConicalSurface) c = Handle(Geom_ConicalSurface)::DownCast(srf);
    if (!c.IsNull()) return c;
    return Handle(Geom_ConicalSurface)::DownCast(untrimmedBasisOf(srf));
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
        const bool isCirc = (kind && std::strcmp(kind, "circ") == 0) ||
                            (!src.IsNull() && src->DynamicType() == STANDARD_TYPE(Geom_Circle));
        const bool isElips = (kind && std::strcmp(kind, "elips") == 0) ||
                             (!src.IsNull() && src->DynamicType() == STANDARD_TYPE(Geom_Ellipse));
        if (isCirc || isElips) {
            try {
                return GeomAPI::To2d(src, pl->Pln());
            } catch (const Standard_Failure&) {
                return {};
            }
        }
        try {
            const gp_Pnt p0 = c3->Value(f);
            const gp_Pnt p1 = c3->Value(l);
            Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0;
            ElSLib::Parameters(pl->Pln(), p0, u0, v0);
            ElSLib::Parameters(pl->Pln(), p1, u1, v1);
            const double df = l - f;
            if (df <= Precision::PConfusion()) return {};
            gp_Vec2d d2(u1 - u0, v1 - v0);
            if (d2.SquareMagnitude() < Precision::PConfusion() * Precision::PConfusion()) {
                return GeomAPI::To2d(src, pl->Pln());
            }
            d2 /= df;
            const gp_Pnt2d loc(u0 - f * d2.X(), v0 - f * d2.Y());
            return new Geom2d_Line(loc, gp_Dir2d(d2));
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
    Handle(Geom_ConicalSurface) cone = basisConeOf(srf);
    if (!cone.IsNull()) {
        // Cone pcurves are exact by construction: OCCT's cone is unit speed in v
        // and angular in u, so a rim circle and a generator are both straight 2d
        // lines. refit_cone_bind builds them and exactMaxAtBind certifies them.
        try {
            const gp_Cone gcone = cone->Cone();
            if (isCirc && !isElips) {
                Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(src);
                if (gc.IsNull()) return {};
                return conePCurveForCircle(gcone, gc->Circ());
            }
            if (isLin) {
                Handle(Geom_Line) gl = Handle(Geom_Line)::DownCast(src);
                if (gl.IsNull()) return {};
                return conePCurveForLine(gcone, gl->Lin(), 0.5 * (f + l));
            }
        } catch (const Standard_Failure&) {
            return {};
        }
        return {};
    }
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

// D-S3-64 / D-S3-70: plane pcurve endpoints from the edge's 3D vertices (the
// shared TShape points), not c3->Value(f/l). SameParameter drift between the
// 3D curve range and the vertex makes two edges that 3D-share a vertex miss
// in UV by tens of microns. Circ/elips stay on To2d (an arc is not a chord).
Handle(Geom2d_Curve) makePlanePCurvePinnedToVertices(const TopoDS_Edge& e,
                                                     const Handle(Geom_Plane)& pl,
                                                     const Handle(Geom_Curve)& c3,
                                                     Standard_Real f, Standard_Real l) {
    if (pl.IsNull() || c3.IsNull() || e.IsNull()) return {};
    try {
        TopoDS_Vertex vf, vl;
        TopExp::Vertices(e, vf, vl, Standard_False);
        const gp_Pnt p0 = vf.IsNull() ? c3->Value(f) : BRep_Tool::Pnt(vf);
        const gp_Pnt p1 = vl.IsNull() ? c3->Value(l) : BRep_Tool::Pnt(vl);
        Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0;
        ElSLib::Parameters(pl->Pln(), p0, u0, v0);
        ElSLib::Parameters(pl->Pln(), p1, u1, v1);
        const double df = l - f;
        if (df <= Precision::PConfusion()) return {};
        gp_Vec2d d2(u1 - u0, v1 - v0);
        if (d2.SquareMagnitude() < Precision::PConfusion() * Precision::PConfusion()) {
            Handle(Geom_Curve) src = basisCurveOf(c3);
            return GeomAPI::To2d(src, pl->Pln());
        }
        d2 /= df;
        const gp_Pnt2d loc(u0 - f * d2.X(), v0 - f * d2.Y());
        return new Geom2d_Line(loc, gp_Dir2d(d2));
    } catch (const Standard_Failure&) {
        return {};
    }
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

// --- TL1: exact-deviation edge tolerance at bindTolFromDev (D-S3-111..113) ---
thread_local int gBindTolUnhandled = 0;
thread_local int gBindTolOverCap = 0;
thread_local int gBindTolOverCapGenuine = 0;
thread_local int gBindTolBrepGtExact = 0;
thread_local int gBindTolBrepLeExact = 0;
thread_local int gCylSideTolChanged = 0;
thread_local int gReconnectFire = 0;
thread_local int gTolRewriteFire = 0;
thread_local std::unordered_set<const void*> gBindTolUnhandledTShapes;
thread_local std::unordered_map<std::string, std::unordered_set<const void*>> gBindTolUnhandledByClass;
thread_local std::unordered_set<const void*> gBindTolOverCapTShapes;
thread_local std::unordered_map<const void*, const char*> gLastTolWriter;
thread_local std::unordered_map<const void*, double> gAllPlateEdgeTol;
thread_local int gTolCheckRows = 0;
thread_local int gTolCheckBrepLeExact = 0;
thread_local int gTolCheckPlaneRows = 0;
thread_local int gTolCheckPlaneBrepLe = 0;
thread_local int gTolCheckCylRows = 0;
thread_local int gTolCheckCylBrepLe = 0;
thread_local int gIncontextPlateFlagged = 0;
thread_local int gIncontextPlateExamined = 0;
thread_local int gIncontextShippedFlagged = 0;
thread_local int gIncontextShippedExamined = 0;

double distToAxis(const gp_Ax1& ax, const gp_Pnt& p) {
    gp_Vec v(ax.Location(), p);
    gp_Vec ad(ax.Direction());
    return (v - ad * v.Dot(ad)).Magnitude();
}

double planeSignedDistBind(const gp_Pln& pl, const gp_Pnt& p) {
    Standard_Real a = 0, b = 0, c = 0, d = 0;
    pl.Coefficients(a, b, c, d);
    return a * p.X() + b * p.Y() + c * p.Z() + d;
}

double curveSurfDevAtBind(const Handle(Geom_Curve)& c3, Standard_Real t,
                          const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                          const TopLoc_Location& loc) {
    const gp_Pnt P = c3->Value(t);
    const gp_Pnt2d uv = c2d->Value(t);
    gp_Pnt Q = srf->Value(uv.X(), uv.Y());
    if (!loc.IsIdentity()) Q.Transform(loc.Transformation());
    return P.Distance(Q);
}

double brepDevMaxOnSurf(const Handle(Geom_Curve)& c3, Standard_Real f3, Standard_Real l3,
                        const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& pc,
                        Standard_Real a, Standard_Real b, const TopLoc_Location& loc);

// D-S3-111 circle-on-plane: |offset| + R·sin(tilt) + R·(1−cos(tilt)) at the same
// extremal parameter as the normal-offset maximum (To2d ellipse in-plane term).
double circleOnPlaneMax(const gp_Pln& pl, const gp_Circ& circ) {
    const double off = std::fabs(planeSignedDistBind(pl, circ.Location()));
    double ang = circ.Axis().Direction().Angle(pl.Axis().Direction());
    ang = std::min(ang, kPi - ang);  // anti-parallel is still coplanar (φ=0, not 2R)
    const double R = circ.Radius();
    return off + R * std::sin(ang) + R * (1.0 - std::cos(ang));
}

// Circle on cylinder: constant radial offset |dist_to_axis − R| when the circle axis is
// parallel to the cylinder axis; otherwise the radial span [d−r, d+r] vs R.
double circleOnCylMax(const Handle(Geom_CylindricalSurface)& cyl, const gp_Circ& circ,
                      const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                      const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                      const TopLoc_Location& loc) {
    const gp_Ax1 cax = cyl->Axis();
    const double Rcyl = cyl->Radius();
    double m = std::max(curveSurfDevAtBind(c3, f, srf, c2d, loc),
                        curveSurfDevAtBind(c3, l, srf, c2d, loc));
    const double d0 = distToAxis(cax, circ.Location());
    const double r = circ.Radius();
    const double parallel = std::fabs(circ.Axis().Direction().Dot(cax.Direction()));
    if (parallel > 1.0 - 1e-9) {
        // Cap circle (center on axis): constant ring distance r. Latitude circle (center on
        // surface): center distance d0 ≈ R. Constant radial offset |d_ring − R|.
        const double dRing = (d0 <= Precision::Confusion()) ? r : d0;
        m = std::max(m, std::fabs(dRing - Rcyl));
    }
    return m;
}

void noteBindTolUnhandled(const TopoDS_Edge& e, const char* cls) {
    const void* ts = diagTShapePtr(e);
    if (!ts) return;
    if (gBindTolUnhandledTShapes.insert(ts).second) {
        gBindTolUnhandled++;
        if (cls) gBindTolUnhandledByClass[cls].insert(ts);
    }
}

// Closed form supremum |C(t)-S(P(t))| at bind site. Returns -1 if unhandled.
double exactMaxAtBind(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                      const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                      const TopLoc_Location& loc, const char** clsOut) {
    if (clsOut) *clsOut = "unhandled";
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return -1.0;
    Handle(Geom_Plane) gpl = Handle(Geom_Plane)::DownCast(srf);
    Handle(Geom_CylindricalSurface) gcyl = basisCylOf(srf);
    // D-130-3: a cone edge is certified by refit_cone_bind, which composes
    // 130-CONE-MATH's surface supremum with the parametrisation term the bind
    // site actually records. It certifies or it refuses; there is no third
    // outcome, and a refusal is tier 2 (counted, per D-130-2).
    if (!basisConeOf(srf).IsNull())
        return coneBindSup(c3, f, l, srf, c2d, loc, clsOut);
    if (gpl.IsNull() && gcyl.IsNull()) {
        if (clsOut) *clsOut = "unhandled-surface";
        return -1.0;
    }
    Handle(Geom_Curve) src = basisCurveOf(c3);
    if (src.IsNull()) src = c3;
    if (src->DynamicType() == STANDARD_TYPE(Geom_Line)) {
        if (clsOut) *clsOut = "affine";
        return std::max(curveSurfDevAtBind(c3, f, srf, c2d, loc),
                        curveSurfDevAtBind(c3, l, srf, c2d, loc));
    }
    if (src->DynamicType() == STANDARD_TYPE(Geom_Circle)) {
        Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(src);
        if (!gc.IsNull()) {
            if (!gpl.IsNull()) {
                gp_Pln pl = gpl->Pln();
                if (!loc.IsIdentity()) pl.Transform(loc.Transformation());
                if (clsOut) *clsOut = "circle-on-plane";
                return circleOnPlaneMax(pl, gc->Circ());
            }
            if (!gcyl.IsNull()) {
                if (clsOut) *clsOut = "circle-on-cylinder";
                return circleOnCylMax(gcyl, gc->Circ(), c3, f, l, srf, c2d, loc);
            }
        }
    }
    Handle(Geom_BSplineCurve) bs = Handle(Geom_BSplineCurve)::DownCast(src);
    if (!bs.IsNull() && bs->Degree() <= 1) {
        if (!gpl.IsNull()) {
            gp_Pln pl = gpl->Pln();
            if (!loc.IsIdentity()) pl.Transform(loc.Transformation());
            double m = std::max(curveSurfDevAtBind(c3, f, srf, c2d, loc),
                                curveSurfDevAtBind(c3, l, srf, c2d, loc));
            for (int i = 1; i <= bs->NbPoles(); i++) {
                const double t = (bs->NbPoles() > 1)
                                     ? f + (l - f) * (double)(i - 1) / (double)(bs->NbPoles() - 1)
                                     : f;
                m = std::max(m, curveSurfDevAtBind(c3, t, srf, c2d, loc));
            }
            if (clsOut) *clsOut = "polyline";
            return m;
        }
        if (!gcyl.IsNull()) {
            if (clsOut) *clsOut = "unhandled-polyline-cyl";
            return -1.0;
        }
    }
    if (clsOut) *clsOut = "unhandled-other";
    return -1.0;
}

double dense4096MaxOnPlane(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                           const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                           const TopLoc_Location& loc) {
    double m = 0.0;
    for (int i = 0; i < 4096; i++) {
        const double u = (double)i / 4095.0;
        const double t = f + u * (l - f);
        m = std::max(m, curveSurfDevAtBind(c3, t, srf, c2d, loc));
    }
    return m;
}

double bindTolFromDevLegacy(double dev, double cap) {
    double t = Precision::Confusion();
    if (std::isfinite(dev) && dev > t) t = dev;
    if (t > cap) t = cap;
    return t;
}

void stampTolWriter(const TopoDS_Edge& e, const char* site) {
    const void* ts = diagTShapePtr(e);
    if (ts && site) gLastTolWriter[ts] = site;
}

void fireTolRewriteEdge(BRep_Builder& B, const TopoDS_Edge& e, double after, const char* site,
                        double meshCap) {
    const double before = BRep_Tool::Tolerance(e);
    if (after <= before + Precision::PConfusion()) return;
    B.UpdateEdge(e, after);
    stampTolWriter(e, site);
    gTolRewriteFire++;
    if (diagP2Enabled())
        std::fprintf(stderr,
                     "DIAG_TOLREWRITE site=%s before=%.9g after=%.9g meshCap=%.9g overCap=%d\n",
                     site ? site : "?", before, after, meshCap,
                     (meshCap > 0.0 && after > meshCap + Precision::PConfusion()) ? 1 : 0);
}

const char* lastTolWriterOf(const TopoDS_Edge& e) {
    const void* ts = diagTShapePtr(e);
    if (!ts) return "unknown";
    auto it = gLastTolWriter.find(ts);
    return (it == gLastTolWriter.end() || !it->second) ? "unknown" : it->second;
}

double exactOrSampledOnPair(const Handle(Geom_Curve)& c3, Standard_Real f, Standard_Real l,
                            const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                            const TopLoc_Location& loc) {
    const char* cls = nullptr;
    const double ex = exactMaxAtBind(c3, f, l, srf, c2d, loc, &cls);
    if (ex >= 0.0 && std::isfinite(ex)) return ex;
    return brepDevMaxOnSurf(c3, f, l, srf, c2d, f, l, loc);
}

// Max exact deviation over every CurveOnSurface the edge already carries, plus extraThis
// (the pcurve about to be bound). Unhandled classes use the 23-sample grid (named residue).
double maxExactAllPcurves(const TopoDS_Edge& e, double extraThis) {
    double m = extraThis;
    try {
        Standard_Real f = 0, l = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
        Handle(BRep_TEdge) te = Handle(BRep_TEdge)::DownCast(e.TShape());
        if (te.IsNull() || c3.IsNull()) return m;
        for (BRep_ListOfCurveRepresentation::Iterator it(te->Curves()); it.More(); it.Next()) {
            Handle(BRep_CurveOnSurface) cos = Handle(BRep_CurveOnSurface)::DownCast(it.Value());
            if (cos.IsNull() || cos->Surface().IsNull() || cos->PCurve().IsNull()) continue;
            const double ex =
                exactOrSampledOnPair(c3, f, l, cos->Surface(), cos->PCurve(), cos->Location());
            if (ex > m) m = ex;
        }
    } catch (const Standard_Failure&) {
    }
    return m;
}

// D-S3-124: one scalar must strictly exceed max over ALL pcurves. Inherited equality
// (existing <= maxAll) is REPLACED. Fat inherited is never lowered. cap is guardrail 3.
bool bindTolFromExact(double exactMax, double sampled, double existing, double cap, double* tolOut,
                      bool* overCapOut) {
    if (overCapOut) *overCapOut = false;
    const double thisDev = (exactMax >= 0.0 && std::isfinite(exactMax)) ? exactMax : sampled;
    if (thisDev > cap) {
        if (overCapOut) *overCapOut = true;
        if (tolOut) *tolOut = existing;
        return false;
    }
    const double want = std::max(Precision::Confusion(), thisDev + Precision::PConfusion());
    double tol;
    if (existing <= thisDev)
        tol = want;
    else
        tol = std::max(existing, want);
    if (tolOut) *tolOut = tol;
    return true;
}

// Apply D-S3-124: declared tol exceeds max over all pcurves. Guardrail 3 applies to
// THIS pcurve's deviation vs this surface's cap — not to a sibling pcurve's exact.
bool applyBindTolAllPcurves(const TopoDS_Edge& e, double exactThis, double sampled, double cap,
                            double* tolOut, bool* overCapOut, double* maxAllOut) {
    const double thisDev = (exactThis >= 0.0 && std::isfinite(exactThis)) ? exactThis : sampled;
    const double existing = BRep_Tool::Tolerance(e);
    const double others = maxExactAllPcurves(e, -1.0);
    const double maxAll = std::max(thisDev, others);
    if (maxAllOut) *maxAllOut = maxAll;
    if (overCapOut) *overCapOut = false;
    if (thisDev > cap) {
        if (overCapOut) *overCapOut = true;
        // Do not raise to thisDev. Still allow D-S3-124 against other pcurves.
        const double wantOther =
            (others > 0.0) ? std::max(Precision::Confusion(), others + Precision::PConfusion())
                            : Precision::Confusion();
        double tol = existing;
        if (existing <= others) tol = std::max(tol, wantOther);
        else tol = std::max(existing, wantOther);
        if (tolOut) *tolOut = tol;
        return false;
    }
    return bindTolFromExact(maxAll, sampled, existing, 1e300, tolOut, overCapOut);
}

double bindTolFromDev(double dev, double cap) {
    return bindTolFromDevLegacy(dev, cap);
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
        Handle(Geom2d_Curve) c2d;
        Handle(Geom_Plane) plBind = Handle(Geom_Plane)::DownCast(srf);
        const bool pinVerts = !plBind.IsNull() && kind &&
                              (std::strcmp(kind, "lin") == 0 || std::strcmp(kind, "poly") == 0);
        if (pinVerts)
            c2d = makePlanePCurvePinnedToVertices(e, plBind, c3, f, l);
        else
            c2d = makePCurveOnSurf(c3, f, l, srf, kind);
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
        const double sampled = std::max(devT, devN);
        const double existing = BRep_Tool::Tolerance(e);
        const char* bindCls = "unhandled";
        double exact = exactMaxAtBind(c3, f, l, srf, c2d, loc, &bindCls);
        const bool isCyl = !Handle(Geom_CylindricalSurface)::DownCast(srf).IsNull();
        if (exact > cap && pinVerts && !plBind.IsNull()) {
            Handle(Geom2d_Curve) projPc = makePCurveOnSurf(c3, f, l, srf, kind);
            if (!projPc.IsNull()) {
                const char* clsP = bindCls;
                const double exactP = exactMaxAtBind(c3, f, l, srf, projPc, loc, &clsP);
                if (exactP >= 0.0 && exactP <= cap) {
                    c2d = projPc;
                    exact = exactP;
                    bindCls = clsP;
                }
            }
        }
        double tol = existing;
        bool overCap = false;
        double maxAll = exact;
        if (exact < 0.0) noteBindTolUnhandled(e, bindCls);
        const bool okTol = applyBindTolAllPcurves(e, exact, sampled, cap, &tol, &overCap, &maxAll);
        if (overCap) {
            gBindTolOverCap++;
            if (const void* ts = diagTShapePtr(e)) gBindTolOverCapTShapes.insert(ts);
            const double brep = brepDevMaxOnSurf(c3, f, l, srf, c2d, f, l, loc);
            const int genuine = (brep > cap) ? 1 : 0;
            if (genuine) gBindTolOverCapGenuine++;
            if (diagP2Enabled() || diagPlatesEnabled())
                std::fprintf(stderr,
                             "DIAG_OVERCAP ci=%d rid=%d class=%s exactMax=%.9f cap=%.9f "
                             "brepMax=%.9f existing=%.9f genuine=%d site=bindAllVariants\n",
                             ci, R.id, bindCls, maxAll, cap, brep, existing, genuine);
            stampTolWriter(e, isCyl ? "cylSide" : "internedPlane-overCap-skip");
            if (!okTol) tol = existing;
        } else if (exact < 0.0) {
            stampTolWriter(e, "legacy");
        } else if (isCyl) {
            stampTolWriter(e, "cylSide");
        } else {
            stampTolWriter(e, "bindTolFromExact");
        }
        if (diagP2Enabled() && (Handle(Geom_Plane)::DownCast(srf) || isCyl))
            std::fprintf(stderr,
                         "DIAG_BINDTOL ci=%d rid=%d class=%s exact=%.9f sampled=%.9f tol=%.9f "
                         "existing=%.9f cap=%.9f overCap=%d\n",
                         ci, R.id, bindCls, exact, sampled, tol, existing, cap, overCap ? 1 : 0);
        if (diagP2Enabled() && isCyl && exact >= 0.0) {
            const double brep = brepDevMaxOnSurf(c3, f, l, srf, c2d, f, l, loc);
            gTolCheckRows++;
            gTolCheckCylRows++;
            if (brep >= 0.0) {
                if (brep <= exact + 1e-12) {
                    gTolCheckBrepLeExact++;
                    gTolCheckCylBrepLe++;
                } else
                    gBindTolBrepGtExact++;
            }
            std::fprintf(stderr,
                         "DIAG_TOLCHECK_CYL ci=%d rid=%d class=%s exactMax=%.9f brepMax=%.9f "
                         "edgeTol=%.9f brepVerdict=%s\n",
                         ci, R.id, bindCls, exact, brep, tol,
                         (brep >= 0.0 && brep > exact + 1e-12) ? "brep>exact" : "brep<=exact");
        }
        if (isCyl && std::fabs(tol - existing) > Precision::PConfusion()) gCylSideTolChanged++;
        BRep_Builder B;
        B.UpdateEdge(e, c2d, srf, loc, tol);
        B.Range(e, srf, loc, f, l);
        acc.nWrite++;
    }
}

// Pcurve birth on the interned plane handle only. Does not enroll gSeamTShapes
// (mesh polylines are not cylinder seams).
void bindEdgePcurveOnInternedPlane(const TopoDS_Edge& e, const Region& plnR, const MeshView& mv,
                                   const char* kind, int ci) {
    if (plnR.type != SurfType::Plane || e.IsNull()) return;
    Handle(Geom_Surface) srf = regionSurf(plnR, SurfVar::Plane);
    if (srf.IsNull()) return;
    Standard_Real f = 0, l = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f, l);
    if (c3.IsNull() || l - f <= Precision::PConfusion()) return;
    const TopLoc_Location loc;
    Standard_Boolean hasPc = Standard_False;
    Standard_Real pf = 0, pl = 0;
    (void)BRep_Tool::CurveOnSurface(e, srf, loc, pf, pl, &hasPc);
    if (hasPc) return;
    Handle(Geom_Plane) plBind = Handle(Geom_Plane)::DownCast(srf);
    Handle(Geom2d_Curve) c2d;
    if (!plBind.IsNull())
        c2d = makePlanePCurvePinnedToVertices(e, plBind, c3, f, l);
    else
        c2d = makePCurveOnSurf(c3, f, l, srf, kind);
    if (c2d.IsNull()) return;
    double cap = meshTolCap(mv, &plnR);
    double devT = 0, devN = 0, dev = 0;
    const char* proj = "closedform";
    bool unproj = false;
    pcurveDevTN(c3, f, l, srf, c2d, devT, devN, &dev, &proj, &unproj);
    const double sampled = std::max(devT, devN);
    const double existing = BRep_Tool::Tolerance(e);
    const char* bindCls = "unhandled";
    double exact = exactMaxAtBind(c3, f, l, srf, c2d, loc, &bindCls);
    if (exact > cap && !plBind.IsNull()) {
        Handle(Geom2d_Curve) projPc = makePCurveOnSurf(c3, f, l, srf, kind);
        if (!projPc.IsNull()) {
            const char* clsP = bindCls;
            const double exactP = exactMaxAtBind(c3, f, l, srf, projPc, loc, &clsP);
            if (exactP >= 0.0 && exactP <= cap) {
                c2d = projPc;
                exact = exactP;
                bindCls = clsP;
            }
        }
    }
    double tol = existing;
    bool overCap = false;
    double maxAll = exact;
    if (exact < 0.0) noteBindTolUnhandled(e, bindCls);
    const bool okTol = applyBindTolAllPcurves(e, exact, sampled, cap, &tol, &overCap, &maxAll);
    if (overCap) {
        gBindTolOverCap++;
        if (const void* ts = diagTShapePtr(e)) gBindTolOverCapTShapes.insert(ts);
        const double brep = brepDevMaxOnSurf(c3, f, l, srf, c2d, f, l, loc);
        const int genuine = (brep > cap) ? 1 : 0;
        if (genuine) gBindTolOverCapGenuine++;
        if (diagP2Enabled() || diagPlatesEnabled())
            std::fprintf(stderr,
                         "DIAG_OVERCAP ci=%d rid=%d class=%s exactMax=%.9f cap=%.9f "
                         "brepMax=%.9f existing=%.9f genuine=%d site=internedPlane\n",
                         ci, plnR.id, bindCls, maxAll, cap, brep, existing, genuine);
        stampTolWriter(e, "internedPlane-overCap-skip");
        if (!okTol) tol = existing;
    } else if (exact < 0.0) {
        stampTolWriter(e, "legacy");
    } else {
        stampTolWriter(e, "bindTolFromExact");
    }
    if (diagP2Enabled())
        std::fprintf(stderr,
                     "DIAG_BINDTOL ci=%d rid=%d class=%s exact=%.9f sampled=%.9f tol=%.9f "
                     "existing=%.9f cap=%.9f overCap=%d\n",
                     ci, plnR.id, bindCls, exact, sampled, tol, existing, cap, overCap ? 1 : 0);
    BRep_Builder B;
    B.UpdateEdge(e, c2d, srf, loc, tol);
    B.Range(e, srf, loc, f, l);
    (void)ci;
}

void dumpPcurveFrames(const Region& r, const TopoDS_Wire& w) {
    if ((!diagP2Enabled() && !diagPlatesEnabled()) || w.IsNull()) return;
    if (r.id < 0 || r.id > 3) return;
    Handle(Geom_Surface) interned = regionSurf(r, SurfVar::Plane);
    const void* internH = interned.IsNull() ? nullptr : interned.get();
    int nE = 0, nMismatch = 0, nStoredInterned = 0, nNoStored = 0, nOtherHandle = 0;
    int ei = 0;
    const TopLoc_Location loc;
    gp_Pnt iOrig;
    gp_Dir iX(1, 0, 0);
    Handle(Geom_Plane) gpl = Handle(Geom_Plane)::DownCast(interned);
    if (!gpl.IsNull()) {
        iOrig = gpl->Location();
        iX = gpl->Pln().XAxis().Direction();
    }
    std::fprintf(stderr,
                 "DIAG_PCFRAME_INTERNED rid=%d interned=%p iOrig=(%.6f,%.6f,%.6f) "
                 "iX=(%.6f,%.6f,%.6f)\n",
                 r.id, internH, iOrig.X(), iOrig.Y(), iOrig.Z(), iX.X(), iX.Y(), iX.Z());
    for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next(), ei++) {
        TopoDS_Edge e = ex.Current();
        nE++;
        Standard_Boolean hasPc = Standard_False;
        Standard_Real pf = 0, pl = 0;
        if (interned) (void)BRep_Tool::CurveOnSurface(e, interned, loc, pf, pl, &hasPc);
        if (!hasPc) nNoStored++;
        else nStoredInterned++;
        int edgeMismatch = 0;
        try {
            Handle(BRep_TEdge) te = Handle(BRep_TEdge)::DownCast(e.TShape());
            if (!te.IsNull()) {
                for (BRep_ListOfCurveRepresentation::Iterator it(te->Curves()); it.More();
                     it.Next()) {
                    Handle(BRep_CurveOnSurface) cos =
                        Handle(BRep_CurveOnSurface)::DownCast(it.Value());
                    if (cos.IsNull()) continue;
                    Handle(Geom_Surface) ps = cos->Surface();
                    const void* ph = ps.IsNull() ? nullptr : static_cast<const void*>(ps.get());
                    if (ph && internH && ph != internH) {
                        edgeMismatch = 1;
                        nOtherHandle++;
                        gp_Pnt sOrig;
                        gp_Dir sX(1, 0, 0);
                        Handle(Geom_Plane) spl = Handle(Geom_Plane)::DownCast(ps);
                        if (!spl.IsNull()) {
                            sOrig = spl->Location();
                            sX = spl->Pln().XAxis().Direction();
                        }
                        std::fprintf(stderr,
                                     "DIAG_PCFRAME rid=%d i=%d interned=%p stored=%p mismatch=1 "
                                     "iOrig=(%.6f,%.6f,%.6f) iX=(%.6f,%.6f,%.6f) "
                                     "sOrig=(%.6f,%.6f,%.6f) sX=(%.6f,%.6f,%.6f)\n",
                                     r.id, ei, internH, ph, iOrig.X(), iOrig.Y(), iOrig.Z(),
                                     iX.X(), iX.Y(), iX.Z(), sOrig.X(), sOrig.Y(), sOrig.Z(),
                                     sX.X(), sX.Y(), sX.Z());
                    }
                }
            }
        } catch (const Standard_Failure&) {
        }
        if (edgeMismatch) nMismatch++;
    }
    std::fprintf(stderr,
                 "DIAG_PCFRAME_SUM rid=%d mismatches=%d/%d storedInterned=%d noStored=%d "
                 "otherHandle=%d examined=%d present=%d\n",
                 r.id, nMismatch, nE, nStoredInterned, nNoStored, nOtherHandle, nE, nE);
    // Oriented UV closure on the interned handle (D-S3-69 functional metric).
    try {
        int nE_iter = 0;
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() == TopAbs_EDGE) nE_iter++;
        }
        std::vector<gp_Pnt2d> uv0s, uv1s;
        std::vector<int> stored;
        int wi = 0;
        for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next(), wi++) {
            TopoDS_Edge e = ex.Current();
            gp_Pnt2d a, b;
            int st = 0;
            if (!interned.IsNull()) {
                Standard_Boolean hasPc = Standard_False;
                Standard_Real pf = 0, pl = 0;
                Handle(Geom2d_Curve) c2 = BRep_Tool::CurveOnSurface(e, interned, loc, pf, pl, &hasPc);
                if (hasPc && !c2.IsNull()) {
                    st = 1;
                    if (e.Orientation() == TopAbs_REVERSED) std::swap(pf, pl);
                    a = c2->Value(pf);
                    b = c2->Value(pl);
                }
            }
            uv0s.push_back(a);
            uv1s.push_back(b);
            stored.push_back(st);
        }
        const int nW = (int)uv0s.size();
        double maxG = 0;
        int gi = -1, gj = -1;
        int nWalkPairs = 0, nWalkStored = 0;
        for (int i = 0; i < nW; i++) {
            if (!stored[(size_t)i] || !stored[(size_t)((i + 1) % nW)]) continue;
            nWalkPairs++;
            nWalkStored++;
            const int j = (i + 1) % nW;
            const double g = uv1s[(size_t)i].Distance(uv0s[(size_t)j]);
            if (g > maxG) {
                maxG = g;
                gi = i;
                gj = j;
            }
        }
        // Iterator-order UV closure (full wire coverage when WireExplorer sees one closed edge).
        double maxIterG = 0;
        int iterGi = -1, iterGj = -1;
        int nIterPairs = 0, nIterStored = 0;
        std::vector<gp_Pnt2d> iuv0, iuv1;
        std::vector<int> istored;
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            const TopoDS_Edge e = TopoDS::Edge(it.Value());
            gp_Pnt2d a, b;
            int st = 0;
            if (!interned.IsNull()) {
                Standard_Boolean hasPc = Standard_False;
                Standard_Real pf = 0, pl = 0;
                Handle(Geom2d_Curve) c2 = BRep_Tool::CurveOnSurface(e, interned, loc, pf, pl, &hasPc);
                if (hasPc && !c2.IsNull()) {
                    st = 1;
                    if (e.Orientation() == TopAbs_REVERSED) std::swap(pf, pl);
                    a = c2->Value(pf);
                    b = c2->Value(pl);
                }
            }
            iuv0.push_back(a);
            iuv1.push_back(b);
            istored.push_back(st);
        }
        const int nIterE = (int)iuv0.size();
        for (int i = 0; i + 1 < nIterE; i++) {
            if (!istored[(size_t)i] || !istored[(size_t)i + 1]) continue;
            nIterPairs++;
            nIterStored++;
            const double g = iuv1[(size_t)i].Distance(iuv0[(size_t)i + 1]);
            if (g > maxIterG) {
                maxIterG = g;
                iterGi = i;
                iterGj = i + 1;
            }
        }
        if (nIterE > 1 && istored[(size_t)nIterE - 1] && istored[0]) {
            nIterPairs++;
            nIterStored++;
            const double g = iuv1[(size_t)nIterE - 1].Distance(iuv0[0]);
            if (g > maxIterG) {
                maxIterG = g;
                iterGi = nIterE - 1;
                iterGj = 0;
            }
        }
        const int closedEdgeWalk = (nW == 1 && nE_iter > 1) ? 1 : 0;
        if (gi >= 0) {
            std::fprintf(stderr,
                         "DIAG_UVWALK rid=%d nE_walk=%d nE_iter=%d closedEdgeWalk=%d "
                         "maxGapUv=%.8f i=%d j=%d walkStored=%d/%d "
                         "uv1i=(%.6f,%.6f) uv0j=(%.6f,%.6f)\n",
                         r.id, nW, nE_iter, closedEdgeWalk, maxG, gi, gj, nWalkStored,
                         std::max(nWalkPairs, 1), uv1s[(size_t)gi].X(), uv1s[(size_t)gi].Y(),
                         uv0s[(size_t)gj].X(), uv0s[(size_t)gj].Y());
        } else {
            std::fprintf(stderr,
                         "DIAG_UVWALK rid=%d nE_walk=%d nE_iter=%d closedEdgeWalk=%d maxGapUv=-1 "
                         "walkStored=%d/%d\n",
                         r.id, nW, nE_iter, closedEdgeWalk, nWalkStored, std::max(nWalkPairs, 1));
        }
        std::fprintf(stderr,
                     "DIAG_UVITER rid=%d nE_iter=%d maxGapUv=%.8f i=%d j=%d iterStored=%d/%d\n",
                     r.id, nIterE, maxIterG, iterGi, iterGj, nIterStored,
                     std::max(nIterPairs, 1));
    } catch (const Standard_Failure&) {
    }
}

void birthMeshEdgePlanePCurves(const RegionSet& rs, const std::vector<TopoDS_Edge>& meshE,
                               const std::vector<char>& edgeOk, const MeshView& mv) {
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        const BoundaryChain& ch = rs.chains[ci];
        const Region* A = regionById(rs, ch.regA);
        const Region* B = regionById(rs, ch.regB);
        const Region* planes[2] = {nullptr, nullptr};
        int np = 0;
        if (A && A->type == SurfType::Plane && isAnalytic(A)) planes[np++] = A;
        if (B && B->type == SurfType::Plane && isAnalytic(B)) planes[np++] = B;
        if (np == 0) continue;
        for (int eid : ch.meshEdges) {
            if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
            if ((size_t)eid < edgeOk.size() && !edgeOk[(size_t)eid]) continue;
            const void* ts = diagTShapePtr(meshE[(size_t)eid]);
            if (ts && gSeamTShapes.count(ts)) continue;
            for (int k = 0; k < np; k++)
                bindEdgePcurveOnInternedPlane(meshE[(size_t)eid], *planes[k], mv, "poly", (int)ci);
        }
    }
}

gp_Pnt probeArcMid(const AnalyticCurve& curve, const BoundaryChain& ch, const MeshView& mv,
                   const Region* A, const Region* B, const std::vector<TopoDS_Vertex>& verts,
                   int ia, int ib) {
    gp_Pnt mid = BRep_Tool::Pnt(verts[(size_t)ia]);
    const Region* cylHint = (A && A->type == SurfType::Cylinder) ? A
                          : (B && B->type == SurfType::Cylinder) ? B
                                                                : nullptr;
    const Region* plnHint = (A && A->type == SurfType::Plane) ? A
                          : (B && B->type == SurfType::Plane) ? B
                                                            : nullptr;
    const bool full = ch.closedLoop;
    const bool partialCylArc =
        cylHint && !cylHint->closed360 && curve.kind == AnalyticCurve::Circ;
    if (partialCylArc && ch.meshVerts.size() < 3)
        return partialCylCapArcMid(*cylHint, plnHint, curve.circ);
    if ((!full && curve.kind == AnalyticCurve::Circ && ch.meshVerts.size() >= 3) ||
        (partialCylArc && ch.meshVerts.size() >= 3)) {
        mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
        // D-130-14: the (uMin+uMax)/2 probe is a PARTIAL cylinder's midpoint.
        // On a closed360 region uMin/uMax are 0 and 2*pi by definition, so it
        // reads pi -- an azimuth an arc chain of that region need not contain,
        // and the "closest vertex" it then returns is an ENDPOINT, which tells
        // makeArc nothing (dm = 0 always takes the short arc). A 288 deg rim
        // chain split by pinch vertices came back as its own 72 deg complement.
        // The chain's own middle vertex is the honest hint there.
        if (cylHint && !(unionBuildOn() && cylHint->closed360)) {
            const double um = 0.5 * (cylHint->uMin + cylHint->uMax);
            double best = 1e300;
            for (int lv : ch.meshVerts) {
                double du = std::fabs(regionU(*cylHint, pntOf(mv, lv)) - um);
                if (du < best) {
                    best = du;
                    mid = pntOf(mv, lv);
                }
            }
        }
        return mid;
    }
    if (!full && curve.kind == AnalyticCurve::Elips && ch.meshVerts.size() >= 3)
        return pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
    if (!full && curve.kind == AnalyticCurve::Elips && ch.meshVerts.size() >= 2) {
        gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)ia]);
        gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)ib]);
        return gp_Pnt(0.5 * (pa.X() + pb.X()), 0.5 * (pa.Y() + pb.Y()),
                      0.5 * (pa.Z() + pb.Z()));
    }
    return mid;
}

TopoDS_Edge probeCollapsedMakeEdge(const ChainGeom& g, const BoundaryChain& ch, const MeshView& mv,
                                   const RegionSet& rs, const std::vector<TopoDS_Vertex>& verts) {
    const AnalyticCurve& curve = g.curve;
    const int ia = g.ia, ib = g.ib;
    if (curve.kind == AnalyticCurve::None || ia < 0 || ib < 0 || (size_t)ia >= verts.size() ||
        (size_t)ib >= verts.size())
        return TopoDS_Edge();
    const Region* A = regionById(rs, g.regA);
    const Region* B = regionById(rs, g.regB);
    const bool full = ch.closedLoop;
    const gp_Pnt mid = probeArcMid(curve, ch, mv, A, B, verts, ia, ib);
    TopoDS_Edge e;
    if (full && curve.kind == AnalyticCurve::Circ) {
        e = makeFullCircle(curve.circ, verts[(size_t)ia]);
        if (e.IsNull()) e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
    } else if (!full && curve.kind == AnalyticCurve::Circ && ch.meshVerts.size() >= 3) {
        e = makeArc(curve.circ, verts[(size_t)ia], verts[(size_t)ib], mid);
        if (e.IsNull()) e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
    } else if (!full && curve.kind == AnalyticCurve::Elips && ch.meshVerts.size() >= 2) {
        e = makeEllipseArc(curve.elips, verts[(size_t)ia], verts[(size_t)ib], mid);
        if (e.IsNull()) e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
    } else {
        e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
    }
    return e;
}

thread_local int gSeamFailFirst = 0;

const char* seamPairKind(AnalyticCurve::Kind a, AnalyticCurve::Kind b) {
    const bool la = a == AnalyticCurve::Lin, lb = b == AnalyticCurve::Lin;
    const bool ca = a == AnalyticCurve::Circ, cb = b == AnalyticCurve::Circ;
    if (la && lb) return "lin-lin";
    if ((la && cb) || (ca && lb)) return "lin-circ";
    if (ca && cb) return "circ-circ";
    if (a == AnalyticCurve::None && la) return "single-lin";
    if (a == AnalyticCurve::None && ca) return "single-circ";
    return "other";
}

int snapSeamVertices(const RegionSet& rs, const MeshView& mv, const std::vector<ChainGeom>& geom,
                     const std::vector<TopoDS_Vertex>& verts, std::unordered_set<int>& seamLocked) {
    int nExamined = 0, nMoved = 0, nNamed = 0, nSkipDeg = 0;
    int examLL = 0, examLC = 0, examCC = 0, examS = 0;
    int movedLL = 0, movedLC = 0, movedCC = 0, movedS = 0;
    int namedLL = 0, namedLC = 0, namedCC = 0, namedS = 0;
    auto bumpKind = [&](const char* pk, int exam, int moved, int named) {
        if (std::strcmp(pk, "lin-lin") == 0) {
            examLL += exam;
            movedLL += moved;
            namedLL += named;
        } else if (std::strcmp(pk, "lin-circ") == 0) {
            examLC += exam;
            movedLC += moved;
            namedLC += named;
        } else if (std::strcmp(pk, "circ-circ") == 0) {
            examCC += exam;
            movedCC += moved;
            namedCC += named;
        } else {
            examS += exam;
            movedS += moved;
            namedS += named;
        }
    };
    std::unordered_map<int, std::vector<int>> atVtx;
    for (size_t ci = 0; ci < geom.size(); ci++) {
        if (geom[ci].curve.kind != AnalyticCurve::Lin && geom[ci].curve.kind != AnalyticCurve::Circ)
            continue;
        if (geom[ci].ia >= 0) atVtx[geom[ci].ia].push_back((int)ci);
        if (geom[ci].ib >= 0 && geom[ci].ib != geom[ci].ia)
            atVtx[geom[ci].ib].push_back((int)ci);
    }
    for (auto& kv : atVtx) {
        const int vid = kv.first;
        std::vector<int> allCh = kv.second;
        std::sort(allCh.begin(), allCh.end());
        allCh.erase(std::unique(allCh.begin(), allCh.end()), allCh.end());
        if (vid < 0 || (size_t)vid >= verts.size() || verts[(size_t)vid].IsNull()) continue;
        if (allCh.empty()) continue;
        struct Pair {
            int c0, c1;
            const char* pk;
        };
        std::vector<Pair> pairs;
        if (allCh.size() == 1) {
            const AnalyticCurve::Kind k = geom[(size_t)allCh[0]].curve.kind;
            pairs.push_back({allCh[0], -1, k == AnalyticCurve::Lin ? "single-lin" : "single-circ"});
        } else {
            for (size_t i = 0; i < allCh.size(); i++) {
                for (size_t j = i + 1; j < allCh.size(); j++) {
                    const char* pk = seamPairKind(geom[(size_t)allCh[i]].curve.kind,
                                                  geom[(size_t)allCh[j]].curve.kind);
                    pairs.push_back({allCh[i], allCh[j], pk});
                }
            }
        }
        const gp_Pnt meshP = BRep_Tool::Pnt(verts[(size_t)vid]);
        std::vector<const Region*> incident;
        auto addReg = [&](int rid) {
            const Region* R = regionById(rs, rid);
            if (!R || !isAnalytic(R)) return;
            for (const Region* x : incident)
                if (x->id == rid) return;
            incident.push_back(R);
        };
        for (int ci : allCh) {
            addReg(geom[(size_t)ci].regA);
            addReg(geom[(size_t)ci].regB);
        }
        for (const BoundaryChain& chV : rs.chains) {
            bool hit = false;
            for (int mvtx : chV.meshVerts) {
                if (mvtx == vid) {
                    hit = true;
                    break;
                }
            }
            if (!hit) continue;
            addReg(chV.regA);
            addReg(chV.regB);
        }
        if (vid >= 0 && (size_t)vid < mv.nVtx && mv.compVtx && mv.tris && mv.compTris) {
            const int gvid = mv.compVtx[vid];
            for (const Region& r : rs.regions) {
                if (!isAnalytic(&r)) continue;
                bool hit = false;
                for (int lt : r.tris) {
                    if (lt < 0 || (size_t)lt >= mv.nTri) continue;
                    const int* T = mv.tris[mv.compTris[lt]];
                    if (T[0] == gvid || T[1] == gvid || T[2] == gvid) {
                        hit = true;
                        break;
                    }
                }
                if (hit) addReg(r.id);
            }
        }
        std::vector<double> ds;
        std::vector<int> boundRids;
        bool devOk = true;
        int failRid = -1;
        for (const Region* R : incident) {
            const double d = regionPointDev(R, meshP);
            const double bound = regionDevBound(mv, R);
            ds.push_back(d);
            boundRids.push_back(R->id);
            if (d > bound && failRid < 0) failRid = R->id;
            if (d > bound) devOk = false;
        }
        int hubRid = -1;
        for (const Region* R : incident)
            if (R->type == SurfType::Plane && R->id >= 0 && R->id <= 3) hubRid = R->id;
        double theta = -1.0;
        if (incident.size() >= 2)
            theta = dihedralBetweenRegions(incident[0], incident[1], meshP);
        for (const Pair& pr : pairs) {
            nExamined++;
            const char* seamKind = (pr.c1 < 0) ? "project" : "junction";
            const char* pk = pr.pk;
            gp_Pnt target;
            bool targetOk = false;
            const ChainGeom& g0 = geom[(size_t)pr.c0];
            const bool closed0 =
                (pr.c0 >= 0 && (size_t)pr.c0 < rs.chains.size() && rs.chains[(size_t)pr.c0].closedLoop) ||
                (g0.ia == g0.ib);
            const gp_Pnt pa0 = BRep_Tool::Pnt(verts[(size_t)g0.ia]);
            const gp_Pnt pb0 = (g0.ib >= 0 && (size_t)g0.ib < verts.size())
                                   ? BRep_Tool::Pnt(verts[(size_t)g0.ib])
                                   : pa0;
            if (pr.c1 < 0) {
                targetOk = projectOntoCurveSpan(g0.curve, meshP, pa0, closed0 ? pa0 : pb0, target);
                if (!targetOk) seamKind = "project-range-fail";
            } else {
                const ChainGeom& g1 = geom[(size_t)pr.c1];
                const bool mixed = (g0.curve.kind == AnalyticCurve::Lin &&
                                    g1.curve.kind == AnalyticCurve::Circ) ||
                                   (g0.curve.kind == AnalyticCurve::Circ &&
                                    g1.curve.kind == AnalyticCurve::Lin);
                seamKind = mixed ? "junction" : "junction-same";
                gp_Pnt p1, p2;
                const int nInt = intersectAnalyticCurves(g0.curve, g1.curve, p1, p2);
                if (nInt <= 0) {
                    seamKind = "no-intersect";
                } else {
                    target = (nInt == 1 || p1.Distance(meshP) <= p2.Distance(meshP)) ? p1 : p2;
                    const gp_Pnt pa1 = BRep_Tool::Pnt(verts[(size_t)g1.ia]);
                    const gp_Pnt pb1 = BRep_Tool::Pnt(verts[(size_t)g1.ib]);
                    const bool closed1 =
                        (pr.c1 >= 0 && (size_t)pr.c1 < rs.chains.size() &&
                         rs.chains[(size_t)pr.c1].closedLoop) ||
                        (g1.ia == g1.ib);
                    const bool on0 = closed0 || pointOnCurveSpan(g0.curve, target, pa0, pb0);
                    const bool on1 = closed1 || pointOnCurveSpan(g1.curve, target, pa1, pb1);
                    targetOk = on0 && on1;
                    if (!targetOk) seamKind = "out-of-range";
                    if (nInt == 2 &&
                        std::fabs(p1.Distance(meshP) - p2.Distance(meshP)) <= Precision::Confusion())
                        targetOk = false;
                }
            }
            const double dist = targetOk ? meshP.Distance(target) : -1.0;
            auto emit = [&](const char* status) {
                if (!(diagP2Enabled() || diagPlatesEnabled())) return;
                std::fprintf(stderr,
                             "DIAG_SEAMVTX rid=%d vid=%d kind=%s d1=%.8f d2=%.8f d3=%.8f "
                             "theta=%.6f dist=%.8f boundRids=%d,%d,%d status=%s ci0=%d ci1=%d "
                             "pair=%s nInc=%zu failRid=%d\n",
                             hubRid, vid, seamKind, ds.size() > 0 ? ds[0] : -1.0,
                             ds.size() > 1 ? ds[1] : -1.0, ds.size() > 2 ? ds[2] : -1.0, theta, dist,
                             boundRids.size() > 0 ? boundRids[0] : -1,
                             boundRids.size() > 1 ? boundRids[1] : -1,
                             boundRids.size() > 2 ? boundRids[2] : -1, status, pr.c0, pr.c1, pk,
                             incident.size(), failRid);
            };
            if (!devOk) {
                nNamed++;
                bumpKind(pk, 1, 0, 1);
                char st[64];
                std::snprintf(st, sizeof(st), "named-dev-%d", failRid);
                emit(st);
                continue;
            }
            if (!targetOk) {
                nNamed++;
                bumpKind(pk, 1, 0, 1);
                emit(seamKind);
                continue;
            }
            // Single-curve: named, not snapped. Enabling MakeVertex on these
            // (10–11 single-circ) sent recover=3 (ok 206→39) and reverted.
            if (pr.c1 < 0) {
                nNamed++;
                bumpKind(pk, 1, 0, 1);
                emit("named-single");
                continue;
            }
            if (dist <= Precision::Confusion() || seamLocked.count(vid)) {
                bumpKind(pk, 1, 0, 0);
                emit("already");
                continue;
            }
            const TopoDS_Vertex saved = verts[(size_t)vid];
            const double t0 = BRep_Tool::Tolerance(saved);
            TopoDS_Vertex nv;
            BRep_Builder Bb;
            Bb.MakeVertex(nv, target, t0 > 0.0 ? t0 : Precision::Confusion());
            const_cast<TopoDS_Vertex&>(verts[(size_t)vid]) = nv;
            bool okProbe = true;
            int probeFailCi = -1;
            for (int ci : allCh) {
                if (ci < 0 || (size_t)ci >= rs.chains.size()) {
                    okProbe = false;
                    probeFailCi = ci;
                    break;
                }
                if (probeCollapsedMakeEdge(geom[(size_t)ci], rs.chains[(size_t)ci], mv, rs, verts)
                        .IsNull()) {
                    okProbe = false;
                    probeFailCi = ci;
                    break;
                }
            }
            if (!okProbe) {
                const_cast<TopoDS_Vertex&>(verts[(size_t)vid]) = saved;
                nNamed++;
                bumpKind(pk, 1, 0, 1);
                char st[64];
                std::snprintf(st, sizeof(st), "makeedge-fail-ci%d", probeFailCi);
                emit(st);
                continue;
            }
            nMoved++;
            seamLocked.insert(vid);
            bumpKind(pk, 1, 1, 0);
            emit("moved");
        }
        (void)nSkipDeg;
    }
    if (diagP2Enabled() || diagPlatesEnabled())
        std::fprintf(stderr,
                     "DIAG_SEAMVTX_SUM examined=%d moved=%d named=%d/%d "
                     "lin-lin=%d/%d/%d lin-circ=%d/%d/%d circ-circ=%d/%d/%d "
                     "single=%d/%d/%d skippedDeg=%d\n",
                     nExamined, nMoved, nNamed, nExamined, examLL, movedLL, namedLL, examLC,
                     movedLC, namedLC, examCC, movedCC, namedCC, examS, movedS, namedS, nSkipDeg);
    return nMoved;
}

void bindPlateMeshPcurves(const RegionSet& rs, const std::vector<char>& collapsed,
                          const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                          const MeshView& mv) {
    for (size_t ci = 0; ci < rs.chains.size(); ci++) {
        if ((size_t)ci < collapsed.size() && collapsed[ci]) continue;
        const BoundaryChain& ch = rs.chains[ci];
        const Region* A = regionById(rs, ch.regA);
        const Region* B = regionById(rs, ch.regB);
        const Region* planes[2] = {nullptr, nullptr};
        int np = 0;
        if (A && A->type == SurfType::Plane && isAnalytic(A)) planes[np++] = A;
        if (B && B->type == SurfType::Plane && isAnalytic(B)) planes[np++] = B;
        if (np == 0) continue;
        if ((A && A->type == SurfType::Cylinder && isAnalytic(A)) ||
            (B && B->type == SurfType::Cylinder && isAnalytic(B)))
            continue;
        for (int eid : ch.meshEdges) {
            if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
            if ((size_t)eid < edgeOk.size() && !edgeOk[(size_t)eid]) continue;
            const void* ts = diagTShapePtr(meshE[(size_t)eid]);
            if (ts && gSeamTShapes.count(ts)) continue;
            for (int k = 0; k < np; k++)
                bindEdgePcurveOnInternedPlane(meshE[(size_t)eid], *planes[k], mv, "poly", (int)ci);
        }
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

    if (!makeFaceKeep(surf, outer, inners, outward, outF)) return false;
    addPcurvesOnFace(outF, sewTol, false);
    return !outF.IsNull();
}

bool hubWireEdgeUv(const TopoDS_Edge& e, const Handle(Geom_Surface)& surf, const TopLoc_Location& loc,
                   gp_Pnt2d& uv0, gp_Pnt2d& uv1) {
    Standard_Boolean hasPc = Standard_False;
    Standard_Real pf = 0, pl = 0;
    Handle(Geom2d_Curve) c2 = BRep_Tool::CurveOnSurface(e, surf, loc, pf, pl, &hasPc);
    if (!hasPc || c2.IsNull()) return false;
    if (e.Orientation() == TopAbs_REVERSED) std::swap(pf, pl);
    uv0 = c2->Value(pf);
    uv1 = c2->Value(pl);
    return true;
}

void dumpDiagHubWire(const Region& r, const TopoDS_Wire& ow, const Handle(Geom_Plane)& gpl,
                     const RegionSet& rs, const MeshView& mv, const std::vector<ChainGeom>& geom,
                     const std::vector<char>& collapsed, const std::vector<TopoDS_Edge>& meshE,
                     const std::vector<TopoDS_Wire>& innerW, bool outward) {
    if (!diagP2Enabled() || ow.IsNull() || gpl.IsNull()) return;
    const double sewTol = mv.sewTol;
    const double tol = std::max(sewTol, Precision::PConfusion());
    const TopLoc_Location loc;
    const Handle(Geom_Surface) surf(gpl);
    const int rid = r.id;

    auto matchChain = [&](const TopoDS_Edge& e, const gp_Pnt& p0, const gp_Pnt& p1) -> int {
        const void* ts = diagTShapePtr(e);
        for (size_t ci = 0; ci < geom.size(); ci++) {
            for (const TopoDS_Edge& ge : geom[ci].edges) {
                if (!ge.IsNull() && (e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts)))
                    return (int)ci;
            }
        }
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            for (int eid : rs.chains[ci].meshEdges) {
                if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
                if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts))
                    return (int)ci;
            }
        }
        int best = -1;
        double bestD = 1e300;
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            const BoundaryChain& ch = rs.chains[ci];
            if (ch.meshVerts.size() < 2) continue;
            gp_Pnt a = pntOf(mv, ch.meshVerts.front());
            gp_Pnt b = pntOf(mv, ch.meshVerts.back());
            const double d = std::min(p0.Distance(a) + p1.Distance(b), p0.Distance(b) + p1.Distance(a));
            if (d < bestD) {
                bestD = d;
                best = (int)ci;
            }
        }
        return (best >= 0 && bestD <= 2.0 * tol) ? best : -1;
    };

    auto edgeInfo = [&](const TopoDS_Edge& e) {
        struct Info {
            const char* cls;
            int ci;
            int partner;
            int fb;
        } info{"fill", -1, -1, 0};
        const void* ts = diagTShapePtr(e);
        const bool isSeamTs = ts && gSeamTShapes.count(ts);
        TopoDS_Vertex v0, v1;
        TopExp::Vertices(e, v0, v1, Standard_True);
        gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
        gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
        if (isSeamTs) {
            // D-130: an analytic edge still belongs to a chain. Returning before
            // matchChain reported ci=-1 for every seam edge and made the tripwire
            // ("analytic edge with no chain index") unmeasurable.
            info.cls = "seam";
            info.ci = matchChain(e, p0, p1);
            if (info.ci >= 0 && (size_t)info.ci < rs.chains.size()) {
                const BoundaryChain& ch = rs.chains[(size_t)info.ci];
                info.partner = (ch.regA == rid) ? ch.regB : ch.regA;
                if (info.partner < 0) info.partner = ch.regA >= 0 ? ch.regA : ch.regB;
                info.fb = gChainSewFallbackCi.count(info.ci) ? 1 : 0;
            }
            return info;
        }
        info.ci = matchChain(e, p0, p1);
        bool isMeshE = false;
        if (info.ci >= 0 && (size_t)info.ci < rs.chains.size()) {
            for (int eid : rs.chains[(size_t)info.ci].meshEdges) {
                if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
                if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts)) {
                    isMeshE = true;
                    break;
                }
            }
            const BoundaryChain& ch = rs.chains[(size_t)info.ci];
            info.partner = (ch.regA == rid) ? ch.regB : ch.regA;
            if (info.partner < 0) info.partner = ch.regA >= 0 ? ch.regA : ch.regB;
            info.fb = gChainSewFallbackCi.count(info.ci) ? 1 : 0;
            if (info.ci >= 0 && (size_t)info.ci < collapsed.size() && collapsed[(size_t)info.ci] &&
                (size_t)info.ci < geom.size() && geom[(size_t)info.ci].collapsed &&
                !geom[(size_t)info.ci].edges.empty()) {
                const Region* rr = regionById(rs, info.partner);
                if (rr && rr->type == SurfType::Cylinder) info.cls = "cyl";
                else if (rr && rr->type == SurfType::Plane) info.cls = "pln";
                else info.cls = "span";
            } else if (isMeshE) {
                info.cls = "meshE";
            } else {
                const Region* rr = regionById(rs, info.partner);
                if (rr && rr->type == SurfType::Cylinder) info.cls = "cyl";
                else if (rr && rr->type == SurfType::Plane) info.cls = "pln";
                else info.cls = "fill";
            }
        }
        return info;
    };

    struct HubE {
        int idx;
        const void* eT;
        const void* v0T;
        const void* v1T;
        const char* cls;
        int ci;
        int partner;
        int fb;
        int storedPc;
        gp_Pnt p0, p1;
        gp_Pnt2d uv0, uv1;
        int haveUv;
        TopoDS_Edge e;
    };

    int r7Hubwire = 0, r7HubwirePresent = 0;
    int r7Share = 0, r7SharePresent = 0;
    int r7Comp = 0, r7CompPresent = 0;
    int r7Pair = 0, r7PairPresent = 0;
    int r7Inner = 0, r7InnerPresent = 0;
    int r7Orient = 0, r7OrientPresent = 0;

    auto collectEdges = [&](const TopoDS_Wire& w, bool explorer) -> std::vector<TopoDS_Edge> {
        std::vector<TopoDS_Edge> out;
        if (explorer) {
            for (BRepTools_WireExplorer ex(w); ex.More(); ex.Next()) out.push_back(ex.Current());
        } else {
            for (TopoDS_Iterator it(w); it.More(); it.Next()) {
                if (it.Value().ShapeType() == TopAbs_EDGE) out.push_back(TopoDS::Edge(it.Value()));
            }
        }
        return out;
    };

    auto dumpWireRow = [&](const TopoDS_Wire& w, int iw, bool explorer) {
        gHubWireRowB = explorer;
        const char* row = hubWireRowLabel();
        const std::vector<TopoDS_Edge> elist = collectEdges(w, explorer);
        std::vector<HubE> edges;
        edges.reserve(elist.size());
        int i = 0;
        for (const TopoDS_Edge& e : elist) {
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(e, v0, v1, Standard_True);
            auto inf = edgeInfo(e);
            HubE he;
            he.idx = i++;
            he.e = e;
            he.eT = diagTShapePtr(e);
            he.v0T = diagTShapePtr(v0);
            he.v1T = diagTShapePtr(v1);
            he.p0 = BRep_Tool::Pnt(v0);
            he.p1 = BRep_Tool::Pnt(v1);
            he.cls = inf.cls;
            he.ci = inf.ci;
            he.partner = inf.partner;
            he.fb = inf.fb;
            he.storedPc = 0;
            he.haveUv = hubWireEdgeUv(e, surf, loc, he.uv0, he.uv1) ? 1 : 0;
            if (he.haveUv) he.storedPc = 1;
            r7Hubwire++;
            r7HubwirePresent++;
            std::fprintf(stderr,
                         "DIAG_HUBWIRE rid=%d row=%s iw=%d i=%d eT=%p v0T=%p v1T=%p class=%s "
                         "ci=%d partner=%d fb=%d storedPc=%d "
                         "p0=(%.4f,%.4f,%.4f) p1=(%.4f,%.4f,%.4f) "
                         "uv0=(%.4f,%.4f) uv1=(%.4f,%.4f) haveUv=%d\n",
                         rid, row, iw, he.idx, he.eT, he.v0T, he.v1T, he.cls, he.ci, he.partner,
                         he.fb, he.storedPc, he.p0.X(), he.p0.Y(), he.p0.Z(), he.p1.X(), he.p1.Y(),
                         he.p1.Z(), he.uv0.X(), he.uv0.Y(), he.uv1.X(), he.uv1.Y(), he.haveUv);
            edges.push_back(he);
        }

        int nGap3dOverTol = 0;
        double maxGap3d = 0.0;
        int nSharedVtx = 0;
        auto gapPair = [&](size_t j, size_t k) {
            if (j >= edges.size() || k >= edges.size()) return;
            const double g3 = edges[j].p1.Distance(edges[k].p0);
            if (g3 > maxGap3d) maxGap3d = g3;
            if (g3 > tol) nGap3dOverTol++;
            const int sameVtx = (edges[j].v1T && edges[j].v1T == edges[k].v0T) ? 1 : 0;
            if (sameVtx) nSharedVtx++;
            double gUv = -1.0;
            if (edges[j].haveUv && edges[k].haveUv) gUv = edges[j].uv1.Distance(edges[k].uv0);
            const bool report = !sameVtx || gUv > tol || g3 > tol ||
                                (rid == 2 && iw == 0 &&
                                 ((int)j == 39 || (int)j == 71 || (int)k == 40 || (int)k == 72));
            r7Pair++;
            if (report) r7PairPresent++;
            if (report) {
                std::fprintf(stderr,
                             "DIAG_HUBWIRE_PAIR rid=%d row=%s iw=%d i=%zu|%zu sameVtx=%d "
                             "gap3d=%.6f gapUv=%.6f ciA=%d ciB=%d clsA=%s clsB=%s fbA=%d fbB=%d\n",
                             rid, row, iw, j, k, sameVtx, g3, gUv, edges[j].ci, edges[k].ci,
                             edges[j].cls, edges[k].cls, edges[j].fb, edges[k].fb);
            }
            std::fprintf(stderr,
                         "DIAG_HUBWIRE rid=%d row=%s gap3d(%zu,%zu)=%.6f gapUv=%.6f "
                         "sameVtxTShape=%d\n",
                         rid, row, j, k, g3, gUv, sameVtx);
        };
        for (size_t j = 0; j + 1 < edges.size(); j++) gapPair(j, j + 1);
        if (!edges.empty() && (w.Closed() || BRep_Tool::IsClosed(w))) gapPair(edges.size() - 1, 0);

        std::unordered_map<const void*, std::vector<int>> vtxEdges;
        for (size_t k = 0; k < edges.size(); k++) {
            if (edges[k].v0T) vtxEdges[edges[k].v0T].push_back((int)k);
            if (edges[k].v1T) vtxEdges[edges[k].v1T].push_back((int)k);
        }
        std::vector<int> parent(edges.size());
        for (size_t k = 0; k < edges.size(); k++) parent[k] = (int)k;
        auto find = [&](int x) {
            while (parent[(size_t)x] != x) {
                parent[(size_t)x] = parent[(size_t)parent[(size_t)x]];
                x = parent[(size_t)x];
            }
            return x;
        };
        auto unite = [&](int a, int b) {
            a = find(a);
            b = find(b);
            if (a != b) parent[(size_t)a] = b;
        };
        for (size_t a = 0; a < edges.size(); a++) {
            for (size_t b = a + 1; b < edges.size(); b++) {
                if (edges[a].v0T && (edges[a].v0T == edges[b].v0T || edges[a].v0T == edges[b].v1T))
                    unite((int)a, (int)b);
                if (edges[a].v1T && (edges[a].v1T == edges[b].v0T || edges[a].v1T == edges[b].v1T))
                    unite((int)a, (int)b);
            }
        }
        std::unordered_map<int, std::vector<int>> comps;
        for (size_t k = 0; k < edges.size(); k++) comps[find((int)k)].push_back((int)k);
        int compIdx = 0;
        for (const auto& cp : comps) {
            const std::vector<int>& el = cp.second;
            std::unordered_set<const void*> vtx;
            for (int ei : el) {
                if (edges[(size_t)ei].v0T) vtx.insert(edges[(size_t)ei].v0T);
                if (edges[(size_t)ei].v1T) vtx.insert(edges[(size_t)ei].v1T);
            }
            int deg1 = 0;
            for (const void* vt : vtx) {
                int deg = 0;
                auto it = vtxEdges.find(vt);
                if (it != vtxEdges.end()) {
                    for (int ei : it->second)
                        if (std::find(el.begin(), el.end(), ei) != el.end()) deg++;
                }
                if (deg == 1) deg1++;
            }
            const char* kind = (deg1 == 0 && !vtx.empty()) ? "cycle" : "path";
            r7Comp++;
            r7CompPresent++;
            std::fprintf(stderr,
                         "DIAG_HUBWIRE_COMP rid=%d row=%s iw=%d comp=%d kind=%s lenE=%zu nVtx=%zu\n",
                         rid, row, iw, compIdx++, kind, el.size(), vtx.size());
        }
        std::fprintf(stderr,
                     "DIAG_HUBWIRE_SUM rid=%d row=%s iw=%d nE=%zu nComp=%zu nGap3dOverTol=%d "
                     "maxGap3d=%.4f nSharedVtx=%d nStoredPc=%zu\n",
                     rid, row, iw, edges.size(), comps.size(), nGap3dOverTol, maxGap3d, nSharedVtx,
                     edges.size());

        if (!explorer && iw == 0) {
            double maxShareUv = 0;
            int si = -1, sj = -1;
            int nPairsShared = 0, nPairsExamined = 0, nPairsOverTol = 0;
            for (size_t a = 0; a < edges.size(); a++) {
                for (size_t b = a + 1; b < edges.size(); b++) {
                    const void* shared = nullptr;
                    if (edges[a].v0T &&
                        (edges[a].v0T == edges[b].v0T || edges[a].v0T == edges[b].v1T))
                        shared = edges[a].v0T;
                    else if (edges[a].v1T &&
                             (edges[a].v1T == edges[b].v0T || edges[a].v1T == edges[b].v1T))
                        shared = edges[a].v1T;
                    if (!shared) continue;
                    nPairsShared++;
                    gp_Pnt2d ua, ub;
                    int oka = 0, okb = 0;
                    if (edges[a].v0T == shared) {
                        ua = edges[a].uv0;
                        oka = edges[a].haveUv;
                    } else if (edges[a].v1T == shared) {
                        ua = edges[a].uv1;
                        oka = edges[a].haveUv;
                    }
                    if (edges[b].v0T == shared) {
                        ub = edges[b].uv0;
                        okb = edges[b].haveUv;
                    } else if (edges[b].v1T == shared) {
                        ub = edges[b].uv1;
                        okb = edges[b].haveUv;
                    }
                    if (!oka || !okb) continue;
                    nPairsExamined++;
                    r7Share++;
                    r7SharePresent++;
                    const double g = ua.Distance(ub);
                    if (g > tol) nPairsOverTol++;
                    if (g > maxShareUv) {
                        maxShareUv = g;
                        si = (int)a;
                        sj = (int)b;
                    }
                }
            }
            std::fprintf(stderr,
                         "DIAG_SHAREDUV_SUM rid=%d maxGapUv=%.8f i=%d j=%d sewTol=%.8f nE=%zu "
                         "pairsShared=%d pairsExamined=%d pairsOverTol=%d/%d\n",
                         rid, maxShareUv, si, sj, sewTol, edges.size(), nPairsShared,
                         nPairsExamined, nPairsOverTol, std::max(nPairsExamined, 1));
        }
    };

    dumpWireRow(ow, 0, false);
    dumpWireRow(ow, 0, true);
    gHubWireRowB = false;

    // try-B / FClass2d / BRepCheck on isolated copies — live wires stay untouched.
    TopoDS_Wire owCopy;
    std::vector<TopoDS_Wire> innerCopy;
    bool copiesOk = false;
    try {
        const TopoDS_Shape owIso = diagIsolatedCopy(ow);
        if (!owIso.IsNull()) {
            owCopy = TopoDS::Wire(owIso);
            copiesOk = true;
            innerCopy.reserve(innerW.size());
            for (const auto& iw : innerW) {
                const TopoDS_Shape iwIso = diagIsolatedCopy(iw);
                if (iwIso.IsNull()) {
                    copiesOk = false;
                    break;
                }
                innerCopy.push_back(TopoDS::Wire(iwIso));
            }
        }
    } catch (const Standard_Failure&) {
        copiesOk = false;
    }
    TopoDS_Face diagF;
    if (copiesOk && makeFaceKeep(gpl, owCopy, innerCopy, outward, diagF) && !diagF.IsNull()) {
        int iiw = 0;
        for (TopExp_Explorer wx(diagF, TopAbs_WIRE); wx.More(); wx.Next(), iiw++) {
            const TopoDS_Wire w = TopoDS::Wire(wx.Current());
            int nEw = 0;
            for (TopoDS_Iterator it(w); it.More(); it.Next())
                if (it.Value().ShapeType() == TopAbs_EDGE) nEw++;
            const char* infPt = "?";
            try {
                BRepTopAdaptor_FClass2d cl(diagF, Precision::PConfusion());
                infPt = topAbsStateName(cl.PerformInfinitePoint());
            } catch (const Standard_Failure&) {
            }
            r7Inner++;
            r7InnerPresent++;
            std::fprintf(stderr, "DIAG_HUBWIRE_INNER rid=%d iw=%d outer=%d infPt=%s nE=%d\n", rid, iiw,
                         w.IsSame(owCopy) ? 1 : 0, infPt, nEw);
        }
        if (rid == 0 || rid == 2 || rid == 3) {
            try {
                BRepCheck_Wire chk(owCopy);
                chk.GeometricControls(Standard_True);
                const BRepCheck_Status c2d = chk.Closed2d(diagF, Standard_False);
                const BRepCheck_Status ori = chk.Orientation(diagF, Standard_False);
                BRepCheck_Face fc(diagF);
                fc.GeometricControls(Standard_True);
                r7Orient++;
                r7OrientPresent++;
                std::fprintf(stderr,
                             "DIAG_HUBORIENT rid=%d closed2d=%s wireOrient=%s intersectWires=%s "
                             "imbrication=%s orientationOfWires=%s\n",
                             rid, brepCheckName((int)c2d), brepCheckName((int)ori),
                             brepCheckName((int)fc.IntersectWires(Standard_False)),
                             brepCheckName((int)fc.ClassifyWires(Standard_False)),
                             brepCheckName((int)fc.OrientationOfWires(Standard_False)));
            } catch (const Standard_Failure&) {
            }
        }
    }

    std::fprintf(stderr,
                 "DIAG_OPEN2D_R7 rid=%d hubwire=%d/%d shareduv=%d/%d comp=%d/%d pair=%d/%d "
                 "inner=%d/%d orient=%d/%d\n",
                 rid, r7HubwirePresent, r7Hubwire, r7SharePresent, r7Share, r7CompPresent, r7Comp,
                 r7PairPresent, r7Pair, r7InnerPresent, r7Inner, r7OrientPresent, r7Orient);
}
void dumpDiagWirePop(const Region& r, const TopoDS_Wire& ow, const RegionSet& rs,
                     const MeshView& mv, const std::vector<ChainGeom>& geom,
                     const std::vector<TopoDS_Edge>& meshE) {
    if (!diagP2Enabled() || ow.IsNull()) return;
    const char* part = sliverCensusPartName((int)mv.nTri);
    const int rid = r.id;
    const double sewTol = mv.sewTol;
    auto chainByTShape = [&](const TopoDS_Edge& e) -> int {
        const void* ts = diagTShapePtr(e);
        for (size_t ci = 0; ci < geom.size(); ci++) {
            for (const TopoDS_Edge& ge : geom[ci].edges) {
                if (!ge.IsNull() && (e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts)))
                    return (int)ci;
            }
        }
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            for (int eid : rs.chains[ci].meshEdges) {
                if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
                if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts))
                    return (int)ci;
            }
        }
        return -1;
    };
    auto plateSide = [&](int ci) -> int {
        if (ci < 0 || (size_t)ci >= rs.chains.size()) return 0;
        const BoundaryChain& ch = rs.chains[(size_t)ci];
        return (ch.regA == rid || ch.regB == rid) ? 1 : 0;
    };

    std::vector<TopoDS_Edge> elist;
    for (TopoDS_Iterator it(ow); it.More(); it.Next()) {
        if (it.Value().ShapeType() == TopAbs_EDGE) elist.push_back(TopoDS::Edge(it.Value()));
    }
    const int nE = (int)elist.size();
    int nChainE = 0, nNonChainE = 0, nSeam = 0;
    int r7E = 0, r7Eok = 0, r7Drop = 0, r7DropOk = 0, r7Cov = 0, r7CovOk = 0;
    std::string sites;
    struct Rec {
        int i, matchCi, originCi, partner, cylRid, closed360, seam, plateSide, fb;
        const void *eT, *v0T, *v1T;
        const char *site, *cls;
        gp_Pnt p0, p1;
        TopoDS_Edge e;
    };
    std::vector<Rec> recs;
    recs.reserve(elist.size());
    std::unordered_set<const void*> chainVtx;
    for (const TopoDS_Edge& e : elist) {
        TopoDS_Vertex v0, v1;
        TopExp::Vertices(e, v0, v1, Standard_True);
        const int mci = chainByTShape(e);
        if (mci >= 0) {
            chainVtx.insert(diagTShapePtr(v0));
            chainVtx.insert(diagTShapePtr(v1));
        }
    }
    for (int i = 0; i < nE; i++) {
        const TopoDS_Edge& e = elist[(size_t)i];
        TopoDS_Vertex v0, v1;
        TopExp::Vertices(e, v0, v1, Standard_True);
        Rec x;
        x.i = i;
        x.e = e;
        x.eT = diagTShapePtr(e);
        x.v0T = diagTShapePtr(v0);
        x.v1T = diagTShapePtr(v1);
        x.p0 = BRep_Tool::Pnt(v0);
        x.p1 = BRep_Tool::Pnt(v1);
        x.matchCi = chainByTShape(e);
        x.originCi = -1;
        x.site = "?";
        if (x.eT && gWirePopCi.count(x.eT)) x.originCi = gWirePopCi[x.eT];
        if (x.eT && gWirePopSite.count(x.eT)) x.site = gWirePopSite[x.eT];
        x.seam = (x.eT && gSeamTShapes.count(x.eT)) ? 1 : 0;
        x.partner = -1;
        x.cylRid = x.seam ? seamCylinderRid(e, rs, geom) : -1;
        x.closed360 = 0;
        x.fb = 0;
        x.plateSide = plateSide(x.matchCi);
        if (x.matchCi >= 0 && (size_t)x.matchCi < rs.chains.size()) {
            const BoundaryChain& ch = rs.chains[(size_t)x.matchCi];
            x.partner = (ch.regA == rid) ? ch.regB : ch.regA;
            x.fb = gChainSewFallbackCi.count(x.matchCi) ? 1 : 0;
        }
        if (x.cylRid >= 0) {
            const Region* cr = regionById(rs, x.cylRid);
            if (cr) x.closed360 = cr->closed360 ? 1 : 0;
        }
        x.cls = x.seam ? "seam" : (x.matchCi >= 0 ? "chain" : "nonchain");
        if (x.matchCi >= 0) nChainE++;
        else nNonChainE++;
        if (x.seam) nSeam++;
        recs.push_back(x);
        if (sites.find(x.site) == std::string::npos) {
            if (!sites.empty()) sites += ",";
            sites += x.site;
        }
        const int vShare = ((x.v0T && chainVtx.count(x.v0T)) || (x.v1T && chainVtx.count(x.v1T)))
                               ? 1
                               : 0;
        if (x.matchCi < 0 || (rid >= 0 && rid <= 3 && x.seam)) {
            r7E++;
            r7Eok++;
            std::fprintf(stderr,
                         "DIAG_WIREPOP_E rid=%d iw=0 i=%d eT=%p v0T=%p v1T=%p cls=%s "
                         "matchCi=%d originCi=%d site=%s partner=%d plateSide=%d seam=%d "
                         "cylRid=%d closed360=%d fb=%d vShareChain=%d "
                         "p0=(%.4f,%.4f,%.4f) p1=(%.4f,%.4f,%.4f)\n",
                         rid, i, x.eT, x.v0T, x.v1T, x.cls, x.matchCi, x.originCi, x.site,
                         x.partner, x.plateSide, x.seam, x.cylRid, x.closed360, x.fb, vShare,
                         x.p0.X(), x.p0.Y(), x.p0.Z(), x.p1.X(), x.p1.Y(), x.p1.Z());
        }
    }
    if (rid == 2 && nE > 40) {
        const Rec& a = recs[39];
        const Rec& b = recs[40];
        const double gOri = a.p1.Distance(b.p0);
        const double dA0B0 = a.p0.Distance(b.p0);
        const double dA0B1 = a.p0.Distance(b.p1);
        const double dA1B0 = a.p1.Distance(b.p0);
        const double dA1B1 = a.p1.Distance(b.p1);
        const int anyShare = (a.v0T == b.v0T || a.v0T == b.v1T || a.v1T == b.v0T || a.v1T == b.v1T)
                                 ? 1
                                 : 0;
        const int seqShare = (a.v1T && a.v1T == b.v0T) ? 1 : 0;
        std::fprintf(stderr,
                     "DIAG_WIREPOP_GAP rid=2 i=39|40 gap3d_oriented=%.6f d_v1_v0=%.6f "
                     "dA0B0=%.6f dA0B1=%.6f dA1B0=%.6f dA1B1=%.6f anyShareVtx=%d seqShare=%d "
                     "matchCi=%d|%d originCi=%d|%d site=%s|%s seam=%d|%d cylRid=%d|%d\n",
                     gOri, dA1B0, dA0B0, dA0B1, dA1B0, dA1B1, anyShare, seqShare, a.matchCi,
                     b.matchCi, a.originCi, b.originCi, a.site, b.site, a.seam, b.seam, a.cylRid,
                     b.cylRid);
    }

    // Removal test: drop TShape-unmatched (true non-chain) edges.
    std::vector<int> kept;
    for (const Rec& x : recs)
        if (x.matchCi >= 0) kept.push_back(x.i);
        else {
            r7Drop++;
            r7DropOk++;
        }
    std::unordered_map<const void*, int> deg;
    std::unordered_map<const void*, const void*> par;
    auto fnd = [&](const void* x) {
        if (!par.count(x)) par[x] = x;
        const void* r0 = x;
        while (par[r0] != r0) r0 = par[r0];
        return r0;
    };
    auto uni = [&](const void* a, const void* b) {
        if (!a || !b) return;
        a = fnd(a);
        b = fnd(b);
        if (a != b) par[a] = b;
    };
    for (int idx : kept) {
        const Rec& x = recs[(size_t)idx];
        if (x.v0T) deg[x.v0T]++;
        if (x.v1T) deg[x.v1T]++;
        uni(x.v0T, x.v1T);
    }
    std::unordered_set<const void*> roots;
    int nDeg1 = 0, nDeg2 = 0, nOdd = 0;
    for (const auto& kv : deg) {
        roots.insert(fnd(kv.first));
        if (kv.second == 1) nDeg1++;
        else if (kv.second == 2) nDeg2++;
        if (kv.second % 2) nOdd++;
    }
    const int nComp = (int)roots.size();
    const int tshapeCycle = (nComp == 1 && nDeg1 == 0 && nOdd == 0 && (int)kept.size() >= 3) ? 1 : 0;
    std::fprintf(stderr,
                 "DIAG_WIREPOP_REMOVE rid=%d nE_before=%d nE_after=%d nDropped=%d nComp=%d "
                 "nDeg1=%d nDeg2=%d tshapeCycle=%d\n",
                 rid, nE, (int)kept.size(), nE - (int)kept.size(), nComp, nDeg1, nDeg2,
                 tshapeCycle);
    if (rid == 2) {
        std::vector<int> keptS;
        std::unordered_map<const void*, int> degS;
        std::unordered_map<const void*, const void*> parS;
        auto fndS = [&](const void* x) {
            if (!x) return x;
            if (!parS.count(x)) parS[x] = x;
            const void* r0 = x;
            while (parS[r0] != r0) r0 = parS[r0];
            return r0;
        };
        auto uniS = [&](const void* a, const void* b) {
            if (!a || !b) return;
            a = fndS(a);
            b = fndS(b);
            if (a != b) parS[a] = b;
        };
        for (const Rec& x : recs) {
            if (x.seam) continue;
            keptS.push_back(x.i);
            if (x.v0T) degS[x.v0T]++;
            if (x.v1T) degS[x.v1T]++;
            uniS(x.v0T, x.v1T);
        }
        int nDeg1S = 0, nCompS = 0;
        std::unordered_set<const void*> rootsS;
        for (const auto& kv : degS) {
            rootsS.insert(fndS(kv.first));
            if (kv.second == 1) nDeg1S++;
        }
        nCompS = (int)rootsS.size();
        const int cycS = (nCompS == 1 && nDeg1S == 0 && (int)keptS.size() >= 3) ? 1 : 0;
        std::fprintf(stderr,
                     "DIAG_WIREPOP_SEAMDROP rid=2 nE_before=%d nSeam=%d nE_after=%d nComp=%d "
                     "nDeg1=%d tshapeCycle=%d\n",
                     nE, nSeam, (int)keptS.size(), nCompS, nDeg1S, cycS);
    }
    if (!tshapeCycle && rid == 2) {
        for (size_t k = 0; k < kept.size(); k++) {
            const Rec& a = recs[(size_t)kept[k]];
            const Rec& b = recs[(size_t)kept[(k + 1) % kept.size()]];
            const int seq = (a.v1T && a.v1T == b.v0T) ? 1 : 0;
            const int any = (a.v0T == b.v0T || a.v0T == b.v1T || a.v1T == b.v0T || a.v1T == b.v1T)
                                ? 1
                                : 0;
            if (seq) continue;
            const double d = a.p1.Distance(b.p0);
            std::fprintf(stderr,
                         "DIAG_WIREPOP_RESIDGAP rid=%d i=%d|%d seqShare=%d anyShare=%d d3d=%.6f "
                         "v1T=%p v0T=%p matchCi=%d|%d cls=%s|%s\n",
                         rid, a.i, b.i, seq, any, d, a.v1T, b.v0T, a.matchCi, b.matchCi, a.cls,
                         b.cls);
        }
        for (const auto& kv : deg) {
            if (kv.second == 2) continue;
            std::fprintf(stderr, "DIAG_WIREPOP_RESIDV rid=%d vT=%p deg=%d\n", rid, kv.first,
                         kv.second);
        }
    }

    // Covering chains for dropped stretches (hub 2).
    if (rid == 2) {
        std::unordered_set<int> inWire;
        for (const Rec& x : recs)
            if (x.matchCi >= 0) inWire.insert(x.matchCi);
        for (const Rec& x : recs) {
            if (x.matchCi >= 0) continue;
            r7Cov++;
            r7CovOk++;
            int nSpan = 0;
            std::string spanList;
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                const BoundaryChain& ch = rs.chains[ci];
                if (ch.regA != rid && ch.regB != rid) continue;
                if (ch.meshVerts.size() < 2) continue;
                gp_Pnt a = pntOf(mv, ch.meshVerts.front());
                gp_Pnt b = pntOf(mv, ch.meshVerts.back());
                const double da = std::min(x.p0.Distance(a), x.p0.Distance(b));
                const double db = std::min(x.p1.Distance(a), x.p1.Distance(b));
                if (da > 2.0 * sewTol && db > 2.0 * sewTol) continue;
                nSpan++;
                if (!spanList.empty()) spanList += ",";
                char buf[64];
                std::snprintf(buf, sizeof(buf), "%d%s", (int)ci, inWire.count((int)ci) ? "*" : "");
                spanList += buf;
            }
            std::fprintf(stderr,
                         "DIAG_WIREPOP_COVER rid=%d i=%d nSpan=%d inWireMarked=* chains=[%s] "
                         "cylRid=%d site=%s\n",
                         rid, x.i, nSpan, spanList.c_str(), x.cylRid, x.site);
        }
    }

    if (sites.empty()) sites = "-";
    std::fprintf(stderr,
                 "DIAG_WIREPOP part=%s rid=%d iw=0 nE=%d nChainE=%d nNonChainE=%d nSeam=%d "
                 "sites=[%s] r7e=%d/%d r7drop=%d/%d r7cov=%d/%d\n",
                 part, rid, nE, nChainE, nNonChainE, nSeam, sites.c_str(), r7Eok, r7E, r7DropOk,
                 r7Drop, r7CovOk, r7Cov);
}

Handle(Geom2d_Curve) basisCurve2dOf(const Handle(Geom2d_Curve)& c) {
    Handle(Geom2d_TrimmedCurve) tr = Handle(Geom2d_TrimmedCurve)::DownCast(c);
    if (!tr.IsNull() && !tr->BasisCurve().IsNull()) return tr->BasisCurve();
    return c;
}

int countPcurvesOnFaceSurf(const TopoDS_Edge& e, const Handle(Geom_Surface)& sF) {
    int n = 0;
    if (e.IsNull() || sF.IsNull()) return 0;
    try {
        Handle(BRep_TEdge) te = Handle(BRep_TEdge)::DownCast(e.TShape());
        if (te.IsNull()) return 0;
        const void* want = static_cast<const void*>(sF.get());
        for (BRep_ListOfCurveRepresentation::Iterator it(te->Curves()); it.More(); it.Next()) {
            Handle(BRep_CurveOnSurface) cos = Handle(BRep_CurveOnSurface)::DownCast(it.Value());
            if (cos.IsNull() || cos->Surface().IsNull()) continue;
            if (static_cast<const void*>(cos->Surface().get()) == want) n++;
        }
    } catch (const Standard_Failure&) {
    }
    return n;
}

double maxDevPcurveVs3d(const TopoDS_Face& f, const TopoDS_Edge& e) {
    double maxD = -1.0;
    try {
        TopLoc_Location loc;
        Handle(Geom_Surface) srf = BRep_Tool::Surface(f, loc);
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        Standard_Real a = 0, b = 0;
        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
        if (c3.IsNull() || pc.IsNull() || srf.IsNull() || l3 <= f3) return maxD;
        const gp_Trsf tr = loc.Transformation();
        for (int i = 0; i < 16; i++) {
            const double u = (double)i / 15.0;
            const gp_Pnt p3 = c3->Value(f3 + u * (l3 - f3));
            const gp_Pnt2d uv = pc->Value(a + u * (b - a));
            gp_Pnt ps = srf->Value(uv.X(), uv.Y());
            if (!loc.IsIdentity()) ps.Transform(tr);
            const double d = p3.Distance(ps);
            if (d > maxD) maxD = d;
        }
    } catch (const Standard_Failure&) {
        return -1.0;
    }
    return maxD;
}

thread_local std::unordered_map<const void*, double> gHpPlateEdgeTol;

void noteHpPlateEdgeTols(const TopoDS_Face& f, int rid) {
    if (f.IsNull() || rid < 0 || rid > 3) return;
    try {
        TopTools_IndexedMapOfShape em;
        TopExp::MapShapes(f, TopAbs_EDGE, em);
        for (int i = 1; i <= em.Extent(); i++) {
            const TopoDS_Edge e = TopoDS::Edge(em(i));
            const void* ts = diagTShapePtr(e);
            if (ts) gHpPlateEdgeTol[ts] = BRep_Tool::Tolerance(e);
        }
    } catch (const Standard_Failure&) {
    }
}

double planeSignedDist(const gp_Pln& pl, const gp_Pnt& p) {
    Standard_Real a = 0, b = 0, c = 0, d = 0;
    pl.Coefficients(a, b, c, d);
    return a * p.X() + b * p.Y() + c * p.Z() + d;
}

bool facePlanePln(const TopoDS_Face& f, gp_Pln& out) {
    TopLoc_Location loc;
    Handle(Geom_Plane) gpl = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(f, loc));
    if (gpl.IsNull()) return false;
    out = gpl->Pln();
    if (!loc.IsIdentity()) out.Transform(loc.Transformation());
    return true;
}

// D-S3-111: affine → |C(t)-S(P(t))| at endpoints (not vertex-to-plane); circle →
// |offset|+R·sin(tilt); polyline → max over poles with pcurve when present.
double exactMaxOnPlane(const TopoDS_Face& f, const TopoDS_Edge& e, const char* c3name) {
    gp_Pln pl;
    if (!facePlanePln(f, pl)) return -1.0;
    Standard_Real f3 = 0, l3 = 0;
    Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
    Standard_Real a = 0, b = 0;
    Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
    TopLoc_Location loc;
    Handle(Geom_Surface) srf = BRep_Tool::Surface(f, loc);
    double m = -1.0;
    if (!c3.IsNull() && !pc.IsNull() && !srf.IsNull() && l3 > f3)
        m = std::max(curveSurfDevAtBind(c3, f3, srf, pc, loc),
                     curveSurfDevAtBind(c3, l3, srf, pc, loc));
    if (c3name && std::strcmp(c3name, "Circle") == 0) {
        try {
            Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(basisCurveOf(c3));
            if (!gc.IsNull()) {
                const double circMax = circleOnPlaneMax(pl, gc->Circ());
                return (m >= 0.0) ? std::max(m, circMax) : circMax;
            }
        } catch (const Standard_Failure&) {
        }
    }
    if (c3name && std::strncmp(c3name, "Polyline", 8) == 0) {
        try {
            Handle(Geom_BSplineCurve) bs =
                Handle(Geom_BSplineCurve)::DownCast(basisCurveOf(c3));
            if (!bs.IsNull()) {
                if (m < 0.0) m = 0.0;
                for (int i = 1; i <= bs->NbPoles(); i++) {
                    if (!pc.IsNull() && !srf.IsNull()) {
                        const double t = (bs->NbPoles() > 1)
                                             ? f3 + (l3 - f3) * (double)(i - 1) /
                                                      (double)(bs->NbPoles() - 1)
                                             : f3;
                        m = std::max(m, curveSurfDevAtBind(c3, t, srf, pc, loc));
                    } else
                        m = std::max(m, std::fabs(planeSignedDist(pl, bs->Pole(i))));
                }
            }
        } catch (const Standard_Failure&) {
        }
    }
    if (m >= 0.0) return m;
    TopoDS_Vertex v0, v1;
    TopExp::Vertices(e, v0, v1, Standard_True);
    const gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
    const gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
    return std::max(std::fabs(planeSignedDist(pl, p0)), std::fabs(planeSignedDist(pl, p1)));
}

// Replicate BRepCheck_Edge::InContext's BRepLib_ValidateEdge sampling.
// NCONTROL in BRepCheck_Edge.cxx is 23; InContext uses ValidateEdge default 22
// unless SetControlPointsNumber is called. We set 23 to match NCONTROL and do
// not SetExitIfToleranceExceeded, so brepMax is the grid maximum not a first-hit.
double brepCheckInContextMax(const TopoDS_Face& f, const TopoDS_Edge& e) {
    double maxD = -1.0;
    try {
        TopLoc_Location loc;
        Handle(Geom_Surface) srf = BRep_Tool::Surface(f, loc);
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        Standard_Real a = 0, b = 0;
        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, f, a, b);
        if (c3.IsNull() || pc.IsNull() || srf.IsNull() || l3 <= f3) return -1.0;
        Handle(Geom_Surface) sb =
            Handle(Geom_Surface)::DownCast(srf->Transformed(loc.Transformation()));
        if (sb.IsNull()) sb = srf;
        const gp_Trsf tr = loc.Transformation();
        const int N = 23;
        for (int i = 0; i < N; i++) {
            const double u = (N > 1) ? (double)i / (double)(N - 1) : 0.0;
            const gp_Pnt p3 = c3->Value(f3 + u * (l3 - f3));
            const gp_Pnt2d uv = pc->Value(a + u * (b - a));
            gp_Pnt ps = sb->Value(uv.X(), uv.Y());
            if (!loc.IsIdentity()) ps.Transform(tr);
            const double d = p3.Distance(ps);
            if (d > maxD) maxD = d;
        }
    } catch (const Standard_Failure&) {
        return -1.0;
    }
    return maxD;
}

double brepDevMaxOnSurf(const Handle(Geom_Curve)& c3, Standard_Real f3, Standard_Real l3,
                        const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& pc,
                        Standard_Real a, Standard_Real b, const TopLoc_Location& loc) {
    double maxD = -1.0;
    try {
        if (c3.IsNull() || pc.IsNull() || srf.IsNull() || l3 <= f3) return -1.0;
        Handle(Geom_Surface) sb =
            Handle(Geom_Surface)::DownCast(srf->Transformed(loc.Transformation()));
        if (sb.IsNull()) sb = srf;
        const gp_Trsf tr = loc.Transformation();
        const int N = 23;
        for (int i = 0; i < N; i++) {
            const double u = (N > 1) ? (double)i / (double)(N - 1) : 0.0;
            const gp_Pnt p3 = c3->Value(f3 + u * (l3 - f3));
            const gp_Pnt2d uv = pc->Value(a + u * (b - a));
            gp_Pnt ps = sb->Value(uv.X(), uv.Y());
            if (!loc.IsIdentity()) ps.Transform(tr);
            const double d = p3.Distance(ps);
            if (d > maxD) maxD = d;
        }
    } catch (const Standard_Failure&) {
        return -1.0;
    }
    return maxD;
}

double uvClosureGapOnFace(const TopoDS_Face& f, const TopoDS_Wire& w) {
    double gap = 0.0;
    gp_Pnt2d firstUv, prevUv;
    bool havePrev = false, haveFirst = false;
    try {
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            const TopoDS_Edge e = TopoDS::Edge(it.Value());
            gp_Pnt2d uv0, uv1;
            int stored = 0;
            if (!edgeUvOnPlane(f, e, uv0, uv1, stored)) continue;
            if (!haveFirst) {
                firstUv = uv0;
                haveFirst = true;
            }
            if (havePrev) {
                const double g = prevUv.Distance(uv0);
                if (g > gap) gap = g;
            }
            prevUv = uv1;
            havePrev = true;
        }
        if (havePrev && haveFirst) {
            const double g = prevUv.Distance(firstUv);
            if (g > gap) gap = g;
        }
    } catch (const Standard_Failure&) {
    }
    return gap;
}

bool inContextHas8or11(const Handle(BRepCheck_Result)& res, const TopoDS_Face& f) {
    auto scan = [](const BRepCheck_ListOfStatus& lst) {
        for (BRepCheck_ListOfStatus::Iterator it(lst); it.More(); it.Next()) {
            if (it.Value() == BRepCheck_InvalidCurveOnSurface ||
                it.Value() == BRepCheck_InvalidSameParameterFlag)
                return true;
        }
        return false;
    };
    if (res.IsNull()) return false;
    try {
        if (res->IsStatusOnShape(f)) return scan(res->StatusOnShape(f));
        res->InitContextIterator();
        for (; res->MoreShapeInContext(); res->NextShapeInContext()) {
            if (res->ContextualShape().IsSame(f) || res->ContextualShape().IsEqual(f))
                return scan(res->StatusOnShape());
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

void edgeExactPlaneAndOther(const TopoDS_Face& f, const TopoDS_Edge& e, double& exactPlane,
                            double& exactOther, double& brepPlane, double& brepOther) {
    exactPlane = exactOther = brepPlane = brepOther = -1.0;
    try {
        TopLoc_Location locF;
        Handle(Geom_Surface) sF = BRep_Tool::Surface(f, locF);
        const void* want = sF.IsNull() ? nullptr : static_cast<const void*>(sF.get());
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        Handle(BRep_TEdge) te = Handle(BRep_TEdge)::DownCast(e.TShape());
        if (te.IsNull() || c3.IsNull()) return;
        for (BRep_ListOfCurveRepresentation::Iterator it(te->Curves()); it.More(); it.Next()) {
            Handle(BRep_CurveOnSurface) cos = Handle(BRep_CurveOnSurface)::DownCast(it.Value());
            if (cos.IsNull() || cos->Surface().IsNull() || cos->PCurve().IsNull()) continue;
            const double ex = exactOrSampledOnPair(c3, f3, l3, cos->Surface(), cos->PCurve(),
                                                    cos->Location());
            const double br = brepDevMaxOnSurf(c3, f3, l3, cos->Surface(), cos->PCurve(), f3, l3,
                                                 cos->Location());
            const void* ph = static_cast<const void*>(cos->Surface().get());
            if (want && ph == want) {
                exactPlane = ex;
                brepPlane = br;
            } else {
                if (ex > exactOther) exactOther = ex;
                if (br > brepOther) brepOther = br;
            }
        }
    } catch (const Standard_Failure&) {
    }
}

// D-S3-99 delta: classify in-context 8/11 edges. Print-only, try A only.
void dumpDiagPcInvalid(int rid, const TopoDS_Face& f, const TopoDS_Wire& ow, const RegionSet& rs,
                       const MeshView& mv, const std::vector<ChainGeom>& geom,
                       const std::vector<char>& collapsed, const std::vector<TopoDS_Edge>& meshE) {
    if (f.IsNull() || rid < 0) return;
    static thread_local std::unordered_set<int> dumped;
    if (!dumped.insert(rid).second) return;
    const double tol = std::max(mv.sewTol, Precision::PConfusion());

    auto matchChain = [&](const TopoDS_Edge& e, const gp_Pnt& p0, const gp_Pnt& p1) -> int {
        const void* ts = diagTShapePtr(e);
        for (size_t ci = 0; ci < geom.size(); ci++) {
            for (const TopoDS_Edge& ge : geom[ci].edges) {
                if (!ge.IsNull() && (e.IsSame(ge) || (ts && diagTShapePtr(ge) == ts)))
                    return (int)ci;
            }
        }
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            for (int eid : rs.chains[ci].meshEdges) {
                if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
                if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts))
                    return (int)ci;
            }
        }
        int best = -1;
        double bestD = 1e300;
        for (size_t ci = 0; ci < rs.chains.size(); ci++) {
            const BoundaryChain& ch = rs.chains[ci];
            if (ch.meshVerts.size() < 2) continue;
            gp_Pnt a = pntOf(mv, ch.meshVerts.front());
            gp_Pnt b = pntOf(mv, ch.meshVerts.back());
            const double d = std::min(p0.Distance(a) + p1.Distance(b), p0.Distance(b) + p1.Distance(a));
            if (d < bestD) {
                bestD = d;
                best = (int)ci;
            }
        }
        return (best >= 0 && bestD <= 2.0 * tol) ? best : -1;
    };

    // D-130 tripwire support: an analytic (bindAllVariants) edge must still name
    // its chain. Identity first -- a collapsed chain owns its edge TShapes -- then
    // the endpoint match. Reporting ci=-1 for every seam edge (as this classifier
    // used to, by returning before ci was computed) made the census unreadable and
    // hid the edges that really have no chain.
    auto chainOfCollapsedEdge = [&](const TopoDS_Edge& e) -> int {
        const void* ts = diagTShapePtr(e);
        if (!ts) return -1;
        for (size_t ci = 0; ci < geom.size(); ci++) {
            if (!geom[ci].collapsed) continue;
            for (const auto& ge : geom[ci].edges)
                if (!ge.IsNull() && diagTShapePtr(ge) == ts) return (int)ci;
        }
        return -1;
    };
    auto classifyEdge = [&](const TopoDS_Edge& e, int& ci, int& partner) -> const char* {
        ci = -1;
        partner = -1;
        const void* ts = diagTShapePtr(e);
        if (ts && gSeamTShapes.count(ts)) {
            ci = chainOfCollapsedEdge(e);
            if (ci < 0) {
                TopoDS_Vertex s0, s1;
                TopExp::Vertices(e, s0, s1, Standard_True);
                ci = matchChain(e, s0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(s0),
                                s1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(s1));
            }
            if (ci >= 0 && (size_t)ci < rs.chains.size()) {
                const BoundaryChain& ch = rs.chains[(size_t)ci];
                partner = (ch.regA == rid) ? ch.regB : ch.regA;
                if (partner < 0) partner = ch.regA >= 0 ? ch.regA : ch.regB;
            }
            return "seam";
        }
        TopoDS_Vertex v0, v1;
        TopExp::Vertices(e, v0, v1, Standard_True);
        gp_Pnt p0 = v0.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v0);
        gp_Pnt p1 = v1.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(v1);
        ci = matchChain(e, p0, p1);
        bool isMeshE = false;
        if (ci >= 0 && (size_t)ci < rs.chains.size()) {
            for (int eid : rs.chains[(size_t)ci].meshEdges) {
                if (eid < 0 || (size_t)eid >= meshE.size() || meshE[(size_t)eid].IsNull()) continue;
                if (e.IsSame(meshE[(size_t)eid]) || (ts && diagTShapePtr(meshE[(size_t)eid]) == ts)) {
                    isMeshE = true;
                    break;
                }
            }
            const BoundaryChain& ch = rs.chains[(size_t)ci];
            partner = (ch.regA == rid) ? ch.regB : ch.regA;
            if (partner < 0) partner = ch.regA >= 0 ? ch.regA : ch.regB;
            if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                !geom[(size_t)ci].edges.empty()) {
                const Region* rr = regionById(rs, partner);
                if (rr && rr->type == SurfType::Cylinder) return "cyl";
                if (rr && rr->type == SurfType::Plane) return "pln";
                return "span";
            }
            if (isMeshE) return "meshE";
            const Region* rr = regionById(rs, partner);
            if (rr && rr->type == SurfType::Cylinder) return "cyl";
            if (rr && rr->type == SurfType::Plane) return "pln";
            return "fill";
        }
        return "fill";
    };

    auto curve3dName = [&](const TopoDS_Edge& e, char* buf, size_t n) {
        Standard_Real f3 = 0, l3 = 0;
        Handle(Geom_Curve) c3 = BRep_Tool::Curve(e, f3, l3);
        if (c3.IsNull()) {
            std::snprintf(buf, n, "other");
            return;
        }
        Handle(Geom_Curve) src = basisCurveOf(c3);
        if (src.IsNull()) src = c3;
        if (src->DynamicType() == STANDARD_TYPE(Geom_Line)) {
            std::snprintf(buf, n, "Line");
            return;
        }
        if (src->DynamicType() == STANDARD_TYPE(Geom_Circle)) {
            std::snprintf(buf, n, "Circle");
            return;
        }
        Handle(Geom_BSplineCurve) bs = Handle(Geom_BSplineCurve)::DownCast(src);
        if (!bs.IsNull()) {
            const int nSeg = std::max(1, bs->NbPoles() - 1);
            if (bs->Degree() <= 1)
                std::snprintf(buf, n, "Polyline:%d", nSeg);
            else
                std::snprintf(buf, n, "BSpline");
            return;
        }
        std::snprintf(buf, n, "other");
    };

    auto pcurveKindName = [&](const TopoDS_Face& face, const TopoDS_Edge& e, char* buf, size_t n) {
        try {
            TopLoc_Location loc;
            Handle(Geom_Plane) pl = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(face, loc));
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, face, a, b);
            if (pc.IsNull()) {
                std::snprintf(buf, n, "none");
                return;
            }
            Handle(Geom2d_Curve) src = basisCurve2dOf(pc);
            if (src.IsNull()) src = pc;
            if (src->DynamicType() == STANDARD_TYPE(Geom2d_Circle) ||
                src->DynamicType() == STANDARD_TYPE(Geom2d_Ellipse)) {
                std::snprintf(buf, n, "Geom2d_Circle/To2d");
                return;
            }
            if (src->DynamicType() == STANDARD_TYPE(Geom2d_BSplineCurve)) {
                std::snprintf(buf, n, "Geom2d_BSpline");
                return;
            }
            if (src->DynamicType() == STANDARD_TYPE(Geom2d_Line)) {
                bool pinned = false;
                if (!pl.IsNull()) {
                    TopoDS_Vertex vf, vl;
                    TopExp::Vertices(e, vf, vl, Standard_False);
                    if (!vf.IsNull() && !vl.IsNull()) {
                        Standard_Real u0 = 0, v0 = 0, u1 = 0, v1 = 0;
                        ElSLib::Parameters(pl->Pln(), BRep_Tool::Pnt(vf), u0, v0);
                        ElSLib::Parameters(pl->Pln(), BRep_Tool::Pnt(vl), u1, v1);
                        const gp_Pnt2d uvA = pc->Value(a), uvB = pc->Value(b);
                        const gp_Pnt2d vtx0(u0, v0), vtx1(u1, v1);
                        const double dFwd = uvA.Distance(vtx0) + uvB.Distance(vtx1);
                        const double dRev = uvA.Distance(vtx1) + uvB.Distance(vtx0);
                        if (std::min(dFwd, dRev) <= 1e-6) pinned = true;
                    }
                }
                std::snprintf(buf, n, pinned ? "Geom2d_Line(pinned)" : "projected");
                return;
            }
            std::snprintf(buf, n, "projected");
        } catch (const Standard_Failure&) {
            std::snprintf(buf, n, "none");
        }
    };

    auto birthSite = [&](const char* cls, int ci) -> const char* {
        if (std::strcmp(cls, "meshE") == 0) return "bindEdgePcurveOnInternedPlane";
        const bool col = (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                          (size_t)ci < geom.size() && geom[(size_t)ci].collapsed);
        if (col) return "bindAllVariants";
        if (std::strcmp(cls, "seam") == 0) return "makeFaceKeep-seam";
        if (std::strcmp(cls, "cyl") == 0 || std::strcmp(cls, "pln") == 0 ||
            std::strcmp(cls, "span") == 0)
            return "bindAllVariants";
        return "unknown";
    };

    try {
        BRepCheck_Analyzer an(f, Standard_True);
        TopLoc_Location loc;
        Handle(Geom_Surface) sF = BRep_Tool::Surface(f, loc);
        TopTools_IndexedMapOfShape emap;
        TopExp::MapShapes(f, TopAbs_EDGE, emap);
        std::map<std::string, int> hist;
        int nFlagged = 0;
        struct Row {
            int ie;
            const char* cls;
            int ci;
            int partner;
            char c3[32];
            char pk[40];
            const char* birth;
            int nPc;
            double dev;
            int sp;
            int sr;
            double etol;
            bool flagged;
            TopoDS_Edge e;
        };
        std::vector<Row> rows;
        rows.reserve((size_t)emap.Extent());
        for (int ie = 1; ie <= emap.Extent(); ie++) {
            const TopoDS_Edge e = TopoDS::Edge(emap(ie));
            Handle(BRepCheck_Result) res = analyzerResultOf(an, e);
            const bool flagged = inContextHas8or11(res, f);
            int ci = -1, partner = -1;
            const char* cls = classifyEdge(e, ci, partner);
            Row r;
            r.ie = ie - 1;
            r.cls = cls;
            r.ci = ci;
            r.partner = partner;
            curve3dName(e, r.c3, sizeof(r.c3));
            pcurveKindName(f, e, r.pk, sizeof(r.pk));
            r.birth = birthSite(cls, ci);
            r.nPc = countPcurvesOnFaceSurf(e, sF);
            r.dev = maxDevPcurveVs3d(f, e);
            r.sp = BRep_Tool::SameParameter(e) ? 1 : 0;
            r.sr = BRep_Tool::SameRange(e) ? 1 : 0;
            r.etol = BRep_Tool::Tolerance(e);
            r.flagged = flagged;
            r.e = e;
            rows.push_back(r);
            if (flagged) {
                nFlagged++;
                char key[128];
                std::snprintf(key, sizeof(key), "%s,%s,%s", cls, r.c3, r.pk);
                hist[key]++;
                std::fprintf(stderr,
                             "DIAG_PCINVALID rid=%d ie=%d class=%s ci=%d partner=%d curve3d=%s "
                             "pcurveKind=%s pcurveBirthSite=%s nPcurvesOnThisFace=%d "
                             "maxDev3d=%.6f sameParam=%d sameRange=%d edgeTol=%.6f\n",
                             rid, r.ie, cls, ci, partner, r.c3, r.pk, r.birth, r.nPc, r.dev, r.sp,
                             r.sr, r.etol);
                double exactPlane = -1, exactOther = -1, brepPlane = -1, brepOther = -1;
                edgeExactPlaneAndOther(f, e, exactPlane, exactOther, brepPlane, brepOther);
                const Region* selfR = regionById(rs, rid);
                const double cap = selfR ? meshTolCap(mv, selfR) : -1.0;
                std::fprintf(stderr,
                             "DIAG_FIRSTKEEP_FLAG rid=%d ie=%d class=%s lastWriter=%s "
                             "tol=%.9g exactPlane=%.9g exactOther=%.9g brepMaxPlane=%.9g "
                             "brepMaxOther=%.9g cap=%.9g\n",
                             rid, r.ie, cls, lastTolWriterOf(e), r.etol, exactPlane, exactOther,
                             brepPlane, brepOther, cap);
            }
        }
        std::fprintf(stderr, "DIAG_PCINVALID_COV rid=%d flagged=%d/%d examined=%d present=%d\n",
                     rid, nFlagged, emap.Extent(), emap.Extent(), emap.Extent());
        for (const auto& kv : hist)
            std::fprintf(stderr, "DIAG_PCINVALID_HIST rid=%d n=%d key=%s\n", rid, kv.second,
                         kv.first.c_str());

        std::unordered_map<std::string, int> ctrlN;
        for (const Row& r : rows) {
            if (r.flagged) continue;
            int& k = ctrlN[r.cls];
            if (k >= 5) continue;
            k++;
            std::fprintf(stderr,
                         "DIAG_PCINVALID_CTRL rid=%d ie=%d class=%s ci=%d partner=%d curve3d=%s "
                         "pcurveKind=%s pcurveBirthSite=%s nPcurvesOnThisFace=%d "
                         "maxDev3d=%.6f sameParam=%d sameRange=%d edgeTol=%.6f\n",
                         rid, r.ie, r.cls, r.ci, r.partner, r.c3, r.pk, r.birth, r.nPc, r.dev, r.sp,
                         r.sr, r.etol);
        }

        noteHpPlateEdgeTols(f, rid);
        try {
            TopTools_IndexedMapOfShape emAll;
            TopExp::MapShapes(f, TopAbs_EDGE, emAll);
            for (int i = 1; i <= emAll.Extent(); i++) {
                const TopoDS_Edge ee = TopoDS::Edge(emAll(i));
                const void* ts = diagTShapePtr(ee);
                if (ts) gAllPlateEdgeTol[ts] = BRep_Tool::Tolerance(ee);
            }
        } catch (const Standard_Failure&) {
        }
        gp_Pln platePln;
        const bool havePln = facePlanePln(f, platePln);
        const double plateTol = std::max(BRep_Tool::Tolerance(f), mv.sewTol);
        const Region* selfR = regionById(rs, rid);

        auto emitTol = [&](const Row& r) {
            const double exact = exactMaxOnPlane(f, r.e, r.c3);
            const double brep = brepCheckInContextMax(f, r.e);
            gTolCheckRows++;
            gTolCheckPlaneRows++;
            if (exact >= 0.0 && brep >= 0.0) {
                if (brep <= exact + 1e-12) {
                    gTolCheckBrepLeExact++;
                    gTolCheckPlaneBrepLe++;
                } else
                    gBindTolBrepGtExact++;
            }
            TopoDS_Vertex va, vb;
            TopExp::Vertices(r.e, va, vb, Standard_True);
            const double vt0 = va.IsNull() ? -1.0 : BRep_Tool::Tolerance(va);
            const double vt1 = vb.IsNull() ? -1.0 : BRep_Tool::Tolerance(vb);
            const double vo0 =
                (!havePln || va.IsNull()) ? -1.0 : std::fabs(planeSignedDist(platePln, BRep_Tool::Pnt(va)));
            const double vo1 =
                (!havePln || vb.IsNull()) ? -1.0 : std::fabs(planeSignedDist(platePln, BRep_Tool::Pnt(vb)));
            const char* verd = (exact > r.etol) ? "exact>tol" : "exact<=tol";
            const char* brepVerd = (exact >= 0.0 && brep >= 0.0 && brep > exact + 1e-12) ? "brep>exact" : "brep<=exact";
            std::fprintf(stderr,
                         "DIAG_TOLCHECK rid=%d ie=%d class=%s curve3d=%s flagged=%d exactMax=%.9f "
                         "edgeTol=%.9f brepMax=%.9f vtxTol0=%.9f vtxOff0=%.9f vtxTol1=%.9f "
                         "vtxOff1=%.9f verdict=%s brepVerdict=%s\n",
                         rid, r.ie, r.cls, r.c3, r.flagged ? 1 : 0, exact, r.etol, brep, vt0, vo0,
                         vt1, vo1, verd, brepVerd);
            if (std::strcmp(r.cls, "seam") == 0) {
                const int cylRid = seamCylinderRid(r.e, rs, geom);
                const Region* cr = regionById(rs, cylRid);
                double angDeg = -1.0;
                int isCirc = (std::strcmp(r.c3, "Circle") == 0) ? 1 : 0;
                int coplanar = 0;
                if (havePln && cr && cr->type == SurfType::Cylinder) {
                    const double ang = cr->ax.Direction().Angle(platePln.Axis().Direction());
                    angDeg = std::min(ang, kPi - ang) * 180.0 / kPi;
                }
                if (isCirc && havePln) {
                    try {
                        Standard_Real f3 = 0, l3 = 0;
                        Handle(Geom_Curve) c3 = BRep_Tool::Curve(r.e, f3, l3);
                        Handle(Geom_Circle) gc = Handle(Geom_Circle)::DownCast(basisCurveOf(c3));
                        if (!gc.IsNull()) {
                            const gp_Circ circ = gc->Circ();
                            const double off = std::fabs(planeSignedDist(platePln, circ.Location()));
                            const double ang =
                                circ.Axis().Direction().Angle(platePln.Axis().Direction());
                            const double tilt = circ.Radius() * std::sin(ang);
                            const double plat = (selfR ? meshTolCap(mv, selfR) : plateTol);
                            if (off <= plat && tilt <= plat) coplanar = 1;
                            if (angDeg < 0.0)
                                angDeg = std::min(ang, kPi - ang) * 180.0 / kPi;
                        }
                    } catch (const Standard_Failure&) {
                    }
                }
                std::fprintf(stderr,
                             "DIAG_SEAMGEOM rid=%d ie=%d cylRid=%d axisNormalAngleDeg=%.4f "
                             "curveIsCircle=%d circlePlaneCoplanarWithPlate=%d exactMax=%.9f\n",
                             rid, r.ie, cylRid, angDeg, isCirc, coplanar, exact);
            }
        };
        std::unordered_map<std::string, int> ctrlN2;
        for (const Row& r : rows) emitTol(r);

        gIncontextPlateExamined++;
        if (nFlagged > 0) gIncontextPlateFlagged++;
        std::fprintf(stderr, "DIAG_INCONTEXT_COV rid=%d flagged=%d/%d examined=%d\n", rid,
                     nFlagged, emap.Extent(), emap.Extent());

        if (rid == 3 && !gHpPlateEdgeTol.empty()) {
            std::vector<double> tols;
            tols.reserve(gHpPlateEdgeTol.size());
            for (const auto& kv : gHpPlateEdgeTol) tols.push_back(kv.second);
            std::sort(tols.begin(), tols.end());
            int nAbove1um = 0;
            for (double t : tols)
                if (t > 0.001) nAbove1um++;
            const double med = tols[tols.size() / 2];
            std::fprintf(stderr,
                         "DIAG_TOLCHECK_DIST n=%zu min=%.9f median=%.9f max=%.9f nAbove1um=%d\n",
                         tols.size(), tols.front(), med, tols.back(), nAbove1um);
        }

        if (rid == 2) {
            TopoDS_Wire outerW;
            int maxNE = -1;
            for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next()) {
                const TopoDS_Wire w = TopoDS::Wire(wx.Current());
                int nE = 0;
                for (TopoDS_Iterator it(w); it.More(); it.Next())
                    if (it.Value().ShapeType() == TopAbs_EDGE) nE++;
                if (nE > maxNE) {
                    maxNE = nE;
                    outerW = w;
                }
            }
            if (outerW.IsNull() && !ow.IsNull()) outerW = ow;
            const char *oOri = "?", *rOri = "?", *oC2 = "?", *rC2 = "?";
            double oGap = -1.0, rGap = -1.0;
            try {
                BRepCheck_Wire chk(outerW);
                chk.GeometricControls(Standard_True);
                oC2 = brepCheckName((int)chk.Closed2d(f, Standard_False));
                oOri = brepCheckName((int)chk.Orientation(f, Standard_False));
                oGap = uvClosureGapOnFace(f, outerW);
            } catch (const Standard_Failure&) {
            }
            try {
                const TopoDS_Wire wr = TopoDS::Wire(outerW.Reversed());
                BRepCheck_Wire chk(wr);
                chk.GeometricControls(Standard_True);
                rC2 = brepCheckName((int)chk.Closed2d(f, Standard_False));
                rOri = brepCheckName((int)chk.Orientation(f, Standard_False));
                rGap = uvClosureGapOnFace(f, wr);
            } catch (const Standard_Failure&) {
            }
            const int clears =
                (std::strcmp(rOri, "NoError") == 0 && std::strcmp(oOri, "NoError") != 0) ? 1 : 0;
            const int gapRemains = (rGap > 1.0) ? 1 : 0;
            std::fprintf(stderr,
                         "DIAG_PCINVALID_WIREREV rid=2 origOrient=%s revOrient=%s origClosed2d=%s "
                         "revClosed2d=%s origGap=%.5f revGap=%.5f clearsSt32=%d gapRemains=%d\n",
                         oOri, rOri, oC2, rC2, oGap, rGap, clears, gapRemains);
        }
    } catch (const Standard_Failure&) {
    }
    (void)ow;
}


// --- OP2p: DIAG_WIREORI2 at first keep-try instant (D-S3-115) ---
std::vector<TopoDS_Edge> wireEdgesIterator2(const TopoDS_Wire& w) {
    std::vector<TopoDS_Edge> out;
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        if (it.Value().ShapeType() == TopAbs_EDGE) out.push_back(TopoDS::Edge(it.Value()));
    }
    return out;
}

std::vector<int> tshapeCycleFromEdge02(const std::vector<TopoDS_Edge>& elist) {
    const int n = (int)elist.size();
    if (n == 0) return {};
    std::unordered_map<const void*, std::vector<std::pair<int, int>>> adj;
    for (int i = 0; i < n; i++) {
        TopoDS_Vertex t, h;
        TopExp::Vertices(elist[(size_t)i], t, h, Standard_True);
        const void* tT = diagTShapePtr(t);
        const void* hT = diagTShapePtr(h);
        if (tT) adj[tT].push_back({i, 1});
        if (hT) adj[hT].push_back({i, 0});
    }
    std::vector<char> used((size_t)n, 0);
    std::vector<int> order;
    order.push_back(0);
    used[0] = 1;
    TopoDS_Vertex t0, h0;
    TopExp::Vertices(elist[0], t0, h0, Standard_True);
    const void* curVtx = diagTShapePtr(h0);
    while ((int)order.size() < n && curVtx) {
        bool found = false;
        auto it = adj.find(curVtx);
        if (it != adj.end()) {
            for (const auto& pr : it->second) {
                if (used[(size_t)pr.first]) continue;
                order.push_back(pr.first);
                used[(size_t)pr.first] = 1;
                TopoDS_Vertex ta, hb;
                TopExp::Vertices(elist[(size_t)pr.first], ta, hb, Standard_True);
                curVtx = pr.second ? diagTShapePtr(hb) : diagTShapePtr(ta);
                found = true;
                break;
            }
        }
        if (!found) break;
    }
    return order;
}

const char* vtxOriChar(TopAbs_Orientation o) {
    if (o == TopAbs_FORWARD) return "F";
    if (o == TopAbs_REVERSED) return "R";
    return "?";
}

std::string analyzerWireStatusList(const BRepCheck_Analyzer& an, const TopoDS_Face& face,
                                   const TopoDS_Wire& w) {
    try {
        Handle(BRepCheck_Result) res = an.Result(w);
        if (res.IsNull()) return "not-determined";
        return formatBRepStatusListAll(res->StatusOnShape(face));
    } catch (const Standard_Failure&) {
        return "not-determined";
    }
}

void dumpDiagWireOri2Edges(int rid, char tryAB, const char* wireLabel, const TopoDS_Face& face,
                           const TopoDS_Wire& w, const RegionSet& rs,
                           const std::vector<ChainGeom>& geom) {
    const std::vector<TopoDS_Edge> elist = wireEdgesIterator2(w);
    const int nE = (int)elist.size();
    const std::vector<int> order = tshapeCycleFromEdge02(elist);
    const void* prevHeadT = nullptr;
    for (size_t ki = 0; ki < order.size(); ki++) {
        const int ie = order[ki];
        if (ie < 0 || ie >= nE) continue;
        const TopoDS_Edge& e = elist[(size_t)ie];
        TopoDS_Vertex tail, head;
        TopExp::Vertices(e, tail, head, Standard_True);
        const void* tailT = diagTShapePtr(tail);
        const void* headT = diagTShapePtr(head);
        const char* junction = "ok";
        if (ki > 0 && prevHeadT && tailT != prevHeadT) {
            if (prevHeadT == headT) junction = "head-head";
            else if (tailT == prevHeadT) junction = "tail-tail";
            else junction = "mismatch";
        }
        int ci = -1;
        const void* eT = diagTShapePtr(e);
        if (eT && gWirePopCi.count(eT)) ci = gWirePopCi[eT];
        const char* site = (eT && gWirePopSite.count(eT)) ? gWirePopSite[eT] : "not-determined";
        int hasPc = 0, sameParam = BRep_Tool::SameParameter(e) ? 1 : 0;
        try {
            Standard_Real a = 0, b = 0;
            Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(e, face, a, b);
            hasPc = pc.IsNull() ? 0 : 1;
        } catch (const Standard_Failure&) {
            hasPc = -1;
        }
        const int nativeDir = (ki > 0 && prevHeadT && tailT != prevHeadT) ? -1 : +1;
        std::fprintf(stderr,
                     "DIAG_WIREORI2_E rid=%d try=%c wire=%s ie=%d ci=%d site=%s tailT=%p headT=%p "
                     "prevHeadT=%p junction=%s vtxOriInEdge=%s,%s pcurveOnFace=%d sameParam=%d "
                     "nativeDir=%+d\n",
                     rid, tryAB, wireLabel, ie, ci, site, tailT, headT, prevHeadT, junction,
                     vtxOriChar(tail.Orientation()), vtxOriChar(head.Orientation()), hasPc,
                     sameParam, nativeDir);
        prevHeadT = headT;
    }
}

void dumpDiagWireOri2StopCause(int rid, char tryAB, const char* wireLabel, const TopoDS_Face& face,
                               const TopoDS_Wire& w, int stopIe) {
    if (stopIe < 0) return;
    const std::vector<TopoDS_Edge> elist = wireEdgesIterator2(w);
    if (stopIe >= (int)elist.size()) return;
    const TopoDS_Edge& stopE = elist[(size_t)stopIe];
    TopoDS_Vertex tail, head;
    TopExp::Vertices(stopE, tail, head, Standard_True);
    const TopoDS_Vertex stopVtx = head;
    const void* stopT = diagTShapePtr(stopVtx);
    const TopoDS_Edge& e0 = elist[0];
    std::fprintf(stderr, "DIAG_WIREORI2_STOP rid=%d try=%c wire=%s stopIe=%d stopVtxT=%p\n", rid,
                 tryAB, wireLabel, stopIe, stopT);
    TopTools_IndexedDataMapOfShapeListOfShape vmap;
    TopExp::MapShapesAndAncestors(face, TopAbs_VERTEX, TopAbs_EDGE, vmap);
    try {
        for (int i = 1; i <= vmap.FindIndex(stopVtx); i++) {
            if (!vmap.FindFromIndex(i).Contains(stopVtx)) continue;
        }
        const TopTools_ListOfShape& adj = vmap.FindFromKey(stopVtx);
        for (TopTools_ListOfShape::Iterator it(adj); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            const TopoDS_Edge ae = TopoDS::Edge(it.Value());
            TopoDS_Vertex ta, hb;
            TopExp::Vertices(ae, ta, hb, Standard_True);
            int hasPc = 0;
            try {
                Standard_Real a = 0, b = 0;
                hasPc = BRep_Tool::CurveOnSurface(ae, face, a, b).IsNull() ? 0 : 1;
            } catch (const Standard_Failure&) {
                hasPc = -1;
            }
            std::fprintf(stderr,
                         "DIAG_WIREORI2_STOPADJ rid=%d try=%c edgeT=%p oriTail=%s oriHead=%s "
                         "pcurveOnFace=%d isSameE0=%d degenerated=%d\n",
                         rid, tryAB, diagTShapePtr(ae), vtxOriChar(ta.Orientation()),
                         vtxOriChar(hb.Orientation()), hasPc, ae.IsSame(e0) ? 1 : 0,
                         BRep_Tool::Degenerated(ae) ? 1 : 0);
        }
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_WIREORI2_STOPADJ rid=%d try=%c adj=not-determined\n", rid,
                     tryAB);
    }
}

void dumpDiagWireOri2One(int rid, char tryAB, const TopoDS_Face& face, const TopoDS_Wire& w,
                         const char* wireLabel, const RegionSet& rs,
                         const std::vector<ChainGeom>& geom, const WireOriTally* tally) {
    if (!diagP2Enabled() || w.IsNull() || face.IsNull()) return;
    const std::vector<TopoDS_Edge> elist = wireEdgesIterator2(w);
    const int nE = (int)elist.size();
    int weVisited = 0, stopIe = -1;
    const char* closed2d = "not-determined";
    const char* orientReal = "not-determined";
    std::string anaWire = "not-determined";
    std::string anaFace = "not-determined";
    const int revBy = tally ? tally->reversedByRule : -1;
    const int notDet = tally ? tally->notDetermined : -1;
    const int seedRev = tally ? tally->seedReversed : -1;
    const int flippedLp = tally ? tally->flipped : -1;
    const char* seedBy = tally && tally->seedOrientedBy ? tally->seedOrientedBy : "none";
    try {
        BRepCheck_Analyzer an(face, Standard_True);
        anaWire = analyzerWireStatusList(an, face, w);
        try {
            Handle(BRepCheck_Result) fres = an.Result(face);
            if (!fres.IsNull()) anaFace = formatBRepStatusListAll(fres->StatusOnShape(face));
        } catch (const Standard_Failure&) {
            anaFace = "not-determined";
        }
        BRepCheck_Wire chk(w);
        chk.GeometricControls(Standard_True);
        closed2d = brepCheckName((int)chk.Closed2d(face, Standard_False));
        orientReal = brepCheckName((int)chk.Orientation(face, Standard_False));
    } catch (const Standard_Failure&) {
    }
    try {
        BRepTools_WireExplorer ex(w, face);
        while (ex.More()) {
            weVisited++;
            ex.Next();
        }
        if (weVisited < nE) stopIe = weVisited - 1;
    } catch (const Standard_Failure&) {
        stopIe = -2;
    }
    const char* stopStr = (stopIe >= 0) ? "ie" : "none";
    char stopBuf[16];
    if (stopIe >= 0) {
        std::snprintf(stopBuf, sizeof(stopBuf), "%d", stopIe);
        stopStr = stopBuf;
    }
    std::fprintf(stderr,
                 "DIAG_WIREORI2 rid=%d try=%c wire=%s nE=%d weVisited=%d/%d stopIe=%s "
                 "closed2dReal=%s orientationReal=%s anaWire=%s anaFace=%s "
                 "reversedByRule=%d notDetermined=%d seedOrientedBy=%s seedReversed=%d "
                 "flipped=%d\n",
                 rid, tryAB, wireLabel, nE, weVisited, nE, stopStr, closed2d, orientReal,
                 anaWire.c_str(), anaFace.c_str(), revBy, notDet, seedBy, seedRev, flippedLp);
    dumpDiagWireOri2Edges(rid, tryAB, wireLabel, face, w, rs, geom);
    if (stopIe >= 0) dumpDiagWireOri2StopCause(rid, tryAB, wireLabel, face, w, stopIe);
}

void dumpDiagWireOri2AtKeepTry(int rid, const TopoDS_Face& face, const TopoDS_Wire& ow,
                               const std::vector<TopoDS_Wire>& innerW, const RegionSet& rs,
                               const std::vector<ChainGeom>& geom) {
    if (!diagP2Enabled() || rid < 0 || face.IsNull()) return;
    static thread_local std::unordered_set<int> printed;
    if (!printed.insert(rid).second) return;
    auto tallyAt = [](size_t i) -> const WireOriTally* {
        if (i >= gPlateWireOri.size()) return nullptr;
        return &gPlateWireOri[i];
    };
    dumpDiagWireOri2One(rid, 'A', face, ow, "OUT", rs, geom, tallyAt(0));
    for (size_t ii = 0; ii < innerW.size(); ii++) {
        char lbl[16];
        std::snprintf(lbl, sizeof(lbl), "IN%zu", ii);
        dumpDiagWireOri2One(rid, 'A', face, innerW[ii], lbl, rs, geom, tallyAt(1 + ii));
    }
    if (rid == 2) {
        try {
            const TopoDS_Wire owB = TopoDS::Wire(ow.Reversed());
            dumpDiagWireOri2One(rid, 'B', face, owB, "OUT", rs, geom, tallyAt(0));
        } catch (const Standard_Failure&) {
        }
    }
}

// D-130 measurement: per-wire sense of a plate face at the instant makeFaceKeep is
// chosen. One line per wire (outer first, then each inner in build order) plus a
// face-level line: BRepCheck_Face needs exactly one wire whose own face classifies
// the infinite point OUT. Diag-only; touches no live shape (probe faces are fresh).
void dumpDiagPlateWireSense(const Region& r, const Handle(Geom_Plane)& gpl, const TopoDS_Wire& ow,
                            const std::vector<TopoDS_Wire>& innerW,
                            const std::vector<char>& innerWireRev, const char* stage) {
    if (!diagP2Enabled() || gpl.IsNull()) return;
    int nOuterLike = 0;
    auto oneWire = [&](const TopoDS_Wire& w, const char* role, int iw, int wireRev) {
        if (w.IsNull()) return;
        int nE = 0;
        TopoDS_Edge only;
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            if (it.Value().ShapeType() != TopAbs_EDGE) continue;
            if (nE == 0) only = TopoDS::Edge(it.Value());
            nE++;
        }
        const int single = (nE == 1 && !only.IsNull() && edgeSpansFullCircle(only)) ? 1 : 0;
        const char* eOri = only.IsNull() ? "-" : topAbsOriName(only.Orientation());
        const double ca = wireCurveSignedAreaUV(w, r);
        const double pa = wireSignedAreaUV(w, r);
        const char* soloInf = "?";
        try {
            BRep_Builder B;
            TopoDS_Face probe;
            B.MakeFace(probe, gpl, Precision::Confusion());
            B.Add(probe, w);
            BRepTopAdaptor_FClass2d cl(probe, Precision::PConfusion());
            const TopAbs_State st = cl.PerformInfinitePoint();
            soloInf = topAbsStateName(st);
            if (st == TopAbs_OUT) nOuterLike++;
        } catch (const Standard_Failure&) {
        }
        std::fprintf(stderr,
                     "DIAG_OP2L_WIRE rid=%d stage=%s role=%s iw=%d nE=%d closed=%d wireOri=%s "
                     "edgeOri=%s singleClosedCircle=%d curveArea=%.6f polyArea=%.6f "
                     "wireRevApplied=%d soloInfPt=%s\n",
                     r.id, stage ? stage : "?", role, iw, nE, w.Closed() ? 1 : 0,
                     topAbsOriName(w.Orientation()), eOri, single, ca, pa, wireRev, soloInf);
    };
    oneWire(ow, "outer", -1, 0);
    for (size_t i = 0; i < innerW.size(); i++)
        oneWire(innerW[i], "inner", (int)i,
                i < innerWireRev.size() ? (int)innerWireRev[i] : 0);
    const char* oriW = "?";
    const char* clsW = "?";
    try {
        BRep_Builder B;
        TopoDS_Face probe;
        B.MakeFace(probe, gpl, Precision::Confusion());
        B.Add(probe, ow);
        for (const auto& iw : innerW)
            if (!iw.IsNull()) B.Add(probe, iw);
        BRepCheck_Face fc(probe);
        fc.GeometricControls(Standard_True);
        oriW = brepCheckName((int)fc.OrientationOfWires(Standard_False));
        clsW = brepCheckName((int)fc.ClassifyWires(Standard_False));
    } catch (const Standard_Failure&) {
    }
    // Subset probe: which wire combination is the one BRepCheck refuses.
    {
        auto probeSubset = [&](const std::vector<int>& useInner, const char* label) {
            const char* verdict = "?";
            std::string own;
            try {
                BRep_Builder B;
                TopoDS_Face probe;
                B.MakeFace(probe, gpl, Precision::Confusion());
                B.Add(probe, ow);
                for (int k : useInner)
                    if (k >= 0 && (size_t)k < innerW.size() && !innerW[(size_t)k].IsNull())
                        B.Add(probe, innerW[(size_t)k]);
                BRepCheck_Analyzer an(probe);
                verdict = an.IsValid() ? "valid" : "invalid";
                const Handle(BRepCheck_Result) res = analyzerResultOf(an, probe);
                own = analyzerStatusOwn(res);
                BRepCheck_Face fc(probe);
                fc.GeometricControls(Standard_True);
                const BRepCheck_Status sInt = fc.IntersectWires(Standard_False);
                const BRepCheck_Status sCls = fc.ClassifyWires(Standard_False);
                const BRepCheck_Status sOri = fc.OrientationOfWires(Standard_False);
                BRepCheck_Wire wc(ow);
                wc.GeometricControls(Standard_True);
                TopoDS_Edge xe1, xe2;
                const BRepCheck_Status wSelf = wc.SelfIntersect(probe, xe1, xe2, Standard_False);
                const BRepCheck_Status wCl2 = wc.Closed2d(probe, Standard_False);
                const BRepCheck_Status wOri = wc.Orientation(probe, Standard_False);
                std::fprintf(stderr,
                             "DIAG_OP2L_SUBDET rid=%d subset=%s intersect=%s classify=%s "
                             "orient=%s wSelfX=%s wClosed2d=%s wOrient=%s\n",
                             r.id, label, brepCheckName((int)sInt), brepCheckName((int)sCls),
                             brepCheckName((int)sOri), brepCheckName((int)wSelf),
                             brepCheckName((int)wCl2), brepCheckName((int)wOri));
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr, "DIAG_OP2L_SUBSET rid=%d subset=%s verdict=%s own=%s\n", r.id,
                         label, verdict, own.empty() ? "-" : own.c_str());
        };
        probeSubset({}, "outer");
        for (size_t k = 0; k < innerW.size(); k++) {
            char lbl[32];
            std::snprintf(lbl, sizeof(lbl), "outer+inner%zu", k);
            probeSubset({(int)k}, lbl);
        }
    }
    // D-130 step 3 localiser: is the refusal the wire sense, or the pcurve
    // residual of an off-plane mesh vertex? Run the same analyzer on a COPY of
    // the outer-only probe whose tolerances have been recomputed from geometry
    // (BRepLib::UpdateTolerances). Copies never touch a live TShape.
    {
        const char* ownRaw = "?";
        const char* ownUpd = "?";
        std::string sRaw, sUpd;
        int unorientableRaw = -1, unorientableUpd = -1;
        double maxEdgeTolUpd = 0.0;
        try {
            BRep_Builder B;
            TopoDS_Face probe;
            B.MakeFace(probe, gpl, Precision::Confusion());
            B.Add(probe, ow);
            {
                BRepCheck_Analyzer an(probe);
                ownRaw = an.IsValid() ? "valid" : "invalid";
                sRaw = analyzerStatusOwn(analyzerResultOf(an, probe));
                Handle(BRepCheck_Face) fr =
                    Handle(BRepCheck_Face)::DownCast(analyzerResultOf(an, probe));
                if (!fr.IsNull()) unorientableRaw = fr->IsUnorientable() ? 1 : 0;
            }
            BRepBuilderAPI_Copy cp(probe);
            TopoDS_Shape cs = cp.Shape();
            if (!cs.IsNull()) {
                BRepLib::UpdateTolerances(cs, Standard_True);
                for (TopExp_Explorer ex(cs, TopAbs_EDGE); ex.More(); ex.Next())
                    maxEdgeTolUpd =
                        std::max(maxEdgeTolUpd, BRep_Tool::Tolerance(TopoDS::Edge(ex.Current())));
                BRepCheck_Analyzer an2(cs);
                ownUpd = an2.IsValid() ? "valid" : "invalid";
                sUpd = analyzerStatusOwn(analyzerResultOf(an2, cs));
                Handle(BRepCheck_Face) fr2 =
                    Handle(BRepCheck_Face)::DownCast(analyzerResultOf(an2, cs));
                if (!fr2.IsNull()) unorientableUpd = fr2->IsUnorientable() ? 1 : 0;
            }
        } catch (const Standard_Failure&) {
        }
        std::fprintf(stderr,
                     "DIAG_OP2L_WHY rid=%d outerOnly raw=%s rawOwn=%s rawUnorientable=%d "
                     "afterTolUpdate=%s updOwn=%s updUnorientable=%d maxEdgeTolUpd=%.6f\n",
                     r.id, ownRaw, sRaw.empty() ? "-" : sRaw.c_str(), unorientableRaw, ownUpd,
                     sUpd.empty() ? "-" : sUpd.c_str(), unorientableUpd, maxEdgeTolUpd);
    // Which wire form does BRepCheck accept? as-built / reversed / MakeWire-rebuilt
        // / rebuilt+reversed, each also with the face flag reversed. Diag-only.
        auto variant = [&](const char* label, const TopoDS_Wire& w, bool revFace) {
            if (w.IsNull()) return;
            const char* v = "?";
            std::string own;
            int uno = -1;
            const char* inf = "?";
            try {
                BRep_Builder B;
                TopoDS_Face probe;
                B.MakeFace(probe, gpl, Precision::Confusion());
                B.Add(probe, w);
                if (revFace) probe.Reverse();
                BRepTopAdaptor_FClass2d cl(probe, Precision::PConfusion());
                inf = topAbsStateName(cl.PerformInfinitePoint());
                BRepCheck_Analyzer an(probe);
                v = an.IsValid() ? "valid" : "invalid";
                own = analyzerStatusOwn(analyzerResultOf(an, probe));
                Handle(BRepCheck_Face) fr =
                    Handle(BRepCheck_Face)::DownCast(analyzerResultOf(an, probe));
                if (!fr.IsNull()) uno = fr->IsUnorientable() ? 1 : 0;
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr,
                         "DIAG_OP2L_VAR rid=%d form=%s revFace=%d verdict=%s own=%s "
                         "unorientable=%d infPt=%s\n",
                         r.id, label, revFace ? 1 : 0, v, own.empty() ? "-" : own.c_str(), uno,
                         inf);
        };
        TopoDS_Wire rebuilt;
        try {
            BRepBuilderAPI_MakeWire mw;
            for (TopoDS_Iterator it(ow); it.More(); it.Next())
                if (it.Value().ShapeType() == TopAbs_EDGE) mw.Add(TopoDS::Edge(it.Value()));
            if (mw.IsDone()) {
                rebuilt = mw.Wire();
                rebuilt.Closed(Standard_True);
            }
        } catch (const Standard_Failure&) {
        }
        variant("asbuilt", ow, false);
        variant("asbuilt", ow, true);
        // Copy + BRepLib::SameParameter: does making every pcurve consistent with
        // its 3d curve (tolerance recomputed on a copy) clear the refusal?
        try {
            BRep_Builder B;
            TopoDS_Face probe;
            B.MakeFace(probe, gpl, Precision::Confusion());
            B.Add(probe, ow);
            BRepBuilderAPI_Copy cp(probe);
            TopoDS_Shape cs = cp.Shape();
            if (!cs.IsNull()) {
                BRepLib::SameParameter(cs, Precision::Confusion(), Standard_True);
                BRepLib::UpdateTolerances(cs, Standard_True);
                double mt = 0.0;
                for (TopExp_Explorer ex(cs, TopAbs_EDGE); ex.More(); ex.Next())
                    mt = std::max(mt, BRep_Tool::Tolerance(TopoDS::Edge(ex.Current())));
                BRepCheck_Analyzer an(cs);
                const std::string own = analyzerStatusOwn(analyzerResultOf(an, cs));
                int uno = -1;
                Handle(BRepCheck_Face) fr =
                    Handle(BRepCheck_Face)::DownCast(analyzerResultOf(an, cs));
                if (!fr.IsNull()) uno = fr->IsUnorientable() ? 1 : 0;
                std::fprintf(stderr,
                             "DIAG_OP2L_VAR rid=%d form=sameparam revFace=0 verdict=%s own=%s "
                             "unorientable=%d infPt=- maxTol=%.6f\n",
                             r.id, an.IsValid() ? "valid" : "invalid",
                             own.empty() ? "-" : own.c_str(), uno, mt);
            }
        } catch (const Standard_Failure&) {
        }
        // Per-edge analyzer status on the outer-only probe (own + in this face).
        try {
            BRep_Builder B;
            TopoDS_Face probe;
            B.MakeFace(probe, gpl, Precision::Confusion());
            B.Add(probe, ow);
            BRepCheck_Analyzer an(probe);
            TopTools_IndexedMapOfShape emap;
            TopExp::MapShapes(probe, TopAbs_EDGE, emap);
            for (int ie = 1; ie <= emap.Extent(); ie++) {
                const TopoDS_Edge e = TopoDS::Edge(emap(ie));
                const Handle(BRepCheck_Result) res = analyzerResultOf(an, e);
                std::fprintf(stderr, "DIAG_OP2L_EDGE rid=%d ie=%d own=%s inFace=%s tol=%.7f\n",
                             r.id, ie - 1, analyzerStatusOwn(res).c_str(),
                             analyzerStatusOnFace(res, probe).c_str(), BRep_Tool::Tolerance(e));
            }
        } catch (const Standard_Failure&) {
        }
        variant("revwire", TopoDS::Wire(ow.Reversed()), false);
        variant("rebuilt", rebuilt, false);
        variant("rebuilt-rev", rebuilt.IsNull() ? rebuilt : TopoDS::Wire(rebuilt.Reversed()),
                false);
        // Per-edge UV chain on the probe face: pcurve presence on the FACE's own
        // surface handle, oriented UV endpoints, and the gap to the previous edge.
        // A 2d gap (or a missing pcurve) is what leaves BRepCheck unable to orient.
        try {
            BRep_Builder B;
            TopoDS_Face probe;
            B.MakeFace(probe, gpl, Precision::Confusion());
            B.Add(probe, ow);
            TopLoc_Location floc;
            Handle(Geom_Surface) fs = BRep_Tool::Surface(probe, floc);
            int nWE = 0, nNoPc = 0;
            double maxGap2d = 0.0;
            gp_Pnt2d prevEnd;
            bool havePrev = false;
            gp_Pnt2d firstStart;
            for (BRepTools_WireExplorer we(ow, probe); we.More(); we.Next()) {
                const TopoDS_Edge e = we.Current();
                nWE++;
                Standard_Real pf = 0, pl = 0;
                Handle(Geom2d_Curve) pc =
                    BRep_Tool::CurveOnSurface(e, fs, floc, pf, pl);
                if (pc.IsNull()) {
                    nNoPc++;
                    std::fprintf(stderr, "DIAG_OP2L_UV rid=%d ie=%d pc=none\n", r.id, nWE - 1);
                    havePrev = false;
                    continue;
                }
                gp_Pnt2d a = pc->Value(pf), b = pc->Value(pl);
                if (e.Orientation() == TopAbs_REVERSED) std::swap(a, b);
                double gap = havePrev ? prevEnd.Distance(a) : 0.0;
                if (havePrev) maxGap2d = std::max(maxGap2d, gap);
                if (nWE == 1) firstStart = a;
                std::fprintf(stderr,
                             "DIAG_OP2L_UV rid=%d ie=%d ori=%s uvA=(%.5f,%.5f) uvB=(%.5f,%.5f) "
                             "gapPrev=%.6f\n",
                             r.id, nWE - 1, topAbsOriName(e.Orientation()), a.X(), a.Y(), b.X(),
                             b.Y(), gap);
                prevEnd = b;
                havePrev = true;
            }
            const double closeGap = havePrev ? prevEnd.Distance(firstStart) : -1.0;
            int nEtot = 0;
            for (TopoDS_Iterator it(ow); it.More(); it.Next())
                if (it.Value().ShapeType() == TopAbs_EDGE) nEtot++;
            std::fprintf(stderr,
                         "DIAG_OP2L_UVSUM rid=%d nEdges=%d wireExplorerVisited=%d noPc=%d "
                         "maxGap2d=%.6f closeGap2d=%.6f\n",
                         r.id, nEtot, nWE, nNoPc, maxGap2d, closeGap);
        } catch (const Standard_Failure&) {
        }
    }
    const char* faceFwd = "?";
    const char* faceOut = "?";
    try {
        BRep_Builder B;
        TopoDS_Face probe;
        B.MakeFace(probe, gpl, Precision::Confusion());
        B.Add(probe, ow);
        for (const auto& iw : innerW)
            if (!iw.IsNull()) B.Add(probe, iw);
        BRepCheck_Analyzer anF(probe);
        faceFwd = anF.IsValid() ? "valid" : "invalid";
        TopoDS_Face probe2 = probe;
        setFaceOutward(probe2, r.outwardNormal);
        BRepCheck_Analyzer anO(probe2);
        faceOut = anO.IsValid() ? "valid" : "invalid";
    } catch (const Standard_Failure&) {
    }
    std::fprintf(stderr,
                 "DIAG_OP2L_FACE rid=%d stage=%s nWires=%zu nSoloOut=%d orientationOfWires=%s "
                 "classifyWires=%s outward=%d probeFwd=%s probeOutward=%s\n",
                 r.id, stage ? stage : "?", innerW.size() + 1, nOuterLike, oriW, clsW,
                 r.outwardNormal ? 1 : 0, faceFwd, faceOut);
}

// DECISION-130 tripwire. An analytic edge on a plate -- one born at
// bindAllVariants, i.e. enrolled in gSeamTShapes -- must name the chain it came
// from. An edge that cannot be resolved to a collapsed chain is an edge no other
// face can reference: the makeFaceKeep-seam / ci=-1 failure mode, which opens
// the shell silently.
//
// It REPORTS, it does not refuse. Refusing cost the two canonical P2 build
// fixtures whose expected builtAs is single/single/single --
// p2buildtest_fillet_strip_2tri and p2buildtest_quarter_round_1tri both explode
// region 1 to facets on an unnamedAnalyticEdges=1 verdict, with openE unchanged
// at 8 either way: the edge those two name IS shared, so the resolver, not the
// face, is what is wrong. Until the resolver is sound a refusal here would be
// the same silent substitution DECISION-130 forbids, only in the other
// direction. STL2STEP_130_TRIPWIRE_ABORT=1 is the hard stop for investigation.
int plateAnalyticEdgesWithoutChain(const Region& r, const TopoDS_Face& f,
                                   const std::vector<ChainGeom>& geom,
                                   TopoDS_Edge* firstBad = nullptr) {
    if (f.IsNull()) return 0;
    int n = 0;
    std::unordered_set<const void*> seen;
    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
        const TopoDS_Edge e = TopoDS::Edge(ex.Current());
        const void* ts = diagTShapePtr(e);
        if (!ts || !gSeamTShapes.count(ts) || !seen.insert(ts).second) continue;
        bool named = false;
        for (size_t ci = 0; ci < geom.size() && !named; ci++)
            for (const TopoDS_Edge& ge : geom[ci].edges)
                if (!ge.IsNull() && diagTShapePtr(ge) == ts) {
                    named = true;
                    break;
                }
        if (named) continue;
        if (n == 0 && firstBad) *firstBad = e;
        n++;
    }
    if (n > 0) {
        gp_Pnt p0, p1;
        if (firstBad && !firstBad->IsNull()) {
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(*firstBad, v0, v1, Standard_True);
            if (!v0.IsNull()) p0 = BRep_Tool::Pnt(v0);
            if (!v1.IsNull()) p1 = BRep_Tool::Pnt(v1);
        }
        std::fprintf(stderr,
                     "DIAG_130_TRIPWIRE rid=%d unnamedAnalyticEdges=%d verdict=reported "
                     "first=(%.4f,%.4f,%.4f)-(%.4f,%.4f,%.4f)\n",
                     r.id, n, p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z());
        if (const char* v = std::getenv("STL2STEP_130_TRIPWIRE_ABORT"); v && v[0] && v[0] != '0')
            std::abort();
    }
    return n;
}

bool buildPlanarFace(const Region& r, const RegionSet& rs, const MeshView& mv,
                     const std::vector<ChainGeom>& geom, const std::vector<char>& collapsed,
                     const std::vector<TopoDS_Edge>& meshE, const std::vector<char>& edgeOk,
                     TopoDS_Face& outF) {
    gPlateWireOri.clear();
    const Loop* outer = nullptr;
    std::vector<const Loop*> inners;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::Outer) outer = &lp;
        else if (lp.role == LoopRole::Inner) inners.push_back(&lp);
    }
    if (!outer) return false;
    TopoDS_Wire ow;
    if (!buildLoopWire(ow, *outer, rs, mv, geom, collapsed, meshE, edgeOk, &r)) return false;
    try {
        Handle(Geom_Plane) gpl = Handle(Geom_Plane)::DownCast(regionSurf(r, SurfVar::Plane));
        std::vector<TopoDS_Wire> innerW;
        std::vector<char> innerRev;
        bool innersOk = true;
        for (const Loop* ip : inners) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, *ip, rs, mv, geom, collapsed, meshE, edgeOk, &r)) {
                innersOk = false;
                break;
            }
            // D-130 step 1: a single closed edge realises its loop sense as that
            // EDGE's own orientation flag (buildLoopWire, D-S3-171 case (c)).
            // Reversing the WIRE on top of that double-negates the sense: the hole
            // then classifies its own infinite point OUT and BRepCheck_Face sees two
            // outer wires (BadOrientationOfSubshape / UnorientableShape). Multi-edge
            // loops keep the existing behaviour. Diag-only assertion: the wire-level
            // reversal must never be the thing that fixes a single-closed-edge loop.
            if (diagP2Enabled() && ip->chainIdx.size() == 1) {
                const int ci = ip->chainIdx[0];
                if (ci >= 0 && (size_t)ci < collapsed.size() && collapsed[(size_t)ci] &&
                    (size_t)ci < geom.size() && geom[(size_t)ci].collapsed &&
                    geom[(size_t)ci].edges.size() == 1 &&
                    edgeSpansFullCircle(geom[(size_t)ci].edges[0])) {
                    const gp_Dir z = r.ax.Direction();
                    const gp_Dir outw = r.outwardNormal ? z : gp_Dir(z.Reversed());
                    const int wantedInner = (z.Dot(outw) >= 0.0) ? -1 : 1;
                    const double ca = wireCurveSignedAreaUV(iw, r);
                    std::fprintf(stderr,
                                 "DIAG_OP2L_SINGLE rid=%d ci=%d wireOri=%s curveArea=%.6f "
                                 "wanted=%+d wireLevelReversalNeeded=%d\n",
                                 r.id, ci, topAbsOriName(iw.Orientation()), ca, wantedInner,
                                 ((ca > 0.0) != (wantedInner > 0)) ? 1 : 0);
                }
            }
            innerW.push_back(iw);
            innerRev.push_back(0);
        }
        if (innersOk) {
            int nBound = 0, nSkippedHasPc = 0, nSkippedNoCurve = 0;
            bindPlanePCurves(ow, gpl, r, 0.0, nBound, nSkippedHasPc, nSkippedNoCurve);
            for (auto& iw : innerW)
                bindPlanePCurves(iw, gpl, r, 0.0, nBound, nSkippedHasPc, nSkippedNoCurve);
            if (diagP2Enabled())
                dumpDiagWirePop(r, ow, rs, mv, geom, meshE);
            if (diagP2Enabled() && r.id >= 0 && r.id <= 3)
                dumpDiagHubWire(r, ow, gpl, rs, mv, geom, collapsed, meshE, innerW, r.outwardNormal);
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_PLANEPCBIND rid=%d nBound=%d nSkippedHasPc=%d nSkippedNoCurve=%d\n",
                             r.id, nBound, nSkippedHasPc, nSkippedNoCurve);
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_PLATETOL rid=%d sewTol=%.7f meshTolCap=%.7f maxVertexDev=%.7f "
                             "plateCap=%.7f\n",
                             r.id, mv.sewTol, meshTolCap(mv, &r), r.maxVertexDev,
                             plateFaceTolCap(mv, r));
            dumpDiagPlateWireSense(r, gpl, ow, innerW, innerRev, "keep-try");
        }
        if (innersOk && makeFaceKeep(gpl, ow, innerW, r.outwardNormal, outF)) {
            if (r.id >= 0 && (diagPlatesEnabled() || diagP2Enabled())) {
                dumpDiagAnalyzerBothTries(r.id, outF, ow);
                dumpDiagPcInvalid(r.id, outF, ow, rs, mv, geom, collapsed, meshE);
            }
            if (diagP2Enabled())
                dumpDiagWireOri2AtKeepTry(r.id, outF, ow, innerW, rs, geom);
            if (faceIsValid(outF) || ensureFaceValid(outF, plateFaceTolCap(mv, r))) {
                if (r.id >= 0 && (diagPlatesEnabled() || diagP2Enabled())) {
                    dumpPcurveFrames(r, ow);
                    dumpDiagPlateLine(r.id, outF);
                }
                TopoDS_Edge tripBad;
                (void)plateAnalyticEdgesWithoutChain(r, outF, geom, &tripBad);
                return !outF.IsNull();
            }
            diagPlateMakeFaceFail(r.id, outF, ow, rs, mv, geom, collapsed, meshE, mv.sewTol);
            if (r.id >= 0 && r.id <= 3) dumpPcurveFrames(r, ow);
        } else if (innersOk) {
            diagPlateMakeFaceFail(r.id, outF, ow, rs, mv, geom, collapsed, meshE, mv.sewTol);
            if (r.id >= 0 && r.id <= 3) dumpPcurveFrames(r, ow);
        }
        if (gpl.IsNull()) return false;
        // Remaining rungs retry makeFaceKeep on reconnect/reverse/flip wires.
        // Never BRepBuilderAPI_MakeFace(asPlane(r), ...): that copies the
        // Geom_Plane Handle and the bound TShapes (fill-boundary free edges).
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
                    TopoDS_Face f2;
                    if (makeFaceKeep(gpl, ow2, innerW, r.outwardNormal, f2) && faceIsValid(f2)) {
                        gReconnectFires++;
                        std::fprintf(stderr, "DIAG_RECONNECT fire=1 rid=%d nTri=%d\n", r.id,
                                     (int)mv.nTri);
                        outF = f2;
                        ow = ow2;
                        gReconnectFire++;
                    }
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF)) {
            // Outer-wire winding vs plane normal: retry reversed wire.
            try {
                TopoDS_Face fR;
                if (makeFaceKeep(gpl, TopoDS::Wire(ow.Reversed()), innerW, r.outwardNormal, fR)) {
                    if (r.id >= 0 && r.id <= 3)
                        diagPlateMakeFaceFail(r.id, fR, TopoDS::Wire(ow.Reversed()), rs, mv, geom,
                                             collapsed, meshE, mv.sewTol);
                    if (faceIsValid(fR)) outF = fR;
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF) && !innerW.empty()) {
            // Retry with inner wires flipped — imbrication/orientation of holes.
            try {
                std::vector<TopoDS_Wire> flipped;
                flipped.reserve(innerW.size());
                for (const auto& iw : innerW) flipped.push_back(TopoDS::Wire(iw.Reversed()));
                TopoDS_Face f2;
                if (makeFaceKeep(gpl, ow, flipped, r.outwardNormal, f2)) {
                    if (!r.outwardNormal) f2.Reverse();
                    if (faceIsValid(f2)) outF = f2;
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!faceIsValid(outF) && inners.empty()) {
            try {
                if (diagP2Enabled())
                    std::fprintf(stderr, "DIAG_130_SFF rid=%d rung=shapefix-face\n", r.id);
                ShapeFix_Face sff(outF);
                sff.FixOrientationMode() = 1;
                sff.FixAddNaturalBoundMode() = 0;
                sff.FixMissingSeamMode() = 0;
                // D-130: this rung exists to fix the WIRE ORIENTATION of a
                // single-wire plate. ShapeFix_Wire (FixWireMode, on by default)
                // closes 2d/3d gaps by MOVING vertices -- and these vertices are
                // shared with every neighbouring face, so the repair opens the
                // shell it was called to close. Measured on
                // linkage_bores_chamfer: 93 of 1601 shared vertices left their
                // snapshot by up to 12.45 mm with tolerances up to 13.71 mm
                // (sewTol 0.0018), which is the same effect as the F7 note on
                // the R1 path ("displaces shared verts[], measured 18.8 mm on
                // S02"). Orientation only; wire geometry is not ours to move.
                sff.FixWireMode() = 0;
                sff.FixSmallAreaWireMode() = 0;
                sff.FixIntersectingWiresMode() = 0;
                sff.FixLoopWiresMode() = 0;
                sff.FixSplitFaceMode() = 0;
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
        ensureFaceValid(outF, plateFaceTolCap(mv, r));
        if (r.id >= 0 && (diagPlatesEnabled() || diagP2Enabled())) {
            dumpPcurveFrames(r, ow);
            dumpDiagPlateLine(r.id, outF);
        }
        TopoDS_Edge tripBad;
        (void)plateAnalyticEdgesWithoutChain(r, outF, geom, &tripBad);
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
    const Loop* outer = nullptr;
    bool cylFailPrinted = false;
    auto printCylFail = [&](const char* stage, const char* reason, const TopoDS_Wire& wire,
                            const Handle(Geom_Surface)& surf, const TopoDS_Face& cand) {
        if (cylFailPrinted) return;
        cylFailPrinted = true;
        dumpDiagCylFail(r, outer, stage, reason, wire, surf, cand);
    };
    if (!isAnalytic(&r)) {
        diagPartial("eprime-admission");
        printCylFail("other", "eprime-admission", TopoDS_Wire(), Handle(Geom_Surface)(),
                     TopoDS_Face());
        return false;
    }
    const char* surfTag = "unset";
    const char* wireTag = "fwd";
    std::vector<TopoDS_Wire> inners;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::Outer) outer = &lp;
        else if (lp.role == LoopRole::Inner) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, lp, rs, mv, geom, collapsed, meshE, edgeOk, nullptr)) {
                printCylFail("keep-wire", "inner-wire", iw, Handle(Geom_Surface)(), TopoDS_Face());
                return false;
            }
            inners.push_back(iw);
        }
    }
    if (!outer) {
        diagPartial("no-outer");
        printCylFail("other", "no-outer", TopoDS_Wire(), Handle(Geom_Surface)(), TopoDS_Face());
        return false;
    }
    TopoDS_Wire ow;
    auto refreshOuterWire = [&]() -> bool {
        if (!buildLoopWire(ow, *outer, rs, mv, geom, collapsed, meshE, edgeOk, nullptr)) return false;
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
        printCylFail("keep-wire", "loop-wire", ow, Handle(Geom_Surface)(), TopoDS_Face());
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
        if (r.id == 100) diagRid100Status(tagbuf, r.id, cand);
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
                stampTolWriter(eW, "other");
            }
        };
        fillMissingPc(ow);
        for (auto& iw : inners) fillMissingPc(iw);
        auto attempt = [&]() -> bool {
            return makeFaceCopy(surf, ow, inners, r.outwardNormal, cand);
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
                Handle(Geom_RectangularTrimmedSurface) trim =
                    Handle(Geom_RectangularTrimmedSurface)::DownCast(regionSurf(r, SurfVar::CylRectTrim));
                if (!trim.IsNull()) {
                    TopoDS_Face cand;
                    bool got = makeFaceBound(trim, cand);
                    if (!got && (mv.nTri < 500 || refreshOuterWire()))
                        got = makeFaceKeep(trim, ow, inners, r.outwardNormal, cand);
                    if (got && finishPartial(cand)) {
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
    if (!cylFailPrinted && diagP2Enabled()) {
        Handle(Geom_Surface) surf = regionSurf(r, SurfVar::CylBase);
        const char* stage = "other";
        const char* reason = "all-failed";
        TopoDS_Face probe;
        if (ow.IsNull() || (!ow.Closed() && !BRep_Tool::IsClosed(ow))) {
            stage = "keep-wire";
            reason = "short";
        } else if (!makeFaceKeep(surf, ow, inners, r.outwardNormal, probe)) {
            stage = "keep-face";
            reason = "makeface-fail";
        } else {
            const int pcm = pcurveMissingOnWireSurf(ow, surf);
            if (pcm > 0) {
                stage = "pcurve";
                reason = "missing-pc";
            } else {
                TopoDS_Face want = probe, other = probe;
                setFaceOutward(want, r.outwardNormal);
                setFaceOutward(other, !r.outwardNormal);
                if (faceIsValid(want))
                    probe = want;
                else if (faceIsValid(other))
                    probe = other;
                else
                    probe = want;
                if (!faceIsValid(probe)) {
                    stage = "brepcheck";
                    reason = "invalid";
                }
            }
        }
        printCylFail(stage, reason, ow, surf, probe);
    }
    if (diagP2Enabled()) {
        if (!outF.IsNull())
            dumpDiagWiresOfFace(r.id, outF, sewTol, "buildPartialCylinder");
        else if (!cylFailPrinted) {
            int nE = 0;
            if (!ow.IsNull()) {
                for (TopoDS_Iterator it(ow); it.More(); it.Next())
                    if (it.Value().ShapeType() == TopAbs_EDGE) nE++;
            }
            std::fprintf(stderr,
                         "DIAG_WIRE rid=%d iw=0 nE=%d ori=? faceOri=? infPt=? areaUV=0 nSelfX=0 "
                         "maxUVgap=0 nDup=0 pcurveMissing=-1 nCylE=0 nFillE=0 nPlnE=0 nSeamFree=0 "
                         "sewTol=%.6e stage=buildPartialCylinder\n",
                         r.id, nE, sewTol);
        }
    }
    return false;
}

TopoDS_Edge makeFullCircle(const gp_Circ& circ, const TopoDS_Vertex& V) {
    if (!V.IsNull()) {
        const gp_Pnt p = BRep_Tool::Pnt(V);
        double t = 0.0;
        try {
            t = ElCLib::Parameter(circ, p);
        } catch (const Standard_Failure&) {
            t = 0.0;
        }
        const double d = ElCLib::Value(t, circ).Distance(p);
        bumpVertexTol(V, std::max(d, Precision::Confusion()));
    }
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
    if (!isAnalytic(&r)) return false;
    const Loop *capL = nullptr, *capH = nullptr;
    std::vector<TopoDS_Wire> inners;
    std::vector<double> innerU;
    for (const Loop& lp : r.loops) {
        if (lp.role == LoopRole::CapLow) capL = &lp;
        else if (lp.role == LoopRole::CapHigh) capH = &lp;
        else if (lp.role == LoopRole::Inner) {
            TopoDS_Wire iw;
            if (!buildLoopWire(iw, lp, rs, mv, geom, collapsed, meshE, edgeOk, nullptr)) return false;
            inners.push_back(std::move(iw));
            if (unionBuildOn()) {
                for (int ci : lp.chainIdx) {
                    if (ci < 0 || (size_t)ci >= rs.chains.size()) continue;
                    for (int lv : rs.chains[(size_t)ci].meshVerts)
                        innerU.push_back(azimuthOf(r, pntOf(mv, lv)));
                }
            }
        }
    }
    if (!capL || !capH || capL->chainIdx.empty() || capH->chainIdx.empty()) return false;

    // D-130-14: a union face carries the enclosed interruptions as inner wires,
    // and the seam is a generator of the SAME face -- it may not run through
    // one. The seam azimuth is therefore the middle of the widest azimuth gap
    // between the interruptions (a measurement of this region, not a constant);
    // with no inner wire this is u = 0 exactly as before and nothing moves.
    double seamU = 0.0;
    if (!innerU.empty()) {
        std::sort(innerU.begin(), innerU.end());
        double gap = innerU.front() + 2.0 * kPi - innerU.back();
        double mid = innerU.back() + 0.5 * gap;
        for (size_t i = 1; i < innerU.size(); i++) {
            const double g = innerU[i] - innerU[i - 1];
            if (g > gap) {
                gap = g;
                mid = innerU[i - 1] + 0.5 * g;
            }
        }
        while (mid > kPi) mid -= 2.0 * kPi;
        while (mid < -kPi) mid += 2.0 * kPi;
        seamU = mid;
    }
    int vL = vertexClosestToUOnLoop(mv, r, *capL, rs, seamU);
    int vH = vertexClosestToUOnLoop(mv, r, *capH, rs, seamU);
    if (!innerU.empty() && vL >= 0) {
        // The seam edge is ONE generator: both cap vertices must lie on it.
        // Nearest-to-seamU on each cap independently can pick two vertices a
        // facet apart on a staggered wall and MakeEdge then refuses. Take the
        // low cap's choice and ask the high cap for the vertex nearest THAT
        // generator -- no new tolerance; the existing snap does the rest.
        const gp_Lin gen(pntOf(mv, vL), r.ax.Direction());
        int bh = -1;
        double bestD = 1e300;
        for (int ci : capH->chainIdx) {
            if (ci < 0 || (size_t)ci >= rs.chains.size()) continue;
            for (int b2 : rs.chains[(size_t)ci].meshVerts) {
                const double d = gen.Distance(pntOf(mv, b2));
                if (d < bestD || (d == bestD && b2 < bh)) {
                    bestD = d;
                    bh = b2;
                }
            }
        }
        if (bh >= 0) vH = bh;
    }
    if (diagP2Enabled() && !innerU.empty())
        std::fprintf(stderr,
                     "DIAG_SEAMPICK rid=%d seamU=%.5f vL=%d uL=%.5f vH=%d uH=%.5f "
                     "innerU=[%.5f..%.5f] n=%zu\n",
                     r.id, seamU, vL, vL >= 0 ? azimuthOf(r, pntOf(mv, vL)) : -1.0, vH,
                     vH >= 0 ? azimuthOf(r, pntOf(mv, vH)) : -1.0, innerU.front(), innerU.back(),
                     innerU.size());
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

    // D-130-16: the face's parameter domain is [u0, u0+2pi] and the SEAM names
    // it. A cylinder is 2pi-periodic in u, so a cap arc bound on another branch
    // is the same curve written as a different number, and the wire then does
    // not close in UV. Measured on handle-pickup's R=4 wall with the rim split
    // in place: the low cap arrived on [2pi, 4pi] while the high cap and the
    // seam were on [0, 2pi], every wire itself clean, and BRepCheck answered
    // UnorientableShape on the face. Shifting a pcurve by a WHOLE PERIOD is
    // exact -- no tolerance is spent, no geometry moves, and the 3d edge and
    // its neighbour face's pcurve are untouched.
    TopoDS_Edge gSeamE;
    double gSeamU0 = 0.0;
    bool gSeamSet = false;
    auto rebranchEdges = [&](const TopoDS_Shape& w, double u0, const TopoDS_Edge& skip) {
        BRep_Builder B;
        for (TopExp_Explorer it(w, TopAbs_EDGE); it.More(); it.Next()) {
            TopoDS_Edge e = TopoDS::Edge(it.Current());
            if (!skip.IsNull() && e.IsSame(skip)) continue;  // the seam pair is deliberate
            Standard_Real f = 0, l = 0;
            Handle(Geom2d_Curve) pc =
                BRep_Tool::CurveOnSurface(e, surf, TopLoc_Location(), f, l);
            if (pc.IsNull()) continue;
            const double uMid = 0.5 * (pc->Value(f).X() + pc->Value(l).X());
            const double k = std::floor((uMid - u0) / (2.0 * kPi));
            if (k == 0.0) continue;
            Handle(Geom2d_Curve) shifted = Handle(Geom2d_Curve)::DownCast(
                pc->Translated(gp_Vec2d(-k * 2.0 * kPi, 0.0)));
            if (shifted.IsNull()) continue;
            B.UpdateEdge(e, shifted, surf, TopLoc_Location(), sewTol);
            B.Range(e, surf, TopLoc_Location(), f, l);
            if (diagP2Enabled())
                std::fprintf(stderr, "DIAG_REBRANCH rid=%d u0=%.5f uMid=%.5f k=%+.0f\n", r.id,
                             u0, uMid, k);
        }
    };

    auto finishFace = [&](TopoDS_Face got, bool publishSimple, int ciL, int ciH,
                          const TopoDS_Edge& eL, const TopoDS_Edge& eH) -> bool {
        if (got.IsNull()) return false;
        setFaceOutward(got, r.outwardNormal);
        addPcurvesOnFace(got, sewTol, true);
        // D-130-16: addPcurvesOnFace is what BIRTHS the inner wire's pcurves,
        // so the branch pass has to run after it too -- an inner loop written
        // on [0, 2pi] under a seam that named [2pi, 4pi] is outside its own
        // face's domain, which is the UnorientableShape.
        if (gSeamSet) rebranchEdges(got, gSeamU0, gSeamE);
        if (gSeamSet && diagP2Enabled()) {
            char st[128] = "[]";
            try {
                BRepCheck_Analyzer an(got, Standard_True);
                formatStatusList(an.Result(got), st, sizeof(st));
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr, "DIAG_REBRANCH_FACE rid=%d valid=%d st=%s\n", r.id,
                         faceIsValid(got) ? 1 : 0, st);
            try {
                TopLoc_Location L;
                Handle(Geom_Surface) sF = BRep_Tool::Surface(got, L);
                int wi = 0;
                for (TopExp_Explorer wx(got, TopAbs_WIRE); wx.More(); wx.Next(), ++wi) {
                    int ei = 0;
                    double uLo = 1e300, uHi = -1e300;
                    for (BRepTools_WireExplorer ex(TopoDS::Wire(wx.Current())); ex.More();
                         ex.Next(), ++ei) {
                        const TopoDS_Edge ee = TopoDS::Edge(ex.Current());
                        Standard_Real a2 = 0, b2 = 0;
                        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(ee, sF, L, a2, b2);
                        if (pc.IsNull()) continue;
                        uLo = std::min({uLo, pc->Value(a2).X(), pc->Value(b2).X()});
                        uHi = std::max({uHi, pc->Value(a2).X(), pc->Value(b2).X()});
                        if (ee.IsSame(gSeamE))
                            std::fprintf(stderr,
                                         "DIAG_RB_SEAM w=%d e=%d ori=%s u0=%.5f u1=%.5f\n", wi,
                                         ei, ee.Orientation() == TopAbs_FORWARD ? "F" : "R",
                                         pc->Value(a2).X(), pc->Value(b2).X());
                    }
                    std::fprintf(stderr,
                                 "DIAG_RB_WIRE rid=%d w=%d ori=%s nE=%d u=[%.5f,%.5f]\n", r.id,
                                 wi, wx.Current().Orientation() == TopAbs_FORWARD ? "F" : "R",
                                 ei, uLo, uHi);
                }
            } catch (const Standard_Failure&) {
            }
        }
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
            if (gSeamSet) rebranchEdges(got, gSeamU0, gSeamE);
        }
        if (!faceIsValid(got) && !ensureFaceValid(got, meshTolCap(mv, &r))) {
            if (diagP2Enabled()) {
                char faceSt[256] = "[]", wireSt[256] = "[]";
                try {
                    BRepCheck_Analyzer an(got, Standard_True);
                    formatStatusList(an.Result(got), faceSt, sizeof(faceSt));
                    size_t used = 0;
                    wireSt[0] = '\0';
                    for (TopExp_Explorer wx(got, TopAbs_WIRE); wx.More(); wx.Next()) {
                        char one[128] = "[]";
                        const TopoDS_Wire ww = TopoDS::Wire(wx.Current());
                        formatStatusList(an.Result(ww), one, sizeof(one));
                        int ne = 0;
                        for (TopExp_Explorer ex(ww, TopAbs_EDGE); ex.More(); ex.Next()) ne++;
                        used += (size_t)std::snprintf(wireSt + used, sizeof(wireSt) - used,
                                                      "w(nE=%d)%s ", ne, one);
                        if (used + 32 >= sizeof(wireSt)) break;
                    }
                } catch (const Standard_Failure&) {
                }
                std::fprintf(stderr,
                             "DIAG_SEAM360FAIL rid=%d R=%.4f inners=%zu faceSt=%s wireSt=%s\n",
                             r.id, r.radius, inners.size(), faceSt, wireSt);
                try {
                    TopLoc_Location L;
                    Handle(Geom_Surface) sF = BRep_Tool::Surface(got, L);
                    int wi = 0;
                    for (TopExp_Explorer wx(got, TopAbs_WIRE); wx.More(); wx.Next(), ++wi) {
                        int ei = 0;
                        for (BRepTools_WireExplorer ex(TopoDS::Wire(wx.Current())); ex.More();
                             ex.Next(), ++ei) {
                            const TopoDS_Edge ee = TopoDS::Edge(ex.Current());
                            Standard_Real a2 = 0, b2 = 0;
                            Handle(Geom2d_Curve) pc =
                                BRep_Tool::CurveOnSurface(ee, sF, L, a2, b2);
                            if (pc.IsNull()) {
                                std::fprintf(stderr, "DIAG_SEAMUV w=%d e=%d NOPC\n", wi, ei);
                                continue;
                            }
                            const gp_Pnt2d q0 = pc->Value(a2), q1 = pc->Value(b2);
                            std::fprintf(stderr,
                                         "DIAG_SEAMUV w=%d e=%d ori=%s uv0=(%.5f,%.5f) "
                                         "uv1=(%.5f,%.5f) rng=[%.5f,%.5f] pc=%s c3=%s\n",
                                         wi, ei,
                                         ee.Orientation() == TopAbs_FORWARD ? "F" : "R", q0.X(),
                                         q0.Y(), q1.X(), q1.Y(), a2, b2,
                                         pc->DynamicType()->Name(),
                                         [&]() -> const char* {
                                             Standard_Real f3 = 0, l3 = 0;
                                             Handle(Geom_Curve) cc = BRep_Tool::Curve(ee, f3, l3);
                                             return cc.IsNull() ? "null"
                                                                : cc->DynamicType()->Name();
                                         }());
                        }
                    }
                } catch (const Standard_Failure&) {
                }
            }
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

    // D-130-14: an inner wire's sense on the face is opposite the outer's, and
    // the wire P1 hands over is oriented for the REGION boundary, not for the
    // face. makeFaceKeep adds what it is given; try the other sense when the
    // first one leaves the face unorientable. No inner wire => unchanged.
    auto buildWithInners = [&](const TopoDS_Wire& ow2, TopoDS_Face& out2) -> bool {
        // makeFaceCopy COPIES the edges (J2) and OCCT recomputes their pcurves
        // on the copy -- which throws away the seam's pcurve PAIR and the
        // monotone cap branch this face is assembled on. A union face carries
        // an inner wire and cannot survive that, so it takes makeFaceKeep
        // first; faces without an inner wire keep the order they had.
        const bool keepFirst = unionBuildOn() && !inners.empty();
        if (keepFirst) {
            if (!makeFaceKeep(surf, ow2, inners, r.outwardNormal, out2) &&
                !makeFaceCopy(surf, ow2, inners, r.outwardNormal, out2))
                return false;
        } else if (!makeFaceCopy(surf, ow2, inners, r.outwardNormal, out2) &&
                   !makeFaceKeep(surf, ow2, inners, r.outwardNormal, out2)) {
            return false;
        }
        if (!unionBuildOn() || inners.empty()) return true;
        // D-130-16: an inner wire's sense on the FACE is opposite the outer's,
        // and the wire P1 hands over is oriented for the REGION boundary. Which
        // sense is right cannot be judged on the bare face: the inner wire's
        // pcurves do not exist until addPcurvesOnFace births them, and until
        // they are on the seam's own branch the face is unorientable either
        // way. So prepare the face fully -- pcurves, branch, outward normal --
        // and let BRepCheck answer; only then try the other sense.
        auto prepare = [&](TopoDS_Face& f) -> bool {
            if (f.IsNull()) return false;
            setFaceOutward(f, r.outwardNormal);
            addPcurvesOnFace(f, sewTol, true);
            if (gSeamSet) rebranchEdges(f, gSeamU0, gSeamE);
            return faceIsValid(f);
        };
        if (prepare(out2)) return true;
        std::vector<TopoDS_Wire> rev;
        rev.reserve(inners.size());
        for (const TopoDS_Wire& iw : inners) rev.push_back(TopoDS::Wire(iw.Reversed()));
        TopoDS_Face alt;
        if ((makeFaceKeep(surf, ow2, rev, r.outwardNormal, alt) ||
             makeFaceCopy(surf, ow2, rev, r.outwardNormal, alt)) &&
            prepare(alt))
            out2 = alt;
        return true;
    };

    TopoDS_Edge eSeam;
    TopoDS_Vertex seamVL = verts[(size_t)vL];
    try {
        gp_Lin lin(pntOf(mv, vL), r.ax.Direction());
        AnalyticCurve acs;
        acs.kind = AnalyticCurve::Lin;
        acs.lin = lin;
        snapVertexToCurve(verts[(size_t)vL], acs, snapCap);
        snapVertexToCurve(verts[(size_t)vH], acs, snapCap);
        BRepBuilderAPI_MakeEdge ms(lin, verts[(size_t)vL], verts[(size_t)vH]);
        if (!ms.IsDone()) {
            if (diagP2Enabled())
                std::fprintf(stderr,
                             "DIAG_SEAMEDGE rid=%d vL=%d vH=%d dist=%.6g snapCap=%.6g "
                             "tolL=%.6g tolH=%.6g\n",
                             r.id, vL, vH, lin.Distance(BRep_Tool::Pnt(verts[(size_t)vH])),
                             snapCap, BRep_Tool::Tolerance(verts[(size_t)vL]),
                             BRep_Tool::Tolerance(verts[(size_t)vH]));
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
            TopoDS_Face got;
            if (!buildWithInners(w, got)) {
                emit(warn, "seamed360: MakeFace not done");
                return false;
            }
            return finishFace(got, true, ciL0, ciH0, eL, eH);
        }

        // Live I8 splits a cap into N chains (S04 boss-top torus junction). Build
        // the actual cap wires and insert the seam at u≈0; do not invent a
        // single circle that would steal neighbour TShapes.
        TopoDS_Wire wL, wH;
        if (!buildLoopWire(wL, *capL, rs, mv, geom, collapsed, meshE, edgeOk, nullptr) ||
            !buildLoopWire(wH, *capH, rs, mv, geom, collapsed, meshE, edgeOk, nullptr)) {
            emit(warn, "seamed360: composite cap wire failed — try TwoHalves");
            return false;
        }
        std::vector<TopoDS_Edge> pathL, pathH;
        if (!collectWireEdges(wH, pathH) || !collectWireEdges(wL, pathL)) {
            emit(warn, "seamed360: composite cap wire failed — try TwoHalves");
            return false;
        }
        if ((unionBuildOn() && !inners.empty()) ||
            !rotateEdgesToVertex(pathH, verts[(size_t)vH]) ||
            !rotateEdgesToVertex(pathL, verts[(size_t)vL])) {
            if (!unionBuildOn()) {
                emit(warn, "seamed360: cap wire does not pass seam vertex — try TwoHalves");
                return false;
            }
            // D-130-14: a union face's caps are chained through pinch vertices,
            // and the mesh vertex nearest u=0 on the LOOP need not survive into
            // the built wire (a collapsed chain publishes one analytic edge and
            // its interior mesh vertices are gone). The seam only needs a
            // generator both caps actually pass through, so ask the wires
            // themselves for the vertex pair nearest u=0 before giving up.
            auto wireVertNearestU = [&](const std::vector<TopoDS_Edge>& path) -> TopoDS_Vertex {
                TopoDS_Vertex best;
                double bestD = 1e300;
                for (const TopoDS_Edge& e : path) {
                    TopoDS_Vertex v1, v2;
                    TopExp::Vertices(e, v1, v2, Standard_False);
                    for (const TopoDS_Vertex* V : {&v1, &v2}) {
                        if (V->IsNull()) continue;
                        double u = azimuthOf(r, BRep_Tool::Pnt(*V));
                        double d = std::fabs(u - seamU);
                        d = std::min(d, 2.0 * kPi - d);
                        if (d < bestD) {
                            bestD = d;
                            best = *V;
                        }
                    }
                }
                return best;
            };
            const TopoDS_Vertex VL = wireVertNearestU(pathL);
            // The high cap must meet the SAME generator, not merely the same
            // azimuth target: on a staggered wall those differ by half a facet.
            TopoDS_Vertex VH;
            if (!VL.IsNull()) {
                const gp_Lin genL(BRep_Tool::Pnt(VL), r.ax.Direction());
                double bd = 1e300;
                for (const TopoDS_Edge& e : pathH) {
                    TopoDS_Vertex a1, b1;
                    TopExp::Vertices(e, a1, b1, Standard_False);
                    for (const TopoDS_Vertex* V : {&a1, &b1}) {
                        if (V->IsNull()) continue;
                        const double d = genL.Distance(BRep_Tool::Pnt(*V));
                        if (d < bd) {
                            bd = d;
                            VH = *V;
                        }
                    }
                }
            }
            seamVL = VL;
            if (VH.IsNull() || VL.IsNull() || !rotateEdgesToVertex(pathH, VH) ||
                !rotateEdgesToVertex(pathL, VL)) {
                emit(warn, "seamed360: cap wire does not pass seam vertex — try TwoHalves");
                return false;
            }
            try {
                gp_Lin lin2(BRep_Tool::Pnt(VL), r.ax.Direction());
                AnalyticCurve acs2;
                acs2.kind = AnalyticCurve::Lin;
                acs2.lin = lin2;
                snapVertexToCurve(VL, acs2, snapCap);
                snapVertexToCurve(VH, acs2, snapCap);
                BRepBuilderAPI_MakeEdge ms2(lin2, VL, VH);
                if (!ms2.IsDone()) {
                    emit(warn, "seamed360: seam MakeEdge failed");
                    return false;
                }
                eSeam = ms2.Edge();
            } catch (const Standard_Failure&) {
                emit(warn, "seamed360: seam threw");
                return false;
            }
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
            // D-130-14: the composite path never bound the seam's SECOND
            // pcurve, so the doubled seam edge read as a self-intersecting wire
            // and every composite 360 face fell through to TwoHalves. A union
            // face cannot: TwoHalves has no inner wire to give the enclosed
            // interruption. Bind the pair here exactly as the single-circle
            // path does, at the seam generator's own azimuth.
            if (unionBuildOn()) {
                BRep_Builder Bs;
                const double u0 = azimuthOf(r, BRep_Tool::Pnt(seamVL));
                Standard_Real fs = 0, ls = 0;
                Handle(Geom_Curve) cs = BRep_Tool::Curve(eSeam, fs, ls);
                Handle(Geom2d_Line) pcS0 =
                    new Geom2d_Line(gp_Pnt2d(u0, r.vMin), gp_Dir2d(0.0, 1.0));
                Handle(Geom2d_Line) pcS1 =
                    new Geom2d_Line(gp_Pnt2d(u0 + 2.0 * kPi, r.vMin), gp_Dir2d(0.0, 1.0));
                Bs.UpdateEdge(eSeam, pcS0, pcS1, surf, TopLoc_Location(), sewTol);
                if (!cs.IsNull()) Bs.Range(eSeam, fs, ls);
                eSeam.Closed(Standard_False);
                gSeamE = eSeam;
                gSeamU0 = u0;
                gSeamSet = true;
                rebranchEdges(ow, u0, eSeam);
                for (TopoDS_Wire& iw : inners) rebranchEdges(iw, u0, eSeam);
            }
            bindCylPCurves(ow, surf, r, sewTol);
            TopoDS_Face got;
            if (!buildWithInners(ow, got)) {
                emit(warn, "seamed360: composite MakeFace not done — try TwoHalves");
                return false;
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
    if (!isAnalytic(&r)) return false;
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
            if (!makeFaceCopy(surf, ow, {}, r.outwardNormal, f) &&
                !makeFaceKeep(surf, ow, {}, r.outwardNormal, f))
                return false;
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

// D-130-9(b) cascade containment. A face that is NOT an analytic region -- an
// island facet or an E' fill triangle, builtRid < 0 -- has no region to
// explode, so `collectFaceCulprits` and `collectShellCulprits` both skip it and
// the culprit set comes back EMPTY. An invalid shell with an empty culprit set
// falls through U0 and U1 in one decide and lands on the U2 blanket, which
// explodes every non-closed360 analytic region on the component: on the plate
// one unorientable facet took out 9 partial cylinders (the R10 cross bore,
// R30 x2, R17 x2, R5 x4) that had each built VALID. That is the escalation this
// function removes. The exploded set for such a face is its own NEIGHBOURHOOD --
// the regions whose built faces share an edge with it -- and nothing else.
// RULE 1.5's 50 % test is NOT changed: it still runs, on this set.
void collectNonRegionNeighbourhood(const TopoDS_Shape& sh, const std::vector<TopoDS_Face>& built,
                                   const std::vector<int>& builtRid,
                                   const std::vector<char>& exploded,
                                   std::vector<CascadeHit>& hits) {
    if (sh.IsNull()) return;
    try {
        BRepCheck_Analyzer an(sh, Standard_True);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(sh, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            const int idx = faceMap.FindIndex(built[i]);
            if (idx <= 0 || i >= builtRid.size()) continue;
            faceRid[(size_t)idx] = builtRid[i];
        }
        TopTools_IndexedDataMapOfShapeListOfShape eanc;
        TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, eanc);
        for (int i = 1; i <= faceMap.Extent(); i++) {
            if (faceRid[(size_t)i] >= 0) continue;  // an analytic region: normal path
            const TopoDS_Face fs = TopoDS::Face(faceMap(i));
            bool bad = brepStatusBad(an.Result(fs));
            if (!bad) {
                BRepCheck_Analyzer fa(fs, Standard_True);
                bad = !fa.IsValid();
            }
            if (!bad) continue;
            std::vector<int> nbr;
            for (TopExp_Explorer ex(fs, TopAbs_EDGE); ex.More(); ex.Next()) {
                const int ei = eanc.FindIndex(ex.Current());
                if (ei <= 0) continue;
                for (TopTools_ListOfShape::Iterator it(eanc(ei)); it.More(); it.Next()) {
                    const int oi = faceMap.FindIndex(it.Value());
                    if (oi <= 0 || oi == i) continue;
                    const int orid = faceRid[(size_t)oi];
                    if (orid < 0 || regionExploded(exploded, orid)) continue;
                    if (std::find(nbr.begin(), nbr.end(), orid) == nbr.end()) nbr.push_back(orid);
                }
            }
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
            if (collapseDiagEnabled()) {
                std::fprintf(stderr, "DIAG_CASCADE nonregion face=%d nNbr=%zu nbrs=[", i,
                             nbr.size());
                for (size_t j = 0; j < nbr.size(); j++)
                    std::fprintf(stderr, "%s%d", j ? "," : "", nbr[j]);
                std::fprintf(stderr, "]\n");
            }
            for (int id : nbr) addCascadeHit(hits, id, "nonregion", true);
        }
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
    std::vector<int> cyls, smallPln, hubs, skipC360;
    std::vector<int> badPlanes;
    for (const CascadeHit& h : culprits) {
        if (h.rid < 0 || regionExploded(exploded, h.rid)) continue;
        const Region* r = regionById(rs, h.rid);
        if (!r) continue;
        if (r->type != SurfType::Plane && r->type != SurfType::Cylinder) continue;
        if (r->type == SurfType::Cylinder) {
            if (r->closed360) continue;
            cyls.push_back(h.rid);
        }
        else if (touchesClosed360Cyl(rs, h.rid))
            skipC360.push_back(h.rid);
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
    uniqueSorted(skipC360);
    if (collapseDiagEnabled()) {
        const int nCyl = countAnalyticCylinders(rs);
        auto emitU0 = [&](int rid, const char* bucket) {
            const char* src = "-";
            int faceValid = 1;
            for (const CascadeHit& h : culprits) {
                if (h.rid != rid) continue;
                src = h.src ? h.src : "-";
                faceValid = h.faceValid ? 1 : 0;
                break;
            }
            int hop = 0;
            for (const BoundaryChain& ch : rs.chains) {
                if (!chainTouches(ch, rid)) continue;
                const Region* o = regionById(rs, otherReg(ch, rid));
                if (o && o->type == SurfType::Cylinder) ++hop;
            }
            std::fprintf(stderr,
                         "DIAG_CASCADE u0 rid=%d bucket=%s hop=%d nCyl=%d hub=%d touchC360=%d "
                         "culpritSrc=%s faceValid=%d\n",
                         rid, bucket, hop, nCyl, isHubPlane(rs, rid) ? 1 : 0,
                         touchesClosed360Cyl(rs, rid) ? 1 : 0, src, faceValid);
            if (rid == 100) {
                char pass[48];
                std::snprintf(pass, sizeof(pass), "u0-%s", bucket);
                if (faceValid)
                    gRid100PrevValid = 1;
                else if (gRid100PrevValid == 1)
                    std::fprintf(stderr, "DIAG_RID100 first-invalid pass=%s\n", pass);
                std::fprintf(stderr, "DIAG_RID100 pass=%s status=faceValid=%d valid=%d\n", pass,
                             faceValid, faceValid);
                if (!faceValid) gRid100PrevValid = 0;
            }
        };
        for (int id : cyls) emitU0(id, "cyls");
        for (int id : skipC360) emitU0(id, "skip-C360");
        for (int id : hubs) emitU0(id, "hubs");
        for (int id : smallPln) emitU0(id, "smallPln");
    }
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

int hopCylCount(const RegionSet& rs, int rid) {
    int hop = 0;
    for (const BoundaryChain& ch : rs.chains) {
        if (!chainTouches(ch, rid)) continue;
        const Region* o = regionById(rs, otherReg(ch, rid));
        if (o && o->type == SurfType::Cylinder) ++hop;
    }
    return hop;
}

void dumpR2Probe(const MeshView& mv, const RegionSet& rs, const std::vector<TopoDS_Face>& built,
                 const std::vector<int>& builtRid, const std::vector<char>& eprimeFill,
                 const std::vector<char>& exploded, double sewTol);

void dumpShippedInContext(const std::vector<TopoDS_Face>& built, const std::vector<int>& builtRid,
                          const RegionSet& rs) {
    if (!diagP2Enabled() && !diagPlatesEnabled()) return;
    int nPlane = 0, nFlagged = 0;
    for (size_t i = 0; i < built.size(); i++) {
        const TopoDS_Face& f = built[i];
        if (f.IsNull()) continue;
        const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
        const Region* r = regionById(rs, rid);
        if (!r || r->type != SurfType::Plane || r->builtAs == BuiltAs::ExplodedToFacets) continue;
        nPlane++;
        try {
            BRepCheck_Analyzer an(f, Standard_True);
            TopTools_IndexedMapOfShape emap;
            TopExp::MapShapes(f, TopAbs_EDGE, emap);
            int flagged = 0;
            for (int ie = 1; ie <= emap.Extent(); ie++) {
                Handle(BRepCheck_Result) res = analyzerResultOf(an, TopoDS::Edge(emap(ie)));
                if (inContextHas8or11(res, f)) flagged++;
            }
            gIncontextShippedExamined++;
            if (flagged > 0) {
                nFlagged++;
                gIncontextShippedFlagged++;
            }
            std::fprintf(stderr, "DIAG_INCONTEXT_SHIPPED rid=%d flagged=%d/%d\n", rid, flagged,
                         emap.Extent());
        } catch (const Standard_Failure&) {
        }
    }
    std::fprintf(stderr,
                 "DIAG_INCONTEXT_SHIPPED_SUM plates=%d flaggedPlates=%d examined=%d flagged=%d\n",
                 nPlane, nFlagged, gIncontextShippedExamined, gIncontextShippedFlagged);
}

void dumpTolBindSummary() {
    if (!diagP2Enabled() && !diagPlatesEnabled()) return;
    std::fprintf(stderr, "DIAG_BINDTOL_SUM unhandled=%d overCap=%d overCapUnique=%zu genuine=%d cylSideChanged=%d\n",
                 gBindTolUnhandled, gBindTolOverCap, gBindTolOverCapTShapes.size(),
                 gBindTolOverCapGenuine, gCylSideTolChanged);
    for (const auto& kv : gBindTolUnhandledByClass)
        std::fprintf(stderr, "DIAG_BINDTOL_UNHANDLED class=%s count=%zu\n", kv.first.c_str(),
                     kv.second.size());
    std::fprintf(stderr,
                 "DIAG_TOLCHECK_SUM rows=%d brepLeExact=%d plane=%d/%d cyl=%d/%d brepGtExact=%d\n",
                 gTolCheckRows, gTolCheckBrepLeExact, gTolCheckPlaneBrepLe, gTolCheckPlaneRows,
                 gTolCheckCylBrepLe, gTolCheckCylRows, gBindTolBrepGtExact);
    std::fprintf(stderr,
                 "DIAG_INCONTEXT_SUM construction=%d/%d shipped=%d/%d\n",
                 gIncontextPlateFlagged, gIncontextPlateExamined, gIncontextShippedFlagged,
                 gIncontextShippedExamined);
    std::fprintf(stderr, "DIAG_REPAIR_SUM reconnect=%d tolrewrite=%d\n", gReconnectFire,
                 gTolRewriteFire);
    if (!gHpPlateEdgeTol.empty()) {
        std::vector<double> tols;
        tols.reserve(gHpPlateEdgeTol.size());
        for (const auto& kv : gHpPlateEdgeTol) tols.push_back(kv.second);
        std::sort(tols.begin(), tols.end());
        int nAbove1um = 0;
        for (double t : tols)
            if (t > 0.001) nAbove1um++;
        const double med = tols[tols.size() / 2];
        std::fprintf(stderr,
                     "DIAG_TOLCHECK_DIST n=%zu min=%.9f median=%.9f max=%.9f nAbove1um=%d\n",
                     tols.size(), tols.front(), med, tols.back(), nAbove1um);
    }
    if (!gAllPlateEdgeTol.empty()) {
        std::vector<double> tols;
        tols.reserve(gAllPlateEdgeTol.size());
        for (const auto& kv : gAllPlateEdgeTol) tols.push_back(kv.second);
        std::sort(tols.begin(), tols.end());
        int nAbove1um = 0;
        for (double t : tols)
            if (t > 0.001) nAbove1um++;
        const double med = tols[tols.size() / 2];
        std::fprintf(stderr,
                     "DIAG_TOLCHECK_DIST_ALL70 n=%zu min=%.9f median=%.9f max=%.9f nAbove1um=%d\n",
                     tols.size(), tols.front(), med, tols.back(), nAbove1um);
    }
}

void dumpR2Probe(const MeshView& mv, const RegionSet& rs, const std::vector<TopoDS_Face>& built,
                 const std::vector<int>& builtRid, const std::vector<char>& eprimeFill,
                 const std::vector<char>& exploded, double sewTol) {
    if (!diagP2Enabled()) return;
    try {
        TopoDS_Shell probe;
        BRep_Builder pb;
        pb.MakeShell(probe);
        for (const auto& f : built) {
            if (!f.IsNull()) pb.Add(probe, f);
        }
        const int closed = BRep_Tool::IsClosed(probe) ? 1 : 0;
        int valid = 0;
        BRepCheck_Analyzer an(probe, Standard_True, Standard_True);
        try {
            valid = an.IsValid() ? 1 : 0;
        } catch (const Standard_Failure&) {
            valid = 0;
        }
        double dVolAbs = 0;
        for (const Region& reg : rs.regions) dVolAbs += std::fabs(reg.dVolPredicted);
        const double meshVol = std::fabs(meshViewVolume(mv));
        const double budget = std::max(1e-4 * meshVol, 3.0 * dVolAbs);
        double shellVol = 0;
        int volPass = 0;
        try {
            GProp_GProps gp;
            BRepGProp::VolumeProperties(probe, gp);
            shellVol = gp.Mass();
            volPass = (std::fabs(shellVol - meshVol) <= budget) ? 1 : 0;
        } catch (const Standard_Failure&) {
            volPass = 0;
        }
        std::fprintf(stderr,
                     "DIAG_R2PROBE comp=%zu nFaces=%d closed=%d valid=%d shellVol=%.6f "
                     "meshVol=%.6f budget=%.6f volPass=%d\n",
                     mv.nTri, (int)built.size(), closed, valid, shellVol, meshVol, budget, volPass);
        dumpM5OnShell(an, built, builtRid, rs);
        if (valid) return;
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(probe, TopAbs_EDGE, TopAbs_FACE, anc);
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(probe, TopAbs_FACE, faceMap);
        std::vector<int> faceRid((size_t)faceMap.Extent() + 1, -1);
        for (size_t i = 0; i < built.size(); i++) {
            int idx = faceMap.FindIndex(built[i]);
            if (idx > 0 && i < builtRid.size()) faceRid[(size_t)idx] = builtRid[i];
        }
        for (size_t i = 0; i < built.size(); i++) {
            const TopoDS_Face& f = built[i];
            if (f.IsNull()) continue;
            Handle(BRepCheck_Result) fr = an.Result(f);
            BRepCheck_Analyzer fan(f, Standard_True);
            bool faceBad = !fan.IsValid();
            TopTools_IndexedMapOfShape subMap;
            TopExp::MapShapes(f, subMap);
            if (!fr.IsNull()) {
                for (BRepCheck_ListOfStatus::Iterator it(fr->Status()); it.More(); it.Next())
                    if (it.Value() != BRepCheck_NoError) faceBad = true;
            }
            if (!faceBad) continue;
            const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
            const Region* rr = regionById(rs, rid);
            const char* kind = "other";
            const char* bas = "NotBuilt";
            if (rid >= 0 && (size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid])
                kind = "fill";
            else if (rr) {
                if (rr->type == SurfType::Plane) kind = "plane";
                else if (rr->type == SurfType::Cylinder) kind = "cyl";
                bas = builtAsName(rr->builtAs);
                if (regionExploded(exploded, rid) || rr->builtAs == BuiltAs::ExplodedToFacets)
                    kind = "facet";
            } else
                kind = "facet";
            char faceSt[128], wireSt[128], edgeSt[128];
            formatStatusList(fr, faceSt, sizeof(faceSt));
            wireSt[0] = '[';
            wireSt[1] = ']';
            wireSt[2] = '\0';
            edgeSt[0] = '[';
            edgeSt[1] = ']';
            edgeSt[2] = '\0';
            try {
                for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next()) {
                    formatStatusList(an.Result(wx.Current()), wireSt, sizeof(wireSt));
                    break;
                }
                for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
                    Handle(BRepCheck_Result) er = an.Result(ex.Current());
                    bool bad = false;
                    if (!er.IsNull()) {
                        for (BRepCheck_ListOfStatus::Iterator it(er->Status()); it.More();
                             it.Next())
                            if (it.Value() != BRepCheck_NoError) bad = true;
                    }
                    if (bad) {
                        formatStatusList(er, edgeSt, sizeof(edgeSt));
                        break;
                    }
                }
            } catch (const Standard_Failure&) {
            }
            std::fprintf(stderr,
                         "DIAG_R2INVALID rid=%d kind=%s builtAs=%s faceSt=%s wireSt=%s edgeSt=%s "
                         "nE=%d pcurveMissing=%d\n",
                         rid, kind, bas, faceSt, wireSt, edgeSt, countFaceEdges(f),
                         pcurveMissingOnFace(f));
            // Name the sub-shape that made this face invalid, and separate a
            // WRONG GEOMETRY from a WRONG FLAG. `an.Result(edge)` above carries
            // only the context-free status list; an edge that is invalid ONLY in
            // the context of this face (InvalidSameParameterFlag,
            // InvalidCurveOnSurface) records its status in the context list and
            // is otherwise invisible -- which is how a bad cone rim reads as an
            // empty faceSt=[] on its neighbour cylinder. `dev` is
            // sup|C3d(t) - S(pc(t))| on a 65-point grid: dev = 0 with sameP = 0
            // is a flag BRepLib could not set (some OTHER face's pcurve on the
            // same edge disagrees), not a gap on this one.
            for (int si = 1; si <= subMap.Extent(); si++) {
                const TopoDS_Shape& sh = subMap(si);
                Handle(BRepCheck_Result) sr = fan.Result(sh);
                if (sr.IsNull()) continue;
                char ctxSt[128];
                ctxSt[0] = '\0';
                size_t cu = 0;
                for (sr->InitContextIterator(); sr->MoreShapeInContext();
                     sr->NextShapeInContext()) {
                    for (BRepCheck_ListOfStatus::Iterator it(sr->StatusOnShape()); it.More();
                         it.Next()) {
                        if (it.Value() == BRepCheck_NoError) continue;
                        const int w = std::snprintf(ctxSt + cu, sizeof(ctxSt) - cu, "%s%s",
                                                    cu ? "," : "", brepCheckName((int)it.Value()));
                        if (w < 0 || (size_t)w >= sizeof(ctxSt) - cu) break;
                        cu += (size_t)w;
                    }
                }
                char ownSt[128];
                formatStatusList(sr, ownSt, sizeof(ownSt));
                const bool ownClean = (ownSt[0] == '\0' || (ownSt[0] == '[' && ownSt[1] == ']'));
                if (ownClean && ctxSt[0] == '\0') continue;
                double dev = -1.0;
                int sameP = -1;
                double etol = 0.0;
                gp_Pnt p0;
                if (sh.ShapeType() == TopAbs_EDGE) {
                    const TopoDS_Edge ee = TopoDS::Edge(sh);
                    sameP = BRep_Tool::SameParameter(ee) ? 1 : 0;
                    etol = BRep_Tool::Tolerance(ee);
                    Standard_Real ef = 0, el = 0, pf = 0, pl = 0;
                    Handle(Geom_Curve) ec = BRep_Tool::Curve(ee, ef, el);
                    Handle(Geom2d_Curve) pc2 = BRep_Tool::CurveOnSurface(ee, f, pf, pl);
                    TopLoc_Location fl;
                    Handle(Geom_Surface) fs2 = BRep_Tool::Surface(f, fl);
                    if (!ec.IsNull()) p0 = ec->Value(ef);
                    if (!ec.IsNull() && !pc2.IsNull() && !fs2.IsNull()) {
                        dev = 0.0;
                        for (int k = 0; k <= 64; k++) {
                            const double t = ef + (el - ef) * (double)k / 64.0;
                            const gp_Pnt2d q = pc2->Value(t);
                            gp_Pnt b3 = fs2->Value(q.X(), q.Y());
                            if (!fl.IsIdentity()) b3.Transform(fl.Transformation());
                            dev = std::max(dev, ec->Value(t).Distance(b3));
                        }
                    }
                } else if (sh.ShapeType() == TopAbs_VERTEX) {
                    p0 = BRep_Tool::Pnt(TopoDS::Vertex(sh));
                    etol = BRep_Tool::Tolerance(TopoDS::Vertex(sh));
                }
                std::fprintf(stderr,
                             "DIAG_R2SUB rid=%d i=%d type=%s own=%s ctx=[%s] dev=%.6g sameP=%d "
                             "tol=%.6g p=(%.5f,%.5f,%.5f)\n",
                             rid, si, sh.ShapeType() == TopAbs_EDGE     ? "EDGE"
                             : sh.ShapeType() == TopAbs_VERTEX ? "VERTEX"
                             : sh.ShapeType() == TopAbs_WIRE   ? "WIRE"
                                                               : "FACE", ownSt, ctxSt, dev, sameP,
                             etol, p0.X(), p0.Y(), p0.Z());
            }
            dumpDiagWiresOfFace(rid, f, sewTol, "r2-invalid", &anc, &faceMap, &faceRid, &eprimeFill,
                                &rs);
        }
    } catch (const Standard_Failure&) {
        std::fprintf(stderr, "DIAG_R2PROBE threw\n");
    }
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
            if (d > BRep_Tool::Tolerance(e)) {
                fireTolRewriteEdge(B, e,
                                   std::min(d * 1.001 + Precision::Confusion(), cap), "other",
                                   cap);
            }
        }
    }
}

bool tryStageP(const MeshView& mv, RegionSet& rs, std::vector<TopoDS_Face>& out);

bool buildFaces(const MeshView& mv, RegionSet& rs, const std::vector<TopoDS_Vertex>& verts,
                std::vector<TopoDS_Face>& out, WarnFn warn) {
    out.clear();
    try {
        gCompNTri = (int)mv.nTri;
        gReconnectFires = 0;
        gOriNotDetermined = 0;
        lastCanonClearPass = -1;
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

        const double sewTol = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
        gRid100PrevValid = -1;
        dumpSliverRegionCensus(mv, rs, sewTol);
        dumpDiagTheta1136(mv, rs);

        // Stage P routing (RULE 5.1a): prismatic -> profiles/prisms/union;
        // decline or build failure falls through byte-identically.
        if (tryStageP(mv, rs, out)) {
            std::vector<ChainGeom> noGeom(rs.chains.size());
            std::vector<char> noCol(rs.chains.size(), 0);
            dumpBoundCredCensus(mv, rs, noGeom, noCol, sewTol);
            return true;
        }

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
        // E′ facet-fill: emit mesh triangles without marking exploded, so the
        // cascade U0/U1/U2 ladder (wouldSaturate) still sees the same analytic
        // set H1 did. exploded[] would jump to U2 + host R2 (D-S3-7 revert).
        std::vector<char> eprimeFill((size_t)std::max(0, maxId) + 1, 0);

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
        birthMeshEdgePlanePCurves(rs, meshE, edgeOk, mv);
        gSeamFailFirst = 0;
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
            birthMeshEdgePlanePCurves(rs, meshE, edgeOk, mv);
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
            int nMix = 0, nNone = 0, nFail = 0, nOk = 0, nChainSewFb = 0;
            (void)tryPlaneLoopCircles(rs, mv, sewTol);
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                collapsed[ci] = 0;
                geom[ci] = ChainGeom{};
                const BoundaryChain& ch = rs.chains[ci];
                const Region* A = regionById(rs, ch.regA);
                const Region* B = regionById(rs, ch.regB);
                if (A && regionExploded(exploded, A->id)) A = nullptr;
                if (B && regionExploded(exploded, B->id)) B = nullptr;
                auto admitAnalytic = [&](const Region* R) {
                    if (!R || !isAnalytic(R)) return false;
                    if (R->id >= 0 && (size_t)R->id < eprimeFill.size() && eprimeFill[(size_t)R->id])
                        return false;
                    return true;
                };
                gp_Ax2 plcAx;
                double plcR = 0.0;
                const bool plc = planeLoopCircleForChain((int)ci, plcAx, plcR) &&
                                 plcR > Precision::Confusion();
                const bool anaPair = admitAnalytic(A) && admitAnalytic(B);
                if (!anaPair) {
                    nMix++;
                    // mixed analytic|faceted, or island|island: polyline verbatim
                    // unless detector B tagged this chain as a plane-loop circle.
                    if ((admitAnalytic(A) && !admitAnalytic(B)) ||
                        (!admitAnalytic(A) && admitAnalytic(B))) {
                        const Region* an = admitAnalytic(A) ? A : B;
                        for (int eid : ch.meshEdges) {
                            if (eid < 0 || (size_t)eid >= mv.nEdge || !edgeOk[(size_t)eid]) continue;
                            const auto& ev = mv.compEdges[eid];
                            double d = mixedEdgeDeviation(mv, verts[(size_t)ev.first],
                                                          verts[(size_t)ev.second], an);
                            rs.stats.maxEdgeTol = std::max(rs.stats.maxEdgeTol, d);
                        }
                    }
                    // D-130: a plane-loop CIRCLE may only replace the polyline
                    // when the OTHER side of the chain is a region that can own
                    // the same edge. Against an island (regB < 0) or a region
                    // already demoted to facets, the circle is an edge nothing
                    // else references: the plate gets one CIRCLE, the facets keep
                    // the N mesh edges under it, and the shell opens by N + 1.
                    // Tier 2 (mesh polyline) is the correct answer there, and it
                    // is shared by construction.
                    const bool filledA =
                        ch.regA >= 0 && (size_t)ch.regA < eprimeFill.size() &&
                        eprimeFill[(size_t)ch.regA];
                    const bool filledB =
                        ch.regB >= 0 && (size_t)ch.regB < eprimeFill.size() &&
                        eprimeFill[(size_t)ch.regB];
                    const bool plcShareable =
                        plc && ch.regA >= 0 && ch.regB >= 0 && !filledA && !filledB &&
                        !regionExploded(exploded, ch.regA) &&
                        !regionExploded(exploded, ch.regB);
                    if (plc && !plcShareable && diagP2Enabled())
                        std::fprintf(stderr,
                                     "DIAG_130_PLCTIER2 ci=%d regA=%d regB=%d filledA=%d "
                                     "filledB=%d R=%.6f\n",
                                     (int)ci, ch.regA, ch.regB, filledA ? 1 : 0, filledB ? 1 : 0,
                                     plcR);
                    if (!plcShareable) continue;
                }
                // analytic | analytic (or mixed plane-loop CIRCLE)
                if (anaPair && recoverPass != lastCanonClearPass) {
                    gRegionSurf.clear();
                    gSurfIdent.clear();
                    gSeamTShapes.clear();
                    gMeshTShapes.clear();
                    gWirePopCi.clear();
                    gWirePopSite.clear();
                    for (const auto& me : meshE) {
                        if (!me.IsNull()) gMeshTShapes.insert(diagTShapePtr(me));
                    }
                    lastCanonClearPass = recoverPass;
                }
                const bool bothCyl = anaPair && A->type == SurfType::Cylinder &&
                                     B->type == SurfType::Cylinder;
                // Banned collapse: keep mesh polyline. No IntAna warning (AC #2)
                // and no chainEdgeFail enrollment (Body11).
                if (anaPair && (collapseBanned[ci] || (bothCyl && fallbackBanned[ci]))) {
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
                if (anaPair) {
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
                }
                // Detector B: a plane-loop CIRCLE answers where IntAna has no
                // answer at all. It never overrides a curve IntAna derived from
                // the two shipped surfaces -- a shared edge has to lie on both,
                // and two distinct planes meet along their Lin and nowhere else,
                // so an arc there is off-surface by construction (on
                // handle-pickup that swap took 42 plane|plane Lin chains, +42
                // STEP edges, +22 faces and +1.26 mm3 of bulge past main).
                // Circ/Elips were never overridden; Lin is now in the same class.
                if (plc && curve.kind == AnalyticCurve::Lin && diagP2Enabled())
                    std::fprintf(stderr,
                                 "DIAG_130_PLCKEEPLIN ci=%d regA=%d regB=%d anaPair=%d "
                                 "tA=%d tB=%d R=%.6f\n",
                                 (int)ci, ch.regA, ch.regB, anaPair ? 1 : 0,
                                 A ? (int)A->type : -1, B ? (int)B->type : -1, plcR);
                if (plc && curve.kind == AnalyticCurve::None) {
                    curve.kind = AnalyticCurve::Circ;
                    curve.circ = gp_Circ(plcAx, plcR);
                }
                // IntAna none: cyl|cyl is legal polyline (IA_CYLCYL_NOGEOM).
                // plane|plane and plane|cyl misses are post-fit edge failures.
                if (curve.kind == AnalyticCurve::None) {
                    nNone++;
                    if (anaPair && A && B &&
                        ((A->type == SurfType::Plane && B->type == SurfType::Plane) ||
                         (A->type == SurfType::Plane && B->type == SurfType::Cylinder) ||
                         (A->type == SurfType::Cylinder && B->type == SurfType::Plane)))
                        chainEdgeFail[ci] = 1;
                    continue;
                }

                // Terminals: open chain = first/last meshVert; closed = seam vertex (or first).
                int ia, ib;
                bool full = ch.closedLoop;
                if (ch.closedLoop) {
                    const Region* cyl = (A && A->type == SurfType::Cylinder) ? A
                                       : (B && B->type == SurfType::Cylinder) ? B
                                                                             : nullptr;
                    const Region* cn = (A && isChamferConeR(*A)) ? A
                                     : (B && isChamferConeR(*B)) ? B
                                                                : nullptr;
                    const Region* axR = cyl ? cyl : cn;
                    int sv = axR ? seamVertexOf(mv, *axR, ch) : ch.meshVerts.front();
                    ia = ib = sv;
                } else {
                    if (ch.meshVerts.size() < 2) continue;
                    ia = ch.meshVerts.front();
                    ib = ch.meshVerts.back();
                }
                if (ia < 0 || ib < 0 || (size_t)ia >= verts.size() || (size_t)ib >= verts.size())
                    continue;
                if (!full && curve.kind == AnalyticCurve::Circ &&
                    pntOf(mv, ia).Distance(pntOf(mv, ib)) <= Precision::Confusion()) {
                    geom[ci] = ChainGeom{};
                    nNone++;
                    continue;
                }
                const Region* cylR = (A && A->type == SurfType::Cylinder) ? A
                                   : (B && B->type == SurfType::Cylinder) ? B
                                                                         : nullptr;
                const Region* plnR = (A && A->type == SurfType::Plane) ? A
                                   : (B && B->type == SurfType::Plane) ? B
                                                                     : nullptr;
                const Region* coneR = (A && isChamferConeR(*A)) ? A
                                    : (B && isChamferConeR(*B)) ? B
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
                    } else if (coneR) {
                        double v = gp_Vec(coneR->ax.Location(), curve.circ.Location())
                                       .Dot(coneR->ax.Direction());
                        if (plnR && planePerpCylinder(*plnR, *coneR))
                            v = planeVOnCylinder(*plnR, *coneR);
                        curve.circ = coneIsoCircle(*coneR, v);
                    }
                }
                double snapCap = analyticSnapCap(mv, A, B);
                // ME_PLNPLN_*_PROJ: terminal-to-line residual (median 0.08 mm,
                // p90 0.45) exceeds sewTol. Bump to residual up to the
                // pickIntAna acceptance (1 mm) on plane|plane only — cyl
                // snaps stay on the I4 floor so extra MakeFace fails do not
                // eat the 127-cylinder small component.
                if (A && B && A->type == SurfType::Plane && B->type == SurfType::Plane) {
                    const double dA = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ia]));
                    const double dB = curveResidual(curve, BRep_Tool::Pnt(verts[(size_t)ib]));
                    const double sew = (mv.sewTol > 0.0) ? mv.sewTol : Precision::Confusion();
                    const double acceptCap = std::max(sew * 50.0, 1.0);
                    snapCap = std::max(snapCap,
                                       std::min(std::max(std::isfinite(dA) ? dA : 0.0,
                                                         std::isfinite(dB) ? dB : 0.0),
                                                acceptCap));
                }
                if (A && B && A->type == SurfType::Cylinder && B->type == SurfType::Cylinder &&
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
                geom[ci].curve = curve;
                geom[ci].ia = ia;
                geom[ci].ib = ib;
                geom[ci].regA = A ? A->id : ch.regA;
                geom[ci].regB = B ? B->id : ch.regB;
                geom[ci].chainSnapCap = snapCap;
                {
                    const double maxDist = chainMaxDistToOwnVerts(curve, mv, ch);
                    const double ratioSew = sewTol > 0.0 ? maxDist / sewTol : maxDist;
                    const char* kind = (curve.kind == AnalyticCurve::Circ)   ? "circ"
                                       : (curve.kind == AnalyticCurve::Lin) ? "lin"
                                                                            : "elips";
                    const int fb = (maxDist > sewTol) ? 1 : 0;
                    if (recoverPass == 0 && rounds == 0 && fallbackGuardPass == 0 &&
                        j6UncollapsePass == 0) {
                        if (fb) noteChainSewFallback((int)ci);
                        std::fprintf(stderr,
                                     "DIAG_CHAINSEW ci=%d kind=%s regA=%d regB=%d ratioSew=%.3f "
                                     "fallback=%d\n",
                                     (int)ci, kind, geom[ci].regA, geom[ci].regB, ratioSew, fb);
                    }
                    // Tagged plane-loop Circs are not exempt: Pratt circles that
                    // miss mesh verts by >> sewTol must demote to polyline like
                    // every other curve. IntAna plane|cyl circs on vertices have
                    // fb=0 (ratioSew≈0) and keep geom.
                    if (fb) {
                        geom[ci] = ChainGeom{};
                        nNone++;
                        nChainSewFb++;
                        continue;
                    }
                }
                continue;
            }
            std::unordered_set<int> seamLocked;
            gSeamFailFirst = 0;
            snapSeamVertices(rs, mv, geom, verts, seamLocked);
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                if (geom[ci].curve.kind == AnalyticCurve::None) continue;
                const int ia = geom[ci].ia, ib = geom[ci].ib;
                if (ia < 0 || ib < 0 || (size_t)ia >= verts.size() || (size_t)ib >= verts.size())
                    continue;
                if (!seamLocked.count(ia))
                    snapVertexToCurve(verts[(size_t)ia], geom[ci].curve, geom[ci].chainSnapCap);
                if (!seamLocked.count(ib))
                    snapVertexToCurve(verts[(size_t)ib], geom[ci].curve, geom[ci].chainSnapCap);
            }
            rebuildMeshEdges();
            birthMeshEdgePlanePCurves(rs, meshE, edgeOk, mv);
            // D-130 step 3: verts[] against its own snapshot, immediately before
            // the chain edges are made. A shared vertex that has left vsnap by
            // more than the sew budget is the thing that makes coincident chain
            // endpoints (and so the "identic" MakeEdge refusals).
            if (diagP2Enabled()) {
                int nOff = 0;
                double maxOff = 0.0, maxTol = 0.0;
                int worst = -1;
                for (size_t i = 0; i < verts.size(); i++) {
                    if (verts[i].IsNull()) continue;
                    const double d = BRep_Tool::Pnt(verts[i]).Distance(vsnap[i].p);
                    const double t = BRep_Tool::Tolerance(verts[i]);
                    if (d > std::max(sewTol, Precision::Confusion())) nOff++;
                    if (d > maxOff) {
                        maxOff = d;
                        worst = (int)i;
                    }
                    maxTol = std::max(maxTol, t);
                }
                std::fprintf(stderr,
                             "DIAG_130_VTXDRIFT pass=%d nVerts=%zu nOffSnapshot=%d maxOff=%.6f "
                             "worstVid=%d maxVtxTol=%.6f sewTol=%.6f\n",
                             recoverPass, verts.size(), nOff, maxOff, worst, maxTol, sewTol);
            }
            for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                if (geom[ci].curve.kind == AnalyticCurve::None) continue;
                const AnalyticCurve curve = geom[ci].curve;
                const int ia = geom[ci].ia;
                const int ib = geom[ci].ib;
                if (ia < 0 || ib < 0 || (size_t)ia >= verts.size() || (size_t)ib >= verts.size())
                    continue;
                const BoundaryChain& ch = rs.chains[ci];
                const Region* A = regionById(rs, geom[ci].regA);
                const Region* B = regionById(rs, geom[ci].regB);
                if (A && regionExploded(exploded, A->id)) A = nullptr;
                if (B && regionExploded(exploded, B->id)) B = nullptr;
                if (!isAnalytic(A) || !isAnalytic(B)) {
                    gp_Ax2 axt;
                    double Rt = 0.0;
                    if (curve.kind != AnalyticCurve::Circ ||
                        !planeLoopCircleForChain((int)ci, axt, Rt))
                        continue;
                    if (!isAnalytic(A) && !isAnalytic(B)) continue;
                }
                if (!A) A = B;
                if (!B) B = A;
                if (!A || !B) continue;
                const bool full = ch.closedLoop;
                WarnFn chainWarn =
                    (recoverPass == 0 && rounds == 0 && fallbackGuardPass == 0 &&
                     j6UncollapsePass == 0)
                        ? warn
                        : nullptr;
                TopoDS_Edge e;
                const Region* cylHint = (A && A->type == SurfType::Cylinder) ? A
                                       : (B && B->type == SurfType::Cylinder) ? B
                                                                             : nullptr;
                const Region* plnHint = (A && A->type == SurfType::Plane) ? A
                                       : (B && B->type == SurfType::Plane) ? B
                                                                         : nullptr;
                gp_Pnt mid = BRep_Tool::Pnt(verts[(size_t)ia]);
                const bool partialCylArcMid =
                    cylHint && !cylHint->closed360 && curve.kind == AnalyticCurve::Circ;
                if (partialCylArcMid && ch.meshVerts.size() < 3)
                    mid = partialCylCapArcMid(*cylHint, plnHint, curve.circ);
                else if ((!full && curve.kind == AnalyticCurve::Circ &&
                          ch.meshVerts.size() >= 3) ||
                         (partialCylArcMid && ch.meshVerts.size() >= 3)) {
                    mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
                    // D-130-14: (uMin+uMax)/2 is a PARTIAL cylinder's midpoint;
                    // on a closed360 region it reads pi, an azimuth the chain
                    // need not contain, and the vertex it returns is then an
                    // endpoint -- makeArc reads dm = 0 and takes the SHORT arc.
                    // A 288 deg rim chain came back as its own 72 deg
                    // complement. The chain's own middle vertex is the hint.
                    if (cylHint && !(unionBuildOn() && cylHint->closed360)) {
                        const double um = 0.5 * (cylHint->uMin + cylHint->uMax);
                        double best = 1e300;
                        for (int lv : ch.meshVerts) {
                            double du = std::fabs(regionU(*cylHint, pntOf(mv, lv)) - um);
                            if (du < best) {
                                best = du;
                                mid = pntOf(mv, lv);
                            }
                        }
                    }
                } else if (!full && curve.kind == AnalyticCurve::Circ &&
                           ch.meshVerts.size() == 2) {
                    const gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)ia]);
                    const gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)ib]);
                    double p1 = ElCLib::Parameter(curve.circ, pa);
                    double p2 = ElCLib::Parameter(curve.circ, pb);
                    double df = p2 - p1;
                    while (df <= 0.0) df += 2.0 * kPi;
                    if (df > kPi) df -= 2.0 * kPi;
                    mid = ElCLib::Value(p1 + 0.5 * df, curve.circ);
                } else if (!full && curve.kind == AnalyticCurve::Elips && ch.meshVerts.size() >= 3)
                    mid = pntOf(mv, ch.meshVerts[ch.meshVerts.size() / 2]);
                else if (!full && curve.kind == AnalyticCurve::Elips && ch.meshVerts.size() >= 2) {
                    gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)ia]);
                    gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)ib]);
                    mid = gp_Pnt(0.5 * (pa.X() + pb.X()), 0.5 * (pa.Y() + pb.Y()),
                                 0.5 * (pa.Z() + pb.Z()));
                }
                if (full && curve.kind == AnalyticCurve::Circ) {
                    e = makeFullCircle(curve.circ, verts[(size_t)ia]);
                    if (e.IsNull())
                        e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
                } else {
                    const bool partialCylArc =
                        cylHint && !cylHint->closed360 && curve.kind == AnalyticCurve::Circ;
                    if (partialCylArc || (!full && curve.kind == AnalyticCurve::Circ &&
                                          ch.meshVerts.size() >= 2)) {
                        e = makeArc(curve.circ, verts[(size_t)ia], verts[(size_t)ib], mid);
                        if (e.IsNull())
                            e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib],
                                                  full);
                    } else if (!full && curve.kind == AnalyticCurve::Elips &&
                               ch.meshVerts.size() >= 2) {
                        e = makeEllipseArc(curve.elips, verts[(size_t)ia], verts[(size_t)ib], mid);
                        if (e.IsNull())
                            e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib],
                                                  full);
                    } else {
                        e = makeEdgeFromCurve(curve, verts[(size_t)ia], verts[(size_t)ib], full);
                    }
                }
                (void)plnHint;
                if (e.IsNull()) {
                    const gp_Pnt paF = BRep_Tool::Pnt(verts[(size_t)ia]);
                    const gp_Pnt pbF = BRep_Tool::Pnt(verts[(size_t)ib]);
                    const bool identicLine =
                        curve.kind == AnalyticCurve::Lin &&
                        (ia == ib || paF.Distance(pbF) <= Precision::Confusion());
                    const bool identicCirc =
                        curve.kind == AnalyticCurve::Circ &&
                        (ia == ib || paF.Distance(pbF) <= Precision::Confusion());
                    if (!identicLine && !identicCirc) chainEdgeFail[ci] = 1;
                    if (!identicLine && !identicCirc) nFail++;
                    static thread_local int gMeFailN = 0;
                    // D-130 step 3: every MakeEdge refusal, not only the
                    // non-degenerate ones -- on this part all 30 were "identic"
                    // (coincident endpoints) and so invisible here.
                    if (gMeFailN < 64) {
                        const char* d130 = std::getenv("STL2STEP_DIAG_130");
                        if (d130 && d130[0] && d130[0] != '0') {
                        gMeFailN++;
                        const gp_Pnt pa = BRep_Tool::Pnt(verts[(size_t)ia]);
                        const gp_Pnt pb = BRep_Tool::Pnt(verts[(size_t)ib]);
                        double p1 = 0, p2 = 0;
                        if (curve.kind == AnalyticCurve::Circ) {
                            p1 = ElCLib::Parameter(curve.circ, pa);
                            p2 = ElCLib::Parameter(curve.circ, pb);
                        }
                        const Region* RA130 = regionById(rs, geom[ci].regA);
                        const Region* RB130 = regionById(rs, geom[ci].regB);
                        auto kindName130 = [](const Region* R) {
                            if (!R) return "none";
                            switch (R->type) {
                            case SurfType::Plane: return "plane";
                            case SurfType::Cylinder: return "cyl";
                            case SurfType::Cone: return "cone";
                            default: return "other";
                            }
                        };
                        std::fprintf(stderr,
                                     "DIAG_130_MEFAIL ci=%d kind=%d full=%d closedLoop=%d nV=%zu "
                                     "identic=%d sameVtxT=%d nMeshE=%zu vFront=%d vBack=%d "
                                     "R=%.6g dA=%.6g dB=%.6g dist=%.6g p1=%.6g p2=%.6g "
                                     "tolA=%.6g tolB=%.6g ia=%d ib=%d ridA=%d ridB=%d "
                                     "kindA=%s kindB=%s snapCap=%.6g "
                                     "vA=(%.5f,%.5f,%.5f) vB=(%.5f,%.5f,%.5f) "
                                     "mA=(%.5f,%.5f,%.5f) mB=(%.5f,%.5f,%.5f) meshLen=%.6g "
                                     "sA=(%.5f,%.5f,%.5f) sB=(%.5f,%.5f,%.5f) sTolA=%.6g\n",
                                     (int)ci, (int)curve.kind, (int)full,
                                     ch.closedLoop ? 1 : 0, ch.meshVerts.size(),
                                     (identicLine || identicCirc) ? 1 : 0,
                                     (diagTShapePtr(verts[(size_t)ia]) ==
                                      diagTShapePtr(verts[(size_t)ib]))
                                         ? 1
                                         : 0,
                                     ch.meshEdges.size(),
                                     ch.meshVerts.empty() ? -1 : ch.meshVerts.front(),
                                     ch.meshVerts.empty() ? -1 : ch.meshVerts.back(),
                                     curve.kind == AnalyticCurve::Circ ? curve.circ.Radius() : 0.0,
                                     curveResidual(curve, pa), curveResidual(curve, pb),
                                     pa.Distance(pb), p1, p2, BRep_Tool::Tolerance(verts[(size_t)ia]),
                                     BRep_Tool::Tolerance(verts[(size_t)ib]), ia, ib,
                                     geom[ci].regA, geom[ci].regB, kindName130(RA130),
                                     kindName130(RB130), geom[ci].chainSnapCap, pa.X(), pa.Y(),
                                     pa.Z(), pb.X(), pb.Y(), pb.Z(), pntOf(mv, ia).X(),
                                     pntOf(mv, ia).Y(), pntOf(mv, ia).Z(), pntOf(mv, ib).X(),
                                     pntOf(mv, ib).Y(), pntOf(mv, ib).Z(),
                                     pntOf(mv, ia).Distance(pntOf(mv, ib)),
                                     vsnap[(size_t)ia].p.X(), vsnap[(size_t)ia].p.Y(),
                                     vsnap[(size_t)ia].p.Z(), vsnap[(size_t)ib].p.X(),
                                     vsnap[(size_t)ib].p.Y(), vsnap[(size_t)ib].p.Z(),
                                     vsnap[(size_t)ia].t);
                        }
                    }
                    if ((diagP2Enabled() || diagPlatesEnabled()) && gSeamFailFirst == 0 &&
                        (seamLocked.count(ia) || seamLocked.count(ib))) {
                        gSeamFailFirst = 1;
                        const Region* RA = regionById(rs, geom[ci].regA);
                        const Region* RB = regionById(rs, geom[ci].regB);
                        std::fprintf(stderr,
                                     "DIAG_SEAMFAIL first=makeedge ci=%d ia=%d ib=%d "
                                     "ridA=%d ridB=%d classA=%s classB=%s kind=%d "
                                     "pA=(%.6f,%.6f,%.6f) pB=(%.6f,%.6f,%.6f)\n",
                                     (int)ci, ia, ib, geom[ci].regA, geom[ci].regB,
                                     RA ? (RA->type == SurfType::Plane ? "plane" : "cyl") : "-",
                                     RB ? (RB->type == SurfType::Plane ? "plane" : "cyl") : "-",
                                     (int)curve.kind, BRep_Tool::Pnt(verts[(size_t)ia]).X(),
                                     BRep_Tool::Pnt(verts[(size_t)ia]).Y(),
                                     BRep_Tool::Pnt(verts[(size_t)ia]).Z(),
                                     BRep_Tool::Pnt(verts[(size_t)ib]).X(),
                                     BRep_Tool::Pnt(verts[(size_t)ib]).Y(),
                                     BRep_Tool::Pnt(verts[(size_t)ib]).Z());
                    }
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
                // D-130-16: the two-arc rim collapse. The chain publishes TWO
                // arcs of its own circle, meeting at its own mesh vertex
                // nearest the seam target, when it is the rim of a closed360
                // cylinder that carries an inner wire and the single arc would
                // otherwise leave the seam nowhere outside that wire to stand.
                std::vector<TopoDS_Edge> chainEdges{e};
                if (unionBuildOn() && !full && curve.kind == AnalyticCurve::Circ) {
                    const int lvSplit = rimSplitVertex(rs, mv, ch, (int)ci, A, B,
                                                       geom[ci].chainSnapCap, curve);
                    if (lvSplit >= 0 && (size_t)lvSplit < verts.size() &&
                        !verts[(size_t)lvSplit].IsNull()) {
                        size_t k = 0;
                        while (k < ch.meshVerts.size() && ch.meshVerts[k] != lvSplit) k++;
                        auto halfMid = [&](size_t i0, size_t i1) -> gp_Pnt {
                            if (i1 - i0 >= 2) return pntOf(mv, ch.meshVerts[(i0 + i1) / 2]);
                            const double p1 =
                                ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[i0]));
                            const double p2 =
                                ElCLib::Parameter(curve.circ, pntOf(mv, ch.meshVerts[i1]));
                            double df = p2 - p1;
                            while (df <= 0.0) df += 2.0 * kPi;
                            if (df > kPi) df -= 2.0 * kPi;
                            return ElCLib::Value(p1 + 0.5 * df, curve.circ);
                        };
                        if (k > 0 && k + 1 < ch.meshVerts.size()) {
                            snapVertexToCurve(verts[(size_t)lvSplit], curve,
                                              geom[ci].chainSnapCap);
                            const TopoDS_Edge e0 =
                                makeArc(curve.circ, verts[(size_t)ia], verts[(size_t)lvSplit],
                                        halfMid(0, k));
                            const TopoDS_Edge e1 =
                                makeArc(curve.circ, verts[(size_t)lvSplit], verts[(size_t)ib],
                                        halfMid(k, ch.meshVerts.size() - 1));
                            if (!e0.IsNull() && !e1.IsNull()) chainEdges = {e0, e1};
                        }
                    }
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
                    for (TopoDS_Edge& ce : chainEdges) {
                        bindAllVariants(ce, *ownerR, kindStr, mv, (int)ci, acc);
                        bindAllVariants(ce, *consumerR, kindStr, mv, (int)ci, acc);
                        bindAllVariants(ce, *ownerR, kindStr, mv, (int)ci, acc);
                    }
                    if (diagP2Enabled()) {
                        std::fprintf(stderr,
                                     "DIAG_SEAMBIND ci=%d nWrite=%d nGuardHit=%d nDoubleWrite=%d\n",
                                     (int)ci, acc.nWrite, acc.nGuardHit, acc.nDoubleWrite);
                    }
                }
                geom[ci].collapsed = true;
                geom[ci].edges = chainEdges;
                for (const TopoDS_Edge& ce : chainEdges) noteWirePop(ce, "rebuildCollapsed", (int)ci);
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
            // Mesh pcurve birth is bindAllVariants (collapsed) + makePCurveOnSurf
            // interned-plane identity. bindPlateMeshPcurves UpdateEdge on shared
            // polylines reverts the component (D-S3-37: fix at birth, not rebind).
            (void)meshE;
            (void)edgeOk;
            (void)nChainSewFb;
            if (recoverPass == 0 && rounds == 0 && fallbackGuardPass == 0 && j6UncollapsePass == 0) {
                dumpChainFitCensus(mv, rs, geom, collapsed, sewTol);
                dumpBoundCredCensus(mv, rs, geom, collapsed, sewTol);
                dumpChainLineCensus(mv, rs, geom);
                std::fprintf(stderr,
                             "DIAG_CHAINSEW_SUM part=%s nFallback=%d nOk=%d nNone=%d total=%zu "
                             "sewTol=%.6f\n",
                             sliverCensusPartName((int)mv.nTri), nChainSewFb, nOk, nNone,
                             rs.chains.size(), sewTol);
            }
            if (diagCollapse) {
                std::fprintf(stderr,
                             "DIAG_COLLAPSE mix=%d none=%d fail=%d ok=%d total=%zu recover=%d "
                             "rounds=%d\n",
                             nMix, nNone, nFail, nOk, rs.chains.size(), recoverPass, rounds);
                int nMixed = 0, nMixedCollapsed = 0;
                for (size_t ci = 0; ci < rs.chains.size(); ci++) {
                    const BoundaryChain& ch = rs.chains[ci];
                    const Region* A = regionById(rs, ch.regA);
                    const Region* B = regionById(rs, ch.regB);
                    if (A && regionExploded(exploded, A->id)) A = nullptr;
                    if (B && regionExploded(exploded, B->id)) B = nullptr;
                    const int anaA = (A && isAnalytic(A)) ? 1 : 0;
                    const int anaB = (B && isAnalytic(B)) ? 1 : 0;
                    const int coll = ((size_t)ci < collapsed.size() && collapsed[ci]) ? 1 : 0;
                    std::fprintf(stderr,
                                 "DIAG_COLLAPSE_CI ci=%zu regA=%d regB=%d anaA=%d anaB=%d "
                                 "collapsed=%d\n",
                                 ci, ch.regA, ch.regB, anaA, anaB, coll);
                    if (anaA != anaB) nMixed++;
                    if (anaA != anaB && coll) nMixedCollapsed++;
                }
                std::fprintf(stderr,
                             "DIAG_COLLAPSE_SUM part=%s nChains=%zu nMixed=%d nMixedCollapsed=%d\n",
                             sliverCensusPartName((int)mv.nTri), rs.chains.size(), nMixed,
                             nMixedCollapsed);
            }
        };

        auto buildOneRegion = [&](Region& r, std::vector<TopoDS_Face>& acc) -> bool {
            if (failRidHits(r.id)) {
                r.reject = Reject::FaceBuildFailed;
                r.builtAs = BuiltAs::NotBuilt;
                return false;
            }
            if (regionExploded(exploded, r.id)) return true;
            if (r.id >= 0 && (size_t)r.id < eprimeFill.size() && eprimeFill[(size_t)r.id]) {
                r.builtAs = BuiltAs::ExplodedToFacets;
                return true;
            }
            r.builtAs = BuiltAs::NotBuilt;
            if (r.type == SurfType::Cone || r.type == SurfType::Sphere || r.type == SurfType::Torus) {
                // U-3: Sphere/Torus stay NYI. Other cones stay ConeNYI.
                // ChamferCone: Geom_ConicalSurface from axis + two radii
                // (radius / v range), MakeFace like cylinder Single.
                if (r.type == SurfType::Cone && r.origin == Origin::ChamferCone) {
                    // Diag only: name the step that refuses a chamfer frustum.
                    auto coneFail = [&](const char* why, double a, double b) {
                        if (diagP2Enabled() || diag130Enabled())
                            std::fprintf(stderr,
                                         "DIAG_CONEFACE rid=%d R=%.4f nTri=%zu why=%s "
                                         "a=%.6g b=%.6g\n",
                                         r.id, r.radius, r.tris.size(), why, a, b);
                        return false;
                    };
                    const Loop *capL = nullptr, *capH = nullptr;
                    std::vector<TopoDS_Wire> inners;
                    for (const Loop& lp : r.loops) {
                        if (lp.role == LoopRole::CapLow) capL = &lp;
                        else if (lp.role == LoopRole::CapHigh) capH = &lp;
                        else if (lp.role == LoopRole::Inner) {
                            TopoDS_Wire iw;
                            if (!buildLoopWire(iw, lp, rs, mv, geom, collapsed, meshE, edgeOk,
                                               nullptr))
                                return coneFail("inner-loop-wire", 0, 0);
                            inners.push_back(iw);
                        }
                    }
                    if (!capL || !capH || capL->chainIdx.empty() || capH->chainIdx.empty())
                        return coneFail("missing-cap-loop", capL ? (double)capL->chainIdx.size() : -1,
                                        capH ? (double)capH->chainIdx.size() : -1);

                    auto circFromLoop = [&](const Loop& lp) -> gp_Circ {
                        double sumV = 0.0, sumR = 0.0;
                        int n = 0;
                        for (int cix : lp.chainIdx) {
                            if (cix < 0 || (size_t)cix >= rs.chains.size()) continue;
                            for (int lv : rs.chains[(size_t)cix].meshVerts) {
                                const gp_Pnt p = pntOf(mv, lv);
                                gp_Vec d(r.ax.Location(), p);
                                const double v = d.Dot(r.ax.Direction());
                                gp_Vec rad = d - gp_Vec(r.ax.Direction()) * v;
                                const double m = rad.Magnitude();
                                if (!(m > Precision::Confusion())) continue;
                                sumV += v;
                                sumR += m;
                                n++;
                            }
                        }
                        const double vm = n ? sumV / (double)n : 0.0;
                        const double Rm = n ? sumR / (double)n : r.radius;
                        const gp_Pnt loc =
                            r.ax.Location().Translated(gp_Vec(r.ax.Direction()) * vm);
                        return gp_Circ(
                            gp_Ax2(loc, r.ax.Direction(), r.ax.XDirection()), Rm);
                    };

                    // ONE surface for the face and the bind site: the region's
                    // own SurfVar::ConeBase. The rim circles the neighbouring
                    // plane and cylinder share are coneIsoCircle() of the same
                    // region fields, so they lie on this surface exactly -- which
                    // is what lets exactMaxAtBind certify them circle-on-cone
                    // rather than record a gap. Re-deriving the surface from the
                    // mesh loops here would put the shared rim off the face.
                    Handle(Geom_ConicalSurface) csurf =
                        Handle(Geom_ConicalSurface)::DownCast(regionSurf(r, SurfVar::ConeBase));
                    if (csurf.IsNull())
                        return coneFail("regionSurf-ConeBase-null", r.radius, chamferRhiOf(r));
                    const double Rlo = r.radius;
                    const double Rhi = chamferRhiOf(r);
                    const double h = chamferHeightOf(r);
                    const double dR = Rhi - Rlo;
                    if (!(std::fabs(h) > Precision::Confusion()) ||
                        !(Rlo > Precision::Confusion()) || dR <= Precision::Confusion())
                        return coneFail("degenerate-frustum", h, dR);
                    const double ang = csurf->SemiAngle();
                    // Which cap loop carries the R_lo rim is decided by
                    // measurement, not by LoopRole: the roles are assigned from a
                    // v ordering a cone region does not carry (vMin/vMax hold the
                    // two RADII).
                    const gp_Circ mL = circFromLoop(*capL);
                    const gp_Circ mH = circFromLoop(*capH);
                    const bool loIsCapL =
                        std::fabs(mL.Radius() - Rlo) <= std::fabs(mH.Radius() - Rlo);
                    const Loop* loopLo = loIsCapL ? capL : capH;
                    const Loop* loopHi = loIsCapL ? capH : capL;
                    const gp_Circ circLo = coneIsoCircle(r, 0.0);
                    const gp_Circ circHi = coneIsoCircle(r, h);

                    int vL = vertexClosestToUOnLoop(mv, r, *loopLo, rs, 0.0);
                    int vH = vertexClosestToUOnLoop(mv, r, *loopHi, rs, 0.0);
                    if (vL < 0 || vH < 0 || (size_t)vL >= verts.size() ||
                        (size_t)vH >= verts.size())
                        return coneFail("no-seam-vertex", (double)vL, (double)vH);
                    AnalyticCurve acL, acH;
                    acL.kind = AnalyticCurve::Circ;
                    acL.circ = circLo;
                    acH.kind = AnalyticCurve::Circ;
                    acH.circ = circHi;
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

                    // The seam is the frustum's GENERATOR, given explicitly, the
                    // way trySeamed360 gives the cylinder's: an edge built from
                    // two vertices alone carries whatever line those two mesh
                    // points happen to define, and its 3d length then disagrees
                    // with the surface's own unit-speed v by the mesh residual
                    // (measured 8.4e-6 here), which is a UV gap at the rim.
                    // Snapping both seam vertices to that generator is the same
                    // snapVertexToCurve the cylinder seam uses.
                    TopoDS_Edge eSeam;
                    try {
                        const gp_Pnt pSeamLo =
                            circLo.Location().Translated(gp_Vec(circLo.XAxis().Direction()) * Rlo);
                        // From the R_lo rim TOWARD the R_hi rim, built from the
                        // frustum's own two numbers so it is right for either
                        // taper sign. Writing it as sin(Ang)*X + cos(Ang)*Z is
                        // the surface's +v direction, which points AWAY from the
                        // R_hi rim when Ang < 0 (h < 0) and gives a seam whose
                        // parameter range runs backwards.
                        const gp_Dir genDir(gp_Vec(circLo.XAxis().Direction()) * dR +
                                            gp_Vec(r.ax.Direction()) * h);
                        const gp_Lin gen(pSeamLo, genDir);
                        AnalyticCurve acs;
                        acs.kind = AnalyticCurve::Lin;
                        acs.lin = gen;
                        snapVertexToCurve(verts[(size_t)vL], acs, snapCap);
                        snapVertexToCurve(verts[(size_t)vH], acs, snapCap);
                        BRepBuilderAPI_MakeEdge ms(gen, verts[(size_t)vL], verts[(size_t)vH]);
                        if (!ms.IsDone()) return coneFail("seam-MakeEdge-failed", 0, 0);
                        eSeam = ms.Edge();
                    } catch (const Standard_Failure&) {
                        return coneFail("seam-MakeEdge-threw", 0, 0);
                    }

                    TopoDS_Edge eL, eH;
                    const int ciL0 = loopLo->chainIdx.front();
                    const int ciH0 = loopHi->chainIdx.front();
                    const bool simple =
                        loopLo->chainIdx.size() == 1 && loopHi->chainIdx.size() == 1 &&
                        takeFullCap(ciL0, acL.circ, verts[(size_t)vL], eL) &&
                        takeFullCap(ciH0, acH.circ, verts[(size_t)vH], eH);

                    // The cone's own pcurves, written where the cylinder path
                    // writes its iso pcurves (trySeamed360::bindIsoPCurves) and
                    // for the same reason: a rim is u-iso and the seam is u-const,
                    // and the seam needs TWO -- one at u0 and one at u0 + 2pi --
                    // which no single-pcurve projection can supply. OCCT's cone is
                    // unit speed in v, so a generator's 3d parameter and its v
                    // agree and the pcurve is a straight 2d line. The tolerance
                    // those edges carry is written by the certified circle-on-cone
                    // / line-on-cone classes at bindAllVariants, not here.
                    auto bindConePCurves = [&](TopoDS_Edge& eCap, double vIso, double& u0Out,
                                               bool writeU0) {
                        BRep_Builder B;
                        Standard_Real cf = 0, cl = 0;
                        Handle(Geom_Curve) cc = BRep_Tool::Curve(eCap, cf, cl);
                        if (cc.IsNull()) cf = 0.0;
                        if (writeU0) u0Out = cf;
                        Handle(Geom2d_Line) pc = new Geom2d_Line(gp_Pnt2d(cf, vIso),
                                                                 gp_Dir2d(1.0, 0.0));
                        B.UpdateEdge(eCap, pc, csurf, TopLoc_Location(), sewTol);
                    };

                    TopoDS_Face f;
                    bool got = false;
                    if (simple) {
                        double u0 = 0.0;
                        bindConePCurves(eL, 0.0, u0, true);
                        double uH = 0.0;
                        bindConePCurves(eH, coneVAtRadius(csurf->Cone(), Rhi), uH, false);
                        {
                            BRep_Builder Bs;
                            Standard_Real fs = 0, ls = 0;
                            Handle(Geom_Curve) cs = BRep_Tool::Curve(eSeam, fs, ls);
                            // The seam edge is parametrised by arc length from
                            // the R_lo rim toward the R_hi rim (genDir above), so
                            // its 3d parameter increases while the surface's v
                            // DECREASES whenever the frustum sits at negative v
                            // (Ang < 0). vSgn is that sign, taken from the
                            // surface itself, not assumed.
                            const double vHiParam = coneVAtRadius(csurf->Cone(), Rhi);
                            const double vSgn = (vHiParam >= 0.0) ? 1.0 : -1.0;
                            Handle(Geom2d_Line) pcS0 =
                                new Geom2d_Line(gp_Pnt2d(u0, 0.0), gp_Dir2d(0.0, vSgn));
                            Handle(Geom2d_Line) pcS1 = new Geom2d_Line(
                                gp_Pnt2d(u0 + 2.0 * kPi, 0.0), gp_Dir2d(0.0, vSgn));
                            Bs.UpdateEdge(eSeam, pcS0, pcS1, csurf, TopLoc_Location(), sewTol);
                            if (!cs.IsNull()) Bs.Range(eSeam, fs, ls);
                            eSeam.Closed(Standard_False);
                        }
                        BRep_Builder B;
                        TopoDS_Wire w;
                        B.MakeWire(w);
                        B.Add(w, eSeam);
                        B.Add(w, eH);
                        B.Add(w, TopoDS::Edge(eSeam.Reversed()));
                        B.Add(w, TopoDS::Edge(eL.Reversed()));
                        w.Closed(Standard_True);
                        got = makeFaceCopy(csurf, w, inners, r.outwardNormal, f) ||
                              makeFaceKeep(csurf, w, inners, r.outwardNormal, f);
                        if (got && ciL0 >= 0 && ciH0 >= 0 && (size_t)ciL0 < geom.size() &&
                            (size_t)ciH0 < geom.size()) {
                            geom[(size_t)ciL0].collapsed = true;
                            geom[(size_t)ciL0].edges = {eL};
                            collapsed[(size_t)ciL0] = 1;
                            geom[(size_t)ciH0].collapsed = true;
                            geom[(size_t)ciH0].edges = {eH};
                            collapsed[(size_t)ciH0] = 1;
                        }
                    }
                    if (!got) {
                        TopoDS_Wire wL, wH;
                        if (buildLoopWire(wL, *capL, rs, mv, geom, collapsed, meshE, edgeOk,
                                          nullptr) &&
                            buildLoopWire(wH, *capH, rs, mv, geom, collapsed, meshE, edgeOk,
                                          nullptr)) {
                            BRep_Builder Bw;
                            TopoDS_Wire ow;
                            Bw.MakeWire(ow);
                            Bw.Add(ow, eSeam);
                            for (BRepTools_WireExplorer ex(wH); ex.More(); ex.Next())
                                Bw.Add(ow, ex.Current());
                            Bw.Add(ow, TopoDS::Edge(eSeam.Reversed()));
                            std::vector<TopoDS_Edge> pathL;
                            for (BRepTools_WireExplorer ex(wL); ex.More(); ex.Next())
                                pathL.push_back(TopoDS::Edge(ex.Current()));
                            for (int i = (int)pathL.size() - 1; i >= 0; --i)
                                Bw.Add(ow, TopoDS::Edge(pathL[(size_t)i].Reversed()));
                            ow.Closed(Standard_True);
                            got = makeFaceCopy(csurf, ow, inners, r.outwardNormal, f) ||
                                  makeFaceKeep(csurf, ow, inners, r.outwardNormal, f);
                        }
                    }
                    if (!got || f.IsNull()) return coneFail("makeFace-failed", simple ? 1 : 0, 0);
                    // Same finish the cylinder Seamed360 face gets
                    // (trySeamed360::finishFace): pcurves completed on the face,
                    // then -- only if it is still invalid -- ShapeFix_Face for the
                    // SEAM and the ORIENTATION and nothing else. Every
                    // vertex-moving mode is pinned off: 130-CORE e7c00a8 measured
                    // that ShapeFix_Wire displaces the shared verts[] and opens
                    // the shell it was called to close (12.4 mm on this part).
                    addPcurvesOnFace(f, sewTol, true);
                    if (!faceIsValid(f)) {
                        try {
                            ShapeFix_Face sff(f);
                            sff.FixMissingSeamMode() = 1;
                            sff.FixAddNaturalBoundMode() = 0;
                            sff.FixOrientationMode() = 1;
                            sff.FixWireMode() = 0;
                            sff.FixSmallAreaWireMode() = 0;
                            sff.FixIntersectingWiresMode() = 0;
                            sff.FixLoopWiresMode() = 0;
                            sff.FixSplitFaceMode() = 0;
                            sff.Perform();
                            TopoDS_Shape res = sff.Result();
                            if (res.IsNull()) res = sff.Face();
                            int nF = 0;
                            TopoDS_Face g2;
                            for (TopExp_Explorer fx(res, TopAbs_FACE); fx.More(); fx.Next()) {
                                nF++;
                                g2 = TopoDS::Face(fx.Current());
                            }
                            if (nF == 1) f = g2;
                        } catch (const Standard_Failure&) {
                        }
                        setFaceOutward(f, r.outwardNormal);
                        addPcurvesOnFace(f, sewTol, true);
                    }
                    if (diagP2Enabled() || diag130Enabled()) {
                        int ie = 0;
                        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next(), ++ie) {
                            const TopoDS_Edge ew = TopoDS::Edge(ex.Current());
                            Standard_Real a = 0, b = 0;
                            Handle(Geom2d_Curve) pc =
                                BRep_Tool::CurveOnSurface(ew, csurf, TopLoc_Location(), a, b);
                            const bool sm = BRep_Tool::IsClosed(ew, f);
                            if (pc.IsNull()) {
                                std::fprintf(stderr,
                                             "DIAG_CONEUV rid=%d ie=%d ori=%s seam=%d pc=null\n",
                                             r.id, ie,
                                             ew.Orientation() == TopAbs_FORWARD ? "F" : "R",
                                             sm ? 1 : 0);
                                continue;
                            }
                            const gp_Pnt2d p0 = pc->Value(a), p1 = pc->Value(b);
                            std::fprintf(stderr,
                                         "DIAG_CONEUV rid=%d ie=%d ori=%s seam=%d rng=[%.6f,%.6f] "
                                         "uv0=(%.6f,%.6f) uv1=(%.6f,%.6f)\n",
                                         r.id, ie, ew.Orientation() == TopAbs_FORWARD ? "F" : "R",
                                         sm ? 1 : 0, a, b, p0.X(), p0.Y(), p1.X(), p1.Y());
                            TopoDS_Vertex va, vb;
                            TopExp::Vertices(ew, va, vb, Standard_False);
                            const gp_Pnt pa = va.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(va);
                            const gp_Pnt pb = vb.IsNull() ? gp_Pnt() : BRep_Tool::Pnt(vb);
                            std::fprintf(stderr,
                                         "DIAG_CONEUVV rid=%d ie=%d sameV=%d "
                                         "pa=(%.5f,%.5f,%.5f) pb=(%.5f,%.5f,%.5f) tol=%.3g\n",
                                         r.id, ie, va.IsSame(vb) ? 1 : 0, pa.X(), pa.Y(), pa.Z(),
                                         pb.X(), pb.Y(), pb.Z(), BRep_Tool::Tolerance(ew));
                        }
                    }
                    const bool fv = ensureFaceValid(f, meshTolCap(mv, &r));
                    if (diagP2Enabled() || diag130Enabled())
                        std::fprintf(stderr,
                                     "DIAG_CONEFACE rid=%d R=%.4f nTri=%zu why=BUILT simple=%d "
                                     "valid=%d Rlo=%.6f Rhi=%.6f ang=%.9f h=%.6f\n",
                                     r.id, r.radius, r.tris.size(), simple ? 1 : 0, fv ? 1 : 0,
                                     Rlo, Rhi, ang, h);
                    r.builtAs = BuiltAs::Seamed360;
                    acc.push_back(f);
                    return true;
                }
                r.reject = (r.type == SurfType::Cone)     ? Reject::ConeNYI
                           : (r.type == SurfType::Sphere) ? Reject::SphereNYI
                                                          : Reject::TorusNYI;
                return false;
            }
            if (r.type == SurfType::Cylinder && !isAnalytic(&r)) {
                if (r.id >= 0 && (size_t)r.id < eprimeFill.size()) eprimeFill[(size_t)r.id] = 1;
                r.builtAs = BuiltAs::ExplodedToFacets;
                return true;
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
                if (diagP2Enabled()) {
                    const Loop* outerL = nullptr;
                    for (const Loop& lp : r.loops)
                        if (lp.role == LoopRole::Outer) outerL = &lp;
                    dumpDiagCylFail(r, outerL, "seam", "seamed360-fail", TopoDS_Wire(),
                                    regionSurf(r, SurfVar::CylBase), TopoDS_Face());
                }
                return false;
            }
            if (r.type == SurfType::Cylinder) {
                TopoDS_Face f;
                if (!buildPartialCylinder(r, rs, mv, geom, collapsed, meshE, edgeOk, sewTol, f))
                    return false;
                r.builtAs = BuiltAs::Single;
                acc.push_back(f);
                if (r.id == 100) diagRid100Status("construction", r.id, f);
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
        int closureHealIter = 0;
        int closureHealTotal = 0;
        int closureHealAdmitted0 = 0;
        bool closureHealCapSet = false;
        int closureHealPendingRid = -1;
        int closureHealPendingFree = -1;
        bool closureHealStop = false;
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
                bool efill = (rid >= 0 && (size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid]);
                if (iid < 0 && !exp && !efill) continue;
                TopoDS_Face f = makeFacet(mv, verts, meshE, edgeOk, k);
                if (!f.IsNull()) {
                    ensureFaceValid(f, meshTolCap(mv, nullptr));
                    if (collapseDiagEnabled()) {
                        BRepCheck_Analyzer fa(f, Standard_True);
                        if (!fa.IsValid()) {
                            const int gt = mv.compTris ? mv.compTris[k] : -1;
                            gp_XYZ a0, b0, c0;
                            if (gt >= 0 && mv.tris && mv.pts) {
                                a0 = mv.pts[mv.tris[gt][0]];
                                b0 = mv.pts[mv.tris[gt][1]];
                                c0 = mv.pts[mv.tris[gt][2]];
                            }
                            const double ar = 0.5 * (b0 - a0).Crossed(c0 - a0).Modulus();
                            std::fprintf(stderr,
                                         "DIAG_FACETBAD k=%zu gt=%d rid=%d iid=%d exp=%d "
                                         "efill=%d area=%.6g A=(%.5f,%.5f,%.5f) "
                                         "B=(%.5f,%.5f,%.5f) C=(%.5f,%.5f,%.5f)\n",
                                         k, gt, rid, iid, exp ? 1 : 0, efill ? 1 : 0, ar, a0.X(),
                                         a0.Y(), a0.Z(), b0.X(), b0.Y(), b0.Z(), c0.X(), c0.Y(),
                                         c0.Z());
                        }
                    }
                    built.push_back(f);
                    builtRid.push_back((exp || efill) ? rid : -1);
                }
            }
        }

        if (built.empty()) return false;

        auto diagStageInvalid = [&](const char* stage) {
            if (!collapseDiagEnabled()) return;
            int nBad = 0;
            int firstI = -1;
            for (size_t i = 0; i < built.size(); i++) {
                if (built[i].IsNull()) continue;
                BRepCheck_Analyzer fa(built[i], Standard_True);
                if (fa.IsValid()) continue;
                if (firstI < 0) firstI = (int)i;
                nBad++;
            }
            std::fprintf(stderr, "DIAG_STAGEBAD stage=%s nBad=%d firstI=%d rid=%d\n", stage, nBad,
                         firstI,
                         (firstI >= 0 && (size_t)firstI < builtRid.size()) ? builtRid[(size_t)firstI]
                                                                          : -2);
        };
        diagStageInvalid("built");

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

        diagStageInvalid("prewalk");
        {
            orientFaceWalk(built);
        }
        diagStageInvalid("postwalk");

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

        diagStageInvalid("preshell");
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
                    stampTolWriter(eW, "other");
                }
            }
        }
        diagStageInvalid("postseampc");
        // D-S3-54: construction-time closure healing (pre-J6). Demote built faces
        // that own free edges via eprimeFill (never exploded[]). Full rebuild
        // with eprimeFill-aware chain admission; no uncollapse (opens neighbours).
        if (!closureHealStop && !unstable && recoverPass == 0 && j6UncollapsePass == 0 &&
            fallbackGuardPass == 0) {
            if (!closureHealCapSet) {
                closureHealAdmitted0 = countAdmittedRegions(rs, eprimeFill);
                closureHealCapSet = true;
            }
            const int freeBeforeHeal = countShellFreeEdges(sh);
            // D-130: the free-edge census that actually decides the revert is the
            // one the healer sees, not the post-explode DIAG_PREJ4 dump. Emit it
            // once, on the first heal decision of the run.
            if (diagP2Enabled() && freeBeforeHeal > 0 && closureHealIter == 0 &&
                closureHealPendingRid < 0) {
                std::fprintf(stderr, "DIAG_HEAL_CENSUS stage=preheal freeEdges=%d\n",
                             freeBeforeHeal);
                dumpDiagFreeEdges(sh, built, builtRid, eprimeFill, rs, mv, meshE, geom, collapsed,
                                  sewTol);
            }
            if (closureHealPendingRid >= 0) {
                if (freeBeforeHeal >= closureHealPendingFree) {
                    if ((size_t)closureHealPendingRid < eprimeFill.size())
                        eprimeFill[(size_t)closureHealPendingRid] = 0;
                    if (diagP2Enabled())
                        std::fprintf(stderr,
                                     "DIAG_HEAL stall rid=%d free=%d reverted=1\n",
                                     closureHealPendingRid, freeBeforeHeal);
                    closureHealPendingRid = -1;
                    closureHealPendingFree = -1;
                    closureHealStop = true;
                    goto try_rebuild;
                }
                closureHealPendingRid = -1;
                closureHealPendingFree = -1;
            }
            if (freeBeforeHeal > 0 && closureHealIter < closureHealAdmitted0) {
                const int admittedBefore = countAdmittedRegions(rs, eprimeFill);
                auto demotions = collectClosureHealDemotions(
                    sh, built, builtRid, eprimeFill, rs, mv, meshE, geom, collapsed, exploded,
                    sewTol);
                if (!demotions.empty()) {
                    std::sort(demotions.begin(), demotions.end(),
                              [&](const ClosureHealAttrib& a, const ClosureHealAttrib& b) {
                                  const Region* ra = regionById(rs, a.rid);
                                  const Region* rb = regionById(rs, b.rid);
                                  const size_t na = ra ? ra->tris.size() : 0;
                                  const size_t nb = rb ? rb->tris.size() : 0;
                                  if (na != nb) return na < nb;
                                  return a.rid < b.rid;
                              });
                    const ClosureHealAttrib& h = demotions.front();
                    if (h.rid >= 0 && (size_t)h.rid < eprimeFill.size() && !eprimeFill[(size_t)h.rid] &&
                        regionClosureHealEligible(h.rid, rs, eprimeFill, exploded)) {
                        eprimeFill[(size_t)h.rid] = 1;
                        if (Region* rr = regionByIdMut(rs, h.rid))
                            rr->builtAs = BuiltAs::ExplodedToFacets;
                        if (diagP2Enabled()) {
                            const Region* r = regionById(rs, h.rid);
                            const char* kind =
                                (r && r->type == SurfType::Plane) ? "plane" : "cyl";
                            const double R =
                                (r && r->type == SurfType::Cylinder) ? r->radius : 0.0;
                            const size_t nTri = r ? r->tris.size() : 0;
                            std::fprintf(stderr,
                                         "DIAG_HEAL iter=%d rid=%d kind=%s R=%.4f nTri=%zu "
                                         "attributedFrom=%s freeEdgesBefore=%d\n",
                                         closureHealIter, h.rid, kind, R, nTri,
                                         h.from ? h.from : "?", freeBeforeHeal);
                        }
                        const int admittedAfter = countAdmittedRegions(rs, eprimeFill);
                        if (admittedAfter < admittedBefore) {
                            closureHealTotal++;
                            closureHealIter++;
                            closureHealPendingRid = h.rid;
                            closureHealPendingFree = freeBeforeHeal;
                            goto try_rebuild;
                        }
                    }
                }
            }
            if (diagP2Enabled()) {
                std::fprintf(stderr,
                             "DIAG_HEAL_SUM part=%s iters=%d healed=%d admittedBefore=%d "
                             "admittedAfter=%d freeEdgesFinal=%d\n",
                             sliverCensusPartName((int)mv.nTri), closureHealIter,
                             closureHealTotal, closureHealAdmitted0,
                             countAdmittedRegions(rs, eprimeFill), countShellFreeEdges(sh));
            }
        }
        if (diagP2Enabled()) {
            dumpFillIdentity(sh, built, builtRid, eprimeFill, rs);
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
            if (gSeamFailFirst == 0 && freeE > 0 && (diagP2Enabled() || diagPlatesEnabled())) {
                gSeamFailFirst = 1;
                std::fprintf(stderr,
                             "DIAG_SEAMFAIL first=freeedge freeEdges=%d closed=%d "
                             "(see DIAG_FREEEDGE for rids/classes/endpoints)\n",
                             freeE, preClosed ? 1 : 0);
            }
            dumpDiagFreeEdges(sh, built, builtRid, eprimeFill, rs, mv, meshE, geom, collapsed,
                              sewTol);
        }
        if (diagWireEnabled() || diagP2Enabled()) {
            try {
                std::unordered_set<int> dumpedIdx;
                TopTools_IndexedDataMapOfShapeListOfShape ancW;
                TopExp::MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, ancW);
                TopTools_IndexedMapOfShape faceMapW;
                TopExp::MapShapes(sh, TopAbs_FACE, faceMapW);
                std::vector<int> faceRidW((size_t)faceMapW.Extent() + 1, -1);
                std::vector<int> faceToBuilt((size_t)faceMapW.Extent() + 1, -1);
                for (size_t i = 0; i < built.size(); i++) {
                    int idx = faceMapW.FindIndex(built[i]);
                    if (idx > 0) {
                        faceToBuilt[(size_t)idx] = (int)i;
                        if (i < builtRid.size()) faceRidW[(size_t)idx] = builtRid[i];
                    }
                }
                auto dumpIdx = [&](size_t i, const char* stage) {
                    if (i >= built.size() || !dumpedIdx.insert((int)i).second) return;
                    const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
                    dumpDiagWiresOfFace(rid, built[i], sewTol, stage, &ancW, &faceMapW, &faceRidW,
                                        &eprimeFill, &rs);
                };
                if (diagPlatesEnabled() || diagP2Enabled()) {
                    for (int hid = 0; hid <= 3; hid++) {
                        int nFaces = 0;
                        for (size_t i = 0; i < builtRid.size(); i++)
                            if (builtRid[i] == hid) nFaces++;
                        if (nFaces == 0 && !regionExploded(exploded, hid)) continue;
                        std::fprintf(stderr,
                                     "DIAG_HUBEXPLODE rid=%d exploded=%d facesAfter=%d\n", hid,
                                     regionExploded(exploded, hid) ? 1 : 0, nFaces);
                    }
                }
                for (size_t i = 0; i < built.size(); i++) {
                    const int rid = (i < builtRid.size()) ? builtRid[i] : -1;
                    if (rid >= 0 && rid <= 3) dumpIdx(i, "hub");
                    const Region* rr = regionById(rs, rid);
                    if (rr && rr->type == SurfType::Plane && (int)rr->tris.size() <= 7 &&
                        rr->builtAs == BuiltAs::Single)
                        dumpIdx(i, "sliver");
                }
                int nNbr = 0;
                for (int i = 1; i <= ancW.Extent() && nNbr < 80; i++) {
                    if (ancW(i).Extent() < 1) continue;
                    int bi[2] = {-1, -1};
                    int nA = 0;
                    bool anyFill = false;
                    for (TopTools_ListOfShape::Iterator it(ancW(i)); it.More() && nA < 2; it.Next()) {
                        const TopoDS_Face fa = TopoDS::Face(it.Value());
                        int ia = faceMapW.FindIndex(fa);
                        if (ia <= 0) continue;
                        int bidx = faceToBuilt[(size_t)ia];
                        if (bidx < 0) continue;
                        bi[nA++] = bidx;
                        const int rid = (size_t)bidx < builtRid.size() ? builtRid[(size_t)bidx] : -1;
                        if (rid >= 0 && (size_t)rid < eprimeFill.size() && eprimeFill[(size_t)rid])
                            anyFill = true;
                    }
                    if (!anyFill) continue;
                    for (int k = 0; k < nA; k++) {
                        dumpIdx((size_t)bi[k], "fill-nbr");
                        nNbr++;
                    }
                }
            } catch (const Standard_Failure&) {
            }
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
                    fireTolRewriteEdge(bb, e, fat, "sameParamRelax", spCap);
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
                    if (sliverCensusEnabled()) {
                        char pass[32];
                        std::snprintf(pass, sizeof(pass), "recover-%d", recoverPass);
                        for (size_t i = 0; i < built.size(); i++)
                            if (builtRid[i] == 100) diagRid100Status(pass, 100, built[i]);
                    }
                    goto try_rebuild;
                }
            } else if (recoverPass < 2 && explodeAll()) {
                recoverPass++;
                if (sliverCensusEnabled()) {
                    char pass[32];
                    std::snprintf(pass, sizeof(pass), "recover-%d", recoverPass);
                    for (size_t i = 0; i < built.size(); i++)
                        if (builtRid[i] == 100) diagRid100Status(pass, 100, built[i]);
                }
                goto try_rebuild;
            }
            restoreShared();
            dumpR2Probe(mv, rs, built, builtRid, eprimeFill, exploded, sewTol);
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
        // The valid first-pass shell has to be printable too, or "valid=1" is a
        // claim about the absence of a log line rather than a measurement.
        if (collapseDiagEnabled() && BRep_Tool::IsClosed(sh) && shValid)
            dumpShellCheck(sh, built, builtRid, rs, sewTol);
        // Site B is closed-but-invalid only. Open shells are site C (untouched).
        if (BRep_Tool::IsClosed(sh) && !shValid) {
            dumpShellCheck(sh, built, builtRid, rs, sewTol);
            std::vector<CascadeHit> culprits;
            collectFaceCulprits(built, builtRid, exploded, culprits);
            collectShellCulprits(sh, built, builtRid, rs, exploded, culprits);
            // D-130-9(b): a non-region face contributes its NEIGHBOURHOOD, not
            // the whole component. Without this the culprit set is empty and
            // RULE 1.4 escalates straight to the U2 blanket.
            collectNonRegionNeighbourhood(sh, built, builtRid, exploded, culprits);
            diagCascadeCulprits(culprits);
            // RULE 1.4: empty culprit set + invalid shell ⇒ escalate, never ship.
            CascadePlan plan = cascadeLadderPlan(cascadeSt, culprits, rs, exploded, false);
            if (applyCascadePlan(plan)) {
                recoverPass++;
                if (sliverCensusEnabled()) {
                    char pass[32];
                    std::snprintf(pass, sizeof(pass), "cascade-recover-%d", recoverPass);
                    for (size_t i = 0; i < built.size(); i++)
                        if (builtRid[i] == 100) diagRid100Status(pass, 100, built[i]);
                }
                goto try_rebuild;
            }
            if (plan.hostR2 || cascadeSt.u2Done || !shValid) {
                restoreShared();
                dumpR2Probe(mv, rs, built, builtRid, eprimeFill, exploded, sewTol);
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
                    dumpR2Probe(mv, rs, built, builtRid, eprimeFill, exploded, sewTol);
                    out.clear();
                    return false;
                }
                // After U2 the host D4.5 probe (RULE 2.0) is the volume authority.
                // Do not R2 a closed+valid U2 shell here — that is census 0.
            }
        }

        fitAnalyticTolerances(sh);

        dumpShippedInContext(built, builtRid, rs);
        dumpTolBindSummary();
        dumpR2Probe(mv, rs, built, builtRid, eprimeFill, exploded, sewTol);
        std::fprintf(stderr, "DIAG_RECONNECT_SUM fires=%d nTri=%d notDetermined=%d\n",
                     gReconnectFires, (int)mv.nTri, gOriNotDetermined);
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
