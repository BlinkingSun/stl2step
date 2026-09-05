#include "grade_step.hpp"

#include <BRep_Tool.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

namespace grade {
namespace {

void silenceOcct() {
    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);
}

SurfParams fromAdaptor(const BRepAdaptor_Surface& ads, SurfClass& cls) {
    SurfParams S;
    switch (ads.GetType()) {
        case GeomAbs_Plane: {
            cls = SurfClass::Plane;
            const gp_Pln pl = ads.Plane();
            const gp_Dir d = pl.Axis().Direction();
            S.cls = cls;
            S.n = normalized(Vec3{d.X(), d.Y(), d.Z()});
            const gp_Pnt loc = pl.Location();
            const Vec3 lp{loc.X(), loc.Y(), loc.Z()};
            S.p0 = S.n * dot(lp, S.n);
            break;
        }
        case GeomAbs_Cylinder: {
            cls = SurfClass::Cylinder;
            const gp_Cylinder c = ads.Cylinder();
            const gp_Ax1 ax = c.Axis();
            S.cls = cls;
            S.n = normalized(Vec3{ax.Direction().X(), ax.Direction().Y(), ax.Direction().Z()});
            S.p0 = Vec3{ax.Location().X(), ax.Location().Y(), ax.Location().Z()};
            S.p0 = footFromOrigin(S.p0, S.n);
            S.R = c.Radius();
            break;
        }
        case GeomAbs_Cone: {
            cls = SurfClass::Cone;
            const gp_Cone c = ads.Cone();
            S.cls = cls;
            S.n = normalized(Vec3{c.Axis().Direction().X(), c.Axis().Direction().Y(),
                                  c.Axis().Direction().Z()});
            const gp_Pnt ap = c.Apex();
            S.apex = Vec3{ap.X(), ap.Y(), ap.Z()};
            S.alpha = std::fabs(c.SemiAngle());
            S.p0 = Vec3{c.Location().X(), c.Location().Y(), c.Location().Z()};
            break;
        }
        case GeomAbs_Sphere: {
            cls = SurfClass::Sphere;
            const gp_Sphere s = ads.Sphere();
            S.cls = cls;
            S.p0 = Vec3{s.Location().X(), s.Location().Y(), s.Location().Z()};
            S.R = s.Radius();
            break;
        }
        case GeomAbs_Torus: {
            cls = SurfClass::Torus;
            const gp_Torus t = ads.Torus();
            S.cls = cls;
            S.n = normalized(Vec3{t.Axis().Direction().X(), t.Axis().Direction().Y(),
                                  t.Axis().Direction().Z()});
            S.p0 = Vec3{t.Location().X(), t.Location().Y(), t.Location().Z()};
            S.p0 = footFromOrigin(S.p0, S.n);
            S.R = t.MajorRadius();
            S.r = t.MinorRadius();
            break;
        }
        case GeomAbs_BSplineSurface:
            cls = SurfClass::Other;
            S.cls = cls;
            break;
        default:
            cls = SurfClass::Other;
            S.cls = cls;
            break;
    }
    return S;
}

}  // namespace

bool loadStep(const std::string& path, const Mesh& mesh, StepModel& out, std::string& err,
              bool computeVolume) {
    silenceOcct();
    out = StepModel{};
    out.cover.assign(mesh.tris.size(), -1);
    STEPControl_Reader r;
    if (r.ReadFile(path.c_str()) != IFSelect_RetDone) {
        err = "STEP ReadFile failed";
        return false;
    }
    if (r.TransferRoots() < 1) {
        err = "STEP TransferRoots produced no shape";
        return false;
    }
    out.shape = r.OneShape();

    BRepCheck_Analyzer ana(out.shape, Standard_True, Standard_False);
    out.valid = ana.IsValid() == Standard_True;

    int nClosed = 0, nShells = 0;
    for (TopExp_Explorer sh(out.shape, TopAbs_SHELL); sh.More(); sh.Next()) {
        ++nShells;
        if (BRep_Tool::IsClosed(sh.Current())) ++nClosed;
    }
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaces;
    TopExp::MapShapesAndAncestors(out.shape, TopAbs_EDGE, TopAbs_FACE, edgeFaces);
    int freeEdges = 0;
    for (int i = 1; i <= edgeFaces.Extent(); ++i) {
        if (edgeFaces(i).Extent() < 2) ++freeEdges;
    }
    // Seamed360 cylinders/spheres carry a seam EDGE with one ancestor
    // FACE; that is not a free edge. D-140-1 watertight is "no closed
    // shell" — use BRep_Tool::IsClosed per shell (census convention).
    (void)freeEdges;
    out.watertight = (nShells > 0 && nClosed == nShells);

    if (computeVolume) {
        GProp_GProps vp;
        BRepGProp::VolumeProperties(out.shape, vp, Precision::Confusion());
        out.volume = vp.Mass();
    }

    TopTools_IndexedMapOfShape fmap;
    TopExp::MapShapes(out.shape, TopAbs_FACE, fmap);
    out.faces.reserve(static_cast<size_t>(fmap.Extent()));
    for (int i = 1; i <= fmap.Extent(); ++i) {
        StepFace sf;
        sf.entity = i;
        sf.face = TopoDS::Face(fmap(i));
        BRepAdaptor_Surface ads(sf.face, Standard_True);
        sf.S = fromAdaptor(ads, sf.cls);
        GProp_GProps sp;
        BRepGProp::SurfaceProperties(sf.face, sp);
        sf.area = sp.Mass();
        if (ads.GetType() == GeomAbs_BSplineSurface) ++out.nBSpline;
        switch (sf.cls) {
            case SurfClass::Plane: ++out.nPlane; break;
            case SurfClass::Cylinder: ++out.nCyl; break;
            case SurfClass::Cone: ++out.nCone; break;
            case SurfClass::Sphere: ++out.nSphere; break;
            case SurfClass::Torus: ++out.nTorus; break;
            default: ++out.nOther; break;
        }
        TopTools_IndexedMapOfShape vmap;
        TopExp::MapShapes(sf.face, TopAbs_VERTEX, vmap);
        bool allMesh = (sf.cls == SurfClass::Plane);
        for (int vi = 1; vi <= vmap.Extent(); ++vi) {
            const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(vmap(vi)));
            const Vec3 pv{p.X(), p.Y(), p.Z()};
            sf.wireVerts.push_back(pv);
            const int mid = mesh.nearestVert(pv, mesh.q);
            sf.meshVerts.push_back(mid);
            if (mid < 0) allMesh = false;
        }
        if (sf.cls == SurfClass::Plane && allMesh && !sf.meshVerts.empty()) {
            std::vector<char> onWire(mesh.verts.size(), 0);
            for (int mid : sf.meshVerts)
                if (mid >= 0) onWire[static_cast<size_t>(mid)] = 1;
            std::vector<int> S;
            double sumA = 0;
            for (int t = 0; t < static_cast<int>(mesh.tris.size()); ++t) {
                const Tri& tr = mesh.tris[static_cast<size_t>(t)];
                if (!onWire[static_cast<size_t>(tr.v[0])] || !onWire[static_cast<size_t>(tr.v[1])] ||
                    !onWire[static_cast<size_t>(tr.v[2])])
                    continue;
                bool onPl = true;
                for (int k = 0; k < 3; ++k) {
                    if (distToSurf(mesh.verts[static_cast<size_t>(tr.v[k])], sf.S) > mesh.tau) {
                        onPl = false;
                        break;
                    }
                }
                if (!onPl) continue;
                S.push_back(t);
                sumA += tr.area;
            }
            // A tessellation facet is one mesh triangle. A design plane covering
            // two or more triangles (S01 cube face, unified verbatim planes) is
            // not a facet face — Step C assigns it as recovered.
            if (S.size() == 1 && std::fabs(sf.area - sumA) <= areaQ(mesh, S)) {
                sf.facet = true;
                ++out.nFacet;
                for (int t : S) out.cover[static_cast<size_t>(t)] = sf.entity;
            }
        }
        out.faces.push_back(std::move(sf));
    }
    return true;
}

}  // namespace grade
