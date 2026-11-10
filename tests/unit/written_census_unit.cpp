// D-gid-7.6: the RESULT census of the written shape (src/written_census.hpp).
// Synthetic shells only: a box of 12 triangles whose provenance is assigned
// by face position, carried through the engine's own guarded unify, plus a
// hand-built history and a synthetic RegionSet for describeShippedFaces.
// SPDX-License-Identifier: MIT

#include "refit.hpp"
#include "written_census.hpp"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools_History.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Pnt.hxx>

#include <cmath>
#include <cstdio>
#include <vector>

namespace stl2step {
namespace detail {
TopoDS_Shape unifySameDomainGuarded(const TopoDS_Shape& shape, double angleDeg,
                                    const std::vector<TopoDS_Face>* keep,
                                    Handle(BRepTools_History)& histOut);
}
}  // namespace stl2step

using namespace stl2step::detail;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static TopoDS_Face triOut(const gp_Pnt& a, const gp_Pnt& b, const gp_Pnt& c,
                          const gp_XYZ& outward) {
    gp_XYZ n = (b.XYZ() - a.XYZ()).Crossed(c.XYZ() - a.XYZ());
    gp_Pnt b2 = b, c2 = c;
    if (n.Dot(outward) < 0) {
        b2 = c;
        c2 = b;
    }
    BRepBuilderAPI_MakePolygon poly(a, b2, c2, Standard_True);
    BRepBuilderAPI_MakeFace mf(poly.Wire(), Standard_True);
    return mf.Face();
}

static TopoDS_Shape sewSolid(const std::vector<TopoDS_Face>& faces) {
    BRepBuilderAPI_Sewing sew(1e-5);
    for (const auto& f : faces) sew.Add(f);
    sew.Perform();
    TopoDS_Shape sewed = sew.SewedShape();
    for (TopExp_Explorer ex(sewed, TopAbs_SHELL); ex.More(); ex.Next()) {
        BRep_Builder B;
        TopoDS_Solid so;
        B.MakeSolid(so);
        B.Add(so, TopoDS::Shell(ex.Current()));
        return so;
    }
    return sewed;
}

// 10 mm cube, two triangles per side.
static std::vector<TopoDS_Face> boxTriangles() {
    const gp_Pnt A(0, 0, 0), B(10, 0, 0), C(10, 10, 0), D(0, 10, 0);
    const gp_Pnt E(0, 0, 10), F(10, 0, 10), G(10, 10, 10), H(0, 10, 10);
    std::vector<TopoDS_Face> f;
    f.push_back(triOut(A, B, C, gp_XYZ(0, 0, -1)));
    f.push_back(triOut(A, C, D, gp_XYZ(0, 0, -1)));
    f.push_back(triOut(A, B, F, gp_XYZ(0, -1, 0)));
    f.push_back(triOut(A, F, E, gp_XYZ(0, -1, 0)));
    f.push_back(triOut(B, C, G, gp_XYZ(1, 0, 0)));
    f.push_back(triOut(B, G, F, gp_XYZ(1, 0, 0)));
    f.push_back(triOut(C, D, H, gp_XYZ(0, 1, 0)));
    f.push_back(triOut(C, H, G, gp_XYZ(0, 1, 0)));
    f.push_back(triOut(D, A, E, gp_XYZ(-1, 0, 0)));
    f.push_back(triOut(D, E, H, gp_XYZ(-1, 0, 0)));
    f.push_back(triOut(E, F, G, gp_XYZ(0, 0, 1)));
    f.push_back(triOut(E, G, H, gp_XYZ(0, 0, 1)));
    return f;
}

static gp_XYZ centroidOf(const TopoDS_Face& f) {
    GProp_GProps g;
    BRepGProp::SurfaceProperties(f, g);
    return g.CentreOfMass().XYZ();
}

int main() {
    // --- 1. identity map, every kind, units --------------------------------
    {
        TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Solid();
        std::vector<TopoDS_Face> faces;
        for (TopExp_Explorer ex(box, TopAbs_FACE); ex.More(); ex.Next())
            faces.push_back(TopoDS::Face(ex.Current()));
        expect(faces.size() == 6, "box-6-faces");
        // two faces of plane region 7, one of plane region 8, one facet,
        // one planar fillet-strip face, one "other" analytic face.
        const std::vector<uint8_t> kind = { kFaceSrcPlane, kFaceSrcPlane, kFaceSrcPlane,
                                            kFaceSrcFacet, kFaceSrcFillet, kFaceSrcOther };
        const std::vector<int> unit = { 7, 7, 8, -1, 9, 10 };
        FaceSources src;
        seedFaceSources(faces, src);
        remapFaceSources(Handle(BRepTools_History)(), src);  // null history: identity
        WrittenFaceCounts c = countWrittenFaces(box, src, kind, unit);
        expect(c.faces == 6, "id-faces-6");
        expect(c.facetFaces == 1, "id-facet-1");
        expect(c.analyticPlanes == 5, "id-analytic-planes-5 (planar fillet + other count)");
        expect(c.filletFaces == 1, "id-fillet-1");
        expect(c.planeUnits == 2, "id-plane-units-2 (region 7 twice, region 8)");
        expect(c.cylinderUnits == 0, "id-cylinder-units-0");
        expect(c.unknownFaces == 0 && c.mixedFaces == 0, "id-unknown-0-mixed-0");

        // An empty map: every written face is unknown, nothing else counts.
        FaceSources none;
        WrittenFaceCounts u = countWrittenFaces(box, none, kind, unit);
        expect(u.unknownFaces == 6 && u.facetFaces == 0 && u.analyticPlanes == 0,
               "empty-map-all-unknown");
    }

    // --- 2. hand-built history: merge, remove, unchanged ----------------------
    {
        TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Solid();
        std::vector<TopoDS_Face> f;
        for (TopExp_Explorer ex(box, TopAbs_FACE); ex.More(); ex.Next())
            f.push_back(TopoDS::Face(ex.Current()));
        FaceSources src;
        seedFaceSources(f, src);
        Handle(BRepTools_History) h = new BRepTools_History();
        h->AddModified(f[0], f[5]);  // f0 merged into f5
        h->Remove(f[1]);             // removed without image
        remapFaceSources(h, src);
        const std::vector<int>* s5 = src.Seek(f[5]);
        expect(s5 && s5->size() == 2, "hist-merge-image-has-2-sources");
        expect(src.Seek(f[0]) == nullptr, "hist-modified-face-gone");
        expect(src.Seek(f[1]) == nullptr, "hist-removed-face-dropped");
        const std::vector<int>* s2 = src.Seek(f[2]);
        expect(s2 && s2->size() == 1 && (*s2)[0] == 2, "hist-untouched-face-kept");
        expect(src.Size() == 4, "hist-map-size-4");
        // orientation is ignored (IsSame)
        expect(src.Seek(f[2].Reversed()) != nullptr, "hist-lookup-ignores-orientation");
    }

    // --- 3. the engine's guarded unify on a 12-triangle cube ----------------
    {
        TopoDS_Shape split = sewSolid(boxTriangles());
        std::vector<TopoDS_Face> built;
        for (TopExp_Explorer ex(split, TopAbs_FACE); ex.More(); ex.Next())
            built.push_back(TopoDS::Face(ex.Current()));
        expect(built.size() == 12, "cube-12-triangles");
        // top (z = 10): both facets -> one merged facet face.
        // bottom (z = 0): one facet + one plane-region triangle -> mixed plane.
        // each side: both triangles of one plane region -> one face, one unit.
        std::vector<uint8_t> kind(built.size(), 0);
        std::vector<int> unit(built.size(), -1);
        bool bottomFacetGiven = false;
        for (size_t i = 0; i < built.size(); ++i) {
            const gp_XYZ c = centroidOf(built[i]);
            if (c.Z() > 9.0) {
                kind[i] = kFaceSrcFacet;
            } else if (c.Z() < 1.0) {
                if (!bottomFacetGiven) {
                    kind[i] = kFaceSrcFacet;
                    bottomFacetGiven = true;
                } else {
                    kind[i] = kFaceSrcPlane;
                    unit[i] = 1;
                }
            } else {
                kind[i] = kFaceSrcPlane;
                unit[i] = c.X() < 1.0 ? 2 : c.X() > 9.0 ? 3 : c.Y() < 1.0 ? 4 : 5;
            }
        }
        FaceSources src;
        seedFaceSources(built, src);
        WrittenFaceCounts before = countWrittenFaces(split, src, kind, unit);
        expect(before.facetFaces == 3 && before.analyticPlanes == 9,
               "cube-before-unify-3-facets-9-planes");
        Handle(BRepTools_History) hist;
        TopoDS_Shape merged = unifySameDomainGuarded(split, 0.001, nullptr, hist);
        expect(!hist.IsNull(), "cube-unify-history-returned");
        remapFaceSources(hist, src);
        WrittenFaceCounts c = countWrittenFaces(merged, src, kind, unit);
        std::fprintf(stderr,
                     "INFO cube after unify faces=%d facet=%d planes=%d units=%d mixed=%d "
                     "unknown=%d\n",
                     c.faces, c.facetFaces, c.analyticPlanes, c.planeUnits, c.mixedFaces,
                     c.unknownFaces);
        expect(c.faces == 6, "cube-unified-6-faces");
        expect(c.unknownFaces == 0, "cube-unified-no-unknown-face");
        expect(c.facetFaces == 1, "cube-unified-1-facet-face (two facets merged)");
        expect(c.analyticPlanes == 5, "cube-unified-5-analytic-planes");
        expect(c.mixedFaces == 1, "cube-unified-1-mixed (facet merged into a plane)");
        expect(c.planeUnits == 5, "cube-unified-5-plane-units");
    }

    // --- 4. describeShippedFaces on a synthetic RegionSet --------------------
    {
        TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Solid();
        TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(2.0, 5.0).Solid();
        std::vector<TopoDS_Face> built;
        for (TopExp_Explorer ex(box, TopAbs_FACE); ex.More(); ex.Next())
            built.push_back(TopoDS::Face(ex.Current()));
        TopoDS_Face lateral;
        for (TopExp_Explorer ex(cyl, TopAbs_FACE); ex.More(); ex.Next()) {
            const TopoDS_Face f = TopoDS::Face(ex.Current());
            if (BRepAdaptor_Surface(f, Standard_False).GetType() == GeomAbs_Cylinder) lateral = f;
        }
        built.push_back(lateral);  // index 6
        built.push_back(lateral);  // index 7
        stl2step::refit::RegionSet rs;
        auto region = [](int id, stl2step::refit::SurfType t, stl2step::refit::Origin o) {
            stl2step::refit::Region r{};
            r.id = id;
            r.type = t;
            r.origin = o;
            return r;
        };
        using stl2step::refit::Origin;
        using stl2step::refit::SurfType;
        rs.regions.push_back(region(0, SurfType::Plane, Origin::PlaneGrow));
        rs.regions.push_back(region(1, SurfType::Cylinder, Origin::CylGrow));
        rs.regions.push_back(region(2, SurfType::Cylinder, Origin::FilletStrip));
        rs.regions.push_back(region(3, SurfType::Torus, Origin::TorusBlend));
        rs.regions.push_back(region(4, SurfType::Plane, Origin::NgonWall));
        // faces: 0 plane r0, 1 facet of exploded r0, 2 island facet, 3 cylinder r1,
        // 4 fillet r2, 5 torus r3, 6 stage-P cylinder, 7 plane r4 (kind by region)
        rs.shipFaceRid = { 0, 0, -1, 1, 2, 3, -1, 4 };
        rs.shipFaceFacet = { 0, 1, 1, 0, 0, 0, 0, 0 };
        std::vector<uint8_t> kind;
        std::vector<int> unit;
        describeShippedFaces(rs, built, kind, unit);
        const std::vector<uint8_t> wantK = { kFaceSrcPlane,  kFaceSrcFacet,    kFaceSrcFacet,
                                             kFaceSrcCylinder, kFaceSrcFillet, kFaceSrcOther,
                                             kFaceSrcCylinder, kFaceSrcPlane };
        const std::vector<int> wantU = { 0, -1, -1, 1, 2, 3, -2 - 6, 4 };
        expect(kind == wantK, "describe-kinds");
        expect(unit == wantU, "describe-units (stage-P face is its own unit)");

        rs.shipFaceFacet.pop_back();  // provenance not parallel to the faces
        describeShippedFaces(rs, built, kind, unit);
        bool allZero = true;
        for (uint8_t k : kind) allZero = allZero && k == 0;
        expect(allZero && kind.size() == built.size(), "describe-mismatch-gives-unknown");
    }

    if (gFail) std::fprintf(stderr, "written_census_unit: %d FAILED\n", gFail);
    else std::fprintf(stderr, "written_census_unit: all passed\n");
    return gFail ? 1 : 0;
}
