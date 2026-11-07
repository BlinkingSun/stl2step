// D-gid-2.1: free-edge guard, dropped-face history predicate, written-body
// retally, and unifySameDomainGuarded on a synthetic shell. A sliver that
// OCCT does not drop is reported; the helpers are still required to pass.
// SPDX-License-Identifier: MIT

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools_History.hxx>
#include <BRep_Builder.hxx>
#include <Precision.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <cstdio>
#include <vector>

namespace stl2step {
namespace detail {
int countFreeEdges(const TopoDS_Shape& shape);
std::vector<TopoDS_Face> facesRemovedWithoutImage(
    const Handle(BRepTools_History)& hist,
    const std::vector<TopoDS_Face>& faces);
TopoDS_Shape keepIfFreeEdgesNotIncreased(const TopoDS_Shape& input,
                                        const TopoDS_Shape& unified);
void tallyWrittenBodies(const TopoDS_Shape& written, int& solids, int& openShells);
TopoDS_Shape unifySameDomainGuarded(const TopoDS_Shape& shape, double angleDeg,
                                    const std::vector<TopoDS_Face>* keep);
}
}

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static int countFaces(const TopoDS_Shape& s) {
    int n = 0;
    for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next()) ++n;
    return n;
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
    TopoDS_Shell sh;
    bool got = false;
    for (TopExp_Explorer ex(sewed, TopAbs_SHELL); ex.More(); ex.Next()) {
        sh = TopoDS::Shell(ex.Current());
        got = true;
        break;
    }
    if (!got) return sewed;
    BRep_Builder B;
    TopoDS_Solid so;
    B.MakeSolid(so);
    B.Add(so, sh);
    return so;
}

static std::vector<TopoDS_Face> boxTriangles(bool sliver) {
    const gp_Pnt A(0, 0, 0), B(10, 0, 0), C(10, 10, 0), D(0, 10, 0);
    const gp_Pnt E(0, 0, 10), F(10, 0, 10), G(10, 10, 10), H(0, 10, 10);
    std::vector<TopoDS_Face> faces;
    auto add = [&](const gp_Pnt& a, const gp_Pnt& b, const gp_Pnt& c, const gp_XYZ& o) {
        faces.push_back(triOut(a, b, c, o));
    };
    add(A, B, C, gp_XYZ(0, 0, -1));
    add(A, C, D, gp_XYZ(0, 0, -1));
    add(A, B, F, gp_XYZ(0, -1, 0));
    add(A, F, E, gp_XYZ(0, -1, 0));
    add(B, C, G, gp_XYZ(1, 0, 0));
    add(B, G, F, gp_XYZ(1, 0, 0));
    add(C, D, H, gp_XYZ(0, 1, 0));
    add(C, H, G, gp_XYZ(0, 1, 0));
    add(D, A, E, gp_XYZ(-1, 0, 0));
    add(D, E, H, gp_XYZ(-1, 0, 0));
    if (!sliver) {
        add(E, F, G, gp_XYZ(0, 0, 1));
        add(E, G, H, gp_XYZ(0, 0, 1));
        return faces;
    }
    // Long thin top triangle, normal ~0.0007 deg off the rest of the lid.
    // 0.000698 deg * pi/180 / 20, from tan(theta) ~= 20*|dz| on a 10 x 0.05 lid.
    const double dz = 6.091e-7;
    const gp_Pnt S(5, 0.05, 10.0 + dz);
    add(E, F, S, gp_XYZ(0, 0, 1));
    add(E, S, H, gp_XYZ(0, 0, 1));
    add(F, G, S, gp_XYZ(0, 0, 1));
    add(S, G, H, gp_XYZ(0, 0, 1));
    return faces;
}

int main() {
    using stl2step::detail::countFreeEdges;
    using stl2step::detail::facesRemovedWithoutImage;
    using stl2step::detail::keepIfFreeEdgesNotIncreased;
    using stl2step::detail::tallyWrittenBodies;
    using stl2step::detail::unifySameDomainGuarded;

    TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Solid();
    expect(countFreeEdges(box) == 0, "box-free-edges-0");
    int solids = -1, openShells = -1;
    tallyWrittenBodies(box, solids, openShells);
    expect(solids == 1 && openShells == 0, "box-tally-1-0");

    std::vector<TopoDS_Face> five;
    int skipped = 0;
    for (TopExp_Explorer ex(box, TopAbs_FACE); ex.More(); ex.Next()) {
        if (skipped++ == 0) continue;
        five.push_back(TopoDS::Face(ex.Current()));
    }
    BRepBuilderAPI_Sewing sew(1e-5);
    for (const auto& f : five) sew.Add(f);
    sew.Perform();
    TopoDS_Shape openSh = sew.SewedShape();
    expect(countFreeEdges(openSh) > 0, "open-shell-has-free-edges");
    tallyWrittenBodies(openSh, solids, openShells);
    expect(solids == 0 && openShells == 1, "open-shell-tally-0-1");

    TopoDS_Shape kept = keepIfFreeEdgesNotIncreased(box, openSh);
    expect(kept.IsSame(box), "guard-keeps-input-when-free-edges-rise");
    TopoDS_Shape same = keepIfFreeEdgesNotIncreased(box, box);
    expect(same.IsSame(box), "guard-keeps-unified-when-free-edges-do-not-rise");

    BRep_Builder B;
    TopoDS_Compound comp;
    B.MakeCompound(comp);
    B.Add(comp, box);
    B.Add(comp, openSh);
    tallyWrittenBodies(comp, solids, openShells);
    expect(solids == 1 && openShells == 1, "compound-tally-1-1");

    Handle(BRepTools_History) hist = new BRepTools_History();
    std::vector<TopoDS_Face> histFaces = five;
    hist->Remove(histFaces[0]);
    hist->AddModified(histFaces[1], histFaces[2]);
    std::vector<TopoDS_Face> dropped = facesRemovedWithoutImage(hist, histFaces);
    expect(dropped.size() == 1 && dropped[0].IsSame(histFaces[0]),
           "history-dropped-is-removed-without-image");
    bool modifiedListed = false;
    for (const auto& f : dropped)
        if (f.IsSame(histFaces[1])) modifiedListed = true;
    expect(!modifiedListed, "history-modified-face-is-not-dropped");

    TopoDS_Shape split = sewSolid(boxTriangles(false));
    expect(countFreeEdges(split) == 0, "split-box-closed");
    expect(countFaces(split) == 12, "split-box-12-faces");
    TopoDS_Shape merged = unifySameDomainGuarded(split, 0.001, nullptr);
    expect(countFreeEdges(merged) == 0, "merged-box-free-edges-0");
    expect(countFaces(merged) == 6, "merged-box-6-faces");
    tallyWrittenBodies(merged, solids, openShells);
    expect(solids == 1 && openShells == 0, "merged-box-tally-1-0");
    BRepCheck_Analyzer mergedOk(merged, Standard_True, Standard_True);
    expect(mergedOk.IsValid() == Standard_True, "merged-box-valid");

    TopoDS_Shape sliverIn = sewSolid(boxTriangles(true));
    const int sliverFreeIn = countFreeEdges(sliverIn);
    const int sliverFacesIn = countFaces(sliverIn);
    expect(sliverFreeIn == 0, "sliver-box-closed");
    std::fprintf(stderr, "INFO sliver input faces=%d free=%d\n", sliverFacesIn, sliverFreeIn);

    ShapeUpgrade_UnifySameDomain usd(sliverIn, Standard_True, Standard_True, Standard_False);
    usd.SetAngularTolerance(0.001 * 3.14159265358979323846 / 180.0);
    usd.SetLinearTolerance(Precision::Confusion());
    std::vector<TopoDS_Face> sliverFaces;
    for (TopExp_Explorer ex(sliverIn, TopAbs_FACE); ex.More(); ex.Next())
        sliverFaces.push_back(TopoDS::Face(ex.Current()));
    usd.Build();
    std::vector<TopoDS_Face> usdDropped =
        facesRemovedWithoutImage(usd.History(), sliverFaces);
    const int usdFree = countFreeEdges(usd.Shape());
    std::fprintf(stderr, "INFO unguarded sliver dropped=%d free=%d faces=%d\n",
                 (int)usdDropped.size(), usdFree, countFaces(usd.Shape()));

    TopoDS_Shape guarded = unifySameDomainGuarded(sliverIn, 0.001, nullptr);
    const int gFree = countFreeEdges(guarded);
    expect(gFree <= sliverFreeIn, "sliver-guard-does-not-add-free-edges");
    tallyWrittenBodies(guarded, solids, openShells);
    expect(solids == 1 && openShells == 0, "sliver-tally-1-0");
    if (!usdDropped.empty()) {
        expect(gFree == 0, "sliver-rebuild-closes");
        expect(countFaces(guarded) > countFaces(usd.Shape()), "sliver-rebuild-keeps-face");
        BRepCheck_Analyzer gOk(guarded, Standard_True, Standard_True);
        expect(gOk.IsValid() == Standard_True, "sliver-rebuild-valid");
        std::fprintf(stderr, "INFO synthetic USD drop reproduced\n");
    } else {
        std::fprintf(stderr,
                     "INFO synthetic USD drop not reproduced; guard + retally helpers covered\n");
    }

    if (gFail) {
        std::fprintf(stderr, "FAILED %d\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "OK\n");
    return 0;
}
