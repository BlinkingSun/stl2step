// Synthetic STL corpus generator for stl2step smooth-feature acceptance gates.
// SPDX-License-Identifier: MIT

#include "corpus_util.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <stdexcept>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Standard_Failure.hxx>

namespace fs = std::filesystem;
using namespace corpus;

namespace {

struct FixtureResult {
    Sidecar sidecar;
    MeshData mesh;
};

Recoverable cylRec(double R, const Vec3& loc, const Vec3& dir, int count, int nSides,
                   bool closed360) {
    Recoverable r;
    r.type = "cylinder";
    r.radius = R;
    r.axis = {loc, normalize(dir)};
    r.count = count;
    r.nSides = nSides;
    r.closed360 = closed360;
    return r;
}

Recoverable planeRec(const Vec3& normal, int count = 1) {
    Recoverable r;
    r.type = "plane";
    r.normal = normalize(normal);
    r.count = count;
    return r;
}

Recoverable coneRec(double R0, double R1, double halfAngleDeg, const Vec3& loc, const Vec3& dir,
                    int count, int nSides) {
    Recoverable r;
    r.type = "cone";
    r.radius = R0;
    r.radius2 = R1;
    r.halfAngleDeg = halfAngleDeg;
    r.axis = {loc, normalize(dir)};
    r.count = count;
    r.nSides = nSides;
    r.closed360 = true;
    return r;
}

// Battery 130 synthetic fixtures: ≥48-side regular polylines (SPEC-130-fixtures).
constexpr int kTessN = 48;

TopoDS_Wire nGonWire(double R, int N, const gp_Ax2& ax, double zAlong = 0.0) {
    BRepBuilderAPI_MakePolygon poly;
    const gp_Pnt loc = ax.Location().Translated(gp_Vec(ax.Direction()) * zAlong);
    for (int i = 0; i < N; ++i) {
        const double ang = 2.0 * M_PI * i / N;
        poly.Add(loc.Translated(gp_Vec(ax.XDirection()) * (R * std::cos(ang)) +
                                gp_Vec(ax.YDirection()) * (R * std::sin(ang))));
    }
    poly.Close();
    return poly.Wire();
}

TopoDS_Shape makeNGonPrism(double R, double H, int N, const gp_Ax2& ax) {
    TopoDS_Face face = BRepBuilderAPI_MakeFace(nGonWire(R, N, ax, 0.0));
    return BRepPrimAPI_MakePrism(face, gp_Vec(ax.Direction()) * H);
}

// Ruled solid frustum (planar N-gon caps + planar trapezoid sides).
TopoDS_Shape makeNGonFrustum(double R0, double R1, double H, int N, const gp_Ax2& ax) {
    BRepOffsetAPI_ThruSections loft(Standard_True, Standard_True);
    loft.AddWire(nGonWire(R0, N, ax, 0.0));
    loft.AddWire(nGonWire(R1, N, ax, H));
    loft.Build();
    return loft.Shape();
}

TopoDS_Shape cutThroughNGon(const TopoDS_Shape& stock, const gp_Pnt& faceCenter, const gp_Dir& dir,
                            double R, double span, int N) {
    gp_Ax2 ax(faceCenter.Translated(gp_Vec(dir) * (-1.0)), dir);
    return BRepAlgoAPI_Cut(stock, makeNGonPrism(R, span + 2.0, N, ax));
}

std::vector<Recoverable> boxSixPlanes() {
    return {planeRec({0, 0, 1}, 1),  planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
            planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1),  planeRec({0, -1, 0}, 1)};
}

TopoDS_Shape makePlateWithHole(double W, double D, double T, double holeR, int holeN,
                               const gp_Pnt& holeCenter) {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(W, D, T);
    gp_Ax2 ax(holeCenter, gp_Dir(0, 0, 1));
    TopoDS_Shape hole = makeNGonPrism(holeR, T + 2.0, holeN, ax);
    gp_Trsf tr;
    tr.SetTranslation(gp_Vec(0, 0, -1.0));
    hole = BRepBuilderAPI_Transform(hole, tr).Shape();
    return BRepAlgoAPI_Cut(plate, hole);
}

// Primitive-count planes: one recoverable entry per connected coplanar patch
// (G2: count is the number of Geom_Plane faces, never tessellation triangles).
void appendMeshPlanes(Sidecar& sc, const MeshData& mesh, double clusterDeg = 8.0) {
    const int n = static_cast<int>(mesh.tris.size());
    if (n == 0) return;
    std::vector<Vec3> nrm(n);
    std::vector<double> area(n);
    for (int i = 0; i < n; ++i) {
        const auto& t = mesh.tris[i];
        nrm[i] = triNormal(mesh.verts[t[0]], mesh.verts[t[1]], mesh.verts[t[2]]);
        area[i] = triArea(mesh.verts[t[0]], mesh.verts[t[1]], mesh.verts[t[2]]);
    }
    auto ek = [](int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
    std::map<std::pair<int, int>, std::vector<int>> e2t;
    for (int i = 0; i < n; ++i) {
        const auto& t = mesh.tris[i];
        e2t[ek(t[0], t[1])].push_back(i);
        e2t[ek(t[1], t[2])].push_back(i);
        e2t[ek(t[2], t[0])].push_back(i);
    }
    const double cth = std::cos(clusterDeg * M_PI / 180.0);
    std::vector<int> seen(n, 0);
    std::vector<Recoverable> planes;
    for (int i = 0; i < n; ++i) {
        if (seen[i]) continue;
        std::vector<int> st = {i};
        seen[i] = 1;
        Vec3 acc = scale(nrm[i], area[i]);
        double aacc = area[i];
        while (!st.empty()) {
            const int u = st.back();
            st.pop_back();
            const auto& t = mesh.tris[u];
            const int vs[3] = {t[0], t[1], t[2]};
            for (int k = 0; k < 3; ++k) {
                const auto it = e2t.find(ek(vs[k], vs[(k + 1) % 3]));
                if (it == e2t.end()) continue;
                for (int v : it->second) {
                    if (seen[v]) continue;
                    if (dot(nrm[u], nrm[v]) < cth) continue;
                    seen[v] = 1;
                    acc = add(acc, scale(nrm[v], area[v]));
                    aacc += area[v];
                    st.push_back(v);
                }
            }
        }
        if (aacc < 1e-12) continue;
        planes.push_back(planeRec(normalize(acc), 1));
    }
    std::sort(planes.begin(), planes.end(), [](const Recoverable& a, const Recoverable& b) {
        if (a.normal.z != b.normal.z) return a.normal.z < b.normal.z;
        if (a.normal.y != b.normal.y) return a.normal.y < b.normal.y;
        return a.normal.x < b.normal.x;
    });
    sc.recoverable.insert(sc.recoverable.end(), planes.begin(), planes.end());
}

// 12 edge-fillet (or 45° chamfer) cylinders of a cube [0,L]^3, radius/setback R.
// vMin/vMax clip out the spherical vertex blends so nSides is the strip ring count.
std::vector<Recoverable> cubeEdgeCylinders(double L, double R) {
    const double I = L - R;
    const double lo = R + 0.45 * R;
    const double hi = L - R - 0.45 * R;
    auto one = [&](Vec3 loc, Vec3 dir) {
        Recoverable r = cylRec(R, loc, dir, 1, 0, false);
        r.vMin = lo;
        r.vMax = hi;
        return r;
    };
    return {
        one({0, R, R}, {1, 0, 0}), one({0, R, I}, {1, 0, 0}),
        one({0, I, R}, {1, 0, 0}), one({0, I, I}, {1, 0, 0}),
        one({R, 0, R}, {0, 1, 0}), one({R, 0, I}, {0, 1, 0}),
        one({I, 0, R}, {0, 1, 0}), one({I, 0, I}, {0, 1, 0}),
        one({R, R, 0}, {0, 0, 1}), one({R, I, 0}, {0, 0, 1}),
        one({I, R, 0}, {0, 0, 1}), one({I, I, 0}, {0, 0, 1}),
    };
}

// Closed box with an exact nBands-row quarter-round at the origin X-edge.
// Built as planar N-gon faces so tessellation cannot invent extra bands.
TopoDS_Shape makeBandedQuarterRoundBox(double L, double R, int nBands) {
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(0, 0, L));
    for (int i = 0; i <= nBands; ++i) {
        const double th = (M_PI / 2.0) * static_cast<double>(i) / static_cast<double>(nBands);
        const double y = R * (1.0 - std::cos(th));
        const double z = R * (1.0 - std::sin(th));
        poly.Add(gp_Pnt(0, y, z));
    }
    poly.Add(gp_Pnt(0, L, 0));
    poly.Add(gp_Pnt(0, L, L));
    poly.Close();
    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire());
    return BRepPrimAPI_MakePrism(face, gp_Vec(L, 0, 0));
}

TopoDS_Shape makeNonConcyclicHexHole(double W, double D, double T, double R,
                                     const gp_Pnt& center) {
    BRepBuilderAPI_MakePolygon poly;
    for (int i = 0; i < 6; ++i) {
        const double ang = 2.0 * M_PI * i / 6.0;
        const double r = (i % 2 == 0) ? R : (1.15 * R);
        poly.Add(center.Translated(gp_Vec(std::cos(ang) * r, std::sin(ang) * r, 0)));
    }
    poly.Close();
    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire());
    TopoDS_Shape hole = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, T + 2.0));
    gp_Trsf tr;
    tr.SetTranslation(gp_Vec(0, 0, -1.0));
    hole = BRepBuilderAPI_Transform(hole, tr).Shape();
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(W, D, T);
    return BRepAlgoAPI_Cut(plate, hole);
}

void fillLiveExpectations(FixtureResult& out);

void fillMeshSidecar(FixtureResult& out) {
    out.sidecar.triangleCount = static_cast<int>(out.mesh.tris.size());
    out.sidecar.meshVolume = meshVolume(out.mesh);
    const EdgeStats es = meshEdgeStats(out.mesh);
    out.sidecar.openEdges = es.open;
    out.sidecar.nonManifoldEdges = es.nonManifold;
    out.sidecar.watertight = (es.open == 0 && es.nonManifold == 0);
    out.sidecar.components = splitMeshComponents(out.mesh);
    for (auto& r : out.sidecar.recoverable) {
        if (r.type == "cylinder")
            r.nSides = measureRecoverableNSides(out.mesh, r);
    }
    fillLiveExpectations(out);
}

void sumRecoverableCensus(const Sidecar& sc, int& planes, int& cylinders, int& cones) {
    planes = cylinders = cones = 0;
    for (const auto& r : sc.recoverable) {
        if (r.type == "plane") planes += r.count;
        else if (r.type == "cylinder") cylinders += r.count;
        else if (r.type == "cone") cones += r.count;
    }
}

void fillLiveExpectations(FixtureResult& out) {
    Sidecar& sc = out.sidecar;
    sc.live.clear();
    if (sc.components.empty()) return;
    int recPlanes = 0, recCyl = 0, recCone = 0;
    sumRecoverableCensus(sc, recPlanes, recCyl, recCone);
    const bool hasFreeform =
        std::any_of(sc.mustRemainFaceted.begin(), sc.mustRemainFaceted.end(),
                    [](const FacetedRegion& f) { return f.type == "freeform"; });
    const bool hasTorus =
        std::any_of(sc.mustRemainFaceted.begin(), sc.mustRemainFaceted.end(),
                    [](const FacetedRegion& f) { return f.type == "torus"; });
    for (const MeshComponentInfo& comp : sc.components) {
        LiveExpectation live;
        live.component = comp.index;
        live.volumeBudgetMM3 = volumeBudgetMM3(std::fabs(comp.meshVolume));
        live.disposition = "PASS";
        if (sc.id == "S02") {
            live.disposition = "ESCALATE";
            live.escalateReason =
                "P1 micro-plane regions on spherical vertex blends (nIslands=0); "
                "target 6 planes + 12 R=2 cylinders after island fix.";
            live.surfaceCensus = {6, 12, 0};
            live.faceCount = 316;
        } else if (sc.id == "S03") {
            live.disposition = "ESCALATE";
            live.escalateReason =
                "P1 shatters drafted cone into 1-tri planes; target seamed360=4, "
                "cylinders=4, faceCount~106 after island fix.";
            live.surfaceCensus = {102, 4, 0};
            live.faceCount = 106;
        } else if (sc.id == "S04") {
            live.disposition = "ESCALATE";
            live.escalateReason =
                "Boss-top torus blend (TorusNYI); P1 must mark blend annulus islands, "
                "not cylinders.";
            live.surfaceCensus = {recPlanes, recCyl, 0};
            live.faceCount = recPlanes + recCyl;
        } else if (sc.id == "S05") {
            live.surfaceCensus = {10, 2, 0};
            live.faceCount = 12;
            live.volumeBudgetMM3 = 21.3;
        } else if (sc.id == "S09" && comp.index == 0) {
            live.surfaceCensus = {6, 0, 0};
            live.faceCount = 6;
        } else if (sc.id == "S09" && comp.index == 1) {
            live.surfaceCensus = {0, 0, comp.triangleCount};
            live.faceCount = comp.triangleCount;
        } else if (sc.id == "S12") {
            live.surfaceCensus = {6, 0, 0};
            live.faceCount = 6;
        } else if (sc.id == "S16-R2-ChainUnstable") {
            live.buildFaces = false;
            live.surfaceCensus = {comp.triangleCount, 0, 0};
            live.faceCount = comp.triangleCount;
        } else if (hasFreeform && comp.index == static_cast<int>(sc.components.size()) - 1) {
            live.surfaceCensus = {0, 0, comp.triangleCount};
            live.faceCount = comp.triangleCount;
        } else {
            live.surfaceCensus = {recPlanes, recCyl, 0, recCone};
            live.faceCount = recPlanes + recCyl + recCone;
            if (!sc.battery.empty()) {
                live.builtCylindersFloor = recCyl;
                live.builtPlanesFloor = recPlanes;
                live.builtConesFloor = recCone;
            }
            if (sc.components.size() == 1) {
                // Single-body fixtures: planes + cylinders (+ cones) is the recoverable census.
            } else if (sc.id == "S12") {
                live.surfaceCensus = {6, 0, 0};
                live.faceCount = 6;
            }
        }
        if (hasTorus && sc.id == "S04") {
            (void)hasTorus;
        }
        sc.live.push_back(live);
    }
}

FixtureResult emitShape(const std::string& id, const std::string& desc, TopoDS_Shape shape,
                        double deflection, double angDeflection, Sidecar sidecar) {
    tessellate(shape, deflection, angDeflection);
    FixtureResult out;
    out.mesh = weldQuantized(quantizeMesh(extractMesh(shape)));
    out.sidecar = std::move(sidecar);
    out.sidecar.id = id;
    out.sidecar.description = desc;
    out.sidecar.deflection = deflection;
    fillMeshSidecar(out);
    return out;
}

// Fixtures whose STL+sidecar are committed (see tests/corpus/.gitignore).
// OCCT tessellation of the S09 loft differs across libm/FMA targets; the
// committed pair is the cross-platform calibration input for P1/P2 gates.
static bool isPinnedCorpusFixture(const std::string& id) {
    return id == "S09";
}

bool writeFixture(const fs::path& dir, const FixtureResult& fx) {
    const fs::path stlPath = dir / (fx.sidecar.id + ".stl");
    const fs::path jsonPath = dir / (fx.sidecar.id + ".expected.json");
    if (isPinnedCorpusFixture(fx.sidecar.id) && fs::exists(stlPath) && fs::exists(jsonPath)) {
        return true;
    }
    const std::string label = ("stl2step corpus " + fx.sidecar.id).substr(0, 79);
    if (!writeBinaryStl(stlPath.string(), fx.mesh, label.c_str())) return false;
    return writeTextFile(jsonPath.string(), writeSidecarJson(fx.sidecar));
}

// ---- S01: 10 mm cube -------------------------------------------------------
FixtureResult buildS01() {
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1)};
    return emitShape("S01", "10 mm cube — six planes", BRepPrimAPI_MakeBox(10, 10, 10), 0.1, 0.5,
                     sc);
}

// ---- S02: filleted cube R=2 --------------------------------------------------
FixtureResult buildS02() {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(20, 20, 20);
    BRepFilletAPI_MakeFillet fillet(box);
    TopExp_Explorer exp(box, TopAbs_EDGE);
    for (; exp.More(); exp.Next()) fillet.Add(2.0, TopoDS::Edge(exp.Current()));
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1)};
    auto cyls = cubeEdgeCylinders(20.0, 2.0);
    sc.recoverable.insert(sc.recoverable.end(), cyls.begin(), cyls.end());
    sc.mustRemainFaceted = {{ "sphere", 8, "vertex blends at cube corners (SphereNYI)" }};
    sc.expectedRejects = {"sphereNYI"};
    return emitShape("S02", "20 mm cube, R=2 mm edge fillets + 8 spherical vertex blends", fillet,
                     0.25, 0.7, sc);
}

// ---- S03: plate with cylindrical holes ---------------------------------------
FixtureResult buildS03() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(100, 60, 10);
    const double holes[][3] = {{20, 15, 5}, {80, 15, 5}, {20, 45, 5}, {80, 45, 5}};
    for (const auto& h : holes) {
        gp_Ax2 ax(gp_Pnt(h[0], h[1], h[2]), gp_Dir(0, 0, 1));
        TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(ax, 4.0, 12.0);
        gp_Trsf tr;
        tr.SetTranslation(gp_Vec(0, 0, -1.0));
        cyl = BRepBuilderAPI_Transform(cyl, tr).Shape();
        plate = BRepAlgoAPI_Cut(plate, cyl);
    }
    {
        gp_Ax2 ax(gp_Pnt(50, 30, 5), gp_Dir(0, 0, 1));
        TopoDS_Shape cone = BRepPrimAPI_MakeCone(ax, 3.0, 2.5, 12.0);
        gp_Trsf tr;
        tr.SetTranslation(gp_Vec(0, 0, -1.0));
        cone = BRepBuilderAPI_Transform(cone, tr).Shape();
        plate = BRepAlgoAPI_Cut(plate, cone);
    }
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1),
                      planeRec({1, 0, 0}, 1), planeRec({-1, 0, 0}, 1),
                      planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1),
                      cylRec(4.0, {20, 15, 5}, {0, 0, 1}, 1, 0, true),
                      cylRec(4.0, {80, 15, 5}, {0, 0, 1}, 1, 0, true),
                      cylRec(4.0, {20, 45, 5}, {0, 0, 1}, 1, 0, true),
                      cylRec(4.0, {80, 45, 5}, {0, 0, 1}, 1, 0, true)};
    // D-130-15(2): `expectedRejects: coneNYI` predates the cone (D-130-3);
    // retired by that ruling. The 5 deg drafted hole is outside detector C's
    // chamfer window, so whatever reject it carries is measured, not expected.
    return emitShape("S03", "100x60x10 plate with four R=4 holes + 5 deg drafted cone hole", plate,
                     0.2, 0.5, sc);
}

// ---- S04: cylinder boss on plate ---------------------------------------------
FixtureResult buildS04() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(50, 50, 10);
    gp_Ax2 ax(gp_Pnt(25, 25, 10), gp_Dir(0, 0, 1));
    TopoDS_Shape boss = BRepPrimAPI_MakeCylinder(ax, 10.0, 15.0);
    BRepFilletAPI_MakeFillet fillet(boss);
    TopExp_Explorer exp(boss, TopAbs_EDGE);
    for (; exp.More(); exp.Next()) {
        TopoDS_Edge e = TopoDS::Edge(exp.Current());
        fillet.Add(3.0, e);
        break;
    }
    TopoDS_Shape fused = BRepAlgoAPI_Fuse(plate, fillet.Shape());
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 2), planeRec({0, 0, -1}, 1),
                      planeRec({1, 0, 0}, 1), planeRec({-1, 0, 0}, 1),
                      planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1),
                      cylRec(10.0, {25, 25, 10}, {0, 0, 1}, 1, 0, true)};
    sc.mustRemainFaceted = {{ "torus", 1, "boss top blend annulus (TorusNYI)" }};
    return emitShape("S04", "50x50 plate + R=10 boss H=15 with torus top blend", fused, 0.2, 0.4,
                     sc);
}

// ---- S05: slot with rounded ends ---------------------------------------------
FixtureResult buildS05() {
    TopoDS_Shape block = BRepPrimAPI_MakeBox(40, 20, 10);
    gp_Ax2 ax1(gp_Pnt(10, 10, 0), gp_Dir(0, 0, 1));
    gp_Ax2 ax2(gp_Pnt(30, 10, 0), gp_Dir(0, 0, 1));
    TopoDS_Shape c1 = BRepPrimAPI_MakeCylinder(ax1, 5.0, 10.0);
    TopoDS_Shape c2 = BRepPrimAPI_MakeCylinder(ax2, 5.0, 10.0);
    gp_Ax2 slotAx(gp_Pnt(10, 5, 0), gp_Dir(0, 0, 1));
    TopoDS_Shape slotTool = BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(slotAx, 20, 10, 10), c1);
    slotTool = BRepAlgoAPI_Fuse(slotTool, c2);
    TopoDS_Shape part = BRepAlgoAPI_Cut(block, slotTool);
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 2), planeRec({0, -1, 0}, 2),
                      cylRec(5.0, {10, 10, 0}, {0, 0, 1}, 1, 0, false),
                      cylRec(5.0, {30, 10, 0}, {0, 0, 1}, 1, 0, false)};
    return emitShape("S05", "40x20x10 block with slot R=5 rounded ends (~180 partials)", part, 0.15,
                     0.5, sc);
}

// ---- S06-S08: exact-N polygonal cylinders ------------------------------------
FixtureResult buildNGonCylinder(const char* id, const char* desc, int N, double defl) {
    gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    TopoDS_Shape cyl = makeNGonPrism(25.0, 100.0, N, ax);
    Sidecar sc;
    sc.recoverable = {cylRec(25.0, {0, 0, 0}, {0, 0, 1}, 1, 0, true),
                      planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1)};
    return emitShape(id, desc, cyl, defl, 45.0, sc);
}

// ---- S09: box + loft freeform ------------------------------------------------
FixtureResult buildS09() {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(40, 40, 20);
    BRepOffsetAPI_ThruSections loft(Standard_True);
    BRepBuilderAPI_MakeWire w1, w2;
    w1.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(40, 10, 20), gp_Pnt(50, 10, 30)).Edge());
    w1.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(50, 10, 30), gp_Pnt(50, 30, 30)).Edge());
    w1.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(50, 30, 30), gp_Pnt(40, 30, 20)).Edge());
    w1.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(40, 30, 20), gp_Pnt(40, 10, 20)).Edge());
    // Inner profile inset in +x so the loft skin does not tessellate a zero-area
    // strip on the shared x=40 seam (adjudicate-p2-real S09c1).
    w2.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(41, 12, 20), gp_Pnt(48, 12, 35)).Edge());
    w2.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(48, 12, 35), gp_Pnt(48, 28, 35)).Edge());
    w2.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(48, 28, 35), gp_Pnt(41, 28, 20)).Edge());
    w2.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(41, 28, 20), gp_Pnt(41, 12, 20)).Edge());
    loft.AddWire(w1.Wire());
    loft.AddWire(w2.Wire());
    loft.Build();
    // Compound (not Fuse): the loft only touches the box along an edge; Fuse
    // tessellates zero-area seam triangles on x=40 that are open==0 yet not
    // removable without opening the shell (adjudicate-p2-real S09c1).
    TopoDS_Compound assembly;
    BRep_Builder builder;
    builder.MakeCompound(assembly);
    builder.Add(assembly, box);
    builder.Add(assembly, loft.Shape());
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, -1}, 1), planeRec({0, 0, 1}, 1), planeRec({-1, 0, 0}, 1),
                      planeRec({0, -1, 0}, 1), planeRec({1, 0, 0}, 1), planeRec({0, 1, 0}, 1)};
    sc.mustRemainFaceted = {{ "freeform", 1, "ThruSections loft patch" }};
    sc.expectedSolids = 2;  // fuse of box+loft does not merge to one body
    return emitShape("S09", "40x40x20 box fused with spline loft (freeform must stay faceted; 2 solids)",
                     assembly, 0.25, 0.5, sc);
}

// ---- S10: thin-wall tube -----------------------------------------------------
FixtureResult buildS10() {
    gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    TopoDS_Shape outer = BRepPrimAPI_MakeCylinder(ax, 20.0, 30.0);
    TopoDS_Shape inner = BRepPrimAPI_MakeCylinder(ax, 17.0, 32.0);
    gp_Trsf tr;
    tr.SetTranslation(gp_Vec(0, 0, -1.0));
    inner = BRepBuilderAPI_Transform(inner, tr).Shape();
    TopoDS_Shape tube = BRepAlgoAPI_Cut(outer, inner);
    Sidecar sc;
    sc.recoverable = {cylRec(20.0, {0, 0, 0}, {0, 0, 1}, 1, 0, true),
                      cylRec(17.0, {0, 0, 0}, {0, 0, 1}, 1, 0, true),
                      planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1)};
    return emitShape("S10", "thin-wall tube outer R=20 inner R=17 H=30", tube, 0.15, 15.0, sc);
}

// ---- S11: asymmetric chamfer -------------------------------------------------
FixtureResult buildS11() {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(30, 30, 30);
    TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
    TopExp::MapShapesAndAncestors(box, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);
    BRepFilletAPI_MakeChamfer chamfer(box);
    // Unequal setbacks on the SAME edge: Add(d1, d2, edge, face), d1/d2 = 1.3.
    const double d1 = 3.9;
    const double d2 = 3.0;
    int added = 0;
    for (int i = 1; i <= edgeToFaces.Extent() && added < 4; ++i) {
        const TopoDS_Edge e = TopoDS::Edge(edgeToFaces.FindKey(i));
        const TopTools_ListOfShape& faces = edgeToFaces(i);
        if (faces.IsEmpty()) continue;
        chamfer.Add(d1, d2, e, TopoDS::Face(faces.First()));
        ++added;
    }
    Sidecar sc;
    sc.expectedRejects = {"filletConsensus"};
    FixtureResult fx = emitShape(
        "S11", "30 mm cube asymmetric chamfer sL/sR=1.3 on same edge (not a fillet strip)", chamfer,
        0.15, 0.5, sc);
    fx.sidecar.recoverable.clear();
    appendMeshPlanes(fx.sidecar, fx.mesh);
    return fx;
}

// ---- S11-b: symmetric 45 deg chamfer -----------------------------------------
FixtureResult buildS11b() {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(25, 25, 25);
    BRepFilletAPI_MakeChamfer chamfer(box);
    TopExp_Explorer exp(box, TopAbs_EDGE);
    for (; exp.More(); exp.Next()) chamfer.Add(2.0, TopoDS::Edge(exp.Current()));
    Sidecar sc;
    sc.expectFillet = true;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1)};
    auto cyls = cubeEdgeCylinders(25.0, 2.0);
    sc.recoverable.insert(sc.recoverable.end(), cyls.begin(), cyls.end());
    return emitShape("S11-b", "25 mm cube symmetric 45 deg chamfer (expect fillet recovery)",
                     chamfer, 0.2, 0.5, sc);
}

// ---- S12: two-body assembly --------------------------------------------------
FixtureResult buildS12() {
    TopoDS_Shape a = BRepPrimAPI_MakeBox(10, 10, 10);
    gp_Trsf tr;
    tr.SetTranslation(gp_Vec(15, 0, 0));
    TopoDS_Shape b = BRepBuilderAPI_Transform(BRepPrimAPI_MakeBox(10, 10, 10), tr).Shape();
    TopoDS_Compound comp;
    BRep_Builder builder;
    builder.MakeCompound(comp);
    builder.Add(comp, a);
    builder.Add(comp, b);
    Sidecar sc;
    sc.recoverable = {planeRec({0, 0, 1}, 2), planeRec({0, 0, -1}, 2), planeRec({1, 0, 0}, 2),
                      planeRec({-1, 0, 0}, 2), planeRec({0, 1, 0}, 2), planeRec({0, -1, 0}, 2)};
    sc.expectedSolids = 2;
    return emitShape("S12", "two 10 mm cubes 5 mm apart (two solids)", comp, 0.1, 0.5, sc);
}

// ---- S12-a: square + pentagonal holes ----------------------------------------
FixtureResult buildS12a() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(80, 50, 8);
    gp_Trsf tr;
    tr.SetTranslation(gp_Vec(0, 0, -1.0));
    {
        gp_Ax2 ax4(gp_Pnt(20, 25, 0), gp_Dir(0, 0, 1));
        TopoDS_Shape h4 = makeNGonPrism(8.0, 10.0, 4, ax4);
        h4 = BRepBuilderAPI_Transform(h4, tr).Shape();
        plate = BRepAlgoAPI_Cut(plate, h4);
        gp_Ax2 ax5(gp_Pnt(60, 25, 0), gp_Dir(0, 0, 1));
        TopoDS_Shape h5 = makeNGonPrism(7.0, 10.0, 5, ax5);
        h5 = BRepBuilderAPI_Transform(h5, tr).Shape();
        plate = BRepAlgoAPI_Cut(plate, h5);
    }
    Sidecar sc;
    FixtureResult fx = emitShape(
        "S12-a",
        "plate with N=4 square + N=5 pentagonal holes (0 cylinders; 4+5 hole-side planes + plate)",
        plate, 0.1, 0.5, sc);
    fx.sidecar.recoverable.clear();
    appendMeshPlanes(fx.sidecar, fx.mesh);
    return fx;
}

// ---- S12-b: non-concyclic hex hole -------------------------------------------
FixtureResult buildS12b() {
    TopoDS_Shape plate = makeNonConcyclicHexHole(70, 50, 8, 10.0, gp_Pnt(35, 25, 0));
    Sidecar sc;
    sc.expectedRejects = {"vertexResidual"};
    FixtureResult fx =
        emitShape("S12-b", "non-concyclic hex hole R/1.15R alternating (G2 negative; 6 hole planes)",
                  plate, 0.1, 0.5, sc);
    fx.sidecar.recoverable.clear();
    appendMeshPlanes(fx.sidecar, fx.mesh);
    return fx;
}

// ---- S13: 1-band (1-quad) quarter-round fillet --------------------------------
FixtureResult buildS13() {
    // Planar 1-band chord between the two tangent generators. OCCT analytic
    // fillet + IncrementalMesh produced 3 bands even at 90° angular deflection.
    TopoDS_Shape shape = makeBandedQuarterRoundBox(15.0, 2.0, 1);
    Sidecar sc;
    sc.expectFillet = true;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1),
                      cylRec(2.0, {0, 2, 2}, {1, 0, 0}, 1, 0, false)};
    return emitShape("S13", "1-band quarter-round fillet strip (D7: vertices on tangent planes)",
                     shape, 0.5, 90.0, sc);
}

// ---- S14: 2-band fillet strip ------------------------------------------------
FixtureResult buildS14() {
    TopoDS_Shape shape = makeBandedQuarterRoundBox(18.0, 2.5, 2);
    Sidecar sc;
    sc.expectFillet = true;
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1),
                      cylRec(2.5, {0, 2.5, 2.5}, {1, 0, 0}, 1, 0, false)};
    return emitShape("S14", "2-band fillet strip (two planar rows spanning 90 deg)", shape, 0.35,
                     45.0, sc);
}

// ---- S15: dirty mesh that sews closed ----------------------------------------
FixtureResult buildS15() {
    FixtureResult clean = buildS01();
    MeshData dirty = clean.mesh;
    // Duplicate one vertex with a gap BELOW auto sewTol (diag*1e-5 ≈ 1.7e-4 mm
    // on a 10 mm cube). 5e-5 mm opens an edge the sew pass must close.
    if (!dirty.tris.empty()) {
        const int dupFrom = dirty.tris[0][1];
        const int dupIdx = static_cast<int>(dirty.verts.size());
        dirty.verts.push_back(dirty.verts[dupFrom]);
        dirty.verts.back().x += 5e-5;
        dirty.verts.back().y += 5e-5;
        dirty.tris[0][1] = dupIdx;
    }
    Sidecar sc = clean.sidecar;
    sc.id = "S15";
    sc.description =
        "dirty 10 mm cube (split shared vertex 5e-5 mm; sews closed under auto sewTol)";
    sc.expectedExit = 0;
    sc.expectedSolids = 1;
    sc.expectedOpenShells = 0;
    FixtureResult out;
    out.mesh = std::move(dirty);
    out.sidecar = std::move(sc);
    fillMeshSidecar(out);
    return out;
}

// Closed wedge plate: top plane slants from zThick at x=0 to zThin (>0) at x=X.
// A vertical N>=6 bore through non-parallel caps is still a closed360 P1 claim
// (45° ridges) but P2 Seamed360 expects iso-v circular caps — slanted elliptical
// caps fail MakeFace → R1 explode. zThin must stay positive so the solid is
// a simple closed trapezoid (negative zThin tore the mesh).
TopoDS_Shape makeWedgePlate(double X, double Y, double zThick, double zThin) {
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(0, 0, 0));
    poly.Add(gp_Pnt(X, 0, 0));
    poly.Add(gp_Pnt(X, 0, zThin));
    poly.Add(gp_Pnt(0, 0, zThick));
    poly.Close();
    return BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(poly.Wire()), gp_Vec(0, Y, 0));
}

TopoDS_Shape cutNGonHole(const TopoDS_Shape& plate, double x, double y, double R, int N,
                         double zLo, double zHi) {
    gp_Ax2 ax(gp_Pnt(x, y, zLo), gp_Dir(0, 0, 1));
    TopoDS_Shape hole = makeNGonPrism(R, zHi - zLo, N, ax);
    return BRepAlgoAPI_Cut(plate, hole);
}

void filletSomeEdges(TopoDS_Shape& shape, double radius, int maxEdges) {
    BRepFilletAPI_MakeFillet fillet(shape);
    TopExp_Explorer exp(shape, TopAbs_EDGE);
    int n = 0;
    for (; exp.More() && n < maxEdges; exp.Next()) {
        fillet.Add(radius, TopoDS::Edge(exp.Current()));
        ++n;
    }
    if (n == 0) return;
    fillet.Build();
    if (fillet.IsDone()) shape = fillet.Shape();
}

// ---- S16 sub-fixtures --------------------------------------------------------
FixtureResult buildS16R1Explode() {
    // N=8 concyclic bore (45° ridges seed a closed360 claim) through a closed
    // wedge. Slanted caps are not iso-v circles → P2 MakeFace fails → R1 explode.
    // Sidecar is the CORRECT STEP outcome: hole faceted, no recoverable
    // closed360 cylinder (G2.5 HARD).
    const double X = 80, Y = 60, zThick = 12, zThin = 3;
    TopoDS_Shape plate = makeWedgePlate(X, Y, zThick, zThin);
    plate = cutNGonHole(plate, X * 0.5, Y * 0.5, 12.0, 8, -2.0, zThick + 2.0);
    Sidecar sc;
    sc.rLadderRole = "R1-explode-success";
    sc.mustRemainFaceted = {{"hole", 1,
                             "N=8 closed360 claim, slanted caps fail Seamed360 MakeFace; R1 explode"}};
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1)};
    return emitShape("S16-R1-explode-success",
                     "N=8 bore through closed wedge: P1 closed360, slanted caps fail MakeFace, R1 explode; no Seamed360 cylinder",
                     plate, 0.15, 45.0, sc);
}

FixtureResult buildS16R1Round2() {
    // Twin N=8 pinched bores sharing the slanted cap + rim fillets. Exploding
    // the first rebuilds the shared plane/fillets; the second then fails
    // MakeFace (round 2) and the cascade stabilises. No recoverable cylinder:
    // both holes explode.
    const double X = 90, Y = 50, zThick = 12, zThin = 3;
    TopoDS_Shape plate = makeWedgePlate(X, Y, zThick, zThin);
    plate = cutNGonHole(plate, X * 0.5, 16.0, 9.0, 8, -2.0, zThick + 2.0);
    plate = cutNGonHole(plate, X * 0.5, 34.0, 9.0, 8, -2.0, zThick + 2.0);
    filletSomeEdges(plate, 1.2, 4);
    Sidecar sc;
    sc.rLadderRole = "R1-round-2";
    sc.mustRemainFaceted = {
        {"hole", 2, "twin N=8 wedge bores + rim fillets: two R1 explodes, then stable"}};
    sc.recoverable = {planeRec({0, 0, 1}, 1), planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1), planeRec({0, -1, 0}, 1)};
    return emitShape("S16-R1-round-2",
                     "twin N=8 wedge bores + fillets: R1 cascade uses both rounds; no recoverable cylinder",
                     plate, 0.15, 0.5, sc);
}

FixtureResult buildS16R2() {
    // Three N=8 pinched bores + fillets: explode cascade exceeds the 2-round
    // cap → ChainUnstable/J6 → whole-component R2. Correct run produces no
    // analytic faces, so recoverable is empty.
    const double X = 80, Y = 50, zThick = 12, zThin = 3;
    TopoDS_Shape plate = makeWedgePlate(X, Y, zThick, zThin);
    plate = cutNGonHole(plate, X * 0.5, 10.0, 7.0, 8, -2.0, zThick + 2.0);
    plate = cutNGonHole(plate, X * 0.5, 25.0, 7.0, 8, -2.0, zThick + 2.0);
    plate = cutNGonHole(plate, X * 0.5, 40.0, 7.0, 8, -2.0, zThick + 2.0);
    filletSomeEdges(plate, 1.0, 6);
    Sidecar sc;
    sc.rLadderRole = "R2-ChainUnstable";
    sc.expectedRejects = {"chainUnstable"};
    sc.mustRemainFaceted = {
        {"component", 1, "R2 reverts the whole component to faceted after ChainUnstable"}};
    FixtureResult fx = emitShape(
        "S16-R2-ChainUnstable",
        "three N=8 wedge bores + fillets: R1 cascade exceeds 2-round cap → ChainUnstable → R2; no cylinder",
        plate, 0.2, 0.5, sc);
    fx.sidecar.recoverable.clear();
    return fx;
}

// ---- 130 battery B2–B6 (synthetic; N=48 regular polylines) -------------------

TopoDS_Edge findBoxEdgeAlongXAtYZ(const TopoDS_Shape& box, double y, double z) {
    TopExp_Explorer exp(box, TopAbs_EDGE);
    for (; exp.More(); exp.Next()) {
        const TopoDS_Edge e = TopoDS::Edge(exp.Current());
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(e, v1, v2);
        const gp_Pnt p1 = BRep_Tool::Pnt(v1);
        const gp_Pnt p2 = BRep_Tool::Pnt(v2);
        const bool atY = std::fabs(p1.Y() - y) < 1e-7 && std::fabs(p2.Y() - y) < 1e-7;
        const bool atZ = std::fabs(p1.Z() - z) < 1e-7 && std::fabs(p2.Z() - z) < 1e-7;
        const bool alongX = std::fabs(p1.X() - p2.X()) > 1.0;
        if (atY && atZ && alongX) return e;
    }
    return {};
}

// External boss: N-gon cylinder then 45° conical chamfer (ruled frustum).
TopoDS_Shape makeNGonBossChamfer(double R, double Hcyl, double chamfer, int N, const gp_Ax2& ax) {
    BRepOffsetAPI_ThruSections loft(Standard_True, Standard_True);
    loft.AddWire(nGonWire(R, N, ax, 0.0));
    loft.AddWire(nGonWire(R, N, ax, Hcyl));
    loft.AddWire(nGonWire(R - chamfer, N, ax, Hcyl + chamfer));
    loft.Build();
    return loft.Shape();
}

void tagBattery130(Sidecar& sc) {
    sc.battery = "130";
    sc.expectedExit = 0;
    sc.expectedSolids = 1;
    sc.expectedOpenShells = 0;
    sc.smoothExpectedExit = -1;  // gate_130 owns the smooth contract
}

// B2: 60×40×30 block, R=8 through Z and R=5 through X, crossing at the centre.
FixtureResult buildCrossBores() {
    TopoDS_Shape block = BRepPrimAPI_MakeBox(60, 40, 30);
    block = cutThroughNGon(block, gp_Pnt(30, 20, 0), gp_Dir(0, 0, 1), 8.0, 30.0, kTessN);
    block = cutThroughNGon(block, gp_Pnt(0, 20, 15), gp_Dir(1, 0, 0), 5.0, 60.0, kTessN);
    Sidecar sc;
    tagBattery130(sc);
    sc.recoverable = boxSixPlanes();
    sc.recoverable.push_back(cylRec(8.0, {30, 20, 0}, {0, 0, 1}, 1, 0, true));
    sc.recoverable.push_back(cylRec(5.0, {0, 20, 15}, {1, 0, 0}, 1, 0, true));
    sc.intersections = {{"cyl R8", "cyl R5", "cylcyl"}};
    return emitShape("cross_bores",
                     "60x40x30 block, through-bore R=8 along Z × R=5 along X at centre (cyl∩cyl)",
                     block, 0.2, 0.5, sc);
}

// B3: 60×60×10 plate, R=10 hole with 45°×2 mm conical entry chamfer, plus
// one straight 45°×3 mm bevel. Straight bevel is a plane (7 planes total).
FixtureResult buildChamferStraightRing() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(60, 60, 10);
    const TopoDS_Edge bevelEdge = findBoxEdgeAlongXAtYZ(plate, 0.0, 10.0);
    if (bevelEdge.IsNull())
        throw std::runtime_error("chamfer_straight_ring: top-front edge not found");
    BRepFilletAPI_MakeChamfer chamfer(plate);
    chamfer.Add(3.0, bevelEdge);
    chamfer.Build();
    if (!chamfer.IsDone())
        throw std::runtime_error("chamfer_straight_ring: MakeChamfer failed");
    plate = chamfer.Shape();
    plate = cutThroughNGon(plate, gp_Pnt(30, 30, 0), gp_Dir(0, 0, 1), 10.0, 10.0, kTessN);
    // 45° conical entry: R 10 at z=8 → R 12 at z=10. Overshoot the top 0.1 mm.
    gp_Ax2 coneAx(gp_Pnt(30, 30, 8), gp_Dir(0, 0, 1));
    plate = BRepAlgoAPI_Cut(plate, makeNGonFrustum(10.0, 12.1, 2.1, kTessN, coneAx));
    Sidecar sc;
    tagBattery130(sc);
    sc.recoverable = boxSixPlanes();
    sc.recoverable.push_back(planeRec(normalize(Vec3(0, -1, 1)), 1));  // 45° straight bevel
    sc.recoverable.push_back(cylRec(10.0, {30, 30, 0}, {0, 0, 1}, 1, 0, true));
    sc.recoverable.push_back(
        coneRec(10.0, 12.0, 45.0, {30, 30, 8}, {0, 0, 1}, 1, kTessN));
    FixtureResult fx = emitShape(
        "chamfer_straight_ring",
        "60x60x10 plate, R=10 hole with 45deg x 2mm conical chamfer + straight 45deg x 3mm bevel",
        plate, 0.2, 0.5, sc);
    // SPEC: write the exact plane count the generator produces (expect 7).
    int planes = 0, cyls = 0, cones = 0;
    sumRecoverableCensus(fx.sidecar, planes, cyls, cones);
    fx.sidecar.description += " (planes=" + std::to_string(planes) + ")";
    return fx;
}

// B4: 60×60×10 plate + external boss R=12 H=15 with 45°×2 mm top chamfer (R 12→10).
FixtureResult buildBossConeChamfer() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(60, 60, 10);
    gp_Ax2 ax(gp_Pnt(30, 30, 10), gp_Dir(0, 0, 1));
    TopoDS_Shape boss = makeNGonBossChamfer(12.0, 13.0, 2.0, kTessN, ax);
    TopoDS_Shape part = BRepAlgoAPI_Fuse(plate, boss);
    Sidecar sc;
    tagBattery130(sc);
    // SPEC GT: plate 6 + boss top 1 = 7 (two +Z patches: plate annulus + boss cap).
    sc.recoverable = {planeRec({0, 0, 1}, 2),  planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1),  planeRec({0, -1, 0}, 1),
                      cylRec(12.0, {30, 30, 10}, {0, 0, 1}, 1, 0, true),
                      coneRec(12.0, 10.0, 45.0, {30, 30, 23}, {0, 0, 1}, 1, kTessN)};
    return emitShape("boss_cone_chamfer",
                     "60x60x10 plate + external boss R=12 H=15 with 45deg x 2mm conical top chamfer "
                     "(R 12→10)",
                     part, 0.2, 0.5, sc);
}

// B5: 60×60×20 plate, through R=6, counterbore R=10 depth 6, 45°×1.5 mm entry chamfer.
FixtureResult buildCounterboreChamfer() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(60, 60, 20);
    plate = cutThroughNGon(plate, gp_Pnt(30, 30, 0), gp_Dir(0, 0, 1), 6.0, 20.0, kTessN);
    // Counterbore from z=14 to above the top.
    gp_Ax2 cbAx(gp_Pnt(30, 30, 14), gp_Dir(0, 0, 1));
    plate = BRepAlgoAPI_Cut(plate, makeNGonPrism(10.0, 7.0, kTessN, cbAx));
    gp_Ax2 coneAx(gp_Pnt(30, 30, 18.5), gp_Dir(0, 0, 1));
    plate = BRepAlgoAPI_Cut(plate, makeNGonFrustum(10.0, 11.6, 1.6, kTessN, coneAx));
    Sidecar sc;
    tagBattery130(sc);
    sc.recoverable = {planeRec({0, 0, 1}, 2),  planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1),  planeRec({0, -1, 0}, 1),
                      cylRec(6.0, {30, 30, 0}, {0, 0, 1}, 1, 0, true),
                      cylRec(10.0, {30, 30, 14}, {0, 0, 1}, 1, 0, true),
                      coneRec(10.0, 11.5, 45.0, {30, 30, 18.5}, {0, 0, 1}, 1, kTessN)};
    return emitShape("counterbore_chamfer",
                     "60x60x20 plate, through R=6, counterbore R=10 depth 6, 45deg x 1.5mm "
                     "conical entry (R 10→11.5)",
                     plate, 0.2, 0.5, sc);
}

// B6: B5 plus a through-bore R=4 along X at z=17 (cuts counterbore wall + chamfer cone).
FixtureResult buildCylMeetsChamfer() {
    TopoDS_Shape plate = BRepPrimAPI_MakeBox(60, 60, 20);
    plate = cutThroughNGon(plate, gp_Pnt(30, 30, 0), gp_Dir(0, 0, 1), 6.0, 20.0, kTessN);
    gp_Ax2 cbAx(gp_Pnt(30, 30, 14), gp_Dir(0, 0, 1));
    plate = BRepAlgoAPI_Cut(plate, makeNGonPrism(10.0, 7.0, kTessN, cbAx));
    gp_Ax2 coneAx(gp_Pnt(30, 30, 18.5), gp_Dir(0, 0, 1));
    plate = BRepAlgoAPI_Cut(plate, makeNGonFrustum(10.0, 11.6, 1.6, kTessN, coneAx));
    plate = cutThroughNGon(plate, gp_Pnt(0, 30, 17), gp_Dir(1, 0, 0), 4.0, 60.0, kTessN);
    Sidecar sc;
    tagBattery130(sc);
    sc.recoverable = {planeRec({0, 0, 1}, 2),  planeRec({0, 0, -1}, 1), planeRec({1, 0, 0}, 1),
                      planeRec({-1, 0, 0}, 1), planeRec({0, 1, 0}, 1),  planeRec({0, -1, 0}, 1),
                      cylRec(6.0, {30, 30, 0}, {0, 0, 1}, 1, 0, true),
                      cylRec(10.0, {30, 30, 14}, {0, 0, 1}, 1, 0, true),
                      cylRec(4.0, {0, 30, 17}, {1, 0, 0}, 1, 0, true),
                      coneRec(10.0, 11.5, 45.0, {30, 30, 18.5}, {0, 0, 1}, 1, kTessN)};
    sc.intersections = {{"cyl R10", "cyl R4", "cylcyl"},
                        {"cyl R6", "cyl R4", "cylcyl"},
                        {"cone R10-11.5", "cyl R4", "conecyl"}};
    return emitShape("cyl_meets_chamfer",
                     "counterbore_chamfer + through-bore R=4 along X at z=17 (cyl∩cyl and cone∩cyl)",
                     plate, 0.2, 0.5, sc);
}

}  // namespace

int main(int argc, char** argv) {
    fs::path outDir = fs::current_path();
    if (argc >= 2) outDir = argv[1];
    fs::create_directories(outDir);

    std::vector<FixtureResult> fixtures;
    int failures = 0;
    const auto run = [&](const char* label, auto&& fn) {
        try {
            fixtures.push_back(fn());
        } catch (const Standard_Failure& e) {
            std::cerr << "OCCT FAIL " << label << ": " << e.GetMessageString() << "\n";
            ++failures;
        } catch (const std::exception& e) {
            std::cerr << "FAIL " << label << ": " << e.what() << "\n";
            ++failures;
        }
    };
    run("S01", [] { return buildS01(); });
    run("S02", [] { return buildS02(); });
    run("S03", [] { return buildS03(); });
    run("S04", [] { return buildS04(); });
    run("S05", [] { return buildS05(); });
    run("S06", [] { return buildNGonCylinder("S06", "N=8 polygonal cylinder D=50 H=100", 8, 0.05); });
    run("S07", [] { return buildNGonCylinder("S07", "N=12 polygonal cylinder D=50 H=100", 12, 0.05); });
    run("S08", [] { return buildNGonCylinder("S08", "N=16 polygonal cylinder D=50 H=100", 16, 0.05); });
    if (!isPinnedCorpusFixture("S09")
        || !fs::exists(outDir / "S09.stl") || !fs::exists(outDir / "S09.expected.json"))
        run("S09", [] { return buildS09(); });
    run("S10", [] { return buildS10(); });
    run("S11", [] { return buildS11(); });
    run("S11-b", [] { return buildS11b(); });
    run("S12", [] { return buildS12(); });
    run("S12-a", [] { return buildS12a(); });
    run("S12-b", [] { return buildS12b(); });
    run("S13", [] { return buildS13(); });
    run("S14", [] { return buildS14(); });
    run("S15", [] { return buildS15(); });
    run("S16-R1-explode-success", [] { return buildS16R1Explode(); });
    run("S16-R1-round-2", [] { return buildS16R1Round2(); });
    run("S16-R2-ChainUnstable", [] { return buildS16R2(); });
    run("cross_bores", [] { return buildCrossBores(); });
    run("chamfer_straight_ring", [] { return buildChamferStraightRing(); });
    run("boss_cone_chamfer", [] { return buildBossConeChamfer(); });
    run("counterbore_chamfer", [] { return buildCounterboreChamfer(); });
    run("cyl_meets_chamfer", [] { return buildCylMeetsChamfer(); });

    for (const auto& fx : fixtures) {
        if (isPinnedCorpusFixture(fx.sidecar.id)
            && fs::exists(outDir / (fx.sidecar.id + ".stl"))
            && fs::exists(outDir / (fx.sidecar.id + ".expected.json"))) {
            continue;
        }
        if (!writeFixture(outDir, fx)) {
            std::cerr << "FAIL write " << fx.sidecar.id << "\n";
            ++failures;
        } else {
            std::cout << "Wrote " << fx.sidecar.id << " tris=" << fx.sidecar.triangleCount
                      << " vol=" << fx.sidecar.meshVolume << "\n";
        }
    }
    return failures ? 1 : 0;
}
