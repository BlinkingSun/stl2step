// Cylinder fixture: prove the census distinguishes Seamed360 from TwoHalves.
//
// 1. BRepPrimAPI_MakeCylinder → one cylindrical face, U-span 2π (Seamed360).
// 2. The same geometry rebuilt as two 180° cylindrical faces (TwoHalves).
// Both are written to STEP and re-read by the independent census — the same
// path G2.5 will use.
//
// SPDX-License-Identifier: MIT

#include "step_census.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

namespace {

int gFails = 0;

void check(bool cond, const char* msg) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", msg);
        gFails++;
    }
}

void silenceOcct() {
    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);
}

double pi() { return std::acos(-1.0); }

bool writeStep(const TopoDS_Shape& shape, const std::string& path, std::string& err) {
    try {
        STEPControl_Writer w;
        const IFSelect_ReturnStatus ts = w.Transfer(shape, STEPControl_AsIs);
        if (ts != IFSelect_RetDone) {
            err = "STEP transfer failed";
            return false;
        }
        const IFSelect_ReturnStatus ws = w.Write(path.c_str());
        if (ws != IFSelect_RetDone) {
            err = "STEP write failed: " + path;
            return false;
        }
        return true;
    } catch (const Standard_Failure& e) {
        const char* m = e.GetMessageString();
        err = m && m[0] ? m : "OCCT exception during STEP write";
        return false;
    }
}

// Full cylinder lateral surface as two 180° faces sharing the same axis/radius.
// This is the TwoHalves rung: G2.5 must NOT accept it as a single hole.
TopoDS_Shape makeTwoHalves(double radius, double height) {
    const gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    const gp_Cylinder cyl(ax, radius);
    const double p = pi();
    const TopoDS_Face lat0 = BRepBuilderAPI_MakeFace(cyl, 0.0, p, 0.0, height).Face();
    const TopoDS_Face lat1 = BRepBuilderAPI_MakeFace(cyl, p, 2.0 * p, 0.0, height).Face();
    TopoDS_Compound comp;
    BRep_Builder b;
    b.MakeCompound(comp);
    b.Add(comp, lat0);
    b.Add(comp, lat1);
    return comp;
}

void assertSeamed360(const step_census::Result& r, double radius) {
    check(r.ok, "seamed360: census ok");
    check(r.cylinder == 1, "seamed360: exactly one cylindrical face");
    check(r.cylinders.size() == 1, "seamed360: cylinders[] length 1");
    if (!r.cylinders.empty()) {
        check(std::fabs(r.cylinders[0].radius - radius) < 1e-6, "seamed360: radius");
        check(r.cylinders[0].fullU, "seamed360: fullU true (U-span covers 2π)");
        check(r.cylinders[0].uSpan + 1e-4 >= 2.0 * pi(), "seamed360: uSpan >= 2π");
    }
    check(r.cylinderGroups.size() == 1, "seamed360: one cylinder group");
    if (!r.cylinderGroups.empty()) {
        check(r.cylinderGroups[0].nFaces == 1, "seamed360: group nFaces == 1");
        check(r.cylinderGroups[0].builtAs == "seamed360",
              "seamed360: builtAs == seamed360");
    }
}

void assertTwoHalves(const step_census::Result& r, double radius) {
    check(r.ok, "twoHalves: census ok");
    check(r.cylinder == 2, "twoHalves: exactly two cylindrical faces");
    check(r.cylinders.size() == 2, "twoHalves: cylinders[] length 2");
    int nFull = 0;
    for (const auto& c : r.cylinders) {
        check(std::fabs(c.radius - radius) < 1e-6, "twoHalves: radius");
        check(!c.fullU, "twoHalves: fullU false");
        check(std::fabs(c.uSpan - pi()) < 0.05, "twoHalves: uSpan ≈ π");
        if (c.fullU) nFull++;
    }
    check(nFull == 0, "twoHalves: no face spans 2π");
    check(r.cylinderGroups.size() == 1, "twoHalves: one cylinder group (same axis)");
    if (!r.cylinderGroups.empty()) {
        check(r.cylinderGroups[0].nFaces == 2, "twoHalves: group nFaces == 2");
        check(r.cylinderGroups[0].builtAs == "twoHalves",
              "twoHalves: builtAs == twoHalves");
    }
}

}  // namespace

int main(int argc, char** argv) {
    silenceOcct();
    const std::string outDir = (argc >= 2) ? argv[1] : ".";
    const double R = 3.0;
    const double H = 10.0;

    TopoDS_Shape seamed;
    try {
        BRepPrimAPI_MakeCylinder mk(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, H);
        seamed = mk.Solid();
        check(!seamed.IsNull(), "MakeCylinder produced a solid");
    } catch (const Standard_Failure& e) {
        std::fprintf(stderr, "FAIL: MakeCylinder threw: %s\n", e.GetMessageString());
        return 1;
    }

    const std::string seamedPath = outDir + "/seamed360.step";
    const std::string halvesPath = outDir + "/twohalves.step";
    std::string err;
    if (!writeStep(seamed, seamedPath, err)) {
        std::fprintf(stderr, "FAIL: write seamed360.step: %s\n", err.c_str());
        return 1;
    }

    const step_census::Result a1 = step_census::censusFile(seamedPath);
    const step_census::Result a2 = step_census::censusFile(seamedPath);
    assertSeamed360(a1, R);
    check(step_census::toJson(a1) == step_census::toJson(a2),
          "determinism: two censuses of seamed360.step are byte-identical");

    TopoDS_Shape halves;
    try {
        halves = makeTwoHalves(R, H);
        check(!halves.IsNull(), "two-halves shape produced");
    } catch (const Standard_Failure& e) {
        std::fprintf(stderr, "FAIL: two-halves construction threw: %s\n", e.GetMessageString());
        return 1;
    }
    if (!writeStep(halves, halvesPath, err)) {
        std::fprintf(stderr, "FAIL: write twohalves.step: %s\n", err.c_str());
        return 1;
    }

    const step_census::Result b = step_census::censusFile(halvesPath);
    assertTwoHalves(b, R);

    check(a1.cylinder != b.cylinder, "discrimination: cylinder face count differs");
    check(!a1.cylinders.empty() && !b.cylinders.empty(),
          "discrimination: both have cylinders[]");
    if (!a1.cylinders.empty() && !b.cylinders.empty()) {
        check(a1.cylinders[0].fullU != b.cylinders[0].fullU,
              "discrimination: fullU differs (2π vs π)");
    }
    check(!a1.cylinderGroups.empty() && !b.cylinderGroups.empty(),
          "discrimination: both have cylinderGroups[]");
    if (!a1.cylinderGroups.empty() && !b.cylinderGroups.empty()) {
        check(a1.cylinderGroups[0].builtAs != b.cylinderGroups[0].builtAs,
              "discrimination: builtAs seamed360 vs twoHalves");
        check(a1.cylinderGroups[0].nFaces != b.cylinderGroups[0].nFaces,
              "discrimination: nFaces 1 vs 2");
    }
    check(step_census::toJson(a1) != step_census::toJson(b),
          "discrimination: JSON of Seamed360 != JSON of TwoHalves");

    const step_census::Result missing = step_census::censusFile("/nonexistent.step");
    check(!missing.ok, "missing file: ok == false");
    check(!missing.error.empty(), "missing file: error string set");
    const std::string missingJson = step_census::toJson(missing);
    check(missingJson.find("\"error\"") != std::string::npos,
          "missing file: JSON contains error");
    check(missingJson.find("\"ok\":false") != std::string::npos,
          "missing file: JSON contains ok:false");

    if (gFails) {
        std::fprintf(stderr, "%d assertion(s) failed\n", gFails);
        std::fprintf(stderr, "seamed360 JSON: %s\n", step_census::toJson(a1).c_str());
        std::fprintf(stderr, "twoHalves JSON: %s\n", step_census::toJson(b).c_str());
        return 1;
    }
    std::printf("ok seamed360.builtAs=%s nFaces=%d twoHalves.builtAs=%s nFaces=%d\n",
                a1.cylinderGroups.empty() ? "?" : a1.cylinderGroups[0].builtAs.c_str(),
                a1.cylinderGroups.empty() ? -1 : a1.cylinderGroups[0].nFaces,
                b.cylinderGroups.empty() ? "?" : b.cylinderGroups[0].builtAs.c_str(),
                b.cylinderGroups.empty() ? -1 : b.cylinderGroups[0].nFaces);
    return 0;
}
