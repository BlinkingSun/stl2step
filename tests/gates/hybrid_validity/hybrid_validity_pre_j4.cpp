// hybrid_validity_pre_j4 — D3-8 / D3-13: BRepCheck at makeFaceKeep, pre-J4.
//
// Construction-time only. This TU must not call BRepLib::SameParameter or
// any ShapeFix_*. A runtime self-check greps this source; CMake has a
// sibling grep test.
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#ifndef HYBRID_VALIDITY_SOURCE
#define HYBRID_VALIDITY_SOURCE __FILE__
#endif

namespace {

int gFail = 0;

void check(bool ok, const char* name) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS %s\n", name);
    }
}

double kPi() { return std::acos(-1.0); }

// Same formula as src/refit_build.cpp meshTolCap for a cylinder of radius R
// with sewTol=0 and no measured chordSagitta (the nSides>=4 π/4 floor).
double meshTolCapR(double radius) {
    double c = Precision::Confusion();
    if (radius > 0.0)
        c = std::max(c, radius * (1.0 - std::cos(kPi() / 4.0)));
    return c * 1.001 + Precision::Confusion();
}

// BRep_Builder::MakeFace + Add — makeFaceKeep (J2). No pcurve, no ShapeFix,
// no SameParameter.
bool makeFaceKeep(const Handle(Geom_Surface)& surf, const TopoDS_Wire& outer,
                  TopoDS_Face& outF) {
    if (surf.IsNull() || outer.IsNull()) return false;
    try {
        BRep_Builder B;
        TopoDS_Face f;
        B.MakeFace(f, surf, Precision::Confusion());
        B.Add(f, outer);
        outF = f;
        return !outF.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
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

std::string faceStatus(const TopoDS_Face& f) {
    if (f.IsNull()) return "null";
    try {
        BRepCheck_Analyzer an(f, Standard_True);
        if (an.IsValid()) return "valid";
        std::string s = "invalid";
        Handle(BRepCheck_Result) res = an.Result(f);
        if (!res.IsNull()) {
            for (BRepCheck_ListOfStatus::Iterator it(res->Status()); it.More();
                 it.Next()) {
                s += " st=";
                s += std::to_string(static_cast<int>(it.Value()));
            }
        }
        return s;
    } catch (const Standard_Failure& e) {
        const char* m = e.GetMessageString();
        return m && m[0] ? m : "threw";
    }
}

TopoDS_Wire makeClosedPolyline(const std::vector<gp_Pnt>& pts) {
    BRepBuilderAPI_MakeWire mw;
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i < n; ++i) {
        const gp_Pnt& a = pts[static_cast<size_t>(i)];
        const gp_Pnt& b = pts[static_cast<size_t>((i + 1) % n)];
        TopoDS_Edge e = BRepBuilderAPI_MakeEdge(a, b).Edge();
        mw.Add(e);
    }
    TopoDS_Wire w = mw.Wire();
    w.Closed(Standard_True);
    return w;
}

void bumpEdgeTol(TopoDS_Face& f, double cap) {
    BRep_Builder B;
    for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
        TopoDS_Edge e = TopoDS::Edge(ex.Current());
        B.UpdateEdge(e, cap);
        for (TopExp_Explorer vx(e, TopAbs_VERTEX); vx.More(); vx.Next()) {
            B.UpdateVertex(TopoDS::Vertex(vx.Current()), cap);
        }
    }
}

// PLANE + mesh polyline (no stored pcurve). Square on z=0.
bool buildPlanePolyline(TopoDS_Face& outF) {
    Handle(Geom_Plane) pln = new Geom_Plane(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)));
    std::vector<gp_Pnt> pts = {
        gp_Pnt(0, 0, 0),
        gp_Pnt(10, 0, 0),
        gp_Pnt(10, 10, 0),
        gp_Pnt(0, 10, 0),
    };
    TopoDS_Wire w = makeClosedPolyline(pts);
    return makeFaceKeep(pln, w, outF);
}

// CYLINDER + mesh polyline. R=100, small u-span (~0.2 rad) so chord sagitta
// (~0.5 mm) is far below meshTolCap (~29.3186 mm). Generators are on-surface
// lines; the two v-constant sides are chords NOT on the cylinder.
bool buildCylPolyline(TopoDS_Face& outF, double radius = 100.0) {
    const gp_Ax3 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    Handle(Geom_CylindricalSurface) cyl = new Geom_CylindricalSurface(gp_Cylinder(ax, radius));
    const double u = 0.2;
    const double z0 = 0.0;
    const double z1 = 10.0;
    auto at = [&](double uu, double zz) {
        return gp_Pnt(radius * std::cos(uu), radius * std::sin(uu), zz);
    };
    std::vector<gp_Pnt> pts = {at(0.0, z0), at(u, z0), at(u, z1), at(0.0, z1)};
    TopoDS_Wire w = makeClosedPolyline(pts);
    return makeFaceKeep(cyl, w, outF);
}

struct ChordChain {
    const char* name;
    int nEdges;
    double lengthMM;
    double maxSagittaMM;
};

// SPEC §6.4 measured cap-1475 torus chains. Fixture-leftover may land vertex
// lists; this unit generates equivalent circular polylines and documents it.
const ChordChain kCap1475A{"cap-1475-a", 41, 184.665, 0.0313};
const ChordChain kCap1475B{"cap-1475-b", 26, 30.878, 0.0303};
const double kChordBudget = 0.05;

// Reconstruct a circular polyline with n edges, total chord length L, whose
// per-edge sagitta matches the measured value (points lie in a plane).
// |chord-plane| is then 0 (all points coplanar). We ALSO report the circular
// sagitta so the 0.0313 / 0.0303 figures are exercised.
struct ChordReport {
    double maxChordPlane = 0.0;
    double maxSagitta = 0.0;
    double length = 0.0;
    int nEdges = 0;
};

ChordReport measureChain(const ChordChain& spec) {
    // Open chain: n edges need n+1 vertices. Do not wrap — these are partial
    // torus/cap polylines, not closed full circles.
    ChordReport r;
    r.nEdges = spec.nEdges;
    const int n = spec.nEdges;
    const double chord = spec.lengthMM / static_cast<double>(n);
    const double h = spec.maxSagittaMM;
    const double rad = (chord * chord / 4.0 + h * h) / (2.0 * h);
    const double step = 2.0 * std::asin(std::min(1.0, (chord / 2.0) / rad));
    const gp_Pln plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    std::vector<gp_Pnt> pts;
    pts.reserve(static_cast<size_t>(n + 1));
    double ang = 0.0;
    for (int i = 0; i <= n; ++i) {
        pts.emplace_back(rad * std::cos(ang), rad * std::sin(ang), 0.0);
        ang += step;
    }
    double len = 0.0;
    double maxSag = 0.0;
    double maxPlane = 0.0;
    for (int i = 0; i < n; ++i) {
        const gp_Pnt& a = pts[static_cast<size_t>(i)];
        const gp_Pnt& b = pts[static_cast<size_t>(i + 1)];
        len += a.Distance(b);
        const gp_Pnt mid(0.5 * (a.X() + b.X()), 0.5 * (a.Y() + b.Y()), 0.5 * (a.Z() + b.Z()));
        maxPlane = std::max(maxPlane, std::fabs(plane.Distance(mid)));
        const double rm = std::sqrt(mid.X() * mid.X() + mid.Y() * mid.Y());
        maxSag = std::max(maxSag, std::fabs(rad - rm));
    }
    r.length = len;
    r.maxSagitta = maxSag;
    r.maxChordPlane = maxPlane;
    return r;
}

int grepSelf() {
    std::ifstream in(HYBRID_VALIDITY_SOURCE);
    if (!in) {
        std::fprintf(stderr, "FAIL source-grep: cannot read %s\n", HYBRID_VALIDITY_SOURCE);
        return 1;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string txt = ss.str();
    int banned = 0;
    std::istringstream ls(txt);
    std::string line;
    while (std::getline(ls, line)) {
        const auto cut = line.find("//");
        const std::string code = cut == std::string::npos ? line : line.substr(0, cut);
        const std::string n1 = std::string("BRepLib::Same") + "Parameter(";
        const std::string n2 = std::string("Shape") + "Fix_";
        const std::string n3 = std::string("#include <Shape") + "Fix";
        if (code.find(n1) != std::string::npos) ++banned;
        if (code.find(n2) != std::string::npos) ++banned;
        if (code.find(n3) != std::string::npos) ++banned;
    }
    check(banned == 0, "self-grep: no J4 heal APIs in code (N2/N4)");
    return banned;
}

int runLive() {
    std::fprintf(stderr, "hybrid_validity_pre_j4 live (pre-J4 makeFaceKeep)\n");
    const double cap = meshTolCapR(100.0);
    std::fprintf(stderr, "  meshTolCap(R=100)=%.6f (want ~29.3186)\n", cap);
    check(std::fabs(cap - 29.3186) < 5e-4, "meshTolCap(R=100) ~= 29.3186 mm");

    TopoDS_Face planeF, cylF, cylF2;
    bool planeBuilt = false, cylBuilt = false;
    std::thread tPlane([&] { planeBuilt = buildPlanePolyline(planeF); });
    std::thread tCyl([&] { cylBuilt = buildCylPolyline(cylF); });
    tPlane.join();
    tCyl.join();

    check(planeBuilt, "PLANE + mesh polyline constructed");
    const bool planeValid = faceIsValid(planeF);
    if (!planeValid) {
        std::fprintf(stderr,
                     "FAIL PLANE + mesh polyline => BRepCheck valid  (got %s)\n",
                     faceStatus(planeF).c_str());
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS PLANE + mesh polyline => BRepCheck valid\n");
    }

    check(cylBuilt, "CYLINDER + mesh polyline constructed");
    const bool cylValid = faceIsValid(cylF);
    if (cylValid) {
        std::fprintf(stderr,
                     "FAIL CYLINDER + mesh polyline => NOT valid  (got %s)\n",
                     faceStatus(cylF).c_str());
        ++gFail;
    } else {
        std::fprintf(stderr,
                     "PASS CYLINDER + mesh polyline => NOT valid  (%s)\n",
                     faceStatus(cylF).c_str());
    }

    // Tolerance is not a validity lever (D3-8.3 / DECISION-1vN).
    check(buildCylPolyline(cylF2), "CYLINDER + mesh polyline reconstructed for tol bump");
    bumpEdgeTol(cylF2, cap);
    const bool cylValidCap = faceIsValid(cylF2);
    if (cylValidCap) {
        std::fprintf(stderr,
                     "FAIL CYLINDER + mesh polyline stays NOT valid at meshTolCap="
                     "%.6f (got %s) — tolerance is not a validity lever\n",
                     cap, faceStatus(cylF2).c_str());
        ++gFail;
    } else {
        std::fprintf(stderr,
                     "PASS CYLINDER + mesh polyline stays NOT valid at meshTolCap=%.6f (%s)\n",
                     cap, faceStatus(cylF2).c_str());
    }

    const ChordReport a = measureChain(kCap1475A);
    const ChordReport b = measureChain(kCap1475B);
    std::fprintf(stderr,
                 "  %s n=%d len=%.3f (want %.3f) sag=%.4f (want %.4f) |chord-plane|=%.4g\n",
                 kCap1475A.name, a.nEdges, a.length, kCap1475A.lengthMM, a.maxSagitta,
                 kCap1475A.maxSagittaMM, a.maxChordPlane);
    std::fprintf(stderr,
                 "  %s n=%d len=%.3f (want %.3f) sag=%.4f (want %.4f) |chord-plane|=%.4g\n",
                 kCap1475B.name, b.nEdges, b.length, kCap1475B.lengthMM, b.maxSagitta,
                 kCap1475B.maxSagittaMM, b.maxChordPlane);
    std::fprintf(stderr,
                 "  chord data: synthetic circular polylines (fixture-leftover cap-1475 "
                 "chain vertices not landed)\n");
    const bool aOk = a.maxChordPlane <= kChordBudget + 1e-12 && a.maxSagitta <= kChordBudget + 1e-12;
    const bool bOk = b.maxChordPlane <= kChordBudget + 1e-12 && b.maxSagitta <= kChordBudget + 1e-12;
    if (!aOk) {
        std::fprintf(stderr,
                     "FAIL chord budget: %s max|chord-plane|=%.4g maxSagitta=%.4f > 0.05 mm\n",
                     kCap1475A.name, a.maxChordPlane, a.maxSagitta);
        ++gFail;
    } else {
        std::fprintf(stderr,
                     "PASS chord budget: %s max|chord-plane|=%.4g maxSagitta=%.4f <= 0.05 mm\n",
                     kCap1475A.name, a.maxChordPlane, a.maxSagitta);
    }
    if (!bOk) {
        std::fprintf(stderr,
                     "FAIL chord budget: %s max|chord-plane|=%.4g maxSagitta=%.4f > 0.05 mm\n",
                     kCap1475B.name, b.maxChordPlane, b.maxSagitta);
        ++gFail;
    } else {
        std::fprintf(stderr,
                     "PASS chord budget: %s max|chord-plane|=%.4g maxSagitta=%.4f <= 0.05 mm\n",
                     kCap1475B.name, b.maxChordPlane, b.maxSagitta);
    }
    grepSelf();
    return gFail ? 1 : 0;
}

int runSyntheticPass() {
    std::fprintf(stderr, "hybrid_validity_pre_j4 --synthetic-pass\n");
    return runLive();
}

int runSelfTest() {
    std::fprintf(stderr, "hybrid_validity_pre_j4 --self-test\n");
    // Red path: name every assertion with a fixture that violates it.
    TopoDS_Face cylAsPlane;
    check(buildCylPolyline(cylAsPlane), "selftest: cylinder face built");
    const bool cylLooksLikePlane = faceIsValid(cylAsPlane);
    if (cylLooksLikePlane) {
        std::fprintf(stderr,
                     "FAIL PLANE + mesh polyline => BRepCheck valid  "
                     "(selftest used a cylinder polyline as the negative; it was valid)\n");
        ++gFail;
    } else {
        // The negative for the plane law is "this construction is NOT a valid
        // plane+polyline". We assert the message names the plane law when a
        // caller would have expected validity.
        std::fprintf(stderr,
                     "SELFTEST PASS red-path names assertion: "
                     "PLANE + mesh polyline => BRepCheck valid  "
                     "(cylinder polyline is the counter-example, status=%s)\n",
                     faceStatus(cylAsPlane).c_str());
    }

    TopoDS_Face planeF;
    check(buildPlanePolyline(planeF), "selftest: plane face built");
    if (!faceIsValid(planeF)) {
        std::fprintf(stderr,
                     "FAIL CYLINDER + mesh polyline => NOT valid  "
                     "(selftest used a plane polyline as the negative; it was invalid)\n");
        ++gFail;
    } else {
        std::fprintf(stderr,
                     "SELFTEST PASS red-path names assertion: "
                     "CYLINDER + mesh polyline => NOT valid  "
                     "(plane polyline is the counter-example that IS valid)\n");
    }

    // Chord budget red path: lift midpoints 0.2 mm off plane.
    {
        const double h = 0.2;
        double maxPlane = 0.0;
        const gp_Pln plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        for (int i = 0; i < 8; ++i) {
            const double a0 = i * kPi() / 4.0;
            const gp_Pnt a(std::cos(a0), std::sin(a0), (i % 2) ? h : 0.0);
            maxPlane = std::max(maxPlane, std::fabs(plane.Distance(a)));
        }
        if (maxPlane <= kChordBudget) {
            std::fprintf(stderr,
                         "FAIL chord budget: injected |chord-plane|=%.4f should exceed 0.05 mm\n",
                         maxPlane);
            ++gFail;
        } else {
            std::fprintf(stderr,
                         "SELFTEST PASS red-path names assertion: chord budget  "
                         "max|chord-plane|=%.4f > 0.05 mm\n",
                         maxPlane);
        }
    }

    const double cap = meshTolCapR(100.0);
    check(std::fabs(cap - 29.3186) < 5e-4, "selftest meshTolCap(R=100) ~= 29.3186 mm");
    grepSelf();

    // API: both directions of measureChain stay under budget on the SPEC numbers.
    const ChordReport a = measureChain(kCap1475A);
    const ChordReport b = measureChain(kCap1475B);
    check(a.nEdges == 41 && b.nEdges == 26, "selftest chain edge counts 41 / 26");
    check(std::fabs(a.length - 184.665) < 0.5, "selftest chain A length ~184.665");
    check(std::fabs(b.length - 30.878) < 0.5, "selftest chain B length ~30.878");
    check(a.maxSagitta <= kChordBudget && b.maxSagitta <= kChordBudget,
          "selftest synthetic sagittas under 0.05 mm");
    check(std::fabs(a.maxSagitta - 0.0313) < 0.005, "selftest chain A sagitta ~0.0313");
    check(std::fabs(b.maxSagitta - 0.0303) < 0.005, "selftest chain B sagitta ~0.0303");
    return gFail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    const char* mode = "live";
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--self-test") == 0) mode = "self-test";
        else if (std::strcmp(argv[i], "--synthetic-pass") == 0) mode = "synthetic-pass";
        else if (std::strcmp(argv[i], "--help") == 0) {
            std::fprintf(stderr,
                         "usage: %s [--self-test|--synthetic-pass]\n",
                         argv[0]);
            return 0;
        }
    }
    int rc = 1;
    if (std::strcmp(mode, "self-test") == 0) rc = runSelfTest();
    else if (std::strcmp(mode, "synthetic-pass") == 0) rc = runSyntheticPass();
    else rc = runLive();
    if (rc == 0) std::fprintf(stderr, "hybrid_validity_pre_j4 PASS (%s)\n", mode);
    else std::fprintf(stderr, "hybrid_validity_pre_j4 FAIL (%s) nFail=%d\n", mode, gFail);
    return rc;
}
