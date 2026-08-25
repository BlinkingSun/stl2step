// stl2step_census — independent STEP reader / witness.
//
// Prints one JSON object describing the B-Rep actually stored in a STEP file:
// per-face surface types, per-cylinder radius/axis/U-span, edge curves and
// pcurves, BRepCheck validity, closedness, volume, continuity, tolerances.
//
// This translation unit includes NO stl2step headers and calls no engine
// function. It reads with STEPControl_Reader. A witness that shares the
// defendant's code is worthless.
//
// SPDX-License-Identifier: MIT

#include "step_census.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_Shape.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

namespace step_census {
namespace {

double kPi() { return std::acos(-1.0); }
double kTwoPi() { return 2.0 * kPi(); }

void silenceOcct() {
    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);
}

std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') {
            o += '\\';
            o += static_cast<char>(c);
        } else if (c >= 0x20) {
            o += static_cast<char>(c);
        }
    }
    return o;
}

std::string jsonNumber(double v) {
    if (!std::isfinite(v)) return "null";
    if (v == 0.0) v = 0.0;  // collapse -0.0
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

int countType(const TopoDS_Shape& s, TopAbs_ShapeEnum t) {
    TopTools_IndexedMapOfShape m;
    TopExp::MapShapes(s, t, m);
    return m.Extent();
}

double axisDist(const Axis& a, const Axis& b) {
    // Distance between two 3d lines (location + unit direction).
    const double dx = a.location.x - b.location.x;
    const double dy = a.location.y - b.location.y;
    const double dz = a.location.z - b.location.z;
    const double cx = a.direction.y * b.direction.z - a.direction.z * b.direction.y;
    const double cy = a.direction.z * b.direction.x - a.direction.x * b.direction.z;
    const double cz = a.direction.x * b.direction.y - a.direction.y * b.direction.x;
    const double cn = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (cn < 1e-12) {
        // Parallel: perpendicular distance from b.location to line a.
        const double fx = b.location.x - a.location.x;
        const double fy = b.location.y - a.location.y;
        const double fz = b.location.z - a.location.z;
        const double px = fy * a.direction.z - fz * a.direction.y;
        const double py = fz * a.direction.x - fx * a.direction.z;
        const double pz = fx * a.direction.y - fy * a.direction.x;
        return std::sqrt(px * px + py * py + pz * pz);
    }
    return std::fabs(dx * cx + dy * cy + dz * cz) / cn;
}

bool dirsParallel(const Vec3& a, const Vec3& b) {
    const double dot = a.x * b.x + a.y * b.y + a.z * b.z;
    return std::fabs(std::fabs(dot) - 1.0) < 1e-8;
}

Axis canonicalAxis(const Axis& in) {
    Axis a = in;
    const double n = std::sqrt(a.direction.x * a.direction.x +
                               a.direction.y * a.direction.y +
                               a.direction.z * a.direction.z);
    if (n > 0) {
        a.direction.x /= n;
        a.direction.y /= n;
        a.direction.z /= n;
    }
    // Flip so the largest-magnitude component is non-negative.
    const double ax = std::fabs(a.direction.x), ay = std::fabs(a.direction.y),
                 az = std::fabs(a.direction.z);
    bool flip = false;
    if (az >= ax && az >= ay) flip = a.direction.z < 0;
    else if (ay >= ax)        flip = a.direction.y < 0;
    else                      flip = a.direction.x < 0;
    if (flip) {
        a.direction.x = -a.direction.x;
        a.direction.y = -a.direction.y;
        a.direction.z = -a.direction.z;
    }
    // Location: point on the axis nearest the origin (stable across faces).
    const double t = a.location.x * a.direction.x + a.location.y * a.direction.y +
                     a.location.z * a.direction.z;
    a.location.x -= t * a.direction.x;
    a.location.y -= t * a.direction.y;
    a.location.z -= t * a.direction.z;
    return a;
}

int cmpD(double a, double b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

int cmpVec(const Vec3& a, const Vec3& b) {
    if (int c = cmpD(a.x, b.x)) return c;
    if (int c = cmpD(a.y, b.y)) return c;
    return cmpD(a.z, b.z);
}

int cmpAxis(const Axis& a, const Axis& b) {
    if (int c = cmpVec(a.location, b.location)) return c;
    return cmpVec(a.direction, b.direction);
}

void bumpContinuity(Result& r, GeomAbs_Shape s) {
    switch (s) {
        case GeomAbs_C0: r.contC0++; break;
        case GeomAbs_G1: r.contG1++; break;
        case GeomAbs_C1: r.contC1++; break;
        case GeomAbs_G2: r.contG2++; break;
        case GeomAbs_C2: r.contC2++; break;
        case GeomAbs_C3: r.contC3++; break;
        case GeomAbs_CN: r.contCN++; break;
        default:         r.contC0++; break;
    }
}

void classifySurface(Result& r, GeomAbs_SurfaceType t) {
    switch (t) {
        case GeomAbs_Plane:           r.plane++; break;
        case GeomAbs_Cylinder:        r.cylinder++; break;
        case GeomAbs_Cone:            r.cone++; break;
        case GeomAbs_Sphere:          r.sphere++; break;
        case GeomAbs_Torus:           r.torus++; break;
        case GeomAbs_BSplineSurface:  r.bspline++; break;
        default:                      r.otherSurf++; break;
    }
}

void classifyCurve(Result& r, GeomAbs_CurveType t) {
    switch (t) {
        case GeomAbs_Line:          r.line++; break;
        case GeomAbs_Circle:        r.circle++; break;
        case GeomAbs_Ellipse:       r.ellipse++; break;
        case GeomAbs_BSplineCurve:  r.bsplineCurve++; break;
        default:                    r.otherCurve++; break;
    }
}

void addCylinderFace(Result& r, const TopoDS_Face& face, const BRepAdaptor_Surface& ads) {
    const gp_Cylinder gc = ads.Cylinder();
    const gp_Ax1 ax = gc.Axis();
    const gp_Pnt loc = ax.Location();
    const gp_Dir dir = ax.Direction();

    Standard_Real umin = 0, umax = 0, vmin = 0, vmax = 0;
    BRepTools::UVBounds(face, umin, umax, vmin, vmax);
    // Face restriction on the adaptor is the fallback if UVBounds is degenerate.
    if (!(std::isfinite(umin) && std::isfinite(umax))) {
        umin = ads.FirstUParameter();
        umax = ads.LastUParameter();
        vmin = ads.FirstVParameter();
        vmax = ads.LastVParameter();
    }

    CylinderFace c;
    c.radius = gc.Radius();
    c.axis.location  = {loc.X(), loc.Y(), loc.Z()};
    c.axis.direction = {dir.X(), dir.Y(), dir.Z()};
    c.uMin = umin;
    c.uMax = umax;
    c.vMin = vmin;
    c.vMax = vmax;
    c.uSpan = std::fabs(umax - umin);
    c.periodicU = ads.IsUPeriodic() == Standard_True;
    c.closedU = ads.IsUClosed() == Standard_True;
    // A full-period face is the Seamed360 signature. Underlying Geom_Cylinder
    // is always U-periodic, so periodicU alone does NOT separate the two
    // constructions — uSpan / fullU does.
    c.fullU = c.uSpan + 1e-4 >= kTwoPi();
    r.cylinders.push_back(c);
}

void groupCylinders(Result& r) {
    const double radTol = Precision::Confusion();
    const double linTol = 1e-6;
    std::vector<int> used(r.cylinders.size(), 0);
    for (size_t i = 0; i < r.cylinders.size(); i++) {
        if (used[i]) continue;
        const Axis ai = canonicalAxis(r.cylinders[i].axis);
        CylinderGroup g;
        g.radius = r.cylinders[i].radius;
        g.axis = ai;
        std::vector<size_t> members;
        for (size_t j = i; j < r.cylinders.size(); j++) {
            if (used[j]) continue;
            if (std::fabs(r.cylinders[j].radius - g.radius) > radTol) continue;
            const Axis aj = canonicalAxis(r.cylinders[j].axis);
            if (!dirsParallel(ai.direction, aj.direction)) continue;
            if (axisDist(ai, aj) > linTol) continue;
            used[j] = 1;
            members.push_back(j);
            g.uSpans.push_back(r.cylinders[j].uSpan);
        }
        g.nFaces = static_cast<int>(members.size());
        std::sort(g.uSpans.begin(), g.uSpans.end());

        int nFull = 0, nHalf = 0;
        for (double s : g.uSpans) {
            if (s + 1e-4 >= kTwoPi()) nFull++;
            else if (std::fabs(s - kPi()) < 0.05) nHalf++;
        }
        if (g.nFaces == 1 && nFull == 1)      g.builtAs = "seamed360";
        else if (g.nFaces == 2 && nHalf == 2 && nFull == 0) g.builtAs = "twoHalves";
        else if (g.nFaces == 1 && nFull == 0) g.builtAs = "partial";
        else                                  g.builtAs = "other";
        r.cylinderGroups.push_back(std::move(g));
    }
    std::sort(r.cylinderGroups.begin(), r.cylinderGroups.end(),
              [](const CylinderGroup& a, const CylinderGroup& b) {
                  if (int c = cmpD(a.radius, b.radius)) return c < 0;
                  if (int c = cmpAxis(a.axis, b.axis)) return c < 0;
                  return a.nFaces < b.nFaces;
              });
}

void censusInto(Result& r, const TopoDS_Shape& shape) {
    r.solids   = countType(shape, TopAbs_SOLID);
    r.shells   = countType(shape, TopAbs_SHELL);
    r.faces    = countType(shape, TopAbs_FACE);
    r.edges    = countType(shape, TopAbs_EDGE);
    r.vertices = countType(shape, TopAbs_VERTEX);

    {
        TopTools_IndexedMapOfShape sm;
        TopExp::MapShapes(shape, TopAbs_SHELL, sm);
        r.shellList.reserve(static_cast<size_t>(sm.Extent()));
        r.openShells = 0;
        for (int i = 1; i <= sm.Extent(); i++) {
            ShellInfo s;
            s.closed = BRep_Tool::IsClosed(sm(i)) == Standard_True;
            if (!s.closed) r.openShells++;
            r.shellList.push_back(s);
        }
        r.closed = (r.shells > 0 && r.openShells == 0);
    }

    {
        BRepCheck_Analyzer ana(shape, Standard_True, Standard_True);
        r.valid = ana.IsValid() == Standard_True;
    }

    {
        GProp_GProps props;
        BRepGProp::VolumeProperties(shape, props);
        r.volume = props.Mass();
    }

    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Current());
        r.maxFaceTol = std::max(r.maxFaceTol, static_cast<double>(BRep_Tool::Tolerance(face)));
        BRepAdaptor_Surface ads(face, Standard_True);
        const GeomAbs_SurfaceType ty = ads.GetType();
        classifySurface(r, ty);
        if (ty == GeomAbs_Cylinder) addCylinderFace(r, face, ads);
    }
    std::sort(r.cylinders.begin(), r.cylinders.end(),
              [](const CylinderFace& a, const CylinderFace& b) {
                  if (int c = cmpD(a.radius, b.radius)) return c < 0;
                  if (int c = cmpAxis(canonicalAxis(a.axis), canonicalAxis(b.axis))) return c < 0;
                  if (int c = cmpD(a.uMin, b.uMin)) return c < 0;
                  if (int c = cmpD(a.uMax, b.uMax)) return c < 0;
                  return cmpD(a.vMin, b.vMin) < 0;
              });
    groupCylinders(r);

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaces;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeFaces);

    TopTools_IndexedMapOfShape em;
    TopExp::MapShapes(shape, TopAbs_EDGE, em);
    for (int i = 1; i <= em.Extent(); i++) {
        const TopoDS_Edge edge = TopoDS::Edge(em(i));
        r.maxEdgeTol = std::max(r.maxEdgeTol, static_cast<double>(BRep_Tool::Tolerance(edge)));

        if (BRep_Tool::Degenerated(edge)) {
            r.otherCurve++;
        } else {
            try {
                BRepAdaptor_Curve ac(edge);
                classifyCurve(r, ac.GetType());
            } catch (const Standard_Failure&) {
                r.otherCurve++;
            }
        }

        const TopTools_ListOfShape* lf = nullptr;
        if (edgeFaces.Contains(edge)) lf = &edgeFaces.FindFromKey(edge);
        const int nAdj = lf ? lf->Extent() : 0;

        if (lf) {
            for (TopTools_ListOfShape::Iterator it(*lf); it.More(); it.Next()) {
                const TopoDS_Face face = TopoDS::Face(it.Value());
                r.pcurveSlots++;
                Standard_Real first = 0, last = 0;
                Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(edge, face, first, last);
                if (pc.IsNull()) r.pcurveMissing++;
                else             r.pcurvePresent++;
            }
        }

        if (nAdj < 2) {
            r.contBoundary++;
        } else if (nAdj > 2) {
            r.contNonManifold++;
        } else {
            const TopoDS_Face f1 = TopoDS::Face(lf->First());
            const TopoDS_Face f2 = TopoDS::Face(lf->Last());
            if (BRep_Tool::HasContinuity(edge, f1, f2)) {
                bumpContinuity(r, BRep_Tool::Continuity(edge, f1, f2));
            } else {
                r.contC0++;
            }
        }
    }

    TopTools_IndexedMapOfShape vm;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vm);
    for (int i = 1; i <= vm.Extent(); i++) {
        r.maxVertexTol = std::max(
            r.maxVertexTol, static_cast<double>(BRep_Tool::Tolerance(TopoDS::Vertex(vm(i)))));
    }
}

Result fail(const std::string& input, const std::string& err) {
    Result r;
    r.ok = false;
    r.input = input;
    r.error = err;
    return r;
}

}  // namespace

Result censusFile(const std::string& path) {
    silenceOcct();
    try {
        STEPControl_Reader reader;
        const IFSelect_ReturnStatus st = reader.ReadFile(path.c_str());
        if (st != IFSelect_RetDone) {
            return fail(path, "not a readable STEP file");
        }
        if (reader.NbRootsForTransfer() < 1) {
            return fail(path, "STEP file has no transferable roots");
        }
        if (reader.TransferRoots() < 1) {
            return fail(path, "STEP transfer produced no shapes");
        }
        const TopoDS_Shape shape = reader.OneShape();
        if (shape.IsNull()) {
            return fail(path, "STEP transfer produced a null shape");
        }
        Result r;
        r.ok = true;
        r.input = path;
        censusInto(r, shape);
        return r;
    } catch (const Standard_Failure& e) {
        const char* m = e.GetMessageString();
        return fail(path, m && m[0] ? std::string("OCCT: ") + m : "OCCT exception");
    } catch (const std::exception& e) {
        return fail(path, e.what());
    } catch (...) {
        return fail(path, "unknown exception");
    }
}

std::string toJson(const Result& r) {
    auto vec3 = [](const Vec3& v) {
        return std::string("[") + jsonNumber(v.x) + "," + jsonNumber(v.y) + "," +
               jsonNumber(v.z) + "]";
    };
    auto axis = [&](const Axis& a) {
        return std::string("{\"location\":") + vec3(a.location) + ",\"direction\":" +
               vec3(a.direction) + "}";
    };
    auto b2s = [](bool b) { return b ? "true" : "false"; };

    if (!r.ok) {
        std::string s = "{\"ok\":false,\"error\":\"";
        s += jsonEscape(r.error);
        s += "\"";
        if (!r.input.empty()) {
            s += ",\"input\":\"";
            s += jsonEscape(r.input);
            s += "\"";
        }
        s += "}";
        return s;
    }

    std::string s;
    s.reserve(2048);
    s += "{\"ok\":true,\"input\":\"";
    s += jsonEscape(r.input);
    s += "\",\"valid\":";
    s += b2s(r.valid);
    s += ",\"solids\":";
    s += std::to_string(r.solids);
    s += ",\"shells\":";
    s += std::to_string(r.shells);
    s += ",\"faces\":";
    s += std::to_string(r.faces);
    s += ",\"edges\":";
    s += std::to_string(r.edges);
    s += ",\"vertices\":";
    s += std::to_string(r.vertices);
    s += ",\"openShells\":";
    s += std::to_string(r.openShells);
    s += ",\"closed\":";
    s += b2s(r.closed);
    s += ",\"volume\":";
    s += jsonNumber(r.volume);

    s += ",\"shellsClosed\":[";
    for (size_t i = 0; i < r.shellList.size(); i++) {
        if (i) s += ",";
        s += b2s(r.shellList[i].closed);
    }
    s += "]";

    s += ",\"surfaces\":{";
    s += "\"plane\":";     s += std::to_string(r.plane);
    s += ",\"cylinder\":"; s += std::to_string(r.cylinder);
    s += ",\"cone\":";     s += std::to_string(r.cone);
    s += ",\"sphere\":";   s += std::to_string(r.sphere);
    s += ",\"torus\":";    s += std::to_string(r.torus);
    s += ",\"bspline\":";  s += std::to_string(r.bspline);
    s += ",\"other\":";    s += std::to_string(r.otherSurf);
    s += "}";

    s += ",\"cylinders\":[";
    for (size_t i = 0; i < r.cylinders.size(); i++) {
        const CylinderFace& c = r.cylinders[i];
        if (i) s += ",";
        s += "{\"radius\":";
        s += jsonNumber(c.radius);
        s += ",\"axis\":";
        s += axis(c.axis);
        s += ",\"uMin\":";  s += jsonNumber(c.uMin);
        s += ",\"uMax\":";  s += jsonNumber(c.uMax);
        s += ",\"vMin\":";  s += jsonNumber(c.vMin);
        s += ",\"vMax\":";  s += jsonNumber(c.vMax);
        s += ",\"uSpan\":"; s += jsonNumber(c.uSpan);
        s += ",\"periodicU\":"; s += b2s(c.periodicU);
        s += ",\"closedU\":";   s += b2s(c.closedU);
        s += ",\"fullU\":";     s += b2s(c.fullU);
        s += "}";
    }
    s += "]";

    s += ",\"cylinderGroups\":[";
    for (size_t i = 0; i < r.cylinderGroups.size(); i++) {
        const CylinderGroup& g = r.cylinderGroups[i];
        if (i) s += ",";
        s += "{\"radius\":";
        s += jsonNumber(g.radius);
        s += ",\"axis\":";
        s += axis(g.axis);
        s += ",\"nFaces\":";
        s += std::to_string(g.nFaces);
        s += ",\"uSpans\":[";
        for (size_t k = 0; k < g.uSpans.size(); k++) {
            if (k) s += ",";
            s += jsonNumber(g.uSpans[k]);
        }
        s += "],\"builtAs\":\"";
        s += jsonEscape(g.builtAs);
        s += "\"}";
    }
    s += "]";

    s += ",\"curves\":{";
    s += "\"line\":";     s += std::to_string(r.line);
    s += ",\"circle\":";  s += std::to_string(r.circle);
    s += ",\"ellipse\":"; s += std::to_string(r.ellipse);
    s += ",\"bspline\":"; s += std::to_string(r.bsplineCurve);
    s += ",\"other\":";   s += std::to_string(r.otherCurve);
    s += "}";

    s += ",\"pcurves\":{";
    s += "\"slots\":";    s += std::to_string(r.pcurveSlots);
    s += ",\"present\":"; s += std::to_string(r.pcurvePresent);
    s += ",\"missing\":"; s += std::to_string(r.pcurveMissing);
    s += "}";

    s += ",\"continuity\":{";
    s += "\"C0\":"; s += std::to_string(r.contC0);
    s += ",\"G1\":"; s += std::to_string(r.contG1);
    s += ",\"C1\":"; s += std::to_string(r.contC1);
    s += ",\"G2\":"; s += std::to_string(r.contG2);
    s += ",\"C2\":"; s += std::to_string(r.contC2);
    s += ",\"C3\":"; s += std::to_string(r.contC3);
    s += ",\"CN\":"; s += std::to_string(r.contCN);
    s += ",\"boundary\":";    s += std::to_string(r.contBoundary);
    s += ",\"nonManifold\":"; s += std::to_string(r.contNonManifold);
    s += "}";

    s += ",\"tolerances\":{";
    s += "\"maxFace\":";   s += jsonNumber(r.maxFaceTol);
    s += ",\"maxEdge\":";  s += jsonNumber(r.maxEdgeTol);
    s += ",\"maxVertex\":"; s += jsonNumber(r.maxVertexTol);
    s += "}}";
    return s;
}

}  // namespace step_census

#ifndef STL2STEP_CENSUS_NO_MAIN
int main(int argc, char** argv) {
    if (argc != 2 || std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        std::fprintf(stderr, "usage: stl2step_census <file.step>\n");
        return argc == 2 ? 0 : 1;
    }
    const step_census::Result r = step_census::censusFile(argv[1]);
    std::puts(step_census::toJson(r).c_str());
    return r.ok ? 0 : 1;
}
#endif
