// stl2step --mesh mode: STEP (B-Rep) -> binary STL + optional Format A edges.
//
// Sibling of convert(). Never called from Converter::run(); never writes
// Result / Result::toJson().
//
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include <Bnd_Box.hxx>
#include <BRep_Tool.hxx>
#include <BRepBndLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <Poly_Polygon3D.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Failure.hxx>
#include <STEPControl_Reader.hxx>
#include <StlAPI_Writer.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

namespace stl2step {
namespace {

namespace fs = std::filesystem;

std::string jsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if ((unsigned char)c >= 0x20) o += c;
    }
    return o;
}

void appendF32LE(std::string& buf, float v) {
    unsigned char ieee[4];
    std::memcpy(ieee, &v, 4);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    std::swap(ieee[0], ieee[3]);
    std::swap(ieee[1], ieee[2]);
#endif
    buf.append(reinterpret_cast<char*>(ieee), 4);
}

struct EdgeDump {
    int drawable = 0;
    int segments = 0;
    std::string bytes;
};

// Unique B-Rep edges via MapShapes. Skip Degenerated. Prefer
// PolygonOnTriangulation index 1 (cylinder seam once); else Polygon3D.
// Apply TopLoc_Location to every node. Consecutive node pairs = 24-byte segments.
EdgeDump extractEdges(const TopoDS_Shape& shape) {
    EdgeDump d;
    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(shape, TopAbs_EDGE, map);
    for (int i = 1; i <= map.Extent(); ++i) {
        const TopoDS_Edge edge = TopoDS::Edge(map(i));
        if (BRep_Tool::Degenerated(edge)) continue;

        std::vector<gp_Pnt> nodes;
        Handle(Poly_PolygonOnTriangulation) poly;
        Handle(Poly_Triangulation) tri;
        TopLoc_Location loc;
        BRep_Tool::PolygonOnTriangulation(edge, poly, tri, loc, 1);
        if (!poly.IsNull() && !tri.IsNull()) {
            const gp_Trsf tr = loc.Transformation();
            const int n = poly->NbNodes();
            nodes.reserve((size_t)n);
            for (int k = 1; k <= n; ++k) {
                gp_Pnt p = tri->Node(poly->Node(k));
                p.Transform(tr);
                nodes.push_back(p);
            }
        } else {
            Handle(Poly_Polygon3D) p3d = BRep_Tool::Polygon3D(edge, loc);
            if (p3d.IsNull()) continue;
            const gp_Trsf tr = loc.Transformation();
            const auto& arr = p3d->Nodes();
            nodes.reserve((size_t)arr.Length());
            for (int k = arr.Lower(); k <= arr.Upper(); ++k) {
                gp_Pnt p = arr(k);
                p.Transform(tr);
                nodes.push_back(p);
            }
        }
        if (nodes.size() < 2) continue;

        int emitted = 0;
        for (size_t k = 0; k + 1 < nodes.size(); ++k) {
            appendF32LE(d.bytes, (float)nodes[k].X());
            appendF32LE(d.bytes, (float)nodes[k].Y());
            appendF32LE(d.bytes, (float)nodes[k].Z());
            appendF32LE(d.bytes, (float)nodes[k + 1].X());
            appendF32LE(d.bytes, (float)nodes[k + 1].Y());
            appendF32LE(d.bytes, (float)nodes[k + 1].Z());
            emitted++;
        }
        if (emitted > 0) {
            d.drawable++;
            d.segments += emitted;
        }
    }
    return d;
}

void removeQuietly(const std::string& path) {
    if (path.empty()) return;
    std::error_code ec;
    fs::remove(path, ec);
}

}  // namespace

std::string MeshResult::toJson() const {
    std::string s = "{\"ok\":";
    s += ok ? "true" : "false";
    s += ",\"faces\":" + std::to_string(faces);
    s += ",\"edges\":" + std::to_string(edges);
    s += ",\"triangles\":" + std::to_string(triangles);
    char sec[64];
    std::snprintf(sec, sizeof sec, "%.2f", seconds);
    s += ",\"seconds\":";
    s += sec;
    s += ",\"input\":\"" + jsonEscape(input) + "\"";
    s += ",\"output\":\"" + jsonEscape(output) + "\"";
    if (!edgesFile.empty())
        s += ",\"edgesFile\":\"" + jsonEscape(edgesFile) + "\"";
    if (!ok)
        s += ",\"error\":\"" + jsonEscape(error) + "\"";
    s += "}";
    return s;
}

MeshResult meshFromStep(const MeshOptions& opt, const LogCallback& log) {
    MeshResult r;
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    auto fail = [&](const std::string& err) -> MeshResult {
        r.ok = false;
        r.exitCode = 1;
        r.error = err;
        r.seconds = elapsed();
        if (log) log(Severity::Error, err);
        return r;
    };
    auto note = [&](const std::string& msg) {
        if (log) log(Severity::Info, msg);
    };

    r.input = opt.input;
    r.output = opt.output;
    if (r.input.empty()) return fail("no input STEP path given");
    const bool outputExplicit = !opt.output.empty();
    if (!outputExplicit) {
        // Default beside the input: <stem>.mesh.stl — never <stem>.stl, which
        // is the user's source mesh after a convert in the same folder.
        fs::path p(r.input);
        p.replace_extension(".mesh.stl");
        r.output = p.string();
        std::error_code ec;
        if (fs::exists(r.output, ec))
            return fail("output exists: " + r.output + " — pass -o");
    }
    // --edges is explicit-only: never invent a default edges path.

    unsigned hw = opt.threads > 0 ? (unsigned)opt.threads
                                  : std::thread::hardware_concurrency();
    if (hw == 0) hw = 4;

    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);

    try {
        note(std::string("stl2step --mesh: ") + r.input + " -> " + r.output);

        STEPControl_Reader reader;
        if (reader.ReadFile(r.input.c_str()) != IFSelect_RetDone)
            return fail("could not read STEP file: " + r.input);
        if (reader.TransferRoots() <= 0)
            return fail("STEP file contains no transferable roots: " + r.input);
        TopoDS_Shape shape = reader.OneShape();
        if (shape.IsNull())
            return fail("STEP file contains no geometry: " + r.input);

        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(shape, TopAbs_FACE, faceMap);
        r.faces = faceMap.Extent();
        if (r.faces < 1)
            return fail("STEP file contains no faces: " + r.input);

        Bnd_Box box;
        BRepBndLib::Add(shape, box, Standard_False);
        if (box.IsVoid())
            return fail("STEP geometry has an empty bounding box");
        const double diag = box.CornerMin().Distance(box.CornerMax());
        const double linDefl = std::max(0.05, 0.001 * diag);

        IMeshTools_Parameters params;
        params.Deflection = linDefl;
        params.Angle = 0.5;
        params.Relative = Standard_False;
        params.InParallel = (hw > 1) ? Standard_True : Standard_False;
        params.AllowQualityDecrease = Standard_False;
        BRepMesh_IncrementalMesh mesher;
        mesher.SetShape(shape);
        mesher.ChangeParameters() = params;
        mesher.Perform();
        if (!mesher.IsDone())
            return fail("tessellation failed");

        StlAPI_Writer writer;
        writer.ASCIIMode() = Standard_False;
        if (!writer.Write(shape, r.output.c_str())) {
            removeQuietly(r.output);
            return fail("could not write binary STL: " + r.output);
        }

        std::error_code ec;
        const auto sz = fs::file_size(r.output, ec);
        if (ec || sz < 84 || ((sz - 84) % 50) != 0) {
            removeQuietly(r.output);
            return fail("STL writer did not produce a binary STL: " + r.output);
        }
        r.triangles = (int)((sz - 84) / 50);
        if (r.triangles < 1) {
            removeQuietly(r.output);
            return fail("tessellation produced no triangles");
        }

        EdgeDump dump = extractEdges(shape);
        r.edges = dump.drawable;

        if (!opt.edgesFile.empty()) {
            std::ofstream out(opt.edgesFile, std::ios::binary | std::ios::trunc);
            if (!out) {
                removeQuietly(r.output);
                return fail("could not write edges file: " + opt.edgesFile);
            }
            if (!dump.bytes.empty())
                out.write(dump.bytes.data(), (std::streamsize)dump.bytes.size());
            out.flush();
            if (!out) {
                removeQuietly(r.output);
                removeQuietly(opt.edgesFile);
                return fail("could not write edges file: " + opt.edgesFile);
            }
            r.edgesFile = opt.edgesFile;
        }

        r.seconds = elapsed();
        r.ok = true;
        r.exitCode = 0;
        note("  mesh      " + std::to_string(r.faces) + " faces, " +
             std::to_string(r.triangles) + " triangles, " +
             std::to_string(r.edges) + " drawable edges");
        return r;

    } catch (const Standard_Failure& f) {
        removeQuietly(r.output);
        if (!opt.edgesFile.empty()) removeQuietly(opt.edgesFile);
        return fail(std::string("OCCT failure: ") +
                    (f.GetMessageString() ? f.GetMessageString() : "?"));
    } catch (const std::exception& e) {
        removeQuietly(r.output);
        if (!opt.edgesFile.empty()) removeQuietly(opt.edgesFile);
        return fail(e.what());
    }
}

}  // namespace stl2step
