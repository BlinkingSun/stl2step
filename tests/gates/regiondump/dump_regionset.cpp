// stl2step_regiondump — weld/split via the mesh harness, run refit::segment(),
// dump each clean component's RegionSet as JSON (frozen regionset.schema.json).
//
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"

#ifndef STL2STEP_HARNESS_HAVE_REFIT
#error "STL2STEP_HARNESS_HAVE_REFIT is not defined; src/refit.hpp must be visible"
#endif

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using stl2step::harness::HarnessComponent;
using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::BoundaryChain;
using stl2step::refit::BuiltAs;
using stl2step::refit::Loop;
using stl2step::refit::LoopRole;
using stl2step::refit::MeshView;
using stl2step::refit::Origin;
using stl2step::refit::RefitStats;
using stl2step::refit::Region;
using stl2step::refit::RegionSet;
using stl2step::refit::Reject;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;

// ---------------------------------------------------------------------------
// JSON primitives. Fixed key order, %.17g floats, no wall-clock.
// ---------------------------------------------------------------------------

static void putIndent(FILE* f, int n) {
    for (int i = 0; i < n; i++) std::fputc(' ', f);
}

static void putFloat(FILE* f, double v) {
    if (!std::isfinite(v)) {
        std::fputs("0", f);
        return;
    }
    if (v == 0.0) v = 0.0;  // collapse -0.0
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    std::fputs(buf, f);
}

static void putIntArray(FILE* f, const std::vector<int>& v) {
    std::fputc('[', f);
    for (size_t i = 0; i < v.size(); i++) {
        if (i) std::fputs(", ", f);
        std::fprintf(f, "%d", v[i]);
    }
    std::fputc(']', f);
}

static void putU8Array(FILE* f, const std::vector<uint8_t>& v) {
    std::fputc('[', f);
    for (size_t i = 0; i < v.size(); i++) {
        if (i) std::fputs(", ", f);
        std::fprintf(f, "%d", (int)v[i]);
    }
    std::fputc(']', f);
}

static void putVec3(FILE* f, double x, double y, double z) {
    std::fputc('[', f);
    putFloat(f, x);
    std::fputs(", ", f);
    putFloat(f, y);
    std::fputs(", ", f);
    putFloat(f, z);
    std::fputc(']', f);
}

static std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') {
            o += '\\';
            o += (char)c;
        } else if (c >= 0x20) {
            o += (char)c;
        }
    }
    return o;
}

static const char* surfTypeName(SurfType t) {
    switch (t) {
        case SurfType::Plane:    return "plane";
        case SurfType::Cylinder: return "cylinder";
        case SurfType::Cone:     return "cone";
        case SurfType::Sphere:   return "sphere";
        case SurfType::Torus:    return "torus";
    }
    return "plane";
}

static const char* originName(Origin o) {
    switch (o) {
        case Origin::PlaneGrow:   return "planeGrow";
        case Origin::CylGrow:     return "cylGrow";
        case Origin::FilletStrip: return "filletStrip";
    }
    return "planeGrow";
}

static const char* builtAsName(BuiltAs b) {
    switch (b) {
        case BuiltAs::NotBuilt:         return "notBuilt";
        case BuiltAs::Single:           return "single";
        case BuiltAs::Seamed360:        return "seamed360";
        case BuiltAs::TwoHalves:        return "twoHalves";
        case BuiltAs::ExplodedToFacets: return "explodedToFacets";
    }
    return "notBuilt";
}

static const char* rejectName(Reject r) {
    switch (r) {
        case Reject::None:               return "none";
        case Reject::GaussPlanarity:     return "gaussPlanarity";
        case Reject::VertexResidual:     return "vertexResidual";
        case Reject::ChordConsistency:   return "chordConsistency";
        case Reject::RadiusSanity:       return "radiusSanity";
        case Reject::Span:               return "span";
        case Reject::FilletConsensus:    return "filletConsensus";
        case Reject::NeighborNotAnalytic:return "neighborNotAnalytic";
        case Reject::StripWidth:         return "stripWidth";
        case Reject::TorusNYI:           return "torusNYI";
        case Reject::ConeNYI:            return "coneNYI";
        case Reject::SphereNYI:          return "sphereNYI";
        case Reject::DirtyComponent:     return "dirtyComponent";
        case Reject::FaceBuildFailed:    return "faceBuildFailed";
        case Reject::ChainUnstable:      return "chainUnstable";
    }
    return "none";
}

static const char* loopRoleName(LoopRole r) {
    switch (r) {
        case LoopRole::Outer:   return "outer";
        case LoopRole::Inner:   return "inner";
        case LoopRole::CapLow:  return "capLow";
        case LoopRole::CapHigh: return "capHigh";
    }
    return "outer";
}

struct Obj {
    FILE* f;
    int ind;
    bool first = true;
    Obj(FILE* fp, int indent) : f(fp), ind(indent) { std::fputc('{', f); }
    void key(const char* k) {
        if (!first) std::fputc(',', f);
        first = false;
        std::fputc('\n', f);
        putIndent(f, ind + 2);
        std::fprintf(f, "\"%s\": ", k);
    }
    void end() {
        if (!first) {
            std::fputc('\n', f);
            putIndent(f, ind);
        }
        std::fputc('}', f);
    }
};

static void writeAx3(FILE* f, const gp_Ax3& ax, int ind) {
    Obj o(f, ind);
    o.key("loc");
    putVec3(f, ax.Location().X(), ax.Location().Y(), ax.Location().Z());
    o.key("dir");
    putVec3(f, ax.Direction().X(), ax.Direction().Y(), ax.Direction().Z());
    o.key("xdir");
    putVec3(f, ax.XDirection().X(), ax.XDirection().Y(), ax.XDirection().Z());
    o.end();
}

static void writeLoop(FILE* f, const Loop& lp, int ind) {
    Obj o(f, ind);
    o.key("chainIdx");
    putIntArray(f, lp.chainIdx);
    o.key("reversed");
    putU8Array(f, lp.reversed);
    o.key("role");
    std::fprintf(f, "\"%s\"", loopRoleName(lp.role));
    o.end();
}

static void writeLoopArray(FILE* f, const std::vector<Loop>& v, int ind) {
    if (v.empty()) {
        std::fputs("[]", f);
        return;
    }
    std::fputc('[', f);
    for (size_t i = 0; i < v.size(); i++) {
        if (i) std::fputc(',', f);
        std::fputc('\n', f);
        putIndent(f, ind + 2);
        writeLoop(f, v[i], ind + 2);
    }
    std::fputc('\n', f);
    putIndent(f, ind);
    std::fputc(']', f);
}

static void writeRegion(FILE* f, const Region& r, int ind) {
    Obj o(f, ind);
    o.key("id");
    std::fprintf(f, "%d", r.id);
    o.key("type");
    std::fprintf(f, "\"%s\"", surfTypeName(r.type));
    o.key("origin");
    std::fprintf(f, "\"%s\"", originName(r.origin));
    o.key("ax");
    writeAx3(f, r.ax, ind + 2);
    o.key("radius");
    putFloat(f, r.radius);
    o.key("uMin");
    putFloat(f, r.uMin);
    o.key("uMax");
    putFloat(f, r.uMax);
    o.key("vMin");
    putFloat(f, r.vMin);
    o.key("vMax");
    putFloat(f, r.vMax);
    o.key("closed360");
    std::fputs(r.closed360 ? "true" : "false", f);
    o.key("outwardNormal");
    std::fputs(r.outwardNormal ? "true" : "false", f);
    o.key("tris");
    putIntArray(f, r.tris);
    o.key("loops");
    writeLoopArray(f, r.loops, ind + 2);
    o.key("maxVertexDev");
    putFloat(f, r.maxVertexDev);
    o.key("rmsVertexDev");
    putFloat(f, r.rmsVertexDev);
    o.key("chordSagitta");
    putFloat(f, r.chordSagitta);
    o.key("nSides");
    std::fprintf(f, "%d", r.nSides);
    o.key("dVolPredicted");
    putFloat(f, r.dVolPredicted);
    o.key("maxVertexSnap");
    putFloat(f, r.maxVertexSnap);
    o.key("reject");
    std::fprintf(f, "\"%s\"", rejectName(r.reject));
    o.key("builtAs");
    std::fprintf(f, "\"%s\"", builtAsName(r.builtAs));
    o.key("filletNbrA");
    std::fprintf(f, "%d", r.filletNbrA);
    o.key("filletNbrB");
    std::fprintf(f, "%d", r.filletNbrB);
    o.end();
}

static void writeRegionArray(FILE* f, const std::vector<Region>& v, int ind) {
    if (v.empty()) {
        std::fputs("[]", f);
        return;
    }
    std::fputc('[', f);
    for (size_t i = 0; i < v.size(); i++) {
        if (i) std::fputc(',', f);
        std::fputc('\n', f);
        putIndent(f, ind + 2);
        writeRegion(f, v[i], ind + 2);
    }
    std::fputc('\n', f);
    putIndent(f, ind);
    std::fputc(']', f);
}

static void writeChain(FILE* f, const BoundaryChain& ch, int ind) {
    Obj o(f, ind);
    o.key("regA");
    std::fprintf(f, "%d", ch.regA);
    o.key("regB");
    std::fprintf(f, "%d", ch.regB);
    o.key("islandA");
    std::fprintf(f, "%d", ch.islandA);
    o.key("islandB");
    std::fprintf(f, "%d", ch.islandB);
    o.key("tangent");
    std::fputs(ch.tangent ? "true" : "false", f);
    o.key("closedLoop");
    std::fputs(ch.closedLoop ? "true" : "false", f);
    o.key("meshEdges");
    putIntArray(f, ch.meshEdges);
    o.key("meshVerts");
    putIntArray(f, ch.meshVerts);
    o.end();
}

static void writeChainArray(FILE* f, const std::vector<BoundaryChain>& v, int ind) {
    if (v.empty()) {
        std::fputs("[]", f);
        return;
    }
    std::fputc('[', f);
    for (size_t i = 0; i < v.size(); i++) {
        if (i) std::fputc(',', f);
        std::fputc('\n', f);
        putIndent(f, ind + 2);
        writeChain(f, v[i], ind + 2);
    }
    std::fputc('\n', f);
    putIndent(f, ind);
    std::fputc(']', f);
}

static void writeStats(FILE* f, const RefitStats& s, int ind) {
    Obj o(f, ind);
    o.key("planes");
    std::fprintf(f, "%d", s.planes);
    o.key("cylinders");
    std::fprintf(f, "%d", s.cylinders);
    o.key("fillets");
    std::fprintf(f, "%d", s.fillets);
    o.key("rejected");
    std::fprintf(f, "%d", s.rejected);
    o.key("facetIslands");
    std::fprintf(f, "%d", s.facetIslands);
    o.key("facetTriangles");
    std::fprintf(f, "%d", s.facetTriangles);
    o.key("distinctRadii");
    std::fprintf(f, "%d", s.distinctRadii);
    o.key("maxVertexDev");
    putFloat(f, s.maxVertexDev);
    o.key("maxEdgeTol");
    putFloat(f, s.maxEdgeTol);
    o.key("dVolPredicted");
    putFloat(f, s.dVolPredicted);
    o.end();
}

static void writeRegionSet(FILE* f, const RegionSet& rs, int ind) {
    Obj o(f, ind);
    o.key("compRoot");
    std::fprintf(f, "%d", rs.compRoot);
    o.key("regions");
    writeRegionArray(f, rs.regions, ind + 2);
    o.key("rejected");
    writeRegionArray(f, rs.rejected, ind + 2);
    o.key("chains");
    writeChainArray(f, rs.chains, ind + 2);
    o.key("triRegion");
    putIntArray(f, rs.triRegion);
    o.key("triIsland");
    putIntArray(f, rs.triIsland);
    o.key("nIslands");
    std::fprintf(f, "%d", rs.nIslands);
    o.key("stats");
    writeStats(f, rs.stats, ind + 2);
    o.end();
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

static void usage(const char* argv0) {
    std::fprintf(stderr,
                 "usage: %s <file.stl> [--scale S] [--weld-tol T] [--sew-tol T]\n"
                 "                    [--component N] [--smooth-tol MM] [--smooth-angle DEG]\n"
                 "                    [--no-fillets] [--out FILE]\n",
                 argv0);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    std::string path;
    std::string outPath;
    double scale = 1.0, weldTol = 0.0, sewTol = 0.0;
    double smoothTolMM = 0.0, smoothAngleDeg = 2.0;
    bool doFillets = true;
    int componentFilter = -1;  // -1 = all

    for (int i = 1; i < argc; i++) {
        auto needArg = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "stl2step_regiondump: %s needs a value\n", flag);
                usage(argv[0]);
                std::exit(1);
            }
            return argv[++i];
        };
        if (std::strcmp(argv[i], "--scale") == 0) {
            scale = std::atof(needArg("--scale"));
        } else if (std::strcmp(argv[i], "--weld-tol") == 0) {
            weldTol = std::atof(needArg("--weld-tol"));
        } else if (std::strcmp(argv[i], "--sew-tol") == 0) {
            sewTol = std::atof(needArg("--sew-tol"));
        } else if (std::strcmp(argv[i], "--component") == 0) {
            componentFilter = std::atoi(needArg("--component"));
        } else if (std::strcmp(argv[i], "--smooth-tol") == 0) {
            smoothTolMM = std::atof(needArg("--smooth-tol"));
        } else if (std::strcmp(argv[i], "--smooth-angle") == 0) {
            smoothAngleDeg = std::atof(needArg("--smooth-angle"));
        } else if (std::strcmp(argv[i], "--no-fillets") == 0) {
            doFillets = false;
        } else if (std::strcmp(argv[i], "--out") == 0) {
            outPath = needArg("--out");
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            return 1;
        } else if (path.empty()) {
            path = argv[i];
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (path.empty()) {
        usage(argv[0]);
        return 1;
    }

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(path, scale, weldTol, sewTol, mesh, err)) {
        std::fprintf(stderr, "stl2step_regiondump: %s\n", err.c_str());
        return 1;
    }

    if (componentFilter >= 0 && componentFilter >= (int)mesh.comps.size()) {
        std::fprintf(stderr,
                     "stl2step_regiondump: --component %d out of range (0..%d)\n",
                     componentFilter, (int)mesh.comps.size() - 1);
        return 1;
    }

    SegmentParams params;  // remaining fields keep refit.hpp defaults
    params.epsPlane = smoothTolMM;
    params.thetaPlaneDeg = smoothAngleDeg;
    params.doFillets = doFillets;

    struct Dumped {
        size_t index;
        const HarnessComponent* c;
        RegionSet rs;
        bool segmentOk;
    };
    std::vector<Dumped> dumped;
    int dirtySkipped = 0;

    for (size_t i = 0; i < mesh.comps.size(); i++) {
        if (componentFilter >= 0 && (int)i != componentFilter) continue;
        const HarnessComponent& c = mesh.comps[i];
        if (!c.clean()) {
            dirtySkipped++;  // counted, never warned (I6)
            continue;
        }
        MeshView mv{};
        toMeshView(mesh, c, mv);
        Dumped d;
        d.index = i;
        d.c = &c;
        d.rs.compRoot = c.root;
        d.segmentOk = stl2step::refit::segment(mv, params, d.rs, nullptr);
        d.rs.compRoot = c.root;  // harness root is the component key
        dumped.push_back(std::move(d));
    }

    FILE* out = stdout;
    if (!outPath.empty()) {
        out = std::fopen(outPath.c_str(), "wb");
        if (!out) {
            std::fprintf(stderr, "stl2step_regiondump: cannot write %s\n",
                         outPath.c_str());
            return 1;
        }
    }

    Obj top(out, 0);
#ifdef STL2STEP_REGIONDUMP_STUB
    // Marker is wrapper-only (not inside RegionSet) so additionalProperties:
    // false still holds. Easy to strip; never emitted on the real path.
    top.key("stub");
    std::fputs("true", out);
#endif
    top.key("file");
    std::fprintf(out, "\"%s\"", jsonEscape(path).c_str());
    top.key("triangles");
    std::fprintf(out, "%zu", mesh.tris.size());
    top.key("vertices");
    std::fprintf(out, "%zu", mesh.pts.size());
    top.key("components");
    std::fprintf(out, "%zu", mesh.comps.size());
    top.key("dirtySkipped");
    std::fprintf(out, "%d", dirtySkipped);
    top.key("dumped");
    std::fprintf(out, "%zu", dumped.size());
    top.key("comps");
    if (dumped.empty()) {
        std::fputs("[]", out);
    } else {
        std::fputc('[', out);
        for (size_t i = 0; i < dumped.size(); i++) {
            if (i) std::fputc(',', out);
            std::fputc('\n', out);
            putIndent(out, 4);
            const Dumped& d = dumped[i];
            Obj co(out, 4);
            co.key("index");
            std::fprintf(out, "%zu", d.index);
            co.key("root");
            std::fprintf(out, "%d", d.c->root);
            co.key("tris");
            std::fprintf(out, "%zu", d.c->compTris.size());
            co.key("vtx");
            std::fprintf(out, "%zu", d.c->compVtx.size());
            co.key("edges");
            std::fprintf(out, "%zu", d.c->compEdges.size());
            co.key("clean");
            std::fputs(d.c->clean() ? "true" : "false", out);
            co.key("segmentOk");
            std::fputs(d.segmentOk ? "true" : "false", out);
            co.key("compVtx");
            putIntArray(out, d.c->compVtx);
            co.key("regionSet");
            writeRegionSet(out, d.rs, 6);
            co.end();
        }
        std::fputc('\n', out);
        putIndent(out, 2);
        std::fputc(']', out);
    }
    top.end();
    std::fputc('\n', out);

    if (out != stdout) std::fclose(out);
    return 0;
}
