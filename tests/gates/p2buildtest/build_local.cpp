// Rework-cycle driver: GRAFT G1 (STL+RegionSet via harness), G2 (J2 by TShape
// identity vs verts[]), G3 (*.expected.json + J6 against open==0).
// SPDX-License-Identifier: MIT

#include "../../harness/mesh_harness.hpp"
#include "../../harness/mesh_harness_refit.hpp"
#include "../../../src/refit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Standard_Failure.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <GProp_GProps.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>

using stl2step::harness::HarnessComponent;
using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;

namespace {

struct JsonValue {
    enum Kind { Null, Bool, Num, Str, Arr, Obj } kind = Null;
    bool b = false;
    double n = 0;
    std::string s;
    std::vector<JsonValue> arr;
    std::unordered_map<std::string, JsonValue> obj;
    const JsonValue& get(const char* k) const {
        static JsonValue none;
        auto it = obj.find(k);
        return it == obj.end() ? none : it->second;
    }
};

class JsonParser {
    const char* p = nullptr;
    void skip() {
        while (*p && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) ++p;
    }
    bool eat(char c) {
        skip();
        if (*p == c) { ++p; return true; }
        return false;
    }
    std::string parseStr() {
        skip();
        if (*p != '"') return "";
        ++p;
        std::string o;
        while (*p && *p != '"') {
            if (*p == '\\' && p[1]) { o += p[1]; p += 2; }
            else { o += *p++; }
        }
        if (*p == '"') ++p;
        return o;
    }
    JsonValue parseVal() {
        skip();
        JsonValue v;
        if (*p == '"') { v.kind = JsonValue::Str; v.s = parseStr(); return v; }
        if (*p == 't' && std::strncmp(p, "true", 4) == 0) {
            p += 4; v.kind = JsonValue::Bool; v.b = true; return v;
        }
        if (*p == 'f' && std::strncmp(p, "false", 5) == 0) {
            p += 5; v.kind = JsonValue::Bool; v.b = false; return v;
        }
        if (*p == 'n' && std::strncmp(p, "null", 4) == 0) {
            p += 4; return v;
        }
        if (*p == '[') {
            ++p; v.kind = JsonValue::Arr; skip();
            if (eat(']')) return v;
            while (true) {
                v.arr.push_back(parseVal());
                skip();
                if (eat(']')) break;
                eat(',');
            }
            return v;
        }
        if (*p == '{') {
            ++p; v.kind = JsonValue::Obj; skip();
            if (eat('}')) return v;
            while (true) {
                std::string k = parseStr();
                eat(':');
                v.obj[k] = parseVal();
                skip();
                if (eat('}')) break;
                eat(',');
            }
            return v;
        }
        v.kind = JsonValue::Num;
        char* end = nullptr;
        v.n = std::strtod(p, &end);
        p = end;
        return v;
    }
public:
    bool parse(const std::string& text, JsonValue& out) {
        p = text.c_str();
        out = parseVal();
        return true;
    }
};

static gp_XYZ readVec3(const JsonValue& a) {
    if (a.kind != JsonValue::Arr || a.arr.size() < 3) return gp_XYZ(0, 0, 0);
    return gp_XYZ(a.arr[0].n, a.arr[1].n, a.arr[2].n);
}
static gp_Ax3 readAx3(const JsonValue& o) {
    return gp_Ax3(gp_Pnt(readVec3(o.obj.at("loc"))), gp_Dir(readVec3(o.obj.at("dir"))),
                  gp_Dir(readVec3(o.obj.at("xdir"))));
}
static stl2step::refit::SurfType parseSurf(const std::string& s) {
    if (s == "cylinder") return stl2step::refit::SurfType::Cylinder;
    if (s == "cone") return stl2step::refit::SurfType::Cone;
    if (s == "sphere") return stl2step::refit::SurfType::Sphere;
    if (s == "torus") return stl2step::refit::SurfType::Torus;
    return stl2step::refit::SurfType::Plane;
}
static stl2step::refit::Origin parseOrigin(const std::string& s) {
    if (s == "cylGrow") return stl2step::refit::Origin::CylGrow;
    if (s == "filletStrip") return stl2step::refit::Origin::FilletStrip;
    return stl2step::refit::Origin::PlaneGrow;
}
static stl2step::refit::LoopRole parseRole(const std::string& s) {
    if (s == "inner") return stl2step::refit::LoopRole::Inner;
    if (s == "capLow") return stl2step::refit::LoopRole::CapLow;
    if (s == "capHigh") return stl2step::refit::LoopRole::CapHigh;
    return stl2step::refit::LoopRole::Outer;
}

static bool loadRegionSet(const JsonValue& root, stl2step::refit::RegionSet& rs) {
    rs.compRoot = static_cast<int>(root.obj.at("compRoot").n);
    rs.nIslands = static_cast<int>(root.obj.at("nIslands").n);
    for (const JsonValue& rj : root.obj.at("regions").arr) {
        stl2step::refit::Region r;
        r.id = static_cast<int>(rj.obj.at("id").n);
        r.type = parseSurf(rj.obj.at("type").s);
        r.origin = parseOrigin(rj.obj.at("origin").s);
        r.ax = readAx3(rj.obj.at("ax"));
        r.radius = rj.obj.at("radius").n;
        r.uMin = rj.obj.at("uMin").n;
        r.uMax = rj.obj.at("uMax").n;
        r.vMin = rj.obj.at("vMin").n;
        r.vMax = rj.obj.at("vMax").n;
        r.closed360 = rj.obj.at("closed360").b;
        r.outwardNormal = rj.obj.at("outwardNormal").b;
        for (const JsonValue& t : rj.obj.at("tris").arr)
            r.tris.push_back(static_cast<int>(t.n));
        for (const JsonValue& lj : rj.obj.at("loops").arr) {
            stl2step::refit::Loop lp;
            for (const JsonValue& ci : lj.obj.at("chainIdx").arr)
                lp.chainIdx.push_back(static_cast<int>(ci.n));
            for (const JsonValue& rv : lj.obj.at("reversed").arr)
                lp.reversed.push_back(static_cast<uint8_t>(rv.n));
            lp.role = parseRole(lj.obj.at("role").s);
            r.loops.push_back(lp);
        }
        r.maxVertexDev = rj.obj.at("maxVertexDev").n;
        r.rmsVertexDev = rj.obj.at("rmsVertexDev").n;
        r.chordSagitta = rj.obj.at("chordSagitta").n;
        r.nSides = static_cast<int>(rj.obj.at("nSides").n);
        r.dVolPredicted = rj.obj.at("dVolPredicted").n;
        r.maxVertexSnap = rj.obj.at("maxVertexSnap").n;
        r.filletNbrA = static_cast<int>(rj.obj.at("filletNbrA").n);
        r.filletNbrB = static_cast<int>(rj.obj.at("filletNbrB").n);
        rs.regions.push_back(r);
    }
    for (const JsonValue& cj : root.obj.at("chains").arr) {
        stl2step::refit::BoundaryChain ch;
        ch.regA = static_cast<int>(cj.obj.at("regA").n);
        ch.regB = static_cast<int>(cj.obj.at("regB").n);
        ch.islandA = static_cast<int>(cj.obj.at("islandA").n);
        ch.islandB = static_cast<int>(cj.obj.at("islandB").n);
        ch.tangent = cj.obj.at("tangent").b;
        ch.closedLoop = cj.obj.at("closedLoop").b;
        for (const JsonValue& e : cj.obj.at("meshEdges").arr)
            ch.meshEdges.push_back(static_cast<int>(e.n));
        for (const JsonValue& v : cj.obj.at("meshVerts").arr)
            ch.meshVerts.push_back(static_cast<int>(v.n));
        rs.chains.push_back(ch);
    }
    for (const JsonValue& t : root.obj.at("triRegion").arr)
        rs.triRegion.push_back(static_cast<int>(t.n));
    for (const JsonValue& t : root.obj.at("triIsland").arr)
        rs.triIsland.push_back(static_cast<int>(t.n));
    return true;
}

static std::string jnum(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

static const char* builtName(stl2step::refit::BuiltAs b) {
    switch (b) {
        case stl2step::refit::BuiltAs::Single: return "single";
        case stl2step::refit::BuiltAs::Seamed360: return "seamed360";
        case stl2step::refit::BuiltAs::TwoHalves: return "twoHalves";
        case stl2step::refit::BuiltAs::ExplodedToFacets: return "explodedToFacets";
        default: return "notBuilt";
    }
}

struct Verdict {
    bool valid = false, closed = false;
    int planes = 0, cylinders = 0, other = 0;
    double volume = 0;
    int j3Missing = 0;
    int slotShared = 0, twins = 0, strays = 0;
    int sameParamFalse = 0;
    double maxVertexTol = 0;
    std::vector<std::string> builtAs;
};

static Verdict runVerdict(const std::vector<TopoDS_Face>& faces,
                          const std::vector<TopoDS_Vertex>& verts,
                          const HarnessMesh& hm,
                          stl2step::refit::RegionSet& rs) {
    Verdict v;
    BRep_Builder B;
    TopoDS_Shell shell;
    B.MakeShell(shell);
    for (const TopoDS_Face& f : faces) B.Add(shell, f);

    try {
        BRepCheck_Analyzer chk(shell, Standard_True);
        v.valid = chk.IsValid() == Standard_True;
    } catch (const Standard_Failure&) {
        v.valid = false;
    }
    v.closed = BRep_Tool::IsClosed(shell) == Standard_True;

    GProp_GProps props;
    BRepGProp::VolumeProperties(shell, props);
    v.volume = props.Mass();

    for (TopExp_Explorer fx(shell, TopAbs_FACE); fx.More(); fx.Next()) {
        const TopoDS_Face& f = TopoDS::Face(fx.Current());
        BRepAdaptor_Surface s(f, Standard_False);
        if (s.GetType() == GeomAbs_Plane) v.planes++;
        else if (s.GetType() == GeomAbs_Cylinder) v.cylinders++;
        else v.other++;
        if (s.GetType() != GeomAbs_Plane) {
            for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
                TopoDS_Edge e = TopoDS::Edge(ex.Current());
                double f0, l0;
                Handle(Geom2d_Curve) c2d = BRep_Tool::CurveOnSurface(e, f, f0, l0);
                if (c2d.IsNull()) v.j3Missing++;
            }
        }
    }
    // F10: SameParameter flag census after buildFaces' forced pass.
    for (TopExp_Explorer ex(shell, TopAbs_EDGE); ex.More(); ex.Next()) {
        if (!BRep_Tool::SameParameter(TopoDS::Edge(ex.Current()))) v.sameParamFalse++;
    }

    // ---- J2: TShape identity of every face vertex against verts[] (GRAFT G2)
    {
        std::set<const void*> slotSet;
        for (const TopoDS_Vertex& sv : verts)
            if (!sv.IsNull()) slotSet.insert(sv.TShape().get());
        std::set<const void*> seen;
        for (const TopoDS_Face& f : faces) {
            for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
                const TopoDS_Vertex fv = TopoDS::Vertex(vx.Current());
                const void* ts = fv.TShape().get();
                if (!seen.insert(ts).second) continue;
                if (slotSet.count(ts)) { v.slotShared++; continue; }
                gp_Pnt P = BRep_Tool::Pnt(fv);
                bool coincident = false;
                for (size_t i = 0; i < verts.size(); ++i) {
                    if (verts[i].IsNull()) continue;
                    if (P.Distance(BRep_Tool::Pnt(verts[i])) <= 1e-7 * (1.0 + hm.diag)) {
                        coincident = true;
                        break;
                    }
                }
                if (coincident) v.twins++;
                else v.strays++;
            }
        }
    }
    for (const auto& sv : verts)
        if (!sv.IsNull())
            v.maxVertexTol = std::max(v.maxVertexTol, BRep_Tool::Tolerance(sv));
    for (const auto& reg : rs.regions) v.builtAs.push_back(builtName(reg.builtAs));
    return v;
}

static std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static std::string sibling(const std::string& regionset, const char* suffix) {
    std::string p = regionset;
    const std::string tag = ".regionset.json";
    auto pos = p.rfind(tag);
    if (pos != std::string::npos) p.replace(pos, tag.size(), suffix);
    else {
        auto dot = p.rfind('.');
        if (dot != std::string::npos) p = p.substr(0, dot) + suffix;
        else p += suffix;
    }
    return p;
}

void splitChainInHalf(stl2step::refit::RegionSet& rs, int ci) {
    if (ci < 0 || (size_t)ci >= rs.chains.size()) return;
    stl2step::refit::BoundaryChain ch = rs.chains[(size_t)ci];
    if (ch.meshEdges.size() < 4 || ch.meshVerts.size() < 3) return;
    const size_t nE = ch.meshEdges.size();
    const size_t mid = nE / 2;
    stl2step::refit::BoundaryChain a = ch, b = ch;
    a.closedLoop = false;
    b.closedLoop = false;
    a.meshEdges.assign(ch.meshEdges.begin(), ch.meshEdges.begin() + (std::ptrdiff_t)mid);
    b.meshEdges.assign(ch.meshEdges.begin() + (std::ptrdiff_t)mid, ch.meshEdges.end());
    a.meshVerts.assign(ch.meshVerts.begin(),
                       ch.meshVerts.begin() + (std::ptrdiff_t)std::min(mid + 1, ch.meshVerts.size()));
    b.meshVerts.assign(ch.meshVerts.begin() + (std::ptrdiff_t)std::min(mid, ch.meshVerts.size()),
                       ch.meshVerts.end());
    if (ch.closedLoop && !ch.meshVerts.empty()) b.meshVerts.push_back(ch.meshVerts.front());
    rs.chains[(size_t)ci] = a;
    const int ni = (int)rs.chains.size();
    rs.chains.push_back(b);
    for (auto& r : rs.regions) {
        for (auto& lp : r.loops) {
            for (size_t i = 0; i < lp.chainIdx.size(); i++) {
                if (lp.chainIdx[i] != ci) continue;
                const uint8_t rev = (i < lp.reversed.size()) ? lp.reversed[i] : 0;
                lp.chainIdx.insert(lp.chainIdx.begin() + (std::ptrdiff_t)i + 1, ni);
                if (i < lp.reversed.size())
                    lp.reversed.insert(lp.reversed.begin() + (std::ptrdiff_t)i + 1, rev);
                break;
            }
        }
    }
}

int runSelfTests(const std::string& rsPath, const std::string& stlPath) {
    JsonValue root;
    JsonParser().parse(slurp(rsPath), root);
    stl2step::refit::RegionSet base;
    if (!loadRegionSet(root, base)) {
        std::fprintf(stderr, "selftest: failed to load RegionSet\n");
        return 1;
    }
    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(stlPath, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        std::fprintf(stderr, "selftest: loadMesh: %s\n", err.c_str());
        return 1;
    }
    const HarnessComponent& comp = mesh.comps[0];
    stl2step::refit::MeshView mv;
    toMeshView(mesh, comp, mv);
    std::vector<TopoDS_Vertex> verts(comp.compVtx.size());
    BRep_Builder bb;
    for (size_t i = 0; i < comp.compVtx.size(); ++i)
        bb.MakeVertex(verts[i], gp_Pnt(mesh.pts[comp.compVtx[i]]), Precision::Confusion());

    int nFail = 0;
    auto call = [&](stl2step::refit::RegionSet rs, const char* tag, bool wantOk) {
        std::vector<TopoDS_Face> faces;
        bool ok = false;
        try {
            ok = stl2step::refit::buildFaces(mv, rs, verts, faces, nullptr);
        } catch (...) {
            std::printf("SELFTEST %s THREW\n", tag);
            nFail++;
            return;
        }
        const bool pass = wantOk ? ok : !ok;
        std::printf("SELFTEST %s buildFaces=%s wantOk=%s %s\n", tag, ok ? "true" : "false",
                    wantOk ? "true" : "false", pass ? "PASS" : "FAIL");
        if (!pass) nFail++;
    };

    {
        auto rs = base;
        if (!rs.regions.empty() && !rs.regions[0].tris.empty())
            rs.regions[0].tris[0] = 1 << 20;
        call(rs, "malformed_tri", false);
    }
    {
        auto rs = base;
        if (!rs.regions.empty() && !rs.regions[0].loops.empty() &&
            !rs.regions[0].loops[0].chainIdx.empty())
            rs.regions[0].loops[0].chainIdx[0] = 1 << 20;
        call(rs, "malformed_loop_chain", false);
    }
    {
        auto rs = base;
        if (!rs.chains.empty() && !rs.chains[0].meshVerts.empty())
            rs.chains[0].meshVerts[0] = 1 << 20;
        call(rs, "malformed_chain_vert", false);
    }

    {
        auto rs = base;
        std::vector<int> caps;
        for (const auto& r : rs.regions) {
            if (!r.closed360 || r.type != stl2step::refit::SurfType::Cylinder) continue;
            for (const auto& lp : r.loops) {
                if (lp.role != stl2step::refit::LoopRole::CapLow &&
                    lp.role != stl2step::refit::LoopRole::CapHigh)
                    continue;
                for (int ci : lp.chainIdx) caps.push_back(ci);
            }
        }
        std::sort(caps.begin(), caps.end());
        caps.erase(std::unique(caps.begin(), caps.end()), caps.end());
        for (int i = (int)caps.size() - 1; i >= 0; --i) splitChainInHalf(rs, caps[(size_t)i]);
        std::vector<TopoDS_Face> faces;
        bool ok = false;
        try {
            ok = stl2step::refit::buildFaces(mv, rs, verts, faces, nullptr);
        } catch (...) {
            std::printf("SELFTEST twoHalves THREW\n");
            return 1;
        }
        bool fired = false;
        for (const auto& r : rs.regions)
            if (r.builtAs == stl2step::refit::BuiltAs::TwoHalves) fired = true;
        std::printf("SELFTEST twoHalves buildFaces=%s fired=%s %s builtAs=[", ok ? "true" : "false",
                    fired ? "true" : "false", fired ? "PASS" : "FAIL");
        for (size_t i = 0; i < rs.regions.size(); i++) {
            if (i) std::printf(",");
            std::printf("\"%s\"", builtName(rs.regions[i].builtAs));
        }
        std::printf("]\n");
        if (!fired) nFail++;
    }

    std::printf("SELFTEST SUMMARY fail=%d\n", nFail);
    return nFail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "--selftest") {
        if (argc < 4) {
            std::fprintf(stderr, "usage: %s --selftest <name.regionset.json> <name.stl>\n",
                         argv[0]);
            return 2;
        }
        return runSelfTests(argv[2], argv[3]);
    }

    int component = 0;
    bool live = false;
    std::string sidecarPath;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--live") {
            live = true;
        } else if (a == "--component") {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "--component needs N\n");
                return 2;
            }
            component = std::atoi(argv[++i]);
        } else if (a == "--sidecar") {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "--sidecar needs a path\n");
                return 2;
            }
            sidecarPath = argv[++i];
        } else if (a[0] == '-') {
            std::fprintf(stderr, "unknown flag %s\n", a.c_str());
            return 2;
        } else {
            pos.push_back(a);
        }
    }
    if (pos.size() < 2) {
        std::fprintf(stderr,
                     "usage: %s <name.regionset.json> <name.stl> [name.expected.json]\n"
                     "       %s --live [--component N] [--sidecar FILE] <rs.json> <stl>\n"
                     "       %s --selftest <name.regionset.json> <name.stl>\n",
                     argv[0], argv[0], argv[0]);
        return 2;
    }

    JsonValue root;
    JsonParser().parse(slurp(pos[0]), root);
    stl2step::refit::RegionSet rs;
    if (!loadRegionSet(root, rs)) {
        std::fprintf(stderr, "failed to load RegionSet\n");
        return 1;
    }

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(pos[1], 1.0, 0.0, 0.0, mesh, err)) {
        std::fprintf(stderr, "loadMesh: %s\n", err.c_str());
        return 1;
    }
    if (mesh.comps.empty()) {
        std::fprintf(stderr, "no components\n");
        return 1;
    }
    if (component < 0 || (size_t)component >= mesh.comps.size()) {
        std::fprintf(stderr, "--component %d out of range (0..%d)\n", component,
                     (int)mesh.comps.size() - 1);
        return 1;
    }
    const HarnessComponent& comp = mesh.comps[(size_t)component];

    stl2step::refit::MeshView mv;
    toMeshView(mesh, comp, mv);

    std::vector<TopoDS_Vertex> verts(comp.compVtx.size());
    BRep_Builder bb;
    for (size_t i = 0; i < comp.compVtx.size(); ++i)
        bb.MakeVertex(verts[i], gp_Pnt(mesh.pts[comp.compVtx[i]]), Precision::Confusion());

    std::vector<std::string> warns;
    auto warnFn = [&](const std::string& m) { warns.push_back(m); };
    std::vector<TopoDS_Face> faces;
    bool built = false;
    try {
        built = stl2step::refit::buildFaces(mv, rs, verts, faces, warnFn);
    } catch (...) {
        built = false;
    }

    Verdict v = runVerdict(faces, verts, mesh, rs);

    // J6: mesh open==0 (not clean()) ⇒ shell closed.
    bool j6 = (comp.open != 0) || v.closed;

    JsonValue sidecar;
    if (!sidecarPath.empty()) {
        std::string txt = slurp(sidecarPath);
        if (!txt.empty()) JsonParser().parse(txt, sidecar);
    }

    JsonValue exp;
    std::string expPath = (pos.size() >= 3) ? pos[2] : (live ? std::string() : sibling(pos[0], ".expected.json"));
    bool haveExp = false;
    bool liveEscalate = false;
    bool expBuildFaces = true;
    bool haveExpBuildFaces = false;
    if (!expPath.empty()) {
        std::string txt = slurp(expPath);
        if (!txt.empty()) {
            JsonParser().parse(txt, exp);
            haveExp = exp.kind == JsonValue::Obj && exp.get("faceCount").kind == JsonValue::Num;
        }
    }
    // Live gate: consume sidecar live[component] faceCount/census/volumeBudget.
    // ESCALATE rows keep the budget but do not hard-fail census (stale faceCount).
    const JsonValue* liveRow = nullptr;
    if (sidecar.kind == JsonValue::Obj) {
        for (const auto& row : sidecar.get("live").arr) {
            if ((int)row.get("component").n == component) {
                liveRow = &row;
                break;
            }
        }
    }
    if (liveRow) {
        liveEscalate = liveRow->get("disposition").s == "ESCALATE";
        if (liveRow->get("buildFaces").kind == JsonValue::Bool) {
            haveExpBuildFaces = true;
            expBuildFaces = liveRow->get("buildFaces").b;
        }
        if (liveRow->get("faceCount").kind == JsonValue::Num) exp = *liveRow;
        if (!liveEscalate)
            haveExp = exp.get("faceCount").kind == JsonValue::Num;
        else
            haveExp = false;
        // Placeholder live[] (faceCount 0 + empty census) is not a hard target.
        if (haveExp && (int)exp.get("faceCount").n == 0) {
            const JsonValue& sc0 = exp.get("surfaceCensus");
            if ((int)sc0.get("plane").n == 0 && (int)sc0.get("cylinder").n == 0 &&
                (int)sc0.get("planar_facets").n == 0)
                haveExp = false;
        }
    }

    int expFaces = haveExp ? (int)exp.get("faceCount").n : -1;
    int expCyl = 0, expPln = 0;
    if (haveExp) {
        const JsonValue& sc = exp.get("surfaceCensus");
        expCyl = (int)sc.get("cylinder").n;
        expPln = (int)sc.get("plane").n + (int)sc.get("planar_facets").n;
    }
    // D4.5: budget is Σ|dVolPredicted| (magnitudes); signed sum is RESULT only.
    double absSum = 0;
    for (const auto& r : rs.regions) absSum += std::fabs(r.dVolPredicted);
    double computedBudget = std::max(1e-4 * std::fabs(comp.vol), 3.0 * absSum);
    double budget = 0;
    if (liveRow && liveRow->get("volumeBudgetMM3").kind == JsonValue::Num)
        budget = liveRow->get("volumeBudgetMM3").n;
    else if (haveExp)
        budget = exp.get("volumeBudgetMM3").n;
    if (budget <= 0)
        budget = computedBudget;
    else
        budget = std::max(budget, computedBudget);
    // Open shells with stepVol=0 are not a pass (vacuous). Closed shells also
    // need |Δ| within D4.5 budget vs |meshVol|.
    bool volNonZero = std::fabs(v.volume) > Precision::Confusion();
    bool volOk = volNonZero;
    if (v.closed)
        volOk = volNonZero && std::fabs(v.volume - std::fabs(comp.vol)) <= budget + 1e-9;

    bool facesOk = !haveExp || (int)faces.size() == expFaces;
    bool censusOk = !haveExp || (v.cylinders == expCyl && v.planes == expPln);
    // All-facet live intent (S09c1): P1 may still emit planar patches; a
    // valid 0-cylinder reducing shell is the contract, not faceCount==nTri.
    if (haveExp && liveRow) {
        const JsonValue& sc = exp.get("surfaceCensus");
        int expFacet = (int)sc.get("planar_facets").n;
        if (expCyl == 0 && expFacet == (int)comp.compTris.size() && v.cylinders == 0 &&
            v.valid && j6) {
            facesOk = true;
            censusOk = true;
        }
        // Sidecar cylinder count is the recoverable-surface contract (S05 2/2).
        // Dump plane/faceCount may differ from the live[] target after islands.
        if (expCyl > 0 && v.cylinders >= (expCyl + 1) / 2 && v.valid && j6) {
            censusOk = true;
            facesOk = true;
            if (!volOk) {
                const double loose = 0.15 * std::fabs(comp.vol) + budget;
                if (std::fabs(v.volume - std::fabs(comp.vol)) <= loose) volOk = true;
            }
        }
        // Recovering extra analytics vs a faceted/R1 sidecar target is a pass.
        if (v.cylinders > expCyl && v.valid && j6 && volOk) {
            censusOk = true;
            facesOk = true;
        }
    }
    bool builtAsOk = true;
    if (haveExp) {
        const JsonValue& regs = exp.get("regions");
        for (size_t i = 0; i < regs.arr.size() && i < rs.regions.size(); i++) {
            std::string want = regs.arr[i].get("builtAs").s;
            if (want.size() && want != builtName(rs.regions[i].builtAs)) builtAsOk = false;
        }
    }

    // J1: face-count reducing, one shell from one component.
    bool j1 = (int)faces.size() <= (int)comp.compTris.size() && !faces.empty();
    int nClosed360 = 0, nSeamed = 0, nHalves = 0, nExploded360 = 0;
    int nExplodedHole = 0, nExplodedCyl = 0, nBuiltCyl = 0;
    for (const auto& r : rs.regions) {
        if (r.type == stl2step::refit::SurfType::Cylinder) {
            if (r.builtAs == stl2step::refit::BuiltAs::ExplodedToFacets) nExplodedCyl++;
            else if (r.builtAs == stl2step::refit::BuiltAs::Single ||
                     r.builtAs == stl2step::refit::BuiltAs::Seamed360 ||
                     r.builtAs == stl2step::refit::BuiltAs::TwoHalves)
                nBuiltCyl++;
        }
        if (!r.closed360) continue;
        nClosed360++;
        if (r.builtAs == stl2step::refit::BuiltAs::Seamed360) nSeamed++;
        else if (r.builtAs == stl2step::refit::BuiltAs::TwoHalves) nHalves++;
        else if (r.builtAs == stl2step::refit::BuiltAs::ExplodedToFacets) {
            nExploded360++;
            if (!r.outwardNormal) nExplodedHole++;
        }
    }

    bool sidecarRecoverable360 = false;
    bool sidecarAllowsExplode = false;
    int nSidecarPartialCyl = 0;
    if (sidecar.kind == JsonValue::Obj) {
        for (const auto& rec : sidecar.get("recoverable").arr) {
            if (rec.get("type").s != "cylinder") continue;
            if (rec.get("closed360").b) sidecarRecoverable360 = true;
            else nSidecarPartialCyl++;
        }
        const JsonValue& role = sidecar.get("rLadderRole");
        if (role.s.find("R1") != std::string::npos || role.s.find("R2") != std::string::npos)
            sidecarAllowsExplode = true;
        for (const auto& mf : sidecar.get("mustRemainFaceted").arr) {
            if (mf.get("type").s == "hole") sidecarAllowsExplode = true;
        }
    }
    bool explodedRecoverableHole = false;
    bool explodedRecoverableSurface = false;
    if (live) {
        if (nExplodedHole > 0 && sidecarRecoverable360) explodedRecoverableHole = true;
        if (nExplodedHole > 0 && sidecarPath.empty()) explodedRecoverableHole = true;
        if (nExploded360 > 0 && sidecarRecoverable360 && !sidecarAllowsExplode)
            explodedRecoverableHole = explodedRecoverableHole || (nExplodedHole > 0);
        if (sidecarAllowsExplode && !sidecarRecoverable360) explodedRecoverableHole = false;
        // Partial recoverable cylinders (S05 slot ends): exploding all of them
        // and emitting facets is a false pass. Fail unless some cylinder remains.
        if (nSidecarPartialCyl > 0 && v.cylinders == 0 && !sidecarAllowsExplode)
            explodedRecoverableSurface = true;
    }
    // Catch 10^3–10^4 mm TShape corruption (adjudication measured 10884 mm
    // on a 20 mm cube). Floor 25 mm still fails those; sagitta stays below.
    const double tolBudget = std::max(0.05 * (mesh.diag > 0 ? mesh.diag : 1.0), 25.0);
    bool tolOk = v.maxVertexTol <= tolBudget + 1e-12;

    bool builtOk = built;
    if (haveExpBuildFaces && !expBuildFaces)
        builtOk = true;  // R2: live[] buildFaces=false; facets still judged below
    bool gate = builtOk && v.valid && j6 && v.twins == 0 && v.j3Missing == 0 && volOk &&
                facesOk && censusOk && builtAsOk;
    if (live) gate = gate && j1 && !explodedRecoverableHole && !explodedRecoverableSurface &&
                     tolOk;

    std::printf("{\"buildFaces\":%s,\"valid\":%s,\"closed\":%s,"
                "\"openE\":%d,\"j6\":%s,\"j1\":%s,\"nTri\":%d,\"faceCount\":%zu,"
                "\"planes\":%d,\"cylinders\":%d,\"volume\":%s,"
                "\"volumeBudget\":%s,\"volumeOk\":%s,\"dVolAbsSum\":%s,"
                "\"slotShared\":%d,\"twins\":%d,\"strays\":%d,"
                "\"j3Missing\":%d,\"sameParamFalse\":%d,"
                "\"closed360\":%d,\"seamed360\":%d,\"twoHalves\":%d,\"exploded360\":%d,"
                "\"explodedRecoverableHole\":%s,\"explodedRecoverableSurface\":%s,"
                "\"maxVertexTol\":%s,\"tolOk\":%s,"
                "\"facesOk\":%s,\"censusOk\":%s,\"builtAsOk\":%s,\"ok\":%s,"
                "\"builtAs\":[",
                built ? "true" : "false", v.valid ? "true" : "false",
                v.closed ? "true" : "false", comp.open, j6 ? "true" : "false",
                j1 ? "true" : "false", (int)comp.compTris.size(),
                faces.size(), v.planes, v.cylinders, jnum(v.volume).c_str(),
                jnum(budget).c_str(), volOk ? "true" : "false", jnum(absSum).c_str(),
                v.slotShared, v.twins, v.strays, v.j3Missing, v.sameParamFalse,
                nClosed360, nSeamed, nHalves, nExploded360,
                explodedRecoverableHole ? "true" : "false",
                explodedRecoverableSurface ? "true" : "false",
                jnum(v.maxVertexTol).c_str(), tolOk ? "true" : "false",
                facesOk ? "true" : "false", censusOk ? "true" : "false",
                builtAsOk ? "true" : "false", gate ? "true" : "false");
    for (size_t i = 0; i < v.builtAs.size(); i++) {
        if (i) std::printf(",");
        std::printf("\"%s\"", v.builtAs[i].c_str());
    }
    std::printf("],\"warnings\":[");
    for (size_t i = 0; i < warns.size(); i++) {
        if (i) std::printf(",");
        std::printf("\"");
        for (char c : warns[i]) {
            if (c == '"' || c == '\\') std::putchar('\\');
            if (c == '\n') { std::printf("\\n"); continue; }
            std::putchar(c);
        }
        std::printf("\"");
    }
    std::printf("]}\n");
    return gate ? 0 : 1;
}
