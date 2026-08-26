// stl2step engine - mesh (STL) -> B-Rep solid (STEP).
//
// Pipeline:
//   read STL (RWStl: binary + ASCII) -> weld duplicate vertices -> split the mesh
//   into MANIFOLD components (connectivity across cleanly-shared edges only, so
//   solids glued along a non-manifold edge separate into clean bodies instead of
//   one dirty component that needs sewing) -> per component: build vertices, then
//   shared edges, then one planar face per triangle -- each phase runs across all
//   cores (independent OCCT objects into pre-assigned slots, no locks); closed,
//   consistently-wound components become oriented solids; dirty components (open
//   edges / flipped facets) fall back to a BRepBuilderAPI_Sewing + ShapeFix_Shell
//   repair pass -> ShapeUpgrade_UnifySameDomain (parallel across bodies) merges
//   coplanar facets into single faces (flats collapse to clean faces; curved
//   surfaces stay faceted -- a mesh carries no analytic geometry to recover) ->
//   vertex/edge tolerances fitted to true planar deviation -> validity check +
//   volume run on a worker thread OVERLAPPED with the STEP write (both only read
//   the finished shape; in the rare invalid case the shape is ShapeFix'd and the
//   file rewritten) -> optional verify by re-reading the file.
//
// This file is the whole engine. The public interface is include/stl2step.hpp.
//
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"
#include "refit.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <RWStl.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_State.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

#if defined(__GNUC__) || defined(__clang__)
#  define STL2STEP_PRINTF(fmtIdx, argIdx) __attribute__((format(printf, fmtIdx, argIdx)))
#else
#  define STL2STEP_PRINTF(fmtIdx, argIdx)
#endif

namespace stl2step {

namespace refit {
void fitAnalyticTolerances(const TopoDS_Shape& shape);
}

namespace {

namespace fs = std::filesystem;

constexpr double kPi = 3.14159265358979323846;  // avoids MSVC's _USE_MATH_DEFINES dance

// ---------------------------------------------------------------- small helpers

std::string fmtInt(long long n) {
    std::string s = std::to_string(n);
    for (int i = (int)s.size() - 3; i > 0; i -= 3) s.insert(i, ",");
    return s;
}

struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    double lap() {
        auto now = std::chrono::steady_clock::now();
        double s = std::chrono::duration<double>(now - t0).count();
        t0 = now;
        return s;
    }
};

std::string jsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if ((unsigned char)c >= 0x20) o += c;
    }
    return o;
}

// Chunked parallel-for over [0, n). Serial below 4096 items (thread spawn costs
// more than tiny workloads); each worker claims 1024-item chunks off an atomic.
template <typename F>
void parallelFor(unsigned threads, size_t n, F&& fn) {
    if (threads <= 1 || n < 4096) {
        for (size_t i = 0; i < n; i++) fn(i);
        return;
    }
    unsigned k = (unsigned)std::min<size_t>(threads, (n + 4095) / 4096);
    std::atomic<size_t> next{ 0 };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < k; t++)
        pool.emplace_back([&]() {
            for (;;) {
                size_t base = next.fetch_add(1024);
                if (base >= n) break;
                size_t end = std::min(n, base + 1024);
                for (size_t i = base; i < end; i++) fn(i);
            }
        });
    for (auto& th : pool) th.join();
}

// Exact-bit vertex key (normalises -0.0 -> +0.0); tolerance welding quantises instead.
struct VKey {
    uint64_t x, y, z;
    bool operator==(const VKey& o) const { return x == o.x && y == o.y && z == o.z; }
};
struct VKeyHash {
    size_t operator()(const VKey& k) const {
        uint64_t h = k.x * 0x9E3779B97F4A7C15ull;
        h ^= k.y + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
        h ^= k.z + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
        return (size_t)h;
    }
};
uint64_t dblBits(double v) {
    if (v == 0.0) v = 0.0;                       // collapse -0.0
    uint64_t u; memcpy(&u, &v, 8); return u;
}

struct UnionFind {
    std::vector<int> p;
    explicit UnionFind(int n) : p(n) { for (int i = 0; i < n; i++) p[i] = i; }
    int find(int a) { while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; } return a; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
};

int countShapes(const TopoDS_Shape& s, TopAbs_ShapeEnum what,
                TopAbs_ShapeEnum avoid = TopAbs_SHAPE) {
    int n = 0;
    for (TopExp_Explorer ex(s, what, avoid); ex.More(); ex.Next()) n++;
    return n;
}

double shapeVolume(const TopoDS_Shape& s) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return props.Mass();
}

long long fileSize(const std::string& path) {
    std::error_code ec;
    auto sz = fs::file_size(fs::path(path), ec);
    return ec ? -1 : (long long)sz;
}

// Raise vertex/edge tolerances to the true measured deviation from each planar
// face. Facet normals in an STL carry float32 noise, so merging coplanar facets
// leaves boundary vertices sitting off the merged plane by more than OCCT's
// default 1e-7 -- the geometry is honest, only the stored tolerance is too
// optimistic. Measuring is linear and a no-op on exact geometry; it replaces a
// full ShapeFix pass. (Edges here are straight lines, so their max deviation
// from a plane is at an endpoint.)
void fitPlanarTolerances(const TopoDS_Shape& shape) {
    BRep_Builder B;
    for (TopExp_Explorer fx(shape, TopAbs_FACE); fx.More(); fx.Next()) {
        const TopoDS_Face& f = TopoDS::Face(fx.Current());
        BRepAdaptor_Surface s(f, Standard_False);
        if (s.GetType() != GeomAbs_Plane) continue;
        gp_Pln pln = s.Plane();
        for (TopExp_Explorer vx(f, TopAbs_VERTEX); vx.More(); vx.Next()) {
            const TopoDS_Vertex& v = TopoDS::Vertex(vx.Current());
            double d = pln.Distance(BRep_Tool::Pnt(v));
            if (d > BRep_Tool::Tolerance(v))
                B.UpdateVertex(v, d * 1.001 + Precision::Confusion());
        }
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) {
            const TopoDS_Edge& e = TopoDS::Edge(ex.Current());
            double d = 0;
            for (TopExp_Explorer vx(e, TopAbs_VERTEX); vx.More(); vx.Next())
                d = std::max(d, pln.Distance(BRep_Tool::Pnt(TopoDS::Vertex(vx.Current()))));
            if (d > BRep_Tool::Tolerance(e))
                B.UpdateEdge(e, d * 1.001 + Precision::Confusion());
        }
    }
}

// ---------------------------------------------------------------- converter

// Holds all per-conversion state (progress sink, collected warnings, io lock),
// so convert() is reentrant and the internal note()/warn() call sites read
// exactly as they did in the original single-shot program.
struct Converter {
    const Options& opt;
    LogCallback    log;
    std::vector<std::string> warnings;
    std::mutex io;                    // guards the log callback + warnings

    explicit Converter(const Options& o, LogCallback l) : opt(o), log(std::move(l)) {}

    void note(const char* fmt, ...) STL2STEP_PRINTF(2, 3) {
        if (!log) return;
        char buf[1400];
        va_list ap; va_start(ap, fmt);
        vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        std::string s(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        std::lock_guard<std::mutex> lock(io);
        log(Severity::Info, s);
    }
    void warn(const std::string& msg) {
        std::lock_guard<std::mutex> lock(io);
        warnings.push_back(msg);
        if (log) log(Severity::Warning, msg);
    }
    // Terminal failure: fill the error fields and hand back the result.
    Result fail(Result& r, const std::string& err) {
        r.ok = false;
        r.exitCode = 1;
        r.error = err;
        r.warnings = warnings;
        if (log) log(Severity::Error, err);
        return r;
    }

    Result run();
};

Result Converter::run() {
    Result r;

    // Resolve paths (std::filesystem keeps this portable across separators).
    const std::string schema = schemaName(opt.schema);
    double scale = opt.scale * (opt.inchInput ? 25.4 : 1.0);
    double weldTol = opt.weldTol, sewTolArg = opt.sewTol, unifyAngleDeg = opt.unifyAngleDeg;
    bool unify = opt.unify, makeSolids = opt.makeSolids, verify = opt.verify, forceSew = opt.forceSew;
    bool smooth = opt.smooth;

    std::string inPath = opt.input;
    std::string outPath = opt.output;
    if (inPath.empty()) return fail(r, "no input path given");
    if (outPath.empty()) {
        fs::path p(inPath);
        p.replace_extension(".step");
        outPath = p.string();
    }
    std::string stem = opt.productName.empty() ? fs::path(outPath).stem().string()
                                               : opt.productName;
    r.input = inPath;
    r.output = outPath;

    unsigned hw = opt.threads > 0 ? (unsigned)opt.threads
                                  : std::thread::hardware_concurrency();
    if (hw == 0) hw = 4;

    // Silence OCCT's own console chatter; convert crashes/signals into C++ exceptions.
    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);

    Timer timer, total;
    int solidsOut = 0, openShellsOut = 0, facesBefore = 0, facesAfter = 0;
    int nTri = 0, nVert = 0, nComp = 0;
    double meshVolume = 0, brepVolume = 0, stepVolume = 0, volDeltaPct = -1;
    bool watertight = true;

    try {
        // ---- 1. read -----------------------------------------------------------
        note("stl2step: %s -> %s\n", inPath.c_str(), outPath.c_str());
        Handle(Poly_Triangulation) mesh = RWStl::ReadFile(inPath.c_str());
        if (mesh.IsNull() || mesh->NbTriangles() < 1)
            return fail(r, "unreadable or empty STL");
        int rawVerts = mesh->NbNodes();
        nTri = mesh->NbTriangles();

        // ---- 2. weld vertices (exact by default, tolerance grid via weldTol) ---
        std::vector<gp_XYZ> pts;             // welded, scaled points
        std::vector<int> remap(rawVerts + 1);
        {
            pts.reserve(rawVerts);
            std::unordered_map<VKey, int, VKeyHash> seen;
            seen.reserve(rawVerts * 2);
            for (int i = 1; i <= rawVerts; i++) {
                gp_XYZ p = mesh->Node(i).XYZ() * scale;
                VKey k;
                if (weldTol > 0) {
                    k = { (uint64_t)(int64_t)llround(p.X() / weldTol),
                          (uint64_t)(int64_t)llround(p.Y() / weldTol),
                          (uint64_t)(int64_t)llround(p.Z() / weldTol) };
                } else {
                    k = { dblBits(p.X()), dblBits(p.Y()), dblBits(p.Z()) };
                }
                auto it = seen.find(k);
                if (it == seen.end()) {
                    seen.emplace(k, (int)pts.size());
                    remap[i] = (int)pts.size();
                    pts.push_back(p);
                } else {
                    remap[i] = it->second;
                }
            }
        }
        nVert = (int)pts.size();

        gp_XYZ bbMin(1e300, 1e300, 1e300), bbMax(-1e300, -1e300, -1e300);
        for (auto& p : pts) {
            bbMin.SetX(std::min(bbMin.X(), p.X())); bbMax.SetX(std::max(bbMax.X(), p.X()));
            bbMin.SetY(std::min(bbMin.Y(), p.Y())); bbMax.SetY(std::max(bbMax.Y(), p.Y()));
            bbMin.SetZ(std::min(bbMin.Z(), p.Z())); bbMax.SetZ(std::max(bbMax.Z(), p.Z()));
        }
        double diag = (bbMax - bbMin).Modulus();
        double sewTol = sewTolArg > 0 ? sewTolArg
                                      : std::min(std::max(1e-6, diag * 1e-5), 0.5);

        note("  read      %s triangles, %s vertices (%s welded), bbox %.2f x %.2f x %.2f mm  [%.2fs]\n",
             fmtInt(nTri).c_str(), fmtInt(nVert).c_str(), fmtInt(rawVerts - nVert).c_str(),
             bbMax.X() - bbMin.X(), bbMax.Y() - bbMin.Y(), bbMax.Z() - bbMin.Z(), timer.lap());
        if (nTri > 500000)
            warn("very large mesh (" + fmtInt(nTri) + " triangles) -- conversion is memory-heavy "
                 "and the STEP will be large; consider decimating the mesh first");

        // ---- 3. edge stats + manifold component split ---------------------------
        struct Tri { int a, b, c; };
        std::vector<Tri> tris;
        tris.reserve(nTri);
        for (int i = 1; i <= nTri; i++) {
            Poly_Triangle t = mesh->Triangle(i);
            int a, b, c; t.Get(a, b, c);
            a = remap[a]; b = remap[b]; c = remap[c];
            if (a == b || b == c || a == c) continue;   // index-degenerate
            tris.push_back({ a, b, c });
        }
        mesh.Nullify();

        auto ekey = [](int a, int b) {
            return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                         : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
        };

        // Pass A: global edge use counts. Pass B: union triangles ONLY across edges
        // used exactly twice, so components split at non-manifold junctions: two
        // solids glued along an edge become two clean manifold bodies (direct path)
        // instead of one dirty component that would need the slow sewing repair.
        UnionFind uf((int)tris.size());
        {
            struct EdgeStat { uint32_t cnt = 0; int t1 = -1, t2 = -1; };
            std::unordered_map<uint64_t, EdgeStat> edges;
            edges.reserve(tris.size() * 2);
            for (int t = 0; t < (int)tris.size(); t++) {
                int v[4] = { tris[t].a, tris[t].b, tris[t].c, tris[t].a };
                for (int s = 0; s < 3; s++) {
                    EdgeStat& e = edges[ekey(v[s], v[s + 1])];
                    e.cnt++;
                    if (e.t1 < 0) e.t1 = t;
                    else if (e.t2 < 0) e.t2 = t;
                }
            }
            for (auto& kv : edges)
                if (kv.second.cnt == 2) uf.unite(kv.second.t1, kv.second.t2);
        }

        // Pass C: per component -- triangle list, local vertex table, edge table
        // with per-triangle edge references (this pre-computed topology is what
        // lets the construction phase run lock-free across all cores), edge
        // health, signed volume.
        struct CompStat {
            int open = 0, conflict = 0, nonManifold = 0;
            double vol = 0;
            std::vector<int> tris;                    // global triangle ids
            std::vector<int> vtx;                     // local vertex -> global index
            std::vector<std::pair<int, int>> edges;   // local edge -> (local vLo, vHi), global-lo first
            std::vector<std::array<int, 3>> triEdges; // per local tri: 3 edge ids
            std::vector<uint8_t> triDirs;             // bit s set = side s runs global-lo -> hi
        };
        std::unordered_map<int, CompStat> comps;
        {
            struct EdgeRec { uint32_t cnt = 0, fwd = 0; };
            std::unordered_map<int, std::unordered_map<uint64_t, int>> edgeIds;
            std::unordered_map<int, std::vector<EdgeRec>> edgeRecs;
            std::unordered_map<int, std::unordered_map<int, int>> vtxIds;
            for (int t = 0; t < (int)tris.size(); t++) {
                int root = uf.find(t);
                CompStat& c = comps[root];
                auto& em = edgeIds[root];
                auto& er = edgeRecs[root];
                auto& vm = vtxIds[root];
                auto localV = [&](int g) {
                    auto it = vm.find(g);
                    if (it != vm.end()) return it->second;
                    int id = (int)c.vtx.size();
                    vm.emplace(g, id);
                    c.vtx.push_back(g);
                    return id;
                };
                int v[4] = { tris[t].a, tris[t].b, tris[t].c, tris[t].a };
                std::array<int, 3> te{};
                uint8_t dirs = 0;
                for (int s = 0; s < 3; s++) {
                    int gu = v[s], gv = v[s + 1];
                    uint64_t k = ekey(gu, gv);
                    auto it = em.find(k);
                    int id;
                    if (it == em.end()) {
                        id = (int)c.edges.size();
                        em.emplace(k, id);
                        c.edges.push_back({ localV(std::min(gu, gv)), localV(std::max(gu, gv)) });
                        er.push_back({});
                    } else {
                        id = it->second;
                    }
                    er[id].cnt++;
                    if (gu < gv) { er[id].fwd++; dirs |= (uint8_t)(1 << s); }
                    te[s] = id;
                }
                c.tris.push_back(t);
                c.triEdges.push_back(te);
                c.triDirs.push_back(dirs);
                const gp_XYZ &A = pts[tris[t].a], &B = pts[tris[t].b], &C = pts[tris[t].c];
                c.vol += (A.X() * (B.Y() * C.Z() - B.Z() * C.Y())
                        - A.Y() * (B.X() * C.Z() - B.Z() * C.X())
                        + A.Z() * (B.X() * C.Y() - B.Y() * C.X())) / 6.0;
            }
            for (auto& kv : edgeRecs) {
                CompStat& c = comps.at(kv.first);
                for (auto& e : kv.second) {
                    if (e.cnt == 1) c.open++;
                    else if (e.cnt > 2) c.nonManifold++;
                    else if (e.fwd != 1) c.conflict++;   // both uses same direction
                }
            }
        }
        nComp = (int)comps.size();
        meshVolume = 0;
        for (auto& kv : comps) {
            if (kv.second.open || kv.second.nonManifold || kv.second.conflict)
                watertight = false;
            meshVolume += std::fabs(kv.second.vol);
        }

        // ---- 4+5. per component: parallel build -> shell -> solid ---------------
        std::vector<int> order;                 // big components first, deterministic
        order.reserve(comps.size());
        for (auto& kv : comps) order.push_back(kv.first);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            size_t sa = comps.at(a).tris.size(), sb = comps.at(b).tris.size();
            return sa != sb ? sa > sb : a < b;
        });
        unsigned nOuter = (unsigned)std::min<size_t>({ (size_t)hw, order.size(), (size_t)10 });
        unsigned nInner = std::max(1u, hw / std::max(1u, nOuter));

        auto isClean = [&](const CompStat& cs) {
            return !forceSew && cs.open == 0 && cs.conflict == 0 && cs.nonManifold == 0;
        };
        auto segmentSummaryStderr = [](int root, const refit::RegionSet& rs) {
            std::fprintf(stderr,
                         "engine segment root=%d regions=%zu rejected=%zu planes=%d "
                         "cylinders=%d fillets=%d facetIslands=%d\n",
                         root, rs.regions.size(), rs.rejected.size(), rs.stats.planes,
                         rs.stats.cylinders, rs.stats.fillets, rs.stats.facetIslands);
            for (const refit::Region& r : rs.regions) {
                const char* ty = "plane";
                if (r.type == refit::SurfType::Cylinder) ty = "cylinder";
                else if (r.type == refit::SurfType::Cone) ty = "cone";
                else if (r.type == refit::SurfType::Sphere) ty = "sphere";
                else if (r.type == refit::SurfType::Torus) ty = "torus";
                std::fprintf(stderr, "  id=%d type=%s tris=%zu radius=%.6g closed360=%d\n",
                             r.id, ty, r.tris.size(), r.radius, r.closed360 ? 1 : 0);
            }
            std::fputs("  cylinder radii:", stderr);
            bool first = true;
            for (const refit::Region& r : rs.regions) {
                if (r.type != refit::SurfType::Cylinder) continue;
                if (!first) std::fputc(',', stderr);
                first = false;
                char buf[64];
                std::snprintf(buf, sizeof buf, "%.17g", r.radius);
                std::fputs(buf, stderr);
            }
            std::fputc('\n', stderr);
        };
        auto fillMeshView = [&](const CompStat& cs, refit::MeshView& mv) {
            mv.pts = pts.data();
            mv.tris = reinterpret_cast<const int(*)[3]>(tris.data());
            mv.compTris = cs.tris.data();
            mv.compVtx = cs.vtx.data();
            mv.compEdges = cs.edges.data();
            mv.triEdges = cs.triEdges.data();
            mv.triDirs = cs.triDirs.data();
            mv.nTri = cs.tris.size();
            mv.nVtx = cs.vtx.size();
            mv.nEdge = cs.edges.size();
            mv.diag = diag;
            mv.weldTol = weldTol;
            mv.sewTol = sewTol;
        };

        std::unordered_map<int, refit::RegionSet> refitPlans;
        refit::RefitStats refitTotals;
        int smoothSkippedComponents = 0;
        double dVolPredSigned = 0, dVolPredAbs = 0;
        if (smooth && !forceSew) {
            refit::SegmentParams segp;
            segp.epsPlane = opt.smoothTolMM;
            segp.thetaPlaneDeg = opt.smoothAngleDeg;
            segp.doFillets = opt.smoothFillets;
            std::mutex planMu;
            parallelFor(hw, order.size(), [&](size_t idx) {
                int root = order[idx];
                const CompStat& cs = comps.at(root);
                if (!isClean(cs)) {
                    std::lock_guard<std::mutex> lk(planMu);
                    smoothSkippedComponents++;
                    return;
                }
                refit::MeshView mv{};
                fillMeshView(cs, mv);
                refit::RegionSet rs;
                rs.compRoot = root;
                refit::segment(mv, segp, rs, nullptr);
                std::lock_guard<std::mutex> lk(planMu);
                if (const char* e = std::getenv("STL2STEP_SEGMENT_SUMMARY");
                    e && e[0] && e[0] != '0')
                    segmentSummaryStderr(root, rs);
                refitPlans.emplace(root, std::move(rs));
            });
            if (smoothSkippedComponents > 0)
                note("  smooth    %d component(s) skipped (dirty mesh, repaired path)\n",
                     smoothSkippedComponents);
        }

        struct CompOut {
            std::vector<TopoDS_Shape> parts;
            std::vector<TopoDS_Face> keepShapes;
            refit::RefitStats refitSt{};
            double refitDVolAbs = 0;
            int solids = 0, openShells = 0;
            long long skipped = 0;
            bool usedRefit = false;
        };
        std::vector<CompOut> outs(order.size());
        std::atomic<size_t> nextComp{ 0 };

        auto buildComponent = [&](int root, CompOut& out, unsigned subThreads) {
            const CompStat& cs = comps.at(root);
            size_t nV = cs.vtx.size(), nE = cs.edges.size(), nT = cs.tris.size();

            // Phase 1: vertices, then edges, then faces -- independent OCCT objects
            // written into pre-assigned slots, so each phase parallelises lock-free.
            std::vector<TopoDS_Vertex> verts(nV);
            parallelFor(subThreads, nV, [&](size_t i) {
                BRep_Builder bb;
                bb.MakeVertex(verts[i], gp_Pnt(pts[cs.vtx[i]]), Precision::Confusion());
            });
            std::vector<TopoDS_Edge> edges(nE);
            std::vector<uint8_t> edgeOk(nE, 0);
            parallelFor(subThreads, nE, [&](size_t i) {
                BRepBuilderAPI_MakeEdge me(verts[cs.edges[i].first], verts[cs.edges[i].second]);
                if (me.IsDone()) { edges[i] = me.Edge(); edgeOk[i] = 1; }
            });
            std::vector<TopoDS_Face> built;
            bool usedRefit = false;
            auto pit = refitPlans.find(root);
            if (pit != refitPlans.end()) {
                refit::MeshView mv{};
                fillMeshView(cs, mv);
                refit::RegionSet rs = pit->second;
                std::vector<TopoDS_Face> rf;
                bool ok = refit::buildFaces(mv, rs, verts, rf,
                                           [&](const std::string& m) { warn(m); })
                    && !rf.empty();
                if (ok) {
                    TopoDS_Shell probe;
                    BRep_Builder pb;
                    pb.MakeShell(probe);
                    for (auto& f : rf) pb.Add(probe, f);
                    ok = BRep_Tool::IsClosed(probe);
                    if (ok) {
                        try {
                            BRepCheck_Analyzer an(probe, Standard_True, Standard_True);
                            ok = an.IsValid();
                        } catch (const Standard_Failure&) { ok = false; }
                    }
                    if (ok) {
                        usedRefit = true;
                        built = std::move(rf);
                        out.keepShapes = built;
                        out.usedRefit = true;
                        out.refitSt = rs.stats;
                        for (const auto& reg : rs.regions)
                            out.refitDVolAbs += std::fabs(reg.dVolPredicted);
                    }
                }
                if (!ok)
                    warn("smooth: analytic rebuild reverted on one component -- kept faceted");
            }
            if (!usedRefit) {
            if (pit != refitPlans.end()) {
                // R2 must be verbatim: P2 UpdateVertex's TShapes in place.
                parallelFor(subThreads, nV, [&](size_t i) {
                    BRep_Builder bb;
                    bb.MakeVertex(verts[i], gp_Pnt(pts[cs.vtx[i]]), Precision::Confusion());
                });
                parallelFor(subThreads, nE, [&](size_t i) {
                    BRepBuilderAPI_MakeEdge me(verts[cs.edges[i].first], verts[cs.edges[i].second]);
                    if (me.IsDone()) { edges[i] = me.Edge(); edgeOk[i] = 1; }
                    else edgeOk[i] = 0;
                });
            }
            std::vector<TopoDS_Face> faces(nT);
            std::vector<uint8_t> faceOk(nT, 0);
            std::atomic<long long> skipped{ 0 };
            parallelFor(subThreads, nT, [&](size_t k) {
                const Tri& T = tris[cs.tris[k]];
                const gp_XYZ &A = pts[T.a], &Bp = pts[T.b], &Cp = pts[T.c];
                gp_XYZ n = (Bp - A).Crossed(Cp - A);
                double mag = n.Modulus();
                double l2 = std::max((Bp - A).SquareModulus(), (Cp - A).SquareModulus());
                if (mag < 1e-12 || mag * mag < l2 * l2 * 1e-20) { skipped++; return; }
                BRep_Builder bb;
                TopoDS_Wire w;
                bb.MakeWire(w);
                for (int s = 0; s < 3; s++) {
                    int id = cs.triEdges[k][s];
                    if (!edgeOk[id]) { skipped++; return; }
                    bool fwd = (cs.triDirs[k] >> s) & 1;
                    bb.Add(w, fwd ? edges[id] : TopoDS::Edge(edges[id].Reversed()));
                }
                w.Closed(Standard_True);
                BRepBuilderAPI_MakeFace mf(gp_Pln(gp_Pnt(A), gp_Dir(n)), w, Standard_True);
                if (!mf.IsDone()) { skipped++; return; }
                faces[k] = mf.Face();
                faceOk[k] = 1;
            });
            out.skipped = skipped;
            built.reserve(nT);
            for (size_t k = 0; k < nT; k++)
                if (faceOk[k]) built.push_back(faces[k]);
            }

            // Phase 2: shell -> solid (serial per component; child-list mutation).
            BRep_Builder B;
            auto orientedSolid = [&](const TopoDS_Shell& sh) -> TopoDS_Shape {
                TopoDS_Solid so;
                B.MakeSolid(so);
                B.Add(so, sh);
                try {
                    BRepClass3d_SolidClassifier cl(so);
                    cl.PerformInfinitePoint(Precision::Confusion());
                    if (cl.State() == TopAbs_IN) {          // inside-out: rebuild reversed
                        TopoDS_Solid so2;
                        B.MakeSolid(so2);
                        B.Add(so2, TopoDS::Shell(sh.Reversed()));
                        so = so2;
                    }
                } catch (const Standard_Failure&) {
                    warn("solid orientation check failed on one component; kept as built");
                }
                return so;
            };
            auto finishShell = [&](TopoDS_Shell sh, bool knownClosed, const char* how) {
                bool closed = knownClosed || BRep_Tool::IsClosed(sh);
                sh.Closed(closed ? Standard_True : Standard_False);
                if (closed && makeSolids) { out.parts.push_back(orientedSolid(sh)); out.solids++; }
                else {
                    if (makeSolids) {
                        out.openShells++;
                        warn(std::string("component is not a closed volume (") + how +
                             ") -- exported as an open shell; CAM may not machine it correctly");
                    }
                    out.parts.push_back(sh);
                }
            };

            bool clean = isClean(cs);
            if (clean) {
                TopoDS_Shell sh;
                B.MakeShell(sh);
                for (auto& f : built) B.Add(sh, f);
                finishShell(sh, usedRefit ? false : true, "direct");
                return;
            }
            // Repair: sew cracks, then fix face orientation.
            Timer rt;
            BRepBuilderAPI_Sewing sew(sewTol);
            for (auto& f : built) sew.Add(f);
            sew.Perform();
            TopoDS_Shape sewed = sew.SewedShape();
            if (sewed.IsNull()) {
                warn("sewing produced nothing for one component -- component dropped");
                return;
            }
            bool anyShell = false;
            for (TopExp_Explorer ex(sewed, TopAbs_SHELL); ex.More(); ex.Next()) {
                anyShell = true;
                TopoDS_Shell sh = TopoDS::Shell(ex.Current());
                try {
                    ShapeFix_Shell fix(sh);
                    fix.Perform();
                    for (TopExp_Explorer fx(fix.Shape(), TopAbs_SHELL); fx.More(); fx.Next())
                        finishShell(TopoDS::Shell(fx.Current()), false, "repaired");
                } catch (const Standard_Failure&) {
                    finishShell(sh, false, "repaired");
                }
            }
            int freeFaces = countShapes(sewed, TopAbs_FACE, TopAbs_SHELL);
            if (freeFaces > 0) {
                TopoDS_Shell leftovers;
                B.MakeShell(leftovers);
                for (TopExp_Explorer ex(sewed, TopAbs_FACE, TopAbs_SHELL); ex.More(); ex.Next())
                    B.Add(leftovers, ex.Current());
                finishShell(leftovers, false, "unsewn leftovers");
            }
            if (!anyShell && freeFaces == 0)
                warn("sewing returned no shells or faces for one component -- component dropped");
            note("  repair    component of %s faces (open %d, conflict %d, non-manifold %d) sewn at %.2g mm  [%.2fs]\n",
                 fmtInt((long long)cs.tris.size()).c_str(),
                 cs.open, cs.conflict, cs.nonManifold, sewTol, rt.lap());
        };

        {
            std::vector<std::thread> pool;
            for (unsigned k = 0; k < nOuter; k++)
                pool.emplace_back([&]() {
                    for (;;) {
                        size_t i = nextComp.fetch_add(1);
                        if (i >= order.size()) break;
                        try {
                            buildComponent(order[i], outs[i], nInner);
                        } catch (const Standard_Failure& f) {
                            warn(std::string("component build failed (") +
                                 (f.GetMessageString() ? f.GetMessageString() : "?") +
                                 ") -- component dropped");
                        } catch (const std::exception& e) {
                            warn(std::string("component build failed (") + e.what() +
                                 ") -- component dropped");
                        }
                    }
                });
            for (auto& th : pool) th.join();
        }

        std::vector<TopoDS_Shape> parts;
        std::vector<const std::vector<TopoDS_Face>*> partKeep;
        long long skipped = 0;
        bool anyRefitUsed = false;
        for (auto& o : outs) {
            for (auto& p : o.parts) {
                parts.push_back(p);
                partKeep.push_back(o.keepShapes.empty() ? nullptr : &o.keepShapes);
            }
            solidsOut += o.solids;
            openShellsOut += o.openShells;
            skipped += o.skipped;
            if (!o.usedRefit) continue;
            anyRefitUsed = true;
            const refit::RefitStats& st = o.refitSt;
            refitTotals.planes += st.planes;
            refitTotals.cylinders += st.cylinders;
            refitTotals.fillets += st.fillets;
            refitTotals.rejected += st.rejected;
            refitTotals.facetIslands += st.facetIslands;
            refitTotals.facetTriangles += st.facetTriangles;
            refitTotals.distinctRadii += st.distinctRadii;
            refitTotals.maxVertexDev = std::max(refitTotals.maxVertexDev, st.maxVertexDev);
            refitTotals.maxEdgeTol = std::max(refitTotals.maxEdgeTol, st.maxEdgeTol);
            dVolPredSigned += st.dVolPredicted;
            dVolPredAbs += o.refitDVolAbs;
        }
        {
            long long idxDegen = (long long)nTri - (long long)tris.size();
            note("  build     %s faces (%s degenerate skipped), %d manifold component%s (%u threads)  [%.2fs]\n",
                 fmtInt((long long)tris.size() - skipped).c_str(),
                 fmtInt(idxDegen + skipped).c_str(), nComp, nComp == 1 ? "" : "s",
                 std::max(nOuter, nInner), timer.lap());
        }

        if (parts.empty())
            return fail(r, "no usable geometry");
        TopoDS_Shape shape;
        if (parts.size() == 1) shape = parts[0];
        else {
            TopoDS_Compound comp;
            BRep_Builder B;
            B.MakeCompound(comp);
            for (auto& p : parts) B.Add(comp, p);
            shape = comp;
        }
        if (makeSolids)
            note("  solid     %d solid%s, %d open shell%s\n",
                 solidsOut, solidsOut == 1 ? "" : "s",
                 openShellsOut, openShellsOut == 1 ? "" : "s");

        // ---- 6. unify coplanar facets (parallel across bodies) ------------------
        facesBefore = countShapes(shape, TopAbs_FACE);
        facesAfter = facesBefore;
        if (unify) {
            auto unifyOne = [&](TopoDS_Shape& s, const std::vector<TopoDS_Face>* keep) {
                ShapeUpgrade_UnifySameDomain usd(s, Standard_True, Standard_True, Standard_False);
                if (unifyAngleDeg > 0) {
                    usd.SetAngularTolerance(unifyAngleDeg * kPi / 180.0);
                    usd.SetLinearTolerance(Precision::Confusion());
                }
                if (keep)
                    for (const auto& f : *keep) usd.KeepShape(f);
                usd.Build();
                s = usd.Shape();
            };
            try {
                if (parts.size() > 1) {
                    unsigned uThreads = (unsigned)std::min<size_t>(
                        { (size_t)hw, parts.size(), (size_t)10 });
                    std::atomic<size_t> nextPart{ 0 };
                    std::atomic<bool> anyFail{ false };
                    std::vector<std::thread> pool;
                    for (unsigned k = 0; k < uThreads; k++)
                        pool.emplace_back([&]() {
                            for (;;) {
                                size_t i = nextPart.fetch_add(1);
                                if (i >= parts.size()) break;
                                try { unifyOne(parts[i], partKeep[i]); }
                                catch (...) { anyFail = true; }   // that body keeps its facets
                            }
                        });
                    for (auto& th : pool) th.join();
                    TopoDS_Compound comp;
                    BRep_Builder B;
                    B.MakeCompound(comp);
                    for (auto& p : parts) B.Add(comp, p);
                    shape = comp;
                    if (anyFail)
                        warn("coplanar merge failed on at least one body -- it keeps per-triangle faces");
                    facesAfter = countShapes(shape, TopAbs_FACE);
                    note("  unify     %s -> %s faces (%u threads)  [%.2fs]\n",
                         fmtInt(facesBefore).c_str(), fmtInt(facesAfter).c_str(),
                         uThreads, timer.lap());
                } else {
                    unifyOne(shape, partKeep.empty() ? nullptr : partKeep[0]);
                    parts[0] = shape;
                    facesAfter = countShapes(shape, TopAbs_FACE);
                    note("  unify     %s -> %s faces  [%.2fs]\n",
                         fmtInt(facesBefore).c_str(), fmtInt(facesAfter).c_str(), timer.lap());
                }
            } catch (const Standard_Failure& f) {
                warn(std::string("coplanar merge failed (") +
                     (f.GetMessageString() ? f.GetMessageString() : "?") +
                     ") -- keeping per-triangle faces");
            }
        }

        int builtPl = 0, builtCy = 0, builtFi = 0, builtCo = 0, revCo = 0;
        if (smooth) {
            size_t pi = 0;
            for (size_t ci = 0; ci < order.size(); ++ci) {
                const CompOut& o = outs[ci];
                int cyls = 0;
                for (size_t k = 0; k < o.parts.size(); ++k) {
                    for (TopExp_Explorer ex(parts[pi], TopAbs_FACE); ex.More(); ex.Next())
                        if (BRepAdaptor_Surface(TopoDS::Face(ex.Current()), Standard_False)
                                .GetType() == GeomAbs_Cylinder)
                            cyls++;
                    ++pi;
                }
                if (!refitPlans.count(order[ci])) continue;
                const refit::RefitStats& st = o.refitSt;
                if (o.usedRefit && (cyls > 0 || (st.cylinders == 0 && st.planes > 0))) {
                    builtCo++;
                    builtPl += st.planes;
                    builtCy += cyls;
                    builtFi += st.fillets;
                } else
                    revCo++;
            }
        }

        // ---- 7. tolerance fit, then validity + volume OVERLAPPED with write ----
        fitPlanarTolerances(shape);
        const int refitFaces = refitTotals.planes + refitTotals.cylinders + refitTotals.fillets;
        if (smooth && anyRefitUsed && refitTotals.cylinders > 0)
            refit::fitAnalyticTolerances(shape);
        if (smooth && anyRefitUsed && refitFaces > 0)
            BRepLib::EncodeRegularity(shape, Precision::Angular());
        note("  fit       tolerances fitted to true planar deviation  [%.2fs]\n", timer.lap());

        bool doCheck = facesAfter <= 50000;
        bool doVol = makeSolids && solidsOut > 0 && openShellsOut == 0;
        struct Analysis {
            bool checkRan = false, valid = true, failed = false, volRan = false;
            double vol = 0, checkS = 0, volS = 0;
        } ana;
        // Both are pure readers of the finished shape (BRepCheck's own parallel
        // mode already reads one shape from many threads), so they can run while
        // the STEP writer serialises the same shape.
        auto runAnalysis = [&]() {
            try {
                if (doCheck) {
                    Timer t;
                    BRepCheck_Analyzer a(shape, Standard_True, Standard_True);
                    ana.valid = a.IsValid();
                    ana.checkRan = true;
                    ana.checkS = t.lap();
                }
                if (doVol) {
                    Timer t;
                    ana.vol = shapeVolume(shape);
                    ana.volRan = true;
                    ana.volS = t.lap();
                }
            } catch (const Standard_Failure&) { ana.failed = true; }
              catch (const std::exception&)   { ana.failed = true; }
        };

        Interface_Static::SetCVal("write.step.schema", schema.c_str());
        Interface_Static::SetCVal("write.step.unit", "MM");
        Interface_Static::SetCVal("write.step.product.name", stem.c_str());
        // 0 = ok, 1 = transfer failed, 2 = write failed
        auto writeStep = [&]() -> int {
            STEPControl_Writer writer;
            if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone) return 1;
            return writer.Write(outPath.c_str()) == IFSelect_RetDone ? 0 : 2;
        };

        int wst;
        double writeS;
        bool overlapped = hw > 1 && (doCheck || doVol);
        {
            Timer wt;
            if (overlapped) {
                std::thread bg(runAnalysis);
                wst = writeStep();
                bg.join();
            } else {
                runAnalysis();
                wst = writeStep();
            }
            writeS = wt.lap();
        }
        if (ana.failed) warn("validity/volume analysis failed with an OCCT exception");
        if (doCheck && ana.checkRan)
            note("  check     %s%s  [%.2fs]\n", ana.valid ? "valid" : "INVALID",
                 overlapped ? " (ran during write)" : "", ana.checkS);
        else if (!doCheck)
            note("  check     skipped (%s faces)\n", fmtInt(facesAfter).c_str());

        if (doCheck && ana.checkRan && !ana.valid && smooth && anyRefitUsed)
            warn("smooth: B-Rep invalid after build -- no ShapeFix rewrite on smooth runs");
        if (doCheck && ana.checkRan && !ana.valid && !(smooth && anyRefitUsed)) {
            // Rare path: repair, then rewrite so the file on disk is the fixed shape.
            try {
                Timer ft;
                ShapeFix_Shape fix(shape);
                fix.Perform();
                shape = fix.Shape();
                BRepCheck_Analyzer a2(shape, Standard_True, Standard_True);
                bool ok2 = a2.IsValid();
                note("  fix       %s after ShapeFix -- rewriting output  [%.2fs]\n",
                     ok2 ? "valid" : "STILL INVALID", ft.lap());
                if (!ok2) warn("shape still reports invalid after ShapeFix -- STEP written anyway");
                Timer rw;
                wst = writeStep();
                if (doVol) ana.vol = shapeVolume(shape);
                writeS += rw.lap();
            } catch (const Standard_Failure&) {
                warn("ShapeFix failed -- STEP left as originally written");
            }
        }
        if (wst == 1) return fail(r, "STEP translation failed");
        if (wst == 2) return fail(r, "could not write output");

        if (doVol && ana.volRan) {
            brepVolume = ana.vol;
            double refVol = watertight ? meshVolume : brepVolume;
            if (watertight && refVol > 0) {
                if (smooth && dVolPredAbs > 0) {
                    double budget = std::max(1e-4 * refVol, 3.0 * dVolPredAbs);
                    if (std::fabs(brepVolume - meshVolume) > budget)
                        warn("B-Rep volume differs from mesh volume beyond refit budget -- inspect the result");
                } else if (std::fabs(brepVolume - meshVolume) / refVol > 1e-4) {
                    warn("B-Rep volume differs from mesh volume by more than 0.01% -- inspect the result");
                }
            }
            note("  volume    B-Rep %.3f mm^3, mesh %.3f mm^3%s  [%.2fs]\n", brepVolume, meshVolume,
                 watertight ? "" : " (mesh not watertight; mesh figure approximate)", ana.volS);
        }
        long long outSize = fileSize(outPath);
        note("  write     %s, %s (%.1f KB)  [%.2fs]\n", schema.c_str(), outPath.c_str(),
             outSize / 1024.0, writeS);
        timer.lap();                                    // reset stage timer for verify

        // ---- 8. verify by re-reading -------------------------------------------
        if (verify) {
            STEPControl_Reader reader;
            if (reader.ReadFile(outPath.c_str()) != IFSelect_RetDone) {
                warn("verification failed: written STEP could not be re-read");
            } else {
                reader.TransferRoots();
                TopoDS_Shape back = reader.OneShape();
                int rSolids = countShapes(back, TopAbs_SOLID);
                int rFaces = countShapes(back, TopAbs_FACE);
                if (rSolids != solidsOut)
                    warn("verification: re-read solid count " + std::to_string(rSolids) +
                         " != written " + std::to_string(solidsOut));
                if (solidsOut > 0 && openShellsOut == 0) {
                    stepVolume = shapeVolume(back);
                    double ref = brepVolume != 0 ? brepVolume : 1;
                    volDeltaPct = std::fabs(stepVolume - brepVolume) / std::fabs(ref) * 100.0;
                    if (volDeltaPct > 0.1)
                        warn("verification: STEP volume deviates " +
                             std::to_string(volDeltaPct) + "% from source");
                }
                note("  verify    re-read %d solid%s, %s faces%s  [%.2fs]\n",
                     rSolids, rSolids == 1 ? "" : "s", fmtInt(rFaces).c_str(),
                     volDeltaPct >= 0
                         ? (", volume delta " + std::to_string(volDeltaPct) + "%").c_str()
                         : "",
                     timer.lap());
            }
        }

        // ---- assemble the result -----------------------------------------------
        r.seconds = total.lap();
        r.triangles = nTri;
        r.vertices = nVert;
        r.components = nComp;
        r.solids = solidsOut;
        r.openShells = openShellsOut;
        r.facesBeforeUnify = facesBefore;
        r.facesAfterUnify = facesAfter;
        r.meshVolumeMM3 = meshVolume;
        r.stepVolumeMM3 = stepVolume;
        r.volumeDeltaPct = volDeltaPct;
        r.watertight = watertight;
        r.warnings = warnings;
        if (smooth) {
            r.smoothPlanes = refitTotals.planes;
            r.smoothCylinders = refitTotals.cylinders;
            r.smoothFillets = refitTotals.fillets;
            r.smoothDistinctRadii = refitTotals.distinctRadii;
            r.smoothRejected = refitTotals.rejected;
            r.smoothFacetFaces = refitTotals.facetTriangles;
            r.facesAfterSmooth = facesAfter;
            r.smoothSkippedComponents = smoothSkippedComponents;
            r.smoothMaxDevMM = refitTotals.maxVertexDev;
            r.smoothMaxEdgeTolMM = refitTotals.maxEdgeTol;
            r.smoothVolPredictedMM3 = dVolPredSigned;
            r.smoothBuiltPlanes = builtPl;
            r.smoothBuiltCylinders = builtCy;
            r.smoothBuiltFillets = builtFi;
            r.smoothBuiltComponents = builtCo;
            r.smoothRevertedComponents = revCo;
        }
        r.ok = true;
        r.exitCode = warnings.empty() ? 0 : 2;
        return r;

    } catch (const Standard_Failure& f) {
        return fail(r, std::string("OCCT failure: ") +
                    (f.GetMessageString() ? f.GetMessageString() : "?"));
    } catch (const std::exception& e) {
        return fail(r, e.what());
    }
}

}  // namespace

// ---------------------------------------------------------------- public API

const char* version() { return STL2STEP_VERSION_STRING; }

const char* schemaName(Schema s) {
    switch (s) {
        case Schema::AP203: return "AP203";
        case Schema::AP242: return "AP242";
        case Schema::AP214: default: return "AP214";
    }
}

bool parseSchema(const std::string& text, Schema& out) {
    std::string u;
    for (char c : text) u += (char)std::toupper((unsigned char)c);
    if (u == "AP203") { out = Schema::AP203; return true; }
    if (u == "AP214") { out = Schema::AP214; return true; }
    if (u == "AP242") { out = Schema::AP242; return true; }
    return false;
}

// toJson() splices smooth* keys iff convert() filled them (Options::smooth).
// Frozen header has no Result::smoothEnabled — splice is member-keyed, not
// thread_local. Hand-filled Results must set facesAfterSmooth to emit keys.
std::string Result::toJson() const {
    if (!ok)
        return std::string("{\"ok\":false,\"error\":\"") + jsonEscape(error) + "\"}";

    std::string wjson;
    for (auto& w : warnings) {
        if (!wjson.empty()) wjson += ",";
        wjson += "\"" + jsonEscape(w) + "\"";
    }
    const char* fmt =
        "{\"ok\":true,\"input\":\"%s\",\"output\":\"%s\",\"triangles\":%d,"
        "\"vertices\":%d,\"components\":%d,\"solids\":%d,\"openShells\":%d,"
        "\"facesBeforeUnify\":%d,\"facesAfterUnify\":%d,\"meshVolumeMM3\":%.6f,"
        "\"stepVolumeMM3\":%.6f,\"volumeDeltaPct\":%.6f,\"watertight\":%s,"
        "\"seconds\":%.2f,\"warnings\":[%s]}";
    const char* fmtSmooth =
        "{\"ok\":true,\"input\":\"%s\",\"output\":\"%s\",\"triangles\":%d,"
        "\"vertices\":%d,\"components\":%d,\"solids\":%d,\"openShells\":%d,"
        "\"facesBeforeUnify\":%d,\"facesAfterUnify\":%d,\"meshVolumeMM3\":%.6f,"
        "\"stepVolumeMM3\":%.6f,\"volumeDeltaPct\":%.6f,\"watertight\":%s,"
        "\"seconds\":%.2f,\"warnings\":[%s],"
        "\"smoothPlanes\":%d,\"smoothCylinders\":%d,\"smoothFillets\":%d,"
        "\"smoothDistinctRadii\":%d,\"smoothRejected\":%d,\"smoothFacetFaces\":%d,"
        "\"facesAfterSmooth\":%d,\"smoothSkippedComponents\":%d,"
        "\"smoothMaxDevMM\":%.6f,\"smoothMaxEdgeTolMM\":%.6f,"
        "\"smoothVolPredictedMM3\":%.6f,"
        "\"smoothBuiltPlanes\":%d,\"smoothBuiltCylinders\":%d,\"smoothBuiltFillets\":%d,"
        "\"smoothBuiltComponents\":%d,\"smoothRevertedComponents\":%d}";
    std::string ei = jsonEscape(input), eo = jsonEscape(output);
    const bool emitSmooth = facesAfterSmooth != 0 || smoothSkippedComponents != 0
        || smoothPlanes != 0 || smoothCylinders != 0 || smoothFillets != 0;
    int n;
    if (emitSmooth) {
        n = std::snprintf(nullptr, 0, fmtSmooth, ei.c_str(), eo.c_str(), triangles, vertices,
                          components, solids, openShells, facesBeforeUnify, facesAfterUnify,
                          meshVolumeMM3, stepVolumeMM3, volumeDeltaPct,
                          watertight ? "true" : "false", seconds, wjson.c_str(),
                          smoothPlanes, smoothCylinders, smoothFillets, smoothDistinctRadii,
                          smoothRejected, smoothFacetFaces, facesAfterSmooth,
                          smoothSkippedComponents, smoothMaxDevMM, smoothMaxEdgeTolMM,
                          smoothVolPredictedMM3, smoothBuiltPlanes, smoothBuiltCylinders,
                          smoothBuiltFillets, smoothBuiltComponents, smoothRevertedComponents);
    } else {
        n = std::snprintf(nullptr, 0, fmt, ei.c_str(), eo.c_str(), triangles, vertices,
                          components, solids, openShells, facesBeforeUnify, facesAfterUnify,
                          meshVolumeMM3, stepVolumeMM3, volumeDeltaPct,
                          watertight ? "true" : "false", seconds, wjson.c_str());
    }
    std::string s((size_t)n + 1, '\0');
    if (emitSmooth) {
        std::snprintf(&s[0], s.size(), fmtSmooth, ei.c_str(), eo.c_str(), triangles, vertices,
                      components, solids, openShells, facesBeforeUnify, facesAfterUnify,
                      meshVolumeMM3, stepVolumeMM3, volumeDeltaPct,
                      watertight ? "true" : "false", seconds, wjson.c_str(),
                      smoothPlanes, smoothCylinders, smoothFillets, smoothDistinctRadii,
                      smoothRejected, smoothFacetFaces, facesAfterSmooth,
                      smoothSkippedComponents, smoothMaxDevMM, smoothMaxEdgeTolMM,
                      smoothVolPredictedMM3, smoothBuiltPlanes, smoothBuiltCylinders,
                      smoothBuiltFillets, smoothBuiltComponents, smoothRevertedComponents);
    } else {
        std::snprintf(&s[0], s.size(), fmt, ei.c_str(), eo.c_str(), triangles, vertices,
                      components, solids, openShells, facesBeforeUnify, facesAfterUnify,
                      meshVolumeMM3, stepVolumeMM3, volumeDeltaPct,
                      watertight ? "true" : "false", seconds, wjson.c_str());
    }
    s.resize((size_t)n);
    return s;
}

Result convert(const Options& opt, const LogCallback& log) {
    Converter c(opt, log);
    return c.run();
}

}  // namespace stl2step
