// Test-side port of Converter::run() stages 2–3 from src/stl2step.cpp.
// Loops are transcribed in source order; do not "improve" them — local ids
// are the I5 sort keys.
//
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <unordered_map>

#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <OSD.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <RWStl.hxx>
#include <Standard_Failure.hxx>

namespace stl2step {
namespace harness {
namespace {

// Exact-bit vertex key (normalises -0.0 -> +0.0); tolerance welding quantises instead.
// Ported from src/stl2step.cpp:150-166.
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

}  // namespace

bool loadMesh(const std::string& stlPath, double scale, double weldTolArg,
              double sewTolArg, HarnessMesh& out, std::string& err) {
    out = HarnessMesh{};
    err.clear();
    try {
        Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
        OSD::SetSignal(Standard_False);

        Handle(Poly_Triangulation) mesh = RWStl::ReadFile(stlPath.c_str());
        if (mesh.IsNull() || mesh->NbTriangles() < 1) {
            err = "unreadable or empty STL";
            return false;
        }
        int rawVerts = mesh->NbNodes();
        int nTri = mesh->NbTriangles();
        double weldTol = weldTolArg;

        // ---- 2. weld vertices (exact by default, tolerance grid via weldTol) ---
        // Ported from src/stl2step.cpp:312-338.
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

        // Ported from src/stl2step.cpp:340-348.
        gp_XYZ bbMin(1e300, 1e300, 1e300), bbMax(-1e300, -1e300, -1e300);
        for (auto& p : pts) {
            bbMin.SetX(std::min(bbMin.X(), p.X())); bbMax.SetX(std::max(bbMax.X(), p.X()));
            bbMin.SetY(std::min(bbMin.Y(), p.Y())); bbMax.SetY(std::max(bbMax.Y(), p.Y()));
            bbMin.SetZ(std::min(bbMin.Z(), p.Z())); bbMax.SetZ(std::max(bbMax.Z(), p.Z()));
        }
        double diag = (bbMax - bbMin).Modulus();
        double sewTol = sewTolArg > 0 ? sewTolArg
                                      : std::min(std::max(1e-6, diag * 1e-5), 0.5);

        // ---- 3. edge stats + manifold component split ---------------------------
        // Ported from src/stl2step.cpp:357-368.
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

        // Ported from src/stl2step.cpp:370-373.
        auto ekey = [](int a, int b) {
            return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                         : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
        };

        // Pass A: global edge use counts. Pass B: union triangles ONLY across edges
        // used exactly twice. Ported from src/stl2step.cpp:379-395.
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

        // Pass C: per component tables. Ported from src/stl2step.cpp:401-466.
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

        out.pts = std::move(pts);
        out.tris.reserve(tris.size());
        for (const Tri& t : tris) out.tris.push_back({ t.a, t.b, t.c });
        out.diag = diag;
        out.weldTol = weldTol;
        out.sewTol = sewTol;

        // Gate-tool order: ascending UnionFind root id. The engine then sorts
        // by triangle-count desc for build scheduling; that is not the contract.
        std::vector<int> roots;
        roots.reserve(comps.size());
        for (auto& kv : comps) roots.push_back(kv.first);
        std::sort(roots.begin(), roots.end());
        out.comps.reserve(roots.size());
        for (int root : roots) {
            CompStat& c = comps.at(root);
            HarnessComponent hc;
            hc.root        = root;
            hc.compTris    = std::move(c.tris);
            hc.compVtx     = std::move(c.vtx);
            hc.compEdges   = std::move(c.edges);
            hc.triEdges    = std::move(c.triEdges);
            hc.triDirs     = std::move(c.triDirs);
            hc.open        = c.open;
            hc.conflict    = c.conflict;
            hc.nonManifold = c.nonManifold;
            hc.vol         = c.vol;
            out.comps.push_back(std::move(hc));
        }
        return true;
    } catch (const Standard_Failure& e) {
        out = HarnessMesh{};
        const char* msg = e.GetMessageString();
        err = (msg && msg[0]) ? msg : "OpenCASCADE failure reading STL";
        return false;
    } catch (const std::exception& e) {
        out = HarnessMesh{};
        err = e.what();
        return false;
    } catch (...) {
        out = HarnessMesh{};
        err = "unknown failure reading STL";
        return false;
    }
}

}  // namespace harness
}  // namespace stl2step
