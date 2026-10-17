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
#include "mesh_heal.hpp"
#include "parallel.hpp"
#include "refit.hpp"
#include "stl_quant.hpp"

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
#include <limits>
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
void prismBindDxfDir(const std::string& d);
void prismBindDxfStem(const std::string& s);
}

namespace detail {

// The faceted builder's existing zero-area test (stl2step.cpp face loop and
// makeFacet in refit_build.cpp). Same expression, both call sites. Not a new gate.
bool triangleAreaDegenerate(const gp_XYZ& A, const gp_XYZ& B, const gp_XYZ& C) {
    const gp_XYZ n = (B - A).Crossed(C - A);
    const double mag = n.Modulus();
    const double l2 = std::max((B - A).SquareModulus(), (C - A).SquareModulus());
    return mag < 1e-12 || mag * mag < l2 * l2 * 1e-20;
}

// Absorb a zero-area T-junction before the edge graph is built.
// A triangle is removed only when it fails triangleAreaDegenerate, one vertex
// lies within tau = 2·quantFloor of the open segment of the opposite edge
// (t strictly inside (0,1)), and that edge is shared by exactly one other
// triangle. The neighbour is split at the vertex, keeping its winding.
// Net triangle count is unchanged. Stable by ascending triangle index.
// quantFloor <= 0 absorbs nothing. Returns the number absorbed.
int absorbZeroAreaTJunctions(std::vector<std::array<int, 3>>& tris,
                             const std::vector<std::array<double, 3>>& pts,
                             double quantFloor, bool diag) {
    const double tau = 2.0 * quantFloor; // D-130-12. The only length threshold.
    auto summary = [&](int n) {
        if (!diag) return;
        std::fprintf(stderr, "DIAG_ABSORB summary n=%d trisNow=%zu\n", n, tris.size());
    };
    if (!(quantFloor > 0.0)) {
        summary(0);
        return 0;
    }

    auto ekey = [](int a, int b) -> uint64_t {
        return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                     : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
    };
    struct EdgeUse { int cnt = 0; int t0 = -1; int t1 = -1; };
    std::unordered_map<uint64_t, EdgeUse> edges;
    auto buildEdges = [&]() {
        edges.clear();
        edges.reserve(tris.size() * 2);
        for (int t = 0; t < (int)tris.size(); ++t) {
            const int v[4] = { tris[t][0], tris[t][1], tris[t][2], tris[t][0] };
            for (int s = 0; s < 3; ++s) {
                EdgeUse& e = edges[ekey(v[s], v[s + 1])];
                if (e.cnt == 0) e.t0 = t;
                else if (e.cnt == 1) e.t1 = t;
                e.cnt++;
            }
        }
    };
    auto at = [&](int i) -> gp_XYZ {
        const std::array<double, 3>& p = pts[(size_t)i];
        return gp_XYZ(p[0], p[1], p[2]);
    };
    auto inRange = [&](int i) {
        return i >= 0 && (size_t)i < pts.size();
    };

    buildEdges();
    int absorbed = 0;
    for (int t = 0; t < (int)tris.size(); ++t) {
        const int v[3] = { tris[t][0], tris[t][1], tris[t][2] };
        if (!inRange(v[0]) || !inRange(v[1]) || !inRange(v[2])) continue;
        if (!triangleAreaDegenerate(at(v[0]), at(v[1]), at(v[2]))) continue;

        bool did = false;
        for (int slot = 0; slot < 3 && !did; ++slot) {
            const int vertex = v[slot];
            const int ea = v[(slot + 1) % 3];
            const int eb = v[(slot + 2) % 3];
            const gp_XYZ P = at(vertex), A = at(ea), B = at(eb);
            const gp_XYZ ab = B - A;
            const double ab2 = ab.SquareModulus();
            if (!(ab2 > 0.0)) continue;
            const double param = (P - A).Dot(ab) / ab2;
            if (!(param > 0.0 && param < 1.0)) continue;
            const gp_XYZ closest(A.X() + ab.X() * param,
                                 A.Y() + ab.Y() * param,
                                 A.Z() + ab.Z() * param);
            const double dist = (P - closest).Modulus();
            if (!(dist <= tau)) continue;

            auto it = edges.find(ekey(ea, eb));
            if (it == edges.end() || it->second.cnt != 2) continue;
            const int other = (it->second.t0 == t) ? it->second.t1 : it->second.t0;
            if (other < 0 || other == t || other >= (int)tris.size()) continue;

            const int nv[3] = { tris[other][0], tris[other][1], tris[other][2] };
            int sFound = -1;
            for (int s = 0; s < 3; ++s) {
                const int x = nv[s], y = nv[(s + 1) % 3];
                if ((x == ea && y == eb) || (x == eb && y == ea)) { sFound = s; break; }
            }
            if (sFound < 0) continue;
            const int vi = nv[sFound];
            const int vj = nv[(sFound + 1) % 3];
            const int apex = nv[(sFound + 2) % 3];
            if (vertex == vi || vertex == vj || vertex == apex) continue;

            // Neighbour winding vi → vj → apex. Insert vertex on vi→vj.
            tris[other] = { vi, vertex, apex };
            tris[t] = { vertex, vj, apex };
            absorbed++;
            did = true;
            if (diag) {
                std::fprintf(stderr,
                    "DIAG_ABSORB degenTri=%d neigh=%d vertex=%d distMM=%.9g "
                    "distTau=%.4f t=%.6f tau=%.9g q=%.9g\n",
                    t, other, vertex, dist, dist / tau, param, tau, quantFloor);
            }
        }
        if (did) buildEdges();
    }
    summary(absorbed);
    return absorbed;
}

// Sorted uniform-grid index over a vertex subset (cell edge h) for the heal's
// neighbour queries. It only bounds which vertices a query looks at; every
// query keeps its exact distance test, so no heal result depends on h. An
// x-sorted sweep is quadratic when many boundary vertices share an x slab (an
// unwelded mesh with a face normal to x); the grid stays near-linear.
class HealCellIndex {
public:
    using Cell = std::array<int64_t, 3>;
    // Indexes positions 0..verts.size()-1 of verts.
    HealCellIndex(const std::vector<std::array<double, 3>>& pts, const std::vector<int>& verts,
                  double h)
        : h_(h) {
        entries_.reserve(verts.size());
        for (int i = 0; i < (int)verts.size(); ++i)
            entries_.push_back({ cellOf(pts[(size_t)verts[(size_t)i]]), i });
        std::sort(entries_.begin(), entries_.end());
    }
    Cell cellOf(const std::array<double, 3>& p) const {
        Cell c;
        for (int k = 0; k < 3; ++k) c[k] = (int64_t)std::floor(p[(size_t)k] / h_);
        return c;
    }
    // f(position) for every indexed position in cell c, ascending position.
    template <class F> void visit(const Cell& c, F&& f) const {
        auto it = std::lower_bound(entries_.begin(), entries_.end(),
                                   std::make_pair(c, std::numeric_limits<int>::min()));
        for (; it != entries_.end() && it->first == c; ++it) f(it->second);
    }

private:
    double h_;
    std::vector<std::pair<Cell, int>> entries_;
};

// Stitch open boundaries before the edge graph (D-fc-5.10). The caller runs it
// on the TrueForm path only, and only when Pass A sees an edge used once.
//   tau       = 2·quantFloor (D-130-12)
//   sMin      = shortest edge longer than tau (an edge <= tau is itself a
//               sub-quantum duplicate)
//   stitchTol = min(sewTol, sMin / 2): the repair path's own sewing tolerance,
//               capped so no edge longer than tau can collapse
// (1) Weld boundary vertices within stitchTol (union-find). A cluster wider
//     than stitchTol is refused whole. Representative = lowest vertex id, one
//     of the file's own points. Interior vertices never move. A triangle that
//     collapses is dropped.
// (2) Split: a boundary vertex within tau of a boundary edge's open segment
//     splits the owning triangle, keeping its winding. The vertex must lie
//     farther than max(stitchTol, tau) from both ends (closer is a duplicate of
//     that end), so every split edge is strictly shorter than its parent and
//     the fixed point terminates. Order: owning triangle, then edge.
// Certificate per stitched cluster (input triangles joined by a shared vertex,
// a weld or a split): the healed cluster has no edge used more than twice, no
// winding conflict, and strictly fewer open edges. A cluster that fails keeps
// its input triangles. No capping: a loop with no partner stays open.
HealStats healOpenBoundaries(std::vector<std::array<int, 3>>& tris,
                             const std::vector<std::array<double, 3>>& pts,
                             double quantFloor, double sewTol, bool diag) {
    HealStats st;
    const double tau = 2.0 * quantFloor; // D-130-12
    if (!(quantFloor > 0.0) || tris.empty()) return st;

    auto ekey = [](int a, int b) -> uint64_t {
        return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                     : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
    };
    auto keyLo = [](uint64_t k) { return (int)(k >> 32); };
    auto keyHi = [](uint64_t k) { return (int)(k & 0xffffffffu); };
    auto dist = [&](int a, int b) {
        const double dx = pts[(size_t)a][0] - pts[(size_t)b][0];
        const double dy = pts[(size_t)a][1] - pts[(size_t)b][1];
        const double dz = pts[(size_t)a][2] - pts[(size_t)b][2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    struct EdgeUse { int cnt = 0, fwd = 0, t0 = -1; };
    std::vector<std::array<int, 3>> work = tris;
    std::vector<char> dead(work.size(), 0);
    auto census = [&](std::unordered_map<uint64_t, EdgeUse>& E) {
        E.clear();
        E.reserve(work.size() * 2);
        for (int t = 0; t < (int)work.size(); ++t) {
            if (dead[(size_t)t]) continue;
            for (int s = 0; s < 3; ++s) {
                const int a = work[(size_t)t][s], b = work[(size_t)t][(s + 1) % 3];
                EdgeUse& e = E[ekey(a, b)];
                if (e.cnt == 0) e.t0 = t;
                e.cnt++;
                if (a < b) e.fwd++;
            }
        }
    };
    std::unordered_map<uint64_t, EdgeUse> E0;
    census(E0);
    for (const auto& kv : E0)
        if (kv.second.cnt == 1) st.boundaryBefore++;
    if (st.boundaryBefore == 0) return st;

    st.sMin = std::numeric_limits<double>::infinity();
    for (const auto& kv : E0) {
        const double d = dist(keyLo(kv.first), keyHi(kv.first));
        if (d > tau) st.sMin = std::min(st.sMin, d);
    }
    if (!std::isfinite(st.sMin)) return st;
    st.stitchTol = std::min(sewTol, 0.5 * st.sMin);
    const double tol = st.stitchTol;
    const double endGap = std::max(tol, tau);

    // Stitched clusters: vertex union-find over the input triangles, joined
    // again by every weld and split. Root = lowest vertex id.
    std::vector<int> cl(pts.size());
    for (size_t i = 0; i < cl.size(); ++i) cl[i] = (int)i;
    auto clFind = [&](int x) {
        while (cl[(size_t)x] != x) { cl[(size_t)x] = cl[(size_t)cl[(size_t)x]]; x = cl[(size_t)x]; }
        return x;
    };
    auto clJoin = [&](int a, int b) {
        a = clFind(a);
        b = clFind(b);
        if (a != b) cl[(size_t)std::max(a, b)] = std::min(a, b);
    };
    for (const auto& t : tris) { clJoin(t[0], t[1]); clJoin(t[1], t[2]); }

    // ---- (1) weld boundary vertices ------------------------------------------
    std::vector<int> bv;
    for (const auto& kv : E0)
        if (kv.second.cnt == 1) { bv.push_back(keyLo(kv.first)); bv.push_back(keyHi(kv.first)); }
    std::sort(bv.begin(), bv.end());
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());
    // Union-find over positions in bv (ascending id, so the root is the lowest id).
    std::vector<int> wp(bv.size());
    for (size_t i = 0; i < wp.size(); ++i) wp[i] = (int)i;
    auto wpFind = [&](int x) {
        while (wp[(size_t)x] != x) { wp[(size_t)x] = wp[(size_t)wp[(size_t)x]]; x = wp[(size_t)x]; }
        return x;
    };
    // Every pair within tol: with cell edge tol such a pair lies in the same or
    // an adjacent cell. The partition (components of the within-tol graph) and
    // its lowest-position roots do not depend on the order pairs are found.
    {
        const HealCellIndex grid(pts, bv, tol);
        for (int i = 0; i < (int)bv.size(); ++i) {
            const int vi = bv[(size_t)i];
            const HealCellIndex::Cell c = grid.cellOf(pts[(size_t)vi]);
            for (int64_t dx = -1; dx <= 1; ++dx)
                for (int64_t dy = -1; dy <= 1; ++dy)
                    for (int64_t dz = -1; dz <= 1; ++dz)
                        grid.visit({ c[0] + dx, c[1] + dy, c[2] + dz }, [&](int j) {
                            if (j <= i || dist(vi, bv[(size_t)j]) > tol) return;
                            const int a = wpFind(i), b = wpFind(j);
                            if (a != b) wp[(size_t)std::max(a, b)] = std::min(a, b);
                        });
        }
    }
    std::vector<std::vector<int>> members(bv.size());
    for (size_t i = 0; i < bv.size(); ++i) members[(size_t)wpFind((int)i)].push_back((int)i);
    struct OpRec { int vertex; double d; };
    std::vector<OpRec> weldRec, splitRec;
    std::vector<int> repOf(pts.size(), -1);
    for (size_t root = 0; root < members.size(); ++root) {
        const std::vector<int>& c = members[root];
        if (c.size() < 2) continue;
        double diam = 0;
        for (size_t i = 0; i < c.size(); ++i)
            for (size_t j = i + 1; j < c.size(); ++j)
                diam = std::max(diam, dist(bv[(size_t)c[i]], bv[(size_t)c[j]]));
        if (diam > tol) { st.refusedWeldClusters++; continue; }
        const int rep = bv[root];
        for (int p : c) {
            const int v = bv[(size_t)p];
            if (v == rep) continue;
            repOf[(size_t)v] = rep;
            weldRec.push_back({ v, dist(v, rep) });
            clJoin(v, rep);
        }
    }
    for (size_t t = 0; t < work.size(); ++t) {
        for (int& x : work[t])
            if (repOf[(size_t)x] >= 0) x = repOf[(size_t)x];
        if (work[t][0] == work[t][1] || work[t][1] == work[t][2] || work[t][0] == work[t][2])
            dead[t] = 1;
    }

    // ---- (2) split boundary edges at boundary vertices on their open segment ----
    std::unordered_map<uint64_t, EdgeUse> E;
    for (;;) {
        census(E);
        std::vector<int> bvx;
        std::vector<std::pair<int, int>> bes;   // directed boundary edge a->b in its triangle
        std::vector<int> beT;
        for (const auto& kv : E) {
            if (kv.second.cnt != 1) continue;
            const int t = kv.second.t0;
            const int lo = keyLo(kv.first), hi = keyHi(kv.first);
            for (int s = 0; s < 3; ++s) {
                const int a = work[(size_t)t][s], b = work[(size_t)t][(s + 1) % 3];
                if ((a == lo && b == hi) || (a == hi && b == lo)) { bes.push_back({ a, b }); beT.push_back(t); }
            }
            bvx.push_back(lo);
            bvx.push_back(hi);
        }
        if (bes.empty()) break;
        std::sort(bvx.begin(), bvx.end());
        bvx.erase(std::unique(bvx.begin(), bvx.end()), bvx.end());
        std::sort(bvx.begin(), bvx.end(), [&](int a, int b) {
            return pts[(size_t)a][0] != pts[(size_t)b][0] ? pts[(size_t)a][0] < pts[(size_t)b][0] : a < b;
        });
        std::vector<size_t> ord(bes.size());
        for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) {
            return beT[x] != beT[y] ? beT[x] < beT[y] : bes[x] < bes[y];
        });
        // A vertex that can split an edge lies within tau of its open segment, so
        // two candidate sources each hold every such vertex: the x window
        // [min x - tau, max x + tau] of bvx (sorted by x, id), and the cells of a
        // grid over bvx positions within reach of samples along the edge (cell
        // edge = this round's median boundary edge length, never below endGap so
        // it is positive and >= tau; samples sampleGap apart leave every point of
        // the edge within sampleGap / 2 of one). Each edge scans the smaller
        // source, in bvx order with a strict minimum, so the winner is the same
        // either way; the window alone is quadratic when many boundary vertices
        // share an x slab, the grid alone is slow on edges much longer than h.
        std::vector<double> lens(bes.size());
        for (size_t i = 0; i < bes.size(); ++i) lens[i] = dist(bes[i].first, bes[i].second);
        std::nth_element(lens.begin(), lens.begin() + (long)(lens.size() / 2), lens.end());
        const double cellH = std::max(lens[lens.size() / 2], endGap);
        const double sampleGap = cellH;
        const HealCellIndex grid(pts, bvx, cellH);
        const int64_t reach = (int64_t)std::ceil((tau + 0.5 * sampleGap) / cellH);
        const int64_t cellsPerSample = (2 * reach + 1) * (2 * reach + 1) * (2 * reach + 1);
        std::vector<HealCellIndex::Cell> cells;
        std::vector<int> cand;
        std::vector<char> touched(work.size(), 0);
        int did = 0;
        for (size_t oi : ord) {
            const int a = bes[oi].first, b = bes[oi].second, t = beT[oi];
            if (touched[(size_t)t]) continue;
            const double abx = pts[(size_t)b][0] - pts[(size_t)a][0];
            const double aby = pts[(size_t)b][1] - pts[(size_t)a][1];
            const double abz = pts[(size_t)b][2] - pts[(size_t)a][2];
            const double L2 = abx * abx + aby * aby + abz * abz, L = std::sqrt(L2);
            const double lo = std::min(pts[(size_t)a][0], pts[(size_t)b][0]) - tau;
            const double hi = std::max(pts[(size_t)a][0], pts[(size_t)b][0]) + tau;
            auto wLo = std::lower_bound(bvx.begin(), bvx.end(), lo,
                                        [&](int v, double x) { return pts[(size_t)v][0] < x; });
            auto wHi = std::upper_bound(wLo, bvx.end(), hi,
                                        [&](double x, int v) { return x < pts[(size_t)v][0]; });
            const int64_t m = std::max<int64_t>(1, (int64_t)std::ceil(L / sampleGap));
            cand.clear();
            if ((int64_t)(wHi - wLo) <= (m + 1) * cellsPerSample) {
                for (auto it = wLo; it != wHi; ++it) cand.push_back((int)(it - bvx.begin()));
            } else {
                cells.clear();
                for (int64_t k = 0; k <= m; ++k) {
                    const double f = (double)k / (double)m;
                    const HealCellIndex::Cell c = grid.cellOf({ pts[(size_t)a][0] + f * abx,
                                                                pts[(size_t)a][1] + f * aby,
                                                                pts[(size_t)a][2] + f * abz });
                    for (int64_t dx = -reach; dx <= reach; ++dx)
                        for (int64_t dy = -reach; dy <= reach; ++dy)
                            for (int64_t dz = -reach; dz <= reach; ++dz)
                                cells.push_back({ c[0] + dx, c[1] + dy, c[2] + dz });
                }
                std::sort(cells.begin(), cells.end());
                cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
                for (const HealCellIndex::Cell& c : cells) grid.visit(c, [&](int p) { cand.push_back(p); });
                std::sort(cand.begin(), cand.end());
            }
            int best = -1;
            double bestT = 2, bestD = 0;
            for (int p : cand) {
                const int v = bvx[(size_t)p];
                if (v == a || v == b) continue;
                const double px = pts[(size_t)v][0] - pts[(size_t)a][0];
                const double py = pts[(size_t)v][1] - pts[(size_t)a][1];
                const double pz = pts[(size_t)v][2] - pts[(size_t)a][2];
                const double tt = (px * abx + py * aby + pz * abz) / L2;
                if (!(tt * L > endGap && (1 - tt) * L > endGap)) continue;
                const double cx = px - tt * abx, cy = py - tt * aby, cz = pz - tt * abz;
                const double d = std::sqrt(cx * cx + cy * cy + cz * cz);
                if (d <= tau && tt < bestT) { best = v; bestT = tt; bestD = d; }
            }
            if (best < 0) continue;
            int s = 0;
            for (; s < 3; ++s)
                if (work[(size_t)t][s] == a && work[(size_t)t][(s + 1) % 3] == b) break;
            const int c = work[(size_t)t][(s + 2) % 3];
            if (best == c) continue;
            work[(size_t)t] = { a, best, c };
            work.push_back({ best, b, c });
            dead.push_back(0);
            touched[(size_t)t] = 1;
            touched.push_back(1);
            clJoin(best, a);
            splitRec.push_back({ best, bestD });
            did++;
        }
        if (!did) break;
    }

    // ---- certificate, per stitched cluster ------------------------------------
    census(E);
    struct ClusterCount { int open0 = 0, nm0 = 0, conf0 = 0, open1 = 0, nm1 = 0, conf1 = 0; };
    std::unordered_map<int, ClusterCount> cc;
    for (const auto& kv : E0) {
        ClusterCount& c = cc[clFind(keyLo(kv.first))];
        if (kv.second.cnt == 1) c.open0++;
        else if (kv.second.cnt > 2) c.nm0++;
        else if (kv.second.fwd != 1) c.conf0++;
    }
    for (const auto& kv : E) {
        ClusterCount& c = cc[clFind(keyLo(kv.first))];
        if (kv.second.cnt == 1) c.open1++;
        else if (kv.second.cnt > 2) c.nm1++;
        else if (kv.second.fwd != 1) c.conf1++;
    }
    std::vector<int> roots;
    for (const auto& kv : cc)
        if (kv.second.open0 > 0) roots.push_back(kv.first);
    std::sort(roots.begin(), roots.end());
    std::unordered_map<int, char> healed;
    for (int root : roots) {
        const ClusterCount& c = cc.at(root);
        const bool ok = c.nm1 == 0 && c.conf1 == 0 && c.open1 < c.open0;
        healed[root] = ok ? 1 : 0;
        if (ok) { st.clustersHealed++; st.boundaryAfter += c.open1; }
        else { st.clustersRefused++; st.boundaryAfter += c.open0; }
    }
    auto isHealed = [&](int v) {
        auto it = healed.find(clFind(v));
        return it != healed.end() && it->second;
    };
    for (const OpRec& w : weldRec)
        if (isHealed(w.vertex)) { st.welds++; st.weldMax = std::max(st.weldMax, w.d); }
    for (const OpRec& s : splitRec)
        if (isHealed(s.vertex)) { st.splits++; st.splitMax = std::max(st.splitMax, s.d); }

    if (st.clustersHealed > 0) {
        // Healed clusters take the stitched triangles; every other cluster keeps
        // its input triangles at their input positions.
        const size_t n0 = tris.size();
        std::vector<std::array<int, 3>> out;
        out.reserve(work.size());
        for (size_t k = 0; k < work.size(); ++k) {
            if (isHealed(k < n0 ? tris[k][0] : work[k][0])) {
                if (!dead[k]) out.push_back(work[k]);
                else st.collapsed++;
            } else if (k < n0) {
                out.push_back(tris[k]);
            }
        }
        tris.swap(out);
    }

    if (diag) {
        std::fprintf(stderr,
                     "DIAG_HEAL boundaryBefore=%d boundaryAfter=%d welds=%d weldMaxTau=%.4f "
                     "collapsed=%d splits=%d splitMaxTau=%.4f refusedWeldClusters=%d "
                     "clustersHealed=%d clustersRefused=%d sMin=%.9g sMinTau=%.4g sewTolTau=%.4g "
                     "stitchTol=%.9g stitchTolTau=%.4g tau=%.9g\n",
                     st.boundaryBefore, st.boundaryAfter, st.welds, st.weldMax / tau, st.collapsed,
                     st.splits, st.splitMax / tau, st.refusedWeldClusters, st.clustersHealed,
                     st.clustersRefused, st.sMin, st.sMin / tau, sewTol / tau, tol, tol / tau, tau);
        for (int root : roots) {
            const ClusterCount& c = cc.at(root);
            std::fprintf(stderr,
                         "DIAG_HEAL_CLUSTER root=%d open=%d->%d nonManifold=%d->%d "
                         "conflicts=%d->%d healed=%d\n",
                         root, c.open0, c.open1, c.nm0, c.nm1, c.conf0, c.conf1,
                         healed.at(root) ? 1 : 0);
        }
    }
    return st;
}

}  // namespace detail

namespace {

namespace fs = std::filesystem;

constexpr double kPi = 3.14159265358979323846;  // avoids MSVC's _USE_MATH_DEFINES dance

// Result has no slot for an additive key (the public header is frozen).
// convert() records the count under the input/output pair; toJson() reads it
// back. A Result that never went through convert() emits no key.
std::mutex gAbsorbMu;
std::unordered_map<std::string, int> gAbsorbCount;

std::string absorbIoKey(const std::string& in, const std::string& out) {
    std::string k;
    k.reserve(in.size() + out.size() + 1);
    k += in;
    k.push_back('\n');
    k += out;
    return k;
}

void noteMeshTJunctionsAbsorbed(const std::string& in, const std::string& out, int n) {
    std::lock_guard<std::mutex> lk(gAbsorbMu);
    gAbsorbCount[absorbIoKey(in, out)] = n;
}

int meshTJunctionsAbsorbedOf(const std::string& in, const std::string& out) {
    std::lock_guard<std::mutex> lk(gAbsorbMu);
    auto it = gAbsorbCount.find(absorbIoKey(in, out));
    return it == gAbsorbCount.end() ? -1 : it->second;
}

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

// Bounded numeric append. Never uses snprintf(nullptr, 0) + string((size_t)n)
// (MSVC: n == -1 wraps the size to 0 and the subsequent write overflows the
// 16-byte SSO into whatever follows, including main's /GS cookie).
void jsonAppendFmt(std::string& o, const char* fmt, ...) STL2STEP_PRINTF(2, 3);
void jsonAppendFmt(std::string& o, const char* fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    if (n < 0) return;
    if ((size_t)n < sizeof buf) {
        o.append(buf, (size_t)n);
        return;
    }
    std::string tmp((size_t)n + 1, '\0');
    va_start(ap, fmt);
    vsnprintf(&tmp[0], tmp.size(), fmt, ap);
    va_end(ap);
    o.append(tmp.c_str(), (size_t)n);
}

// Chunked parallel-for over [0, n). Serial below 4096 items (thread spawn costs
// more than tiny workloads); each worker claims 1024-item chunks off an atomic.
// Construction goes through detail::runPool (D-142-1): threads<=1 => zero std::thread.
template <typename F>
void parallelFor(unsigned threads, size_t n, F&& fn) {
    if (threads <= 1 || n < 4096) {
        for (size_t i = 0; i < n; i++) fn(i);
        return;
    }
    unsigned k = (unsigned)std::min<size_t>(threads, (n + 4095) / 4096);
    std::atomic<size_t> next{ 0 };
    detail::runPool(k, [&](unsigned) {
        for (;;) {
            size_t base = next.fetch_add(1024);
            if (base >= n) break;
            size_t end = std::min(n, base + 1024);
            for (size_t i = base; i < end; i++) fn(i);
        }
    });
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
    // Adaptive 2D Gauss (Eps overload). The no-Eps VolumeProperties uses a
    // fixed scheme that under-integrates a cylindrical face whose inner/cap
    // wire is a many-span polyline: B2's R8 wall with two 68-segment windows
    // reports 62410.25 mm³ against a true B-Rep volume of 62445.89 (37 mm³ /
    // 0.0598 % artifact). Adaptive areas of those faces match the exact
    // Steinmetz construction to 0.3 mm²; the residual 1.7 mm³ is the D-130-2
    // tier-2 chord, counted not absorbed. Precision::Confusion() is OCCT's
    // existing linear precision (the adaptive switch is Eps <= 0.001), not a
    // new geometry threshold. (D-130-20(2))
    BRepGProp::VolumeProperties(s, props, Precision::Confusion());
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
        std::sort(r.warnings.begin(), r.warnings.end());  // D-130-10(2)
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

    if (!opt.dxfDir.empty()) {
        std::error_code ec;
        fs::create_directories(opt.dxfDir, ec);
    }

    unsigned hw = detail::resolveThreadCount(opt.threads);
    detail::ThreadBudgetScope threadBudget(opt.threads);

    // Silence OCCT's own console chatter; convert crashes/signals into C++ exceptions.
    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
    OSD::SetSignal(Standard_False);

    Timer timer, total;
    int solidsOut = 0, openShellsOut = 0, facesBefore = 0, facesAfter = 0;
    int nTri = 0, nVert = 0, nComp = 0;
    int meshTJunctionsAbsorbed = 0;
    long long healTriDelta = 0;   // stitched triangles minus input (D-fc-5.10)
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

        // D-130-12: the mesh file's own coordinate quantization, read from the
        // FILE (the loaded points are doubles and no longer say what format they
        // arrived in), then carried into the same mm the welded points live in.
        const StlQuantFloor qf = stlQuantFloor(inPath);
        const double quantFloor = qf.q * std::fabs(scale);
        if (const char* e = std::getenv("STL2STEP_LAWBAND_DIAG"); e && e[0] && e[0] != '0')
            std::fprintf(stderr,
                         "DIAG_QUANTFLOOR %s maxAbs=%.9g q=%.9g fmt=%s res=%.9g "
                         "sig=%d scale=%.9g ok=%d\n",
                         inPath.c_str(), qf.maxAbs * std::fabs(scale), quantFloor,
                         qf.ascii ? "ascii" : "binary", qf.res * std::fabs(scale),
                         qf.sigDigits, scale, qf.ok ? 1 : 0);

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

        // T-junction absorb, before ekey / UnionFind / chartsA1. Heals the mesh
        // the segmenter and the faceted builder both see.
        {
            std::vector<std::array<int, 3>> raw;
            raw.reserve(tris.size());
            for (const Tri& t : tris) raw.push_back({ t.a, t.b, t.c });
            std::vector<std::array<double, 3>> xyz;
            xyz.reserve(pts.size());
            for (const gp_XYZ& p : pts) xyz.push_back({ p.X(), p.Y(), p.Z() });
            const char* envAbsorb = std::getenv("STL2STEP_DIAG_ABSORB");
            const bool absorbDiag = envAbsorb && envAbsorb[0] && envAbsorb[0] != '0';
            meshTJunctionsAbsorbed =
                detail::absorbZeroAreaTJunctions(raw, xyz, quantFloor, absorbDiag);
            tris.clear();
            tris.reserve(raw.size());
            for (const auto& t : raw) tris.push_back({ t[0], t[1], t[2] });
        }

        auto ekey = [](int a, int b) {
            return a < b ? ((uint64_t)(uint32_t)a << 32) | (uint32_t)b
                         : ((uint64_t)(uint32_t)b << 32) | (uint32_t)a;
        };

        // Pass A: global edge use counts. Pass B: union triangles ONLY across edges
        // used exactly twice, so components split at non-manifold junctions: two
        // solids glued along an edge become two clean manifold bodies (direct path)
        // instead of one dirty component that would need the slow sewing repair.
        UnionFind uf(0);
        {
            struct EdgeStat { uint32_t cnt = 0; int t1 = -1, t2 = -1; };
            std::unordered_map<uint64_t, EdgeStat> edges;
            auto passA = [&]() {
                edges = std::unordered_map<uint64_t, EdgeStat>();
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
            };
            passA();

            // Stitch open boundaries (D-fc-5.10), TrueForm path only: Verbatim
            // stays byte-identical and --force-sew stays a diagnostic. Runs only
            // when Pass A sees an edge used once, so a closed mesh does no more work.
            if (smooth && !forceSew
                && std::any_of(edges.begin(), edges.end(),
                               [](const auto& kv) { return kv.second.cnt == 1; })) {
                std::vector<std::array<int, 3>> raw;
                raw.reserve(tris.size());
                for (const Tri& t : tris) raw.push_back({ t.a, t.b, t.c });
                std::vector<std::array<double, 3>> xyz;
                xyz.reserve(pts.size());
                for (const gp_XYZ& p : pts) xyz.push_back({ p.X(), p.Y(), p.Z() });
                const char* envRel = std::getenv("STL2STEP_DIAG_REL");
                const bool healDiag = envRel && envRel[0] && envRel[0] != '0';
                const detail::HealStats hs =
                    detail::healOpenBoundaries(raw, xyz, quantFloor, sewTol, healDiag);
                if (hs.clustersHealed > 0) {
                    healTriDelta = (long long)raw.size() - (long long)tris.size();
                    tris.clear();
                    tris.reserve(raw.size());
                    for (const auto& t : raw) tris.push_back({ t[0], t[1], t[2] });
                    passA();
                    note("  heal      %d open edges -> %d (%d welds, %d splits, stitch %.2g mm)\n",
                         hs.boundaryBefore, hs.boundaryAfter, hs.welds, hs.splits, hs.stitchTol);
                }
            }

            uf = UnionFind((int)tris.size());
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
            mv.quantFloor = quantFloor;
        };

        std::unordered_map<int, refit::RegionSet> refitPlans;
        refit::RefitStats refitTotals;
        int smoothSkippedComponents = 0;
        double dVolPredSigned = 0, dVolPredAbs = 0;
        refit::SegmentParams segp;
        segp.epsPlane = opt.smoothTolMM;
        segp.thetaPlaneDeg = opt.smoothAngleDeg;
        segp.doFillets = opt.smoothFillets;
        segp.requestedThreads = opt.threads;
        if (smooth && !forceSew) {
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
            refit::prismBindDxfDir(opt.dxfDir);
            refit::prismBindDxfStem(fs::path(opt.input).stem().string());
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
                std::vector<std::string> torusProbeWarn;
                auto tryBuild = [&](refit::RegionSet& rset, bool capture) -> bool {
                    rf.clear();
                    return refit::buildFaces(
                               mv, rset, verts, rf,
                               [&](const std::string& m) {
                                   if (capture) torusProbeWarn.push_back(m);
                                   else warn(m);
                               })
                        && !rf.empty();
                };
                const bool torusProbe = rs.torusRevertValid && !rs.torusRevertApplied;
                torusProbeWarn.clear();
                // The shipped shell's own test: closed, valid, inside the volume
                // budget. U-3 covers every torus (D-fc-1.3), so a probe whose
                // faces built but whose shell fails this test takes the pre-T
                // snapshot too, before the verbatim revert (D-fc-4.6).
                auto shellOk = [&]() -> bool {
                    TopoDS_Shell probe;
                    BRep_Builder pb;
                    pb.MakeShell(probe);
                    for (auto& f : rf) pb.Add(probe, f);
                    bool good = BRep_Tool::IsClosed(probe);
                    if (good) {
                        try {
                            BRepCheck_Analyzer an(probe, Standard_True, Standard_True);
                            good = an.IsValid();
                        } catch (const Standard_Failure&) {
                            good = false;
                        }
                    }
                    if (good) {
                        double dVolAbs = 0;
                        for (const auto& reg : rs.regions)
                            dVolAbs += std::fabs(reg.dVolPredicted);
                        const double meshVol = std::fabs(cs.vol);
                        const double budget = std::max(1e-4 * meshVol, 3.0 * dVolAbs);
                        try {
                            const double shellVol = shapeVolume(probe);
                            if (std::fabs(shellVol - meshVol) > budget) good = false;
                        } catch (const Standard_Failure&) {
                            good = false;
                        }
                    }
                    return good;
                };
                bool ok = tryBuild(rs, torusProbe);
                bool shellChecked = false;
                if (ok && torusProbe) {
                    ok = shellOk();
                    shellChecked = true;
                }
                // U-3 per torus (D-fc-4.11 (a)): the tori the probe could not ship
                // are withheld -- the stages from T on gather them but leave their
                // triangles and atoms in their pre-T state -- and the component is
                // built again with every other torus. Each pass withholds at least
                // one more torus; when none is left, or a failed torus cannot be
                // withheld, the pre-T snapshot below decides as before.
                std::vector<int> withheld;
                auto heldTorus = [&](const refit::Region& reg) {
                    for (int t : reg.tris)
                        if (std::binary_search(withheld.begin(), withheld.end(), t)) return true;
                    return false;
                };
                while (!ok && torusProbe && !rs.torusFailedTris.empty()) {
                    withheld.insert(withheld.end(), rs.torusFailedTris.begin(),
                                    rs.torusFailedTris.end());
                    std::sort(withheld.begin(), withheld.end());
                    withheld.erase(std::unique(withheld.begin(), withheld.end()),
                                   withheld.end());
                    // Lossless prune: a pass that withholds every torus of this set
                    // leaves none (or keeps one it cannot withhold), so it ends here.
                    bool anyFree = false;
                    for (const auto& reg : rs.regions)
                        if (reg.origin == refit::Origin::TorusBlend && !heldTorus(reg))
                            anyFree = true;
                    if (!anyFree) break;
                    refit::SegmentParams heldp = segp;
                    heldp.torusWithheldTris = withheld;
                    refit::RegionSet held;
                    if (!refit::resegmentTorus(mv, heldp, rs, held, nullptr)) break;
                    bool anyTorus = false, stillHeld = false;
                    for (const auto& reg : held.regions) {
                        if (reg.origin != refit::Origin::TorusBlend) continue;
                        anyTorus = true;
                        if (heldTorus(reg)) stillHeld = true;
                    }
                    if (!anyTorus || stillHeld || !held.torusRevertValid) break;
                    parallelFor(subThreads, nV, [&](size_t i) {
                        BRep_Builder bb;
                        bb.MakeVertex(verts[i], gp_Pnt(pts[cs.vtx[i]]), Precision::Confusion());
                    });
                    held.compRoot = rs.compRoot;
                    rs = std::move(held);
                    torusProbeWarn.clear();
                    ok = tryBuild(rs, true);
                    shellChecked = false;
                    if (ok) {
                        ok = shellOk();
                        shellChecked = true;
                    }
                }
                if (!ok && torusProbe) {
                    shellChecked = false;
                    torusProbeWarn.clear();
                    parallelFor(subThreads, nV, [&](size_t i) {
                        BRep_Builder bb;
                        bb.MakeVertex(verts[i], gp_Pnt(pts[cs.vtx[i]]), Precision::Confusion());
                    });
                    refit::RegionSet pre;
                    pre.compRoot = rs.compRoot;
                    pre.regions = rs.torusRevert.regions;
                    pre.chains = rs.torusRevert.chains;
                    pre.triRegion = rs.torusRevert.triRegion;
                    pre.triIsland = rs.torusRevert.triIsland;
                    pre.nIslands = rs.torusRevert.nIslands;
                    pre.rejected = rs.torusRevert.rejected;
                    pre.torusRevertApplied = true;
                    // Discard the probe's construction census the same way
                    // torusProbeWarn is discarded: RESULT must count the
                    // shipped (pre-T) pass, not the U-3 attempt that reverted.
                    pre.stats = rs.torusRevert.stats;
                    pre.stats.tori = 0;
                    rs = pre;
                    ok = tryBuild(rs, false);
                } else if (ok && torusProbe) {
                    for (const auto& m : torusProbeWarn) warn(m);
                }
                if (ok) {
                    if (!shellChecked) ok = shellOk();
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
            std::atomic<long long> degenSkipped{ 0 };
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
                if (detail::triangleAreaDegenerate(A, Bp, Cp)) {
                    skipped++;
                    degenSkipped++;
                    return;
                }
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
                // A skipped degenerate is not a face, so the shell is open.
                // knownClosed would wrap it as a solid and solidsOut would
                // count a body ShapeFix rewrites as an OPEN_SHELL.
                const bool knownClosed = !usedRefit && degenSkipped.load() == 0;
                finishShell(sh, knownClosed, "direct");
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
            detail::runPool(nOuter, [&](unsigned) {
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
            refitTotals.edgeAnalytic += st.edgeAnalytic;
            refitTotals.edgePolylineTier2 += st.edgePolylineTier2;
            refitTotals.edgeUnhandled += st.edgeUnhandled;
            refitTotals.edgeOverTol += st.edgeOverTol;
            refitTotals.edgeOverCap += st.edgeOverCap;
            refitTotals.radiusN += st.radiusN;
            refitTotals.radiusMaxAbs = std::max(refitTotals.radiusMaxAbs, st.radiusMaxAbs);
            refitTotals.radiusMaxRel = std::max(refitTotals.radiusMaxRel, st.radiusMaxRel);
            dVolPredSigned += st.dVolPredicted;
            dVolPredAbs += o.refitDVolAbs;
        }
        {
            long long idxDegen = (long long)nTri + healTriDelta - (long long)tris.size();
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
            auto unifyOne = [&](TopoDS_Shape& s, const std::vector<TopoDS_Face>* keep,
                                double angleDeg) {
                ShapeUpgrade_UnifySameDomain usd(s, Standard_True, Standard_True, Standard_False);
                if (angleDeg > 0) {
                    usd.SetAngularTolerance(angleDeg * kPi / 180.0);
                    usd.SetLinearTolerance(Precision::Confusion());
                }
                if (keep)
                    for (const auto& f : *keep) usd.KeepShape(f);
                usd.Build();
                s = usd.Shape();
            };
            auto unifyAll = [&](double angleDeg, const char* label) {
                const int before = countShapes(shape, TopAbs_FACE);
                if (parts.size() > 1) {
                    unsigned uThreads = (unsigned)std::min<size_t>(
                        { (size_t)hw, parts.size(), (size_t)10 });
                    std::atomic<size_t> nextPart{ 0 };
                    std::atomic<bool> anyFail{ false };
                    detail::runPool(uThreads, [&](unsigned) {
                        for (;;) {
                            size_t i = nextPart.fetch_add(1);
                            if (i >= parts.size()) break;
                            try { unifyOne(parts[i], partKeep[i], angleDeg); }
                            catch (...) { anyFail = true; }   // that body keeps its facets
                        }
                    });
                    TopoDS_Compound comp;
                    BRep_Builder B;
                    B.MakeCompound(comp);
                    for (auto& p : parts) B.Add(comp, p);
                    shape = comp;
                    if (anyFail)
                        warn("coplanar merge failed on at least one body -- it keeps per-triangle faces");
                    facesAfter = countShapes(shape, TopAbs_FACE);
                    note("  %s     %s -> %s faces (%u threads)  [%.2fs]\n",
                         label, fmtInt(before).c_str(), fmtInt(facesAfter).c_str(),
                         uThreads, timer.lap());
                } else {
                    unifyOne(shape, partKeep.empty() ? nullptr : partKeep[0], angleDeg);
                    parts[0] = shape;
                    facesAfter = countShapes(shape, TopAbs_FACE);
                    note("  %s     %s -> %s faces  [%.2fs]\n",
                         label, fmtInt(before).c_str(), fmtInt(facesAfter).c_str(), timer.lap());
                }
            };
            try {
                unifyAll(unifyAngleDeg, "unify");
                // TrueForm-only: mesh-jitter-aware flat consolidation on faceted
                // islands (revert fallback and analytic+facet mixes). Verbatim
                // defaults stay at unifyAngleDeg (G0.1 byte-identity).
                // Coarse Fusion/STLB band only (500–1200 tris; sync with
                // coarseFusionBand in refit_internal.hpp). Large real-CAD meshes
                // (e.g. Body11) must not get a second unify — it breaks closure.
                if (smooth && !forceSew && nTri >= 500 && nTri <= 1200) {
                    double segMaxDev = refitTotals.maxVertexDev;
                    for (const auto& kv : refitPlans)
                        segMaxDev = std::max(segMaxDev, kv.second.stats.maxVertexDev);
                    const double epsMesh = std::max({weldTol, 1e-4 * diag, 1e-3});
                    const double dev = std::max(segMaxDev, epsMesh);
                    const double L = std::max(diag, dev);
                    const double smoothFlatDeg =
                        std::max(0.01, dev / L * (180.0 / kPi));
                    if (smoothFlatDeg > unifyAngleDeg)
                        unifyAll(smoothFlatDeg, "smooth-flat");
                }
            } catch (const Standard_Failure& f) {
                warn(std::string("coplanar merge failed (") +
                     (f.GetMessageString() ? f.GetMessageString() : "?") +
                     ") -- keeping per-triangle faces");
            }
        }

        int builtPl = 0, builtCy = 0, builtCn = 0, builtTo = 0, builtFi = 0, builtCo = 0, revCo = 0;
        if (smooth) {
            size_t pi = 0;
            for (size_t ci = 0; ci < order.size(); ++ci) {
                const CompOut& o = outs[ci];
                int cyls = 0, cones = 0, tori = 0;
                for (size_t k = 0; k < o.parts.size(); ++k) {
                    for (TopExp_Explorer ex(parts[pi], TopAbs_FACE); ex.More(); ex.Next()) {
                        // D-130-10(1): the shipped cone faces are counted the
                        // same way the shipped cylinders are -- off the written
                        // shape, not off a construction counter.
                        const GeomAbs_SurfaceType t =
                            BRepAdaptor_Surface(TopoDS::Face(ex.Current()), Standard_False)
                                .GetType();
                        if (t == GeomAbs_Cylinder) cyls++;
                        else if (t == GeomAbs_Cone) cones++;
                        else if (t == GeomAbs_Torus) tori++;
                    }
                    ++pi;
                }
                if (!refitPlans.count(order[ci])) continue;
                const refit::RefitStats& st = o.refitSt;
                if (o.usedRefit && (cyls > 0 || (st.cylinders == 0 && st.planes > 0))) {
                    builtCo++;
                    builtPl += st.planes;
                    builtCy += cyls;
                    builtCn += cones;
                    builtTo += tori;
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
                detail::overlap(opt.threads, runAnalysis, [&]() { wst = writeStep(); });
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
        r.threadsUsed = (int)hw;
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
        // D-130-10(2): component builds run on a thread pool, so the arrival
        // order of warnings is scheduler state, not a property of the mesh
        // (S03 was measured emitting the same 34 strings in six orders with
        // identical STEP). Sorting makes RESULT deterministic by construction.
        r.warnings = warnings;
        std::sort(r.warnings.begin(), r.warnings.end());
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
            r.smoothBuiltCones = builtCn;
            r.smoothBuiltTori = builtTo;
            r.smoothBuiltFillets = builtFi;
            r.smoothBuiltComponents = builtCo;
            r.smoothRevertedComponents = revCo;
            r.edgeClassAnalytic = refitTotals.edgeAnalytic;
            r.edgeClassPolylineTier2 = refitTotals.edgePolylineTier2;
            r.edgeClassUnhandled = refitTotals.edgeUnhandled;
            r.edgeClassOverTol = refitTotals.edgeOverTol;
            r.edgeClassOverCap = refitTotals.edgeOverCap;
            r.radiusDriftN = refitTotals.radiusN;
            r.radiusDriftMaxAbs = refitTotals.radiusMaxAbs;
            r.radiusDriftMaxRel = refitTotals.radiusMaxRel;
        }
        noteMeshTJunctionsAbsorbed(r.input, r.output, meshTJunctionsAbsorbed);
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
//
// Built by append (same shape as MeshResult::toJson). The previous two-pass
// snprintf(nullptr,0) + string((size_t)n+1) writer overflowed MSVC's 16-byte
// SSO when n < 0 (size wraps to 0) and smashed main's /GS cookie next to the
// hidden Result (D-I-4, GitHub #7).
std::string Result::toJson() const {
    if (!ok)
        return std::string("{\"ok\":false,\"error\":\"") + jsonEscape(error) + "\"}";

    std::string o;
    o.reserve(512);
    o += "{\"ok\":true,\"input\":\"";
    o += jsonEscape(input);
    o += "\",\"output\":\"";
    o += jsonEscape(output);
    jsonAppendFmt(o,
        "\",\"triangles\":%d,\"vertices\":%d,\"components\":%d,\"solids\":%d,"
        "\"openShells\":%d,\"facesBeforeUnify\":%d,\"facesAfterUnify\":%d,"
        "\"meshVolumeMM3\":%.6f,\"stepVolumeMM3\":%.6f,\"volumeDeltaPct\":%.6f,"
        "\"watertight\":%s,\"seconds\":%.2f,\"warnings\":[",
        triangles, vertices, components, solids, openShells, facesBeforeUnify,
        facesAfterUnify, meshVolumeMM3, stepVolumeMM3, volumeDeltaPct,
        watertight ? "true" : "false", seconds);
    bool first = true;
    for (const auto& w : warnings) {
        if (!first) o += ',';
        first = false;
        o += '"';
        o += jsonEscape(w);
        o += '"';
    }
    o += ']';
    const bool emitSmooth = facesAfterSmooth != 0 || smoothSkippedComponents != 0
        || smoothPlanes != 0 || smoothCylinders != 0 || smoothFillets != 0;
    if (emitSmooth) {
        jsonAppendFmt(o,
            ",\"smoothPlanes\":%d,\"smoothCylinders\":%d,\"smoothFillets\":%d,"
            "\"smoothDistinctRadii\":%d,\"smoothRejected\":%d,\"smoothFacetFaces\":%d,"
            "\"facesAfterSmooth\":%d,\"smoothSkippedComponents\":%d,"
            "\"smoothMaxDevMM\":%.6f,\"smoothMaxEdgeTolMM\":%.6f,"
            "\"smoothVolPredictedMM3\":%.6f,"
            "\"smoothBuiltPlanes\":%d,\"smoothBuiltCylinders\":%d,\"smoothBuiltCones\":%d,"
            "\"smoothBuiltTori\":%d,\"smoothBuiltFillets\":%d,"
            "\"smoothBuiltComponents\":%d,\"smoothRevertedComponents\":%d,"
            "\"edgeClasses\":{\"analytic\":%d,\"polylineTier2\":%d,\"unhandled\":%d,"
            "\"overTol\":%d,\"overCap\":%d},"
            "\"radiusDrift\":{\"n\":%d,\"maxAbs\":%.9f,\"maxRel\":%.9f}",
            smoothPlanes, smoothCylinders, smoothFillets, smoothDistinctRadii,
            smoothRejected, smoothFacetFaces, facesAfterSmooth,
            smoothSkippedComponents, smoothMaxDevMM, smoothMaxEdgeTolMM,
            smoothVolPredictedMM3, smoothBuiltPlanes, smoothBuiltCylinders,
            smoothBuiltCones, smoothBuiltTori, smoothBuiltFillets,
            smoothBuiltComponents, smoothRevertedComponents,
            edgeClassAnalytic, edgeClassPolylineTier2,
            edgeClassUnhandled, edgeClassOverTol, edgeClassOverCap, radiusDriftN,
            radiusDriftMaxAbs, radiusDriftMaxRel);
    }
    const int absorbed = meshTJunctionsAbsorbedOf(input, output);
    // Off-path RESULT stays character-identical to 1.0.0 (AGENTS.md §9).
    // The absorb always runs; the key is a feature result, so emit it only
    // when at least one triangle was absorbed.
    if (absorbed > 0)
        jsonAppendFmt(o, ",\"meshTJunctionsAbsorbed\":%d", absorbed);
    o += '}';
    return o;
}

Result convert(const Options& opt, const LogCallback& log) {
    Converter c(opt, log);
    return c.run();
}

}  // namespace stl2step
