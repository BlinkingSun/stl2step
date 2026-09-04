// stl2step P1 stage D — boundary chains, I8 split, loops, ids, stats.
// Include allowlist: D5.3 (gp_ value types + listed std headers only).
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kDeg3 = 3.0 * M_PI / 180.0;

struct Part {
    int reg = -1;
    int isl = -1;
};

inline bool partEq(Part a, Part b) { return a.reg == b.reg && a.isl == b.isl; }

inline bool partVoid(Part p) { return p.reg < 0 && p.isl < 0; }

// Regions before islands before void; then by id. Used only as a local
// comparison (not a published sort key).
inline bool partLess(Part a, Part b) {
    const int ka = a.reg >= 0 ? 0 : (a.isl >= 0 ? 1 : 2);
    const int kb = b.reg >= 0 ? 0 : (b.isl >= 0 ? 1 : 2);
    if (ka != kb) return ka < kb;
    if (a.reg >= 0) return a.reg < b.reg;
    if (a.isl >= 0) return a.isl < b.isl;
    return false;
}

bool unionDiagOn() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_LAWBAND_DIAG");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

// D-130-14 / D-130-16 -- the union's landing stage (see refit_grow.cpp).
// With it off the loop stitch and the seam generator are exactly what the
// branch tip produced.
bool unionOn() {
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("STL2STEP_UNION");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached != 0;
}

inline int otherVert(const std::pair<int, int>& e, int v) {
    return e.first == v ? e.second : e.first;
}

void localVertsOfTri(const MeshView& mv, int t, int lv[3]) {
    for (int s = 0; s < 3; ++s) {
        const auto& e = mv.compEdges[mv.triEdges[t][s]];
        lv[s] = ((mv.triDirs[t] >> s) & 1) ? e.first : e.second;
    }
}

gp_XYZ triNormal(const MeshView& mv, int t) {
    const int g = mv.compTris[t];
    const gp_XYZ& A = mv.pts[mv.tris[g][0]];
    const gp_XYZ& B = mv.pts[mv.tris[g][1]];
    const gp_XYZ& C = mv.pts[mv.tris[g][2]];
    return (B - A).Crossed(C - A);
}

double triArea(const MeshView& mv, int t) {
    return 0.5 * triNormal(mv, t).Modulus();
}

double regionArea(const MeshView& mv, const std::vector<int>& tris) {
    // I5: accumulate in ascending local triangle id.
    std::vector<int> order = tris;
    std::sort(order.begin(), order.end());
    double a = 0;
    for (int t : order) a += triArea(mv, t);
    return a;
}

int minTriId(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    int m = tris[0];
    for (int t : tris) if (t < m) m = t;
    return m;
}

gp_XYZ localPnt(const MeshView& mv, int lv) {
    return mv.pts[mv.compVtx[lv]];
}

Part labelOfTri(int t, const std::vector<int>& triRegion, const std::vector<int>& triIsland) {
    Part p;
    if (t < 0) return p;
    p.reg = triRegion[t];
    p.isl = triIsland[t];
    return p;
}

// Directed edge vFrom -> vTo on mesh edge e: which incident triangle is on
// the LEFT of the walk (triangle winding). Returns local tri id or -1.
int leftTriangle(const MeshView& mv, int e, int vFrom, int vTo,
                 const std::vector<std::vector<int>>& edgeTris) {
    for (int t : edgeTris[e]) {
        int lv[3];
        localVertsOfTri(mv, t, lv);
        for (int s = 0; s < 3; ++s) {
            const int a = lv[s];
            const int b = lv[(s + 1) % 3];
            if (a == vFrom && b == vTo) return t;
        }
    }
    return -1;
}

int rightTriangle(const MeshView& mv, int e, int vFrom, int vTo,
                  const std::vector<std::vector<int>>& edgeTris) {
    return leftTriangle(mv, e, vTo, vFrom, edgeTris);
}

double distPointPlane(const gp_Ax3& ax, const gp_XYZ& p) {
    const gp_XYZ n = ax.Direction().XYZ();
    return std::fabs(n.Dot(p - ax.Location().XYZ()));
}

double distPointCylinder(const Region& r, const gp_XYZ& p) {
    const gp_XYZ a = r.ax.Direction().XYZ();
    const gp_XYZ d = p - r.ax.Location().XYZ();
    const gp_XYZ radial = d.Crossed(a);
    return std::fabs(radial.Modulus() - r.radius);
}

double distLinePlane(const gp_Ax3& plane, const gp_Ax3& cyl) {
    return distPointPlane(plane, cyl.Location().XYZ());
}

double distAxes(const gp_Ax3& a, const gp_Ax3& b) {
    const gp_XYZ u = a.Direction().XYZ();
    const gp_XYZ v = b.Direction().XYZ();
    const gp_XYZ w = b.Location().XYZ() - a.Location().XYZ();
    const gp_XYZ cr = u.Crossed(v);
    const double crn = cr.Modulus();
    if (crn <= std::sin(kDeg3) + 1e-15) {
        // parallel: |(w × u)|
        return w.Crossed(u).Modulus();
    }
    return std::fabs(w.Dot(cr)) / crn;
}

bool axesParallel3(const gp_Dir& a, const gp_Dir& b) {
    return a.XYZ().Crossed(b.XYZ()).Modulus() <= std::sin(kDeg3) + 1e-15;
}

// Coarse Fusion exports: fitted cylinder axes can be several degrees off
// true parallel; widen the tangent gate so cyl|cyl generator construction
// fires (CYLEDGES lane).
bool axesNearParallel8(const gp_Dir& a, const gp_Dir& b) {
    return a.XYZ().Crossed(b.XYZ()).Modulus() <= std::sin(8.0 * M_PI / 180.0) + 1e-15;
}

bool g1Tangent(const Region& A, const Region& B, double epsPlane) {
    const bool aPln = A.type == SurfType::Plane;
    const bool bPln = B.type == SurfType::Plane;
    const bool aCyl = A.type == SurfType::Cylinder;
    const bool bCyl = B.type == SurfType::Cylinder;

    if (aPln && bPln) {
        return A.ax.Direction().Dot(B.ax.Direction()) >= std::cos(kDeg3) - 1e-15;
    }
    if ((aPln && bCyl) || (aCyl && bPln)) {
        const Region& pl = aPln ? A : B;
        const Region& cy = aCyl ? A : B;
        const double adn = std::fabs(cy.ax.Direction().Dot(pl.ax.Direction()));
        if (adn > std::sin(kDeg3) + 1e-15) return false;
        // Mesh fit residual can exceed the epsPlane floor on coarse exports
        // (handle-lock reg55|17: |dist-R|=0.022 mm vs epsPlane=0.02).
        const double tol =
            std::max(epsPlane, std::max(pl.maxVertexDev, cy.maxVertexDev));
        return std::fabs(distLinePlane(pl.ax, cy.ax) - cy.radius) <= tol;
    }
    if (aCyl && bCyl) {
        if (!axesNearParallel8(A.ax.Direction(), B.ax.Direction())) return false;
        const double tol =
            std::max(epsPlane, std::max(A.maxVertexDev, B.maxVertexDev));
        const double d = distAxes(A.ax, B.ax);
        const double ext = std::fabs(A.radius + B.radius);
        const double inn = std::fabs(std::fabs(A.radius - B.radius));
        return std::fabs(d - ext) <= tol || std::fabs(d - inn) <= tol;
    }
    return false;
}

double chordDevToRegion(const Region& r, const gp_XYZ& p, const gp_XYZ& q) {
    const gp_XYZ m((p.X() + q.X()) * 0.5, (p.Y() + q.Y()) * 0.5, (p.Z() + q.Z()) * 0.5);
    if (r.type == SurfType::Plane) return distPointPlane(r.ax, m);
    if (r.type == SurfType::Cylinder) return distPointCylinder(r, m);
    return 0;
}

struct IslandAcc {
    std::vector<int> tris;
    double area = 0;
    int minTri = 0;
};

struct ChainWalk {
    std::vector<int> edges;
    std::vector<int> verts;
    bool closed = false;
};

// Rotate a closed chain so the lowest local vertex id is first. Direction
// (I3) is preserved. First-edge-id tie-break is applied by the caller
// before calling this (they pick the walk direction).
void rotateClosedToMinVert(ChainWalk& w) {
    if (w.verts.empty()) return;
    int best = 0;
    for (int i = 1; i < (int)w.verts.size(); ++i) {
        if (w.verts[i] < w.verts[best]) best = i;
    }
    if (best == 0) return;
    std::vector<int> ve, ed;
    ve.reserve(w.verts.size());
    ed.reserve(w.edges.size());
    for (int i = 0; i < (int)w.verts.size(); ++i) {
        ve.push_back(w.verts[(best + i) % (int)w.verts.size()]);
        ed.push_back(w.edges[(best + i) % (int)w.edges.size()]);
    }
    w.verts.swap(ve);
    w.edges.swap(ed);
}

void reverseOpen(ChainWalk& w) {
    std::reverse(w.verts.begin(), w.verts.end());
    std::reverse(w.edges.begin(), w.edges.end());
}

// Reverse a closed chain, preserving verts[i] --edges[i]--> verts[(i+1)%n]
// (including the wrap edge) and keeping the start vertex in slot 0.
void reverseClosedKeepStart(ChainWalk& w) {
    const int n = (int)w.verts.size();
    if (n < 2 || (int)w.edges.size() != n) return;
    std::vector<int> ve(n), ed(n);
    ve[0] = w.verts[0];
    for (int i = 1; i < n; ++i) ve[i] = w.verts[n - i];
    ed[0] = w.edges[n - 1];
    for (int i = 1; i < n; ++i) ed[i] = w.edges[n - 1 - i];
    w.verts.swap(ve);
    w.edges.swap(ed);
}

std::pair<double, double> projectPlane(const gp_Ax3& ax, const gp_XYZ& p) {
    const gp_XYZ d = p - ax.Location().XYZ();
    return {d.Dot(ax.XDirection().XYZ()), d.Dot(ax.YDirection().XYZ())};
}

std::pair<double, double> projectCylDev(const Region& r, const gp_XYZ& p, double& prevU, bool first) {
    const gp_XYZ a = r.ax.Direction().XYZ();
    const gp_XYZ d = p - r.ax.Location().XYZ();
    const double v = d.Dot(a);
    const gp_XYZ rho = d - a.Multiplied(v);
    double u = std::atan2(rho.Dot(r.ax.YDirection().XYZ()),
                          rho.Dot(r.ax.XDirection().XYZ()));
    if (!first) {
        while (u - prevU > M_PI) u -= 2.0 * M_PI;
        while (prevU - u > M_PI) u += 2.0 * M_PI;
    }
    prevU = u;
    const double R = r.radius > 0 ? r.radius : 1.0;
    return {R * u, v};
}

double shoelace(const std::vector<std::pair<double, double>>& uv) {
    const int n = (int)uv.size();
    if (n < 3) return 0;
    double a = 0;
    for (int i = 0; i < n; ++i) {
        const auto& p = uv[i];
        const auto& q = uv[(i + 1) % n];
        a += p.first * q.second - q.first * p.second;
    }
    return 0.5 * a;
}

int chainStartVert(const BoundaryChain& c, bool rev) {
    if (c.meshVerts.empty()) return -1;
    if (c.closedLoop) return c.meshVerts.front();
    return rev ? c.meshVerts.back() : c.meshVerts.front();
}

int chainEndVert(const BoundaryChain& c, bool rev) {
    if (c.meshVerts.empty()) return -1;
    if (c.closedLoop) return c.meshVerts.front();
    return rev ? c.meshVerts.front() : c.meshVerts.back();
}

bool loopUsesChain(const Loop& lp, int ci) {
    for (int x : lp.chainIdx) if (x == ci) return true;
    return false;
}

}  // namespace

bool buildTopologyD(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                    SegmentWork& work, RegionSet& out) {
    (void)p;
    const int nTri = (int)mv.nTri;
    const int nVtx = (int)mv.nVtx;
    const int nEdge = (int)mv.nEdge;

    const int savedRoot = out.compRoot;
    out = RegionSet{};
    out.compRoot = savedRoot;
    out.triRegion.assign(nTri, -1);
    out.triIsland.assign(nTri, -1);

    if (nTri < 0 || nVtx < 0 || nEdge < 0) return false;
    if (nTri == 0) {
        out.rejected = work.rejected;
        out.stats.rejected = (int)work.rejected.size();
        return true;
    }

    // --- accepted regions: sort (-area, minLocalTriId), dense ids ----------
    struct Acc {
        Region r;
        double area = 0;
        int minTri = 0;
        int inIdx = 0;
    };
    std::vector<Acc> acc;
    acc.reserve(work.accepted.size());
    for (int i = 0; i < (int)work.accepted.size(); ++i) {
        Acc a;
        a.r = work.accepted[i];
        std::sort(a.r.tris.begin(), a.r.tris.end());
        for (int t : a.r.tris) {
            if (t < 0 || t >= nTri) return false;
        }
        a.area = regionArea(mv, a.r.tris);
        a.minTri = minTriId(a.r.tris);
        a.inIdx = i;
        a.r.loops.clear();
        acc.push_back(std::move(a));
    }
    std::sort(acc.begin(), acc.end(), [](const Acc& a, const Acc& b) {
        if (a.area != b.area) return a.area > b.area;           // -area ascending
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        return a.inIdx < b.inIdx;
    });
    out.regions.resize(acc.size());
    std::vector<int> accRemap(work.accepted.size(), -1);
    for (int i = 0; i < (int)acc.size(); ++i) {
        acc[i].r.id = i;
        if (acc[i].inIdx >= 0 && acc[i].inIdx < (int)accRemap.size())
            accRemap[acc[i].inIdx] = i;
        out.regions[i] = std::move(acc[i].r);
        for (int t : out.regions[i].tris) {
            if (out.triRegion[t] >= 0) return false;            // overlap
            out.triRegion[t] = i;
        }
    }
    std::vector<int> provToDense((int)work.provisionals.size(), -1);
    for (int i = 0; i < (int)work.provisionals.size(); ++i) {
        for (int t : work.provisionals[i].tris) {
            if (t >= 0 && t < nTri && out.triRegion[t] >= 0) {
                provToDense[i] = out.triRegion[t];
                break;
            }
        }
    }
    auto remapFilletNbr = [&](int nb) -> int {
        if (nb < 0) return -1;
        if (nb < (int)provToDense.size() && provToDense[nb] >= 0)
            return provToDense[nb];
        if (nb < (int)accRemap.size() && accRemap[nb] >= 0)
            return accRemap[nb];
        return -1;
    };
    for (Region& r : out.regions) {
        r.filletNbrA = remapFilletNbr(r.filletNbrA);
        r.filletNbrB = remapFilletNbr(r.filletNbrB);
    }

    // --- islands: maximal connected unclaimed, sort (-area, minLocalTriId)
    std::vector<int> uf(nTri);
    for (int i = 0; i < nTri; ++i) uf[i] = i;
    auto find = [&](auto&& self, int x) -> int {
        return uf[x] == x ? x : uf[x] = self(self, uf[x]);
    };
    auto unite = [&](int a, int b) {
        a = find(find, a);
        b = find(find, b);
        if (a == b) return;
        if (a > b) std::swap(a, b);
        uf[b] = a;
    };

    std::vector<std::vector<int>> edgeTris((size_t)nEdge);
    std::vector<int> edgeCount((size_t)nEdge, 0);
    for (int t = 0; t < nTri; ++t) {
        for (int s = 0; s < 3; ++s) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || e >= nEdge) return false;
            edgeTris[e].push_back(t);
            edgeCount[e]++;
        }
    }
    // Unique-and-sort incident tris per edge (a tri listed once per side).
    for (int e = 0; e < nEdge; ++e) {
        auto& v = edgeTris[e];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        edgeCount[e] = (int)v.size();
    }

    for (int e = 0; e < nEdge; ++e) {
        if (edgeCount[e] != 2) continue;
        const int t0 = edgeTris[e][0], t1 = edgeTris[e][1];
        if (out.triRegion[t0] < 0 && out.triRegion[t1] < 0) unite(t0, t1);
    }

    std::vector<int> rootOf(nTri, -1);
    std::vector<IslandAcc> islands;
    std::vector<int> rootIndex;  // parallel to islands, the uf root
    for (int t = 0; t < nTri; ++t) {
        if (out.triRegion[t] >= 0) continue;
        const int r = find(find, t);
        int slot = -1;
        for (int i = 0; i < (int)rootIndex.size(); ++i) {
            if (rootIndex[i] == r) { slot = i; break; }
        }
        if (slot < 0) {
            slot = (int)islands.size();
            rootIndex.push_back(r);
            islands.push_back({});
        }
        islands[slot].tris.push_back(t);
    }
    for (auto& is : islands) {
        std::sort(is.tris.begin(), is.tris.end());
        is.area = 0;
        for (int t : is.tris) is.area += triArea(mv, t);        // already ascending
        is.minTri = is.tris.empty() ? std::numeric_limits<int>::max() : is.tris.front();
    }
    std::vector<int> islOrd(islands.size());
    for (int i = 0; i < (int)islOrd.size(); ++i) islOrd[i] = i;
    std::sort(islOrd.begin(), islOrd.end(), [&](int a, int b) {
        if (islands[a].area != islands[b].area) return islands[a].area > islands[b].area;
        if (islands[a].minTri != islands[b].minTri) return islands[a].minTri < islands[b].minTri;
        return a < b;
    });
    out.nIslands = (int)islands.size();
    for (int id = 0; id < (int)islOrd.size(); ++id) {
        for (int t : islands[islOrd[id]].tris) out.triIsland[t] = id;
    }

    // I1
    for (int t = 0; t < nTri; ++t) {
        const bool inR = out.triRegion[t] >= 0;
        const bool inI = out.triIsland[t] >= 0;
        if (inR == inI) return false;
    }

    // --- I8 split vertices ------------------------------------------------
    std::vector<char> isSplit((size_t)nVtx, 0);
    std::vector<std::vector<int>> vtxTris((size_t)nVtx);
    std::vector<std::vector<int>> vtxEdges((size_t)nVtx);
    for (int e = 0; e < nEdge; ++e) {
        const int a = mv.compEdges[e].first;
        const int b = mv.compEdges[e].second;
        if (a < 0 || a >= nVtx || b < 0 || b >= nVtx) return false;
        vtxEdges[a].push_back(e);
        vtxEdges[b].push_back(e);
        const bool open = edgeCount[e] == 1;
        const bool nm = edgeCount[e] >= 3;
        if (open || nm) {
            isSplit[a] = 1;
            isSplit[b] = 1;
        }
    }
    for (int t = 0; t < nTri; ++t) {
        int lv[3];
        localVertsOfTri(mv, t, lv);
        for (int k = 0; k < 3; ++k) {
            if (lv[k] < 0 || lv[k] >= nVtx) return false;
            vtxTris[lv[k]].push_back(t);
        }
    }
    for (int v = 0; v < nVtx; ++v) {
        auto& ts = vtxTris[v];
        std::sort(ts.begin(), ts.end());
        ts.erase(std::unique(ts.begin(), ts.end()), ts.end());
        // unique Parts at this vertex
        std::vector<Part> parts;
        for (int t : ts) {
            Part lab = labelOfTri(t, out.triRegion, out.triIsland);
            bool seen = false;
            for (const Part& q : parts) if (partEq(q, lab)) { seen = true; break; }
            if (!seen) parts.push_back(lab);
        }
        if ((int)parts.size() >= 3) isSplit[v] = 1;
    }

    auto isBoundaryEdge = [&](int e) -> bool {
        if (edgeCount[e] == 1) return true;                     // open
        if (edgeCount[e] != 2) return false;                    // NM: vertices split, edge not chained
        const Part a = labelOfTri(edgeTris[e][0], out.triRegion, out.triIsland);
        const Part b = labelOfTri(edgeTris[e][1], out.triRegion, out.triIsland);
        return !partEq(a, b);
    };

    std::vector<char> boundary((size_t)nEdge, 0);
    std::vector<int> boundAtV((size_t)nVtx, 0);
    for (int e = 0; e < nEdge; ++e) {
        if (!isBoundaryEdge(e)) continue;
        boundary[e] = 1;
        boundAtV[mv.compEdges[e].first]++;
        boundAtV[mv.compEdges[e].second]++;
    }
    // A vertex with != 2 boundary incidences cannot be walked through.
    for (int v = 0; v < nVtx; ++v) {
        if (boundAtV[v] != 2) isSplit[v] = 1;
    }

    std::vector<char> used((size_t)nEdge, 0);

    auto nextBoundEdge = [&](int v, int incoming) -> int {
        int found = -1;
        for (int e : vtxEdges[v]) {
            if (!boundary[e] || used[e] || e == incoming) continue;
            if (found < 0 || e < found) found = e;              // lowest id if several
        }
        return found;
    };

    auto walkFrom = [&](int startE, int startV) -> ChainWalk {
        ChainWalk w;
        int e = startE;
        int v = startV;
        w.verts.push_back(v);
        while (e >= 0 && !used[e]) {
            used[e] = 1;
            w.edges.push_back(e);
            const int nv = otherVert(mv.compEdges[e], v);
            w.verts.push_back(nv);
            v = nv;
            if (isSplit[v]) break;
            const int ne = nextBoundEdge(v, e);
            if (ne < 0) break;
            e = ne;
        }
        return w;
    };

    std::vector<ChainWalk> walks;

    // Open chains: start at every unused boundary edge leaving a split vertex.
    std::vector<int> splitVerts;
    for (int v = 0; v < nVtx; ++v) if (isSplit[v]) splitVerts.push_back(v);
    for (int v : splitVerts) {
        // Collect unused boundary edges at v, then walk each (sorted by edge id).
        std::vector<int> leaves;
        for (int e : vtxEdges[v]) {
            if (boundary[e] && !used[e]) leaves.push_back(e);
        }
        std::sort(leaves.begin(), leaves.end());
        for (int e : leaves) {
            if (used[e]) continue;
            ChainWalk w = walkFrom(e, v);
            // A walk that leaves a split vertex and comes back to THE SAME
            // vertex has traversed a complete cycle: walkFrom stops at the
            // first split vertex it reaches, so front == back can only mean it
            // closed on itself. Emitting it as an open chain is what makes a
            // region's whole outer loop one closedLoop=false chain when its
            // only neighbour is an island (I7: loops are complete and closed).
            // It is a cycle, so it is recorded as one; the closed branch then
            // applies the same canonical form every other cycle gets (rotate
            // to the min vertex, first step the lower-id edge). Nothing about
            // the mesh changes -- only the chain's own description of itself.
            w.closed = (w.edges.size() >= 3 && w.verts.size() == w.edges.size() + 1 &&
                        w.verts.front() == w.verts.back());
            if (w.closed) w.verts.pop_back();
            if (!w.edges.empty()) walks.push_back(std::move(w));
        }
    }

    // Remaining unused boundary edges are closed cycles with no split vertex.
    for (int e0 = 0; e0 < nEdge; ++e0) {
        if (!boundary[e0] || used[e0]) continue;
        const int a = mv.compEdges[e0].first;
        const int b = mv.compEdges[e0].second;
        const int startV = a < b ? a : b;
        // Direction: first step along the lower-id incident boundary edge,
        // applied after we extract the cycle by walking from e0.
        ChainWalk w;
        int e = e0;
        int v = startV;
        // If startV is not on e0 (shouldn't happen), use endpoint a.
        if (v != a && v != b) v = a;
        // Prefer to begin at a vertex of e0; we will rotate later.
        v = a;  // walk, then rotate to min vertex + first-edge rule
        w.verts.push_back(v);
        const int guard = nEdge + 2;
        int steps = 0;
        while (e >= 0 && !used[e] && steps++ < guard) {
            used[e] = 1;
            w.edges.push_back(e);
            const int nv = otherVert(mv.compEdges[e], v);
            v = nv;
            if (v == w.verts.front() && (int)w.edges.size() > 1) break;
            w.verts.push_back(nv);
            const int ne = nextBoundEdge(v, e);
            e = ne;
            if (ne < 0) break;
            if (ne == e0 && v == w.verts.front()) break;
        }
        // Drop the duplicated close vertex if we pushed it.
        if (!w.verts.empty() && w.verts.size() == w.edges.size() + 1
            && w.verts.back() == w.verts.front()) {
            w.verts.pop_back();
        }
        // Closed: verts == edges (no repeated start).
        if ((int)w.verts.size() == (int)w.edges.size() && !w.edges.empty()) {
            w.closed = true;
            walks.push_back(std::move(w));
        } else if (!w.edges.empty()) {
            w.closed = false;
            walks.push_back(std::move(w));
        }
    }

    // --- materialise BoundaryChain: I3 orientation, D5.6 start vertex -----
    out.chains.reserve(walks.size());
    for (ChainWalk& w : walks) {
        if (w.edges.empty() || w.verts.size() < 2) continue;

        if (!w.closed) {
            // D5.6: start at the lowest terminal vertex id.
            if (w.verts.front() > w.verts.back()) reverseOpen(w);
        } else {
            if ((int)w.verts.size() != (int)w.edges.size()) {
                if ((int)w.verts.size() == (int)w.edges.size() + 1
                    && w.verts.back() == w.verts.front()) {
                    w.verts.pop_back();
                }
            }
            if ((int)w.verts.size() != (int)w.edges.size()) return false;
            rotateClosedToMinVert(w);
            // From the min vertex, the first step is the lower-id incident edge.
            if (w.edges.size() >= 2 && w.edges.back() < w.edges.front()) {
                reverseClosedKeepStart(w);
            }
        }

        const int vFrom = w.verts[0];
        const int vTo = w.verts[1];
        const int e0 = w.edges[0];
        int tL = leftTriangle(mv, e0, vFrom, vTo, edgeTris);
        int tR = rightTriangle(mv, e0, vFrom, vTo, edgeTris);
        // Two-manifold sides are the two incident triangles. I3 puts the
        // winding-left one on A; the other triangle is B even when the mesh
        // winding is conflicting (both "left").
        if ((int)edgeTris[e0].size() == 2) {
            const int t0 = edgeTris[e0][0], t1 = edgeTris[e0][1];
            if (tL < 0 && tR >= 0) tL = (tR == t0) ? t1 : t0;
            else if (tR < 0 && tL >= 0) tR = (tL == t0) ? t1 : t0;
            else if (tL < 0 && tR < 0) { tL = t0; tR = t1; }
            else if (tL >= 0 && tR == tL) tR = (tL == t0) ? t1 : t0;
        }
        const Part L = labelOfTri(tL, out.triRegion, out.triIsland);
        const Part R = labelOfTri(tR, out.triRegion, out.triIsland);

        BoundaryChain ch;
        ch.regA = L.reg;
        ch.regB = R.reg;
        ch.islandA = L.isl;
        ch.islandB = R.isl;
        ch.closedLoop = w.closed;
        ch.meshEdges = w.edges;
        ch.meshVerts = w.verts;
        if (w.closed) {
            if ((int)ch.meshVerts.size() != (int)ch.meshEdges.size()) return false;
        } else {
            if ((int)ch.meshVerts.size() != (int)ch.meshEdges.size() + 1) return false;
        }

        // tangent: geometry, never Origin (D1 §1.6)
        if (ch.regA >= 0 && ch.regB >= 0) {
            ch.tangent = g1Tangent(out.regions[ch.regA], out.regions[ch.regB], tol.epsPlane);
        } else {
            ch.tangent = false;
        }
        out.chains.push_back(std::move(ch));
    }

    // Sort chains by minLocalMeshEdgeId (I2 => unique).
    std::sort(out.chains.begin(), out.chains.end(),
              [](const BoundaryChain& a, const BoundaryChain& b) {
                  int ma = a.meshEdges.empty() ? std::numeric_limits<int>::max() : a.meshEdges[0];
                  for (int e : a.meshEdges) if (e < ma) ma = e;
                  int mb = b.meshEdges.empty() ? std::numeric_limits<int>::max() : b.meshEdges[0];
                  for (int e : b.meshEdges) if (e < mb) mb = e;
                  if (ma != mb) return ma < mb;
                  return false;  // I2: each edge in ≤1 chain, so minEdge is unique
              });

    // --- loop assembly per region -----------------------------------------
    const double areaTie = 1e-12 * (mv.diag > 0 ? mv.diag * mv.diag : 1.0);

    auto classifyAndPush = [&](Region& reg, std::vector<Loop>& loops) -> bool {
        if (loops.empty()) return false;

        struct Scored {
            Loop lp;
            double area = 0;          // signed shoelace (non-360) or unused
            double meanV = 0;         // closed360
            int minCh = 0;
        };
        std::vector<Scored> sc;
        sc.reserve(loops.size());
        for (Loop& lp : loops) {
            Scored s;
            s.lp = std::move(lp);
            s.minCh = s.lp.chainIdx.empty() ? std::numeric_limits<int>::max()
                                            : s.lp.chainIdx[0];
            for (int ci : s.lp.chainIdx) if (ci < s.minCh) s.minCh = ci;

            std::vector<int> vs;
            for (int k = 0; k < (int)s.lp.chainIdx.size(); ++k) {
                const BoundaryChain& c = out.chains[s.lp.chainIdx[k]];
                const bool rev = s.lp.reversed[k] != 0;
                if (c.closedLoop) {
                    if (rev) {
                        vs.push_back(c.meshVerts.front());
                        for (int i = (int)c.meshVerts.size() - 1; i >= 1; --i)
                            vs.push_back(c.meshVerts[i]);
                    } else {
                        vs.insert(vs.end(), c.meshVerts.begin(), c.meshVerts.end());
                    }
                } else {
                    if (!rev) {
                        const int begin = vs.empty() ? 0 : 1;
                        for (int i = begin; i < (int)c.meshVerts.size(); ++i)
                            vs.push_back(c.meshVerts[i]);
                        if (vs.empty()) vs.push_back(c.meshVerts.front());
                    } else {
                        const int last = (int)c.meshVerts.size() - 1;
                        const int begin = vs.empty() ? last : last - 1;
                        for (int i = begin; i >= 0; --i) vs.push_back(c.meshVerts[i]);
                        if (vs.empty()) vs.push_back(c.meshVerts.back());
                    }
                }
            }
            // drop duplicated close
            if (vs.size() >= 2 && vs.front() == vs.back()) vs.pop_back();

            if (reg.closed360) {
                double sum = 0;
                const gp_XYZ a = reg.ax.Direction().XYZ();
                const gp_XYZ o = reg.ax.Location().XYZ();
                int n = 0;
                for (int lv : vs) {
                    sum += (localPnt(mv, lv) - o).Dot(a);
                    ++n;
                }
                s.meanV = n ? sum / n : 0;
            } else {
                std::vector<std::pair<double, double>> uv;
                uv.reserve(vs.size());
                double prevU = 0;
                bool first = true;
                for (int lv : vs) {
                    const gp_XYZ p = localPnt(mv, lv);
                    if (reg.type == SurfType::Plane) {
                        uv.push_back(projectPlane(reg.ax, p));
                    } else {
                        uv.push_back(projectCylDev(reg, p, prevU, first));
                        first = false;
                    }
                }
                s.area = shoelace(uv);
            }
            sc.push_back(std::move(s));
        }

        if (reg.closed360) {
            if (sc.size() < 2) return false;
            int iLow = 0, iHigh = 0;
            for (int i = 1; i < (int)sc.size(); ++i) {
                if (sc[i].meanV < sc[iLow].meanV) iLow = i;
                if (sc[i].meanV > sc[iHigh].meanV) iHigh = i;
            }
            if (iLow == iHigh) return false;
            for (int i = 0; i < (int)sc.size(); ++i) {
                if (i == iLow) sc[i].lp.role = LoopRole::CapLow;
                else if (i == iHigh) sc[i].lp.role = LoopRole::CapHigh;
                else sc[i].lp.role = LoopRole::Inner;
            }
        } else {
            int iOuter = 0;
            double bestAbs = -1;
            for (int i = 0; i < (int)sc.size(); ++i) {
                const double aa = std::fabs(sc[i].area);
                if (aa > bestAbs + areaTie
                    || (std::fabs(aa - bestAbs) <= areaTie && sc[i].minCh < sc[iOuter].minCh)
                    || bestAbs < 0) {
                    bestAbs = aa;
                    iOuter = i;
                }
            }
            // D3: Outer sign should agree with outwardNormal. A mismatch is
            // ChainUnstable on the region, not R3 — ax frame is upstream
            // geometry; we still emit a complete loops vector for P2.
            for (int i = 0; i < (int)sc.size(); ++i)
                sc[i].lp.role = (i == iOuter) ? LoopRole::Outer : LoopRole::Inner;
        }

        std::sort(sc.begin(), sc.end(), [](const Scored& a, const Scored& b) {
            const int ra = (int)a.lp.role, rb = (int)b.lp.role;
            if (ra != rb) return ra < rb;                       // (role, minChainIdx)
            return a.minCh < b.minCh;
        });
        reg.loops.clear();
        for (auto& s : sc) reg.loops.push_back(std::move(s.lp));
        return true;
    };

    // D-130-14 -- THE LOOP-LEVEL UNION, loop half.
    //
    // A region claimed on ONE certified surface may be cut by interruptions
    // that it does not own (the junction facets round a crossing bore, the
    // sawtooth of refused triangles at a mouth). Where two of its pieces meet a
    // vertex from opposite sides, that vertex carries FOUR boundary edges, not
    // two, so `isSplit` breaks the chains there and the stitch below has a
    // choice to make. Picking the lowest chain id closes each piece off into
    // its own loop -- the per-piece answer `partial_recovery_gate` refused.
    //
    // The right pairing is geometric and needs no tolerance: at that vertex the
    // fan alternates the region's own sectors with the interruption's, and the
    // two boundary edges that bound ONE interruption sector are the two ends of
    // ONE wire of the face. Rotating from the incoming edge through the
    // non-member triangles to the next boundary edge names it exactly. Every
    // interruption enclosed by the region then becomes one inner wire, the
    // outer/cap wires stay whole, and an interruption that actually cuts the
    // region leaves the pieces in different regions upstream (they share no
    // vertex) -- separate faces, as D-130-14 requires. At an ordinary boundary
    // vertex (two boundary edges) this returns what the old rule returned.
    auto pairAcrossInterruption = [&](int v, int eIn, const Region& reg) -> int {
        if (v < 0 || v >= nVtx || eIn < 0 || eIn >= nEdge) return -1;
        int t = -1;
        for (int cand : edgeTris[eIn]) {
            const Part lab = labelOfTri(cand, out.triRegion, out.triIsland);
            if (lab.reg != reg.id) {
                t = cand;
                break;
            }
        }
        if (t < 0) return -1;                       // open edge: nothing to rotate through
        int e = eIn;
        const int guard = (int)vtxEdges[v].size() + 2;
        for (int step = 0; step < guard; ++step) {
            int eNext = -1;
            for (int k = 0; k < 3; ++k) {
                const int ee = mv.triEdges[t][k];
                if (ee == e) continue;
                const auto& pe = mv.compEdges[ee];
                if (pe.first == v || pe.second == v) {
                    eNext = ee;
                    break;
                }
            }
            if (eNext < 0) return -1;
            if (boundary[eNext]) return eNext;
            if ((int)edgeTris[eNext].size() != 2) return -1;
            t = (edgeTris[eNext][0] == t) ? edgeTris[eNext][1] : edgeTris[eNext][0];
            e = eNext;
        }
        return -1;
    };
    auto chainEdgeAt = [&](const BoundaryChain& c, bool rev, bool atEnd) -> int {
        if (c.meshEdges.empty()) return -1;
        const bool back = (rev != atEnd);            // rev XOR atEnd == false => back
        return back ? c.meshEdges.front() : c.meshEdges.back();
    };

    for (Region& reg : out.regions) {
        std::vector<char> usedCh(out.chains.size(), 0);
        std::vector<Loop> loops;

        auto revFor = [&](const BoundaryChain& c) -> int {
            if (c.regA == reg.id) return 0;
            if (c.regB == reg.id) return 1;
            return -1;
        };

        // Closed single-chain loops first (each is already a cycle).
        for (int ci = 0; ci < (int)out.chains.size(); ++ci) {
            const BoundaryChain& c = out.chains[ci];
            const int rv = revFor(c);
            if (rv < 0) continue;
            if (!c.closedLoop) continue;
            Loop lp;
            lp.chainIdx.push_back(ci);
            lp.reversed.push_back((uint8_t)rv);
            lp.role = LoopRole::Inner;  // filled in classify
            loops.push_back(std::move(lp));
            usedCh[ci] = 1;
        }

        // Stitch open chains into cycles.
        while (true) {
            int seed = -1;
            for (int ci = 0; ci < (int)out.chains.size(); ++ci) {
                if (usedCh[ci]) continue;
                if (revFor(out.chains[ci]) < 0) continue;
                seed = ci;
                break;
            }
            if (seed < 0) break;

            Loop lp;
            int ci = seed;
            int rv = revFor(out.chains[ci]);
            const int startV = chainStartVert(out.chains[ci], rv != 0);
            const int guard = (int)out.chains.size() + 2;
            int steps = 0;
            bool closedOk = false;
            while (ci >= 0 && !usedCh[ci] && steps++ < guard) {
                usedCh[ci] = 1;
                lp.chainIdx.push_back(ci);
                lp.reversed.push_back((uint8_t)rv);
                const int endV = chainEndVert(out.chains[ci], rv != 0);
                if (endV == startV && (int)lp.chainIdx.size() >= 1) {
                    // A single open chain cannot close unless its terminals match,
                    // which would mean it was actually closed — handled above.
                    if ((int)lp.chainIdx.size() >= 2 || out.chains[ci].closedLoop) {
                        closedOk = true;
                        break;
                    }
                    // 1-chain open with start==end: treat as closed.
                    if (startV == endV) { closedOk = true; break; }
                }
                // Find next unused chain of this region that starts at endV.
                // When several do (a pinch vertex, D-130-14), the successor is
                // the one whose first edge bounds the SAME interruption sector
                // the incoming edge does.
                const int eIn = chainEdgeAt(out.chains[ci], rv != 0, true);
                const int eWant = unionOn() ? pairAcrossInterruption(endV, eIn, reg) : -1;
                int next = -1, nextRev = -1;
                for (int k = 0; k < (int)out.chains.size(); ++k) {
                    if (usedCh[k]) continue;
                    const int rvk = revFor(out.chains[k]);
                    if (rvk < 0) continue;
                    if (chainStartVert(out.chains[k], rvk != 0) != endV) continue;
                    if (eWant >= 0 && chainEdgeAt(out.chains[k], rvk != 0, false) != eWant)
                        continue;
                    if (next < 0 || k < next) {
                        next = k;
                        nextRev = rvk;
                    }
                }
                if (next < 0 && eWant >= 0) {
                    for (int k = 0; k < (int)out.chains.size(); ++k) {
                        if (usedCh[k]) continue;
                        const int rvk = revFor(out.chains[k]);
                        if (rvk < 0) continue;
                        if (chainStartVert(out.chains[k], rvk != 0) != endV) continue;
                        if (next < 0 || k < next) {
                            next = k;
                            nextRev = rvk;
                        }
                    }
                }
                if (next < 0) break;
                ci = next;
                rv = nextRev;
                if (chainEndVert(out.chains[ci], rv != 0) == startV) {
                    usedCh[ci] = 1;
                    lp.chainIdx.push_back(ci);
                    lp.reversed.push_back((uint8_t)rv);
                    closedOk = true;
                    break;
                }
            }
            if (!closedOk || lp.chainIdx.empty()) return false;
            loops.push_back(std::move(lp));
        }

        // Every chain touching this region must appear in exactly one loop (I7/I7b).
        for (int ci = 0; ci < (int)out.chains.size(); ++ci) {
            const int rv = revFor(out.chains[ci]);
            if (rv < 0) continue;
            int hits = 0;
            for (const Loop& lp : loops) if (loopUsesChain(lp, ci)) ++hits;
            if (hits != 1) return false;
        }

        if (!classifyAndPush(reg, loops)) return false;

        // D-130-14: a union face's SEAM is a generator of that same face, so it
        // may not run through an interruption the face carries as an inner
        // wire. The 360 deg seam sits at u = 0 by construction, so u = 0 is
        // moved: XDirection becomes the region's own mesh-vertex generator
        // nearest the middle of the widest azimuth gap between the inner wires.
        // It stays a REAL MESH VERTEX azimuth (the existing rule -- the seam
        // lands on a facet generator, never bisecting a facet); only the choice
        // among them changes, and only for a region that HAS an inner wire, so
        // every region built before this lane keeps the azimuth it had.
        if (unionOn() && reg.closed360 && reg.type == SurfType::Cylinder) {
            std::vector<double> innerU;
            const gp_XYZ o = reg.ax.Location().XYZ();
            const gp_XYZ a = reg.ax.Direction().XYZ();
            const gp_XYZ xd = reg.ax.XDirection().XYZ();
            const gp_XYZ yd = a.Crossed(xd);
            auto azim = [&](const gp_XYZ& p) {
                const gp_XYZ d = p - o;
                return std::atan2(d.Dot(yd), d.Dot(xd));
            };
            for (const Loop& lp : reg.loops) {
                if (lp.role != LoopRole::Inner) continue;
                for (int ci : lp.chainIdx)
                    for (int lv : out.chains[(size_t)ci].meshVerts)
                        innerU.push_back(azim(localPnt(mv, lv)));
            }
            if (!innerU.empty()) {
                std::sort(innerU.begin(), innerU.end());
                double gap = innerU.front() + 2.0 * M_PI - innerU.back();
                double mid = innerU.back() + 0.5 * gap;
                for (size_t i = 1; i < innerU.size(); ++i) {
                    const double g = innerU[i] - innerU[i - 1];
                    if (g > gap) {
                        gap = g;
                        mid = innerU[i - 1] + 0.5 * g;
                    }
                }
                int bestV = -1;
                double bestD = 1e300;
                for (int t : reg.tris) {
                    int lv3[3];
                    localVertsOfTri(mv, t, lv3);
                    for (int k = 0; k < 3; ++k) {
                        double d = std::fabs(azim(localPnt(mv, lv3[k])) - mid);
                        while (d > M_PI) d = std::fabs(d - 2.0 * M_PI);
                        if (d < bestD - 1e-15 || (std::fabs(d - bestD) <= 1e-15 &&
                                                  (bestV < 0 || lv3[k] < bestV))) {
                            bestD = d;
                            bestV = lv3[k];
                        }
                    }
                }
                if (bestV >= 0) {
                    const gp_XYZ d = localPnt(mv, bestV) - o;
                    gp_XYZ rad = d - a * d.Dot(a);
                    if (rad.Modulus() > 0.0) {
                        rad.Normalize();
                        reg.ax = gp_Ax3(reg.ax.Location(),
                                        gp_Dir(a.X(), a.Y(), a.Z()),
                                        gp_Dir(rad.X(), rad.Y(), rad.Z()));
                    }
                }
            }
        }

        // D-130-14 loop census: what the union actually shipped, per region --
        // one cap/outer set plus one inner wire per enclosed interruption.
        if (unionDiagOn() && reg.type == SurfType::Cylinder) {
            std::fprintf(stderr, "DIAG_130_LOOPS rid=%d R=%.6f nTri=%zu closed360=%d loops=%zu",
                         reg.id, reg.radius, reg.tris.size(), reg.closed360 ? 1 : 0,
                         reg.loops.size());
            for (const Loop& lp : reg.loops) {
                size_t nv = 0;
                for (int ci : lp.chainIdx)
                    nv += out.chains[(size_t)ci].meshVerts.size();
                const char* rn = lp.role == LoopRole::Outer     ? "outer"
                                 : lp.role == LoopRole::Inner   ? "inner"
                                 : lp.role == LoopRole::CapLow  ? "capLo"
                                                                : "capHi";
                std::fprintf(stderr, " [%s ch=%zu v=%zu", rn, lp.chainIdx.size(), nv);
                for (int ci : lp.chainIdx) {
                    const BoundaryChain& c = out.chains[(size_t)ci];
                    const gp_XYZ o2 = reg.ax.Location().XYZ();
                    const gp_XYZ a2 = reg.ax.Direction().XYZ();
                    const gp_XYZ x2 = reg.ax.XDirection().XYZ();
                    const gp_XYZ y2 = a2.Crossed(x2);
                    double lo = 1e300, hi = -1e300;
                    for (int lv : c.meshVerts) {
                        const gp_XYZ d2 = localPnt(mv, lv) - o2;
                        const double uu = std::atan2(d2.Dot(y2), d2.Dot(x2));
                        lo = std::min(lo, uu);
                        hi = std::max(hi, uu);
                    }
                    std::fprintf(stderr, " (ci=%d nV=%zu u=[%.1f,%.1f] A=%d B=%d iA=%d iB=%d)", ci,
                                 c.meshVerts.size(), lo * 180.0 / M_PI, hi * 180.0 / M_PI, c.regA,
                                 c.regB, c.islandA, c.islandB);
                }
                std::fprintf(stderr, "]");
            }
            std::fprintf(stderr, "\n");
        }

        if (reg.closed360) {
            int nLow = 0, nHigh = 0, nOuter = 0;
            for (const Loop& lp : reg.loops) {
                if (lp.role == LoopRole::CapLow) ++nLow;
                else if (lp.role == LoopRole::CapHigh) ++nHigh;
                else if (lp.role == LoopRole::Outer) ++nOuter;
            }
            if (nLow != 1 || nHigh != 1 || nOuter != 0) return false;
        } else {
            int nOuter = 0, nCap = 0;
            for (const Loop& lp : reg.loops) {
                if (lp.role == LoopRole::Outer) ++nOuter;
                if (lp.role == LoopRole::CapLow || lp.role == LoopRole::CapHigh) ++nCap;
            }
            if (nOuter != 1 || nCap != 0) return false;
        }
    }

    // FINDING 2: after the area-sort dense ids and chain walk, filletNbrA/B
    // are the two tangent-adjacent regions, not the C1 provisional indices.
    for (Region& r : out.regions) {
        if (r.origin != Origin::FilletStrip) continue;
        std::vector<int> tn;
        for (const BoundaryChain& ch : out.chains) {
            if (!ch.tangent) continue;
            if (ch.regA == r.id && ch.regB >= 0) tn.push_back(ch.regB);
            if (ch.regB == r.id && ch.regA >= 0) tn.push_back(ch.regA);
        }
        std::sort(tn.begin(), tn.end());
        tn.erase(std::unique(tn.begin(), tn.end()), tn.end());
        if (tn.size() >= 2) {
            r.filletNbrA = tn[0];
            r.filletNbrB = tn[1];
        }
    }

    // --- rejected[]: D1.3-A6b / D5.6-A1 ----------------------------------
    // tris = full pre-peel claim set, ascending, non-empty (clearing FORBIDDEN).
    // loops empty by contract. id = rejectOrdinal assigned at record time,
    // unique dense over {0..n-1}, NOT reassigned after the sort.
    // Sort key: (-area, minLocalTriId, id).
    struct Rej {
        Region r;
        double area = 0;
        int minTri = 0;
    };
    std::vector<Rej> rej;
    rej.reserve(work.rejected.size());
    for (int i = 0; i < (int)work.rejected.size(); ++i) {
        Rej x;
        x.r = work.rejected[i];
        std::sort(x.r.tris.begin(), x.r.tris.end());
        x.r.tris.erase(std::unique(x.r.tris.begin(), x.r.tris.end()), x.r.tris.end());
        x.r.loops.clear();
        x.area = regionArea(mv, x.r.tris);
        x.minTri = minTriId(x.r.tris);
        if (x.r.id < 0) x.r.id = i;
        rej.push_back(std::move(x));
    }
    std::sort(rej.begin(), rej.end(), [](const Rej& a, const Rej& b) {
        if (a.area != b.area) return a.area > b.area;
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        return a.r.id < b.r.id;
    });
    out.rejected.resize(rej.size());
    for (int i = 0; i < (int)rej.size(); ++i) {
        Region& rr = rej[i].r;
        rr.filletNbrA = remapFilletNbr(rr.filletNbrA);
        rr.filletNbrB = remapFilletNbr(rr.filletNbrB);
        out.rejected[i] = std::move(rr);
    }

    // --- RefitStats -------------------------------------------------------
    // dVolPredicted here is the SIGNED sum (D4.5 / RESULT). The budget that
    // sums magnitudes is not a RegionSet field; P2/P3 compute it from the
    // per-region values. We do not store dVolPredAbs.
    RefitStats st;
    st.rejected = (int)out.rejected.size();
    st.facetIslands = out.nIslands;
    int facetTris = 0;
    for (int t = 0; t < nTri; ++t) if (out.triIsland[t] >= 0) ++facetTris;
    st.facetTriangles = facetTris;

    std::vector<double> radii;
    for (const Region& r : out.regions) {
        switch (r.origin) {
            case Origin::FilletStrip:
                ++st.fillets;
                break;
            case Origin::CylGrow:
            case Origin::NgonWall:
                if (r.type == SurfType::Plane) ++st.planes;
                else if (r.type == SurfType::Cylinder) ++st.cylinders;
                break;
            case Origin::ChamferCone:
                break;  // not a FilletStrip
            case Origin::PlaneGrow:
                if (r.type == SurfType::Plane) ++st.planes;
                else if (r.type == SurfType::Cylinder) ++st.cylinders;
                break;
        }

        if (r.maxVertexDev > st.maxVertexDev) st.maxVertexDev = r.maxVertexDev;
        st.dVolPredicted += r.dVolPredicted;                    // signed

        bool takeR = false;
        switch (r.origin) {
            case Origin::ChamferCone:
                break;  // two-radius cone, not a fillet/cyl sample
            case Origin::FilletStrip:
            case Origin::CylGrow:
            case Origin::NgonWall:
                takeR = r.radius > 0;
                break;
            case Origin::PlaneGrow:
                takeR = r.type == SurfType::Cylinder && r.radius > 0;
                break;
        }
        if (takeR) radii.push_back(r.radius);
    }
    std::sort(radii.begin(), radii.end());
    if (!radii.empty()) {
        st.distinctRadii = 1;
        for (int i = 1; i < (int)radii.size(); ++i) {
            const double eps = std::max(tol.epsMesh, 1e-9);
            if (std::fabs(radii[i] - radii[i - 1]) > eps) ++st.distinctRadii;
        }
    }

    double maxEdge = 0;
    for (const BoundaryChain& c : out.chains) {
        for (int i = 0; i < (int)c.meshEdges.size(); ++i) {
            int a, b;
            if (c.closedLoop) {
                a = c.meshVerts[i];
                b = c.meshVerts[(i + 1) % (int)c.meshVerts.size()];
            } else {
                a = c.meshVerts[i];
                b = c.meshVerts[i + 1];
            }
            const gp_XYZ pa = localPnt(mv, a);
            const gp_XYZ pb = localPnt(mv, b);
            if (c.regA >= 0) {
                const double d = chordDevToRegion(out.regions[c.regA], pa, pb);
                if (d > maxEdge) maxEdge = d;
            }
            if (c.regB >= 0) {
                const double d = chordDevToRegion(out.regions[c.regB], pa, pb);
                if (d > maxEdge) maxEdge = d;
            }
        }
    }
    st.maxEdgeTol = maxEdge;
    out.stats = st;
    return true;
}

}  // namespace refit
}  // namespace stl2step
