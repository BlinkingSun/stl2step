// D-fc-6.1(d). Bodies for refit_island_absorb.hpp. One definition, linked by
// refit_chains.cpp (which cannot include the header; D5.3) and by the unit test.
// SPDX-License-Identifier: MIT

#include "refit_island_absorb.hpp"

#include <cmath>
#include <unordered_set>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

gp_XYZ islandTriCorner(const MeshView& mv, int t, int k) {
    if (!mv.tris || !mv.compTris || t < 0 || (size_t)t >= mv.nTri) return gp_XYZ(0, 0, 0);
    const int g = mv.compTris[t];
    return mv.pts[mv.tris[g][k]];
}

gp_XYZ islandTriCentroid(const MeshView& mv, int t) {
    return (islandTriCorner(mv, t, 0) + islandTriCorner(mv, t, 1) + islandTriCorner(mv, t, 2)) *
           (1.0 / 3.0);
}

gp_XYZ islandTriNormal(const MeshView& mv, int t) {
    const gp_XYZ A = islandTriCorner(mv, t, 0);
    const gp_XYZ B = islandTriCorner(mv, t, 1);
    const gp_XYZ C = islandTriCorner(mv, t, 2);
    return (B - A).Crossed(C - A);
}

// Local vertex ids (component edge endpoints), the same ids the chain builder uses.
void islandTriLocalVerts(const MeshView& mv, int t, int lv[3]) {
    lv[0] = lv[1] = lv[2] = -1;
    if (!mv.triEdges || !mv.compEdges || !mv.triDirs || t < 0 || (size_t)t >= mv.nTri) return;
    for (int s = 0; s < 3; ++s) {
        const int e = mv.triEdges[t][s];
        if (e < 0 || (size_t)e >= mv.nEdge) return;
        const std::pair<int, int>& ce = mv.compEdges[e];
        lv[s] = ((mv.triDirs[t] >> s) & 1) ? ce.first : ce.second;
    }
}

double radialSagitta(const gp_Ax3& ax, double radius, const gp_XYZ& c) {
    const gp_XYZ a = ax.Direction().XYZ();
    const gp_XYZ d = c - ax.Location().XYZ();
    const double radial = (d - a * d.Dot(a)).Modulus();
    return std::fabs(radial - radius);
}

bool f1EdgeAdjacent(const MeshView& mv, const Region& S, const std::vector<int>& island) {
    if (!mv.triEdges) return false;
    std::unordered_set<int> edges;
    edges.reserve(S.tris.size() * 3 + 1);
    for (int t : S.tris) {
        if (t < 0 || (size_t)t >= mv.nTri) continue;
        for (int k = 0; k < 3; ++k) edges.insert(mv.triEdges[t][k]);
    }
    for (int t : island) {
        if (t < 0 || (size_t)t >= mv.nTri) continue;
        for (int k = 0; k < 3; ++k)
            if (edges.count(mv.triEdges[t][k])) return true;
    }
    return false;
}

}  // namespace

double hostSagitta(const MeshView& mv, const Region& S) {
    double m = 0.0;
    if (!(S.radius > 0.0)) return 0.0;
    for (int t : S.tris) {
        if (t < 0 || (size_t)t >= mv.nTri) continue;
        m = std::max(m, radialSagitta(S.ax, S.radius, islandTriCentroid(mv, t)));
    }
    return m;
}

bool f1WouldAbsorb(const MeshView& mv, const Region& S, const std::vector<int>& island,
                   double lateralSin, double* sagOut) {
    if (sagOut) *sagOut = 0.0;
    if (S.type != SurfType::Cylinder || !S.closed360 || !(S.radius > 0.0)) return false;
    if (island.empty() || !(lateralSin >= 0.0) || !std::isfinite(lateralSin)) return false;
    std::unordered_set<int> verts;
    verts.reserve(S.tris.size() * 3 + 1);
    for (int t : S.tris) {
        int lv[3];
        islandTriLocalVerts(mv, t, lv);
        for (int k = 0; k < 3; ++k)
            if (lv[k] >= 0) verts.insert(lv[k]);
    }
    if (verts.empty()) return false;
    for (int t : island) {
        int lv[3];
        islandTriLocalVerts(mv, t, lv);
        for (int k = 0; k < 3; ++k)
            if (lv[k] < 0 || !verts.count(lv[k])) return false;
    }
    const double cap = hostSagitta(mv, S);
    double sag = 0.0;
    const gp_XYZ axis = S.ax.Direction().XYZ();
    for (int t : island) {
        const gp_XYZ n = islandTriNormal(mv, t);
        const double len = n.Modulus();
        if (!(len > 0.0) || !std::isfinite(len)) return false;
        if (std::fabs(n.Dot(axis) / len) > lateralSin) return false;
        sag = std::max(sag, radialSagitta(S.ax, S.radius, islandTriCentroid(mv, t)));
    }
    if (sagOut) *sagOut = sag;
    if (!(sag <= cap)) return false;
    return f1EdgeAdjacent(mv, S, island);
}

bool uniqueF1Host(const MeshView& mv, const std::vector<Region>& regions,
                  const std::vector<int>& island, double lateralSin, int& host, double& sag,
                  double& hostSag) {
    host = -1;
    sag = 0.0;
    hostSag = 0.0;
    int n = 0;
    int found = -1;
    double foundSag = 0.0;
    for (int i = 0; i < (int)regions.size(); ++i) {
        double s = 0.0;
        if (!f1WouldAbsorb(mv, regions[(size_t)i], island, lateralSin, &s)) continue;
        ++n;
        found = i;
        foundSag = s;
        if (n > 1) return false;
    }
    if (n != 1) return false;
    host = found;
    sag = foundSag;
    hostSag = hostSagitta(mv, regions[(size_t)found]);
    return true;
}

}  // namespace refit
}  // namespace stl2step
