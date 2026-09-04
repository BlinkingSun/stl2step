// stl2step P1 stage A — closed N-gon cylinder claim (docs/1.3.0-linkage-detectors.md).
// Unclaimed A2 planar walls, N>=6, prismatic generators. R = circumradius
// (mean vertex distance to axis), not Eberly. Origin::NgonWall.
// Commit only a through-hole or outer cylindrical boss: closed360, N in [6,256],
// max|ρ−R| ≤ max(sewTol, chordSagitta). Extra reject tests (tryAcceptCycle):
//   N>256 noise; axis-turn blend (internal gens or shared outside ridge);
//   zero cap-plane rims; N<=8 without both rims in cap planes (6-facet bump);
//   one cap end without a chamfer ring on the other end.
// Prefer two-cap-plane holes (collectCapPairCycles) and larger N.
// gp_ / stdlib only. Never throws. Returns true even if zero claims.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kTiny = 1e-30;
constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
constexpr int kNMin = 6;
constexpr int kNMax = 256;  // N>256 is tessellation noise, not a designed hole
constexpr int kMaxCycles = 4096;
constexpr double kCapAreaMul = 4.0;
// Chamfer-cone α window (same as refit_chamfer_cone.cpp) for the other-end ring.
constexpr double kSinChamLo = 0.42261826174069944;  // sin(25°)
constexpr double kSinChamHi = 0.9063077870366499;   // sin(65°)

// Measurement only (STL2STEP_DIAG_130): why a candidate cycle was refused.
// Never changes a decision -- every caller passes the reason it already took.
bool ngonDiagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_130");
    return e && e[0] != '\0' && e[0] != '0';
}

gp_XYZ localTriVert(const MeshView& mv, int lt, int corner) {
    const int gt = mv.compTris[lt];
    const int gv = mv.tris[gt][corner];
    return mv.pts[gv];
}

double triAreaLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    return 0.5 * (b - a).Crossed(c - a).Modulus();
}

gp_XYZ triCentroidLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    return (a + b + c) / 3.0;
}

gp_Dir triNormalLocal(const MeshView& mv, int lt) {
    const gp_XYZ a = localTriVert(mv, lt, 0);
    const gp_XYZ b = localTriVert(mv, lt, 1);
    const gp_XYZ c = localTriVert(mv, lt, 2);
    gp_XYZ n = (b - a).Crossed(c - a);
    const double m = n.Modulus();
    if (m < kTiny) return gp_Dir(0, 0, 1);
    n /= m;
    return gp_Dir(n.X(), n.Y(), n.Z());
}

int minTriId(const Provisional& p) {
    return p.tris.empty() ? std::numeric_limits<int>::max() : p.tris.front();
}

int minTriOf(const std::vector<int>& tris) {
    return tris.empty() ? std::numeric_limits<int>::max() : tris.front();
}

gp_Dir normalizeXYZ(const gp_XYZ& v) {
    const double m = v.Modulus();
    if (m < kTiny) return gp_Dir(1, 0, 0);
    return gp_Dir(v.X() / m, v.Y() / m, v.Z() / m);
}

double eigenAxisSign(const gp_XYZ& w) {
    if (std::abs(w.X()) > 1e-9) return w.X() > 0.0 ? 1.0 : -1.0;
    if (std::abs(w.Y()) > 1e-9) return w.Y() > 0.0 ? 1.0 : -1.0;
    if (std::abs(w.Z()) > 1e-9) return w.Z() > 0.0 ? 1.0 : -1.0;
    return 1.0;
}

gp_Dir canonicalAxis(const gp_XYZ& w) {
    return normalizeXYZ(w * eigenAxisSign(w));
}

bool provEligible(const Provisional& p) {
    return p.claim == ProvClaim::Unclaimed && !p.tris.empty();
}

struct EdgeAdj {
    std::vector<std::array<int, 2>> tri;
};

EdgeAdj buildEdgeAdj(const MeshView& mv) {
    EdgeAdj ea;
    ea.tri.resize(mv.nEdge, { -1, -1 });
    for (int lt = 0; lt < static_cast<int>(mv.nTri); lt++) {
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[lt][s];
            if (e < 0 || static_cast<size_t>(e) >= mv.nEdge) continue;
            if (ea.tri[static_cast<size_t>(e)][0] < 0)
                ea.tri[static_cast<size_t>(e)][0] = lt;
            else
                ea.tri[static_cast<size_t>(e)][1] = lt;
        }
    }
    return ea;
}

gp_XYZ edgeDir(const MeshView& mv, int e) {
    const int v0 = mv.compEdges[e].first;
    const int v1 = mv.compEdges[e].second;
    const gp_XYZ p0 = mv.pts[mv.compVtx[v0]];
    const gp_XYZ p1 = mv.pts[mv.compVtx[v1]];
    return p1 - p0;
}

bool planePlaneIntersect(const gp_Ax3& A, const gp_Ax3& B, gp_XYZ& p, gp_XYZ& d) {
    const gp_XYZ n1 = A.Direction().XYZ();
    const gp_XYZ n2 = B.Direction().XYZ();
    d = n1.Crossed(n2);
    const double det = d.Dot(d);
    if (det < 1e-20) return false;
    const double k1 = n1.Dot(A.Location().XYZ());
    const double k2 = n2.Dot(B.Location().XYZ());
    p = (n2.Crossed(d) * k1 + d.Crossed(n1) * k2) / det;
    d /= std::sqrt(det);
    return true;
}

bool dirsParallel(const gp_XYZ& a, const gp_XYZ& b, double sinTol) {
    const double ma = a.Modulus();
    const double mb = b.Modulus();
    if (ma < kTiny || mb < kTiny) return false;
    gp_XYZ ua = a / ma;
    gp_XYZ ub = b / mb;
    if (ua.Dot(ub) < 0.0) ub.Reverse();
    return ua.Crossed(ub).Modulus() <= sinTol;
}

bool edgeParallelTo(const MeshView& mv, int e, const gp_XYZ& axisUnit, double sewTol) {
    const gp_XYZ d = edgeDir(mv, e);
    return d.Crossed(axisUnit).Modulus() <= sewTol;
}

std::vector<int> buildG2L(const MeshView& mv) {
    if (mv.nVtx == 0 || !mv.compVtx) return {};
    int mx = 0;
    for (size_t i = 0; i < mv.nVtx; i++) mx = std::max(mx, mv.compVtx[i]);
    std::vector<int> g2l(static_cast<size_t>(mx) + 1, -1);
    for (size_t i = 0; i < mv.nVtx; i++)
        g2l[static_cast<size_t>(mv.compVtx[i])] = static_cast<int>(i);
    return g2l;
}

std::vector<int> uniqueLocalVerts(const MeshView& mv, const std::vector<int>& g2l,
                                  const std::vector<int>& tris) {
    std::vector<int> lvs;
    for (int lt : tris) {
        if (lt < 0 || static_cast<size_t>(lt) >= mv.nTri) continue;
        const int gt = mv.compTris[lt];
        for (int k = 0; k < 3; k++) {
            const int gv = mv.tris[gt][k];
            if (gv < 0 || static_cast<size_t>(gv) >= g2l.size()) continue;
            const int lv = g2l[static_cast<size_t>(gv)];
            if (lv >= 0) lvs.push_back(lv);
        }
    }
    std::sort(lvs.begin(), lvs.end());
    lvs.erase(std::unique(lvs.begin(), lvs.end()), lvs.end());
    return lvs;
}

gp_XYZ localPt(const MeshView& mv, int lv) {
    return mv.pts[mv.compVtx[lv]];
}

std::vector<int> mergeCycleTris(const std::vector<Provisional>& provs,
                                const std::vector<int>& members) {
    std::vector<int> tris;
    for (int m : members) {
        if (m < 0 || static_cast<size_t>(m) >= provs.size()) continue;
        for (int t : provs[static_cast<size_t>(m)].tris) tris.push_back(t);
    }
    std::sort(tris.begin(), tris.end());
    tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
    return tris;
}

double medianOf(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

bool computeOutwardCylinder(const MeshView& mv, const std::vector<int>& tris,
                            const gp_Dir& axis, const gp_Pnt& loc) {
    const gp_XYZ aw(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ lxyz(loc.X(), loc.Y(), loc.Z());
    double sigma = 0.0;
    for (int lt : tris) {
        const double a = triAreaLocal(mv, lt);
        const gp_XYZ cent = triCentroidLocal(mv, lt);
        const gp_XYZ d = cent - lxyz;
        const double ad = d.Dot(aw);
        gp_XYZ radial = d - aw * ad;
        const double rm = radial.Modulus();
        if (rm < kTiny) continue;
        radial /= rm;
        sigma += a * triNormalLocal(mv, lt).Dot(gp_Dir(radial));
    }
    return sigma > 0.0;
}

void sortRegions(std::vector<Region>& regs) {
    std::sort(regs.begin(), regs.end(), [](const Region& a, const Region& b) {
        return minTriOf(a.tris) < minTriOf(b.tris);
    });
}

struct CycleCand {
    std::vector<int> order;
    int minTri = std::numeric_limits<int>::max();
    bool fromCapPair = false;  // both rims already known to lie in cap planes
};

void canonicalizeCycle(CycleCand& c, const std::vector<Provisional>& provs) {
    if (c.order.size() < 2) return;
    size_t best = 0;
    int bestTri = minTriId(provs[static_cast<size_t>(c.order[0])]);
    for (size_t i = 1; i < c.order.size(); i++) {
        const int mt = minTriId(provs[static_cast<size_t>(c.order[i])]);
        if (mt < bestTri) {
            bestTri = mt;
            best = i;
        }
    }
    std::rotate(c.order.begin(), c.order.begin() + static_cast<std::ptrdiff_t>(best),
                c.order.end());
    const int t1 = minTriId(provs[static_cast<size_t>(c.order[1])]);
    const int tN = minTriId(provs[static_cast<size_t>(c.order.back())]);
    if (t1 > tN)
        std::reverse(c.order.begin() + 1, c.order.end());
    c.minTri = bestTri;
}

std::vector<int> memberKey(const std::vector<int>& order) {
    std::vector<int> k = order;
    std::sort(k.begin(), k.end());
    return k;
}

void addCycle(std::vector<CycleCand>& out, std::set<std::vector<int>>& seen,
              CycleCand c, const std::vector<Provisional>& provs) {
    if (c.order.size() < static_cast<size_t>(kNMin) ||
        c.order.size() > static_cast<size_t>(kNMax))
        return;
    if (static_cast<int>(out.size()) >= kMaxCycles) return;
    canonicalizeCycle(c, provs);
    const std::vector<int> key = memberKey(c.order);
    if (!seen.insert(key).second) return;
    out.push_back(std::move(c));
}

// --- dual graph of unclaimed provisionals (any shared mesh edge) -------------

struct ProvLink {
    int other = -1;
    std::vector<int> edges;  // local mesh-edge ids, I5-sorted
};

using ProvAdj = std::vector<std::vector<ProvLink>>;

ProvAdj buildUnclaimedAdj(const MeshView& mv, const EdgeAdj& ea,
                          const std::vector<int>& triToProv, int nProv) {
    struct Raw {
        int lo, hi, e;
    };
    std::vector<Raw> raw;
    for (int e = 0; e < static_cast<int>(mv.nEdge); e++) {
        const int t0 = ea.tri[static_cast<size_t>(e)][0];
        const int t1 = ea.tri[static_cast<size_t>(e)][1];
        if (t0 < 0 || t1 < 0) continue;
        const int p0 = triToProv[static_cast<size_t>(t0)];
        const int p1 = triToProv[static_cast<size_t>(t1)];
        if (p0 < 0 || p1 < 0 || p0 == p1) continue;
        raw.push_back({ std::min(p0, p1), std::max(p0, p1), e });
    }
    std::sort(raw.begin(), raw.end(), [](const Raw& a, const Raw& b) {
        if (a.lo != b.lo) return a.lo < b.lo;
        if (a.hi != b.hi) return a.hi < b.hi;
        return a.e < b.e;
    });

    ProvAdj adj(static_cast<size_t>(nProv));
    size_t i = 0;
    while (i < raw.size()) {
        const int lo = raw[i].lo;
        const int hi = raw[i].hi;
        std::vector<int> eds;
        while (i < raw.size() && raw[i].lo == lo && raw[i].hi == hi) {
            eds.push_back(raw[i].e);
            i++;
        }
        ProvLink pa;
        pa.other = hi;
        pa.edges = eds;
        adj[static_cast<size_t>(lo)].push_back(pa);
        pa.other = lo;
        adj[static_cast<size_t>(hi)].push_back(pa);
    }
    for (int p = 0; p < nProv; p++) {
        std::sort(adj[static_cast<size_t>(p)].begin(), adj[static_cast<size_t>(p)].end(),
                  [](const ProvLink& a, const ProvLink& b) { return a.other < b.other; });
    }
    return adj;
}

const ProvLink* findLink(const ProvAdj& adj, int a, int b) {
    if (a < 0 || static_cast<size_t>(a) >= adj.size()) return nullptr;
    for (const ProvLink& L : adj[static_cast<size_t>(a)]) {
        if (L.other == b) return &L;
    }
    return nullptr;
}

bool adjacent(const ProvAdj& adj, int a, int b) {
    return findLink(adj, a, b) != nullptr;
}

// Shared mesh edges of a wall-wall link are generators iff they run parallel
// to the plane-plane intersection (n_a × n_b) within sewTol.
bool isGeneratorRidge(const MeshView& mv, const std::vector<Provisional>& provs,
                      const ProvLink& L, int self, double sewTol) {
    if (L.edges.empty()) return false;
    const Provisional& A = provs[static_cast<size_t>(self)];
    const Provisional& B = provs[static_cast<size_t>(L.other)];
    gp_XYZ p, d;
    if (!planePlaneIntersect(A.plane, B.plane, p, d)) return false;
    gp_XYZ mean(0, 0, 0);
    for (int e : L.edges) {
        if (!edgeParallelTo(mv, e, d, sewTol)) return false;
        gp_XYZ ed = edgeDir(mv, e);
        if (ed.Dot(d) < 0.0) ed.Reverse();
        const double m = ed.Modulus();
        if (m > kTiny) mean += ed / m;
    }
    return mean.SquareModulus() > kTiny;
}

using GenAdj = std::vector<std::vector<int>>;

// Cluster generator ridges whose intersection dirs are parallel, then walk
// 2-regular components. Same-axis holes land in one cluster as disjoint cycles.
void collectGeneratorCycles(const MeshView& mv, const std::vector<Provisional>& provs,
                            const ProvAdj& adj, const DerivedTols& tol, double sewTol,
                            std::vector<CycleCand>& out, std::set<std::vector<int>>& seen) {
    struct Ridge {
        int a, b;
        gp_XYZ dir;
        int minTri;
        int maxTri;
    };
    const double sinSharp = std::sin(tol.thetaSharp);
    std::vector<Ridge> ridges;
    for (size_t i = 0; i < adj.size(); i++) {
        if (!provEligible(provs[i])) continue;
        for (const ProvLink& L : adj[i]) {
            if (L.other < static_cast<int>(i)) continue;
            if (!provEligible(provs[static_cast<size_t>(L.other)])) continue;
            if (!isGeneratorRidge(mv, provs, L, static_cast<int>(i), sewTol)) continue;
            gp_XYZ p, d;
            if (!planePlaneIntersect(provs[i].plane,
                                     provs[static_cast<size_t>(L.other)].plane, p, d))
                continue;
            d = canonicalAxis(d).XYZ();
            Ridge r;
            r.a = static_cast<int>(i);
            r.b = L.other;
            r.dir = d;
            const int ta = minTriId(provs[i]);
            const int tb = minTriId(provs[static_cast<size_t>(L.other)]);
            r.minTri = std::min(ta, tb);
            r.maxTri = std::max(ta, tb);
            ridges.push_back(r);
        }
    }
    std::sort(ridges.begin(), ridges.end(), [](const Ridge& a, const Ridge& b) {
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        if (a.maxTri != b.maxTri) return a.maxTri < b.maxTri;
        if (a.a != b.a) return a.a < b.a;
        return a.b < b.b;
    });

    std::vector<int> clusterOf(ridges.size(), -1);
    std::vector<gp_XYZ> cdir;
    for (size_t i = 0; i < ridges.size(); i++) {
        gp_XYZ d = ridges[i].dir;
        int hit = -1;
        for (size_t c = 0; c < cdir.size(); c++) {
            gp_XYZ u = d;
            if (u.Dot(cdir[c]) < 0.0) u.Reverse();
            if (u.Crossed(cdir[c]).Modulus() <= sinSharp) {
                hit = static_cast<int>(c);
                break;
            }
        }
        if (hit < 0) {
            hit = static_cast<int>(cdir.size());
            cdir.push_back(d);
        }
        clusterOf[i] = hit;
    }

    for (size_t c = 0; c < cdir.size(); c++) {
        GenAdj g(adj.size());
        int nR = 0;
        for (size_t i = 0; i < ridges.size(); i++) {
            if (clusterOf[i] != static_cast<int>(c)) continue;
            nR++;
            g[static_cast<size_t>(ridges[i].a)].push_back(ridges[i].b);
            g[static_cast<size_t>(ridges[i].b)].push_back(ridges[i].a);
        }
        if (nR < kNMin) continue;
        for (size_t i = 0; i < g.size(); i++) {
            std::sort(g[i].begin(), g[i].end());
            g[i].erase(std::unique(g[i].begin(), g[i].end()), g[i].end());
        }

        std::vector<char> seenV(g.size(), 0);
        std::vector<int> starts;
        for (size_t i = 0; i < g.size(); i++) {
            if (g[i].size() >= 2 && provEligible(provs[i]))
                starts.push_back(static_cast<int>(i));
        }
        std::sort(starts.begin(), starts.end(), [&](int a, int b) {
            return minTriId(provs[static_cast<size_t>(a)]) <
                   minTriId(provs[static_cast<size_t>(b)]);
        });

        for (int v : starts) {
            if (seenV[static_cast<size_t>(v)]) continue;
            if (g[static_cast<size_t>(v)].size() != 2) continue;
            std::vector<int> cyc;
            std::vector<char> walk(g.size(), 0);
            int prev = -1;
            int cur = v;
            bool closed = false;
            while (!walk[static_cast<size_t>(cur)]) {
                walk[static_cast<size_t>(cur)] = 1;
                cyc.push_back(cur);
                const std::vector<int>& nbrs = g[static_cast<size_t>(cur)];
                int nxt = -1;
                for (int n : nbrs) {
                    if (n != prev) {
                        nxt = n;
                        break;
                    }
                }
                if (nxt < 0) break;
                if (nxt == v) {
                    closed = true;
                    break;
                }
                if (g[static_cast<size_t>(nxt)].size() != 2) break;
                prev = cur;
                cur = nxt;
            }
            if (closed && cyc.size() >= static_cast<size_t>(kNMin)) {
                for (int m : cyc) seenV[static_cast<size_t>(m)] = 1;
                CycleCand cc;
                cc.order = std::move(cyc);
                addCycle(out, seen, std::move(cc), provs);
            }
        }

        // Non-2-regular leftovers: bounded DFS from unused degree>2 nodes.
        auto dfs = [&](auto&& self, int start, int cur, int prev, std::vector<int>& path,
                       std::vector<char>& on) -> void {
            if (static_cast<int>(out.size()) >= kMaxCycles) return;
            if (path.size() >= static_cast<size_t>(kNMax)) return;
            for (int nxt : g[static_cast<size_t>(cur)]) {
                if (nxt == prev) continue;
                if (nxt == start && path.size() >= static_cast<size_t>(kNMin)) {
                    CycleCand cc;
                    cc.order = path;
                    addCycle(out, seen, std::move(cc), provs);
                    continue;
                }
                if (on[static_cast<size_t>(nxt)]) continue;
                path.push_back(nxt);
                on[static_cast<size_t>(nxt)] = 1;
                self(self, start, nxt, cur, path, on);
                on[static_cast<size_t>(nxt)] = 0;
                path.pop_back();
            }
        };
        for (int v : starts) {
            if (static_cast<int>(out.size()) >= kMaxCycles) break;
            if (seenV[static_cast<size_t>(v)]) continue;
            if (g[static_cast<size_t>(v)].size() <= 2) continue;
            std::vector<int> path = { v };
            std::vector<char> on(g.size(), 0);
            on[static_cast<size_t>(v)] = 1;
            dfs(dfs, v, v, -1, path, on);
        }
    }
}

// Through-hole: two large cap planes (area >> wall) with (anti)parallel
// normals. Axis = cap normal. Walls = unclaimed laterals adjacent to a cap;
// order by azimuth and require a closed dual cycle.
void collectCapPairCycles(const std::vector<Provisional>& provs, const ProvAdj& adj,
                          const DerivedTols& tol, std::vector<CycleCand>& out,
                          std::set<std::vector<int>>& seen) {
    const double sinSharp = std::sin(tol.thetaSharp);
    const double cosSharp = std::cos(tol.thetaSharp);
    std::vector<double> areas;
    std::vector<int> uncl;
    for (size_t i = 0; i < provs.size(); i++) {
        if (!provEligible(provs[i])) continue;
        uncl.push_back(static_cast<int>(i));
        areas.push_back(provs[i].area);
    }
    if (uncl.size() < static_cast<size_t>(kNMin + 2)) return;
    const double med = medianOf(areas);
    if (!(med > 0.0)) return;

    std::vector<int> caps;
    for (int i : uncl) {
        if (provs[static_cast<size_t>(i)].area >= kCapAreaMul * med)
            caps.push_back(i);
    }
    std::sort(caps.begin(), caps.end(), [&](int a, int b) {
        if (provs[static_cast<size_t>(a)].area != provs[static_cast<size_t>(b)].area)
            return provs[static_cast<size_t>(a)].area > provs[static_cast<size_t>(b)].area;
        return minTriId(provs[static_cast<size_t>(a)]) < minTriId(provs[static_cast<size_t>(b)]);
    });

    auto wallsNextTo = [&](int cap, const gp_Dir& axis) {
        std::vector<int> w;
        for (const ProvLink& L : adj[static_cast<size_t>(cap)]) {
            if (!provEligible(provs[static_cast<size_t>(L.other)])) continue;
            if (std::abs(provs[static_cast<size_t>(L.other)].plane.Direction().Dot(axis)) >=
                sinSharp)
                continue;
            if (provs[static_cast<size_t>(L.other)].area >=
                provs[static_cast<size_t>(cap)].area / kCapAreaMul)
                continue;
            w.push_back(L.other);
        }
        std::sort(w.begin(), w.end());
        w.erase(std::unique(w.begin(), w.end()), w.end());
        return w;
    };

    for (size_t i = 0; i < caps.size(); i++) {
        for (size_t j = i + 1; j < caps.size(); j++) {
            const gp_Dir n0 = provs[static_cast<size_t>(caps[i])].plane.Direction();
            const gp_Dir n1 = provs[static_cast<size_t>(caps[j])].plane.Direction();
            if (std::abs(n0.Dot(n1)) < cosSharp) continue;
            const gp_Dir axis = canonicalAxis(n0.XYZ());

            std::vector<int> wA = wallsNextTo(caps[i], axis);
            std::vector<int> wB = wallsNextTo(caps[j], axis);
            std::vector<int> both;
            std::set_intersection(wA.begin(), wA.end(), wB.begin(), wB.end(),
                                  std::back_inserter(both));
            if (static_cast<int>(both.size()) < kNMin) continue;

            std::vector<char> isW(provs.size(), 0);
            for (int w : both) isW[static_cast<size_t>(w)] = 1;
            std::vector<int> rest = both;
            std::vector<char> vis(provs.size(), 0);
            while (!rest.empty()) {
                std::sort(rest.begin(), rest.end(), [&](int a, int b) {
                    return minTriId(provs[static_cast<size_t>(a)]) <
                           minTriId(provs[static_cast<size_t>(b)]);
                });
                int seed = -1;
                for (int r : rest) {
                    if (!vis[static_cast<size_t>(r)]) {
                        seed = r;
                        break;
                    }
                }
                if (seed < 0) break;
                std::vector<int> comp;
                std::vector<int> st = { seed };
                vis[static_cast<size_t>(seed)] = 1;
                while (!st.empty()) {
                    const int u = st.back();
                    st.pop_back();
                    comp.push_back(u);
                    for (const ProvLink& L : adj[static_cast<size_t>(u)]) {
                        if (!isW[static_cast<size_t>(L.other)]) continue;
                        if (vis[static_cast<size_t>(L.other)]) continue;
                        vis[static_cast<size_t>(L.other)] = 1;
                        st.push_back(L.other);
                    }
                }
                rest.erase(std::remove_if(rest.begin(), rest.end(),
                                          [&](int x) { return vis[static_cast<size_t>(x)] != 0; }),
                           rest.end());
                if (static_cast<int>(comp.size()) < kNMin) continue;

                struct Az {
                    int idx;
                    double psi;
                    int minTri;
                };
                std::vector<Az> az;
                az.reserve(comp.size());
                gp_XYZ tmp = (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
                const gp_XYZ xD = axis.XYZ().Crossed(tmp);
                const gp_Dir xDir = normalizeXYZ(xD);
                const gp_XYZ yD = axis.XYZ().Crossed(xDir.XYZ());
                const gp_Dir yDir = normalizeXYZ(yD);
                for (int w : comp) {
                    const gp_XYZ n = provs[static_cast<size_t>(w)].plane.Direction().XYZ();
                    Az z;
                    z.idx = w;
                    z.psi = std::atan2(n.Dot(yDir.XYZ()), n.Dot(xDir.XYZ()));
                    z.minTri = minTriId(provs[static_cast<size_t>(w)]);
                    az.push_back(z);
                }
                std::sort(az.begin(), az.end(), [](const Az& a, const Az& b) {
                    if (a.psi != b.psi) return a.psi < b.psi;
                    return a.minTri < b.minTri;
                });
                bool ring = true;
                for (size_t k = 0; k < az.size(); k++) {
                    const int a = az[k].idx;
                    const int b = az[(k + 1) % az.size()].idx;
                    if (!adjacent(adj, a, b)) {
                        ring = false;
                        break;
                    }
                }
                if (!ring) continue;
                CycleCand cc;
                cc.fromCapPair = true;
                cc.order.reserve(az.size());
                for (const Az& z : az) cc.order.push_back(z.idx);
                addCycle(out, seen, std::move(cc), provs);
            }
        }
    }
}

void fillNgonCylinder(const MeshView& mv, const std::vector<int>& tris,
                      const gp_Dir& axis, const gp_Pnt& center, double radius, int nSides,
                      Region& reg) {
    double areaReg = 0.0;
    for (int lt : tris) areaReg += triAreaLocal(mv, lt);

    const gp_XYZ aw(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ cxyz(center.X(), center.Y(), center.Z());

    gp_XYZ centroid(0, 0, 0);
    double totalArea = 0.0;
    for (int lt : tris) {
        const double ar = triAreaLocal(mv, lt);
        centroid += triCentroidLocal(mv, lt) * ar;
        totalArea += ar;
    }
    if (totalArea > kTiny) centroid /= totalArea;
    const double tLoc = (centroid - cxyz).Dot(aw);
    const gp_XYZ loc = cxyz + aw * tLoc;
    const gp_Pnt Loc(loc);

    const std::vector<int> g2l = buildG2L(mv);
    const std::vector<int> lvs = uniqueLocalVerts(mv, g2l, tris);
    int minLv = std::numeric_limits<int>::max();
    gp_XYZ ps(0, 0, 0);
    for (int lv : lvs) {
        if (lv < minLv) {
            minLv = lv;
            ps = localPt(mv, lv);
        }
    }
    const gp_XYZ dps = ps - loc;
    const double ads = dps.Dot(aw);
    gp_XYZ xDir = dps - aw * ads;
    if (xDir.Modulus() < 1e-12) {
        gp_XYZ tmp = (std::abs(aw.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
        xDir = aw.Crossed(tmp);
    }

    // Field assignments match fillCylinderRegion (refit_grow.cpp); origin and
    // closed360 are the N-gon overrides. Cannot call that static helper.
    reg.type = SurfType::Cylinder;
    reg.origin = Origin::NgonWall;
    reg.tris = tris;
    reg.ax = gp_Ax3(Loc, axis, normalizeXYZ(xDir));
    reg.radius = radius;
    reg.uMin = 0.0;
    reg.uMax = kTwoPi;
    double vMin = std::numeric_limits<double>::infinity();
    double vMax = -std::numeric_limits<double>::infinity();
    for (int lv : lvs) {
        const double v = (localPt(mv, lv) - loc).Dot(aw);
        vMin = std::min(vMin, v);
        vMax = std::max(vMax, v);
    }
    if (!std::isfinite(vMin)) {
        vMin = 0.0;
        vMax = 0.0;
    }
    reg.vMin = vMin;
    reg.vMax = vMax;
    reg.closed360 = true;
    reg.nSides = nSides;
    reg.chordSagitta = chordSagitta(radius, nSides);
    reg.outwardNormal = computeOutwardCylinder(mv, tris, axis, Loc);
    reg.dVolPredicted = dVolCylinderSector(areaReg, radius, nSides, reg.outwardNormal);

    double maxDev = 0.0;
    double sumSq = 0.0;
    int nV = 0;
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            const gp_XYZ p = localTriVert(mv, lt, k);
            const double rr = aw.Crossed(p - loc).Modulus();
            const double d = rr - radius;
            maxDev = std::max(maxDev, std::abs(d));
            sumSq += d * d;
            nV++;
        }
    }
    reg.maxVertexDev = maxDev;
    reg.rmsVertexDev = nV > 0 ? std::sqrt(sumSq / static_cast<double>(nV)) : 0.0;
}

// Reject: cycle shares a generator ridge with an outside wall whose
// plane-plane intersection is not ∥ the cycle axis (turning-axis blend).
bool sharesAxisTurnBlend(const MeshView& mv, const std::vector<Provisional>& provs,
                         const ProvAdj& adj, const CycleCand& cyc, const gp_XYZ& axyz,
                         double sewTol, double sinSharp) {
    std::vector<char> inCyc(provs.size(), 0);
    for (int m : cyc.order) {
        if (m >= 0 && static_cast<size_t>(m) < inCyc.size())
            inCyc[static_cast<size_t>(m)] = 1;
    }
    for (int w : cyc.order) {
        if (w < 0 || static_cast<size_t>(w) >= adj.size()) continue;
        for (const ProvLink& L : adj[static_cast<size_t>(w)]) {
            if (L.other < 0 || static_cast<size_t>(L.other) >= provs.size()) continue;
            if (inCyc[static_cast<size_t>(L.other)]) continue;
            if (!provEligible(provs[static_cast<size_t>(L.other)])) continue;
            if (!isGeneratorRidge(mv, provs, L, w, sewTol)) continue;
            gp_XYZ p, d;
            if (!planePlaneIntersect(provs[static_cast<size_t>(w)].plane,
                                     provs[static_cast<size_t>(L.other)].plane, p, d))
                continue;
            if (!dirsParallel(d, axyz, sinSharp)) return true;
        }
    }
    return false;
}

// Cap-plane rims at the two axial ends. A through-hole / outer boss has both
// ends in cap planes; a chamfered through-hole has one cap and a slope ring.
struct CycleEnds {
    int nCapEnds = 0;
    bool otherEndChamfer = false;
};

CycleEnds classifyCycleEnds(const MeshView& mv, const std::vector<Provisional>& provs,
                            const ProvAdj& adj, const CycleCand& cyc, const gp_Dir& axis,
                            const std::vector<int>& lvs, double sewTol, double cosSharp) {
    CycleEnds out;
    const int N = static_cast<int>(cyc.order.size());
    if (N <= 0 || lvs.empty()) return out;

    const gp_XYZ axyz(axis.X(), axis.Y(), axis.Z());
    std::vector<char> inCyc(provs.size(), 0);
    for (int m : cyc.order) {
        if (m >= 0 && static_cast<size_t>(m) < inCyc.size())
            inCyc[static_cast<size_t>(m)] = 1;
    }

    double vMin = std::numeric_limits<double>::infinity();
    double vMax = -std::numeric_limits<double>::infinity();
    for (int lv : lvs) {
        const double v = localPt(mv, lv).Dot(axyz);
        vMin = std::min(vMin, v);
        vMax = std::max(vMax, v);
    }
    if (!std::isfinite(vMin) || !(vMax > vMin)) return out;
    const double mid = 0.5 * (vMin + vMax);

    int capLo = 0, capHi = 0, chamLo = 0, chamHi = 0;
    for (int w : cyc.order) {
        if (w < 0 || static_cast<size_t>(w) >= adj.size()) continue;
        bool wCapLo = false, wCapHi = false, wChamLo = false, wChamHi = false;
        for (const ProvLink& L : adj[static_cast<size_t>(w)]) {
            if (L.other < 0 || static_cast<size_t>(L.other) >= provs.size()) continue;
            if (inCyc[static_cast<size_t>(L.other)]) continue;
            if (!provEligible(provs[static_cast<size_t>(L.other)])) continue;
            const double nd =
                std::abs(provs[static_cast<size_t>(L.other)].plane.Direction().Dot(axis));
            int nRim = 0;
            double vSum = 0.0;
            for (int e : L.edges) {
                if (e < 0 || static_cast<size_t>(e) >= mv.nEdge) continue;
                if (!mv.compEdges) continue;
                if (edgeParallelTo(mv, e, axyz, sewTol)) continue;
                const gp_XYZ midE = 0.5 * (localPt(mv, mv.compEdges[e].first) +
                                           localPt(mv, mv.compEdges[e].second));
                vSum += midE.Dot(axyz);
                nRim++;
            }
            if (nRim == 0) continue;
            const double v = vSum / static_cast<double>(nRim);
            const bool lo = v < mid;
            if (nd >= cosSharp) {
                if (lo) wCapLo = true;
                else wCapHi = true;
            } else if (nd >= kSinChamLo && nd <= kSinChamHi) {
                if (lo) wChamLo = true;
                else wChamHi = true;
            }
        }
        if (wCapLo) capLo++;
        if (wCapHi) capHi++;
        if (wChamLo) chamLo++;
        if (wChamHi) chamHi++;
    }

    const int need = (N + 1) / 2;
    const bool loCap = capLo >= need;
    const bool hiCap = capHi >= need;
    const bool loCham = chamLo >= need;
    const bool hiCham = chamHi >= need;
    out.nCapEnds = (loCap ? 1 : 0) + (hiCap ? 1 : 0);
    out.otherEndChamfer = (loCap && hiCham) || (hiCap && loCham);
    return out;
}

// Through-hole or outer boss only. Two cap-plane rims, or (N>8) one cap plus
// a chamfer ring. N<=8 without both caps is a 6-facet bump, not a hole.
bool isThroughHoleOrOuterBoss(const CycleEnds& ends, int N, bool fromCapPair) {
    if (fromCapPair) return true;  // collectCapPairCycles already required both caps
    if (ends.nCapEnds == 2) return true;
    if (N > 8 && ends.nCapEnds == 1 && ends.otherEndChamfer) return true;
    return false;
}

bool tryAcceptCycle(const MeshView& mv, const DerivedTols& tol, SegmentWork& work,
                    const ProvAdj& adj, const CycleCand& cyc, const std::vector<int>& g2l,
                    const char** why = nullptr, double* devOut = nullptr,
                    double* tolOut = nullptr, double* rOut = nullptr) {
    const char* sink = nullptr;
    if (!why) why = &sink;
    *why = "accepted";
    auto no = [&](const char* r) {
        *why = r;
        return false;
    };
    const std::vector<Provisional>& provs = work.provisionals;
    const int N = static_cast<int>(cyc.order.size());
    if (N < kNMin || N > kNMax) return no("N-out-of-range");
    for (int m : cyc.order) {
        if (m < 0 || static_cast<size_t>(m) >= provs.size()) return no("bad-member");
        if (!provEligible(provs[static_cast<size_t>(m)])) return no("member-claimed");
    }

    const double sinSharp = std::sin(tol.thetaSharp);
    const double sewTol = mv.sewTol;

    gp_XYZ axisSum(0, 0, 0);
    std::vector<gp_XYZ> genDirs;
    genDirs.reserve(static_cast<size_t>(N));
    for (int i = 0; i < N; i++) {
        const int a = cyc.order[static_cast<size_t>(i)];
        const int b = cyc.order[static_cast<size_t>((i + 1) % N)];
        gp_XYZ p, d;
        if (!planePlaneIntersect(provs[static_cast<size_t>(a)].plane,
                                 provs[static_cast<size_t>(b)].plane, p, d))
            return no("gen-plane-parallel");
        genDirs.push_back(d);
        if (axisSum.SquareModulus() > kTiny && d.Dot(axisSum) < 0.0) d.Reverse();
        axisSum += d;
    }

    gp_Dir axis;
    bool haveAxis = axisSum.SquareModulus() > kTiny;
    if (haveAxis) axis = canonicalAxis(axisSum);

    // Fallback: two cap planes with area >> wall, normals ~ axis.
    if (!haveAxis) {
        const double meanWall = [&]() {
            double s = 0.0;
            for (int m : cyc.order) s += provs[static_cast<size_t>(m)].area;
            return s / static_cast<double>(N);
        }();
        std::vector<int> caps;
        std::vector<char> inCyc(provs.size(), 0);
        for (int m : cyc.order) inCyc[static_cast<size_t>(m)] = 1;
        for (size_t pi = 0; pi < provs.size(); pi++) {
            if (!provEligible(provs[pi]) || inCyc[pi]) continue;
            if (!(provs[pi].area >= kCapAreaMul * meanWall)) continue;
            bool nbr = false;
            for (int m : cyc.order) {
                if (adjacent(adj, static_cast<int>(pi), m)) {
                    nbr = true;
                    break;
                }
            }
            if (nbr) caps.push_back(static_cast<int>(pi));
        }
        std::sort(caps.begin(), caps.end(), [&](int a, int b) {
            if (provs[static_cast<size_t>(a)].area != provs[static_cast<size_t>(b)].area)
                return provs[static_cast<size_t>(a)].area > provs[static_cast<size_t>(b)].area;
            return minTriId(provs[static_cast<size_t>(a)]) < minTriId(provs[static_cast<size_t>(b)]);
        });
        if (caps.size() < 2) return no("axis-fallback-no-caps");
        const gp_Dir n0 = provs[static_cast<size_t>(caps[0])].plane.Direction();
        const gp_Dir n1 = provs[static_cast<size_t>(caps[1])].plane.Direction();
        if (std::abs(n0.Dot(n1)) < std::cos(tol.thetaSharp)) return no("axis-fallback-caps-skew");
        axis = canonicalAxis(n0.XYZ());
        haveAxis = true;
    }
    if (!haveAxis) return no("no-axis");

    const gp_XYZ axyz(axis.X(), axis.Y(), axis.Z());

    // Skip turning-axis blends: successive generators must stay ∥ mean axis.
    for (const gp_XYZ& d : genDirs) {
        if (!dirsParallel(d, axyz, sinSharp)) return no("gen-not-parallel-axis");
    }
    // Reject cycles that share an axis-turn (blend) with an outside wall.
    if (sharesAxisTurnBlend(mv, provs, adj, cyc, axyz, sewTol, sinSharp))
        return no("axis-turn-blend");

    // Skip cones: |n_i · axis| must be < sin(thetaSharp) on every wall.
    for (int m : cyc.order) {
        const gp_Dir n = provs[static_cast<size_t>(m)].plane.Direction();
        if (std::abs(n.Dot(axis)) >= sinSharp) return no("wall-is-cone");
    }

    // Lateral (shared wall-wall) edges ∥ axis within sewTol.
    for (int i = 0; i < N; i++) {
        const int a = cyc.order[static_cast<size_t>(i)];
        const int b = cyc.order[static_cast<size_t>((i + 1) % N)];
        const ProvLink* L = findLink(adj, a, b);
        if (!L || L->edges.empty()) return no("no-lateral-link");
        for (int e : L->edges) {
            if (!edgeParallelTo(mv, e, axyz, sewTol)) return no("lateral-edge-not-parallel");
        }
    }

    const std::vector<int> tris = mergeCycleTris(provs, cyc.order);
    if (tris.empty()) return no("no-tris");
    const std::vector<int> lvs = uniqueLocalVerts(mv, g2l, tris);
    if (lvs.size() < 3) return no("no-verts");

    const double cosSharp = std::cos(tol.thetaSharp);
    const CycleEnds ends = classifyCycleEnds(mv, provs, adj, cyc, axis, lvs, sewTol, cosSharp);
    if (!isThroughHoleOrOuterBoss(ends, N, cyc.fromCapPair)) {
        if (ngonDiagOn())
            std::fprintf(stderr, "  DIAG_A_ENDS N=%d capEnds=%d otherEndChamfer=%d\n", N,
                         ends.nCapEnds, ends.otherEndChamfer ? 1 : 0);
        return no("not-through-hole-or-boss");
    }

    // Axis point = I5 centroid of unique verts in the plane ⊥ axis.
    gp_XYZ sumPerp(0, 0, 0);
    for (int lv : lvs) {
        const gp_XYZ p = localPt(mv, lv);
        sumPerp += p - axyz * p.Dot(axyz);
    }
    const gp_XYZ locPerp = sumPerp / static_cast<double>(lvs.size());

    // R = mean vertex distance to axis (circumradius), I5 vertex-id order.
    double sumR = 0.0;
    for (int lv : lvs) {
        const gp_XYZ p = localPt(mv, lv);
        sumR += axyz.Crossed(p - locPerp).Modulus();
    }
    const double R = sumR / static_cast<double>(lvs.size());
    if (rOut) *rOut = R;
    if (!(R > 0.0) || R >= 2.0 * mv.diag) return no("bad-radius");

    double maxDev = 0.0;
    for (int lv : lvs) {
        const gp_XYZ p = localPt(mv, lv);
        const double rho = axyz.Crossed(p - locPerp).Modulus();
        maxDev = std::max(maxDev, std::abs(rho - R));
    }
    const double acceptTol = std::max(sewTol, chordSagitta(R, N));
    if (devOut) *devOut = maxDev;
    if (tolOut) *tolOut = acceptTol;
    if (maxDev > acceptTol) return no("residual-over-accept-tol");

    // closed360: vertex azimuths around the axis, D2-style max-gap test.
    gp_XYZ x0 = (std::abs(axyz.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    gp_XYZ xD = axyz.Crossed(x0);
    if (xD.Modulus() < kTiny) return no("degenerate-basis");
    xD /= xD.Modulus();
    const gp_XYZ yD = axyz.Crossed(xD);
    std::vector<double> chi;
    chi.reserve(lvs.size());
    for (int lv : lvs) {
        const gp_XYZ d = localPt(mv, lv) - locPerp;
        chi.push_back(std::atan2(d.Dot(yD), d.Dot(xD)));
    }
    std::sort(chi.begin(), chi.end());
    double maxGap = 0.0;
    for (size_t i = 1; i < chi.size(); i++)
        maxGap = std::max(maxGap, chi[i] - chi[i - 1]);
    if (!chi.empty())
        maxGap = std::max(maxGap, kTwoPi - chi.back() + chi.front());
    const double bandArc = kTwoPi / static_cast<double>(N);
    if (!(maxGap <= 1.5 * bandArc)) return no("not-closed360");

    const gp_Pnt center(locPerp.X(), locPerp.Y(), locPerp.Z());
    Region reg;
    fillNgonCylinder(mv, tris, axis, center, R, N, reg);
    work.accepted.push_back(reg);
    for (int m : cyc.order) {
        work.provisionals[static_cast<size_t>(m)].claim = ProvClaim::ConsumedCylinder;
        work.provisionals[static_cast<size_t>(m)].tris.clear();
    }
    return true;
}

}  // namespace

bool claimNgonWallsA(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work) {
    (void)p;
    if (work.provisionals.empty() || mv.nTri == 0) return true;

    const int nProv = static_cast<int>(work.provisionals.size());
    std::vector<int> triToProv(mv.nTri, -1);
    for (int pi = 0; pi < nProv; pi++) {
        if (!provEligible(work.provisionals[static_cast<size_t>(pi)])) continue;
        for (int t : work.provisionals[static_cast<size_t>(pi)].tris) {
            if (t >= 0 && static_cast<size_t>(t) < mv.nTri)
                triToProv[static_cast<size_t>(t)] = pi;
        }
    }

    const EdgeAdj ea = buildEdgeAdj(mv);
    const ProvAdj adj = buildUnclaimedAdj(mv, ea, triToProv, nProv);
    const std::vector<int> g2l = buildG2L(mv);

    std::vector<CycleCand> cands;
    std::set<std::vector<int>> seen;
    collectCapPairCycles(work.provisionals, adj, tol, cands, seen);
    const size_t nCapPairCands = cands.size();
    collectGeneratorCycles(mv, work.provisionals, adj, tol, mv.sewTol, cands, seen);
    if (ngonDiagOn()) {
        int nElig = 0;
        int byClaim[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        int emptyTris = 0;
        for (const Provisional& pr : work.provisionals) {
            if (provEligible(pr)) nElig++;
            byClaim[static_cast<int>(pr.claim) & 7]++;
            if (pr.tris.empty()) emptyTris++;
        }
        std::fprintf(stderr,
                     "DIAG_A_CANDS nProv=%d eligible=%d emptyTris=%d claim[U/inCyl/cyl/"
                     "inFil/fil/plane]=%d/%d/%d/%d/%d/%d cands=%zu capPair=%zu gen=%zu\n",
                     nProv, nElig, emptyTris, byClaim[0], byClaim[1], byClaim[2],
                     byClaim[3], byClaim[4], byClaim[5], cands.size(), nCapPairCands,
                     cands.size() - nCapPairCands);
        // Per-provisional ridge degree in the unclaimed dual graph: a ring of a
        // closed N-gon bore has exactly two generator ridges per wall.
        std::vector<int> degHist(8, 0);
        int nRidge = 0;
        for (size_t i = 0; i < adj.size(); i++) {
            if (!provEligible(work.provisionals[i])) continue;
            int deg = 0;
            for (const ProvLink& L : adj[i]) {
                if (!provEligible(work.provisionals[static_cast<size_t>(L.other)])) continue;
                if (isGeneratorRidge(mv, work.provisionals, L, static_cast<int>(i), mv.sewTol))
                    deg++;
            }
            nRidge += deg;
            degHist[static_cast<size_t>(std::min(deg, 7))]++;
        }
        std::fprintf(stderr,
                     "DIAG_A_RIDGEDEG halfEdges=%d deg0=%d deg1=%d deg2=%d deg3=%d "
                     "deg4=%d deg5=%d deg6=%d deg7+=%d\n",
                     nRidge, degHist[0], degHist[1], degHist[2], degHist[3], degHist[4],
                     degHist[5], degHist[6], degHist[7]);
        for (const Region& r : work.accepted)
            std::fprintf(stderr,
                         "DIAG_A_PRIOR origin=%d type=%d R=%.5f nTri=%zu nSides=%d "
                         "closed360=%d minTri=%d\n",
                         static_cast<int>(r.origin), static_cast<int>(r.type), r.radius,
                         r.tris.size(), r.nSides, r.closed360 ? 1 : 0, minTriOf(r.tris));
    }

    // Prefer holes whose two rims lie in cap planes, then larger N (real
    // bores over 6-facet bumps), then I5 min-tri / member order.
    std::sort(cands.begin(), cands.end(), [](const CycleCand& a, const CycleCand& b) {
        if (a.fromCapPair != b.fromCapPair) return a.fromCapPair && !b.fromCapPair;
        if (a.order.size() != b.order.size()) return a.order.size() > b.order.size();
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        return a.order < b.order;
    });

    for (const CycleCand& c : cands) {
        const char* why = nullptr;
        double dev = -1.0, atol = -1.0, R = -1.0;
        const bool ok = tryAcceptCycle(mv, tol, work, adj, c, g2l, &why, &dev, &atol, &R);
        if (ngonDiagOn())
            std::fprintf(stderr,
                         "DIAG_A_CYC N=%zu minTri=%d capPair=%d R=%.6f dev=%.6g "
                         "acceptTol=%.6g accepted=%d why=%s\n",
                         c.order.size(), c.minTri, c.fromCapPair ? 1 : 0, R, dev, atol,
                         ok ? 1 : 0, why ? why : "?");
    }

    sortRegions(work.accepted);
    return true;
}

}}  // namespace stl2step::refit
