// stl2step P1 stage C1 — fillet strip claim (D1 §1.5, D7 radius).
// gp_ / stdlib only. Never throws. Silent on every reject / NYI (D5.5).
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kCoplanarG = 1e-6;
constexpr double kDeg30 = 30.0;
constexpr double kDeg180 = 180.0;

struct Nbr {
    enum Kind : uint8_t { None = 0, Prov = 1, Acc = 2 };
    Kind kind = None;
    int  idx  = -1;   // provisional index, or work.accepted index for Acc.
                      // Stage D remaps filletNbrA/B from these to Region::id.

    bool operator==(const Nbr& o) const { return kind == o.kind && idx == o.idx; }
    bool operator<(const Nbr& o) const {
        if (kind != o.kind) return kind < o.kind;
        return idx < o.idx;
    }
    bool valid() const { return kind != None && idx >= 0; }
};

struct NbrTally {
    Nbr n;
    int count = 0;
};

double medianInPlace(std::vector<double>& v) {
    const std::size_t n = v.size();
    if (n == 0) return 0.0;
    std::sort(v.begin(), v.end());
    if (n & 1u) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
}

int minTriOfProv(const Provisional& p) {
    return minTriOf(p.tris);
}

double triAreaNormal(const MeshView& mv, int t, gp_XYZ& nOut) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    const gp_XYZ a = mv.pts[mv.tris[gt][0]];
    const gp_XYZ b = mv.pts[mv.tris[gt][1]];
    const gp_XYZ c = mv.pts[mv.tris[gt][2]];
    const gp_XYZ n = (b - a).Crossed(c - a);
    const double mag = n.Modulus();
    if (mag > 0.0) nOut = n / mag;
    else nOut = gp_XYZ(0, 0, 1);
    return 0.5 * mag;
}

double triArea(const MeshView& mv, int t) {
    gp_XYZ n;
    return triAreaNormal(mv, t, n);
}

double regionArea(const MeshView& mv, const std::vector<int>& tris) {
    std::vector<int> ord = tris;
    std::sort(ord.begin(), ord.end());
    double a = 0.0;
    for (int t : ord) a += triArea(mv, t);
    return a;
}

gp_XYZ triCentroid(const MeshView& mv, int t) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return (mv.pts[mv.tris[gt][0]] + mv.pts[mv.tris[gt][1]] + mv.pts[mv.tris[gt][2]])
           / 3.0;
}

// local vertex id of triangle corner k; MeshView.tris stores global point ids.
int cornerLocal(const std::vector<int>& g2l, const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    const int g = mv.tris[gt][k];
    if (g < 0 || g >= (int)g2l.size()) return -1;
    return g2l[g];
}

gp_XYZ localPt(const MeshView& mv, int localV) {
    return mv.pts[mv.compVtx ? mv.compVtx[localV] : localV];
}

int edgeOtherTri(const std::vector<std::array<int, 2>>& edgeTris, int e, int t) {
    if (e < 0 || e >= (int)edgeTris.size()) return -1;
    const int a = edgeTris[e][0], b = edgeTris[e][1];
    if (a == t) return b;
    if (b == t) return a;
    return -1;
}

void tallyAdd(std::vector<NbrTally>& tallies, const Nbr& n) {
    for (NbrTally& x : tallies) {
        if (x.n == n) {
            x.count++;
            return;
        }
    }
    NbrTally x;
    x.n = n;
    x.count = 1;
    tallies.push_back(x);
}

double d7Radius(double dL, double dR, double g, bool& ok) {
    ok = false;
    const double omg = 1.0 - g;
    if (!(omg >= kCoplanarG)) return 0.0;
    if (dL < 0.0) dL = 0.0;
    if (dR < 0.0) dR = 0.0;
    double disc = 2.0 * (1.0 + g) * dL * dR;
    if (disc < 0.0) disc = 0.0;
    const double Ri = (dL + dR + std::sqrt(disc)) / omg;
    if (!std::isfinite(Ri) || Ri < 0.0) return 0.0;
    ok = true;
    return Ri;
}

// D2.2 facet-normal clustering. Ported from p1-fillet-r3 (GRAFT-2).
// psi = atan2(n·v0, n·u0) in (-pi, pi]; thetaBin in radians.
int nBandsFromNormals(const std::vector<double>& psi, double thetaBin) {
    if (psi.empty()) return 0;
    std::vector<double> s = psi;
    std::sort(s.begin(), s.end());
    const int n = (int)s.size();
    if (n == 1) return 1;
    std::vector<double> gaps((std::size_t)n);
    for (int i = 0; i < n - 1; ++i)
        gaps[(std::size_t)i] = s[(std::size_t)i + 1] - s[(std::size_t)i];
    gaps[(std::size_t)n - 1] = (s[0] + 2.0 * kPi) - s[(std::size_t)n - 1];
    int imax = 0;
    for (int i = 1; i < n; ++i)
        if (gaps[(std::size_t)i] > gaps[(std::size_t)imax]) imax = i;
    std::vector<double> rest;
    rest.reserve((std::size_t)n - 1);
    for (int i = 0; i < n; ++i)
        if (i != imax) rest.push_back(gaps[(std::size_t)i]);
    double p75 = 0.0;
    if (!rest.empty()) {
        std::sort(rest.begin(), rest.end());
        const std::size_t idx = (rest.size() * 3u) / 4u;
        p75 = rest[idx < rest.size() ? idx : rest.size() - 1u];
    }
    const double th = std::max(thetaBin, 0.5 * p75);
    const int start = (imax + 1) % n;
    int bands = 1;
    double prev = s[(std::size_t)start];
    for (int k = 1; k < n; ++k) {
        const int idx = (start + k) % n;
        const double cur = s[(std::size_t)idx];
        double gap = cur - prev;
        if (gap < 0.0) gap += 2.0 * kPi;
        if (gap > th) ++bands;
        prev = cur;
    }
    return bands;
}

// GRAFT-1: sign a band normal to the wedge bisector (mL+mR), not mL alone.
// dot(mL) flips any band more than 90° from L and drops fillets with span
// ≳108° (~40% of D7.3's sanctioned [30°, 180°] window).
void signTowardBisector(gp_XYZ& nn, const gp_XYZ& mL, const gp_XYZ& mR) {
    if (nn.Dot(mL + mR) < 0.0) nn.Reverse();
}

gp_XYZ planeNormal(const gp_Ax3& ax) {
    return ax.Direction().XYZ();
}

gp_XYZ planePoint(const gp_Ax3& ax) {
    return ax.Location().XYZ();
}

bool planeOfNbr(const Nbr& n, const SegmentWork& work, gp_Ax3& ax) {
    if (n.kind == Nbr::Prov) {
        if (n.idx < 0 || n.idx >= (int)work.provisionals.size()) return false;
        ax = work.provisionals[n.idx].plane;
        return true;
    }
    if (n.kind == Nbr::Acc) {
        if (n.idx < 0 || n.idx >= (int)work.accepted.size()) return false;
        ax = work.accepted[n.idx].ax;
        return true;
    }
    return false;
}

bool nbrIsCylinder(const Nbr& n, const SegmentWork& work) {
    if (n.kind == Nbr::Acc) {
        if (n.idx < 0 || n.idx >= (int)work.accepted.size()) return false;
        const Region& r = work.accepted[n.idx];
        // ChamferCone is not a FilletStrip and not a cylinder flank.
        if (r.origin == Origin::ChamferCone) return false;
        // NgonWall is a committed cylinder, same as CylGrow.
        return r.type == SurfType::Cylinder;
    }
    if (n.kind == Nbr::Prov) {
        if (n.idx < 0 || n.idx >= (int)work.provisionals.size()) return false;
        const ProvClaim c = work.provisionals[n.idx].claim;
        // ConsumedCylinder (B1 / NgonWall / ChamferCone) is ineligible as a
        // C1 member (Unclaimed-only). As a neighbour it is cylinder-like, so
        // the strip becomes TorusNYI rather than a stolen FilletStrip.
        return c == ProvClaim::ConsumedCylinder || c == ProvClaim::InCylinderClaim;
    }
    return false;
}

bool nbrIsPlaneLike(const Nbr& n, const SegmentWork& work) {
    if (n.kind == Nbr::Prov) {
        if (n.idx < 0 || n.idx >= (int)work.provisionals.size()) return false;
        return work.provisionals[n.idx].claim == ProvClaim::Unclaimed;
    }
    if (n.kind == Nbr::Acc) {
        if (n.idx < 0 || n.idx >= (int)work.accepted.size()) return false;
        const Origin o = work.accepted[n.idx].origin;
        if (o == Origin::NgonWall || o == Origin::CylGrow || o == Origin::ChamferCone)
            return false;
        return work.accepted[n.idx].type == SurfType::Plane;
    }
    return false;
}

gp_XYZ stripExtrusionAxis(const MeshView& mv, const std::vector<int>& g2l,
                          const std::vector<int>& tris) {
    gp_XYZ mn(1e300, 1e300, 1e300), mx(-1e300, -1e300, -1e300);
    for (int t : tris) {
        for (int k = 0; k < 3; ++k) {
            const int lv = cornerLocal(g2l, mv, t, k);
            if (lv < 0) continue;
            const gp_XYZ p = localPt(mv, lv);
            mn.SetCoord(std::min(mn.X(), p.X()), std::min(mn.Y(), p.Y()),
                        std::min(mn.Z(), p.Z()));
            mx.SetCoord(std::max(mx.X(), p.X()), std::max(mx.Y(), p.Y()),
                        std::max(mx.Z(), p.Z()));
        }
    }
    const gp_XYZ span(mx.X() - mn.X(), mx.Y() - mn.Y(), mx.Z() - mn.Z());
    if (span.X() >= span.Y() && span.X() >= span.Z() && span.X() > 1e-12)
        return gp_XYZ(1, 0, 0);
    if (span.Y() >= span.Z() && span.Y() > 1e-12) return gp_XYZ(0, 1, 0);
    if (span.Z() > 1e-12) return gp_XYZ(0, 0, 1);
    return gp_XYZ(1, 0, 0);
}

bool provIsExtrusionCap(const MeshView& mv, const std::vector<int>& g2l,
                        const Provisional& pr) {
    gp_XYZ n = planeNormal(pr.plane);
    const double nm = n.Modulus();
    if (nm < 1e-15) return false;
    n /= nm;
    const double ax = std::abs(n.X());
    const double ay = std::abs(n.Y());
    const double az = std::abs(n.Z());
    const double am = std::max(ax, std::max(ay, az));
    if (am < 0.95) return false;

    gp_XYZ mn(1e300, 1e300, 1e300), mx(-1e300, -1e300, -1e300);
    for (int t : pr.tris) {
        for (int k = 0; k < 3; ++k) {
            const int lv = cornerLocal(g2l, mv, t, k);
            if (lv < 0) continue;
            const gp_XYZ p = localPt(mv, lv);
            mn.SetCoord(std::min(mn.X(), p.X()), std::min(mn.Y(), p.Y()),
                        std::min(mn.Z(), p.Z()));
            mx.SetCoord(std::max(mx.X(), p.X()), std::max(mx.Y(), p.Y()),
                        std::max(mx.Z(), p.Z()));
        }
    }
    const gp_XYZ span(mx.X() - mn.X(), mx.Y() - mn.Y(), mx.Z() - mn.Z());
    auto isCap = [](double thin, double u, double v) {
        const double inMax = std::max(u, v);
        const double inMin = std::min(u, v);
        if (inMax < 1e-12) return false;
        // Real extrusion caps are thin along the normal and have two
        // comparable in-plane spans. Fillet slivers are long in one
        // direction and fail the in-plane aspect test.
        if (thin > inMax * 0.05 + 1e-9) return false;
        return inMin > 0.25 * inMax;
    };
    if (ax == am && isCap(span.X(), span.Y(), span.Z())) return true;
    if (ay == am && isCap(span.Y(), span.X(), span.Z())) return true;
    if (az == am && isCap(span.Z(), span.X(), span.Y())) return true;
    return false;
}

bool samePlaneLikeNbr(const Nbr& a, const Nbr& b, const SegmentWork& work,
                      const DerivedTols& tol) {
    if (!nbrIsPlaneLike(a, work) || !nbrIsPlaneLike(b, work)) return false;
    gp_Ax3 axA, axB;
    if (!planeOfNbr(a, work, axA) || !planeOfNbr(b, work, axB)) return false;
    const gp_XYZ nA = planeNormal(axA);
    const gp_XYZ nB = planeNormal(axB);
    return nA.Dot(nB) >= std::cos(tol.thetaPlane) - 1e-12;
}

void mergeCoplanarTallies(std::vector<NbrTally>& tallies, const SegmentWork& work,
                          const DerivedTols& tol) {
    std::vector<NbrTally> merged;
    merged.reserve(tallies.size());
    for (const NbrTally& t : tallies) {
        bool placed = false;
        for (NbrTally& m : merged) {
            if (!samePlaneLikeNbr(t.n, m.n, work, tol)) continue;
            m.count += t.count;
            if (t.n < m.n) m.n = t.n;   // deterministic representative (I5)
            placed = true;
            break;
        }
        if (!placed) merged.push_back(t);
    }
    tallies.swap(merged);
}

int nbrId(const Nbr& n) {
    return n.idx;
}

void canonicalizeDir(gp_XYZ& a) {
    if (a.Z() < -1e-15
        || (std::abs(a.Z()) <= 1e-15 && a.Y() < -1e-15)
        || (std::abs(a.Z()) <= 1e-15 && std::abs(a.Y()) <= 1e-15 && a.X() < 0.0)) {
        a.Reverse();
    }
}

double angleFrom(const gp_XYZ& n, const gp_XYZ& mL, const gp_XYZ& axis) {
    const gp_XYZ v = axis.Crossed(mL);   // toward mR in the {mL,mR} plane
    return std::atan2(n.Dot(v), n.Dot(mL));
}

void uniqueSortedVerts(const MeshView& mv, const std::vector<int>& g2l,
                       const std::vector<int>& tris, std::vector<int>& out) {
    out.clear();
    std::vector<int> ord = tris;
    std::sort(ord.begin(), ord.end());
    for (int t : ord) {
        for (int k = 0; k < 3; ++k) {
            const int lv = cornerLocal(g2l, mv, t, k);
            if (lv >= 0) out.push_back(lv);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

void sortAcceptedRejected(const MeshView& mv, SegmentWork& work) {
    auto keyLess = [&](const Region& a, const Region& b) {
        const double aa = regionArea(mv, a.tris);
        const double ba = regionArea(mv, b.tris);
        if (aa != ba) return aa > ba;           // -area ascending
        const int am = minTriOf(a.tris);
        const int bm = minTriOf(b.tris);
        if (am != bm) return am < bm;
        return a.origin < b.origin;             // trailing discriminator
    };
    std::sort(work.accepted.begin(), work.accepted.end(), keyLess);
    std::sort(work.rejected.begin(), work.rejected.end(), keyLess);
}

}  // namespace

bool claimFilletsC1(const MeshView& mv, const SegmentParams& p,
                    const DerivedTols& tol, SegmentWork& work) {
    try {
        if (!p.doFillets) return true;
        if (!mv.pts || !mv.tris || mv.nTri == 0) return true;
        if (work.provisionals.empty()) return true;

        const int nProv = (int)work.provisionals.size();
        const int nTri  = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const int nVtx  = (int)mv.nVtx;

        // --- local ids -------------------------------------------------------
        std::vector<int> g2l;
        for (int v = 0; v < nVtx; ++v) {
            const int g = mv.compVtx ? mv.compVtx[v] : v;
            if (g >= (int)g2l.size()) g2l.resize((std::size_t)g + 1, -1);
            g2l[g] = v;
        }

        std::vector<int> triProv((std::size_t)nTri, -1);
        for (int i = 0; i < nProv; ++i) {
            std::sort(work.provisionals[i].tris.begin(),
                      work.provisionals[i].tris.end());
            for (int t : work.provisionals[i].tris) {
                if (t >= 0 && t < nTri) triProv[t] = i;
            }
        }

        std::vector<int> triAcc((std::size_t)nTri, -1);
        for (int i = 0; i < (int)work.accepted.size(); ++i) {
            std::sort(work.accepted[i].tris.begin(), work.accepted[i].tris.end());
            for (int t : work.accepted[i].tris) {
                if (t >= 0 && t < nTri) triAcc[t] = i;
            }
        }

        // edge -> up to two incident triangles
        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{-1, -1});
        if (mv.triEdges) {
            for (int t = 0; t < nTri; ++t) {
                for (int s = 0; s < 3; ++s) {
                    const int e = mv.triEdges[t][s];
                    if (e < 0 || e >= nEdge) continue;
                    if (edgeTris[e][0] < 0) edgeTris[e][0] = t;
                    else if (edgeTris[e][1] < 0) edgeTris[e][1] = t;
                }
            }
        }

        auto nbrOfTri = [&](int nt) -> Nbr {
            if (nt < 0 || nt >= nTri) return {};
            if (triAcc[nt] >= 0) {
                const Origin o = work.accepted[triAcc[nt]].origin;
                // ChamferCone is not a FilletStrip neighbour (ineligible).
                // NgonWall falls through as Acc cylinder, same as CylGrow.
                if (o == Origin::ChamferCone) return {};
                Nbr n;
                n.kind = Nbr::Acc;
                n.idx = triAcc[nt];
                return n;
            }
            const int pi = triProv[nt];
            if (pi < 0) return {};
            const ProvClaim cl = work.provisionals[pi].claim;
            if (cl == ProvClaim::ConsumedCylinder || cl == ProvClaim::InCylinderClaim) {
                // B1 / NgonWall / ChamferCone consumed this band. Enough to
                // keep it out of Unclaimed members; Acc-miss still a cyl nbr.
                Nbr n;
                n.kind = Nbr::Prov;
                n.idx = pi;
                return n;
            }
            if (cl == ProvClaim::Unclaimed) {
                Nbr n;
                n.kind = Nbr::Prov;
                n.idx = pi;
                return n;
            }
            return {};
        };

        // unclaimed adjacency (shared mesh edge)
        std::vector<std::vector<int>> adj((std::size_t)nProv);
        if (mv.triEdges) {
            for (int e = 0; e < nEdge; ++e) {
                const int ta = edgeTris[e][0], tb = edgeTris[e][1];
                if (ta < 0 || tb < 0) continue;
                const int pa = triProv[ta], pb = triProv[tb];
                if (pa < 0 || pb < 0 || pa == pb) continue;
                if (work.provisionals[pa].claim != ProvClaim::Unclaimed) continue;
                if (work.provisionals[pb].claim != ProvClaim::Unclaimed) continue;
                adj[pa].push_back(pb);
                adj[pb].push_back(pa);
            }
        }
        for (int i = 0; i < nProv; ++i) {
            auto& v = adj[i];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }

        std::vector<int> unclaimed;
        unclaimed.reserve((std::size_t)nProv);
        for (int i = 0; i < nProv; ++i) {
            if (work.provisionals[i].claim == ProvClaim::Unclaimed
                && !work.provisionals[i].tris.empty()) {
                unclaimed.push_back(i);
            }
        }

        struct Cand {
            std::vector<int> mem;
            int n = 0;
            int minTri = std::numeric_limits<int>::max();
            int minProv = std::numeric_limits<int>::max();
        };

        std::vector<Cand> raw;
        auto finalizeCand = [&](Cand& cd) {
            std::sort(cd.mem.begin(), cd.mem.end());
            cd.mem.erase(std::unique(cd.mem.begin(), cd.mem.end()), cd.mem.end());
            while (!cd.mem.empty() && cd.mem.front() < 0) cd.mem.erase(cd.mem.begin());
            cd.n = (int)cd.mem.size();
            cd.minProv = cd.n ? cd.mem[0] : std::numeric_limits<int>::max();
            cd.minTri = std::numeric_limits<int>::max();
            for (int m : cd.mem)
                cd.minTri = std::min(cd.minTri, minTriOfProv(work.provisionals[m]));
        };
        auto pushCand = [&](int a, int b, int c) {
            Cand cd;
            if (a >= 0) cd.mem.push_back(a);
            if (b >= 0) cd.mem.push_back(b);
            if (c >= 0) cd.mem.push_back(c);
            finalizeCand(cd);
            if (cd.n > 0) raw.push_back(std::move(cd));
        };
        auto pushCandMem = [&](std::vector<int> mem) {
            Cand cd;
            cd.mem = std::move(mem);
            finalizeCand(cd);
            if (cd.n > 0) raw.push_back(std::move(cd));
        };

        auto isAdj = [&](int a, int b) -> bool {
            if (a < 0 || b < 0 || a >= nProv || b >= nProv) return false;
            const auto& v = adj[a];
            return std::binary_search(v.begin(), v.end(), b);
        };

        // Emit each |S|<=3 subset once from its min-index vertex. Replaces the
        // O(|raw|^2) already() scan: size-1/2 are unique by i / (i<j); size-3
        // is unique because the third vertex is taken from adj[i] (k>j) or
        // from adj[j] excluding adj[i] (so triangles are not double-counted).
        std::vector<char> seenK((std::size_t)nProv, 0);
        for (int i : unclaimed) {
            if (provIsExtrusionCap(mv, g2l, work.provisionals[i])) continue;
            pushCand(i, -1, -1);
            for (int j : adj[i]) {
                if (j <= i) continue;
                if (provIsExtrusionCap(mv, g2l, work.provisionals[j])) continue;
                pushCand(i, j, -1);
                if (nProv > 0) std::fill(seenK.begin(), seenK.end(), 0);
                for (int k : adj[i]) {
                    if (k <= j || k == i) continue;
                    if (provIsExtrusionCap(mv, g2l, work.provisionals[k])) continue;
                    if (seenK[(std::size_t)k]) continue;
                    seenK[(std::size_t)k] = 1;
                    pushCand(i, j, k);
                }
                for (int k : adj[j]) {
                    if (k == i || k <= i) continue;
                    if (isAdj(i, k)) continue;
                    if (provIsExtrusionCap(mv, g2l, work.provisionals[k])) continue;
                    if (seenK[(std::size_t)k]) continue;
                    seenK[(std::size_t)k] = 1;
                    pushCand(i, j, k);
                }
            }
        }

        // Tessellated edge rounds (S02 nSides=20) are longer than |S|<=3.
        // Flood small unclaimed bands from each large-plane seed; stop at a
        // third large plane (vertex blend / cube corner).
        {
            double maxA = 0.0;
            for (int i : unclaimed)
                maxA = std::max(maxA, work.provisionals[i].area);
            const double stripCut = 0.25 * maxA;
            std::vector<char> isLarge((std::size_t)nProv, 0);
            for (int i : unclaimed) {
                if (work.provisionals[i].area > stripCut + 1e-15)
                    isLarge[i] = 1;
            }
            auto largeNbrsOf = [&](int pi, std::vector<int>& out) {
                out.clear();
                for (int j : adj[pi]) if (isLarge[j]) out.push_back(j);
                std::sort(out.begin(), out.end());
                out.erase(std::unique(out.begin(), out.end()), out.end());
            };
            auto axisId = [&](int pi) {
                gp_XYZ mn(1e300, 1e300, 1e300), mx(-1e300, -1e300, -1e300);
                for (int t : work.provisionals[pi].tris) {
                    for (int k = 0; k < 3; ++k) {
                        const int lv = cornerLocal(g2l, mv, t, k);
                        if (lv < 0) continue;
                        const gp_XYZ p = localPt(mv, lv);
                        mn.SetCoord(std::min(mn.X(), p.X()), std::min(mn.Y(), p.Y()),
                                    std::min(mn.Z(), p.Z()));
                        mx.SetCoord(std::max(mx.X(), p.X()), std::max(mx.Y(), p.Y()),
                                    std::max(mx.Z(), p.Z()));
                    }
                }
                const double ax = std::abs(mx.X() - mn.X());
                const double ay = std::abs(mx.Y() - mn.Y());
                const double az = std::abs(mx.Z() - mn.Z());
                if (ax >= ay && ax >= az) return 0;
                if (ay >= az) return 1;
                return 2;
            };
            std::vector<int> largeSeeds;
            for (int i : unclaimed) if (isLarge[i]) largeSeeds.push_back(i);
            std::sort(largeSeeds.begin(), largeSeeds.end(), [&](int a, int b) {
                const int ma = minTriOfProv(work.provisionals[a]);
                const int mb = minTriOfProv(work.provisionals[b]);
                if (ma != mb) return ma < mb;
                return a < b;
            });
            std::vector<int> ln, lnV;
            for (int P : largeSeeds) {
                std::vector<int> seeds = adj[P];
                std::sort(seeds.begin(), seeds.end(), [&](int a, int b) {
                    const int ma = minTriOfProv(work.provisionals[a]);
                    const int mb = minTriOfProv(work.provisionals[b]);
                    if (ma != mb) return ma < mb;
                    return a < b;
                });
                for (int seed : seeds) {
                    if (isLarge[seed]) continue;
                    if (provIsExtrusionCap(mv, g2l, work.provisionals[seed])) continue;
                    std::vector<char> in((std::size_t)nProv, 0);
                    std::vector<int> stack = {seed};
                    std::vector<int> comp;
                    std::vector<int> flanks = {P};
                    const int ax0 = axisId(seed);
                    in[seed] = 1;
                    while (!stack.empty() && (int)comp.size() < 32) {
                        const int u = stack.back();
                        stack.pop_back();
                        if (axisId(u) != ax0) continue;
                        comp.push_back(u);
                        largeNbrsOf(u, ln);
                        for (int L : ln) {
                            if (std::find(flanks.begin(), flanks.end(), L) == flanks.end())
                                flanks.push_back(L);
                        }
                        for (int v : adj[u]) {
                            if (in[v] || isLarge[v]) continue;
                            if (axisId(v) != ax0) continue;
                            if (provIsExtrusionCap(mv, g2l, work.provisionals[v])) continue;
                            largeNbrsOf(v, lnV);
                            std::vector<int> f2 = flanks;
                            for (int L : lnV) {
                                if (std::find(f2.begin(), f2.end(), L) == f2.end())
                                    f2.push_back(L);
                            }
                            if ((int)f2.size() > 2) continue;
                            if ((int)f2.size() == 2) {
                                gp_XYZ nA = planeNormal(work.provisionals[f2[0]].plane);
                                gp_XYZ nB = planeNormal(work.provisionals[f2[1]].plane);
                                const double mA = nA.Modulus(), mB = nB.Modulus();
                                if (mA > 1e-15 && mB > 1e-15
                                    && nA.Dot(nB) < -0.95 * mA * mB)
                                    continue;  // through-hole top/bottom, not a fillet
                                gp_XYZ a = nA.Crossed(nB);
                                const double am = a.Modulus();
                                if (am > 1e-12) {
                                    a /= am;
                                    gp_XYZ nn = planeNormal(work.provisionals[v].plane);
                                    const double nm = nn.Modulus();
                                    if (nm > 1e-15) {
                                        nn /= nm;
                                        if (std::abs(nn.Dot(a)) > tol.gaussAxisTiltSin() + 1e-15)
                                            continue;
                                    }
                                }
                            }
                            in[v] = 1;
                            flanks.swap(f2);
                            stack.push_back(v);
                        }
                    }
                    if ((int)comp.size() >= 2 && (int)flanks.size() >= 1)
                        pushCandMem(std::move(comp));
                }
            }
        }

        auto inS = [](const Cand& c, int pidx) {
            for (int k = 0; k < c.n; ++k)
                if (c.mem[k] == pidx) return true;
            return false;
        };

        auto memberSet = [](const Cand& c, std::vector<char>& mark, int nProv_) {
            mark.assign((std::size_t)nProv_, 0);
            for (int k = 0; k < c.n; ++k)
                if (c.mem[k] >= 0) mark[c.mem[k]] = 1;
        };

        struct Ready {
            Cand cand;
            Nbr L, R;
            gp_Ax3 axL, axR;
            gp_XYZ mL, mR, p0L, p0R;
            double g = 0;
            double wedge = 0;          // radians, angle(mL, mR)
            bool cylNbr = false;
        };

        std::vector<Ready> valid;
        int skTally = 0, skCover = 0, skArea = 0, skNa = 0, skPar = 0, skMono = 0;
        int skF1 = 0, skF2 = 0, nF1pair = 0, nF2pair = 0;

        auto diagOn = []() {
            const char* e = std::getenv("STL2STEP_P1_DIAG");
            return e && e[0] && e[0] != '0';
        };

        for (Cand c : raw) {
            std::vector<char> mark;
            memberSet(c, mark, nProv);

            std::vector<int> sTrisCand;
            for (int k = 0; k < c.n; ++k) {
                const auto& ts = work.provisionals[c.mem[k]].tris;
                sTrisCand.insert(sTrisCand.end(), ts.begin(), ts.end());
            }
            std::sort(sTrisCand.begin(), sTrisCand.end());
            const gp_XYZ stripRun = stripExtrusionAxis(mv, g2l, sTrisCand);

            std::vector<NbrTally> tallies;
            int nOut = 0;
            for (int k = 0; k < c.n; ++k) {
                const int pi = c.mem[k];
                for (int t : work.provisionals[pi].tris) {
                    if (!mv.triEdges) continue;
                    for (int s = 0; s < 3; ++s) {
                        const int e = mv.triEdges[t][s];
                        const int ot = edgeOtherTri(edgeTris, e, t);
                        if (ot < 0) continue;                       // open: skip
                        const int op = (ot < nTri) ? triProv[ot] : -1;
                        if (op >= 0 && op < nProv && mark[op]) continue;  // interior
                        Nbr nb = nbrOfTri(ot);
                        if (!nb.valid()) continue;
                        if (nb.kind == Nbr::Prov && nb.idx >= 0 && nb.idx < nProv
                            && mark[nb.idx])
                            continue;
                        if (nbrIsPlaneLike(nb, work)) {
                            gp_Ax3 axN;
                            if (planeOfNbr(nb, work, axN)
                                && std::abs(planeNormal(axN).Dot(stripRun)) > 0.95) {
                                continue;
                            }
                        }
                        nOut++;
                        tallyAdd(tallies, nb);
                    }
                }
            }
            if (nOut <= 0 || tallies.size() < 2) {
                ++skTally;
                if (diagOn() && c.n >= 4)
                    std::fprintf(stderr, "  skip n=%d nOut=%d nTally=%zu\n",
                                 c.n, nOut, tallies.size());
                continue;
            }

            mergeCoplanarTallies(tallies, work, tol);
            if (tallies.size() < 2) { ++skTally; continue; }

            auto nbrArea = [&](const Nbr& n) -> double {
                if (n.kind == Nbr::Prov && n.idx >= 0 && n.idx < nProv)
                    return work.provisionals[n.idx].area;
                if (n.kind == Nbr::Acc && n.idx >= 0
                    && n.idx < (int)work.accepted.size())
                    return regionArea(mv, work.accepted[n.idx].tris);
                return 0.0;
            };
            // D7.5: rank by (count desc, nbrArea desc, Nbr id asc); cap at 8.
            std::sort(tallies.begin(), tallies.end(),
                      [&](const NbrTally& a, const NbrTally& b) {
                          if (a.count != b.count) return a.count > b.count;
                          const double aa = nbrArea(a.n), ba = nbrArea(b.n);
                          if (aa != ba) return aa > ba;
                          return a.n < b.n;
                      });
            if (tallies.size() > 8) {
                std::fprintf(stderr,
                             "C1 D7.5: truncating neighbour tallies %zu -> 8 "
                             "(n=%d minTri=%d)\n",
                             tallies.size(), c.n, c.minTri);
                tallies.resize(8);
            }

            double aS = 0.0;
            for (int k = 0; k < c.n; ++k)
                aS += work.provisionals[c.mem[k]].area;

            struct PairPick {
                bool ok = false;
                Ready rd;
                double spread = 1e300;
                int countSum = 0;
                double areaSum = 0.0;
                Nbr nLo, nHi;
                double Rref = 0.0;
            };
            PairPick best;
            bool anyF1 = false, anyF2 = false, anyPairOk = false;

            const int nT = (int)tallies.size();
            for (int ia = 0; ia < nT; ++ia) {
                for (int ib = ia + 1; ib < nT; ++ib) {
                    Nbr n0 = tallies[ia].n;
                    Nbr n1 = tallies[ib].n;
                    if (n0 == n1) continue;
                    Cand workC = c;

                    const double a0 = nbrArea(n0), a1 = nbrArea(n1);
                    if (a0 < 2.0 * aS && a1 < 2.0 * aS) continue;

                    const bool n0ok = nbrIsPlaneLike(n0, work) || nbrIsCylinder(n0, work);
                    const bool n1ok = nbrIsPlaneLike(n1, work) || nbrIsCylinder(n1, work);
                    if (!n0ok || !n1ok) continue;
                    if (n0.kind == Nbr::Prov && inS(workC, n0.idx)) continue;
                    if (n1.kind == Nbr::Prov && inS(workC, n1.idx)) continue;

                    gp_Ax3 ax0, ax1;
                    if (!planeOfNbr(n0, work, ax0) || !planeOfNbr(n1, work, ax1))
                        continue;

                    gp_XYZ flank0 = planeNormal(ax0);
                    gp_XYZ flank1 = planeNormal(ax1);
                    const double f0m = flank0.Modulus();
                    const double f1m = flank1.Modulus();
                    if (f0m < 1e-15 || f1m < 1e-15) continue;
                    flank0 /= f0m;
                    flank1 /= f1m;
                    if (workC.n > 1 && workC.n < 4) {
                        const double cosPlane = std::cos(tol.thetaPlane);
                        bool hasFlankPlaneMember = false;
                        for (int k = 0; k < workC.n; ++k) {
                            gp_XYZ pn = planeNormal(work.provisionals[workC.mem[k]].plane);
                            const double pnm = pn.Modulus();
                            if (pnm < 1e-15) {
                                hasFlankPlaneMember = true;
                                break;
                            }
                            pn /= pnm;
                            if (std::abs(pn.Dot(flank0)) >= cosPlane - 1e-12
                                || std::abs(pn.Dot(flank1)) >= cosPlane - 1e-12) {
                                hasFlankPlaneMember = true;
                                break;
                            }
                        }
                        if (hasFlankPlaneMember) continue;
                    }

                    const bool cyl = nbrIsCylinder(n0, work) || nbrIsCylinder(n1, work);

                    std::vector<int> sTris;
                    for (int k = 0; k < workC.n; ++k)
                        sTris.insert(sTris.end(),
                                     work.provisionals[workC.mem[k]].tris.begin(),
                                     work.provisionals[workC.mem[k]].tris.end());
                    std::sort(sTris.begin(), sTris.end());
                    std::vector<int> sVerts;
                    uniqueSortedVerts(mv, g2l, sTris, sVerts);
                    if (sVerts.empty()) continue;

                    auto signedNormal = [&](const gp_Ax3& ax) {
                        const gp_XYZ n = planeNormal(ax);
                        const gp_XYZ p0 = planePoint(ax);
                        double acc = 0.0;
                        for (int lv : sVerts) acc += n.Dot(localPt(mv, lv) - p0);
                        const double sgn = (acc >= 0.0) ? 1.0 : -1.0;
                        return n * sgn;
                    };

                    gp_XYZ m0 = signedNormal(ax0);
                    gp_XYZ m1 = signedNormal(ax1);
                    const double mag0 = m0.Modulus();
                    const double mag1 = m1.Modulus();
                    if (mag0 < 1e-15 || mag1 < 1e-15) continue;
                    m0 /= mag0;
                    m1 /= mag1;

                    gp_XYZ axis = m0.Crossed(m1);
                    double axisMag = axis.Modulus();
                    if (axisMag >= 1e-12 && workC.n >= 4) {
                        const gp_XYZ aHat = axis / axisMag;
                        const double sin3 = tol.gaussAxisTiltSin();
                        std::vector<int> kept;
                        kept.reserve(workC.mem.size());
                        for (int pi : workC.mem) {
                            gp_XYZ nn = planeNormal(work.provisionals[pi].plane);
                            const double nm = nn.Modulus();
                            if (nm < 1e-15) continue;
                            nn /= nm;
                            if (std::abs(nn.Dot(aHat)) <= sin3 + 1e-15)
                                kept.push_back(pi);
                        }
                        if (kept.empty()) continue;
                        if (kept.size() != workC.mem.size()) {
                            workC.mem.swap(kept);
                            workC.n = (int)workC.mem.size();
                            sTris.clear();
                            for (int pi : workC.mem) {
                                sTris.insert(sTris.end(),
                                             work.provisionals[pi].tris.begin(),
                                             work.provisionals[pi].tris.end());
                            }
                            std::sort(sTris.begin(), sTris.end());
                            sTris.erase(std::unique(sTris.begin(), sTris.end()),
                                        sTris.end());
                            uniqueSortedVerts(mv, g2l, sTris, sVerts);
                            if (sVerts.empty()) continue;
                        }
                    }
                    bool haveAxis = axisMag >= 1e-12;
                    if (haveAxis) axis /= axisMag;

                    const double g = std::max(-1.0, std::min(1.0, m0.Dot(m1)));
                    const double wedge = std::acos(g);
                    if (std::abs(g) > 0.90) continue;

                    auto monotoneFrom = [&](const gp_XYZ& mL, const gp_XYZ& ax,
                                            const gp_XYZ& mR, double wedgeRad) -> bool {
                        if (!haveAxis) return wedgeRad < 1e-6;
                        std::vector<std::pair<double, int>> ang;
                        ang.reserve((std::size_t)workC.n);
                        for (int k = 0; k < workC.n; ++k) {
                            gp_XYZ nn = planeNormal(
                                work.provisionals[workC.mem[k]].plane);
                            const double nm = nn.Modulus();
                            if (nm < 1e-15) return false;
                            nn /= nm;
                            signTowardBisector(nn, mL, mR);
                            const double th = angleFrom(nn, mL, ax);
                            ang.push_back({th, workC.mem[k]});
                        }
                        std::sort(ang.begin(), ang.end(),
                                  [](const std::pair<double, int>& a,
                                     const std::pair<double, int>& b) {
                                      if (a.first != b.first) return a.first < b.first;
                                      return a.second < b.second;
                                  });
                        const double hi = wedgeRad + 1e-6;
                        for (std::size_t i = 0; i < ang.size(); ++i) {
                            if (!(ang[i].first > -0.25 && ang[i].first < hi + 0.25))
                                return false;
                        }
                        if (workC.n < 4) {
                            for (std::size_t i = 1; i < ang.size(); ++i) {
                                const int a = ang[i - 1].second, b = ang[i].second;
                                bool hit = false;
                                for (int nb : adj[a])
                                    if (nb == b) {
                                        hit = true;
                                        break;
                                    }
                                if (!hit && ang.size() > 1) return false;
                            }
                        }
                        return true;
                    };

                    Nbr L = n0, R = n1;
                    gp_XYZ mL = m0, mR = m1;
                    gp_Ax3 axL = ax0, axR = ax1;
                    gp_XYZ axisUse = axis;

                    bool mono = haveAxis && monotoneFrom(mL, axisUse, mR, wedge);
                    if (!mono && haveAxis) {
                        L = n1;
                        R = n0;
                        mL = m1;
                        mR = m0;
                        axL = ax1;
                        axR = ax0;
                        axisUse = mL.Crossed(mR);
                        const double am = axisUse.Modulus();
                        if (am >= 1e-12) {
                            axisUse /= am;
                            mono = monotoneFrom(mL, axisUse, mR, wedge);
                        }
                    }
                    if (!haveAxis) {
                        L = n0;
                        R = n1;
                        mL = m0;
                        mR = m1;
                        axL = ax0;
                        axR = ax1;
                        mono = true;
                    }
                    if (!mono) continue;

                    double minTh = 1e300, maxTh = -1e300;
                    if (haveAxis) {
                        for (int k = 0; k < workC.n; ++k) {
                            gp_XYZ nn = planeNormal(
                                work.provisionals[workC.mem[k]].plane);
                            const double nm = nn.Modulus();
                            if (nm < 1e-15) continue;
                            nn /= nm;
                            signTowardBisector(nn, mL, mR);
                            const double th = angleFrom(nn, mL, axisUse);
                            minTh = std::min(minTh, th);
                            maxTh = std::max(maxTh, th);
                        }
                    }
                    const double stripSpan = (maxTh >= minTh) ? (maxTh - minTh) : 0.0;
                    const double oneBand = (workC.n >= 1) ? wedge / (double)workC.n : wedge;
                    const double spanPlus = haveAxis ? (stripSpan + oneBand) : wedge;
                    const double spanDeg = spanPlus * 180.0 / kPi;
                    const bool spanOk =
                        (spanDeg + 1e-9 >= kDeg30 && spanDeg - 1e-9 <= kDeg180)
                        || (!haveAxis);
                    if (!spanOk) continue;

                    const gp_XYZ p0L = planePoint(axL);
                    const gp_XYZ p0R = planePoint(axR);
                    const double planeLim = std::max(tol.epsMesh, 1e-12 * mv.diag);
                    std::vector<int> chainL, chainR;
                    for (int lv : sVerts) {
                        const gp_XYZ v = localPt(mv, lv);
                        if (std::abs(mL.Dot(v - p0L)) <= planeLim) chainL.push_back(lv);
                        if (std::abs(mR.Dot(v - p0R)) <= planeLim) chainR.push_back(lv);
                    }
                    // D7.6 (f1): two vertex-disjoint flank chains of size >= 2.
                    bool overlap = false;
                    for (int a : chainL) {
                        for (int b : chainR) {
                            if (a == b) {
                                overlap = true;
                                break;
                            }
                        }
                        if (overlap) break;
                    }
                    if (overlap || (int)chainL.size() < 2 || (int)chainR.size() < 2) {
                        ++nF1pair;
                        anyF1 = true;
                        continue;
                    }

                    auto riOf = [&](int lv, bool& ok) {
                        const gp_XYZ v = localPt(mv, lv);
                        double dL = mL.Dot(v - p0L);
                        double dR = mR.Dot(v - p0R);
                        return d7Radius(dL, dR, g, ok);
                    };
                    std::vector<double> chainRi;
                    bool rOk = true;
                    for (int lv : chainL) {
                        bool ok = false;
                        const double Ri = riOf(lv, ok);
                        if (!ok) { rOk = false; break; }
                        chainRi.push_back(Ri);
                    }
                    for (int lv : chainR) {
                        bool ok = false;
                        const double Ri = riOf(lv, ok);
                        if (!ok) { rOk = false; break; }
                        chainRi.push_back(Ri);
                    }
                    if (!rOk || chainRi.empty()) {
                        ++nF2pair;
                        anyF2 = true;
                        continue;
                    }
                    const double Rref = medianInPlace(chainRi);
                    const double runBound = std::max(tol.epsMesh, 0.25 * Rref);
                    bool f2ok = true;
                    std::vector<double> allRi;
                    allRi.reserve(sVerts.size());
                    for (int lv : sVerts) {
                        bool ok = false;
                        const double Ri = riOf(lv, ok);
                        if (!ok) { f2ok = false; break; }
                        allRi.push_back(Ri);
                        if (std::abs(Ri - Rref) > runBound) {
                            f2ok = false;
                            break;
                        }
                    }
                    if (!f2ok) {
                        ++nF2pair;
                        anyF2 = true;
                        continue;
                    }
                    const double Rmed = medianInPlace(allRi);
                    double maxDev = 0.0;
                    for (double Ri : allRi)
                        maxDev = std::max(maxDev, std::abs(Ri - Rmed));
                    const double denom = std::max(Rmed, tol.epsMesh);
                    const double spread = (denom > 0.0) ? (maxDev / denom) : maxDev;

                    const int countSum = tallies[ia].count + tallies[ib].count;
                    const double areaSum = a0 + a1;
                    Nbr nLo = n0, nHi = n1;
                    if (nHi < nLo) std::swap(nLo, nHi);

                    bool better = !best.ok;
                    if (!better) {
                        if (spread != best.spread) better = spread < best.spread;
                        else if (countSum != best.countSum)
                            better = countSum > best.countSum;
                        else if (areaSum != best.areaSum) better = areaSum > best.areaSum;
                        else if (nLo < best.nLo) better = true;
                        else if (!(best.nLo < nLo) && nHi < best.nHi) better = true;
                    }
                    if (!better) continue;

                    anyPairOk = true;
                    best.ok = true;
                    best.spread = spread;
                    best.countSum = countSum;
                    best.areaSum = areaSum;
                    best.nLo = nLo;
                    best.nHi = nHi;
                    best.Rref = Rref;
                    best.rd.cand = workC;
                    best.rd.L = L;
                    best.rd.R = R;
                    best.rd.axL = axL;
                    best.rd.axR = axR;
                    best.rd.mL = mL;
                    best.rd.mR = mR;
                    best.rd.p0L = p0L;
                    best.rd.p0R = p0R;
                    best.rd.g = g;
                    best.rd.wedge = wedge;
                    best.rd.cylNbr = cyl;
                    if (diagOn() && workC.n >= 1) {
                        std::fprintf(stderr,
                                     "  D7.5 pair n=%d minTri=%d spread=%.4g "
                                     "Rref=%.5g count=%d area=%.4g\n",
                                     workC.n, workC.minTri, spread, Rref,
                                     countSum, areaSum);
                    }
                }
            }

            if (!best.ok) {
                if (!anyPairOk && anyF1 && !anyF2) ++skF1;
                else if (!anyPairOk && anyF2) ++skF2;
                continue;
            }
            const int cover = best.countSum;
            if (c.n == 1 && cover * 2 < nOut) { ++skCover; continue; }
            if (c.n >= 2 && c.n < 4 && cover * 5 < nOut * 4) { ++skCover; continue; }
            valid.push_back(best.rd);
        }

        // Maximal: drop S if a valid S' properly contains it.
        std::vector<char> keep(valid.size(), 1);
        for (std::size_t i = 0; i < valid.size(); ++i) {
            for (std::size_t j = 0; j < valid.size(); ++j) {
                if (i == j) continue;
                if (valid[j].cand.n <= valid[i].cand.n) continue;
                bool subset = true;
                for (int k = 0; k < valid[i].cand.n; ++k) {
                    if (!inS(valid[j].cand, valid[i].cand.mem[k])) {
                        subset = false;
                        break;
                    }
                }
                if (subset) {
                    keep[i] = 0;
                    break;
                }
            }
        }
        std::vector<Ready> maximal;
        for (std::size_t i = 0; i < valid.size(); ++i)
            if (keep[i]) maximal.push_back(valid[i]);

        // D1 §1.6: visit in ascending minTriId(S); key ends in an id.
        std::sort(maximal.begin(), maximal.end(),
                  [](const Ready& a, const Ready& b) {
                      if (a.cand.minTri != b.cand.minTri)
                          return a.cand.minTri < b.cand.minTri;
                      if (a.cand.minProv != b.cand.minProv)
                          return a.cand.minProv < b.cand.minProv;
                      return a.cand.n < b.cand.n;
                  });

        std::vector<char> used((std::size_t)nProv, 0);

        auto emitRejected = [&](const Ready& rd, Reject code) {
            Region r;
            r.id = (int)work.rejected.size();
            r.type = SurfType::Cylinder;
            // Rejected strips are not committed; leave origin at the default so
            // FilletConsensus fixtures have Origin::FilletStrip nowhere in output.
            r.origin = Origin::PlaneGrow;
            r.reject = code;
            r.builtAs = BuiltAs::NotBuilt;
            // Provisional indices (not Region::id — those do not exist until
            // stage D). buildTopologyD remaps filletNbrA/B onto dense ids.
            r.filletNbrA = nbrId(rd.L);
            r.filletNbrB = nbrId(rd.R);
            r.closed360 = false;
            for (int k = 0; k < rd.cand.n; ++k) {
                const auto& ts = work.provisionals[rd.cand.mem[k]].tris;
                r.tris.insert(r.tris.end(), ts.begin(), ts.end());
            }
            std::sort(r.tris.begin(), r.tris.end());
            r.tris.erase(std::unique(r.tris.begin(), r.tris.end()), r.tris.end());
            work.rejected.push_back(std::move(r));
        };

        for (const Ready& rd : maximal) {
            bool overlap = false;
            for (int k = 0; k < rd.cand.n; ++k) {
                if (used[rd.cand.mem[k]]) {
                    overlap = true;
                    break;
                }
                if (work.provisionals[rd.cand.mem[k]].claim != ProvClaim::Unclaimed) {
                    overlap = true;
                    break;
                }
            }
            if (overlap) continue;

            if (rd.cylNbr) {
                emitRejected(rd, Reject::TorusNYI);
                continue;   // leave Unclaimed for A3
            }

            // Both plane: D7 radius.
            std::vector<int> sTris;
            for (int k = 0; k < rd.cand.n; ++k) {
                work.provisionals[rd.cand.mem[k]].claim = ProvClaim::InFilletClaim;
                const auto& ts = work.provisionals[rd.cand.mem[k]].tris;
                sTris.insert(sTris.end(), ts.begin(), ts.end());
            }
            std::sort(sTris.begin(), sTris.end());
            sTris.erase(std::unique(sTris.begin(), sTris.end()), sTris.end());

            std::vector<int> sVerts;
            uniqueSortedVerts(mv, g2l, sTris, sVerts);

            const gp_XYZ mL = rd.mL, mR = rd.mR;
            const gp_XYZ p0L = rd.p0L, p0R = rd.p0R;
            const double g = rd.g;

            auto RiOf = [&](int lv, bool& ok) {
                const gp_XYZ v = localPt(mv, lv);
                double dL = mL.Dot(v - p0L);
                double dR = mR.Dot(v - p0R);
                return d7Radius(dL, dR, g, ok);
            };

            bool guardHit = !(1.0 - g >= kCoplanarG);
            std::vector<double> allRi;
            allRi.reserve(sVerts.size());
            bool anyOk = false;
            if (!guardHit) {
                for (int lv : sVerts) {
                    bool ok = false;
                    const double Ri = RiOf(lv, ok);
                    if (!ok || !std::isfinite(Ri)) {
                        guardHit = true;
                        break;
                    }
                    allRi.push_back(Ri);
                    anyOk = true;
                }
            }

            auto rollback = [&]() {
                for (int k = 0; k < rd.cand.n; ++k)
                    work.provisionals[rd.cand.mem[k]].claim = ProvClaim::Unclaimed;
            };

            if (guardHit || !anyOk || allRi.empty()) {
                rollback();
                emitRejected(rd, Reject::FilletConsensus);
                continue;
            }

            const double Rmed = medianInPlace(allRi);
            // re-compute allRi (medianInPlace sorted a copy... it sorted allRi)
            // rebuild for maxdev
            std::vector<double> allRi2;
            allRi2.reserve(sVerts.size());
            double maxDev = 0.0;
            for (int lv : sVerts) {
                bool ok = false;
                const double Ri = RiOf(lv, ok);
                if (!ok) {
                    maxDev = 1e300;
                    break;
                }
                allRi2.push_back(Ri);
                maxDev = std::max(maxDev, std::abs(Ri - Rmed));
            }
            const double R = Rmed;
            const double epsF = tol.epsFilletTol(R);
            if (const char* e = std::getenv("STL2STEP_P1_DIAG"); e && e[0] && e[0] != '0') {
                std::fprintf(stderr,
                             "  C1 strip n=%d minTri=%d R=%.6g maxDev=%.6g "
                             "wedge=%.3f L=%d/%d R=%d/%d verts=%zu\n",
                             rd.cand.n, rd.cand.minTri, R, maxDev,
                             rd.wedge * 180.0 / kPi,
                             (int)rd.L.kind, rd.L.idx, (int)rd.R.kind, rd.R.idx,
                             sVerts.size());
                int shown = 0;
                for (int lv : sVerts) {
                    const gp_XYZ v = localPt(mv, lv);
                    const double dL = mL.Dot(v - p0L);
                    const double dR = mR.Dot(v - p0R);
                    bool ok = false;
                    const double Ri = d7Radius(dL, dR, g, ok);
                    std::fprintf(stderr,
                                 "    v=%d dL=%.5g dR=%.5g Ri=%.5g ok=%d\n",
                                 lv, dL, dR, Ri, (int)ok);
                    if (++shown >= 12) break;
                }
            }

            // Long chains: vertices of dual edges S — L and S — R.
            std::vector<char> onL((std::size_t)std::max(nVtx, 0), 0);
            std::vector<char> onR((std::size_t)std::max(nVtx, 0), 0);
            std::vector<char> markP((std::size_t)nProv, 0);
            for (int k = 0; k < rd.cand.n; ++k) markP[rd.cand.mem[k]] = 1;

            auto markEdgeVerts = [&](int t, int e, std::vector<char>& on) {
                if (!mv.compEdges || e < 0 || e >= nEdge) {
                    // fall back: all corners of t
                    for (int k = 0; k < 3; ++k) {
                        const int lv = cornerLocal(g2l, mv, t, k);
                        if (lv >= 0 && lv < nVtx) on[lv] = 1;
                    }
                    return;
                }
                const int lo = mv.compEdges[e].first;
                const int hi = mv.compEdges[e].second;
                if (lo >= 0 && lo < nVtx) on[lo] = 1;
                if (hi >= 0 && hi < nVtx) on[hi] = 1;
            };

            for (int t : sTris) {
                if (!mv.triEdges) continue;
                for (int s = 0; s < 3; ++s) {
                    const int e = mv.triEdges[t][s];
                    const int ot = edgeOtherTri(edgeTris, e, t);
                    if (ot < 0) continue;
                    const Nbr nb = nbrOfTri(ot);
                    if (nb == rd.L) markEdgeVerts(t, e, onL);
                    if (nb == rd.R) markEdgeVerts(t, e, onR);
                }
            }

            auto chainMedian = [&](const std::vector<char>& on, bool& okOut) {
                okOut = false;
                std::vector<double> rs;
                for (int lv = 0; lv < nVtx; ++lv) {
                    if (!on[lv]) continue;
                    bool ok = false;
                    const double Ri = RiOf(lv, ok);
                    if (!ok) return 0.0;
                    rs.push_back(Ri);
                }
                if (rs.empty()) return 0.0;
                okOut = true;
                return medianInPlace(rs);
            };

            bool cLok = false, cRok = false;
            const double medL = chainMedian(onL, cLok);
            const double medR = chainMedian(onR, cRok);

            const double spanDeg = rd.wedge * 180.0 / kPi;
            const bool passSpread = maxDev <= epsF + 1e-15;
            const bool passG4 = (R > 2.0 * tol.epsPlane) && (R < 2.0 * mv.diag);
            const bool passG5 = (spanDeg + 1e-9 >= kDeg30 && spanDeg - 1e-9 <= kDeg180);
            const bool passChains = cLok && cRok
                                    && std::abs(medL - R) <= epsF + 1e-15
                                    && std::abs(medR - R) <= epsF + 1e-15;
            if (!passG5) {
                rollback();   // D5.5: G5 decline is silent
                continue;
            }
            if (!(passSpread && passG4 && passChains)) {
                if (const char* e = std::getenv("STL2STEP_P1_DIAG"); e && e[0] && e[0] != '0') {
                    std::fprintf(stderr,
                                 "  C1 D7 fail n=%d R=%.5g maxDev=%.5g epsF=%.5g "
                                 "spread=%d G4=%d chains=%d medL=%.5g medR=%.5g wedge=%.3f\n",
                                 rd.cand.n, R, maxDev, epsF,
                                 (int)passSpread, (int)passG4, (int)passChains,
                                 medL, medR, spanDeg);
                }
                rollback();
                emitRejected(rd, Reject::FilletConsensus);
                continue;
            }
            // ---- commit Origin::FilletStrip --------------------------------
            gp_XYZ axis = mL.Crossed(mR);
            double am = axis.Modulus();
            if (am < 1e-15) {
                rollback();
                emitRejected(rd, Reject::FilletConsensus);
                continue;
            }
            axis /= am;
            canonicalizeDir(axis);

            // Intersection of offset planes p0 + R m.
            const double d1 = mL.Dot(p0L) + R;
            const double d2 = mR.Dot(p0R) + R;
            const double den = 1.0 - g * g;
            gp_XYZ loc0;
            if (std::abs(den) < 1e-12) {
                loc0 = (p0L + mL * R + p0R + mR * R) * 0.5;
            } else {
                const double alpha = (d1 - g * d2) / den;
                const double beta  = (d2 - g * d1) / den;
                loc0 = mL * alpha + mR * beta;
            }

            // Area-weighted centroid, ascending tri-id (I5).
            gp_XYZ cent(0, 0, 0);
            double areaSum = 0.0;
            for (int t : sTris) {
                const double at = triArea(mv, t);
                areaSum += at;
                cent += triCentroid(mv, t) * at;
            }
            if (areaSum > 0.0) cent /= areaSum;
            gp_XYZ loc = loc0 + axis * (cent - loc0).Dot(axis);

            // XDirection: lowest local vertex id of the strip.
            gp_XYZ xdir(0, 0, 0);
            bool gotX = false;
            for (int lv : sVerts) {
                gp_XYZ rad = localPt(mv, lv) - loc;
                rad -= axis * rad.Dot(axis);
                if (rad.Modulus() > 1e-12) {
                    xdir = rad / rad.Modulus();
                    gotX = true;
                    break;
                }
            }
            if (!gotX) {
                gp_XYZ w = (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0)
                                                      : gp_XYZ(0, 1, 0);
                xdir = axis.Crossed(w);
                const double xm = xdir.Modulus();
                if (xm < 1e-15) {
                    rollback();
                    emitRejected(rd, Reject::FilletConsensus);
                    continue;
                }
                xdir /= xm;
            }

            // nBands: D2.2 clustering of facet-normal azimuths about `axis`,
            // floored at |S| so a designed N-row fillet cannot under-split.
            int wk = 0;
            double wmin = std::abs(axis.X());
            if (std::abs(axis.Y()) < wmin) { wmin = std::abs(axis.Y()); wk = 1; }
            if (std::abs(axis.Z()) < wmin) { wk = 2; }
            const gp_XYZ ww(wk == 0 ? 1.0 : 0.0,
                            wk == 1 ? 1.0 : 0.0,
                            wk == 2 ? 1.0 : 0.0);
            gp_XYZ u0 = axis.Crossed(ww);
            const double u0m = u0.Modulus();
            if (u0m > 1e-15) u0 /= u0m;
            else u0 = gp_XYZ(1.0, 0.0, 0.0);
            const gp_XYZ v0 = axis.Crossed(u0);
            std::vector<double> psi;
            psi.reserve(sTris.size());
            for (int t : sTris) {
                gp_XYZ nt;
                triAreaNormal(mv, t, nt);
                psi.push_back(std::atan2(nt.Dot(v0), nt.Dot(u0)));
            }
            int nBands = nBandsFromNormals(psi, tol.thetaBin);
            if (nBands < 1) nBands = rd.cand.n;
            if (rd.cand.n > nBands) nBands = rd.cand.n;

            const double arcSpan = rd.wedge;    // n_L to n_R, radians
            if (arcSpan <= 0.0 || nBands < 1) {
                rollback();
                emitRejected(rd, Reject::FilletConsensus);
                continue;
            }
            int nSides = (int)std::llround(2.0 * kPi * (double)nBands / arcSpan);
            if (nSides < 1) nSides = 1;

            // Vertex azimuths for uMin/uMax / vMin/vMax. closed360 = false.
            const gp_XYZ vdir = axis.Crossed(xdir);
            std::vector<double> chi;
            chi.reserve(sVerts.size());
            double vMin = 1e300, vMax = -1e300;
            for (int lv : sVerts) {
                const gp_XYZ d = localPt(mv, lv) - loc;
                vMin = std::min(vMin, d.Dot(axis));
                vMax = std::max(vMax, d.Dot(axis));
                const double ur = d.Dot(xdir);
                const double vr = d.Dot(vdir);
                if (ur * ur + vr * vr < 1e-24) continue;
                chi.push_back(std::atan2(vr, ur));
            }
            std::sort(chi.begin(), chi.end());
            double spanV = arcSpan;
            double uMin = 0.0;
            if (chi.size() >= 2) {
                double maxGap = 0.0;
                std::size_t gapAt = 0;
                for (std::size_t i = 0; i + 1 < chi.size(); ++i) {
                    const double gp = chi[i + 1] - chi[i];
                    if (gp > maxGap) {
                        maxGap = gp;
                        gapAt = i;
                    }
                }
                const double wrap = (chi.front() + 2.0 * kPi) - chi.back();
                if (wrap > maxGap) {
                    maxGap = wrap;
                    gapAt = chi.size() - 1;
                }
                spanV = 2.0 * kPi - maxGap;
                if (spanV < 0.0) spanV = 0.0;
                const std::size_t start = (gapAt + 1) % chi.size();
                uMin = chi[start];
                // wrap to (-pi, pi]
                if (uMin > kPi) uMin -= 2.0 * kPi;
                if (uMin <= -kPi) uMin += 2.0 * kPi;
            }
            const double uMax = uMin + spanV;

            // outwardNormal: material on the axis side (D4.4).
            double sigmaAcc = 0.0;
            for (int t : sTris) {
                gp_XYZ nt;
                const double at = triAreaNormal(mv, t, nt);
                gp_XYZ rho = triCentroid(mv, t) - loc;
                rho -= axis * rho.Dot(axis);
                const double rm = rho.Modulus();
                if (rm < 1e-15) continue;
                rho /= rm;
                sigmaAcc += at * nt.Dot(rho);
            }
            const bool outward = sigmaAcc > 0.0;

            // vertex residuals vs cylinder, unique verts ascending id.
            double maxDevV = 0.0, sumSq = 0.0;
            int nDev = 0;
            for (int lv : sVerts) {
                const gp_XYZ d = localPt(mv, lv) - loc;
                const double rad = axis.Crossed(d).Modulus();
                const double dev = std::abs(rad - R);
                maxDevV = std::max(maxDevV, dev);
                sumSq += dev * dev;
                nDev++;
            }
            const double rms = (nDev > 0) ? std::sqrt(sumSq / (double)nDev) : 0.0;

            Region r;
            r.id = 0;                                  // stage D
            r.type = SurfType::Cylinder;
            r.origin = Origin::FilletStrip;
            r.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(axis), gp_Dir(xdir));
            r.radius = R;
            r.uMin = uMin;
            r.uMax = uMax;
            r.vMin = (vMin > vMax) ? 0.0 : vMin;
            r.vMax = (vMin > vMax) ? 0.0 : vMax;
            r.closed360 = false;
            r.outwardNormal = outward;
            r.tris = sTris;
            r.maxVertexDev = maxDevV;
            r.rmsVertexDev = rms;
            r.nSides = nSides;
            r.chordSagitta = chordSagitta(R, nSides);
            r.dVolPredicted = dVolCylinderSector(areaSum, R, nSides, outward);
            r.maxVertexSnap = 0.0;
            r.reject = Reject::None;
            r.builtAs = BuiltAs::NotBuilt;
            // Provisional indices. Region ids are assigned in stage D, which
            // remaps filletNbrA/B onto the committed neighbour Region::id.
            r.filletNbrA = nbrId(rd.L);
            r.filletNbrB = nbrId(rd.R);

            for (int k = 0; k < rd.cand.n; ++k)
                work.provisionals[rd.cand.mem[k]].claim = ProvClaim::ConsumedFillet;

            for (int k = 0; k < rd.cand.n; ++k)
                used[rd.cand.mem[k]] = 1;

            work.accepted.push_back(std::move(r));
        }

        if (const char* e = std::getenv("STL2STEP_P1_DIAG"); e && e[0] && e[0] != '0') {
            int nLong = 0;
            for (const Cand& c : raw)
                if (c.n >= 4) ++nLong;
            int nFil = 0;
            for (const Region& r : work.accepted)
                if (r.origin == Origin::FilletStrip) ++nFil;
            int nFc = 0;
            for (const Region& r : work.rejected)
                if (r.reject == Reject::FilletConsensus) ++nFc;
            std::fprintf(stderr,
                         "C1: raw=%zu long>=4=%d valid=%zu maximal=%zu "
                         "filletStrip=%d filletConsensus=%d unclaimed=%d "
                         "skip tally=%d cover=%d area=%d na=%d par=%d "
                         "f1=%d f2=%d f1pair=%d f2pair=%d\n",
                         raw.size(), nLong, valid.size(), maximal.size(),
                         nFil, nFc, (int)unclaimed.size(),
                         skTally, skCover, skArea, skNa, skPar,
                         skF1, skF2, nF1pair, nF2pair);
        }

        sortAcceptedRejected(mv, work);
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
