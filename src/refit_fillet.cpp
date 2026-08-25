// stl2step P1 stage C1 — fillet strip claim (D1 §1.5, D7 radius).
// gp_ / stdlib only. Never throws. Silent on every reject / NYI (D5.5).
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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
        return work.accepted[n.idx].type == SurfType::Cylinder;
    }
    if (n.kind == Nbr::Prov) {
        if (n.idx < 0 || n.idx >= (int)work.provisionals.size()) return false;
        const ProvClaim c = work.provisionals[n.idx].claim;
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
        return work.accepted[n.idx].type == SurfType::Plane;
    }
    return false;
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
                Nbr n;
                n.kind = Nbr::Acc;
                n.idx = triAcc[nt];
                return n;
            }
            const int pi = triProv[nt];
            if (pi < 0) return {};
            const ProvClaim cl = work.provisionals[pi].claim;
            if (cl == ProvClaim::ConsumedCylinder || cl == ProvClaim::InCylinderClaim) {
                // B1 consumed this band but the accepted-region map missed it.
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
            std::array<int, 3> mem{{-1, -1, -1}};
            int n = 0;
            int minTri = std::numeric_limits<int>::max();
            int minProv = std::numeric_limits<int>::max();
        };

        auto pack = [](int a, int b, int c) {
            std::array<int, 3> m = {a, b, c};
            std::sort(m.begin(), m.end(), [](int x, int y) {
                // push -1 to the back so size-1/2 packs compare cleanly
                if (x < 0) return false;
                if (y < 0) return true;
                return x < y;
            });
            return m;
        };

        std::vector<Cand> raw;
        auto pushCand = [&](int a, int b, int c) {
            const std::array<int, 3> m = pack(a, b, c);
            Cand cd;
            cd.mem = m;
            cd.n = (m[0] >= 0) + (m[1] >= 0) + (m[2] >= 0);
            cd.minProv = m[0];
            cd.minTri = std::numeric_limits<int>::max();
            for (int k = 0; k < cd.n; ++k)
                cd.minTri = std::min(cd.minTri, minTriOfProv(work.provisionals[m[k]]));
            raw.push_back(cd);
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
            pushCand(i, -1, -1);
            for (int j : adj[i]) {
                if (j <= i) continue;
                pushCand(i, j, -1);
                if (nProv > 0) std::fill(seenK.begin(), seenK.end(), 0);
                for (int k : adj[i]) {
                    if (k <= j || k == i) continue;
                    if (seenK[(std::size_t)k]) continue;
                    seenK[(std::size_t)k] = 1;
                    pushCand(i, j, k);
                }
                for (int k : adj[j]) {
                    if (k == i || k <= i) continue;
                    if (isAdj(i, k)) continue;
                    if (seenK[(std::size_t)k]) continue;
                    seenK[(std::size_t)k] = 1;
                    pushCand(i, j, k);
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

        for (const Cand& c : raw) {
            std::vector<char> mark;
            memberSet(c, mark, nProv);

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
                        // a neighbour that is itself a member (via Acc vs Prov) skip
                        if (nb.kind == Nbr::Prov && nb.idx >= 0 && nb.idx < nProv
                            && mark[nb.idx])
                            continue;
                        nOut++;
                        tallyAdd(tallies, nb);
                    }
                }
            }
            if (nOut <= 0 || tallies.size() < 2) continue;

            std::sort(tallies.begin(), tallies.end(),
                      [](const NbrTally& a, const NbrTally& b) {
                          if (a.count != b.count) return a.count > b.count;
                          return a.n < b.n;
                      });
            const int cover = tallies[0].count + tallies[1].count;
            if (cover * 5 < nOut * 4) continue;                     // < 80 %
            if (tallies[0].n == tallies[1].n) continue;

            Nbr n0 = tallies[0].n;
            Nbr n1 = tallies[1].n;

            const bool n0ok = nbrIsPlaneLike(n0, work) || nbrIsCylinder(n0, work);
            const bool n1ok = nbrIsPlaneLike(n1, work) || nbrIsCylinder(n1, work);
            if (!n0ok || !n1ok) continue;
            // A neighbour in S is not a "neighbour region".
            if (n0.kind == Nbr::Prov && inS(c, n0.idx)) continue;
            if (n1.kind == Nbr::Prov && inS(c, n1.idx)) continue;

            gp_Ax3 ax0, ax1;
            if (!planeOfNbr(n0, work, ax0) || !planeOfNbr(n1, work, ax1)) continue;

            const bool cyl = nbrIsCylinder(n0, work) || nbrIsCylinder(n1, work);

            // Sign both plane normals toward the strip (D7: toward the fillet axis).
            std::vector<int> sTris;
            for (int k = 0; k < c.n; ++k)
                sTris.insert(sTris.end(),
                             work.provisionals[c.mem[k]].tris.begin(),
                             work.provisionals[c.mem[k]].tris.end());
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
            // Near-parallel: still a candidate for the coplanar guard; skip
            // monotonicity that needs a well-defined axis.
            bool haveAxis = axisMag >= 1e-12;
            if (haveAxis) axis /= axisMag;

            const double g = std::max(-1.0, std::min(1.0, m0.Dot(m1)));
            const double wedge = std::acos(g);

            auto monotoneFrom = [&](const gp_XYZ& mL, const gp_XYZ& ax,
                                    const gp_XYZ& mR, double wedgeRad,
                                    Nbr /*L*/, Nbr /*R*/) -> bool {
                if (!haveAxis) return wedgeRad < 1e-6;     // degenerate
                std::vector<std::pair<double, int>> ang;
                ang.reserve((std::size_t)c.n);
                for (int k = 0; k < c.n; ++k) {
                    gp_XYZ nn = planeNormal(work.provisionals[c.mem[k]].plane);
                    const double nm = nn.Modulus();
                    if (nm < 1e-15) return false;
                    nn /= nm;
                    signTowardBisector(nn, mL, mR);
                    const double th = angleFrom(nn, mL, ax);
                    ang.push_back({th, c.mem[k]});
                }
                std::sort(ang.begin(), ang.end(),
                          [](const std::pair<double, int>& a,
                             const std::pair<double, int>& b) {
                              if (a.first != b.first) return a.first < b.first;
                              return a.second < b.second;
                          });
                // Monotonicity is the angular window plus the adjacency walk
                // in angle order. Re-checking that the just-sorted sequence is
                // sorted is tautological and is not the rule.
                const double hi = wedgeRad + 1e-6;
                for (std::size_t i = 0; i < ang.size(); ++i) {
                    if (!(ang[i].first > -0.25 && ang[i].first < hi + 0.25))
                        return false;
                }
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
                return true;
            };

            Nbr L = n0, R = n1;
            gp_XYZ mL = m0, mR = m1;
            gp_Ax3 axL = ax0, axR = ax1;
            gp_XYZ axisUse = axis;

            bool mono = haveAxis && monotoneFrom(mL, axisUse, mR, wedge, L, R);
            if (!mono && haveAxis) {
                // swap L/R
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
                    mono = monotoneFrom(mL, axisUse, mR, wedge, L, R);
                }
            }
            if (!haveAxis) {
                // near-coplanar: treat as a candidate so the D7 guard fires
                L = n0;
                R = n1;
                mL = m0;
                mR = m1;
                axL = ax0;
                axR = ax1;
                mono = true;
            }
            if (!mono) continue;

            // G5-fillet / candidate rule 4: strip-normal span + one band.
            double minTh = 1e300, maxTh = -1e300;
            if (haveAxis) {
                for (int k = 0; k < c.n; ++k) {
                    gp_XYZ nn = planeNormal(work.provisionals[c.mem[k]].plane);
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
            const double oneBand = (c.n >= 1) ? wedge / (double)c.n : wedge;
            const double spanPlus = haveAxis ? (stripSpan + oneBand) : wedge;
            const double spanDeg = spanPlus * 180.0 / kPi;
            // Near-coplanar (wedge ~ 0) is admitted so the 1-g guard can fire.
            const bool spanOk = (spanDeg + 1e-9 >= kDeg30 && spanDeg - 1e-9 <= kDeg180)
                                || (!haveAxis);
            if (!spanOk) continue;

            Ready rd;
            rd.cand = c;
            rd.L = L;
            rd.R = R;
            rd.axL = axL;
            rd.axR = axR;
            rd.mL = mL;
            rd.mR = mR;
            rd.p0L = planePoint(axL);
            rd.p0R = planePoint(axR);
            rd.g = g;
            rd.wedge = wedge;
            rd.cylNbr = cyl;
            valid.push_back(rd);
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
            r.id = 0;
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

            for (int k = 0; k < rd.cand.n; ++k)
                used[rd.cand.mem[k]] = 1;

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
            if (!(passSpread && passG4 && passG5 && passChains)) {
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

            work.accepted.push_back(std::move(r));
        }

        sortAcceptedRejected(mv, work);
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
