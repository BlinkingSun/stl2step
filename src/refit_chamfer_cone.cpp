// stl2step P1 stage C — chamfer-cone claim (1.3.0 detector C).
// After B1 / NgonWall, before C1. gp_ / stdlib only. Never throws.
// Always returns true (stage ok). Silent on every skip.
//
// Predicates (all required):
//   1. Unclaimed provisional ring (flood of constant-slope provisionals).
//   2. N ≥ 6 (shared-edge verts with the cylinder, else azimuth clusters).
//   3. Adjacent to one accepted closed360 Cylinder (NgonWall or CylGrow).
//   4. Adjacent to one Unclaimed or already-accepted Plane cap (|n·axis|≈1).
//   5. Shared axis with that cylinder (inherited; gaussMapAxis rejects cones).
//   6. R(v) linear (least-squares ρ = r0 + k v, residual ≤ sagTol).
//   7. Two distinct end radii; one equals cyl R ± sagTol.
//   8. Constant slope: |n·axis| ≈ sin α, α typically ~45° (window 25°–65°).
//   9. Dihedral to the cylinder neighbour is not a round:
//        acute crease φ < thetaCylLo → G1 / seed-band round, skip;
//        |n·axis| in B1's lateral cylinder band → leftover cyl strip, skip.
//        A 45° chamfer crease may land in [thetaCylLo, thetaCylHi]; the
//        |n·axis| ≈ sin α window is what lets it through.
//  10. Not a turning-axis round: in-ring n_i×n_j is ∥ cylinder axis (cone),
//      not meridional (torus / |S| fillet path).
//
// Do not claim: S02-class rounds, S11-b plane–plane 45° cube chamfers
// (no closed360 cylinder neighbour → C1 keeps 12 FilletStrip cylinders).
//
// On accept: SurfType::Cone, Origin::ChamferCone. Provisionals marked
// ConsumedCylinder so C1 (fillet.cpp) treats them as ineligible.
//
// Region fields:
//   ax            cylinder axis; Location at the R_lo ring (cyl-match if that
//                 end is the smaller R), XDir = lowest-id mesh vertex azimuth
//   radius        R_lo (smaller of the two end radii; OCCT cone R at Location)
//   vMin, vMax    the two radii R_lo, R_hi  (Region has no second radius field)
//   maxVertexDev  max |ρ − (r0+k v)| residual (I4-honest)
//   rmsVertexDev  rms of that residual
//   uMin, uMax    0, 2π  (closed ring)
//   closed360     true
//   nSides        N
//   chordSagitta  chordSagitta(R_lo, N)
//   filletNbrA    work.accepted index of the neighbour cylinder (pre-D)
//   outwardNormal cylinder-style σ (material on the axis side)
//   maxVertexSnap frustum height along the axis (mm), SIGNED: the axial offset
//                 from the R_lo rim to the R_hi rim measured along +Direction,
//                 so it is negative when the taper runs against the canonical
//                 axis. The axis itself is never flipped -- see the comment at
//                 the axOut/hSigned site below.
//
// SPDX-License-Identifier: MIT

#include "refit_cone_math.hpp"
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

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTiny = 1e-30;

// α window around a typical 45° hole-mouth chamfer.
constexpr double kSinALo = 0.42261826174069944;  // sin(25°)
constexpr double kSinAHi = 0.9063077870366499;   // sin(65°)
constexpr double kCapDot = 0.9659258262890683;   // cos(15°)
constexpr double kTurnDot = 0.7;                 // |(n×n)·axis| / |n×n| for cone

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
}

// Median of samples weighted by how much of the measured quantity each one
// actually carries. The turning-axis test reads the DIRECTION of n_i x n_j; for
// unit normals |n_i x n_j| is the sine of the fold, so a pair with a small
// modulus has measured almost nothing and its direction is correspondingly
// uncertain (the relative error of that direction goes as normal-noise/modulus).
// Weighting each sample by its own modulus is therefore weighting by inverse
// angular uncertainty -- a measurement statement, not a threshold. No sample is
// discarded and there is no cut-off to tune.
double weightedMedianInPlace(std::vector<std::pair<double, double>>& vw) {
    if (vw.empty()) return 0.0;
    std::sort(vw.begin(), vw.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                  if (a.second != b.second) return a.second < b.second;
                  return a.first < b.first;
              });
    double total = 0.0;
    for (const auto& e : vw) total += e.first;
    if (!(total > 0.0)) return vw[vw.size() / 2].second;
    double acc = 0.0;
    for (const auto& e : vw) {
        acc += e.first;
        if (acc >= 0.5 * total) return e.second;
    }
    return vw.back().second;
}

double medianInPlace(std::vector<double>& v) {
    const std::size_t n = v.size();
    if (n == 0) return 0.0;
    std::sort(v.begin(), v.end());
    if (n & 1u) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

gp_XYZ triVert(const MeshView& mv, int lt, int k) {
    const int gt = mv.compTris ? mv.compTris[lt] : lt;
    return mv.pts[mv.tris[gt][k]];
}

double triAreaNormal(const MeshView& mv, int lt, gp_XYZ& nOut) {
    const gp_XYZ a = triVert(mv, lt, 0);
    const gp_XYZ b = triVert(mv, lt, 1);
    const gp_XYZ c = triVert(mv, lt, 2);
    const gp_XYZ n = (b - a).Crossed(c - a);
    const double mag = n.Modulus();
    if (mag > kTiny) nOut = n / mag;
    else nOut = gp_XYZ(0, 0, 1);
    return 0.5 * mag;
}

gp_XYZ triCentroid(const MeshView& mv, int lt) {
    return (triVert(mv, lt, 0) + triVert(mv, lt, 1) + triVert(mv, lt, 2)) / 3.0;
}

int edgeOtherTri(const std::vector<std::array<int, 2>>& edgeTris, int e, int t) {
    if (e < 0 || e >= (int)edgeTris.size()) return -1;
    const int a = edgeTris[e][0], b = edgeTris[e][1];
    if (a == t) return b;
    if (b == t) return a;
    return -1;
}

int cornerLocal(const std::vector<int>& g2l, const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    const int g = mv.tris[gt][k];
    if (g < 0 || g >= (int)g2l.size()) return -1;
    return g2l[g];
}

gp_XYZ localPt(const MeshView& mv, int lv) {
    return mv.pts[mv.compVtx ? mv.compVtx[lv] : lv];
}

void canonicalizeDir(gp_XYZ& a) {
    if (a.Z() < -1e-15
        || (std::abs(a.Z()) <= 1e-15 && a.Y() < -1e-15)
        || (std::abs(a.Z()) <= 1e-15 && std::abs(a.Y()) <= 1e-15 && a.X() < 0.0)) {
        a.Reverse();
    }
}

double rhoOf(const gp_XYZ& p, const gp_XYZ& loc, const gp_XYZ& axis) {
    return axis.Crossed(p - loc).Modulus();
}

double vOf(const gp_XYZ& p, const gp_XYZ& loc, const gp_XYZ& axis) {
    return (p - loc).Dot(axis);
}

bool isClosedCyl(const Region& r) {
    if (r.type != SurfType::Cylinder || !r.closed360) return false;
    return r.origin == Origin::NgonWall || r.origin == Origin::CylGrow;
}

double sagTolOf(const Region& cyl, int nRing, const DerivedTols& tol, const MeshView& mv) {
    const int nUse = std::max(cyl.nSides, std::max(nRing, 3));
    const double sag = std::max(cyl.chordSagitta, chordSagitta(cyl.radius, nUse));
    return std::max({tol.epsMesh, tol.epsPlane, mv.sewTol, sag});
}

int nBandsAzimuth(std::vector<double>& psi, double thetaBin) {
    if (psi.empty()) return 0;
    std::sort(psi.begin(), psi.end());
    const int n = (int)psi.size();
    if (n == 1) return 1;
    std::vector<double> gaps((std::size_t)n);
    for (int i = 0; i < n - 1; ++i)
        gaps[(std::size_t)i] = psi[(std::size_t)i + 1] - psi[(std::size_t)i];
    gaps[(std::size_t)n - 1] = (psi[0] + kTwoPi) - psi[(std::size_t)n - 1];
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
    double prev = psi[(std::size_t)start];
    for (int k = 1; k < n; ++k) {
        const int idx = (start + k) % n;
        const double cur = psi[(std::size_t)idx];
        double gap = cur - prev;
        if (gap < 0.0) gap += kTwoPi;
        if (gap > th) ++bands;
        prev = cur;
    }
    return bands;
}

void sortAccepted(const MeshView& mv, SegmentWork& work) {
    auto areaOf = [&](const Region& r) {
        double a = 0.0;
        std::vector<int> ord = r.tris;
        std::sort(ord.begin(), ord.end());
        for (int t : ord) {
            gp_XYZ n;
            a += triAreaNormal(mv, t, n);
        }
        return a;
    };
    std::sort(work.accepted.begin(), work.accepted.end(),
              [&](const Region& a, const Region& b) {
                  const double aa = areaOf(a), ba = areaOf(b);
                  if (aa != ba) return aa > ba;
                  const int am = minTriOf(a.tris), bm = minTriOf(b.tris);
                  if (am != bm) return am < bm;
                  return a.origin < b.origin;
              });
}

}  // namespace

bool claimChamferConesC(const MeshView& mv, const SegmentParams& /*p*/,
                        const DerivedTols& tol, SegmentWork& work) {
    try {
        if (!mv.pts || !mv.tris || mv.nTri == 0) return true;
        if (work.provisionals.empty() || work.accepted.empty()) return true;
        if (!mv.triEdges) return true;

        const int nProv = (int)work.provisionals.size();
        const int nTri = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const int nVtx = (int)mv.nVtx;

        std::vector<int> cylIds;
        for (int i = 0; i < (int)work.accepted.size(); ++i) {
            if (isClosedCyl(work.accepted[(std::size_t)i])) cylIds.push_back(i);
        }
        // Diag only (STL2STEP_DIAG_130): NAME the predicate that refuses each
        // candidate ring. Prints nothing and decides nothing when off; no
        // threshold, window or predicate is changed by any of it.
        const char* d130 = std::getenv("STL2STEP_DIAG_130");
        const bool cDiag = d130 && d130[0] && d130[0] != '0';
        if (cDiag) {
            int nClosed = 0;
            for (const Region& a : work.accepted)
                if (isClosedCyl(a)) nClosed++;
            std::fprintf(stderr,
                         "DIAG_C_ENTRY nProv=%d nAcc=%d closed360Cyls=%d\n",
                         (int)work.provisionals.size(), (int)work.accepted.size(), nClosed);
        }
        if (cylIds.empty()) {
            if (cDiag)
                std::fprintf(stderr, "DIAG_C_NONE why=no-closed360-cylinder-neighbour\n");
            return true;
        }
        std::sort(cylIds.begin(), cylIds.end(), [&](int a, int b) {
            const int ma = minTriOf(work.accepted[(std::size_t)a].tris);
            const int mb = minTriOf(work.accepted[(std::size_t)b].tris);
            if (ma != mb) return ma < mb;
            return a < b;
        });

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
                if (t >= 0 && t < nTri) triProv[(std::size_t)t] = i;
            }
        }

        std::vector<int> triAcc((std::size_t)nTri, -1);
        for (int i = 0; i < (int)work.accepted.size(); ++i) {
            std::sort(work.accepted[i].tris.begin(), work.accepted[i].tris.end());
            for (int t : work.accepted[i].tris) {
                if (t >= 0 && t < nTri) triAcc[(std::size_t)t] = i;
            }
        }

        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{-1, -1});
        for (int t = 0; t < nTri; ++t) {
            for (int s = 0; s < 3; ++s) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || e >= nEdge) continue;
                if (edgeTris[(std::size_t)e][0] < 0) edgeTris[(std::size_t)e][0] = t;
                else if (edgeTris[(std::size_t)e][1] < 0) edgeTris[(std::size_t)e][1] = t;
            }
        }

        std::vector<std::vector<int>> adj((std::size_t)nProv);
        for (int e = 0; e < nEdge; ++e) {
            const int ta = edgeTris[(std::size_t)e][0];
            const int tb = edgeTris[(std::size_t)e][1];
            if (ta < 0 || tb < 0) continue;
            const int pa = triProv[(std::size_t)ta], pb = triProv[(std::size_t)tb];
            if (pa < 0 || pb < 0 || pa == pb) continue;
            adj[(std::size_t)pa].push_back(pb);
            adj[(std::size_t)pb].push_back(pa);
        }
        for (int i = 0; i < nProv; ++i) {
            auto& v = adj[(std::size_t)i];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }

        auto provSlopeAbs = [&](int pi, const gp_XYZ& axis) -> double {
            const Provisional& pr = work.provisionals[(std::size_t)pi];
            gp_XYZ acc(0, 0, 0);
            double area = 0.0;
            for (int t : pr.tris) {
                gp_XYZ n;
                const double a = triAreaNormal(mv, t, n);
                acc += n * a;
                area += a;
            }
            if (area <= kTiny) {
                const gp_XYZ n = pr.plane.Direction().XYZ();
                return std::abs(n.Dot(axis));
            }
            const double m = acc.Modulus();
            if (m <= kTiny) return 0.0;
            return std::abs((acc / m).Dot(axis));
        };

        auto inSlopeWindow = [&](double s) {
            return s >= kSinALo - 1e-12 && s <= kSinAHi + 1e-12;
        };

        bool anyCommit = false;

        for (int cylI : cylIds) {
            const Region& cyl = work.accepted[(std::size_t)cylI];
            gp_XYZ axis = cyl.ax.Direction().XYZ();
            const double am = axis.Modulus();
            if (am < 1e-15) continue;
            axis /= am;
            canonicalizeDir(axis);
            const gp_XYZ cylLoc = cyl.ax.Location().XYZ();

            std::vector<char> touchCyl((std::size_t)nProv, 0);
            for (int t : cyl.tris) {
                if (t < 0 || t >= nTri) continue;
                for (int s = 0; s < 3; ++s) {
                    const int e = mv.triEdges[t][s];
                    const int ot = edgeOtherTri(edgeTris, e, t);
                    if (ot < 0 || ot >= nTri) continue;
                    const int pi = triProv[(std::size_t)ot];
                    if (pi < 0) continue;
                    if (work.provisionals[(std::size_t)pi].claim != ProvClaim::Unclaimed)
                        continue;
                    touchCyl[(std::size_t)pi] = 1;
                }
            }

            // Diag only: the ring's provisionals as they REACH C -- claim state,
            // slope, and the accepted region that consumed them. This is the
            // measurement that separates "A2 never linked them" from "an earlier
            // stage took them away". Decides nothing.
            if (cDiag) {
                for (int i = 0; i < nProv; ++i) {
                    const Provisional& pr = work.provisionals[(std::size_t)i];
                    if (pr.tris.empty()) continue;
                    const double sl = provSlopeAbs(i, axis);
                    if (!inSlopeWindow(sl)) continue;
                    int owner = -1;
                    for (int t : pr.tris) {
                        if (t >= 0 && t < nTri && triAcc[(std::size_t)t] >= 0) {
                            owner = triAcc[(std::size_t)t];
                            break;
                        }
                    }
                    std::fprintf(stderr,
                                 "DIAG_C_PROV cyl=%d p=%d nTri=%d minTri=%d claim=%d "
                                 "slope=%.6f touch=%d owner=%d\n",
                                 cylI, i, (int)pr.tris.size(), minTriOf(pr.tris),
                                 (int)pr.claim, sl, touchCyl[(std::size_t)i] ? 1 : 0,
                                 owner);
                }
            }

            std::vector<int> seeds;
            for (int i = 0; i < nProv; ++i) {
                if (!touchCyl[(std::size_t)i]) continue;
                if (work.provisionals[(std::size_t)i].claim != ProvClaim::Unclaimed)
                    continue;
                if (work.provisionals[(std::size_t)i].tris.empty()) continue;
                if (!inSlopeWindow(provSlopeAbs(i, axis))) continue;
                seeds.push_back(i);
            }
            if (seeds.empty()) continue;
            std::sort(seeds.begin(), seeds.end(), [&](int a, int b) {
                const int ma = minTriOf(work.provisionals[(std::size_t)a].tris);
                const int mb = minTriOf(work.provisionals[(std::size_t)b].tris);
                if (ma != mb) return ma < mb;
                return a < b;
            });

            std::vector<char> used((std::size_t)nProv, 0);

            // Diag only: name which predicate refused a candidate ring.
            auto crej = [&](int sd, const char* why, double a, double b) {
                if (cDiag)
                    std::fprintf(stderr,
                                 "DIAG_C_RING cyl=%d cylR=%.4f closed=%d seed=%d why=%s "
                                 "a=%.6g b=%.6g\n",
                                 cylI, cyl.radius, cyl.closed360 ? 1 : 0, sd, why, a, b);
            };
            if (cDiag)
                std::fprintf(stderr, "DIAG_C_CYL cyl=%d R=%.4f nTri=%d nSides=%d seeds=%d\n",
                             cylI, cyl.radius, (int)cyl.tris.size(), cyl.nSides,
                             (int)seeds.size());

            for (int seed : seeds) {
                if (used[(std::size_t)seed]) continue;
                if (work.provisionals[(std::size_t)seed].claim != ProvClaim::Unclaimed)
                    continue;

                std::vector<int> mem;
                std::vector<int> stack = {seed};
                used[(std::size_t)seed] = 1;
                while (!stack.empty()) {
                    const int u = stack.back();
                    stack.pop_back();
                    mem.push_back(u);
                    for (int v : adj[(std::size_t)u]) {
                        // Diag only: name why a shared-edge neighbour of a ring
                        // provisional is not linked into the flood.
                        auto lrej = [&](const char* why, double a, double b) {
                            if (cDiag)
                                std::fprintf(stderr,
                                             "DIAG_C_LINK cyl=%d seed=%d u=%d v=%d why=%s "
                                             "a=%.6g b=%.6g\n",
                                             cylI, seed, u, v, why, a, b);
                        };
                        if (used[(std::size_t)v]) continue;
                        if (work.provisionals[(std::size_t)v].claim != ProvClaim::Unclaimed) {
                            lrej("not-unclaimed",
                                 (double)(int)work.provisionals[(std::size_t)v].claim, 0.0);
                            continue;
                        }
                        if (work.provisionals[(std::size_t)v].tris.empty()) {
                            lrej("empty-provisional", 0.0, 0.0);
                            continue;
                        }
                        if (!inSlopeWindow(provSlopeAbs(v, axis))) {
                            lrej("slope-window", provSlopeAbs(v, axis), kSinALo);
                            continue;
                        }
                        used[(std::size_t)v] = 1;
                        stack.push_back(v);
                    }
                }
                std::sort(mem.begin(), mem.end());
                mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
                if (mem.empty()) { crej(seed, "empty-flood", 0, 0); continue; }

                std::vector<char> inMem((std::size_t)nProv, 0);
                for (int m : mem) inMem[(std::size_t)m] = 1;

                std::vector<int> sTris;
                for (int m : mem) {
                    const auto& ts = work.provisionals[(std::size_t)m].tris;
                    sTris.insert(sTris.end(), ts.begin(), ts.end());
                }
                std::sort(sTris.begin(), sTris.end());
                sTris.erase(std::unique(sTris.begin(), sTris.end()), sTris.end());
                if (sTris.size() < 2) { crej(seed, "under-2-triangles", (double)sTris.size(), 0); continue; }

                // Shared-edge vertices with the cylinder → N.
                std::vector<int> juncVerts;
                bool touchesThisCyl = false;
                for (int t : sTris) {
                    for (int s = 0; s < 3; ++s) {
                        const int e = mv.triEdges[t][s];
                        const int ot = edgeOtherTri(edgeTris, e, t);
                        if (ot < 0 || ot >= nTri) continue;
                        if (triAcc[(std::size_t)ot] != cylI) continue;
                        touchesThisCyl = true;
                        if (mv.compEdges && e >= 0 && e < nEdge) {
                            const int lo = mv.compEdges[e].first;
                            const int hi = mv.compEdges[e].second;
                            if (lo >= 0) juncVerts.push_back(lo);
                            if (hi >= 0) juncVerts.push_back(hi);
                        } else {
                            for (int k = 0; k < 3; ++k) {
                                const int lv = cornerLocal(g2l, mv, t, k);
                                if (lv >= 0) juncVerts.push_back(lv);
                            }
                        }
                    }
                }
                if (!touchesThisCyl) { crej(seed, "no-shared-edge-with-this-cylinder", (double)sTris.size(), 0); continue; }
                std::sort(juncVerts.begin(), juncVerts.end());
                juncVerts.erase(std::unique(juncVerts.begin(), juncVerts.end()),
                                juncVerts.end());

                // Cap: Unclaimed provisional or accepted Plane, |n·axis| ≈ 1.
                bool hasCap = false;
                for (int t : sTris) {
                    for (int s = 0; s < 3; ++s) {
                        const int e = mv.triEdges[t][s];
                        const int ot = edgeOtherTri(edgeTris, e, t);
                        if (ot < 0 || ot >= nTri) continue;
                        const int op = triProv[(std::size_t)ot];
                        if (op >= 0 && inMem[(std::size_t)op]) continue;
                        gp_XYZ nCap(0, 0, 0);
                        bool capCand = false;
                        if (op >= 0
                            && work.provisionals[(std::size_t)op].claim
                                   == ProvClaim::Unclaimed
                            && !inMem[(std::size_t)op]) {
                            nCap = work.provisionals[(std::size_t)op].plane.Direction().XYZ();
                            capCand = true;
                        } else {
                            const int ai = triAcc[(std::size_t)ot];
                            if (ai >= 0 && ai != cylI
                                && work.accepted[(std::size_t)ai].type == SurfType::Plane) {
                                nCap = work.accepted[(std::size_t)ai].ax.Direction().XYZ();
                                capCand = true;
                            }
                        }
                        if (!capCand) continue;
                        const double nm = nCap.Modulus();
                        if (nm < 1e-15) continue;
                        nCap /= nm;
                        if (std::abs(nCap.Dot(axis)) >= kCapDot) {
                            hasCap = true;
                            break;
                        }
                    }
                    if (hasCap) break;
                }
                if (!hasCap) { crej(seed, "no-cap-plane-neighbour", (double)sTris.size(), kCapDot); continue; }

                // Unique verts, ascending local id (I5).
                std::vector<std::pair<int, gp_XYZ>> uniq;
                for (int t : sTris) {
                    for (int k = 0; k < 3; ++k) {
                        const int lv = cornerLocal(g2l, mv, t, k);
                        if (lv < 0) continue;
                        uniq.emplace_back(lv, localPt(mv, lv));
                    }
                }
                std::sort(uniq.begin(), uniq.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });
                uniq.erase(std::unique(uniq.begin(), uniq.end(),
                                       [](const auto& a, const auto& b) {
                                           return a.first == b.first;
                                       }),
                           uniq.end());
                if (uniq.size() < 6) { crej(seed, "under-6-unique-vertices", (double)uniq.size(), 0); continue; }

                // N: junction verts preferred, else azimuth clusters, else |mem|.
                std::vector<double> chi;
                chi.reserve(uniq.size());
                {
                    const gp_XYZ tmp =
                        (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
                    gp_XYZ xD = axis.Crossed(tmp);
                    const double xm = xD.Modulus();
                    if (xm > 1e-15) xD /= xm;
                    const gp_XYZ yD = axis.Crossed(xD);
                    for (const auto& uv : uniq) {
                        const gp_XYZ d = uv.second - cylLoc;
                        const gp_XYZ rad = d - axis * d.Dot(axis);
                        if (rad.SquareModulus() < 1e-24) continue;
                        chi.push_back(std::atan2(rad.Dot(yD), rad.Dot(xD)));
                    }
                }
                const int nAz = nBandsAzimuth(chi, tol.thetaBin);
                int nSides = (int)juncVerts.size();
                if (nSides < DerivedTols::kG5NSidesMin) nSides = nAz;
                if (nSides < DerivedTols::kG5NSidesMin) nSides = (int)mem.size();
                if (nSides < DerivedTols::kG5NSidesMin) { crej(seed, "nSides-under-floor", (double)nSides, (double)DerivedTols::kG5NSidesMin); continue; }

                // Closed ring: largest azimuth gap vs 2π/N. chi is sorted.
                if (chi.size() < 3) { crej(seed, "under-3-azimuths", (double)chi.size(), 0); continue; }
                {
                    double maxGap = 0.0;
                    for (std::size_t i = 0; i + 1 < chi.size(); ++i)
                        maxGap = std::max(maxGap, chi[i + 1] - chi[i]);
                    maxGap = std::max(maxGap, (chi.front() + kTwoPi) - chi.back());
                    const double bandArc = kTwoPi / (double)std::max(nSides, 3);
                    if (maxGap > 1.5 * bandArc + 1e-12) { crej(seed, "ring-not-closed", maxGap, 1.5 * bandArc); continue; }
                }

                const double sagTol = sagTolOf(cyl, nSides, tol, mv);

                // Constant slope |n·axis| ≈ sin α.
                std::vector<double> slopes;
                slopes.reserve(sTris.size());
                double slopeMin = 1e300, slopeMax = -1e300;
                for (int t : sTris) {
                    gp_XYZ n;
                    triAreaNormal(mv, t, n);
                    const double s = std::abs(n.Dot(axis));
                    slopes.push_back(s);
                    slopeMin = std::min(slopeMin, s);
                    slopeMax = std::max(slopeMax, s);
                }
                const double slopeMed = medianInPlace(slopes);
                if (!inSlopeWindow(slopeMed)) { crej(seed, "slope-outside-25-65deg-window", slopeMed, 0); continue; }
                if (slopeMax - slopeMin > 0.12) { crej(seed, "slope-not-constant", slopeMax - slopeMin, 0.12); continue; }  // torus / round

                // Dihedral to cylinder neighbour is not a round.
                // Acute crease φ < thetaCylLo: G1 seed-band round.
                // |n·axis| in B1 lateral band (≲ sin 3° / thetaCylLo): leftover cyl.
                {
                    std::vector<double> localPhi;
                    for (int t : sTris) {
                        gp_XYZ nR;
                        triAreaNormal(mv, t, nR);
                        for (int s = 0; s < 3; ++s) {
                            const int e = mv.triEdges[t][s];
                            const int ot = edgeOtherTri(edgeTris, e, t);
                            if (ot < 0 || ot >= nTri || triAcc[(std::size_t)ot] != cylI)
                                continue;
                            gp_XYZ nC;
                            triAreaNormal(mv, ot, nC);
                            const double d = std::min(1.0, std::max(-1.0, nC.Dot(nR)));
                            double phi = std::acos(d);
                            if (phi > kPi * 0.5) phi = kPi - phi;
                            localPhi.push_back(phi);
                        }
                    }
                    if (localPhi.empty()) { crej(seed, "no-dihedral-to-cylinder", 0, 0); continue; }
                    const double phiMed = medianInPlace(localPhi);
                    const double loEps = std::max(1e-12, 8.0 * std::numeric_limits<double>::epsilon()
                                                             * std::max(1.0, tol.thetaCylLo));
                    if (phiMed < tol.thetaCylLo - loEps) { crej(seed, "crease-below-thetaCylLo-tangent-round", phiMed, tol.thetaCylLo); continue; }
                    const double lat = std::max(tol.gaussAxisTiltSin(),
                                                std::sin(tol.thetaCylLo));
                    if (slopeMed <= lat + 1e-12) { crej(seed, "slope-inside-B1-lateral-band", slopeMed, lat); continue; }
                }

                // Turning-axis round: in-ring n_i×n_j must be ∥ cylinder axis.
                {
                    // Each in-ring pair contributes its alignment WEIGHTED by
                    // |n_i x n_j|, the sine of the fold it realises. Half the
                    // adjacent pairs of a quad-tessellated ring are the two
                    // triangles of ONE planar quad: they measure no fold at all
                    // and their cross product is float32 round-off with an
                    // arbitrary direction. Measured on the plate's chamfer ring
                    // the two populations do not overlap -- 92 real folds at
                    // |cross| = 0.0483 with alignment 0.7071026..0.7071100 (the
                    // exact 1/sqrt(2) a 45 deg cone must give), and 92 round-off
                    // pairs at |cross| in [3.7e-8, 2.7e-5] carrying 1.0e-4 of the
                    // total weight and an arbitrary 0.649. An unweighted median
                    // of the two is 0.6779 and refuses a perfect cone; the
                    // weighted median is 0.707106781. kTurnDot is NOT touched.
                    std::vector<std::pair<double, double>> align;
                    for (int t : sTris) {
                        gp_XYZ nA;
                        triAreaNormal(mv, t, nA);
                        for (int s = 0; s < 3; ++s) {
                            const int e = mv.triEdges[t][s];
                            const int ot = edgeOtherTri(edgeTris, e, t);
                            if (ot < 0 || ot >= nTri) continue;
                            const int op = triProv[(std::size_t)ot];
                            if (op < 0 || !inMem[(std::size_t)op]) continue;
                            if (ot <= t) continue;
                            gp_XYZ nB;
                            triAreaNormal(mv, ot, nB);
                            const gp_XYZ cr = nA.Crossed(nB);
                            const double cm = cr.Modulus();
                            if (cm < 1e-8) continue;
                            align.emplace_back(cm, std::abs(cr.Dot(axis)) / cm);
                        }
                    }
                    if (!align.empty()) {
                        const double al = weightedMedianInPlace(align);
                        if (al < kTurnDot) { crej(seed, "turning-axis-round", al, kTurnDot); continue; }
                    }
                }

                // R(v) linear about the cylinder axis.
                double vMin = 1e300, vMax = -1e300;
                double sumV = 0, sumR = 0, sumVV = 0, sumVR = 0;
                const int nU = (int)uniq.size();
                for (const auto& uv : uniq) {
                    const double v = vOf(uv.second, cylLoc, axis);
                    const double r = rhoOf(uv.second, cylLoc, axis);
                    vMin = std::min(vMin, v);
                    vMax = std::max(vMax, v);
                    sumV += v;
                    sumR += r;
                    sumVV += v * v;
                    sumVR += v * r;
                }
                const double det = (double)nU * sumVV - sumV * sumV;
                if (!(std::abs(det) > 1e-18) || !(vMax > vMin + 1e-15)) { crej(seed, "degenerate-Rv-system", det, vMax - vMin); continue; }
                const double kSlope = ((double)nU * sumVR - sumV * sumR) / det;
                const double r0 = (sumR - kSlope * sumV) / (double)nU;
                double maxLin = 0.0, sumSq = 0.0;
                for (const auto& uv : uniq) {
                    const double v = vOf(uv.second, cylLoc, axis);
                    const double r = rhoOf(uv.second, cylLoc, axis);
                    const double pred = r0 + kSlope * v;
                    const double d = std::abs(r - pred);
                    maxLin = std::max(maxLin, d);
                    sumSq += d * d;
                }
                if (maxLin > sagTol + 1e-12) { crej(seed, "Rv-not-linear", maxLin, sagTol); continue; }  // circular R(v) = round

                double R1 = r0 + kSlope * vMin;
                double R2 = r0 + kSlope * vMax;
                if (!(R1 > 0.0) || !(R2 > 0.0)) { crej(seed, "non-positive-end-radius", R1, R2); continue; }
                if (R1 > R2) std::swap(R1, R2);
                if (!(R2 - R1 > std::max(4.0 * sagTol, 3.0 * tol.epsPlane))) { crej(seed, "end-radii-not-distinct", R2 - R1, std::max(4.0 * sagTol, 3.0 * tol.epsPlane)); continue; }

                const double d1 = std::abs(R1 - cyl.radius);
                const double d2 = std::abs(R2 - cyl.radius);
                if (std::min(d1, d2) > sagTol + 1e-12) { crej(seed, "neither-end-matches-cylinder-R", std::min(d1, d2), sagTol); continue; }

                if (!(2.0 * tol.epsPlane < R1 && R2 < 2.0 * mv.diag)) { crej(seed, "radii-out-of-scale", R1, R2); continue; }

                // Place Location at the R_lo ring; XDir = lowest-id vertex.
                const bool loAtVmin = std::abs((r0 + kSlope * vMin) - R1)
                                      <= std::abs((r0 + kSlope * vMax) - R1);
                const double vLo = loAtVmin ? vMin : vMax;
                const gp_XYZ loc = cylLoc + axis * vLo;
                const double height = std::abs(vMax - vMin);

                int minLv = std::numeric_limits<int>::max();
                gp_XYZ ps(0, 0, 0);
                for (const auto& uv : uniq) {
                    if (uv.first < minLv) {
                        minLv = uv.first;
                        ps = uv.second;
                    }
                }
                // XDirection is INHERITED from the neighbour cylinder, not taken
                // from this ring's own lowest-id vertex. A frustum's seam edge has
                // to join a vertex on each rim, and the two rims' seam vertices
                // are chosen (seamVertexOf) in the frame of whichever region the
                // chain machinery names as the axis region -- the CYLINDER for
                // the rim they share, the CONE for the rim that meets a cap
                // plane. Two different frames put the two seam vertices at two
                // different azimuths and the seam is then not a generator at all:
                // measured on the plate, the low rim's vertex sat at azimuth 0
                // (the cylinder's XDirection) and the high rim's at 3.913 deg
                // (this ring's own lowest-id vertex), the wire came apart
                // (NotConnected / UnorientableShape) and the component reverted.
                // Sharing the cylinder's XDirection makes the two frames one, and
                // it is just as deterministic -- it is the cylinder's own.
                gp_XYZ xDir = cyl.ax.XDirection().XYZ();
                xDir -= axis * xDir.Dot(axis);
                if (xDir.Modulus() < 1e-12)
                    xDir = (ps - loc) - axis * (ps - loc).Dot(axis);
                if (xDir.Modulus() < 1e-12) {
                    const gp_XYZ tmp =
                        (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
                    xDir = axis.Crossed(tmp);
                }
                const double xdm = xDir.Modulus();
                if (xdm < 1e-15) { crej(seed, "degenerate-xdir", xdm, 0); continue; }
                xDir /= xdm;

                double sigma = 0.0;
                double areaSum = 0.0;
                for (int t : sTris) {
                    gp_XYZ n;
                    const double a = triAreaNormal(mv, t, n);
                    areaSum += a;
                    gp_XYZ rho = triCentroid(mv, t) - loc;
                    rho -= axis * rho.Dot(axis);
                    const double rm = rho.Modulus();
                    if (rm < 1e-15) continue;
                    rho /= rm;
                    sigma += a * n.Dot(rho);
                }

                // ONE FRAME FOR BOTH TAPER SIGNS. The cylinder axis this cone
                // inherits has been through canonicalizeDir and carries no
                // information about which way the taper runs, while every
                // consumer evaluates the frustum as
                // rho(v) = R_lo + (v/h)*(R_hi - R_lo) from Location along
                // +Direction. That sign is loAtVmin, and it is carried HERE, in
                // the SIGNED height, never by flipping the axis.
                //
                // Flipping the axis was measured to be wrong twice over. It
                // fixes the radius (the plane|cone rim came out at R 7.5 instead
                // of 9.5, 2.0 mm out, ratioSew 1114, 93 free edges) but it puts
                // the cone's frame ANTI-PARALLEL to the neighbour cylinder's:
                // gp_Ax2/gp_Ax3 are right-handed, so negating Direction while
                // keeping XDirection negates YDirection too, and the cone's u
                // then sweeps the opposite way round the shared rim circle --
                // which is built in the cylinder's frame. The cone's rim pcurve
                // is a 2d line of slope +1 in u either way, so on the flipped
                // copy it disagrees with the 3d circle by up to 2R = 17 mm.
                // BRepLib::SameParameter then cannot set the flag on an edge
                // whose OTHER pcurve is exact, and the neighbour cylinder reads
                // as invalid with no status of its own:
                //   DIAG_R2SUB rid=18 type=EDGE own=[] ctx=[InvalidSameParameterFlag]
                //              dev=0 sameP=0 p=(144.50000,-7.50000,-44.00000)
                // (the plate's LOWER chamfer rim, the copy whose taper runs
                // against +Z). Measured: with the axis flipped, cone rid 30
                // reaches faceSt=[27:UnorientableShape] and the component
                // reverts; with the signed height it reaches valid=1.
                //
                // The signed height is exact: gp_Cone accepts a negative
                // SemiAngle (its guard is on |Ang|), so
                // Ang = atan((R_hi - R_lo)/h) < 0 describes a frustum whose
                // radius SHRINKS along +Direction with no change of frame, and
                // rho(v) = R_lo + (v/h)*(R_hi - R_lo) is unchanged for both
                // signs of h.
                const gp_XYZ axOut = axis;
                const double hSigned = loAtVmin ? height : -height;

                // D-130-11(2): a cone predicts its facet-to-surface volume the
                // way a cylinder does. dVolPredicted = 0 on a built ChamferCone
                // made the cascade budget and D4's guard read the face as
                // zero-area (I9 FAIL region30/31 on the plate).
                //
                // Closed form, summed over the triangles this region CLAIMS --
                // not over the whole frustum -- because a detector-C ring may
                // hold only part of the skin. The density is exact per triangle:
                // kappa is affine in the axial coordinate and the axial
                // coordinate is affine over a planar facet, so
                // area(t) * kappa(rho at centroid(t)) IS the integral over t,
                // and over the complete skin the sum reproduces
                // coneFrustumChordVolume. The sign is I9's sigma, the same
                // outwardNormal that signs dVolCylinderSector.
                double dVolSum = 0.0;
                for (int t : sTris) {
                    gp_XYZ nT;
                    const double a = triAreaNormal(mv, t, nT);
                    if (!(a > 0.0)) continue;
                    const double u = (triCentroid(mv, t) - loc).Dot(axis);
                    const double rhoC = R1 + (u / hSigned) * (R2 - R1);
                    dVolSum += a * coneFrustumDVolDensity(rhoC, R1, R2, hSigned, nSides);
                }
                if (!std::isfinite(dVolSum)) dVolSum = 0.0;

                Region r;
                r.id = 0;
                r.type = SurfType::Cone;
                r.origin = Origin::ChamferCone;
                r.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(axOut), gp_Dir(xDir));
                r.radius = R1;
                r.uMin = 0.0;
                r.uMax = kTwoPi;
                r.vMin = R1;   // two radii (no second radius field)
                r.vMax = R2;
                r.closed360 = true;
                r.outwardNormal = sigma > 0.0;
                r.tris = sTris;
                r.maxVertexDev = maxLin;
                r.rmsVertexDev = std::sqrt(sumSq / (double)nU);
                r.chordSagitta = chordSagitta(R1, nSides);
                r.nSides = nSides;
                r.dVolPredicted = (sigma > 0.0 ? 1.0 : -1.0) * dVolSum;
                r.maxVertexSnap = hSigned;
                r.reject = Reject::None;
                r.builtAs = BuiltAs::NotBuilt;
                r.filletNbrA = cylI;
                r.filletNbrB = -1;
                (void)areaSum;

                for (int m : mem)
                    work.provisionals[(std::size_t)m].claim = ProvClaim::ConsumedCylinder;

                if (cDiag)
                    std::fprintf(stderr,
                                 "DIAG_C_ACCEPT cyl=%d cylR=%.4f seed=%d R1=%.6f R2=%.6f "
                                 "nSides=%d nTri=%d maxLin=%.3g height=%.6f\n",
                                 cylI, cyl.radius, seed, R1, R2, nSides, (int)sTris.size(),
                                 maxLin, height);
                work.accepted.push_back(std::move(r));
                anyCommit = true;
            }
        }

        if (anyCommit) sortAccepted(mv, work);
        return true;
    } catch (...) {
        return true;
    }
}

}  // namespace refit
}  // namespace stl2step
