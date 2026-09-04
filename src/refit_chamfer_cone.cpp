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
//   maxVertexSnap frustum height along the axis (mm)
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
        if (cylIds.empty()) return true;
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
                        if (used[(std::size_t)v]) continue;
                        if (work.provisionals[(std::size_t)v].claim != ProvClaim::Unclaimed)
                            continue;
                        if (work.provisionals[(std::size_t)v].tris.empty()) continue;
                        if (!inSlopeWindow(provSlopeAbs(v, axis))) continue;
                        used[(std::size_t)v] = 1;
                        stack.push_back(v);
                    }
                }
                std::sort(mem.begin(), mem.end());
                mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
                if (mem.empty()) continue;

                std::vector<char> inMem((std::size_t)nProv, 0);
                for (int m : mem) inMem[(std::size_t)m] = 1;

                std::vector<int> sTris;
                for (int m : mem) {
                    const auto& ts = work.provisionals[(std::size_t)m].tris;
                    sTris.insert(sTris.end(), ts.begin(), ts.end());
                }
                std::sort(sTris.begin(), sTris.end());
                sTris.erase(std::unique(sTris.begin(), sTris.end()), sTris.end());
                if (sTris.size() < 2) continue;

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
                if (!touchesThisCyl) continue;
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
                if (!hasCap) continue;

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
                if (uniq.size() < 6) continue;

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
                if (nSides < DerivedTols::kG5NSidesMin) continue;

                // Closed ring: largest azimuth gap vs 2π/N. chi is sorted.
                if (chi.size() < 3) continue;
                {
                    double maxGap = 0.0;
                    for (std::size_t i = 0; i + 1 < chi.size(); ++i)
                        maxGap = std::max(maxGap, chi[i + 1] - chi[i]);
                    maxGap = std::max(maxGap, (chi.front() + kTwoPi) - chi.back());
                    const double bandArc = kTwoPi / (double)std::max(nSides, 3);
                    if (maxGap > 1.5 * bandArc + 1e-12) continue;
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
                if (!inSlopeWindow(slopeMed)) continue;
                if (slopeMax - slopeMin > 0.12) continue;  // not constant (torus / round)

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
                    if (localPhi.empty()) continue;
                    const double phiMed = medianInPlace(localPhi);
                    const double loEps = std::max(1e-12, 8.0 * std::numeric_limits<double>::epsilon()
                                                             * std::max(1.0, tol.thetaCylLo));
                    if (phiMed < tol.thetaCylLo - loEps) continue;  // tangent round
                    const double lat = std::max(tol.gaussAxisTiltSin(),
                                                std::sin(tol.thetaCylLo));
                    if (slopeMed <= lat + 1e-12) continue;  // cylinder-like band
                }

                // Turning-axis round: in-ring n_i×n_j must be ∥ cylinder axis.
                {
                    std::vector<double> align;
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
                            align.push_back(std::abs(cr.Dot(axis)) / cm);
                        }
                    }
                    if (!align.empty()) {
                        const double al = medianInPlace(align);
                        if (al < kTurnDot) continue;  // meridional turning axis
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
                if (!(std::abs(det) > 1e-18) || !(vMax > vMin + 1e-15)) continue;
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
                if (maxLin > sagTol + 1e-12) continue;  // circular R(v) = round

                double R1 = r0 + kSlope * vMin;
                double R2 = r0 + kSlope * vMax;
                if (!(R1 > 0.0) || !(R2 > 0.0)) continue;
                if (R1 > R2) std::swap(R1, R2);
                if (!(R2 - R1 > std::max(4.0 * sagTol, 3.0 * tol.epsPlane))) continue;

                const double d1 = std::abs(R1 - cyl.radius);
                const double d2 = std::abs(R2 - cyl.radius);
                if (std::min(d1, d2) > sagTol + 1e-12) continue;

                if (!(2.0 * tol.epsPlane < R1 && R2 < 2.0 * mv.diag)) continue;

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
                gp_XYZ xDir = (ps - loc) - axis * (ps - loc).Dot(axis);
                if (xDir.Modulus() < 1e-12) {
                    const gp_XYZ tmp =
                        (std::abs(axis.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
                    xDir = axis.Crossed(tmp);
                }
                const double xdm = xDir.Modulus();
                if (xdm < 1e-15) continue;
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

                Region r;
                r.id = 0;
                r.type = SurfType::Cone;
                r.origin = Origin::ChamferCone;
                r.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(axis), gp_Dir(xDir));
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
                r.dVolPredicted = 0.0;
                r.maxVertexSnap = height;
                r.reject = Reject::None;
                r.builtAs = BuiltAs::NotBuilt;
                r.filletNbrA = cylI;
                r.filletNbrB = -1;
                (void)areaSum;

                for (int m : mem)
                    work.provisionals[(std::size_t)m].claim = ProvClaim::ConsumedCylinder;

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
