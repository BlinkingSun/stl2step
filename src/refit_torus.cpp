// stl2step P1 stage T — toroidal-round claim (D-140-6 §1/§3(7), SPEC-torus-2).
//
// Consumes ProvClaim::Unclaimed provisionals and eligible Origin::CylGrow
// regions (maxVertexDev > tau) in Phase C's position: after C1, before A3.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

constexpr int kFitIters = 16;

bool torusDiagEnabled() {
    const char* a = std::getenv("STL2STEP_DIAG_130");
    const char* b = std::getenv("STL2STEP_DIAG_P2");
    return (a && a[0] && a[0] != '0') || (b && b[0] && b[0] != '0');
}

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
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
    if (mag > gp::Resolution()) nOut = n / mag;
    else nOut = gp_XYZ(0, 0, 1);
    return 0.5 * mag;
}

gp_XYZ triCentroid(const MeshView& mv, int lt) {
    return (triVert(mv, lt, 0) + triVert(mv, lt, 1) + triVert(mv, lt, 2)) / 3.0;
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

int edgeOtherTri(const std::vector<std::array<int, 2>>& edgeTris, int e, int t) {
    if (e < 0 || e >= (int)edgeTris.size()) return -1;
    const int a = edgeTris[(std::size_t)e][0], b = edgeTris[(std::size_t)e][1];
    if (a == t) return b;
    if (b == t) return a;
    return -1;
}

bool isClosedCyl(const Region& r) {
    if (r.type != SurfType::Cylinder || !r.closed360) return false;
    return r.origin == Origin::NgonWall || r.origin == Origin::CylGrow;
}

bool eligibleCylGrow(const Region& r, double tau) {
    return r.origin == Origin::CylGrow && r.maxVertexDev > tau;
}

bool axisDonor(const Region& r, double tau) {
    return isClosedCyl(r) && r.maxVertexDev <= tau;
}

struct ProfPt {
    double rho = 0, z = 0, psi = 0;
    int lv = -1;
};

struct TorusCert {
    bool ok = false;
    bool determinate = false;
    const char* why = "none";
    double Rmaj = 0, z0 = 0, Rmin = 0;
    double maxDev = 0;
    double lineRes = 0;
    int levels = 0, cols = 0, profPts = 0;
    double rmsDev = 0;
};

int countColumns(std::vector<double> psi, const std::vector<double>& rho, double tau) {
    const std::size_t n = psi.size();
    if (n == 0) return 0;
    double rhoMax = 0.0;
    for (double r : rho) rhoMax = std::max(rhoMax, r);
    if (!(rhoMax > gp::Resolution())) return 1;
    const double angTol = tau / rhoMax;
    std::sort(psi.begin(), psi.end());
    int cols = 1;
    double ref = psi[0];
    for (std::size_t i = 1; i < n; ++i) {
        if (psi[i] - ref > angTol) {
            ++cols;
            ref = psi[i];
        }
    }
    if (cols > 1 && (psi[0] + kTwoPi) - ref <= angTol) --cols;
    return cols;
}

int countLevels(std::vector<double> z, double tau) {
    if (z.empty()) return 0;
    std::sort(z.begin(), z.end());
    int lv = 1;
    double ref = z[0];
    for (std::size_t i = 1; i < z.size(); ++i) {
        if (z[i] - ref > tau) {
            ++lv;
            ref = z[i];
        }
    }
    return lv;
}

int countProfPts(const std::vector<ProfPt>& P, double tau) {
    std::vector<std::pair<double, double>> reps;
    for (const ProfPt& p : P) {
        bool found = false;
        for (const auto& q : reps)
            if (std::fabs(p.rho - q.first) <= tau && std::fabs(p.z - q.second) <= tau) {
                found = true;
                break;
            }
        if (!found) reps.emplace_back(p.rho, p.z);
    }
    return (int)reps.size();
}

double profileLineMaxResid(const std::vector<ProfPt>& P) {
    const std::size_t n = P.size();
    if (n < 2) return 0.0;
    double sx = 0, sy = 0;
    for (const ProfPt& p : P) { sx += p.rho; sy += p.z; }
    const double mx = sx / (double)n, my = sy / (double)n;
    double sxx = 0, sxy = 0, syy = 0;
    for (const ProfPt& p : P) {
        const double dx = p.rho - mx, dy = p.z - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    const double tr = sxx + syy, det = sxx * syy - sxy * sxy;
    const double disc = std::max(0.0, 0.25 * tr * tr - det);
    const double lmin = 0.5 * tr - std::sqrt(disc);
    double nx = sxy, ny = lmin - sxx;
    double nm = std::hypot(nx, ny);
    if (!(nm > gp::Resolution())) { nx = lmin - syy; ny = sxy; nm = std::hypot(nx, ny); }
    if (!(nm > gp::Resolution())) return 0.0;
    nx /= nm; ny /= nm;
    double m = 0.0;
    for (const ProfPt& p : P)
        m = std::max(m, std::fabs((p.rho - mx) * nx + (p.z - my) * ny));
    return m;
}

bool fitProfileCircle(const std::vector<ProfPt>& P, double tau, double& Rmaj, double& z0,
                      double& Rmin) {
    const std::size_t n = P.size();
    if (n < 3) return false;
    double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sz = 0, Sxz = 0, Syz = 0;
    for (const ProfPt& p : P) {
        const double x = p.rho, y = p.z, zq = x * x + y * y;
        Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y;
        Sz += zq; Sxz += x * zq; Syz += y * zq;
    }
    const double N = (double)n;
    const double a11 = 2.0 * (Sxx - Sx * Sx / N);
    const double a12 = 2.0 * (Sxy - Sx * Sy / N);
    const double a22 = 2.0 * (Syy - Sy * Sy / N);
    const double b1 = Sxz - Sx * Sz / N;
    const double b2 = Syz - Sy * Sz / N;
    const double det = a11 * a22 - a12 * a12;
    if (!(std::fabs(det) > gp::Resolution())) return false;
    double a = (b1 * a22 - b2 * a12) / det;
    double b = (a11 * b2 - a12 * b1) / det;
    double r = 0.0;
    for (const ProfPt& p : P) r += std::hypot(p.rho - a, p.z - b);
    r /= N;
    for (int it = 0; it < kFitIters; ++it) {
        double h11 = 0, h12 = 0, h13 = 0, h22 = 0, h23 = 0, h33 = 0;
        double g1 = 0, g2 = 0, g3 = 0;
        for (const ProfPt& p : P) {
            const double dx = p.rho - a, dy = p.z - b;
            const double d = std::hypot(dx, dy);
            if (!(d > gp::Resolution())) continue;
            const double j1 = -dx / d, j2 = -dy / d, j3 = -1.0;
            const double res = d - r;
            h11 += j1 * j1; h12 += j1 * j2; h13 += j1 * j3;
            h22 += j2 * j2; h23 += j2 * j3; h33 += j3 * j3;
            g1 += j1 * res; g2 += j2 * res; g3 += j3 * res;
        }
        const double m[3][3] = {{h11, h12, h13}, {h12, h22, h23}, {h13, h23, h33}};
        const double d0 = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                          m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                          m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if (!(std::fabs(d0) > gp::Resolution())) break;
        const double rhs[3] = {-g1, -g2, -g3};
        auto solveCol = [&](int col) {
            double c[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) c[i][j] = (j == col) ? rhs[i] : m[i][j];
            return (c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) -
                    c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) +
                    c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])) / d0;
        };
        const double da = solveCol(0), db = solveCol(1), dr = solveCol(2);
        if (!std::isfinite(da) || !std::isfinite(db) || !std::isfinite(dr)) break;
        a += da; b += db; r += dr;
        if (std::hypot(std::hypot(da, db), dr) <= tau) break;
    }
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(r)) return false;
    Rmaj = a; z0 = b; Rmin = r;
    return true;
}

TorusCert certifyTorus(const std::vector<ProfPt>& P, double tau) {
    TorusCert c;
    if (P.size() < 3) { c.why = "profile-underdetermined"; return c; }
    std::vector<double> zs, ps, rs;
    zs.reserve(P.size()); ps.reserve(P.size()); rs.reserve(P.size());
    for (const ProfPt& p : P) { zs.push_back(p.z); ps.push_back(p.psi); rs.push_back(p.rho); }
    c.levels = countLevels(zs, tau);
    c.cols = countColumns(ps, rs, tau);
    c.profPts = countProfPts(P, tau);
    if (c.levels < 3 || c.cols < 2 || c.profPts < 5) {
        c.why = "profile-underdetermined";
        return c;
    }
    c.determinate = true;
    c.lineRes = profileLineMaxResid(P);
    if (c.lineRes <= tau) { c.why = "cone-or-cylinder-or-plane"; return c; }
    if (!fitProfileCircle(P, tau, c.Rmaj, c.z0, c.Rmin)) { c.why = "fit-singular"; return c; }
    if (!(c.Rmaj > tau)) { c.why = "sphere"; return c; }
    if (!(c.Rmin > tau)) { c.why = "degenerate"; return c; }
    if (!(c.Rmaj > c.Rmin)) { c.why = "spindle-torus"; return c; }
    double sum2 = 0.0;
    for (const ProfPt& p : P) {
        const double d = std::fabs(std::hypot(p.rho - c.Rmaj, p.z - c.z0) - c.Rmin);
        c.maxDev = std::max(c.maxDev, d);
        sum2 += d * d;
    }
    c.rmsDev = std::sqrt(sum2 / (double)P.size());
    if (c.maxDev > tau) { c.why = "vertex-off-surface"; return c; }
    c.ok = true;
    c.why = "accept";
    return c;
}

void arcExtent(std::vector<double> a, double& lo, double& hi, double& maxGap) {
    lo = hi = maxGap = 0.0;
    if (a.empty()) return;
    std::sort(a.begin(), a.end());
    if (a.size() == 1) { lo = hi = a.front(); maxGap = kTwoPi; return; }
    std::size_t ig = a.size() - 1;
    maxGap = (a.front() + kTwoPi) - a.back();
    for (std::size_t i = 0; i + 1 < a.size(); ++i) {
        const double g = a[i + 1] - a[i];
        if (g > maxGap) { maxGap = g; ig = i; }
    }
    lo = a[(ig + 1) % a.size()];
    hi = a[ig];
    if (hi < lo) hi += kTwoPi;
}

enum class MemKind : uint8_t { Prov, CylGrow };

struct MemKey {
    MemKind kind = MemKind::Prov;
    int idx = -1;
    bool operator<(const MemKey& o) const {
        if (kind != o.kind) return kind < o.kind;
        return idx < o.idx;
    }
    bool operator==(const MemKey& o) const { return kind == o.kind && idx == o.idx; }
};

}  // namespace

bool claimToriT(const MeshView& mv, const SegmentParams& p, const DerivedTols& /*tol*/,
                SegmentWork& work) {
    try {
        const bool tDiag = torusDiagEnabled();
        int nCand = 0, nAdmit = 0, nRefProfile = 0, nRefLine = 0, nRefSphere = 0,
            nRefSpindle = 0, nRefVertex = 0;

        auto summary = [&]() {
            if (tDiag)
                std::fprintf(stderr,
                             "DIAG_TORUS_SUM candidates=%d admitted=%d refusedProfile=%d "
                             "refusedLine=%d refusedSphere=%d refusedSpindle=%d refusedVertex=%d\n",
                             nCand, nAdmit, nRefProfile, nRefLine, nRefSphere, nRefSpindle,
                             nRefVertex);
        };

        if (!mv.pts || !mv.tris || mv.nTri == 0) { summary(); return true; }
        if (!mv.triEdges) { summary(); return true; }
        const double q = mv.quantFloor;
        if (!(q > 0.0)) { summary(); return true; }
        if (work.accepted.empty()) { summary(); return true; }

        const int nTri = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const int nVtx = (int)mv.nVtx;
        const int nProv = (int)work.provisionals.size();
        const double tau = 2.0 * q;

        std::vector<int> g2l;
        for (int v = 0; v < nVtx; ++v) {
            const int g = mv.compVtx ? mv.compVtx[v] : v;
            if (g >= (int)g2l.size()) g2l.resize((std::size_t)g + 1, -1);
            g2l[g] = v;
        }

        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{-1, -1});
        for (int t = 0; t < nTri; ++t)
            for (int s = 0; s < 3; ++s) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || e >= nEdge) continue;
                if (edgeTris[(std::size_t)e][0] < 0) edgeTris[(std::size_t)e][0] = t;
                else if (edgeTris[(std::size_t)e][1] < 0) edgeTris[(std::size_t)e][1] = t;
            }

        std::vector<char> provUsed((std::size_t)nProv, 0);
        std::vector<char> cylReclaimed((std::size_t)work.accepted.size(), 0);

        bool committed = true;
        while (committed) {
            committed = false;

            std::vector<int> triProv((std::size_t)nTri, -1);
            for (int i = 0; i < nProv; ++i)
                for (int t : work.provisionals[(std::size_t)i].tris)
                    if (t >= 0 && t < nTri) triProv[(std::size_t)t] = i;

            std::vector<int> triCylGrow((std::size_t)nTri, -1);
            std::vector<int> axisIds, eligibleIds;
            for (int i = 0; i < (int)work.accepted.size(); ++i) {
                const Region& r = work.accepted[(std::size_t)i];
                if (cylReclaimed[(std::size_t)i]) continue;
                if (axisDonor(r, tau)) axisIds.push_back(i);
                if (eligibleCylGrow(r, tau)) {
                    eligibleIds.push_back(i);
                    for (int t : r.tris)
                        if (t >= 0 && t < nTri) triCylGrow[(std::size_t)t] = i;
                }
            }
            if (axisIds.empty()) break;

            std::sort(axisIds.begin(), axisIds.end(), [&](int a, int b) {
                const int ma = minTriOf(work.accepted[(std::size_t)a].tris);
                const int mb = minTriOf(work.accepted[(std::size_t)b].tris);
                if (ma != mb) return ma < mb;
                return a < b;
            });

            auto memOfTri = [&](int t) -> MemKey {
                MemKey m;
                m.idx = -1;
                if (t < 0 || t >= nTri) return m;
                const int pi = triProv[(std::size_t)t];
                if (pi >= 0 && !provUsed[(std::size_t)pi] &&
                    work.provisionals[(std::size_t)pi].claim == ProvClaim::Unclaimed) {
                    m.kind = MemKind::Prov;
                    m.idx = pi;
                    return m;
                }
                const int ci = triCylGrow[(std::size_t)t];
                if (ci >= 0 && !cylReclaimed[(std::size_t)ci]) {
                    m.kind = MemKind::CylGrow;
                    m.idx = ci;
                    return m;
                }
                return m;
            };

            std::vector<MemKey> allMem;
            for (int t = 0; t < nTri; ++t) {
                const MemKey mk = memOfTri(t);
                if (mk.idx < 0) continue;
                bool seen = false;
                for (const MemKey& x : allMem)
                    if (x == mk) { seen = true; break; }
                if (!seen) allMem.push_back(mk);
            }
            std::sort(allMem.begin(), allMem.end());
            std::vector<std::vector<MemKey>> memAdj(allMem.size());
            auto findMem = [&](const MemKey& k) -> int {
                for (int i = 0; i < (int)allMem.size(); ++i)
                    if (allMem[(std::size_t)i] == k) return i;
                return -1;
            };
            for (int e = 0; e < nEdge; ++e) {
                const int ta = edgeTris[(std::size_t)e][0], tb = edgeTris[(std::size_t)e][1];
                if (ta < 0 || tb < 0) continue;
                const MemKey ma = memOfTri(ta), mb = memOfTri(tb);
                if (ma.idx < 0 || mb.idx < 0 || ma == mb) continue;
                const int ia = findMem(ma), ib = findMem(mb);
                if (ia < 0 || ib < 0) continue;
                memAdj[(std::size_t)ia].push_back(mb);
                memAdj[(std::size_t)ib].push_back(ma);
            }
            for (auto& v : memAdj) {
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
            }

            auto trisOf = [&](const MemKey& mk) -> const std::vector<int>& {
                if (mk.kind == MemKind::Prov) return work.provisionals[(std::size_t)mk.idx].tris;
                return work.accepted[(std::size_t)mk.idx].tris;
            };

            for (int cylI : axisIds) {
                const Region& cyl = work.accepted[(std::size_t)cylI];
                gp_XYZ axis = cyl.ax.Direction().XYZ();
                const double am = axis.Modulus();
                if (!(am > gp::Resolution())) continue;
                axis /= am;
                const gp_XYZ axLoc = cyl.ax.Location().XYZ();
                gp_XYZ e1 = cyl.ax.XDirection().XYZ();
                e1 -= axis * e1.Dot(axis);
                if (!(e1.Modulus() > gp::Resolution())) continue;
                e1 /= e1.Modulus();
                const gp_XYZ e2 = axis.Crossed(e1);

                auto profOfLv = [&](int lv, ProfPt& out) {
                    const gp_XYZ d = localPt(mv, lv) - axLoc;
                    out.lv = lv;
                    out.z = d.Dot(axis);
                    const double x = d.Dot(e1), y = d.Dot(e2);
                    out.rho = std::hypot(x, y);
                    out.psi = std::atan2(y, x);
                    if (out.psi < 0.0) out.psi += kTwoPi;
                };
                // A z=const disk (the mouth-round's planar cap) has every
                // vertex on one v-iso, so certifiesTorus still passes; it is
                // a plane, not the blend. A3 must keep it.
                auto isAxialDisk = [&](const std::vector<int>& tris) -> bool {
                    double zLo = 1e300, zHi = -1e300, rLo = 1e300, rHi = -1e300;
                    int n = 0;
                    for (int t : tris) {
                        for (int k = 0; k < 3; ++k) {
                            const int lv = cornerLocal(g2l, mv, t, k);
                            if (lv < 0) continue;
                            ProfPt pp;
                            profOfLv(lv, pp);
                            zLo = std::min(zLo, pp.z);
                            zHi = std::max(zHi, pp.z);
                            rLo = std::min(rLo, pp.rho);
                            rHi = std::max(rHi, pp.rho);
                            ++n;
                        }
                    }
                    if (n < 3) return false;
                    return (zHi - zLo) <= tau;
                };

                std::vector<MemKey> seeds;
                for (int t : cyl.tris) {
                    if (t < 0 || t >= nTri) continue;
                    for (int s = 0; s < 3; ++s) {
                        const int ot = edgeOtherTri(edgeTris, mv.triEdges[t][s], t);
                        if (ot < 0 || ot >= nTri) continue;
                        const MemKey mk = memOfTri(ot);
                        if (mk.idx >= 0) seeds.push_back(mk);
                    }
                }
                std::sort(seeds.begin(), seeds.end());
                seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
                std::sort(seeds.begin(), seeds.end(), [&](const MemKey& a, const MemKey& b) {
                    const int ma = minTriOf(trisOf(a));
                    const int mb = minTriOf(trisOf(b));
                    if (ma != mb) return ma < mb;
                    return a < b;
                });

                for (const MemKey& seed : seeds) {
                    if (isAxialDisk(trisOf(seed))) continue;
                    nCand++;
                    std::vector<MemKey> members = {seed};
                    std::vector<char> inMem(allMem.size(), 0);
                    const int seedPos = findMem(seed);
                    if (seedPos < 0) continue;
                    inMem[(std::size_t)seedPos] = 1;

                    std::vector<char> vtxIn((std::size_t)nVtx, 0);
                    std::vector<ProfPt> curPts;
                    auto foldMem = [&](const MemKey& mk) {
                        for (int t : trisOf(mk))
                            for (int k = 0; k < 3; ++k) {
                                const int lv = cornerLocal(g2l, mv, t, k);
                                if (lv < 0 || vtxIn[(std::size_t)lv]) continue;
                                vtxIn[(std::size_t)lv] = 1;
                                ProfPt pp;
                                profOfLv(lv, pp);
                                curPts.push_back(pp);
                            }
                        std::sort(curPts.begin(), curPts.end(),
                                  [](const ProfPt& a, const ProfPt& b) { return a.lv < b.lv; });
                    };
                    foldMem(seed);

                    bool changed = true;
                    while (changed) {
                        changed = false;
                        std::vector<MemKey> frontier;
                        for (const MemKey& m : members) {
                            const int mi = findMem(m);
                            if (mi < 0) continue;
                            for (const MemKey& nb : memAdj[(std::size_t)mi]) {
                                if (inMem[(std::size_t)findMem(nb)]) continue;
                                frontier.push_back(nb);
                            }
                        }
                        std::sort(frontier.begin(), frontier.end(),
                                  [&](const MemKey& a, const MemKey& b) {
                                      const int ma = minTriOf(trisOf(a));
                                      const int mb = minTriOf(trisOf(b));
                                      if (ma != mb) return ma < mb;
                                      return a < b;
                                  });
                        frontier.erase(std::unique(frontier.begin(), frontier.end()),
                                       frontier.end());

                        for (const MemKey& u : frontier) {
                            const int ui = findMem(u);
                            if (ui < 0 || inMem[(std::size_t)ui]) continue;
                            if (isAxialDisk(trisOf(u))) continue;
                            std::vector<ProfPt> tent = curPts;
                            for (int t : trisOf(u))
                                for (int k = 0; k < 3; ++k) {
                                    const int lv = cornerLocal(g2l, mv, t, k);
                                    if (lv < 0) continue;
                                    bool have = false;
                                    for (const ProfPt& pp : tent)
                                        if (pp.lv == lv) { have = true; break; }
                                    if (have) continue;
                                    ProfPt pp;
                                    profOfLv(lv, pp);
                                    tent.push_back(pp);
                                }
                            std::sort(tent.begin(), tent.end(),
                                      [](const ProfPt& a, const ProfPt& b) { return a.lv < b.lv; });
                            TorusCert cert = certifyTorus(tent, tau);
                            const bool admit = !cert.determinate || cert.ok;
                            if (!admit) continue;
                            inMem[(std::size_t)ui] = 1;
                            members.push_back(u);
                            curPts = tent;
                            changed = true;
                        }
                    }

                    // Absorb every unclaimed provisional sharing an edge with the
                    // current band (mesh islands between the band and plane/cyl
                    // neighbours are not in memAdj but must be in the region).
                    {
                        std::vector<char> inTriBand((std::size_t)nTri, 0);
                        std::vector<char> provAbsorbed((std::size_t)nProv, 0);
                        for (const MemKey& m : members)
                            for (int t : trisOf(m)) inTriBand[(std::size_t)t] = 1;
                        bool touchGrew = true;
                        while (touchGrew) {
                            touchGrew = false;
                            std::vector<int> addProv;
                            for (int pi = 0; pi < nProv; ++pi) {
                                if (provUsed[(std::size_t)pi] || provAbsorbed[(std::size_t)pi]) continue;
                                if (work.provisionals[(std::size_t)pi].claim != ProvClaim::Unclaimed)
                                    continue;
                                bool touches = false;
                                for (int t : work.provisionals[(std::size_t)pi].tris) {
                                    if (t < 0 || t >= nTri) continue;
                                    for (int s = 0; s < 3 && !touches; ++s) {
                                        const int ot =
                                            edgeOtherTri(edgeTris, mv.triEdges[t][s], t);
                                        if (ot >= 0 && inTriBand[(std::size_t)ot]) touches = true;
                                    }
                                }
                                if (!touches) continue;
                                if (isAxialDisk(work.provisionals[(std::size_t)pi].tris)) continue;
                                std::vector<ProfPt> tent = curPts;
                                for (int t : work.provisionals[(std::size_t)pi].tris)
                                    for (int k = 0; k < 3; ++k) {
                                        const int lv = cornerLocal(g2l, mv, t, k);
                                        if (lv < 0) continue;
                                        bool have = false;
                                        for (const ProfPt& pp : tent)
                                            if (pp.lv == lv) { have = true; break; }
                                        if (have) continue;
                                        ProfPt pp;
                                        profOfLv(lv, pp);
                                        tent.push_back(pp);
                                    }
                                std::sort(tent.begin(), tent.end(),
                                          [](const ProfPt& a, const ProfPt& b) {
                                              return a.lv < b.lv;
                                          });
                                TorusCert cert = certifyTorus(tent, tau);
                                const bool admit = !cert.determinate || cert.ok;
                                if (!admit) continue;
                                addProv.push_back(pi);
                            }
                            for (int pi : addProv) {
                                provAbsorbed[(std::size_t)pi] = 1;
                                MemKey mk;
                                mk.kind = MemKind::Prov;
                                mk.idx = pi;
                                const int ui = findMem(mk);
                                if (ui >= 0) inMem[(std::size_t)ui] = 1;
                                members.push_back(mk);
                                for (int t : work.provisionals[(std::size_t)pi].tris)
                                    inTriBand[(std::size_t)t] = 1;
                                for (int t : work.provisionals[(std::size_t)pi].tris)
                                    for (int k = 0; k < 3; ++k) {
                                        const int lv = cornerLocal(g2l, mv, t, k);
                                        if (lv < 0) continue;
                                        bool have = false;
                                        for (const ProfPt& pp : curPts)
                                            if (pp.lv == lv) { have = true; break; }
                                        if (have) continue;
                                        ProfPt pp;
                                        profOfLv(lv, pp);
                                        curPts.push_back(pp);
                                    }
                                std::sort(curPts.begin(), curPts.end(),
                                          [](const ProfPt& a, const ProfPt& b) {
                                              return a.lv < b.lv;
                                          });
                                touchGrew = true;
                            }
                        }
                    }

                    std::sort(members.begin(), members.end());
                    std::vector<int> allTris;
                    for (const MemKey& m : members)
                        for (int t : trisOf(m)) allTris.push_back(t);
                    std::sort(allTris.begin(), allTris.end());
                    allTris.erase(std::unique(allTris.begin(), allTris.end()), allTris.end());
                    if (allTris.size() < 2) continue;

                    TorusCert final = certifyTorus(curPts, tau);
                    if (tDiag)
                        std::fprintf(stderr,
                                     "DIAG_TORUS rid=%d nTri=%zu levels=%d cols=%d profPts=%d "
                                     "Rmaj=%.6f Rmin=%.6f z0=%.6f maxDev=%.6e q=%.6e lineRes=%.6e "
                                     "admit=%d why=%s\n",
                                     cylI, allTris.size(), final.levels, final.cols, final.profPts,
                                     final.Rmaj, final.Rmin, final.z0, final.maxDev, q, final.lineRes,
                                     final.ok ? 1 : 0, final.why);

                    if (!final.ok) {
                        const char* w = final.why;
                        if (std::strcmp(w, "profile-underdetermined") == 0 ||
                            std::strcmp(w, "fit-singular") == 0)
                            ++nRefProfile;
                        else if (std::strcmp(w, "cone-or-cylinder-or-plane") == 0)
                            ++nRefLine;
                        else if (std::strcmp(w, "sphere") == 0 || std::strcmp(w, "degenerate") == 0)
                            ++nRefSphere;
                        else if (std::strcmp(w, "spindle-torus") == 0)
                            ++nRefSpindle;
                        else
                            ++nRefVertex;
                        continue;
                    }

                    bool reclaimOk = true;
                    std::vector<int> reclaimCyl;
                    for (const MemKey& m : members) {
                        if (m.kind != MemKind::CylGrow) continue;
                        const Region& cr = work.accepted[(std::size_t)m.idx];
                        const bool taken = final.maxDev < cr.maxVertexDev;
                        if (tDiag)
                            std::fprintf(stderr,
                                         "DIAG_TORUS_RECLAIM rid=%d cylR=%.4f nTri=%zu cylDev=%.6e "
                                         "cylDevQ=%.1f torusDevQ=%.3f taken=%d\n",
                                         m.idx, cr.radius, cr.tris.size(), cr.maxVertexDev,
                                         cr.maxVertexDev / q, final.maxDev / q, taken ? 1 : 0);
                        if (!taken) reclaimOk = false;
                        else reclaimCyl.push_back(m.idx);
                    }
                    if (!reclaimOk) continue;

                    std::vector<double> psis, vs;
                    for (const ProfPt& pp : curPts) {
                        psis.push_back(pp.psi);
                        double v = std::atan2(pp.z - final.z0, pp.rho - final.Rmaj);
                        if (v < 0.0) v += kTwoPi;
                        vs.push_back(v);
                    }
                    double uLo = 0, uHi = 0, uGap = 0, vLo = 0, vHi = 0, vGap = 0;
                    arcExtent(psis, uLo, uHi, uGap);
                    arcExtent(vs, vLo, vHi, vGap);
                    if (vHi - vLo >= kPi) continue;
                    while (vLo > kPi) vLo -= kTwoPi;
                    while (vLo <= -kPi) vLo += kTwoPi;
                    while (vHi > kPi) vHi -= kTwoPi;
                    while (vHi <= -kPi) vHi += kTwoPi;
                    if (vHi < vLo) vHi += kTwoPi;

                    int lvMin = std::numeric_limits<int>::max();
                    for (const ProfPt& pp : curPts) lvMin = std::min(lvMin, pp.lv);
                    gp_XYZ xdir(0, 0, 0);
                    {
                        const gp_XYZ d = localPt(mv, lvMin) - axLoc;
                        xdir = d - axis * d.Dot(axis);
                        if (!(xdir.Modulus() > gp::Resolution())) xdir = e1;
                        else xdir /= xdir.Modulus();
                    }

                    const gp_XYZ centre = axLoc + axis * final.z0;
                    const int nSides = std::max(final.cols, 1);
                    const double bandArc = kTwoPi / (double)nSides;
                    double hMin = 0.0;
                    for (int t : allTris) {
                        const gp_XYZ a = triVert(mv, t, 0);
                        const gp_XYZ b = triVert(mv, t, 1);
                        const gp_XYZ c = triVert(mv, t, 2);
                        const gp_XYZ ab = b - a, bc = c - b, ca = a - c;
                        const double la = bc.Modulus(), lb = ca.Modulus(), lc = ab.Modulus();
                        const double aa = 0.5 * ab.Crossed(ca).Modulus();
                        if (!(aa > gp::Resolution())) continue;
                        const double h = 2.0 * aa / std::max({la, lb, lc});
                        if (h > 0.0 && (hMin <= 0.0 || h < hMin)) hMin = h;
                    }
                    const double uSlop =
                        (hMin > gp::Resolution()) ? std::asin(std::min(1.0, tau / hMin)) : 0.0;
                    const bool closed360 =
                        final.cols >= 3 && uGap <= bandArc + uSlop;

                    double sigma = 0.0;
                    const double gamma = kTwoPi / (double)nSides;
                    const double kSeg = (gamma - std::sin(gamma)) / (4.0 * std::sin(0.5 * gamma));
                    double dVolSum = 0.0;
                    for (int t : allTris) {
                        gp_XYZ n;
                        const double a2 = triAreaNormal(mv, t, n);
                        const gp_XYZ d = triCentroid(mv, t) - axLoc;
                        const double zc = d.Dot(axis);
                        gp_XYZ rad = d - axis * zc;
                        const double rc = rad.Modulus();
                        if (!(rc > gp::Resolution())) continue;
                        rad /= rc;
                        gp_XYZ mhat = rad * (rc - final.Rmaj) + axis * (zc - final.z0);
                        const double mm2 = mhat.Modulus();
                        if (mm2 > gp::Resolution()) sigma += a2 * n.Dot(mhat / mm2);
                        dVolSum += a2 * rc * kSeg * std::fabs(n.Dot(rad));
                    }

                    Region r;
                    r.type = SurfType::Torus;
                    r.origin = Origin::TorusBlend;
                    r.ax = gp_Ax3(gp_Pnt(centre), gp_Dir(axis), gp_Dir(xdir));
                    r.radius = final.Rmaj;
                    r.radius2 = final.Rmin;
                    r.uMin = closed360 ? 0.0 : uLo;
                    r.uMax = closed360 ? kTwoPi : uHi;
                    r.vMin = vLo;
                    r.vMax = vHi;
                    r.closed360 = closed360;
                    r.outwardNormal = sigma >= 0.0;
                    r.tris = allTris;
                    r.maxVertexDev = final.maxDev;
                    r.rmsVertexDev = final.rmsDev;
                    r.chordSagitta = chordSagitta(final.Rmaj, nSides);
                    r.nSides = nSides;
                    r.dVolPredicted = (sigma >= 0.0 ? 1.0 : -1.0) * dVolSum;
                    r.reject = Reject::None;
                    r.builtAs = BuiltAs::NotBuilt;
                    r.filletNbrA = cylI;
                    r.filletNbrB = -1;

                    for (const MemKey& m : members) {
                        if (m.kind == MemKind::Prov)
                            work.provisionals[(std::size_t)m.idx].claim = ProvClaim::ConsumedTorus;
                    }

                    TorusPendingReclaim pend;
                    for (const MemKey& m : members) {
                        if (m.kind == MemKind::Prov) pend.provIdx.push_back(m.idx);
                    }
                    for (int ci : reclaimCyl) {
                        pend.stashedCylinders.push_back(work.accepted[(std::size_t)ci]);
                        pend.reclaimAccIdx.push_back(ci);
                    }
                    const int packIdx = (int)work.torusPending.size();
                    work.torusPending.push_back(std::move(pend));
                    r.torusPending = packIdx;

                    for (int ci : reclaimCyl) cylReclaimed[(std::size_t)ci] = 1;

                    std::vector<Region> kept;
                    for (int i = 0; i < (int)work.accepted.size(); ++i) {
                        if (cylReclaimed[(std::size_t)i] &&
                            eligibleCylGrow(work.accepted[(std::size_t)i], tau))
                            continue;
                        kept.push_back(work.accepted[(std::size_t)i]);
                    }
                    kept.push_back(std::move(r));
                    work.accepted.swap(kept);

                    std::vector<char> newReclaimed(work.accepted.size(), 0);
                    cylReclaimed.swap(newReclaimed);
                    work.torusAdmitted = true;

                    for (const MemKey& m : members)
                        if (m.kind == MemKind::Prov) provUsed[(std::size_t)m.idx] = 1;

                    ++nAdmit;
                    committed = true;
                    break;
                }
                if (committed) break;
            }
        }

        if (tDiag) {
            for (int i = 0; i < (int)work.accepted.size(); ++i) {
                const Region& cr = work.accepted[(std::size_t)i];
                if (!eligibleCylGrow(cr, tau)) continue;
                std::fprintf(stderr,
                             "DIAG_TORUS_RECLAIM rid=%d cylR=%.4f nTri=%zu cylDev=%.6e "
                             "cylDevQ=%.1f torusDevQ=0.000 taken=0\n",
                             i, cr.radius, cr.tris.size(), cr.maxVertexDev, cr.maxVertexDev / q);
            }
        }

        summary();
        return true;
    } catch (...) {
        return true;
    }
}

}  // namespace refit
}  // namespace stl2step
