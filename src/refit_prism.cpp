// RULE 5.1 — prismaticity predicate + cap-level extraction.
// Tolerances are SELF-COMPUTED (RULE 4.2a). Never throws.
// SPDX-License-Identifier: MIT

#include "refit_prism.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <utility>
#include <vector>

#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {

void derivePrismTols(const MeshView& mv, const RegionSet& rs, PrismTols& t) {
    const double weld = mv.weldTol > 0.0 ? mv.weldTol : 0.0;
    const double diag = std::isfinite(mv.diag) && mv.diag > 0.0 ? mv.diag : 0.0;
    const double tauSurf = std::max(5e-5, std::max(4.0 * weld, 1e-6 * diag));
    t.tauSurf = tauSurf;
    t.tauLvl = tauSurf;
    t.tauFit = tauSurf;
    double hMin = 0.0;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder) continue;
        const double h = std::abs(r.vMax - r.vMin);
        if (!(h > 0.0) || !std::isfinite(h)) continue;
        if (hMin <= 0.0 || h < hMin) hMin = h;
    }
    // Floor of the formula when the component has no usable cylinder height.
    t.tauAx = 1e-6;
    if (hMin > 0.0)
        t.tauAx = std::max(1e-6, 2.0 * tauSurf / hMin);
}

namespace {

bool prismDiagOn() {
    const char* e = std::getenv("STL2STEP_PRISM_DIAG");
    return e && e[0] && e[0] != '0';
}

int compIndex(const RegionSet& rs) {
    return rs.compRoot >= 0 ? rs.compRoot : 0;
}

void emitDiag(const RegionSet& rs, const PrismLevels& lv, const PrismTols& t,
              int nCyl, int nPlane, int nCap, int nLat, int nOblique) {
    if (!prismDiagOn()) return;
    const int c = compIndex(rs);
    std::fprintf(stderr,
                 "DIAG_PRISM comp=%d ok=%d failedCond=%d nCyl=%d nPlane=%d "
                 "nCap=%d nLat=%d nOblique=%d tauAx=%.3e tauLvl=%.3e\n",
                 c, lv.ok ? 1 : 0, lv.failedCond, nCyl, nPlane, nCap, nLat, nOblique,
                 t.tauAx, t.tauLvl);
    for (size_t i = 0; i < lv.y.size(); ++i) {
        const int rid = i < lv.capRegion.size() ? lv.capRegion[i] : -1;
        std::fprintf(stderr, "DIAG_PRISMLVL comp=%d i=%d y=%.6f capRegion=%d\n",
                     c, static_cast<int>(i), lv.y[i], rid);
    }
}

template <class Fn>
void parallelFor(size_t n, Fn fn) {
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    if (n < 2 || hw < 2) {
        for (size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    const unsigned w = static_cast<unsigned>(std::min<size_t>(hw, n));
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(w);
    for (unsigned k = 0; k < w; ++k) {
        pool.emplace_back([&]() {
            for (;;) {
                const size_t i = next.fetch_add(1);
                if (i >= n) break;
                fn(i);
            }
        });
    }
    for (auto& th : pool) th.join();
}

double triArea(const MeshView& mv, int localTri) {
    if (!mv.pts || !mv.tris || !mv.compTris) return 0.0;
    if (localTri < 0 || static_cast<size_t>(localTri) >= mv.nTri) return 0.0;
    const int g = mv.compTris[localTri];
    if (g < 0) return 0.0;
    const int ia = mv.tris[g][0];
    const int ib = mv.tris[g][1];
    const int ic = mv.tris[g][2];
    const gp_XYZ a = mv.pts[ia];
    const gp_XYZ ab = mv.pts[ib] - a;
    const gp_XYZ ac = mv.pts[ic] - a;
    return 0.5 * ab.Crossed(ac).Modulus();
}

double regionArea(const MeshView& mv, const Region& r) {
    const size_t n = r.tris.size();
    if (n == 0) return 0.0;
    if (n < 8) {
        double s = 0.0;
        for (int t : r.tris) s += triArea(mv, t);
        return s;
    }
    std::vector<double> part(n, 0.0);
    parallelFor(n, [&](size_t i) { part[i] = triArea(mv, r.tris[i]); });
    double s = 0.0;
    for (double v : part) s += v;
    return s;
}

double chainLength(const MeshView& mv, const BoundaryChain& ch) {
    if (!mv.pts || !mv.compVtx || ch.meshVerts.size() < 2) return 0.0;
    double len = 0.0;
    auto pt = [&](int loc) -> gp_XYZ {
        if (loc < 0 || static_cast<size_t>(loc) >= mv.nVtx) return gp_XYZ(0, 0, 0);
        return mv.pts[mv.compVtx[loc]];
    };
    for (size_t i = 1; i < ch.meshVerts.size(); ++i)
        len += (pt(ch.meshVerts[i]) - pt(ch.meshVerts[i - 1])).Modulus();
    if (ch.closedLoop && ch.meshVerts.size() > 2)
        len += (pt(ch.meshVerts.front()) - pt(ch.meshVerts.back())).Modulus();
    return len;
}

double capPerimeter(const MeshView& mv, const RegionSet& rs, const Region& r) {
    double peri = 0.0;
    for (const Loop& lp : r.loops) {
        for (int ci : lp.chainIdx) {
            if (ci < 0 || static_cast<size_t>(ci) >= rs.chains.size()) continue;
            peri += chainLength(mv, rs.chains[static_cast<size_t>(ci)]);
        }
    }
    if (peri > 0.0) return peri;
    const double a = regionArea(mv, r);
    return a > 0.0 ? 4.0 * std::sqrt(a) : 0.0;
}

gp_XYZ unitOrZero(const gp_Dir& d) {
    gp_XYZ x = d.XYZ();
    const double m = x.Modulus();
    if (!(m > 0.0) || !std::isfinite(m)) return gp_XYZ(0, 0, 0);
    x /= m;
    return x;
}

}  // namespace

PrismLevels detectPrismatic(const MeshView& mv, const RegionSet& rs, const PrismTols& tIn) {
    PrismLevels out;
    out.ok = false;
    out.failedCond = 1;
    PrismTols t = tIn;
    int nCyl = 0, nPlane = 0, nCap = 0, nLat = 0, nOblique = 0;
    try {
        std::vector<const Region*> cyls;
        std::vector<const Region*> planes;
        cyls.reserve(rs.regions.size());
        planes.reserve(rs.regions.size());
        for (const Region& r : rs.regions) {
            if (r.type == SurfType::Cylinder) cyls.push_back(&r);
            else if (r.type == SurfType::Plane) planes.push_back(&r);
        }
        nCyl = static_cast<int>(cyls.size());
        nPlane = static_cast<int>(planes.size());

        if (t.tauSurf <= 0.0 || t.tauAx <= 0.0)
            derivePrismTols(mv, rs, t);

        // --- 1: at least two recognized cylinders --------------------------------
        if (nCyl < 2) {
            out.failedCond = 1;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }

        // --- 2: common axis ------------------------------------------------------
        std::vector<gp_XYZ> axes(static_cast<size_t>(nCyl));
        for (int i = 0; i < nCyl; ++i)
            axes[static_cast<size_t>(i)] = unitOrZero(cyls[static_cast<size_t>(i)]->ax.Direction());

        const size_t pairs = static_cast<size_t>(nCyl) * static_cast<size_t>(nCyl);
        std::vector<double> pairSin(pairs, 0.0);
        parallelFor(pairs, [&](size_t k) {
            const int i = static_cast<int>(k / static_cast<size_t>(nCyl));
            const int j = static_cast<int>(k % static_cast<size_t>(nCyl));
            if (j <= i) return;
            pairSin[k] = axes[static_cast<size_t>(i)]
                             .Crossed(axes[static_cast<size_t>(j)])
                             .Modulus();
        });
        double maxSin = 0.0;
        for (double s : pairSin)
            if (s > maxSin) maxSin = s;
        if (!(maxSin < t.tauAx)) {
            out.failedCond = 2;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }

        gp_XYZ ref = axes[0];
        if (ref.Modulus() <= 0.0) {
            out.failedCond = 2;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }
        gp_XYZ sum(0, 0, 0);
        for (const gp_XYZ& a : axes) {
            gp_XYZ u = a;
            if (u.Dot(ref) < 0.0) u.Reverse();
            sum += u;
        }
        const double sm = sum.Modulus();
        if (!(sm > 0.0)) {
            out.failedCond = 2;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }
        sum /= sm;
        out.axis = gp_Dir(sum);
        const gp_XYZ ahat = sum;

        // --- 3: no oblique planes ------------------------------------------------
        std::vector<int> kind(static_cast<size_t>(nPlane), 2);
        parallelFor(static_cast<size_t>(nPlane), [&](size_t i) {
            const gp_XYZ n = unitOrZero(planes[i]->ax.Direction());
            const double nd = std::abs(n.Dot(ahat));
            if (nd > 1.0 - t.tauAx) kind[i] = 0;
            else if (nd < t.tauAx) kind[i] = 1;
            else kind[i] = 2;
        });
        std::vector<const Region*> caps;
        for (size_t i = 0; i < kind.size(); ++i) {
            if (kind[i] == 0) {
                ++nCap;
                caps.push_back(planes[i]);
            } else if (kind[i] == 1) {
                ++nLat;
            } else {
                ++nOblique;
            }
        }
        if (nOblique > 0) {
            out.failedCond = 3;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }

        // --- 4: >=2 distinct cap levels (cluster n.p along â) --------------------
        struct CapOff {
            double y = 0.0;
            int id = -1;
            const Region* r = nullptr;
        };
        std::vector<CapOff> offs;
        offs.reserve(caps.size());
        for (const Region* r : caps) {
            CapOff o;
            o.y = ahat.Dot(r->ax.Location().XYZ());
            o.id = r->id;
            o.r = r;
            offs.push_back(o);
        }
        std::sort(offs.begin(), offs.end(), [](const CapOff& a, const CapOff& b) {
            if (a.y < b.y) return true;
            if (a.y > b.y) return false;
            return a.id < b.id;
        });
        std::vector<std::vector<CapOff>> clusters;
        for (const CapOff& o : offs) {
            if (clusters.empty() || (o.y - clusters.back().back().y) > t.tauLvl)
                clusters.emplace_back();
            clusters.back().push_back(o);
        }
        out.y.clear();
        out.capRegion.clear();
        for (const auto& cl : clusters) {
            double mean = 0.0;
            int bestId = cl[0].id;
            double bestA = -1.0;
            for (const CapOff& o : cl) {
                mean += o.y;
                const double ar = o.r ? regionArea(mv, *o.r) : 0.0;
                if (ar > bestA || (ar == bestA && o.id < bestId)) {
                    bestA = ar;
                    bestId = o.id;
                }
            }
            mean /= static_cast<double>(cl.size());
            out.y.push_back(mean);
            out.capRegion.push_back(bestId);
        }
        if (static_cast<int>(out.y.size()) < 2) {
            out.failedCond = 4;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }

        // --- 5: every cylinder spans a contiguous run of levels ------------------
        const size_t L = out.y.size();
        auto nearest = [&](double yv) -> int {
            int best = -1;
            double bestD = t.tauLvl;
            for (size_t k = 0; k < L; ++k) {
                const double d = std::abs(yv - out.y[k]);
                if (d <= t.tauLvl && (best < 0 || d < bestD)) {
                    best = static_cast<int>(k);
                    bestD = d;
                }
            }
            return best;
        };
        std::vector<int> spanOk(static_cast<size_t>(nCyl), 0);
        parallelFor(static_cast<size_t>(nCyl), [&](size_t i) {
            const Region& r = *cyls[i];
            const gp_XYZ c = r.ax.Location().XYZ();
            const gp_XYZ ad = unitOrZero(r.ax.Direction());
            const double y0 = ahat.Dot(c + ad * r.vMin);
            const double y1 = ahat.Dot(c + ad * r.vMax);
            int ia = nearest(y0);
            int ib = nearest(y1);
            if (ia < 0 || ib < 0) return;
            if (ia > ib) std::swap(ia, ib);
            if (ib > ia) spanOk[i] = 1;
        });
        for (int v : spanOk) {
            if (!v) {
                out.failedCond = 5;
                emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
                return out;
            }
        }

        // --- 6: signed cap-area closure (per-level flux sums to ~0) --------------
        std::vector<double> signedA(caps.size(), 0.0);
        std::vector<double> peri(caps.size(), 0.0);
        parallelFor(caps.size(), [&](size_t i) {
            const Region& r = *caps[i];
            const gp_XYZ n = unitOrZero(r.ax.Direction());
            signedA[i] = regionArea(mv, r) * n.Dot(ahat);
            peri[i] = capPerimeter(mv, rs, r);
        });
        double flux = 0.0;
        double periScale = 0.0;
        for (size_t i = 0; i < caps.size(); ++i) {
            flux += signedA[i];
            periScale += peri[i];
        }
        if (periScale <= 0.0) periScale = 1.0;
        if (!(std::abs(flux) < t.tauFit * periScale)) {
            out.failedCond = 6;
            emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
            return out;
        }

        out.ok = true;
        out.failedCond = 0;
        emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
        return out;
    } catch (...) {
        out.ok = false;
        if (out.failedCond <= 0) out.failedCond = 1;
        emitDiag(rs, out, t, nCyl, nPlane, nCap, nLat, nOblique);
        return out;
    }
}

}  // namespace refit
}  // namespace stl2step
