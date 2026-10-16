// Relational reclaim (D-fc-1.11, D-fc-4). Line mode belongs to lane b1-reclaim.
// Torus mode (D-fc-4.6): atoms join the tori stages T and H1 certified (F2),
// then every torus takes its closure from its own rails (F4). Never throws.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
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

bool torusDiagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

gp_XYZ cornerOf(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.pts[mv.tris[gt][k]];
}

int cornerId(const MeshView& mv, int t, int k) {
    const int gt = mv.compTris ? mv.compTris[t] : t;
    return mv.tris[gt][k];
}

// Stage T's extent rule (refit_torus.cpp): the complement of the largest gap.
void arcExtentOf(std::vector<double> a, double& lo, double& hi, double& maxGap) {
    lo = hi = maxGap = 0.0;
    if (a.empty()) return;
    std::sort(a.begin(), a.end());
    if (a.size() == 1) {
        lo = hi = a.front();
        maxGap = kTwoPi;
        return;
    }
    std::size_t ig = a.size() - 1;
    maxGap = (a.front() + kTwoPi) - a.back();
    for (std::size_t i = 0; i + 1 < a.size(); ++i) {
        const double g = a[i + 1] - a[i];
        if (g > maxGap) {
            maxGap = g;
            ig = i;
        }
    }
    lo = a[(ig + 1) % a.size()];
    hi = a[ig];
    if (hi < lo) hi += kTwoPi;
}

// The torus keeps its certified surface; the union's own vertices give its
// deviation, extents and predicted volume (stage T's formulas).
void refreshTorusRegion(const MeshView& mv, Region& r) {
    std::sort(r.tris.begin(), r.tris.end());
    r.tris.erase(std::unique(r.tris.begin(), r.tris.end()), r.tris.end());
    const gp_XYZ o = r.ax.Location().XYZ();
    const gp_XYZ a = r.ax.Direction().XYZ();
    const gp_XYZ ex = r.ax.XDirection().XYZ();
    const gp_XYZ ey = a.Crossed(ex);
    std::vector<int> ids;
    ids.reserve(r.tris.size() * 3);
    for (int t : r.tris)
        for (int k = 0; k < 3; ++k) ids.push_back(cornerId(mv, t, k));
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    double maxD = 0.0, sum2 = 0.0;
    std::vector<double> us, vs;
    us.reserve(ids.size());
    vs.reserve(ids.size());
    for (int g : ids) {
        const gp_XYZ d = mv.pts[g] - o;
        const double z = d.Dot(a);
        const gp_XYZ rad = d - a * z;
        const double rho = rad.Modulus();
        const double dev = std::fabs(std::hypot(rho - r.radius, z) - r.radius2);
        maxD = std::max(maxD, dev);
        sum2 += dev * dev;
        double u = std::atan2(rad.Dot(ey), rad.Dot(ex));
        if (u < 0.0) u += kTwoPi;
        us.push_back(u);
        double v = std::atan2(z, rho - r.radius);
        if (v < 0.0) v += kTwoPi;
        vs.push_back(v);
    }
    r.maxVertexDev = maxD;
    r.rmsVertexDev = ids.empty() ? 0.0 : std::sqrt(sum2 / (double)ids.size());
    double uLo = 0, uHi = 0, uGap = 0, vLo = 0, vHi = 0, vGap = 0;
    arcExtentOf(us, uLo, uHi, uGap);
    arcExtentOf(vs, vLo, vHi, vGap);
    while (vLo > kPi) vLo -= kTwoPi;
    while (vLo <= -kPi) vLo += kTwoPi;
    while (vHi > kPi) vHi -= kTwoPi;
    while (vHi <= -kPi) vHi += kTwoPi;
    if (vHi < vLo) vHi += kTwoPi;
    r.vMin = vLo;
    r.vMax = vHi;
    if (!r.closed360) {
        r.uMin = uLo;
        r.uMax = uHi;
    }
    const int nSides = std::max(r.nSides, 1);
    const double gamma = kTwoPi / (double)nSides;
    const double sHalf = std::sin(0.5 * gamma);
    const double kSeg = (sHalf > gp::Resolution()) ? (gamma - std::sin(gamma)) / (4.0 * sHalf) : 0.0;
    double dVolSum = 0.0;
    for (int t : r.tris) {
        const gp_XYZ p0 = cornerOf(mv, t, 0), p1 = cornerOf(mv, t, 1), p2 = cornerOf(mv, t, 2);
        gp_XYZ n = (p1 - p0).Crossed(p2 - p0);
        const double mag = n.Modulus();
        if (!(mag > gp::Resolution())) continue;
        n /= mag;
        const double a2 = 0.5 * mag;
        const gp_XYZ d = (p0 + p1 + p2) / 3.0 - o;
        const double zc = d.Dot(a);
        gp_XYZ rad = d - a * zc;
        const double rc = rad.Modulus();
        if (!(rc > gp::Resolution())) continue;
        rad /= rc;
        dVolSum += a2 * rc * kSeg * std::fabs(n.Dot(rad));
    }
    r.dVolPredicted = (r.outwardNormal ? 1.0 : -1.0) * dVolSum;
}

bool reclaimTorusMode(const MeshView& mv, SegmentWork& work) {
    const double tau = tauOf(mv);
    if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;
    bool anyTorus = false;
    for (const Region& r : work.accepted)
        if (r.origin == Origin::TorusBlend) anyTorus = true;
    if (!anyTorus) return true;
    const bool diag = torusDiagOn();
    const int nTri = (int)mv.nTri;
    const int nEdge = (int)mv.nEdge;
    std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                             std::array<int, 2>{-1, -1});
    for (int t = 0; t < nTri; ++t)
        for (int k = 0; k < 3; ++k) {
            const int e = mv.triEdges[t][k];
            if (e < 0 || e >= nEdge) continue;
            if (edgeTris[(std::size_t)e][0] < 0) edgeTris[(std::size_t)e][0] = t;
            else if (edgeTris[(std::size_t)e][1] < 0) edgeTris[(std::size_t)e][1] = t;
        }
    const int nAcc = (int)work.accepted.size();
    std::vector<int> owner((std::size_t)nTri, -1);
    for (int i = 0; i < nAcc; ++i)
        for (int t : work.accepted[(std::size_t)i].tris)
            if (t >= 0 && t < nTri) owner[(std::size_t)t] = i;
    std::vector<char> dead((std::size_t)nAcc, 0);
    int nDead = 0;

    // F2: an atom edge-adjacent to a certified torus joins it whole when every
    // vertex lies within tau. The torus band was certified without it, so the
    // leave-one-region-out count holds; the union only grows to a fixed point.
    for (int s = 0; s < nAcc; ++s) {
        if (work.accepted[(std::size_t)s].origin != Origin::TorusBlend) continue;
        int nAbs = 0, nAbsTri = 0, nRefused = 0;
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<int> cand;
            for (int t : work.accepted[(std::size_t)s].tris) {
                if (t < 0 || t >= nTri) continue;
                for (int k = 0; k < 3; ++k) {
                    const int e = mv.triEdges[t][k];
                    if (e < 0 || e >= nEdge) continue;
                    const int o = edgeTris[(std::size_t)e][0] == t ? edgeTris[(std::size_t)e][1]
                                                                    : edgeTris[(std::size_t)e][0];
                    if (o < 0) continue;
                    const int ri = owner[(std::size_t)o];
                    if (ri >= 0 && ri != s && !dead[(std::size_t)ri]) cand.push_back(ri);
                }
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
            for (int ri : cand) {
                const Region& R = work.accepted[(std::size_t)ri];
                if (!torusAtom(R) || ri == work.accepted[(std::size_t)s].filletNbrA) continue;
                if (!(regionTorusDev(mv, work.accepted[(std::size_t)s], R) <= tau)) {
                    ++nRefused;
                    continue;
                }
                Region& S = work.accepted[(std::size_t)s];
                S.tris.insert(S.tris.end(), R.tris.begin(), R.tris.end());
                for (int t : R.tris)
                    if (t >= 0 && t < nTri) owner[(std::size_t)t] = s;
                if (S.torusPending >= 0 && (std::size_t)S.torusPending < work.torusPending.size()) {
                    TorusPendingReclaim& pack = work.torusPending[(std::size_t)S.torusPending];
                    pack.stashedCylinders.push_back(R);
                    pack.reclaimAccIdx.push_back(ri);
                }
                nAbsTri += (int)R.tris.size();
                ++nAbs;
                dead[(std::size_t)ri] = 1;
                ++nDead;
                changed = true;
            }
        }
        if (nAbs > 0) refreshTorusRegion(mv, work.accepted[(std::size_t)s]);
        if (diag && (nAbs > 0 || nRefused > 0)) {
            const Region& S = work.accepted[(std::size_t)s];
            std::fprintf(stderr,
                         "DIAG_REL_TORUS_ABSORB acc=%d Rmaj=%.6f Rmin=%.6f absorbed=%d absTri=%d "
                         "refusedAdj=%d nTri=%zu maxDevTau=%.4f\n",
                         s, S.radius, S.radius2, nAbs, nAbsTri, nRefused, S.tris.size(),
                         S.maxVertexDev / tau);
        }
    }

    // F4: closure by the rails, the measure the builder calls fullRing.
    for (int s = 0; s < nAcc; ++s) {
        Region& r = work.accepted[(std::size_t)s];
        if (dead[(std::size_t)s] || r.origin != Origin::TorusBlend || r.closed360) continue;
        const TorusRails rails = torusRails(mv, r);
        if (diag)
            std::fprintf(stderr,
                         "DIAG_REL_TORUS_RAILS acc=%d Rmaj=%.6f Rmin=%.6f nTri=%zu gap=%.6f "
                         "maxStep=%.6f full=%d\n",
                         s, r.radius, r.radius2, r.tris.size(), rails.gap, rails.maxStep,
                         rails.full ? 1 : 0);
        if (!rails.full) continue;
        r.closed360 = true;
        r.uMin = 0.0;
        r.uMax = kTwoPi;
    }

    if (nDead > 0) {
        std::vector<int> remap((std::size_t)nAcc, -1);
        std::vector<Region> kept;
        kept.reserve((std::size_t)(nAcc - nDead));
        for (int i = 0; i < nAcc; ++i) {
            if (dead[(std::size_t)i]) continue;
            remap[(std::size_t)i] = (int)kept.size();
            kept.push_back(std::move(work.accepted[(std::size_t)i]));
        }
        for (Region& r : kept)
            if (r.origin == Origin::TorusBlend && r.filletNbrA >= 0 && r.filletNbrA < nAcc)
                r.filletNbrA = remap[(std::size_t)r.filletNbrA];
        work.accepted.swap(kept);
    }
    return true;
}

}  // namespace

bool torusAtom(const Region& r) {
    return r.origin == Origin::CylGrow && r.type == SurfType::Cylinder && !r.closed360 &&
           !r.tris.empty();
}

double torusSurfaceDev(const Region& torus, const gp_XYZ& p) {
    const gp_XYZ o = torus.ax.Location().XYZ();
    const gp_XYZ a = torus.ax.Direction().XYZ();
    const gp_XYZ d = p - o;
    const double z = d.Dot(a);
    const double rho = (d - a * z).Modulus();
    return std::fabs(std::hypot(rho - torus.radius, z) - torus.radius2);
}

double regionTorusDev(const MeshView& mv, const Region& torus, const Region& region) {
    double mx = 0.0;
    for (int t : region.tris) {
        if (t < 0 || (std::size_t)t >= mv.nTri) return std::numeric_limits<double>::infinity();
        for (int k = 0; k < 3; ++k) mx = std::max(mx, torusSurfaceDev(torus, cornerOf(mv, t, k)));
    }
    return mx;
}

bool leaveOneOutDeterminate(const std::vector<std::array<ProfPt, 3>>& triPts,
                            const std::vector<int>& triAtom, int nAtoms, double tau) {
    if (nAtoms <= 0) return true;
    if (triAtom.size() != triPts.size()) return false;
    std::vector<ProfPt> pts;
    pts.reserve(triPts.size() * 3);
    for (int k = 0; k < nAtoms; ++k) {
        pts.clear();
        for (std::size_t i = 0; i < triPts.size(); ++i) {
            if (triAtom[i] == k) continue;
            for (const ProfPt& p : triPts[i]) pts.push_back(p);
        }
        std::sort(pts.begin(), pts.end(), [](const ProfPt& a, const ProfPt& b) { return a.lv < b.lv; });
        pts.erase(std::unique(pts.begin(), pts.end(),
                              [](const ProfPt& a, const ProfPt& b) { return a.lv == b.lv; }),
                  pts.end());
        if (!certifyTorus(pts, tau).determinate) return false;
    }
    return true;
}

TorusRails torusRails(const MeshView& mv, const Region& torus) {
    TorusRails out;
    out.gap = kTwoPi;
    const int nTri = (int)mv.nTri;
    const int nEdge = (int)mv.nEdge;
    if (!mv.pts || !mv.tris || !mv.triEdges || !mv.compEdges || nTri == 0) return out;
    const gp_XYZ o = torus.ax.Location().XYZ();
    const gp_XYZ a = torus.ax.Direction().XYZ();
    const gp_XYZ ex = torus.ax.XDirection().XYZ();
    const gp_XYZ ey = a.Crossed(ex);
    auto az = [&](const gp_XYZ& p) {
        const gp_XYZ d = p - o;
        return std::atan2(d.Dot(ey), d.Dot(ex));
    };
    std::vector<char> in((std::size_t)nTri, 0);
    for (int t : torus.tris)
        if (t >= 0 && t < nTri) in[(std::size_t)t] = 1;
    std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                             std::array<int, 2>{-1, -1});
    for (int t : torus.tris) {
        if (t < 0 || t >= nTri) continue;
        for (int k = 0; k < 3; ++k) {
            const int e = mv.triEdges[t][k];
            if (e < 0 || e >= nEdge) continue;
            if (edgeTris[(std::size_t)e][0] < 0) edgeTris[(std::size_t)e][0] = t;
            else if (edgeTris[(std::size_t)e][1] < 0) edgeTris[(std::size_t)e][1] = t;
        }
    }
    std::vector<int> ids;
    ids.reserve(torus.tris.size() * 3);
    double maxStep = 0.0;
    for (int t : torus.tris) {
        if (t < 0 || t >= nTri) continue;
        for (int k = 0; k < 3; ++k) {
            ids.push_back(cornerId(mv, t, k));
            const int e = mv.triEdges[t][k];
            if (e < 0 || e >= nEdge) continue;
            // A boundary edge has one triangle of this region: the mesh's other
            // triangle across it, if any, belongs to a neighbour.
            if (edgeTris[(std::size_t)e][1] >= 0) continue;
            const std::pair<int, int>& ce = mv.compEdges[e];
            const int ga = mv.compVtx ? mv.compVtx[ce.first] : ce.first;
            const int gb = mv.compVtx ? mv.compVtx[ce.second] : ce.second;
            double dd = std::fabs(az(mv.pts[ga]) - az(mv.pts[gb]));
            if (dd > kPi) dd = kTwoPi - dd;
            maxStep = std::max(maxStep, dd);
        }
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    std::vector<double> ps;
    ps.reserve(ids.size());
    for (int g : ids) ps.push_back(az(mv.pts[g]));
    std::sort(ps.begin(), ps.end());
    double gap = ps.empty() ? kTwoPi : (kTwoPi - ps.back() + ps.front());
    for (std::size_t i = 1; i < ps.size(); ++i) gap = std::max(gap, ps[i] - ps[i - 1]);
    out.gap = gap;
    out.maxStep = maxStep;
    out.full = maxStep > 0.0 && gap <= maxStep;
    return out;
}

// Line mode: lane b1-reclaim (D-fc-1.11); no-op here.
bool reclaimRelational(const MeshView& mv, const SegmentParams&, const DerivedTols&,
                       SegmentWork& work, RelRevolveMode mode) {
    try {
        if (mode == RelRevolveMode::Torus) return reclaimTorusMode(mv, work);
        return true;
    } catch (...) {
        return true;
    }
}

}  // namespace refit
}  // namespace stl2step
