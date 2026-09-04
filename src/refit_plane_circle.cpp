// Detector B — CIRCLE on already-analytic plane loops (1.3.0).
//
// Builder (agent 06) calls tryPlaneLoopCircles during collapse, then
// loopIsCircle / planeLoopCircleOf / planeLoopCircleForChain to emit
// Geom_Circle. This TU does not rewrite STEP edges. IntAna plane|plane is
// line-only and cannot consume these tags; the builder must ask here when
// a plane loop has no cylinder neighbor (or when IntAna misses).
//
// Accept iff max 3D vertex residual ≤ max(sewTol, Precision::Confusion()).
// chordSagitta(R,N) is not the cap: at N=6 sagitta is ~R(1−cos(π/6)) and
// garbage Pratt circles pass, then DIAG_CHAINSEW as ratioSew thousands.
// Stadium / slot mixed loops fail that residual gate and stay polylines.
// Regular hole ngons with verts already on a circle (R std ~0) still pass.
// Never force-polyline: a false from loopIsCircle is not an instruction to
// discard an IntAna curve the collapse already accepted.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <Precision.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
// prattCircleFit gates spanRad ∈ (0, π) for short arcs. A closed hole is 2π;
// the algebraic Pratt fit does not use spanRad after that check, so a dummy
// in-range value lets us reuse the in-TU solver for full loops.
constexpr double kPrattFullSpan = 1.0;

struct PlaneCircleTag {
    int rid = -1;
    int loopIx = -1;
    gp_Ax2 ax;
    double R = 0.0;
    std::vector<int> chains;
};

// Filled by tryPlaneLoopCircles; collapse reads it on the same thread.
// convert() fans out components, so this must not be process-global.
thread_local std::vector<PlaneCircleTag> gPlaneCircles;

gp_Pnt pntOf(const MeshView& mv, int lv) {
    if (!mv.pts || !mv.compVtx || lv < 0 || static_cast<size_t>(lv) >= mv.nVtx)
        return gp_Pnt(0.0, 0.0, 0.0);
    return gp_Pnt(mv.pts[mv.compVtx[lv]]);
}

bool sameLoop(const Loop& a, const Loop& b) {
    return a.role == b.role && a.chainIdx == b.chainIdx && a.reversed == b.reversed;
}

const Region* planeOwnerOf(const RegionSet& rs, const Loop& loop) {
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane) continue;
        for (const Loop& lp : r.loops) {
            if (sameLoop(lp, loop)) return &r;
        }
    }
    return nullptr;
}

// Ordered local verts around a loop, matching refit_chains loop assembly
// (drop the duplicated close). Sequential duplicates at chain joins skipped.
std::vector<int> loopLocalVerts(const RegionSet& rs, const Loop& lp) {
    std::vector<int> vs;
    const int nCh = static_cast<int>(std::min(lp.chainIdx.size(), lp.reversed.size()));
    for (int k = 0; k < nCh; ++k) {
        const int ci = lp.chainIdx[static_cast<size_t>(k)];
        if (ci < 0 || static_cast<size_t>(ci) >= rs.chains.size()) continue;
        const BoundaryChain& c = rs.chains[static_cast<size_t>(ci)];
        if (c.meshVerts.empty()) continue;
        const bool rev = lp.reversed[static_cast<size_t>(k)] != 0;
        if (c.closedLoop) {
            if (rev) {
                vs.push_back(c.meshVerts.front());
                for (int i = static_cast<int>(c.meshVerts.size()) - 1; i >= 1; --i)
                    vs.push_back(c.meshVerts[static_cast<size_t>(i)]);
            } else {
                vs.insert(vs.end(), c.meshVerts.begin(), c.meshVerts.end());
            }
        } else if (!rev) {
            const int begin = vs.empty() ? 0 : 1;
            for (int i = begin; i < static_cast<int>(c.meshVerts.size()); ++i)
                vs.push_back(c.meshVerts[static_cast<size_t>(i)]);
            if (vs.empty()) vs.push_back(c.meshVerts.front());
        } else {
            const int last = static_cast<int>(c.meshVerts.size()) - 1;
            const int begin = vs.empty() ? last : last - 1;
            for (int i = begin; i >= 0; --i) vs.push_back(c.meshVerts[static_cast<size_t>(i)]);
            if (vs.empty()) vs.push_back(c.meshVerts.back());
        }
    }
    if (vs.size() >= 2 && vs.front() == vs.back()) vs.pop_back();
    return vs;
}

// Unique by global welded point, first-occurrence order. N is this count.
std::vector<int> uniqueLoopVerts(const MeshView& mv, const std::vector<int>& vs) {
    std::vector<int> out;
    if (!mv.compVtx || mv.nVtx == 0) return out;
    std::vector<char> seenG;
    int maxG = -1;
    for (int lv : vs) {
        if (lv < 0 || static_cast<size_t>(lv) >= mv.nVtx) continue;
        const int g = mv.compVtx[lv];
        if (g > maxG) maxG = g;
    }
    if (maxG < 0) return out;
    seenG.assign(static_cast<size_t>(maxG) + 1, 0);
    out.reserve(vs.size());
    for (int lv : vs) {
        if (lv < 0 || static_cast<size_t>(lv) >= mv.nVtx) continue;
        const int g = mv.compVtx[lv];
        if (g < 0 || seenG[static_cast<size_t>(g)]) continue;
        seenG[static_cast<size_t>(g)] = 1;
        out.push_back(lv);
    }
    return out;
}

gp_Pnt projectOnPlane(const gp_Ax3& ax, const gp_Pnt& p) {
    const gp_Vec n(ax.Direction());
    const gp_Vec d(ax.Location(), p);
    return p.Translated(n * (-d.Dot(n)));
}

double circleResidual(const gp_Pnt& center, const gp_Dir& nrm, double R, const gp_Pnt& p) {
    const gp_Vec n(nrm);
    const gp_Vec d(center, p);
    const double off = std::fabs(d.Dot(n));
    const gp_Vec rho = d - n * d.Dot(n);
    const double r = rho.Magnitude();
    const double dr = std::fabs(r - R);
    return std::hypot(dr, off);
}

bool fitPlaneLoopCircle(const RegionSet& rs, const Loop& loop, const MeshView& mv,
                        const Region& plane, double sewTol, gp_Ax2& ax, double& ROut) {
    const std::vector<int> uniq = uniqueLoopVerts(mv, loopLocalVerts(rs, loop));
    const int N = static_cast<int>(uniq.size());
    if (N < 6) return false;
    if (!mv.pts || !mv.compVtx) return false;

    std::vector<gp_Pnt> pts;
    pts.reserve(uniq.size());
    for (int lv : uniq) pts.push_back(projectOnPlane(plane.ax, pntOf(mv, lv)));

    gp_Pnt center;
    double R = 0.0;
    if (!prattCircleFit(pts.data(), pts.size(), kPrattFullSpan, center, R)) return false;
    if (!(R > 0.0) || !std::isfinite(R)) return false;

    center = projectOnPlane(plane.ax, center);
    const gp_Dir nrm = plane.ax.Direction();

    double maxRes = 0.0;
    for (int lv : uniq) {
        const double res = circleResidual(center, nrm, R, pntOf(mv, lv));
        if (!std::isfinite(res)) return false;
        if (res > maxRes) maxRes = res;
    }

    const double accept = std::max(sewTol, Precision::Confusion());
    if (!(maxRes <= accept)) return false;

    gp_Dir xdir = plane.ax.XDirection();
    const gp_Vec xtry(center, pts.front());
    const gp_Vec xin = xtry - gp_Vec(nrm) * xtry.Dot(gp_Vec(nrm));
    if (xin.Magnitude() > 1e-15) xdir = gp_Dir(xin);

    ax = gp_Ax2(center, nrm, xdir);
    ROut = R;
    return true;
}

}  // namespace

bool loopIsCircle(const RegionSet& rs, const Loop& loop, const MeshView& mv, gp_Ax2& ax,
                 double& R, double sewTol) {
    const Region* plane = planeOwnerOf(rs, loop);
    if (!plane) return false;
    return fitPlaneLoopCircle(rs, loop, mv, *plane, sewTol, ax, R);
}

bool tryPlaneLoopCircles(RegionSet& rs, const MeshView& mv, double sewTol) {
    gPlaneCircles.clear();
    if (!mv.pts || !mv.compVtx) return false;

    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane) continue;
        for (int i = 0; i < static_cast<int>(r.loops.size()); ++i) {
            const Loop& lp = r.loops[static_cast<size_t>(i)];
            gp_Ax2 ax;
            double R = 0.0;
            if (!fitPlaneLoopCircle(rs, lp, mv, r, sewTol, ax, R)) continue;
            PlaneCircleTag tag;
            tag.rid = r.id;
            tag.loopIx = i;
            tag.ax = ax;
            tag.R = R;
            tag.chains = lp.chainIdx;
            gPlaneCircles.push_back(tag);
        }
    }
    // Tagging pass completed. Individual stadium/slot loops stay untagged.
    // Region has no per-loop circle field; the builder reads the thread_local
    // via planeLoopCircleOf / planeLoopCircleForChain, or refits with loopIsCircle.
    return true;
}

bool planeLoopCircleOf(int regionId, int loopIndex, gp_Ax2& ax, double& R) {
    for (const PlaneCircleTag& t : gPlaneCircles) {
        if (t.rid == regionId && t.loopIx == loopIndex) {
            ax = t.ax;
            R = t.R;
            return true;
        }
    }
    return false;
}

bool planeLoopCircleForChain(int chainIndex, gp_Ax2& ax, double& R) {
    if (chainIndex < 0) return false;
    for (const PlaneCircleTag& t : gPlaneCircles) {
        for (int ci : t.chains) {
            if (ci == chainIndex) {
                ax = t.ax;
                R = t.R;
                return true;
            }
        }
    }
    return false;
}

}  // namespace refit
}  // namespace stl2step
