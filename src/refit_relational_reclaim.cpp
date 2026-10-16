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

    // F8: one axis per revolution axis. A closed torus touching a ChamferCone (the
    // H0 donor frame) takes that cone's axis line, and the closed cylinders it
    // touches follow; each only while its own vertices stay within tau. Circle-
    // on-torus certifies exactly coaxial rims only (D-140-6 §3(4)).
    {
        std::vector<int> own((std::size_t)nTri, -1);
        for (int i = 0; i < nAcc; ++i)
            if (!dead[(std::size_t)i])
                for (int t : work.accepted[(std::size_t)i].tris)
                    if (t >= 0 && t < nTri) own[(std::size_t)t] = i;
        auto touching = [&](int s) {
            std::vector<int> out;
            for (int t : work.accepted[(std::size_t)s].tris) {
                if (t < 0 || t >= nTri) continue;
                for (int k = 0; k < 3; ++k) {
                    const int e = mv.triEdges[t][k];
                    if (e < 0 || e >= nEdge) continue;
                    const int o = edgeTris[(std::size_t)e][0] == t ? edgeTris[(std::size_t)e][1]
                                                                    : edgeTris[(std::size_t)e][0];
                    if (o >= 0 && own[(std::size_t)o] >= 0 && own[(std::size_t)o] != s)
                        out.push_back(own[(std::size_t)o]);
                }
            }
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
            return out;
        };
        auto onLine = [](const gp_Ax3& a, const gp_XYZ& lineP, const gp_XYZ& lineD) {
            gp_XYZ d = a.Direction().XYZ();
            gp_XYZ nd = lineD;
            if (d.Dot(nd) < 0.0) nd.Reverse();
            const gp_XYZ q = a.Location().XYZ();
            const gp_XYZ loc = lineP + nd * (q - lineP).Dot(nd);
            gp_XYZ x = a.XDirection().XYZ();
            x -= nd * x.Dot(nd);
            if (!(x.Modulus() > gp::Resolution())) return a;
            return gp_Ax3(gp_Pnt(loc), gp_Dir(nd), gp_Dir(x));
        };
        auto cylDev = [&](const Region& c, const gp_Ax3& ax) {
            double mx = 0.0;
            const gp_XYZ o = ax.Location().XYZ(), d = ax.Direction().XYZ();
            for (int t : c.tris)
                for (int k = 0; k < 3; ++k) {
                    const gp_XYZ v = cornerOf(mv, t, k) - o;
                    mx = std::max(mx, std::fabs((v - d * v.Dot(d)).Modulus() - c.radius));
                }
            return mx;
        };
        for (int s = 0; s < nAcc; ++s) {
            Region& T = work.accepted[(std::size_t)s];
            if (dead[(std::size_t)s] || T.origin != Origin::TorusBlend || !T.closed360) continue;
            const std::vector<int> nb = touching(s);
            int cone = -1;
            for (int i : nb) {
                const Region& c = work.accepted[(std::size_t)i];
                if (c.origin != Origin::ChamferCone) continue;
                const gp_XYZ cd = c.ax.Direction().XYZ();
                const gp_XYZ dd = T.ax.Location().XYZ() - c.ax.Location().XYZ();
                const double off = (dd - cd * dd.Dot(cd)).Modulus();
                const double ang = cd.Crossed(T.ax.Direction().XYZ()).Modulus();
                if (off <= tau && ang * (T.radius + T.radius2) <= tau) {
                    cone = i;
                    break;
                }
            }
            if (cone < 0) continue;
            const gp_XYZ lp = work.accepted[(std::size_t)cone].ax.Location().XYZ();
            const gp_XYZ ld = work.accepted[(std::size_t)cone].ax.Direction().XYZ();
            Region moved = T;
            moved.ax = onLine(T.ax, lp, ld);
            if (!(regionTorusDev(mv, moved, moved) <= tau)) continue;
            const bool changed = !(moved.ax.Location().IsEqual(T.ax.Location(), 0.0) &&
                                   moved.ax.Direction().IsEqual(T.ax.Direction(), 0.0));
            T.ax = moved.ax;
            for (int i : nb) {
                Region& c = work.accepted[(std::size_t)i];
                if (c.type != SurfType::Cylinder || !c.closed360) continue;
                const gp_Ax3 ca = onLine(c.ax, lp, ld);
                if (cylDev(c, ca) <= tau) c.ax = ca;
            }
            if (diag && changed)
                std::fprintf(stderr, "DIAG_REL_TORUS_AXIS acc=%d cone=%d\n", s, cone);
        }
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

JointAxisFit jointAxisAbout(const std::vector<std::vector<gp_XYZ>>& groups, const gp_XYZ& origin,
                            const gp_XYZ& axis) {
    JointAxisFit f;
    f.origin = origin;
    f.axis = axis;
    const double am = axis.Modulus();
    if (am > gp::Resolution()) f.axis /= am;
    for (const std::vector<gp_XYZ>& g : groups) {
        double sum = 0;
        std::vector<double> rho;
        rho.reserve(g.size());
        for (const gp_XYZ& p : g) {
            const gp_XYZ d = p - f.origin;
            rho.push_back((d - f.axis * d.Dot(f.axis)).Modulus());
            sum += rho.back();
        }
        const double R = g.empty() ? 0.0 : sum / (double)g.size();
        double dev = 0;
        for (double r : rho) dev = std::max(dev, std::fabs(r - R));
        f.radius.push_back(R);
        f.maxDev.push_back(dev);
    }
    return f;
}

bool jointAxisCertifies(const JointAxisFit& f, double tau) {
    if (!(tau > 0.0) || f.maxDev.empty()) return false;
    for (std::size_t g = 0; g < f.maxDev.size(); ++g)
        if (!(f.maxDev[g] <= tau) || !(f.radius[g] > tau)) return false;
    return true;
}

// Levenberg-Marquardt on (axis tilt along U, V; axis point along U, V), each
// donor's radius projected out as the mean of its own rho (variable projection,
// the same step cylinderFitLS takes for one radius).
bool jointAxisFit(const std::vector<std::vector<gp_XYZ>>& groups, const gp_XYZ& originSeed,
                  const gp_XYZ& axisSeed, int maxIters, JointAxisFit& out) {
    std::vector<gp_XYZ> P;
    std::vector<int> gOf;
    for (int g = 0; g < (int)groups.size(); ++g)
        for (const gp_XYZ& p : groups[(std::size_t)g]) {
            if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
                return false;
            P.push_back(p);
            gOf.push_back(g);
        }
    const int nG = (int)groups.size();
    if (nG == 0 || P.size() < 6) return false;
    gp_XYZ mean(0, 0, 0);
    for (const gp_XYZ& p : P) mean += p;
    mean /= (double)P.size();
    gp_XYZ W = axisSeed;
    if (!(W.Modulus() > gp::Resolution())) return false;
    W /= W.Modulus();
    const gp_XYZ Wseed = W;
    gp_XYZ C = originSeed;
    std::vector<double> rho(P.size()), Rg((std::size_t)nG), cnt((std::size_t)nG, 0.0);
    for (int g : gOf) cnt[(std::size_t)g] += 1.0;
    auto evaluate = [&](const gp_XYZ& Wc, const gp_XYZ& Cc, double& E) {
        std::fill(Rg.begin(), Rg.end(), 0.0);
        for (std::size_t i = 0; i < P.size(); ++i) {
            const gp_XYZ d = P[i] - Cc;
            const double m = (d - Wc * d.Dot(Wc)).Modulus();
            if (!std::isfinite(m) || m <= gp::Resolution()) return false;
            rho[i] = m;
            Rg[(std::size_t)gOf[i]] += m;
        }
        for (int g = 0; g < nG; ++g)
            if (cnt[(std::size_t)g] > 0) Rg[(std::size_t)g] /= cnt[(std::size_t)g];
        E = 0;
        for (std::size_t i = 0; i < P.size(); ++i) {
            const double f = rho[i] - Rg[(std::size_t)gOf[i]];
            E += f * f;
        }
        return std::isfinite(E);
    };
    double E = 0;
    if (!evaluate(W, C, E)) return false;
    double lambda = 1e-6;
    const int iterCap = maxIters > 0 ? maxIters : 16;
    for (int it = 0; it < iterCap; ++it) {
        const gp_XYZ tmp = (std::fabs(W.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
        gp_XYZ U = W.Crossed(tmp);
        U /= U.Modulus();
        const gp_XYZ V = W.Crossed(U);
        std::vector<std::array<double, 4>> J(P.size());
        std::vector<std::array<double, 4>> colMean((std::size_t)nG, {{0, 0, 0, 0}});
        for (std::size_t i = 0; i < P.size(); ++i) {
            const gp_XYZ d = P[i] - C;
            const double wd = d.Dot(W);
            const gp_XYZ rhat = (d - W * wd) / rho[i];
            J[i] = {{-wd * d.Dot(U) / rho[i], -wd * d.Dot(V) / rho[i], -rhat.Dot(U), -rhat.Dot(V)}};
            std::array<double, 4>& cm = colMean[(std::size_t)gOf[i]];
            for (int k = 0; k < 4; ++k) cm[(std::size_t)k] += J[i][(std::size_t)k];
        }
        for (int g = 0; g < nG; ++g)
            for (int k = 0; k < 4; ++k)
                if (cnt[(std::size_t)g] > 0)
                    colMean[(std::size_t)g][(std::size_t)k] /= cnt[(std::size_t)g];
        double A[4][4] = {{0}};
        double b[4] = {0, 0, 0, 0};
        for (std::size_t i = 0; i < P.size(); ++i) {
            double gr[4];
            const std::array<double, 4>& cm = colMean[(std::size_t)gOf[i]];
            for (int k = 0; k < 4; ++k) gr[k] = J[i][(std::size_t)k] - cm[(std::size_t)k];
            const double f = rho[i] - Rg[(std::size_t)gOf[i]];
            for (int k = 0; k < 4; ++k) {
                b[k] -= gr[k] * f;
                for (int l = 0; l < 4; ++l) A[k][l] += gr[k] * gr[l];
            }
        }
        bool stepped = false;
        for (int trial = 0; trial < 8 && !stepped; ++trial) {
            double M[4][5];
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) M[k][l] = A[k][l];
                M[k][k] = M[k][k] * (1.0 + lambda) + lambda;
                M[k][4] = b[k];
            }
            bool singular = false;
            for (int col = 0; col < 4 && !singular; ++col) {
                int piv = col;
                for (int row = col + 1; row < 4; ++row)
                    if (std::fabs(M[row][col]) > std::fabs(M[piv][col])) piv = row;
                if (!(std::fabs(M[piv][col]) > gp::Resolution())) {
                    singular = true;
                    break;
                }
                if (piv != col)
                    for (int l = 0; l < 5; ++l) std::swap(M[piv][l], M[col][l]);
                for (int row = col + 1; row < 4; ++row) {
                    const double fac = M[row][col] / M[col][col];
                    for (int l = col; l < 5; ++l) M[row][l] -= fac * M[col][l];
                }
            }
            if (singular) break;
            double dx[4];
            for (int row = 3; row >= 0; --row) {
                double acc = M[row][4];
                for (int l = row + 1; l < 4; ++l) acc -= M[row][l] * dx[l];
                dx[row] = acc / M[row][row];
            }
            if (!std::isfinite(dx[0]) || !std::isfinite(dx[1]) || !std::isfinite(dx[2]) ||
                !std::isfinite(dx[3]))
                break;
            gp_XYZ Wn = W + U * dx[0] + V * dx[1];
            if (!(Wn.Modulus() > gp::Resolution())) break;
            Wn /= Wn.Modulus();
            const gp_XYZ Cn = C + U * dx[2] + V * dx[3];
            double En = 0;
            const std::vector<double> rhoKeep = rho, RgKeep = Rg;
            if (evaluate(Wn, Cn, En) && En <= E) {
                W = Wn;
                C = Cn;
                E = En;
                lambda = std::max(1e-12, lambda / 3.0);
                stepped = true;
            } else {
                rho = rhoKeep;
                Rg = RgKeep;
                lambda *= 5.0;
                if (lambda > 1e6) break;
            }
        }
        if (!stepped) break;
    }
    if (W.Dot(Wseed) < 0.0) W.Reverse();
    C = C + W * ((mean - C).Dot(W));
    out = jointAxisAbout(groups, C, W);
    return std::isfinite(C.X()) && std::isfinite(C.Y()) && std::isfinite(C.Z());
}

std::vector<CoaxialFamily> coaxialFamilies(const std::vector<AxisDonor>& donors,
                                           const std::vector<int>& order, double tau,
                                           int maxIters) {
    std::vector<CoaxialFamily> fams;
    if (!(tau > 0.0)) return fams;
    std::vector<char> taken(donors.size(), 0);
    for (std::size_t oi = 0; oi < order.size(); ++oi) {
        const int i = order[oi];
        if (i < 0 || (std::size_t)i >= donors.size() || taken[(std::size_t)i]) continue;
        taken[(std::size_t)i] = 1;
        const AxisDonor& lead = donors[(std::size_t)i];
        CoaxialFamily fam;
        fam.members.push_back(i);
        std::vector<std::vector<gp_XYZ>> groups{lead.pts};
        fam.fit = jointAxisAbout(groups, lead.origin, lead.axis);
        double famLen = lead.length, famRmax = lead.radius;
        for (std::size_t oj = oi + 1; oj < order.size(); ++oj) {
            const int j = order[oj];
            if (j < 0 || (std::size_t)j >= donors.size() || taken[(std::size_t)j]) continue;
            const AxisDonor& e = donors[(std::size_t)j];
            gp_XYZ ea = e.axis;
            if (!(ea.Modulus() > gp::Resolution())) continue;
            ea /= ea.Modulus();
            const double sinA = ea.Crossed(fam.fit.axis).Modulus();
            if (sinA > tanAngTol(tau, tau, std::max(std::min(e.length, famLen), tau))) continue;
            const gp_XYZ d = e.origin - fam.fit.origin;
            const double lineDist = (d - fam.fit.axis * d.Dot(fam.fit.axis)).Modulus();
            if (lineDist > e.radius + famRmax) continue;
            groups.push_back(e.pts);
            JointAxisFit f = jointAxisAbout(groups, fam.fit.origin, fam.fit.axis);
            bool refit = false;
            if (!jointAxisCertifies(f, tau)) {
                JointAxisFit g;
                if (jointAxisFit(groups, fam.fit.origin, fam.fit.axis, maxIters, g) &&
                    jointAxisCertifies(g, tau)) {
                    f = g;
                    refit = true;
                } else {
                    groups.pop_back();
                    continue;
                }
            }
            fam.fit = f;
            fam.refit = fam.refit || refit;
            fam.members.push_back(j);
            taken[(std::size_t)j] = 1;
            famLen = std::max(famLen, e.length);
            famRmax = std::max(famRmax, e.radius);
        }
        fams.push_back(std::move(fam));
    }
    return fams;
}

bool torusChainSupersetFails(const std::vector<ProfPt>& candidatePts, std::size_t nCandidates,
                             std::size_t nBand, double tau) {
    if (nBand >= nCandidates) return false;
    return !certifyTorus(candidatePts, tau).ok;
}

}  // namespace refit
}  // namespace stl2step
