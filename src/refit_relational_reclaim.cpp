// Relational reclaim (D-fc-1.11(ii), D-fc-4.4, D-fc-4.5). Line mode runs at H0,
// right after claimRelRevolve(Line); Torus mode at H1 belongs to the tori lane.
// SPDX-License-Identifier: MIT
//
// Line mode, in order:
//   1. Narrow donor heal (D-fc-4.5). An accepted CylGrow region that the revolve
//      donor test refuses on its vertices (own LS residual > tau) is re-certified
//      on its tau-inliers I; every other triangle of it must be chord-bounded
//      cover. It then grows over free triangles from I, refitting and regrowing
//      to a fixed point, and takes a same-surface partial reached only through
//      the healed gap.
//   2. Same-surface closure (D-fc-4.4). Certified partial CylGrow regions that
//      are one surface at tau merge greedily to a fixed point; the Unclaimed
//      I/Cover triangles inside the group's extent fill the gap (D-fc-3.2.7).
//      The group commits only if the result is ONE edge-component and a disk.
// Every threshold is tau = 2q or a quantity measured on the certified set.

#include "refit_profile_cert.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTiny = 1e-30;

// D-fc-1.6(2): a free fit needs >= 2 axial levels.
constexpr int kMinAxialLevels = 2;
// D-fc-4.5: a cover triangle has >= 2 vertices within tau (the third within sDev).
constexpr int kCoverOnVertices = 2;

bool diagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_REL");
    return e && e[0] && e[0] != '0';
}

int minTriOf(const std::vector<int>& tris) {
    if (tris.empty()) return std::numeric_limits<int>::max();
    return *std::min_element(tris.begin(), tris.end());
}

gp_XYZ vtxOf(const MeshView& mv, int lt, int k) {
    const int gt = mv.compTris ? mv.compTris[lt] : lt;
    return mv.pts[mv.tris[gt][k]];
}

double triArea(const MeshView& mv, int lt) {
    const gp_XYZ a = vtxOf(mv, lt, 0), b = vtxOf(mv, lt, 1), c = vtxOf(mv, lt, 2);
    return 0.5 * (b - a).Crossed(c - a).Modulus();
}

gp_XYZ triCentroid(const MeshView& mv, int lt) {
    return (vtxOf(mv, lt, 0) + vtxOf(mv, lt, 1) + vtxOf(mv, lt, 2)) / 3.0;
}

gp_Dir triNormal(const MeshView& mv, int lt) {
    const gp_XYZ a = vtxOf(mv, lt, 0), b = vtxOf(mv, lt, 1), c = vtxOf(mv, lt, 2);
    gp_XYZ n = (b - a).Crossed(c - a);
    const double m = n.Modulus();
    if (m < kTiny) return gp_Dir(0, 0, 1);
    n /= m;
    return gp_Dir(n.X(), n.Y(), n.Z());
}

double wrapPi(double t) {
    while (t <= -kPi) t += kTwoPi;
    while (t > kPi) t -= kTwoPi;
    return t;
}

// A cylinder surface: axis point C, unit axis A, radius R.
struct Cyl {
    gp_XYZ C, A;
    double R = 0;
};

Cyl cylOf(const Region& r) {
    return {r.ax.Location().XYZ(), r.ax.Direction().XYZ(), r.radius};
}

double radialDev(const gp_XYZ& p, const Cyl& s) {
    return std::abs(s.A.Crossed(p - s.C).Modulus() - s.R);
}

double triDev(const MeshView& mv, int t, const Cyl& s) {
    double m = 0;
    for (int k = 0; k < 3; ++k) m = std::max(m, radialDev(vtxOf(mv, t, k), s));
    return m;
}

double setDev(const MeshView& mv, const std::vector<int>& ts, const Cyl& s) {
    double m = 0;
    for (int t : ts) m = std::max(m, triDev(mv, t, s));
    return m;
}

void axialRange(const MeshView& mv, const std::vector<int>& ts, const Cyl& s, double& lo,
                double& hi) {
    lo = std::numeric_limits<double>::infinity();
    hi = -std::numeric_limits<double>::infinity();
    for (int t : ts)
        for (int k = 0; k < 3; ++k) {
            const double v = (vtxOf(mv, t, k) - s.C).Dot(s.A);
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
}

// Axial levels of one triangle's three vertices at tau (D-fc-3.2.1 section ear).
int triLevels(const MeshView& mv, int t, const Cyl& s, double tau) {
    double v[3];
    for (int k = 0; k < 3; ++k) v[k] = (vtxOf(mv, t, k) - s.C).Dot(s.A);
    std::sort(v, v + 3);
    int n = 1;
    for (int k = 1; k < 3; ++k)
        if (v[k] - v[k - 1] > tau) ++n;
    return n;
}

// Mesh adjacency shared by both passes.
struct RMesh {
    const MeshView* mv = nullptr;
    std::vector<std::array<int, 2>> edgeTri;  // local edge -> two local tris
    std::vector<int> g2l;                     // global vertex -> local vertex

    int nbr(int t, int s) const {
        const int e = mv->triEdges[t][s];
        if (e < 0 || (std::size_t)e >= edgeTri.size()) return -1;
        return edgeTri[(std::size_t)e][0] == t ? edgeTri[(std::size_t)e][1]
                                               : edgeTri[(std::size_t)e][0];
    }
    int localVtx(int t, int k) const {
        const int gt = mv->compTris ? mv->compTris[t] : t;
        const int g = mv->tris[gt][k];
        return (g >= 0 && (std::size_t)g < g2l.size()) ? g2l[(std::size_t)g] : -1;
    }
};

RMesh makeRMesh(const MeshView& mv) {
    RMesh rm;
    rm.mv = &mv;
    const int nTri = (int)mv.nTri;
    rm.edgeTri.assign(mv.nEdge, std::array<int, 2>{{-1, -1}});
    for (int t = 0; t < nTri; ++t)
        for (int s = 0; s < 3; ++s) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || (std::size_t)e >= mv.nEdge) continue;
            auto& et = rm.edgeTri[(std::size_t)e];
            if (et[0] < 0) et[0] = t;
            else if (et[0] != t) et[1] = t;
        }
    if (mv.compVtx) {
        int gMax = -1;
        for (std::size_t i = 0; i < mv.nVtx; ++i) gMax = std::max(gMax, mv.compVtx[i]);
        rm.g2l.assign((std::size_t)(gMax + 1), -1);
        for (std::size_t i = 0; i < mv.nVtx; ++i) {
            const int g = mv.compVtx[i];
            if (g >= 0 && rm.g2l[(std::size_t)g] < 0) rm.g2l[(std::size_t)g] = (int)i;
        }
    }
    return rm;
}

// Distinct vertex points of a triangle set, ascending local vertex id.
std::vector<std::pair<int, gp_XYZ>> uniqueVerts(const RMesh& rm, const std::vector<int>& tris) {
    std::vector<std::pair<int, gp_XYZ>> u;
    u.reserve(tris.size() * 3);
    for (int t : tris)
        for (int k = 0; k < 3; ++k) {
            const int lv = rm.localVtx(t, k);
            if (lv < 0) continue;
            u.emplace_back(lv, vtxOf(*rm.mv, t, k));
        }
    std::sort(u.begin(), u.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    u.erase(std::unique(u.begin(), u.end(),
                        [](const auto& a, const auto& b) { return a.first == b.first; }),
            u.end());
    return u;
}

// D-fc-1.4 over-determination of a partial cylinder: azimuth columns and axial
// levels of the set's vertices, merged at tau (the revolve donor test's counts).
void stationsOf(const RMesh& rm, const std::vector<int>& tris, const Cyl& s, double tau,
                int& cols, int& levels) {
    const auto u = uniqueVerts(rm, tris);
    gp_XYZ e1 = std::abs(s.A.X()) < 0.9 ? gp_XYZ(1, 0, 0).Crossed(s.A)
                                       : gp_XYZ(0, 1, 0).Crossed(s.A);
    e1 /= e1.Modulus();
    const gp_XYZ e2 = s.A.Crossed(e1);
    std::vector<double> ps, rs, zs;
    for (const auto& v : u) {
        const gp_XYZ d = v.second - s.C;
        const double x = d.Dot(e1), y = d.Dot(e2);
        double psi = std::atan2(y, x);
        if (psi < 0.0) psi += kTwoPi;
        ps.push_back(psi);
        rs.push_back(std::hypot(x, y));
        zs.push_back(d.Dot(s.A));
    }
    cols = countColumns(ps, rs, tau);
    levels = countLevels(zs, tau);
}

bool overDetermined(int cols, int levels) {
    return cols >= minStations(kParamsCircle) && levels >= kMinAxialLevels;
}

// Region fields by the existing rules (D-fc-1.6(13)): B1's D2 metrics (normal-azimuth
// bands, vertex-azimuth span, closed360, nSides), the sigma of computeOutwardCylinder,
// chordSagitta and dVolCylinderSector over the region's own triangles.
struct CylFields {
    gp_Ax3 ax;
    double uMin = 0, uMax = 0, vMin = 0, vMax = 0;
    bool closed360 = false;
    int nSides = 0;
    bool outward = true;
    double area = 0;
    double maxDev = 0, rmsDev = 0;
};

double percentile75(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[(std::size_t)std::floor(0.75 * (double)(v.size() - 1))];
}

bool cylFields(const RMesh& rm, const DerivedTols& tol, const std::vector<int>& tris,
               const Cyl& s, CylFields& f) {
    const MeshView& mv = *rm.mv;
    if (tris.empty()) return false;
    const gp_XYZ aw = s.A;
    gp_XYZ centroid(0, 0, 0);
    double totalArea = 0.0;
    for (int lt : tris) {
        const double ar = triArea(mv, lt);
        centroid += triCentroid(mv, lt) * ar;
        totalArea += ar;
    }
    if (totalArea > kTiny) centroid /= totalArea;
    const gp_XYZ loc = s.C + aw * (centroid - s.C).Dot(aw);

    const auto uniq = uniqueVerts(rm, tris);
    if (uniq.size() < 2) return false;
    // XDirection: the region vertex with the lowest local id (refit.hpp Region::ax).
    gp_XYZ xDir = uniq.front().second - loc;
    xDir -= aw * xDir.Dot(aw);
    if (xDir.Modulus() < 1e-12)
        xDir = aw.Crossed((std::abs(aw.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0));
    xDir /= xDir.Modulus();
    f.ax = gp_Ax3(gp_Pnt(loc), gp_Dir(aw), gp_Dir(xDir));
    const gp_XYZ xAx = f.ax.XDirection().XYZ(), yAx = f.ax.YDirection().XYZ();

    std::vector<double> psi;
    psi.reserve(tris.size());
    for (int lt : tris) {
        const gp_Dir n = triNormal(mv, lt);
        psi.push_back(wrapPi(std::atan2(n.XYZ().Dot(yAx), n.XYZ().Dot(xAx))));
    }
    std::sort(psi.begin(), psi.end());
    int nBands = 1;
    if (psi.size() >= 2) {
        std::vector<double> gaps;
        for (std::size_t i = 1; i < psi.size(); ++i) gaps.push_back(psi[i] - psi[i - 1]);
        gaps.push_back(kTwoPi - psi.back() + psi.front());
        gaps.erase(std::max_element(gaps.begin(), gaps.end()));
        const double thetaBin = std::max(tol.thetaBin, 0.5 * percentile75(gaps));
        for (std::size_t i = 1; i < psi.size(); ++i)
            if (psi[i] - psi[i - 1] > thetaBin) ++nBands;
    }

    std::vector<double> chi;
    chi.reserve(uniq.size());
    for (const auto& uv : uniq) {
        const gp_XYZ d = uv.second - loc;
        chi.push_back(wrapPi(std::atan2(d.Dot(yAx), d.Dot(xAx))));
    }
    std::sort(chi.begin(), chi.end());
    const std::size_t m = chi.size();
    int jmax = 0;
    double maxGap = chi[1] - chi[0];
    for (std::size_t j = 1; j + 1 < m; ++j) {
        const double g = chi[j + 1] - chi[j];
        if (g > maxGap) {
            maxGap = g;
            jmax = (int)j;
        }
    }
    const double wrapGap = kTwoPi - chi.back() + chi.front();
    if (wrapGap > maxGap) {
        maxGap = wrapGap;
        jmax = (int)m - 1;
    }
    const double bandArc = kTwoPi / (double)nBands;
    f.closed360 = nBands >= 3 && maxGap <= 1.5 * bandArc;
    const double span = f.closed360 ? kTwoPi : (kTwoPi - maxGap);
    if (!(span > 0)) return false;
    f.nSides = (int)std::llround(kTwoPi * nBands / span);
    if (f.closed360) {
        f.uMin = 0.0;
        f.uMax = kTwoPi;
    } else {
        f.uMin = wrapPi(chi[(std::size_t)(jmax + 1) % m]);
        f.uMax = f.uMin + span;
    }
    f.vMin = std::numeric_limits<double>::infinity();
    f.vMax = -std::numeric_limits<double>::infinity();
    double mx = 0, s2 = 0;
    for (const auto& uv : uniq) {
        const double v = (uv.second - loc).Dot(aw);
        f.vMin = std::min(f.vMin, v);
        f.vMax = std::max(f.vMax, v);
    }
    int nv = 0;
    for (int lt : tris)
        for (int k = 0; k < 3; ++k) {
            const double d = radialDev(vtxOf(mv, lt, k), s);
            mx = std::max(mx, d);
            s2 += d * d;
            ++nv;
        }
    f.maxDev = mx;
    f.rmsDev = nv ? std::sqrt(s2 / nv) : 0.0;

    double sigma = 0.0;
    for (int lt : tris) {
        const gp_XYZ d = triCentroid(mv, lt) - loc;
        gp_XYZ radial = d - aw * d.Dot(aw);
        const double rmod = radial.Modulus();
        if (rmod < kTiny) continue;
        radial /= rmod;
        sigma += triArea(mv, lt) * triNormal(mv, lt).XYZ().Dot(radial);
    }
    f.outward = sigma > 0.0;
    f.area = totalArea;
    return true;
}

void applyCylFields(Region& r, const Cyl& s, const CylFields& f) {
    r.type = SurfType::Cylinder;
    r.ax = f.ax;
    r.radius = s.R;
    r.uMin = f.uMin;
    r.uMax = f.uMax;
    r.vMin = f.vMin;
    r.vMax = f.vMax;
    r.closed360 = f.closed360;
    r.nSides = f.nSides;
    r.chordSagitta = chordSagitta(s.R, f.nSides);
    r.outwardNormal = f.outward;
    r.dVolPredicted = dVolCylinderSector(f.area, s.R, f.nSides, f.outward);
    r.maxVertexDev = f.maxDev;
    r.rmsVertexDev = f.rmsDev;
}

std::vector<std::vector<int>> edgeComponents(const RMesh& rm, const std::vector<int>& ts);

// A2 plane fields of a provisional from its own triangles (PCA plane, area, vertex
// residuals), as the straddle trim of the other relational stages does.
void refitProvisional(const RMesh& rm, Provisional& pr) {
    const MeshView& mv = *rm.mv;
    std::sort(pr.tris.begin(), pr.tris.end());
    gp_Ax3 pl = pr.plane;
    if (pcaPlane(mv, pr.tris, pl)) pr.plane = pl;
    double area = 0;
    for (int t : pr.tris) area += triArea(mv, t);
    pr.area = area;
    const gp_XYZ n = pr.plane.Direction().XYZ(), o = pr.plane.Location().XYZ();
    double mx = 0, s2 = 0;
    const auto u = uniqueVerts(rm, pr.tris);
    for (const auto& v : u) {
        const double d = std::abs((v.second - o).Dot(n));
        mx = std::max(mx, d);
        s2 += d * d;
    }
    pr.maxVertexDev = mx;
    pr.rmsVertexDev = u.empty() ? 0.0 : std::sqrt(s2 / (double)u.size());
}

// D-fc-1.6(12) straddle trim: claimed triangles leave the Unclaimed provisionals.
// An emptied provisional is consumed; a remainder is split into edge-connected
// components (the first keeps the provisional, the rest are new Unclaimed
// provisionals), each refit on its own triangles. Returns the straddle count.
int trimClaimed(const RMesh& rm, SegmentWork& work, const std::vector<char>& inC) {
    int straddle = 0;
    const int nProv = (int)work.provisionals.size();
    for (int i = 0; i < nProv; ++i) {
        Provisional& pr = work.provisionals[(std::size_t)i];
        if (pr.claim != ProvClaim::Unclaimed) continue;
        std::vector<int> keep;
        bool drop = false;
        for (int t : pr.tris) {
            if (t >= 0 && (std::size_t)t < inC.size() && inC[(std::size_t)t]) drop = true;
            else keep.push_back(t);
        }
        if (!drop) continue;
        if (keep.empty()) {
            pr.claim = ProvClaim::ConsumedCylinder;
            continue;
        }
        ++straddle;
        auto comps = edgeComponents(rm, keep);
        std::sort(comps.begin(), comps.end(), [](const std::vector<int>& a,
                                                 const std::vector<int>& b) {
            return minTriOf(a) < minTriOf(b);
        });
        pr.tris = comps[0];
        refitProvisional(rm, pr);
        for (std::size_t c = 1; c < comps.size(); ++c) {
            Provisional np = work.provisionals[(std::size_t)i];
            np.tris = comps[c];
            np.claim = ProvClaim::Unclaimed;
            refitProvisional(rm, np);
            work.provisionals.push_back(std::move(np));
        }
    }
    return straddle;
}

std::vector<std::vector<int>> edgeComponents(const RMesh& rm, const std::vector<int>& ts) {
    const int nTri = (int)rm.mv->nTri;
    std::vector<char> in((std::size_t)nTri, 0), seen((std::size_t)nTri, 0);
    for (int t : ts) in[(std::size_t)t] = 1;
    std::vector<std::vector<int>> out;
    for (int t0 : ts) {
        if (seen[(std::size_t)t0]) continue;
        std::vector<int> c, st{t0};
        seen[(std::size_t)t0] = 1;
        while (!st.empty()) {
            const int x = st.back();
            st.pop_back();
            c.push_back(x);
            for (int s = 0; s < 3; ++s) {
                const int u = rm.nbr(x, s);
                if (u < 0 || !in[(std::size_t)u] || seen[(std::size_t)u]) continue;
                seen[(std::size_t)u] = 1;
                st.push_back(u);
            }
        }
        std::sort(c.begin(), c.end());
        out.push_back(std::move(c));
    }
    return out;
}

// Connected components of the boundary-edge graph of a triangle set; -1 when the
// mesh view carries no edge endpoints.
int boundaryLoops(const RMesh& rm, const std::vector<int>& ts) {
    const MeshView& mv = *rm.mv;
    if (!mv.compEdges) return -1;
    std::map<int, int> cnt;
    for (int t : ts)
        for (int s = 0; s < 3; ++s) ++cnt[mv.triEdges[t][s]];
    std::map<int, std::vector<int>> adj;
    for (const auto& kv : cnt) {
        if (kv.second != 1) continue;
        const auto& e = mv.compEdges[kv.first];
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    std::set<int> seen;
    int loops = 0;
    for (const auto& kv : adj) {
        if (seen.count(kv.first)) continue;
        ++loops;
        std::vector<int> st{kv.first};
        seen.insert(kv.first);
        while (!st.empty()) {
            const int x = st.back();
            st.pop_back();
            for (int y : adj[x])
                if (seen.insert(y).second) st.push_back(y);
        }
    }
    return loops;
}

// Partial-arc membership with the largest-gap convention of Region::uMin/uMax.
bool inArc(double chi, double uMin, double uMax, double slack) {
    double x = chi;
    while (x < uMin - slack) x += kTwoPi;
    while (x > uMin - slack + kTwoPi) x -= kTwoPi;
    return x <= uMax + slack;
}

double azimuthOf(const gp_XYZ& p, const Cyl& s, const gp_XYZ& xd) {
    const gp_XYZ d = p - s.C;
    return std::atan2(d.Dot(s.A.Crossed(xd)), d.Dot(xd));
}

struct BandDiag {
    const char* mode = "";
    int seed = -1, nTri = 0, nV = 0, stations = 0, levels = 0, members = 0;
    int nI = 0, nCover = 0, comps = 0, loops = 0, nSides = 0, straddle = 0;
    double dVol = 0;
    double radius = 0, maxDev = 0, maxDevTau = 0, certDevTau = 0, sTau = 0;
    gp_XYZ axis{0, 0, 0};
    const char* why = "";
};

void printBand(const BandDiag& d) {
    if (!diagOn()) return;
    std::fprintf(stderr,
                 "DIAG_REL_BAND stage=reclaim mode=%s seed=%d nTri=%d nV=%d stations=%d "
                 "levels=%d params=%d fit=free radii=%.9g axis=%.6g,%.6g,%.6g maxDev=%.6g "
                 "maxDevTau=%.4f certDevTau=%.4f sTau=%.4f members=%d nI=%d nCover=%d "
                 "comps=%d loops=%d nSides=%d dVol=%.6g straddle=%d why=%s\n",
                 d.mode, d.seed, d.nTri, d.nV, d.stations, d.levels, kParamsCircle + kParamsAxisDir,
                 d.radius, d.axis.X(), d.axis.Y(), d.axis.Z(), d.maxDev, d.maxDevTau,
                 d.certDevTau, d.sTau, d.members, d.nI, d.nCover, d.comps, d.loops, d.nSides,
                 d.dVol, d.straddle, d.why);
}

// Erase the accepted regions this stage absorbed, highest index first.
void eraseRegions(SegmentWork& work, std::vector<int> idx) {
    std::sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
    for (auto it = idx.rbegin(); it != idx.rend(); ++it)
        work.accepted.erase(work.accepted.begin() + *it);
}

bool isCylGrowCylinder(const Region& r) {
    return r.type == SurfType::Cylinder && r.origin == Origin::CylGrow && !r.tris.empty();
}

// ---- D-fc-4.5: narrow donor heal ------------------------------------------------
int healRefusedDonors(const RMesh& rm, const SegmentParams& p, const DerivedTols& tol,
                      double tau, SegmentWork& work) {
    const MeshView& mv = *rm.mv;
    const int nTri = (int)mv.nTri;
    const int maxIt = std::max(1, p.maxRefineIters);
    std::vector<int> owner((std::size_t)nTri, -1), triProv((std::size_t)nTri, -1);
    for (int i = 0; i < (int)work.accepted.size(); ++i)
        for (int t : work.accepted[(std::size_t)i].tris)
            if (t >= 0 && t < nTri) owner[(std::size_t)t] = i;
    for (int i = 0; i < (int)work.provisionals.size(); ++i)
        for (int t : work.provisionals[(std::size_t)i].tris)
            if (t >= 0 && t < nTri) triProv[(std::size_t)t] = i;
    auto isFree = [&](int t) {
        if (t < 0 || t >= nTri || owner[(std::size_t)t] >= 0) return false;
        const int pi = triProv[(std::size_t)t];
        return pi >= 0 && work.provisionals[(std::size_t)pi].claim == ProvClaim::Unclaimed;
    };
    std::vector<int> order;
    for (int i = 0; i < (int)work.accepted.size(); ++i)
        if (isCylGrowCylinder(work.accepted[(std::size_t)i])) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return minTriOf(work.accepted[(std::size_t)a].tris) <
               minTriOf(work.accepted[(std::size_t)b].tris);
    });

    int healed = 0;
    std::vector<int> absorbedRegions;
    std::vector<char> inSet((std::size_t)nTri, 0);
    for (int ri : order) {
        Region& R = work.accepted[(std::size_t)ri];
        if (R.tris.empty()) continue;
        gp_Dir ax;
        gp_Pnt c;
        double rad = 0, mres = 0;
        if (!cylinderFitLS(mv, R.tris, R.ax.Direction(), ax, c, rad, mres, maxIt)) continue;
        // Narrow trigger: the donor test refuses this region on its vertices.
        // A region that certifies is never touched (the broad variant is refused).
        if (mres <= tau) continue;

        BandDiag dg;
        dg.mode = "heal";
        dg.seed = minTriOf(R.tris);
        Cyl S{c.XYZ(), ax.XYZ(), rad};
        auto inliersOf = [&](const Cyl& s, double& mI) {
            std::vector<int> I;
            mI = 0;
            for (int t : R.tris) {
                const double d = triDev(mv, t, s);
                if (d <= tau) {
                    I.push_back(t);
                    mI = std::max(mI, d);
                }
            }
            return I;
        };
        double mI = 0;
        std::vector<int> I = inliersOf(S, mI);
        if (I.empty() || !cylinderFitLS(mv, I, ax, ax, c, rad, mres, maxIt)) {
            dg.why = "underdetermined";
            printBand(dg);
            continue;
        }
        S = {c.XYZ(), ax.XYZ(), rad};
        I = inliersOf(S, mI);
        int cols = 0, levels = 0;
        stationsOf(rm, I, S, tau, cols, levels);
        dg.stations = cols;
        dg.levels = levels;
        dg.radius = rad;
        dg.axis = S.A;
        dg.nI = (int)I.size();
        dg.certDevTau = mI / tau;
        if (I.empty() || !overDetermined(cols, levels)) {
            dg.why = "underdetermined";
            printBand(dg);
            continue;
        }
        // Chord bound of the certified band: the max mid-edge radial deviation of
        // I's own edges (D-fc-3.2.7).
        auto chordBound = [&](const Cyl& s, const std::vector<int>& set) {
            double sd = 0;
            for (int t : set)
                for (int k = 0; k < 3; ++k)
                    sd = std::max(sd, radialDev(0.5 * (vtxOf(mv, t, k) + vtxOf(mv, t, (k + 1) % 3)), s));
            return sd;
        };
        const double sDev = chordBound(S, I);
        dg.sTau = sDev / tau;
        auto coverOk = [&](int t, const Cyl& s) {
            int on = 0;
            for (int k = 0; k < 3; ++k) {
                const double d = radialDev(vtxOf(mv, t, k), s);
                if (d <= tau) ++on;
                if (d > std::max(sDev, tau)) return false;
            }
            return on >= kCoverOnVertices;
        };
        std::fill(inSet.begin(), inSet.end(), 0);
        for (int t : I) inSet[(std::size_t)t] = 1;
        int nCov = 0;
        bool covered = true;
        for (int t : R.tris) {
            if (inSet[(std::size_t)t]) continue;
            if (!coverOk(t, S)) {
                covered = false;
                break;
            }
            ++nCov;
        }
        dg.nCover = nCov;
        if (!covered) {
            dg.why = "vertex-off";
            printBand(dg);
            continue;
        }

        // Grow over free triangles from I (inliers chain, chord-bounded cover is
        // traversed); refit on I + grown inliers and regrow until the set is stable.
        std::vector<int> absorbed, cov;
        std::vector<char> seen((std::size_t)nTri, 0);
        bool stable = false;
        for (int it = 0; it < maxIt; ++it) {
            std::fill(seen.begin(), seen.end(), 0);
            for (int t : R.tris) seen[(std::size_t)t] = 1;
            std::vector<int> stack = I, ab, cv;
            while (!stack.empty()) {
                const int u = stack.back();
                stack.pop_back();
                for (int s = 0; s < 3; ++s) {
                    const int w = rm.nbr(u, s);
                    if (w < 0 || seen[(std::size_t)w] || !isFree(w)) continue;
                    if (triDev(mv, w, S) <= tau) {
                        seen[(std::size_t)w] = 1;
                        ab.push_back(w);
                        stack.push_back(w);
                    } else if (coverOk(w, S)) {
                        seen[(std::size_t)w] = 1;
                        cv.push_back(w);
                        stack.push_back(w);
                    }
                }
            }
            std::sort(ab.begin(), ab.end());
            std::sort(cv.begin(), cv.end());
            const bool same = (ab == absorbed && cv == cov);
            absorbed.swap(ab);
            cov.swap(cv);
            if (same || absorbed.empty()) {
                stable = true;
                break;
            }
            std::vector<int> fitSet(I);
            fitSet.insert(fitSet.end(), absorbed.begin(), absorbed.end());
            gp_Dir ax2;
            gp_Pnt c2;
            double r2 = 0, m2 = 0;
            if (!cylinderFitLS(mv, fitSet, ax, ax2, c2, r2, m2, maxIt)) {
                stable = true;
                break;
            }
            const Cyl S2{c2.XYZ(), ax2.XYZ(), r2};
            // The refit must certify the whole fit set, else the current fit stands.
            if (!(setDev(mv, fitSet, S2) <= tau)) {
                stable = true;
                break;
            }
            ax = ax2;
            c = c2;
            rad = r2;
            S = S2;
            I = inliersOf(S, mI);
        }
        if (!stable) {
            dg.why = "unstable";
            printBand(dg);
            continue;
        }
        if (absorbed.empty()) {
            dg.why = "not-joined";
            printBand(dg);
            continue;
        }

        // A same-surface partial reached only through the healed gap joins (D-130-14).
        std::fill(inSet.begin(), inSet.end(), 0);
        for (int t : R.tris)
            for (int s = 0; s < 3; ++s) {
                const int w = rm.nbr(t, s);
                if (w >= 0) inSet[(std::size_t)w] = 1;
            }
        std::vector<int> partners;
        for (int u : absorbed)
            for (int s = 0; s < 3; ++s) {
                const int w = rm.nbr(u, s);
                if (w < 0 || owner[(std::size_t)w] < 0 || owner[(std::size_t)w] == ri) continue;
                const int qi = owner[(std::size_t)w];
                const Region& Q = work.accepted[(std::size_t)qi];
                if (!isCylGrowCylinder(Q)) continue;
                if (std::find(partners.begin(), partners.end(), qi) != partners.end()) continue;
                bool adjBefore = false;
                for (int t : Q.tris)
                    if (inSet[(std::size_t)t]) {
                        adjBefore = true;
                        break;
                    }
                if (adjBefore) continue;
                std::vector<int> qIn;
                bool ok = true;
                for (int t : Q.tris) {
                    if (triDev(mv, t, S) <= tau) qIn.push_back(t);
                    else if (!coverOk(t, S)) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;
                int qc = 0, ql = 0;
                stationsOf(rm, qIn, S, tau, qc, ql);
                if (overDetermined(qc, ql)) partners.push_back(qi);
            }

        std::vector<int> grown = R.tris;
        grown.insert(grown.end(), absorbed.begin(), absorbed.end());
        grown.insert(grown.end(), cov.begin(), cov.end());
        int nPartTris = 0;
        for (int qi : partners) {
            const auto& qt = work.accepted[(std::size_t)qi].tris;
            nPartTris += (int)qt.size();
            grown.insert(grown.end(), qt.begin(), qt.end());
        }
        std::sort(grown.begin(), grown.end());
        grown.erase(std::unique(grown.begin(), grown.end()), grown.end());
        CylFields f;
        if (!cylFields(rm, tol, grown, S, f)) {
            dg.why = "degenerate";
            printBand(dg);
            continue;
        }
        std::vector<char> inAdd((std::size_t)nTri, 0);
        for (int t : absorbed) inAdd[(std::size_t)t] = 1;
        for (int t : cov) inAdd[(std::size_t)t] = 1;
        dg.straddle = trimClaimed(rm, work, inAdd);
        for (int qi : partners) {
            work.accepted[(std::size_t)qi].tris.clear();
            absorbedRegions.push_back(qi);
        }
        for (int t : grown) owner[(std::size_t)t] = ri;
        R.tris = std::move(grown);
        applyCylFields(R, S, f);
        ++healed;
        dg.nTri = (int)R.tris.size();
        dg.nV = (int)uniqueVerts(rm, R.tris).size();
        dg.members = 1 + (int)partners.size();
        dg.nI = (int)I.size() + (int)absorbed.size();
        dg.nCover = nCov + (int)cov.size();
        dg.maxDev = R.maxVertexDev;
        dg.maxDevTau = R.maxVertexDev / tau;
        dg.nSides = R.nSides;
        dg.dVol = R.dVolPredicted;
        dg.certDevTau = mI / tau;
        dg.radius = rad;
        dg.axis = S.A;
        dg.comps = (int)edgeComponents(rm, R.tris).size();
        dg.loops = boundaryLoops(rm, R.tris);
        dg.why = "accept";
        printBand(dg);
        if (diagOn())
            std::fprintf(stderr, "DIAG_REL_HEAL seed=%d absorbed=%d absorbCover=%d partners=%d "
                                 "partnerTris=%d closed360=%d\n",
                         dg.seed, (int)absorbed.size(), (int)cov.size(), (int)partners.size(),
                         nPartTris, (int)R.closed360);
    }
    eraseRegions(work, absorbedRegions);
    return healed;
}

// ---- D-fc-4.4: same-surface closure ---------------------------------------------
int closeSameSurface(const RMesh& rm, const SegmentParams& p, const DerivedTols& tol,
                     double tau, SegmentWork& work, int& deferred) {
    const MeshView& mv = *rm.mv;
    const int nTri = (int)mv.nTri;
    const int maxIt = std::max(1, p.maxRefineIters);
    std::vector<int> cand;
    for (int i = 0; i < (int)work.accepted.size(); ++i) {
        const Region& r = work.accepted[(std::size_t)i];
        if (!isCylGrowCylinder(r) || r.closed360 || !(r.maxVertexDev <= tau)) continue;
        int cols = 0, levels = 0;
        stationsOf(rm, r.tris, cylOf(r), tau, cols, levels);
        if (overDetermined(cols, levels)) cand.push_back(i);
    }
    struct Group {
        std::vector<int> members;  // accepted indices
        std::vector<int> tris;
        Cyl s;
    };
    std::vector<Group> groups;
    for (int i : cand) groups.push_back({{i}, work.accepted[(std::size_t)i].tris,
                                         cylOf(work.accepted[(std::size_t)i])});
    // Greedy in ascending order to a fixed point: mutual vertex containment at tau,
    // strictly overlapping axial extents (stacked coaxial walls stay apart), joint
    // LS residual at tau.
    bool merged = true;
    while (merged) {
        merged = false;
        for (std::size_t i = 0; i < groups.size() && !merged; ++i) {
            if (groups[i].members.empty()) continue;
            for (std::size_t j = i + 1; j < groups.size(); ++j) {
                if (groups[j].members.empty()) continue;
                Group& gi = groups[i];
                Group& gj = groups[j];
                if (setDev(mv, gj.tris, gi.s) > tau || setDev(mv, gi.tris, gj.s) > tau) continue;
                double li, hi, lj, hj;
                axialRange(mv, gi.tris, gi.s, li, hi);
                axialRange(mv, gj.tris, gi.s, lj, hj);
                if (!(li < hj && lj < hi)) continue;
                std::vector<int> trial = gi.tris;
                trial.insert(trial.end(), gj.tris.begin(), gj.tris.end());
                std::sort(trial.begin(), trial.end());
                gp_Dir fa;
                gp_Pnt fc;
                double fr = 0, fres = 0;
                if (!cylinderFitLS(mv, trial, gp_Dir(gi.s.A), fa, fc, fr, fres, maxIt) ||
                    fres > tau)
                    continue;
                gi.members.insert(gi.members.end(), gj.members.begin(), gj.members.end());
                gi.tris = std::move(trial);
                gi.s = {fc.XYZ(), fa.XYZ(), fr};
                gj.members.clear();
                gj.tris.clear();
                merged = true;
                break;
            }
        }
    }

    std::vector<char> unclaimed((std::size_t)nTri, 0);
    for (const Provisional& pr : work.provisionals)
        if (pr.claim == ProvClaim::Unclaimed)
            for (int t : pr.tris)
                if (t >= 0 && t < nTri) unclaimed[(std::size_t)t] = 1;
    for (const Region& r : work.accepted)
        for (int t : r.tris)
            if (t >= 0 && t < nTri) unclaimed[(std::size_t)t] = 0;

    int committed = 0;
    std::vector<int> erase;
    for (Group& g : groups) {
        if (g.members.size() < 2) continue;  // one surface => one face triggers the fill
        const Cyl S = g.s;
        const Region& R0 = work.accepted[(std::size_t)g.members.front()];
        gp_XYZ xd = R0.ax.XDirection().XYZ();
        xd -= S.A * xd.Dot(S.A);
        xd /= xd.Modulus();
        // The group's own angular extent (largest-gap convention) and axial extent.
        std::vector<double> chi;
        double vlo, vhi;
        axialRange(mv, g.tris, S, vlo, vhi);
        for (int t : g.tris)
            for (int k = 0; k < 3; ++k) chi.push_back(azimuthOf(vtxOf(mv, t, k), S, xd));
        std::sort(chi.begin(), chi.end());
        double maxGap = 0;
        std::size_t gapAt = 0;
        for (std::size_t i = 0; i + 1 < chi.size(); ++i)
            if (chi[i + 1] - chi[i] > maxGap) {
                maxGap = chi[i + 1] - chi[i];
                gapAt = i;
            }
        const double wrap = chi.front() + kTwoPi - chi.back();
        if (wrap > maxGap) {
            maxGap = wrap;
            gapAt = chi.size() - 1;
        }
        const double u0 = wrapPi(chi[(gapAt + 1) % chi.size()]);
        const double u1 = u0 + (kTwoPi - maxGap);
        const double s = claimChordSagitta(mv, g.tris, gp_Ax1(gp_Pnt(S.C), gp_Dir(S.A)), S.R);

        std::vector<int> I, Cv;
        for (int t = 0; t < nTri; ++t) {
            if (!unclaimed[(std::size_t)t]) continue;
            int nin = 0;
            double md = 0;
            bool inside = true;
            for (int k = 0; k < 3 && inside; ++k) {
                const gp_XYZ pt = vtxOf(mv, t, k);
                const double d = radialDev(pt, S);
                md = std::max(md, d);
                if (d <= tau) ++nin;
                const double v = (pt - S.C).Dot(S.A);
                if (v < vlo - tau || v > vhi + tau) inside = false;
                if (!inArc(azimuthOf(pt, S, xd), u0, u1, tau / S.R)) inside = false;
            }
            if (!inside) continue;
            if (nin == 3 && triLevels(mv, t, S, tau) >= kMinAxialLevels) I.push_back(t);
            else if (nin >= 1 && md <= s) Cv.push_back(t);
        }
        std::vector<int> C = g.tris;
        C.insert(C.end(), I.begin(), I.end());
        C.insert(C.end(), Cv.begin(), Cv.end());
        std::sort(C.begin(), C.end());
        C.erase(std::unique(C.begin(), C.end()), C.end());
        std::vector<int> cert = g.tris;
        cert.insert(cert.end(), I.begin(), I.end());
        std::sort(cert.begin(), cert.end());
        gp_Dir fa;
        gp_Pnt fc;
        double fr = 0, fres = 0;
        const bool fitOk = cylinderFitLS(mv, cert, gp_Dir(S.A), fa, fc, fr, fres, maxIt);
        const auto cc = edgeComponents(rm, C);
        const int loops = cc.size() == 1 ? boundaryLoops(rm, C) : -1;

        BandDiag dg;
        dg.mode = "closure";
        dg.seed = minTriOf(g.tris);
        dg.members = (int)g.members.size();
        dg.nTri = (int)C.size();
        dg.nV = (int)uniqueVerts(rm, C).size();
        dg.nI = (int)I.size();
        dg.nCover = (int)Cv.size();
        dg.comps = (int)cc.size();
        dg.loops = loops;
        dg.sTau = s / tau;
        dg.radius = fitOk ? fr : S.R;
        dg.axis = fitOk ? fa.XYZ() : S.A;
        dg.certDevTau = fres / tau;
        int cols = 0, levels = 0;
        stationsOf(rm, cert, fitOk ? Cyl{fc.XYZ(), fa.XYZ(), fr} : S, tau, cols, levels);
        dg.stations = cols;
        dg.levels = levels;
        if (!fitOk || fres > tau) {
            dg.why = "vertex-off";
            printBand(dg);
            ++deferred;
            continue;
        }
        // A result whose pieces the fill does not join stays separate faces (D-130-14);
        // a result with an inner loop is the Seamed360/inner-wire build (owner seam).
        if (cc.size() != 1) {
            dg.why = "not-joined";
            printBand(dg);
            ++deferred;
            continue;
        }
        if (loops != 1) {
            dg.why = "domain-not-disk";
            printBand(dg);
            ++deferred;
            continue;
        }
        const Cyl N{fc.XYZ(), fa.XYZ(), fr};
        CylFields f;
        if (!cylFields(rm, tol, C, N, f)) {
            dg.why = "degenerate";
            printBand(dg);
            ++deferred;
            continue;
        }
        Region nr = R0;
        nr.tris = C;
        applyCylFields(nr, N, f);
        std::vector<char> inC((std::size_t)nTri, 0);
        for (int t : C) {
            inC[(std::size_t)t] = 1;
            unclaimed[(std::size_t)t] = 0;
        }
        dg.straddle = trimClaimed(rm, work, inC);
        dg.maxDev = nr.maxVertexDev;
        dg.maxDevTau = nr.maxVertexDev / tau;
        dg.nSides = nr.nSides;
        dg.dVol = nr.dVolPredicted;
        if (diagOn())
            for (int m : g.members)
                std::fprintf(stderr, "DIAG_REL_CLOSURE_MEMBER seed=%d member=%d nTri=%zu nSides=%d dVol=%.6g\n",
                             dg.seed, minTriOf(work.accepted[(std::size_t)m].tris),
                             work.accepted[(std::size_t)m].tris.size(),
                             work.accepted[(std::size_t)m].nSides,
                             work.accepted[(std::size_t)m].dVolPredicted);
        work.accepted[(std::size_t)g.members.front()] = std::move(nr);
        for (std::size_t k = 1; k < g.members.size(); ++k) erase.push_back(g.members[k]);
        ++committed;
        dg.why = "accept";
        printBand(dg);
    }
    eraseRegions(work, erase);
    return committed;
}

bool reclaimLine(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                 SegmentWork& work) {
    const double tau = tauOf(mv);
    if (!(tau > 0.0) || mv.nTri == 0 || !mv.triEdges) return true;
    const RMesh rm = makeRMesh(mv);
    const int healed = healRefusedDonors(rm, p, tol, tau, work);
    int deferred = 0;
    const int closed = closeSameSurface(rm, p, tol, tau, work, deferred);
    if (diagOn())
        std::fprintf(stderr, "DIAG_REL_SUM stage=reclaim mode=line healed=%d closed=%d deferred=%d\n",
                     healed, closed, deferred);
    return true;
}

}  // namespace

bool reclaimRelational(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                       SegmentWork& work, RelRevolveMode mode) {
    if (mode == RelRevolveMode::Line) return reclaimLine(mv, p, tol, work);
    return true;  // Torus mode (H1): the tori lane.
}

}}  // namespace stl2step::refit
