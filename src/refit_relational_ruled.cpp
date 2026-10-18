// RULED QUAD between certified rails (D-fc-5.9). Stage claimRuledQuads runs at
// the tail of H0: after reclaimRelational(Line), before the pre-T snapshot, so
// it sees only what B1/C1 and the relational stages left, and U-3 never
// reverts it.
//
// Certificate order (D-fc-1.6; no constant but tau = 2q and parameter counts):
//   1. domain = one A1 chart whose triangles are all in Unclaimed provisionals;
//      its boundary is one loop (a disk) in four runs, a run being a maximal
//      stretch against one neighbouring chart;
//   2. caps: an opposite pair of runs, each lying in ONE plane anchor
//      (anchorCertified, a3WouldCommit); no such pair => no-anchor;
//   3. every run line-first certified (>= minStations(kParamsLine) stations,
//      every vertex within tau of its LS line);
//   4. corners = rail line ∩ cap plane -- zero free parameters;
//   5. line first on the surface: corners coplanar within tau => line (A3);
//   6. evidence: >= minStations(kParamsRuled) interior vertices;
//   7. EVERY domain vertex within tau of the face (max, not RMS);
//   8. same surface: an adjacent chart all of whose vertices are within tau
//      of S joins the domain, and the certificate is re-run on the union
//      (D-130-14); a failing union is refused whole.
// The region takes whole charts. A2 grows provisionals inside a chart, so a
// chart is a union of whole provisionals and no provisional straddles it.
// SPDX-License-Identifier: MIT

#include "refit_relational_ruled.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {

namespace {

bool diagRelOn() {
    const char* v = std::getenv("STL2STEP_DIAG_REL");
    return v && v[0] && v[0] != '0';
}

gp_XYZ vtxPt(const MeshView& mv, int lv) { return mv.pts[mv.compVtx[lv]]; }

// Local vertex ids of triangle t in its winding order. Side s runs
// compEdges[e].first -> .second when bit s of triDirs is set.
std::array<int, 3> triLocalVerts(const MeshView& mv, int t) {
    std::array<int, 3> out{};
    for (int s = 0; s < 3; ++s) {
        const int e = mv.triEdges[t][s];
        const bool fwd = (mv.triDirs[t] >> s) & 1;
        out[(std::size_t)s] = fwd ? mv.compEdges[e].first : mv.compEdges[e].second;
    }
    return out;
}

struct Lsq3 {
    // Covariance of a point set about its centroid; principal axes ascending.
    gp_XYZ c;
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    bool ok = false;
};

Lsq3 lsq3(const std::vector<gp_XYZ>& pts) {
    Lsq3 L;
    if (pts.empty()) return L;
    gp_XYZ c(0, 0, 0);
    for (const gp_XYZ& p : pts) c += p;
    c /= (double)pts.size();
    std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (const gp_XYZ& p : pts) {
        const double d[3] = {p.X() - c.X(), p.Y() - c.Y(), p.Z() - c.Z()};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) cov[(std::size_t)i][(std::size_t)j] += d[i] * d[j];
    }
    L.c = c;
    L.ok = jacobiEigenSymmetric3(cov, L.eval, L.evec);
    return L;
}

}  // namespace

ChainLineCert certifyChainLine(const std::vector<gp_XYZ>& pts, double tau) {
    ChainLineCert C;
    const Lsq3 L = lsq3(pts);
    if (!L.ok) return C;
    gp_XYZ d(L.evec[2][0], L.evec[2][1], L.evec[2][2]);
    const double dm = d.Modulus();
    if (!(dm > 0)) return C;
    d /= dm;
    C.p = L.c;
    C.d = d;
    std::vector<double> s;
    s.reserve(pts.size());
    for (const gp_XYZ& p : pts) {
        const gp_XYZ w = p - L.c;
        const double t = w.Dot(d);
        s.push_back(t);
        C.maxDev = std::max(C.maxDev, (w - d * t).Modulus());
    }
    std::sort(s.begin(), s.end());
    for (std::size_t i = 0; i < s.size(); ++i)
        if (i == 0 || s[i] - s[i - 1] > tau) ++C.stations;
    if (C.stations < minStations(kParamsLine)) {
        C.why = "underdetermined";
        return C;
    }
    if (!(C.maxDev <= tau)) {
        C.why = "vertex-off";
        return C;
    }
    C.ok = true;
    C.why = "accept";
    return C;
}

bool ruledCapAnchor(bool unclaimed, bool a3Commits, double maxVertexDev, double tau) {
    return unclaimed && a3Commits && anchorCertified(tau, maxVertexDev);
}

bool railCapCorner(const gp_XYZ& railP, const gp_XYZ& railD, const gp_XYZ& planeP,
                   const gp_XYZ& planeN, gp_XYZ& out) {
    const double dn = railD.Dot(planeN);
    if (dn == 0.0) return false;
    const double t = (planeP - railP).Dot(planeN) / dn;
    if (!std::isfinite(t)) return false;
    out = railP + railD * t;
    return true;
}

double ruledCornerPlanarity(const RuledPoles& P) {
    const Lsq3 L = lsq3(std::vector<gp_XYZ>(P.begin(), P.end()));
    if (!L.ok) return 0.0;
    const gp_XYZ n(L.evec[0][0], L.evec[0][1], L.evec[0][2]);
    const double nm = n.Modulus();
    if (!(nm > 0)) return 0.0;
    double m = 0;
    for (const gp_XYZ& p : P) m = std::max(m, std::fabs((p - L.c).Dot(n) / nm));
    return m;
}

gp_XYZ ruledPoint(const RuledPoles& P, double u, double v) {
    return P[0] * ((1 - u) * (1 - v)) + P[1] * (u * (1 - v)) + P[2] * ((1 - u) * v) +
           P[3] * (u * v);
}

gp_XYZ ruledNormal(const RuledPoles& P, double u, double v) {
    const gp_XYZ Su = (P[1] - P[0]) * (1 - v) + (P[3] - P[2]) * v;
    const gp_XYZ Sv = (P[2] - P[0]) * (1 - u) + (P[3] - P[1]) * u;
    return Su.Crossed(Sv);
}

namespace {

// The start of the closest-point iteration: the closed-form inverse of the
// bilinear map in the plane normal to S at the patch centre (the quadratic in
// v of (e - v b) x (a + v c) = 0), the root nearest [0, 1].
void ruledInverseStart(const RuledPoles& P, const gp_XYZ& x, double& u, double& v) {
    const gp_XYZ a = P[1] - P[0], b = P[2] - P[0], c = P[3] - P[1] - P[2] + P[0];
    u = 0.5;
    v = 0.5;
    gp_XYZ n = ruledNormal(P, 0.5, 0.5);
    const double nm = n.Modulus();
    const double am = a.Modulus();
    if (!(nm > 0) || !(am > 0)) return;
    n /= nm;
    gp_XYZ t1 = a / am;
    t1 -= n * t1.Dot(n);
    const double t1m = t1.Modulus();
    if (!(t1m > 0)) return;
    t1 /= t1m;
    const gp_XYZ t2 = n.Crossed(t1);
    auto cr = [&](const gp_XYZ& p, const gp_XYZ& q) {
        return p.Dot(t1) * q.Dot(t2) - p.Dot(t2) * q.Dot(t1);
    };
    const gp_XYZ e = x - P[0];
    const double A = -cr(b, c), B = cr(e, c) - cr(b, a), C = cr(e, a);
    auto dist01 = [](double t) { return t < 0 ? -t : (t > 1 ? t - 1 : 0.0); };
    double vr = 0.5;
    if (A == 0.0) {
        if (B != 0.0) vr = -C / B;
    } else {
        const double disc = B * B - 4 * A * C;
        if (disc >= 0) {
            const double q = -0.5 * (B + std::copysign(std::sqrt(disc), B));
            const double r1 = q / A;
            const double r2 = q != 0.0 ? C / q : r1;
            vr = dist01(r1) <= dist01(r2) ? r1 : r2;
        } else {
            vr = -B / (2 * A);
        }
    }
    if (!std::isfinite(vr)) return;
    const gp_XYZ g = a + c * vr;
    const double gg = g.Dot(g);
    if (!(gg > 0)) return;
    u = (e - b * vr).Dot(g) / gg;
    v = vr;
    if (!std::isfinite(u)) {
        u = 0.5;
        v = 0.5;
    }
}

}  // namespace

double ruledDistance(const RuledPoles& P, const gp_XYZ& x, int maxIters, bool clampDomain,
                     double* uOut, double* vOut, bool* converged) {
    const gp_XYZ a = P[1] - P[0], b = P[2] - P[0], c = P[3] - P[1] - P[2] + P[0];
    double u = 0.5, v = 0.5;
    ruledInverseStart(P, x, u, v);
    // Gauss-Newton on |S(u, v) - x|^2: J^T J is positive semi-definite, so a
    // step never heads for the saddle the full Hessian can point at.
    // Converged means a critical point: the gradient is the zero vector, or
    // the step does not change (u, v) in double. Both are exact. A singular
    // Jacobian or an exhausted iteration cap is non-convergence (D-fc-5.18 F-2).
    bool ok = false;
    for (int it = 0; it < std::max(1, maxIters); ++it) {
        const gp_XYZ S = P[0] + a * u + b * v + c * (u * v);
        const gp_XYZ Su = a + c * v, Sv = b + c * u;
        const gp_XYZ r = S - x;
        const double fu = r.Dot(Su), fv = r.Dot(Sv);
        if (fu == 0.0 && fv == 0.0) {
            ok = true;
            break;
        }
        const double juu = Su.Dot(Su), jvv = Sv.Dot(Sv), juv = Su.Dot(Sv);
        const double det = juu * jvv - juv * juv;
        if (!(det > 0) || !std::isfinite(det)) break;
        const double du = (jvv * fu - juv * fv) / det;
        const double dv = (juu * fv - juv * fu) / det;
        if (!std::isfinite(du) || !std::isfinite(dv)) break;
        if (u - du == u && v - dv == v) {
            ok = true;
            break;
        }
        u -= du;
        v -= dv;
    }
    if (clampDomain) {
        u = std::min(1.0, std::max(0.0, u));
        v = std::min(1.0, std::max(0.0, v));
    }
    if (uOut) *uOut = u;
    if (vOut) *vOut = v;
    if (converged) *converged = ok;
    return (ruledPoint(P, u, v) - x).Modulus();
}

RuledVertexCert certifyRuledVertices(const RuledPoles& P, const std::vector<gp_XYZ>& pts,
                                     double tau, int maxIters) {
    RuledVertexCert C;
    double ss = 0;
    for (const gp_XYZ& x : pts) {
        const double d = ruledDistance(P, x, maxIters, true);
        C.maxDev = std::max(C.maxDev, d);
        ss += d * d;
    }
    C.rmsDev = pts.empty() ? 0.0 : std::sqrt(ss / (double)pts.size());
    C.ok = !pts.empty() && C.maxDev <= tau;
    return C;
}

bool ruledSameSurface(const RuledPoles& P, const std::vector<gp_XYZ>& pts, double tau,
                      int maxIters) {
    if (pts.empty()) return false;
    for (const gp_XYZ& x : pts) {
        bool conv = false;
        const double d = ruledDistance(P, x, maxIters, false, nullptr, nullptr, &conv);
        // D-fc-5.18 F-2: non-convergence refuses, even when the upper bound
        // happens to sit within tau.
        if (!ruledDistanceCertifies(d, conv, tau)) return false;
    }
    return true;
}

double dVolRuledTriangle(const RuledPoles& P, const gp_XYZ& a, const gp_XYZ& b,
                         const gp_XYZ& c, int maxIters, double* sagOut) {
    gp_XYZ n = (b - a).Crossed(c - a);
    const double twiceA = n.Modulus();
    if (sagOut) *sagOut = 0;
    if (!(twiceA > 0)) return 0.0;
    n /= twiceA;
    // 7-point rule, exact for cubics: centroid 27/60, vertices 3/60 each,
    // edge midpoints 8/60 each.
    const gp_XYZ q[7] = {(a + b + c) / 3.0, a, b, c, (a + b) * 0.5, (b + c) * 0.5,
                         (c + a) * 0.5};
    const double w[7] = {27.0 / 60, 3.0 / 60, 3.0 / 60, 3.0 / 60, 8.0 / 60, 8.0 / 60, 8.0 / 60};
    double acc = 0, sag = 0;
    for (int k = 0; k < 7; ++k) {
        double u = 0, v = 0;
        (void)ruledDistance(P, q[k], maxIters, false, &u, &v);
        const gp_XYZ d = ruledPoint(P, u, v) - q[k];
        acc += w[k] * d.Dot(n);
        sag = std::max(sag, d.Modulus());
    }
    if (sagOut) *sagOut = sag;
    return 0.5 * twiceA * acc;
}

namespace {

struct BEdge {
    int from = -1, to = -1;
    int acrossTri = -1, acrossChart = -1;
};

struct QuadTry {
    const char* why = "not-quad";
    RuledPoles P{};
    std::array<ChainLineCert, 4> line{};
    std::array<int, 4> nbrChart{{-1, -1, -1, -1}};
    int capsPair = -1;  // 0: runs (0,2) are the caps; 1: runs (1,3)
    const char* failAt = "-";  // rail | cap | surface: where the certificate stopped
    double planarity = 0;
    int interior = 0;
    RuledVertexCert vc;
    std::vector<int> verts;  // domain vertices, ascending
};

}  // namespace

bool claimRuledQuads(const MeshView& mv, const SegmentParams& p, const DerivedTols& tol,
                     SegmentWork& work) {
    const double tau = tauOf(mv);
    const int nTri = (int)mv.nTri;
    if (!(tau > 0) || nTri == 0 || work.nCharts <= 0 || (int)work.triChart.size() != nTri)
        return true;
    const bool diag = diagRelOn();
    const auto t0 = std::chrono::steady_clock::now();
    const int maxIters = p.maxRefineIters;

    std::vector<int> triProv((std::size_t)nTri, -1);
    for (int pi = 0; pi < (int)work.provisionals.size(); ++pi)
        for (int t : work.provisionals[(std::size_t)pi].tris)
            if (t >= 0 && t < nTri) triProv[(std::size_t)t] = pi;
    // A3's own exemption for fillet flanks (a3WouldCommit).
    std::vector<char> flank(work.provisionals.size(), 0);
    for (const Region& r : work.accepted) {
        if (r.origin != Origin::FilletStrip) continue;
        if (r.filletNbrA >= 0 && r.filletNbrA < (int)flank.size()) flank[(std::size_t)r.filletNbrA] = 1;
        if (r.filletNbrB >= 0 && r.filletNbrB < (int)flank.size()) flank[(std::size_t)r.filletNbrB] = 1;
    }
    std::vector<std::array<int, 2>> edgeTri(mv.nEdge, {{-1, -1}});
    for (int t = 0; t < nTri; ++t)
        for (int s = 0; s < 3; ++s) {
            const int e = mv.triEdges[t][s];
            auto& et = edgeTri[(std::size_t)e];
            if (et[0] < 0) et[0] = t;
            else if (et[1] < 0) et[1] = t;
        }
    std::vector<std::vector<int>> chartTris((std::size_t)work.nCharts);
    for (int t = 0; t < nTri; ++t) {
        const int c = work.triChart[(std::size_t)t];
        if (c >= 0 && c < work.nCharts) chartTris[(std::size_t)c].push_back(t);
    }
    // A candidate chart is wholly unclaimed: no earlier stage owns any of it.
    std::vector<char> eligible((std::size_t)work.nCharts, 0);
    for (int c = 0; c < work.nCharts; ++c) {
        bool ok = !chartTris[(std::size_t)c].empty();
        for (int t : chartTris[(std::size_t)c]) {
            const int pi = triProv[(std::size_t)t];
            if (pi < 0 || work.provisionals[(std::size_t)pi].claim != ProvClaim::Unclaimed) {
                ok = false;
                break;
            }
        }
        eligible[(std::size_t)c] = ok ? 1 : 0;
    }

    std::vector<char> inDom((std::size_t)nTri, 0);
    // The certificate of one domain (a set of whole charts).
    auto tryQuad = [&](const std::vector<int>& dom) -> QuadTry {
        QuadTry Q;
        for (int t : dom) inDom[(std::size_t)t] = 1;
        struct Clear {
            std::vector<char>& m;
            const std::vector<int>& d;
            ~Clear() {
                for (int t : d) m[(std::size_t)t] = 0;
            }
        } clear{inDom, dom};
        // Boundary edges in the triangles' winding (the loop keeps the domain on
        // its left), and the disk test V - E + F = 1.
        std::vector<BEdge> be;
        std::vector<int> vs, es;
        for (int t : dom) {
            const std::array<int, 3> lv = triLocalVerts(mv, t);
            for (int s = 0; s < 3; ++s) {
                vs.push_back(lv[(std::size_t)s]);
                const int e = mv.triEdges[t][s];
                es.push_back(e);
                const auto& et = edgeTri[(std::size_t)e];
                const int o = et[0] == t ? et[1] : et[0];
                if (o < 0) return Q;  // open edge: not a clean disk
                if (inDom[(std::size_t)o]) continue;
                BEdge b;
                b.from = lv[(std::size_t)s];
                b.to = lv[(std::size_t)((s + 1) % 3)];
                b.acrossTri = o;
                b.acrossChart = work.triChart[(std::size_t)o];
                be.push_back(b);
            }
        }
        std::sort(vs.begin(), vs.end());
        vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
        std::sort(es.begin(), es.end());
        es.erase(std::unique(es.begin(), es.end()), es.end());
        if ((long)vs.size() - (long)es.size() + (long)dom.size() != 1 || be.size() < 4) return Q;
        // One loop through every boundary edge; a pinch vertex is not a disk.
        std::vector<std::pair<int, int>> out;  // (from, edge index)
        for (int i = 0; i < (int)be.size(); ++i) out.emplace_back(be[(std::size_t)i].from, i);
        std::sort(out.begin(), out.end());
        for (std::size_t i = 1; i < out.size(); ++i)
            if (out[i].first == out[i - 1].first) return Q;
        auto next = [&](int v) -> int {
            auto it = std::lower_bound(out.begin(), out.end(), std::make_pair(v, -1));
            return (it != out.end() && it->first == v) ? it->second : -1;
        };
        std::vector<int> loop;
        int cur = out.front().second;
        for (std::size_t k = 0; k < be.size(); ++k) {
            loop.push_back(cur);
            cur = next(be[(std::size_t)cur].to);
            if (cur < 0) return Q;
            if (cur == loop.front()) break;
        }
        if (loop.size() != be.size() || cur != loop.front()) return Q;
        // Runs by neighbouring chart.
        const int L = (int)loop.size();
        int start = -1;
        for (int k = 0; k < L; ++k) {
            const int a = be[(std::size_t)loop[(std::size_t)((k + L - 1) % L)]].acrossChart;
            const int b = be[(std::size_t)loop[(std::size_t)k]].acrossChart;
            if (a != b) {
                start = k;
                break;
            }
        }
        if (start < 0) return Q;
        std::vector<std::vector<int>> runs;  // edge indices
        for (int k = 0; k < L; ++k) {
            const int ei = loop[(std::size_t)((start + k) % L)];
            if (runs.empty() || be[(std::size_t)runs.back().front()].acrossChart !=
                                    be[(std::size_t)ei].acrossChart)
                runs.emplace_back();
            runs.back().push_back(ei);
        }
        if (runs.size() != 4) return Q;
        // Canonical first run: the one starting at the lowest local vertex id.
        int r0 = 0;
        for (int k = 1; k < 4; ++k)
            if (be[(std::size_t)runs[(std::size_t)k].front()].from <
                be[(std::size_t)runs[(std::size_t)r0].front()].from)
                r0 = k;
        std::rotate(runs.begin(), runs.begin() + r0, runs.end());
        // Caps: one plane anchor holds every triangle across the run.
        std::array<int, 4> capProv{{-1, -1, -1, -1}};
        for (int k = 0; k < 4; ++k) {
            Q.nbrChart[(std::size_t)k] = be[(std::size_t)runs[(std::size_t)k].front()].acrossChart;
            int pi = -2;
            for (int ei : runs[(std::size_t)k]) {
                const int q = triProv[(std::size_t)be[(std::size_t)ei].acrossTri];
                if (pi == -2) pi = q;
                else if (pi != q) pi = -1;
            }
            if (pi < 0) continue;
            const Provisional& pr = work.provisionals[(std::size_t)pi];
            if (ruledCapAnchor(pr.claim == ProvClaim::Unclaimed,
                               a3WouldCommit(mv, tol, pr, flank[(std::size_t)pi] != 0),
                               pr.maxVertexDev, tau))
                capProv[(std::size_t)k] = pi;
        }
        if (capProv[0] >= 0 && capProv[2] >= 0) Q.capsPair = 0;
        else if (capProv[1] >= 0 && capProv[3] >= 0) Q.capsPair = 1;
        else {
            Q.why = "no-anchor";
            return Q;
        }
        // Line first on every boundary chain, points ascending by vertex id.
        for (int k = 0; k < 4; ++k) {
            std::vector<int> rv;
            for (int ei : runs[(std::size_t)k]) {
                rv.push_back(be[(std::size_t)ei].from);
                rv.push_back(be[(std::size_t)ei].to);
            }
            std::sort(rv.begin(), rv.end());
            rv.erase(std::unique(rv.begin(), rv.end()), rv.end());
            std::vector<gp_XYZ> pts;
            for (int v : rv) pts.push_back(vtxPt(mv, v));
            Q.line[(std::size_t)k] = certifyChainLine(pts, tau);
        }
        for (int k = 0; k < 4; ++k)
            if (!Q.line[(std::size_t)k].ok) {
                Q.why = Q.line[(std::size_t)k].why;
                Q.failAt = (k % 2 == Q.capsPair) ? "cap" : "rail";
                return Q;
            }
        // Walk order: cap0 = run c0 runs from railA to railB; cap1 = run c0 + 2.
        const int c0 = Q.capsPair, rB = c0 + 1, c1 = c0 + 2, rA = (c0 + 3) % 4;
        const Provisional& cap0 = work.provisionals[(std::size_t)capProv[(std::size_t)c0]];
        const Provisional& cap1 = work.provisionals[(std::size_t)capProv[(std::size_t)c1]];
        const gp_XYZ n0 = cap0.plane.Direction().XYZ(), n1 = cap1.plane.Direction().XYZ();
        const gp_XYZ o0 = cap0.plane.Location().XYZ(), o1 = cap1.plane.Location().XYZ();
        const ChainLineCert& la = Q.line[(std::size_t)rA];
        const ChainLineCert& lb = Q.line[(std::size_t)rB];
        if (!railCapCorner(la.p, la.d, o0, n0, Q.P[0]) || !railCapCorner(lb.p, lb.d, o0, n0, Q.P[1]) ||
            !railCapCorner(la.p, la.d, o1, n1, Q.P[2]) || !railCapCorner(lb.p, lb.d, o1, n1, Q.P[3])) {
            Q.why = "degenerate";
            return Q;
        }
        Q.failAt = "surface";
        Q.planarity = ruledCornerPlanarity(Q.P);
        if (Q.planarity <= tau) {
            Q.why = "line";
            return Q;
        }
        {
            std::vector<int> bv;
            for (const BEdge& b : be) bv.push_back(b.from);
            std::sort(bv.begin(), bv.end());
            for (int v : vs)
                if (!std::binary_search(bv.begin(), bv.end(), v)) ++Q.interior;
        }
        if (Q.interior < minStations(kParamsRuled)) {
            Q.why = "underdetermined";
            return Q;
        }
        Q.verts = vs;
        std::vector<gp_XYZ> pts;
        pts.reserve(vs.size());
        for (int v : vs) pts.push_back(vtxPt(mv, v));
        Q.vc = certifyRuledVertices(Q.P, pts, tau, maxIters);
        Q.why = Q.vc.ok ? "accept" : "vertex-off";
        if (Q.vc.ok) Q.failAt = "-";
        return Q;
    };

    std::vector<char> done((std::size_t)work.nCharts, 0);
    int nQuad = 0, nAccept = 0, nMerged = 0;
    int nLine = 0, nOff = 0, nAnchor = 0, nUnder = 0, nSame = 0, nDegen = 0;
    for (int c = 0; c < work.nCharts; ++c) {
        if (!eligible[(std::size_t)c] || done[(std::size_t)c]) continue;
        std::vector<int> charts{c};
        std::vector<int> dom = chartTris[(std::size_t)c];
        QuadTry Q = tryQuad(dom);
        if (std::string(Q.why) == "not-quad") continue;
        ++nQuad;
        // Same surface: an adjacent chart on S joins; the union is certified
        // again, from its own rails and caps. A union that does not certify is
        // refused whole -- the chart alone is a passing subset of it.
        const char* why = Q.why;
        for (int round = 0; std::string(why) == "accept" && round < work.nCharts; ++round) {
            std::vector<int> join;
            bool foreign = false;
            for (int k = 0; k < 4; ++k) {
                const int d = Q.nbrChart[(std::size_t)k];
                if (d < 0 || std::find(charts.begin(), charts.end(), d) != charts.end() ||
                    std::find(join.begin(), join.end(), d) != join.end())
                    continue;
                std::vector<int> dv;
                for (int t : chartTris[(std::size_t)d])
                    for (int v : triLocalVerts(mv, t)) dv.push_back(v);
                std::sort(dv.begin(), dv.end());
                dv.erase(std::unique(dv.begin(), dv.end()), dv.end());
                std::vector<gp_XYZ> pts;
                for (int v : dv) pts.push_back(vtxPt(mv, v));
                if (!ruledSameSurface(Q.P, pts, tau, maxIters)) continue;
                if (!eligible[(std::size_t)d] || done[(std::size_t)d]) foreign = true;
                else join.push_back(d);
            }
            if (foreign) {
                why = "same-surface";
                break;
            }
            if (join.empty()) break;
            std::sort(join.begin(), join.end());
            for (int d : join) {
                charts.push_back(d);
                dom.insert(dom.end(), chartTris[(std::size_t)d].begin(),
                           chartTris[(std::size_t)d].end());
            }
            std::sort(dom.begin(), dom.end());
            nMerged += (int)join.size();
            Q = tryQuad(dom);
            why = Q.why;
            if (std::string(why) == "not-quad") why = "same-surface";
        }
        for (int d : charts) done[(std::size_t)d] = 1;
        const std::string w(why);
        if (w == "line") ++nLine;
        else if (w == "vertex-off") ++nOff;
        else if (w == "no-anchor") ++nAnchor;
        else if (w == "underdetermined") ++nUnder;
        else if (w == "same-surface") ++nSame;
        else if (w == "degenerate") ++nDegen;

        Region r;
        double dVol = 0, sag = 0, sigma = 0, area = 0;
        if (w == "accept") {
            for (int t : dom) {
                const std::array<int, 3> lv = triLocalVerts(mv, t);
                const gp_XYZ a = vtxPt(mv, lv[0]), b = vtxPt(mv, lv[1]), cc = vtxPt(mv, lv[2]);
                double s = 0;
                dVol += dVolRuledTriangle(Q.P, a, b, cc, maxIters, &s);
                sag = std::max(sag, s);
                const gp_XYZ nt = (b - a).Crossed(cc - a);
                double u = 0, v = 0;
                (void)ruledDistance(Q.P, (a + b + cc) / 3.0, maxIters, true, &u, &v);
                const gp_XYZ ns = ruledNormal(Q.P, u, v);
                const double nm = ns.Modulus();
                if (nm > 0) sigma += nt.Dot(ns) / nm;
                area += 0.5 * nt.Modulus();
            }
        }
        if (diag) {
            const ChainLineCert* L = Q.line.data();
            double lineRes = 0;
            for (int k = 0; k < 4; ++k) lineRes = std::max(lineRes, L[k].maxDev);
            std::fprintf(stderr,
                         "DIAG_REL_BAND stage=ruled mode=quad seed=%d charts=%zu nTri=%zu nV=%zu "
                         "stations=%d/%d/%d/%d caps=%d params=%d fit=anchored interior=%d "
                         "planarity=%.6g planarityTau=%.4g lineResTau=%.4f maxDev=%.6g "
                         "maxDevTau=%.4f rmsTau=%.4f dVol=%.6f sag=%.6g area=%.6f at=%s why=%s\n",
                         dom.empty() ? -1 : dom.front(), charts.size(), dom.size(),
                         Q.verts.size(), L[0].stations, L[1].stations, L[2].stations,
                         L[3].stations, Q.capsPair, kParamsRuled, Q.interior, Q.planarity,
                         Q.planarity / tau, lineRes / tau, Q.vc.maxDev, Q.vc.maxDev / tau,
                         Q.vc.rmsDev / tau, dVol, sag, area, Q.failAt, why);
        }
        if (w != "accept") continue;
        ++nAccept;
        const gp_XYZ nMid = ruledNormal(Q.P, 0.5, 0.5);
        gp_XYZ xd = Q.P[1] - Q.P[0];
        const double nmm = nMid.Modulus();
        if (!(nmm > 0) || !(xd.Modulus() > 0)) continue;
        const gp_XYZ zn = nMid / nmm;
        xd -= zn * xd.Dot(zn);
        if (!(xd.Modulus() > 0)) continue;
        r.id = 0;
        r.type = SurfType::BSpline;
        r.origin = Origin::RuledQuad;
        r.ax = gp_Ax3(gp_Pnt(Q.P[0]), gp_Dir(zn), gp_Dir(xd));
        r.ruledPoles = Q.P;
        r.uMin = 0;
        r.uMax = 1;
        r.vMin = 0;
        r.vMax = 1;
        r.closed360 = false;
        // The material side from the mesh: the outward triangle normals against
        // S_u x S_v (the fit decides, never both sides; D-fc-1.6.8).
        r.outwardNormal = sigma > 0;
        r.tris = dom;
        r.maxVertexDev = Q.vc.maxDev;
        r.rmsVertexDev = Q.vc.rmsDev;
        r.nSides = 0;
        r.chordSagitta = sag;
        r.dVolPredicted = dVol;
        r.reject = Reject::None;
        r.builtAs = BuiltAs::NotBuilt;
        work.accepted.push_back(r);
        for (int t : dom) {
            const int pi = triProv[(std::size_t)t];
            if (pi >= 0) work.provisionals[(std::size_t)pi].claim = ProvClaim::ConsumedRuled;
        }
    }
    if (diag)
        std::fprintf(stderr,
                     "DIAG_REL_SUM stage=ruled charts=%d quads=%d accept=%d merged=%d line=%d "
                     "vertex-off=%d no-anchor=%d underdetermined=%d same-surface=%d "
                     "degenerate=%d seconds=%.6f\n",
                     work.nCharts, nQuad, nAccept, nMerged, nLine, nOff, nAnchor, nUnder, nSame,
                     nDegen,
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return true;
}

}  // namespace refit
}  // namespace stl2step
