// Tessellation-law refutation of NgonWall cylinders (D-fc-1.9).
// A monotone chord/angle law never gives a larger circle a coarser step than
// a smaller one. An NgonWall whose step violates that bound is removed and
// its walls return to Unclaimed so A3 commits them as planes.
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

bool tesslawDiagOn() {
    const char* e = std::getenv("STL2STEP_DIAG_TESSLAW");
    return e && e[0] != '\0' && e[0] != '0';
}

bool meshReadable(const MeshView& mv) {
    return mv.pts && mv.tris && mv.compTris && mv.nTri > 0;
}

gp_XYZ triVert(const MeshView& mv, int lt, int corner) {
    const int gt = mv.compTris[lt];
    return mv.pts[mv.tris[gt][corner]];
}

double triArea(const MeshView& mv, int lt) {
    const gp_XYZ a = triVert(mv, lt, 0);
    const gp_XYZ b = triVert(mv, lt, 1);
    const gp_XYZ c = triVert(mv, lt, 2);
    return 0.5 * (b - a).Crossed(c - a).Modulus();
}

double pointPlaneDist(const gp_Ax3& plane, const gp_XYZ& p) {
    return std::abs((p - plane.Location().XYZ()).Dot(plane.Direction().XYZ()));
}

double triPlaneDist(const MeshView& mv, int lt, const gp_Ax3& plane) {
    double m = 0.0;
    for (int k = 0; k < 3; ++k)
        m = std::max(m, pointPlaneDist(plane, triVert(mv, lt, k)));
    return m;
}

int minTriOf(const std::vector<int>& tris) {
    int m = 0;
    bool any = false;
    for (int t : tris) {
        if (!any || t < m) {
            m = t;
            any = true;
        }
    }
    return any ? m : -1;
}

bool isLawCylinder(const Region& r) {
    if (r.type != SurfType::Cylinder || r.tris.empty()) return false;
    return r.origin == Origin::CylGrow || r.origin == Origin::FilletStrip
        || r.origin == Origin::NgonWall;
}

// Azimuth stations of the region's own vertices about the refit axis.
// Vertices within τ of the axis have no azimuth. Stations merge at angle τ/r
// (arc length τ). The walk starts after the largest gap so a cluster that
// crosses angle 0 stays one station. One station (every gap ≤ τ/r) returns
// a single angle.
bool azimuthStations(const MeshView& mv, const std::vector<int>& tris, const gp_Dir& axis,
                     const gp_Pnt& center, double radius, double tau,
                     std::vector<double>& stationAng) {
    stationAng.clear();
    if (!(radius > 0.0) || !(tau > 0.0)) return false;
    const gp_XYZ W(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ C = center.XYZ();
    const gp_XYZ tmp = (std::abs(W.X()) < 0.9) ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 1, 0);
    gp_XYZ U = W.Crossed(tmp);
    const double um = U.Modulus();
    if (!(um > 0.0)) return false;
    U.Divide(um);
    gp_XYZ V = W.Crossed(U);
    const double vm = V.Modulus();
    if (!(vm > 0.0)) return false;
    V.Divide(vm);

    std::vector<int> gids;
    for (int t : tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        const int gt = mv.compTris[t];
        if (gt < 0) continue;
        for (int k = 0; k < 3; ++k) gids.push_back(mv.tris[gt][k]);
    }
    std::sort(gids.begin(), gids.end());
    gids.erase(std::unique(gids.begin(), gids.end()), gids.end());

    std::vector<double> angs;
    angs.reserve(gids.size());
    for (int g : gids) {
        if (g < 0) continue;
        const gp_XYZ rad = (mv.pts[g] - C);
        const gp_XYZ perp = rad - W * rad.Dot(W);
        if (!(perp.Modulus() > tau)) continue;
        double a = std::atan2(perp.Dot(V), perp.Dot(U));
        if (a < 0.0) a += kTwoPi;
        if (a >= kTwoPi) a = 0.0;
        angs.push_back(a);
    }
    if (angs.empty()) return false;
    std::sort(angs.begin(), angs.end());

    const double delta = tau / radius;
    const int n = static_cast<int>(angs.size());
    int splitAfter = n - 1;
    double maxGap = -1.0;
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        double g = (j == 0) ? (angs[0] + kTwoPi - angs[i]) : (angs[j] - angs[i]);
        if (g > maxGap) {
            maxGap = g;
            splitAfter = i;
        }
    }
    if (!(maxGap > delta)) {
        stationAng.push_back(angs[0]);
        return true;
    }

    const int start = (splitAfter + 1) % n;
    double sx = 0.0, sy = 0.0;
    int cnt = 0;
    auto flush = [&]() {
        if (cnt <= 0) return;
        double a = std::atan2(sx, sy);
        if (a < 0.0) a += kTwoPi;
        if (a >= kTwoPi) a = 0.0;
        stationAng.push_back(a);
        sx = 0.0;
        sy = 0.0;
        cnt = 0;
    };
    for (int s = 0; s < n; ++s) {
        const int idx = (start + s) % n;
        const double a = angs[idx];
        if (cnt == 0) {
            sx = std::sin(a);
            sy = std::cos(a);
            cnt = 1;
            continue;
        }
        const int prev = (start + s - 1) % n;
        double g = angs[idx] - angs[prev];
        if (g < 0.0) g += kTwoPi;
        if (g <= delta) {
            sx += std::sin(a);
            sy += std::cos(a);
            ++cnt;
        } else {
            flush();
            sx = std::sin(a);
            sy = std::cos(a);
            cnt = 1;
        }
    }
    flush();
    return true;
}

// k steps on the occupied arc, θ_obs = largest of those steps.
// A closed cylinder counts every gap, including the seam. An open arc drops
// the single exterior gap (the walk was started across it).
bool arcStep(const std::vector<double>& st, bool closed360, int& k, double& thetaObs) {
    const int m = static_cast<int>(st.size());
    if (m < 2) return false;
    std::vector<double> gaps(static_cast<size_t>(m), 0.0);
    for (int i = 0; i < m; ++i) {
        const int j = (i + 1) % m;
        double g = st[static_cast<size_t>(j)] - st[static_cast<size_t>(i)];
        if (g <= 0.0) g += kTwoPi;
        gaps[static_cast<size_t>(i)] = g;
    }
    if (closed360) {
        k = m;
        thetaObs = *std::max_element(gaps.begin(), gaps.end());
        return k >= 2 && thetaObs > 0.0 && std::isfinite(thetaObs);
    }
    if (m < 3) return false;
    thetaObs = 0.0;
    for (int i = 0; i < m - 1; ++i)
        thetaObs = std::max(thetaObs, gaps[static_cast<size_t>(i)]);
    k = m - 1;
    return k >= 2 && thetaObs > 0.0 && std::isfinite(thetaObs);
}

struct Meas {
    int accIdx = -1;
    int minTri = -1;
    bool isNgon = false;
    bool certified = false;
    double rRefit = 0.0;
    double rStored = 0.0;
    int N = 0;
    int stations = 0;
    int k = 0;
    double thetaObs = 0.0;
};

bool measureRef(const MeshView& mv, const SegmentParams& p, double tau, int accIdx,
                const Region& reg, Meas& m) {
    m = Meas{};
    m.accIdx = accIdx;
    m.minTri = minTriOf(reg.tris);
    m.isNgon = reg.origin == Origin::NgonWall;
    m.rStored = reg.radius;
    m.N = reg.nSides;
    if (!isLawCylinder(reg)) return false;

    gp_Dir axis;
    gp_Pnt center;
    double radius = 0.0, maxResid = 0.0;
    const int iters = p.maxRefineIters > 0 ? p.maxRefineIters : 16;
    if (!cylinderFitLS(mv, reg.tris, reg.ax.Direction(), axis, center, radius, maxResid, iters))
        return true;
    if (!(radius > 0.0) || !std::isfinite(radius) || !std::isfinite(maxResid)) return true;
    m.rRefit = radius;
    if (!(maxResid <= tau)) return true;

    std::vector<double> stations;
    if (!azimuthStations(mv, reg.tris, axis, center, radius, tau, stations)) return true;
    m.stations = static_cast<int>(stations.size());
    if (m.stations < minStations(kParamsCircle)) return true;
    int k = 0;
    double theta = 0.0;
    if (!arcStep(stations, reg.closed360, k, theta)) return true;
    m.k = k;
    m.thetaObs = theta;
    m.certified = true;
    return true;
}

// Witness reference: a certified cylinder, not this N-gon, not a refuted
// N-gon, with r ≤ R_N. `hits` when 2π/N is at least the law bound.
struct Witness {
    bool any = false;
    bool hits = false;
    int minTri = -1;
    double r = 0.0;
    double theta = 0.0;
    int k = 0;
    double bound = 0.0;
};

bool refutes(const Meas& ref, const Meas& ngon) {
    if (!ref.certified || ref.k < 2 || ref.accIdx == ngon.accIdx) return false;
    if (!(ref.rRefit <= ngon.rStored) || !(ngon.rStored > 0.0) || ngon.N < 3) return false;
    const double bound = ref.thetaObs * static_cast<double>(ref.k) / static_cast<double>(ref.k - 1);
    const double step = kTwoPi / static_cast<double>(ngon.N);
    return std::isfinite(bound) && step >= bound;
}

void consider(const Meas& ref, const Meas& ngon, Witness& w) {
    if (!ref.certified || ref.k < 2 || ref.accIdx == ngon.accIdx) return;
    if (!(ref.rRefit <= ngon.rStored) || !(ngon.rStored > 0.0) || ngon.N < 3) return;
    const double bound = ref.thetaObs * static_cast<double>(ref.k) / static_cast<double>(ref.k - 1);
    if (!std::isfinite(bound)) return;
    const double step = kTwoPi / static_cast<double>(ngon.N);
    const bool hit = step >= bound;
    const bool betterHit = hit && (!w.hits || ref.rRefit < w.r
                                   || (ref.rRefit == w.r && ref.minTri < w.minTri));
    const bool betterKeep = !hit && !w.hits && (!w.any || bound < w.bound
                                                || (bound == w.bound && ref.minTri < w.minTri));
    if (!w.any || betterHit || betterKeep) {
        w.any = true;
        w.hits = hit;
        w.minTri = ref.minTri;
        w.r = ref.rRefit;
        w.theta = ref.thetaObs;
        w.k = ref.k;
        w.bound = bound;
    } else if (hit && w.hits && ref.rRefit == w.r && ref.minTri < w.minTri) {
        w.minTri = ref.minTri;
        w.theta = ref.thetaObs;
        w.k = ref.k;
        w.bound = bound;
    }
}

void recomputeProv(const MeshView& mv, Provisional& prov) {
    prov.area = 0.0;
    double maxD = 0.0, sumSq = 0.0;
    int nV = 0;
    for (int t : prov.tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        prov.area += triArea(mv, t);
        for (int k = 0; k < 3; ++k) {
            const double d = pointPlaneDist(prov.plane, triVert(mv, t, k));
            maxD = std::max(maxD, d);
            sumSq += d * d;
            ++nV;
        }
    }
    prov.maxVertexDev = maxD;
    prov.rmsVertexDev = nV > 0 ? std::sqrt(sumSq / static_cast<double>(nV)) : 0.0;
}

gp_Ax3 planeFromTri(const MeshView& mv, int lt) {
    const gp_XYZ a = triVert(mv, lt, 0);
    const gp_XYZ b = triVert(mv, lt, 1);
    const gp_XYZ c = triVert(mv, lt, 2);
    gp_XYZ n = (b - a).Crossed(c - a);
    const double m = n.Modulus();
    const gp_XYZ g = (a + b + c) / 3.0;
    if (!(m > 0.0)) return gp_Ax3(gp_Pnt(g.X(), g.Y(), g.Z()), gp_Dir(0, 0, 1));
    n.Divide(m);
    return gp_Ax3(gp_Pnt(g.X(), g.Y(), g.Z()), gp_Dir(n.X(), n.Y(), n.Z()));
}

bool coplanarPair(const MeshView& mv, int t0, int t1, double tau) {
    const gp_Ax3 p0 = planeFromTri(mv, t0);
    const gp_Ax3 p1 = planeFromTri(mv, t1);
    return triPlaneDist(mv, t1, p0) <= tau && triPlaneDist(mv, t0, p1) <= tau;
}

// Put every triangle of a refuted N-gon back on an Unclaimed provisional.
// Production claimNgonWallsA clears the member provisional's triangle list
// and keeps the A2 plane, so membership is recovered by that plane within τ.
// Triangles that match none are grouped into τ-coplanar edge components and
// appended. Returns false if a triangle could not be homed.
bool releaseToUnclaimed(const MeshView& mv, double tau, const Region& region, SegmentWork& work) {
    std::vector<char> inReg(mv.nTri, 0);
    int nNeed = 0;
    for (int t : region.tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        if (!inReg[static_cast<size_t>(t)]) {
            inReg[static_cast<size_t>(t)] = 1;
            ++nNeed;
        }
    }
    if (nNeed == 0) return true;

    std::vector<char> homed(mv.nTri, 0);
    int nHomed = 0;
    const int nProv = static_cast<int>(work.provisionals.size());

    for (int pi = 0; pi < nProv; ++pi) {
        Provisional& prov = work.provisionals[static_cast<size_t>(pi)];
        if (prov.tris.empty()) continue;
        bool subset = true;
        for (int t : prov.tris) {
            if (t < 0 || static_cast<size_t>(t) >= mv.nTri || !inReg[static_cast<size_t>(t)]) {
                subset = false;
                break;
            }
        }
        if (!subset) continue;
        prov.claim = ProvClaim::Unclaimed;
        std::sort(prov.tris.begin(), prov.tris.end());
        prov.tris.erase(std::unique(prov.tris.begin(), prov.tris.end()), prov.tris.end());
        recomputeProv(mv, prov);
        for (int t : prov.tris) {
            if (!homed[static_cast<size_t>(t)]) {
                homed[static_cast<size_t>(t)] = 1;
                ++nHomed;
            }
        }
    }

    std::vector<std::vector<int>> assign(static_cast<size_t>(nProv));
    std::vector<int> leftover;
    std::vector<int> need;
    need.reserve(static_cast<size_t>(nNeed));
    for (int t = 0; t < static_cast<int>(mv.nTri); ++t)
        if (inReg[static_cast<size_t>(t)] && !homed[static_cast<size_t>(t)])
            need.push_back(t);

    for (int t : need) {
        int best = -1;
        double bestScore = 0.0;
        for (int pi = 0; pi < nProv; ++pi) {
            const Provisional& prov = work.provisionals[static_cast<size_t>(pi)];
            if (prov.claim != ProvClaim::ConsumedCylinder || !prov.tris.empty()) continue;
            const double score = triPlaneDist(mv, t, prov.plane);
            if (!(score <= tau)) continue;
            if (best < 0 || score < bestScore) {
                best = pi;
                bestScore = score;
            }
        }
        if (best >= 0) assign[static_cast<size_t>(best)].push_back(t);
        else leftover.push_back(t);
    }
    for (int pi = 0; pi < nProv; ++pi) {
        if (assign[static_cast<size_t>(pi)].empty()) continue;
        Provisional& prov = work.provisionals[static_cast<size_t>(pi)];
        prov.tris = std::move(assign[static_cast<size_t>(pi)]);
        std::sort(prov.tris.begin(), prov.tris.end());
        prov.claim = ProvClaim::Unclaimed;
        recomputeProv(mv, prov);
        for (int t : prov.tris) {
            if (!homed[static_cast<size_t>(t)]) {
                homed[static_cast<size_t>(t)] = 1;
                ++nHomed;
            }
        }
    }

    if (!leftover.empty()) {
        const int nL = static_cast<int>(leftover.size());
        std::vector<int> parent(static_cast<size_t>(nL));
        for (int i = 0; i < nL; ++i) parent[static_cast<size_t>(i)] = i;
        auto find = [&](int x) {
            int r = x;
            while (parent[static_cast<size_t>(r)] != r) r = parent[static_cast<size_t>(r)];
            while (parent[static_cast<size_t>(x)] != r) {
                const int n = parent[static_cast<size_t>(x)];
                parent[static_cast<size_t>(x)] = r;
                x = n;
            }
            return r;
        };
        auto unite = [&](int a, int b) {
            a = find(a);
            b = find(b);
            if (a == b) return;
            if (a > b) std::swap(a, b);
            parent[static_cast<size_t>(b)] = a;
        };
        std::vector<std::pair<std::pair<int, int>, int>> edges;
        for (int i = 0; i < nL; ++i) {
            const int t = leftover[static_cast<size_t>(i)];
            int v[3] = {mv.tris[mv.compTris[t]][0], mv.tris[mv.compTris[t]][1],
                        mv.tris[mv.compTris[t]][2]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[(s + 1) % 3]);
                const int hi = std::max(v[s], v[(s + 1) % 3]);
                edges.push_back({{lo, hi}, i});
            }
        }
        std::sort(edges.begin(), edges.end());
        for (size_t i = 0; i < edges.size();) {
            size_t j = i + 1;
            while (j < edges.size() && edges[j].first == edges[i].first) ++j;
            for (size_t a = i; a < j; ++a) {
                for (size_t b = a + 1; b < j; ++b) {
                    const int ia = edges[a].second;
                    const int ib = edges[b].second;
                    if (coplanarPair(mv, leftover[static_cast<size_t>(ia)],
                                     leftover[static_cast<size_t>(ib)], tau))
                        unite(ia, ib);
                }
            }
            i = j;
        }
        std::vector<std::vector<int>> comps(static_cast<size_t>(nL));
        for (int i = 0; i < nL; ++i)
            comps[static_cast<size_t>(find(i))].push_back(leftover[static_cast<size_t>(i)]);
        std::vector<std::vector<int>> groups;
        for (auto& c : comps)
            if (!c.empty()) groups.push_back(std::move(c));
        std::sort(groups.begin(), groups.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
            return minTriOf(a) < minTriOf(b);
        });
        for (auto& g : groups) {
            std::sort(g.begin(), g.end());
            Provisional np;
            np.claim = ProvClaim::Unclaimed;
            np.chartId = -1;
            np.tris = g;
            if (!pcaPlane(mv, np.tris, np.plane))
                np.plane = planeFromTri(mv, np.tris.front());
            recomputeProv(mv, np);
            for (int t : np.tris) {
                if (!homed[static_cast<size_t>(t)]) {
                    homed[static_cast<size_t>(t)] = 1;
                    ++nHomed;
                }
            }
            work.provisionals.push_back(std::move(np));
        }
    }
    return nHomed == nNeed;
}

struct DiagLine {
    int minTri = -1;
    int accIdx = -1;
    int N = 0;
    double R = 0.0;
    double step = 0.0;
    int refRid = -1;
    double refR = 0.0;
    double refTheta = 0.0;
    int refK = 0;
    double bound = 0.0;
    const char* verdict = "no-ref";
};

}  // namespace

bool refuteNgonByLaw(const MeshView& mv, const SegmentParams& p, const DerivedTols&,
                     SegmentWork& work) {
    const double tau = tauOf(mv);
    const bool diag = tesslawDiagOn();
    if (!(tau > 0.0) || !meshReadable(mv)) {
        if (diag)
            std::fprintf(stderr, "DIAG_TESSLAW summary n=0 refute=0 keep=0 noref=0\n");
        return true;
    }

    bool anyNgon = false;
    for (const Region& r : work.accepted) {
        if (r.origin == Origin::NgonWall && r.type == SurfType::Cylinder && !r.tris.empty()) {
            anyNgon = true;
            break;
        }
    }
    if (!anyNgon) {
        if (diag)
            std::fprintf(stderr, "DIAG_TESSLAW summary n=0 refute=0 keep=0 noref=0\n");
        return true;
    }

    std::vector<Meas> meas;
    meas.reserve(work.accepted.size());
    for (int i = 0; i < static_cast<int>(work.accepted.size()); ++i) {
        const Region& r = work.accepted[static_cast<size_t>(i)];
        if (!isLawCylinder(r) && !(r.origin == Origin::NgonWall && r.type == SurfType::Cylinder))
            continue;
        Meas m;
        measureRef(mv, p, tau, i, r, m);
        if (m.isNgon || m.certified) meas.push_back(m);
    }

    std::vector<int> nonNgon;
    std::vector<int> ngons;
    for (int i = 0; i < static_cast<int>(meas.size()); ++i) {
        if (!meas[static_cast<size_t>(i)].isNgon) {
            if (meas[static_cast<size_t>(i)].certified) nonNgon.push_back(i);
        } else {
            ngons.push_back(i);
        }
    }

    // Least refutation: a refuted N-gon is not a reference (D-fc-1.9 / SPEC §5).
    // Non-N-gon certificates seed the set. An N-gon stays a reference only when
    // no still-kept reference refutes it, so two N-gons that only refute each
    // other are both kept.
    std::vector<char> isRef(meas.size(), 0);
    for (int gi : ngons) {
        for (int u : nonNgon) {
            if (refutes(meas[static_cast<size_t>(u)], meas[static_cast<size_t>(gi)])) {
                isRef[static_cast<size_t>(gi)] = 1;
                break;
            }
        }
    }
    std::vector<char> grounded(meas.size(), 0);
    bool changed = true;
    for (int guard = 0; changed && guard < static_cast<int>(ngons.size()) + 2; ++guard) {
        changed = false;
        std::fill(grounded.begin(), grounded.end(), 0);
        for (int gi : ngons) {
            if (!meas[static_cast<size_t>(gi)].certified || isRef[static_cast<size_t>(gi)])
                continue;
            bool attacked = false;
            for (int hj : ngons) {
                if (hj == gi || isRef[static_cast<size_t>(hj)]) continue;
                if (!meas[static_cast<size_t>(hj)].certified) continue;
                if (refutes(meas[static_cast<size_t>(hj)], meas[static_cast<size_t>(gi)])) {
                    attacked = true;
                    break;
                }
            }
            if (!attacked) grounded[static_cast<size_t>(gi)] = 1;
        }
        for (int gi : ngons) {
            if (isRef[static_cast<size_t>(gi)]) continue;
            for (int hj : ngons) {
                if (!grounded[static_cast<size_t>(hj)] || hj == gi) continue;
                if (refutes(meas[static_cast<size_t>(hj)], meas[static_cast<size_t>(gi)])) {
                    isRef[static_cast<size_t>(gi)] = 1;
                    changed = true;
                    break;
                }
            }
        }
    }
    // Grounded is the set of N-gon references that stand without a circular
    // witness. A clique that only refutes itself has no grounded member and
    // is kept (a refuted N-gon is not a reference).
    std::fill(grounded.begin(), grounded.end(), 0);
    for (int gi : ngons) {
        if (!meas[static_cast<size_t>(gi)].certified || isRef[static_cast<size_t>(gi)]) continue;
        bool attacked = false;
        for (int hj : ngons) {
            if (hj == gi || isRef[static_cast<size_t>(hj)]) continue;
            if (!meas[static_cast<size_t>(hj)].certified) continue;
            if (refutes(meas[static_cast<size_t>(hj)], meas[static_cast<size_t>(gi)])) {
                attacked = true;
                break;
            }
        }
        if (!attacked) grounded[static_cast<size_t>(gi)] = 1;
    }

    std::vector<DiagLine> lines;
    std::vector<int> drop;
    int nRefute = 0, nKeep = 0, nNoref = 0;
    for (int gi : ngons) {
        const Meas& g = meas[static_cast<size_t>(gi)];
        Witness w;
        for (int u : nonNgon)
            consider(meas[static_cast<size_t>(u)], g, w);
        for (int hj : ngons) {
            if (!grounded[static_cast<size_t>(hj)] || hj == gi) continue;
            consider(meas[static_cast<size_t>(hj)], g, w);
        }
        DiagLine line;
        line.minTri = g.minTri;
        line.accIdx = g.accIdx;
        line.N = g.N;
        line.R = g.rStored;
        line.step = g.N >= 3 ? kTwoPi / static_cast<double>(g.N) : 0.0;
        if (w.hits) {
            line.verdict = "refute";
            line.refRid = w.minTri;
            line.refR = w.r;
            line.refTheta = w.theta;
            line.refK = w.k;
            line.bound = w.bound;
            ++nRefute;
            drop.push_back(g.accIdx);
        } else if (w.any) {
            line.verdict = "keep";
            line.refRid = w.minTri;
            line.refR = w.r;
            line.refTheta = w.theta;
            line.refK = w.k;
            line.bound = w.bound;
            ++nKeep;
        } else {
            line.verdict = "no-ref";
            ++nNoref;
        }
        lines.push_back(line);
    }

    std::sort(lines.begin(), lines.end(), [](const DiagLine& a, const DiagLine& b) {
        if (a.minTri != b.minTri) return a.minTri < b.minTri;
        return a.accIdx < b.accIdx;
    });
    std::sort(drop.begin(), drop.end());
    drop.erase(std::unique(drop.begin(), drop.end()), drop.end());

    // Home triangles before erasing, highest accepted index first.
    std::vector<int> erased;
    for (int i = static_cast<int>(drop.size()) - 1; i >= 0; --i) {
        const int idx = drop[static_cast<size_t>(i)];
        if (idx < 0 || static_cast<size_t>(idx) >= work.accepted.size()) continue;
        if (!releaseToUnclaimed(mv, tau, work.accepted[static_cast<size_t>(idx)], work))
            continue;
        work.accepted.erase(work.accepted.begin() + idx);
        erased.push_back(idx);
    }
    if (static_cast<int>(erased.size()) != static_cast<int>(drop.size())) {
        // A wall that could not be homed stays accepted. Correct the diag so
        // the line matches the region that is still in the set.
        for (int idx : drop) {
            if (std::find(erased.begin(), erased.end(), idx) != erased.end()) continue;
            for (DiagLine& line : lines) {
                if (line.accIdx == idx && line.verdict[0] == 'r') {
                    line.verdict = "keep";
                    --nRefute;
                    ++nKeep;
                }
            }
        }
    }

    if (diag) {
        for (const DiagLine& line : lines) {
            std::fprintf(stderr,
                         "DIAG_TESSLAW rid=%d N=%d R=%.8g step=%.8g refRid=%d refR=%.8g "
                         "refThetaObs=%.8g refK=%d bound=%.8g verdict=%s\n",
                         line.minTri, line.N, line.R, line.step, line.refRid, line.refR,
                         line.refTheta, line.refK, line.bound, line.verdict);
        }
        std::fprintf(stderr, "DIAG_TESSLAW summary n=%d refute=%d keep=%d noref=%d\n",
                     nRefute + nKeep + nNoref, nRefute, nKeep, nNoref);
    }
    return true;
}

}}  // namespace stl2step::refit
