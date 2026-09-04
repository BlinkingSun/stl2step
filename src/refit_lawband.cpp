// Law-band recognizer core (AC2-L1). Pure functions; not wired into runStages.
// Inverse R = median w_i/(2 sin(theta_i/2)). Parameter-free Tier-1 accept.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>
#include <vector>

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kCvThetaMax = 1e-3;
constexpr double kRelRMax = 5e-4;
constexpr double kCvRMax = 1e-6;
constexpr double kTauSurfFloor = 5e-5;

bool finite3(double x, double y, double z) {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

double wrapPi(double t) {
    while (t <= -kPi) t += kTwoPi;
    while (t > kPi) t -= kTwoPi;
    return t;
}

double wrapTwoPi(double t) {
    t = std::fmod(t, kTwoPi);
    if (t < 0.0) t += kTwoPi;
    return t;
}

double tauSurf(const MeshView& mv) {
    return std::max({kTauSurfFloor, 4.0 * mv.weldTol, 1e-6 * mv.diag});
}

double linTol(const MeshView& mv, const DerivedTols& tol) {
    return std::max({tol.epsMesh, 4.0 * mv.weldTol, 1e-6 * mv.diag, kTauSurfFloor});
}

bool diagOn() {
    const char* e = std::getenv("STL2STEP_LAWBAND_DIAG");
    return e && e[0] == '1' && e[1] == '\0';
}

// Every virtual-generator ATTEMPT and the named reason it was refused. Its own
// switch: the seed sweep tries thousands of triples, so this is a probe, not
// part of STL2STEP_LAWBAND_DIAG's per-band census.
bool vgenDiagOn() {
    const char* e = std::getenv("STL2STEP_LAWVGEN_DIAG");
    return e && e[0] == '1' && e[1] == '\0';
}

std::mutex& diagMu() {
    static std::mutex m;
    return m;
}

double popMean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double popCV(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double m = popMean(v);
    if (!(std::abs(m) > 0.0)) return 0.0;
    double ss = 0.0;
    for (double x : v) {
        const double d = x - m;
        ss += d * d;
    }
    const double sd = std::sqrt(ss / static_cast<double>(v.size()));
    return sd / std::abs(m);
}

double medianOf(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n & 1u) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

void signNormalize(std::array<double, 3>& v) {
    int imax = 0;
    double amax = std::abs(v[0]);
    for (int i = 1; i < 3; ++i) {
        const double ai = std::abs(v[i]);
        if (ai > amax) {
            amax = ai;
            imax = i;
        }
    }
    if (v[imax] < 0.0) {
        v[0] = -v[0];
        v[1] = -v[1];
        v[2] = -v[2];
    }
}

bool axisFrame(const gp_XYZ& aUnit, gp_XYZ& u, gp_XYZ& v) {
    const double nx = std::abs(aUnit.X());
    const double ny = std::abs(aUnit.Y());
    const double nz = std::abs(aUnit.Z());
    gp_XYZ w;
    if (nx <= ny && nx <= nz)
        w = gp_XYZ(1.0, 0.0, 0.0);
    else if (ny <= nz)
        w = gp_XYZ(0.0, 1.0, 0.0);
    else
        w = gp_XYZ(0.0, 0.0, 1.0);
    u = aUnit.Crossed(w);
    const double um = u.Modulus();
    if (um < 1e-15) return false;
    u.Divide(um);
    v = aUnit.Crossed(u);
    const double vm = v.Modulus();
    if (vm < 1e-15) return false;
    v.Divide(vm);
    return true;
}

bool triCorners(const MeshView& mv, int localTri, gp_XYZ& a, gp_XYZ& b, gp_XYZ& c) {
    if (!mv.pts || !mv.tris || !mv.compTris) return false;
    if (localTri < 0 || static_cast<size_t>(localTri) >= mv.nTri) return false;
    const int g = mv.compTris[localTri];
    a = mv.pts[mv.tris[g][0]];
    b = mv.pts[mv.tris[g][1]];
    c = mv.pts[mv.tris[g][2]];
    return finite3(a.X(), a.Y(), a.Z()) && finite3(b.X(), b.Y(), b.Z()) &&
           finite3(c.X(), c.Y(), c.Z());
}

bool unitTriNormal(const MeshView& mv, int lt, gp_XYZ& nOut, double& areaOut) {
    gp_XYZ a, b, c;
    if (!triCorners(mv, lt, a, b, c)) return false;
    const gp_XYZ nU = (b - a).Crossed(c - a);
    areaOut = 0.5 * nU.Modulus();
    if (!(areaOut > 0.0)) return false;
    const double nm = nU.Modulus();
    if (nm < 1e-15) return false;
    nOut = nU;
    nOut.Divide(nm);
    return true;
}

void sortedUnique(const std::vector<int>& in, std::vector<int>& out) {
    out = in;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

bool recoverAxisDir(const MeshView& mv, const std::vector<int>& ids, gp_XYZ& axis) {
    gp_XYZ nbar(0.0, 0.0, 0.0);
    double areaSum = 0.0;
    std::vector<gp_XYZ> ns;
    std::vector<double> areas;
    ns.reserve(ids.size());
    areas.reserve(ids.size());
    for (int t : ids) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, t, n, area)) continue;
        areaSum += area;
        nbar += n * area;
        ns.push_back(n);
        areas.push_back(area);
    }
    if (!(areaSum > 0.0) || ns.size() < 2) return false;
    nbar.Divide(areaSum);

    std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (size_t i = 0; i < ns.size(); ++i) {
        const double dx = ns[i].X() - nbar.X();
        const double dy = ns[i].Y() - nbar.Y();
        const double dz = ns[i].Z() - nbar.Z();
        const double w = areas[i];
        cov[0][0] += w * dx * dx;
        cov[0][1] += w * dx * dy;
        cov[0][2] += w * dx * dz;
        cov[1][1] += w * dy * dy;
        cov[1][2] += w * dy * dz;
        cov[2][2] += w * dz * dz;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return false;
    axis = gp_XYZ(evec[0][0], evec[0][1], evec[0][2]);
    const double am = axis.Modulus();
    if (am < 1e-15) return false;
    axis.Divide(am);
    std::array<double, 3> av{axis.X(), axis.Y(), axis.Z()};
    signNormalize(av);
    axis = gp_XYZ(av[0], av[1], av[2]);
    return true;
}

void collectUniqueVerts(const MeshView& mv, const std::vector<int>& ids,
                        std::vector<gp_XYZ>& out) {
    out.clear();
    if (!mv.pts || !mv.tris || !mv.compTris) return;
    std::vector<int> gids;
    gids.reserve(ids.size() * 3);
    for (int t : ids) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        const int g = mv.compTris[t];
        gids.push_back(mv.tris[g][0]);
        gids.push_back(mv.tris[g][1]);
        gids.push_back(mv.tris[g][2]);
    }
    std::sort(gids.begin(), gids.end());
    gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
    out.reserve(gids.size());
    for (int g : gids) out.push_back(mv.pts[g]);
}

void cluster1d(std::vector<double> vals, double tol, std::vector<double>& modes) {
    modes.clear();
    if (vals.empty()) return;
    std::sort(vals.begin(), vals.end());
    std::vector<double> cur{vals[0]};
    for (size_t i = 1; i < vals.size(); ++i) {
        if (std::abs(vals[i] - cur.back()) <= tol)
            cur.push_back(vals[i]);
        else {
            modes.push_back(popMean(cur));
            cur = {vals[i]};
        }
    }
    modes.push_back(popMean(cur));
}

void uniquePositions(const std::vector<gp_XYZ>& in, double merge, std::vector<gp_XYZ>& out) {
    out.clear();
    for (const gp_XYZ& p : in) {
        bool hit = false;
        for (gp_XYZ& q : out) {
            if ((p - q).Modulus() <= merge) {
                q = (q + p) * 0.5;
                hit = true;
                break;
            }
        }
        if (!hit) out.push_back(p);
    }
}

bool lsBisectorCenter(const std::vector<double>& x, const std::vector<double>& y, double& cx,
                      double& cy) {
    const size_t n = x.size();
    if (n < 3) return false;
    double a00 = 0, a01 = 0, a11 = 0, b0 = 0, b1 = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double dx = x[j] - x[i];
            const double dy = y[j] - y[i];
            const double L = std::sqrt(dx * dx + dy * dy);
            if (L < 1e-12) continue;
            const double nx = dx / L;
            const double ny = dy / L;
            const double mx = 0.5 * (x[i] + x[j]);
            const double my = 0.5 * (y[i] + y[j]);
            const double w = L * L;
            a00 += w * nx * nx;
            a01 += w * nx * ny;
            a11 += w * ny * ny;
            b0 += w * nx * (nx * mx + ny * my);
            b1 += w * ny * (nx * mx + ny * my);
        }
    }
    const double det = a00 * a11 - a01 * a01;
    if (!(std::abs(det) > 1e-18)) return false;
    cx = (a11 * b0 - a01 * b1) / det;
    cy = (a00 * b1 - a01 * b0) / det;
    return std::isfinite(cx) && std::isfinite(cy);
}

bool circumcenter2(double ax, double ay, double bx, double by, double cx, double cy, double& ox,
                   double& oy) {
    const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (!(std::abs(d) > 1e-18)) return false;
    const double a2 = ax * ax + ay * ay;
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    ox = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    oy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    return std::isfinite(ox) && std::isfinite(oy);
}

bool circumMedianCenter(const std::vector<double>& x, const std::vector<double>& y, double& cx,
                        double& cy) {
    const size_t n = x.size();
    if (n < 3) return false;
    std::vector<double> xs, ys, areas;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            for (size_t k = j + 1; k < n; ++k) {
                const double area = 0.5 * std::abs((x[j] - x[i]) * (y[k] - y[i]) -
                                                   (y[j] - y[i]) * (x[k] - x[i]));
                if (area < 1e-10) continue;
                double ox = 0, oy = 0;
                if (!circumcenter2(x[i], y[i], x[j], y[j], x[k], y[k], ox, oy)) continue;
                xs.push_back(ox);
                ys.push_back(oy);
                areas.push_back(area);
            }
        }
    }
    if (xs.empty()) return false;
    std::vector<size_t> ord(xs.size());
    for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return areas[a] > areas[b]; });
    const size_t take = std::max<size_t>(1, ord.size() / 4);
    std::vector<double> tx, ty;
    tx.reserve(take);
    ty.reserve(take);
    for (size_t i = 0; i < take; ++i) {
        tx.push_back(xs[ord[i]]);
        ty.push_back(ys[ord[i]]);
    }
    cx = medianOf(tx);
    cy = medianOf(ty);
    return true;
}

double azimuth(const gp_XYZ& p, const gp_XYZ& origin, const gp_XYZ& axis, const gp_XYZ& u,
               const gp_XYZ& v) {
    const gp_XYZ d = p - origin;
    const gp_XYZ rad = d - axis * d.Dot(axis);
    return std::atan2(rad.Dot(v), rad.Dot(u));
}

double rhoOf(const gp_XYZ& p, const gp_XYZ& origin, const gp_XYZ& axis) {
    const gp_XYZ d = p - origin;
    const gp_XYZ rad = d - axis * d.Dot(axis);
    return rad.Modulus();
}

void clusterAngles(std::vector<double> angs, double tol, std::vector<double>& gens) {
    gens.clear();
    if (angs.empty()) return;
    for (double& a : angs) a = wrapTwoPi(a);
    std::sort(angs.begin(), angs.end());
    std::vector<std::vector<double>> cl;
    cl.push_back({angs[0]});
    for (size_t i = 1; i < angs.size(); ++i) {
        if (std::abs(wrapPi(angs[i] - cl.back().back())) <= tol)
            cl.back().push_back(angs[i]);
        else
            cl.push_back({angs[i]});
    }
    if (cl.size() > 1 && std::abs(wrapPi(cl.front().front() - cl.back().back())) <= tol) {
        cl.front().insert(cl.front().begin(), cl.back().begin(), cl.back().end());
        cl.pop_back();
    }
    gens.reserve(cl.size());
    for (const auto& c : cl) {
        // circular mean via vector average
        double sx = 0, sy = 0;
        for (double a : c) {
            sx += std::cos(a);
            sy += std::sin(a);
        }
        gens.push_back(wrapTwoPi(std::atan2(sy, sx)));
    }
    std::sort(gens.begin(), gens.end());
}

struct StripPair {
    double span = 0;
    double a0 = 0;
    double a1 = 0;
};

void makeIntervals(const std::vector<double>& gens, std::vector<StripPair>& pairs, bool& closed) {
    pairs.clear();
    closed = false;
    const size_t nG = gens.size();
    if (nG < 2) return;
    for (size_t i = 0; i < nG; ++i) {
        const double a0 = gens[i];
        const double a1 = gens[(i + 1) % nG];
        const double span = wrapTwoPi(a1 - a0);
        if (span > 1e-12) pairs.push_back({span, a0, a1});
    }
    if (pairs.empty()) return;
    double mx = pairs[0].span;
    for (const auto& p : pairs) mx = std::max(mx, p.span);
    std::vector<double> others;
    for (const auto& p : pairs)
        if (p.span < mx - 1e-12) others.push_back(p.span);
    if (others.empty()) others.push_back(mx);
    const double med = medianOf(others);
    closed = (nG >= 3) && (mx < 1.8 * med);
    if (!closed && pairs.size() >= 2) {
        std::sort(pairs.begin(), pairs.end(),
                  [](const StripPair& a, const StripPair& b) { return a.span > b.span; });
        pairs.erase(pairs.begin());  // drop unused circular sector
        std::sort(pairs.begin(), pairs.end(),
                  [](const StripPair& a, const StripPair& b) { return a.a0 < b.a0; });
    }
}

gp_XYZ nearestAtAngle(const std::vector<gp_XYZ>& pts, double ang, const gp_XYZ& origin,
                      const gp_XYZ& axis, const gp_XYZ& u, const gp_XYZ& v) {
    gp_XYZ best = pts.front();
    double bd = 1e300;
    for (const gp_XYZ& p : pts) {
        const double d = std::abs(wrapPi(azimuth(p, origin, axis, u, v) - ang));
        if (d < bd) {
            bd = d;
            best = p;
        }
    }
    return best;
}

double axisLineSep(const gp_Ax1& a, const gp_Ax1& b) {
    const gp_XYZ d = b.Location().XYZ() - a.Location().XYZ();
    const gp_XYZ da = a.Direction().XYZ();
    const gp_XYZ db = b.Direction().XYZ();
    const gp_XYZ cr = da.Crossed(db);
    const double cm = cr.Modulus();
    if (cm < 1e-12) return da.Crossed(d).Modulus();
    return std::abs(d.Dot(cr)) / cm;
}

double tauQuant(const MeshView& mv);

// ---------------------------------------------------------------------------
// D-130-8 -- VIRTUAL GENERATORS.
//
// On a tessellated cylinder every facet is a chord plane, and the line where
// two adjacent chord planes meet is the inscribed polygon's own corner line:
// it lies exactly ON the surface and is exactly parallel to the axis. That
// line is a generator whether or not the mesh drew an edge along it, so it is
// the generator of a STAGGERED strip -- one whose two rings sit at alternating
// azimuths, where at most azimuths there is no axis-parallel mesh edge to pair
// on and the end-ring pairing above stops at the first regular run (the
// plate's R=8.5 bore, its four R=3 slot ends).
//
// What the intersection lines resolve is the DIRECTION: they are parallel to
// the axis by construction, so their mean carries it to ~1e-8 rad where the
// normal-covariance PCA of a 39-facet half-wall leaves 1e-4 rad -- and 1e-4 rad
// across a 15 mm band is a 1.5e-3 mm radial error, an order of magnitude past
// any certificate. The AZIMUTHS then come from the band's own vertices: every
// vertex of an inscribed tessellation lies on one of those same corner lines,
// and unlike the mesh edges there is one at every azimuth, staggered or not.
//
// PARALLELISM is the scale-relative test testCommonAxis already uses for a mesh
// generator edge -- the line's lateral drift across the band's own axial span
// must stay inside the band's linear tolerance. No degree gate (D-130-8).
struct VGenAxis {
    std::vector<gp_XYZ> pts;  // one point on each kept generator line
    gp_XYZ dir = gp_XYZ(0, 0, 1);
    int nKept = 0;
    int nLines = 0;
    int nCand = 0;
    double tiltMax = 0.0;
    double spreadRad = 0.0;
    bool ok = false;
};

VGenAxis virtualGeneratorAxis(const MeshView& mv, const std::vector<int>& ids,
                              const gp_XYZ& axisSeed, double span, double lim, double qFloor) {
    VGenAxis vg;
    if (!mv.triEdges || ids.size() < 2 || !(span > 0.0)) return vg;

    // Interior edges of the SET: an edge id seen by exactly two of its triangles.
    std::vector<std::pair<int, int>> et;
    et.reserve(ids.size() * 3);
    for (int t : ids) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        for (int s = 0; s < 3; ++s) et.push_back({mv.triEdges[t][s], t});
    }
    std::sort(et.begin(), et.end());

    // A line is a generator only if its lateral drift over the band's own axial
    // span stays inside the band's linear tolerance.
    const double tiltMax = lim / span;
    std::vector<gp_XYZ> dirs;
    dirs.reserve(et.size() / 2);
    vg.pts.reserve(et.size() / 2);
    for (size_t i = 0; i + 1 < et.size(); ++i) {
        if (et[i].first != et[i + 1].first) continue;
        if (i + 2 < et.size() && et[i + 2].first == et[i].first) {  // non-manifold
            while (i + 1 < et.size() && et[i + 1].first == et[i].first) ++i;
            continue;
        }
        const int t0 = et[i].second;
        const int t1 = et[i + 1].second;
        ++i;
        gp_XYZ n0, n1;
        double a0 = 0.0, a1 = 0.0;
        if (!unitTriNormal(mv, t0, n0, a0) || !unitTriNormal(mv, t1, n1, a1)) continue;
        gp_XYZ d = n0.Crossed(n1);
        const double m = d.Modulus();
        if (m < 1e-12) continue;  // coplanar pair: one facet split by a diagonal
        // Two facets are ONE facet unless the mesh resolves the fold: the far
        // corner of each triangle must stand off the other's plane by more than
        // the file's own coordinate quantization. Below that the cross product
        // is round-off and its direction is noise -- which is exactly what put
        // two junk lines beside every real one on a four-facet seed and lost
        // the plate's slot ends. Measured, not an angle.
        gp_XYZ a0v, b0v, c0v, a1v, b1v, c1v;
        if (!triCorners(mv, t0, a0v, b0v, c0v) || !triCorners(mv, t1, a1v, b1v, c1v)) continue;
        double fold = 0.0;
        const gp_XYZ w0[3] = {a1v, b1v, c1v};
        const gp_XYZ w1[3] = {a0v, b0v, c0v};
        for (int k = 0; k < 3; ++k) {
            fold = std::max(fold, std::abs((w0[k] - a0v).Dot(n0)));
            fold = std::max(fold, std::abs((w1[k] - a1v).Dot(n1)));
        }
        if (!(fold > qFloor)) continue;
        d.Divide(m);
        ++vg.nLines;
        if (d.Dot(axisSeed) < 0.0) d.Reverse();
        // The point of that line nearest the shared edge's midpoint: solve
        // n0.x = n0.p0, n1.x = n1.p1, d.x = d.m. It lies ON the surface, so it
        // is both the azimuth and the radius the chain law needs.
        const int va = mv.compEdges ? mv.compEdges[et[i].first].first : -1;
        const int vb2 = mv.compEdges ? mv.compEdges[et[i].first].second : -1;
        if (va < 0 || vb2 < 0 || !mv.compVtx) continue;
        const gp_XYZ mid = (mv.pts[mv.compVtx[va]] + mv.pts[mv.compVtx[vb2]]) * 0.5;
        const gp_XYZ rows[3] = {n0, n1, d};
        const double rhs[3] = {n0.Dot(a0v), n1.Dot(a1v), d.Dot(mid)};
        auto det3 = [](const gp_XYZ& r0, const gp_XYZ& r1, const gp_XYZ& r2) {
            return r0.X() * (r1.Y() * r2.Z() - r1.Z() * r2.Y()) -
                   r0.Y() * (r1.X() * r2.Z() - r1.Z() * r2.X()) +
                   r0.Z() * (r1.X() * r2.Y() - r1.Y() * r2.X());
        };
        const double det = det3(rows[0], rows[1], rows[2]);
        if (!(std::abs(det) > 1e-15)) continue;
        auto colSub = [&](int col) {
            gp_XYZ r[3] = {rows[0], rows[1], rows[2]};
            for (int k = 0; k < 3; ++k) r[k].SetCoord(col + 1, rhs[k]);
            return det3(r[0], r[1], r[2]) / det;
        };
        const double x = colSub(0);
        const double y = colSub(1);
        const double z = colSub(2);
        if (!finite3(x, y, z)) continue;
        dirs.push_back(d);
        vg.pts.push_back(gp_XYZ(x, y, z));
    }
    vg.nCand = static_cast<int>(dirs.size());
    vg.tiltMax = tiltMax;
    if (dirs.size() < 3) return vg;

    // "the running axis common to every fold" (D-130-8). The lines find the
    // direction themselves: a generator is parallel to every other generator,
    // a chord (the diagonal of a staggered quad, the ridge onto a neighbouring
    // wall) is parallel to nothing. So take the direction the MOST lines agree
    // with inside tiltMax -- the normal-covariance PCA of a three-facet seed is
    // 1e-2 rad off and cannot be the referee here, which is exactly why a
    // staggered seed found no generator when it was.
    size_t bestI = 0;
    int bestN = -1;
    for (size_t i = 0; i < dirs.size(); ++i) {
        int c = 0;
        for (size_t j = 0; j < dirs.size(); ++j)
            if (dirs[i].Crossed(dirs[j]).Modulus() <= tiltMax) ++c;
        if (c > bestN) {
            bestN = c;
            bestI = i;
        }
    }
    if (bestN < 3) return vg;

    gp_XYZ acc(0, 0, 0);
    std::vector<gp_XYZ> keptPts;
    keptPts.reserve(static_cast<size_t>(bestN));
    for (size_t i = 0; i < dirs.size(); ++i) {
        if (dirs[i].Crossed(dirs[bestI]).Modulus() > tiltMax) continue;
        acc += dirs[i];
        keptPts.push_back(vg.pts[i]);
    }
    const double am = acc.Modulus();
    if (!(am > 1e-15)) return vg;
    acc.Divide(am);
    std::array<double, 3> av{acc.X(), acc.Y(), acc.Z()};
    signNormalize(av);
    vg.dir = gp_XYZ(av[0], av[1], av[2]);

    for (size_t i = 0; i < dirs.size(); ++i) {
        if (dirs[i].Crossed(dirs[bestI]).Modulus() > tiltMax) continue;
        vg.spreadRad = std::max(vg.spreadRad, dirs[i].Crossed(vg.dir).Modulus());
    }
    vg.pts.swap(keptPts);
    vg.nKept = static_cast<int>(vg.pts.size());
    vg.ok = vg.spreadRad <= tiltMax;
    return vg;
}

// The chain of a staggered strip: axis from the facet-plane intersection lines,
// azimuths from every band vertex. Everything downstream of the generator set
// (clusterAngles / makeIntervals / w = 2 R sin(theta/2)) is the same machinery
// extractChain uses -- only the generators differ, which is the whole of
// D-130-8's construction.
bool extractChainVirtual(const MeshView& mv, const std::vector<int>& ids, const DerivedTols& tol,
                         LawBand& out, VGenAxis& vgOut) {
    out = LawBand{};
    out.tris = ids;

    gp_XYZ axisSeed;
    if (!recoverAxisDir(mv, ids, axisSeed)) return false;

    std::vector<gp_XYZ> verts;
    collectUniqueVerts(mv, ids, verts);
    if (verts.size() < 3) return false;

    gp_XYZ centroid(0, 0, 0);
    for (const gp_XYZ& p : verts) centroid += p;
    centroid.Divide(static_cast<double>(verts.size()));

    double aLo = 1e300, aHi = -1e300;
    for (const gp_XYZ& p : verts) {
        const double x = (p - centroid).Dot(axisSeed);
        aLo = std::min(aLo, x);
        aHi = std::max(aHi, x);
    }
    const double span = aHi - aLo;
    if (!(span > 0.0)) return false;

    vgOut = virtualGeneratorAxis(mv, ids, axisSeed, span, linTol(mv, tol), tauQuant(mv));
    if (!vgOut.ok || vgOut.pts.size() < 3) return false;
    const gp_XYZ axis = vgOut.dir;
    const std::vector<gp_XYZ>& gpts = vgOut.pts;

    gp_XYZ u, v;
    if (!axisFrame(axis, u, v)) return false;

    gp_XYZ gcen(0, 0, 0);
    for (const gp_XYZ& p : gpts) gcen += p;
    gcen.Divide(static_cast<double>(gpts.size()));

    std::vector<double> xs, ys;
    xs.reserve(gpts.size());
    ys.reserve(gpts.size());
    for (const gp_XYZ& p : gpts) {
        const gp_XYZ d = p - gcen;
        xs.push_back(d.Dot(u));
        ys.push_back(d.Dot(v));
    }
    double cx = 0, cy = 0;
    if (!lsBisectorCenter(xs, ys, cx, cy)) {
        if (!circumMedianCenter(xs, ys, cx, cy)) return false;
    }
    const gp_XYZ origin = gcen + u * cx + v * cy;

    std::vector<double> angs;
    angs.reserve(gpts.size());
    for (const gp_XYZ& p : gpts) angs.push_back(azimuth(p, origin, axis, u, v));
    std::vector<double> gens;
    clusterAngles(angs, 5e-4, gens);
    if (gens.size() < 3) return false;

    std::vector<StripPair> pairs;
    bool closed = false;
    makeIntervals(gens, pairs, closed);
    if (pairs.empty()) return false;

    std::vector<double> Ri;
    Ri.reserve(pairs.size());
    for (const StripPair& p : pairs) {
        const gp_XYZ p0 = nearestAtAngle(gpts, p.a0, origin, axis, u, v);
        const gp_XYZ p1 = nearestAtAngle(gpts, p.a1, origin, axis, u, v);
        const double rho = 0.5 * (rhoOf(p0, origin, axis) + rhoOf(p1, origin, axis));
        const double s = std::sin(0.5 * p.span);
        if (!(s > 1e-15)) continue;
        const double w = 2.0 * rho * s;
        if (!(w > 0.0)) continue;
        out.theta.push_back(p.span);
        out.w.push_back(w);
        Ri.push_back(w / (2.0 * s));
    }
    out.axis = gp_Ax1(gp_Pnt(origin.X(), origin.Y(), origin.Z()),
                      gp_Dir(axis.X(), axis.Y(), axis.Z()));
    if (out.theta.size() < 2) {
        out.N = static_cast<int>(out.theta.size());
        if (!Ri.empty()) out.R = medianOf(Ri);
        return false;
    }
    out.R = medianOf(Ri);
    out.N = static_cast<int>(out.theta.size());
    out.closed360 = closed;
    out.phi = closed ? kTwoPi : 0.0;
    if (!closed)
        for (double t : out.theta) out.phi += t;
    out.cvTheta = popCV(out.theta);
    out.cvR = popCV(Ri);

    double maxRes = 0.0;
    for (const gp_XYZ& p : verts)
        maxRes = std::max(maxRes, std::abs(rhoOf(p, origin, axis) - out.R));
    out.maxVertResid = maxRes;
    return true;
}

bool extractChain(const MeshView& mv, const std::vector<int>& ids, const DerivedTols& tol,
                  LawBand& out) {
    out = LawBand{};
    out.tris = ids;

    gp_XYZ axis;
    if (!recoverAxisDir(mv, ids, axis)) return false;
    out.axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(axis.X(), axis.Y(), axis.Z()));

    std::vector<gp_XYZ> verts;
    collectUniqueVerts(mv, ids, verts);
    if (verts.size() < 3) return false;

    gp_XYZ centroid(0, 0, 0);
    for (const gp_XYZ& p : verts) centroid += p;
    centroid.Divide(static_cast<double>(verts.size()));

    std::vector<double> axial;
    axial.reserve(verts.size());
    for (const gp_XYZ& p : verts) axial.push_back((p - centroid).Dot(axis));
    const double amin = *std::min_element(axial.begin(), axial.end());
    const double amax = *std::max_element(axial.begin(), axial.end());
    const double span = amax - amin;
    const double modeTol = std::max(1e-4, 1e-3 * std::max(span, 1e-9));
    std::vector<double> modes;
    cluster1d(axial, modeTol, modes);
    const double lo = modes.front();
    const double hi = modes.back();
    const double ringTol = std::max(5e-4, linTol(mv, tol));

    std::vector<gp_XYZ> endPts;
    for (size_t i = 0; i < verts.size(); ++i) {
        if (std::abs(axial[i] - lo) <= ringTol || std::abs(axial[i] - hi) <= ringTol)
            endPts.push_back(verts[i]);
    }
    if (endPts.size() < 3) endPts = verts;

    const double merge = std::max({4.0 * mv.weldTol, 1e-4, 1e-6 * mv.diag});
    std::vector<gp_XYZ> gens3;
    uniquePositions(endPts, merge, gens3);
    if (gens3.size() < 3) return false;

    gp_XYZ u, v;
    if (!axisFrame(axis, u, v)) return false;

    std::vector<double> xs, ys;
    xs.reserve(gens3.size());
    ys.reserve(gens3.size());
    for (const gp_XYZ& p : gens3) {
        const gp_XYZ d = p - centroid;
        xs.push_back(d.Dot(u));
        ys.push_back(d.Dot(v));
    }

    double cx = 0, cy = 0;
    if (!lsBisectorCenter(xs, ys, cx, cy)) {
        if (!circumMedianCenter(xs, ys, cx, cy)) return false;
    }
    gp_XYZ origin = centroid + u * cx + v * cy;

    // Refine axis from end-ring generator pairs (same azimuth, two rings).
    {
        gp_XYZ acc(0, 0, 0);
        int np = 0;
        for (const gp_XYZ& p : gens3) {
            const double az = azimuth(p, origin, axis, u, v);
            const double vax = (p - centroid).Dot(axis);
            for (const gp_XYZ& q : gens3) {
                if ((&q) == (&p)) continue;
                if (std::abs(wrapPi(azimuth(q, origin, axis, u, v) - az)) > 5e-3) continue;
                const double vbx = (q - centroid).Dot(axis);
                if (std::abs(vbx - vax) < 1e-3) continue;
                gp_XYZ d = q - p;
                const double m = d.Modulus();
                if (m < 1e-9) continue;
                d.Divide(m);
                if (d.Dot(axis) < 0.0) d.Reverse();
                acc += d;
                ++np;
            }
        }
        if (np >= 2 && acc.Modulus() > 1e-15) {
            acc.Divide(acc.Modulus());
            std::array<double, 3> av{acc.X(), acc.Y(), acc.Z()};
            signNormalize(av);
            axis = gp_XYZ(av[0], av[1], av[2]);
            if (axisFrame(axis, u, v)) {
                xs.clear();
                ys.clear();
                for (const gp_XYZ& p : gens3) {
                    const gp_XYZ d = p - centroid;
                    xs.push_back(d.Dot(u));
                    ys.push_back(d.Dot(v));
                }
                if (lsBisectorCenter(xs, ys, cx, cy) || circumMedianCenter(xs, ys, cx, cy))
                    origin = centroid + u * cx + v * cy;
            }
        }
    }

    std::vector<double> angs;
    angs.reserve(gens3.size());
    for (const gp_XYZ& p : gens3) angs.push_back(azimuth(p, origin, axis, u, v));
    std::vector<double> gens;
    clusterAngles(angs, 5e-4, gens);
    if (gens.size() < 2) return false;

    std::vector<StripPair> pairs;
    bool closed = false;
    makeIntervals(gens, pairs, closed);
    if (pairs.empty()) return false;

    out.theta.reserve(pairs.size());
    out.w.reserve(pairs.size());
    std::vector<double> Ri;
    Ri.reserve(pairs.size());
    for (const StripPair& p : pairs) {
        const gp_XYZ p0 = nearestAtAngle(gens3, p.a0, origin, axis, u, v);
        const gp_XYZ p1 = nearestAtAngle(gens3, p.a1, origin, axis, u, v);
        const gp_XYZ d = p1 - p0;
        const double w3 = (d - axis * d.Dot(axis)).Modulus();
        const double rho = 0.5 * (rhoOf(p0, origin, axis) + rhoOf(p1, origin, axis));
        const double s = std::sin(0.5 * p.span);
        if (!(s > 1e-15)) continue;
        // Circumferential chord of the two generators about the recovered axis.
        // Polar form is the same chord on the inscribed circle (RULE 4.1d).
        const double w = (rho > 0.0) ? (2.0 * rho * s) : w3;
        if (!(w > 0.0)) continue;
        out.theta.push_back(p.span);
        out.w.push_back(w);
        Ri.push_back(w / (2.0 * s));
    }
    if (out.theta.size() < 2) {
        out.N = static_cast<int>(out.theta.size());
        if (!Ri.empty()) out.R = medianOf(Ri);
        out.axis = gp_Ax1(gp_Pnt(origin.X(), origin.Y(), origin.Z()),
                          gp_Dir(axis.X(), axis.Y(), axis.Z()));
        return false;
    }

    // Equal-θ inverse (RULE 4.1d): R_i from measured (w_i, θ_i).
    out.R = medianOf(Ri);
    out.N = static_cast<int>(out.theta.size());
    out.closed360 = closed;
    out.phi = closed ? kTwoPi : 0.0;
    if (!closed)
        for (double t : out.theta) out.phi += t;
    out.cvTheta = popCV(out.theta);
    out.cvR = popCV(Ri);
    out.axis = gp_Ax1(gp_Pnt(origin.X(), origin.Y(), origin.Z()),
                      gp_Dir(axis.X(), axis.Y(), axis.Z()));

    const double tau = tauSurf(mv);
    double maxRes = 0.0;
    const unsigned nV = static_cast<unsigned>(verts.size());
    unsigned nTh = std::thread::hardware_concurrency();
    if (nTh == 0) nTh = 2;
    if (nV < 48) nTh = 1;
    nTh = std::min(nTh, nV);
    if (nTh <= 1) {
        for (const gp_XYZ& p : verts)
            maxRes = std::max(maxRes, std::abs(rhoOf(p, origin, axis) - out.R));
    } else {
        std::vector<double> part(nTh, 0.0);
        std::vector<std::thread> pool;
        pool.reserve(nTh);
        for (unsigned t = 0; t < nTh; ++t) {
            pool.emplace_back([&, t]() {
                double m = 0.0;
                for (unsigned i = t; i < nV; i += nTh)
                    m = std::max(m, std::abs(rhoOf(verts[i], origin, axis) - out.R));
                part[t] = m;
            });
        }
        for (auto& th : pool) th.join();
        for (double m : part) maxRes = std::max(maxRes, m);
    }
    out.maxVertResid = maxRes;
    (void)tau;
    return true;
}

bool testEqualTheta(const LawBand& b) { return b.cvTheta < kCvThetaMax; }

bool testRcons(const LawBand& b) {
    if (!(b.R > 0.0) || b.theta.size() != b.w.size() || b.theta.empty()) return false;
    double mx = 0.0;
    std::vector<double> Ri;
    Ri.reserve(b.theta.size());
    for (size_t i = 0; i < b.theta.size(); ++i) {
        const double s = std::sin(0.5 * b.theta[i]);
        if (!(s > 1e-15)) return false;
        const double r = b.w[i] / (2.0 * s);
        Ri.push_back(r);
        mx = std::max(mx, std::abs(r - b.R) / b.R);
    }
    // CV(R)<1e-6 is the tautological (w=2 R sin(θ/2)) measurement. On mesh
    // chords, equal-θ (test 1) already implies one R; keep the 5e-4 envelope.
    const double cv = popCV(Ri);
    return mx < kRelRMax && (cv < kCvRMax || b.cvTheta < kCvThetaMax);
}

bool testCommonAxis(const MeshView& mv, const std::vector<int>& ids, const LawBand& b,
                    const DerivedTols& tol) {
    if (!mv.compEdges || !mv.triEdges || !mv.compVtx || !mv.pts) return true;
    const gp_XYZ a = b.axis.Direction().XYZ();
    const gp_XYZ o = b.axis.Location().XYZ();
    const double lim = linTol(mv, tol);
    int nGen = 0;
    int nOk = 0;
    std::vector<char> in(mv.nTri, 0);
    for (int t : ids)
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) in[static_cast<size_t>(t)] = 1;
    for (int t : ids) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        for (int s = 0; s < 3; ++s) {
            const int e = mv.triEdges[t][s];
            const int v0 = mv.compEdges[e].first;
            const int v1 = mv.compEdges[e].second;
            const gp_XYZ p0 = mv.pts[mv.compVtx[v0]];
            const gp_XYZ p1 = mv.pts[mv.compVtx[v1]];
            const gp_XYZ d = p1 - p0;
            const double len = d.Modulus();
            if (!(len > 1e-12)) continue;
            const gp_XYZ ehat = d / len;
            const double par = std::abs(ehat.Dot(a));
            // generator: nearly parallel to axis (scale-relative tilt, no degree gate)
            const double tilt = ehat.Crossed(a).Modulus();
            if (tilt * len > lim * 8.0) continue;
            if (par < 0.98) continue;
            ++nGen;
            const double midRho = rhoOf((p0 + p1) * 0.5, o, a);
            if (std::abs(midRho - b.R) < std::max(lim, 5e-3 * std::max(b.R, 1.0))) ++nOk;
        }
    }
    if (nGen < 2) return true;  // no generator edges resolved; axis from PCA stands
    return nOk * 2 >= nGen;     // majority of generators coaxial
}

bool testOnSurface(const LawBand& b, const MeshView& mv) {
    return b.maxVertResid < tauSurf(mv);
}

// D-130-12: the certificate a VIRTUAL-generator band is judged on is the mesh
// FILE's own coordinate quantization -- the radius of the ball a vertex may sit
// anywhere inside purely from having been written to this file. It is 13x
// tighter than tauSurf's part-scale budget on the plate (1.32e-05 against
// 1.795e-04), never looser; when the file could not answer (quantFloor == 0) it
// certifies nothing and the part-scale budget stands.
double tauQuant(const MeshView& mv) {
    if (std::isfinite(mv.quantFloor) && mv.quantFloor > 0.0) return mv.quantFloor;
    return tauSurf(mv);
}

// D-130-8: equal-theta is how a REGULAR tessellation is recognised, not what
// certifies a cylinder -- a staggered strip whose generators the mesh drew at
// only some azimuths reads theta and 2*theta and is no less a cylinder for it.
// What must hold at EVERY generator is the chain law itself: w = 2 R sin(t/2)
// with one R. So the disjunction testRcons allows (cvR OR equal-theta) is not
// available here; the band must carry the law on its own.
// A cylinder claim must be one the mesh can tell from a PLANE. If the band's
// vertices lie on one plane within the same certificate that admitted them to
// the cylinder, the mesh resolves no curvature there and no radius fitted to
// them is certified -- the plate's flat slot side walls otherwise read as
// R = 69.696 from four facets whose deviation from their own plane is zero.
// Parameter-free: the same tau on both sides, plane residual against cylinder
// residual, nothing angular.
double maxPlaneResid(const MeshView& mv, const std::vector<int>& ids) {
    std::vector<gp_XYZ> verts;
    collectUniqueVerts(mv, ids, verts);
    if (verts.size() < 3) return 0.0;
    gp_XYZ c(0, 0, 0);
    for (const gp_XYZ& p : verts) c += p;
    c.Divide(static_cast<double>(verts.size()));
    std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (const gp_XYZ& p : verts) {
        const double dx = p.X() - c.X(), dy = p.Y() - c.Y(), dz = p.Z() - c.Z();
        cov[0][0] += dx * dx;
        cov[0][1] += dx * dy;
        cov[0][2] += dx * dz;
        cov[1][1] += dy * dy;
        cov[1][2] += dy * dz;
        cov[2][2] += dz * dz;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return 0.0;
    gp_XYZ n(evec[0][0], evec[0][1], evec[0][2]);
    const double nm = n.Modulus();
    if (nm < 1e-15) return 0.0;
    n.Divide(nm);
    double mx = 0.0;
    for (const gp_XYZ& p : verts) mx = std::max(mx, std::abs((p - c).Dot(n)));
    return mx;
}

// D-130-8: equal-theta generalised to the lattice a tessellation actually has.
// A tessellated cylinder carries ONE generator pitch; where the mesh drew no
// ridge the chain reads an exact INTEGER MULTIPLE of it (the plate's R=3 slot
// ends read 12.857 deg and 25.714 deg, both multiples of the wall's 6.4286 deg
// 56-gon pitch). So the test is not that the thetas are equal but that their
// reduced pitches are, inside the same kCvThetaMax envelope testEqualTheta
// uses. A symmetric bevel -- a flat chord between two equal folds -- fits a
// cylinder exactly and is refused here, because its two pitches are not
// commensurate.
bool testLatticeTheta(const LawBand& b) {
    if (b.theta.size() < 2) return false;
    double base = b.theta.front();
    for (double t : b.theta) base = std::min(base, t);
    if (!(base > 0.0)) return false;
    std::vector<double> pitch;
    pitch.reserve(b.theta.size());
    for (double t : b.theta) {
        const double k = std::floor(t / base + 0.5);
        if (!(k >= 1.0) || !std::isfinite(k)) return false;
        pitch.push_back(t / k);
    }
    return popCV(pitch) < kCvThetaMax;
}

bool testRconsStrict(const LawBand& b) {
    if (!(b.R > 0.0) || b.theta.size() != b.w.size() || b.theta.empty()) return false;
    double mx = 0.0;
    std::vector<double> Ri;
    Ri.reserve(b.theta.size());
    for (size_t i = 0; i < b.theta.size(); ++i) {
        const double s = std::sin(0.5 * b.theta[i]);
        if (!(s > 1e-15)) return false;
        const double r = b.w[i] / (2.0 * s);
        Ri.push_back(r);
        mx = std::max(mx, std::abs(r - b.R) / b.R);
    }
    return mx < kRelRMax && popCV(Ri) < kCvRMax;
}

bool capCircleOk(const MeshView& mv, const std::vector<int>& ids, const LawBand& b) {
    std::vector<gp_XYZ> verts;
    collectUniqueVerts(mv, ids, verts);
    if (verts.size() < 3 || !(b.R > 0.0)) return false;
    const gp_XYZ o = b.axis.Location().XYZ();
    const gp_XYZ a = b.axis.Direction().XYZ();
    std::vector<double> rhos;
    rhos.reserve(verts.size());
    for (const gp_XYZ& p : verts) rhos.push_back(rhoOf(p, o, a));
    const double med = medianOf(rhos);
    if (!(med > 0.0)) return false;
    double mx = 0.0;
    for (double r : rhos) mx = std::max(mx, std::abs(r - b.R) / b.R);
    return mx < kRelRMax && popCV(rhos) < 1e-3;
}

void emitLawband(const LawBand& b, bool accept) {
    if (!diagOn()) return;
    int rid = b.tris.empty() ? -1 : b.tris[0];
    for (int t : b.tris)
        if (t < rid) rid = t;
    std::lock_guard<std::mutex> g(diagMu());
    std::fprintf(stderr,
                 "DIAG_LAWBAND rid=%d n=%d N=%d R=%.6f phi=%.6f cvT=%.3g cvR=%.3g resid=%.3g "
                 "accept=%d\n",
                 rid, static_cast<int>(b.tris.size()), b.N, b.R, b.phi, b.cvTheta, b.cvR,
                 b.maxVertResid, accept ? 1 : 0);
}

void emitLawcal(const TessLawInterval& li) {
    if (!diagOn()) return;
    const double rad2deg = 180.0 / kPi;
    std::lock_guard<std::mutex> g(diagMu());
    std::fprintf(stderr, "DIAG_LAWCAL dLo=%.6f dHi=%.6f aLoDeg=%.4f aHiDeg=%.4f nD=%d nA=%d empty=%d\n",
                 li.dLo, li.dHi, li.alphaLo * rad2deg, li.alphaHi * rad2deg, li.nDLimited,
                 li.nAlphaLimited, li.empty ? 1 : 0);
}

struct DInt {
    double lo = 0, hi = 0;
    bool ok = false;
};

DInt dFeasible(const LawBand& b) {
    DInt o;
    if (b.N < 2 || !(b.R > 0.0) || !(b.phi > 0.0)) return o;
    const double thLo = b.phi / static_cast<double>(b.N);
    const double thHi = b.phi / static_cast<double>(b.N - 1);
    if (!(thLo > 0.0) || thLo >= kPi) return o;
    auto dFrom = [&](double th) { return b.R * (1.0 - std::cos(0.5 * th)); };
    o.lo = dFrom(thLo);
    o.hi = (thHi >= kPi) ? (2.0 * b.R) : dFrom(thHi);
    o.ok = o.hi > o.lo && std::isfinite(o.lo) && std::isfinite(o.hi);
    return o;
}

double thetaSurf(double R, double d) {
    if (!(R > 0.0) || !(d > 0.0) || d >= 2.0 * R) return 0.0;
    const double x = 1.0 - d / R;
    if (x <= -1.0) return kPi;
    if (x >= 1.0) return 0.0;
    return 2.0 * std::acos(x);
}

}  // namespace

bool lawChainAccept(const MeshView& mv, const std::vector<int>& tris, const DerivedTols& tol,
                    LawBand& out) {
    out = LawBand{};
    std::vector<int> ids;
    sortedUnique(tris, ids);
    if (ids.size() < 2) {
        emitLawband(out, false);
        return false;
    }
    if (!extractChain(mv, ids, tol, out)) {
        emitLawband(out, false);
        return false;
    }

    const bool t1 = testEqualTheta(out);
    const bool t2 = testRcons(out);
    const bool t3 = testCommonAxis(mv, ids, out, tol);
    const bool t4 = testOnSurface(out, mv);

    bool accept = false;
    out.lowConfidence = false;
    if (out.N <= 1) {
        accept = false;
    } else if (out.N == 2) {
        // RULE 4.2d: N==2 needs tests 2–4 + cap-circle. A 4-tri plane/cyl
        // chimera (rid=5) can look like two strips — require a real chain.
        const bool cap = capCircleOk(mv, ids, out);
        accept = t1 && t2 && t3 && t4 && cap && ids.size() >= 6;
        out.lowConfidence = accept;
    } else {
        accept = t1 && t2 && t3 && t4;
    }
    // D-130-8 -- virtual generators, tried ONLY where the end-ring pairing did
    // not certify a band. A set the mesh-edge chain accepts is left exactly as
    // it was, band and all; only what stage L refuses today can reach here.
    if (!accept) {
        LawBand vb;
        VGenAxis vg;
        const bool got = extractChainVirtual(mv, ids, tol, vb, vg);
        const char* why = "ok";
        if (!got) why = vg.nLines == 0 ? "nolines" : (!vg.ok ? "spread" : "chain");
        else if (vb.N < 2) why = "shortN";
        else if (!testLatticeTheta(vb)) why = "lattice";
        else if (!testRconsStrict(vb)) why = "rcons";
        else if (!testCommonAxis(mv, ids, vb, tol)) why = "axis";
        else if (!(vb.maxVertResid < tauQuant(mv))) why = "resid";
        else if (!(maxPlaneResid(mv, ids) > tauQuant(mv))) why = "planar";
        else if (vb.N == 2 && !(capCircleOk(mv, ids, vb) && ids.size() >= 6)) why = "cap";
        if (vgenDiagOn()) {
            std::lock_guard<std::mutex> g(diagMu());
            std::fprintf(stderr,
                         "DIAG_LAWVGENX rid=%d n=%zu why=%s lines=%d cand=%d kept=%d "
                         "tiltMax=%.3g spread=%.3g N=%d R=%.9f cvT=%.3g cvR=%.3g resid=%.3g "
                         "tauQ=%.3g\n",
                         ids.empty() ? -1 : ids.front(), ids.size(), why, vg.nLines, vg.nCand,
                         vg.nKept, vg.tiltMax, vg.spreadRad, vb.N, vb.R, vb.cvTheta, vb.cvR,
                         vb.maxVertResid, tauQuant(mv));
        }
        if (got && std::strcmp(why, "ok") == 0) {
            {
                if (diagOn()) {
                    std::lock_guard<std::mutex> g(diagMu());
                    std::fprintf(stderr,
                                 "DIAG_LAWVGEN rid=%d n=%zu lines=%d kept=%d spread=%.3g "
                                 "N=%d R=%.9f cvT=%.3g cvR=%.3g resid=%.3g tauQ=%.3g "
                                 "closed=%d\n",
                                 ids.empty() ? -1 : ids.front(), ids.size(), vg.nLines,
                                 vg.nKept, vg.spreadRad, vb.N, vb.R, vb.cvTheta, vb.cvR,
                                 vb.maxVertResid, tauQuant(mv), vb.closed360 ? 1 : 0);
                }
                vb.lowConfidence = (vb.N == 2);
                vb.virtualGen = true;
                out = std::move(vb);
                accept = true;
            }
        }
    }
    emitLawband(out, accept);
    return accept;
}

// The folds the mesh RESOLVES in a triangle set (a fold whose far corner
// stands off its neighbour's plane by more than the file's own quantization),
// and how many of them are mutually parallel. These are the two numbers a seed
// bootstrap needs: a virtual-generator chain cannot exist below three parallel
// folds, and a set that gains no new fold when it grows will never reach them.
void lawVirtualFoldCount(const MeshView& mv, const std::vector<int>& tris,
                         const DerivedTols& tol, int& nFolds, int& nParallel) {
    nFolds = 0;
    nParallel = 0;
    std::vector<int> ids;
    sortedUnique(tris, ids);
    if (ids.size() < 2) return;
    gp_XYZ axisSeed;
    if (!recoverAxisDir(mv, ids, axisSeed)) return;
    std::vector<gp_XYZ> verts;
    collectUniqueVerts(mv, ids, verts);
    if (verts.size() < 3) return;
    gp_XYZ c(0, 0, 0);
    for (const gp_XYZ& p : verts) c += p;
    c.Divide(static_cast<double>(verts.size()));
    double lo = 1e300, hi = -1e300;
    for (const gp_XYZ& p : verts) {
        const double x = (p - c).Dot(axisSeed);
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    if (!(hi > lo)) return;
    const VGenAxis vg =
        virtualGeneratorAxis(mv, ids, axisSeed, hi - lo, linTol(mv, tol), tauQuant(mv));
    nFolds = vg.nLines;
    nParallel = vg.nKept;
}

TessLawInterval lawCalibrate(const std::vector<LawBand>& bands) {
    TessLawInterval li;
    std::vector<DInt> iv;
    std::vector<int> idx;
    iv.reserve(bands.size());
    for (size_t i = 0; i < bands.size(); ++i) {
        DInt d = dFeasible(bands[i]);
        iv.push_back(d);
        if (d.ok) idx.push_back(static_cast<int>(i));
    }
    if (idx.empty()) {
        li.empty = true;
        emitLawcal(li);
        return li;
    }

    // Maximum subset with nonempty common intersection = max interval stabbing.
    struct Ev {
        double x;
        int s;
        int id;
    };
    std::vector<Ev> ev;
    for (int i : idx) {
        ev.push_back({iv[static_cast<size_t>(i)].lo, +1, i});
        ev.push_back({iv[static_cast<size_t>(i)].hi, -1, i});
    }
    std::sort(ev.begin(), ev.end(), [](const Ev& a, const Ev& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.s > b.s;
    });
    int cov = 0, best = 0;
    double bestX = ev.front().x;
    for (const Ev& e : ev) {
        cov += e.s;
        if (cov > best) {
            best = cov;
            bestX = e.x;
        }
    }
    std::vector<int> dset, leftover;
    for (int i : idx) {
        if (iv[static_cast<size_t>(i)].lo <= bestX && bestX < iv[static_cast<size_t>(i)].hi)
            dset.push_back(i);
        else
            leftover.push_back(i);
    }
    if (dset.empty()) {
        li.empty = true;
        emitLawcal(li);
        return li;
    }

    double dLo = iv[static_cast<size_t>(dset[0])].lo;
    double dHi = iv[static_cast<size_t>(dset[0])].hi;
    for (int i : dset) {
        dLo = std::max(dLo, iv[static_cast<size_t>(i)].lo);
        dHi = std::min(dHi, iv[static_cast<size_t>(i)].hi);
    }

    const double dMid = 0.5 * (dLo + dHi);
    std::vector<int> aset;
    bool mixed = false;
    for (int i : leftover) {
        const LawBand& b = bands[static_cast<size_t>(i)];
        const double thEq = (b.N > 0) ? (b.phi / static_cast<double>(b.N)) : 0.0;
        const double ths = thetaSurf(b.R, dMid);
        if (ths > thEq && b.N >= 2)
            aset.push_back(i);
        else
            mixed = true;
    }
    if (mixed) {
        li.empty = true;
        li.nDLimited = static_cast<int>(dset.size());
        li.nAlphaLimited = static_cast<int>(aset.size());
        emitLawcal(li);
        return li;
    }

    li.dLo = dLo;
    li.dHi = dHi;
    li.nDLimited = static_cast<int>(dset.size());
    li.nAlphaLimited = static_cast<int>(aset.size());
    li.empty = !(dHi > dLo);

    if (!aset.empty()) {
        double aLo = 0.0, aHi = 1e300;
        bool first = true;
        for (int i : aset) {
            const LawBand& b = bands[static_cast<size_t>(i)];
            if (b.N < 2) continue;
            const double lo = b.phi / static_cast<double>(b.N);
            const double hi = b.phi / static_cast<double>(b.N - 1);
            if (first) {
                aLo = lo;
                aHi = hi;
                first = false;
            } else {
                aLo = std::max(aLo, lo);
                aHi = std::min(aHi, hi);
            }
        }
        if (!first && aHi > aLo) {
            li.alphaLo = aLo;
            li.alphaHi = aHi;
        }
    }
    emitLawcal(li);
    return li;
}

bool lawNConsistent(const LawBand& b, const TessLawInterval& li) {
    if (li.empty || !(li.dHi > li.dLo)) return false;
    if (b.N < 2 || !(b.R > 0.0) || !(b.phi > 0.0)) return false;

    auto thMax = [&](double d) {
        double th = thetaSurf(b.R, d);
        if (li.alphaHi > li.alphaLo) th = std::min(th, li.alphaLo);
        return th;
    };
    auto nFrom = [&](double th) -> int {
        if (!(th > 0.0)) return 0;
        return static_cast<int>(std::ceil(b.phi / th - 1e-15));
    };

    const int nA = nFrom(thMax(li.dLo + (li.dHi - li.dLo) * 1e-9));
    const int nB = nFrom(thMax(li.dHi));
    const int nLo = std::min(nA, nB);
    const int nHi = std::max(nA, nB);
    if (b.N >= nLo && b.N <= nHi) return true;

    // Invert: d-range that yields this N, then overlap.
    const double thLo = b.phi / static_cast<double>(b.N);
    const double thHi = b.phi / static_cast<double>(b.N - 1);
    auto dFrom = [&](double th) {
        if (!(b.R > 0.0) || th <= 0.0) return 0.0;
        return b.R * (1.0 - std::cos(0.5 * std::min(th, kPi)));
    };
    double lo = dFrom(thLo);
    double hi = dFrom(thHi);
    if (hi < lo) std::swap(lo, hi);
    return std::max(lo, li.dLo) < std::min(hi, li.dHi);
}

bool lawBandsMergeable(const LawBand& a, const LawBand& b, const DerivedTols& tol) {
    const bool shortFrag = (a.tris.size() <= 4 || b.tris.size() <= 4);
    const bool haveR = (a.R > 0.0 && b.R > 0.0);
    const double rAvg = haveR ? 0.5 * (a.R + b.R) : std::max(a.R, b.R);
    if (haveR && !shortFrag && std::abs(a.R - b.R) / std::max(rAvg, 1e-12) >= kRelRMax)
        return false;
    if (!haveR && !shortFrag && a.R <= 0.0 && b.R <= 0.0) return false;

    const double cdir = std::abs(a.axis.Direction().Dot(b.axis.Direction()));
    if (!shortFrag && cdir < 0.9999995) return false;  // ~1e-6 rad
    if (!shortFrag) {
        const double sep = axisLineSep(a.axis, b.axis);
        const double lim = std::max({tol.epsMesh, tol.epsPlane, 4.0 * kTauSurfFloor});
        if (sep > std::max(lim, kRelRMax * std::max(rAvg, 1.0))) return false;
    }

    // Shared generator proxy: triangle-id ranges overlap or abut (same-face split).
    if (a.tris.empty() || b.tris.empty()) return false;
    int a0 = a.tris[0], a1 = a.tris[0], b0 = b.tris[0], b1 = b.tris[0];
    for (int t : a.tris) {
        a0 = std::min(a0, t);
        a1 = std::max(a1, t);
    }
    for (int t : b.tris) {
        b0 = std::min(b0, t);
        b1 = std::max(b1, t);
    }
    const bool shareTri = !(a1 < b0 - 1 || b1 < a0 - 1);
    if (!shareTri) {
        for (int t : a.tris)
            for (int u : b.tris)
                if (t == u) goto shared;
        return false;
    }
shared:
    // Short leftover (N<=1 / nTri<=4) cannot carry a stable R; RULE 4.1b still
    // merges it into the abutting equal-θ host.
    const size_t na = a.tris.size(), nb = b.tris.size();
    const LawBand* host = (na >= nb) ? &a : &b;
    const LawBand* frag = (na >= nb) ? &b : &a;
    if (frag->tris.size() <= 4 && host->N >= 3 && host->R > 0.0) return true;
    std::vector<double> cat = a.theta;
    cat.insert(cat.end(), b.theta.begin(), b.theta.end());
    return popCV(cat) < kCvThetaMax;
}

}  // namespace refit
}  // namespace stl2step
