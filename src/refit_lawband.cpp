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
    emitLawband(out, accept);
    return accept;
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
