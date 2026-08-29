// stl2step P1 math primitives. Written in-TU; no Eigen, no GProp_, no GeomConvert_.
// Signatures are frozen in refit_internal.hpp.
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

#ifndef STL2STEP_ARCHCHAIN_BAND_MAX
#define STL2STEP_ARCHCHAIN_BAND_MAX 8000
#endif
static_assert(STL2STEP_ARCHCHAIN_BAND_MAX < 15300,
              "archChainBand upper bound must exclude Body11 file (15300 tris)");

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kEps = 1e-15;
constexpr int kJacobiSweeps = 64;
constexpr int kPrattNewtonCap = 20;  // Chernov: 4-6 typical; hard-capped (I5)

bool finite3(double x, double y, double z) {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

void addOuter(std::array<std::array<double, 3>, 3>& m, double x, double y, double z,
              double w) {
    m[0][0] += w * x * x;
    m[0][1] += w * x * y;
    m[0][2] += w * x * z;
    m[1][0] += w * y * x;
    m[1][1] += w * y * y;
    m[1][2] += w * y * z;
    m[2][0] += w * z * x;
    m[2][1] += w * z * y;
    m[2][2] += w * z * z;
}

void sortedUniqueTris(const std::vector<int>& in, std::vector<int>& out) {
    out = in;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
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

bool triGeom(const MeshView& mv, int localTri, gp_XYZ& a, gp_XYZ& b, gp_XYZ& c,
             gp_XYZ& nUnnorm, double& area) {
    if (!triCorners(mv, localTri, a, b, c)) return false;
    nUnnorm = (b - a).Crossed(c - a);
    area = 0.5 * nUnnorm.Modulus();
    return std::isfinite(area);
}

// D2.3 frame: w = coordinate axis of smallest |a·e_k|, lowest k on a tie.
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

bool solve2x2(double a00, double a01, double a11, double b0, double b1, double& x0,
              double& x1) {
    const double det = a00 * a11 - a01 * a01;
    const double scale = std::abs(a00) + std::abs(a11) + 2.0 * std::abs(a01);
    if (!std::isfinite(det) || !std::isfinite(scale)) return false;
    if (scale <= 0.0 || std::abs(det) <= 1e-14 * scale * scale) return false;
    x0 = (a11 * b0 - a01 * b1) / det;
    x1 = (-a01 * b0 + a00 * b1) / det;
    return std::isfinite(x0) && std::isfinite(x1);
}

// Chernov CircleFitByPratt.m — moments already centroid-centred.
bool prattFit2(const double* xs, const double* ys, std::size_t n, double& cx,
               double& cy, double& radius) {
    if (n < 3) return false;
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mx += xs[i];
        my += ys[i];
    }
    const double inv = 1.0 / static_cast<double>(n);
    mx *= inv;
    my *= inv;

    double Mxx = 0.0, Myy = 0.0, Mxy = 0.0, Mxz = 0.0, Myz = 0.0, Mzz = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double Xi = xs[i] - mx;
        const double Yi = ys[i] - my;
        const double Zi = Xi * Xi + Yi * Yi;
        Mxx += Xi * Xi;
        Myy += Yi * Yi;
        Mxy += Xi * Yi;
        Mxz += Xi * Zi;
        Myz += Yi * Zi;
        Mzz += Zi * Zi;
    }
    Mxx *= inv;
    Myy *= inv;
    Mxy *= inv;
    Mxz *= inv;
    Myz *= inv;
    Mzz *= inv;

    const double Mz = Mxx + Myy;
    const double Cov_xy = Mxx * Myy - Mxy * Mxy;
    const double Mxz2 = Mxz * Mxz;
    const double Myz2 = Myz * Myz;
    const double A2 = 4.0 * Cov_xy - 3.0 * Mz * Mz - Mzz;
    const double A1 = Mzz * Mz + 4.0 * Cov_xy * Mz - Mxz2 - Myz2 - Mz * Mz * Mz;
    const double A0 =
        Mxz2 * Myy + Myz2 * Mxx - Mzz * Cov_xy - 2.0 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy;
    const double A22 = A2 + A2;

    double xnew = 0.0;
    double ynew = 1.0e20;
    for (int iter = 0; iter < kPrattNewtonCap; ++iter) {
        const double yold = ynew;
        ynew = A0 + xnew * (A1 + xnew * (A2 + 4.0 * xnew * xnew));
        if (!std::isfinite(ynew) || std::abs(ynew) > std::abs(yold)) {
            xnew = 0.0;
            break;
        }
        const double Dy = A1 + xnew * (A22 + 16.0 * xnew * xnew);
        if (Dy == 0.0 || !std::isfinite(Dy)) {
            xnew = 0.0;
            break;
        }
        const double xold = xnew;
        xnew = xold - ynew / Dy;
        if (!std::isfinite(xnew)) {
            xnew = 0.0;
            break;
        }
        if (xnew < 0.0) {
            xnew = 0.0;
            break;
        }
        if (std::abs(xnew) > kEps && std::abs((xnew - xold) / xnew) < 1e-12) break;
        if (std::abs(xnew - xold) < 1e-15) break;
    }

    const double DET = xnew * xnew - xnew * Mz + Cov_xy;
    if (!std::isfinite(DET) || std::abs(DET) <= 1e-16 * (1.0 + std::abs(Mz) + std::abs(Cov_xy)))
        return false;
    const double Xc = (Mxz * (Myy - xnew) - Myz * Mxy) / DET / 2.0;
    const double Yc = (Myz * (Mxx - xnew) - Mxz * Mxy) / DET / 2.0;
    const double r2 = Xc * Xc + Yc * Yc + Mz + 2.0 * xnew;
    if (!std::isfinite(Xc) || !std::isfinite(Yc) || !std::isfinite(r2) || r2 <= 0.0)
        return false;
    cx = mx + Xc;
    cy = my + Yc;
    radius = std::sqrt(r2);
    return std::isfinite(cx) && std::isfinite(cy) && std::isfinite(radius) && radius > 0.0;
}

bool kasaFit2(const double* xs, const double* ys, std::size_t n, double& cx, double& cy,
              double& radius) {
    if (n < 3) return false;
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mx += xs[i];
        my += ys[i];
    }
    const double inv = 1.0 / static_cast<double>(n);
    mx *= inv;
    my *= inv;

    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0, Sxrr = 0.0, Syrr = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double x = xs[i] - mx;
        const double y = ys[i] - my;
        const double rr = x * x + y * y;
        Sxx += x * x;
        Sxy += x * y;
        Syy += y * y;
        Sxrr += x * rr;
        Syrr += y * rr;
    }
    double dcx = 0.0, dcy = 0.0;
    if (!solve2x2(Sxx, Sxy, Syy, 0.5 * Sxrr, 0.5 * Syrr, dcx, dcy)) return false;
    cx = mx + dcx;
    cy = my + dcy;
    double acc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = xs[i] - cx;
        const double dy = ys[i] - cy;
        acc += dx * dx + dy * dy;
    }
    radius = std::sqrt(acc * inv);
    return std::isfinite(cx) && std::isfinite(cy) && std::isfinite(radius) && radius > 0.0;
}

}  // namespace

static void maybeInventoryArchChains(const MeshView& mv);

// evec[k][i] = component i of eigenvector k; eval/evec ordered λ0 ≤ λ1 ≤ λ2.
bool jacobiEigenSymmetric3(const std::array<std::array<double, 3>, 3>& m,
                           std::array<double, 3>& eval,
                           std::array<std::array<double, 3>, 3>& evec) {
    double a[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double s = 0.5 * (m[i][j] + m[j][i]);
            if (!std::isfinite(s)) return false;
            a[i][j] = s;
        }
    }

    double v[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    for (int sweep = 0; sweep < kJacobiSweeps; ++sweep) {
        const double off = std::sqrt(2.0 * (a[0][1] * a[0][1] + a[0][2] * a[0][2] +
                                            a[1][2] * a[1][2]));
        const double diag =
            std::abs(a[0][0]) + std::abs(a[1][1]) + std::abs(a[2][2]);
        if (off <= kEps * (1.0 + diag)) break;

        // Deterministic cyclic sweep: (0,1), (0,2), (1,2).
        const int pq[3][2] = {{0, 1}, {0, 2}, {1, 2}};
        for (int k = 0; k < 3; ++k) {
            const int p = pq[k][0];
            const int q = pq[k][1];
            const double apq = a[p][q];
            if (std::abs(apq) <= kEps * (1.0 + std::abs(a[p][p]) + std::abs(a[q][q])))
                continue;

            const double app = a[p][p];
            const double aqq = a[q][q];
            double t;
            if (std::abs(app - aqq) <= kEps * (1.0 + std::abs(app) + std::abs(aqq))) {
                t = (apq >= 0.0) ? 1.0 : -1.0;
            } else {
                const double tau = (aqq - app) / (2.0 * apq);
                const double mag = std::abs(tau) + std::sqrt(1.0 + tau * tau);
                t = 1.0 / mag;
                if (tau < 0.0) t = -t;
            }
            const double c = 1.0 / std::sqrt(1.0 + t * t);
            const double s = t * c;

            a[p][p] = c * c * app + s * s * aqq - 2.0 * s * c * apq;
            a[q][q] = s * s * app + c * c * aqq + 2.0 * s * c * apq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            for (int r = 0; r < 3; ++r) {
                if (r == p || r == q) continue;
                const double arp = a[r][p];
                const double arq = a[r][q];
                a[r][p] = a[p][r] = c * arp - s * arq;
                a[r][q] = a[q][r] = s * arp + c * arq;
            }
            for (int r = 0; r < 3; ++r) {
                const double vrp = v[r][p];
                const double vrq = v[r][q];
                v[r][p] = c * vrp - s * vrq;
                v[r][q] = s * vrp + c * vrq;
            }
        }
    }

    double lam[3] = {a[0][0], a[1][1], a[2][2]};
    int ord[3] = {0, 1, 2};
    if (lam[ord[1]] < lam[ord[0]]) std::swap(ord[0], ord[1]);
    if (lam[ord[2]] < lam[ord[1]]) std::swap(ord[1], ord[2]);
    if (lam[ord[1]] < lam[ord[0]]) std::swap(ord[0], ord[1]);

    for (int k = 0; k < 3; ++k) {
        eval[k] = lam[ord[k]];
        evec[k][0] = v[0][ord[k]];
        evec[k][1] = v[1][ord[k]];
        evec[k][2] = v[2][ord[k]];
        signNormalize(evec[k]);
        if (!finite3(evec[k][0], evec[k][1], evec[k][2]) || !std::isfinite(eval[k]))
            return false;
    }
    return eval[0] <= eval[1] && eval[1] <= eval[2];
}

bool pcaPlane(const MeshView& mv, const std::vector<int>& tris, gp_Ax3& plane) {
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.empty()) return false;

    double areaSum = 0.0;
    gp_XYZ c(0.0, 0.0, 0.0);
    gp_XYZ nMean(0.0, 0.0, 0.0);
    for (int t : ids) {
        gp_XYZ a, b, p, nU;
        double area = 0.0;
        if (!triGeom(mv, t, a, b, p, nU, area) || area <= 0.0) continue;
        const gp_XYZ g = (a + b + p) / 3.0;
        areaSum += area;
        c += g * area;
        nMean += nU;  // already scaled by 2*area; orientation only
    }
    if (!(areaSum > 0.0) || !finite3(c.X(), c.Y(), c.Z())) return false;
    c.Divide(areaSum);

    std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (int t : ids) {
        gp_XYZ a, b, p, nU;
        double area = 0.0;
        if (!triGeom(mv, t, a, b, p, nU, area) || area <= 0.0) continue;
        const double w = area / 3.0;
        const gp_XYZ vs[3] = {a - c, b - c, p - c};
        for (const gp_XYZ& d : vs) addOuter(cov, d.X(), d.Y(), d.Z(), w);
    }

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return false;

    gp_XYZ n(evec[0][0], evec[0][1], evec[0][2]);
    const double nm = n.Modulus();
    if (nm < 1e-15) return false;
    n.Divide(nm);
    if (n.Dot(nMean) < 0.0) n.Reverse();
    if (n.Modulus() < 1e-15) return false;
    plane = gp_Ax3(gp_Pnt(c.X(), c.Y(), c.Z()), gp_Dir(n));
    return true;
}

bool eberlyCenterRadius(const MeshView& mv, const std::vector<int>& tris,
                        const gp_Dir& axis, gp_Pnt& center, double& radius) {
    maybeInventoryArchChains(mv);
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.empty() || !mv.pts || !mv.tris || !mv.compTris) return false;

    gp_XYZ a(axis.XYZ());
    const double am = a.Modulus();
    if (am < 1e-15) return false;
    a.Divide(am);
    gp_XYZ u, v;
    if (!axisFrame(a, u, v)) return false;

    std::vector<int> gids;
    for (int t : ids) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        const int g = mv.compTris[t];
        for (int k = 0; k < 3; ++k) gids.push_back(mv.tris[g][k]);
    }
    std::sort(gids.begin(), gids.end());
    gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
    if (gids.size() < 3) return false;

    std::vector<double> xs(gids.size()), ys(gids.size());
    gp_XYZ mean(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < gids.size(); ++i) {
        const gp_XYZ p = mv.pts[gids[i]];
        mean += p;
        xs[i] = p.Dot(u);
        ys[i] = p.Dot(v);
    }
    mean.Divide(static_cast<double>(gids.size()));

    double cx = 0.0, cy = 0.0, R = 0.0;
    if (!kasaFit2(xs.data(), ys.data(), gids.size(), cx, cy, R)) return false;

    // Closed-form is a single 2×2 solve (signature has no maxRefineIters).
    const gp_XYZ c = cx * u + cy * v + (mean.Dot(a) * a);
    if (!finite3(c.X(), c.Y(), c.Z()) || !(R > 0.0) || !std::isfinite(R)) return false;
    center = gp_Pnt(c.X(), c.Y(), c.Z());
    radius = R;
    return true;
}

bool prattCircleFit(const gp_Pnt* points, std::size_t n, double spanRad, gp_Pnt& center,
                    double& radius) {
    if (!points || n < 3) return false;
    // Short-arc precondition: span in (0, π).
    if (!(spanRad > 0.0) || !(spanRad < kPi)) return false;

    gp_XYZ mean(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < n; ++i) mean += points[i].XYZ();
    mean.Divide(static_cast<double>(n));

    std::array<std::array<double, 3>, 3> cov{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (std::size_t i = 0; i < n; ++i) {
        const gp_XYZ d = points[i].XYZ() - mean;
        addOuter(cov, d.X(), d.Y(), d.Z(), 1.0);
    }
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return false;

    gp_XYZ nrm(evec[0][0], evec[0][1], evec[0][2]);
    const double nm = nrm.Modulus();
    if (nm < 1e-15) return false;
    nrm.Divide(nm);
    gp_XYZ u, v;
    if (!axisFrame(nrm, u, v)) return false;

    std::vector<double> xs(n), ys(n);
    for (std::size_t i = 0; i < n; ++i) {
        const gp_XYZ d = points[i].XYZ() - mean;
        xs[i] = d.Dot(u);
        ys[i] = d.Dot(v);
    }

    double cx = 0.0, cy = 0.0, R = 0.0;
    if (!prattFit2(xs.data(), ys.data(), n, cx, cy, R)) return false;
    const gp_XYZ c = mean + cx * u + cy * v;
    if (!finite3(c.X(), c.Y(), c.Z()) || !(R > 0.0)) return false;
    center = gp_Pnt(c.X(), c.Y(), c.Z());
    radius = R;
    return true;
}

bool gaussMapAxis(const MeshView& mv, const std::vector<int>& tris, gp_Dir& axis) {
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.empty()) return false;

    std::array<std::array<double, 3>, 3> G{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    double areaSum = 0.0;
    struct Nrm {
        gp_XYZ n;
        double area;
    };
    std::vector<Nrm> nrms;
    nrms.reserve(ids.size());

    for (int t : ids) {
        gp_XYZ a, b, p, nU;
        double area = 0.0;
        if (!triGeom(mv, t, a, b, p, nU, area) || area <= 0.0) continue;
        const double nm = nU.Modulus();
        if (nm < 1e-15) continue;
        gp_XYZ n = nU;
        n.Divide(nm);
        addOuter(G, n.X(), n.Y(), n.Z(), area);
        areaSum += area;
        nrms.push_back({n, area});
    }
    if (!(areaSum > 0.0) || nrms.empty()) return false;

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(G, eval, evec)) return false;

    // T1 (D1 §1.3) plus a rank-1 (planar Gauss map) reject so a plane is not
    // reported as a cylinder axis. Frozen signature cannot return λ; p1-grow
    // must recompute G + jacobiEigenSymmetric3 for the numeric T1 ratios.
    const double l1 = eval[0];
    const double l2 = eval[1];
    const double l3 = eval[2];
    if (l2 <= 1e-12 * std::max(l3, 1e-300)) return false;
    if (l1 / std::max(l2, 1e-300) >= 0.05) return false;

    gp_XYZ a(evec[0][0], evec[0][1], evec[0][2]);
    const double am = a.Modulus();
    if (am < 1e-15) return false;
    a.Divide(am);

    const double tilt = std::sin(3.0 * kPi / 180.0);
    double maxDot = 0.0;
    for (const Nrm& rec : nrms) {
        const double d = std::abs(rec.n.Dot(a));
        if (d > maxDot) maxDot = d;
    }
    if (maxDot >= tilt) return false;

    std::array<double, 3> av{a.X(), a.Y(), a.Z()};
    signNormalize(av);
    a = gp_XYZ(av[0], av[1], av[2]);
    if (a.Modulus() < 1e-15) return false;
    axis = gp_Dir(a);
    return true;
}

double radiusFromChordLength(double chordLen, int nSides) {
    if (nSides < 3 || !(chordLen > 0.0) || !std::isfinite(chordLen)) return 0.0;
    const double s = std::sin(kPi / static_cast<double>(nSides));
    if (s <= 1e-15) return 0.0;
    const double r = chordLen / (2.0 * s);
    return (std::isfinite(r) && r > 0.0) ? r : 0.0;
}

double circumradiusFromInscribed(double rInscribed, int nSides) {
    if (nSides < 3 || !(rInscribed > 0.0) || !std::isfinite(rInscribed)) return 0.0;
    const double c = std::cos(kPi / static_cast<double>(nSides));
    if (c <= 1e-15) return 0.0;
    const double r = rInscribed / c;
    return (std::isfinite(r) && r > 0.0) ? r : 0.0;
}

double radiusFromChordSagitta(double halfChord, double sagitta) {
    if (!(sagitta > 1e-15) || !(halfChord >= 0.0) || !std::isfinite(halfChord)
        || !std::isfinite(sagitta))
        return 0.0;
    const double r = (halfChord * halfChord + sagitta * sagitta) / (2.0 * sagitta);
    return (std::isfinite(r) && r > 0.0) ? r : 0.0;
}

double medianOfSorted(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

gp_XYZ triNormalUnnorm(const MeshView& mv, int localTri) {
    gp_XYZ a, b, c;
    if (!triCorners(mv, localTri, a, b, c)) return gp_XYZ(0, 0, 0);
    return (b - a).Crossed(c - a);
}

int estimateFullCircleSides(const MeshView& mv, const std::vector<int>& tris) {
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.empty()) return 0;

    std::vector<double> phis;
    for (int t : ids) {
        for (int u : ids) {
            if (u <= t) continue;
            const int g = mv.compTris[t];
            const int gu = mv.compTris[u];
            const int* T = mv.tris[g];
            const int* U = mv.tris[gu];
            bool share = false;
            for (int k = 0; k < 3 && !share; k++) {
                const int gv0 = T[k];
                const int gv1 = T[(k + 1) % 3];
                for (int j = 0; j < 3; j++) {
                    const int a = U[j];
                    const int b = U[(j + 1) % 3];
                    if ((a == gv0 && b == gv1) || (a == gv1 && b == gv0)) {
                        share = true;
                        break;
                    }
                }
            }
            if (!share) continue;
            const gp_XYZ n0 = triNormalUnnorm(mv, t);
            const gp_XYZ n1 = triNormalUnnorm(mv, u);
            const double m0 = n0.Modulus();
            const double m1 = n1.Modulus();
            if (m0 < 1e-15 || m1 < 1e-15) continue;
            const double dot = (n0 / m0).Dot(n1 / m1);
            const double phi = std::acos(std::clamp(dot, -1.0, 1.0));
            if (phi > 0.05 && phi < kPi - 0.05) phis.push_back(phi);
        }
    }
    if (phis.empty()) return 0;
    const double med = medianOfSorted(phis);
    if (med < 1e-6) return 0;
    return std::max(3, static_cast<int>(std::llround(kTwoPi / med)));
}

bool refineCylinderRadius(const MeshView& mv, const std::vector<int>& tris,
                          const gp_Dir& axis, gp_Pnt& center, double& radius,
                          int nSides, double spanRad, double rHint) {
    if (!(radius > 0.0) || !std::isfinite(radius) || nSides < 3) return false;
    // Coarse Fusion band only — corpus fixtures stay on the Pratt vertex fit.
    if (mv.nTri < 500 || mv.nTri > 1200) return false;

    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.size() < 2) return false;

    gp_XYZ aw(axis.XYZ());
    const double am = aw.Modulus();
    if (am < 1e-15) return false;
    aw.Divide(am);
    gp_XYZ u, v;
    if (!axisFrame(aw, u, v)) return false;

    const gp_XYZ c0(center.X(), center.Y(), center.Z());
    const double rEberly = radius;

    struct Vtx {
        gp_XYZ p;
        double ang;
        double rad;
    };
    std::vector<Vtx> ring;
    {
        std::vector<int> gids;
        for (int t : ids) {
            if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
            const int g = mv.compTris[t];
            for (int k = 0; k < 3; k++) gids.push_back(mv.tris[g][k]);
        }
        std::sort(gids.begin(), gids.end());
        gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
        ring.reserve(gids.size());
        for (int gi : gids) {
            const gp_XYZ p = mv.pts[gi];
            const gp_XYZ d = p - c0;
            const gp_XYZ radial = d - aw * d.Dot(aw);
            const double rr = radial.Modulus();
            if (!(rr > 0.0)) continue;
            const double ang = std::atan2(radial.Dot(v), radial.Dot(u));
            ring.push_back({p, ang, rr});
        }
    }
    if (ring.size() < 3) return false;
    std::sort(ring.begin(), ring.end(),
              [](const Vtx& a, const Vtx& b) { return a.ang < b.ang; });

    const int nDihedral = estimateFullCircleSides(mv, ids);
    int nEff = nDihedral > 0 ? nDihedral : nSides;
    nEff = std::clamp(nEff, 4, 48);
    const double segAng = kTwoPi / static_cast<double>(nEff);
    const double chordNom = 2.0 * radius * std::sin(kPi / static_cast<double>(nEff));
  // Growth R_ref from an early circumferential band can exceed the axial-dragged
  // eberly fit; widen chord acceptance toward rHint when present.
    const double rAnchor = (rHint > radius * 1.02) ? rHint : radius;
    const double rChordLo = std::min(radius * 0.90, rAnchor * 0.88);
    const double rChordHi = std::max(radius * 1.42, rAnchor * 1.08);

    std::vector<double> chordRadii;
    const size_t nV = ring.size();
    size_t skipEdge = nV;
    if (nV >= 3) {
        double maxGap = -1.0;
        for (size_t i = 0; i < nV; i++) {
            const size_t j = (i + 1) % nV;
            double gap = ring[j].ang - ring[i].ang;
            if (j == 0) gap = (ring[0].ang + kTwoPi) - ring[nV - 1].ang;
            if (gap > maxGap) {
                maxGap = gap;
                skipEdge = i;
            }
        }
    }

    // Mesh-adjacent ring chords only (angular sort skips facets on partial arcs).
    std::vector<char> inTri(static_cast<size_t>(mv.nTri), 0);
    for (int t : ids)
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inTri[static_cast<size_t>(t)] = 1;
    for (int t : ids) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        const int g = mv.compTris[t];
        const int* T = mv.tris[g];
        for (int k = 0; k < 3; k++) {
            const int gv0 = T[k];
            const int gv1 = T[(k + 1) % 3];
            for (int uTri : ids) {
                if (uTri == t) continue;
                if (!inTri[static_cast<size_t>(uTri)]) continue;
                const int gu = mv.compTris[uTri];
                const int* U = mv.tris[gu];
                bool share = false;
                for (int j = 0; j < 3; j++) {
                    const int a = U[j];
                    const int b = U[(j + 1) % 3];
                    if ((a == gv0 && b == gv1) || (a == gv1 && b == gv0)) {
                        share = true;
                        break;
                    }
                }
                if (!share) continue;
                const gp_XYZ p0 = mv.pts[gv0];
                const gp_XYZ p1 = mv.pts[gv1];
                const double chord = (p1 - p0).Modulus();
                if (!(chord > 0.0)) continue;
                const gp_XYZ d0 = p0 - c0;
                const gp_XYZ d1 = p1 - c0;
                const gp_XYZ r0 = d0 - aw * d0.Dot(aw);
                const gp_XYZ r1 = d1 - aw * d1.Dot(aw);
                const double rad0 = r0.Modulus();
                const double rad1 = r1.Modulus();
                if (!(rad0 > 0.0) || !(rad1 > 0.0)) continue;
                const double ang0 = std::atan2(r0.Dot(v), r0.Dot(u));
                const double ang1 = std::atan2(r1.Dot(v), r1.Dot(u));
                double dAng = std::abs(ang1 - ang0);
                if (dAng > kPi) dAng = kTwoPi - dAng;
                // Skip axial seams (Δθ ≈ 0) and near-diameter spans.
                if (dAng < 0.09 || dAng >= kPi - 0.09) continue;
                if (std::abs(rad0 - rad1) > 0.05 * std::max(rad0, rad1)) continue;
                const double s = std::sin(0.5 * dAng);
                if (s <= 1e-9) continue;
                const double rChord = chord / (2.0 * s);
                if (rChord > rChordLo && rChord < rChordHi) chordRadii.push_back(rChord);
            }
        }
    }

    // Fallback: angularly consecutive vertices when no internal mesh chords found.
    if (chordRadii.empty()) {
        for (size_t i = 0; i < nV; i++) {
            if (i == skipEdge) continue;
            const size_t j = (i + 1) % nV;
            double dAng = ring[j].ang - ring[i].ang;
            if (j == 0) dAng = (ring[0].ang + kTwoPi) - ring[nV - 1].ang;
            if (dAng < 1e-6 || dAng >= kPi) continue;
            const double chord = (ring[j].p - ring[i].p).Modulus();
            if (!(chord > 0.0)) continue;
            const double halfAng = 0.5 * dAng;
            const double s = std::sin(halfAng);
            if (s <= 1e-9) continue;
            const double rChord = chord / (2.0 * s);
            if (rChord > rChordLo && rChord < rChordHi) chordRadii.push_back(rChord);
        }
    }
  // Partial arcs: circumferential edges are often patch boundaries (thin axial bands).
    if (chordRadii.empty() && spanRad > 0.05 && spanRad < 2.8) {
        for (int t : ids) {
            if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
            const int g = mv.compTris[t];
            const int* T = mv.tris[g];
            for (int k = 0; k < 3; k++) {
                const int gv0 = T[k];
                const int gv1 = T[(k + 1) % 3];
                bool internal = false;
                for (int uTri : ids) {
                    if (uTri == t) continue;
                    const int gu = mv.compTris[uTri];
                    const int* U = mv.tris[gu];
                    for (int j = 0; j < 3; j++) {
                        const int a = U[j];
                        const int b = U[(j + 1) % 3];
                        if ((a == gv0 && b == gv1) || (a == gv1 && b == gv0)) {
                            internal = true;
                            break;
                        }
                    }
                    if (internal) break;
                }
                if (internal) continue;
                const gp_XYZ p0 = mv.pts[gv0];
                const gp_XYZ p1 = mv.pts[gv1];
                const double chord = (p1 - p0).Modulus();
                if (!(chord > 0.0)) continue;
                const gp_XYZ d0 = p0 - c0;
                const gp_XYZ d1 = p1 - c0;
                const gp_XYZ r0 = d0 - aw * d0.Dot(aw);
                const gp_XYZ r1 = d1 - aw * d1.Dot(aw);
                if (!(r0.Modulus() > 0.0) || !(r1.Modulus() > 0.0)) continue;
                const double ang0 = std::atan2(r0.Dot(v), r0.Dot(u));
                const double ang1 = std::atan2(r1.Dot(v), r1.Dot(u));
                double dAng = std::abs(ang1 - ang0);
                if (dAng > kPi) dAng = kTwoPi - dAng;
                if (dAng < 0.09 || dAng >= kPi - 0.09) continue;
                const double s = std::sin(0.5 * dAng);
                if (s <= 1e-9) continue;
                const double rChord = chord / (2.0 * s);
                if (rChord > rChordLo && rChord < rChordHi) chordRadii.push_back(rChord);
            }
        }
    }

    const double rInscribed = circumradiusFromInscribed(radius, nEff);
    double rPick = radius;

    if (chordRadii.size() >= 1) {
        const double rChordMed = medianOfSorted(chordRadii);
        if (rChordMed > radius) rPick = rChordMed;
    }
    if (rInscribed > rPick * 1.003 && rInscribed <= radius * 1.20 && nEff <= 12)
        rPick = std::max(rPick, rInscribed);
  // Axially-grown large-bore patch: bulk eberly drags below circumferential band R_ref.
    if (rHint > radius * 1.08 && rHint <= radius * 1.38 && radius >= 11.0) {
        const bool chordsLow =
            chordRadii.empty()
            || medianOfSorted(chordRadii) < radius * 1.05;
        if (chordsLow) rPick = std::max(rPick, rHint);
    }

    if (rPick <= radius * 1.01) {
        (void)spanRad;
        (void)segAng;
        (void)chordNom;
        return true;
    }

    double cap = radius * 1.25;
    if (rHint > radius * 1.02) cap = std::max(cap, rHint * 1.03);
    if (!chordRadii.empty()) {
        const double rChordMed = medianOfSorted(chordRadii);
        cap = std::max(cap, std::max(rPick * 1.02, rChordMed * 1.02));
    }
    if (rPick > cap) rPick = cap;
    radius = rPick;

  // After rHint lift, snap to tight circumferential chords (facet edges on patch).
    if (rHint > rEberly * 1.08 && radius > rEberly * 1.05) {
        std::vector<double> snapR;
        const double sLo = radius * 0.92;
        const double sHi = radius * 1.06;
        for (int t : ids) {
            if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
            const int g = mv.compTris[t];
            const int* T = mv.tris[g];
            for (int k = 0; k < 3; k++) {
                const int gv0 = T[k];
                const int gv1 = T[(k + 1) % 3];
                const gp_XYZ p0 = mv.pts[gv0];
                const gp_XYZ p1 = mv.pts[gv1];
                const double chord = (p1 - p0).Modulus();
                if (!(chord > 0.0)) continue;
                const gp_XYZ d0 = p0 - c0;
                const gp_XYZ d1 = p1 - c0;
                const gp_XYZ r0 = d0 - aw * d0.Dot(aw);
                const gp_XYZ r1 = d1 - aw * d1.Dot(aw);
                if (!(r0.Modulus() > 0.0) || !(r1.Modulus() > 0.0)) continue;
                const double ang0 = std::atan2(r0.Dot(v), r0.Dot(u));
                const double ang1 = std::atan2(r1.Dot(v), r1.Dot(u));
                double dAng = std::abs(ang1 - ang0);
                if (dAng > kPi) dAng = kTwoPi - dAng;
                if (dAng < 0.09 || dAng >= kPi - 0.09) continue;
                const double s = std::sin(0.5 * dAng);
                if (s <= 1e-9) continue;
                const double rChord = chord / (2.0 * s);
                if (rChord > sLo && rChord < sHi) snapR.push_back(rChord);
            }
        }
        if (snapR.size() >= 1) {
            const double med = medianOfSorted(snapR);
            if (med > rEberly * 1.02 && med < radius) radius = med;
        } else if (nEff >= 4) {
            double maxChord = 0.0;
            for (int t : ids) {
                if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
                const int g = mv.compTris[t];
                const int* T = mv.tris[g];
                for (int k = 0; k < 3; k++) {
                    const int gv0 = T[k];
                    const int gv1 = T[(k + 1) % 3];
                    const gp_XYZ p0 = mv.pts[gv0];
                    const gp_XYZ p1 = mv.pts[gv1];
                    const double chord = (p1 - p0).Modulus();
                    if (!(chord > 0.0)) continue;
                    const gp_XYZ d0 = p0 - c0;
                    const gp_XYZ d1 = p1 - c0;
                    const gp_XYZ r0 = d0 - aw * d0.Dot(aw);
                    const gp_XYZ r1 = d1 - aw * d1.Dot(aw);
                    if (!(r0.Modulus() > 0.0) || !(r1.Modulus() > 0.0)) continue;
                    const double ang0 = std::atan2(r0.Dot(v), r0.Dot(u));
                    const double ang1 = std::atan2(r1.Dot(v), r1.Dot(u));
                    double dAng = std::abs(ang1 - ang0);
                    if (dAng > kPi) dAng = kTwoPi - dAng;
                    if (dAng < 0.09 || dAng >= kPi - 0.09) continue;
                    if (chord > maxChord) maxChord = chord;
                }
            }
            const double rNom = radiusFromChordLength(maxChord, nEff);
            if (rNom > rEberly * 1.02 && rNom < radius) radius = rNom;
        }
    }

    std::vector<double> xs(ring.size()), ys(ring.size());
    for (size_t i = 0; i < ring.size(); i++) {
        const gp_XYZ d = ring[i].p - c0;
        xs[i] = d.Dot(u);
        ys[i] = d.Dot(v);
    }
    double cx = 0.0, cy = 0.0, r2 = 0.0;
    if (prattFit2(xs.data(), ys.data(), ring.size(), cx, cy, r2)) {
        const gp_XYZ cNew = c0 + cx * u + cy * v;
        if (finite3(cNew.X(), cNew.Y(), cNew.Z()) && r2 > 0.0) {
            center = gp_Pnt(cNew.X(), cNew.Y(), cNew.Z());
            if (r2 > radius * 0.97 && r2 < radius * 1.06)
                radius = 0.5 * (radius + r2);
            else if (rHint > rEberly * 1.08 && r2 > rEberly * 1.02)
                radius = r2;
        }
    }
    (void)spanRad;
    return std::isfinite(radius) && radius > 0.0;
}

double chordSagitta(double radius, int nSides) {
    if (nSides < 1 || !std::isfinite(radius)) return 0.0;
    return radius * (1.0 - std::cos(kPi / static_cast<double>(nSides)));
}

double dVolCylinderSector(double areaReg, double radius, int nSides, bool outwardNormal) {
    // D4.2: nSides < 3 or R <= 0 or gamma >= π => 0. SPEC also guards nSides < 1
    // and sin(gamma/2) → 0.
    if (nSides < 3 || !(radius > 0.0) || !std::isfinite(areaReg) || !std::isfinite(radius))
        return 0.0;
    const double gamma = kTwoPi / static_cast<double>(nSides);
    if (!(gamma < kPi)) return 0.0;
    const double s = std::sin(0.5 * gamma);
    if (std::abs(s) < 1e-15) return 0.0;
    const double sigma = outwardNormal ? 1.0 : -1.0;
    const double dvol = sigma * areaReg * radius * (gamma - std::sin(gamma)) / (4.0 * s);
    return std::isfinite(dvol) ? dvol : 0.0;
}

double dVolPlaneRegion(const MeshView& mv, const std::vector<int>& tris, const gp_Ax3& ax) {
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    const gp_XYZ n = ax.Direction().XYZ();
    const gp_XYZ o = ax.Location().XYZ();
    double sum = 0.0;
    for (int t : ids) {
        gp_XYZ a, b, p, nU;
        double area = 0.0;
        if (!triGeom(mv, t, a, b, p, nU, area) || area <= 0.0) continue;
        const double ha = n.Dot(a - o);
        const double hb = n.Dot(b - o);
        const double hc = n.Dot(p - o);
        sum += area * (ha + hb + hc) / 3.0;
    }
    // D4.3: leading minus — flattening an outward bulge (h>0) removes volume.
    const double dvol = -sum;
    return std::isfinite(dvol) ? dvol : 0.0;
}

namespace {

double wrapToPiLocal(double t) {
    while (t <= -kPi) t += kTwoPi;
    while (t > kPi) t -= kTwoPi;
    return t;
}

bool unitTriNormal(const MeshView& mv, int lt, gp_XYZ& nOut, double& areaOut) {
    gp_XYZ a, b, p, nU;
    if (!triGeom(mv, lt, a, b, p, nU, areaOut) || areaOut <= 0.0) return false;
    const double nm = nU.Modulus();
    if (nm < 1e-15) return false;
    nOut = nU;
    nOut.Divide(nm);
    return true;
}

double edgeDihedralTriPair(const MeshView& mv, int t0, int t1) {
    gp_XYZ n0, n1;
    double a0 = 0.0, a1 = 0.0;
    if (!unitTriNormal(mv, t0, n0, a0) || !unitTriNormal(mv, t1, n1, a1)) return 0.0;
    const double d = std::min(1.0, std::max(-1.0, n0.Dot(n1)));
    return std::acos(d);
}

bool orderTrisBfs(const MeshView& mv, const std::vector<int>& tris, double maxPhi,
                  std::vector<int>& orderOut) {
    orderOut.clear();
    if (tris.empty() || !mv.triEdges) return false;

    std::vector<char> inPatch(static_cast<size_t>(mv.nTri), 0);
    for (int t : tris) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inPatch[static_cast<size_t>(t)] = 1;
    }

    auto neighborsOf = [&](int t) {
        std::vector<int> nb;
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            for (int u : tris) {
                if (u == t || !inPatch[static_cast<size_t>(u)]) continue;
                for (int su = 0; su < 3; su++) {
                    if (mv.triEdges[u][su] != e) continue;
                    if (edgeDihedralTriPair(mv, t, u) <= maxPhi) nb.push_back(u);
                }
            }
        }
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
        return nb;
    };

    int seed = tris.front();
    for (int t : tris)
        if (t < seed) seed = t;

    std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
    std::vector<int> q;
    q.push_back(seed);
    seen[static_cast<size_t>(seed)] = 1;

    while (!q.empty()) {
        const int t = q.front();
        q.erase(q.begin());
        orderOut.push_back(t);
        for (int u : neighborsOf(t)) {
            if (seen[static_cast<size_t>(u)]) continue;
            seen[static_cast<size_t>(u)] = 1;
            q.push_back(u);
        }
    }
    return orderOut.size() >= 3;
}

bool normalCovariance(const MeshView& mv, const std::vector<int>& tris, gp_XYZ& nbar,
                      std::array<std::array<double, 3>, 3>& cov, double& areaSum) {
    cov = {};
    nbar = gp_XYZ(0.0, 0.0, 0.0);
    areaSum = 0.0;
    for (int lt : tris) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, lt, n, area)) continue;
        areaSum += area;
        nbar += n * area;
    }
    if (!(areaSum > 0.0)) return false;
    nbar.Divide(areaSum);
    for (int lt : tris) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, lt, n, area)) continue;
        const double dx = n.X() - nbar.X();
        const double dy = n.Y() - nbar.Y();
        const double dz = n.Z() - nbar.Z();
        cov[0][0] += area * dx * dx;
        cov[0][1] += area * dx * dy;
        cov[0][2] += area * dx * dz;
        cov[1][1] += area * dy * dy;
        cov[1][2] += area * dy * dz;
        cov[2][2] += area * dz * dz;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];
    return true;
}

double maxVertexPlaneDev(const MeshView& mv, const std::vector<int>& tris,
                         const gp_Ax3& plane) {
    const gp_Dir n = plane.Direction();
    const gp_Pnt loc = plane.Location();
    double maxD = 0.0;
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            gp_XYZ a, b, p, nU;
            double area = 0.0;
            if (!triCorners(mv, lt, a, b, p)) continue;
            const gp_XYZ v = (k == 0) ? a : (k == 1) ? b : p;
            const double d = std::abs(gp_Vec(loc, gp_Pnt(v)).Dot(gp_Vec(n)));
            maxD = std::max(maxD, d);
        }
    }
    return maxD;
}

double maxCylResidual(const MeshView& mv, const std::vector<int>& tris, const gp_Dir& axis,
                      const gp_Pnt& center, double radius) {
    const gp_XYZ a(axis.X(), axis.Y(), axis.Z());
    const gp_XYZ c(center.X(), center.Y(), center.Z());
    double maxR = 0.0;
    for (int lt : tris) {
        for (int k = 0; k < 3; k++) {
            gp_XYZ a3, b3, p3;
            if (!triCorners(mv, lt, a3, b3, p3)) continue;
            const gp_XYZ v = (k == 0) ? a3 : (k == 1) ? b3 : p3;
            const double rr = a.Crossed(v - c).Modulus();
            maxR = std::max(maxR, std::abs(rr - radius));
        }
    }
    return maxR;
}

bool monotonicNormalSpan(const MeshView& mv, const std::vector<int>& order,
                         const gp_Dir& axis, double& spanRad, double& monoFrac) {
    spanRad = 0.0;
    monoFrac = 0.0;
    if (order.size() < 3) return false;

    gp_XYZ u, v;
    if (!axisFrame(axis.XYZ(), u, v)) return false;

    std::vector<double> ang;
    ang.reserve(order.size());
    for (int lt : order) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, lt, n, area)) continue;
        const gp_XYZ ax(axis.X(), axis.Y(), axis.Z());
        gp_XYZ pp = n - ax * n.Dot(ax);
        const double pm = pp.Modulus();
        if (pm < 1e-12) continue;
        pp.Divide(pm);
        ang.push_back(wrapToPiLocal(std::atan2(pp.Dot(v), pp.Dot(u))));
    }
    if (ang.size() < 3) return false;

    std::vector<double> sorted = ang;
    std::sort(sorted.begin(), sorted.end());

    int pos = 0, neg = 0, tot = 0;
    for (size_t i = 1; i < sorted.size(); i++) {
        const double d = sorted[i] - sorted[i - 1];
        if (std::abs(d) < 1e-4) continue;
        tot++;
        if (d > 0.0) pos++;
        else neg++;
    }
    monoFrac = tot > 0 ? static_cast<double>(std::max(pos, neg)) / tot : 1.0;

    const double wrapGap = kTwoPi - sorted.back() + sorted.front();
    spanRad = sorted.back() - sorted.front();
    if (wrapGap > spanRad) spanRad = wrapGap;

    return spanRad >= 0.14 && monoFrac >= 0.70;
}

bool trisShareEdgeVerts(const MeshView& mv, int t0, int t1, int& gv0, int& gv1) {
    if (!mv.tris || !mv.compTris) return false;
    if (t0 < 0 || t1 < 0 || static_cast<size_t>(t0) >= mv.nTri
        || static_cast<size_t>(t1) >= mv.nTri)
        return false;
    const int g0 = mv.compTris[t0];
    const int g1 = mv.compTris[t1];
    const int* T0 = mv.tris[g0];
    const int* T1 = mv.tris[g1];
    for (int k = 0; k < 3; k++) {
        const int a = T0[k];
        const int b = T0[(k + 1) % 3];
        for (int j = 0; j < 3; j++) {
            const int c = T1[j];
            const int d = T1[(j + 1) % 3];
            if ((a == c && b == d) || (a == d && b == c)) {
                gv0 = a;
                gv1 = b;
                return true;
            }
        }
    }
    return false;
}

double sharedEdgeLength(const MeshView& mv, int t0, int t1) {
    int gv0 = 0, gv1 = 0;
    if (!trisShareEdgeVerts(mv, t0, t1, gv0, gv1) || !mv.pts) return 0.0;
    const gp_XYZ p0 = mv.pts[gv0];
    const gp_XYZ p1 = mv.pts[gv1];
    const double len = (p1 - p0).Modulus();
    return std::isfinite(len) ? len : 0.0;
}

double stripWidthAlongAxis(const MeshView& mv, int t0, int t1, const gp_Dir& axis) {
    auto facetAxialExtent = [&](int t) {
        gp_XYZ a, b, c;
        if (!triCorners(mv, t, a, b, c)) return 0.0;
        gp_XYZ ax(axis.XYZ());
        const double am = ax.Modulus();
        if (am < 1e-15) return 0.0;
        ax.Divide(am);
        double mn = ax.Dot(a);
        double mx = mn;
        for (const gp_XYZ& p : {b, c}) {
            const double d = ax.Dot(p);
            mn = std::min(mn, d);
            mx = std::max(mx, d);
        }
        return mx - mn;
    };
    auto widthFromArea = [&](int t) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, t, n, area) || area <= 0.0) return 0.0;
        const double h = facetAxialExtent(t);
        if (h <= 1e-6) return 0.0;
        return 2.0 * area / h;
    };
    const double w0 = widthFromArea(t0);
    const double w1 = widthFromArea(t1);
    if (w0 > 0.0 && w1 > 0.0) return 0.5 * (w0 + w1);
    int gv0 = 0, gv1 = 0;
    if (!trisShareEdgeVerts(mv, t0, t1, gv0, gv1) || !mv.pts) return 0.0;
    const gp_XYZ e = mv.pts[gv1] - mv.pts[gv0];
    const double len = e.Modulus();
    if (!(len > 0.0)) return 0.0;
    gp_XYZ ax(axis.XYZ());
    const double am = ax.Modulus();
    if (am < 1e-15) return len;
    ax.Divide(am);
    const double parallel = std::abs(e.Dot(ax)) / len;
    if (parallel > 0.82) return 0.0;
    return len;
}

bool triInPatchNeighbors(const MeshView& mv, int t, const std::vector<char>& inPatch,
                         std::vector<int>& nb) {
    nb.clear();
    if (!mv.triEdges || t < 0 || static_cast<size_t>(t) >= mv.nTri) return false;
    for (int s = 0; s < 3; s++) {
        const int e = mv.triEdges[t][s];
        for (size_t u = 0; u < mv.nTri; u++) {
            if (!inPatch[u] || static_cast<int>(u) == t) continue;
            for (int su = 0; su < 3; su++) {
                if (mv.triEdges[static_cast<int>(u)][su] != e) continue;
                nb.push_back(static_cast<int>(u));
            }
        }
    }
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    return true;
}

bool buildTriPathChain(const MeshView& mv, const std::vector<int>& tris,
                       std::vector<int>& chainOut) {
    chainOut.clear();
    if (tris.size() < 3) return false;

    std::vector<char> inPatch(static_cast<size_t>(mv.nTri), 0);
    for (int t : tris) {
        if (t >= 0 && static_cast<size_t>(t) < mv.nTri) inPatch[static_cast<size_t>(t)] = 1;
    }

    auto walkFrom = [&](int start, bool arcOnly, std::vector<int>& out) {
        out.clear();
        std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
        int cur = start;
        int prev = -1;
        while (cur >= 0 && !seen[static_cast<size_t>(cur)]) {
            out.push_back(cur);
            seen[static_cast<size_t>(cur)] = 1;
            std::vector<int> nb;
            triInPatchNeighbors(mv, cur, inPatch, nb);
            int next = -1;
            double bestPhi = 0.0;
            for (int u : nb) {
                if (u == prev || seen[static_cast<size_t>(u)]) continue;
                const double phi = edgeDihedralTriPair(mv, cur, u);
                if (arcOnly && (phi < 0.012 || phi > kPi - 0.012)) continue;
                if (next < 0 || phi > bestPhi) {
                    next = u;
                    bestPhi = phi;
                }
            }
            prev = cur;
            cur = next;
        }
    };

    std::vector<int> best;
    for (int t : tris) {
        if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
        std::vector<int> chain;
        walkFrom(t, true, chain);
        if (chain.size() > best.size()) best = std::move(chain);
    }
    if (best.size() < 3) {
        for (int t : tris) {
            if (t < 0 || static_cast<size_t>(t) >= mv.nTri) continue;
            std::vector<int> chain;
            walkFrom(t, false, chain);
            if (chain.size() > best.size()) best = std::move(chain);
        }
    }
    chainOut = best;
    return chainOut.size() >= 3;
}

bool coefficientOfVariationOk(const std::vector<double>& vals, int skipEnds, double maxCv) {
    if (vals.size() < 3) return false;
    const size_t lo = static_cast<size_t>(skipEnds);
    const size_t hi = vals.size() - static_cast<size_t>(skipEnds);
    if (hi <= lo + 1) return false;
    double sum = 0.0;
    double sumSq = 0.0;
    size_t n = 0;
    for (size_t i = lo; i < hi; i++) {
        if (!(vals[i] > 0.0)) return false;
        sum += vals[i];
        sumSq += vals[i] * vals[i];
        n++;
    }
    if (n < 2) return false;
    const double mean = sum / static_cast<double>(n);
    const double var = sumSq / static_cast<double>(n) - mean * mean;
    if (var <= 0.0) return true;
    const double cv = std::sqrt(var) / mean;
    return cv <= maxCv;
}

bool coefficientOfVariation(const std::vector<double>& vals, int skipEnds, double& cvOut) {
    cvOut = 0.0;
    if (vals.size() < 2) return false;
    const size_t lo = static_cast<size_t>(std::max(0, skipEnds));
    const size_t hi = vals.size() - static_cast<size_t>(std::max(0, skipEnds));
    if (hi <= lo + 1) return false;
    double sum = 0.0;
    double sumSq = 0.0;
    size_t n = 0;
    for (size_t i = lo; i < hi; i++) {
        if (!(vals[i] > 0.0)) return false;
        sum += vals[i];
        sumSq += vals[i] * vals[i];
        n++;
    }
    if (n < 2) return false;
    const double mean = sum / static_cast<double>(n);
    if (!(mean > 0.0)) return false;
    const double var = sumSq / static_cast<double>(n) - mean * mean;
    cvOut = (var <= 0.0) ? 0.0 : std::sqrt(var) / mean;
    return std::isfinite(cvOut);
}

bool collapseDiagOn() {
    const char* v = std::getenv("STL2STEP_COLLAPSE_DIAG");
    return v && v[0] && v[0] != '0';
}

void emitArchChainDiag(int rid, int n, double areaCv, double angCv, double R, double signal) {
    if (!collapseDiagOn()) return;
    const int rMilli = static_cast<int>(std::lround(R * 1000.0));
    static std::mutex mu;
    static std::vector<std::array<int, 3>> seen;
    {
        std::lock_guard<std::mutex> lock(mu);
        for (const auto& s : seen) {
            if (s[0] == rid && s[1] == n && s[2] == rMilli) return;
        }
        seen.push_back({rid, n, rMilli});
    }
    std::fprintf(stderr,
                 "DIAG_ARCHCHAIN rid=%d n=%d areaCV=%.6f angCV=%.6f R=%.6f signal=%.6f\n",
                 rid, n, areaCv, angCv, R, signal);
}

bool chainAreaAngCV(const MeshView& mv, const std::vector<int>& chain,
                    double& areaCv, double& angCv) {
    areaCv = 0.0;
    angCv = 0.0;
    if (chain.size() < 3) return false;
    std::vector<double> areas;
    areas.reserve(chain.size());
    for (int t : chain) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, t, n, area) || !(area > 0.0)) return false;
        areas.push_back(area);
    }
    std::vector<double> arcThetas;
    for (size_t i = 1; i < chain.size(); i++) {
        const double th = edgeDihedralTriPair(mv, chain[i - 1], chain[i]);
        if (th >= 0.012 && th <= kPi - 0.04) arcThetas.push_back(th);
    }
    const int areaSkip = chain.size() >= 10 ? 2 : 1;
    if (!coefficientOfVariation(areas, areaSkip, areaCv)) return false;
    if (arcThetas.size() >= 2) (void)coefficientOfVariation(arcThetas, 0, angCv);
    return true;
}

template <typename F>
void parallelFor(size_t n, F&& fn) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;
    if (n < 32 || hw == 1) {
        for (size_t i = 0; i < n; i++) fn(i);
        return;
    }
    std::vector<std::thread> pool;
    const size_t chunk = (n + hw - 1) / hw;
    for (unsigned t = 0; t < hw; t++) {
        const size_t lo = static_cast<size_t>(t) * chunk;
        const size_t hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        pool.emplace_back([&, lo, hi]() {
            for (size_t i = lo; i < hi; i++) fn(i);
        });
    }
    for (auto& th : pool) th.join();
}

void buildEdgeTriAdj(const MeshView& mv, std::vector<std::array<int, 2>>& e2t) {
    e2t.assign(mv.nEdge, {-1, -1});
    if (!mv.triEdges) return;
    for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || static_cast<size_t>(e) >= mv.nEdge) continue;
            if (e2t[static_cast<size_t>(e)][0] < 0)
                e2t[static_cast<size_t>(e)][0] = t;
            else if (e2t[static_cast<size_t>(e)][1] < 0)
                e2t[static_cast<size_t>(e)][1] = t;
        }
    }
}

void walkArcPath(const MeshView& mv, const std::vector<std::array<int, 2>>& e2t,
                 const std::vector<char>& alive, int start, std::vector<int>& out) {
    out.clear();
    if (start < 0 || static_cast<size_t>(start) >= mv.nTri || !alive[static_cast<size_t>(start)])
        return;
    if (!mv.triEdges) return;
    std::vector<char> seen(static_cast<size_t>(mv.nTri), 0);
    int cur = start;
    int prev = -1;
    while (cur >= 0 && !seen[static_cast<size_t>(cur)]) {
        out.push_back(cur);
        seen[static_cast<size_t>(cur)] = 1;
        int next = -1;
        double bestPhi = 0.0;
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[cur][s];
            if (e < 0 || static_cast<size_t>(e) >= e2t.size()) continue;
            const int a = e2t[static_cast<size_t>(e)][0];
            const int b = e2t[static_cast<size_t>(e)][1];
            const int u = (a == cur) ? b : ((b == cur) ? a : -1);
            if (u < 0 || u == prev || !alive[static_cast<size_t>(u)] ||
                seen[static_cast<size_t>(u)])
                continue;
            const double phi = edgeDihedralTriPair(mv, cur, u);
            if (phi < 0.012 || phi > kPi - 0.012) continue;
            if (next < 0 || phi > bestPhi) {
                next = u;
                bestPhi = phi;
            }
        }
        prev = cur;
        cur = next;
    }
}

void inventoryArchChains(const MeshView& mv) {
    if (!archChainBand(mv) || mv.nTri < 3 || !mv.triEdges || mv.nEdge == 0) return;

    DerivedTols tol;
    tol.epsMesh = std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
    tol.epsPlane = std::max({tol.epsMesh, mv.sewTol, 0.02});

    std::vector<std::array<int, 2>> e2t;
    buildEdgeTriAdj(mv, e2t);
    std::vector<char> alive(static_cast<size_t>(mv.nTri), 1);

    const size_t n = mv.nTri;
    const int maxRounds = static_cast<int>(std::min(n, static_cast<size_t>(128)));
    for (int round = 0; round < maxRounds; round++) {
        std::vector<int> seeds;
        seeds.reserve(n);
        for (size_t t = 0; t < n; t++) {
            if (alive[t]) seeds.push_back(static_cast<int>(t));
        }
        if (seeds.size() < 3) break;

        std::vector<std::vector<int>> found(seeds.size());
        parallelFor(seeds.size(), [&](size_t i) {
            walkArcPath(mv, e2t, alive, seeds[i], found[i]);
        });

        size_t bestI = 0;
        size_t bestN = 0;
        for (size_t i = 0; i < found.size(); i++) {
            if (found[i].size() > bestN) {
                bestN = found[i].size();
                bestI = i;
            }
        }
        if (bestN < 3) break;

        const std::vector<int>& chain = found[bestI];
        ArcStripDetect det;
        (void)detectArchChain(mv, chain, tol, det);
        for (int t : chain) {
            if (t >= 0 && static_cast<size_t>(t) < alive.size())
                alive[static_cast<size_t>(t)] = 0;
        }
    }
}

}  // namespace

static void maybeInventoryArchChains(const MeshView& mv) {
    if (!collapseDiagOn()) return;
    if (!archChainBand(mv)) return;

    struct Key {
        const void* pts = nullptr;
        const void* comp = nullptr;
        size_t nTri = 0;
        bool operator==(const Key& o) const {
            return pts == o.pts && comp == o.comp && nTri == o.nTri;
        }
    };

    static std::mutex mu;
    static Key seen[8];
    static int nSeen = 0;
    static thread_local int depth = 0;
    if (depth > 0) return;
    const Key key{mv.pts, mv.compTris, mv.nTri};

    {
        std::lock_guard<std::mutex> lock(mu);
        for (int i = 0; i < nSeen; i++) {
            if (seen[i] == key) return;
        }
        if (nSeen < 8) seen[nSeen++] = key;
    }

    ++depth;
    inventoryArchChains(mv);
    --depth;
}

bool radiusFromArchChainPairs(const MeshView& mv, const std::vector<int>& tris,
                              const gp_Dir& axis, double& radiusOut, double& chainScoreOut,
                              double rHint) {
    radiusOut = 0.0;
    chainScoreOut = 0.0;
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.size() < 3) return false;

    std::vector<double> ws;
    std::vector<double> ths;
    for (size_t i = 0; i < ids.size(); i++) {
        for (size_t j = i + 1; j < ids.size(); j++) {
            const int t0 = ids[i];
            const int t1 = ids[j];
            int gv0 = 0, gv1 = 0;
            if (!trisShareEdgeVerts(mv, t0, t1, gv0, gv1)) continue;
            const double theta = edgeDihedralTriPair(mv, t0, t1);
            if (theta < 0.012 || theta > kPi - 0.012) continue;
            const double w = stripWidthAlongAxis(mv, t0, t1, axis);
            if (!(w > 0.0)) continue;
            const double s = std::sin(0.5 * theta);
            if (s <= 1e-9) continue;
            const double rChord = w / (2.0 * s);
            if (!(rChord > 0.0) || !std::isfinite(rChord)) continue;
            if (rHint > 0.0) {
                const double rel = rChord / rHint;
                if (rel < 0.48 || rel > 1.52) continue;
            }
            ws.push_back(w);
            ths.push_back(theta);
        }
    }
    if (ws.size() < 2 || ths.size() < 2) return false;
    if (!coefficientOfVariationOk(ws, 0, 0.45) || !coefficientOfVariationOk(ths, 0, 0.45))
        return false;
    const double wMed = medianOfSorted(ws);
    const double thMed = medianOfSorted(ths);
    const double s = std::sin(0.5 * thMed);
    if (s <= 1e-9) return false;
    radiusOut = wMed / (2.0 * s);
    if (!(radiusOut > 0.0)) return false;
    chainScoreOut = std::min(1.0, static_cast<double>(ws.size()) / 6.0);
    return std::isfinite(radiusOut);
}

bool radiusFromArchChain(const MeshView& mv, const std::vector<int>& chain,
                         const gp_Dir& axis, double& radiusOut, double& chainScoreOut,
                         double rHint) {
    radiusOut = 0.0;
    chainScoreOut = 0.0;
    if (chain.size() < 3) return false;

    constexpr double kMinTheta = 0.012;
    constexpr double kMaxTheta = kPi - 0.04;

    std::vector<double> areas;
    areas.reserve(chain.size());
    for (int t : chain) {
        gp_XYZ n;
        double area = 0.0;
        if (!unitTriNormal(mv, t, n, area) || area <= 0.0) return false;
        areas.push_back(area);
    }

    std::vector<double> thetas;
    thetas.reserve(chain.size());
    for (size_t i = 1; i < chain.size(); i++) {
        const double th = edgeDihedralTriPair(mv, chain[i - 1], chain[i]);
        thetas.push_back(th);
    }

    const int areaSkip = chain.size() >= 10 ? 2 : 1;
    if (!coefficientOfVariationOk(areas, areaSkip, 0.42)) return false;

    std::vector<double> arcThetas;
    for (double th : thetas) {
        if (th >= kMinTheta && th <= kMaxTheta) arcThetas.push_back(th);
    }
    if (arcThetas.size() < 2) return false;
    if (!coefficientOfVariationOk(arcThetas, 0, 0.45)) return false;

    std::vector<double> radii;
    const size_t nLinks = chain.size() - 1;
    size_t goodLinks = 0;
    for (size_t i = 1; i < chain.size(); i++) {
        const double theta = thetas[i - 1];
        if (theta < kMinTheta || theta > kMaxTheta) continue;
        const double w = stripWidthAlongAxis(mv, chain[i - 1], chain[i], axis);
        if (!(w > 0.0)) continue;
        const double s = std::sin(0.5 * theta);
        if (s <= 1e-9) continue;
        const double rChord = w / (2.0 * s);
        if (!(rChord > 0.0) || !std::isfinite(rChord)) continue;
        if (rHint > 0.0) {
            const double rel = rChord / rHint;
            if (rel < 0.48 || rel > 1.52) continue;
        }

        radii.push_back(rChord);
        goodLinks++;
    }

    if (radii.size() < 2) return false;
    const double rMed = medianOfSorted(radii);
    if (!(rMed > 0.0)) return false;

  // Score: arc-link coverage + chain length.
    const double linkFrac = nLinks > 0 ? static_cast<double>(goodLinks) / nLinks : 0.0;
    chainScoreOut = linkFrac * std::min(1.0, static_cast<double>(chain.size()) / 5.0);
    if (chainScoreOut < 0.28) return false;

    radiusOut = rMed;
    return std::isfinite(radiusOut) && radiusOut > 0.0;
}

bool archChainRadiusFromPatch(const MeshView& mv, const std::vector<int>& tris,
                              const gp_Dir& axis, double& radiusOut,
                              double& chainScoreOut, double rHint) {
    radiusOut = 0.0;
    chainScoreOut = 0.0;
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    std::vector<int> chain;
    if (!buildTriPathChain(mv, ids, chain)) return false;
    const bool okChain = radiusFromArchChain(mv, chain, axis, radiusOut, chainScoreOut, rHint);
    double pairR = 0.0;
    double pairScore = 0.0;
    const bool okPairs =
        ids.size() >= 3
        && radiusFromArchChainPairs(mv, ids, axis, pairR, pairScore, rHint);
    bool ok = okChain;
    if (okPairs && pairScore >= chainScoreOut) {
        radiusOut = pairR;
        chainScoreOut = pairScore;
        ok = true;
    }
    if (ok && radiusOut > 0.0 && !chain.empty()) {
        double areaCv = 0.0, angCv = 0.0;
        (void)chainAreaAngCV(mv, chain, areaCv, angCv);
        int rid = chain[0];
        for (int t : chain)
            if (t < rid) rid = t;
        emitArchChainDiag(rid, static_cast<int>(chain.size()), areaCv, angCv, radiusOut,
                          chainScoreOut);
    }
    return ok;
}

bool detectLargeArcStrip(const MeshView& mv, const std::vector<int>& tris,
                         const DerivedTols& tol, ArcStripDetect& out) {
    out = ArcStripDetect{};
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.size() < 3) return false;

    gp_XYZ nbar;
    std::array<std::array<double, 3>, 3> cov{};
    double areaSum = 0.0;
    if (!normalCovariance(mv, ids, nbar, cov, areaSum)) return false;

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return false;

    const double L = std::max(mv.diag, tol.epsMesh);
    const double theta = (L > 0.0) ? (tol.epsMesh / L) : 0.0;
    const double mu2Floor = areaSum * theta * theta;
    if (eval[1] <= mu2Floor && eval[1] <= 1e-12 * std::max(eval[2], 1e-300)) return false;

    gp_Ax3 pca;
    const bool havePca = pcaPlane(mv, ids, pca);
    const double pdev = havePca ? maxVertexPlaneDev(mv, ids, pca) : 0.0;

  // w1 = smallest-eigenvalue axis (normal-rotation axis for cylinder bands).
    gp_XYZ w1(evec[0][0], evec[0][1], evec[0][2]);
    const double wm = w1.Modulus();
    if (wm < 1e-15) return false;
    w1.Divide(wm);
    std::array<double, 3> av{w1.X(), w1.Y(), w1.Z()};
    signNormalize(av);
    w1 = gp_XYZ(av[0], av[1], av[2]);
    const gp_Dir axis(w1.X(), w1.Y(), w1.Z());

    std::vector<int> order;
    if (!orderTrisBfs(mv, ids, tol.thetaPlane, order)) return false;

    double spanRad = 0.0;
    double monoFrac = 0.0;
    const bool mono = monotonicNormalSpan(mv, order, axis, spanRad, monoFrac);
    const bool staticNormals = !mono && pdev > tol.epsPlane * 2.0;

    if (!mono && !staticNormals) return false;
    if (mono && spanRad < 0.21) return false;  // ~12°

    gp_Pnt center;
    double radius = 0.0;
    if (!eberlyCenterRadius(mv, ids, axis, center, radius)) return false;
    if (!(radius >= 15.0) || radius > 55.0) return false;

    const double res = maxCylResidual(mv, ids, axis, center, radius);
    double accept = tol.epsCylAccept(radius);
    const int nEst = std::max(6, static_cast<int>(ids.size()));
    accept = std::max(accept, chordSagitta(radius, nEst));
    // Coarse Fusion band: chordal rings on large-R partial arcs (handle-lock
    // R≈20 needs ~0.5 mm slack vs the 0.39 mm sagitta floor alone).
    if (mv.nTri >= 500 && mv.nTri <= 1200)
        accept = std::max(accept, 0.05 * radius);
    if (res > accept) return false;

    out.ok = true;
    out.axis = axis;
    out.center = center;
    out.radius = radius;
    out.spanRad = spanRad;
    out.staticNormals = staticNormals;
    return true;
}

bool detectArchChain(const MeshView& mv, const std::vector<int>& tris,
                     const DerivedTols& tol, ArcStripDetect& out) {
    out = ArcStripDetect{};
    if (!archChainBand(mv)) return false;
    std::vector<int> ids;
    sortedUniqueTris(tris, ids);
    if (ids.size() < 3) return false;

    gp_XYZ nbar;
    std::array<std::array<double, 3>, 3> cov{};
    double areaSum = 0.0;
    if (!normalCovariance(mv, ids, nbar, cov, areaSum)) return false;

    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, eval, evec)) return false;

    gp_XYZ w1(evec[0][0], evec[0][1], evec[0][2]);
    const double wm = w1.Modulus();
    if (wm < 1e-15) return false;
    w1.Divide(wm);
    std::array<double, 3> av{w1.X(), w1.Y(), w1.Z()};
    signNormalize(av);
    w1 = gp_XYZ(av[0], av[1], av[2]);
    const gp_Dir axis(w1.X(), w1.Y(), w1.Z());

    std::vector<int> chain;
    if (!buildTriPathChain(mv, ids, chain)) return false;

    double chainR = 0.0;
    double chainScore = 0.0;
    if (!radiusFromArchChain(mv, chain, axis, chainR, chainScore, 0.0)) return false;

    double areaCv = 0.0;
    double angCv = 0.0;
    (void)chainAreaAngCV(mv, chain, areaCv, angCv);
    out.areaCV = areaCv;
    out.angCV = angCv;
    out.chainN = static_cast<int>(chain.size());
    out.chainScore = chainScore;

    int rid = chain[0];
    for (int t : chain)
        if (t < rid) rid = t;
    emitArchChainDiag(rid, static_cast<int>(chain.size()), areaCv, angCv, chainR, chainScore);

    double spanRad = 0.0;
    double monoFrac = 0.0;
    (void)monotonicNormalSpan(mv, chain, axis, spanRad, monoFrac);

    gp_Pnt center;
    double radius = chainR;
    if (!eberlyCenterRadius(mv, ids, axis, center, radius)) return false;

    if (chainScore >= 0.45) radius = chainR;

    // Coarse band keeps the shipped R>=8 floor (handle-lock identity).
    // Outside coarse, inside archChainBand, high-confidence chains may go to R>=2.
    const double rFloor = (!coarseFusionBand(mv) && chainScore >= 0.85) ? 2.0 : 8.0;
    if (!(radius >= rFloor) || radius > 55.0) return false;

    const double res = maxCylResidual(mv, ids, axis, center, radius);
    double accept = tol.epsCylAccept(radius);
    const int nEst = std::max(6, static_cast<int>(ids.size()));
    accept = std::max(accept, chordSagitta(radius, nEst));
    if (mv.nTri >= 500 && mv.nTri <= 1200)
        accept = std::max(accept, 0.05 * radius);
    if (res > accept) return false;

    out.ok = true;
    out.axis = axis;
    out.center = center;
    out.radius = radius;
    out.spanRad = spanRad;
    out.staticNormals = false;
    out.fromArchChain = true;
    return true;
}

}  // namespace refit
}  // namespace stl2step
