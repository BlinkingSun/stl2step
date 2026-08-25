// stl2step P1 math primitives. Written in-TU; no Eigen, no GProp_, no GeomConvert_.
// Signatures are frozen in refit_internal.hpp.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

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

}  // namespace refit
}  // namespace stl2step
