// Shared profile certificate (D-fc-1.10). Bodies moved verbatim from refit_torus.cpp.
// SPDX-License-Identifier: MIT

#include "refit_profile_cert.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <gp.hxx>

namespace stl2step {
namespace refit {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

constexpr int kFitIters = 16;

int countColumns(std::vector<double> psi, const std::vector<double>& rho, double tau) {
    const std::size_t n = psi.size();
    if (n == 0) return 0;
    double rhoMax = 0.0;
    for (double r : rho) rhoMax = std::max(rhoMax, r);
    if (!(rhoMax > gp::Resolution())) return 1;
    const double angTol = tau / rhoMax;
    std::sort(psi.begin(), psi.end());
    int cols = 1;
    double ref = psi[0];
    for (std::size_t i = 1; i < n; ++i) {
        if (psi[i] - ref > angTol) {
            ++cols;
            ref = psi[i];
        }
    }
    if (cols > 1 && (psi[0] + kTwoPi) - ref <= angTol) --cols;
    return cols;
}

int countLevels(std::vector<double> z, double tau) {
    if (z.empty()) return 0;
    std::sort(z.begin(), z.end());
    int lv = 1;
    double ref = z[0];
    for (std::size_t i = 1; i < z.size(); ++i) {
        if (z[i] - ref > tau) {
            ++lv;
            ref = z[i];
        }
    }
    return lv;
}

int countProfPts(const std::vector<ProfPt>& P, double tau) {
    std::vector<std::pair<double, double>> reps;
    for (const ProfPt& p : P) {
        bool found = false;
        for (const auto& q : reps)
            if (std::fabs(p.rho - q.first) <= tau && std::fabs(p.z - q.second) <= tau) {
                found = true;
                break;
            }
        if (!found) reps.emplace_back(p.rho, p.z);
    }
    return (int)reps.size();
}

double profileLineMaxResid(const std::vector<ProfPt>& P) {
    const std::size_t n = P.size();
    if (n < 2) return 0.0;
    double sx = 0, sy = 0;
    for (const ProfPt& p : P) { sx += p.rho; sy += p.z; }
    const double mx = sx / (double)n, my = sy / (double)n;
    double sxx = 0, sxy = 0, syy = 0;
    for (const ProfPt& p : P) {
        const double dx = p.rho - mx, dy = p.z - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    const double tr = sxx + syy, det = sxx * syy - sxy * sxy;
    const double disc = std::max(0.0, 0.25 * tr * tr - det);
    const double lmin = 0.5 * tr - std::sqrt(disc);
    double nx = sxy, ny = lmin - sxx;
    double nm = std::hypot(nx, ny);
    if (!(nm > gp::Resolution())) { nx = lmin - syy; ny = sxy; nm = std::hypot(nx, ny); }
    if (!(nm > gp::Resolution())) return 0.0;
    nx /= nm; ny /= nm;
    double m = 0.0;
    for (const ProfPt& p : P)
        m = std::max(m, std::fabs((p.rho - mx) * nx + (p.z - my) * ny));
    return m;
}

bool fitProfileCircle(const std::vector<ProfPt>& P, double tau, double& Rmaj, double& z0,
                      double& Rmin) {
    const std::size_t n = P.size();
    if (n < 3) return false;
    double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sz = 0, Sxz = 0, Syz = 0;
    for (const ProfPt& p : P) {
        const double x = p.rho, y = p.z, zq = x * x + y * y;
        Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y;
        Sz += zq; Sxz += x * zq; Syz += y * zq;
    }
    const double N = (double)n;
    const double a11 = 2.0 * (Sxx - Sx * Sx / N);
    const double a12 = 2.0 * (Sxy - Sx * Sy / N);
    const double a22 = 2.0 * (Syy - Sy * Sy / N);
    const double b1 = Sxz - Sx * Sz / N;
    const double b2 = Syz - Sy * Sz / N;
    const double det = a11 * a22 - a12 * a12;
    if (!(std::fabs(det) > gp::Resolution())) return false;
    double a = (b1 * a22 - b2 * a12) / det;
    double b = (a11 * b2 - a12 * b1) / det;
    double r = 0.0;
    for (const ProfPt& p : P) r += std::hypot(p.rho - a, p.z - b);
    r /= N;
    for (int it = 0; it < kFitIters; ++it) {
        double h11 = 0, h12 = 0, h13 = 0, h22 = 0, h23 = 0, h33 = 0;
        double g1 = 0, g2 = 0, g3 = 0;
        for (const ProfPt& p : P) {
            const double dx = p.rho - a, dy = p.z - b;
            const double d = std::hypot(dx, dy);
            if (!(d > gp::Resolution())) continue;
            const double j1 = -dx / d, j2 = -dy / d, j3 = -1.0;
            const double res = d - r;
            h11 += j1 * j1; h12 += j1 * j2; h13 += j1 * j3;
            h22 += j2 * j2; h23 += j2 * j3; h33 += j3 * j3;
            g1 += j1 * res; g2 += j2 * res; g3 += j3 * res;
        }
        const double m[3][3] = {{h11, h12, h13}, {h12, h22, h23}, {h13, h23, h33}};
        const double d0 = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                          m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                          m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if (!(std::fabs(d0) > gp::Resolution())) break;
        const double rhs[3] = {-g1, -g2, -g3};
        auto solveCol = [&](int col) {
            double c[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) c[i][j] = (j == col) ? rhs[i] : m[i][j];
            return (c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) -
                    c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) +
                    c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])) / d0;
        };
        const double da = solveCol(0), db = solveCol(1), dr = solveCol(2);
        if (!std::isfinite(da) || !std::isfinite(db) || !std::isfinite(dr)) break;
        a += da; b += db; r += dr;
        if (std::hypot(std::hypot(da, db), dr) <= tau) break;
    }
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(r)) return false;
    Rmaj = a; z0 = b; Rmin = r;
    return true;
}

TorusCert certifyTorus(const std::vector<ProfPt>& P, double tau) {
    TorusCert c;
    if (P.size() < 3) { c.why = "profile-underdetermined"; return c; }
    std::vector<double> zs, ps, rs;
    zs.reserve(P.size()); ps.reserve(P.size()); rs.reserve(P.size());
    for (const ProfPt& p : P) { zs.push_back(p.z); ps.push_back(p.psi); rs.push_back(p.rho); }
    c.levels = countLevels(zs, tau);
    c.cols = countColumns(ps, rs, tau);
    c.profPts = countProfPts(P, tau);
    if (c.levels < 3 || c.cols < 2 || c.profPts < 5) {
        c.why = "profile-underdetermined";
        return c;
    }
    c.determinate = true;
    c.lineRes = profileLineMaxResid(P);
    if (c.lineRes <= tau) { c.why = "cone-or-cylinder-or-plane"; return c; }
    if (!fitProfileCircle(P, tau, c.Rmaj, c.z0, c.Rmin)) { c.why = "fit-singular"; return c; }
    if (!(c.Rmaj > tau)) { c.why = "sphere"; return c; }
    if (!(c.Rmin > tau)) { c.why = "degenerate"; return c; }
    if (!(c.Rmaj > c.Rmin)) { c.why = "spindle-torus"; return c; }
    double sum2 = 0.0;
    for (const ProfPt& p : P) {
        const double d = std::fabs(std::hypot(p.rho - c.Rmaj, p.z - c.z0) - c.Rmin);
        c.maxDev = std::max(c.maxDev, d);
        sum2 += d * d;
    }
    c.rmsDev = std::sqrt(sum2 / (double)P.size());
    if (c.maxDev > tau) { c.why = "vertex-off-surface"; return c; }
    c.ok = true;
    c.why = "accept";
    return c;
}

}}  // namespace stl2step::refit
