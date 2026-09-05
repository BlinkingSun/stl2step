#include "grade_oracle.hpp"

#include "grade_sha256.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <unordered_set>

namespace grade {
namespace {

constexpr int kMaxIters = 50;

struct Stamp {
    std::vector<int> gen;
    int g = 1;
    explicit Stamp(int n) : gen(static_cast<size_t>(n), 0) {}
    void next() {
        ++g;
        if (g == 0) {
            std::fill(gen.begin(), gen.end(), 0);
            g = 1;
        }
    }
    bool mark(int i) {
        if (gen[static_cast<size_t>(i)] == g) return false;
        gen[static_cast<size_t>(i)] = g;
        return true;
    }
    bool has(int i) const { return gen[static_cast<size_t>(i)] == g; }
};

void uniqueVerts(const Mesh& m, const std::vector<int>& region, std::vector<int>& out, Stamp& st) {
    out.clear();
    st.next();
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        for (int k = 0; k < 3; ++k) {
            if (st.mark(tr.v[k])) out.push_back(tr.v[k]);
        }
    }
}

bool fitPlane(const Mesh& m, const std::vector<int>& region, SurfParams& S) {
    double wsum = 0;
    Vec3 c{};
    Vec3 nmean{};
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        wsum += tr.area;
        c = c + tr.centroid * tr.area;
        nmean = nmean + tr.n * tr.area;
    }
    if (!(wsum > 0.0)) return false;
    c = c * (1.0 / wsum);
    nmean = normalized(nmean);
    double Cov[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        for (int k = 0; k < 3; ++k) {
            const Vec3 d = m.verts[static_cast<size_t>(tr.v[k])] - c;
            const double w = tr.area;
            Cov[0][0] += w * d.x * d.x;
            Cov[0][1] += w * d.x * d.y;
            Cov[0][2] += w * d.x * d.z;
            Cov[1][1] += w * d.y * d.y;
            Cov[1][2] += w * d.y * d.z;
            Cov[2][2] += w * d.z * d.z;
        }
    }
    Cov[1][0] = Cov[0][1];
    Cov[2][0] = Cov[0][2];
    Cov[2][1] = Cov[1][2];
    double eval[3];
    double evec[3][3];
    eigen3(Cov, eval, evec);
    int imin = 0;
    if (eval[1] < eval[imin]) imin = 1;
    if (eval[2] < eval[imin]) imin = 2;
    Vec3 n{evec[0][imin], evec[1][imin], evec[2][imin]};
    n = normalized(n);
    if (dot(n, nmean) < 0.0) n = -n;
    S.cls = SurfClass::Plane;
    S.n = n;
    S.p0 = n * dot(c, n);  // foot of the perpendicular from the origin onto the plane
    return true;
}

bool kasaCircle(const std::vector<Vec3>& pts, const Vec3& origin, const Vec3& u, const Vec3& v,
                Vec3& centre, double& R) {
    const int n = static_cast<int>(pts.size());
    if (n < 3) return false;
    double AtA[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double Atb[3] = {0, 0, 0};
    for (const Vec3& p : pts) {
        const Vec3 d = p - origin;
        const double x = dot(d, u), y = dot(d, v);
        const double z = x * x + y * y;
        const double row[3] = {x, y, 1.0};
        for (int i = 0; i < 3; ++i) {
            Atb[i] += row[i] * z;
            for (int j = 0; j < 3; ++j) AtA[i][j] += row[i] * row[j];
        }
    }
    double A[9], b[3], x[3];
    for (int i = 0; i < 3; ++i) {
        b[i] = Atb[i];
        for (int j = 0; j < 3; ++j) A[i * 3 + j] = AtA[i][j];
    }
    if (!solveN(3, A, b, x)) return false;
    const double cx = x[0] * 0.5, cy = x[1] * 0.5;
    const double rad2 = cx * cx + cy * cy + x[2];
    if (rad2 <= 0.0) return false;
    R = std::sqrt(rad2);
    centre = origin + u * cx + v * cy;
    return true;
}

void frameFromAxis(const Vec3& a, Vec3& u, Vec3& v) {
    const Vec3 an = normalized(a);
    Vec3 tmp = (std::fabs(an.z) < 0.9) ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
    u = normalized(cross(tmp, an));
    v = cross(an, u);
}

bool gaussNewton(const Mesh& m, const std::vector<int>& verts, SurfClass c, SurfParams& S) {
    const int nV = static_cast<int>(verts.size());
    if (nV < paramCount(c) + 1) return false;
    const double q = m.q;
    const double extent = [&]() {
        Vec3 mn = m.verts[static_cast<size_t>(verts[0])];
        Vec3 mx = mn;
        for (int vi : verts) {
            const Vec3& p = m.verts[static_cast<size_t>(vi)];
            mn.x = std::min(mn.x, p.x);
            mn.y = std::min(mn.y, p.y);
            mn.z = std::min(mn.z, p.z);
            mx.x = std::max(mx.x, p.x);
            mx.y = std::max(mx.y, p.y);
            mx.z = std::max(mx.z, p.z);
        }
        return std::max(dist(mn, mx), q);
    }();
    const double dirTol = q / extent;

    auto pack = [&](std::vector<double>& p) {
        p.clear();
        switch (c) {
            case SurfClass::Cylinder:
                p = {S.p0.x, S.p0.y, S.p0.z, S.n.x, S.n.y, S.n.z, S.R};
                break;
            case SurfClass::Cone:
                p = {S.apex.x, S.apex.y, S.apex.z, S.n.x, S.n.y, S.n.z, S.alpha};
                break;
            case SurfClass::Sphere:
                p = {S.p0.x, S.p0.y, S.p0.z, S.R};
                break;
            case SurfClass::Torus:
                p = {S.p0.x, S.p0.y, S.p0.z, S.n.x, S.n.y, S.n.z, S.R, S.r};
                break;
            default:
                break;
        }
    };
    auto unpack = [&](const std::vector<double>& p) {
        switch (c) {
            case SurfClass::Cylinder:
                S.p0 = {p[0], p[1], p[2]};
                S.n = normalized(Vec3{p[3], p[4], p[5]});
                S.p0 = footFromOrigin(S.p0, S.n);
                S.R = std::fabs(p[6]);
                break;
            case SurfClass::Cone:
                S.apex = {p[0], p[1], p[2]};
                S.n = normalized(Vec3{p[3], p[4], p[5]});
                S.alpha = p[6];
                break;
            case SurfClass::Sphere:
                S.p0 = {p[0], p[1], p[2]};
                S.R = std::fabs(p[3]);
                break;
            case SurfClass::Torus:
                S.p0 = {p[0], p[1], p[2]};
                S.n = normalized(Vec3{p[3], p[4], p[5]});
                S.R = std::fabs(p[6]);
                S.r = std::fabs(p[7]);
                break;
            default:
                break;
        }
    };

    std::vector<double> p;
    pack(p);
    const int np = static_cast<int>(p.size());
    if (np == 0) return true;
    std::vector<double> r(static_cast<size_t>(nV)), J(static_cast<size_t>(nV * np));
    std::vector<double> p2 = p;
    for (int it = 0; it < kMaxIters; ++it) {
        unpack(p);
        for (int i = 0; i < nV; ++i)
            r[static_cast<size_t>(i)] = distToSurf(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], S);
        for (int j = 0; j < np; ++j) {
            p2 = p;
            p2[static_cast<size_t>(j)] += q;
            unpack(p2);
            for (int i = 0; i < nV; ++i) {
                const double rp =
                    distToSurf(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], S);
                J[static_cast<size_t>(i * np + j)] = (rp - r[static_cast<size_t>(i)]) / q;
            }
        }
        unpack(p);
        std::vector<double> JtJ(static_cast<size_t>(np * np), 0.0), Jtr(static_cast<size_t>(np), 0.0);
        for (int i = 0; i < nV; ++i) {
            for (int j = 0; j < np; ++j) {
                Jtr[static_cast<size_t>(j)] += J[static_cast<size_t>(i * np + j)] * r[static_cast<size_t>(i)];
                for (int k = 0; k < np; ++k)
                    JtJ[static_cast<size_t>(j * np + k)] +=
                        J[static_cast<size_t>(i * np + j)] * J[static_cast<size_t>(i * np + k)];
            }
        }
        for (int j = 0; j < np; ++j) Jtr[static_cast<size_t>(j)] = -Jtr[static_cast<size_t>(j)];
        std::vector<double> dp(static_cast<size_t>(np), 0.0);
        if (!solveN(np, JtJ.data(), Jtr.data(), dp.data())) break;
        double maxStep = 0;
        for (int j = 0; j < np; ++j) {
            p[static_cast<size_t>(j)] += dp[static_cast<size_t>(j)];
            const bool dirP = (c != SurfClass::Sphere && (j == 3 || j == 4 || j == 5) && c != SurfClass::Plane);
            const double tol = dirP ? dirTol : q;
            maxStep = std::max(maxStep, std::fabs(dp[static_cast<size_t>(j)]) / std::max(tol, q));
        }
        if (maxStep <= 1.0) break;
    }
    unpack(p);
    return true;
}

bool fitCylinder(const Mesh& m, const std::vector<int>& region, const std::vector<int>& verts,
                 SurfParams& S, bool refine) {
    double NNT[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int t : region) {
        const Vec3& n = m.tris[static_cast<size_t>(t)].n;
        NNT[0][0] += n.x * n.x;
        NNT[0][1] += n.x * n.y;
        NNT[0][2] += n.x * n.z;
        NNT[1][1] += n.y * n.y;
        NNT[1][2] += n.y * n.z;
        NNT[2][2] += n.z * n.z;
    }
    NNT[1][0] = NNT[0][1];
    NNT[2][0] = NNT[0][2];
    NNT[2][1] = NNT[1][2];
    double eval[3];
    double evec[3][3];
    eigen3(NNT, eval, evec);
    int imin = 0;
    if (eval[1] < eval[imin]) imin = 1;
    if (eval[2] < eval[imin]) imin = 2;
    Vec3 a{evec[0][imin], evec[1][imin], evec[2][imin]};
    a = normalized(a);
    Vec3 u, v;
    frameFromAxis(a, u, v);
    std::vector<Vec3> pts;
    pts.reserve(verts.size());
    for (int vi : verts) pts.push_back(m.verts[static_cast<size_t>(vi)]);
    Vec3 c{};
    double R = 0;
    if (!kasaCircle(pts, Vec3{}, u, v, c, R)) return false;
    S.cls = SurfClass::Cylinder;
    S.n = a;
    S.p0 = footFromOrigin(c, a);
    S.R = R;
    if (refine) gaussNewton(m, verts, SurfClass::Cylinder, S);
    return S.R > 0.0;
}

bool fitSphere(const Mesh& m, const std::vector<int>& verts, SurfParams& S, bool refine) {
    const int n = static_cast<int>(verts.size());
    if (n < 4) return false;
    double AtA[4][4] = {};
    double Atb[4] = {};
    for (int vi : verts) {
        const Vec3& p = m.verts[static_cast<size_t>(vi)];
        const double row[4] = {2.0 * p.x, 2.0 * p.y, 2.0 * p.z, 1.0};
        const double rhs = norm2(p);
        for (int i = 0; i < 4; ++i) {
            Atb[i] += row[i] * rhs;
            for (int j = 0; j < 4; ++j) AtA[i][j] += row[i] * row[j];
        }
    }
    double A[16], b[4], x[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = Atb[i];
        for (int j = 0; j < 4; ++j) A[i * 4 + j] = AtA[i][j];
    }
    if (!solveN(4, A, b, x)) return false;
    S.cls = SurfClass::Sphere;
    S.p0 = {x[0], x[1], x[2]};
    const double r2 = norm2(S.p0) + x[3];
    if (r2 <= 0.0) return false;
    S.R = std::sqrt(r2);
    if (refine) gaussNewton(m, verts, SurfClass::Sphere, S);
    return S.R > 0.0;
}

bool fitCone(const Mesh& m, const std::vector<int>& region, const std::vector<int>& verts,
             SurfParams& S, bool refine) {
    // Cone of normals: n · â = sin(α) (constant). Axis is the smallest
    // eigenvector of the covariance of normals about their mean — not the
    // stacked-NNT cylinder seed (that axis is ⟂ every n, i.e. a cylinder).
    Vec3 nmean{};
    double wsum = 0;
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        nmean = nmean + tr.n * tr.area;
        wsum += tr.area;
    }
    if (!(wsum > 0.0)) return false;
    nmean = normalized(nmean);
    double Cov[3][3] = {};
    for (int t : region) {
        const Vec3 d = m.tris[static_cast<size_t>(t)].n - nmean;
        Cov[0][0] += d.x * d.x;
        Cov[0][1] += d.x * d.y;
        Cov[0][2] += d.x * d.z;
        Cov[1][1] += d.y * d.y;
        Cov[1][2] += d.y * d.z;
        Cov[2][2] += d.z * d.z;
    }
    Cov[1][0] = Cov[0][1];
    Cov[2][0] = Cov[0][2];
    Cov[2][1] = Cov[1][2];
    double NNT[3][3] = {};
    for (int t : region) {
        const Vec3& n = m.tris[static_cast<size_t>(t)].n;
        NNT[0][0] += n.x * n.x;
        NNT[0][1] += n.x * n.y;
        NNT[0][2] += n.x * n.z;
        NNT[1][1] += n.y * n.y;
        NNT[1][2] += n.y * n.z;
        NNT[2][2] += n.z * n.z;
    }
    NNT[1][0] = NNT[0][1];
    NNT[2][0] = NNT[0][2];
    NNT[2][1] = NNT[1][2];
    auto axisOf = [](const double A[3][3]) {
        double eval[3];
        double evec[3][3];
        eigen3(A, eval, evec);
        int imin = 0;
        if (eval[1] < eval[imin]) imin = 1;
        if (eval[2] < eval[imin]) imin = 2;
        return normalized(Vec3{evec[0][imin], evec[1][imin], evec[2][imin]});
    };
    Vec3 c{};
    for (int vi : verts) c = c + m.verts[static_cast<size_t>(vi)];
    c = c * (1.0 / static_cast<double>(verts.size()));
    auto seedAxis = [&](Vec3 a, SurfParams& T) -> double {
        a = normalized(a);
        if (!(norm2(a) > 0.0)) return 1e300;
        if (dot(a, nmean) < 0.0) a = -a;
        double meanAbsDot = 0;
        for (int t : region)
            meanAbsDot += std::fabs(dot(m.tris[static_cast<size_t>(t)].n, a));
        meanAbsDot /= static_cast<double>(region.size());
        double alpha = std::asin(clamp1(meanAbsDot));
        if (alpha < 0.0) alpha = -alpha;
        double meanAx = 0, meanRho = 0;
        for (int vi : verts) {
            const Vec3 u = m.verts[static_cast<size_t>(vi)] - c;
            const double ax = dot(u, a);
            meanAx += ax;
            meanRho += norm(u - a * ax);
        }
        meanAx /= static_cast<double>(verts.size());
        meanRho /= static_cast<double>(verts.size());
        const double ta = std::tan(std::max(alpha, m.q));
        const double along = (ta > 0.0) ? (meanRho / ta) : 0.0;
        auto evalApex = [&](const Vec3& apex) {
            SurfParams U;
            U.cls = SurfClass::Cone;
            U.n = a;
            U.alpha = alpha;
            U.apex = apex;
            double mx = 0;
            for (int vi : verts) mx = std::max(mx, distToSurf(m.verts[static_cast<size_t>(vi)], U));
            return std::make_pair(mx, U);
        };
        const auto p1 = evalApex(c + a * meanAx - a * along);
        const auto p2 = evalApex(c + a * meanAx + a * along);
        auto p = (p1.first <= p2.first) ? p1 : p2;
        // 1-D search of apex along the axis (extent-scaled; not a tolerance).
        Vec3 mn{1e300, 1e300, 1e300}, mxb{-1e300, -1e300, -1e300};
        for (int vi : verts) {
            const Vec3& q = m.verts[static_cast<size_t>(vi)];
            mn.x = std::min(mn.x, q.x);
            mn.y = std::min(mn.y, q.y);
            mn.z = std::min(mn.z, q.z);
            mxb.x = std::max(mxb.x, q.x);
            mxb.y = std::max(mxb.y, q.y);
            mxb.z = std::max(mxb.z, q.z);
        }
        const double ext = std::max(dist(mn, mxb), m.q);
        const Vec3 A0 = c + a * meanAx;
        double bestT = dot(p.second.apex - A0, a);
        double bestMx = p.first;
        const double span = std::max(ext * 20.0, along * 4.0);
        for (int k = 0; k <= 40; ++k) {
            const double t = -span + (2.0 * span) * (static_cast<double>(k) / 40.0);
            const auto pk = evalApex(A0 + a * t);
            if (pk.first < bestMx) {
                bestMx = pk.first;
                bestT = t;
            }
        }
        T = evalApex(A0 + a * bestT).second;
        return bestMx;
    };
    const Vec3 cands[5] = {axisOf(Cov), axisOf(NNT), Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
    bool any = false;
    double best = 1e300;
    for (const Vec3& ax : cands) {
        SurfParams T;
        const double mx = seedAxis(ax, T);
        if (mx < best) {
            best = mx;
            S = T;
            any = true;
        }
    }
    if (!any) return false;
    if (refine) {
        auto coneMax = [&](const SurfParams& T) {
            double mx = 0;
            for (int vi : verts) mx = std::max(mx, distToSurf(m.verts[static_cast<size_t>(vi)], T));
            return mx;
        };
        const SurfParams seed = S;
        const double r0 = coneMax(S);
        gaussNewton(m, verts, SurfClass::Cone, S);
        S.alpha = std::fabs(S.alpha);
        if (coneMax(S) > r0) S = seed;
    }
    S.alpha = std::fabs(S.alpha);
    return S.alpha > 0.0 && S.alpha < M_PI * 0.5;
}

Vec3 smallestEigenvec(const double A[3][3]) {
    double eval[3];
    double evec[3][3];
    eigen3(A, eval, evec);
    int imin = 0;
    if (eval[1] < eval[imin]) imin = 1;
    if (eval[2] < eval[imin]) imin = 2;
    return normalized(Vec3{evec[0][imin], evec[1][imin], evec[2][imin]});
}

bool seedTorusAxis(const Mesh& m, const std::vector<int>& verts, const Vec3& aIn, SurfParams& S) {
    Vec3 a = normalized(aIn);
    if (!(norm2(a) > 0.0)) return false;
    Vec3 cen{};
    for (int vi : verts) cen = cen + m.verts[static_cast<size_t>(vi)];
    cen = cen * (1.0 / static_cast<double>(verts.size()));
    // Meridional Kåsa: (ρ, axial) is a circle of radius R_min about (R_maj, z0).
    // A 90° rim-blend is a quarter of that circle; Kåsa on the arc is stable,
    // unlike Kåsa on the XY annulus.
    std::vector<Vec3> mer;
    mer.reserve(verts.size());
    for (int vi : verts) {
        const Vec3 w = m.verts[static_cast<size_t>(vi)] - cen;
        const double ax = dot(w, a);
        const double rho = norm(w - a * ax);
        mer.push_back(Vec3{rho, ax, 0});
    }
    Vec3 c2{};
    double rmin = 0;
    if (!kasaCircle(mer, Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, c2, rmin) || !(rmin > 0.0)) {
        double meanRho = 0;
        for (const Vec3& p : mer) meanRho += p.x;
        meanRho /= static_cast<double>(mer.size());
        double rsum = 0;
        for (const Vec3& p : mer)
            rsum += std::sqrt((p.x - meanRho) * (p.x - meanRho) + p.y * p.y);
        S.cls = SurfClass::Torus;
        S.n = a;
        S.p0 = cen;
        S.R = meanRho;
        S.r = rsum / static_cast<double>(mer.size());
        return S.R > 0.0 && S.r > 0.0;
    }
    S.cls = SurfClass::Torus;
    S.n = a;
    S.p0 = cen + a * c2.y;
    S.R = std::fabs(c2.x);
    S.r = rmin;
    return S.R > 0.0 && S.r > 0.0;
}

bool fitTorus(const Mesh& m, const std::vector<int>& region, const std::vector<int>& verts,
              SurfParams& S, bool refine) {
    // Two algebraic axis seeds: (1) smallest eigenvector of the point
    // covariance (a torus ring is planar in the major-circle plane — this is
    // the axis), (2) stacked-normals smallest (cylinder-like, weaker on a
    // 90° fillet whose normals include +axis). Keep the seed with the smaller
    // max vertex residual, then Gauss–Newton.
    Vec3 cen{};
    for (int vi : verts) cen = cen + m.verts[static_cast<size_t>(vi)];
    cen = cen * (1.0 / static_cast<double>(verts.size()));
    double Pcov[3][3] = {};
    for (int vi : verts) {
        const Vec3 d = m.verts[static_cast<size_t>(vi)] - cen;
        Pcov[0][0] += d.x * d.x;
        Pcov[0][1] += d.x * d.y;
        Pcov[0][2] += d.x * d.z;
        Pcov[1][1] += d.y * d.y;
        Pcov[1][2] += d.y * d.z;
        Pcov[2][2] += d.z * d.z;
    }
    Pcov[1][0] = Pcov[0][1];
    Pcov[2][0] = Pcov[0][2];
    Pcov[2][1] = Pcov[1][2];
    double NNT[3][3] = {};
    for (int t : region) {
        const Vec3& n = m.tris[static_cast<size_t>(t)].n;
        NNT[0][0] += n.x * n.x;
        NNT[0][1] += n.x * n.y;
        NNT[0][2] += n.x * n.z;
        NNT[1][1] += n.y * n.y;
        NNT[1][2] += n.y * n.z;
        NNT[2][2] += n.z * n.z;
    }
    NNT[1][0] = NNT[0][1];
    NNT[2][0] = NNT[0][2];
    NNT[2][1] = NNT[1][2];
    const Vec3 aPts = smallestEigenvec(Pcov);
    const Vec3 aNrm = smallestEigenvec(NNT);
    const Vec3 cands[5] = {aPts, aNrm, Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
    auto maxR = [&](const SurfParams& T) {
        double mx = 0;
        for (int vi : verts) mx = std::max(mx, distToSurf(m.verts[static_cast<size_t>(vi)], T));
        return mx;
    };
    bool any = false;
    double best = 1e300;
    for (const Vec3& ax : cands) {
        SurfParams T;
        if (!seedTorusAxis(m, verts, ax, T)) continue;
        const double mx = maxR(T);
        if (!any || mx < best) {
            any = true;
            best = mx;
            S = T;
        }
    }
    if (!any) return false;
    if (refine) {
        const SurfParams seed = S;
        const double r0 = maxR(S);
        gaussNewton(m, verts, SurfClass::Torus, S);
        if (maxR(S) > r0) S = seed;
    }
    return S.R > 0.0 && S.r > 0.0;
}

}  // namespace

double distToSurf(const Vec3& v, const SurfParams& S) {
    switch (S.cls) {
        case SurfClass::Plane:
            return std::fabs(dot(v - S.p0, S.n));
        case SurfClass::Cylinder: {
            const Vec3 w = v - S.p0;
            const double ax = dot(w, S.n);
            const double rho = norm(w - S.n * ax);
            return std::fabs(rho - S.R);
        }
        case SurfClass::Cone: {
            const double a = std::fabs(S.alpha);
            const Vec3 u = v - S.apex;
            const double ax = dot(u, S.n);
            const double r = norm(u - S.n * ax);
            return std::fabs(r * std::cos(a) - ax * std::sin(a));
        }
        case SurfClass::Sphere:
            return std::fabs(norm(v - S.p0) - S.R);
        case SurfClass::Torus: {
            const Vec3 u = v - S.p0;
            const double ax = dot(u, S.n);
            const double r = norm(u - S.n * ax);
            return std::fabs(std::sqrt((r - S.R) * (r - S.R) + ax * ax) - S.r);
        }
        default:
            return 1e300;
    }
}

Vec3 normalAt(const SurfParams& S, const Vec3& p) {
    switch (S.cls) {
        case SurfClass::Plane:
            return S.n;
        case SurfClass::Cylinder: {
            const Vec3 w = p - S.p0;
            const double ax = dot(w, S.n);
            return normalized(w - S.n * ax);
        }
        case SurfClass::Cone: {
            const Vec3 u = p - S.apex;
            const double ax = dot(u, S.n);
            const Vec3 rad = u - S.n * ax;
            const Vec3 rn = normalized(rad);
            const double a = std::fabs(S.alpha);
            return normalized(rn * std::cos(a) - S.n * std::sin(a));
        }
        case SurfClass::Sphere:
            return normalized(p - S.p0);
        case SurfClass::Torus: {
            const Vec3 u = p - S.p0;
            const double ax = dot(u, S.n);
            const Vec3 rad = u - S.n * ax;
            const Vec3 onMaj = S.p0 + normalized(rad) * S.R;
            return normalized(p - onMaj);
        }
        default:
            return {0, 0, 1};
    }
}

bool fitClassEx(const Mesh& m, const std::vector<int>& region, SurfClass c, SurfParams& S,
                bool refine) {
    if (region.empty()) return false;
    Stamp st(static_cast<int>(m.verts.size()));
    std::vector<int> verts;
    uniqueVerts(m, region, verts, st);
    S = SurfParams{};
    S.cls = c;
    switch (c) {
        case SurfClass::Plane:
            return fitPlane(m, region, S);
        case SurfClass::Cylinder:
            return fitCylinder(m, region, verts, S, refine);
        case SurfClass::Cone:
            return fitCone(m, region, verts, S, refine);
        case SurfClass::Sphere:
            return fitSphere(m, verts, S, refine);
        case SurfClass::Torus:
            return fitTorus(m, region, verts, S, refine);
        default:
            return false;
    }
}

bool fitClass(const Mesh& m, const std::vector<int>& region, SurfClass c, SurfParams& S) {
    return fitClassEx(m, region, c, S, true);
}

int distinctNormalClusters(const Mesh& m, const std::vector<int>& triIds) {
    std::vector<Vec3> reps;
    reps.reserve(triIds.size());
    for (int t : triIds) {
        const Vec3 n = m.tris[static_cast<size_t>(t)].n;
        bool found = false;
        for (const Vec3& r : reps) {
            if (angleUnit(n, r) <= m.tris[static_cast<size_t>(t)].thetaQ) {
                found = true;
                break;
            }
        }
        if (!found) reps.push_back(n);
    }
    return static_cast<int>(reps.size());
}

bool surroundsAxis(const Mesh& m, const std::vector<int>& triIds, Vec3 axis) {
    axis = normalized(axis);
    Vec3 uu, vv;
    frameFromAxis(axis, uu, vv);
    std::vector<double> ang;
    ang.reserve(triIds.size());
    for (int t : triIds) {
        const Vec3 n = m.tris[static_cast<size_t>(t)].n;
        const Vec3 p = n - axis * dot(n, axis);
        if (!(norm2(p) > 0.0)) continue;
        const Vec3 pn = normalized(p);
        ang.push_back(std::atan2(dot(pn, vv), dot(pn, uu)));
    }
    if (static_cast<int>(ang.size()) < 3) return false;
    std::sort(ang.begin(), ang.end());
    double maxGap = 0;
    for (size_t i = 1; i < ang.size(); ++i) maxGap = std::max(maxGap, ang[i] - ang[i - 1]);
    maxGap = std::max(maxGap, (ang.front() + 2.0 * M_PI) - ang.back());
    return maxGap <= M_PI;
}

bool cylinderNormalsOk(const Mesh& m, const std::vector<int>& triIds, const SurfParams& S) {
    for (int t : triIds) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        if (std::fabs(dot(tr.n, S.n)) > std::sin(tr.thetaQ)) return false;
    }
    // Closed N<6 prisms (cube walls wrap an axis with 4 face-normals).
    // Applied during growth too — a 2-tri seed does not surround, so fillets
    // can still start. Partial strips stay.
    if (surroundsAxis(m, triIds, S.n) &&
        distinctNormalClusters(m, triIds) < paramCount(SurfClass::Cylinder) + 1)
        return false;
    return true;
}

bool admits(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S) {
    // Growth admission: vertex-on-surface (and the plane normal clause). Size
    // gates are applied only when a region is emitted (certifies).
    Stamp st(static_cast<int>(m.verts.size()));
    std::vector<int> verts;
    uniqueVerts(m, region, verts, st);
    for (int vi : verts) {
        if (distToSurf(m.verts[static_cast<size_t>(vi)], S) > m.tau) return false;
    }
    if (c == SurfClass::Plane) {
        for (int t : region) {
            const Tri& tr = m.tris[static_cast<size_t>(t)];
            if (angleUnit(tr.n, S.n) > tr.thetaQ) return false;
        }
    } else if (c == SurfClass::Cylinder) {
        if (!cylinderNormalsOk(m, region, S)) return false;
    }
    return true;
}

bool certifies(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S,
               double* maxResidOut, double* maxNDevOut) {
    if (static_cast<int>(region.size()) < 2) return false;
    Stamp st(static_cast<int>(m.verts.size()));
    std::vector<int> verts;
    uniqueVerts(m, region, verts, st);
    if (static_cast<int>(verts.size()) < paramCount(c) + 1) return false;
    if (!admits(m, region, c, S)) return false;
    double maxR = 0;
    for (int vi : verts) maxR = std::max(maxR, distToSurf(m.verts[static_cast<size_t>(vi)], S));
    double maxNd = 0;
    if (c != SurfClass::Plane) {
        for (int t : region) {
            const Tri& tr = m.tris[static_cast<size_t>(t)];
            const Vec3 ns = normalAt(S, tr.centroid);
            maxNd = std::max(maxNd, angleUnit(tr.n, ns) - tr.thetaQ);
        }
    }
    if (maxResidOut) *maxResidOut = maxR;
    if (maxNDevOut) *maxNDevOut = maxNd;
    return true;
}

namespace {

void bboxOf(const Mesh& m, const std::vector<int>& region, Vec3& mn, Vec3& mx) {
    mn = {1e300, 1e300, 1e300};
    mx = {-1e300, -1e300, -1e300};
    for (int t : region) {
        for (int k = 0; k < 3; ++k) {
            const Vec3& p = m.verts[static_cast<size_t>(m.tris[static_cast<size_t>(t)].v[k])];
            mn.x = std::min(mn.x, p.x);
            mn.y = std::min(mn.y, p.y);
            mn.z = std::min(mn.z, p.z);
            mx.x = std::max(mx.x, p.x);
            mx.y = std::max(mx.y, p.y);
            mx.z = std::max(mx.z, p.z);
        }
    }
}

int minVertOf(const Mesh& m, const std::vector<int>& region) {
    int mv = m.verts.empty() ? 0 : static_cast<int>(m.verts.size());
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        mv = std::min(mv, std::min(tr.v[0], std::min(tr.v[1], tr.v[2])));
    }
    return mv;
}

void fillOracleStats(const Mesh& m, Oracle& o) {
    Stamp st(static_cast<int>(m.verts.size()));
    uniqueVerts(m, o.tris, o.verts, st);
    o.w = regionArea(m, o.tris);
    certifies(m, o.tris, o.cls, o.S, &o.maxResid, &o.maxNormalDev);
    bboxOf(m, o.tris, o.bboxMin, o.bboxMax);
    o.minVertIndex = minVertOf(m, o.tris);
    // edge-connected pieces
    std::vector<char> in(m.tris.size(), 0);
    for (int t : o.tris) in[static_cast<size_t>(t)] = 1;
    std::vector<char> seen(m.tris.size(), 0);
    std::vector<std::vector<int>> pieces;
    for (int t : o.tris) {
        if (seen[static_cast<size_t>(t)]) continue;
        std::vector<int> stack{t}, piece;
        seen[static_cast<size_t>(t)] = 1;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            piece.push_back(u);
            for (int n : m.adj[static_cast<size_t>(u)]) {
                if (in[static_cast<size_t>(n)] && !seen[static_cast<size_t>(n)]) {
                    seen[static_cast<size_t>(n)] = 1;
                    stack.push_back(n);
                }
            }
        }
        pieces.push_back(std::move(piece));
    }
    o.oraclePieces = static_cast<int>(pieces.size());
    o.oracleMinSeparation = 0;
    if (pieces.size() >= 2) {
        std::sort(pieces.begin(), pieces.end(),
                  [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });
        // BFS from piece 0 through non-member tris to piece 1
        std::vector<int> dist(m.tris.size(), -1);
        std::vector<int> q;
        for (int t : pieces[0]) {
            dist[static_cast<size_t>(t)] = 0;
            q.push_back(t);
        }
        std::unordered_set<int> goal(pieces[1].begin(), pieces[1].end());
        int found = -1;
        size_t qh = 0;
        while (qh < q.size()) {
            const int u = q[qh++];
            if (goal.count(u) && dist[static_cast<size_t>(u)] > 0) {
                found = dist[static_cast<size_t>(u)];
                break;
            }
            for (int n : m.adj[static_cast<size_t>(u)]) {
                if (dist[static_cast<size_t>(n)] >= 0) continue;
                const int step = in[static_cast<size_t>(n)] ? 0 : 1;
                dist[static_cast<size_t>(n)] = dist[static_cast<size_t>(u)] + step;
                q.push_back(n);
            }
        }
        o.oracleMinSeparation = found < 0 ? 0 : found;
    }
}

bool evenOdd(const std::vector<Vec3>& poly, double x, double y) {
    bool in = false;
    const int n = static_cast<int>(poly.size());
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const double xi = poly[static_cast<size_t>(i)].x, yi = poly[static_cast<size_t>(i)].y;
        const double xj = poly[static_cast<size_t>(j)].x, yj = poly[static_cast<size_t>(j)].y;
        const bool hit = ((yi > y) != (yj > y)) &&
                         (x < (xj - xi) * (y - yi) / ((yj - yi) == 0.0 ? 1.0 : (yj - yi)) + xi);
        if (hit) in = !in;
    }
    return in;
}

void partitionPlanes(const Mesh& m, std::vector<Oracle>& oracles) {
    std::vector<Oracle> extra;
    for (Oracle& o : oracles) {
        if (o.cls != SurfClass::Plane) continue;
        Vec3 n = o.S.n;
        Vec3 u, v;
        frameFromAxis(n, u, v);
        std::vector<char> in(m.tris.size(), 0);
        for (int t : o.tris) in[static_cast<size_t>(t)] = 1;
        std::vector<char> seen(m.tris.size(), 0);
        std::vector<std::vector<int>> pieces;
        for (int t : o.tris) {
            if (seen[static_cast<size_t>(t)]) continue;
            std::vector<int> stack{t}, piece;
            seen[static_cast<size_t>(t)] = 1;
            while (!stack.empty()) {
                const int cur = stack.back();
                stack.pop_back();
                piece.push_back(cur);
                for (int nb : m.adj[static_cast<size_t>(cur)]) {
                    if (in[static_cast<size_t>(nb)] && !seen[static_cast<size_t>(nb)]) {
                        seen[static_cast<size_t>(nb)] = 1;
                        stack.push_back(nb);
                    }
                }
            }
            pieces.push_back(std::move(piece));
        }
        if (pieces.size() <= 1) continue;
        // Re-unite Q into P iff Q's centroid is inside a P inner loop.
        // Approximate inner loop: for each piece, collect unused-once projected edges,
        // walk cycles, largest |signed area| is outer, others inner.
        struct PieceInfo {
            std::vector<int> tris;
            Vec3 c2{};
            std::vector<std::vector<Vec3>> inners;
            bool taken = false;
        };
        std::vector<PieceInfo> info;
        for (auto& piece : pieces) {
            PieceInfo pi;
            pi.tris = piece;
            double w = 0;
            Vec3 c{};
            for (int t : piece) {
                const Tri& tr = m.tris[static_cast<size_t>(t)];
                const Vec3 p = tr.centroid;
                const Vec3 d = p - o.S.p0;
                const Vec3 p2{dot(d, u), dot(d, v), 0};
                c = c + p2 * tr.area;
                w += tr.area;
            }
            pi.c2 = w > 0.0 ? c * (1.0 / w) : Vec3{};
            // boundary loops
            std::unordered_map<uint64_t, int> use;
            auto ek = [](int a, int b) {
                if (a > b) std::swap(a, b);
                return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
            };
            for (int t : piece) {
                const Tri& tr = m.tris[static_cast<size_t>(t)];
                for (int s = 0; s < 3; ++s) ++use[ek(tr.v[s], tr.v[(s + 1) % 3])];
            }
            std::unordered_map<int, std::vector<int>> nexts;
            for (int t : piece) {
                const Tri& tr = m.tris[static_cast<size_t>(t)];
                for (int s = 0; s < 3; ++s) {
                    const int a = tr.v[s], b = tr.v[(s + 1) % 3];
                    if (use[ek(a, b)] == 1) nexts[a].push_back(b);
                }
            }
            std::unordered_set<uint64_t> walked;
            std::vector<std::pair<double, std::vector<Vec3>>> loops;
            for (auto& kv : nexts) {
                int start = kv.first;
                if (kv.second.empty()) continue;
                std::vector<Vec3> loop;
                int cur = start;
                int guard = 0;
                while (guard++ < 10000) {
                    const Vec3 d = m.verts[static_cast<size_t>(cur)] - o.S.p0;
                    loop.push_back(Vec3{dot(d, u), dot(d, v), 0});
                    auto it = nexts.find(cur);
                    if (it == nexts.end() || it->second.empty()) break;
                    int nxt = -1;
                    for (int cand : it->second) {
                        const uint64_t e = ek(cur, cand);
                        if (!walked.count(e)) {
                            nxt = cand;
                            walked.insert(e);
                            break;
                        }
                    }
                    if (nxt < 0) break;
                    cur = nxt;
                    if (cur == start) break;
                }
                if (loop.size() < 3) continue;
                double area = 0;
                for (size_t i = 0; i < loop.size(); ++i) {
                    const Vec3& a = loop[i];
                    const Vec3& b = loop[(i + 1) % loop.size()];
                    area += a.x * b.y - b.x * a.y;
                }
                loops.push_back({area * 0.5, std::move(loop)});
            }
            if (!loops.empty()) {
                std::sort(loops.begin(), loops.end(),
                          [](const auto& A, const auto& B) { return std::fabs(A.first) > std::fabs(B.first); });
                for (size_t i = 1; i < loops.size(); ++i) pi.inners.push_back(std::move(loops[i].second));
            }
            info.push_back(std::move(pi));
        }
        std::vector<int> parent(info.size(), -1);
        for (size_t q = 0; q < info.size(); ++q) {
            for (size_t p = 0; p < info.size(); ++p) {
                if (p == q) continue;
                for (const auto& inner : info[p].inners) {
                    if (evenOdd(inner, info[q].c2.x, info[q].c2.y)) {
                        parent[q] = static_cast<int>(p);
                        break;
                    }
                }
                if (parent[q] >= 0) break;
            }
        }
        std::vector<std::vector<int>> groups(info.size());
        std::vector<int> roots;
        for (size_t i = 0; i < info.size(); ++i) {
            if (parent[i] < 0) {
                roots.push_back(static_cast<int>(i));
                groups[i] = info[i].tris;
            }
        }
        for (size_t i = 0; i < info.size(); ++i) {
            if (parent[i] >= 0) {
                int r = parent[i];
                while (parent[static_cast<size_t>(r)] >= 0) r = parent[static_cast<size_t>(r)];
                groups[static_cast<size_t>(r)].insert(groups[static_cast<size_t>(r)].end(),
                                                      info[i].tris.begin(), info[i].tris.end());
            }
        }
        bool first = true;
        for (int r : roots) {
            if (first) {
                o.tris = groups[static_cast<size_t>(r)];
                fitClass(m, o.tris, SurfClass::Plane, o.S);
                first = false;
            } else {
                Oracle nO = o;
                nO.tris = groups[static_cast<size_t>(r)];
                fitClass(m, nO.tris, SurfClass::Plane, nO.S);
                extra.push_back(std::move(nO));
            }
        }
    }
    oracles.insert(oracles.end(), extra.begin(), extra.end());
}

std::string canonString(const Oracle& o, double q) {
    const int k = printDigits(q);
    char buf[512];
    auto pr = [&](double x) {
        const double r = roundToQ(x, q);
        char b[64];
        std::snprintf(b, sizeof b, "%.*f", k, r);
        return std::string(b);
    };
    Vec3 dir = o.S.n;
    canonDir(dir);
    switch (o.cls) {
        case SurfClass::Plane: {
            Vec3 p0 = o.S.p0;
            std::snprintf(buf, sizeof buf, "plane|%s|%s|%s|%s|%s|%s", pr(dir.x).c_str(), pr(dir.y).c_str(),
                          pr(dir.z).c_str(), pr(p0.x).c_str(), pr(p0.y).c_str(), pr(p0.z).c_str());
            break;
        }
        case SurfClass::Cylinder: {
            Vec3 c = footFromOrigin(o.S.p0, dir);
            std::snprintf(buf, sizeof buf, "cylinder|%s|%s|%s|%s|%s|%s|%s", pr(dir.x).c_str(),
                          pr(dir.y).c_str(), pr(dir.z).c_str(), pr(c.x).c_str(), pr(c.y).c_str(),
                          pr(c.z).c_str(), pr(o.S.R).c_str());
            break;
        }
        case SurfClass::Cone: {
            std::snprintf(buf, sizeof buf, "cone|%s|%s|%s|%s|%s|%s|%s", pr(dir.x).c_str(), pr(dir.y).c_str(),
                          pr(dir.z).c_str(), pr(o.S.apex.x).c_str(), pr(o.S.apex.y).c_str(),
                          pr(o.S.apex.z).c_str(), pr(std::fabs(o.S.alpha)).c_str());
            break;
        }
        case SurfClass::Sphere: {
            std::snprintf(buf, sizeof buf, "sphere|%s|%s|%s|%s", pr(o.S.p0.x).c_str(), pr(o.S.p0.y).c_str(),
                          pr(o.S.p0.z).c_str(), pr(o.S.R).c_str());
            break;
        }
        case SurfClass::Torus: {
            std::snprintf(buf, sizeof buf, "torus|%s|%s|%s|%s|%s|%s|%s|%s", pr(dir.x).c_str(),
                          pr(dir.y).c_str(), pr(dir.z).c_str(), pr(o.S.p0.x).c_str(), pr(o.S.p0.y).c_str(),
                          pr(o.S.p0.z).c_str(), pr(o.S.R).c_str(), pr(o.S.r).c_str());
            break;
        }
        default:
            std::snprintf(buf, sizeof buf, "other");
            break;
    }
    return buf;
}

bool growOnce(const Mesh& m, std::vector<char>& claimed, SurfClass c, std::vector<int>& R,
              SurfParams& S, bool reverse) {
    bool changed = false;
    std::vector<int> cand;
    std::unordered_set<int> inR(R.begin(), R.end());
    for (int t : R) {
        for (int n : m.adj[static_cast<size_t>(t)]) {
            if (claimed[static_cast<size_t>(n)] || inR.count(n)) continue;
            cand.push_back(n);
        }
    }
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    if (reverse) std::reverse(cand.begin(), cand.end());
    for (int u : cand) {
        std::vector<int> R2 = R;
        R2.push_back(u);
        SurfParams S2;
        if (!fitClassEx(m, R2, c, S2, false)) continue;
        if (admits(m, R2, c, S2)) {
            R.swap(R2);
            S = S2;
            changed = true;
        }
    }
    return changed;
}

bool tryGrow(const Mesh& m, std::vector<char>& claimed, SurfClass c, int seed, bool reverse,
             Oracle& out, std::vector<char>& skipComp) {
    if (claimed[static_cast<size_t>(seed)]) return false;
    if (skipComp[static_cast<size_t>(seed)]) return false;

    auto finish = [&](std::vector<int>& R, SurfParams& S) -> bool {
        while (growOnce(m, claimed, c, R, S, reverse)) {
        }
        if (!fitClassEx(m, R, c, S, true) || !certifies(m, R, c, S)) return false;
        if (c != SurfClass::Plane) {
            SurfParams P;
            if (fitPlane(m, R, P) && certifies(m, R, SurfClass::Plane, P)) return false;
        }
        // N>=6 is the cylinder/prism law (AGENTS.md): a closed prism with
        // fewer than d_C+1 triangles cannot be a cylinder. Cone is gated
        // only by certifies (|V| >= d_C+1, |R| >= 2).
        if (c == SurfClass::Cylinder) {
            if (static_cast<int>(R.size()) < paramCount(SurfClass::Cylinder) + 1)
                return false;
            // S02 cube-face circumcylinders have R larger than the whole
            // mesh (R=53 on a 20 mm cube). A cylinder the mesh can round
            // cannot exceed the mesh's own bounding diagonal.
            if (!m.verts.empty()) {
                Vec3 mn = m.verts[0], mx = m.verts[0];
                for (const Vec3& p : m.verts) {
                    mn.x = std::min(mn.x, p.x);
                    mn.y = std::min(mn.y, p.y);
                    mn.z = std::min(mn.z, p.z);
                    mx.x = std::max(mx.x, p.x);
                    mx.y = std::max(mx.y, p.y);
                    mx.z = std::max(mx.z, p.z);
                }
                if (S.R > dist(mn, mx)) return false;
            }
        }
        if (c == SurfClass::Torus && !m.verts.empty()) {
            Vec3 mn = m.verts[0], mx = m.verts[0];
            for (const Vec3& p : m.verts) {
                mn.x = std::min(mn.x, p.x);
                mn.y = std::min(mn.y, p.y);
                mn.z = std::min(mn.z, p.z);
                mx.x = std::max(mx.x, p.x);
                mx.y = std::max(mx.y, p.y);
                mx.z = std::max(mx.z, p.z);
            }
            const double diag = dist(mn, mx);
            if (S.R > diag || S.r > diag) return false;
        }
        // Torus is tried before sphere so a rim-blend is not eaten as
        // osculating-sphere patches, but a true sphere (S02 corners)
        // also certifies as a degenerate torus — refuse, leave for sphere.
        if (c == SurfClass::Torus) {
            SurfParams Sph;
            if (fitClassEx(m, R, SurfClass::Sphere, Sph, true) &&
                certifies(m, R, SurfClass::Sphere, Sph))
                return false;
        }
        out.cls = c;
        out.S = S;
        out.tris = std::move(R);
        return true;
    };

    if (c == SurfClass::Plane) {
        std::vector<int> R{seed};
        SurfParams S;
        fitClass(m, R, c, S);
        return finish(R, S);
    }

    // Curved classes: a single triangle cannot determine the surface. Pair the
    // seed with each unclaimed neighbour (ascending, or reversed).
    std::vector<int> nbrs = m.adj[static_cast<size_t>(seed)];
    std::sort(nbrs.begin(), nbrs.end());
    if (reverse) std::reverse(nbrs.begin(), nbrs.end());
    for (int u : nbrs) {
        if (claimed[static_cast<size_t>(u)]) continue;
        if (c == SurfClass::Cylinder || c == SurfClass::Cone || c == SurfClass::Torus) {
            const double ang = angleUnit(m.tris[static_cast<size_t>(seed)].n,
                                         m.tris[static_cast<size_t>(u)].n);
            const double th = std::max(m.tris[static_cast<size_t>(seed)].thetaQ,
                                       m.tris[static_cast<size_t>(u)].thetaQ);
            if (ang <= th) continue;  // coplanar pair: a plane, not a curve seed
        }
        std::vector<int> R{seed, u};
        SurfParams S;
        if (!fitClass(m, R, c, S) || !admits(m, R, c, S)) continue;
        if (finish(R, S)) return true;
    }
    // Cone / torus: a local patch is under-determined and looks like a
    // cylinder/sphere. Take the whole unclaimed connected component (the
    // S04 512-tri band, the S03 drafted hole), then fit. finish() still
    // requires certifies — a freeform island fails and stays unclaimed.
    if (c == SurfClass::Cone || c == SurfClass::Torus) {
        std::vector<int> R{seed};
        std::unordered_set<int> inR{seed};
        std::vector<int> stack{seed};
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int n : m.adj[static_cast<size_t>(t)]) {
                if (claimed[static_cast<size_t>(n)] || inR.count(n)) continue;
                inR.insert(n);
                R.push_back(n);
                stack.push_back(n);
            }
        }
        std::sort(R.begin(), R.end());
        SurfParams S;
        const bool fitOk = fitClass(m, R, c, S);
        const bool adOk = fitOk && admits(m, R, c, S);
        if (fitOk && adOk && finish(R, S)) return true;
        for (int t : R) skipComp[static_cast<size_t>(t)] = 1;
    }
    // Sphere: also try the seed plus any two neighbours (four+ vertices).
    if (c == SurfClass::Sphere) {
        std::vector<int> R{seed};
        SurfParams S;
        fitClass(m, R, c, S);
        if (finish(R, S)) return true;
    }
    return false;
}

}  // namespace

void assignFeatureIds(const Mesh& m, OracleSet& set) {
    std::unordered_map<std::string, std::vector<int>> coll;
    for (int i = 0; i < static_cast<int>(set.oracles.size()); ++i) {
        Oracle& o = set.oracles[static_cast<size_t>(i)];
        Sha256 sha;
        const std::string canon = canonString(o, m.q);
        sha.update(canon);
        o.featureId = std::string(className(o.cls)) + ":" + sha.hex12();
        coll[o.featureId].push_back(i);
    }
    for (auto& kv : coll) {
        if (kv.second.size() < 2) continue;
        std::vector<int>& ids = kv.second;
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            const Oracle& A = set.oracles[static_cast<size_t>(a)];
            const Oracle& B = set.oracles[static_cast<size_t>(b)];
            if (A.w != B.w) return A.w > B.w;
            return A.minVertIndex < B.minVertIndex;
        });
        for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
            set.oracles[static_cast<size_t>(ids[static_cast<size_t>(i)])].featureId +=
                "#" + std::to_string(i);
        }
    }
}

void buildOracle(const Mesh& m, OracleSet& out, bool reverseSeeds) {
    out = OracleSet{};
    out.owner.assign(m.tris.size(), -1);
    std::vector<char> claimed(m.tris.size(), 0);
    // Sphere is after torus: a torus band's osculating-sphere patches (S04)
    // must not consume the band before the torus class sees it. A true
    // sphere does not certify as a torus (degenerate Rmaj/Rmin) and is
    // still claimed next. Plane remains first (S04 boss-top class order).
    const SurfClass order[5] = {SurfClass::Plane, SurfClass::Cylinder, SurfClass::Cone,
                                SurfClass::Torus, SurfClass::Sphere};

    std::vector<Oracle> heldPlanes;  // |R|==2 tessellation quads; curved classes may claim them

    auto seedList = [&]() {
        std::vector<int> s;
        s.reserve(m.tris.size());
        if (!reverseSeeds) {
            for (int i = 0; i < static_cast<int>(m.tris.size()); ++i) s.push_back(i);
        } else {
            for (int i = static_cast<int>(m.tris.size()) - 1; i >= 0; --i) s.push_back(i);
        }
        return s;
    };

    for (SurfClass c : order) {
        std::vector<char> skipComp(m.tris.size(), 0);
        const auto seeds = seedList();
        for (int seed : seeds) {
            if (claimed[static_cast<size_t>(seed)]) continue;
            Oracle o;
            if (!tryGrow(m, claimed, c, seed, reverseSeeds, o, skipComp)) continue;
            if (c == SurfClass::Plane && static_cast<int>(o.tris.size()) == 2) {
                // Tessellated cylinders/spheres are planar quads. Hold a 2-tri
                // plane whose neighbour dihedral is the tessellation step
                // (≤ π / (d_C+1) would be a hex; use the cylinder over-
                // determination count). Cube faces against a fillet also look
                // shallow — cylinderNormalsOk keeps those off the cylinder,
                // and they commit after the curved pass.
                bool sharp = true;
                std::unordered_set<int> inR(o.tris.begin(), o.tris.end());
                const Vec3 pn = o.S.n;
                for (int t : o.tris) {
                    for (int n : m.adj[static_cast<size_t>(t)]) {
                        if (inR.count(n)) continue;
                        if (angleUnit(pn, m.tris[static_cast<size_t>(n)].n) <=
                            2.0 * M_PI / static_cast<double>(paramCount(SurfClass::Cylinder) + 1))
                            sharp = false;
                    }
                }
                if (sharp) {
                    for (int t : o.tris) {
                        claimed[static_cast<size_t>(t)] = 1;
                        out.owner[static_cast<size_t>(t)] = static_cast<int>(out.oracles.size());
                    }
                    o.id = static_cast<int>(out.oracles.size());
                    out.oracles.push_back(std::move(o));
                    continue;
                }
                heldPlanes.push_back(std::move(o));
                continue;
            }
            for (int t : o.tris) {
                claimed[static_cast<size_t>(t)] = 1;
                out.owner[static_cast<size_t>(t)] = static_cast<int>(out.oracles.size());
            }
            o.id = static_cast<int>(out.oracles.size());
            out.oracles.push_back(std::move(o));
        }
        // Claim pass (islands)
        bool claimChanged = true;
        while (claimChanged) {
            claimChanged = false;
            for (int u : seedList()) {
                if (claimed[static_cast<size_t>(u)]) continue;
                int best = -1;
                double bestResid = 1e300;
                SurfParams bestS;
                for (int ri = 0; ri < static_cast<int>(out.oracles.size()); ++ri) {
                    Oracle& Rk = out.oracles[static_cast<size_t>(ri)];
                    if (Rk.cls != c) continue;
                    std::vector<int> R2 = Rk.tris;
                    R2.push_back(u);
                    SurfParams S2;
                    if (!fitClass(m, R2, c, S2)) continue;
                    double maxR = 0;
                    if (!certifies(m, R2, c, S2, &maxR)) continue;
                    if (maxR < bestResid || (maxR == bestResid && (best < 0 || ri < best))) {
                        bestResid = maxR;
                        best = ri;
                        bestS = S2;
                    }
                }
                if (best >= 0) {
                    Oracle& Rk = out.oracles[static_cast<size_t>(best)];
                    Rk.tris.push_back(u);
                    Rk.S = bestS;
                    claimed[static_cast<size_t>(u)] = 1;
                    out.owner[static_cast<size_t>(u)] = best;
                    claimChanged = true;
                }
            }
        }
        // Merge pass
        bool merged = true;
        while (merged) {
            merged = false;
            for (int i = 0; i < static_cast<int>(out.oracles.size()) && !merged; ++i) {
                if (out.oracles[static_cast<size_t>(i)].cls != c) continue;
                for (int j = i + 1; j < static_cast<int>(out.oracles.size()); ++j) {
                    if (out.oracles[static_cast<size_t>(j)].cls != c) continue;
                    std::vector<int> U = out.oracles[static_cast<size_t>(i)].tris;
                    U.insert(U.end(), out.oracles[static_cast<size_t>(j)].tris.begin(),
                             out.oracles[static_cast<size_t>(j)].tris.end());
                    SurfParams S;
                    if (!fitClass(m, U, c, S)) continue;
                    if (!certifies(m, U, c, S)) continue;
                    out.oracles[static_cast<size_t>(i)].tris.swap(U);
                    out.oracles[static_cast<size_t>(i)].S = S;
                    for (int t : out.oracles[static_cast<size_t>(j)].tris)
                        out.owner[static_cast<size_t>(t)] = i;
                    out.oracles.erase(out.oracles.begin() + j);
                    for (int t = 0; t < static_cast<int>(m.tris.size()); ++t) {
                        if (out.owner[static_cast<size_t>(t)] > j) --out.owner[static_cast<size_t>(t)];
                    }
                    for (int k = 0; k < static_cast<int>(out.oracles.size()); ++k)
                        out.oracles[static_cast<size_t>(k)].id = k;
                    merged = true;
                    break;
                }
            }
        }
    }

    // Commit held 2-triangle planes that curved classes did not take.
    for (Oracle& hp : heldPlanes) {
        std::vector<int> left;
        for (int t : hp.tris)
            if (!claimed[static_cast<size_t>(t)]) left.push_back(t);
        if (static_cast<int>(left.size()) < 2) continue;
        SurfParams S;
        if (!fitClass(m, left, SurfClass::Plane, S) || !certifies(m, left, SurfClass::Plane, S))
            continue;
        Oracle o;
        o.cls = SurfClass::Plane;
        o.S = S;
        o.tris = std::move(left);
        for (int t : o.tris) {
            claimed[static_cast<size_t>(t)] = 1;
            out.owner[static_cast<size_t>(t)] = static_cast<int>(out.oracles.size());
        }
        o.id = static_cast<int>(out.oracles.size());
        out.oracles.push_back(std::move(o));
    }

    partitionPlanes(m, out.oracles);
    // Re-number owners after partition
    out.owner.assign(m.tris.size(), -1);
    for (int i = 0; i < static_cast<int>(out.oracles.size()); ++i) {
        out.oracles[static_cast<size_t>(i)].id = i;
        std::sort(out.oracles[static_cast<size_t>(i)].tris.begin(),
                  out.oracles[static_cast<size_t>(i)].tris.end());
        fitClass(m, out.oracles[static_cast<size_t>(i)].tris, out.oracles[static_cast<size_t>(i)].cls,
                 out.oracles[static_cast<size_t>(i)].S);
        for (int t : out.oracles[static_cast<size_t>(i)].tris)
            out.owner[static_cast<size_t>(t)] = i;
        fillOracleStats(m, out.oracles[static_cast<size_t>(i)]);
    }

    // Residue
    std::vector<int> rest;
    for (int t = 0; t < static_cast<int>(m.tris.size()); ++t)
        if (out.owner[static_cast<size_t>(t)] < 0) rest.push_back(t);
    out.residueTris = static_cast<int>(rest.size());
    std::vector<char> seen(m.tris.size(), 0);
    for (int t : rest) {
        if (seen[static_cast<size_t>(t)]) continue;
        ResidueComp rc;
        std::vector<int> stack{t};
        seen[static_cast<size_t>(t)] = 1;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            rc.tris.push_back(u);
            rc.area += m.tris[static_cast<size_t>(u)].area;
            for (int n : m.adj[static_cast<size_t>(u)]) {
                if (out.owner[static_cast<size_t>(n)] < 0 && !seen[static_cast<size_t>(n)]) {
                    seen[static_cast<size_t>(n)] = 1;
                    stack.push_back(n);
                }
            }
        }
        bboxOf(m, rc.tris, rc.bboxMin, rc.bboxMax);
        // ruled: smallest singular vector of stacked normals, then verify
        double NNT[3][3] = {};
        for (int u : rc.tris) {
            const Vec3& n = m.tris[static_cast<size_t>(u)].n;
            NNT[0][0] += n.x * n.x;
            NNT[0][1] += n.x * n.y;
            NNT[0][2] += n.x * n.z;
            NNT[1][1] += n.y * n.y;
            NNT[1][2] += n.y * n.z;
            NNT[2][2] += n.z * n.z;
        }
        NNT[1][0] = NNT[0][1];
        NNT[2][0] = NNT[0][2];
        NNT[2][1] = NNT[1][2];
        double eval[3];
        double evec[3][3];
        eigen3(NNT, eval, evec);
        int imin = 0;
        if (eval[1] < eval[imin]) imin = 1;
        if (eval[2] < eval[imin]) imin = 2;
        Vec3 d{evec[0][imin], evec[1][imin], evec[2][imin]};
        d = normalized(d);
        bool ok = !rc.tris.empty();
        for (int u : rc.tris) {
            const Tri& tr = m.tris[static_cast<size_t>(u)];
            if (std::fabs(dot(tr.n, d)) > std::sin(tr.thetaQ)) {
                ok = false;
                break;
            }
        }
        rc.ruled = ok;
        if (ok) {
            canonDir(d);
            rc.direction = d;
        }
        out.residueArea += rc.area;
        out.residue.push_back(std::move(rc));
    }
    std::sort(out.residue.begin(), out.residue.end(), [](const ResidueComp& a, const ResidueComp& b) {
        if (a.area != b.area) return a.area > b.area;
        const int ia = a.tris.empty() ? 0 : *std::min_element(a.tris.begin(), a.tris.end());
        const int ib = b.tris.empty() ? 0 : *std::min_element(b.tris.begin(), b.tris.end());
        return ia < ib;
    });

    // unprovable singletons: leftover triangles that would certify a class alone
    // if |R|>=2 / |V|>=d+1 were dropped — counted, never oracles.
    for (int t : rest) {
        bool any = false;
        for (SurfClass c : order) {
            std::vector<int> R{t};
            SurfParams S;
            if (!fitClass(m, R, c, S)) continue;
            Stamp st(static_cast<int>(m.verts.size()));
            std::vector<int> verts;
            uniqueVerts(m, R, verts, st);
            bool vertsOk = true;
            for (int vi : verts) {
                if (distToSurf(m.verts[static_cast<size_t>(vi)], S) > m.tau) {
                    vertsOk = false;
                    break;
                }
            }
            if (!vertsOk) continue;
            if (c == SurfClass::Plane) {
                const Tri& tr = m.tris[static_cast<size_t>(t)];
                if (angleUnit(tr.n, S.n) > tr.thetaQ) continue;
            }
            any = true;
            break;
        }
        if (any) ++out.unprovableSingletons;
    }

    assignFeatureIds(m, out);
}

}  // namespace grade
