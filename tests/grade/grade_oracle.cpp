#include "grade_oracle.hpp"

#include "grade_sha256.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <unordered_set>

namespace grade {
namespace {

constexpr int kMaxIters = 50;

// Measurement instrument, default OFF (STL2STEP_GRADE_DIAG). It reads the
// oracle state and prints; it never changes a partition, a tolerance or a
// certificate. Every graded artefact in this report was produced with it off.
bool gdiag() {
    static const bool on = [] {
        const char* e = std::getenv("STL2STEP_GRADE_DIAG");
        return e && e[0] && std::strcmp(e, "0") != 0;
    }();
    return on;
}

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

// SPEC §9 data structure: certifies() and admits() are the innermost loop of
// the oracle (100 % of the grip's sampled stacks). A fresh `Stamp` there is an
// O(|V|) zero-fill plus two heap allocations per certificate — 19 % of the
// grip's samples sat in __bzero under admits(). One generation-stamped scratch
// per call site, grown but never re-zeroed, removes both. The grader is
// single-threaded (§7.4); thread_local keeps it re-entrant anyway. No vertex
// is skipped: the scratch holds exactly the same marks the local Stamp did.
struct Scratch {
    Stamp st{0};
    std::vector<int> verts;
    void ensure(size_t n) {
        if (st.gen.size() < n) {
            st.gen.assign(n, 0);
            st.g = 1;
        }
    }
};

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

double signedDist(const Vec3& v, const SurfParams& S) {
    switch (S.cls) {
        case SurfClass::Plane:
            return dot(v - S.p0, S.n);
        case SurfClass::Cylinder: {
            const Vec3 w = v - S.p0;
            const double ax = dot(w, S.n);
            return norm(w - S.n * ax) - S.R;
        }
        case SurfClass::Cone: {
            const double a = S.alpha;
            const Vec3 u = v - S.apex;
            const double ax = dot(u, S.n);
            const double r = norm(u - S.n * ax);
            return r * std::cos(a) - ax * std::sin(a);
        }
        case SurfClass::Sphere:
            return norm(v - S.p0) - S.R;
        case SurfClass::Torus: {
            const Vec3 u = v - S.p0;
            const double ax = dot(u, S.n);
            const double r = norm(u - S.n * ax);
            return std::sqrt((r - S.R) * (r - S.R) + ax * ax) - S.r;
        }
        default:
            return 1e300;
    }
}

double rssSigned(const Mesh& m, const std::vector<int>& verts, const SurfParams& S) {
    double s = 0;
    for (int vi : verts) {
        const double r = signedDist(m.verts[static_cast<size_t>(vi)], S);
        s += r * r;
    }
    return s;
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

    auto pack = [&](std::vector<double>& p, const SurfParams& T) {
        p.clear();
        switch (c) {
            case SurfClass::Cylinder:
                p = {T.p0.x, T.p0.y, T.p0.z, T.n.x, T.n.y, T.n.z, T.R};
                break;
            case SurfClass::Cone:
                p = {T.apex.x, T.apex.y, T.apex.z, T.n.x, T.n.y, T.n.z, T.alpha};
                break;
            case SurfClass::Sphere:
                p = {T.p0.x, T.p0.y, T.p0.z, T.R};
                break;
            case SurfClass::Torus:
                p = {T.p0.x, T.p0.y, T.p0.z, T.n.x, T.n.y, T.n.z, T.R, T.r};
                break;
            default:
                break;
        }
    };
    auto unpack = [&](const std::vector<double>& p, SurfParams& T) {
        T.cls = c;
        switch (c) {
            case SurfClass::Cylinder:
                T.p0 = {p[0], p[1], p[2]};
                T.n = normalized(Vec3{p[3], p[4], p[5]});
                T.p0 = footFromOrigin(T.p0, T.n);
                T.R = std::fabs(p[6]);
                break;
            case SurfClass::Cone:
                T.apex = {p[0], p[1], p[2]};
                T.n = normalized(Vec3{p[3], p[4], p[5]});
                T.alpha = p[6];
                break;
            case SurfClass::Sphere:
                T.p0 = {p[0], p[1], p[2]};
                T.R = std::fabs(p[3]);
                break;
            case SurfClass::Torus:
                T.p0 = {p[0], p[1], p[2]};
                T.n = normalized(Vec3{p[3], p[4], p[5]});
                T.R = std::fabs(p[6]);
                T.r = std::fabs(p[7]);
                break;
            default:
                break;
        }
    };

    std::vector<double> p;
    pack(p, S);
    const int np = static_cast<int>(p.size());
    if (np == 0) return true;
    std::vector<double> r(static_cast<size_t>(nV)), J(static_cast<size_t>(nV * np));
    std::vector<double> p2 = p;
    double lam = 1e-3;
    double bestRss = rssSigned(m, verts, S);
    SurfParams best = S;
    for (int it = 0; it < kMaxIters; ++it) {
        SurfParams cur;
        unpack(p, cur);
        for (int i = 0; i < nV; ++i)
            r[static_cast<size_t>(i)] =
                signedDist(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], cur);
        const double fd = std::max(q, extent * 1e-8);
        for (int j = 0; j < np; ++j) {
            p2 = p;
            p2[static_cast<size_t>(j)] += fd;
            SurfParams T;
            unpack(p2, T);
            for (int i = 0; i < nV; ++i) {
                const double rp =
                    signedDist(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], T);
                J[static_cast<size_t>(i * np + j)] = (rp - r[static_cast<size_t>(i)]) / fd;
            }
        }
        std::vector<double> JtJ(static_cast<size_t>(np * np), 0.0), Jtr(static_cast<size_t>(np), 0.0);
        for (int i = 0; i < nV; ++i) {
            for (int j = 0; j < np; ++j) {
                Jtr[static_cast<size_t>(j)] += J[static_cast<size_t>(i * np + j)] * r[static_cast<size_t>(i)];
                for (int k = 0; k < np; ++k)
                    JtJ[static_cast<size_t>(j * np + k)] +=
                        J[static_cast<size_t>(i * np + j)] * J[static_cast<size_t>(i * np + k)];
            }
        }
        for (int j = 0; j < np; ++j) {
            Jtr[static_cast<size_t>(j)] = -Jtr[static_cast<size_t>(j)];
            JtJ[static_cast<size_t>(j * np + j)] *= (1.0 + lam);
        }
        std::vector<double> dp(static_cast<size_t>(np), 0.0);
        if (!solveN(np, JtJ.data(), Jtr.data(), dp.data())) {
            lam = std::min(lam * 10.0, 1e12);
            if (lam >= 1e12) break;
            continue;
        }
        std::vector<double> pTry = p;
        double maxStep = 0;
        for (int j = 0; j < np; ++j) {
            pTry[static_cast<size_t>(j)] += dp[static_cast<size_t>(j)];
            const bool dirP =
                (c != SurfClass::Sphere && (j == 3 || j == 4 || j == 5) && c != SurfClass::Plane);
            const double tol = dirP ? dirTol : q;
            maxStep = std::max(maxStep, std::fabs(dp[static_cast<size_t>(j)]) / std::max(tol, q));
        }
        SurfParams T;
        unpack(pTry, T);
        const double rss = rssSigned(m, verts, T);
        if (rss <= bestRss) {
            p.swap(pTry);
            bestRss = rss;
            best = T;
            lam = std::max(lam * 0.1, 1e-12);
            if (maxStep <= 1.0) break;
        } else {
            lam = std::min(lam * 10.0, 1e12);
            if (lam >= 1e12) break;
        }
    }
    S = best;
    if (c == SurfClass::Cone) S.alpha = std::fabs(S.alpha);
    return true;
}

// SPEC §5.2 cone: 6 parameters (apex 3 + axis 2-DOF + α). A 3-component unit
// axis is rank-deficient; the local (u,v) chart is the 2-DOF parameterisation.
bool gaussNewtonCone(const Mesh& m, const std::vector<int>& verts, SurfParams& S) {
    const int nV = static_cast<int>(verts.size());
    if (nV < paramCount(SurfClass::Cone) + 1) return false;
    const double q = m.q;
    Vec3 mn = m.verts[static_cast<size_t>(verts[0])], mx = mn;
    for (int vi : verts) {
        const Vec3& p = m.verts[static_cast<size_t>(vi)];
        mn.x = std::min(mn.x, p.x);
        mn.y = std::min(mn.y, p.y);
        mn.z = std::min(mn.z, p.z);
        mx.x = std::max(mx.x, p.x);
        mx.y = std::max(mx.y, p.y);
        mx.z = std::max(mx.z, p.z);
    }
    const double extent = std::max(dist(mn, mx), q);
    const double dirTol = q / extent;
    Vec3 a = normalized(S.n);
    Vec3 u, v;
    frameFromAxis(a, u, v);
    double p[6] = {S.apex.x, S.apex.y, S.apex.z, 0, 0, S.alpha};
    auto apply = [&](const double* pp, SurfParams& T) {
        T.cls = SurfClass::Cone;
        T.apex = {pp[0], pp[1], pp[2]};
        T.n = normalized(a + u * pp[3] + v * pp[4]);
        T.alpha = pp[5];
    };
    auto rechart = [&]() {
        a = normalized(S.n);
        frameFromAxis(a, u, v);
        p[0] = S.apex.x;
        p[1] = S.apex.y;
        p[2] = S.apex.z;
        p[3] = 0;
        p[4] = 0;
        p[5] = S.alpha;
    };
    rechart();
    double lam = 1e-3;
    double bestRss = rssSigned(m, verts, S);
    SurfParams best = S;
    const int np = 6;
    std::vector<double> r(static_cast<size_t>(nV)), J(static_cast<size_t>(nV * np));
    for (int it = 0; it < kMaxIters; ++it) {
        SurfParams cur;
        apply(p, cur);
        for (int i = 0; i < nV; ++i)
            r[static_cast<size_t>(i)] =
                signedDist(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], cur);
        const double fd = std::max(q, extent * 1e-8);
        for (int j = 0; j < np; ++j) {
            double p2[6];
            for (int k = 0; k < np; ++k) p2[k] = p[k];
            p2[j] += fd;
            SurfParams T;
            apply(p2, T);
            for (int i = 0; i < nV; ++i) {
                const double rp =
                    signedDist(m.verts[static_cast<size_t>(verts[static_cast<size_t>(i)])], T);
                J[static_cast<size_t>(i * np + j)] = (rp - r[static_cast<size_t>(i)]) / fd;
            }
        }
        std::vector<double> JtJ(static_cast<size_t>(np * np), 0.0), Jtr(static_cast<size_t>(np), 0.0);
        for (int i = 0; i < nV; ++i) {
            for (int j = 0; j < np; ++j) {
                Jtr[static_cast<size_t>(j)] += J[static_cast<size_t>(i * np + j)] * r[static_cast<size_t>(i)];
                for (int k = 0; k < np; ++k)
                    JtJ[static_cast<size_t>(j * np + k)] +=
                        J[static_cast<size_t>(i * np + j)] * J[static_cast<size_t>(i * np + k)];
            }
        }
        for (int j = 0; j < np; ++j) {
            Jtr[static_cast<size_t>(j)] = -Jtr[static_cast<size_t>(j)];
            JtJ[static_cast<size_t>(j * np + j)] *= (1.0 + lam);
        }
        std::vector<double> dp(static_cast<size_t>(np), 0.0);
        if (!solveN(np, JtJ.data(), Jtr.data(), dp.data())) {
            lam = std::min(lam * 10.0, 1e12);
            if (lam >= 1e12) break;
            continue;
        }
        double pTry[6];
        double maxStep = 0;
        for (int j = 0; j < np; ++j) {
            pTry[j] = p[j] + dp[static_cast<size_t>(j)];
            const double tol = (j == 3 || j == 4 || j == 5) ? dirTol : q;
            maxStep = std::max(maxStep, std::fabs(dp[static_cast<size_t>(j)]) / std::max(tol, q));
        }
        SurfParams T;
        apply(pTry, T);
        T.alpha = std::fabs(T.alpha);
        const double rss = rssSigned(m, verts, T);
        if (rss <= bestRss) {
            for (int j = 0; j < np; ++j) p[j] = pTry[j];
            bestRss = rss;
            best = T;
            S = T;
            rechart();
            lam = std::max(lam * 0.1, 1e-12);
            if (maxStep <= 1.0) break;
        } else {
            lam = std::min(lam * 10.0, 1e12);
            if (lam >= 1e12) break;
        }
    }
    S = best;
    S.alpha = std::fabs(S.alpha);
    S.n = normalized(S.n);
    return S.alpha > 0.0 && S.alpha < M_PI * 0.5;
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

bool seedConeFromAxis(const Mesh& m, const std::vector<int>& verts, Vec3 a, SurfParams& T) {
    a = normalized(a);
    if (!(norm2(a) > 0.0) || verts.empty()) return false;
    Vec3 c{};
    for (int vi : verts) c = c + m.verts[static_cast<size_t>(vi)];
    c = c * (1.0 / static_cast<double>(verts.size()));
    // SPEC §5.2: apex from least squares on ρ = (apex-offset)·tanα, i.e. a
    // line in the (axial, radial) plane of the vertices — not facet normals.
    double AtA00 = 0, AtA01 = 0, AtA11 = 0, Atb0 = 0, Atb1 = 0;
    int n = 0;
    for (int vi : verts) {
        const Vec3 w = m.verts[static_cast<size_t>(vi)] - c;
        const double ax = dot(w, a);
        const double rho = norm(w - a * ax);
        AtA00 += 1.0;
        AtA01 += ax;
        AtA11 += ax * ax;
        Atb0 += rho;
        Atb1 += rho * ax;
        ++n;
    }
    if (n < 3) return false;
    const double det = AtA00 * AtA11 - AtA01 * AtA01;
    if (!(std::fabs(det) > 0.0)) return false;
    const double p0 = (AtA11 * Atb0 - AtA01 * Atb1) / det;
    const double p1 = (AtA00 * Atb1 - AtA01 * Atb0) / det;
    if (!(std::fabs(p1) > 0.0)) return false;
    T.cls = SurfClass::Cone;
    T.n = a;
    T.alpha = std::atan(std::fabs(p1));
    const double tApex = -p0 / p1;
    T.apex = c + a * tApex;
    return T.alpha > 0.0 && T.alpha < M_PI * 0.5;
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
    auto maxAbs = [&](const SurfParams& T) {
        double mx = 0;
        for (int vi : verts) mx = std::max(mx, distToSurf(m.verts[static_cast<size_t>(vi)], T));
        return mx;
    };
    const Vec3 cands[5] = {axisOf(Cov), axisOf(NNT), Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
    bool any = false;
    double best = 1e300;
    for (const Vec3& ax0 : cands) {
        for (int sgn = 0; sgn < 2; ++sgn) {
            Vec3 ax = (sgn == 0) ? ax0 : (ax0 * -1.0);
            if (dot(ax, nmean) < 0.0) ax = ax * -1.0;
            SurfParams T;
            if (!seedConeFromAxis(m, verts, ax, T)) continue;
            const double mx = maxAbs(T);
            if (!any || mx < best) {
                any = true;
                best = mx;
                S = T;
            }
        }
    }
    if (!any) return false;
    if (refine) {
        const SurfParams seed = S;
        const double r0 = maxAbs(S);
        gaussNewtonCone(m, verts, S);
        S.alpha = std::fabs(S.alpha);
        if (!(S.alpha > 0.0 && S.alpha < M_PI * 0.5) || maxAbs(S) > r0) S = seed;
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

// Centroids wrap an axis (max azimuth gap ≤ π). A torus latitude ring
// around the boss wraps; an S02 octant in a cube corner does not.
bool wrapsAxis(const Mesh& m, const std::vector<int>& triIds, const Vec3& origin, Vec3 axis) {
    axis = normalized(axis);
    if (!(norm2(axis) > 0.0)) return false;
    Vec3 uu, vv;
    frameFromAxis(axis, uu, vv);
    const int bins = paramCount(SurfClass::Cylinder) + 1;
    unsigned seen = 0;
    int nAng = 0;
    for (int t : triIds) {
        const Vec3 w = m.tris[static_cast<size_t>(t)].centroid - origin;
        const Vec3 p = w - axis * dot(w, axis);
        if (!(norm2(p) > 0.0)) continue;
        const Vec3 pn = normalized(p);
        const double ang = std::atan2(dot(pn, vv), dot(pn, uu));
        int b = static_cast<int>(std::floor((ang + M_PI) / (2.0 * M_PI) * bins));
        if (b < 0) b = 0;
        if (b >= bins) b = bins - 1;
        seen |= 1u << b;
        ++nAng;
    }
    if (nAng < paramCount(SurfClass::Plane)) return false;  // 3 azimuths minimum
    return seen == (1u << bins) - 1u;
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

// SPEC §5.3 analog for spheres: the region's normals must span ℝ³. A fillet
// band or plane cluster is rank ≤ 2. Noise floor is Σ area·sin²(θ_q) — no
// size constant. Emit-only (a 2-tri seed is rank 2).
// Radial Gauss-map (port 49345f6): n(t) must track (centroid − centre)
// within the patch's own Gauss radius + θ_q. Handle-lock leftovers were
// rank-3 with maxNormalDev ≈ 0.75 rad — larger than their own span.
bool sphereNormalsSpan(const Mesh& m, const std::vector<int>& triIds, const SurfParams& S) {
    double G[3][3] = {};
    double noise = 0;
    double wsum = 0;
    for (int t : triIds) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        const double w = tr.area;
        wsum += w;
        const Vec3& n = tr.n;
        G[0][0] += w * n.x * n.x;
        G[0][1] += w * n.x * n.y;
        G[0][2] += w * n.x * n.z;
        G[1][1] += w * n.y * n.y;
        G[1][2] += w * n.y * n.z;
        G[2][2] += w * n.z * n.z;
        const double s = std::sin(tr.thetaQ);
        noise += w * s * s;
    }
    G[1][0] = G[0][1];
    G[2][0] = G[0][2];
    G[2][1] = G[1][2];
    if (!(wsum > 0.0)) return false;
    double eval[3];
    double evec[3][3];
    eigen3(G, eval, evec);
    double lmin = eval[0];
    if (eval[1] < lmin) lmin = eval[1];
    if (eval[2] < lmin) lmin = eval[2];
    if (!(lmin > noise)) return false;
    if (!(S.R > 0.0)) return false;
    Vec3 nmean{};
    double thMax = 0;
    for (int t : triIds) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        nmean = nmean + tr.n;
        thMax = std::max(thMax, tr.thetaQ);
    }
    nmean = normalized(nmean);
    if (!(norm2(nmean) > 0.0)) return false;
    double gaussRad = 0;
    for (int t : triIds)
        gaussRad = std::max(gaussRad, angleUnit(m.tris[static_cast<size_t>(t)].n, nmean));
    for (int t : triIds) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        const Vec3 rad = normalized(tr.centroid - S.p0);
        if (!(norm2(rad) > 0.0)) return false;
        if (angleUnit(tr.n, rad) > gaussRad + thMax) return false;
    }
    return true;
}

// The facet-normal clause for the ruled/curved classes (cone, torus).
//
// SPEC §5.3 names the normal clause "THE discriminator" and applies it to the
// plane with theta_q(t), because a facet inscribed in a plane carries the
// plane's normal exactly. On a curved surface it does not: the facet normal is
// the normal at some interior point, so the exact supremum of its deviation
// from the normal at the centroid is the FITTED SURFACE'S OWN normal spread
// over that facet's vertices — a closed form in (S, the mesh's own vertices),
// with theta_q(t) added for quantization. Nothing is chosen: a coarser
// tessellation widens the bound by exactly the amount its own facet widens the
// Gauss image, and a facet that does not lie on the surface fails it.
//
// Without this clause cone/torus growth is unbounded: `admits` was only the
// vertex test, and an under-determined 2- or 3-triangle cone can be re-solved
// to swallow any neighbour, so every seed runs off its feature and the
// over-grown set then fails `certifies`. Measured on the linkage plate: 0 of
// 511 seeds in the cross-bore component grew a certifying cone, although 200
// of those triangles are a 45.000000-degree frustum at 0.107 q (D-140-6 M2).
double normalSpread(const Mesh& m, int t, const SurfParams& S) {
    const Tri& tr = m.tris[static_cast<size_t>(t)];
    const Vec3 nc = normalAt(S, tr.centroid);
    if (!(norm2(nc) > 0.0)) return 0.0;
    double sp = 0;
    for (int k = 0; k < 3; ++k) {
        const Vec3 nv = normalAt(S, m.verts[static_cast<size_t>(tr.v[k])]);
        if (!(norm2(nv) > 0.0)) continue;
        sp = std::max(sp, angleUnit(nc, nv));
    }
    return sp;
}

bool curvedNormalsOk(const Mesh& m, const std::vector<int>& region, const SurfParams& S) {
    // Orientation of the fitted surface against the mesh is a property of the
    // whole region (a bore's facet normals point at the axis), taken once from
    // the region's own agreement, never per triangle.
    double sgn = 0;
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        sgn += dot(tr.n, normalAt(S, tr.centroid));
    }
    const double s = (sgn < 0.0) ? -1.0 : 1.0;
    for (int t : region) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        const Vec3 nc = normalAt(S, tr.centroid) * s;
        if (!(norm2(nc) > 0.0)) return false;
        if (angleUnit(tr.n, nc) > normalSpread(m, t, S) + tr.thetaQ) return false;
    }
    return true;
}

// D-140-6 §1(2): a surface of revolution about `a` through `p0` is a plane
// curve in (rho, z); a torus is a CIRCLE there and dist_torus is identically
// the profile-circle distance. The whole torus certificate is therefore a 2-D
// statement, and over-determination must be counted on the profile.
struct ProfileCensus {
    int levels = 0;      // distinct z, collapsed within tau
    int columns = 0;     // distinct azimuth columns, collapsed within tau of arc
    double lineRes = 0;  // max perpendicular residual of the best profile line
};

ProfileCensus profileCensus(const Mesh& m, const std::vector<int>& verts, const SurfParams& S) {
    ProfileCensus pc;
    if (verts.empty()) return pc;
    const Vec3 a = normalized(S.n);
    if (!(norm2(a) > 0.0)) return pc;
    Vec3 e1{1, 0, 0};
    if (std::fabs(a.x) > std::fabs(a.y)) e1 = Vec3{0, 1, 0};
    e1 = normalized(e1 - a * dot(e1, a));
    const Vec3 e2 = cross(a, e1);
    std::vector<double> zs, phis, rhos;
    zs.reserve(verts.size());
    phis.reserve(verts.size());
    rhos.reserve(verts.size());
    for (int vi : verts) {
        const Vec3 w = m.verts[static_cast<size_t>(vi)] - S.p0;
        const double z = dot(w, a);
        const Vec3 rad = w - a * z;
        zs.push_back(z);
        rhos.push_back(norm(rad));
        phis.push_back(std::atan2(dot(rad, e2), dot(rad, e1)));
    }
    // levels: distinct z within tau
    std::vector<double> zz = zs;
    std::sort(zz.begin(), zz.end());
    pc.levels = 1;
    double last = zz.front();
    for (double v : zz) {
        if (v - last > m.tau) {
            ++pc.levels;
            last = v;
        }
    }
    // columns: distinct azimuth, separated by more than the arc that tau
    // subtends at the region's own largest radius (so the split is the mesh's,
    // not a chosen angle).
    double rhoMax = 0;
    for (double r : rhos) rhoMax = std::max(rhoMax, r);
    const double dphi = (rhoMax > 0.0) ? (m.tau / rhoMax) : 0.0;
    std::vector<double> pp = phis;
    std::sort(pp.begin(), pp.end());
    pc.columns = 1;
    double lastp = pp.front();
    for (double v : pp) {
        if (v - lastp > dphi) {
            ++pc.columns;
            lastp = v;
        }
    }
    // best line through (rho, z) by total least squares, and its max
    // perpendicular residual. D-140-6 §1(4): a straight profile is a cone,
    // a cylinder or a plane, never a torus.
    double mr = 0, mz = 0;
    const double n = static_cast<double>(verts.size());
    for (size_t i = 0; i < rhos.size(); ++i) {
        mr += rhos[i];
        mz += zs[i];
    }
    mr /= n;
    mz /= n;
    double srr = 0, srz = 0, szz = 0;
    for (size_t i = 0; i < rhos.size(); ++i) {
        const double dr = rhos[i] - mr, dz = zs[i] - mz;
        srr += dr * dr;
        srz += dr * dz;
        szz += dz * dz;
    }
    const double tr2 = srr + szz;
    const double det = srr * szz - srz * srz;
    const double disc = std::sqrt(std::max(0.0, tr2 * tr2 * 0.25 - det));
    const double lmin = tr2 * 0.5 - disc;
    double nx = srz, ny = lmin - srr;
    if (!(nx * nx + ny * ny > 0.0)) {
        nx = lmin - szz;
        ny = srz;
    }
    if (!(nx * nx + ny * ny > 0.0)) {
        nx = 1.0;
        ny = 0.0;
    }
    const double nl = std::sqrt(nx * nx + ny * ny);
    nx /= nl;
    ny /= nl;
    for (size_t i = 0; i < rhos.size(); ++i)
        pc.lineRes = std::max(pc.lineRes, std::fabs((rhos[i] - mr) * nx + (zs[i] - mz) * ny));
    return pc;
}

bool admits(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S) {
    // Growth admission: vertex-on-surface (and the plane normal clause). Size
    // gates are applied only when a region is emitted (certifies).
    thread_local Scratch sc;
    sc.ensure(m.verts.size());
    std::vector<int>& verts = sc.verts;
    uniqueVerts(m, region, verts, sc.st);
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
    } else if (c == SurfClass::Cone || c == SurfClass::Torus) {
        if (!curvedNormalsOk(m, region, S)) return false;
    }
    return true;
}

bool certifies(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S,
               double* maxResidOut, double* maxNDevOut) {
    if (static_cast<int>(region.size()) < 2) return false;
    // Own scratch: admits() below reuses its own, so `verts` here stays valid.
    thread_local Scratch sc;
    sc.ensure(m.verts.size());
    std::vector<int>& verts = sc.verts;
    uniqueVerts(m, region, verts, sc.st);
    if (static_cast<int>(verts.size()) < paramCount(c) + 1) return false;
    if (!admits(m, region, c, S)) return false;
    if (c == SurfClass::Sphere && !sphereNormalsSpan(m, region, S)) return false;
    if ((c == SurfClass::Cylinder || c == SurfClass::Sphere) && S.R > m.meshDiag)
        return false;
    if (c == SurfClass::Torus) {
        if (S.R > m.meshDiag || S.r > m.meshDiag) return false;
        // D-140-6 §1, clauses (3) and (4) — the certificate, not a bound.
        const ProfileCensus pc = profileCensus(m, verts, S);
        // (3) PROFILE over-determination. Two levels make the profile a chord:
        // every circle through those two rings fits at zero residual and the
        // fit is a one-parameter family (measured: the plate's 200-triangle
        // band admits (10,2), (12,2) AND a 45-degree cone, all within 0.16 q).
        // The floor is on profile rows, never on |V| in 3-space.
        // A cylinder is an axis line plus a radius, so an axis line carries
        // paramCount(Cylinder) - 1 degrees of freedom; a torus is that same
        // line plus the profile circle, so d_profile = paramCount(Torus) -
        // (paramCount(Cylinder) - 1) = 3, and the floor is d_profile levels.
        // Two columns is what makes the set a surface of revolution at all.
        const int axisDof = paramCount(SurfClass::Cylinder) - 1;
        const int dProfile = paramCount(SurfClass::Torus) - axisDof;
        if (pc.levels < dProfile || pc.columns < axisDof - paramCount(SurfClass::Plane) + 1)
            return false;
        // (4) fewest parameters first, both degeneracies, refused by
        // measurement: a straight profile is a cone / cylinder / plane, and a
        // vanishing major radius is a sphere (d = 4 < 7).
        if (pc.lineRes <= m.tau) return false;
        if (!(S.R > m.tau) || !(S.r > m.tau)) return false;
        if (!(S.R > S.r)) return false;
    }
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
    (void)reverse;  // neighbour order is always ascending (SPEC §5.4)
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

// D-140-6 §1(1)+(2): THE AXIS IS INHERITED, NOT FITTED FREE. Given an axis
// line taken from an already-certified neighbour, every vertex reduces to its
// profile point (rho, z) and the class becomes a 2-D fit there: a cone is a
// LINE in the profile, a torus is a CIRCLE. Two or three triangles cannot
// determine a free 6- or 7-parameter surface, which is why a free seed either
// runs away or never admits; about an inherited axis the same two triangles
// determine the profile exactly.
bool fitAboutAxis(const Mesh& m, const std::vector<int>& verts, SurfClass c, const Vec3& aDirIn,
                  const Vec3& aLoc, SurfParams& S) {
    if (verts.size() < 2) return false;
    Vec3 a = normalized(aDirIn);
    if (!(norm2(a) > 0.0)) return false;
    std::vector<double> rho, zz;
    rho.reserve(verts.size());
    zz.reserve(verts.size());
    for (int vi : verts) {
        const Vec3 w = m.verts[static_cast<size_t>(vi)] - aLoc;
        const double z = dot(w, a);
        rho.push_back(norm(w - a * z));
        zz.push_back(z);
    }
    const double n = static_cast<double>(verts.size());
    double mr = 0, mz = 0;
    for (size_t i = 0; i < rho.size(); ++i) {
        mr += rho[i];
        mz += zz[i];
    }
    mr /= n;
    mz /= n;
    S = SurfParams{};
    S.cls = c;
    if (c == SurfClass::Cone) {
        double srr = 0, srz = 0, szz = 0;
        for (size_t i = 0; i < rho.size(); ++i) {
            const double dr = rho[i] - mr, dz = zz[i] - mz;
            srr += dr * dr;
            srz += dr * dz;
            szz += dz * dz;
        }
        // Total least squares line in (rho, z): the profile of a cone.
        const double tr2 = srr + szz;
        const double disc = std::sqrt(std::max(0.0, tr2 * tr2 * 0.25 - (srr * szz - srz * srz)));
        const double lmax = tr2 * 0.5 + disc;
        double dx = srz, dy = lmax - srr;  // direction along the line, (d rho, d z)
        if (!(dx * dx + dy * dy > 0.0)) {
            dx = lmax - szz;
            dy = srz;
        }
        if (!(dx * dx + dy * dy > 0.0)) return false;
        if (!(std::fabs(dy) > 0.0)) return false;  // rho independent of z: a plane
        const double slope = dx / dy;              // d rho / d z = tan(alpha)
        if (!(std::fabs(slope) > 0.0)) return false;  // constant rho: a cylinder
        const double zApex = mz - mr / slope;
        // The axis must point from the apex into the region, so ax > 0 there.
        if (mz - zApex < 0.0) {
            a = -a;
            S.apex = aLoc + (-a) * zApex;
        } else {
            S.apex = aLoc + a * zApex;
        }
        S.n = a;
        S.alpha = std::atan(std::fabs(slope));
        return S.alpha > 0.0 && S.alpha < M_PI * 0.5;
    }
    if (c == SurfClass::Torus) {
        std::vector<Vec3> prof;
        prof.reserve(verts.size());
        for (size_t i = 0; i < rho.size(); ++i) prof.push_back(Vec3{rho[i], zz[i], 0});
        Vec3 c2{};
        double rmin = 0;
        if (!kasaCircle(prof, Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, c2, rmin)) return false;
        if (!(rmin > 0.0) || !(c2.x > 0.0)) return false;
        S.n = a;
        S.p0 = aLoc + a * c2.y;
        S.R = c2.x;
        S.r = rmin;
        return true;
    }
    return false;
}

// Measurement only (gdiag): where a seed's growth stopped and which clause
// refused it. Never read by the algorithm.
struct GrowTrace {
    int pairs = 0, pairFit = 0, pairAdmit = 0;
    size_t bestGrown = 0;
    int failFit = 0, failCert = 0, failPlane = 0, failSize = 0;
};

// The already-certified neighbours a seed may inherit an axis from (D-140-6
// §1(1)): the axis line of any certified cylinder that owns a triangle sharing
// an edge with the seed. `prior`/`owner` are the oracle state at the phase's
// start; nothing else about them is read.
struct AxisLine {
    Vec3 dir{}, loc{};
};

void inheritedAxes(const Mesh& m, int seed, const std::vector<Oracle>* prior,
                   const std::vector<int>* owner, std::vector<AxisLine>& out) {
    out.clear();
    if (!prior || !owner) return;
    for (int nb : m.adj[static_cast<size_t>(seed)]) {
        const int oi = (*owner)[static_cast<size_t>(nb)];
        if (oi < 0 || oi >= static_cast<int>(prior->size())) continue;
        const Oracle& o = (*prior)[static_cast<size_t>(oi)];
        if (o.cls != SurfClass::Cylinder) continue;
        const AxisLine L{normalized(o.S.n), o.S.p0};
        // D-140-6 §2: the blend class 1.4 recognises is plane-perpendicular-
        // to-cylinder-axis x cylinder. The seed must therefore ALSO neighbour
        // a certified plane whose normal is parallel to that axis, within that
        // plane's own theta_q — the mesh's angle, never a chosen one. An
        // oblique plane (blend.oblique-plane-cyl) and a cylinder-cylinder
        // junction (blend.canal-cylcyl) are both deferred there. The test is
        // on the GROWN REGION, not on the seed: a mouth round touches the
        // plane along one rim and the cylinder along the other, so no single
        // triangle of it neighbours both. See `rimOk` at the emit site.
        bool dup = false;
        for (const AxisLine& e : out) {
            if (std::fabs(std::fabs(dot(e.dir, L.dir)) - 1.0) > 0.0) continue;
            if (norm(cross(L.loc - e.loc, e.dir)) > 0.0) continue;
            dup = true;
        }
        if (!dup) out.push_back(L);
    }
}

bool tryGrow(const Mesh& m, std::vector<char>& claimed, SurfClass c, int seed, bool reverse,
             Oracle& out, std::vector<char>& skipComp, GrowTrace* tr = nullptr,
             const std::vector<Oracle>* prior = nullptr,
             const std::vector<int>* owner = nullptr) {
    if (claimed[static_cast<size_t>(seed)]) return false;
    // For cone/torus `skipComp` records that the WHOLE-COMPONENT fallback has
    // already been tried and refused for this component; it must not retire
    // the component's per-seed growth as well. Measured on the linkage plate:
    // it did, and one failing seed cost the other 510 triangles of the
    // cross-bore component their entire cone phase (pairs=1 for 511 seeds).
    if (skipComp[static_cast<size_t>(seed)] && c != SurfClass::Cone && c != SurfClass::Torus)
        return false;

    auto finish = [&](std::vector<int>& R, SurfParams& S, bool regrow = true) -> bool {
        if (regrow) {
            while (growOnce(m, claimed, c, R, S, reverse)) {
            }
        }
        if (tr) tr->bestGrown = std::max(tr->bestGrown, R.size());
        if (regrow && !fitClassEx(m, R, c, S, true)) {
            if (tr) ++tr->failFit;
            return false;
        }
        if (!certifies(m, R, c, S)) {
            if (tr) ++tr->failCert;
            return false;
        }
        // D-140-1 §3 / D-140-6 §1(4): fewest parameters first, refused by
        // measurement and not by a bound. A region a cheaper class certifies
        // IS that class — a near-zero half-angle cone is a cylinder with an
        // apex 4.9e8 mm away, and a torus whose profile circle degenerates is
        // a cone, a cylinder or a sphere. Each cheaper class is tried in
        // paramCount order; nothing here is a threshold.
        if (c != SurfClass::Plane) {
            SurfParams P;
            if (fitPlane(m, R, P) && certifies(m, R, SurfClass::Plane, P)) {
                if (tr) ++tr->failPlane;
                return false;
            }
        }
        if (c == SurfClass::Cone) {
            // Addendum C's instrument, applied to the cone's own degeneracy:
            // a half-angle that vanishes puts the apex outside the mesh, and a
            // surface whose defining point the mesh cannot reach is a cylinder
            // in this mesh, not a cone (d = 5 < 6). Measured: an apex at
            // x = 4.9e8 mm on a 152 mm part. meshDiag is the same instrument
            // the cylinder/sphere/torus radii are already bounded by.
            Vec3 mn, mx;
            bboxOf(m, R, mn, mx);
            const Vec3 mid = (mn + mx) * 0.5;
            if (norm(S.apex - mid) > m.meshDiag) return false;
        }
        if (c == SurfClass::Torus) {
            for (SurfClass cheaper : {SurfClass::Cylinder, SurfClass::Cone}) {
                SurfParams Q;
                if (fitClassEx(m, R, cheaper, Q, true) && certifies(m, R, cheaper, Q)) {
                    if (tr) ++tr->failPlane;
                    return false;
                }
            }
        }
        if (c == SurfClass::Sphere) {
            const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            bool wrap = false;
            for (const Vec3& ax : axes) {
                if (wrapsAxis(m, R, S.p0, ax)) {
                    wrap = true;
                    break;
                }
            }
            if (wrap) return false;
        }
        // N>=6 is the cylinder/prism law (AGENTS.md): a closed prism with
        // fewer than d_C+1 triangles cannot be a cylinder. Cone is gated
        // only by certifies (|V| >= d_C+1, |R| >= 2).
        if (c == SurfClass::Cylinder) {
            if (static_cast<int>(R.size()) < paramCount(SurfClass::Cylinder) + 1)
                return false;
        }
        // SPEC §5.4 tries sphere before torus. A true sphere (S02) is
        // already claimed; refuse a leftover patch that still certifies
        // as a sphere so it is not emitted as a degenerate torus.
        if (c == SurfClass::Torus) {
            if (!(S.R > S.r) || !(S.r > m.tau)) return false;
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
    (void)reverse;
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
        if (tr) ++tr->pairs;
        if (!fitClass(m, R, c, S)) continue;
        if (tr) ++tr->pairFit;
        if (!admits(m, R, c, S)) continue;
        if (tr) ++tr->pairAdmit;
        if (finish(R, S)) return true;
        if (c == SurfClass::Sphere &&
            static_cast<int>(R.size()) >= paramCount(SurfClass::Sphere) + 1) {
            for (int t : R) skipComp[static_cast<size_t>(t)] = 1;
            break;
        }
    }
    // D-140-6 §1(1): grow about an axis inherited from an already-certified
    // neighbour. The profile fit is exact from the seed pair onward, so a
    // feature that is a strict subset of its connected component (a mouth
    // round inside a cross-bore component) is reachable without ever fitting
    // 6 or 7 free parameters. The free re-solve is then admitted only as a
    // refinement, and only if it still certifies.
    if (c == SurfClass::Cone || c == SurfClass::Torus) {
        std::vector<AxisLine> axes;
        inheritedAxes(m, seed, prior, owner, axes);
        thread_local Scratch axSc;
        axSc.ensure(m.verts.size());
        for (const AxisLine& L : axes) {
            std::vector<int> R{seed};
            std::vector<int> vs;
            SurfParams S;
            bool grew = true;
            while (grew) {
                grew = false;
                std::vector<int> cand;
                std::unordered_set<int> inR(R.begin(), R.end());
                for (int t : R)
                    for (int nb : m.adj[static_cast<size_t>(t)]) {
                        if (claimed[static_cast<size_t>(nb)] || inR.count(nb)) continue;
                        cand.push_back(nb);
                    }
                std::sort(cand.begin(), cand.end());
                cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
                for (int u : cand) {
                    std::vector<int> R2 = R;
                    R2.push_back(u);
                    uniqueVerts(m, R2, vs, axSc.st);
                    SurfParams S2;
                    if (!fitAboutAxis(m, vs, c, L.dir, L.loc, S2)) continue;
                    if (!admits(m, R2, c, S2)) continue;
                    R.swap(R2);
                    S = S2;
                    grew = true;
                }
            }
            if (static_cast<int>(R.size()) < 2) continue;
            uniqueVerts(m, R, vs, axSc.st);
            if (!fitAboutAxis(m, vs, c, L.dir, L.loc, S)) continue;
            // D-140-6 §2 scope, evaluated over the region: some triangle of R
            // must neighbour a certified plane whose normal is parallel to the
            // inherited axis within that plane's own theta_q.
            bool rimOk = false;
            {
                std::unordered_set<int> inR2(R.begin(), R.end());
                for (int t : R) {
                    for (int nb2 : m.adj[static_cast<size_t>(t)]) {
                        if (inR2.count(nb2)) continue;
                        const int oj = (*owner)[static_cast<size_t>(nb2)];
                        if (oj < 0 || oj >= static_cast<int>(prior->size())) continue;
                        const Oracle& pl = (*prior)[static_cast<size_t>(oj)];
                        if (pl.cls != SurfClass::Plane) continue;
                        double th = 0;
                        for (int pt : pl.tris)
                            th = std::max(th, m.tris[static_cast<size_t>(pt)].thetaQ);
                        const double ang = angleUnit(normalized(pl.S.n), L.dir);
                        if (std::min(ang, M_PI - ang) <= th) rimOk = true;
                    }
                }
            }
            if (!rimOk) continue;
            std::vector<int> Rk = R;
            SurfParams Sk = S;
            if (finish(Rk, Sk, false)) return true;
        }
    }
    // Cone / torus: a local patch is under-determined and looks like a
    // cylinder/sphere. Take the whole unclaimed connected component (the
    // S04 512-tri band, the S03 drafted hole), then fit. finish() still
    // requires certifies — a freeform island fails and stays unclaimed.
    if (c == SurfClass::Cone || c == SurfClass::Torus) {
        if (skipComp[static_cast<size_t>(seed)]) return false;
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

void buildOracle(const Mesh& m, OracleSet& out, bool reverseSeeds, int seedOrder) {
    out = OracleSet{};
    out.owner.assign(m.tris.size(), -1);
    std::vector<char> claimed(m.tris.size(), 0);
    // SPEC §5.4: fewest parameters first {plane, cylinder, cone, sphere, torus}.
    // S04 boss-top stays plane (d=3 before d=7). If the torus band then
    // falls to sphere, that is sphereNormalsSpan's defect — never reorder.
    const SurfClass order[5] = {SurfClass::Plane, SurfClass::Cylinder, SurfClass::Cone,
                                SurfClass::Sphere, SurfClass::Torus};

    std::vector<Oracle> heldPlanes;  // |R|==2 tessellation quads; curved classes may claim them

    auto seedList = [&]() {
        std::vector<int> s;
        s.reserve(m.tris.size());
        if (!reverseSeeds) {
            for (int i = 0; i < static_cast<int>(m.tris.size()); ++i) s.push_back(i);
        } else {
            for (int i = static_cast<int>(m.tris.size()) - 1; i >= 0; --i) s.push_back(i);
        }
        if (seedOrder != 0) {
            // §8 case 9 third permutation, test-only: enter the mesh in
            // triangle-centroid order instead of file order — a permutation
            // taken from the geometry itself, so it introduces no constant.
            // seedOrder > 0 ascending, < 0 descending; ties by index, so the
            // order is total and platform-independent.
            const bool desc = seedOrder < 0;
            std::stable_sort(s.begin(), s.end(), [&](int a, int b) {
                const Vec3& ca = m.tris[static_cast<size_t>(a)].centroid;
                const Vec3& cb = m.tris[static_cast<size_t>(b)].centroid;
                if (ca.x != cb.x) return desc ? (ca.x > cb.x) : (ca.x < cb.x);
                if (ca.y != cb.y) return desc ? (ca.y > cb.y) : (ca.y < cb.y);
                if (ca.z != cb.z) return desc ? (ca.z > cb.z) : (ca.z < cb.z);
                return a < b;
            });
        }
        return s;
    };

    auto restabilize = [&]() {
        for (Oracle& o : out.oracles) {
            o.w = regionArea(m, o.tris);
            o.minVertIndex = minVertOf(m, o.tris);
            std::sort(o.tris.begin(), o.tris.end());
        }
        std::stable_sort(out.oracles.begin(), out.oracles.end(),
                         [](const Oracle& A, const Oracle& B) {
                             if (A.cls != B.cls)
                                 return static_cast<int>(A.cls) < static_cast<int>(B.cls);
                             if (A.w != B.w) return A.w > B.w;
                             return A.minVertIndex < B.minVertIndex;
                         });
        out.owner.assign(m.tris.size(), -1);
        for (int i = 0; i < static_cast<int>(out.oracles.size()); ++i) {
            out.oracles[static_cast<size_t>(i)].id = i;
            for (int t : out.oracles[static_cast<size_t>(i)].tris)
                out.owner[static_cast<size_t>(t)] = i;
        }
    };

    // One pass of SPEC §5.4's growth loop over `seeds`: committed regions go to
    // `sink`, deferred 2-triangle planes to `held`, and `claimedV` is marked.
    auto walk = [&](SurfClass c, const std::vector<int>& seeds, std::vector<char>& claimedV,
                    std::vector<Oracle>& sink, std::vector<Oracle>& held) {
        const std::vector<Oracle> prior = out.oracles;
        const std::vector<int> priorOwner = out.owner;
        std::vector<char> skipComp(m.tris.size(), 0);
        for (int seed : seeds) {
            if (claimedV[static_cast<size_t>(seed)]) continue;
            Oracle o;
            if (!tryGrow(m, claimedV, c, seed, false, o, skipComp, nullptr, &prior, &priorOwner))
                continue;
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
                    for (int t : o.tris) claimedV[static_cast<size_t>(t)] = 1;
                    sink.push_back(std::move(o));
                    continue;
                }
                held.push_back(std::move(o));
                continue;
            }
            for (int t : o.tris) claimedV[static_cast<size_t>(t)] = 1;
            sink.push_back(std::move(o));
        }
    };

    // Stage A of each class phase: WHICH COMPONENTS CARRY THIS CLASS.
    // `ambient` is the unclaimed set at the phase's start and is frozen for
    // the whole pass, so no seed consumes another seed's triangles and R(t) —
    // the maximal region grown from t — is a function of (mesh, ambient,
    // class, t). Write CORE(c) = { t : R(t) certifies }; it names no seed
    // order. The pass answers ONE BOOLEAN PER COMPONENT of the unclaimed set,
    // `carries[K]`, and never materialises CORE(c) or a coverage bitmap:
    //   * carries[K] is set only when some seed t in K grew a certifying
    //     region, i.e. only when K meets CORE(c);
    //   * conversely, if K meets CORE(c), take any t in K ∩ CORE(c). seedList()
    //     enumerates every triangle, so t is reached; either carries[K] was
    //     already set by an earlier seed of K, or the early-out below does not
    //     fire and R(t) certifies, setting it. Either way carries[K] = 1.
    // So carries[K] ⟺ K ∩ CORE(c) ≠ ∅, with no reference to the seed order —
    // and THAT is the quantity the partition is taken over (see `mine` below).
    // The early-out (skip a component once it has answered yes) is what makes
    // this affordable: one certifying grow retires a whole component instead
    // of leaving the rest of it to be re-grown. Which triangle inside a
    // component a stage-A grow happens to cover still depends on the walk
    // (growth stops after a fixed run, so a seed entering a tessellated band
    // mid-way chops it differently); which components carry the class does
    // not, and stage B re-derives the partition from the lowest index.
    // `skipComp` is tryGrow's own negative cache, not state this pass reads:
    // it retires triangles a grow already tried and failed on (for cone/torus
    // that is the whole component, so the retirement is itself per-component).
    auto coverage = [&](SurfClass c, const std::vector<int>& seeds,
                        const std::vector<char>& ambient, const std::vector<int>& compId,
                        int nComp) {
        std::vector<char> carries(static_cast<size_t>(nComp) + 1, 0);
        std::vector<char> amb = ambient;  // tryGrow reads it; the pass never writes it
        const std::vector<Oracle> prior = out.oracles;
        const std::vector<int> priorOwner = out.owner;
        std::vector<char> skipComp(m.tris.size(), 0);
        for (int seed : seeds) {
            if (amb[static_cast<size_t>(seed)]) continue;
            const int k = compId[static_cast<size_t>(seed)];
            if (carries[static_cast<size_t>(k)]) continue;  // already answered for this component
            Oracle o;
            if (!tryGrow(m, amb, c, seed, false, o, skipComp, nullptr, &prior, &priorOwner))
                continue;
            carries[static_cast<size_t>(k)] = 1;
        }
        return carries;
    };

    // Edge-connected components of the unclaimed set. Derived from the mesh's
    // own adjacency; no tolerance, no constant.
    auto componentsOf = [&](const std::vector<char>& ambient, std::vector<int>& compId) {
        compId.assign(m.tris.size(), -1);
        int n = 0;
        for (int t = 0; t < static_cast<int>(m.tris.size()); ++t) {
            if (ambient[static_cast<size_t>(t)] || compId[static_cast<size_t>(t)] >= 0) continue;
            compId[static_cast<size_t>(t)] = n;
            std::vector<int> stack{t};
            while (!stack.empty()) {
                const int x = stack.back();
                stack.pop_back();
                for (int nb : m.adj[static_cast<size_t>(x)]) {
                    if (ambient[static_cast<size_t>(nb)] || compId[static_cast<size_t>(nb)] >= 0)
                        continue;
                    compId[static_cast<size_t>(nb)] = n;
                    stack.push_back(nb);
                }
            }
            ++n;
        }
        return n;
    };

    for (SurfClass c : order) {
        std::vector<int> compId;
        const int nComp = componentsOf(claimed, compId);
        if (gdiag()) {
            int unc = 0;
            for (int t = 0; t < static_cast<int>(m.tris.size()); ++t)
                if (!claimed[static_cast<size_t>(t)]) ++unc;
            std::fprintf(stderr, "GRADE_PHASE class=%s unclaimed=%d comps=%d q=%.9g\n",
                         className(c), unc, nComp, m.q);
            if (c == SurfClass::Cone || c == SurfClass::Torus) {
                std::vector<std::vector<int>> comps(static_cast<size_t>(nComp));
                for (int t = 0; t < static_cast<int>(m.tris.size()); ++t)
                    if (!claimed[static_cast<size_t>(t)])
                        comps[static_cast<size_t>(compId[static_cast<size_t>(t)])].push_back(t);
                std::vector<int> ord(static_cast<size_t>(nComp));
                for (int i = 0; i < nComp; ++i) ord[static_cast<size_t>(i)] = i;
                std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
                    return regionArea(m, comps[static_cast<size_t>(a)]) >
                           regionArea(m, comps[static_cast<size_t>(b)]);
                });
                for (int k = 0; k < nComp; ++k) {
                    const std::vector<int>& K = comps[static_cast<size_t>(ord[static_cast<size_t>(k)])];
                    if (K.size() < 2) continue;
                    Vec3 mn, mx;
                    bboxOf(m, K, mn, mx);
                    SurfParams S;
                    const bool fo = fitClass(m, K, c, S);
                    double mr = 0;
                    const bool ce = fo && certifies(m, K, c, S, &mr);
                    // What the growth walk actually finds inside this
                    // component: the largest certifying region over a frozen
                    // ambient (no consumption), and how many seeds certify.
                    {
                        std::vector<char> amb(m.tris.size(), 1);
                        for (int t : K) amb[static_cast<size_t>(t)] = 0;
                        std::vector<char> skipC(m.tris.size(), 0);
                        int nOk = 0;
                        size_t bestN = 0;
                        double bestA = 0;
                        Oracle bo;
                        GrowTrace tr;
                        for (int seed : K) {
                            Oracle o;
                            if (!tryGrow(m, amb, c, seed, false, o, skipC, &tr, &out.oracles,
                                         &out.owner))
                                continue;
                            ++nOk;
                            const double a2 = regionArea(m, o.tris);
                            if (o.tris.size() > bestN) {
                                bestN = o.tris.size();
                                bestA = a2;
                                bo = o;
                            }
                        }
                        std::fprintf(stderr,
                                     "GRADE_GROW class=%s comp=%d seedsOk=%d bestN=%zu bestArea=%.6f "
                                     "R=%.6f r=%.6f alpha=%.6f pairs=%d pairFit=%d pairAdmit=%d "
                                     "bestGrown=%zu failFit=%d failCert=%d failPlane=%d\n",
                                     className(c), ord[static_cast<size_t>(k)], nOk, bestN, bestA,
                                     bo.S.R, bo.S.r, bo.S.alpha, tr.pairs, tr.pairFit, tr.pairAdmit,
                                     tr.bestGrown, tr.failFit, tr.failCert, tr.failPlane);
                    }
                    std::fprintf(stderr,
                                 "GRADE_CAND class=%s comp=%d n=%zu area=%.6f fit=%d cert=%d "
                                 "R=%.6f r=%.6f alpha=%.6f maxResid=%.6eq bbox=[%.4f %.4f %.4f]-[%.4f %.4f %.4f]\n",
                                 className(c), ord[static_cast<size_t>(k)], K.size(),
                                 regionArea(m, K), fo ? 1 : 0, ce ? 1 : 0, S.R, S.r, S.alpha,
                                 m.q > 0 ? mr / m.q : 0.0, mn.x, mn.y, mn.z, mx.x, mx.y, mx.z);
                }
            }
        }
        // Stage A — the growth-order probe runs here, over seedList().
        const std::vector<char> carries = coverage(c, seedList(), claimed, compId, nComp);
        // Stage B — CANONICAL PARTITION. The domain is every unclaimed
        // triangle of every component the class was certified in — a set the
        // seed direction cannot change (proof above). Re-run §5.4's walk over
        // it, entered at the lowest triangle index: same certificate,
        // canonical entry point. Triangles the canonical walk does not take
        // go back to unclaimed, and a component in which nothing certified is
        // never entered at all.
        std::vector<int> mine;
        for (int t = 0; t < static_cast<int>(m.tris.size()); ++t)
            if (!claimed[static_cast<size_t>(t)] &&
                carries[static_cast<size_t>(compId[static_cast<size_t>(t)])])
                mine.push_back(t);
        std::vector<char> mask(m.tris.size(), 1);
        for (int t : mine) mask[static_cast<size_t>(t)] = 0;
        std::vector<Oracle> canon, canonHeld;
        walk(c, mine, mask, canon, canonHeld);
        for (int t : mine) claimed[static_cast<size_t>(t)] = mask[static_cast<size_t>(t)];
        for (Oracle& o : canon) {
            o.id = static_cast<int>(out.oracles.size());
            out.oracles.push_back(std::move(o));
        }
        for (Oracle& o : canonHeld) heldPlanes.push_back(std::move(o));
        if (gdiag())
            std::fprintf(stderr, "GRADE_PHASE_OUT class=%s domain=%zu emitted=%zu held=%zu\n",
                         className(c), mine.size(), canon.size(), canonHeld.size());
        restabilize();
        // Claim pass (islands). Triangle order is always ascending — seed
        // reversal is a growth-order probe; SPEC §5.4's claim/merge must be
        // a canonical maximal partition (area desc, then min welded-vertex
        // index — never input index / region id).
        bool claimChanged = true;
        while (claimChanged) {
            claimChanged = false;
            for (int u = 0; u < static_cast<int>(m.tris.size()); ++u) {
                if (claimed[static_cast<size_t>(u)]) continue;
                int best = -1;
                double bestResid = 1e300;
                double bestArea = -1;
                int bestMinV = 0;
                SurfParams bestS;
                for (int ri = 0; ri < static_cast<int>(out.oracles.size()); ++ri) {
                    Oracle& Rk = out.oracles[static_cast<size_t>(ri)];
                    if (Rk.cls != c) continue;
                    std::vector<int> R2 = Rk.tris;
                    R2.push_back(u);
                    SurfParams S2;
                    if (!fitClassEx(m, R2, c, S2, false)) continue;
                    double maxR = 0;
                    if (!certifies(m, R2, c, S2, &maxR)) continue;
                    const double area = regionArea(m, Rk.tris);
                    const int mv = minVertOf(m, Rk.tris);
                    const bool better =
                        (maxR < bestResid) ||
                        (maxR == bestResid && area > bestArea) ||
                        (maxR == bestResid && area == bestArea && (best < 0 || mv < bestMinV));
                    if (better) {
                        bestResid = maxR;
                        bestArea = area;
                        bestMinV = mv;
                        best = ri;
                        bestS = S2;
                    }
                }
                if (best >= 0) {
                    Oracle& Rk = out.oracles[static_cast<size_t>(best)];
                    Rk.tris.push_back(u);
                    Rk.S = bestS;
                    Rk.w = regionArea(m, Rk.tris);
                    Rk.minVertIndex = minVertOf(m, Rk.tris);
                    claimed[static_cast<size_t>(u)] = 1;
                    out.owner[static_cast<size_t>(u)] = best;
                    claimChanged = true;
                }
            }
        }
        // Merge pass: oracles of this class in canonical order (area desc,
        // min welded-vertex), then first certifying pair. Deterministic
        // without an all-pairs scan each round (F12).
        bool merged = true;
        while (merged) {
            merged = false;
            std::vector<int> idx;
            for (int i = 0; i < static_cast<int>(out.oracles.size()); ++i)
                if (out.oracles[static_cast<size_t>(i)].cls == c) idx.push_back(i);
            std::sort(idx.begin(), idx.end(), [&](int a, int b) {
                const double wa = regionArea(m, out.oracles[static_cast<size_t>(a)].tris);
                const double wb = regionArea(m, out.oracles[static_cast<size_t>(b)].tris);
                if (wa != wb) return wa > wb;
                return minVertOf(m, out.oracles[static_cast<size_t>(a)].tris) <
                       minVertOf(m, out.oracles[static_cast<size_t>(b)].tris);
            });
            for (size_t ii = 0; ii < idx.size() && !merged; ++ii) {
                const int i = idx[ii];
                for (size_t jj = ii + 1; jj < idx.size(); ++jj) {
                    const int j = idx[jj];
                    std::vector<int> U = out.oracles[static_cast<size_t>(i)].tris;
                    U.insert(U.end(), out.oracles[static_cast<size_t>(j)].tris.begin(),
                             out.oracles[static_cast<size_t>(j)].tris.end());
                    SurfParams S;
                    if (!fitClassEx(m, U, c, S, false)) continue;
                    if (!certifies(m, U, c, S)) continue;
                    const int lo = std::min(i, j), hi = std::max(i, j);
                    out.oracles[static_cast<size_t>(lo)].tris.swap(U);
                    out.oracles[static_cast<size_t>(lo)].S = S;
                    for (int t : out.oracles[static_cast<size_t>(hi)].tris)
                        out.owner[static_cast<size_t>(t)] = lo;
                    out.oracles.erase(out.oracles.begin() + hi);
                    for (int t = 0; t < static_cast<int>(m.tris.size()); ++t) {
                        if (out.owner[static_cast<size_t>(t)] > hi)
                            --out.owner[static_cast<size_t>(t)];
                    }
                    for (int k = 0; k < static_cast<int>(out.oracles.size()); ++k)
                        out.oracles[static_cast<size_t>(k)].id = k;
                    merged = true;
                    break;
                }
            }
        }
        restabilize();
    }

    // Commit held 2-triangle planes that curved classes did not take.
    std::sort(heldPlanes.begin(), heldPlanes.end(), [](const Oracle& A, const Oracle& B) {
        const int ma = A.tris.empty() ? 0 : *std::min_element(A.tris.begin(), A.tris.end());
        const int mb = B.tris.empty() ? 0 : *std::min_element(B.tris.begin(), B.tris.end());
        if (ma != mb) return ma < mb;
        return A.w > B.w;
    });
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
        // D-130-12 / D-140-6 §1(5): a re-solve is admitted only if the region
        // still certifies under it; otherwise roll back to the parameters the
        // region was certified with. Without this a curved oracle can report
        // parameters no clause ever passed (measured: a 6-triangle pickup
        // torus reporting Rmin > Rmaj).
        {
            Oracle& oi = out.oracles[static_cast<size_t>(i)];
            const SurfParams keep = oi.S;
            if (!fitClass(m, oi.tris, oi.cls, oi.S) || !certifies(m, oi.tris, oi.cls, oi.S))
                oi.S = keep;
        }
        for (int t : out.oracles[static_cast<size_t>(i)].tris)
            out.owner[static_cast<size_t>(t)] = i;
        fillOracleStats(m, out.oracles[static_cast<size_t>(i)]);
    }

    if (const char* op = std::getenv("STL2STEP_GRADE_OWNERS")) {
        if (op[0]) {
            if (FILE* f = std::fopen(op, "w")) {
                for (int t = 0; t < static_cast<int>(m.tris.size()); ++t) {
                    const int oi = out.owner[static_cast<size_t>(t)];
                    const Vec3& cn = m.tris[static_cast<size_t>(t)].centroid;
                    std::fprintf(f, "%d %d %s %.9g %.9g %.9g\n", t, oi,
                                 oi < 0 ? "residue"
                                        : className(out.oracles[static_cast<size_t>(oi)].cls),
                                 cn.x, cn.y, cn.z);
                }
                std::fclose(f);
            }
        }
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
    // Canonical oracle order: area desc, then featureId (SPEC §5.4 / §7.2).
    // Seed order must not leak into assignment ties.
    std::sort(out.oracles.begin(), out.oracles.end(), [](const Oracle& a, const Oracle& b) {
        if (a.w != b.w) return a.w > b.w;
        if (a.featureId != b.featureId) return a.featureId < b.featureId;
        return a.minVertIndex < b.minVertIndex;
    });
    out.owner.assign(m.tris.size(), -1);
    for (int i = 0; i < static_cast<int>(out.oracles.size()); ++i) {
        out.oracles[static_cast<size_t>(i)].id = i;
        for (int t : out.oracles[static_cast<size_t>(i)].tris)
            out.owner[static_cast<size_t>(t)] = i;
    }
}

}  // namespace grade
