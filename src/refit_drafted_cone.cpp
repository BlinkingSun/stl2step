// Drafted wall about a certified plane normal (D-gid-8.18).
//
// A moulded wall is drafted: it is a cone whose axis is the pull direction,
// and the pull direction is the normal of the certified planes the part is
// pulled along. Taking the direction from such a plane leaves the cone four
// free parameters -- the circle it cuts in the plane (centre 2, radius 1) and
// the draft (1) -- which every vertex of the wall over-determines. A cone with
// no certified plane normal for its axis is a free cone and is never claimed
// here (D-gid-3.6).
//
// Stage: runTorusStages, after stage T, H1, the anchored rounds and the sphere
// corners -- every family a certified round or donor anchors -- and before A3
// (refit_segment.cpp). Band: whole atoms -- an Unclaimed provisional
// that is not itself a certified plane, or an accepted region that does not
// certify its surface (D-fc-4.3 robbableUncertified) -- grown from a seed of
// two edge-adjacent atoms. An atom joins only whole and only while the refit
// band keeps the certificate (certifyDraftedCone).
//
// Every threshold is tau = 2q or a count of parameters (D-fc-1.6).
// SPDX-License-Identifier: MIT

#include "refit_drafted_cone.hpp"

#include "refit_internal.hpp"
#include "refit_profile_cert.hpp"
#include "refit_band_volume.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

namespace stl2step {
namespace refit {

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;

// I5 determinism cap of the least-squares fits (fitProfileCircle's own), not
// an acceptance threshold.
constexpr int kDraftFitIters = 16;

bool diagOn() {
    static const bool on = [] {
        const char* e = std::getenv("STL2STEP_DIAG_DRAFT");
        return e && e[0] && e[0] != '0';
    }();
    return on;
}

void diagPrintf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
}

// Gaussian elimination with partial pivoting on an n x n system (n <= 4).
// False when a pivot vanishes.
template <int N>
bool solveN(std::array<std::array<double, N>, N> a, std::array<double, N> b,
            std::array<double, N>& x) {
    for (int c = 0; c < N; ++c) {
        int piv = c;
        for (int r = c + 1; r < N; ++r)
            if (std::fabs(a[r][c]) > std::fabs(a[piv][c])) piv = r;
        if (!(std::fabs(a[piv][c]) > gp::Resolution())) return false;
        std::swap(a[c], a[piv]);
        std::swap(b[c], b[piv]);
        for (int r = c + 1; r < N; ++r) {
            const double f = a[r][c] / a[c][c];
            for (int k = c; k < N; ++k) a[r][k] -= f * a[c][k];
            b[r] -= f * b[c];
        }
    }
    for (int c = N - 1; c >= 0; --c) {
        double s = b[c];
        for (int k = c + 1; k < N; ++k) s -= a[c][k] * x[k];
        x[c] = s / a[c][c];
        if (!std::isfinite(x[c])) return false;
    }
    return true;
}

// Unit vectors e1, e2 completing the unit axis a to a right-handed frame.
void frameOf(const gp_XYZ& a, gp_XYZ& e1, gp_XYZ& e2) {
    const double ax = std::fabs(a.X()), ay = std::fabs(a.Y()), az = std::fabs(a.Z());
    const gp_XYZ t = (ax <= ay && ax <= az) ? gp_XYZ(1, 0, 0)
                     : (ay <= az)           ? gp_XYZ(0, 1, 0)
                                            : gp_XYZ(0, 0, 1);
    e1 = a.Crossed(t);
    e1 /= e1.Modulus();
    e2 = a.Crossed(e1);
}

// Max distance of the points to their least-squares plane (three parameters).
double planeResidual(const std::vector<gp_XYZ>& pts) {
    if (pts.size() < 3) return 0.0;
    gp_XYZ m(0, 0, 0);
    for (const gp_XYZ& p : pts) m += p;
    m /= (double)pts.size();
    std::array<std::array<double, 3>, 3> cov{};
    for (const gp_XYZ& p : pts) {
        const double d[3] = {p.X() - m.X(), p.Y() - m.Y(), p.Z() - m.Z()};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) cov[i][j] += d[i] * d[j];
    }
    std::array<double, 3> ev{};
    std::array<std::array<double, 3>, 3> evec{};
    if (!jacobiEigenSymmetric3(cov, ev, evec)) return 1e300;
    const gp_XYZ n(evec[0][0], evec[0][1], evec[0][2]);
    const double nm = n.Modulus();
    if (!(nm > gp::Resolution())) return 1e300;
    double mx = 0.0;
    for (const gp_XYZ& p : pts) mx = std::max(mx, std::fabs((p - m).Dot(n) / nm));
    return mx;
}

// A cylinder about the axis direction (three parameters): the Kasa circle of
// the section points, then Gauss-Newton on (cx, cy, R). fitted is false when
// the Kasa system is singular.
struct AxisCylinder {
    bool fitted = false;
    double cx = 0, cy = 0, R = 0, dev = 1e300;
};

AxisCylinder fitAxisCylinder(const std::vector<double>& X, const std::vector<double>& Y,
                             double tau) {
    AxisCylinder c;
    const std::size_t N = std::min(X.size(), Y.size());
    std::array<std::array<double, 3>, 3> A{};
    std::array<double, 3> b{};
    for (std::size_t i = 0; i < N; ++i) {
        const double row[3] = {X[i], Y[i], 1.0};
        const double u = X[i] * X[i] + Y[i] * Y[i];
        for (int r = 0; r < 3; ++r) {
            for (int k = 0; k < 3; ++k) A[r][k] += row[r] * row[k];
            b[r] += row[r] * u;
        }
    }
    std::array<double, 3> s{};
    if (!solveN<3>(A, b, s)) return c;
    double cx = 0.5 * s[0], cy = 0.5 * s[1];
    double R = std::sqrt(std::max(0.0, s[2] + cx * cx + cy * cy));
    for (int it = 0; it < kDraftFitIters; ++it) {
        std::array<std::array<double, 3>, 3> JJ{};
        std::array<double, 3> Jr{};
        for (std::size_t i = 0; i < N; ++i) {
            const double dx = X[i] - cx, dy = Y[i] - cy;
            const double rho = std::hypot(dx, dy);
            if (!(rho > gp::Resolution())) continue;
            const double j[3] = {-dx / rho, -dy / rho, -1.0};
            const double res = rho - R;
            for (int r = 0; r < 3; ++r) {
                for (int k = 0; k < 3; ++k) JJ[r][k] += j[r] * j[k];
                Jr[r] -= j[r] * res;
            }
        }
        std::array<double, 3> d{};
        if (!solveN<3>(JJ, Jr, d)) break;
        cx += d[0];
        cy += d[1];
        R += d[2];
        if (std::hypot(std::hypot(d[0], d[1]), d[2]) <= tau) break;
    }
    double dev = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        dev = std::max(dev, std::fabs(std::hypot(X[i] - cx, Y[i] - cy) - R));
    c.fitted = true;
    c.cx = cx;
    c.cy = cy;
    c.R = R;
    c.dev = std::isfinite(dev) ? dev : 1e300;
    return c;
}

// Max distance of the points (section coordinates X, Y and height H) to their
// least-squares sphere (four parameters): the algebraic sphere, then
// Gauss-Newton on (centre, radius), iteration cap kDraftFitIters. 1e300 when
// the algebraic system is singular.
double sphereResidual(const std::vector<double>& X, const std::vector<double>& Y,
                      const std::vector<double>& H, double tau) {
    const std::size_t N = std::min(std::min(X.size(), Y.size()), H.size());
    std::array<std::array<double, 4>, 4> A{};
    std::array<double, 4> b{};
    for (std::size_t i = 0; i < N; ++i) {
        const double row[4] = {X[i], Y[i], H[i], 1.0};
        const double u = X[i] * X[i] + Y[i] * Y[i] + H[i] * H[i];
        for (int r = 0; r < 4; ++r) {
            for (int k = 0; k < 4; ++k) A[r][k] += row[r] * row[k];
            b[r] += row[r] * u;
        }
    }
    std::array<double, 4> q{};
    if (!solveN<4>(A, b, q)) return 1e300;
    double sx = 0.5 * q[0], sy = 0.5 * q[1], sz = 0.5 * q[2];
    double sr = std::sqrt(std::max(0.0, q[3] + sx * sx + sy * sy + sz * sz));
    for (int it = 0; it < kDraftFitIters; ++it) {
        std::array<std::array<double, 4>, 4> JJ{};
        std::array<double, 4> Jr{};
        for (std::size_t i = 0; i < N; ++i) {
            const double dx = X[i] - sx, dy = Y[i] - sy, dz = H[i] - sz;
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (!(d > gp::Resolution())) continue;
            const double j[4] = {-dx / d, -dy / d, -dz / d, -1.0};
            const double res = d - sr;
            for (int r = 0; r < 4; ++r) {
                for (int k = 0; k < 4; ++k) JJ[r][k] += j[r] * j[k];
                Jr[r] -= j[r] * res;
            }
        }
        std::array<double, 4> dd{};
        if (!solveN<4>(JJ, Jr, dd)) break;
        sx += dd[0];
        sy += dd[1];
        sz += dd[2];
        sr += dd[3];
        if (std::hypot(std::hypot(dd[0], dd[1]), std::hypot(dd[2], dd[3])) <= tau) break;
    }
    double dev = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double dx = X[i] - sx, dy = Y[i] - sy, dz = H[i] - sz;
        dev = std::max(dev, std::fabs(std::sqrt(dx * dx + dy * dy + dz * dz) - sr));
    }
    return std::isfinite(dev) ? dev : 1e300;
}

// D-gid-8.29. Fewest parameters first, at a tie too: a band that certifies at
// tau as another surface family with no more parameters than the drafted cone
// (kParamsDraftCone) is ambiguous -- the mesh does not determine its surface --
// and the cone, the newer anchored hypothesis, does not claim it; that
// family's own stage keeps its certificate. The families: a plane (3), a
// cylinder about the axis direction (3: no resolvable draft) and a sphere (4).
// Returns the family that certifies, or nullptr; f receives the residuals of
// the families tried.
// D-gid-8.30: the fits run in band-centred coordinates (the centroid of the
// band's points), so the refusal does not depend on where the frame's origin
// lies: an algebraic fit far from its origin is ill-conditioned. The plane is
// PCA about the points' mean already.
const char* ambiguousWithCheaperFamily(const std::vector<gp_XYZ>& pts,
                                       const std::vector<double>& X,
                                       const std::vector<double>& Y,
                                       const std::vector<double>& H, double tau,
                                       DraftedConeFit& f) {
    f.planeDev = planeResidual(pts);
    if (f.planeDev <= tau) return "plane";
    const std::size_t N = std::min(std::min(X.size(), Y.size()), H.size());
    double mx = 0.0, my = 0.0, mh = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        mx += X[i];
        my += Y[i];
        mh += H[i];
    }
    if (N > 0) {
        mx /= (double)N;
        my /= (double)N;
        mh /= (double)N;
    }
    std::vector<double> Xc(N), Yc(N), Hc(N);
    for (std::size_t i = 0; i < N; ++i) {
        Xc[i] = X[i] - mx;
        Yc[i] = Y[i] - my;
        Hc[i] = H[i] - mh;
    }
    const AxisCylinder cyl = fitAxisCylinder(Xc, Yc, tau);
    f.cylDev = cyl.dev;
    if (cyl.fitted && cyl.dev <= tau) return "cylinder";
    f.sphereDev = sphereResidual(Xc, Yc, Hc, tau);
    if (f.sphereDev <= tau) return "sphere";
    return nullptr;
}

}  // namespace

// D-gid-8.25 (a). Two heights are one level unless more than levelSep(tau) =
// 2 tau apart: a point is certified anywhere in a slab of width 2 tau about the
// surface, so heights within that diameter are one level to the certificate.
// The levels are the clusters of the sorted heights, a new one starting at a
// gap wider than levelSep (any two heights in different levels are then more
// than levelSep apart). Each entry of the result is one level's point indices,
// ascending height.
std::vector<std::vector<std::size_t>> levelClusters(const std::vector<double>& h, double tau) {
    std::vector<std::size_t> ord(h.size());
    for (std::size_t i = 0; i < h.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
        if (h[a] != h[b]) return h[a] < h[b];
        return a < b;
    });
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t i = 0; i < ord.size(); ++i) {
        if (i == 0 || h[ord[i]] - h[ord[i - 1]] > levelSep(tau)) out.emplace_back();
        out.back().push_back(ord[i]);
    }
    return out;
}

int countLevelClusters(const std::vector<double>& h, double tau) {
    return (int)levelClusters(h, tau).size();
}

int maxLevelStations(const std::vector<double>& h, const std::vector<ProfPt>& sec, double tau) {
    int best = 0;
    for (const std::vector<std::size_t>& lv : levelClusters(h, tau)) {
        std::vector<ProfPt> p;
        p.reserve(lv.size());
        for (std::size_t i : lv)
            if (i < sec.size()) p.push_back(sec[i]);
        best = std::max(best, countProfPts(p, tau));
    }
    return best;
}

DraftedConeFit certifyDraftedCone(const std::vector<gp_XYZ>& pts, const gp_XYZ& n,
                                  const gp_XYZ& p0, double tau) {
    DraftedConeFit f;
    const double nm = n.Modulus();
    if (!(tau > 0.0) || !(nm > gp::Resolution()) || pts.empty()) return f;
    f.axis = n / nm;
    gp_XYZ e1, e2;
    frameOf(f.axis, e1, e2);
    const std::size_t N = pts.size();
    std::vector<double> X(N), Y(N), H(N);
    std::vector<ProfPt> sec(N);
    for (std::size_t i = 0; i < N; ++i) {
        const gp_XYZ d = pts[i] - p0;
        H[i] = d.Dot(f.axis);
        X[i] = d.Dot(e1);
        Y[i] = d.Dot(e2);
        sec[i].rho = X[i];
        sec[i].z = Y[i];
    }
    f.hLo = *std::min_element(H.begin(), H.end());
    f.hHi = *std::max_element(H.begin(), H.end());
    f.stations = countProfPts(sec, tau);
    f.levels = countLevelClusters(H, tau);
    f.levelStations = maxLevelStations(H, sec, tau);
    if (f.stations < minStations(kParamsDraftCone) || f.levels < minStations(kParamsLine) ||
        f.levelStations < kParamsCircle) {
        f.why = "underdetermined";
        return f;
    }
    f.determinate = true;
    // D-gid-8.29: a band another family with no more parameters certifies is
    // ambiguous; the drafted cone does not claim it.
    if (const char* family = ambiguousWithCheaperFamily(pts, X, Y, H, tau, f)) {
        f.why = family;
        return f;
    }
    // The cone fit's starting point: the axis cylinder in the certificate's
    // frame (a start, not a certificate; the cone's own residual decides).
    const AxisCylinder cyl = fitAxisCylinder(X, Y, tau);
    if (!cyl.fitted) {
        f.why = "fit-singular";
        return f;
    }
    double cx = cyl.cx, cy = cyl.cy;
    const double R = cyl.R;
    // The cone about n: rho = rhoP + k h. Start from the cylinder's centre and
    // the regression of rho on h; Gauss-Newton on the orthogonal distance.
    double rhoP = R, k = 0.0;
    {
        double sh = 0, shh = 0, sr = 0, shr = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const double rho = std::hypot(X[i] - cx, Y[i] - cy);
            sh += H[i];
            shh += H[i] * H[i];
            sr += rho;
            shr += H[i] * rho;
        }
        const double det = (double)N * shh - sh * sh;
        if (std::fabs(det) > gp::Resolution()) {
            k = ((double)N * shr - sh * sr) / det;
            rhoP = (sr - k * sh) / (double)N;
        }
    }
    const double lever = std::max(std::fabs(f.hLo), std::fabs(f.hHi));
    for (int it = 0; it < kDraftFitIters; ++it) {
        const double s = std::sqrt(1.0 + k * k);
        std::array<std::array<double, 4>, 4> JJ{};
        std::array<double, 4> Jr{};
        for (std::size_t i = 0; i < N; ++i) {
            const double dx = X[i] - cx, dy = Y[i] - cy;
            const double rho = std::hypot(dx, dy);
            if (!(rho > gp::Resolution())) continue;
            const double res = (rho - rhoP - k * H[i]) / s;
            const double j[4] = {-dx / rho / s, -dy / rho / s, -1.0 / s,
                                 -H[i] / s - res * k / (s * s)};
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) JJ[r][c] += j[r] * j[c];
                Jr[r] -= j[r] * res;
            }
        }
        std::array<double, 4> d{};
        if (!solveN<4>(JJ, Jr, d)) break;
        cx += d[0];
        cy += d[1];
        rhoP += d[2];
        k += d[3];
        if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(rhoP) || !std::isfinite(k))
            break;
        if (std::hypot(std::hypot(d[0], d[1]), std::hypot(d[2], d[3] * lever)) <= tau) break;
    }
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(rhoP) || !std::isfinite(k)) {
        f.why = "fit-singular";
        return f;
    }
    f.rhoP = rhoP;
    f.k = k;
    f.c0 = p0 + e1 * cx + e2 * cy;
    const double s = std::sqrt(1.0 + k * k);
    // D-gid-8.25 (b). Fewest parameters first: the draft is resolvable only
    // when the cone's radius changes over the band by more than the
    // certificate's diameter; otherwise the band is a cylinder to the
    // certificate, not a cone.
    if (!(std::fabs(k) * (f.hHi - f.hLo) / s > levelSep(tau))) {
        f.why = "cylinder";
        return f;
    }
    if (!((rhoP + k * f.hLo) / s > tau) || !((rhoP + k * f.hHi) / s > tau)) {
        f.why = "apex";
        return f;
    }
    double sum2 = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double d = std::fabs(std::hypot(X[i] - cx, Y[i] - cy) - rhoP - k * H[i]) / s;
        f.maxDev = std::max(f.maxDev, d);
        sum2 += d * d;
    }
    f.rmsDev = std::sqrt(sum2 / (double)N);
    if (!(f.maxDev <= tau)) {
        f.why = "vertex-off-surface";
        return f;
    }
    f.ok = true;
    f.why = "accept";
    return f;
}

// D-gid-8.18. After the sphere corners (refit_segment.cpp). See the file comment.
bool claimDraftedCones(const MeshView& mv, const SegmentParams& /*p*/, const DerivedTols& tol,
                       SegmentWork& work) {
    try {
        const double tau = tauOf(mv);
        if (!(tau > 0.0) || !mv.pts || !mv.tris || !mv.triEdges || mv.nTri == 0) return true;
        const int nTri = (int)mv.nTri;
        const int nEdge = (int)mv.nEdge;
        const bool dg = diagOn();
        auto gidOf = [&](int t, int k) {
            const int gt = mv.compTris ? mv.compTris[t] : t;
            return mv.tris[gt][k];
        };
        auto triArea = [&](int t) {
            const gp_XYZ a = mv.pts[gidOf(t, 0)], b = mv.pts[gidOf(t, 1)], c = mv.pts[gidOf(t, 2)];
            return 0.5 * (b - a).Crossed(c - a).Modulus();
        };

        std::vector<std::array<int, 2>> edgeTris((std::size_t)std::max(nEdge, 0),
                                                 std::array<int, 2>{{-1, -1}});
        for (int t = 0; t < nTri; ++t)
            for (int s = 0; s < 3; ++s) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || e >= nEdge) continue;
                auto& et = edgeTris[(std::size_t)e];
                if (et[0] < 0) et[0] = t;
                else if (et[0] != t) et[1] = t;
            }
        auto otherTri = [&](int e, int t) -> int {
            if (e < 0 || e >= nEdge) return -1;
            const auto& et = edgeTris[(std::size_t)e];
            if (et[0] == t) return et[1];
            if (et[1] == t) return et[0];
            return -1;
        };

        // The pull directions: the normals of the certified planes (an
        // Unclaimed provisional A3 would commit whose own residual is within
        // tau, or an accepted plane within tau), largest first. Two normals
        // closer than tau / diag move no point of the part by tau: one
        // direction.
        struct PullDir {
            gp_XYZ n, p0;
            double area = 0;
            int minTri = 0;
            int src = -1;  // provisional index, or -(accepted index + 2)
        };
        std::vector<char> provPlane(work.provisionals.size(), 0);
        std::vector<PullDir> cand;
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
            const Provisional& pr = work.provisionals[(std::size_t)pi];
            if (pr.tris.empty() || !a3WouldCommit(mv, tol, pr, false) ||
                !anchorCertified(tau, pr.maxVertexDev))
                continue;
            provPlane[(std::size_t)pi] = 1;
            PullDir d;
            d.n = pr.plane.Direction().XYZ();
            d.p0 = pr.plane.Location().XYZ();
            d.area = pr.area;
            d.minTri = *std::min_element(pr.tris.begin(), pr.tris.end());
            d.src = pi;
            cand.push_back(d);
        }
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai) {
            const Region& A = work.accepted[(std::size_t)ai];
            if (A.type != SurfType::Plane || A.tris.empty() ||
                !anchorCertified(tau, A.maxVertexDev))
                continue;
            PullDir d;
            d.n = A.ax.Direction().XYZ();
            d.p0 = A.ax.Location().XYZ();
            for (int t : A.tris)
                if (t >= 0 && t < nTri) d.area += triArea(t);
            d.minTri = *std::min_element(A.tris.begin(), A.tris.end());
            d.src = -(ai + 2);
            cand.push_back(d);
        }
        std::sort(cand.begin(), cand.end(), [](const PullDir& a, const PullDir& b) {
            if (a.area != b.area) return a.area > b.area;
            return a.minTri < b.minTri;
        });
        std::vector<PullDir> dirs;
        const double dirTol = (mv.diag > 0.0) ? tau / mv.diag : 0.0;
        for (const PullDir& d : cand) {
            const double m = d.n.Modulus();
            if (!(m > gp::Resolution())) continue;
            bool dup = false;
            for (const PullDir& q : dirs)
                if ((d.n / m).Crossed(q.n).Modulus() <= dirTol) {
                    dup = true;
                    break;
                }
            if (dup) continue;
            PullDir u = d;
            u.n = d.n / m;
            dirs.push_back(u);
        }

        // Atoms: whole provisionals (Unclaimed, not a certified plane) and
        // whole accepted regions that do not certify their surface.
        struct Atom {
            int prov = -1, acc = -1;
            int minTri = 0;
            std::vector<int> tris;
            std::vector<int> nbrs;  // atom ranks, ascending
        };
        std::vector<Atom> atoms;
        for (int pi = 0; pi < (int)work.provisionals.size(); ++pi) {
            const Provisional& pr = work.provisionals[(std::size_t)pi];
            if (pr.tris.empty() || pr.claim != ProvClaim::Unclaimed || provPlane[(std::size_t)pi])
                continue;
            Atom a;
            a.prov = pi;
            a.tris = pr.tris;
            atoms.push_back(std::move(a));
        }
        std::vector<char> triAccepted((std::size_t)nTri, 0);
        for (const Region& A : work.accepted)
            for (int t : A.tris)
                if (t >= 0 && t < nTri) triAccepted[(std::size_t)t] = 1;
        for (int ai = 0; ai < (int)work.accepted.size(); ++ai) {
            const Region& A = work.accepted[(std::size_t)ai];
            if (A.tris.empty() || !robbableUncertified(A, tau)) continue;
            Atom a;
            a.acc = ai;
            a.tris = A.tris;
            atoms.push_back(std::move(a));
        }
        std::vector<int> atomOf((std::size_t)nTri, -1);
        for (Atom& a : atoms) {
            std::vector<int> keep;
            for (int t : a.tris)
                if (t >= 0 && t < nTri) keep.push_back(t);
            std::sort(keep.begin(), keep.end());
            keep.erase(std::unique(keep.begin(), keep.end()), keep.end());
            a.tris.swap(keep);
            a.minTri = a.tris.empty() ? nTri : a.tris.front();
        }
        // A provisional triangle an accepted region already holds is not free.
        atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
                                   [&](const Atom& a) {
                                       if (a.tris.empty()) return true;
                                       if (a.prov < 0) return false;
                                       for (int t : a.tris)
                                           if (triAccepted[(std::size_t)t]) return true;
                                       return false;
                                   }),
                    atoms.end());
        std::sort(atoms.begin(), atoms.end(), [](const Atom& a, const Atom& b) {
            if (a.minTri != b.minTri) return a.minTri < b.minTri;
            return a.prov < b.prov;
        });
        for (int ri = 0; ri < (int)atoms.size(); ++ri)
            for (int t : atoms[(std::size_t)ri].tris) atomOf[(std::size_t)t] = ri;
        for (int ri = 0; ri < (int)atoms.size(); ++ri) {
            Atom& a = atoms[(std::size_t)ri];
            for (int t : a.tris)
                for (int s = 0; s < 3; ++s) {
                    const int o = otherTri(mv.triEdges[t][s], t);
                    if (o < 0) continue;
                    const int b = atomOf[(std::size_t)o];
                    if (b >= 0 && b != ri) a.nbrs.push_back(b);
                }
            std::sort(a.nbrs.begin(), a.nbrs.end());
            a.nbrs.erase(std::unique(a.nbrs.begin(), a.nbrs.end()), a.nbrs.end());
        }

        auto pointsOf = [&](const std::vector<int>& ranks) {
            std::vector<int> g;
            for (int r : ranks)
                for (int t : atoms[(std::size_t)r].tris)
                    for (int k = 0; k < 3; ++k) g.push_back(gidOf(t, k));
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            std::vector<gp_XYZ> P;
            P.reserve(g.size());
            for (int v : g) P.push_back(mv.pts[v]);
            return P;
        };

        // D-gid-8.18. The claim. Region convention of a partial coaxial cone
        // (D-fc-5.5, refit_relational_revolve.cpp emitPartialCone): Location on
        // the R_lo rim's axis point, Direction the pull direction, radius R_lo,
        // vMin / vMax = R_lo / R_hi, maxVertexSnap the signed height from the
        // R_lo rim to the R_hi rim, XDirection at the start of the admitted
        // arc (u in [0, span]). One annulus about the axis (two boundary
        // loops, each winding once, D-fc-5.7 O-1) is a closed frustum.
        int nClaimed = 0;
        std::vector<int> dissolved;
        auto emitCone = [&](const std::vector<int>& band, const std::vector<int>& tris,
                            const DraftedConeFit& fit, const PullDir& D) -> const char* {
            const gp_XYZ a = D.n;
            gp_XYZ e1, e2;
            frameOf(a, e1, e2);
            const double rA = fit.rhoP + fit.k * fit.hLo, rB = fit.rhoP + fit.k * fit.hHi;
            const bool loAtMin = rA <= rB;
            const double Rlo = loAtMin ? rA : rB, Rhi = loAtMin ? rB : rA;
            const double zLo = loAtMin ? fit.hLo : fit.hHi;
            const double height = fit.hHi - fit.hLo;
            const double hSigned = loAtMin ? height : -height;
            if (!(Rhi - Rlo > tau) || !(Rlo > tau)) return "degenerate";
            // The axis point at h = 0 is c0; heights are measured from D.p0.
            const gp_XYZ axisAt0 = fit.c0;
            auto psiOf = [&](const gp_XYZ& x) {
                const gp_XYZ d = x - axisAt0;
                double u = std::atan2(d.Dot(e2), d.Dot(e1));
                if (u < 0.0) u += 2.0 * kPi;
                return u;
            };
            // Boundary cycles of the band and their windings about the axis.
            std::vector<char> inBand((std::size_t)nTri, 0);
            for (int t : tris) inBand[(std::size_t)t] = 1;
            std::vector<std::pair<int, int>> half;  // directed boundary edges (ga -> gb)
            for (int t : tris)
                for (int sd = 0; sd < 3; ++sd) {
                    const int o = otherTri(mv.triEdges[t][sd], t);
                    if (o >= 0 && inBand[(std::size_t)o]) continue;
                    const int e = mv.triEdges[t][sd];
                    if (e < 0 || e >= nEdge || !mv.compEdges) return "no-edge-table";
                    const std::pair<int, int>& ce = mv.compEdges[e];
                    const int ga = mv.compVtx ? mv.compVtx[ce.first] : ce.first;
                    const int gb = mv.compVtx ? mv.compVtx[ce.second] : ce.second;
                    int ka = -1, kb = -1;
                    for (int k = 0; k < 3; ++k) {
                        if (gidOf(t, k) == ga) ka = k;
                        if (gidOf(t, k) == gb) kb = k;
                    }
                    if (ka < 0 || kb < 0) return "no-edge-table";
                    if ((ka + 1) % 3 == kb) half.emplace_back(ga, gb);
                    else half.emplace_back(gb, ga);
                }
            std::sort(half.begin(), half.end());
            for (std::size_t i = 1; i < half.size(); ++i)
                if (half[i].first == half[i - 1].first) return "pinch";
            std::vector<int> winds;
            {
                std::vector<char> used(half.size(), 0);
                for (std::size_t i0 = 0; i0 < half.size(); ++i0) {
                    if (used[i0]) continue;
                    std::vector<double> az;
                    std::size_t i = i0;
                    bool closedLoop = false;
                    for (std::size_t guard = 0; guard <= half.size(); ++guard) {
                        used[i] = 1;
                        const gp_XYZ d = mv.pts[half[i].first] - axisAt0;
                        az.push_back(std::atan2(d.Dot(e2), d.Dot(e1)));
                        const int nx = half[i].second;
                        auto it = std::lower_bound(half.begin(), half.end(),
                                                   std::make_pair(nx, -1));
                        if (it == half.end() || it->first != nx) break;
                        i = (std::size_t)(it - half.begin());
                        if (i == i0) {
                            closedLoop = true;
                            break;
                        }
                    }
                    if (!closedLoop) return "open-boundary";
                    winds.push_back(windingNumber(az));
                }
            }
            const bool closed = annulusAboutAxis(1, winds);
            if (!closed)
                for (int w : winds)
                    if (w != 0) return "not-one-face";
            // The admitted arc: the complement of the largest vertex-azimuth gap.
            std::vector<int> gv;
            for (int t : tris)
                for (int k = 0; k < 3; ++k) gv.push_back(gidOf(t, k));
            std::sort(gv.begin(), gv.end());
            gv.erase(std::unique(gv.begin(), gv.end()), gv.end());
            double uLo = 0.0, span = 2.0 * kPi;
            if (!closed) {
                std::vector<double> u;
                u.reserve(gv.size());
                for (int g : gv) u.push_back(psiOf(mv.pts[g]));
                std::sort(u.begin(), u.end());
                std::size_t ig = u.size() - 1;
                double gap = (u.front() + 2.0 * kPi) - u.back();
                for (std::size_t i = 0; i + 1 < u.size(); ++i)
                    if (u[i + 1] - u[i] > gap) {
                        gap = u[i + 1] - u[i];
                        ig = i;
                    }
                uLo = u[(ig + 1) % u.size()];
                span = 2.0 * kPi - gap;
                if (!(span > 0.0)) return "degenerate";
            } else {
                // A real mesh vertex azimuth: a corner of the band's lowest triangle.
                uLo = psiOf(mv.pts[gidOf(tris.front(), 0)]);
            }
            const gp_XYZ xDir = e1 * std::cos(uLo) + e2 * std::sin(uLo);
            const gp_XYZ yDir = a.Crossed(xDir);
            const gp_XYZ loc = axisAt0 + a * zLo;
            // Orientation (I9), the exact chord volume (revolvedChordTerm, the
            // apex and partial cones' own) and the side count of the widest
            // facet the band holds.
            double sigma = 0.0, dVol = 0.0, maxStep = 0.0;
            for (int t : tris) {
                const gp_XYZ p0 = mv.pts[gidOf(t, 0)], p1 = mv.pts[gidOf(t, 1)],
                             p2 = mv.pts[gidOf(t, 2)];
                const gp_XYZ nn = (p1 - p0).Crossed(p2 - p0);
                gp_XYZ rad = (p0 + p1 + p2) / 3.0 - loc;
                rad -= a * rad.Dot(a);
                const double rm = rad.Modulus();
                if (rm > gp::Resolution()) sigma += 0.5 * nn.Dot(rad / rm);
                dVol += revolvedChordTerm(p0, p1, p2, loc, a, xDir, yDir, Rlo,
                                          (Rhi - Rlo) / hSigned);
                // The band's oblique rim edges (its open ends and a rim on a
                // plane not normal to the axis): coneRimTerm.
                const int sg = revolvedParamSign(p0, p1, p2, loc, a, xDir, yDir);
                for (int sd = 0; sd < 3; ++sd) {
                    const int o = otherTri(mv.triEdges[t][sd], t);
                    if (o >= 0 && inBand[(std::size_t)o]) continue;
                    const std::pair<int, int>& ce = mv.compEdges[mv.triEdges[t][sd]];
                    const int ga = mv.compVtx ? mv.compVtx[ce.first] : ce.first;
                    const int gb = mv.compVtx ? mv.compVtx[ce.second] : ce.second;
                    int ka = -1, kb = -1;
                    for (int k = 0; k < 3; ++k) {
                        if (gidOf(t, k) == ga) ka = k;
                        if (gidOf(t, k) == gb) kb = k;
                    }
                    const bool fwd = (ka + 1) % 3 == kb;
                    const gp_XYZ& pa = mv.pts[fwd ? ga : gb];
                    const gp_XYZ& pb = mv.pts[fwd ? gb : ga];
                    dVol += (double)sg * coneRimTerm(pa, pb, loc, a, xDir, yDir, Rlo,
                                                     (Rhi - Rlo) / hSigned);
                }
                const gp_XYZ P[3] = {p0, p1, p2};
                for (int k = 0; k < 3; ++k) {
                    double du = std::fabs(psiOf(P[k]) - psiOf(P[(k + 1) % 3]));
                    if (du > kPi) du = 2.0 * kPi - du;
                    maxStep = std::max(maxStep, du);
                }
            }
            const int nSides =
                (maxStep > 0.0) ? std::max(3, (int)std::ceil(2.0 * kPi / maxStep)) : 3;
            double sum2 = 0.0;
            {
                const double s = std::sqrt(1.0 + fit.k * fit.k);
                for (int g : gv) {
                    const gp_XYZ d = mv.pts[g] - D.p0;
                    const double h = d.Dot(a);
                    gp_XYZ r = mv.pts[g] - axisAt0;
                    r -= a * r.Dot(a);
                    const double dv = (r.Modulus() - fit.rhoP - fit.k * h) / s;
                    sum2 += dv * dv;
                }
            }
            Region r;
            r.id = 0;
            r.type = SurfType::Cone;
            r.origin = Origin::ChamferCone;
            r.ax = gp_Ax3(gp_Pnt(loc.X(), loc.Y(), loc.Z()), gp_Dir(a.X(), a.Y(), a.Z()),
                          gp_Dir(xDir.X(), xDir.Y(), xDir.Z()));
            r.radius = Rlo;
            r.uMin = 0.0;
            r.uMax = span;
            r.vMin = Rlo;
            r.vMax = Rhi;
            r.closed360 = closed;
            r.outwardNormal = sigma > 0.0;
            r.tris = tris;
            r.maxVertexDev = fit.maxDev;
            r.rmsVertexDev = gv.empty() ? 0.0 : std::sqrt(sum2 / (double)gv.size());
            r.chordSagitta = chordSagitta(Rlo, nSides);
            r.nSides = nSides;
            r.dVolPredicted = (sigma > 0.0 ? 1.0 : -1.0) * dVol;
            r.maxVertexSnap = hSigned;
            r.reject = Reject::None;
            r.builtAs = BuiltAs::NotBuilt;
            r.filletNbrA = -1;
            r.filletNbrB = -1;
            // Every atom is whole inside the band: its provisional is consumed;
            // an uncertified region it robs dissolves (D-fc-4.3).
            for (int q : band) {
                const Atom& at = atoms[(std::size_t)q];
                if (at.prov >= 0)
                    work.provisionals[(std::size_t)at.prov].claim = ProvClaim::ConsumedCylinder;
                else if (at.acc >= 0) {
                    if (dg)
                        diagPrintf("DIAG_DRAFT_ROB region=%d origin=%d radius=%.9g devTau=%.4f "
                                   "nTri=%zu\n",
                                   at.acc, (int)work.accepted[(std::size_t)at.acc].origin,
                                   work.accepted[(std::size_t)at.acc].radius,
                                   work.accepted[(std::size_t)at.acc].maxVertexDev / tau,
                                   work.accepted[(std::size_t)at.acc].tris.size());
                    work.accepted[(std::size_t)at.acc].tris.clear();
                    dissolved.push_back(at.acc);
                }
            }
            if (dg)
                diagPrintf("DIAG_DRAFT_REGION seed=%d closed=%d Rlo=%.9g Rhi=%.9g h=%.9g "
                           "span=%.6f nSides=%d outward=%d dVol=%.9g\n",
                           atoms[(std::size_t)band.front()].minTri, closed ? 1 : 0, Rlo, Rhi,
                           hSigned, span, nSides, sigma > 0.0 ? 1 : 0, r.dVolPredicted);
            work.accepted.push_back(std::move(r));
            return "accept";
        };

        // Seed outcomes by why (diag census only).
        const char* kWhys[] = {"underdetermined", "plane", "cylinder", "apex",
                               "vertex-off-surface", "fit-singular", "accept", "sphere"};
        long whyCount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        auto tallyWhy = [&](const char* w) {
            for (int i = 0; i < 8; ++i)
                if (std::strcmp(w, kWhys[i]) == 0) {
                    ++whyCount[i];
                    return;
                }
        };
        // D-gid-8.25 (d). Early-out: the directions that can anchor a band of
        // each atom component (a band grows only inside its component). Under
        // a direction, a band's level is a subset of one of the component's
        // levels (levelClusters: a gap chain in the band is one in the
        // component), so when no component level holds kParamsCircle vertices,
        // no band in the component holds kParamsCircle stations on one level
        // and every certificate there is "underdetermined". Exact: the
        // directions dropped are the ones that cannot certify.
        std::vector<int> compOf(atoms.size(), -1);
        int nComp = 0;
        for (int a0 = 0; a0 < (int)atoms.size(); ++a0) {
            if (compOf[(std::size_t)a0] >= 0) continue;
            std::vector<int> st = {a0};
            compOf[(std::size_t)a0] = nComp;
            while (!st.empty()) {
                const int u = st.back();
                st.pop_back();
                for (int o : atoms[(std::size_t)u].nbrs)
                    if (compOf[(std::size_t)o] < 0) {
                        compOf[(std::size_t)o] = nComp;
                        st.push_back(o);
                    }
            }
            ++nComp;
        }
        std::vector<std::vector<int>> compDirs((std::size_t)nComp);
        int nCompSkipped = 0;
        {
            std::vector<std::vector<int>> compAtoms((std::size_t)nComp);
            for (int a0 = 0; a0 < (int)atoms.size(); ++a0)
                compAtoms[(std::size_t)compOf[(std::size_t)a0]].push_back(a0);
            for (int c = 0; c < nComp; ++c) {
                const std::vector<gp_XYZ> P = pointsOf(compAtoms[(std::size_t)c]);
                if (P.size() < (std::size_t)minStations(kParamsDraftCone)) {
                    ++nCompSkipped;
                    continue;
                }
                std::vector<double> h(P.size());
                for (int k = 0; k < (int)dirs.size(); ++k) {
                    for (std::size_t i = 0; i < P.size(); ++i)
                        h[i] = (P[i] - dirs[(std::size_t)k].p0).Dot(dirs[(std::size_t)k].n);
                    std::size_t most = 0;
                    for (const auto& lv : levelClusters(h, tau)) most = std::max(most, lv.size());
                    if (most >= (std::size_t)kParamsCircle) compDirs[(std::size_t)c].push_back(k);
                }
                if (compDirs[(std::size_t)c].empty()) ++nCompSkipped;
            }
        }
        std::vector<char> taken(atoms.size(), 0);
        long nFits = 0;
        int nSeeds = 0, nBands = 0, nCert = 0;
        for (int a0 = 0; a0 < (int)atoms.size(); ++a0) {
            if (taken[(std::size_t)a0]) continue;
            const std::vector<int>& cdirs = compDirs[(std::size_t)compOf[(std::size_t)a0]];
            if (cdirs.empty()) continue;
            for (int b0 : atoms[(std::size_t)a0].nbrs) {
                if (taken[(std::size_t)a0]) break;
                if (taken[(std::size_t)b0]) continue;
                std::vector<int> band = {std::min(a0, b0), std::max(a0, b0)};
                const std::vector<gp_XYZ> seedPts = pointsOf(band);
                if (seedPts.size() < (std::size_t)minStations(kParamsDraftCone)) continue;
                if (planeResidual(seedPts) <= tau) continue;
                ++nSeeds;
                int di = -1;
                DraftedConeFit fit;
                for (int k : cdirs) {
                    ++nFits;
                    fit = certifyDraftedCone(seedPts, dirs[(std::size_t)k].n,
                                             dirs[(std::size_t)k].p0, tau);
                    tallyWhy(fit.why);
                    if (fit.ok) {
                        di = k;
                        break;
                    }
                }
                if (di < 0) continue;
                // Grow: an adjacent atom joins whole while the refit band
                // keeps the certificate. Ascending rank, to a fixed point.
                const PullDir& D = dirs[(std::size_t)di];
                for (bool grew = true; grew;) {
                    grew = false;
                    std::vector<int> front;
                    for (int r : band)
                        for (int o : atoms[(std::size_t)r].nbrs)
                            if (!taken[(std::size_t)o] &&
                                !std::binary_search(band.begin(), band.end(), o))
                                front.push_back(o);
                    std::sort(front.begin(), front.end());
                    front.erase(std::unique(front.begin(), front.end()), front.end());
                    for (int o : front) {
                        std::vector<int> trial = band;
                        trial.insert(std::upper_bound(trial.begin(), trial.end(), o), o);
                        ++nFits;
                        const DraftedConeFit f2 = certifyDraftedCone(pointsOf(trial), D.n, D.p0, tau);
                        if (!f2.ok) continue;
                        band.swap(trial);
                        fit = f2;
                        grew = true;
                    }
                }
                ++nBands;
                // D-fc-1.11 (ii)(a): without any one atom the band still meets
                // the station and level counts.
                const char* why = "accept";
                {
                    gp_XYZ e1, e2;
                    frameOf(D.n, e1, e2);
                    for (int r : band) {
                        if (band.size() < 2) break;
                        std::vector<int> rest;
                        for (int q : band)
                            if (q != r) rest.push_back(q);
                        const std::vector<gp_XYZ> P = pointsOf(rest);
                        std::vector<ProfPt> sec(P.size());
                        std::vector<double> H(P.size());
                        for (std::size_t i = 0; i < P.size(); ++i) {
                            const gp_XYZ d = P[i] - D.p0;
                            sec[i].rho = d.Dot(e1);
                            sec[i].z = d.Dot(e2);
                            H[i] = d.Dot(D.n);
                        }
                        if (countProfPts(sec, tau) < minStations(kParamsDraftCone) ||
                            countLevelClusters(H, tau) < minStations(kParamsLine) ||
                            maxLevelStations(H, sec, tau) < kParamsCircle) {
                            why = "underdetermined";
                            break;
                        }
                    }
                }
                // One side of the wall: every triangle faces the same way about
                // the axis.
                std::vector<int> tris;
                for (int r : band)
                    tris.insert(tris.end(), atoms[(std::size_t)r].tris.begin(),
                                atoms[(std::size_t)r].tris.end());
                std::sort(tris.begin(), tris.end());
                int nPos = 0, nNeg = 0, nRob = 0;
                double area = 0.0;
                for (int r : band)
                    if (atoms[(std::size_t)r].acc >= 0) ++nRob;
                for (int t : tris) {
                    const gp_XYZ a = mv.pts[gidOf(t, 0)], b = mv.pts[gidOf(t, 1)],
                                 c = mv.pts[gidOf(t, 2)];
                    const gp_XYZ nn = (b - a).Crossed(c - a);
                    area += 0.5 * nn.Modulus();
                    gp_XYZ rad = (a + b + c) / 3.0 - fit.c0;
                    rad -= D.n * rad.Dot(D.n);
                    const double s = nn.Dot(rad);
                    if (s > 0.0) ++nPos;
                    else if (s < 0.0) ++nNeg;
                }
                if (nPos > 0 && nNeg > 0 && std::strcmp(why, "accept") == 0) why = "not-one-side";
                if (std::strcmp(why, "accept") == 0) ++nCert;
                if (dg) {
                    const double s = std::sqrt(1.0 + fit.k * fit.k);
                    diagPrintf("DIAG_DRAFT_BAND seed=%d dir=%d plane=%d nAtoms=%zu robbed=%d "
                               "nTri=%zu nV=%zu stations=%d levels=%d levelSt=%d rhoLo=%.9g rhoHi=%.9g "
                               "hLo=%.9g hHi=%.9g draftDeg=%.6f axis=%.9f,%.9f,%.9f "
                               "c0=%.9f,%.9f,%.9f maxDevTau=%.4f planeDevTau=%.1f cylDevTau=%.1f "
                               "area=%.6g why=%s\n",
                               atoms[(std::size_t)band.front()].minTri, di, D.src, band.size(),
                               nRob, tris.size(), pointsOf(band).size(), fit.stations,
                               fit.levels, fit.levelStations, fit.rhoP + fit.k * fit.hLo, fit.rhoP + fit.k * fit.hHi,
                               fit.hLo, fit.hHi, std::atan(fit.k) * 180.0 / kPi, D.n.X(),
                               D.n.Y(), D.n.Z(), fit.c0.X(), fit.c0.Y(), fit.c0.Z(),
                               fit.maxDev / tau, fit.planeDev / tau, fit.cylDev / tau, area, why);
                    std::fprintf(stderr, "DIAG_DRAFT_TRIS seed=%d tris=",
                                 atoms[(std::size_t)band.front()].minTri);
                    for (int t : tris) std::fprintf(stderr, "%d,", t);
                    std::fprintf(stderr, "\n");
                    std::fprintf(stderr, "DIAG_DRAFT_PTS seed=%d pts=",
                                 atoms[(std::size_t)band.front()].minTri);
                    for (const gp_XYZ& q : pointsOf(band))
                        std::fprintf(stderr, "%.9f,%.9f,%.9f;", q.X(), q.Y(), q.Z());
                    std::fprintf(stderr, "\n");
                    (void)s;
                }
                // A refused band's atoms are not seeds again (B1's rule); a
                // certified band is claimed as one cone face.
                for (int r : band) taken[(std::size_t)r] = 1;
                if (std::strcmp(why, "accept") == 0) {
                    const char* ew = emitCone(band, tris, fit, D);
                    if (dg)
                        diagPrintf("DIAG_DRAFT_CLAIM seed=%d nTri=%zu why=%s\n",
                                   atoms[(std::size_t)band.front()].minTri, tris.size(), ew);
                    if (std::strcmp(ew, "accept") == 0) ++nClaimed;
                }
                break;
            }
        }
        if (dg)
            diagPrintf("DIAG_DRAFT_SEEDWHY under=%ld plane=%ld cyl=%ld sphere=%ld apex=%ld "
                       "off=%ld singular=%ld accept=%ld\n",
                       whyCount[0], whyCount[1], whyCount[2], whyCount[7], whyCount[3],
                       whyCount[4], whyCount[5], whyCount[6]);
        if (dg)
            diagPrintf("DIAG_DRAFT_COMPS comps=%d skipped=%d\n", nComp, nCompSkipped);
        if (dg)
            diagPrintf("DIAG_DRAFT_SUM nTri=%d tau=%.9g dirs=%zu atoms=%zu seeds=%d fits=%ld "
                       "bands=%d certified=%d claimed=%d\n",
                       nTri, tau, dirs.size(), atoms.size(), nSeeds, nFits, nBands, nCert,
                       nClaimed);
        compactDissolved(work, dissolved);
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
