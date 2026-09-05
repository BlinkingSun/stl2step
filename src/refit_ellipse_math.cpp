// stl2step ellipse math — the closed forms behind the 1.4.0 ELLIPSE bind class
// (lane ellipse-math, SPEC-ellipse PART A, D-140-9 §2).  Signatures and the
// derivations are published in src/refit_ellipse_math.hpp; the proofs and the
// unit table are in _team/reports/ellipse-math.md.
//
// This is deliberately NOT part of src/refit_math.cpp, for the reason
// refit_cone_math.cpp (0d67eba) and refit_torus_math.cpp (048598a) are not:
// that file is one of the five P1 sources the D5.3 include allowlist covers,
// and the allowlist — a gate instrument, not a style rule — admits no project
// header beyond refit.hpp / refit_internal.hpp and no bare gp.hxx.  Nothing
// here touches MeshView, Region or any engine state: it is pure gp_ geometry,
// which is also what lets tests/unit/ellipse_math_test.cpp compile this one
// file and certify the math without linking the engine.
//
// Everything below is closed form.  The plane class is an EXACT supremum (the
// quartic critical set plus the endpoints); the cylinder class is a certified
// UPPER BOUND whose every term is finite and fixed by the pcurve's own degree
// (D-130-2 door (b), "certified-sample").  Nothing samples the surface to find
// a maximum, nothing iterates to convergence, and no tolerance is a literal:
// Precision::Angular / Confusion / PConfusion and gp::Resolution only.  The
// numeric literals are the algebra's own integers (Ferrari's 3/8, 27, 256, the
// eight equispaced shift samples) and pi.
//
// SPDX-License-Identifier: MIT

#include "refit_ellipse_math.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <BSplCLib.hxx>
#include <gp.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {

const char* const kEllipseClausePlaneParam = "unhandled-ellipse-plane-param";
const char* const kEllipseClausePlaneDegenerate = "unhandled-ellipse-plane-degenerate";
const char* const kEllipseClauseCylOnCyl = "unhandled-ellipse-cyl-oncyl";
const char* const kEllipseClauseCylUAffine = "unhandled-ellipse-cyl-uaffine";
const char* const kEllipseClauseCylSpan = "unhandled-ellipse-cyl-span";
const char* const kEllipseClauseCylDegenerate = "unhandled-ellipse-cyl-degenerate";
const char* const kEllipseClauseGrazing = "unhandled-ellipse-grazing";

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr int kShiftSamples = 8;  // equispaced samples of h that pick the tan(t/2) shift

bool finite3(const gp_XYZ& v) {
    return std::isfinite(v.X()) && std::isfinite(v.Y()) && std::isfinite(v.Z());
}

bool tolsOk(double angTol, double linTol) {
    return std::isfinite(angTol) && angTol >= 0.0 && std::isfinite(linTol) && linTol >= 0.0;
}

void setRefused(EllipseDevClass* clsOut, const char** clauseOut, const char* clause) {
    if (clsOut) *clsOut = EllipseDevClass::Unhandled;
    if (clauseOut) *clauseOut = clause;
}

void setAdmitted(EllipseDevClass* clsOut, const char** clauseOut, EllipseDevClass c) {
    if (clsOut) *clsOut = c;
    if (clauseOut) *clauseOut = ellipseDevClassName(c);
}

// ---------------------------------------------------------------------------
// the trig polynomial  g(t) = |W + A cos t + B sin t|^2  and h(t) = g'(t)/2
// ---------------------------------------------------------------------------
struct TrigForm {
    gp_XYZ W, A, B;
    double P, Q, Ra, Rb;  // |B|^2 - |A|^2, A.B, W.A, W.B

    TrigForm(const gp_XYZ& w, const gp_XYZ& a, const gp_XYZ& b) : W(w), A(a), B(b) {
        P = B.SquareModulus() - A.SquareModulus();
        Q = A.Dot(B);
        Ra = W.Dot(A);
        Rb = W.Dot(B);
    }
    double g(double t) const {
        const gp_XYZ r = W + A * std::cos(t) + B * std::sin(t);
        return r.SquareModulus();
    }
    double h(double t) const {
        const double c = std::cos(t), s = std::sin(t);
        return P * s * c + Q * (c * c - s * s) - Ra * s + Rb * c;
    }
    double hPrime(double t) const {
        const double c = std::cos(t), s = std::sin(t);
        return P * (c * c - s * s) - 4.0 * Q * s * c - Ra * c - Rb * s;
    }
};

// ---------------------------------------------------------------------------
// closed-form real roots: quadratic, cubic (largest real root), quartic
// ---------------------------------------------------------------------------

// y^2 + beta y + gamma = 0.  A negative discriminant (a complex pair, or a real
// double root lost to rounding) contributes the real part -beta/2 as a
// CANDIDATE: an extra candidate never changes the maximum of a set that
// already contains every critical point, and if the pair was real-and-merged
// the Newton polish recovers it from its midpoint.
int quadraticCandidates(double beta, double gamma, double* out) {
    const double disc = beta * beta - 4.0 * gamma;
    if (!(disc >= 0.0)) {
        out[0] = -0.5 * beta;
        return 1;
    }
    const double s = std::sqrt(disc);
    out[0] = 0.5 * (-beta + s);
    out[1] = 0.5 * (-beta - s);
    return 2;
}

// Largest real root of m^3 + c2 m^2 + c1 m + c0 = 0 (Cardano / trigonometric).
double cubicLargestRealRoot(double c2, double c1, double c0) {
    const double sh = c2 / 3.0;                        // m = z - sh
    const double P = c1 - c2 * c2 / 3.0;
    const double Q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    const double D = Q * Q / 4.0 + P * P * P / 27.0;
    double z;
    if (D <= 0.0) {
        // three real roots; the largest is k = 0 of the trigonometric form
        const double rho = std::sqrt(std::max(0.0, -P / 3.0));
        if (rho <= 0.0) {
            z = 0.0;
        } else {
            double arg = (Q / 2.0) / (rho * rho * rho);  // = 3Q/(2P) sqrt(-3/P), sign-safe form
            arg = std::max(-1.0, std::min(1.0, -arg));
            z = 2.0 * rho * std::cos(std::acos(arg) / 3.0);
        }
    } else {
        const double sq = std::sqrt(D);
        const double u = std::cbrt(-Q / 2.0 - (Q >= 0.0 ? sq : -sq));
        z = (u == 0.0) ? 0.0 : u - P / (3.0 * u);
    }
    return z - sh;
}

// Real-root CANDIDATES of the monic quartic x^4 + a x^3 + b x^2 + c x + d
// (Ferrari).  Every real root is among the candidates; a candidate that is not
// a root is harmless (see quadraticCandidates).
int quarticCandidates(double a, double b, double c, double d, double* out) {
    const double a2 = a * a;
    const double p = b - 3.0 * a2 / 8.0;
    const double q = c - 0.5 * a * b + a2 * a / 8.0;
    const double r = d - 0.25 * a * c + a2 * b / 16.0 - 3.0 * a2 * a2 / 256.0;
    double ys[4];
    int ny = 0;
    bool biquadratic = (q == 0.0);
    if (!biquadratic) {
        // resolvent  m^3 + p m^2 + (p^2/4 - r) m - q^2/8 = 0  has a root m > 0
        const double m = cubicLargestRealRoot(p, p * p / 4.0 - r, -q * q / 8.0);
        if (m > 0.0 && std::isfinite(m)) {
            const double s = std::sqrt(2.0 * m);
            const double k = 0.5 * p + m;
            ny += quadraticCandidates(-s, k + q / (2.0 * s), ys + ny);
            ny += quadraticCandidates(s, k - q / (2.0 * s), ys + ny);
        } else {
            biquadratic = true;  // q^2 below rounding: y^4 + p y^2 + r
        }
    }
    if (biquadratic) {
        double y2[2];
        const int n2 = quadraticCandidates(p, r, y2);
        for (int i = 0; i < n2; ++i) {
            if (y2[i] >= 0.0) {
                const double s = std::sqrt(y2[i]);
                ys[ny++] = s;
                ys[ny++] = -s;
            } else {
                ys[ny++] = 0.0;  // real part of the complex pair
            }
        }
    }
    int n = 0;
    for (int i = 0; i < ny; ++i)
        if (std::isfinite(ys[i])) out[n++] = ys[i] - a / 4.0;
    return n;
}

// sup over [f, l] of sqrt(g) — the EXACT supremum: endpoints plus every
// critical point (the tau-quartic's real roots, translated by 2 pi into range).
double trigSupOnRange(const TrigForm& tf, double f, double l) {
    const double twoPi = 2.0 * kPi;
    // g is 2 pi-periodic: over any interval of length >= 2 pi its value set is
    // the whole range, so the supremum over [f, l] is the one over [f, f + 2 pi]
    l = std::min(l, f + twoPi);
    double gMax = std::max(tf.g(f), tf.g(l));

    // Threshold-free shift: put the lost point tau = inf (t0 + pi) where |h| is
    // largest among eight equispaced samples.  A degree-2 trig polynomial with
    // eight zeros in a period is identically zero, so |h| == 0 there means g is
    // constant and the endpoints already carry the supremum.
    double tStar = 0.0, hStar = -1.0;
    for (int k = 0; k < kShiftSamples; ++k) {
        const double ts = (double)k * twoPi / (double)kShiftSamples;
        const double hk = std::fabs(tf.h(ts));
        if (hk > hStar) {
            hStar = hk;
            tStar = ts;
        }
    }
    if (!(hStar > 0.0)) return std::sqrt(gMax);
    const double t0 = tStar - kPi;

    // A cos t + B sin t = A' cos th + B' sin th,  th = t - t0
    const double c0 = std::cos(t0), s0 = std::sin(t0);
    const gp_XYZ Ap = tf.A * c0 + tf.B * s0;
    const gp_XYZ Bp = tf.A * (-s0) + tf.B * c0;
    const TrigForm sh(tf.W, Ap, Bp);
    // (1 + tau^2)^2 h(th) = c4 tau^4 + c3 tau^3 + c2 tau^2 + c1 tau + c0
    const double q4 = sh.Q - sh.Rb;
    const double q3 = -2.0 * (sh.P + sh.Ra);
    const double q2 = -6.0 * sh.Q;
    const double q1 = 2.0 * (sh.P - sh.Ra);
    const double q0 = sh.Q + sh.Rb;
    if (q4 == 0.0 || !std::isfinite(q4)) return std::sqrt(gMax);  // |h| below the floor: constant g

    double taus[4];
    const int nr = quarticCandidates(q3 / q4, q2 / q4, q1 / q4, q0 / q4, taus);

    // candidates in t: the roots, the lost point, and the shift samples
    double cand[4 + 1 + kShiftSamples];
    int nc = 0;
    for (int i = 0; i < nr; ++i) cand[nc++] = t0 + 2.0 * std::atan(taus[i]);
    cand[nc++] = t0 + kPi;
    for (int k = 0; k < kShiftSamples; ++k) cand[nc++] = (double)k * twoPi / (double)kShiftSamples;

    for (int i = 0; i < nc; ++i) {
        double t = cand[i];
        if (!std::isfinite(t)) continue;
        // three fixed Newton steps on h (removes Ferrari's cancellation; the
        // supremum is second-order insensitive to the root's position).  A step
        // longer than half a period is not a polish: the seed is kept as is.
        for (int it = 0; it < 3; ++it) {
            const double hp = tf.hPrime(t);
            if (hp == 0.0) break;
            const double tn = t - tf.h(t) / hp;
            if (!std::isfinite(tn) || std::fabs(tn - t) > kPi) break;
            t = tn;
        }
        // every translate t + 2 pi k inside [f, l]  (at most two, l - f <= 2 pi)
        const double kLo = std::ceil((f - t) / twoPi);
        const double kHi = std::floor((l - t) / twoPi);
        if (!std::isfinite(kLo) || !std::isfinite(kHi) || kHi < kLo) continue;
        const int nk = std::min(2, (int)(kHi - kLo) + 1);
        for (int j = 0; j < nk; ++j) {
            const double tk = t + (kLo + (double)j) * twoPi;
            if (tk < f || tk > l) continue;
            gMax = std::max(gMax, tf.g(tk));
        }
    }
    return std::sqrt(gMax);
}

// ---------------------------------------------------------------------------
// cylinder frame helpers
// ---------------------------------------------------------------------------
struct CylFrame {
    gp_XYZ O, X, Y, D;
    double R;
    explicit CylFrame(const gp_Cylinder& cyl)
        : O(cyl.Location().XYZ()),
          X(cyl.Position().XDirection().XYZ()),
          Y(cyl.Position().YDirection().XYZ()),
          D(cyl.Position().Direction().XYZ()),
          R(cyl.Radius()) {}
    bool ok() const { return finite3(O) && finite3(X) && finite3(Y) && finite3(D) && std::isfinite(R); }
};

struct SectionData {
    double theta0 = 0.0, sigma = 1.0, z0 = 0.0, lambda = 0.0, psi = 0.0;
    double zc = 0.0, zs = 0.0;  // z(t) = z0 + zc cos t + zs sin t
    double rhoRes = 0.0;        // |centre off axis| + sigma_max(M - M0): |rho - R| <= rhoRes
};

// Step 1 (SPEC A.2).  False when E is not a section of cyl: centre and radii
// within linTol, Uperp/Vperp orthogonal within angTol (cosine of their angle).
bool sectionOf(const gp_Cylinder& cyl, const gp_Elips& el, double angTol, double linTol, SectionData& sd) {
    const CylFrame cf(cyl);
    if (!cf.ok() || !(cf.R > gp::Resolution())) return false;
    const double a = el.MajorRadius(), b = el.MinorRadius();
    if (!std::isfinite(a) || !std::isfinite(b) || !(a > gp::Resolution()) || !(b > gp::Resolution()))
        return false;
    const gp_XYZ C = el.Location().XYZ();
    const gp_XYZ U = el.XAxis().Direction().XYZ();
    const gp_XYZ V = el.YAxis().Direction().XYZ();
    if (!finite3(C) || !finite3(U) || !finite3(V)) return false;

    const gp_XYZ w = C - cf.O;
    sd.z0 = w.Dot(cf.D);
    const gp_XYZ wPerp = w - cf.D * sd.z0;
    const double dCentre = wPerp.Modulus();
    if (!(dCentre <= linTol)) return false;

    const double ud = U.Dot(cf.D), vd = V.Dot(cf.D);
    const gp_XYZ uPerp = U - cf.D * ud;
    const gp_XYZ vPerp = V - cf.D * vd;
    const double au = a * uPerp.Modulus(), bv = b * vPerp.Modulus();
    if (!(au > gp::Resolution()) || !(bv > gp::Resolution())) return false;
    if (!(std::fabs(au - cf.R) <= linTol)) return false;
    if (!(std::fabs(bv - cf.R) <= linTol)) return false;
    if (!(std::fabs(uPerp.Dot(vPerp)) <= angTol * uPerp.Modulus() * vPerp.Modulus())) return false;

    sd.theta0 = std::atan2(uPerp.Dot(cf.Y), uPerp.Dot(cf.X));
    sd.sigma = (uPerp.Crossed(vPerp).Dot(cf.D) >= 0.0) ? 1.0 : -1.0;
    sd.zc = a * ud;
    sd.zs = b * vd;
    sd.lambda = std::hypot(sd.zc, sd.zs);
    sd.psi = std::atan2(sd.zs, sd.zc);

    // the identity's residual, exact in closed form (header §3)
    const double ct = std::cos(sd.theta0), st = std::sin(sd.theta0);
    const gp_XYZ col1(a * uPerp.Dot(cf.X) - cf.R * ct, a * uPerp.Dot(cf.Y) - cf.R * st, 0.0);
    const gp_XYZ col2(b * vPerp.Dot(cf.X) + sd.sigma * cf.R * st,
                      b * vPerp.Dot(cf.Y) - sd.sigma * cf.R * ct, 0.0);
    sd.rhoRes = dCentre + ellipseTrigSupZeroOffset(gp_Vec(col1), gp_Vec(col2));
    return std::isfinite(sd.rhoRes);
}

// de Boor on the full (expanded) knot vector; degree n, poles p[0..nPoles-1].
double deBoor(int n, int nPoles, const double* K, const double* p, double t) {
    // span index s with K[s] <= t < K[s+1], clamped to the valid range [n, nPoles-1]
    int s = n;
    while (s < nPoles - 1 && t >= K[s + 1]) ++s;
    std::vector<double> d((size_t)n + 1);
    for (int j = 0; j <= n; ++j) d[(size_t)j] = p[s - n + j];
    for (int r = 1; r <= n; ++r) {
        for (int j = n; j >= r; --j) {
            const int i = s - n + j;
            const double den = K[i + n - r + 1] - K[i];
            const double alpha = (den > 0.0) ? (t - K[i]) / den : 0.0;
            // d[j-1] + alpha (d[j] - d[j-1]): exact when the two coefficients are
            // equal, so a constant v (a circle section) evaluates to its own bits
            d[(size_t)j] = d[(size_t)j - 1] + alpha * (d[(size_t)j] - d[(size_t)j - 1]);
        }
    }
    return d[(size_t)n];
}

}  // namespace

// ---------------------------------------------------------------------------
// classes
// ---------------------------------------------------------------------------

const char* ellipseDevClassName(EllipseDevClass c) {
    switch (c) {
        case EllipseDevClass::OnPlane:         return "ellipse-on-plane";
        case EllipseDevClass::OnCylinder:      return "ellipse-on-cylinder";
        case EllipseDevClass::OnCylinderBound: return "ellipse-on-cylinder-bound";
        case EllipseDevClass::Unhandled:       break;
    }
    return "unhandled-ellipse";
}

bool ellipseDevClassIsExact(EllipseDevClass c) {
    return c == EllipseDevClass::OnPlane || c == EllipseDevClass::OnCylinder;
}

// ---------------------------------------------------------------------------
// the plane class
// ---------------------------------------------------------------------------

double ellipseTrigSupZeroOffset(const gp_Vec& A, const gp_Vec& B) {
    const double aa = A.SquareMagnitude(), bb = B.SquareMagnitude(), ab = A.Dot(B);
    const double diff = aa - bb;
    const double disc = std::sqrt(diff * diff + 4.0 * ab * ab);
    return std::sqrt(0.5 * (aa + bb + disc));
}

double ellipseOnPlaneMax(const gp_Pln& pl, const gp_Elips& el, double f, double l,
                         const gp_Pnt& p0, const gp_Dir& u2, const gp_Dir& v2,
                         double a2, double b2,
                         EllipseDevClass* clsOut, const char** clauseOut) {
    setRefused(clsOut, clauseOut, kEllipseClausePlaneDegenerate);
    const double a = el.MajorRadius(), b = el.MinorRadius();
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(a2) || !std::isfinite(b2) ||
        !std::isfinite(f) || !std::isfinite(l))
        return -1.0;
    if (!(a > gp::Resolution()) || !(b > gp::Resolution())) return -1.0;
    if (!(a2 > gp::Resolution()) || !(b2 > gp::Resolution())) return -1.0;
    if (!(l - f > Precision::PConfusion())) return -1.0;
    const gp_XYZ C = el.Location().XYZ();
    const gp_XYZ U = el.XAxis().Direction().XYZ();
    const gp_XYZ V = el.YAxis().Direction().XYZ();
    const gp_XYZ P0 = p0.XYZ(), U2 = u2.XYZ(), V2 = v2.XYZ();
    const gp_XYZ O = pl.Location().XYZ(), N = pl.Axis().Direction().XYZ();
    if (!finite3(C) || !finite3(U) || !finite3(V) || !finite3(P0) || !finite3(U2) || !finite3(V2) ||
        !finite3(O) || !finite3(N))
        return -1.0;

    // "plane-param": the lifted To2d frame must lie in the plane and be a frame
    setRefused(clsOut, clauseOut, kEllipseClausePlaneParam);
    if (!(std::fabs((P0 - O).Dot(N)) <= Precision::Confusion())) return -1.0;
    if (!(std::fabs(U2.Dot(N)) <= Precision::Angular())) return -1.0;
    if (!(std::fabs(V2.Dot(N)) <= Precision::Angular())) return -1.0;
    if (!(std::fabs(U2.Dot(V2)) <= Precision::Angular())) return -1.0;

    const gp_XYZ W = C - P0;
    const gp_XYZ A = U * a - U2 * a2;
    const gp_XYZ B = V * b - V2 * b2;
    setAdmitted(clsOut, clauseOut, EllipseDevClass::OnPlane);

    // W = 0 on a range covering a half period: |A cos t + B sin t| has period
    // pi, so the supremum is sigma_max([A B]) — the published degenerate branch.
    if (W.X() == 0.0 && W.Y() == 0.0 && W.Z() == 0.0 && l - f >= kPi)
        return ellipseTrigSupZeroOffset(gp_Vec(A), gp_Vec(B));

    return trigSupOnRange(TrigForm(W, A, B), f, l);
}

// ---------------------------------------------------------------------------
// the cylinder section identity
// ---------------------------------------------------------------------------

bool ellipseIsCylinderSection(const gp_Cylinder& cyl, const gp_Elips& el,
                              double& theta0, double& sigma, double& z0,
                              double& lambda, double& psi, double linTol) {
    if (!std::isfinite(linTol) || linTol < 0.0) return false;
    SectionData sd;
    if (!sectionOf(cyl, el, Precision::Angular(), linTol, sd)) return false;
    theta0 = sd.theta0;
    sigma = sd.sigma;
    z0 = sd.z0;
    lambda = sd.lambda;
    psi = sd.psi;
    return true;
}

// ---------------------------------------------------------------------------
// the grazing predicate
// ---------------------------------------------------------------------------

double ellipseCylGrazingSin(const gp_Cylinder& cyl, const gp_Pln& pl) {
    const gp_XYZ d = cyl.Position().Direction().XYZ();
    const gp_XYZ n = pl.Axis().Direction().XYZ();
    if (!finite3(d) || !finite3(n)) return -1.0;
    return std::min(1.0, std::fabs(n.Dot(d)));
}

double ellipseVertexBound(double deltaPlane, double deltaCyl, double sinGamma,
                          double cosGamma, double rhoMinCurv) {
    if (!std::isfinite(deltaPlane) || !std::isfinite(deltaCyl) || !std::isfinite(sinGamma) ||
        !std::isfinite(cosGamma) || !std::isfinite(rhoMinCurv))
        return -1.0;
    if (deltaPlane < 0.0 || deltaCyl < 0.0 || !(rhoMinCurv > 0.0)) return -1.0;
    if (!(sinGamma > 0.0)) return std::numeric_limits<double>::infinity();
    const double h = deltaPlane;
    const double s = (deltaCyl + deltaPlane * std::fabs(cosGamma)) / sinGamma;
    return std::hypot(s, h) + s * s / (2.0 * rhoMinCurv);
}

// ---------------------------------------------------------------------------
// the cylinder class
// ---------------------------------------------------------------------------

double ellipseChebyshevLebesgue(int degree) {
    if (degree < 0) return -1.0;
    return (2.0 / kPi) * std::log((double)degree + 1.0) + 1.0;
}

double ellipseChebyshevRemainder(double lambda, double H, int degree) {
    if (!std::isfinite(lambda) || !std::isfinite(H) || lambda < 0.0 || H < 0.0) return -1.0;
    if (degree < 0 || degree > BSplCLib::MaxDegree()) return -1.0;
    const int n = degree;
    double fact = 1.0;  // (n + 1)!
    for (int k = 2; k <= n + 1; ++k) fact *= (double)k;
    return lambda * std::pow(H, (double)(n + 1)) / (std::ldexp(1.0, 2 * n + 1) * fact);
}

double ellipseOnCylMax(const gp_Cylinder& cyl, const gp_Elips& el, double f, double l,
                       int degree, int nPoles, const double* poleU, const double* poleV,
                       int nKnots, const double* knots, const int* mults,
                       EllipseDevClass* clsOut, double angTol, double linTol,
                       const char** clauseOut) {
    // -degenerate
    setRefused(clsOut, clauseOut, kEllipseClauseCylDegenerate);
    if (!tolsOk(angTol, linTol)) return -1.0;
    if (!std::isfinite(f) || !std::isfinite(l) || !(l - f > Precision::PConfusion())) return -1.0;
    const CylFrame cf(cyl);
    if (!cf.ok() || !(cf.R > gp::Resolution())) return -1.0;
    const double a = el.MajorRadius(), b = el.MinorRadius();
    if (!std::isfinite(a) || !std::isfinite(b) || !(a > gp::Resolution()) || !(b > gp::Resolution()))
        return -1.0;
    if (!finite3(el.Location().XYZ())) return -1.0;
    if (!poleU || !poleV || !knots || !mults) return -1.0;
    for (int i = 0; i < nPoles; ++i)
        if (!std::isfinite(poleU[i]) || !std::isfinite(poleV[i])) return -1.0;
    for (int i = 0; i < nKnots; ++i)
        if (!std::isfinite(knots[i])) return -1.0;

    // -span: a clamped 2-D B-spline (a Bézier is the one-span case) whose
    // parameter range IS [f, l] within PConfusion (the D-140-9 §6 range clause)
    setRefused(clsOut, clauseOut, kEllipseClauseCylSpan);
    if (degree < 1 || degree > BSplCLib::MaxDegree()) return -1.0;
    if (nKnots < 2 || nPoles < degree + 1) return -1.0;
    int sumMults = 0;
    for (int i = 0; i < nKnots; ++i) {
        if (mults[i] < 1) return -1.0;
        const bool end = (i == 0 || i == nKnots - 1);
        if (end && mults[i] != degree + 1) return -1.0;
        if (!end && mults[i] > degree) return -1.0;
        if (i > 0 && !(knots[i] > knots[i - 1])) return -1.0;
        sumMults += mults[i];
    }
    if (sumMults != nPoles + degree + 1) return -1.0;
    if (!(std::fabs(knots[0] - f) <= Precision::PConfusion())) return -1.0;
    if (!(std::fabs(knots[nKnots - 1] - l) <= Precision::PConfusion())) return -1.0;

    // -oncyl: Step 1
    setRefused(clsOut, clauseOut, kEllipseClauseCylOnCyl);
    SectionData sd;
    if (!sectionOf(cyl, el, angTol, linTol, sd)) return -1.0;
    if (!(sd.rhoRes < cf.R)) return -1.0;  // azimuth of the section would be unbounded
    const double eta = std::asin(std::min(1.0, sd.rhoRes / cf.R));  // |theta - affine| <= eta

    // full knot vector and Greville abscissae
    std::vector<double> K;
    K.reserve((size_t)(nPoles + degree + 1));
    for (int i = 0; i < nKnots; ++i)
        for (int m = 0; m < mults[i]; ++m) K.push_back(knots[i]);

    // -uaffine: Step 2, eps_u exact from the poles (convex hull), theta0 on the
    // poles' own 2 pi branch
    setRefused(clsOut, clauseOut, kEllipseClauseCylUAffine);
    const double twoPi = 2.0 * kPi;
    double theta0 = sd.theta0;
    {
        double xi0 = 0.0;
        for (int j = 1; j <= degree; ++j) xi0 += K[j];
        xi0 /= (double)degree;
        const double branch = std::round((poleU[0] - (theta0 + sd.sigma * xi0)) / twoPi);
        theta0 += branch * twoPi;
    }
    double epsU = 0.0;
    for (int i = 0; i < nPoles; ++i) {
        double xi = 0.0;
        for (int j = 1; j <= degree; ++j) xi += K[i + j];
        xi /= (double)degree;
        epsU = std::max(epsU, std::fabs(poleU[i] - (theta0 + sd.sigma * xi)));
    }
    // The derived margin: the azimuthal term 4 R (R + d) sin^2((eps_u + eta)/2)
    // is monotone in eps_u only while eps_u + eta < pi; beyond it the composition
    // is no longer a bound.  Below it eps_u is CARRIED EXACTLY, never dropped —
    // measured on GeomProjLib::Curve2d output (unit case 12b) eps_u ranges from
    // 1e-15 to ~1e-9 depending on the arc, so a fixed angular tolerance here
    // would be a hand constant, and none is used.
    if (!(epsU + eta < kPi)) return -1.0;

    // Step 3: per span, Chebyshev nodes; L_n max|v - z| + Lambda H^(n+1)/(2^(2n+1)(n+1)!)
    const double Ln = ellipseChebyshevLebesgue(degree);
    double dzMax = 0.0;
    for (int sIdx = 0; sIdx + 1 < nKnots; ++sIdx) {
        const double lo = knots[sIdx], hi = knots[sIdx + 1];
        const double H = hi - lo, mid = 0.5 * (lo + hi);
        double rMax = 0.0;
        for (int j = 0; j <= degree; ++j) {
            const double x = mid + 0.5 * H * std::cos((2.0 * (double)j + 1.0) * kPi / (2.0 * ((double)degree + 1.0)));
            const double v = deBoor(degree, nPoles, K.data(), poleV, x);
            const double z = sd.z0 + sd.zc * std::cos(x) + sd.zs * std::sin(x);
            rMax = std::max(rMax, std::fabs(v - z));
        }
        const double rem = ellipseChebyshevRemainder(sd.lambda, H, degree);
        if (rem < 0.0) return -1.0;
        dzMax = std::max(dzMax, Ln * rMax + rem);
    }

    // Step 4: compose
    const double ang = std::min(kPi, epsU + eta);
    const double sh = std::sin(0.5 * ang);
    const double dev2 = 4.0 * cf.R * (cf.R + sd.rhoRes) * sh * sh + sd.rhoRes * sd.rhoRes + dzMax * dzMax;
    if (!std::isfinite(dev2)) {
        setRefused(clsOut, clauseOut, kEllipseClauseCylDegenerate);
        return -1.0;
    }
    const bool exact = (epsU == 0.0 && sd.rhoRes == 0.0 && dzMax == 0.0);
    setAdmitted(clsOut, clauseOut, exact ? EllipseDevClass::OnCylinder : EllipseDevClass::OnCylinderBound);
    return exact ? 0.0 : std::sqrt(dev2);
}

}  // namespace refit
}  // namespace stl2step
