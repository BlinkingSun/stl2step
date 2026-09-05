// ellipse-math — certifies the closed forms in src/refit_ellipse_math.cpp against
// independent ground truth (SPEC-ellipse PART A.4, D-140-9 §2, D-130-2, D-S3-111).
//
// Two oracles, neither of which shares a line of derivation with the code under
// test:
//   (O1) the direct 3-D distance |E(t) - S(pc(t))| at >= 10^4 uniform parameters
//        (kSamples), E from ElCLib, pc from OCCT's own Geom2d curve, S from
//        ElSLib — then the best sample is refined by a 120-step ternary search so
//        the oracle's maximum is the true supremum to ~1e-15, not to the grid;
//   (O2) GeomAPI_ProjectPointOnSurf on Geom_Plane / Geom_CylindricalSurface for
//        the surface-distance cases.
//
// The plane class is asserted EXACT: equal to (O1) and an upper bound at every
// sample (`over == 0`).  The cylinder class is asserted a BOUND: >= (O1) at every
// sample, with the bound/truth ratio reported.  The pcurves are OCCT's own —
// GeomAPI::To2d on the plane (parameter preservation asserted, case 1) and
// GeomProjLib::Curve2d on the cylinder (the engine's projector, case 12b) —
// plus S11-b's measured degree-8 Bézier poles verbatim (case 7).
//
// SPDX-License-Identifier: MIT

#include "refit_ellipse_math.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Plane.hxx>
#include <GeomProjLib.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <gp.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Elips2d.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

using namespace stl2step::refit;

static const double kPi = 3.14159265358979323846264338327950288;
static const int kSamples = 10000;   // SPEC-ellipse A.4 (O1): >= 10^4 uniform parameters
static const double kOverSlack = 1e-12;  // rounding allowance for "upper bound at every sample"

static int gPass = 0;
static int gFail = 0;
static double gMaxAbsErr = 0.0;   // max |closed form - refined sampled supremum| over exact cases
static long gOverSamples = 0;     // samples exceeding the closed form / bound (must stay 0)
static double gBoundRatio = 0.0;  // S11-b bound / measured truth

static void check(bool ok, const char* name) {
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

static void checkNear(double got, double want, double tol, const char* name) {
    const bool ok = std::isfinite(got) && std::isfinite(want) && std::fabs(got - want) <= tol;
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s (%.17g vs %.17g)\n", name, got, want);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s got %.17g want %.17g tol %.3g\n", name, got, want, tol);
    }
}

static void checkLE(double got, double bound, const char* name) {
    const bool ok = std::isfinite(got) && std::isfinite(bound) && got <= bound;
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s (%.17g <= %.17g)\n", name, got, bound);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s got %.17g not <= %.17g\n", name, got, bound);
    }
}

static void checkSameBits(double got, double want, const char* name) {
    const bool ok = std::isfinite(got) && std::isfinite(want) && got == want;
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s (%.17g == %.17g, identical ulp)\n", name, got, want);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s got %.17g want %.17g (differ by %.3g)\n", name, got, want,
                     got - want);
    }
}

static void checkStr(const char* got, const char* want, const char* name) {
    const bool ok = got && want && std::strcmp(got, want) == 0;
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s (\"%s\")\n", name, got);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s got \"%s\" want \"%s\"\n", name, got ? got : "(null)", want);
    }
}

// deterministic LCG (SPEC A.4 case 12)
struct Lcg {
    std::uint64_t s;
    explicit Lcg(std::uint64_t seed) : s(seed) {}
    double next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double)(s >> 11) * (1.0 / 9007199254740992.0);
    }
    double range(double lo, double hi) { return lo + (hi - lo) * next(); }
};

static gp_Ax3 randomFrame(Lcg& g) {
    const gp_Pnt o(g.range(-50.0, 50.0), g.range(-50.0, 50.0), g.range(-50.0, 50.0));
    gp_Dir n(g.range(-1.0, 1.0), g.range(-1.0, 1.0), g.range(-1.0, 1.0));
    gp_Vec x(g.range(-1.0, 1.0), g.range(-1.0, 1.0), g.range(-1.0, 1.0));
    x -= gp_Vec(n) * x.Dot(gp_Vec(n));
    if (x.Magnitude() < 1e-6) x = gp_Vec(n).Crossed(gp_Vec(0, 0, 1));
    return gp_Ax3(o, n, gp_Dir(x));
}

// ---------------------------------------------------------------------------
// oracle (O1): sampled supremum with ternary refinement of the best bracket
// ---------------------------------------------------------------------------
struct Sup {
    double max = -1.0;   // refined supremum
    double gridMax = -1.0;
    double tMax = 0.0;
    int over = 0;        // samples > closed + slack
};

static Sup sampledSup(const std::function<double(double)>& dev, double f, double l, double closed,
                      double slack = kOverSlack) {
    Sup s;
    std::vector<double> d((size_t)kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
        const double t = f + (l - f) * (double)i / (double)kSamples;
        d[(size_t)i] = dev(t);
        if (d[(size_t)i] > s.gridMax) {
            s.gridMax = d[(size_t)i];
            s.tMax = t;
        }
        if (closed >= 0.0 && d[(size_t)i] > closed + slack) ++s.over;
    }
    s.max = s.gridMax;
    // Refine every local maximum of the grid (the endpoints included) — a
    // degree-2 trig polynomial has at most two per period, but a maximum
    // that nearly merges with a minimum can hide a narrow peak between two
    // coarse samples, so each bracket gets a fine grid before the ternary
    // search.  At most the 8 highest local maxima are refined (a constant
    // deviation makes every sample one).
    std::vector<std::pair<double, int>> locals;
    for (int i = 0; i <= kSamples; ++i) {
        const bool leftOk = (i == 0) || d[(size_t)i] >= d[(size_t)i - 1];
        const bool rightOk = (i == kSamples) || d[(size_t)i] >= d[(size_t)i + 1];
        if (leftOk && rightOk) locals.emplace_back(-d[(size_t)i], i);
    }
    std::sort(locals.begin(), locals.end());
    const double h = (l - f) / (double)kSamples;
    for (size_t n = 0; n < locals.size() && n < 8; ++n) {
        const int best = locals[n].second;
        double lo = std::max(f, f + h * (double)(best - 1));
        double hi = std::min(l, f + h * (double)(best + 1));
        const int fine = 2000;
        double fBest = -1.0, tBest = lo;
        for (int j = 0; j <= fine; ++j) {
            const double t = lo + (hi - lo) * (double)j / (double)fine;
            const double v = dev(t);
            if (v > fBest) {
                fBest = v;
                tBest = t;
            }
        }
        const double fh = (hi - lo) / (double)fine;
        lo = std::max(f, tBest - fh);
        hi = std::min(l, tBest + fh);
        for (int it = 0; it < 100; ++it) {
            const double m1 = lo + (hi - lo) / 3.0;
            const double m2 = hi - (hi - lo) / 3.0;
            if (dev(m1) >= dev(m2))
                hi = m2;
            else
                lo = m1;
        }
        const double tm = 0.5 * (lo + hi);
        const double vm = std::max(fBest, dev(tm));
        if (vm > s.max) {
            s.max = vm;
            s.tMax = tm;
        }
    }
    return s;
}

// ---------------------------------------------------------------------------
// plane-side helpers
// ---------------------------------------------------------------------------
struct Lifted {
    bool ok = false;
    gp_Pnt p0;
    gp_Dir u2, v2;
    double a2 = 0.0, b2 = 0.0;
};

// The 2-D ellipse To2d wrote, lifted through the plane's own frame.
static Lifted liftTo2d(const gp_Pln& pl, const Handle(Geom2d_Curve)& c2d) {
    Lifted L;
    Handle(Geom2d_Curve) c = c2d;
    Handle(Geom2d_TrimmedCurve) tr = Handle(Geom2d_TrimmedCurve)::DownCast(c);
    if (!tr.IsNull()) c = tr->BasisCurve();
    Handle(Geom2d_Ellipse) e2 = Handle(Geom2d_Ellipse)::DownCast(c);
    if (e2.IsNull()) return L;
    const gp_Elips2d el2 = e2->Elips2d();
    const gp_Ax3& fr = pl.Position();
    const gp_XYZ X = fr.XDirection().XYZ(), Y = fr.YDirection().XYZ(), O = fr.Location().XYZ();
    const gp_Pnt2d c0 = el2.Location();
    const gp_Dir2d xd = el2.XAxis().Direction();
    const gp_Dir2d yd = el2.YAxis().Direction();
    L.p0 = gp_Pnt(O + X * c0.X() + Y * c0.Y());
    L.u2 = gp_Dir(X * xd.X() + Y * xd.Y());
    L.v2 = gp_Dir(X * yd.X() + Y * yd.Y());
    L.a2 = el2.MajorRadius();
    L.b2 = el2.MinorRadius();
    L.ok = true;
    return L;
}

static double planeDevAt(const gp_Pln& pl, const gp_Elips& el, const Handle(Geom2d_Curve)& c2d,
                         double t) {
    const gp_Pnt2d uv = c2d->Value(t);
    return ElCLib::Value(t, el).Distance(ElSLib::Value(uv.X(), uv.Y(), pl));
}

// Ellipse lying IN the plane: centre (cx, cy) in plane coordinates, major axis
// rotated by phi from X; flipN puts the ellipse normal anti-parallel to N.
static gp_Elips inPlaneEllipse(const gp_Pln& pl, double a, double b, double phi, double cx, double cy,
                               bool flipN = false) {
    const gp_Ax3& fr = pl.Position();
    const gp_XYZ X = fr.XDirection().XYZ(), Y = fr.YDirection().XYZ(), N = fr.Direction().XYZ();
    const gp_Pnt C(fr.Location().XYZ() + X * cx + Y * cy);
    const gp_Dir U(X * std::cos(phi) + Y * std::sin(phi));
    const gp_Dir n(flipN ? N * (-1.0) : N);
    return gp_Elips(gp_Ax2(C, n, U), a, b);
}

// General ellipse near the plane: the in-plane one rotated by `tilt` about its
// own major axis and lifted by `off` along N.
static gp_Elips tiltedEllipse(const gp_Pln& pl, double a, double b, double phi, double cx, double cy,
                              double tilt, double off) {
    gp_Elips e = inPlaneEllipse(pl, a, b, phi, cx, cy);
    e.Rotate(e.XAxis(), tilt);
    e.Translate(gp_Vec(pl.Position().Direction()) * off);
    return e;
}

// refit_build.cpp:7361 verbatim (D-S3-111 circle-on-plane triangle bound) — the
// sibling class the a = b case is compared against; it is not the code under test.
static double circleOnPlaneMaxRef(const gp_Pln& pl, const gp_Elips& circ) {
    Standard_Real A = 0, B = 0, C = 0, D = 0;
    pl.Coefficients(A, B, C, D);
    const gp_Pnt c = circ.Location();
    const double off = std::fabs(A * c.X() + B * c.Y() + C * c.Z() + D);
    double ang = circ.Axis().Direction().Angle(pl.Axis().Direction());
    ang = std::min(ang, kPi - ang);
    const double R = circ.MajorRadius();
    return off + R * std::sin(ang) + R * (1.0 - std::cos(ang));
}

// Runs the plane class on OCCT's To2d image and compares with (O1). Returns the
// closed form (or -1).  `exactCase` feeds gMaxAbsErr.
struct PlaneRun {
    double closed = -1.0;
    EllipseDevClass cls = EllipseDevClass::Unhandled;
    const char* clause = nullptr;
    Sup sup;
    Lifted lift;
};

static PlaneRun runPlane(const gp_Pln& pl, const gp_Elips& el, double f, double l) {
    PlaneRun r;
    Handle(Geom_Ellipse) ge = new Geom_Ellipse(el);
    Handle(Geom2d_Curve) c2d;
    try {
        c2d = GeomAPI::To2d(ge, pl);
    } catch (const Standard_Failure&) {
        return r;
    }
    if (c2d.IsNull()) return r;
    r.lift = liftTo2d(pl, c2d);
    if (!r.lift.ok) return r;
    r.closed = ellipseOnPlaneMax(pl, el, f, l, r.lift.p0, r.lift.u2, r.lift.v2, r.lift.a2, r.lift.b2,
                                 &r.cls, &r.clause);
    r.sup = sampledSup([&](double t) { return planeDevAt(pl, el, c2d, t); }, f, l, r.closed);
    return r;
}

// ---------------------------------------------------------------------------
// cylinder-side helpers
// ---------------------------------------------------------------------------

// The plane∩cylinder section: plane normal at angle alpha from the axis, its
// in-plane-of-cross-section direction e at azimuth thetaE, cutting the axis at
// axial z = zc.  Major axis U = sin(alpha) d - cos(alpha) e (a = R / cos alpha),
// minor axis V = d x e (b = R).  theta(t) = thetaE + pi - t, z = zc + R tan(alpha) cos t.
struct Section {
    gp_Elips el;
    gp_Pln pl;
    double alpha, thetaE, zc;
};

static Section sectionOf(const gp_Cylinder& cyl, double alpha, double thetaE, double zc) {
    const gp_Ax3& pos = cyl.Position();
    const gp_XYZ O = pos.Location().XYZ(), X = pos.XDirection().XYZ(), Y = pos.YDirection().XYZ(),
                 D = pos.Direction().XYZ();
    const double R = cyl.Radius();
    const gp_XYZ e = X * std::cos(thetaE) + Y * std::sin(thetaE);
    const gp_XYZ N = D * std::cos(alpha) + e * std::sin(alpha);
    const gp_XYZ U = D * std::sin(alpha) - e * std::cos(alpha);
    const gp_XYZ V = D.Crossed(e);
    const gp_Pnt C(O + D * zc);
    Section s;
    s.alpha = alpha;
    s.thetaE = thetaE;
    s.zc = zc;
    s.pl = gp_Pln(gp_Ax3(C, gp_Dir(N), gp_Dir(U)));
    s.el = gp_Elips(gp_Ax2(C, gp_Dir(U.Crossed(V)), gp_Dir(U)), R / std::cos(alpha), R);
    return s;
}

static double cylDevAt(const gp_Cylinder& cyl, const gp_Elips& el, const Handle(Geom2d_Curve)& c2d,
                       double t) {
    const gp_Pnt2d uv = c2d->Value(t);
    return ElCLib::Value(t, el).Distance(ElSLib::Value(uv.X(), uv.Y(), cyl));
}

struct Poles2d {
    bool ok = false;
    int degree = 0;
    std::vector<double> u, v, knots;
    std::vector<int> mults;
    double first = 0.0, last = 0.0;
};

static Poles2d extractBSpline(const Handle(Geom2d_Curve)& c2d) {
    Poles2d p;
    Handle(Geom2d_Curve) c = c2d;
    Handle(Geom2d_TrimmedCurve) tr = Handle(Geom2d_TrimmedCurve)::DownCast(c);
    if (!tr.IsNull()) c = tr->BasisCurve();
    Handle(Geom2d_BSplineCurve) bs = Handle(Geom2d_BSplineCurve)::DownCast(c);
    if (bs.IsNull() || bs->IsPeriodic()) return p;
    p.degree = bs->Degree();
    for (int i = 1; i <= bs->NbPoles(); ++i) {
        p.u.push_back(bs->Pole(i).X());
        p.v.push_back(bs->Pole(i).Y());
    }
    for (int i = 1; i <= bs->NbKnots(); ++i) {
        p.knots.push_back(bs->Knot(i));
        p.mults.push_back(bs->Multiplicity(i));
    }
    p.first = bs->FirstParameter();
    p.last = bs->LastParameter();
    p.ok = true;
    return p;
}

static Handle(Geom2d_BSplineCurve) makeBSpline(const Poles2d& p) {
    TColgp_Array1OfPnt2d poles(1, (int)p.u.size());
    for (int i = 0; i < (int)p.u.size(); ++i) poles.SetValue(i + 1, gp_Pnt2d(p.u[(size_t)i], p.v[(size_t)i]));
    TColStd_Array1OfReal knots(1, (int)p.knots.size());
    TColStd_Array1OfInteger mults(1, (int)p.mults.size());
    for (int i = 0; i < (int)p.knots.size(); ++i) {
        knots.SetValue(i + 1, p.knots[(size_t)i]);
        mults.SetValue(i + 1, p.mults[(size_t)i]);
    }
    return new Geom2d_BSplineCurve(poles, knots, mults, p.degree);
}

struct CylRun {
    double bound = -1.0;
    EllipseDevClass cls = EllipseDevClass::Unhandled;
    const char* clause = nullptr;
    Sup sup;
};

static CylRun runCyl(const gp_Cylinder& cyl, const gp_Elips& el, double f, double l, const Poles2d& p,
                     const Handle(Geom2d_Curve)& c2d, double angTol = Precision::Angular(),
                     double linTol = Precision::Confusion()) {
    CylRun r;
    r.bound = ellipseOnCylMax(cyl, el, f, l, p.degree, (int)p.u.size(), p.u.data(), p.v.data(),
                              (int)p.knots.size(), p.knots.data(), p.mults.data(), &r.cls, angTol,
                              linTol, &r.clause);
    if (!c2d.IsNull())
        r.sup = sampledSup([&](double t) { return cylDevAt(cyl, el, c2d, t); }, f, l, r.bound);
    return r;
}

// S11-b ci=4 (rid 6 cylinder | rid 18 corner plane): the shipped pcurve verbatim
// (SPEC-ellipse PART 0).  Cylinder R = 2, 90° arc, degree-8 Bézier.
static Poles2d s11bPoles() {
    Poles2d p;
    p.ok = true;
    p.degree = 8;
    const double U[9] = {1.570796326795, 1.374446785946, 1.178097245096, 0.981747704245,
                         0.785398163397, 0.589048622547, 0.392699081699, 0.196349540849,
                         1.110223024625E-15};
    const double V[9] = {-10.5, -10.1073009183, -9.802723295323, -9.609337384929, -9.542964997599,
                         -9.609337384931, -9.802723295323, -10.1073009183, -10.5};
    p.u.assign(U, U + 9);
    p.v.assign(V, V + 9);
    p.knots = {5.497787143782, 7.068583470577};
    p.mults = {9, 9};
    p.first = p.knots[0];
    p.last = p.knots[1];
    return p;
}

// ---------------------------------------------------------------------------
// 0. class names
// ---------------------------------------------------------------------------
static void testNames() {
    checkStr(ellipseDevClassName(EllipseDevClass::OnPlane), "ellipse-on-plane", "name OnPlane");
    checkStr(ellipseDevClassName(EllipseDevClass::OnCylinder), "ellipse-on-cylinder", "name OnCylinder");
    checkStr(ellipseDevClassName(EllipseDevClass::OnCylinderBound), "ellipse-on-cylinder-bound",
             "name OnCylinderBound");
    checkStr(ellipseDevClassName(EllipseDevClass::Unhandled), "unhandled-ellipse", "name Unhandled");
    check(ellipseDevClassIsExact(EllipseDevClass::OnPlane) && ellipseDevClassIsExact(EllipseDevClass::OnCylinder) &&
              !ellipseDevClassIsExact(EllipseDevClass::OnCylinderBound) &&
              !ellipseDevClassIsExact(EllipseDevClass::Unhandled),
          "ellipseDevClassIsExact: OnPlane, OnCylinder only");
    checkStr(kEllipseClauseGrazing, "unhandled-ellipse-grazing", "grazing clause name");
}

// ---------------------------------------------------------------------------
// 1. To2d parameter preservation — a certificate clause, asserted against OCCT
// ---------------------------------------------------------------------------
static void testTo2dParamPreservation() {
    Lcg g(0x140E11);
    double worst = 0.0;
    int nOk = 0;
    for (int k = 0; k < 13; ++k) {
        const gp_Pln pl(randomFrame(g));
        const double a = g.range(0.5, 20.0), b = a * g.range(0.1, 1.0);
        const gp_Elips el = inPlaneEllipse(pl, a, b, g.range(0.0, 2.0 * kPi), g.range(-30.0, 30.0),
                                           g.range(-30.0, 30.0), (k % 3) == 1);
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(el);
        Handle(Geom2d_Curve) c2d = GeomAPI::To2d(ge, pl);
        Lifted L = liftTo2d(pl, c2d);
        if (!L.ok) {
            check(false, "To2d of an ellipse is a Geom2d_Ellipse");
            continue;
        }
        // the 2-D point at parameter t IS the image of the 3-D point at t
        double w = 0.0;
        for (int i = 0; i <= kSamples; ++i) {
            const double t = 2.0 * kPi * (double)i / (double)kSamples;
            w = std::max(w, planeDevAt(pl, el, c2d, t));
        }
        worst = std::max(worst, w);
        if (L.a2 == a && L.b2 == b) ++nOk;
    }
    checkNear(worst, 0.0, 1e-12, "To2d is parameter-preserving on 13 in-plane placements (incl. flipped normal)");
    check(nOk == 13, "To2d keeps MajorRadius / MinorRadius bitwise on all 13");
}

// ---------------------------------------------------------------------------
// 2. the shipped case: exact plane section -> 0.0 bitwise
// ---------------------------------------------------------------------------
static void testShippedCase() {
    const gp_Pln tilted(gp_Ax3(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3)));
    const gp_Pln located(gp_Ax3(gp_Pnt(12.0, -40.0, 3.0), gp_Dir(0, 0, 1), gp_Dir(0.6, 0.8, 0)));
    const gp_Pln planes[2] = {tilted, located};
    const char* names[2] = {"tilted frame", "located frame"};
    for (int k = 0; k < 2; ++k) {
        const gp_Pln& pl = planes[k];
        // S11-b's radii: a = 2 sqrt(3), b = 2
        const gp_Elips el = inPlaneEllipse(pl, 3.4641016151377544, 2.0, 0.7, 4.0, -2.5);
        // (i) bitwise-consistent inputs: the lifted image IS the ellipse
        EllipseDevClass cls;
        const char* clause = nullptr;
        const double v = ellipseOnPlaneMax(pl, el, 5.497787143782, 7.068583470577, el.Location(),
                                           el.XAxis().Direction(), el.YAxis().Direction(),
                                           el.MajorRadius(), el.MinorRadius(), &cls, &clause);
        char nm[160];
        std::snprintf(nm, sizeof nm, "shipped case (%s): A = B = W = 0 -> 0.0 bitwise", names[k]);
        checkSameBits(v, 0.0, nm);
        std::snprintf(nm, sizeof nm, "shipped case (%s): class ellipse-on-plane", names[k]);
        checkStr(clause, "ellipse-on-plane", nm);
        check(cls == EllipseDevClass::OnPlane && ellipseDevClassIsExact(cls), "shipped case: exact class");
        // full period, same inputs
        const double v2 = ellipseOnPlaneMax(pl, el, 0.0, 2.0 * kPi, el.Location(), el.XAxis().Direction(),
                                            el.YAxis().Direction(), el.MajorRadius(), el.MinorRadius());
        checkSameBits(v2, 0.0, "shipped case, full period: 0.0 bitwise (W = 0 branch)");
        // (ii) the real To2d round trip: the DIAG_PCDEV 1.3e-15 measurement reproduced
        PlaneRun r = runPlane(pl, el, 5.497787143782, 7.068583470577);
        std::snprintf(nm, sizeof nm, "shipped case (%s): To2d round trip <= 1e-13 (measured 1.29e-15)", names[k]);
        checkLE(r.closed, 1e-13, nm);
        check(r.sup.over == 0, "shipped case: upper bound at every sample");
        gOverSamples += r.sup.over;
        gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(r.closed - r.sup.max));
    }
}

// ---------------------------------------------------------------------------
// 3. tilted / offset ellipse: the general quartic path vs (O1)
// ---------------------------------------------------------------------------
static void testGeneralQuartic() {
    Lcg g(0x3A11);
    struct Case { double a, b, phi, cx, cy, tilt, off, f, l; };
    const Case cases[] = {
        {5.0, 3.0, 0.4, 2.0, -1.0, 0.20, 0.30, 0.0, 2.0 * kPi},
        {5.0, 3.0, 0.4, 2.0, -1.0, 0.20, 0.30, 0.7, 2.9},
        {10.0, 2.0, 1.9, -8.0, 3.0, 0.05, -0.02, 5.497787143782, 7.068583470577},
        {1.0, 0.999, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 1.0},
        {7.5, 7.5, 2.2, 1.0, 1.0, 0.35, 0.10, -1.0, 4.0},
        {3.0, 1.0, 0.0, 0.0, 0.0, 1.20, 0.00, 0.0, 2.0 * kPi},
        {3.0, 1.0, 0.0, 0.0, 0.0, 0.00, 0.50, 0.0, 2.0 * kPi},   // pure offset (A = B = 0)
        {20.0, 4.0, 0.9, 12.0, -7.0, 0.002, 0.0001, 1.0, 1.0 + 3e-3},  // tiny range: endpoint sup
    };
    int idx = 0;
    for (const Case& c : cases) {
        const gp_Pln pl(randomFrame(g));
        const gp_Elips el = tiltedEllipse(pl, c.a, c.b, c.phi, c.cx, c.cy, c.tilt, c.off);
        PlaneRun r = runPlane(pl, el, c.f, c.l);
        char nm[200];
        std::snprintf(nm, sizeof nm, "general quartic case %d: closed == (O1) refined (<= 1e-10)", idx);
        checkNear(r.closed, r.sup.max, 1e-10, nm);
        std::snprintf(nm, sizeof nm, "general quartic case %d: upper bound at every sample (over == 0)", idx);
        check(r.sup.over == 0, nm);
        std::snprintf(nm, sizeof nm, "general quartic case %d: class ellipse-on-plane", idx);
        checkStr(r.clause, "ellipse-on-plane", nm);
        gOverSamples += r.sup.over;
        gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(r.closed - r.sup.max));
        ++idx;
    }
}

// ---------------------------------------------------------------------------
// 4. W = 0: sigma_max branch bitwise identical to the general path's dispatch,
//    and the quartic path reproduces sigma_max on a sub-period range
// ---------------------------------------------------------------------------
static void testZeroOffsetBranch() {
    Lcg g(0x77);
    double worstQuartic = 0.0;
    for (int k = 0; k < 8; ++k) {
        const gp_Pln pl(randomFrame(g));
        const double a = g.range(1.0, 15.0), b = a * g.range(0.2, 1.0);
        gp_Elips el = inPlaneEllipse(pl, a, b, g.range(0.0, 2.0 * kPi), g.range(-10.0, 10.0), g.range(-10.0, 10.0));
        el.Rotate(el.XAxis(), g.range(-0.6, 0.6));  // tilted about its own major axis, centre kept
        // lifted image with p0 == C exactly (W = 0): To2d's directions and radii,
        // the centre taken from the ellipse itself
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(el);
        Handle(Geom2d_Curve) c2d = GeomAPI::To2d(ge, pl);
        Lifted L = liftTo2d(pl, c2d);
        // p0 = C is on the plane (centre kept), so the param clause admits it
        const gp_Vec A = gp_Vec(el.XAxis().Direction()) * a - gp_Vec(L.u2) * L.a2;
        const gp_Vec B = gp_Vec(el.YAxis().Direction()) * b - gp_Vec(L.v2) * L.b2;
        const double sig = ellipseTrigSupZeroOffset(A, B);
        EllipseDevClass cls;
        const double v = ellipseOnPlaneMax(pl, el, 0.0, 2.0 * kPi, el.Location(), L.u2, L.v2, L.a2, L.b2, &cls);
        char nm[160];
        std::snprintf(nm, sizeof nm, "W = 0 placement %d: ellipseOnPlaneMax == sigma_max([A B]) bitwise", k);
        checkSameBits(v, sig, nm);
        // the quartic path on [t* - 1, t* + 1] (length 2 < pi, so no dispatch)
        const double tStar = 0.5 * std::atan2(2.0 * A.Dot(B), A.SquareMagnitude() - B.SquareMagnitude());
        const double vq = ellipseOnPlaneMax(pl, el, tStar - 1.0, tStar + 1.0, el.Location(), L.u2, L.v2, L.a2, L.b2);
        worstQuartic = std::max(worstQuartic, std::fabs(vq - sig) / std::max(1.0, sig));
        // and (O1) on the same W = 0 image, lifted back with p0 = C
        const gp_Pnt P0 = el.Location();
        const gp_XYZ U2 = L.u2.XYZ(), V2 = L.v2.XYZ();
        Sup s = sampledSup(
            [&](double t) {
                const gp_XYZ S = P0.XYZ() + U2 * (L.a2 * std::cos(t)) + V2 * (L.b2 * std::sin(t));
                return (ElCLib::Value(t, el).XYZ() - S).Modulus();
            },
            0.0, 2.0 * kPi, v);
        std::snprintf(nm, sizeof nm, "W = 0 placement %d: sigma_max == (O1) refined", k);
        checkNear(v, s.max, 1e-10, nm);
        check(s.over == 0, "W = 0: upper bound at every sample");
        gOverSamples += s.over;
        gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(v - s.max));
    }
    checkNear(worstQuartic, 0.0, 1e-12, "W = 0: the quartic path reproduces sigma_max on a sub-period range (rel)");
    // A = B = 0 exactly
    checkSameBits(ellipseTrigSupZeroOffset(gp_Vec(0, 0, 0), gp_Vec(0, 0, 0)), 0.0, "sigma_max(0, 0) == 0.0");
    checkNear(ellipseTrigSupZeroOffset(gp_Vec(3, 0, 0), gp_Vec(0, 4, 0)), 4.0, 1e-15, "sigma_max of orthogonal (3, 4) = 4");
    checkNear(ellipseTrigSupZeroOffset(gp_Vec(1, 0, 0), gp_Vec(1, 0, 0)), std::sqrt(2.0), 1e-15,
              "sigma_max of parallel unit pair = sqrt 2");
}

// ---------------------------------------------------------------------------
// 5. a = b: the sanity relation with circleOnPlaneMax (D-S3-111 triangle bound)
// ---------------------------------------------------------------------------
static void testCircleRelation() {
    Lcg g(0xC1C);
    const double tilts[] = {0.0, 0.05, 0.35, 0.8, 1.3};
    const double offs[] = {0.0, 0.01, 0.7};
    int idx = 0;
    for (double tilt : tilts) {
        for (double off : offs) {
            const gp_Pln pl(randomFrame(g));
            const double R = g.range(0.5, 12.0);
            const gp_Elips circ = tiltedEllipse(pl, R, R, g.range(0.0, 2.0 * kPi), g.range(-5.0, 5.0),
                                                g.range(-5.0, 5.0), tilt, off);
            PlaneRun r = runPlane(pl, circ, 0.0, 2.0 * kPi);
            const double tri = circleOnPlaneMaxRef(pl, circ);
            char nm[200];
            std::snprintf(nm, sizeof nm, "a = b case %d (tilt %.2f off %.2f): ellipseOnPlaneMax <= circleOnPlaneMax", idx, tilt, off);
            checkLE(r.closed, tri + 1e-12, nm);
            std::snprintf(nm, sizeof nm, "a = b case %d: ellipseOnPlaneMax >= (O1) sampled max", idx);
            checkLE(r.sup.gridMax, r.closed + kOverSlack, nm);
            std::snprintf(nm, sizeof nm, "a = b case %d: exact (== (O1) refined)", idx);
            checkNear(r.closed, r.sup.max, 1e-10, nm);
            gOverSamples += r.sup.over;
            gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(r.closed - r.sup.max));
            ++idx;
        }
    }
}

// ---------------------------------------------------------------------------
// 6. ellipseIsCylinderSection on corpus-scale geometries, and its refusals
// ---------------------------------------------------------------------------
static void testSectionIdentity() {
    struct Geo { const char* name; double R, alpha, thetaE, zc; gp_Ax3 frame; };
    const gp_Ax3 zFrame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Ax3 tiltedFrame(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3));
    const std::vector<Geo> geos = {
        {"S11-b fillet R=2, corner plane alpha=acos(1/sqrt3)", 2.0, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5, zFrame},
        {"small fillet R=0.5, alpha=30deg, tilted frame", 0.5, kPi / 6.0, 0.4, 3.0, tiltedFrame},
        {"bore R=5, alpha=60deg", 5.0, kPi / 3.0, 2.0, -7.0, zFrame},
        {"large R=10, alpha=80deg (a = 57.6), tilted frame", 10.0, 80.0 * kPi / 180.0, 5.5, 0.25, tiltedFrame},
    };
    for (const Geo& gg : geos) {
        const gp_Cylinder cyl(gg.frame, gg.R);
        const Section s = sectionOf(cyl, gg.alpha, gg.thetaE, gg.zc);
        double th0 = 0, sg = 0, z0 = 0, lam = 0, psi = 0;
        const bool ok = ellipseIsCylinderSection(cyl, s.el, th0, sg, z0, lam, psi);
        char nm[220];
        std::snprintf(nm, sizeof nm, "section identity: %s", gg.name);
        check(ok, nm);
        if (!ok) continue;
        checkNear(sg, -1.0, 0.0, "  sigma = -1 (theta = thetaE + pi - t)");
        checkNear(z0, gg.zc, 1e-12, "  z0 = zc");
        checkNear(lam, gg.R * std::tan(gg.alpha), 1e-9, "  Lambda = R tan(alpha)");
        double dth = std::fmod(th0 - (gg.thetaE + kPi), 2.0 * kPi);
        if (dth < -kPi) dth += 2.0 * kPi;
        if (dth > kPi) dth -= 2.0 * kPi;
        checkNear(dth, 0.0, 1e-12, "  theta0 = thetaE + pi (mod 2 pi)");
        checkNear(std::fmod(std::fabs(psi), 2.0 * kPi), 0.0, 1e-9, "  psi = 0 (z = zc + Lambda cos t)");
        // every point of the ellipse is ON the cylinder: (O2) GeomAPI_ProjectPointOnSurf
        Handle(Geom_CylindricalSurface) gs = new Geom_CylindricalSurface(cyl);
        double worst = 0.0;
        for (int i = 0; i < 360; ++i) {
            const double t = 2.0 * kPi * (double)i / 360.0;
            GeomAPI_ProjectPointOnSurf proj(ElCLib::Value(t, s.el), gs);
            if (proj.NbPoints() > 0) worst = std::max(worst, proj.LowerDistance());
            // and the identity's own claims
            const gp_XYZ w = ElCLib::Value(t, s.el).XYZ() - cyl.Location().XYZ();
            const double z = w.Dot(cyl.Position().Direction().XYZ());
            const gp_XYZ wp = w - cyl.Position().Direction().XYZ() * z;
            worst = std::max(worst, std::fabs(wp.Modulus() - gg.R));
            worst = std::max(worst, std::fabs(z - (z0 + lam * std::cos(t - psi))));
        }
        checkNear(worst, 0.0, 1e-9, "  rho == R, z == z0 + Lambda cos(t - psi) at 360 points; (O2) on-surface");
        // grazing sin for this section = cos(alpha) = b/a
        checkNear(ellipseCylGrazingSin(cyl, s.pl), std::cos(gg.alpha), 1e-12, "  ellipseCylGrazingSin = cos(alpha)");
        checkNear(ellipseCylGrazingSin(cyl, s.pl), s.el.MinorRadius() / s.el.MajorRadius(), 1e-12,
                  "  ellipseCylGrazingSin = b / a");
    }
    // refusals
    {
        const gp_Cylinder cyl(zFrame, 2.0);
        const Section s = sectionOf(cyl, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5);
        double th0, sg, z0, lam, psi;
        gp_Elips tilted = s.el;
        tilted.Rotate(tilted.XAxis(), 0.01);
        check(!ellipseIsCylinderSection(cyl, tilted, th0, sg, z0, lam, psi), "section identity refuses a tilted ellipse (0.01 rad)");
        gp_Elips off = s.el;
        off.Translate(gp_Vec(0.5, 0.0, 0.0));
        check(!ellipseIsCylinderSection(cyl, off, th0, sg, z0, lam, psi), "section identity refuses an off-axis centre (0.5 mm)");
        gp_Elips wrongB = s.el;
        wrongB.SetMinorRadius(2.001);
        check(!ellipseIsCylinderSection(cyl, wrongB, th0, sg, z0, lam, psi), "section identity refuses wrong b (2.001 on R = 2)");
        gp_Elips wrongA = s.el;
        wrongA.SetMajorRadius(3.5);
        check(!ellipseIsCylinderSection(cyl, wrongA, th0, sg, z0, lam, psi), "section identity refuses wrong a (3.5 vs 2 sqrt 3)");
        gp_Elips axial = s.el;
        axial.Translate(gp_Vec(0.0, 0.0, 3.0));
        check(ellipseIsCylinderSection(cyl, axial, th0, sg, z0, lam, psi) && std::fabs(z0 - (-9.5)) < 1e-12,
              "section identity admits an axial translation (z0 follows)");
    }
}

// ---------------------------------------------------------------------------
// 7. S11-b's own numbers — the acceptance case for Step 3's constant
// ---------------------------------------------------------------------------
static void testS11b() {
    const gp_Ax3 zFrame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Cylinder cyl(zFrame, 2.0);
    const Section s = sectionOf(cyl, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5);
    checkNear(s.el.MajorRadius(), 3.464101615138, 1e-11, "S11-b: a = 3.464101615138 (the shipped 2-D ELLIPSE)");
    checkNear(s.el.MinorRadius(), 2.0, 0.0, "S11-b: b = 2");
    const Poles2d p = s11bPoles();
    Handle(Geom2d_BSplineCurve) bs = makeBSpline(p);
    const double f = p.knots[0], l = p.knots[1];
    // the measured shape: u affine with du = (pi/2)/8
    double epsU = 0.0;
    for (int i = 0; i < 9; ++i)
        epsU = std::max(epsU, std::fabs(p.u[(size_t)i] - (0.5 * kPi - (double)i * (0.5 * kPi / 8.0))));
    checkLE(epsU, 2e-12, "S11-b: u-poles in exact arithmetic progression, eps_u <= 2e-12 (12-digit STEP print)");
    // The poles here are the STEP file's 12-significant-digit print, so eps_u is
    // 1.9e-12 (the print's own rounding, above Precision::Angular()); the class
    // carries it exactly into the azimuthal term (4 R^2 sin^2(eps_u/2) ~ 1.5e-23)
    // rather than gating on an angular tolerance.
    std::fprintf(stderr, "INFO S11-b 12-digit print: eps_u = %.4e (carried exactly, not gated)\n", epsU);
    CylRun r = runCyl(cyl, s.el, f, l, p, bs);
    checkStr(r.clause, "ellipse-on-cylinder-bound", "S11-b: class ellipse-on-cylinder-bound");
    check(r.cls == EllipseDevClass::OnCylinderBound && !ellipseDevClassIsExact(r.cls), "S11-b: a bound, not exact");
    checkLE(4.05358e-10, r.bound, "S11-b: bound >= measured DIAG_PCDEV 4.05358e-10");
    checkLE(r.bound, 1e-8, "S11-b: bound <= 1e-8");
    checkLE(r.sup.max, r.bound, "S11-b: bound >= (O1) refined supremum");
    check(r.sup.over == 0, "S11-b: bound >= every sample");
    gOverSamples += r.sup.over;
    gBoundRatio = r.bound / r.sup.max;
    std::fprintf(stderr, "INFO S11-b bound=%.6e truth(O1)=%.6e ratio=%.3f remainder(Lambda=%.6f,H=%.6f,n=8)=%.6e L8=%.4f\n",
                 r.bound, r.sup.max, gBoundRatio, 2.0 * std::sqrt(2.0), l - f,
                 ellipseChebyshevRemainder(2.0 * std::sqrt(2.0), l - f, 8), ellipseChebyshevLebesgue(8));
    checkNear(r.sup.max, 4.05358e-10, 5e-11, "S11-b: (O1) reproduces the measured 4.05358e-10 (poles printed to 12 digits)");
    // the midpoint-Taylor remainder the spec REFUSES, for the record: 2^n looser
    const double taylor = 2.0 * std::sqrt(2.0) * std::pow(0.5 * (l - f), 9.0) / 362880.0;
    const double cheb = ellipseChebyshevRemainder(2.0 * std::sqrt(2.0), l - f, 8);
    checkNear(taylor / cheb, 256.0, 1e-9, "S11-b: midpoint-Taylor / Chebyshev remainder = 2^8 (refused form is 256x looser)");
    // the same edge under a tilted, located cylinder frame: bound invariant
    const gp_Ax3 tiltedFrame(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3));
    const gp_Cylinder cyl2(tiltedFrame, 2.0);
    const Section s2 = sectionOf(cyl2, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5);
    CylRun r2 = runCyl(cyl2, s2.el, f, l, p, bs);
    checkNear(r2.bound, r.bound, 1e-12, "S11-b under a tilted/located frame: same bound");
    check(r2.sup.over == 0, "S11-b tilted frame: bound >= every sample");
    gOverSamples += r2.sup.over;
}

// ---------------------------------------------------------------------------
// 8. eps_u exactness: a perturbed u-pole
// ---------------------------------------------------------------------------
static void testEpsUExactness() {
    const gp_Ax3 zFrame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Cylinder cyl(zFrame, 2.0);
    const Section s = sectionOf(cyl, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5);
    Poles2d p = s11bPoles();
    const double delta = 0.01;
    p.u[3] += delta;  // one interior pole
    Handle(Geom2d_BSplineCurve) bs = makeBSpline(p);
    const double f = p.knots[0], l = p.knots[1];
    // admitted at the default: eps_u is carried exactly, |dTheta(t)| <= eps_u at
    // every sample, and the bound holds
    CylRun r = runCyl(cyl, s.el, f, l, p, bs);
    checkStr(r.clause, "ellipse-on-cylinder-bound", "perturbed u-pole (0.01 rad): admitted as a bound, eps_u carried");
    double worstDth = 0.0;
    for (int i = 0; i <= kSamples; ++i) {
        const double t = f + (l - f) * (double)i / (double)kSamples;
        const double u = bs->Value(t).X();
        const double affine = 2.25 * kPi - t;  // theta0 + sigma t for S11-b
        worstDth = std::max(worstDth, std::fabs(u - affine));
    }
    checkLE(worstDth, delta, "perturbed u-pole: |dTheta(t)| <= eps_u = delta at every sample (convex hull)");
    checkLE(r.sup.max, r.bound, "perturbed u-pole: bound >= (O1)");
    check(r.sup.over == 0, "perturbed u-pole: bound >= every sample");
    checkLE(2.0 * 2.0 * std::sin(0.5 * (delta - 1e-11)), r.bound, "perturbed u-pole: bound >= 2 R sin(eps_u / 2) (eps_u = delta - print rounding)");
    gOverSamples += r.sup.over;
    // the derived margin: one pole 3.2 rad off makes eps_u >= pi -> -uaffine by name
    {
        Poles2d bad = s11bPoles();
        bad.u[5] += 3.2;
        Handle(Geom2d_BSplineCurve) bs2 = makeBSpline(bad);
        CylRun rb = runCyl(cyl, s.el, f, l, bad, bs2);
        checkStr(rb.clause, "unhandled-ellipse-cyl-uaffine", "eps_u >= pi (one pole 3.2 rad off): refused -uaffine (the derived margin)");
        check(rb.bound == -1.0, "eps_u >= pi: -1.0");
    }
    // a pcurve running the WRONG way round (u poles reversed): eps_u = pi/2 at the
    // ends, admitted, and the bound 2 R sin(eps_u/2) is still >= the truth
    {
        Poles2d rev = s11bPoles();
        std::reverse(rev.u.begin(), rev.u.end());
        Handle(Geom2d_BSplineCurve) bs2 = makeBSpline(rev);
        CylRun rr = runCyl(cyl, s.el, f, l, rev, bs2);
        checkStr(rr.clause, "ellipse-on-cylinder-bound", "reversed u-poles: admitted as a (large) bound");
        checkLE(rr.sup.max, rr.bound, "reversed u-poles: bound >= (O1) refined (truth ~ 2 R sin(pi/4))");
        check(rr.sup.over == 0, "reversed u-poles: bound >= every sample");
        checkNear(rr.bound, 2.0 * 2.0 * std::sin(0.25 * kPi), 1e-8, "reversed u-poles: bound = 2 R sin(pi/4) + axial (tight)");
        gOverSamples += rr.sup.over;
    }
}

// ---------------------------------------------------------------------------
// 9. Chebyshev remainder: monotone in n, exact 2^(n+1) under halving H
// ---------------------------------------------------------------------------
static void testChebyshevRemainder() {
    const double lam = 2.0 * std::sqrt(2.0), H = 0.5 * kPi;
    bool mono = true;
    for (int n = 1; n < 20; ++n)
        if (!(ellipseChebyshevRemainder(lam, H, n + 1) < ellipseChebyshevRemainder(lam, H, n))) mono = false;
    check(mono, "Chebyshev remainder strictly decreasing in n (1..20) at H = pi/2");
    bool halving = true;
    for (int n = 1; n <= 12; ++n) {
        const double ratio = ellipseChebyshevRemainder(lam, H, n) / ellipseChebyshevRemainder(lam, 0.5 * H, n);
        if (std::fabs(ratio - std::ldexp(1.0, n + 1)) > 1e-9 * std::ldexp(1.0, n + 1)) halving = false;
    }
    check(halving, "halving the span divides the remainder by exactly 2^(n+1) (n = 1..12)");
    checkNear(ellipseChebyshevRemainder(1.0, 4.0, 0), 2.0, 1e-15, "n = 0: Lambda H / 2 (= 2 (H/4)^1)");
    checkNear(ellipseChebyshevRemainder(1.0, 1.0, 1), 1.0 / 16.0, 1e-15, "n = 1: H^2 / (2^3 2!) = 1/16");
    check(ellipseChebyshevRemainder(1.0, 1.0, -1) == -1.0 && ellipseChebyshevRemainder(1.0, 1.0, 26) == -1.0,
          "remainder refuses n < 0 and n > BSplCLib::MaxDegree()");
    checkNear(ellipseChebyshevLebesgue(8), (2.0 / kPi) * std::log(9.0) + 1.0, 1e-15, "L_8 bound = (2/pi) ln 9 + 1");
    // the bound's monotonicity in the class itself: the S11-b pcurve split into two
    // Bézier spans by knot insertion has a smaller remainder term (exact 2^9)
}

// ---------------------------------------------------------------------------
// 10. grazing sweep — the open predicate
// ---------------------------------------------------------------------------
static void testGrazing() {
    const double q = 1.321449896e-05;  // the plate's measured q (D-140-8 U1)
    const double tau = 2.0 * q;
    const double R = 2.0;
    const double cap = 0.000433546;  // S11-b's census meshTolCap
    auto boundAt = [&](double sinG, double dPl, double dCy) {
        const double cosG = std::sqrt(std::max(0.0, 1.0 - sinG * sinG));
        return ellipseVertexBound(dPl, dCy, sinG, cosG, R * sinG);  // rho_min = R cos(alpha) = R sin(gamma)
    };
    const int N = 200;
    std::vector<double> sinG(N + 1), bnd(N + 1);
    bool finiteAll = true, mono = true, smooth = true;
    int firstRefused = -1;
    for (int k = 0; k <= N; ++k) {
        sinG[(size_t)k] = std::pow(10.0, -6.0 * (double)k / (double)N);  // 1 -> 1e-6
        bnd[(size_t)k] = boundAt(sinG[(size_t)k], tau, tau);
        if (!std::isfinite(bnd[(size_t)k]) || bnd[(size_t)k] < 0.0) finiteAll = false;
        if (k > 0 && bnd[(size_t)k] < bnd[(size_t)k - 1]) mono = false;
        if (bnd[(size_t)k] > cap && firstRefused < 0) firstRefused = k;
    }
    check(finiteAll, "grazing sweep: bound finite and >= 0 for sin gamma in [1e-6, 1]");
    check(mono, "grazing sweep: bound non-increasing in sin gamma (refusal set is one interval at the grazing end)");
    check(firstRefused > 0 && firstRefused < N, "grazing sweep: refusal begins strictly inside the sweep");
    for (int k = 1; k < firstRefused; ++k) {
        const double rel = std::fabs(bnd[(size_t)k] / bnd[(size_t)k - 1] - 1.0);
        if (rel > 0.5) smooth = false;  // sweep step 10^(6/200) = 7.2 %; cos(gamma) = sqrt(1 - s^2) rises to 0.36 in the first step (+25 %)
    }
    check(smooth, "grazing sweep: no discontinuity on the admitted side (relative step <= 50 %)");
    // the derived margin sin gamma*: bound(sin gamma*) = cap, by bisection in the test
    auto crossing = [&](double dPl, double dCy, double c) {
        double lo = 1e-9, hi = 1.0;
        for (int it = 0; it < 200; ++it) {
            const double m = 0.5 * (lo + hi);
            if (boundAt(m, dPl, dCy) > c) lo = m; else hi = m;
        }
        return 0.5 * (lo + hi);
    };
    const double sStar = crossing(tau, tau, cap);
    std::fprintf(stderr, "INFO grazing: q=%.6e tau=%.6e R=%g cap=%.6e -> sin gamma* = %.6e (gamma* = %.4f deg)\n",
                 q, tau, R, cap, sStar, std::asin(sStar) * 180.0 / kPi);
    check(bnd[(size_t)firstRefused] > cap && bnd[(size_t)firstRefused - 1] <= cap &&
              sinG[(size_t)firstRefused] < sStar && sStar <= sinG[(size_t)firstRefused - 1],
          "grazing sweep: the class refuses exactly where ellipseVertexBound > meshTolCap");
    check(crossing(2.0 * tau, 2.0 * tau, cap) > sStar, "margin derives from q: doubling tau raises sin gamma*");
    check(crossing(tau, tau, 2.0 * cap) < sStar, "margin derives from meshTolCap: doubling the cap lowers sin gamma*");
    checkNear(crossing(2.0 * tau, 2.0 * tau, cap) / sStar, 2.0, 0.05, "small-angle regime: sin gamma* ~ 2 tau / cap (doubles with tau)");
    // the bound's own algebra
    checkNear(ellipseVertexBound(0.0, 0.0, 0.5, std::sqrt(0.75), 1.0), 0.0, 0.0, "vertex bound: zero deviations -> 0");
    checkNear(ellipseVertexBound(1e-3, 0.0, 1.0, 0.0, 1.0), 1e-3, 1e-18, "vertex bound: gamma = 90deg, only h -> deltaPlane");
    {
        const double s = (2e-3 + 1e-3 * 0.6) / 0.8;
        checkNear(ellipseVertexBound(1e-3, 2e-3, 0.8, 0.6, 5.0), std::hypot(s, 1e-3) + s * s / 10.0, 1e-18,
                  "vertex bound: sqrt(s^2 + h^2) + s^2 / (2 rho_min) verbatim");
    }
    check(std::isinf(ellipseVertexBound(1e-3, 1e-3, 0.0, 1.0, 1.0)), "vertex bound: sin gamma = 0 -> +inf (tangential cut)");
    check(ellipseVertexBound(std::nan(""), 1e-3, 0.5, 0.5, 1.0) == -1.0 && ellipseVertexBound(1e-3, 1e-3, 0.5, 0.5, 0.0) == -1.0,
          "vertex bound: NaN / rho_min = 0 -> -1.0");
    // ellipseCylGrazingSin extremes
    const gp_Ax3 zFrame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Cylinder cyl(zFrame, 2.0);
    checkSameBits(ellipseCylGrazingSin(cyl, gp_Pln(gp_Pnt(0, 0, 5), gp_Dir(0, 0, 1))), 1.0, "grazing sin: plane perp axis (a circle) = 1");
    checkSameBits(ellipseCylGrazingSin(cyl, gp_Pln(gp_Pnt(2, 0, 0), gp_Dir(1, 0, 0))), 0.0, "grazing sin: plane tangent to the cylinder = 0");
    checkSameBits(ellipseCylGrazingSin(cyl, gp_Pln(gp_Pnt(0, 0, 5), gp_Dir(0, 0, -1))), 1.0, "grazing sin: anti-parallel normal = 1");
}

// ---------------------------------------------------------------------------
// 11. degenerate inputs, each refused by NAME
// ---------------------------------------------------------------------------
static void testDegenerate() {
    const gp_Pln pl(gp_Ax3(gp_Pnt(1, 2, 3), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)));
    const gp_Elips el = inPlaneEllipse(pl, 3.0, 2.0, 0.3, 1.0, 1.0);
    const gp_Pnt p0 = el.Location();
    const gp_Dir u2 = el.XAxis().Direction(), v2 = el.YAxis().Direction();
    const double nan = std::nan("");
    EllipseDevClass cls;
    const char* clause = nullptr;
    auto plane = [&](const gp_Elips& e, double f, double l, const gp_Pnt& q, const gp_Dir& a, const gp_Dir& b,
                     double a2, double b2) { return ellipseOnPlaneMax(pl, e, f, l, q, a, b, a2, b2, &cls, &clause); };
    check(plane(el, nan, 1.0, p0, u2, v2, 3.0, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
          "plane: non-finite f -> -1.0, -plane-degenerate");
    check(plane(el, 0.0, 1.0, gp_Pnt(nan, 0, 0), u2, v2, 3.0, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
          "plane: non-finite p0 -> -plane-degenerate");
    check(plane(el, 0.0, 1.0, p0, u2, v2, nan, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
          "plane: non-finite a2 -> -plane-degenerate");
    {
        gp_Elips tiny(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.0, 1.0);
        tiny.SetMinorRadius(0.0);
        check(plane(tiny, 0.0, 1.0, p0, u2, v2, 1.0, 1.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
              "plane: b = 0 (<= gp::Resolution()) -> -plane-degenerate");
        tiny.SetMajorRadius(0.0);
        check(plane(tiny, 0.0, 1.0, p0, u2, v2, 1.0, 1.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
              "plane: a = 0 -> -plane-degenerate");
    }
    check(plane(el, 1.0, 1.0 + 0.5 * Precision::PConfusion(), p0, u2, v2, 3.0, 2.0) == -1.0 &&
              std::strcmp(clause, kEllipseClausePlaneDegenerate) == 0,
          "plane: l - f <= PConfusion -> -plane-degenerate");
    check(plane(el, 0.0, 1.0, gp_Pnt(1, 2, 3.5), u2, v2, 3.0, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneParam) == 0,
          "plane: lifted centre off the plane -> -plane-param");
    check(plane(el, 0.0, 1.0, p0, gp_Dir(0, 0, 1), v2, 3.0, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneParam) == 0,
          "plane: lifted u2 out of the plane -> -plane-param");
    check(plane(el, 0.0, 1.0, p0, u2, u2, 3.0, 2.0) == -1.0 && std::strcmp(clause, kEllipseClausePlaneParam) == 0,
          "plane: u2 not orthogonal to v2 -> -plane-param");
    check(plane(el, 0.0, 1.0, p0, u2, v2, 3.0, 2.0) >= 0.0 && std::strcmp(clause, "ellipse-on-plane") == 0,
          "plane: the same inputs, well-formed -> admitted");

    // cylinder
    const gp_Ax3 zFrame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    const gp_Cylinder cyl(zFrame, 2.0);
    const Section s = sectionOf(cyl, std::acos(1.0 / std::sqrt(3.0)), 1.25 * kPi, -12.5);
    const Poles2d p = s11bPoles();
    const double f = p.knots[0], l = p.knots[1];
    auto cy = [&](const gp_Elips& e, double ff, double ll, const Poles2d& pp) {
        return ellipseOnCylMax(cyl, e, ff, ll, pp.degree, (int)pp.u.size(), pp.u.data(), pp.v.data(),
                               (int)pp.knots.size(), pp.knots.data(), pp.mults.data(), &cls,
                               Precision::Angular(), Precision::Confusion(), &clause);
    };
    check(cy(s.el, nan, l, p) == -1.0 && std::strcmp(clause, kEllipseClauseCylDegenerate) == 0, "cyl: non-finite f -> -cyl-degenerate");
    check(cy(s.el, f, f + 0.5 * Precision::PConfusion(), p) == -1.0 && std::strcmp(clause, kEllipseClauseCylDegenerate) == 0,
          "cyl: l - f <= PConfusion -> -cyl-degenerate");
    {
        Poles2d bad = p;
        bad.v[2] = nan;
        check(cy(s.el, f, l, bad) == -1.0 && std::strcmp(clause, kEllipseClauseCylDegenerate) == 0, "cyl: non-finite pole -> -cyl-degenerate");
        gp_Elips e0 = s.el;
        e0.SetMinorRadius(0.0);
        check(cy(e0, f, l, p) == -1.0 && std::strcmp(clause, kEllipseClauseCylDegenerate) == 0, "cyl: b = 0 -> -cyl-degenerate");
        const gp_Cylinder cyl0(zFrame, 0.0);
        check(ellipseOnCylMax(cyl0, s.el, f, l, p.degree, 9, p.u.data(), p.v.data(), 2, p.knots.data(), p.mults.data(), &cls,
                              Precision::Angular(), Precision::Confusion(), &clause) == -1.0 &&
                  std::strcmp(clause, kEllipseClauseCylDegenerate) == 0,
              "cyl: R = 0 -> -cyl-degenerate");
    }
    {
        Poles2d bad = p;
        bad.u.pop_back();
        bad.v.pop_back();  // 8 poles for degree 8
        check(cy(s.el, f, l, bad) == -1.0 && std::strcmp(clause, kEllipseClauseCylSpan) == 0, "cyl: nPoles != degree + 1 -> -cyl-span");
        Poles2d shifted = p;
        shifted.knots[0] += 1e-3;
        shifted.knots[1] += 1e-3;
        check(cy(s.el, f, l, shifted) == -1.0 && std::strcmp(clause, kEllipseClauseCylSpan) == 0, "cyl: knot range != [f, l] by 1e-3 -> -cyl-span (D-140-9 §6 range clause)");
        check(cy(s.el, f, l + 2.0 * Precision::PConfusion(), p) == -1.0 && std::strcmp(clause, kEllipseClauseCylSpan) == 0,
              "cyl: range off by 2 PConfusion -> -cyl-span");
        check(cy(s.el, f, l + 0.5 * Precision::PConfusion(), p) >= 0.0, "cyl: range off by PConfusion/2 -> admitted");
        Poles2d unclamped = p;
        unclamped.mults = {8, 8};
        unclamped.u.pop_back();
        unclamped.v.pop_back();
        unclamped.u.pop_back();
        unclamped.v.pop_back();  // 7 poles: sum mults 16 = 7 + 8 + 1 but ends not degree+1
        check(cy(s.el, f, l, unclamped) == -1.0 && std::strcmp(clause, kEllipseClauseCylSpan) == 0, "cyl: unclamped ends -> -cyl-span");
        Poles2d deg0 = p;
        deg0.degree = 0;
        check(cy(s.el, f, l, deg0) == -1.0 && std::strcmp(clause, kEllipseClauseCylSpan) == 0, "cyl: degree 0 -> -cyl-span");
    }
    {
        gp_Elips off = s.el;
        off.Translate(gp_Vec(0.5, 0, 0));
        check(cy(off, f, l, p) == -1.0 && std::strcmp(clause, kEllipseClauseCylOnCyl) == 0, "cyl: ellipse not a section -> -cyl-oncyl");
    }
    check(cy(s.el, f, l, p) >= 0.0 && std::strcmp(clause, "ellipse-on-cylinder-bound") == 0, "cyl: S11-b well-formed -> admitted");
    // an exact circle section with an exact pcurve -> ellipse-on-cylinder, 0.0
    {
        // alpha = 0 (a circle, a = b = R) in the cylinder's own frame: theta0 = 0,
        // sigma = +1, z = 4, u(t) = t, v = 4 — every residual is exactly 0.0 in double
        const gp_Elips c0(gp_Ax2(gp_Pnt(0, 0, 4), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), 2.0, 2.0);
        Poles2d lin;
        lin.ok = true;
        lin.degree = 1;
        const double ff = 0.2, ll = 1.7;
        lin.u = {ff, ll};
        lin.v = {4.0, 4.0};
        lin.knots = {ff, ll};
        lin.mults = {2, 2};
        const double v0 = ellipseOnCylMax(cyl, c0, ff, ll, 1, 2, lin.u.data(), lin.v.data(), 2, lin.knots.data(),
                                          lin.mults.data(), &cls, Precision::Angular(), Precision::Confusion(), &clause);
        checkSameBits(v0, 0.0, "cyl: circle section (alpha = 0) with its exact degree-1 pcurve -> 0.0 bitwise");
        checkStr(clause, "ellipse-on-cylinder", "cyl: ... class ellipse-on-cylinder (exact)");
        check(cls == EllipseDevClass::OnCylinder && ellipseDevClassIsExact(cls), "cyl: ... IsExact");
    }
}

// ---------------------------------------------------------------------------
// 12. random deterministic placements
// ---------------------------------------------------------------------------
static void testRandomPlane() {
    Lcg g(0x140140);
    int admitted = 0, over = 0;
    double worst = 0.0;
    const int N = 2000;
    for (int k = 0; k < N; ++k) {
        const gp_Pln pl(randomFrame(g));
        const double a = g.range(0.2, 30.0), b = a * g.range(0.05, 1.0);
        const int kind = k % 4;  // 0 in-plane, 1 tilted, 2 offset, 3 both
        const double tilt = (kind == 1 || kind == 3) ? g.range(-1.4, 1.4) : 0.0;
        const double off = (kind == 2 || kind == 3) ? g.range(-2.0, 2.0) : 0.0;
        const gp_Elips el = tiltedEllipse(pl, a, b, g.range(0.0, 2.0 * kPi), g.range(-40.0, 40.0), g.range(-40.0, 40.0), tilt, off);
        double f = g.range(-2.0 * kPi, 2.0 * kPi);
        double l = f + ((k % 5 == 0) ? 2.0 * kPi : g.range(0.05, 2.0 * kPi));
        PlaneRun r = runPlane(pl, el, f, l);
        if (r.closed < 0.0) continue;
        ++admitted;
        over += r.sup.over;
        const double err = std::fabs(r.closed - r.sup.max);
        if (err > 1e-10)
            std::fprintf(stderr, "DEBUG random plane k=%d a=%.6g b=%.6g tilt=%.6g off=%.6g f=%.9g l=%.9g closed=%.12g grid=%.12g refined=%.12g tMax=%.9g over=%d\n",
                         k, a, b, tilt, off, f, l, r.closed, r.sup.gridMax, r.sup.max, r.sup.tMax, r.sup.over);
        worst = std::max(worst, err);
    }
    std::fprintf(stderr, "INFO random plane: admitted=%d/%d over=%d maxAbsErr=%.3e\n", admitted, N, over, worst);
    check(admitted == N, "random plane: every one of 2000 placements admitted as ellipse-on-plane");
    check(over == 0, "random plane: upper bound at every sample of every placement (over == 0)");
    checkLE(worst, 1e-10, "random plane: |closed - (O1) refined| <= 1e-10 on all 2000");
    gOverSamples += over;
    gMaxAbsErr = std::max(gMaxAbsErr, worst);
}

static void testRandomCylinder() {
    Lcg g(0xCC140);
    int admitted = 0, over = 0, projected = 0, notBSpline = 0;
    double worstRatio = 0.0, bestRatio = std::numeric_limits<double>::max(), maxBound = 0.0;
    std::map<std::string, int> shapes, clauses;
    std::vector<double> epsUs;
    const int N = 300;
    for (int k = 0; k < N; ++k) {
        const gp_Ax3 fr = randomFrame(g);
        const double R = g.range(0.5, 20.0);
        const gp_Cylinder cyl(fr, R);
        const double alpha = g.range(0.0, 75.0 * kPi / 180.0);
        const double L = g.range(0.3, 3.0);
        const double f = g.range(-3.0, 3.0), l = f + L;
        // keep u = thetaE + pi - t inside (0.05, 2 pi - 0.05): no seam crossing
        const double thetaE = g.range(l - kPi + 0.05, kPi - 0.05 + f);
        const Section s = sectionOf(cyl, alpha, thetaE, g.range(-20.0, 20.0));
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(s.el);
        Handle(Geom_CylindricalSurface) gs = new Geom_CylindricalSurface(cyl);
        Handle(Geom2d_Curve) c2d;
        try {
            c2d = GeomProjLib::Curve2d(ge, f, l, gs);
        } catch (const Standard_Failure&) {
        }
        if (c2d.IsNull()) continue;
        ++projected;
        Poles2d p = extractBSpline(c2d);
        if (!p.ok) {
            ++notBSpline;
            continue;
        }
        char sh[64];
        std::snprintf(sh, sizeof sh, "deg%d/spans%d", p.degree, (int)p.knots.size() - 1);
        shapes[sh]++;
        CylRun r = runCyl(cyl, s.el, f, l, p, c2d);
        clauses[r.clause ? r.clause : "(null)"]++;
        // test-side eps_u: poles vs the affine map at the Greville abscissae
        {
            std::vector<double> K;
            for (size_t i = 0; i < p.knots.size(); ++i)
                for (int m = 0; m < p.mults[i]; ++m) K.push_back(p.knots[i]);
            double e = 0.0;
            for (size_t i = 0; i < p.u.size(); ++i) {
                double xi = 0.0;
                for (int j = 1; j <= p.degree; ++j) xi += K[i + (size_t)j];
                xi /= (double)p.degree;
                double aff = thetaE + kPi - xi;
                double du = p.u[i] - aff;
                du -= 2.0 * kPi * std::round(du / (2.0 * kPi));
                e = std::max(e, std::fabs(du));
            }
            epsUs.push_back(e);
        }
        if (r.bound < 0.0) continue;
        ++admitted;
        over += r.sup.over;
        maxBound = std::max(maxBound, r.bound);
        if (r.sup.max > 0.0) {
            worstRatio = std::max(worstRatio, r.bound / r.sup.max);
            bestRatio = std::min(bestRatio, r.bound / r.sup.max);
        }
    }
    std::fprintf(stderr, "INFO random cylinder (GeomProjLib::Curve2d): projected=%d notBSpline=%d admitted=%d over=%d maxBound=%.3e ratio[min,max]=[%.2f, %.2f]\n",
                 projected, notBSpline, admitted, over, maxBound, bestRatio, worstRatio);
    for (const auto& kv : shapes) std::fprintf(stderr, "INFO   pcurve shape %s: %d\n", kv.first.c_str(), kv.second);
    for (const auto& kv : clauses) std::fprintf(stderr, "INFO   class %s: %d\n", kv.first.c_str(), kv.second);
    if (!epsUs.empty()) {
        std::sort(epsUs.begin(), epsUs.end());
        int above1e12 = 0;
        for (double e : epsUs) if (e > 1e-12) ++above1e12;
        std::fprintf(stderr, "INFO   eps_u of the projector's u-poles: min=%.2e median=%.2e max=%.2e; > Precision::Angular(): %d/%d\n",
                     epsUs.front(), epsUs[epsUs.size() / 2], epsUs.back(), above1e12, (int)epsUs.size());
    }
    check(projected >= N * 9 / 10, "random cylinder: the projector produced a pcurve on >= 90 %");
    check(admitted == projected - notBSpline, "random cylinder: the class admits every one of the projector's own pcurves");
    check(over == 0, "random cylinder: the bound is >= every sample of every admitted pcurve (over == 0)");
    gOverSamples += over;
}

int main() {
    testNames();
    testTo2dParamPreservation();
    testShippedCase();
    testGeneralQuartic();
    testZeroOffsetBranch();
    testCircleRelation();
    testSectionIdentity();
    testS11b();
    testEpsUExactness();
    testChebyshevRemainder();
    testGrazing();
    testDegenerate();
    testRandomPlane();
    testRandomCylinder();
    std::fprintf(stderr, "\nellipse_math_unit: %d passed, %d failed\n", gPass, gFail);
    std::fprintf(stderr, "ELLIPSE-MATH-UNIT exactClass=ellipse-on-plane boundClass=ellipse-on-cylinder-bound unit=%d/%d maxAbsErr=%.3e boundRatio=%.3f overSamples=%ld\n",
                 gPass, gPass + gFail, gMaxAbsErr, gBoundRatio, gOverSamples);
    return gFail == 0 ? 0 : 1;
}
