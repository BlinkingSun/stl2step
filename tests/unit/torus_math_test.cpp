// torus-math — certifies the closed-form torus math in src/refit_torus_math.cpp
// against independent ground truth (SPEC-torus PART B, D-140-6 §3(4), D-130-2).
//
// Two independent oracles, neither of which shares a line of derivation with
// the code under test:
//   * numPointTorusDist — minimises, over the azimuth u, the 3-D distance from
//     the point to the MERIDIAN CIRCLE at u, built straight from the OCCT
//     parametric equation S(u,v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z.
//     It never uses the (rho, z) reduction the implementation is built on.
//   * GeomAPI_ProjectPointOnSurf on a Geom_ToroidalSurface — OCCT's own answer.
//
// The supremum is checked against >= 10^4 sampled circle points per circle
// (SPEC-torus B.2), asserting (i) the closed form is an UPPER bound at every
// sample, (ii) it equals the sampled maximum, and (iii) every point of a coaxial
// circle is equidistant (max - min over the samples ~ 0), which is the proof.
//
// Numbers under test are the fixtures' (SPEC-torus A.2 / A.5): S04 (Rmaj 7,
// Rmin 3, centre (25,25,22), R_cyl 10, tangency rho=10@z=22 and rho=7@z=25),
// S19 torus A (12, 2, centre z=18) and torus B (11, 1, centre z=1), R_cyl 10.
//
// SPDX-License-Identifier: MIT

#include "refit_torus_math.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Torus.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

using stl2step::refit::TorusDevClass;
using stl2step::refit::circleOnTorusMax;
using stl2step::refit::pointTorusDist;
using stl2step::refit::torusAxialCoord;
using stl2step::refit::torusDevClassIsExact;
using stl2step::refit::torusDevClassName;
using stl2step::refit::torusIsRing;
using stl2step::refit::torusRadialCoord;
using stl2step::refit::torusVIsoCircle;
using stl2step::refit::torusVOfProfilePoint;

static const double kPi = 3.14159265358979323846264338327950288;
static const int kSamples = 10000;  // SPEC-torus B.2: >= 10^4 sampled circle points

static int gPass = 0;
static int gFail = 0;
static double gMaxAbsErr = 0.0;  // max |closed form - sampled supremum| over every exact case

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

// Bitwise equality of two doubles (the "identical-ulp" assertion).
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

// ---------------------------------------------------------------------------
// oracle 1: distance to the torus from the OCCT parametric equation directly
// ---------------------------------------------------------------------------

// Distance from p to the MERIDIAN CIRCLE at azimuth u: centre O + R e_r(u),
// radius r, normal e_theta(u). Every point of that circle is S(u, v) for some
// v, so this is a distance to a genuine surface point. Point-to-circle in 3-D:
// split w = p - C into its normal component and its in-plane component.
static double meridianCircleDist(const gp_Torus& t, const gp_Pnt& p, double u) {
    const gp_Ax3& ax = t.Position();
    const gp_XYZ O = ax.Location().XYZ();
    const gp_XYZ X = ax.XDirection().XYZ();
    const gp_XYZ Y = ax.YDirection().XYZ();
    const gp_XYZ er = X * std::cos(u) + Y * std::sin(u);
    const gp_XYZ n = X * (-std::sin(u)) + Y * std::cos(u);  // unit normal of the meridian plane
    const gp_XYZ C = O + er * t.MajorRadius();
    const gp_XYZ w = p.XYZ() - C;
    const double wn = w.Dot(n);
    const gp_XYZ win = w - n * wn;
    return std::hypot(win.Modulus() - t.MinorRadius(), wn);
}

static double numPointTorusDist(const gp_Torus& t, const gp_Pnt& p, int nU = 256) {
    double bestU = 0.0;
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < nU; ++i) {
        const double u = 2.0 * kPi * (double)i / (double)nU;
        const double d = meridianCircleDist(t, p, u);
        if (d < best) {
            best = d;
            bestU = u;
        }
    }
    const double h = 2.0 * kPi / (double)nU;
    double lo = bestU - h;
    double hi = bestU + h;
    // 120 iterations, not 60: at a circle that lies ON the surface the distance
    // is V-shaped in u (hypot(0, x) = |x|), and a ternary search converges only
    // linearly, so 60 halvings-by-thirds left a 3e-12 floor. 120 reach the
    // double-precision floor of the bracket, ~1e-15 * r_c.
    for (int it = 0; it < 120; ++it) {
        const double m1 = lo + (hi - lo) / 3.0;
        const double m2 = hi - (hi - lo) / 3.0;
        if (meridianCircleDist(t, p, m1) <= meridianCircleDist(t, p, m2))
            hi = m2;
        else
            lo = m1;
    }
    return std::min(best, meridianCircleDist(t, p, 0.5 * (lo + hi)));
}

// The numeric supremum over >= 10^4 sampled circle points, with the statistics
// the assertions need: the sampled max and min, and how many samples exceed the
// closed form by more than `slack` (must be 0: the closed form is an upper
// bound at every sample).
struct SampleStats {
    double max = -1.0;
    double min = std::numeric_limits<double>::max();
    int over = 0;
};

static SampleStats numCircleTorusSup(const gp_Torus& t, const gp_Circ& q, double closed,
                                     double slack, int nT = kSamples) {
    SampleStats s;
    for (int i = 0; i < nT; ++i) {
        const double tt = 2.0 * kPi * (double)i / (double)nT;
        const double d = numPointTorusDist(t, ElCLib::Value(tt, q));
        s.max = std::max(s.max, d);
        s.min = std::min(s.min, d);
        if (d > closed + slack) ++s.over;
    }
    return s;
}

// ---------------------------------------------------------------------------
// fixtures (SPEC-torus A.2 / A.5), in their own canonical +Z frames
// ---------------------------------------------------------------------------

struct Fixture {
    const char* name;
    double Rmaj, Rmin, Rcyl;
    gp_Pnt centre;    // profile-circle centre, on the axis
    double planeSide; // +1: the plane-side tangency is at h = +Rmin; -1: at h = -Rmin
};

static std::vector<Fixture> fixtures() {
    return {{"S04 boss-top round", 7.0, 3.0, 10.0, gp_Pnt(25.0, 25.0, 22.0), +1.0},
            {"S19 torus A (top rim r=2)", 12.0, 2.0, 10.0, gp_Pnt(30.0, 30.0, 18.0), +1.0},
            {"S19 torus B (bottom rim r=1)", 11.0, 1.0, 10.0, gp_Pnt(30.0, 30.0, 1.0), -1.0}};
}

static gp_Torus canonicalTorus(const Fixture& f) {
    return gp_Torus(gp_Ax3(f.centre, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), f.Rmaj, f.Rmin);
}

// A coaxial circle of radius rc at axial offset h in the torus's own frame.
static gp_Circ coaxial(const gp_Torus& t, double rc, double h) {
    const gp_Dir za = t.Axis().Direction();
    const gp_Pnt c = t.Location().Translated(gp_Vec(za) * h);
    return gp_Circ(gp_Ax2(c, za, t.Position().XDirection()), rc);
}

// A deliberately off-origin, obliquely-oriented placement: nothing here may
// depend on the torus sitting at the world origin along Z.
static gp_Ax3 tiltedPlacement() {
    return gp_Ax3(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3));
}

// ---------------------------------------------------------------------------
// 1. parametrisation, frame helpers, pointTorusDist
// ---------------------------------------------------------------------------

static void testParametrisationAndPointDist() {
    const gp_Torus t(tiltedPlacement(), 7.0, 3.0);
    check(torusIsRing(t), "S04 torus is a ring");

    // torusVIsoCircle(v) at parameter u IS S(u, v) of OCCT's own parametrisation.
    double worst = 0.0;
    for (int i = 0; i < 13; ++i) {
        const double v = -kPi + 2.0 * kPi * (double)i / 13.0;
        gp_Circ c;
        if (!torusVIsoCircle(t, v, c)) {
            check(false, "torusVIsoCircle on a ring");
            return;
        }
        for (int j = 0; j < 17; ++j) {
            const double u = 2.0 * kPi * (double)j / 17.0;
            worst = std::max(worst, ElCLib::Value(u, c).Distance(ElSLib::Value(u, v, t)));
        }
        // and its (radius, offset) round-trips through torusVOfProfilePoint
        const double rc = c.Radius();
        const double h = torusAxialCoord(t, c.Location());
        double dv = std::fabs(torusVOfProfilePoint(t, rc, h) - v);
        dv = std::min(dv, std::fabs(dv - 2.0 * kPi));
        if (dv > 1e-12) check(false, "torusVOfProfilePoint round-trips torusVIsoCircle");
        // every point of a v-iso circle is ON the surface
        for (int j = 0; j < 17; ++j) {
            const double u = 2.0 * kPi * (double)j / 17.0;
            worst = std::max(worst, std::fabs(pointTorusDist(t, ElCLib::Value(u, c))));
        }
    }
    checkNear(worst, 0.0, 1e-12, "torusVIsoCircle == ElSLib::Value(u, v, torus), and on-surface");

    // frame helpers against the parametrisation
    {
        const gp_Pnt s = ElSLib::Value(1.1, 0.4, t);
        checkNear(torusAxialCoord(t, s), 3.0 * std::sin(0.4), 1e-12, "torusAxialCoord = r sin v");
        checkNear(torusRadialCoord(t, s), 7.0 + 3.0 * std::cos(0.4), 1e-12,
                  "torusRadialCoord = R + r cos v");
    }

    // pointTorusDist against OCCT's own projector, at deterministic points that
    // sweep inside the tube, outside it, near the axis and far away.
    Handle(Geom_ToroidalSurface) gs = new Geom_ToroidalSurface(t.Position(), 7.0, 3.0);
    double worstOcc = 0.0;
    double worstNum = 0.0;
    int nProj = 0;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 7; ++j)
            for (int k = 0; k < 5; ++k) {
                const double rho = 0.25 + 1.75 * (double)i;   // 0.25 .. 14.25
                const double phi = 2.0 * kPi * (double)j / 7.0;
                const double z = -4.0 + 2.0 * (double)k;      // -4 .. 4
                const gp_Ax3& ax = t.Position();
                const gp_Pnt p(ax.Location().XYZ() + ax.XDirection().XYZ() * (rho * std::cos(phi)) +
                               ax.YDirection().XYZ() * (rho * std::sin(phi)) +
                               ax.Direction().XYZ() * z);
                const double mine = pointTorusDist(t, p);
                worstNum = std::max(worstNum, std::fabs(mine - numPointTorusDist(t, p)));
                try {
                    GeomAPI_ProjectPointOnSurf proj(p, gs);
                    if (proj.NbPoints() > 0) {
                        worstOcc = std::max(worstOcc, std::fabs(mine - proj.LowerDistance()));
                        ++nProj;
                    }
                } catch (const Standard_Failure&) {
                }
            }
    checkNear(worstNum, 0.0, 1e-10, "pointTorusDist == meridian-circle oracle (315 points)");
    checkNear(worstOcc, 0.0, 1e-7, "pointTorusDist == GeomAPI_ProjectPointOnSurf");
    check(nProj >= 300, "GeomAPI_ProjectPointOnSurf answered on >= 300 of 315 points");
}

// ---------------------------------------------------------------------------
// 2. the two tangency circles of a mouth round: EXACT zeros (SPEC-torus B.2(a))
// ---------------------------------------------------------------------------

static void testTangencyCirclesExactZero() {
    char nm[256];
    for (const Fixture& f : fixtures()) {
        const gp_Torus t = canonicalTorus(f);
        TorusDevClass cls = TorusDevClass::Unhandled;

        // cylinder-side tangency: radius R_cyl (= Rmaj +- Rmin) at h = 0
        const gp_Circ cyl = coaxial(t, f.Rcyl, 0.0);
        const double dCyl = circleOnTorusMax(t, cyl, &cls);
        std::snprintf(nm, sizeof nm, "%s: cylinder-side circle rho=%g h=0 -> exactly 0", f.name,
                      f.Rcyl);
        check(dCyl == 0.0, nm);
        std::snprintf(nm, sizeof nm, "%s: cylinder-side class is circle-on-torus", f.name);
        check(cls == TorusDevClass::Coaxial && torusDevClassIsExact(cls) &&
                  std::string(torusDevClassName(cls)) == "circle-on-torus",
              nm);

        // plane-side tangency: radius Rmaj at h = +-Rmin
        const gp_Circ pln = coaxial(t, f.Rmaj, f.planeSide * f.Rmin);
        const double dPln = circleOnTorusMax(t, pln, &cls);
        std::snprintf(nm, sizeof nm, "%s: plane-side circle rho=%g h=%+g -> exactly 0", f.name,
                      f.Rmaj, f.planeSide * f.Rmin);
        check(dPln == 0.0, nm);
        std::snprintf(nm, sizeof nm, "%s: plane-side class is circle-on-torus", f.name);
        check(cls == TorusDevClass::Coaxial, nm);

        // the numeric oracle agrees: >= 10^4 samples, every one within 1e-12 of 0
        for (int side = 0; side < 2; ++side) {
            const gp_Circ& q = side ? pln : cyl;
            const SampleStats s = numCircleTorusSup(t, q, 0.0, 1e-12);
            gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(s.max - 0.0));
            std::snprintf(nm, sizeof nm, "%s: %s tangency circle, oracle max over %d samples", f.name,
                          side ? "plane-side" : "cylinder-side", kSamples);
            checkNear(s.max, 0.0, 1e-12, nm);
            std::snprintf(nm, sizeof nm, "%s: %s tangency circle, closed form >= every sample", f.name,
                          side ? "plane-side" : "cylinder-side");
            check(s.over == 0, nm);
        }

        // the four cardinal v-iso circles (v = 0, PI/2, PI, 3PI/2) are the same
        // two pairs seen through torusVIsoCircle; cos/sin of PI/2 are not exact
        // in binary so these are zero to the arithmetic, not bitwise.
        double worstIso = 0.0;
        for (int k = 0; k < 4; ++k) {
            gp_Circ c;
            if (!torusVIsoCircle(t, 0.5 * kPi * (double)k, c)) {
                check(false, "torusVIsoCircle");
                continue;
            }
            worstIso = std::max(worstIso, circleOnTorusMax(t, c, &cls));
        }
        std::snprintf(nm, sizeof nm, "%s: the four cardinal v-iso circles deviate by 0", f.name);
        checkNear(worstIso, 0.0, 1e-14, nm);
    }

    // The same tangency circles under a tilted placement: zero to
    // Precision::Confusion() (SPEC-torus B.2(a)), and in fact to 1e-13.
    for (const Fixture& f : fixtures()) {
        const gp_Torus t(tiltedPlacement(), f.Rmaj, f.Rmin);
        TorusDevClass cls = TorusDevClass::Unhandled;
        const double a = circleOnTorusMax(t, coaxial(t, f.Rcyl, 0.0), &cls);
        const bool ca = cls == TorusDevClass::Coaxial;
        const double b = circleOnTorusMax(t, coaxial(t, f.Rmaj, f.planeSide * f.Rmin), &cls);
        const bool cb = cls == TorusDevClass::Coaxial;
        std::snprintf(nm, sizeof nm, "%s tilted: both tangency circles <= Precision::Confusion()",
                      f.name);
        check(a >= 0.0 && b >= 0.0 && a <= Precision::Confusion() && b <= Precision::Confusion(), nm);
        std::snprintf(nm, sizeof nm, "%s tilted: both tangency circles zero to 1e-13", f.name);
        checkNear(std::max(a, b), 0.0, 1e-13, nm);
        std::snprintf(nm, sizeof nm, "%s tilted: both class circle-on-torus", f.name);
        check(ca && cb, nm);
    }
}

// ---------------------------------------------------------------------------
// 3. general coaxial circles: identical-ulp closed form, oracle, equidistance
//    (SPEC-torus B.2(b))
// ---------------------------------------------------------------------------

struct RcH { double rc, h; const char* what; };

static const RcH kGeneralCases[] = {
    {10.0, 1.5, "outside the tube, above the equator"},
    {8.5, 0.0, "equatorial plane (h = 0), inside the tube"},
    {7.0, 4.0, "on the core cylinder rho = Rmaj, past the top"},
    {12.0, -0.7, "outside, below the equator"},
    {7.5, 0.5, "inside the tube, near the core circle"},
    {0.0, 2.25, "r_c = 0: a point on the axis"},
    {7.0, 0.0, "r_c = Rmaj, h = 0: the core circle itself (distance Rmin)"},
    {4.1, -2.9, "inside the hole, below"},
};

static void testGeneralCoaxial() {
    char nm[256];
    const double R = 7.0, r = 3.0;

    // (i) canonical frame at the origin: the value is BITWISE the closed form
    //     |hypot(r_c - R, h) - r| written independently here.
    {
        const gp_Torus t(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), R, r);
        for (const RcH& c : kGeneralCases) {
            TorusDevClass cls = TorusDevClass::Unhandled;
            const double got = circleOnTorusMax(t, coaxial(t, c.rc, c.h), &cls);
            const double want = std::fabs(std::hypot(c.rc - R, c.h) - r);
            std::snprintf(nm, sizeof nm, "origin frame (%g, %g) %s: identical ulp", c.rc, c.h, c.what);
            checkSameBits(got, want, nm);
            std::snprintf(nm, sizeof nm, "origin frame (%g, %g): class circle-on-torus", c.rc, c.h);
            check(cls == TorusDevClass::Coaxial, nm);
        }
    }

    // (ii) S04's placed frame and a tilted placement: against the closed form
    //      to 1e-12, against the >= 10^4-sample oracle, an upper bound at every
    //      sample, and every point equidistant.
    const gp_Ax3 frames[] = {gp_Ax3(gp_Pnt(25.0, 25.0, 22.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)),
                             tiltedPlacement()};
    const char* frameNames[] = {"S04 frame", "tilted frame"};
    for (int fi = 0; fi < 2; ++fi) {
        const gp_Torus t(frames[fi], R, r);
        for (const RcH& c : kGeneralCases) {
            TorusDevClass cls = TorusDevClass::Unhandled;
            const gp_Circ q = coaxial(t, c.rc, c.h);
            const double got = circleOnTorusMax(t, q, &cls);
            const double want = std::fabs(std::hypot(c.rc - R, c.h) - r);
            std::snprintf(nm, sizeof nm, "%s (%g, %g): closed form", frameNames[fi], c.rc, c.h);
            checkNear(got, want, 1e-12, nm);

            const SampleStats s = numCircleTorusSup(t, q, got, 1e-12);
            gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(got - s.max));
            std::snprintf(nm, sizeof nm, "%s (%g, %g): == oracle max over %d samples", frameNames[fi],
                          c.rc, c.h, kSamples);
            checkNear(got, s.max, 1e-10, nm);
            std::snprintf(nm, sizeof nm, "%s (%g, %g): upper bound at every sample", frameNames[fi],
                          c.rc, c.h);
            check(s.over == 0, nm);
            std::snprintf(nm, sizeof nm, "%s (%g, %g): every point equidistant (max-min)",
                          frameNames[fi], c.rc, c.h);
            checkNear(s.max - s.min, 0.0, 1e-10, nm);
        }
    }

    // (iii) the sign convention: the modulus drops inside/outside; the two
    //       circles at equal |hypot - r| on either side of the tube agree.
    {
        const gp_Torus t(frames[0], R, r);
        const double inside = circleOnTorusMax(t, coaxial(t, R + 1.0, 0.0));   // hypot 1 -> r - 1 = 2
        const double outside = circleOnTorusMax(t, coaxial(t, R + 5.0, 0.0));  // hypot 5 -> 5 - r = 2
        checkNear(inside, outside, 0.0, "inside (r - g) and outside (g - r) both read a distance");
        checkNear(inside, 2.0, 0.0, "the distance is 2 on both");
    }

    // (iv) the two S19 tori, placed, one non-zero case each, against the oracle.
    for (const Fixture& f : fixtures()) {
        const gp_Torus t = canonicalTorus(f);
        const gp_Circ q = coaxial(t, f.Rcyl, 0.5 * f.Rmin);  // half-way up the cylinder side
        TorusDevClass cls = TorusDevClass::Unhandled;
        const double got = circleOnTorusMax(t, q, &cls);
        const double want = std::fabs(std::hypot(f.Rcyl - f.Rmaj, 0.5 * f.Rmin) - f.Rmin);
        std::snprintf(nm, sizeof nm, "%s: (R_cyl, Rmin/2) closed form", f.name);
        checkNear(got, want, 1e-12, nm);
        const SampleStats s = numCircleTorusSup(t, q, got, 1e-12);
        gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(got - s.max));
        std::snprintf(nm, sizeof nm, "%s: (R_cyl, Rmin/2) == oracle, upper bound everywhere", f.name);
        check(std::fabs(got - s.max) <= 1e-10 && s.over == 0, nm);
    }
}

// ---------------------------------------------------------------------------
// 4. Rmaj -> 0: the sphere limit (SPEC-torus B.2(c)) and the spindle refusal
// ---------------------------------------------------------------------------

static void testSphereLimitAndSpindle() {
    char nm[256];
    const double r = 3.0;
    const gp_Ax3 fr(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));

    // Admitted while Rmaj > Rmin: exact and an upper bound at every sample, down
    // to a near-horn torus Rmaj - Rmin = 1e-6.
    const double admitted[] = {9.0, 6.0, 3.5, 3.000001};
    for (double R : admitted) {
        const gp_Torus t(fr, R, r);
        std::snprintf(nm, sizeof nm, "Rmaj=%g Rmin=3: is a ring", R);
        check(torusIsRing(t), nm);
        gp_Circ iso;
        torusVIsoCircle(t, 2.3, iso);  // a v-iso circle: on the surface
        TorusDevClass cls = TorusDevClass::Unhandled;
        const double d0 = circleOnTorusMax(t, iso, &cls);
        std::snprintf(nm, sizeof nm, "Rmaj=%g: v-iso circle deviates by 0", R);
        checkNear(d0, 0.0, 1e-14, nm);
        const gp_Circ q = coaxial(t, R + 1.0, 1.0);
        const double got = circleOnTorusMax(t, q, &cls);
        const SampleStats s = numCircleTorusSup(t, q, got, 1e-12);
        gMaxAbsErr = std::max(gMaxAbsErr, std::fabs(got - s.max));
        std::snprintf(nm, sizeof nm, "Rmaj=%g: (R+1, 1) == oracle over %d samples, upper bound", R,
                      kSamples);
        check(cls == TorusDevClass::Coaxial && std::fabs(got - s.max) <= 1e-10 && s.over == 0, nm);
    }

    // Refused at and below the horn: Rmaj - Rmin <= gp::Resolution() is
    // gp_Torus::SetMajorRadius's own guard; Rmaj -> 0 is a sphere (D-140-6 §1(4)).
    const double refused[] = {3.0, 2.0, 1.0, 0.5, 0.0};
    for (double R : refused) {
        const gp_Torus t(fr, R, r);
        TorusDevClass cls = TorusDevClass::Coaxial;
        const double d = circleOnTorusMax(t, coaxial(t, R + 1.0, 1.0), &cls);
        std::snprintf(nm, sizeof nm, "Rmaj=%g Rmin=3 (%s): refused -> -1 / Unhandled", R,
                      R == 3.0 ? "horn" : (R == 0.0 ? "sphere" : "spindle"));
        check(d == -1.0 && cls == TorusDevClass::Unhandled && !torusIsRing(t), nm);
        std::snprintf(nm, sizeof nm, "Rmaj=%g: pointTorusDist refuses too", R);
        check(pointTorusDist(t, gp_Pnt(1, 0, 0)) == -1.0, nm);
        gp_Circ c;
        std::snprintf(nm, sizeof nm, "Rmaj=%g: torusVIsoCircle refuses too", R);
        check(!torusVIsoCircle(t, 0.0, c) && std::isnan(torusVOfProfilePoint(t, 1.0, 0.0)), nm);
    }

    // WHY the spindle is refused rather than evaluated: on Rmaj=1, Rmin=3 the
    // circle rho = Rmin - Rmaj = 2 at h = 0 is the spindle's MIRRORED inner
    // equator — it lies ON the untrimmed surface (the oracle reads 0) — while
    // the coaxial closed form would read |hypot(2 - 1, 0) - 3| = 2. The form is
    // an upper bound there but not the supremum, and D-130-2 records only exact
    // classes. Measured, not asserted from the derivation.
    {
        const gp_Torus spindle(fr, 1.0, 3.0);
        const gp_Circ q = coaxial(spindle, 2.0, 0.0);
        const double wouldBe = std::fabs(std::hypot(2.0 - 1.0, 0.0) - 3.0);
        const SampleStats s = numCircleTorusSup(spindle, q, wouldBe, 1e-12, 2000);
        std::fprintf(stderr, "spindle witness: coaxial form %.17g, oracle max %.3g min %.3g\n",
                     wouldBe, s.max, s.min);
        checkNear(s.max, 0.0, 1e-10, "spindle: the mirrored equator lies on the surface (oracle 0)");
        check(wouldBe > 1.0, "spindle: the coaxial form would read 2 there -- not the supremum");
    }
}

// ---------------------------------------------------------------------------
// 5. refusals: non-coaxial circles, degenerate input, bad tolerances
// ---------------------------------------------------------------------------

static void testRefusals() {
    const gp_Torus t(tiltedPlacement(), 7.0, 3.0);
    const gp_Ax3& ax = t.Position();
    const gp_Dir za = ax.Direction();
    const gp_Dir xa = ax.XDirection();
    TorusDevClass cls = TorusDevClass::Coaxial;

    // tilted normal (5 degrees): a general plane section -> refused
    {
        gp_Dir n = za;
        n.Rotate(gp_Ax1(ax.Location(), xa), 5.0 * kPi / 180.0);
        const gp_Circ q(gp_Ax2(ax.Location(), n, xa), 10.0);
        const double d = circleOnTorusMax(t, q, &cls);
        check(d == -1.0 && cls == TorusDevClass::Unhandled &&
                  std::string(torusDevClassName(cls)) == "unhandled-circle-on-torus",
              "tilted circle (5 deg) -> refused, unhandled-circle-on-torus");
    }
    // perpendicular plane but centre off the axis -> refused
    {
        const gp_Pnt c = ax.Location().Translated(gp_Vec(xa) * 0.5);
        const gp_Circ q(gp_Ax2(c, za, xa), 10.0);
        const double d = circleOnTorusMax(t, q, &cls);
        check(d == -1.0 && cls == TorusDevClass::Unhandled, "off-axis circle (0.5 mm) -> refused");
        // ... and just inside Precision::Confusion() it is admitted, with a value
        // that is still >= the truth (the box), within Confusion of the coaxial value.
        const gp_Pnt c2 = ax.Location().Translated(gp_Vec(xa) * (0.5 * Precision::Confusion()));
        const gp_Circ q2(gp_Ax2(c2, za, xa), 10.0);
        const double d2 = circleOnTorusMax(t, q2, &cls);
        const double coax = circleOnTorusMax(t, coaxial(t, 10.0, 0.0));
        check(cls == TorusDevClass::Coaxial && d2 >= coax && d2 <= coax + Precision::Confusion(),
              "centre within Confusion of the axis: admitted, value >= coaxial, <= +Confusion");
        const SampleStats s = numCircleTorusSup(t, q2, d2, 1e-12, 2000);
        check(s.over == 0, "centre within Confusion: still an upper bound at every sample");
    }
    // anti-parallel normal: only the axis LINE matters -> admitted, same value
    {
        const gp_Circ up(gp_Ax2(ax.Location(), za, xa), 10.0);
        const gp_Circ down(gp_Ax2(ax.Location(), za.Reversed(), xa), 10.0);
        TorusDevClass c1 = TorusDevClass::Unhandled, c2 = TorusDevClass::Unhandled;
        const double a = circleOnTorusMax(t, up, &c1);
        const double b = circleOnTorusMax(t, down, &c2);
        check(c1 == TorusDevClass::Coaxial && c2 == TorusDevClass::Coaxial && a == b,
              "anti-parallel normal is still coaxial, identical value");
    }
    // negative tolerance -> refused; degenerate torus (no tube) -> refused
    {
        const double neg = circleOnTorusMax(t, coaxial(t, 10.0, 0.0), &cls, -1.0);
        check(neg == -1.0 && cls == TorusDevClass::Unhandled, "negative tolerance -> refused");
        const gp_Torus noTube(tiltedPlacement(), 7.0, 0.0);
        check(circleOnTorusMax(noTube, coaxial(noTube, 7.0, 0.0), &cls) == -1.0 &&
                  cls == TorusDevClass::Unhandled && !torusIsRing(noTube),
              "Rmin = 0 (no tube) -> refused");
        gp_Torus inf = t;
        inf.SetMajorRadius(std::numeric_limits<double>::infinity());
        check(circleOnTorusMax(inf, coaxial(t, 10.0, 0.0), &cls) == -1.0, "non-finite Rmaj -> refused");
        check(std::string(torusDevClassName(TorusDevClass::Unhandled)) == "unhandled-circle-on-torus" &&
                  !torusDevClassIsExact(TorusDevClass::Unhandled),
              "Unhandled class name / not exact");
    }
}

// ---------------------------------------------------------------------------
// 6. placement invariance (the cone lane's measured regression)
// ---------------------------------------------------------------------------

// Deterministic LCG — the battery must reproduce byte-identically (I5).
static uint64_t gSeed = 0xD1B54A32D192ED03ull;
static double rnd(double lo, double hi) {
    gSeed = gSeed * 6364136223846793005ull + 1442695040888963407ull;
    const double u = (double)((gSeed >> 11) & ((1ull << 53) - 1)) / (double)(1ull << 53);
    return lo + (hi - lo) * u;
}

static void testPlacementInvariance() {
    int nNotExact = 0;
    double worst = 0.0;
    for (int i = 0; i < 2000; ++i) {
        gp_XYZ z(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
        if (z.Modulus() < 0.1) z = gp_XYZ(0.3, -0.5, 0.81);
        gp_XYZ x(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
        x -= z * (x.Dot(z) / z.SquareModulus());
        if (x.Modulus() < 0.1) x = z.Crossed(gp_XYZ(1, 0, 0));
        const gp_Ax3 pl(gp_Pnt(rnd(-50, 50), rnd(-50, 50), rnd(-50, 50)), gp_Dir(z), gp_Dir(x));
        const gp_Torus t(pl, 7.0, 3.0);
        for (int k = 0; k < 4; ++k) {
            gp_Circ c;
            torusVIsoCircle(t, 0.5 * kPi * (double)k, c);
            TorusDevClass cls = TorusDevClass::Unhandled;
            const double d = circleOnTorusMax(t, c, &cls);
            if (cls != TorusDevClass::Coaxial) ++nNotExact;
            worst = std::max(worst, d);
        }
    }
    check(nNotExact == 0, "2000 random placements x 4 rims: every rim classes circle-on-torus");
    checkNear(worst, 0.0, 1e-13, "2000 random placements x 4 rims: worst rim deviation");
}

int main() {
    testParametrisationAndPointDist();
    testTangencyCirclesExactZero();
    testGeneralCoaxial();
    testSphereLimitAndSpindle();
    testRefusals();
    testPlacementInvariance();

    std::fprintf(stderr, "torus_math_unit: maxAbsErr(closed form vs sampled supremum) = %.3g\n",
                 gMaxAbsErr);
    std::fprintf(stderr, "torus_math_unit: %d/%d PASS\n", gPass, gPass + gFail);
    return gFail ? 1 : 0;
}
