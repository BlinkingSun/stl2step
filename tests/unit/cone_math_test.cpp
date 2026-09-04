// 130-CONE-MATH — certifies the closed-form cone math in src/refit_cone_math.cpp
// against independent ground truth (SPEC-130-cone-math, D-130-2, D-130-3).
//
// Two independent oracles, neither of which shares a line of derivation with
// the code under test:
//   * numPointConeDist  — minimises the point-to-line distance over the cone's
//     meridian LINES, built straight from the OCCT parametric equation
//     S(u,v) = O + (R + v sinA)(cos u X + sin u Y) + v cosA Z. It never uses
//     the (rho, z) reduction the implementation is built on.
//   * GeomAPI_ProjectPointOnSurf on a Geom_ConicalSurface, and
//     IntAna_QuadQuadGeo for the intersection cases — OCCT's own answers, and
//     IntAna is the very code the bind site (pickIntAna) runs.
//
// Numbers under test are the battery's: B1 (45.017 deg, R 8.5 -> 9.5) and the
// synthetic B3/B4/B5 (45 deg, R 10 -> 12, R 12 -> 10, R 10 -> 11.5).
//
// SPDX-License-Identifier: MIT

#include "refit_cone_math.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <ElCLib.hxx>
#include <Geom_ConicalSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

using stl2step::refit::ConeDevClass;
using stl2step::refit::circleOnConeMax;
using stl2step::refit::coneApexV;
using stl2step::refit::coneApexZ;
using stl2step::refit::coneAxialCoord;
using stl2step::refit::coneDevClassIsExact;
using stl2step::refit::coneDevClassName;
using stl2step::refit::coneFromRims;
using stl2step::refit::coneFrustumChordVolume;
using stl2step::refit::coneFrustumDVolDensity;
using stl2step::refit::coneFrustumVRange;
using stl2step::refit::coneFrustumVRangeFromRadii;
using stl2step::refit::coneIntConeCircle;
using stl2step::refit::coneIntCylCircle;
using stl2step::refit::coneIntPlaneCircle;
using stl2step::refit::coneRadialCoord;
using stl2step::refit::coneRadiusAt;
using stl2step::refit::coneVAtRadius;
using stl2step::refit::coneVAtZ;
using stl2step::refit::pointConeDist;

static const double kPi = 3.14159265358979323846264338327950288;

static int gPass = 0;
static int gFail = 0;

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

// ---------------------------------------------------------------------------
// oracle 1: distance to the cone from the OCCT parametric equation directly
// ---------------------------------------------------------------------------

// Distance from p to the FULL meridian line at parameter u. Every point of that
// line is S(u, v) for some v, so this is a distance to a genuine surface point.
static double meridianLineDist(const gp_Cone& c, const gp_Pnt& p, double u) {
    const gp_Ax3& ax = c.Position();
    const gp_XYZ O = ax.Location().XYZ();
    const gp_XYZ er = ax.XDirection().XYZ() * std::cos(u) + ax.YDirection().XYZ() * std::sin(u);
    const double A = c.SemiAngle();
    const gp_XYZ P0 = O + er * c.RefRadius();
    const gp_XYZ d = er * std::sin(A) + ax.Direction().XYZ() * std::cos(A);  // unit
    gp_XYZ w = p.XYZ() - P0;
    w -= d * w.Dot(d);
    return w.Modulus();
}

static double numPointConeDist(const gp_Cone& c, const gp_Pnt& p, int nU = 512) {
    double bestU = 0.0;
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < nU; ++i) {
        const double u = 2.0 * kPi * (double)i / (double)nU;
        const double d = meridianLineDist(c, p, u);
        if (d < best) {
            best = d;
            bestU = u;
        }
    }
    const double h = 2.0 * kPi / (double)nU;
    double lo = bestU - h;
    double hi = bestU + h;
    for (int it = 0; it < 60; ++it) {
        const double m1 = lo + (hi - lo) / 3.0;
        const double m2 = hi - (hi - lo) / 3.0;
        if (meridianLineDist(c, p, m1) <= meridianLineDist(c, p, m2))
            hi = m2;
        else
            lo = m1;
    }
    return std::min(best, meridianLineDist(c, p, 0.5 * (lo + hi)));
}

// Supremum over the circle of the distance to the cone, by sampling + a ternary
// refinement of the best bracket. A sample is always a LOWER bound on the true
// supremum, which is the direction the bound assertions need.
static double numCircleConeSup(const gp_Cone& c, const gp_Circ& q, int nT = 512, int nU = 512) {
    double bestT = 0.0;
    double best = -1.0;
    for (int i = 0; i < nT; ++i) {
        const double t = 2.0 * kPi * (double)i / (double)nT;
        const double d = numPointConeDist(c, ElCLib::Value(t, q), nU);
        if (d > best) {
            best = d;
            bestT = t;
        }
    }
    const double h = 2.0 * kPi / (double)nT;
    double lo = bestT - h;
    double hi = bestT + h;
    for (int it = 0; it < 40; ++it) {
        const double m1 = lo + (hi - lo) / 3.0;
        const double m2 = hi - (hi - lo) / 3.0;
        if (numPointConeDist(c, ElCLib::Value(m1, q), nU) >=
            numPointConeDist(c, ElCLib::Value(m2, q), nU))
            hi = m2;
        else
            lo = m1;
    }
    return std::max(best, numPointConeDist(c, ElCLib::Value(0.5 * (lo + hi), q), nU));
}

// ---------------------------------------------------------------------------
// battery frusta
// ---------------------------------------------------------------------------

struct Frustum {
    const char* name;
    double halfAngleDeg;
    double rhoA;  // rim in the cone's reference plane
    double rhoB;  // rim `height` further along the axis
    double height;
};

// height is derived, never assumed: |rhoB - rhoA| / tan(halfAngle). The SIGN of
// the semi-angle then falls out of (rhoB - rhoA), so B4 (12 -> 10) exercises the
// negative-semi-angle branch that OCCT allows and gp_Cone::Apex() flips on.
static Frustum mkFrustum(const char* name, double halfAngleDeg, double rhoA, double rhoB) {
    Frustum f;
    f.name = name;
    f.halfAngleDeg = halfAngleDeg;
    f.rhoA = rhoA;
    f.rhoB = rhoB;
    f.height = std::fabs(rhoB - rhoA) / std::tan(halfAngleDeg * kPi / 180.0);
    return f;
}

static std::vector<Frustum> battery() {
    return {mkFrustum("B1 linkage chamfer", 45.017, 8.5, 9.5),
            mkFrustum("B3 chamfer_straight_ring", 45.0, 10.0, 12.0),
            mkFrustum("B4 boss_cone_chamfer", 45.0, 12.0, 10.0),
            mkFrustum("B5 counterbore_chamfer", 45.0, 10.0, 11.5)};
}

// A deliberately off-origin, obliquely-oriented placement: nothing here may
// depend on the cone sitting at the world origin along Z.
static gp_Ax3 tiltedPlacement() {
    return gp_Ax3(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3));
}

// ---------------------------------------------------------------------------
// 1. parametrisation + frustum helpers
// ---------------------------------------------------------------------------

static void testParametrisation() {
    const gp_Ax3 pl = tiltedPlacement();
    for (const Frustum& f : battery()) {
        gp_Cone cone;
        char nm[256];
        std::snprintf(nm, sizeof nm, "%s: coneFromRims", f.name);
        check(coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone), nm);

        const double signedAng = std::atan((f.rhoB - f.rhoA) / f.height);
        std::snprintf(nm, sizeof nm, "%s: semi-angle magnitude", f.name);
        checkNear(std::fabs(cone.SemiAngle()) * 180.0 / kPi, f.halfAngleDeg, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "%s: semi-angle sign follows the rims", f.name);
        checkNear(cone.SemiAngle(), signedAng, 1e-15, nm);

        std::snprintf(nm, sizeof nm, "%s: coneRadiusAt(0) == rhoA", f.name);
        checkNear(coneRadiusAt(cone, 0.0), f.rhoA, 1e-13, nm);
        std::snprintf(nm, sizeof nm, "%s: coneRadiusAt(height) == rhoB", f.name);
        checkNear(coneRadiusAt(cone, f.height), f.rhoB, 1e-12, nm);

        // The parametric equation itself: S(u,v) must have radial coordinate
        // |R + v sinA| and axial coordinate v cosA, for arbitrary (u, v).
        Handle(Geom_ConicalSurface) gs = new Geom_ConicalSurface(cone);
        double worstRad = 0.0;
        double worstAx = 0.0;
        double worstRz = 0.0;
        for (int i = 0; i < 37; ++i) {
            const double u = 2.0 * kPi * (double)i / 37.0;
            for (int j = -6; j <= 6; ++j) {
                const double v = 1.7 * (double)j;
                const gp_Pnt P = gs->Value(u, v);
                const double z = coneAxialCoord(cone, P);
                const double rho = coneRadialCoord(cone, P);
                worstRad = std::max(worstRad,
                                    std::fabs(rho - std::fabs(cone.RefRadius() +
                                                              v * std::sin(cone.SemiAngle()))));
                worstAx = std::max(worstAx, std::fabs(z - v * std::cos(cone.SemiAngle())));
                // and the identity every later formula rests on: rho == |rho_cone(z)|
                worstRz = std::max(worstRz, std::fabs(rho - std::fabs(coneRadiusAt(cone, z))));
            }
        }
        std::snprintf(nm, sizeof nm, "%s: S(u,v) radial == |R + v sinA|", f.name);
        checkNear(worstRad, 0.0, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "%s: S(u,v) axial == v cosA", f.name);
        checkNear(worstAx, 0.0, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "%s: surface trace is rho == |rho_cone(z)|", f.name);
        checkNear(worstRz, 0.0, 1e-12, nm);

        // v <-> z <-> rho round trips.
        std::snprintf(nm, sizeof nm, "%s: coneVAtZ(height) == height/cosA", f.name);
        checkNear(coneVAtZ(cone, f.height), f.height / std::cos(cone.SemiAngle()), 1e-12, nm);
        std::snprintf(nm, sizeof nm, "%s: coneVAtRadius(rhoB) == coneVAtZ(height)", f.name);
        checkNear(coneVAtRadius(cone, f.rhoB), coneVAtZ(cone, f.height), 1e-10, nm);
    }
}

static void testApexAndVRange() {
    const gp_Ax3 pl = tiltedPlacement();
    for (const Frustum& f : battery()) {
        gp_Cone cone;
        if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
            check(false, "apex: coneFromRims");
            continue;
        }
        char nm[256];

        const double zApex = coneApexZ(cone);
        std::snprintf(nm, sizeof nm, "%s: coneApexZ agrees with gp_Cone::Apex()", f.name);
        checkNear(zApex, coneAxialCoord(cone, cone.Apex()), 1e-9, nm);
        std::snprintf(nm, sizeof nm, "%s: coneRadiusAt(apex) == 0", f.name);
        checkNear(coneRadiusAt(cone, zApex), 0.0, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "%s: coneApexV == coneVAtZ(apexZ)", f.name);
        checkNear(coneApexV(cone), coneVAtZ(cone, zApex), 1e-9, nm);
        std::snprintf(nm, sizeof nm, "%s: apex lies outside [0, height]", f.name);
        check(zApex < 0.0 ? (zApex < -1e-9) : (zApex > f.height + 1e-9), nm);

        double vLo = 0.0;
        double vHi = 0.0;
        std::snprintf(nm, sizeof nm, "%s: coneFrustumVRange accepts", f.name);
        check(coneFrustumVRange(cone, 0.0, f.height, vLo, vHi), nm);
        std::snprintf(nm, sizeof nm, "%s: v-range ascending", f.name);
        check(vLo < vHi, nm);
        std::snprintf(nm, sizeof nm, "%s: v-range EXCLUDES the apex", f.name);
        check(coneApexV(cone) < vLo || coneApexV(cone) > vHi, nm);
        std::snprintf(nm, sizeof nm, "%s: v-range spans the two rim radii", f.name);
        checkNear(vHi - vLo, std::fabs(f.rhoB - f.rhoA) / std::fabs(std::sin(cone.SemiAngle())),
                  1e-10, nm);

        double wLo = 0.0;
        double wHi = 0.0;
        std::snprintf(nm, sizeof nm, "%s: FromRadii accepts", f.name);
        check(coneFrustumVRangeFromRadii(cone, f.rhoA, f.rhoB, wLo, wHi), nm);
        std::snprintf(nm, sizeof nm, "%s: FromRadii == FromZ", f.name);
        check(std::fabs(wLo - vLo) <= 1e-9 && std::fabs(wHi - vHi) <= 1e-9, nm);

        // The rims are exactly on the surface at those v. Which END carries which
        // rim depends on the SIGN of the semi-angle (v ascending shrinks the
        // radius when Ang < 0), so this asserts the pair, not an order.
        Handle(Geom_ConicalSurface) gs = new Geom_ConicalSurface(cone);
        std::snprintf(nm, sizeof nm, "%s: vLo/vHi carry the two rim radii", f.name);
        const double rLo = coneRadialCoord(cone, gs->Value(0.0, vLo));
        const double rHi = coneRadialCoord(cone, gs->Value(0.0, vHi));
        check(std::fabs(std::min(rLo, rHi) - std::min(f.rhoA, f.rhoB)) <= 1e-10 &&
                  std::fabs(std::max(rLo, rHi) - std::max(f.rhoA, f.rhoB)) <= 1e-10,
              nm);
        std::snprintf(nm, sizeof nm, "%s: v ascends with the radius iff Ang > 0", f.name);
        check((cone.SemiAngle() > 0.0) == (rHi > rLo), nm);

        // Refusals: zero height, and a rim at/behind the apex.
        std::snprintf(nm, sizeof nm, "%s: v-range refuses zero height", f.name);
        check(!coneFrustumVRange(cone, f.height, f.height, vLo, vHi), nm);
        std::snprintf(nm, sizeof nm, "%s: v-range refuses a rim past the apex", f.name);
        check(!coneFrustumVRange(cone, zApex + (zApex < 0.0 ? -1.0 : 1.0), 0.0, vLo, vHi), nm);
        std::snprintf(nm, sizeof nm, "%s: FromRadii refuses a non-positive rim", f.name);
        check(!coneFrustumVRangeFromRadii(cone, 0.0, f.rhoB, wLo, wHi), nm);
        std::snprintf(nm, sizeof nm, "%s: FromRadii refuses equal rims", f.name);
        check(!coneFrustumVRangeFromRadii(cone, f.rhoA, f.rhoA, wLo, wHi), nm);
    }

    gp_Cone tmp;
    check(!coneFromRims(tiltedPlacement(), 10.0, 10.0, 2.0, tmp),
          "coneFromRims refuses equal rims (cylinder, not a cone)");
    check(!coneFromRims(tiltedPlacement(), 10.0, 12.0, 0.0, tmp),
          "coneFromRims refuses zero height");
    check(!coneFromRims(tiltedPlacement(), -1.0, 12.0, 2.0, tmp),
          "coneFromRims refuses a negative rim");
}

// ---------------------------------------------------------------------------
// 2. exact point-to-surface distance
// ---------------------------------------------------------------------------

static void testPointConeDist() {
    const gp_Ax3 pl = tiltedPlacement();
    const Frustum f = mkFrustum("B1", 45.017, 8.5, 9.5);
    gp_Cone cone;
    if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
        check(false, "pointConeDist: coneFromRims");
        return;
    }
    Handle(Geom_ConicalSurface) gs = new Geom_ConicalSurface(cone);

    // Points strewn over the whole surface neighbourhood, INCLUDING past the
    // apex (v very negative), where the naive "distance to the near ray" answer
    // is wrong and the double-nappe formula is right.
    double worstNum = 0.0;
    double worstOcc = 0.0;
    int n = 0;
    for (int i = 0; i < 13; ++i) {
        const double u = 2.0 * kPi * (double)i / 13.0;
        for (int j = -8; j <= 8; ++j) {
            const double v = 3.1 * (double)j;
            for (int k = -2; k <= 2; ++k) {
                const double off = 0.85 * (double)k;
                // push the surface point off along the surface normal-ish radial
                gp_Pnt S = gs->Value(u, v);
                const gp_XYZ a = cone.Axis().Direction().XYZ();
                gp_XYZ rad = S.XYZ() - cone.Location().XYZ();
                rad -= a * rad.Dot(a);
                const double m = rad.Modulus();
                if (m < 1e-12) continue;
                rad /= m;
                const gp_Pnt P(S.XYZ() + rad * off + a * (0.4 * (double)k));
                const double got = pointConeDist(cone, P);
                worstNum = std::max(worstNum, std::fabs(got - numPointConeDist(cone, P)));
                GeomAPI_ProjectPointOnSurf pr(P, gs);
                if (pr.IsDone() && pr.NbPoints() > 0)
                    worstOcc = std::max(worstOcc, std::fabs(got - pr.LowerDistance()));
                ++n;
            }
        }
    }
    std::fprintf(stderr, "pointConeDist: %d probes\n", n);
    checkNear(worstNum, 0.0, 1e-9, "pointConeDist == meridian-line oracle");
    checkNear(worstOcc, 0.0, 1e-7, "pointConeDist == GeomAPI_ProjectPointOnSurf");

    // On the surface it is exactly zero.
    double worstOn = 0.0;
    for (int i = 0; i < 11; ++i)
        for (int j = -5; j <= 5; ++j)
            worstOn = std::max(worstOn, std::fabs(pointConeDist(
                                            cone, gs->Value(2.0 * kPi * (double)i / 11.0,
                                                            2.3 * (double)j))));
    checkNear(worstOn, 0.0, 1e-12, "pointConeDist == 0 on the surface");
}

// ---------------------------------------------------------------------------
// 3. circleOnConeMax
// ---------------------------------------------------------------------------

// The rim circle of the frustum at axial height z (radius = the cone radius
// there), optionally displaced radially and axially.
static gp_Circ rimCircle(const gp_Cone& cone, double z, double dRadius, double dAxial) {
    const gp_Dir za = cone.Axis().Direction();
    const gp_Pnt c = cone.Location().Translated(gp_Vec(za) * (z + dAxial));
    const double rho = std::fabs(coneRadiusAt(cone, z)) + dRadius;
    return gp_Circ(gp_Ax2(c, za, cone.Position().XDirection()), rho);
}

static void testCircleOnConeCoaxial() {
    const gp_Ax3 pl = tiltedPlacement();
    for (const Frustum& f : battery()) {
        gp_Cone cone;
        if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
            check(false, "coaxial: coneFromRims");
            continue;
        }
        const double ca = std::cos(cone.SemiAngle());
        const double sa = std::fabs(std::sin(cone.SemiAngle()));
        char nm[256];

        for (int r = 0; r < 2; ++r) {
            const double z = (r == 0) ? 0.0 : f.height;
            ConeDevClass cls = ConeDevClass::Unhandled;

            // (i) the rim itself: deviation 0, class "circle-on-cone".
            const double d0 = circleOnConeMax(cone, rimCircle(cone, z, 0.0, 0.0), &cls);
            std::snprintf(nm, sizeof nm, "%s rim%d: exact rim deviates by 0", f.name, r);
            checkNear(d0, 0.0, 1e-12, nm);
            std::snprintf(nm, sizeof nm, "%s rim%d: class is circle-on-cone", f.name, r);
            check(cls == ConeDevClass::CoaxialRim && coneDevClassIsExact(cls) &&
                      std::string(coneDevClassName(cls)) == "circle-on-cone",
                  nm);

            // (ii) radial displacement delta -> delta * cos(halfAngle), exactly.
            const double dr = 0.25;
            const double dRad = circleOnConeMax(cone, rimCircle(cone, z, dr, 0.0), &cls);
            std::snprintf(nm, sizeof nm, "%s rim%d: radial %.2f -> delta*cosA", f.name, r, dr);
            checkNear(dRad, dr * ca, 1e-12, nm);

            // (iii) axial displacement delta -> delta * |sin(halfAngle)|, exactly.
            const double dz = 0.4;
            const double dAx = circleOnConeMax(cone, rimCircle(cone, z, 0.0, dz), &cls);
            std::snprintf(nm, sizeof nm, "%s rim%d: axial %.2f -> delta*sinA", f.name, r, dz);
            checkNear(dAx, dz * sa, 1e-12, nm);

            // (iv) both, against the independent oracle.
            const gp_Circ q = rimCircle(cone, z, -0.6, 0.35);
            const double got = circleOnConeMax(cone, q, &cls);
            const double want = numCircleConeSup(cone, q, 64, 512);
            std::snprintf(nm, sizeof nm, "%s rim%d: displaced rim == oracle", f.name, r);
            checkNear(got, want, 1e-9, nm);
        }

        // A coaxial circle far past the apex still gets the exact double-nappe
        // answer (the near-ray reading would be wrong here).
        ConeDevClass cls = ConeDevClass::Unhandled;
        const double zPast = coneApexZ(cone) + (coneApexZ(cone) < 0.0 ? -20.0 : 20.0);
        const gp_Circ q = rimCircle(cone, zPast, 1.3, 0.0);
        const double got = circleOnConeMax(cone, q, &cls);
        std::snprintf(nm, sizeof nm, "%s: coaxial circle past the apex == oracle", f.name);
        checkNear(got, numCircleConeSup(cone, q, 32, 512), 1e-9, nm);
        std::snprintf(nm, sizeof nm, "%s: past-apex circle is still exact", f.name);
        check(coneDevClassIsExact(cls), nm);
    }
}

static void testCircleOnConeOffAxisPerp() {
    const gp_Ax3 pl = tiltedPlacement();
    const Frustum f = mkFrustum("B5", 45.0, 10.0, 11.5);
    gp_Cone cone;
    if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
        check(false, "perp: coneFromRims");
        return;
    }
    const gp_Dir za = cone.Axis().Direction();
    const gp_Dir xa = cone.Position().XDirection();

    // Centre off the axis, plane still perpendicular: EXACT, and it must equal
    // the oracle for a centre both inside and outside the cone radius.
    const double offs[] = {0.75, 3.0, 9.0, 14.0};
    const double radii[] = {1.5, 6.0, 11.0};
    double worst = 0.0;
    bool allPerp = true;
    for (double o : offs) {
        for (double rr : radii) {
            const gp_Pnt c = cone.Location().Translated(gp_Vec(za) * (0.5 * f.height) +
                                                        gp_Vec(xa) * o);
            const gp_Circ q(gp_Ax2(c, za, xa), rr);
            ConeDevClass cls = ConeDevClass::Unhandled;
            const double got = circleOnConeMax(cone, q, &cls);
            allPerp = allPerp && cls == ConeDevClass::PerpPlane && coneDevClassIsExact(cls);
            worst = std::max(worst, std::fabs(got - numCircleConeSup(cone, q, 256, 256)));
        }
    }
    checkNear(worst, 0.0, 1e-8, "off-axis perpendicular circle == oracle (exact)");
    check(allPerp, "off-axis perpendicular circle classes as circle-on-cone-perp");
}

// Deterministic LCG — the battery must reproduce byte-identically (I5).
static uint64_t gSeed = 0x9E3779B97F4A7C15ull;
static double rnd(double lo, double hi) {
    gSeed = gSeed * 6364136223846793005ull + 1442695040888963407ull;
    const double u = (double)((gSeed >> 11) & ((1ull << 53) - 1)) / (double)(1ull << 53);
    return lo + (hi - lo) * u;
}

// The class must not depend on where the cone sits in world space: a rim whose
// normal IS the cone axis is exact for every placement. It is not automatic —
// sin(beta) from sqrt(1 - dot*dot) reports a spurious ~1e-8 tilt whenever the
// axis coordinates dot to 1 -+ 1ulp, which demoted a third of all placements to
// GeneralBound and inflated the deviation past Precision::Confusion(). This is
// the regression test for that; the cross-product form is exactly zero.
// Its own stream, so adding this test cannot shift the circles
// testCircleOnConeGeneralBound draws (I5: the battery reproduces exactly).
static uint64_t gPlaceSeed = 0xD1B54A32D192ED03ull;
static double rndPlace(double lo, double hi) {
    gPlaceSeed = gPlaceSeed * 6364136223846793005ull + 1442695040888963407ull;
    const double u = (double)((gPlaceSeed >> 11) & ((1ull << 53) - 1)) / (double)(1ull << 53);
    return lo + (hi - lo) * u;
}

static void testCoaxialClassIsPlacementInvariant() {
    int nCases = 0;
    int nNotExact = 0;
    int nNotCoaxial = 0;
    double worstDev = 0.0;
    double worstRimDev = 0.0;
    for (int i = 0; i < 400; ++i) {
        gp_Ax3 pl;
        try {
            const gp_Dir za(rndPlace(-1.0, 1.0), rndPlace(-1.0, 1.0), rndPlace(-1.0, 1.0));
            gp_Vec v(rndPlace(-1.0, 1.0), rndPlace(-1.0, 1.0), rndPlace(-1.0, 1.0));
            v -= gp_Vec(za) * v.Dot(gp_Vec(za));
            if (v.Magnitude() < 1e-3) continue;
            pl = gp_Ax3(gp_Pnt(rndPlace(-40.0, 40.0), rndPlace(-40.0, 40.0), rndPlace(-40.0, 40.0)), za,
                        gp_Dir(v));
        } catch (const Standard_Failure&) {
            continue;
        }
        gp_Cone cone;
        if (!coneFromRims(pl, 8.5, 9.5, 1.0, cone)) continue;
        ++nCases;

        // (i) the exact rim, taken from the tier-1 intersection itself.
        const gp_Dir za = cone.Axis().Direction();
        gp_Circ rim;
        if (!coneIntPlaneCircle(cone, gp_Pln(cone.Location(), za), rim)) {
            ++nNotCoaxial;
            continue;
        }
        ConeDevClass cls = ConeDevClass::Unhandled;
        const double dev = circleOnConeMax(cone, rim, &cls);
        if (cls != ConeDevClass::CoaxialRim) ++nNotCoaxial;
        worstRimDev = std::max(worstRimDev, std::fabs(dev));

        // (ii) an off-axis circle in the same perpendicular plane stays exact,
        // and matches the closed form cosA * max(d0 + r - xi, xi - |d0 - r|)
        // written out here from the derivation, not from the implementation.
        const double d0 = 2.5;
        const double rr = 4.0;
        const double xi = 8.5;  // = rho_cone(0) = rhoA
        const gp_Pnt c = cone.Location().Translated(gp_Vec(cone.Position().XDirection()) * d0);
        const gp_Circ off(gp_Ax2(c, za, cone.Position().XDirection()), rr);
        ConeDevClass cls2 = ConeDevClass::Unhandled;
        const double d2 = circleOnConeMax(cone, off, &cls2);
        if (cls2 != ConeDevClass::PerpPlane) ++nNotExact;
        const double want = std::cos(cone.SemiAngle()) *
                            std::max(d0 + rr - xi, xi - std::fabs(d0 - rr));
        worstDev = std::max(worstDev, std::fabs(d2 - want));
    }
    std::fprintf(stderr,
                 "placement invariance: %d placements, %d not CoaxialRim, %d not PerpPlane, "
                 "worst rim dev %.3g, worst perp err %.3g\n",
                 nCases, nNotCoaxial, nNotExact, worstRimDev, worstDev);
    check(nCases >= 300, "placement sweep produced enough cones");
    check(nNotCoaxial == 0, "the exact rim is CoaxialRim for EVERY cone placement");
    check(nNotExact == 0, "the off-axis perpendicular circle is PerpPlane for EVERY placement");
    checkNear(worstRimDev, 0.0, 1e-13, "the exact rim deviates by 0 for every placement");
    checkNear(worstDev, 0.0, 1e-12, "the perpendicular closed form holds for every placement");
}

static void testCircleOnConeGeneralBound() {
    const gp_Ax3 pl = tiltedPlacement();
    const Frustum f = mkFrustum("B3", 45.0, 10.0, 12.0);
    gp_Cone cone;
    if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
        check(false, "bound: coneFromRims");
        return;
    }

    int nBound = 0;
    int nViolations = 0;
    int nUnhandled = 0;
    double worstSlack = 0.0;
    double tightest = 1e300;
    for (int i = 0; i < 150; ++i) {
        const double nx = rnd(-1.0, 1.0);
        const double ny = rnd(-1.0, 1.0);
        const double nz = rnd(-1.0, 1.0);
        const gp_Pnt c(rnd(-20.0, 20.0), rnd(-20.0, 20.0), rnd(-20.0, 20.0));
        const double r = rnd(0.3, 12.0);
        gp_Circ q;
        try {
            q = gp_Circ(gp_Ax2(c, gp_Dir(nx, ny, nz)), r);
        } catch (const Standard_Failure&) {
            continue;
        }
        ConeDevClass cls = ConeDevClass::Unhandled;
        const double got = circleOnConeMax(cone, q, &cls);
        if (cls == ConeDevClass::Unhandled) {
            ++nUnhandled;
            continue;
        }
        if (cls == ConeDevClass::GeneralBound) ++nBound;
        const double truth = numCircleConeSup(cone, q, 128, 128);
        // A sample of the supremum is a LOWER bound on it, so a valid upper
        // bound must never fall below it.
        if (got < truth - 1e-9) ++nViolations;
        worstSlack = std::max(worstSlack, got - truth);
        tightest = std::min(tightest, got - truth);
    }
    std::fprintf(stderr,
                 "general bound: %d circles bounded, %d unhandled, worst slack %.6g, "
                 "tightest %.6g\n",
                 nBound, nUnhandled, worstSlack, tightest);
    check(nBound >= 140, "random tilted circles class as circle-on-cone-bound");
    check(nUnhandled == 0, "no random circle is unhandled");
    check(nViolations == 0, "every bound is >= the true supremum");
    check(tightest >= -1e-9, "no bound falls below the sampled supremum");

    // A tilted circle whose centre IS on the axis: the bound must stay tight,
    // because a bound of order r there would refuse every near-perfect rim.
    const gp_Dir za = cone.Axis().Direction();
    const gp_Dir xa = cone.Position().XDirection();
    const double tiltDeg = 0.05;
    const double t = tiltDeg * kPi / 180.0;
    const gp_Dir tilted(za.XYZ() * std::cos(t) + xa.XYZ() * std::sin(t));
    const gp_Pnt c = cone.Location().Translated(gp_Vec(za) * (0.5 * f.height));
    const gp_Circ q(gp_Ax2(c, tilted, xa), std::fabs(coneRadiusAt(cone, 0.5 * f.height)));
    ConeDevClass cls = ConeDevClass::Unhandled;
    const double got = circleOnConeMax(cone, q, &cls);
    const double truth = numCircleConeSup(cone, q, 512, 512);
    check(cls == ConeDevClass::GeneralBound, "0.05 deg tilt classes as a bound");
    check(got >= truth - 1e-12, "0.05 deg tilt bound >= supremum");
    std::fprintf(stderr, "0.05 deg tilt: bound %.6g truth %.6g\n", got, truth);
    check(got <= 0.05, "0.05 deg tilt bound stays small (not O(r))");
}

static void testCircleOnConeUnhandled() {
    gp_Cone bad;
    bad.SetRadius(std::numeric_limits<double>::infinity());
    ConeDevClass cls = ConeDevClass::CoaxialRim;
    const double d = circleOnConeMax(bad, gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 5.0),
                                     &cls);
    check(d < 0.0 && cls == ConeDevClass::Unhandled, "degenerate cone -> -1 / Unhandled");
    check(std::string(coneDevClassName(ConeDevClass::Unhandled)) == "unhandled-circle-on-cone",
          "Unhandled class name");
    check(!coneDevClassIsExact(ConeDevClass::GeneralBound), "GeneralBound is not exact");
    check(std::string(coneDevClassName(ConeDevClass::GeneralBound)) == "circle-on-cone-bound",
          "GeneralBound class name");
    check(std::string(coneDevClassName(ConeDevClass::PerpPlane)) == "circle-on-cone-perp",
          "PerpPlane class name");

    const gp_Ax3 pl = tiltedPlacement();
    gp_Cone cone;
    coneFromRims(pl, 8.5, 9.5, 1.0, cone);
    const double neg =
        circleOnConeMax(cone, gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 5.0), &cls, -1.0);
    check(neg < 0.0 && cls == ConeDevClass::Unhandled, "negative tolerance -> Unhandled");
    check(pointConeDist(bad, gp_Pnt(0, 0, 0)) < 0.0, "pointConeDist rejects a degenerate cone");
    check(!std::isfinite(coneRadiusAt(bad, 1.0)), "coneRadiusAt is NaN on a degenerate cone");
}

// ---------------------------------------------------------------------------
// 4. tier-1 IntAna cone cases
// ---------------------------------------------------------------------------

static bool circMatches(const gp_Circ& a, const gp_Circ& b, double tol) {
    return std::fabs(a.Radius() - b.Radius()) <= tol &&
           a.Location().Distance(b.Location()) <= tol &&
           a.Axis().Direction().IsParallel(b.Axis().Direction(), 1e-7);
}

// Does IntAna_QuadQuadGeo report a circle equal to `want`?
static bool intAnaHasCircle(const IntAna_QuadQuadGeo& iq, const gp_Circ& want, double tol) {
    if (!iq.IsDone()) return false;
    const IntAna_ResultType ty = iq.TypeInter();
    if (ty != IntAna_Circle && ty != IntAna_PointAndCircle) return false;
    for (int i = 1; i <= iq.NbSolutions(); ++i) {
        try {
            if (circMatches(iq.Circle(i), want, tol)) return true;
        } catch (const Standard_Failure&) {
        }
    }
    return false;
}

static void testIntPlane() {
    const gp_Ax3 pl = tiltedPlacement();
    for (const Frustum& f : battery()) {
        gp_Cone cone;
        if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) continue;
        const gp_Dir za = cone.Axis().Direction();
        char nm[256];

        for (int r = 0; r < 2; ++r) {
            const double z = (r == 0) ? 0.0 : f.height;
            const double want = (r == 0) ? f.rhoA : f.rhoB;
            const gp_Pln plane(cone.Location().Translated(gp_Vec(za) * z), za);
            gp_Circ got;
            std::snprintf(nm, sizeof nm, "%s rim%d: cone^perp-plane accepts", f.name, r);
            check(coneIntPlaneCircle(cone, plane, got), nm);
            std::snprintf(nm, sizeof nm, "%s rim%d: circle radius == rim radius", f.name, r);
            checkNear(got.Radius(), want, 1e-12, nm);
            std::snprintf(nm, sizeof nm, "%s rim%d: circle lies ON the cone", f.name, r);
            checkNear(pointConeDist(cone, ElCLib::Value(0.7, got)), 0.0, 1e-12, nm);
            std::snprintf(nm, sizeof nm, "%s rim%d: circle lies IN the plane", f.name, r);
            checkNear(plane.Distance(ElCLib::Value(2.4, got)), 0.0, 1e-12, nm);
            std::snprintf(nm, sizeof nm, "%s rim%d: supremum on the cone is 0", f.name, r);
            checkNear(circleOnConeMax(cone, got), 0.0, 1e-12, nm);

            IntAna_QuadQuadGeo iq(plane, cone, Precision::Angular(), Precision::Confusion());
            std::snprintf(nm, sizeof nm, "%s rim%d: agrees with IntAna_QuadQuadGeo", f.name, r);
            check(intAnaHasCircle(iq, got, 1e-7), nm);
        }

        // Refusals: oblique plane, and the plane through the apex.
        gp_Circ got;
        const gp_Dir oblique(za.XYZ() * std::cos(0.2) +
                             cone.Position().XDirection().XYZ() * std::sin(0.2));
        std::snprintf(nm, sizeof nm, "%s: oblique plane -> none (tier 2)", f.name);
        check(!coneIntPlaneCircle(cone, gp_Pln(cone.Location(), oblique), got), nm);
        std::snprintf(nm, sizeof nm, "%s: apex plane -> none (a point, not a circle)", f.name);
        check(!coneIntPlaneCircle(cone, gp_Pln(cone.Apex(), za), got), nm);
    }
}

static void testIntCylinder() {
    const gp_Ax3 pl = tiltedPlacement();
    const Frustum f = mkFrustum("B1", 45.017, 8.5, 9.5);
    gp_Cone cone;
    if (!coneFromRims(pl, f.rhoA, f.rhoB, f.height, cone)) {
        check(false, "cyl: coneFromRims");
        return;
    }
    const gp_Dir za = cone.Axis().Direction();
    const gp_Dir xa = cone.Position().XDirection();

    // The B1 bore: a coaxial cylinder of radius 8.5 meets the chamfer cone at
    // its lower rim; radius 9.0 meets it mid-frustum.
    const double rcs[] = {8.5, 9.0, 9.5};
    for (double rc : rcs) {
        const gp_Cylinder cyl(gp_Ax3(cone.Location().Translated(gp_Vec(za) * 5.0), za, xa), rc);
        gp_Circ got;
        char nm[256];
        std::snprintf(nm, sizeof nm, "B1: cone^coaxial-cyl R=%.2f accepts", rc);
        check(coneIntCylCircle(cone, cyl, got), nm);
        std::snprintf(nm, sizeof nm, "B1: R=%.2f circle radius == R_cyl", rc);
        checkNear(got.Radius(), rc, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "B1: R=%.2f circle lies ON the cone", rc);
        checkNear(pointConeDist(cone, ElCLib::Value(1.1, got)), 0.0, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "B1: R=%.2f circle lies ON the cylinder", rc);
        checkNear(coneRadialCoord(cone, ElCLib::Value(1.1, got)), rc, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "B1: R=%.2f height solves rho_cone(z) = R_cyl", rc);
        checkNear(coneRadiusAt(cone, coneAxialCoord(cone, got.Location())), rc, 1e-12, nm);
        std::snprintf(nm, sizeof nm, "B1: R=%.2f supremum on the cone is 0", rc);
        checkNear(circleOnConeMax(cone, got), 0.0, 1e-12, nm);

        IntAna_QuadQuadGeo iq(cyl, cone, Precision::Confusion());
        std::snprintf(nm, sizeof nm, "B1: R=%.2f agrees with IntAna_QuadQuadGeo", rc);
        check(intAnaHasCircle(iq, got, 1e-7), nm);
    }

    // Refusals: axis offset off the cone axis, and a skew axis.
    gp_Circ got;
    const gp_Cylinder offset(
        gp_Ax3(cone.Location().Translated(gp_Vec(xa) * 0.5), za, xa), 9.0);
    check(!coneIntCylCircle(cone, offset, got), "B1: offset cylinder -> none (tier 2)");
    const gp_Dir skew(za.XYZ() * std::cos(0.3) + xa.XYZ() * std::sin(0.3));
    check(!coneIntCylCircle(cone, gp_Cylinder(gp_Ax3(cone.Location(), skew), 9.0), got),
          "B1: skew cylinder -> none (tier 2)");
}

static void testIntCone() {
    const gp_Ax3 pl = tiltedPlacement();
    const gp_Dir za = pl.Direction();
    const gp_Dir xa = pl.XDirection();

    gp_Cone a;
    gp_Cone b;
    if (!coneFromRims(pl, 10.0, 12.0, 2.0, a)) {  // +45 deg
        check(false, "cone^cone: coneFromRims a");
        return;
    }
    // A counter-chamfer: -45 deg, reference radius 13 at z = 1 in a's frame, so
    // rho_b(z) = 13 - (z - 1). rho_a(z) = 10 + z. They cross at z = 2, rho = 12.
    const gp_Ax3 pb(pl.Location().Translated(gp_Vec(za) * 1.0), za, xa);
    if (!coneFromRims(pb, 13.0, 11.0, 2.0, b)) {
        check(false, "cone^cone: coneFromRims b");
        return;
    }
    gp_Circ got;
    check(coneIntConeCircle(a, b, got), "cone^coaxial-cone accepts");
    checkNear(got.Radius(), 12.0, 1e-11, "cone^cone: crossing radius == 12");
    checkNear(coneAxialCoord(a, got.Location()), 2.0, 1e-11, "cone^cone: crossing at z == 2");
    checkNear(pointConeDist(a, ElCLib::Value(0.9, got)), 0.0, 1e-11, "cone^cone: circle on a");
    checkNear(pointConeDist(b, ElCLib::Value(0.9, got)), 0.0, 1e-11, "cone^cone: circle on b");
    checkNear(circleOnConeMax(a, got), 0.0, 1e-11, "cone^cone: supremum on a is 0");
    checkNear(circleOnConeMax(b, got), 0.0, 1e-11, "cone^cone: supremum on b is 0");
    {
        IntAna_QuadQuadGeo iq(a, b, Precision::Confusion());
        check(intAnaHasCircle(iq, got, 1e-7), "cone^cone agrees with IntAna_QuadQuadGeo");
    }

    // An ANTIPARALLEL axis for b must give the same circle: coaxiality is about
    // the axis LINE, not its sense.
    gp_Cone bFlip;
    const gp_Ax3 pbFlip(pl.Location().Translated(gp_Vec(za) * 3.0), gp_Dir(za.XYZ() * -1.0), xa);
    if (coneFromRims(pbFlip, 11.0, 13.0, 2.0, bFlip)) {
        gp_Circ gotFlip;
        check(coneIntConeCircle(a, bFlip, gotFlip), "cone^cone accepts an antiparallel axis");
        check(circMatches(got, gotFlip, 1e-9), "antiparallel axis gives the same circle");
    } else {
        check(false, "cone^cone: coneFromRims bFlip");
    }

    // Refusals.
    gp_Cone parallelGen;
    const gp_Ax3 pp(pl.Location().Translated(gp_Vec(za) * 0.0 + gp_Vec(xa) * 0.0), za, xa);
    if (coneFromRims(pp, 14.0, 16.0, 2.0, parallelGen))  // same +45 deg, different radius
        check(!coneIntConeCircle(a, parallelGen, got), "cone^cone: parallel generators -> none");
    else
        check(false, "cone^cone: coneFromRims parallelGen");

    gp_Cone offAxis;
    const gp_Ax3 po(pl.Location().Translated(gp_Vec(xa) * 0.5 + gp_Vec(za) * 1.0), za, xa);
    if (coneFromRims(po, 13.0, 11.0, 2.0, offAxis))
        check(!coneIntConeCircle(a, offAxis, got), "cone^cone: off-axis -> none (tier 2)");
    else
        check(false, "cone^cone: coneFromRims offAxis");

    gp_Cone oblique;
    const gp_Dir skew(za.XYZ() * std::cos(0.25) + xa.XYZ() * std::sin(0.25));
    if (coneFromRims(gp_Ax3(pl.Location(), skew), 13.0, 11.0, 2.0, oblique))
        check(!coneIntConeCircle(a, oblique, got), "cone^cone: oblique axis -> none (tier 2)");
    else
        check(false, "cone^cone: coneFromRims oblique");

    // A crossing where BOTH signed radii are negative is a real circle on the
    // far nappe of both double cones — outside every frustum's v-range, so it is
    // refused rather than shipped. rho_a(z) = 10 + z against a 60 deg cone
    // rho_c(z) = 20 + z*tan(60) cross at z = -13.66, rho = -3.66.
    gp_Cone farNappe;
    const double h60 = 2.0 / std::tan(60.0 * kPi / 180.0);
    if (coneFromRims(pl, 20.0, 22.0, h60, farNappe)) {
        checkNear(std::fabs(farNappe.SemiAngle()) * 180.0 / kPi, 60.0, 1e-12,
                  "cone^cone: far-nappe probe is a 60 deg cone");
        check(!coneIntConeCircle(a, farNappe, got),
              "cone^cone: far-nappe crossing (negative radius) -> none");
    } else {
        check(false, "cone^cone: coneFromRims farNappe");
    }
}

// ---------------------------------------------------------------------------
// oracle 4: the facet-to-frustum gap volume by numeric integration
// ---------------------------------------------------------------------------
//
// D-130-11(2). Independent of the closed form under test: it slices the frustum
// into nZ axial slabs, and in each slab it integrates the AREA BETWEEN the
// inscribed N-gon and the circle by a 2D sweep -- polygon area from the
// shoelace formula, circle area from pi r^2 -- then multiplies by the slab
// height. Nothing here knows the circular-segment identity
// (r^2/2)(gamma - sin gamma) that coneFrustumChordVolume is built on; it only
// knows what a regular polygon inscribed in a circle looks like.
static double numFrustumGapVolume(double R1, double R2, double height, int nSides, int nZ) {
    const double h = std::fabs(height);
    double v = 0.0;
    for (int i = 0; i < nZ; ++i) {
        // Midpoint in z of slab i. rho is linear in z, so the slab's exact
        // contribution is integral rho^2 dz over the slab; the midpoint rule is
        // NOT exact for rho^2 and Simpson is, so the slab uses Simpson in z and
        // the polygon/circle difference in the plane.
        const double z0 = h * (double)i / (double)nZ;
        const double z1 = h * (double)(i + 1) / (double)nZ;
        const auto gapAt = [&](double z) {
            const double rho = R1 + (z / h) * (R2 - R1);
            // Shoelace area of the regular N-gon inscribed in radius rho, built
            // vertex by vertex, against the circle's pi rho^2.
            double poly = 0.0;
            for (int k = 0; k < nSides; ++k) {
                const double a0 = 2.0 * kPi * (double)k / (double)nSides;
                const double a1 = 2.0 * kPi * (double)(k + 1) / (double)nSides;
                const double x0 = rho * std::cos(a0), y0 = rho * std::sin(a0);
                const double x1 = rho * std::cos(a1), y1 = rho * std::sin(a1);
                poly += x0 * y1 - x1 * y0;
            }
            poly *= 0.5;
            return kPi * rho * rho - poly;
        };
        v += (z1 - z0) / 6.0 * (gapAt(z0) + 4.0 * gapAt(0.5 * (z0 + z1)) + gapAt(z1));
    }
    return v;
}

// The plate's frustum, triangulated the way the mesh is: N quad facets, each
// split into two triangles by the rim-to-rim diagonal. Sums
// area(t) * kappa(rho at centroid(t)) over the triangles -- exactly what
// detector C does over the triangles a ChamferCone region claims.
static double sumDensityOverSkin(double R1, double R2, double height, int nSides) {
    const double h = std::fabs(height);
    double sum = 0.0;
    for (int k = 0; k < nSides; ++k) {
        const double a0 = 2.0 * kPi * (double)k / (double)nSides;
        const double a1 = 2.0 * kPi * (double)(k + 1) / (double)nSides;
        const gp_XYZ p00(R1 * std::cos(a0), R1 * std::sin(a0), 0.0);
        const gp_XYZ p10(R1 * std::cos(a1), R1 * std::sin(a1), 0.0);
        const gp_XYZ p01(R2 * std::cos(a0), R2 * std::sin(a0), h);
        const gp_XYZ p11(R2 * std::cos(a1), R2 * std::sin(a1), h);
        const gp_XYZ tri[2][3] = {{p00, p10, p11}, {p00, p11, p01}};
        for (int t = 0; t < 2; ++t) {
            const gp_XYZ e1 = tri[t][1] - tri[t][0];
            const gp_XYZ e2 = tri[t][2] - tri[t][0];
            const double area = 0.5 * e1.Crossed(e2).Modulus();
            const double zc = (tri[t][0].Z() + tri[t][1].Z() + tri[t][2].Z()) / 3.0;
            const double rho = R1 + (zc / h) * (R2 - R1);
            sum += area * coneFrustumDVolDensity(rho, R1, R2, h, nSides);
        }
    }
    return sum;
}

static void testFrustumChordVolume() {
    // D-130-11(2), the plate's own frustum: R 8.5 -> 9.5, 45 deg, 92 sides.
    const double R1 = 8.5, R2 = 9.5, H = 1.0;
    const int N = 92;
    const double closed = coneFrustumChordVolume(R1, R2, H, N);
    const double num = numFrustumGapVolume(R1, R2, H, N, 4096);
    checkNear(closed, num, 1e-9 * std::fabs(num) + 1e-12,
              "plate frustum: closed form == numeric integral (R 8.5->9.5, 45 deg, 92 sides)");
    // 130-CONE-ARCS 1.3 attributed -0.197976 mm^3 per plate frustum in closed
    // form, independently of this code. The two must be the same number.
    checkNear(closed, 0.197976, 5e-7, "plate frustum: 0.197976 mm^3 (130-CONE-ARCS 1.3)");

    // Summing the density over the triangulated skin IS the total, to the last
    // bit the arithmetic carries: kappa is affine over a planar facet.
    checkNear(sumDensityOverSkin(R1, R2, H, N), closed, 1e-12 * closed,
              "plate frustum: sum over claimed triangles == closed form");

    // The sign of the height must not change the volume -- the frame convention
    // 30f67dd landed carries the taper sign in h, and a volume is not signed by
    // a frame.
    checkNear(coneFrustumChordVolume(R1, R2, -H, N), closed, 0.0,
              "negative signed height gives the same volume");
    checkNear(sumDensityOverSkin(R1, R2, -H, N), closed, 1e-12 * closed,
              "negative signed height: density sum unchanged");

    // The battery's synthetic frustums (B3/B4/B5), against the same oracle.
    struct Case { double R1, R2, H; int N; const char* name; };
    const Case cases[] = {
        {10.0, 12.0, 2.0, 64, "B3 45 deg R 10->12"},
        {10.0, 12.0, 2.0, 24, "coarse N=24 R 10->12"},
        {12.0, 10.0, 2.0, 64, "B4 45 deg R 12->10 (taper reversed)"},
        {10.0, 11.5, 1.5, 128, "B5 45 deg R 10->11.5"},
        {0.75, 1.25, 0.5, 16, "small radii, N=16"},
        {30.0, 30.5, 8.0, 360, "shallow taper, N=360"},
    };
    for (const Case& c : cases) {
        const double cf = coneFrustumChordVolume(c.R1, c.R2, c.H, c.N);
        const double nf = numFrustumGapVolume(c.R1, c.R2, c.H, c.N, 4096);
        checkNear(cf, nf, 1e-9 * std::fabs(nf) + 1e-12, c.name);
        checkNear(sumDensityOverSkin(c.R1, c.R2, c.H, c.N), cf, 1e-11 * cf + 1e-15,
                  "density sum matches closed form");
    }

    // The cylinder limit: kappa collapses to dVolCylinderSector's density
    // R*(gamma - sin gamma)/(4 sin(gamma/2)). Stated here as the independent
    // expression, not by calling the engine (this TU does not link it).
    for (int N2 : {6, 24, 92, 360}) {
        const double R = 8.5;
        const double g = 2.0 * kPi / (double)N2;
        const double want = R * (g - std::sin(g)) / (4.0 * std::sin(0.5 * g));
        // R1 != R2 is required by a frustum, so approach the limit: a taper of
        // 1e-9 mm over a height of 1 mm.
        const double got = coneFrustumDVolDensity(R, R, R + 1e-9, 1.0, N2);
        checkNear(got, want, 1e-9 * want, "kappa -> dVolCylinderSector density");
    }

    // Degenerate input returns 0, never a plausible-looking number.
    check(coneFrustumChordVolume(8.5, 9.5, 1.0, 2) == 0.0, "N < 3 -> 0");
    check(coneFrustumChordVolume(-1.0, 9.5, 1.0, 92) == 0.0, "negative rim radius -> 0");
    check(coneFrustumChordVolume(8.5, 9.5, 0.0, 92) == 0.0, "zero height -> 0");
    check(coneFrustumDVolDensity(8.5, 8.5, 9.5, 0.0, 92) == 0.0, "density: zero height -> 0");
    check(coneFrustumDVolDensity(std::numeric_limits<double>::quiet_NaN(), 8.5, 9.5, 1.0, 92)
              == 0.0,
          "density: non-finite rho -> 0");
}

int main() {
    testParametrisation();
    testApexAndVRange();
    testPointConeDist();
    testCircleOnConeCoaxial();
    testCircleOnConeOffAxisPerp();
    testCoaxialClassIsPlacementInvariant();
    testCircleOnConeGeneralBound();
    testCircleOnConeUnhandled();
    testIntPlane();
    testIntCylinder();
    testIntCone();
    testFrustumChordVolume();

    std::fprintf(stderr, "cone_math_unit: %d/%d PASS\n", gPass, gPass + gFail);
    return gFail ? 1 : 0;
}
