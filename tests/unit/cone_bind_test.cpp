// 130-CONE-BUILD — certifies the `circle-on-cone` / `line-on-cone` BIND classes
// in src/refit_cone_bind.cpp (D-130-3's hard gate, D-130-2's invariant).
//
// The quantity under test is NOT the circle-to-surface distance 130-CONE-MATH
// certifies; it is what a bind site records,
//
//     sup_t | C3d(t) - S(pcurve(t)) |,
//
// and the two differ by exactly the parametrisation.  The oracle is therefore
// direct: evaluate that expression on a dense parameter grid with OCCT's own
// Geom_ConicalSurface / Geom_Circle / Geom2d_Line evaluators and compare against
// the closed form the code returns.  A grid sample is a LOWER bound on a
// supremum, which is the direction every assertion below needs — a returned
// value that is >= every sample and equal to the largest one is the supremum.
//
// Everything runs on an off-origin, obliquely-placed cone: nothing may depend on
// the cone sitting at the world origin along Z.  Frustum numbers are the
// battery's — B1 is the plate's own measured chamfer (45.000004 deg, R 8.5 ->
// 9.5, the vertex fit; the 45.017 in 130-arcs is that tessellation's FACET
// NORMAL angle, atan(1/cos(pi/92)), not the surface semi-angle), and B3/B4/B5
// are the synthetic fixtures (45 deg, 10 -> 12, 12 -> 10, 10 -> 11.5).
//
// SPDX-License-Identifier: MIT

#include "refit_cone_bind.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Lin.hxx>
#include <gp_Lin2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include "refit_cone_math.hpp"

using stl2step::refit::ConeDevClass;
using stl2step::refit::circleOnConeMax;
using stl2step::refit::coneBindSup;
using stl2step::refit::coneCircleBindSup;
using stl2step::refit::coneFromRims;
using stl2step::refit::coneFrustumVRangeFromRadii;
using stl2step::refit::coneLineBindSup;
using stl2step::refit::conePCurveForCircle;
using stl2step::refit::conePCurveForLine;
using stl2step::refit::coneRadiusAt;
using stl2step::refit::coneVAtRadius;
using stl2step::refit::coneVAtZ;

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

static void checkStr(const char* got, const char* want, const char* name) {
    const bool ok = got && want && std::string(got) == std::string(want);
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s (%s)\n", name, got ? got : "(null)");
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s got %s want %s\n", name, got ? got : "(null)", want);
    }
}

// ---------------------------------------------------------------------------
// the oracle: |C3d(t) - S(pc(t))| sampled directly, sharing no derivation with
// the code under test
// ---------------------------------------------------------------------------

struct GridStat {
    double maxGap = 0.0;
    double minGap = 0.0;
    bool valid = false;
};

static GridStat sampleGap(const Handle(Geom_Curve)& c3, double f, double l,
                          const Handle(Geom_Surface)& srf, const Handle(Geom2d_Curve)& c2d,
                          int n = 257) {
    GridStat s;
    if (c3.IsNull() || srf.IsNull() || c2d.IsNull()) return s;
    s.maxGap = -1.0;
    s.minGap = 1e300;
    for (int i = 0; i < n; ++i) {
        const double t = f + (l - f) * (double)i / (double)(n - 1);
        const gp_Pnt P = c3->Value(t);
        const gp_Pnt2d uv = c2d->Value(t);
        const gp_Pnt Q = srf->Value(uv.X(), uv.Y());
        const double d = P.Distance(Q);
        s.maxGap = std::max(s.maxGap, d);
        s.minGap = std::min(s.minGap, d);
    }
    s.valid = true;
    return s;
}

// ---------------------------------------------------------------------------
// the battery's frusta, all on ONE off-origin oblique placement
// ---------------------------------------------------------------------------

static gp_Ax3 obliquePlacement() {
    return gp_Ax3(gp_Pnt(-3.25, 7.5, 11.0), gp_Dir(0.3, -0.5, 0.81), gp_Dir(0.81, 0.0, -0.3));
}

struct Frustum {
    std::string name;
    gp_Cone cone;
    double rLo = 0.0, rHi = 0.0;
    double vLo = 0.0, vHi = 0.0;
};

// height DERIVED from the rims and the half-angle, never assumed: |dR|/tan(A).
static Frustum mkFrustum(const char* name, double halfDeg, double rA, double rB) {
    Frustum fr;
    fr.name = name;
    const double h = std::fabs(rB - rA) / std::tan(halfDeg * kPi / 180.0);
    if (!coneFromRims(obliquePlacement(), rA, rB, h, fr.cone)) return fr;
    fr.rLo = std::min(rA, rB);
    fr.rHi = std::max(rA, rB);
    if (!coneFrustumVRangeFromRadii(fr.cone, rA, rB, fr.vLo, fr.vHi)) fr.vLo = fr.vHi = 0.0;
    return fr;
}

static std::vector<Frustum> battery() {
    // B1 is the plate's measured chamfer; B3/B4/B5 the synthetic fixtures.
    return {mkFrustum("B1 plate chamfer", 45.000004, 8.5, 9.5), mkFrustum("B3", 45.0, 10.0, 12.0),
            mkFrustum("B4", 45.0, 12.0, 10.0), mkFrustum("B5", 45.0, 10.0, 11.5)};
}

// The exact rim circle at surface parameter v: centre on the axis at that
// height, radius the cone's own radius there, axis the cone axis. Built from the
// parametrisation, not from the code under test.
static gp_Circ rimAt(const gp_Cone& c, double v) {
    const gp_Ax3& ax = c.Position();
    const double A = c.SemiAngle();
    const double z = v * std::cos(A);
    const double rho = coneRadiusAt(c, z);
    const gp_Pnt ctr = ax.Location().Translated(gp_Vec(ax.Direction()) * z);
    return gp_Circ(gp_Ax2(ctr, ax.Direction(), ax.XDirection()), std::fabs(rho));
}

static Handle(Geom_Surface) surfOf(const gp_Cone& c) { return new Geom_ConicalSurface(c); }

// ---------------------------------------------------------------------------
// 1. an exact rim binds at zero and certifies as circle-on-cone
// ---------------------------------------------------------------------------

static void testExactRim() {
    const TopLoc_Location id;
    for (const Frustum& fr : battery()) {
        if (fr.vHi <= fr.vLo) {
            check(false, (fr.name + ": frustum v-range").c_str());
            continue;
        }
        Handle(Geom_Surface) S = surfOf(fr.cone);
        for (int end = 0; end < 2; ++end) {
            const double v = end ? fr.vHi : fr.vLo;
            const gp_Circ rim = rimAt(fr.cone, v);
            Handle(Geom_Curve) c3 = new Geom_Circle(rim);
            Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, rim);
            if (pc.IsNull()) {
                check(false, (fr.name + ": rim pcurve built").c_str());
                continue;
            }
            const char* cls = nullptr;
            const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, id, &cls);
            checkStr(cls, "circle-on-cone", (fr.name + ": exact rim class").c_str());
            checkNear(sup, 0.0, 1e-11, (fr.name + ": exact rim binds at 0").c_str());
            const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, pc);
            checkNear(g.maxGap, 0.0, 1e-11, (fr.name + ": oracle agrees the rim is on").c_str());
        }
    }
}

// ---------------------------------------------------------------------------
// 2. THE POINT OF THIS FILE: the bind supremum carries the parametrisation
//    A radial displacement delta sits delta*cos(A) from the SURFACE but delta
//    from the point its own pcurve names. Recording the surface figure would
//    understate a 45 deg chamfer's edge tolerance by a factor of sqrt(2).
// ---------------------------------------------------------------------------

static void testParametrisationTerm() {
    const TopLoc_Location id;
    for (const Frustum& fr : battery()) {
        Handle(Geom_Surface) S = surfOf(fr.cone);
        const double A = fr.cone.SemiAngle();
        const double v = 0.5 * (fr.vLo + fr.vHi);
        const gp_Circ rim = rimAt(fr.cone, v);
        for (const double delta : {1e-3, 5e-2}) {
            // radial: same centre and axis, radius + delta
            const gp_Circ off(rim.Position(), rim.Radius() + delta);
            Handle(Geom_Curve) c3 = new Geom_Circle(off);
            Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, off);
            const char* cls = nullptr;
            const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, id, &cls);
            checkStr(cls, "circle-on-cone", (fr.name + ": radial-offset class").c_str());
            checkNear(sup, delta, 1e-11, (fr.name + ": radial delta binds at delta").c_str());
            ConeDevClass mc = ConeDevClass::Unhandled;
            const double surfSup = circleOnConeMax(fr.cone, off, &mc);
            checkNear(surfSup, delta * std::cos(A), 1e-11,
                      (fr.name + ": surface supremum is delta*cos(A)").c_str());
            check(sup > surfSup + 1e-9,
                  (fr.name + ": bind supremum strictly exceeds the surface one").c_str());
            const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, pc);
            checkNear(g.maxGap, sup, 1e-11, (fr.name + ": oracle max == returned").c_str());
            checkNear(g.minGap, g.maxGap, 1e-11,
                      (fr.name + ": the certificate's constancy holds").c_str());

            // axial: same radius, moved delta along the axis. The pcurve pins v
            // to the NEW height, so the surface names a circle of radius
            // rho + delta*tan(A) and the gap is delta*|tan(A)|.
            const gp_Ax2 moved(rim.Location().Translated(gp_Vec(fr.cone.Position().Direction()) *
                                                         delta),
                               rim.Axis().Direction(), rim.XAxis().Direction());
            const gp_Circ ax(moved, rim.Radius());
            Handle(Geom_Curve) c3a = new Geom_Circle(ax);
            Handle(Geom2d_Curve) pca = conePCurveForCircle(fr.cone, ax);
            const char* clsA = nullptr;
            const double supA = coneCircleBindSup(c3a, 0.0, 2.0 * kPi, S, pca, id, &clsA);
            checkStr(clsA, "circle-on-cone", (fr.name + ": axial-offset class").c_str());
            checkNear(supA, delta * std::fabs(std::tan(A)), 1e-11,
                      (fr.name + ": axial delta binds at delta*|tan A|").c_str());
            const GridStat ga = sampleGap(c3a, 0.0, 2.0 * kPi, S, pca);
            checkNear(ga.maxGap, supA, 1e-11, (fr.name + ": oracle max == returned (axial)").c_str());
        }
    }
}

// ---------------------------------------------------------------------------
// 3. the certificate is invariant under rigid placement (the defect
//    130-CONE-MATH 3.1 found in the surface term must not reappear here)
// ---------------------------------------------------------------------------

static void testPlacementInvariance() {
    const TopLoc_Location id;
    std::uint64_t st = 0x5851f42d4c957f2dULL;
    auto rnd = [&]() {
        st = st * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double)((st >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
    };
    int bad = 0, badVal = 0;
    double worst = 0.0;
    const int N = 400;
    for (int i = 0; i < N; ++i) {
        const gp_Dir zd(2.0 * rnd() - 1.0 + 1e-9, 2.0 * rnd() - 1.0, 2.0 * rnd() - 1.0 + 1e-9);
        gp_Dir xd(2.0 * rnd() - 1.0, 2.0 * rnd() - 1.0 + 1e-9, 2.0 * rnd() - 1.0);
        const gp_Vec proj = gp_Vec(xd) - gp_Vec(zd) * gp_Vec(xd).Dot(gp_Vec(zd));
        if (proj.Magnitude() < 1e-6) continue;
        xd = gp_Dir(proj);
        const gp_Ax3 pos(gp_Pnt(20.0 * rnd() - 10.0, 20.0 * rnd() - 10.0, 20.0 * rnd() - 10.0), zd,
                         xd);
        gp_Cone c;
        if (!coneFromRims(pos, 8.5, 9.5, 1.0, c)) continue;
        double vLo = 0, vHi = 0;
        if (!coneFrustumVRangeFromRadii(c, 8.5, 9.5, vLo, vHi)) continue;
        Handle(Geom_Surface) S = surfOf(c);
        const gp_Circ rim = rimAt(c, vHi);
        Handle(Geom_Curve) c3 = new Geom_Circle(rim);
        Handle(Geom2d_Curve) pc = conePCurveForCircle(c, rim);
        const char* cls = nullptr;
        const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, id, &cls);
        if (!cls || std::string(cls) != "circle-on-cone") ++bad;
        if (!(sup >= 0.0 && sup < 1e-10)) ++badVal;
        worst = std::max(worst, sup);
    }
    std::fprintf(stderr, "placement invariance: %d placements, %d not circle-on-cone, %d not 0, worst %.3g\n",
                 N, bad, badVal, worst);
    check(bad == 0, "placement invariance: every exact rim classes circle-on-cone");
    check(badVal == 0, "placement invariance: every exact rim binds at 0");
}

// ---------------------------------------------------------------------------
// 4. refusals — and each refusal is WARRANTED: the oracle shows the gap really
//    is non-constant, so the endpoint value really would not be the supremum
// ---------------------------------------------------------------------------

static void testRefusals() {
    const TopLoc_Location id;
    const Frustum fr = mkFrustum("B1", 45.000004, 8.5, 9.5);
    Handle(Geom_Surface) S = surfOf(fr.cone);
    const gp_Ax3& pos = fr.cone.Position();
    const gp_Circ rim = rimAt(fr.cone, fr.vHi);

    // (a) tilted rim: axis 0.5 deg off the cone axis.
    {
        const gp_Vec tilt(gp_Vec(pos.Direction()) * std::cos(0.5 * kPi / 180.0) +
                          gp_Vec(pos.XDirection()) * std::sin(0.5 * kPi / 180.0));
        const gp_Circ tilted(gp_Ax2(rim.Location(), gp_Dir(tilt), pos.YDirection()), rim.Radius());
        Handle(Geom_Curve) c3 = new Geom_Circle(tilted);
        Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, tilted);
        check(pc.IsNull(), "tilted rim: no pcurve is offered");
        // bind against the untilted rim's pcurve to prove the refusal is not vacuous
        Handle(Geom2d_Curve) pcRim = conePCurveForCircle(fr.cone, rim);
        const char* cls = nullptr;
        const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pcRim, id, &cls);
        checkStr(cls, "unhandled-circle-on-cone", "tilted rim: refused");
        checkNear(sup, -1.0, 0.0, "tilted rim: returns -1");
        const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, pcRim);
        check(g.maxGap - g.minGap > 1e-6, "tilted rim: the gap really is non-constant");
    }

    // (b) off-axis centre, plane still perpendicular. Two sub-cases, because
    //     they say different things about the certificate:
    //     (b1) radius != the cone radius at that height -- the gap genuinely
    //          sweeps, so certificate (3) is load-bearing;
    //     (b2) radius == the cone radius -- the offset happens to be constant,
    //          and the class is refused anyway. D-130-2 certifies what is
    //          PROVED, not what happens to be true for one placement; a
    //          conservative refusal is tier 2 (counted), never a wrong number.
    {
        const gp_Pnt c2 = rim.Location().Translated(gp_Vec(pos.XDirection()) * 0.75);
        const gp_Circ offAxis(gp_Ax2(c2, pos.Direction(), pos.XDirection()), rim.Radius() + 0.3);
        Handle(Geom_Curve) c3 = new Geom_Circle(offAxis);
        Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, offAxis);
        const char* cls = nullptr;
        const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, id, &cls);
        checkStr(cls, "unhandled-circle-on-cone", "off-axis rim: refused");
        checkNear(sup, -1.0, 0.0, "off-axis rim: returns -1");
        const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, pc);
        check(g.maxGap - g.minGap > 1e-6, "off-axis rim: the gap really is non-constant");

        const gp_Circ same(gp_Ax2(c2, pos.Direction(), pos.XDirection()), rim.Radius());
        Handle(Geom_Curve) c3s = new Geom_Circle(same);
        Handle(Geom2d_Curve) pcs = conePCurveForCircle(fr.cone, same);
        const char* clsS = nullptr;
        checkNear(coneCircleBindSup(c3s, 0.0, 2.0 * kPi, S, pcs, id, &clsS), -1.0, 0.0,
                  "off-axis rim, equal radius: still returns -1");
        checkStr(clsS, "unhandled-circle-on-cone", "off-axis rim, equal radius: still refused");
        const GridStat gs = sampleGap(c3s, 0.0, 2.0 * kPi, S, pcs);
        checkNear(gs.maxGap, 0.75, 1e-9,
                  "off-axis rim, equal radius: the gap is the offset (refusal is conservative)");
    }

    // (c) reversed traversal sense: a pcurve that sweeps u the wrong way. The
    //     endpoints still agree; the interior does not.
    {
        Handle(Geom2d_Curve) good = conePCurveForCircle(fr.cone, rim);
        Handle(Geom2d_Line) gl = Handle(Geom2d_Line)::DownCast(good);
        check(!gl.IsNull(), "reversed pcurve: base pcurve is a 2d line");
        if (!gl.IsNull()) {
            const gp_Pnt2d p0 = gl->Lin2d().Location();
            const gp_Dir2d d0 = gl->Lin2d().Direction();
            Handle(Geom2d_Curve) bad = new Geom2d_Line(p0, gp_Dir2d(-d0.X(), -d0.Y()));
            Handle(Geom_Curve) c3 = new Geom_Circle(rim);
            const char* cls = nullptr;
            const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, bad, id, &cls);
            checkStr(cls, "unhandled-circle-on-cone", "reversed pcurve: refused");
            checkNear(sup, -1.0, 0.0, "reversed pcurve: returns -1");
            const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, bad);
            check(g.maxGap - g.minGap > 1.0,
                  "reversed pcurve: the gap sweeps the whole diameter");
        }
    }

    // (d) a pcurve that also moves in v (a projected / reparametrised curve).
    {
        Handle(Geom_Curve) c3 = new Geom_Circle(rim);
        Handle(Geom2d_Curve) bad = new Geom2d_Line(gp_Pnt2d(0.0, fr.vHi), gp_Dir2d(1.0, 0.05));
        const char* cls = nullptr;
        const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, bad, id, &cls);
        checkStr(cls, "unhandled-circle-on-cone", "v-varying pcurve: refused");
        checkNear(sup, -1.0, 0.0, "v-varying pcurve: returns -1");
    }

    // (e) a cylinder is not a cone: coneBindSup must refuse, not guess.
    {
        Handle(Geom_Curve) c3 = new Geom_Circle(rim);
        Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, rim);
        Handle(Geom_Surface) notCone;
        const char* cls = nullptr;
        checkNear(coneBindSup(c3, 0.0, 2.0 * kPi, notCone, pc, id, &cls), -1.0, 0.0,
                  "null surface: returns -1");
        checkStr(cls, "unhandled-circle-on-cone", "null surface: refused");
    }

    // (f) a curve that is neither circle nor line.
    {
        Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, rim);
        Handle(Geom_Curve) none;
        const char* cls = nullptr;
        checkNear(coneBindSup(none, 0.0, 1.0, S, pc, id, &cls), -1.0, 0.0,
                  "null curve: returns -1");
        checkStr(cls, "unhandled-curve-on-cone", "null curve: refused");
    }
}

// ---------------------------------------------------------------------------
// 5. the seam generator
// ---------------------------------------------------------------------------

static void testSeamLine() {
    const TopLoc_Location id;
    for (const Frustum& fr : battery()) {
        Handle(Geom_Surface) S = surfOf(fr.cone);
        const gp_Ax3& pos = fr.cone.Position();
        const double A = fr.cone.SemiAngle();
        // the generator at u = 0, unit speed in v, starting at v = vLo
        const gp_Vec rhoHat(pos.XDirection());
        const gp_Vec gen = rhoHat * std::sin(A) + gp_Vec(pos.Direction()) * std::cos(A);
        const gp_Pnt p0 = pos.Location()
                              .Translated(rhoHat * (fr.cone.RefRadius() + fr.vLo * std::sin(A)))
                              .Translated(gp_Vec(pos.Direction()) * (fr.vLo * std::cos(A)));
        const gp_Lin lin(p0, gp_Dir(gen));
        Handle(Geom_Curve) c3 = new Geom_Line(lin);
        const double len = fr.vHi - fr.vLo;
        Handle(Geom2d_Curve) pc = conePCurveForLine(fr.cone, lin, 0.5 * len);
        if (pc.IsNull()) {
            check(false, (fr.name + ": seam pcurve built").c_str());
            continue;
        }
        const char* cls = nullptr;
        const double sup = coneLineBindSup(c3, 0.0, len, S, pc, id, &cls);
        checkStr(cls, "line-on-cone", (fr.name + ": seam class").c_str());
        checkNear(sup, 0.0, 1e-11, (fr.name + ": exact generator binds at 0").c_str());
        const GridStat g = sampleGap(c3, 0.0, len, S, pc);
        checkNear(g.maxGap, 0.0, 1e-11, (fr.name + ": oracle agrees the seam is on").c_str());

        // A line displaced off the generator is affine, so the supremum really
        // is at an endpoint — assert against the dense grid, not by assumption.
        const gp_Lin moved(p0.Translated(rhoHat * 0.01), gp_Dir(gen));
        Handle(Geom_Curve) c3m = new Geom_Line(moved);
        const char* clsM = nullptr;
        const double supM = coneLineBindSup(c3m, 0.0, len, S, pc, id, &clsM);
        checkStr(clsM, "line-on-cone", (fr.name + ": displaced seam class").c_str());
        const GridStat gm = sampleGap(c3m, 0.0, len, S, pc);
        check(supM >= gm.maxGap - 1e-11,
              (fr.name + ": displaced seam supremum dominates the grid").c_str());
        checkNear(supM, gm.maxGap, 1e-9,
                  (fr.name + ": displaced seam supremum IS the grid max").c_str());

        // A line that is not a generator gets no pcurve.
        const gp_Lin skew(p0, pos.Direction());
        check(conePCurveForLine(fr.cone, skew, 0.5 * len).IsNull(),
              (fr.name + ": a non-generator line gets no pcurve").c_str());
    }
}

// ---------------------------------------------------------------------------
// 6. a trimmed surface and a trimmed curve must behave as their bases do
// ---------------------------------------------------------------------------

static void testTrimmedForms() {
    const TopLoc_Location id;
    const Frustum fr = mkFrustum("B1", 45.000004, 8.5, 9.5);
    Handle(Geom_ConicalSurface) cs = new Geom_ConicalSurface(fr.cone);
    Handle(Geom_Surface) trimmed =
        new Geom_RectangularTrimmedSurface(cs, 0.0, 2.0 * kPi, fr.vLo, fr.vHi);
    const gp_Circ rim = rimAt(fr.cone, fr.vHi);
    Handle(Geom_Curve) c3 = new Geom_TrimmedCurve(new Geom_Circle(rim), 0.0, 2.0 * kPi);
    Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, rim);
    Handle(Geom2d_Curve) pcT = new Geom2d_TrimmedCurve(pc, 0.0, 2.0 * kPi);
    const char* cls = nullptr;
    const double sup = coneBindSup(c3, 0.0, 2.0 * kPi, trimmed, pcT, id, &cls);
    checkStr(cls, "circle-on-cone", "trimmed surface + trimmed curves: class");
    checkNear(sup, 0.0, 1e-11, "trimmed surface + trimmed curves: binds at 0");
}

// ---------------------------------------------------------------------------
// 7. a placed (TopLoc_Location) face: the cone must be taken to world space
// ---------------------------------------------------------------------------

static void testLocation() {
    const Frustum fr = mkFrustum("B1", 45.000004, 8.5, 9.5);
    Handle(Geom_Surface) S = surfOf(fr.cone);
    gp_Trsf T;
    T.SetTranslation(gp_Vec(1.5, -2.25, 3.75));
    const TopLoc_Location loc(T);
    gp_Cone world = fr.cone;
    world.Transform(T);
    const gp_Circ rim = rimAt(world, fr.vHi);
    Handle(Geom_Curve) c3 = new Geom_Circle(rim);
    // The pcurve lives in the surface's OWN parameter space, so it is built from
    // the untransformed cone against the untransformed rim.
    Handle(Geom2d_Curve) pc = conePCurveForCircle(fr.cone, rimAt(fr.cone, fr.vHi));
    const char* cls = nullptr;
    const double sup = coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, loc, &cls);
    checkStr(cls, "circle-on-cone", "located face: class");
    checkNear(sup, 0.0, 1e-11, "located face: binds at 0");
}

// ---------------------------------------------------------------------------
// 8. BOTH TAPER SIGNS IN ONE FRAME  (130-CONE-ARCS)
//
// A chamfer frustum's radius may grow either way along the axis its neighbour
// cylinder carries (that axis has been through canonicalizeDir and says nothing
// about the taper). Detector C carries that sign in a SIGNED HEIGHT and never
// flips the axis. This group asserts both halves of that: the signed height
// describes either taper exactly, and flipping the axis instead does not.
// ---------------------------------------------------------------------------

static void testTaperSigns() {
    const gp_Ax3 P = obliquePlacement();  // Location = the R_lo rim, Direction = the cylinder's
    const double rLo = 8.5, rHi = 9.5;    // the plate's chamfer

    for (int sgn = 0; sgn < 2; ++sgn) {
        // +1: the R_hi rim is one unit along +Direction (the plate's UPPER
        // chamfer). -1: it is one unit along -Direction (the LOWER chamfer,
        // whose taper runs against the canonical axis).
        const double h = (sgn == 0) ? 1.0 : -1.0;
        const std::string nm = (sgn == 0) ? "taper +h" : "taper -h";
        gp_Cone cone;
        const bool built = coneFromRims(P, rLo, rHi, h, cone);
        check(built, (nm + ": coneFromRims builds").c_str());
        if (!built) continue;
        checkNear(cone.SemiAngle(), (h > 0.0 ? 1.0 : -1.0) * kPi / 4.0, 1e-12,
                  (nm + ": SemiAngle carries the sign").c_str());
        checkNear(cone.RefRadius(), rLo, 1e-12, (nm + ": RefRadius is R_lo").c_str());
        // ONE FRAME: Direction, XDirection and YDirection are the neighbour
        // cylinder's for both signs, so every circle either side builds in it.
        check(cone.Position().Direction().IsEqual(P.Direction(), 1e-12) &&
                  cone.Position().XDirection().IsEqual(P.XDirection(), 1e-12) &&
                  cone.Position().YDirection().IsEqual(P.YDirection(), 1e-12),
              (nm + ": frame is the neighbour cylinder's, unchanged").c_str());

        for (int k = 0; k < 2; ++k) {
            // Exactly what coneIsoCircle builds: the axial offset is a fraction
            // of the SIGNED height and the circle takes the region's own
            // Direction and XDirection.
            const double vAx = (k == 0) ? 0.0 : h;
            const double R = rLo + (vAx / h) * (rHi - rLo);
            const std::string rn = nm + (k == 0 ? " R_lo rim" : " R_hi rim");
            checkNear(R, k == 0 ? rLo : rHi, 1e-12, (rn + ": radius").c_str());
            const gp_Pnt ctr = P.Location().Translated(gp_Vec(P.Direction()) * vAx);
            const gp_Circ rim(gp_Ax2(ctr, P.Direction(), P.XDirection()), R);
            // 130-CONE-MATH: the rim is ON the surface, both signs.
            ConeDevClass mc = ConeDevClass::Unhandled;
            checkNear(circleOnConeMax(cone, rim, &mc), 0.0, 1e-12,
                      (rn + ": circleOnConeMax = 0").c_str());
            // and it binds at 0 through the certified class.
            Handle(Geom_Surface) S = surfOf(cone);
            Handle(Geom2d_Curve) pc = conePCurveForCircle(cone, rim);
            check(!pc.IsNull(), (rn + ": certified pcurve exists").c_str());
            if (pc.IsNull()) continue;
            Handle(Geom_Curve) c3 = new Geom_Circle(rim);
            const char* cls = nullptr;
            const double sup =
                coneCircleBindSup(c3, 0.0, 2.0 * kPi, S, pc, TopLoc_Location(), &cls);
            checkStr(cls, "circle-on-cone", (rn + ": class").c_str());
            checkNear(sup, 0.0, 1e-11, (rn + ": binds at 0").c_str());
            const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, pc);
            checkNear(g.maxGap, 0.0, 1e-11, (rn + ": oracle grid max = 0").c_str());
            // u sweeps FORWARD in this frame for both taper signs -- which is
            // what lets the face builder write a slope +1 rim pcurve at all.
            const gp_Pnt2d a2 = pc->Value(0.0), b2 = pc->Value(1.0);
            check(b2.X() - a2.X() > 0.0, (rn + ": u sweeps forward").c_str());
        }

        // The seam generator runs from the R_lo rim TOWARD the R_hi rim. Built
        // from the frustum's own two numbers it is right for either sign;
        // sin(Ang)*X + cos(Ang)*Z is the surface's +v and points AWAY from the
        // R_hi rim when Ang < 0.
        const gp_Dir gen(gp_Vec(P.XDirection()) * (rHi - rLo) + gp_Vec(P.Direction()) * h);
        const gp_Pnt pLo = P.Location().Translated(gp_Vec(P.XDirection()) * rLo);
        const gp_Pnt pHi = P.Location()
                               .Translated(gp_Vec(P.Direction()) * h)
                               .Translated(gp_Vec(P.XDirection()) * rHi);
        const double len = std::fabs(coneVAtRadius(cone, rHi));
        checkNear(pLo.Translated(gp_Vec(gen) * len).Distance(pHi), 0.0, 1e-12,
                  (nm + ": generator reaches the R_hi seam vertex").c_str());
    }

    // The defect the signed height replaces. Flipping Direction and keeping
    // XDirection gives a RIGHT-handed frame with YDirection negated, so the
    // cone's u runs the other way round the rim circle the neighbour cylinder
    // built -- and a slope +1 rim pcurve there is off by the full diameter.
    {
        const gp_Ax3 F(P.Location(), P.Direction().Reversed(), P.XDirection());
        gp_Cone flipped;
        check(coneFromRims(F, rLo, rHi, 1.0, flipped), "flipped axis: coneFromRims builds");
        check(flipped.Position().YDirection().IsOpposite(P.YDirection(), 1e-9),
              "flipped axis: YDirection is negated (gp_Ax3 is right-handed)");
        const gp_Circ rim(gp_Ax2(P.Location(), P.Direction(), P.XDirection()), rLo);
        Handle(Geom_Surface) S = surfOf(flipped);
        Handle(Geom2d_Curve) pc = conePCurveForCircle(flipped, rim);
        check(!pc.IsNull(), "flipped axis: certified pcurve exists");
        if (!pc.IsNull()) {
            const gp_Pnt2d a2 = pc->Value(0.0), b2 = pc->Value(1.0);
            check(b2.X() - a2.X() < 0.0,
                  "flipped axis: the certified pcurve must sweep u BACKWARD");
        }
        Handle(Geom_Curve) c3 = new Geom_Circle(rim);
        Handle(Geom2d_Curve) slopePlus1 =
            new Geom2d_Line(gp_Pnt2d(0.0, coneVAtZ(flipped, 0.0)), gp_Dir2d(1.0, 0.0));
        const GridStat g = sampleGap(c3, 0.0, 2.0 * kPi, S, slopePlus1);
        checkNear(g.maxGap, 2.0 * rLo, 1e-9,
                  "flipped axis + slope +1 rim pcurve: gap is the full diameter");
    }
}

int main() {
    testExactRim();
    testParametrisationTerm();
    testPlacementInvariance();
    testRefusals();
    testSeamLine();
    testTrimmedForms();
    testLocation();
    testTaperSigns();
    std::fprintf(stderr, "cone_bind_unit: %d/%d PASS\n", gPass, gPass + gFail);
    return gFail ? 1 : 0;
}
