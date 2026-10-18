// cyl|cyl contact predicates (D-fc-5.17). Pure: no mesh, no edge.
// A Steinmetz section is the plane cut of cylinder A by a bisector of the two
// axes. It is admissible only when the axes are not parallel and
// |Ra-Rb| + axisSeparation <= τ, τ = tauOf = 2q. The published chord is a
// separate certificate: polylineOnCylSup is the same supremum exactMaxAtBind
// records for a degree-1 span on a cylinder.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_CYLCYL_CONTACT_HPP
#define STL2STEP_REFIT_CYLCYL_CONTACT_HPP

#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>

#include <algorithm>
#include <cmath>

namespace stl2step {
namespace refit {

struct CylAxis {
    gp_Ax3 ax;
    double radius = 0;
};

struct AxisSep {
    double distance = 0;
    gp_Pnt onA;
    gp_Pnt onB;
    bool parallel = false;
};

struct SteinmetzEllipse {
    gp_Elips elips;
};

struct SteinmetzPair {
    int n = 0;
    SteinmetzEllipse e0;
    SteinmetzEllipse e1;
};

inline double cylcylPi() { return std::acos(-1.0); }

// Shortest distance between the two axis lines, and a closest point on each.
inline AxisSep cylinderAxisSeparation(const gp_Ax3& a, const gp_Ax3& b) {
    const gp_Vec d1(a.Direction());
    const gp_Vec d2(b.Direction());
    const gp_Vec w(a.Location(), b.Location());
    const double c = d1.Dot(d2);
    const double s2 = std::max(0.0, 1.0 - c * c);
    const double ang = Precision::Angular();
    AxisSep o;
    o.parallel = !(s2 > ang * ang);
    if (o.parallel) {
        const double t = w.Dot(d1);
        o.onA = a.Location().Translated(d1 * t);
        o.onB = b.Location();
        o.distance = d1.Crossed(gp_Vec(o.onA, o.onB)).Magnitude();
        return o;
    }
    const gp_Vec n = d1.Crossed(d2);
    const double n2 = n.SquareMagnitude();
    o.distance = std::fabs(w.Dot(n)) / std::sqrt(n2);
    const double t = w.Crossed(d2).Dot(n) / n2;
    const double s = w.Crossed(d1).Dot(n) / n2;
    o.onA = a.Location().Translated(d1 * t);
    o.onB = b.Location().Translated(d2 * s);
    return o;
}

// |Ra-Rb| + axisSeparation. The 1-Lipschitz budget that admits a section of A
// as a contact candidate. Parallel axes are not a Steinmetz pair.
inline double steinmetzGap(double radiusA, double radiusB, const AxisSep& sep) {
    return std::fabs(radiusA - radiusB) + sep.distance;
}

inline bool steinmetzAdmissible(double radiusA, double radiusB, const AxisSep& sep, double tau) {
    if (!(tau > 0.0) || sep.parallel) return false;
    if (!(radiusA > 0.0) || !(radiusB > 0.0)) return false;
    if (!std::isfinite(sep.distance)) return false;
    return steinmetzGap(radiusA, radiusB, sep) <= tau;
}

inline double distToCylinder(const gp_Ax3& ax, double radius, const gp_Pnt& p) {
    const gp_Vec rho(ax.Location(), p);
    const gp_Vec axial = gp_Vec(ax.Direction()) * rho.Dot(ax.Direction());
    return std::fabs((rho - axial).Magnitude() - radius);
}

inline bool withinTauOfBoth(const gp_Pnt& p, const CylAxis& a, const CylAxis& b, double tau) {
    if (!(tau > 0.0)) return false;
    return distToCylinder(a.ax, a.radius, p) <= tau && distToCylinder(b.ax, b.radius, p) <= tau;
}

// g(h) from exactMaxAtBind's polyline-on-cylinder term. h is the half azimuth.
inline double polylineAzimuthGap(double h) {
    const double halfPi = 0.5 * cylcylPi();
    if (!(h > 0.0) || !(h < halfPi)) return 0.0;
    const double th = std::tan(h);
    if (!(th > 0.0)) return 0.0;
    const double q2 = th / h - 1.0;
    const double sStar = q2 > 0.0 ? std::sqrt(q2) / th : 0.0;
    return std::fabs(std::atan(th * sStar) - h * sStar);
}

// Distance from m to the segment p0–p1. A closed intersection meets itself,
// so endpoint distance alone does not prove the parameter span is the chord.
inline double segmentGap(const gp_Pnt& p0, const gp_Pnt& p1, const gp_Pnt& m) {
    const gp_Vec ab(p0, p1);
    const double ab2 = ab.SquareMagnitude();
    if (!(ab2 > gp::Resolution() * gp::Resolution())) return m.Distance(p0);
    const double u = std::max(0.0, std::min(1.0, gp_Vec(p0, m).Dot(ab) / ab2));
    return m.Distance(p0.Translated(ab * u));
}

// |chord - cylinder| supremum for one degree-1 span. Same expression as
// exactMaxAtBind (endpoint radial error + R(1-cos h) + R g(h)).
inline double polylineOnCylSup(const gp_Pnt& p0, const gp_Pnt& p1, const gp_Ax3& pos, double radius) {
    const gp_Vec d(pos.Direction());
    auto radial = [&](const gp_Pnt& P) {
        const gp_Vec rho(pos.Location(), P);
        return rho - d * rho.Dot(d);
    };
    const gp_Vec w0 = radial(p0);
    const gp_Vec w1 = radial(p1);
    const double r0 = w0.Magnitude();
    const double r1 = w1.Magnitude();
    const double e = std::max(std::fabs(r0 - radius), std::fabs(r1 - radius));
    double h = 0.0;
    if (r0 > gp::Resolution() && r1 > gp::Resolution()) {
        const double cs = std::max(-1.0, std::min(1.0, w0.Dot(w1) / (r0 * r1)));
        h = 0.5 * std::acos(cs);
    }
    return e + radius * (1.0 - std::cos(h)) + radius * polylineAzimuthGap(h);
}

// How far a contact of this chain may sit from the mesh vertices. Each mesh
// edge already leaves both cylinders by polylineOnCylSup; the sum is that
// edge's certificate. The chain's band is the worst edge. A branch of the
// infinite-cylinder intersection farther than the band is not this chain.
inline double contactChordBand(const gp_Pnt* pts, int n, bool closed, const gp_Ax3& axA, double rA,
                               const gp_Ax3& axB, double rB) {
    if (pts == nullptr || n < 2) return 0.0;
    double band = 0.0;
    auto accum = [&](const gp_Pnt& a, const gp_Pnt& b) {
        const double sA = polylineOnCylSup(a, b, axA, rA);
        const double sB = polylineOnCylSup(a, b, axB, rB);
        if (std::isfinite(sA) && std::isfinite(sB)) band = std::max(band, sA + sB);
    };
    for (int i = 1; i < n; ++i) accum(pts[i - 1], pts[i]);
    if (closed) accum(pts[n - 1], pts[0]);
    return band;
}

// Largest half-azimuth whose sag term R(1-cos h)+R g(h) is within budget.
// Zero when no positive h fits.
inline double polylineHalfAzimuth(double radius, double budget) {
    if (!(radius > 0.0) || !(budget > 0.0)) return 0.0;
    const double hiLimit = 0.5 * cylcylPi() - Precision::Angular();
    auto sag = [&](double h) {
        return radius * (1.0 - std::cos(h)) + radius * polylineAzimuthGap(h);
    };
    if (sag(hiLimit) <= budget) return hiLimit;
    double lo = 0.0;
    double hi = hiLimit;
    for (;;) {
        const double mid = 0.5 * (lo + hi);
        if (!(mid > lo) || !(mid < hi)) break;
        if (sag(mid) <= budget) lo = mid;
        else hi = mid;
    }
    return lo;
}

// Plane through onA (so the section lies on cylinder A) with normal N.
inline bool steinmetzSection(const gp_Ax3& axisA, double radiusA, const gp_Dir& normal, gp_Elips& out) {
    const gp_Vec d(axisA.Direction());
    const double dn = std::fabs(d.Dot(normal));
    if (!(dn > std::sin(Precision::Angular()))) return false;
    const double major = radiusA / dn;
    const double minor = radiusA;
    if (!(major >= minor)) return false;
    gp_Vec maj = d - gp_Vec(normal) * d.Dot(normal);
    if (!(maj.Magnitude() > gp::Resolution())) return false;
    maj.Normalize();
    try {
        out = gp_Elips(gp_Ax2(axisA.Location(), normal, gp_Dir(maj)), major, minor);
    } catch (const Standard_Failure&) {
        return false;
    }
    return major > minor;
}

inline SteinmetzPair steinmetzSections(const CylAxis& a, const CylAxis& b, double tau) {
    SteinmetzPair o;
    const AxisSep sep = cylinderAxisSeparation(a.ax, b.ax);
    if (!steinmetzAdmissible(a.radius, b.radius, sep, tau)) return o;
    gp_Ax3 axisA = a.ax;
    axisA.SetLocation(sep.onA);
    const gp_Vec d1(axisA.Direction());
    const gp_Vec d2(b.ax.Direction());
    const gp_Vec arms[2] = {d1 - d2, d1 + d2};
    for (const gp_Vec& arm : arms) {
        if (!(arm.Magnitude() > gp::Resolution())) continue;
        gp_Elips el;
        if (!steinmetzSection(axisA, a.radius, gp_Dir(arm), el)) continue;
        if (o.n == 0) o.e0.elips = el;
        else o.e1.elips = el;
        o.n++;
        if (o.n == 2) break;
    }
    return o;
}

}  // namespace refit
}  // namespace stl2step

#endif
