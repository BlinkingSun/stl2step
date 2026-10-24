// Synthetic sphere-corner cases (D-fc-5.11 F8): the predicates of
// refit_relational_sphere.hpp and claimSphereCorners on a box corner of three
// quarter-cylinder rounds. No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_internal.hpp"
#include "refit_relational.hpp"
#include "refit_relational_sphere.hpp"

#include <ElSLib.hxx>
#include <Precision.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kQ = 1e-6;
static constexpr double kTau = 2 * kQ;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;

    int add(const gp_XYZ& p) {
        pts.push_back(p);
        return (int)pts.size() - 1;
    }
    int addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
        return (int)(triFlat.size() / 3) - 1;
    }
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2],
                              triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]);
                const int hi = std::max(v[s], v[s + 1]);
                auto it = ids.find({lo, hi});
                int id;
                if (it == ids.end()) {
                    id = (int)edges.size();
                    ids.emplace(std::make_pair(lo, hi), id);
                    edges.emplace_back(lo, hi);
                } else {
                    id = it->second;
                }
                triEdges[t][(std::size_t)s] = id;
                if (v[s] < v[s + 1]) triDirs[t] = (uint8_t)(triDirs[t] | (uint8_t)(1u << s));
            }
        }
        compVtx.resize(pts.size());
        for (std::size_t i = 0; i < pts.size(); ++i) compVtx[i] = (int)i;
    }
    MeshView view() const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int (*)[3]>(triFlat.data());
        mv.compTris = compTris.data();
        mv.compVtx = compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = triEdges.data();
        mv.triDirs = triDirs.data();
        mv.nTri = triFlat.size() / 3;
        mv.nVtx = pts.size();
        mv.nEdge = edges.size();
        mv.diag = 100.0;
        mv.quantFloor = kQ;
        return mv;
    }
};

// --- predicates ---------------------------------------------------------

static void casePredicates() {
    RoundSpine line;
    line.o = gp_XYZ(1, 2, 3);
    line.a = gp_XYZ(0, 0, 1);
    line.r = 2;
    expect(std::fabs(spineDistance(line, gp_XYZ(4, 6, -7)) - 5.0) < 1e-12, "line spine distance");
    RoundSpine circ;
    circ.circle = true;
    circ.o = gp_XYZ(0, 0, 0);
    circ.a = gp_XYZ(0, 0, 1);
    circ.R = 10;
    circ.r = 1;
    expect(std::fabs(spineDistance(circ, gp_XYZ(13, 0, 4)) - 5.0) < 1e-12,
           "circle spine distance");
    gp_XYZ t;
    expect(spineTangentAt(circ, gp_XYZ(10, 0, 0), t) && std::fabs(t.Y() - 1.0) < 1e-12,
           "circle spine tangent");
    expect(!spineTangentAt(circ, gp_XYZ(0, 0, 5), t), "no tangent on the circle spine's axis");

    // Two lines meeting at (1, 1, 1); a parallel pair has no point.
    RoundSpine lx, ly, lz;
    lx.o = gp_XYZ(-5, 1, 1);
    lx.a = gp_XYZ(1, 0, 0);
    ly.o = gp_XYZ(1, 7, 1);
    ly.a = gp_XYZ(0, 1, 0);
    lz.o = gp_XYZ(1, 1, 4);
    lz.a = gp_XYZ(0, 0, 1);
    gp_XYZ c;
    expect(lineSpinesPoint({lx, ly}, c) && (c - gp_XYZ(1, 1, 1)).Modulus() < 1e-12,
           "two meeting lines give their point");
    RoundSpine lx2 = lx;
    lx2.o = gp_XYZ(0, 3, 1);
    expect(!lineSpinesPoint({lx, lx2}, c), "parallel lines give no point");
    expect(!lineSpinesPoint({lx}, c), "one line gives no point");

    // sphereOfRounds: the anchor centre and the anchor radius. The shipped
    // radius is fitted later, on the corner's vertices.
    lx.r = ly.r = lz.r = 2.0;
    SphereCornerFit f = sphereOfRounds({lx, ly, lz}, kTau);
    expect(f.ok && std::strcmp(f.why, "accept") == 0 && std::fabs(f.r - 2.0) < 1e-12 &&
               (f.c - gp_XYZ(1, 1, 1)).Modulus() < 1e-12,
           "three rounds of one radius meeting: accept");
    lz.r = 2.0 + 1.5 * kTau;
    f = sphereOfRounds({lx, ly, lz}, kTau);
    expect(!f.ok && std::strcmp(f.why, "not-blend") == 0, "radii apart by more than tau: refuse");
    lz.r = 2.0 + 0.5 * kTau;
    f = sphereOfRounds({lx, ly, lz}, kTau);
    expect(f.ok, "radii within tau: accept");
    lz.r = 2.0;
    // Least squares over the three lines splits a 3 tau offset: each spine
    // ends 1.5 tau from the centre.
    lz.o = gp_XYZ(1, 1 + 3.0 * kTau, 4);
    f = sphereOfRounds({lx, ly, lz}, kTau);
    expect(!f.ok && std::strcmp(f.why, "not-blend") == 0, "a spine off the centre by > tau: refuse");
    lz.o = gp_XYZ(1, 1 + 1.5 * kTau, 4);
    f = sphereOfRounds({lx, ly, lz}, kTau);
    expect(f.ok && f.spineDev <= kTau, "a 1.5 tau offset shared by the fit: within tau");
    lz.o = gp_XYZ(1, 1, 4);
    RoundSpine c1 = circ, c2 = circ;
    f = sphereOfRounds({lx, c1, c2}, kTau);
    expect(!f.ok && std::strcmp(f.why, "underdetermined") == 0,
           "fewer than two line spines: underdetermined");
    // A torus spine through the corner counts like a line.
    RoundSpine tor;
    tor.circle = true;
    tor.o = gp_XYZ(1, 11, 1);
    tor.a = gp_XYZ(0, 0, 1);
    tor.R = 10;
    tor.r = 2;
    f = sphereOfRounds({lx, lz, tor}, kTau);
    expect(f.ok, "two lines and a torus spine through the corner: accept");

    // Least-squares radius with the centre held: the mean of the distances.
    {
        const gp_XYZ o(0, 0, 0);
        const std::vector<gp_XYZ> pts = {gp_XYZ(3, 0, 0), gp_XYZ(0, 4, 0), gp_XYZ(0, 0, 5)};
        expect(std::fabs(sphereRadiusLeastSquares(o, pts) - 4.0) < 1e-12,
               "least-squares radius is the mean distance");
        const std::vector<gp_XYZ> none;
        expect(!std::isfinite(sphereRadiusLeastSquares(o, none)), "no vertices: no radius");
        const std::vector<gp_XYZ> bad = {gp_XYZ(1, 0, 0),
                                         gp_XYZ(std::numeric_limits<double>::quiet_NaN(), 0, 0)};
        expect(!std::isfinite(sphereRadiusLeastSquares(o, bad)), "a non-finite vertex: no radius");
    }
    // Contact at tau, not a wider band. A cylinder on the axis and a torus
    // whose tube circle holds the centre.
    {
        RoundSpine cyl;
        cyl.o = gp_XYZ(1, 1, 4);
        cyl.a = gp_XYZ(0, 0, 1);
        cyl.r = 2.0;
        expect(sphereRoundContact(gp_XYZ(1, 1, 1), 2.0 + 0.5 * kTau, cyl, kTau),
               "cylinder contact within tau");
        expect(!sphereRoundContact(gp_XYZ(1, 1, 1), 2.0 + 1.5 * kTau, cyl, kTau),
               "cylinder contact beyond tau: refuse");
        RoundSpine torC;
        torC.circle = true;
        torC.o = gp_XYZ(0, 0, 0);
        torC.a = gp_XYZ(0, 0, 1);
        torC.R = 10.0;
        torC.r = 3.0;
        expect(sphereRoundContact(gp_XYZ(10, 0, 0), 3.0 + 0.5 * kTau, torC, kTau),
               "torus contact within tau");
        expect(!sphereRoundContact(gp_XYZ(10, 0, 0), 3.0 + 1.5 * kTau, torC, kTau),
               "torus contact beyond tau: refuse");
    }

    // Plane first.
    std::vector<gp_XYZ> pl = {gp_XYZ(0, 0, 0), gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0),
                              gp_XYZ(1, 1, 0.5 * kTau), gp_XYZ(2, 3, 0)};
    expect(pointsPlanar(pl, kTau), "points within tau of a plane: planar");
    expect(planeGateRefuses(pl, kTau), "a solved plane within tau refuses");
    pl.push_back(gp_XYZ(0.5, 0.5, 3 * kTau));
    expect(!pointsPlanar(pl, kTau), "a point 3 tau off: not planar");
    expect(!planeGateRefuses(pl, kTau), "a solved plane beyond tau does not refuse");
    // D-fc-5.20. Coordinates whose covariance overflows make the eigensolve
    // return false. That is a refusal, not a skip of the plane gate. The
    // magnitude is only large enough to overflow a product of doubles. It is
    // not a measured length from a part.
    const double big = 1e200;
    const std::vector<gp_XYZ> boom = {gp_XYZ(0, 0, 0), gp_XYZ(big, 0, 0), gp_XYZ(0, big, 0),
                                      gp_XYZ(0, 0, big), gp_XYZ(big, big, 0)};
    expect(!pointsPlanar(boom, kTau), "an unsolved plane is not a measured plane");
    expect(planeGateRefuses(boom, kTau), "a failed eigensolve refuses the plane gate");
    expect(planeResidual(boom) == HUGE_VAL, "a failed eigensolve is not a finite residual");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<gp_XYZ> nanPts = {gp_XYZ(0, 0, 0), gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0),
                                        gp_XYZ(0, 0, 1), gp_XYZ(nan, 0, 0)};
    expect(planeGateRefuses(nanPts, kTau), "a non-finite cloud refuses the plane gate");
    expect(planeGateRefuses(pl, 0.0), "tau 0 refuses the plane gate");
    expect(planeGateRefuses(pl, nan), "a non-finite tau refuses the plane gate");

    // pinVertices. A gap at the cap used to be stored past the cap by the pad.
    const double cap = 1.0;
    expect(pinnedVertexTol(cap, cap) == cap, "a vertex at the cap is stored at the cap");
    const double mid = pinnedVertexTol(0.5 * cap, cap);
    expect(mid > 0.5 * cap && mid < cap, "a vertex at half the cap keeps the pad inside the cap");
    const double near = cap / 1.001;
    expect(near < cap && pinnedVertexTol(near, cap) == cap,
           "a pad that would pass the cap is clamped to the cap");
    expect(pinnedVertexTol(0.0, cap) == Precision::Confusion(),
           "a vertex on the pcurve end stores only Confusion, under the cap");

    // I9 cap volume: the flat octant triangle of a unit sphere bounds pi/6 - 1/6.
    const double v = sphereCapVolume(gp_XYZ(0, 0, 0), 1.0, gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0),
                                     gp_XYZ(0, 0, 1));
    expect(std::fabs(v - (kPi / 6.0 - 1.0 / 6.0)) < 1e-12, "octant cap volume pi/6 - 1/6");

    // Iso bound: an exact meridian / parallel is 0; a circle tilted by e is ~ r e.
    const gp_XYZ X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1), O(0, 0, 0);
    SphereIsoLine L = sphereIsoLineOf(O, 2.0, X, Y, Z, O, 2.0, X, Z, 0.7);
    expect(L.meridian && L.bound < 1e-12, "great circle through the poles: exact meridian");
    L = sphereIsoLineOf(O, 2.0, X, Y, Z, O, 2.0, X, Y, 0.7);
    expect(!L.meridian && L.bound < 1e-12, "the equator: exact parallel");
    L = sphereIsoLineOf(O, 2.0, X, Y, Z, gp_XYZ(0, 0, 1), std::sqrt(3.0), X, Y, 2.0);
    expect(!L.meridian && L.bound < 1e-12 && std::fabs(L.v0 - kPi / 6.0) < 1e-12,
           "a latitude circle: exact parallel at its latitude");
    const double e = 1e-4;
    const gp_XYZ Zt(0, -std::sin(e), std::cos(e));
    L = sphereIsoLineOf(O, 2.0, X, Y, Z, O, 2.0, X, Zt, 0.7);
    expect(L.meridian && L.bound > 0.5 * 2.0 * e && L.bound < 4.0 * 2.0 * e,
           "a meridian tilted by e: bound of order r e");

    // Contact circle of a cylinder round.
    double vv = 0;
    expect(sphereCylinderContact(gp_XYZ(1, 1, 1), 2.0, gp_XYZ(1, 1, 4), Z, 2.0 + 0.5 * kTau,
                                 2 * kTau, vv) &&
               std::fabs(vv + 3.0) < 1e-12,
           "cylinder contact: the v-iso through the centre");
    expect(!sphereCylinderContact(gp_XYZ(1, 1 + 3 * kTau, 1), 2.0, gp_XYZ(1, 1, 4), Z, 2.0,
                                  2 * kTau, vv),
           "cylinder contact: centre off the axis by more than tol: none");

    // A torus spine within tanDistTol still meets; past it does not.
    tor.R = 10.0 - 1.5 * kTau;
    f = sphereOfRounds({lx, lz, tor}, kTau);
    expect(f.ok, "torus spine within tanDistTol of the centre: accept");
    tor.R = 10.0 - (tanDistTol(kTau, kTau) + kTau);
    f = sphereOfRounds({lx, lz, tor}, kTau);
    expect(!f.ok && std::strcmp(f.why, "not-blend") == 0,
           "torus spine beyond tanDistTol of the centre: refuse");
    tor.R = 10.0;

    // Torus contact: the centre lies on the major circle. A point on the tube
    // surface is one minor radius off that circle and is not the contact.
    const gp_Ax3 fr(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    double u = -1.0;
    const double tolT = tanDistTol(kTau, kTau);
    expect(sphereTorusContact(gp_XYZ(10, 0, 0), 3.0, fr, 10.0, 3.0, tolT, u) &&
               std::fabs(u) < 1e-12,
           "torus contact: centre on the tube circle at u=0");
    u = -1.0;
    const bool atHalf = sphereTorusContact(gp_XYZ(0, 10, 0), 3.0, fr, 10.0, 3.0, tolT, u) &&
                        std::fabs(u - 0.5 * kPi) < 1e-12;
    expect(atHalf, "torus contact: centre on the tube circle at u=pi/2");
    if (atHalf) {
        const gp_Circ iso = ElSLib::TorusUIso(fr, 10.0, 3.0, u);
        expect(std::fabs(iso.Radius() - 3.0) < 1e-12 &&
                   iso.Location().Distance(gp_Pnt(0, 10, 0)) < 1e-9,
               "torus contact azimuth is the u-iso through the centre");
    } else {
        expect(false, "torus contact azimuth is the u-iso through the centre");
    }
    expect(!sphereTorusContact(gp_XYZ(13, 0, 0), 3.0, fr, 10.0, 3.0, tolT, u),
           "torus contact: a tube-surface point is not the centre");
    expect(!sphereTorusContact(gp_XYZ(10, 0, 0), 3.0 + 3.0 * kTau, fr, 10.0, 3.0, tolT, u),
           "torus contact: minor radius disagrees by more than tol");
    expect(!sphereTorusContact(gp_XYZ(2, 0, 0), 3.0, fr, 2.0, 3.0, tolT, u),
           "torus contact: a spindle is not a ring");
    expect(!sphereTorusContact(gp_XYZ(10, 0, 0) + gp_XYZ(0, 0, tolT + kTau), 3.0, fr, 10.0, 3.0,
                               tolT, u),
           "torus contact: centre off the tube circle by more than tol");
}

// --- round certificate --------------------------------------------------

static void caseRoundCert() {
    RoundSpine s;
    s.o = gp_XYZ(0, 0, 0);
    s.a = gp_XYZ(0, 0, 1);
    s.r = 2.0;
    std::vector<std::array<gp_XYZ, 3>> tris;
    const int n = 6;
    for (int i = 0; i < n; ++i) {
        const double a0 = 0.5 * kPi * i / n, a1 = 0.5 * kPi * (i + 1) / n;
        const gp_XYZ p00(2 * std::cos(a0), 2 * std::sin(a0), 0), p10(2 * std::cos(a1), 2 * std::sin(a1), 0);
        const gp_XYZ p01(p00.X(), p00.Y(), 5), p11(p10.X(), p10.Y(), 5);
        tris.push_back({p00, p10, p01});
        tris.push_back({p01, p10, p11});
    }
    RoundCert rc = certifyCylinderRound(s, tris, kTau);
    expect(rc.ok && rc.nI == 2 * n && rc.nCover == 0 && rc.stations >= minStations(kParamsCircle) &&
               rc.levels == 2,
           "exact quarter cylinder certifies on its own triangles");
    // A fill triangle with one vertex on the round and the others inside the
    // chord bound s of I is cover.
    std::vector<std::array<gp_XYZ, 3>> t2 = tris;
    const gp_XYZ a(2, 0, 2), b(2 * std::cos(0.1), 2 * std::sin(0.1), 2);
    const gp_XYZ in = (a + b) * 0.5;  // a chord midpoint: s off the cylinder
    t2.push_back({a, b, in});
    rc = certifyCylinderRound(s, t2, kTau);
    expect(rc.ok && rc.nCover == 1, "a chord-bounded fill is cover");
    // A triangle far off the round is neither: not certified.
    std::vector<std::array<gp_XYZ, 3>> t3 = tris;
    t3.push_back({a, b, gp_XYZ(1.0, 0.0, 2.0)});
    rc = certifyCylinderRound(s, t3, kTau);
    expect(!rc.ok && rc.nOther == 1, "a triangle beyond the chord bound: not certified");
    // One column of triangles cannot over-determine the circle.
    std::vector<std::array<gp_XYZ, 3>> t4(tris.begin(), tris.begin() + 2);
    rc = certifyCylinderRound(s, t4, kTau);
    expect(!rc.ok, "two stations: under-determined, not certified");
}

// --- the stage on a box corner -----------------------------------------

struct Corner {
    Synth s;
    std::vector<int> sphereTris;
    std::vector<std::vector<int>> roundTris;  // X, Y, Z
};

// Sphere octant x, y, z >= 0 about the origin (radius r) by an n-level
// barycentric grid projected onto the sphere, and the three quarter-cylinder
// rounds along the axes, each ending on the octant's boundary circle and
// running length L away from the corner. Outward winding.
static Corner makeCorner(double r, int n, double L, double rZ, double bump) {
    Corner k;
    std::map<std::array<int, 3>, int> id;
    auto vid = [&](int i, int j, int m) {
        const std::array<int, 3> key = {i, j, m};
        auto it = id.find(key);
        if (it != id.end()) return it->second;
        gp_XYZ p((double)i, (double)j, (double)m);
        p /= p.Modulus();
        p *= r;
        const int v = k.s.add(p);
        id[key] = v;
        return v;
    };
    for (int i = 0; i < n; ++i)
        for (int j = 0; i + j < n; ++j) {
            const int m = n - i - j;
            k.sphereTris.push_back(k.s.addTri(vid(i + 1, j, m - 1), vid(i, j + 1, m - 1),
                                              vid(i, j, m)));
            if (i + j < n - 1)
                k.sphereTris.push_back(k.s.addTri(vid(i + 1, j, m - 1), vid(i + 1, j + 1, m - 2),
                                                  vid(i, j + 1, m - 1)));
        }
    if (bump != 0.0) {
        // Move one interior vertex radially (it must exist for n >= 3).
        const int v = vid(1, 1, n - 2);
        gp_XYZ p = k.s.pts[(std::size_t)v];
        k.s.pts[(std::size_t)v] = p * ((p.Modulus() + bump) / p.Modulus());
    }
    // Rounds: axis a, quadrant spanned by (e1, e2); the boundary ring is the
    // grid edge with the a-coordinate 0.
    auto round = [&](int axis, double rad) {
        std::vector<int> tris;
        std::vector<int> ring, far;
        for (int t = 0; t <= n; ++t) {
            std::array<int, 3> g = {0, 0, 0};
            const int e1 = (axis + 1) % 3, e2 = (axis + 2) % 3;
            g[(std::size_t)e1] = n - t;
            g[(std::size_t)e2] = t;
            const int v = vid(g[0], g[1], g[2]);
            ring.push_back(v);
            gp_XYZ p = k.s.pts[(std::size_t)v] * (rad / r);
            gp_XYZ q = p;
            q.SetCoord(axis + 1, -L);
            far.push_back(k.s.add(q));
        }
        for (int t = 0; t < n; ++t) {
            tris.push_back(k.s.addTri(ring[(std::size_t)t], far[(std::size_t)t],
                                      ring[(std::size_t)t + 1]));
            tris.push_back(k.s.addTri(ring[(std::size_t)t + 1], far[(std::size_t)t],
                                      far[(std::size_t)t + 1]));
        }
        return tris;
    };
    k.roundTris.push_back(round(0, r));
    k.roundTris.push_back(round(1, r));
    k.roundTris.push_back(round(2, rZ));
    k.s.finish();
    return k;
}

static Region roundRegion(int axis, double rad, const std::vector<int>& tris, Origin o,
                          double dev) {
    Region g;
    g.type = SurfType::Cylinder;
    g.origin = o;
    gp_Dir d(axis == 0 ? 1 : 0, axis == 1 ? 1 : 0, axis == 2 ? 1 : 0);
    gp_Dir x(axis == 0 ? 0 : 1, axis == 0 ? 1 : 0, 0);
    g.ax = gp_Ax3(gp_Pnt(0, 0, 0), d, x);
    g.radius = rad;
    g.tris = tris;
    g.maxVertexDev = dev;
    g.closed360 = false;
    return g;
}

static int runCorner(const Corner& k, double rZ, Origin o, double devZ, Region* out,
                     double* dvol) {
    SegmentWork w;
    for (int q = 0; q < 3; ++q)
        w.accepted.push_back(roundRegion(q, q == 2 ? rZ : 2.0, k.roundTris[(std::size_t)q], o,
                                         q == 2 ? devZ : 0.0));
    Provisional pr;
    pr.tris = k.sphereTris;
    std::sort(pr.tris.begin(), pr.tris.end());
    pr.claim = ProvClaim::Unclaimed;
    w.provisionals.push_back(pr);
    const MeshView mv = k.s.view();
    SegmentParams p;
    DerivedTols tol;
    if (!claimSphereCorners(mv, p, tol, w)) return -1;
    int nSph = 0;
    for (const Region& g : w.accepted)
        if (g.type == SurfType::Sphere && g.origin == Origin::SphereCorner) {
            ++nSph;
            if (out) *out = g;
            if (dvol) *dvol = g.dVolPredicted;
        }
    return nSph;
}

static void caseCorner() {
    const double r = 2.0;
    const int n = 6;
    {
        const Corner k = makeCorner(r, n, 5.0, r, 0.0);
        Region g;
        double dv = 0;
        const int ns = runCorner(k, r, Origin::FilletStrip, 0.0, &g, &dv);
        expect(ns == 1, "box corner of three R rounds: one sphere");
        expect(ns == 1 && (int)g.tris.size() == n * n, "the claim is the whole octant");
        expect(ns == 1 && std::fabs(g.radius - r) < 1e-12 &&
                   g.ax.Location().XYZ().Modulus() < 1e-12,
               "centre is the rounds' and the fitted radius matches them");
        expect(ns == 1 && g.outwardNormal, "convex corner: outward");
        // Polar frame: one round perpendicular to the other two gives the axis.
        const gp_XYZ zd = g.ax.Direction().XYZ();
        expect(ns == 1 && std::fabs(std::fabs(zd.X()) + std::fabs(zd.Y()) + std::fabs(zd.Z()) - 1.0) <
                              1e-12,
               "perpendicular rounds: the frame axis is a round's axis");
        // dVol: the octant's sector minus the facets' cones from the centre.
        double poly = 0;
        const MeshView mv = k.s.view();
        for (int t : k.sphereTris) {
            const gp_XYZ a = mv.pts[mv.tris[t][0]], b = mv.pts[mv.tris[t][1]],
                         c = mv.pts[mv.tris[t][2]];
            poly += std::fabs(a.Dot(b.Crossed(c))) / 6.0;
        }
        const double want = kPi * r * r * r / 6.0 - poly;
        expect(ns == 1 && std::fabs(dv - want) < 1e-9, "dVolPredicted = sector - facet cones");
    }
    {
        const Corner k = makeCorner(r, n, 5.0, r + 3 * kTau, 0.0);
        expect(runCorner(k, r + 3 * kTau, Origin::FilletStrip, 0.0, nullptr, nullptr) == 0,
               "one round's radius 3 tau off: no sphere");
    }
    {
        const Corner k = makeCorner(r, n, 5.0, r, 2 * kTau);
        expect(runCorner(k, r, Origin::FilletStrip, 0.0, nullptr, nullptr) == 0,
               "one interior vertex 2 tau off the sphere: refused whole");
    }
    {
        const Corner k = makeCorner(r, n, 5.0, r, 0.0);
        expect(runCorner(k, r, Origin::CylGrow, 3 * kTau, nullptr, nullptr) == 0,
               "an uncertified round (CylGrow, 3 tau): no sphere");
    }
    {
        // One interior vertex 0.8 tau off the anchor radius. The least-squares
        // radius moves with it and every vertex stays within tau, so the
        // corner ships. The shipped radius is that fit, not the anchor mean.
        const Corner k = makeCorner(r, n, 5.0, r, 0.8 * kTau);
        Region g;
        const int ns = runCorner(k, r, Origin::FilletStrip, 0.0, &g, nullptr);
        const MeshView mv = k.s.view();
        double sum = 0.0;
        int nv = 0;
        std::vector<char> seen(k.s.pts.size(), 0);
        if (ns == 1) {
            for (int t : g.tris) {
                const int gt = mv.compTris ? mv.compTris[t] : t;
                for (int q = 0; q < 3; ++q) {
                    const int vi = mv.tris[gt][q];
                    if (vi < 0 || (std::size_t)vi >= seen.size() || seen[(std::size_t)vi]) continue;
                    seen[(std::size_t)vi] = 1;
                    sum += (k.s.pts[(std::size_t)vi] - g.ax.Location().XYZ()).Modulus();
                    ++nv;
                }
            }
        }
        const double want = nv > 0 ? sum / (double)nv : -1.0;
        expect(ns == 1 && nv > 0 && std::fabs(g.radius - want) < 1e-9,
               "shipped radius is the least-squares mean of the corner vertices");
        expect(ns == 1 && std::fabs(g.radius - r) > 1e-9 && g.maxVertexDev <= kTau,
               "the fitted radius leaves the anchor mean and still certifies");
    }
}

// A dissolved uncertified region leaves a hole in the accepted list. A torus
// whose neighbour index named that hole is cleared. A torus whose neighbour
// sat above the hole takes the dense index, not the stale one.
// The Z round of the box corner, retargeted onto a ring torus whose tube-centre
// circle passes through the origin. Axis Y, centre (-Rmaj, 0, 0), X along +X,
// so u = 0 is the origin and the contact meridian is the XY quarter circle the
// octant already shares. The far ring sits at azimuth u1, on the z < 0 side.
static void retargetFarRingToTorus(Synth& s, double r, double Rmaj, double u1, double L) {
    const gp_XYZ X(1, 0, 0), Yd(0, 0, -1), Za(0, 1, 0);
    const gp_XYZ O(-Rmaj, 0, 0);
    for (gp_XYZ& p : s.pts) {
        if (!(p.Z() < -0.5 * L)) continue;
        const double v = std::atan2(p.Y(), p.X());
        const gp_XYZ radial = X * std::cos(u1) + Yd * std::sin(u1);
        p = O + radial * (Rmaj + r * std::cos(v)) + Za * (r * std::sin(v));
    }
}

static Region torusRoundRegion(const gp_Ax3& fr, double Rmaj, double Rmin,
                               const std::vector<int>& tris, double dev) {
    Region g;
    g.type = SurfType::Torus;
    g.origin = Origin::TorusBlend;
    g.ax = fr;
    g.radius = Rmaj;
    g.radius2 = Rmin;
    g.tris = tris;
    g.maxVertexDev = dev;
    g.closed360 = false;
    return g;
}

static int runTorusCorner(const Corner& k, const gp_Ax3& fr, double Rmaj, double dev,
                          Region* out) {
    const double r = 2.0;
    SegmentWork w;
    w.accepted.push_back(roundRegion(0, r, k.roundTris[0], Origin::FilletStrip, 0.0));
    w.accepted.push_back(roundRegion(1, r, k.roundTris[1], Origin::FilletStrip, 0.0));
    w.accepted.push_back(torusRoundRegion(fr, Rmaj, r, k.roundTris[2], dev));
    Provisional pr;
    pr.tris = k.sphereTris;
    std::sort(pr.tris.begin(), pr.tris.end());
    pr.claim = ProvClaim::Unclaimed;
    w.provisionals.push_back(pr);
    const MeshView mv = k.s.view();
    SegmentParams p;
    DerivedTols tol;
    if (!claimSphereCorners(mv, p, tol, w)) return -1;
    int nSph = 0;
    for (const Region& g : w.accepted)
        if (g.type == SurfType::Sphere && g.origin == Origin::SphereCorner) {
            ++nSph;
            if (out) *out = g;
        }
    return nSph;
}

static void caseTorusCorner() {
    const double r = 2.0;
    const double Rmaj = 10.0;
    const double L = 5.0;
    const double u1 = 0.4;
    const int n = 6;
    const gp_Ax3 fr(gp_Pnt(-Rmaj, 0, 0), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
    Corner k = makeCorner(r, n, L, r, 0.0);
    retargetFarRingToTorus(k.s, r, Rmaj, u1, L);
    {
        Region g;
        const int ns = runTorusCorner(k, fr, Rmaj, 0.0, &g);
        expect(ns == 1, "two cylinders and one torus of equal radius: one sphere");
        expect(ns == 1 && std::fabs(g.radius - r) < 1e-9 && g.ax.Location().XYZ().Modulus() < 1e-8,
               "the sphere centre and radius are the rounds'");
        double u = -1.0;
        const bool hit = ns == 1 && sphereTorusContact(g.ax.Location().XYZ(), g.radius, fr, Rmaj, r,
                                                       tanDistTol(kTau, kTau), u);
        const gp_Circ iso = ElSLib::TorusUIso(fr, Rmaj, r, hit ? u : 0.0);
        const gp_Pnt ring(r, 0, 0);
        expect(hit && std::fabs(u) < 1e-8 && std::fabs(iso.Radius() - r) < 1e-9 &&
                   iso.Distance(ring) < 1e-8,
               "the torus contact is the u-iso through the corner");
    }
    {
        // Shift the tube circle so the corner centre is past tanDistTol.
        const double off = tanDistTol(kTau, kTau) + kTau;
        const gp_Ax3 far(gp_Pnt(-Rmaj + off, 0, 0), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
        expect(runTorusCorner(k, far, Rmaj, 0.0, nullptr) == 0,
               "corner centre off the tube circle by more than tanDistTol: refused");
    }
    {
        expect(runTorusCorner(k, fr, Rmaj, 3.0 * kTau, nullptr) == 0,
               "a torus whose maxVertexDev exceeds tau is not a round");
    }
}

// A CylGrow that owns the corner and one triangle off the sphere is not one
// preimage. F-D leaves it whole: no sphere, and the region is not dissolved.
static void casePartialRob() {
    const double r = 2.0;
    const int n = 6;
    Corner k = makeCorner(r, n, 5.0, r, 0.0);
    const int a = k.s.add(gp_XYZ(8, 0, 0));
    const int b = k.s.add(gp_XYZ(8, 1, 0));
    const int c = k.s.add(gp_XYZ(8, 0, 1));
    const int far = k.s.addTri(a, b, c);
    k.s.finish();
    SegmentWork w;
    for (int q = 0; q < 3; ++q)
        w.accepted.push_back(roundRegion(q, r, k.roundTris[(std::size_t)q], Origin::FilletStrip, 0.0));
    std::vector<int> held = k.sphereTris;
    held.push_back(far);
    w.accepted.push_back(roundRegion(0, r, held, Origin::CylGrow, 3 * kTau));
    const int nBefore = (int)held.size();
    const MeshView mv = k.s.view();
    SegmentParams p;
    DerivedTols tol;
    expect(claimSphereCorners(mv, p, tol, w), "partial rob returns");
    int nSph = 0, heldLeft = -1;
    for (const Region& g : w.accepted) {
        if (g.type == SurfType::Sphere && g.origin == Origin::SphereCorner) ++nSph;
        if (g.origin == Origin::CylGrow && g.maxVertexDev > kTau) heldLeft = (int)g.tris.size();
    }
    expect(nSph == 0, "a partial CylGrow is not a sphere preimage");
    expect(heldLeft == nBefore, "a partial CylGrow is not dissolved");
}

static void caseDissolveNbr() {
    const double r = 2.0;
    const int n = 6;
    const Corner k = makeCorner(r, n, 5.0, r, 0.0);
    SegmentWork w;
    for (int q = 0; q < 3; ++q)
        w.accepted.push_back(roundRegion(q, r, k.roundTris[(std::size_t)q], Origin::FilletStrip, 0.0));
    // Index 3: uncertified CylGrow holding the octant. The corner robs it.
    w.accepted.push_back(roundRegion(0, r, k.sphereTris, Origin::CylGrow, 3 * kTau));
    Region torGone;
    torGone.type = SurfType::Torus;
    torGone.origin = Origin::TorusBlend;
    torGone.radius = 4.0;
    torGone.radius2 = 1.0;
    torGone.filletNbrA = 3;  // the region that dissolves
    w.accepted.push_back(torGone);
    Region donor;
    donor.type = SurfType::Cylinder;
    donor.origin = Origin::CylGrow;
    donor.radius = r;
    donor.maxVertexDev = 0.0;
    w.accepted.push_back(donor);  // index 5, survives
    Region torLive;
    torLive.type = SurfType::Torus;
    torLive.origin = Origin::TorusBlend;
    torLive.radius = 4.0;
    torLive.radius2 = 2.0;
    torLive.filletNbrA = 5;  // the surviving donor, above the hole
    w.accepted.push_back(torLive);
    Region cone;
    cone.type = SurfType::Cone;
    cone.origin = Origin::ChamferCone;
    cone.radius = r;
    cone.filletNbrA = 5;
    w.accepted.push_back(cone);
    const MeshView mv = k.s.view();
    SegmentParams p;
    DerivedTols tol;
    expect(claimSphereCorners(mv, p, tol, w), "dissolve claim returns");
    int nSph = 0, donorAt = -1, torGoneNbr = 99, torLiveNbr = 99, coneNbr = 99;
    for (int i = 0; i < (int)w.accepted.size(); ++i) {
        const Region& g = w.accepted[(std::size_t)i];
        if (g.type == SurfType::Sphere && g.origin == Origin::SphereCorner) ++nSph;
        if (g.origin == Origin::TorusBlend && g.radius2 == 1.0) torGoneNbr = g.filletNbrA;
        if (g.origin == Origin::TorusBlend && g.radius2 == 2.0) torLiveNbr = g.filletNbrA;
        if (g.origin == Origin::ChamferCone) coneNbr = g.filletNbrA;
        if (g.origin == Origin::CylGrow && g.tris.empty() && g.maxVertexDev == 0.0) donorAt = i;
    }
    expect(nSph == 1, "the corner still claims one sphere");
    expect(torGoneNbr == -1, "a torus naming the dissolved region is cleared");
    expect(donorAt >= 0 && torLiveNbr == donorAt,
           "a torus naming a later donor takes the dense index");
    expect(donorAt >= 0 && coneNbr == donorAt, "a chamfer cone's neighbour is remapped too");
    expect(torLiveNbr >= 0 && torLiveNbr < (int)w.accepted.size(),
           "the remapped torus index is inside the kept list");
}

int main() {
    casePredicates();
    caseRoundCert();
    caseCorner();
    caseTorusCorner();
    casePartialRob();
    caseDissolveNbr();
    std::fprintf(stderr, "relational_sphere_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS",
                 gFail);
    return gFail ? 1 : 0;
}
