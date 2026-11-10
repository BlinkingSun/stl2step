// D-gid-8.8: a tau-certified closed apex cone (a drill point) takes the exact
// chord volume of its fan (apexConeChordTerm via dVolExactOfRegion).
// D-gid-8.6: a ring-certified CylGrow wall takes the exact chord volume
// whatever the chord sag in its maxVertexDev. Synthetic meshes only; every
// truth is analytic: cone minus pyramid, h (pi R^2 - A) / 3, with A the rim
// polygon's area (an independent shoelace), and L (pi R^2 - A) for a bore.
// SPDX-License-Identifier: MIT

#include "refit_band_volume.hpp"
#include "refit_internal.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;
static constexpr double kPi = 3.14159265358979323846264338327950288;

static void expect(bool ok, const char* name, double got, double want) {
    std::fprintf(stderr, "%s %s got=%.15g want=%.15g diff=%.3g\n", ok ? "PASS" : "FAIL", name, got,
                 want, got - want);
    if (!ok) ++gFail;
}

static bool close(double got, double want, double rel) {
    return std::fabs(got - want) <= rel * std::fabs(want);
}

struct Mesh {
    std::vector<gp_XYZ> pts;
    std::vector<int> tri;
    std::vector<int> comp;
    int add(const gp_XYZ& p) {
        pts.push_back(p);
        return (int)pts.size() - 1;
    }
    void face(int a, int b, int c) {
        tri.push_back(a);
        tri.push_back(b);
        tri.push_back(c);
        comp.push_back((int)comp.size());
    }
    MeshView view(double q) const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int (*)[3]>(tri.data());
        mv.compTris = comp.data();
        mv.nTri = comp.size();
        mv.nVtx = pts.size();
        mv.diag = 100.0;
        mv.quantFloor = q;  // tau = 2q
        return mv;
    }
    std::vector<int> all() const {
        std::vector<int> t(comp.size());
        for (std::size_t i = 0; i < t.size(); ++i) t[i] = (int)i;
        return t;
    }
};

// A general frame (nothing axis-aligned).
struct Frame3 {
    gp_XYZ o, a, e1, e2;
};
static Frame3 frame() {
    Frame3 f;
    f.o = gp_XYZ(3.0, -2.0, 7.0);
    f.a = gp_XYZ(1.0, 2.0, 2.0) / 3.0;
    gp_XYZ e1 = gp_XYZ(2.0, -1.0, 0.0);
    f.e1 = e1 / e1.Modulus();
    f.e2 = f.a.Crossed(f.e1);
    return f;
}
static gp_XYZ at(const Frame3& f, double rho, double th, double u) {
    return f.o + f.e1 * (rho * std::cos(th)) + f.e2 * (rho * std::sin(th)) + f.a * u;
}

// Shoelace area of a rim polygon in the plane perpendicular to the axis.
static double rimArea(const Frame3& f, const std::vector<gp_XYZ>& rim) {
    double a = 0.0;
    for (std::size_t i = 0; i < rim.size(); ++i) {
        const gp_XYZ p = rim[i] - f.o, q = rim[(i + 1) % rim.size()] - f.o;
        a += p.Dot(f.e1) * q.Dot(f.e2) - q.Dot(f.e1) * p.Dot(f.e2);
    }
    return 0.5 * std::fabs(a);
}

static Region apexRegion(const Frame3& f, double R, double hSigned, bool outward,
                         const std::vector<int>& tris) {
    Region r;
    r.type = SurfType::Cone;
    r.origin = Origin::ChamferCone;
    r.ax = gp_Ax3(gp_Pnt(f.o.X(), f.o.Y(), f.o.Z()), gp_Dir(f.a.X(), f.a.Y(), f.a.Z()),
                  gp_Dir(f.e1.X(), f.e1.Y(), f.e1.Z()));
    r.radius = 0.0;
    r.vMin = 0.0;
    r.vMax = R;
    r.closed360 = true;
    r.maxVertexSnap = hSigned;
    r.outwardNormal = outward;
    r.maxVertexDev = 0.0;
    r.tris = tris;
    return r;
}

// Fan about an apex vertex; normalsOut = winding normal away from the axis.
static Mesh fan(const Frame3& f, const std::vector<gp_XYZ>& rim, const gp_XYZ& apexV,
                bool normalsOut) {
    Mesh m;
    const int A = m.add(apexV);
    std::vector<int> ids;
    for (const gp_XYZ& p : rim) ids.push_back(m.add(p));
    for (std::size_t i = 0; i < ids.size(); ++i) {
        int b = ids[i], c = ids[(i + 1) % ids.size()];
        const gp_XYZ pa = m.pts[A], pb = m.pts[b], pc = m.pts[c];
        const gp_XYZ cen = (pa + pb + pc) / 3.0 - f.o;
        const gp_XYZ n = (pb - pa).Crossed(pc - pa);
        const bool out = (cen - f.a * cen.Dot(f.a)).Dot(n) > 0.0;
        if (out != normalsOut) std::swap(b, c);
        m.face(A, b, c);
    }
    return m;
}

// --- A. regular 16-gon drill point (material outside the cone) --------------
// Truth: h (pi R^2 - (N/2) R^2 sin(2 pi / N)) / 3, negative (the cone removes
// material: outwardNormal false, winding normals toward the axis).
static void caseRegular() {
    const Frame3 f = frame();
    const double R = 1.15, h = 0.62;
    const int N = 16;
    std::vector<gp_XYZ> rim;
    for (int i = 0; i < N; ++i) rim.push_back(at(f, R, 2.0 * kPi * i / N + 0.3, h));
    Mesh m = fan(f, rim, f.o, false);
    const double truth = h * (kPi * R * R - 0.5 * N * R * R * std::sin(2.0 * kPi / N)) / 3.0;
    const MeshView mv = m.view(1e-6);
    double got = 0.0, phys = 0.0;
    const bool ok = dVolExactOfRegion(mv, apexRegion(f, R, h, false, m.all()), got, &phys);
    expect(ok && close(got, -truth, 1e-12), "apex: regular 16-gon drill point = -(cone - pyramid)",
           got, -truth);
    expect(close(phys, -truth, 1e-12), "apex: physical volume carries the winding sign", phys,
           -truth);
    // The stage formula it replaces (revolvedChordTerm through the apex vertex).
    double old = 0.0;
    for (int t : m.all()) {
        const gp_XYZ a = m.pts[m.tri[3 * t]], b = m.pts[m.tri[3 * t + 1]], c = m.pts[m.tri[3 * t + 2]];
        old += revolvedChordTerm(a, b, c, f.o, f.a, f.e1, f.e2, 0.0, R / h);
    }
    std::fprintf(stderr, "INFO apex: revolvedChordTerm through the apex = %.9g, truth %.9g\n", old,
                 truth);
    expect(!close(old, truth, 0.05), "apex: the straight-footprint term misreads the apex vertex",
           old, truth);
}

// --- B. irregular fan: unequal azimuths, a T-junction chord vertex on one rim
// chord, the apex vertex off the axis by a sub-tau amount, the mouth BELOW the
// apex (negative height), winding normals away from the axis. Truth: the cone
// holds pi R^2 |h| / 3 over the rim plane, the fan's pyramid A |h| / 3 whatever
// the apex's offset in that plane.
static void caseIrregular() {
    const Frame3 f = frame();
    const double R = 2.5, h = -1.7;
    const double w[] = {1, 3, 2, 1, 4, 2, 2, 1, 3, 1, 2, 2};
    double s = 0.0;
    for (double x : w) s += x;
    std::vector<gp_XYZ> rim;
    double acc = 0.1;
    for (int i = 0; i < 12; ++i) {
        rim.push_back(at(f, R, acc, h));
        if (i == 4) {  // chord vertex: on the chord to the next corner, inside the circle
            const gp_XYZ nx = at(f, R, acc + 2.0 * kPi * w[i] / s, h);
            rim.push_back(rim.back() * 0.6 + nx * 0.4);
        }
        acc += 2.0 * kPi * w[i] / s;
    }
    const gp_XYZ apexV = f.o + f.e1 * 3e-7 - f.e2 * 2e-7;  // sub-tau off the axis (tau = 2e-6)
    Mesh m = fan(f, rim, apexV, true);
    const double A = rimArea(f, rim);
    const double truth = std::fabs(h) * (kPi * R * R - A) / 3.0;
    const MeshView mv = m.view(1e-6);
    double got = 0.0, phys = 0.0;
    const bool ok = dVolExactOfRegion(mv, apexRegion(f, R, h, true, m.all()), got, &phys);
    expect(ok && close(got, truth, 1e-11),
           "apex: irregular fan, chord vertex, apex off axis, mouth below = cone - pyramid", got,
           truth);
    expect(close(phys, truth, 1e-11), "apex: physical = + (the cone lies beyond outward facets)",
           phys, truth);
    // Gates: the vertex certificate, the apex, the closure.
    Region r = apexRegion(f, R, h, true, m.all());
    r.maxVertexDev = 2.5e-6;
    expect(!dVolExactOfRegion(mv, r, got), "apex: over tau keeps its stage estimate", 0, 0);
    r = apexRegion(f, R, h, true, m.all());
    r.radius = 0.4;
    expect(!dVolExactOfRegion(mv, r, got), "apex: a frustum (radius > 0) is not this rule", 0, 0);
    r = apexRegion(f, R, h, true, m.all());
    r.closed360 = false;
    expect(!dVolExactOfRegion(mv, r, got), "apex: an open cone is not this rule", 0, 0);
}

// --- C. D-gid-8.6: a ring-certified CylGrow bore over tau ---------------------
// Truth: L (pi R^2 - A), A the irregular rim polygon.
static void caseRingCertifiedBore() {
    const Frame3 f = frame();
    const double R = 1.5, L = 4.0;
    const double w[] = {2, 1, 3, 1, 2, 2, 1, 3, 2, 1};
    double s = 0.0;
    for (double x : w) s += x;
    Mesh m;
    std::vector<int> lo, hi;
    std::vector<gp_XYZ> rim;
    double acc = 0.0;
    for (double x : w) {
        rim.push_back(at(f, R, acc, 0.0));
        lo.push_back(m.add(at(f, R, acc, 0.0)));
        hi.push_back(m.add(at(f, R, acc, L)));
        acc += 2.0 * kPi * x / s;
    }
    const std::size_t n = lo.size();
    for (std::size_t i = 0; i < n; ++i) {  // winding normals toward the axis (a hole)
        const std::size_t j = (i + 1) % n;
        m.face(lo[i], hi[i], lo[j]);
        m.face(lo[j], hi[i], hi[j]);
    }
    const double truth = L * (kPi * R * R - rimArea(f, rim));
    Region r;
    r.type = SurfType::Cylinder;
    r.origin = Origin::CylGrow;
    r.ax = gp_Ax3(gp_Pnt(f.o.X(), f.o.Y(), f.o.Z()), gp_Dir(f.a.X(), f.a.Y(), f.a.Z()),
                  gp_Dir(f.e1.X(), f.e1.Y(), f.e1.Z()));
    r.radius = R;
    r.closed360 = true;
    r.outwardNormal = false;
    r.tris = m.all();
    r.maxVertexDev = 1e-2;  // the chord sag the ring certificate excuses, far over tau
    const MeshView mv = m.view(1e-6);
    double got = 0.0;
    r.ringCertified = false;
    expect(!dVolExactOfRegion(mv, r, got), "bore: CylGrow over tau, no ring certificate: stage estimate",
           0, 0);
    r.ringCertified = true;
    const bool ok = dVolExactOfRegion(mv, r, got);
    expect(ok && close(got, -truth, 1e-12), "bore: ring-certified CylGrow over tau takes the exact volume",
           got, -truth);
}

int main() {
    caseRegular();
    caseIrregular();
    caseRingCertifiedBore();
    std::fprintf(stderr, "%s apex_cone_dvol_unit (%d failure(s))\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
