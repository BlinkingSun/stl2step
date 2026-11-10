// D-gid-7.6: dVolPredicted of FilletStrip cylinders and tori is the exact chord
// volume over the region's own triangles (dVolExactOfRegion, torusChordTerm,
// torusRimTerm). Synthetic meshes only; every truth is analytic or an
// independent fine quadrature. No training mesh.
// SPDX-License-Identifier: MIT

#include "refit_band_volume.hpp"
#include "refit_internal.hpp"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <map>
#include <utility>
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
    MeshView view() const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int (*)[3]>(tri.data());
        mv.compTris = comp.data();
        mv.nTri = comp.size();
        mv.nVtx = pts.size();
        mv.diag = 100.0;
        return mv;
    }
    std::vector<int> all() const {
        std::vector<int> t(comp.size());
        for (std::size_t i = 0; i < t.size(); ++i) t[i] = (int)i;
        return t;
    }
};

static double tetSum(const std::vector<gp_XYZ>& P, const std::vector<int>& T) {
    double v = 0.0;
    for (std::size_t i = 0; i + 2 < T.size(); i += 3)
        v += P[T[i]].Dot(P[T[i + 1]].Crossed(P[T[i + 2]])) / 6.0;
    return v;
}

// Orient a closed triangle surface consistently (propagate from the first
// triangle across shared edges); every edge must be shared by exactly two.
static void orientClosed(std::vector<int>& T) {
    const std::size_t n = T.size() / 3;
    std::map<std::pair<int, int>, std::vector<std::size_t>> use;
    for (std::size_t f = 0; f < n; ++f)
        for (int e = 0; e < 3; ++e) {
            const int a = T[3 * f + e], b = T[3 * f + (e + 1) % 3];
            use[{std::min(a, b), std::max(a, b)}].push_back(f);
        }
    std::vector<int> done(n, 0);
    std::vector<std::size_t> q{0};
    done[0] = 1;
    auto hasDir = [&](std::size_t f, int a, int b) {
        for (int e = 0; e < 3; ++e)
            if (T[3 * f + e] == a && T[3 * f + (e + 1) % 3] == b) return true;
        return false;
    };
    while (!q.empty()) {
        const std::size_t f = q.back();
        q.pop_back();
        for (int e = 0; e < 3; ++e) {
            const int a = T[3 * f + e], b = T[3 * f + (e + 1) % 3];
            for (std::size_t g : use[{std::min(a, b), std::max(a, b)}]) {
                if (g == f || done[g]) continue;
                if (hasDir(g, a, b)) std::swap(T[3 * g + 1], T[3 * g + 2]);
                done[g] = 1;
                q.push_back(g);
            }
        }
    }
}

// Irregular steps that still sum to the span: weights w_i.
static std::vector<double> steps(double a, double b, const std::vector<double>& w) {
    double s = 0.0;
    for (double x : w) s += x;
    std::vector<double> out{a};
    double acc = a;
    for (double x : w) {
        acc += (b - a) * x / s;
        out.push_back(acc);
    }
    out.back() = b;
    return out;
}

static gp_XYZ torusPt(double R, double r, double phi, double psi) {
    const double rho = R + r * std::cos(psi);
    return gp_XYZ(rho * std::cos(phi), rho * std::sin(phi), r * std::sin(psi));
}

// --- A. convex quarter-round FilletStrip with unequal rows ------------------
// Truth: L * sum_i R^2 (g_i - sin g_i) / 2 (the circular segments of the rows).
static void caseStrip() {
    const double R = 2.25, L = 7.0;
    const std::vector<double> th = steps(0.0, 0.5 * kPi, {1, 2, 1, 3, 1, 2, 2});
    Mesh m;
    std::vector<int> lo, hi;
    for (double t : th) {
        lo.push_back(m.add(gp_XYZ(R * std::cos(t), R * std::sin(t), 0.0)));
        hi.push_back(m.add(gp_XYZ(R * std::cos(t), R * std::sin(t), L)));
    }
    double truth = 0.0;
    for (std::size_t i = 0; i + 1 < th.size(); ++i) {
        const double g = th[i + 1] - th[i];
        truth += L * R * R * (g - std::sin(g)) / 2.0;
        // outward winding (normal away from the axis), alternating diagonals
        if (i % 2 == 0) {
            m.face(lo[i], lo[i + 1], hi[i + 1]);
            m.face(lo[i], hi[i + 1], hi[i]);
        } else {
            m.face(lo[i], lo[i + 1], hi[i]);
            m.face(lo[i + 1], hi[i + 1], hi[i]);
        }
    }
    const MeshView mv = m.view();
    Region r;
    r.type = SurfType::Cylinder;
    r.origin = Origin::FilletStrip;
    r.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.radius = R;
    r.outwardNormal = true;
    r.tris = m.all();
    double got = 0.0;
    const bool ok = dVolExactOfRegion(mv, r, got);
    expect(ok && close(got, truth, 1e-12), "strip: convex quarter round, unequal rows = sum of segments",
           got, truth);
    r.outwardNormal = false;
    dVolExactOfRegion(mv, r, got);
    expect(close(got, -truth, 1e-12), "strip: concave sign (material off the axis side)", got, -truth);
    // The regular-N density it replaces: with the true row count it misses the
    // unequal rows; with the side count doubled (two normals per row) it is ~4x low.
    double area = 0.0;
    for (int t : r.tris) {
        const gp_XYZ a = m.pts[m.tri[3 * t]], b = m.pts[m.tri[3 * t + 1]], c = m.pts[m.tri[3 * t + 2]];
        area += 0.5 * (b - a).Crossed(c - a).Modulus();
    }
    const int nRows = (int)th.size() - 1;
    const double reg = dVolCylinderSector(area, R, 4 * nRows, true);
    const double reg2 = dVolCylinderSector(area, R, 8 * nRows, true);
    std::fprintf(stderr, "INFO strip regularN(nSides=%d)=%.9g regularN(nSides=%d)=%.9g exact=%.9g\n",
                 4 * nRows, reg, 8 * nRows, reg2, truth);
    expect(reg < truth * 0.95 && reg2 < truth * 0.3,
           "strip: the regular-N estimate under-counts unequal / over-split rows", reg, truth);
    // D-gid-7.15 / D-gid-7.17: a grown cylinder (CylGrow) takes the same exact
    // volume when its own certificate holds (maxVertexDev <= tau = 2q), and
    // keeps its stage estimate when it does not.
    r.origin = Origin::CylGrow;
    r.outwardNormal = true;
    Mesh mq = m;
    MeshView mvq = mq.view();
    mvq.quantFloor = 1e-6;  // tau = 2e-6
    r.maxVertexDev = 2e-6;
    expect(dVolExactOfRegion(mvq, r, got) && close(got, truth, 1e-12),
           "cylinder: CylGrow at maxVertexDev = tau takes the exact volume", got, truth);
    r.maxVertexDev = 2.5e-6;
    expect(!dVolExactOfRegion(mvq, r, got),
           "cylinder: CylGrow over tau keeps its stage estimate (D-gid-7.17)", 0, 0);
    r.maxVertexDev = 2e-6;
    r.origin = Origin::FilletStrip;
    expect(dVolExactOfRegion(mvq, r, got), "cylinder: FilletStrip is not gated by tau", got, truth);
    // Any other origin keeps its own prediction.
    r.origin = Origin::NgonWall;
    expect(!dVolExactOfRegion(mv, r, got), "cylinder: NgonWall is not rewritten", 0, 0);
}

// --- A2. closed bore (CylGrow, closed360) with unequal azimuth steps ----------
// Truth: L * (pi R^2 - polygon area), the polygon of the irregular N-gon. A
// regular N-gon reproduces dVolCylinderSector.
static void caseBore() {
    const double R = 2.75, L = 11.0;
    for (int irregular = 0; irregular < 2; ++irregular) {
        std::vector<double> w;
        for (int i = 0; i < 40; ++i) w.push_back(irregular ? 1.0 + 0.6 * ((i * 7) % 5) / 4.0 : 1.0);
        const std::vector<double> th = steps(0.0, 2 * kPi, w);
        Mesh m;
        std::vector<int> lo, hi;
        for (std::size_t i = 0; i + 1 < th.size(); ++i) {
            lo.push_back(m.add(gp_XYZ(R * std::cos(th[i]), R * std::sin(th[i]), 0.0)));
            hi.push_back(m.add(gp_XYZ(R * std::cos(th[i]), R * std::sin(th[i]), L)));
        }
        const int n = (int)lo.size();
        double poly = 0.0;
        for (int i = 0; i < n; ++i) {
            const int j = (i + 1) % n;
            poly += 0.5 * R * R * std::sin(th[i + 1] - th[i]);
            // hole: the mesh normal points to the axis (material outside)
            m.face(lo[i], hi[j], lo[j]);
            m.face(lo[i], hi[i], hi[j]);
        }
        const double truth = -L * (kPi * R * R - poly);
        Region r;
        r.type = SurfType::Cylinder;
        r.origin = Origin::CylGrow;
        r.closed360 = true;
        r.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
        r.radius = R;
        r.outwardNormal = false;
        r.tris = m.all();
        double got = 0.0;
        const bool ok = dVolExactOfRegion(m.view(), r, got);
        expect(ok && close(got, truth, 1e-12),
               irregular ? "bore: closed CylGrow, unequal steps = L (pi R^2 - polygon)"
                         : "bore: closed CylGrow, regular 40-gon = L (pi R^2 - polygon)",
               got, truth);
        if (!irregular) {
            double area = 0.0;
            for (int t : r.tris) {
                const gp_XYZ a = m.pts[m.tri[3 * t]], b = m.pts[m.tri[3 * t + 1]],
                             c = m.pts[m.tri[3 * t + 2]];
                area += 0.5 * (b - a).Crossed(c - a).Modulus();
            }
            const double reg = dVolCylinderSector(area, R, 40, false);
            expect(close(reg, truth, 1e-12), "bore: regular 40-gon equals dVolCylinderSector", reg,
                   truth);
        }
    }
}

// Grid on a torus between phi0..phi1 and psi0..psi1, outward (phi x psi) winding.
static void torusGrid(Mesh& m, double R, double r, const std::vector<double>& ph,
                      const std::vector<double>& ps, bool wrapPhi, bool wrapPsi,
                      std::vector<std::vector<int>>& id) {
    const std::size_t nf = wrapPhi ? ph.size() - 1 : ph.size();
    const std::size_t ns = wrapPsi ? ps.size() - 1 : ps.size();
    id.assign(ph.size(), std::vector<int>(ps.size(), -1));
    for (std::size_t i = 0; i < nf; ++i)
        for (std::size_t j = 0; j < ns; ++j) id[i][j] = m.add(torusPt(R, r, ph[i], ps[j]));
    for (std::size_t j = 0; j < ps.size(); ++j)
        if (wrapPhi) id[ph.size() - 1][j] = id[0][j];
    for (std::size_t i = 0; i < ph.size(); ++i)
        if (wrapPsi) id[i][ps.size() - 1] = id[i][0];
    for (std::size_t i = 0; i + 1 < ph.size(); ++i)
        for (std::size_t j = 0; j + 1 < ps.size(); ++j) {
            const int a = id[i][j], b = id[i + 1][j], c = id[i + 1][j + 1], d = id[i][j + 1];
            if ((i + j) % 2 == 0) {
                m.face(a, b, c);
                m.face(a, c, d);
            } else {
                m.face(a, b, d);
                m.face(b, c, d);
            }
        }
}

static Region torusRegion(double R, double r, bool outward, const std::vector<int>& tris) {
    Region t;
    t.type = SurfType::Torus;
    t.origin = Origin::TorusBlend;
    t.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    t.radius = R;
    t.radius2 = r;
    t.outwardNormal = outward;
    t.tris = tris;
    return t;
}

// --- B. full torus: truth 2 pi^2 R r^2 - polyhedron -------------------------
static void caseFullTorus() {
    const double R = 6.0, r = 1.75;
    const std::vector<double> ph = steps(0.0, 2 * kPi, {1, 1.3, 0.7, 1, 2, 1, 1.5, 0.5, 1, 1, 1.2, 0.8,
                                                         1, 1, 1, 1, 1.1, 0.9, 1, 1, 1, 1, 1, 1});
    // one very thin row exercises the series branch of S0/S1
    const std::vector<double> ps = steps(-kPi, kPi, {1, 1, 0.002, 1.5, 0.5, 1, 1, 2, 1, 1, 0.7, 1.3});
    Mesh m;
    std::vector<std::vector<int>> id;
    torusGrid(m, R, r, ph, ps, true, true, id);
    const double truth = 2 * kPi * kPi * R * r * r - tetSum(m.pts, m.tri);
    double got = 0.0;
    const bool ok = dVolExactOfRegion(m.view(), torusRegion(R, r, true, m.all()), got);
    expect(ok && close(got, truth, 1e-10), "torus: closed ring (saddle side included)", got, truth);
}

// --- C. corner patch closed by flat faces ------------------------------------
// The band, two meridian polygons (to the axis) and two fans in the planes
// v = vLo, v = vHi close a polyhedron P. The solid bounded by the torus patch
// and the same flat faces is  Phi * integral rho_S(v)^2 / 2 dv  (Pappus).
static void cornerPatch(bool outer, bool coarsePhi, const char* name) {
    const double R = 5.0, r = 1.25;
    const double Phi = 0.5 * kPi;
    const double psLo = outer ? 0.0 : kPi, psHi = outer ? 0.5 * kPi : 1.5 * kPi;
    const std::vector<double> ph = coarsePhi
                                       ? steps(0.0, Phi, {1, 1.3, 0.7})
                                       : steps(0.0, Phi, {1, 1.4, 0.6, 1, 1, 2, 1, 1, 0.5, 1.5, 1, 1});
    const std::vector<double> ps = steps(psLo, psHi, {1, 2, 1, 1.5, 0.5, 1, 1});
    Mesh m;
    std::vector<std::vector<int>> id;
    torusGrid(m, R, r, ph, ps, false, false, id);
    const std::vector<int> band = m.all();
    // Closing faces. Axis points at the two rim heights; each meridian polygon
    // (axis, profile, axis) is a fan from its low axis point.
    const double vA = r * std::sin(psLo), vB = r * std::sin(psHi);
    const int axA = m.add(gp_XYZ(0, 0, vA)), axB = m.add(gp_XYZ(0, 0, vB));
    const std::size_t I = ph.size() - 1, J = ps.size() - 1;
    std::vector<int> closed = m.tri;
    auto cf = [&](int a, int b, int c) {
        closed.push_back(a);
        closed.push_back(b);
        closed.push_back(c);
    };
    for (std::size_t i = 0; i < I; ++i) {
        cf(axA, id[i][0], id[i + 1][0]);
        cf(axB, id[i][J], id[i + 1][J]);
    }
    for (std::size_t j = 0; j < J; ++j) {
        cf(axA, id[0][j], id[0][j + 1]);
        cf(axA, id[I][j], id[I][j + 1]);
    }
    cf(axA, id[0][J], axB);
    cf(axA, id[I][J], axB);
    orientClosed(closed);
    double vPoly = std::fabs(tetSum(m.pts, closed));
    // V_true = Phi/2 * integral (R +- sqrt(r^2 - v^2))^2 dv over [min(vA,vB), max]
    const double sgn = outer ? 1.0 : -1.0;
    auto F = [&](double v) {
        const double q = std::sqrt(std::max(0.0, (r - v) * (r + v)));
        const double A = 0.5 * (v * q + r * r * std::asin(std::max(-1.0, std::min(1.0, v / r))));
        return (R * R + r * r) * v - v * v * v / 3.0 + 2.0 * sgn * R * A;
    };
    const double lo = std::min(vA, vB), hi = std::max(vA, vB);
    const double vTrue = 0.5 * Phi * (F(hi) - F(lo));
    const double truth = vTrue - vPoly;
    // Material: inside the tube for the outer corner (mesh normal along the tube
    // normal); for the inner corner the material is on the axis side, outside
    // the tube.
    double got = 0.0, phys = 0.0;
    bool ok = dVolExactOfRegion(m.view(), torusRegion(R, r, outer, band), got, &phys);
    // The inner band was wound phi x psi (= outward tube normal); its material
    // is on the other side, so the mesh is re-wound for that region.
    if (!outer) {
        Mesh w = m;
        for (std::size_t k = 0; k < band.size(); ++k) std::swap(w.tri[3 * k + 1], w.tri[3 * k + 2]);
        ok = dVolExactOfRegion(w.view(), torusRegion(R, r, false, band), got, &phys);
    }
    // physicalOut is the net material the surface adds; dVolPredicted carries
    // I9's sign (outwardNormal) with the same magnitude.
    expect(ok && close(phys, truth, 1e-10), name, phys, truth);
    const double sigma = outer ? 1.0 : -1.0;
    expect(close(got, sigma * std::fabs(truth), 1e-10) && (got == 0.0 || (got > 0) == outer),
           "torus: dVolPredicted = sigma * |net| (I9 sign, exact magnitude)", got,
           sigma * std::fabs(truth));
}

// --- D. ragged rim: one triangle off the band ---------------------------------
// Additivity (the oblique rim edges cancel between the two pieces) and the
// single triangle against an independent fine quadrature of the volume between
// it and the torus along the rays perpendicular to the axis.
static double fineRadialVolume(const gp_XYZ& p0, const gp_XYZ& p1, const gp_XYZ& p2, double R,
                               double r, int n) {
    // V = integral_T (rho_S(v)^2 - rho^2) / 2 |d(theta, v) / d(s, t)| ds dt
    static const double w7[7] = {9.0 / 40, 0.13239415278850619, 0.13239415278850619,
                                 0.13239415278850619, 0.12593918054482714, 0.12593918054482714,
                                 0.12593918054482714};
    static const double a7[7][2] = {{1.0 / 3, 1.0 / 3},
                                    {0.059715871789769820, 0.47014206410511509},
                                    {0.47014206410511509, 0.059715871789769820},
                                    {0.47014206410511509, 0.47014206410511509},
                                    {0.79742698535308732, 0.10128650732345634},
                                    {0.10128650732345634, 0.79742698535308732},
                                    {0.10128650732345634, 0.10128650732345634}};
    const gp_XYZ e1 = p1 - p0, e2 = p2 - p0;
    auto f = [&](double s, double t) {
        const gp_XYZ x = p0 + e1 * s + e2 * t;
        const double rho2 = x.X() * x.X() + x.Y() * x.Y();
        const double v = x.Z();
        const double rs = R + std::sqrt(std::max(0.0, (r - v) * (r + v)));
        // theta_s = (x e_y - y e_x) / rho^2, v_s = e_z
        const double ts = (x.X() * e1.Y() - x.Y() * e1.X()) / rho2;
        const double tt = (x.X() * e2.Y() - x.Y() * e2.X()) / rho2;
        const double J = std::fabs(ts * e2.Z() - tt * e1.Z());
        return 0.5 * (rs * rs - rho2) * J;
    };
    double v = 0.0;
    const double h = 1.0 / n;
    for (int i = 0; i < n; ++i)
        for (int j = 0; i + j < n; ++j) {
            for (int up = 0; up < 2; ++up) {
                if (up && i + j + 1 >= n) continue;
                // sub-triangle corners in (s,t)
                double cs[3][2];
                if (!up) {
                    cs[0][0] = i * h; cs[0][1] = j * h;
                    cs[1][0] = (i + 1) * h; cs[1][1] = j * h;
                    cs[2][0] = i * h; cs[2][1] = (j + 1) * h;
                } else {
                    cs[0][0] = (i + 1) * h; cs[0][1] = j * h;
                    cs[1][0] = (i + 1) * h; cs[1][1] = (j + 1) * h;
                    cs[2][0] = i * h; cs[2][1] = (j + 1) * h;
                }
                double acc = 0.0;
                for (int k = 0; k < 7; ++k) {
                    const double l1 = a7[k][0], l2 = a7[k][1], l0 = 1 - l1 - l2;
                    const double s = l0 * cs[0][0] + l1 * cs[1][0] + l2 * cs[2][0];
                    const double t = l0 * cs[0][1] + l1 * cs[1][1] + l2 * cs[2][1];
                    acc += w7[k] * f(s, t);
                }
                v += acc * 0.5 * h * h;
            }
        }
    return v;
}

static void caseRaggedRim() {
    const double R = 5.0, r = 1.25;
    const std::vector<double> ph = steps(0.0, 0.5 * kPi, {1, 1, 1, 1, 1, 1, 1, 1});
    const std::vector<double> ps = steps(0.0, 0.5 * kPi, {1, 1, 1, 1, 1, 1});
    Mesh m;
    std::vector<std::vector<int>> id;
    torusGrid(m, R, r, ph, ps, false, false, id);
    // Notch: a triangle on the phi = 0 rim at mid-profile, and one on the equator rim.
    const int notchA = 2 * (0 * 6 + 3);      // quad (i=0, j=3), first triangle
    const int notchB = 2 * (4 * 6 + 0);      // quad (i=4, j=0), first triangle
    std::vector<int> rest, onlyA{notchA}, onlyB{notchB};
    for (int t : m.all())
        if (t != notchA && t != notchB) rest.push_back(t);
    const MeshView mv = m.view();
    double full = 0, vr = 0, va = 0, vb = 0;
    dVolExactOfRegion(mv, torusRegion(R, r, true, m.all()), full);
    dVolExactOfRegion(mv, torusRegion(R, r, true, rest), vr);
    dVolExactOfRegion(mv, torusRegion(R, r, true, onlyA), va);
    dVolExactOfRegion(mv, torusRegion(R, r, true, onlyB), vb);
    expect(close(vr + va + vb, full, 1e-12), "torus: ragged rim pieces add up to the band",
           vr + va + vb, full);
    for (int t : {notchA, notchB}) {
        const gp_XYZ a = m.pts[m.tri[3 * t]], b = m.pts[m.tri[3 * t + 1]], c = m.pts[m.tri[3 * t + 2]];
        const double ref = fineRadialVolume(a, b, c, R, r, 96);
        const double got = t == notchA ? va : vb;
        expect(close(got, ref, 1e-9),
               t == notchA ? "torus: oblique rim triangle (meridian notch) = fine quadrature"
                           : "torus: oblique rim triangle (equator notch) = fine quadrature",
               got, ref);
    }
    // Without the rim term the oblique edges are taken straight in (phi, psi):
    // the same triangle then misses by far more than the quadrature's error.
    {
        const int t = notchA;
        const gp_XYZ a = m.pts[m.tri[3 * t]], b = m.pts[m.tri[3 * t + 1]], c = m.pts[m.tri[3 * t + 2]];
        const gp_XYZ o(0, 0, 0), ax(0, 0, 1), e1(1, 0, 0), e2(0, 1, 0);
        const double base = torusChordTerm(a, b, c, o, ax, e1, e2, R, r).vol;
        const double ref = fineRadialVolume(a, b, c, R, r, 96);
        std::fprintf(stderr, "INFO notch straight-(phi,psi) only=%.12g with rim=%.12g ref=%.12g\n", base,
                     va, ref);
        expect(std::fabs(base - ref) > 1e3 * std::fabs(va - ref),
               "torus: the rim term is what makes the oblique rim exact", base, ref);
    }
}

int main() {
    caseStrip();
    caseBore();
    caseFullTorus();
    cornerPatch(true, false, "torus: outer corner quarter round (convex) = Pappus - polyhedron");
    cornerPatch(false, false, "torus: inner corner quarter (saddle) = Pappus - polyhedron");
    // Few columns: the saddle's major-circle chords dominate and the net material
    // the surface adds is positive although outwardNormal is false.
    cornerPatch(false, true, "torus: inner corner, coarse azimuth (net opposes sigma) = Pappus - polyhedron");
    caseRaggedRim();
    std::fprintf(stderr, "%s relational_dvol_exact_unit (%d failed)\n", gFail ? "FAIL" : "PASS", gFail);
    return gFail ? 1 : 0;
}
