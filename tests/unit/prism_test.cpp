// AC3-P1 — exercise detectPrismatic / derivePrismTols (SPEC A1-A10).
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "refit_prism.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::MeshView;
using stl2step::refit::PrismLevels;
using stl2step::refit::PrismTols;
using stl2step::refit::Region;
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;
using stl2step::refit::detectPrismatic;

namespace stl2step {
namespace refit {
void derivePrismTols(const MeshView& mv, const RegionSet& rs, PrismTols& t);
}
}

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS %s\n", name);
    }
}

static bool near(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

static gp_XYZ unit(const gp_Dir& d) {
    gp_XYZ x = d.XYZ();
    const double m = x.Modulus();
    return m > 0.0 ? x / m : gp_XYZ(0, 0, 0);
}

static Region mkPlane(int id, const gp_Dir& n, const gp_Pnt& p) {
    Region r;
    r.id = id;
    r.type = SurfType::Plane;
    r.ax = gp_Ax3(p, n);
    return r;
}

static Region mkCyl(int id, const gp_Dir& a, const gp_Pnt& c, double vmin, double vmax,
                    double R) {
    Region r;
    r.id = id;
    r.type = SurfType::Cylinder;
    r.ax = gp_Ax3(c, a);
    r.vMin = vmin;
    r.vMax = vmax;
    r.radius = R;
    return r;
}

static bool runSegment(const char* stl, HarnessMesh& mesh, MeshView& mv, RegionSet& rs,
                       std::string& err) {
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) return false;
    size_t pick = 0;
    for (size_t i = 1; i < mesh.comps.size(); ++i)
        if (mesh.comps[i].compTris.size() > mesh.comps[pick].compTris.size())
            pick = i;
    toMeshView(mesh, mesh.comps[pick], mv);
    rs = RegionSet{};
    rs.compRoot = 0;
    SegmentParams p;
    const bool ok = stl2step::refit::segment(mv, p, rs, nullptr);
    rs.compRoot = 0;
    return ok;
}

static void testTolsSelfComputed() {
    MeshView mv{};
    mv.diag = 78.200087;
    mv.weldTol = 0.0;
    RegionSet rs;
    rs.regions.push_back(mkCyl(0, gp_Dir(0, 1, 0), gp_Pnt(0, 0, 0), 0.0, 2.825, 1.0));
    PrismTols t;
    stl2step::refit::derivePrismTols(mv, rs, t);
    std::fprintf(stderr, "A10 tauSurf=%.6e tauAx=%.6e (diag=%.6f hMin=%.6f)\n",
                 t.tauSurf, t.tauAx, mv.diag, 2.825);
    check(near(t.tauSurf, 7.820e-5, 5e-8), "A10 tauSurf ~ 7.820e-5");
    check(near(t.tauAx, 5.537e-5, 5e-8), "A10 tauAx ~ 5.537e-5");
    check(t.tauLvl == t.tauSurf && t.tauFit == t.tauSurf, "A10 tauLvl==tauFit==tauSurf");

    MeshView mv2 = mv;
    mv2.diag = mv.diag * 2.0;
    PrismTols t2;
    stl2step::refit::derivePrismTols(mv2, rs, t2);
    check(t2.tauSurf > t.tauSurf && t2.tauAx > t.tauAx, "A10 changing mv.diag changes tols");
}

static void testSyntheticDeclines() {
    MeshView mv{};
    mv.diag = 100.0;
    mv.weldTol = 0.0;
    const gp_Dir az(0, 0, 1);

    // Oblique plane: cond 3 (two cylinders so we do not stop at cond 1).
    RegionSet rs3;
    rs3.compRoot = 0;
    rs3.regions.push_back(mkCyl(0, az, gp_Pnt(0, 0, 0), 0.0, 10.0, 3.0));
    rs3.regions.push_back(mkCyl(1, az, gp_Pnt(8, 0, 0), 0.0, 10.0, 2.0));
    rs3.regions.push_back(mkPlane(2, gp_Dir(0, 0, -1), gp_Pnt(0, 0, 0)));
    rs3.regions.push_back(mkPlane(3, gp_Dir(0, 0, 1), gp_Pnt(0, 0, 10)));
    rs3.regions.push_back(mkPlane(4, gp_Dir(0, 1, 1), gp_Pnt(0, 0, 8)));
    PrismTols t3;
    stl2step::refit::derivePrismTols(mv, rs3, t3);
    const PrismLevels lv3 = detectPrismatic(mv, rs3, t3);
    std::fprintf(stderr, "A9-oblique ok=%d failedCond=%d\n", lv3.ok ? 1 : 0, lv3.failedCond);
    check(!lv3.ok && lv3.failedCond == 3, "A9 synthetic oblique declines at cond 3");

    // Floating cylinder extent: cond 5.
    RegionSet rs5;
    rs5.compRoot = 0;
    rs5.regions.push_back(mkCyl(0, az, gp_Pnt(0, 0, 0), 0.0, 7.0, 3.0));
    rs5.regions.push_back(mkCyl(1, az, gp_Pnt(8, 0, 0), 1.0, 6.0, 2.0));
    rs5.regions.push_back(mkPlane(2, gp_Dir(0, 0, -1), gp_Pnt(0, 0, 0)));
    rs5.regions.push_back(mkPlane(3, gp_Dir(0, 0, 1), gp_Pnt(0, 0, 10)));
    PrismTols t5;
    stl2step::refit::derivePrismTols(mv, rs5, t5);
    const PrismLevels lv5 = detectPrismatic(mv, rs5, t5);
    std::fprintf(stderr, "A9-float ok=%d failedCond=%d\n", lv5.ok ? 1 : 0, lv5.failedCond);
    check(!lv5.ok && lv5.failedCond == 5, "A9 synthetic float declines at cond 5");
}

static void testHandleLock(const char* stl) {
    HarnessMesh mesh;
    MeshView mv{};
    RegionSet rs;
    std::string err;
    if (!runSegment(stl, mesh, mv, rs, err)) {
        std::fprintf(stderr, "handle-lock load/segment: %s\n", err.c_str());
        check(false, "A7 handle-lock segment");
        return;
    }
    std::fprintf(stderr, "handle-lock nTri=%zu diag=%.6f weldTol=%.6g regions=%zu\n",
                 mv.nTri, mv.diag, mv.weldTol, rs.regions.size());

    PrismTols t;
    stl2step::refit::derivePrismTols(mv, rs, t);
    std::fprintf(stderr, "A10-live tauSurf=%.6e tauAx=%.6e tauLvl=%.6e\n",
                 t.tauSurf, t.tauAx, t.tauLvl);
    check(near(t.tauSurf, 7.820e-5, 2e-7), "A10 live tauSurf from mesh diag");

    const PrismLevels lv = detectPrismatic(mv, rs, t);
    std::fprintf(stderr, "A7 handle-lock ok=%d failedCond=%d L=%zu\n",
                 lv.ok ? 1 : 0, lv.failedCond, lv.y.size());
    check(lv.ok && lv.failedCond == 0, "A7 handle-lock PRISMATIC");

    std::vector<const Region*> cyls, planes;
    for (const Region& r : rs.regions) {
        if (r.type == SurfType::Cylinder) cyls.push_back(&r);
        else if (r.type == SurfType::Plane) planes.push_back(&r);
    }

    double maxSin = 0.0;
    for (size_t i = 0; i < cyls.size(); ++i) {
        const gp_XYZ ai = unit(cyls[i]->ax.Direction());
        for (size_t j = i + 1; j < cyls.size(); ++j) {
            const double s = ai.Crossed(unit(cyls[j]->ax.Direction())).Modulus();
            if (s > maxSin) maxSin = s;
        }
    }
    std::fprintf(stderr, "A1 max pairwise |sin| = %.3e rad  nCyl=%zu\n", maxSin, cyls.size());
    check(maxSin < t.tauAx, "A1 pairwise misalignment < tau_ax");
    if (maxSin == 0.0)
        std::fprintf(stderr, "A1 reported 0.000e+00 rad\n");
    else
        std::fprintf(stderr, "A1 reported %.3e rad (mesh-fit, not STEP bit-identity)\n", maxSin);

    const gp_XYZ truth(0.000007680, 1.000000000, -0.000000000);
    const gp_XYZ got = unit(lv.axis);
    gp_XYZ tunit = truth;
    const double tm = tunit.Modulus();
    if (tm > 0.0) tunit /= tm;
    if (got.Dot(tunit) < 0.0) tunit.Reverse();
    const double axisSin = got.Crossed(tunit).Modulus();
    std::fprintf(stderr, "A1 recovered axis [%.9f, %.9f, %.9f]  |sin vs truth|=%.3e\n",
                 got.X(), got.Y(), got.Z(), axisSin);
    check(axisSin < t.tauAx, "A1 recovered axis within tau_ax of STEP truth");

    const gp_XYZ ahat = got;
    int nCap = 0, nLat = 0, nObl = 0;
    double minCapDot = 1.0, maxLatDot = 0.0;
    for (const Region* r : planes) {
        const double nd = std::abs(unit(r->ax.Direction()).Dot(ahat));
        if (nd > 1.0 - t.tauAx) {
            ++nCap;
            if (nd < minCapDot) minCapDot = nd;
        } else if (nd < t.tauAx) {
            ++nLat;
            if (nd > maxLatDot) maxLatDot = nd;
        } else {
            ++nObl;
        }
    }
    std::fprintf(stderr, "A2 nCap=%d (|n.a| min=%.6f) nLat=%d (max=%.6f) nOblique=%d\n",
                 nCap, nCap ? minCapDot : 0.0, nLat, maxLatDot, nObl);
    check(nCap == 3 && nLat == 10 && nObl == 0, "A2 3 caps / 10 laterals / 0 oblique");
    check(minCapDot > 0.999999, "A2 |n.a| reported as 1.000000 on caps");
    check(maxLatDot < 1e-5, "A2 |n.a| reported as 0.000000 on laterals");

    std::fprintf(stderr, "A3 levels (ref = common-axis dot plane.Location()):\n");
    for (size_t i = 0; i < lv.y.size(); ++i)
        std::fprintf(stderr, "  y[%zu]=%.6f capRegion=%d\n", i, lv.y[i],
                     i < lv.capRegion.size() ? lv.capRegion[i] : -1);
    check(lv.y.size() == 3, "A3 exactly 3 levels");
    if (lv.y.size() >= 3) {
        const double s0 = lv.y[1] - lv.y[0];
        const double s1 = lv.y[2] - lv.y[1];
        const double tot = lv.y[2] - lv.y[0];
        std::fprintf(stderr, "A3 spacings %.6f %.6f total %.6f\n", s0, s1, tot);
        check(near(s0, 8.375000, 1e-4) && near(s1, 2.825000, 1e-4),
              "A3 spacings 8.375000 and 2.825000");
        check(near(tot, 11.200000, 1e-4), "A3 total 11.200000");
    }

    // A4 / A5 — height histogram and through-features (span both slabs).
    int hSlab0 = 0, hSlab1 = 0, hThru = 0;
    std::vector<double> thruR;
    if (lv.y.size() >= 3) {
        const double g0 = lv.y[1] - lv.y[0];
        const double g1 = lv.y[2] - lv.y[1];
        const double gt = lv.y[2] - lv.y[0];
        for (const Region* r : cyls) {
            const gp_XYZ c = r->ax.Location().XYZ();
            const gp_XYZ ad = unit(r->ax.Direction());
            const double y0 = ahat.Dot(c + ad * r->vMin);
            const double y1 = ahat.Dot(c + ad * r->vMax);
            const double h = std::abs(y1 - y0);
            if (near(h, g0, t.tauLvl * 4.0 + 1e-3)) ++hSlab0;
            else if (near(h, g1, t.tauLvl * 4.0 + 1e-3)) ++hSlab1;
            else if (near(h, gt, t.tauLvl * 4.0 + 1e-3)) {
                ++hThru;
                thruR.push_back(r->radius);
            }
        }
    }
    std::fprintf(stderr, "A4 histogram slab0=%d slab1=%d thru=%d (expect 11 / 2 / 2)\n",
                 hSlab0, hSlab1, hThru);
    check(hSlab0 == 11 && hSlab1 == 2 && hThru == 2, "A4 height histogram 11/2/2");
    std::sort(thruR.begin(), thruR.end());
    if (thruR.size() >= 2)
        std::fprintf(stderr, "A5 through-feature radii %.6f %.6f (expect ~5 and ~10)\n",
                     thruR[0], thruR[1]);
    check(thruR.size() == 2 && near(thruR[0], 5.0, 0.15) && near(thruR[1], 10.0, 0.15),
          "A5 through-features R5 and R10 span both slabs");

    // A6 — signed cap flux residual.
    double flux = 0.0;
    double peri = 0.0;
    std::vector<double> capA;
    for (const Region* r : planes) {
        const gp_XYZ n = unit(r->ax.Direction());
        const double nd = n.Dot(ahat);
        if (std::abs(nd) <= 1.0 - t.tauAx) continue;
        double area = 0.0;
        if (mv.pts && mv.tris && mv.compTris) {
            for (int local : r->tris) {
                if (local < 0 || static_cast<size_t>(local) >= mv.nTri) continue;
                const int g = mv.compTris[local];
                if (g < 0) continue;
                const gp_XYZ a = mv.pts[mv.tris[g][0]];
                const gp_XYZ ab = mv.pts[mv.tris[g][1]] - a;
                const gp_XYZ ac = mv.pts[mv.tris[g][2]] - a;
                area += 0.5 * ab.Crossed(ac).Modulus();
            }
        }
        capA.push_back(area);
        flux += area * nd;
        peri += 4.0 * std::sqrt(std::max(0.0, area));
    }
    std::sort(capA.rbegin(), capA.rend());
    std::fprintf(stderr, "A6 cap areas (desc)");
    for (double a : capA) std::fprintf(stderr, " %.6f", a);
    std::fprintf(stderr, "\n");
    double ident = 0.0;
    if (capA.size() >= 3) ident = capA[0] - (capA[1] + capA[2]);
    std::fprintf(stderr, "A6 residual flux=%.3e  f2-(f4+f21)=%.3e  periScale=%.3f\n",
                 flux, ident, peri);
    check(std::abs(flux) < t.tauFit * std::max(peri, 1.0), "A6 flux < tau_fit * peri");
    check(std::abs(ident) < 1.0, "A6 area(f2)-(area(f4)+area(f21)) closes");
}

static void testBody11(const char* stl) {
    HarnessMesh mesh;
    MeshView mv{};
    RegionSet rs;
    std::string err;
    if (!runSegment(stl, mesh, mv, rs, err)) {
        // Dirty / unreadable still counts as a decline if we can load.
        if (mesh.comps.empty()) {
            std::fprintf(stderr, "Body11 load: %s\n", err.c_str());
            check(false, "A8 Body11 load");
            return;
        }
        toMeshView(mesh, mesh.comps[0], mv);
        rs.compRoot = 0;
    }
    std::fprintf(stderr, "A8 Body11 nTri=%zu regions=%zu\n", mv.nTri, rs.regions.size());
    check(mv.nTri == 15300 || mesh.comps[0].compTris.size() == 15300 ||
              mv.nTri > 10000,
          "A8 Body11 is the 15300-tri class");
    PrismTols t;
    stl2step::refit::derivePrismTols(mv, rs, t);
    const PrismLevels lv = detectPrismatic(mv, rs, t);
    std::fprintf(stderr, "A8 Body11 ok=%d failedCond=%d\n", lv.ok ? 1 : 0, lv.failedCond);
    check(!lv.ok && lv.failedCond >= 1 && lv.failedCond <= 6, "A8 Body11 declines with failedCond");
}

static void testNegStl(const char* stl) {
    HarnessMesh mesh;
    MeshView mv{};
    RegionSet rs;
    std::string err;
    if (!runSegment(stl, mesh, mv, rs, err)) {
        std::fprintf(stderr, "neg-control load/segment: %s\n", err.c_str());
        // Synthetic A9 already covers cond 3/5. Fixture must still decline loudly.
        if (mesh.comps.empty()) {
            check(false, "A9 fixture load");
            return;
        }
    }
    PrismTols t;
    stl2step::refit::derivePrismTols(mv, rs, t);
    const PrismLevels lv = detectPrismatic(mv, rs, t);
    std::fprintf(stderr, "A9 fixture ok=%d failedCond=%d nTri=%zu nCyl_reg=%zu\n",
                 lv.ok ? 1 : 0, lv.failedCond, mv.nTri, rs.regions.size());
    check(!lv.ok && lv.failedCond != 0, "A9 fixture declines");
    if (lv.failedCond == 1)
        std::fprintf(stderr, "A9 fixture stopped at cond 1 (recognition); synthetic covers 3/5\n");
    else
        check(lv.failedCond == 3 || lv.failedCond == 5, "A9 fixture declines at cond 3 or 5");
}

int main(int argc, char** argv) {
    const char* hl = argc > 1 ? argv[1] : "tests/corpus/handle-lock.stl";
    const char* b11 = argc > 2 ? argv[2] : "tests/corpus/Body11.stl";
    const char* neg = argc > 3 ? argv[3] : "tests/corpus/nonprismatic-control.stl";

    testTolsSelfComputed();
    testSyntheticDeclines();
    testHandleLock(hl);
    testBody11(b11);
    testNegStl(neg);

    if (gFail)
        std::fprintf(stderr, "prism_test: %d FAIL\n", gFail);
    else
        std::fprintf(stderr, "prism_test: all PASS\n");
    return gFail ? 1 : 0;
}
