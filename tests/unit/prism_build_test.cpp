// AC3-P3 — exercise buildPrismSolid / Stage-P routing (SPEC C1-C10, C16).
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "refit_prism.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GProp_GProps.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_XYZ.hxx>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::MeshView;
using stl2step::refit::PrismLevels;
using stl2step::refit::PrismTols;
using stl2step::refit::ProfLoop;
using stl2step::refit::ProfSeg;
using stl2step::refit::Profile;
using stl2step::refit::Region;
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;
using stl2step::refit::buildFaces;
using stl2step::refit::buildPrismSolid;
using stl2step::refit::detectPrismatic;
using stl2step::refit::fitProfile;
using stl2step::refit::segment;
using stl2step::refit::sliceProfiles;

namespace stl2step {
namespace refit {
void derivePrismTols(const MeshView& mv, const RegionSet& rs, PrismTols& t);
void prismBindSketchOrigin(const gp_XYZ& o);
void prismResetStageFlags();
bool prismStagePUsed();
int prismPlatePathHits();
int prismLastReverted();
}  // namespace refit
}  // namespace stl2step

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
    return std::fabs(a - b) <= tol;
}

struct Census {
    int faces = 0, planes = 0, cyls = 0, lines = 0, circles = 0, otherE = 0;
    double vol = 0.0;
    bool closed = false;
    bool valid = false;
    std::vector<double> radii;
    std::vector<double> capAreas;
    std::vector<double> cylHeights;
};

static Census censusOf(const TopoDS_Shape& s) {
    Census c;
    if (s.IsNull()) return c;
    try {
        GProp_GProps gp;
        BRepGProp::VolumeProperties(s, gp);
        c.vol = gp.Mass();
    } catch (...) {
    }
    try {
        BRepCheck_Analyzer an(s, Standard_True);
        c.valid = an.IsValid() == Standard_True;
    } catch (...) {
        c.valid = false;
    }
    BRep_Builder B;
    TopoDS_Shell sh;
    B.MakeShell(sh);
    for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) {
        const TopoDS_Face f = TopoDS::Face(fx.Current());
        ++c.faces;
        B.Add(sh, f);
        try {
            BRepAdaptor_Surface sa(f, Standard_False);
            if (sa.GetType() == GeomAbs_Plane) {
                ++c.planes;
                GProp_GProps sp;
                BRepGProp::SurfaceProperties(f, sp);
                const gp_Dir n = sa.Plane().Axis().Direction();
                if (std::fabs(n.Y()) > 0.9) c.capAreas.push_back(sp.Mass());
            } else if (sa.GetType() == GeomAbs_Cylinder) {
                ++c.cyls;
                c.radii.push_back(sa.Cylinder().Radius());
                double umin = 0, umax = 0, vmin = 0, vmax = 0;
                BRepTools::UVBounds(f, umin, umax, vmin, vmax);
                c.cylHeights.push_back(std::fabs(vmax - vmin));
            }
        } catch (...) {
        }
    }
    c.closed = BRep_Tool::IsClosed(sh);
    for (TopExp_Explorer ex(s, TopAbs_EDGE); ex.More(); ex.Next()) {
        try {
            BRepAdaptor_Curve ca(TopoDS::Edge(ex.Current()));
            if (ca.GetType() == GeomAbs_Line) ++c.lines;
            else if (ca.GetType() == GeomAbs_Circle) ++c.circles;
            else ++c.otherE;
        } catch (...) {
            ++c.otherE;
        }
    }
    std::sort(c.capAreas.begin(), c.capAreas.end(), std::greater<double>());
    std::sort(c.radii.begin(), c.radii.end());
    return c;
}

static bool hasRadius(const std::vector<double>& rs, double R, double rel) {
    for (double v : rs)
        if (R > 0.0 && std::fabs(v - R) / R <= rel) return true;
    return false;
}

static bool hasHeight(const std::vector<double>& hs, double h, double tol) {
    for (double v : hs)
        if (std::fabs(v - h) <= tol) return true;
    return false;
}

static Profile squareProfile(int slab, double x0, double y0, double x1, double y1) {
    Profile p;
    p.slab = slab;
    ProfLoop lp;
    lp.outer = true;
    const double xs[4] = {x0, x1, x1, x0};
    const double ys[4] = {y0, y0, y1, y1};
    for (int i = 0; i < 4; ++i) {
        ProfSeg s;
        s.a = gp_Pnt2d(xs[i], ys[i]);
        s.b = gp_Pnt2d(xs[(i + 1) % 4], ys[(i + 1) % 4]);
        lp.segs.push_back(s);
    }
    lp.area = std::fabs((x1 - x0) * (y1 - y0));
    p.loops.push_back(lp);
    return p;
}

static void addCircleHole(Profile& p, double cx, double cy, double R) {
    ProfLoop lp;
    lp.outer = false;
    lp.area = std::acos(-1.0) * R * R;
    ProfSeg s;
    s.isArc = true;
    s.center = gp_Pnt2d(cx, cy);
    s.R = R;
    s.phi = 2.0 * std::acos(-1.0);
    s.ccw = true;
    s.a = gp_Pnt2d(cx + R, cy);
    s.b = s.a;
    lp.segs.push_back(s);
    p.loops.push_back(lp);
}

static bool runSegment(const char* stl, HarnessMesh& mesh, MeshView& mv, RegionSet& rs,
                       std::string& err) {
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) return false;
    size_t pick = 0;
    for (size_t i = 1; i < mesh.comps.size(); ++i)
        if (mesh.comps[i].compTris.size() > mesh.comps[pick].compTris.size()) pick = i;
    toMeshView(mesh, mesh.comps[pick], mv);
    rs = RegionSet{};
    rs.compRoot = 0;
    SegmentParams sp;
    return segment(mv, sp, rs, nullptr);
}

int main(int argc, char** argv) {
    const char* stl = argc > 1 ? argv[1] : "tests/corpus/handle-lock.stl";

    // Synthetic API: two stacked squares (union of axially disjoint prisms).
    {
        stl2step::refit::prismResetStageFlags();
        stl2step::refit::prismBindSketchOrigin(gp_XYZ(0, 0, 0));
        PrismLevels lv;
        lv.ok = true;
        lv.axis = gp_Dir(0, 0, 1);
        lv.y = {0.0, 4.0, 6.0};
        lv.capRegion = {0, 1, 2};
        std::vector<Profile> profs;
        profs.push_back(squareProfile(0, 0.0, 0.0, 10.0, 10.0));
        profs.push_back(squareProfile(1, 2.0, 2.0, 8.0, 8.0));
        addCircleHole(profs[0], 5.0, 5.0, 1.0);
        addCircleHole(profs[1], 5.0, 5.0, 1.0);
        TopoDS_Shape sh;
        check(buildPrismSolid(profs, lv, sh), "api_build_two_slabs");
        const Census c = censusOf(sh);
        check(c.valid, "api_valid");
        check(c.closed, "api_closed");
        const double hole = std::acos(-1.0) * 1.0 * 1.0 * 6.0;
        const double expect = 10.0 * 10.0 * 4.0 + 6.0 * 6.0 * 2.0 - hole;
        check(near(c.vol, expect, 0.05), "api_volume");
        check(c.cyls >= 1, "api_has_cyl");
        std::fprintf(stderr, "api census faces=%d planes=%d cyls=%d vol=%.6f expect=%.6f\n",
                     c.faces, c.planes, c.cyls, c.vol, expect);
    }

    // Bad profile: missing outer -> refuse.
    {
        PrismLevels lv;
        lv.ok = true;
        lv.axis = gp_Dir(0, 0, 1);
        lv.y = {0.0, 2.0};
        Profile bad;
        bad.slab = 0;
        TopoDS_Shape sh;
        check(!buildPrismSolid({bad}, lv, sh), "api_refuse_empty");
    }

    HarnessMesh mesh;
    MeshView mv{};
    RegionSet rs;
    std::string err;
    if (!runSegment(stl, mesh, mv, rs, err)) {
        std::fprintf(stderr, "segment failed: %s\n", err.c_str());
        return 2;
    }

    PrismTols tol;
    stl2step::refit::derivePrismTols(mv, rs, tol);
    const PrismLevels lv = detectPrismatic(mv, rs, tol);
    check(lv.ok, "detect_ok");
    check(lv.y.size() == 3, "C3_three_levels");
    if (lv.y.size() >= 3) {
        const double h0 = lv.y[1] - lv.y[0];
        const double h1 = lv.y[2] - lv.y[1];
        const double ht = lv.y[2] - lv.y[0];
        check(near(h0, 8.375, 1e-3), "C3_spacing0");
        check(near(h1, 2.825, 1e-3), "C3_spacing1");
        check(near(ht, 11.2, 1e-3), "C3_total");
        std::fprintf(stderr, "C3 levels y=%.6f %.6f %.6f  h=%.6f %.6f tot=%.6f\n", lv.y[0],
                     lv.y[1], lv.y[2], h0, h1, ht);
    }

    std::vector<Profile> profs;
    check(sliceProfiles(mv, rs, lv, tol, profs), "sliceProfiles");
    int nDecl = 0;
    for (Profile& p : profs) {
        int d = 0;
        check(fitProfile(mv, tol, p, d), "fitProfile");
        nDecl += d;
    }
    check(nDecl == 0, "fit_nDecl");

    gp_XYZ origin(0, 0, 0);
    if (!lv.capRegion.empty()) {
        for (const Region& r : rs.regions)
            if (r.id == lv.capRegion[0]) origin = r.ax.Location().XYZ();
    }
    stl2step::refit::prismBindSketchOrigin(origin);

    for (const Profile& p : profs) {
        std::fprintf(stderr, "prof slab=%d loops=%zu\n", p.slab, p.loops.size());
        for (size_t i = 0; i < p.loops.size(); ++i) {
            const ProfLoop& lp = p.loops[i];
            int nA = 0, nL = 0;
            for (const ProfSeg& s : lp.segs) {
                if (s.isArc) ++nA;
                else ++nL;
            }
            std::fprintf(stderr, "  loop %zu outer=%d segs=%zu line=%d arc=%d area=%.3f\n", i,
                         (int)lp.outer, lp.segs.size(), nL, nA, lp.area);
        }
    }

    TopoDS_Shape solid;
    {
        double dVolAbs = 0.0;
        for (const Region& r : rs.regions) dVolAbs += std::fabs(r.dVolPredicted);
        std::fprintf(stderr, "dVolAbs=%.6f budget=%.6f\n", dVolAbs,
                     std::max(1e-4 * 15868.884516, 3.0 * dVolAbs));
    }
    check(buildPrismSolid(profs, lv, solid), "buildPrismSolid");
    const Census c = censusOf(solid);
    {
        std::vector<std::array<double, 4>> planes;
        for (TopExp_Explorer fx(solid, TopAbs_FACE); fx.More(); fx.Next()) {
            const TopoDS_Face f = TopoDS::Face(fx.Current());
            try {
                BRepAdaptor_Surface sa(f, Standard_False);
                if (sa.GetType() != GeomAbs_Plane) continue;
                const gp_Pln pl = sa.Plane();
                const gp_Dir n = pl.Axis().Direction();
                const double d = pl.Location().XYZ().Dot(n.XYZ());
                planes.push_back({n.X(), n.Y(), n.Z(), d});
            } catch (...) {
            }
        }
        int uniq = 0;
        std::vector<char> used(planes.size(), 0);
        for (size_t i = 0; i < planes.size(); ++i) {
            if (used[i]) continue;
            ++uniq;
            for (size_t j = i + 1; j < planes.size(); ++j) {
                if (used[j]) continue;
                const double dn = std::fabs(planes[i][0] * planes[j][0] +
                                            planes[i][1] * planes[j][1] +
                                            planes[i][2] * planes[j][2]);
                const double dd = std::fabs(planes[i][3] - planes[j][3]);
                if (dn > 0.999 && dd < 0.05) used[j] = 1;
            }
        }
        std::fprintf(stderr, "C1 unique-planes~%d of %d\n", uniq, (int)planes.size());
    }
    std::fprintf(stderr,
                 "C1 census faces=%d planes=%d cyls=%d lines=%d circles=%d otherE=%d "
                 "vol=%.6f closed=%d valid=%d\n",
                 c.faces, c.planes, c.cyls, c.lines, c.circles, c.otherE, c.vol,
                 (int)c.closed, (int)c.valid);
    check(c.faces == 28 && c.planes == 13 && c.cyls == 15, "C1_face_census");
    check(c.faces != 26, "C1_not_26");

    const double truthR[] = {0.5, 3.7872, 5.0, 5.75, 6.0, 9.0, 10.0, 15.0, 16.0, 20.0, 30.0};
    int nHit = 0;
    for (double R : truthR)
        if (hasRadius(c.radii, R, 0.003)) ++nHit;
    check(nHit == 11, "C2_distinct_radii");
    check(hasRadius(c.radii, 16.0, 0.003), "C2_R16");
    check(hasRadius(c.radii, 20.0, 0.003), "C2_R20");
    check(hasRadius(c.radii, 30.0, 0.003), "C2_R30");
    std::fprintf(stderr, "C2 radii:");
    for (double r : c.radii) std::fprintf(stderr, " %.4f", r);
    std::fprintf(stderr, "\nC2 heights:");
    for (double h : c.cylHeights) std::fprintf(stderr, " %.4f", h);
    std::fprintf(stderr, "\n");

    if (c.capAreas.size() >= 3) {
        const double a0 = c.capAreas[0], a1 = c.capAreas[1], a2 = c.capAreas[2];
        const double resid = a0 - (a1 + a2);
        std::fprintf(stderr, "C4 caps %.6f %.6f %.6f resid=%.6e\n", a0, a1, a2, resid);
        check(std::fabs(resid) < 1e-3, "C4_cap_closure");
        check(near(a0, 1751.434630, 1.0), "C4_cap0");
        check(near(a1, 1326.457013, 1.0) || near(a1, 424.977617, 1.0), "C4_cap1");
    } else {
        check(false, "C4_three_caps");
    }

    check(hasHeight(c.cylHeights, 11.2, 0.05), "C5_through_height");
    check(hasRadius(c.radii, 5.0, 0.003), "C5_f5_R5");
    check(hasRadius(c.radii, 10.0, 0.003), "C5_f13_R10");

    const double meshVol = 15868.884516;
    const double dPct = meshVol > 0.0 ? 100.0 * std::fabs(c.vol - meshVol) / meshVol : -1.0;
    std::fprintf(stderr, "C6 vol=%.6f mesh=%.6f deltaPct=%.6f\n", c.vol, meshVol, dPct);
    check(dPct >= 0.0 && dPct < 1.0, "C6_volume");

    check(c.closed && c.valid, "C7_watertight_valid");
    check(c.otherE == 0, "C8_alphabet_only");
    check(c.lines == 48 && c.circles == 30, "C8_edge_counts");

    // C9 — buildFaces Stage-P path does not set BuiltAs (plate/cascade not entered).
    {
        stl2step::refit::prismResetStageFlags();
        std::vector<TopoDS_Vertex> verts(mv.nVtx);
        BRep_Builder bb;
        for (size_t i = 0; i < mv.nVtx; ++i) {
            const gp_XYZ p = mv.pts[mv.compVtx[i]];
            bb.MakeVertex(verts[i], gp_Pnt(p), 1e-7);
        }
        RegionSet rs2 = rs;
        std::vector<TopoDS_Face> faces;
        const bool ok = buildFaces(mv, rs2, verts, faces, nullptr);
        check(ok && !faces.empty(), "C9_buildFaces");
        check(stl2step::refit::prismStagePUsed(), "C9_stageP_used");
        check(stl2step::refit::prismPlatePathHits() == 0, "C9_plate_not_called");
        check(stl2step::refit::prismLastReverted() == 0, "C9_not_reverted");
        int nBuilt = 0;
        for (const Region& r : rs2.regions)
            if (r.builtAs != stl2step::refit::BuiltAs::NotBuilt) ++nBuilt;
        check(nBuilt == 0, "C9_builtAs_untouched");
        std::fprintf(stderr, "C9 faces=%zu stageP=%d plate=%d builtAs=%d\n", faces.size(),
                     (int)stl2step::refit::prismStagePUsed(),
                     stl2step::refit::prismPlatePathHits(), nBuilt);
    }

    // C10 — injected bad profile reverts (Stage P not consumed).
    {
        stl2step::refit::prismResetStageFlags();
        ::setenv("STL2STEP_PRISM_INJECT_BAD", "1", 1);
        std::vector<TopoDS_Vertex> verts(mv.nVtx);
        BRep_Builder bb;
        for (size_t i = 0; i < mv.nVtx; ++i) {
            const gp_XYZ p = mv.pts[mv.compVtx[i]];
            bb.MakeVertex(verts[i], gp_Pnt(p), 1e-7);
        }
        RegionSet rs3 = rs;
        std::vector<TopoDS_Face> faces;
        const bool ok = buildFaces(mv, rs3, verts, faces, nullptr);
        ::unsetenv("STL2STEP_PRISM_INJECT_BAD");
        check(ok || !ok, "C10_legacy_ran");
        check(!stl2step::refit::prismStagePUsed(), "C10_stageP_not_used");
        check(stl2step::refit::prismLastReverted() == 1, "C10_reverted_flag");
        std::fprintf(stderr, "C10 ok=%d stageP=%d reverted=%d faces=%zu\n", (int)ok,
                     (int)stl2step::refit::prismStagePUsed(),
                     stl2step::refit::prismLastReverted(), faces.size());
    }

    // Public convert() API on the fixture.
    {
        stl2step::Options opt;
        opt.input = stl;
        opt.output = "/tmp/p3-r3-api.step";
        opt.smooth = true;
        opt.verify = true;
        const stl2step::Result r = stl2step::convert(opt);
        std::fprintf(stderr,
                     "convert ok=%d wt=%d open=%d volΔ=%.6f mesh=%.6f step=%.6f warn=%zu\n",
                     (int)r.ok, (int)r.watertight, r.openShells, r.volumeDeltaPct,
                     r.meshVolumeMM3, r.stepVolumeMM3, r.warnings.size());
        check(r.ok, "api_convert_ok");
        check(r.watertight && r.openShells == 0, "api_convert_closed");
    }

    std::fprintf(stderr, "SUMMARY fail=%d\n", gFail);
    return gFail ? 4 : 0;
}
