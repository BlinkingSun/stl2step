// D-gid-8.13. A cone/cylinder pair between certified coaxial rims shares its
// junction circle and one seam half-plane. The stagger predicates
// (rimSeamAzimuth, coneRimsStagger), then a synthetic plate with a chamfered
// through-hole whose chamfer rims carry a vertex count and a phase different from
// the bore's rims, converted end to end: the bore and both cones ship, every
// vertex and edge on them is within tau, the STEP volume equals the mesh volume
// plus the predicted change within (hole area) * tau, and both cones mesh. A
// control hole whose chamfer rims align with the bore's ships the same way.
// Synthetic geometry only.
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include "posix_compat.hpp"
#include "refit_axis_seam.hpp"
#include "stl_quant.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

// ---- predicates ------------------------------------------------------------------

static void casePredicates() {
    const gp_Ax1 axis(gp_Pnt(1.0, 2.0, -0.5), gp_Dir(0.0, 0.0, 1.0));
    const gp_Dir refX(1.0, 0.0, 0.0);
    auto at = [&](double r, double a, double z) {
        return gp_Pnt(1.0 + r * std::cos(a), 2.0 + r * std::sin(a), z);
    };
    std::vector<gp_Pnt> rim = {at(5.0, 0.3, 1.0), at(5.0, -0.1, 1.0), at(5.0, 2.0, 1.0)};
    double phi = 99.0;
    expect(rimSeamAzimuth(axis, refX, rim, phi) && std::fabs(phi + 0.1) < 1e-12,
           "azimuth_nearest_signed");
    const gp_Dir turned(std::cos(0.2), std::sin(0.2), 0.0);
    expect(rimSeamAzimuth(axis, turned, rim, phi) && std::fabs(phi - 0.1) < 1e-12,
           "azimuth_from_half_plane");
    std::vector<gp_Pnt> onAxis = {gp_Pnt(1.0, 2.0, 3.0)};
    expect(!rimSeamAzimuth(axis, refX, onAxis, phi), "azimuth_refuses_axis_point");
    expect(!rimSeamAzimuth(axis, refX, {}, phi), "azimuth_refuses_empty_rim");

    const double tau = 2.0e-5, rLo = 4.0, lim = tau / rLo;
    expect(!coneRimsStagger(0.1, 0.1 + 0.5 * lim, rLo, tau), "stagger_within_rim_start_tol");
    expect(coneRimsStagger(0.1, 0.1 + 2.0 * lim, rLo, tau), "stagger_beyond_rim_start_tol");
    expect(!coneRimsStagger(kPi - 0.25 * lim, -kPi + 0.25 * lim, rLo, tau),
           "stagger_wraps_across_pi");
    expect(coneRimsStagger(kPi - lim, -kPi + lim, rLo, tau), "stagger_wrapped_beyond");
    expect(!coneRimsStagger(0.0, 1.0, rLo, 0.0), "stagger_needs_tau");
    expect(!coneRimsStagger(0.0, 1.0, 0.0, tau), "stagger_needs_radius");
}

// ---- synthetic chamfered through-hole ------------------------------------------

struct P3 {
    double x, y, z, a;
};
using Tri = std::array<P3, 3>;

// A square plate (side 2*half, height H) whose centre (cx, cy) carries a bore of
// radius Rb with nb sides and two 45 deg chamfers of size C whose outer rims carry
// nc vertices starting at phase phDeg.
static std::vector<Tri> chamferedHole(int nb, int nc, double phDeg, double cx, double cy) {
    const double H = 6.0, Rb = 3.0, C = 0.75, half = 12.0, Rc = Rb + C;
    const double twoPi = 2.0 * kPi;
    auto wrap = [&](double a) {
        a = std::fmod(a, twoPi);
        return a < 0.0 ? a + twoPi : a;
    };
    auto ring = [&](int n, double r, double z, double ph) {
        std::vector<P3> v;
        for (int k = 0; k < n; ++k) {
            const double a = ph + twoPi * k / n;
            v.push_back({cx + r * std::cos(a), cy + r * std::sin(a), z, wrap(a)});
        }
        std::sort(v.begin(), v.end(), [](const P3& p, const P3& q) { return p.a < q.a; });
        return v;
    };
    auto sq = [&](double a, double z) {
        const double c = std::cos(a), s = std::sin(a);
        const double t = half / std::max(std::fabs(c), std::fabs(s));
        return P3{cx + t * c, cy + t * s, z, wrap(a)};
    };
    std::vector<Tri> tris;
    auto tri = [&](const P3& a, const P3& b, const P3& c, bool flip) {
        if (flip) tris.push_back({a, c, b});
        else tris.push_back({a, b, c});
    };
    // Band between two closed rings sorted by azimuth, advancing on the smaller next
    // azimuth. flip = false gives a normal toward the axis side of the outer ring
    // turned down (-z for a horizontal band).
    auto zipper = [&](const std::vector<P3>& in, const std::vector<P3>& out, bool flip) {
        const size_t ni = in.size(), no = out.size();
        const double base = std::min(in[0].a, out[0].a);
        auto unroll = [&](const std::vector<P3>& r) {
            std::vector<double> u;
            for (const P3& p : r) {
                double a = p.a;
                while (a < base - 1e-12) a += twoPi;
                if (!u.empty())
                    while (a < u.back()) a += twoPi;
                u.push_back(a);
            }
            u.push_back(u.front() + twoPi);
            return u;
        };
        const std::vector<double> ia = unroll(in), oa = unroll(out);
        size_t i = 0, j = 0;
        while (i < ni || j < no) {
            if (j >= no || (i < ni && ia[i + 1] <= oa[j + 1])) {
                tri(in[i % ni], in[(i + 1) % ni], out[j % no], flip);
                ++i;
            } else {
                tri(in[i % ni], out[(j + 1) % no], out[j % no], flip);
                ++j;
            }
        }
    };
    const double ph = phDeg * kPi / 180.0;
    const std::vector<P3> topC = ring(nc, Rc, H, ph), botC = ring(nc, Rc, 0.0, ph);
    const std::vector<P3> topB = ring(nb, Rb, H - C, 0.0), botB = ring(nb, Rb, C, 0.0);
    std::vector<double> outA;
    for (const P3& p : topC) outA.push_back(p.a);
    for (double d : {45.0, 135.0, 225.0, 315.0}) outA.push_back(d * kPi / 180.0);
    std::sort(outA.begin(), outA.end());
    outA.erase(std::unique(outA.begin(), outA.end(),
                           [](double a, double b) { return std::fabs(a - b) < 1e-12; }),
               outA.end());
    std::vector<P3> topS, botS;
    for (double a : outA) {
        topS.push_back(sq(a, H));
        botS.push_back(sq(a, 0.0));
    }
    zipper(topC, topS, true);
    zipper(botC, botS, false);
    for (size_t k = 0; k < topS.size(); ++k) {
        const size_t k1 = (k + 1) % topS.size();
        tri(botS[k], botS[k1], topS[k1], false);
        tri(botS[k], topS[k1], topS[k], false);
    }
    zipper(topB, topC, true);
    zipper(botB, botC, false);
    for (int k = 0; k < nb; ++k) {
        const int k1 = (k + 1) % nb;
        tri(botB[(size_t)k], topB[(size_t)k1], botB[(size_t)k1], false);
        tri(botB[(size_t)k], topB[(size_t)k], topB[(size_t)k1], false);
    }
    return tris;
}

static bool writeStl(const std::string& path, const std::vector<Tri>& tris) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    char hdr[80];
    std::memset(hdr, 0, sizeof(hdr));
    std::snprintf(hdr, sizeof(hdr), "synthetic chamfered hole");
    std::fwrite(hdr, 1, 80, f);
    const std::uint32_t n = (std::uint32_t)tris.size();
    std::fwrite(&n, 4, 1, f);
    for (const Tri& t : tris) {
        const double ux = t[1].x - t[0].x, uy = t[1].y - t[0].y, uz = t[1].z - t[0].z;
        const double vx = t[2].x - t[0].x, vy = t[2].y - t[0].y, vz = t[2].z - t[0].z;
        double nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
        const double m = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (m > 0.0) {
            nx /= m;
            ny /= m;
            nz /= m;
        }
        float rec[12] = {(float)nx, (float)ny, (float)nz};
        for (int k = 0; k < 3; ++k) {
            rec[3 + 3 * k] = (float)t[(size_t)k].x;
            rec[4 + 3 * k] = (float)t[(size_t)k].y;
            rec[5 + 3 * k] = (float)t[(size_t)k].z;
        }
        std::fwrite(rec, 4, 12, f);
        const std::uint16_t attr = 0;
        std::fwrite(&attr, 2, 1, f);
    }
    std::fclose(f);
    return true;
}

static void caseHole(const char* tag, int nb, int nc, double phDeg) {
    const std::string stl = std::string("cone-bore-seam-") + tag + ".stl";
    const std::string step = std::string("cone-bore-seam-") + tag + ".step";
    const bool wrote = writeStl(stl, chamferedHole(nb, nc, phDeg, 7.5, -3.25));
    expect(wrote, (std::string(tag) + "_stl_written").c_str());
    if (!wrote) return;
    const double tau = 2.0 * stl2step::stlQuantFloor(stl).q;
    expect(tau > 0.0, (std::string(tag) + "_tau").c_str());

    stl2step::Options opt;
    opt.input = stl;
    opt.output = step;
    opt.smooth = true;
    opt.verify = true;
    opt.threads = 1;
    const stl2step::Result r = stl2step::convert(opt);
    const std::string t(tag);
    expect(r.ok && r.solids == 1 && r.smoothRevertedComponents == 0, (t + "_convert").c_str());
    expect(r.smoothBuiltCylinders == 1 && r.smoothBuiltCones == 2 && r.smoothFacetFaces == 0,
           (t + "_bore_and_two_cones").c_str());

    STEPControl_Reader rd;
    const bool read = rd.ReadFile(step.c_str()) == IFSelect_RetDone;
    expect(read, (t + "_step_read").c_str());
    if (!read) {
        std::remove(stl.c_str());
        return;
    }
    rd.TransferRoots();
    const TopoDS_Shape s = rd.OneShape();
    int nSolid = 0;
    for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next()) ++nSolid;
    expect(nSolid == 1 && BRepCheck_Analyzer(s).IsValid(), (t + "_reread_valid").c_str());

    TopTools_IndexedMapOfShape fm;
    TopExp::MapShapes(s, TopAbs_FACE, fm);
    double worstV = 0.0, worstE = 0.0, holeArea = 0.0;
    int nCone = 0, nCyl = 0;
    std::vector<int> cones;
    for (int i = 1; i <= fm.Extent(); ++i) {
        const TopoDS_Face f = TopoDS::Face(fm(i));
        const GeomAbs_SurfaceType ty = BRepAdaptor_Surface(f, Standard_False).GetType();
        if (ty != GeomAbs_Cone && ty != GeomAbs_Cylinder) continue;
        if (ty == GeomAbs_Cone) {
            ++nCone;
            cones.push_back(i);
        } else {
            ++nCyl;
        }
        GProp_GProps sp;
        BRepGProp::SurfaceProperties(f, sp);
        holeArea += sp.Mass();
        for (TopExp_Explorer ev(f, TopAbs_VERTEX); ev.More(); ev.Next())
            worstV = std::max(worstV, BRep_Tool::Tolerance(TopoDS::Vertex(ev.Current())));
        for (TopExp_Explorer ee(f, TopAbs_EDGE); ee.More(); ee.Next())
            worstE = std::max(worstE, BRep_Tool::Tolerance(TopoDS::Edge(ee.Current())));
    }
    std::fprintf(stderr, "INFO %s tau=%.6g worstVertexTol=%.6g worstEdgeTol=%.6g\n", tag, tau,
                 worstV, worstE);
    expect(nCone == 2 && nCyl == 1, (t + "_step_faces").c_str());
    expect(worstV <= tau, (t + "_vertex_tol_within_tau").c_str());
    expect(worstE <= tau, (t + "_edge_tol_within_tau").c_str());
    const double gap = r.stepVolumeMM3 - (r.meshVolumeMM3 + r.smoothVolPredictedMM3);
    std::fprintf(stderr, "INFO %s step=%.6f mesh=%.6f pred=%.6f gap=%.3g bound=%.3g\n", tag,
                 r.stepVolumeMM3, r.meshVolumeMM3, r.smoothVolPredictedMM3, gap, holeArea * tau);
    expect(std::fabs(gap) <= holeArea * tau, (t + "_volume_matches_prediction").c_str());
    for (double defl : {1e-2, 1e-3, 1e-4}) {
        STEPControl_Reader rm;
        rm.ReadFile(step.c_str());
        rm.TransferRoots();
        TopoDS_Shape m = rm.OneShape();
        BRepMesh_IncrementalMesh inc(m, defl, Standard_False, 0.5, Standard_False);
        TopTools_IndexedMapOfShape fm2;
        TopExp::MapShapes(m, TopAbs_FACE, fm2);
        int bad = 0;
        for (int i = 1; i <= fm2.Extent(); ++i) {
            TopLoc_Location L;
            Handle(Poly_Triangulation) T = BRep_Tool::Triangulation(TopoDS::Face(fm2(i)), L);
            if (T.IsNull() || T->NbTriangles() == 0) ++bad;
        }
        char nm[96];
        std::snprintf(nm, sizeof(nm), "%s_meshes_at_%g", tag, defl);
        expect(bad == 0, nm);
    }
    std::remove(stl.c_str());
    std::remove(step.c_str());
}

int main() {
    casePredicates();
    // Chamfer rims of 24 vertices at a 4 deg phase never meet a 40-vertex bore
    // column: every seam corner would stagger.
    caseHole("staggered", 40, 24, 4.0);
    // Control: chamfer rims on the bore's own columns.
    caseHole("aligned", 40, 40, 0.0);
    std::fprintf(stderr, "%s cone_bore_seam_unit (%d failure%s)\n", gFail ? "FAIL" : "PASS", gFail,
                 gFail == 1 ? "" : "s");
    return gFail ? 1 : 0;
}
