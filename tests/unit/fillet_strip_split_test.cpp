// D-gid-4.1. A mixed FilletStrip (one loop, two mesh-disjoint pieces on the
// strip's existing cylinder) builds one analytic face per piece that certifies.
// A piece whose vertices are outside τ facets alone. Synthetic mesh only.
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "posix_compat.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Builder.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::BuiltAs;
using stl2step::refit::Loop;
using stl2step::refit::LoopRole;
using stl2step::refit::MeshView;
using stl2step::refit::Origin;
using stl2step::refit::Region;
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;
using stl2step::refit::buildFaces;

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS %s\n", name);
    }
}

struct V {
    double x, y, z;
};

static V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V cross(V a, V b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// Bar along +Z. The +X/+Y corner is a quarter-cylinder of radius R on two
// axial pieces and a sharp corner between them, so the two fillets share no
// mesh edge. Vertices of each fillet lie on that cylinder.
static bool writeMesh(const std::string& path) {
    const double W = 20.0, H = 16.0, R = 4.0;
    const int nAng = 4;
    const double z1[] = {0, 5, 10, 15};
    const double zS[] = {15, 22};
    const double z2[] = {22, 27, 32, 37};
    const gp_XYZ center(W - R, H - R, 0);

    std::vector<V> fillet;
    fillet.push_back({0, 0, 0});
    fillet.push_back({W, 0, 0});
    for (int k = 0; k <= nAng; ++k) {
        const double th = (0.5 * std::acos(-1.0) * k) / nAng;
        fillet.push_back({center.X() + R * std::cos(th), center.Y() + R * std::sin(th), 0});
    }
    fillet.push_back({0, H, 0});

    std::vector<V> sharp;
    sharp.push_back({0, 0, 0});
    sharp.push_back({W, 0, 0});
    sharp.push_back({W, H - R, 0});
    sharp.push_back({W, H, 0});
    sharp.push_back({W - R, H, 0});
    sharp.push_back({0, H, 0});

    std::vector<std::array<V, 3>> tris;
    auto emit = [&](V a, V b, V c, V pref) {
        if (dot(cross(b - a, c - a), pref) < 0) std::swap(b, c);
        tris.push_back({a, b, c});
    };
    auto at = [](V p, double z) {
        p.z = z;
        return p;
    };
    auto stitch = [&](const std::vector<V>& ring, double z0, double z1) {
        const int n = (int)ring.size();
        for (int k = 0; k < n; ++k) {
            const V a = at(ring[k], z0);
            const V b = at(ring[(k + 1) % n], z0);
            const V c = at(ring[(k + 1) % n], z1);
            const V d = at(ring[k], z1);
            const V mid{(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, 0};
            emit(a, b, c, {mid.x - 0.5 * W, mid.y - 0.5 * H, 0});
            emit(a, c, d, {mid.x - 0.5 * W, mid.y - 0.5 * H, 0});
        }
    };
    auto cap = [&](const std::vector<V>& ring, double z, double nz) {
        V c{0, 0, z};
        for (const V& p : ring) {
            c.x += p.x;
            c.y += p.y;
        }
        c.x /= (double)ring.size();
        c.y /= (double)ring.size();
        const int n = (int)ring.size();
        for (int k = 0; k < n; ++k)
            emit(at(ring[k], z), at(ring[(k + 1) % n], z), c, {0, 0, nz});
    };
    auto lid = [&](double z, double nz) {
        std::vector<V> arc;
        for (int k = 0; k <= nAng; ++k) {
            const double th = (0.5 * std::acos(-1.0) * k) / nAng;
            arc.push_back({center.X() + R * std::cos(th), center.Y() + R * std::sin(th), z});
        }
        const V s{W, H, z};
        for (int k = 0; k < nAng; ++k) emit(arc[k], arc[k + 1], s, {0, 0, nz});
    };

    cap(fillet, z1[0], -1);
    for (int i = 0; i + 1 < 4; ++i) stitch(fillet, z1[i], z1[i + 1]);
    lid(z1[3], -1);
    stitch(sharp, zS[0], zS[1]);
    lid(z2[0], +1);
    for (int i = 0; i + 1 < 4; ++i) stitch(fillet, z2[i], z2[i + 1]);
    cap(fillet, z2[3], +1);

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "solid mixedstrip\n");
    for (const auto& t : tris) {
        const V n = cross(t[1] - t[0], t[2] - t[0]);
        std::fprintf(f, "facet normal %.4f %.4f %.4f\n  outer loop\n", n.x, n.y, n.z);
        for (int k = 0; k < 3; ++k)
            std::fprintf(f, "    vertex %.4f %.4f %.4f\n", t[k].x, t[k].y, t[k].z);
        std::fprintf(f, "  endloop\nendfacet\n");
    }
    std::fprintf(f, "endsolid mixedstrip\n");
    std::fclose(f);
    std::fprintf(stderr, "mesh tris=%zu path=%s\n", tris.size(), path.c_str());
    return true;
}

struct FaceCensus {
    int cyl = 0, plane = 0, facet3 = 0, other = 0;
};

static FaceCensus censusOf(const std::vector<TopoDS_Face>& faces) {
    FaceCensus c;
    for (const TopoDS_Face& f : faces) {
        if (f.IsNull()) continue;
        int nE = 0;
        for (TopExp_Explorer ex(f, TopAbs_EDGE); ex.More(); ex.Next()) ++nE;
        BRepAdaptor_Surface as(f, Standard_False);
        if (as.GetType() == GeomAbs_Cylinder) ++c.cyl;
        else if (as.GetType() == GeomAbs_Plane && nE == 3) ++c.facet3;
        else if (as.GetType() == GeomAbs_Plane) ++c.plane;
        else ++c.other;
    }
    return c;
}

static const char* originName(Origin o) {
    switch (o) {
    case Origin::PlaneGrow: return "plane";
    case Origin::CylGrow: return "cyl";
    case Origin::FilletStrip: return "fillet";
    case Origin::NgonWall: return "ngon";
    case Origin::ChamferCone: return "chamfer";
    case Origin::TorusBlend: return "torus";
    case Origin::RuledQuad: return "ruled";
    case Origin::SphereCorner: return "sphere";
    }
    return "?";
}

static void dumpRegions(const RegionSet& rs) {
    for (const Region& r : rs.regions) {
        std::fprintf(stderr,
                     "  rid=%d origin=%s type=%d R=%.6f nTri=%zu loops=%zu dev=%.6g "
                     "built=%d\n",
                     r.id, originName(r.origin), (int)r.type, r.radius, r.tris.size(),
                     r.loops.size(), r.maxVertexDev, (int)r.builtAs);
    }
}

static bool prepare(const std::string& path, HarnessMesh& mesh, MeshView& mv, RegionSet& rs,
                    std::string& err) {
    if (!loadMesh(path, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) return false;
    if (mesh.comps.size() != 1 || !mesh.comps[0].clean()) {
        std::fprintf(stderr, "mesh comps=%zu open=%d conflict=%d nm=%d\n", mesh.comps.size(),
                     mesh.comps[0].open, mesh.comps[0].conflict, mesh.comps[0].nonManifold);
        return false;
    }
    toMeshView(mesh, mesh.comps[0], mv);
    rs = RegionSet{};
    SegmentParams p;
    if (!stl2step::refit::segment(mv, p, rs, nullptr)) return false;
    std::fprintf(stderr, "segment regions=%zu q=%.6g\n", rs.regions.size(), mv.quantFloor);
    dumpRegions(rs);
    return true;
}

// Concatenate B's outer-loop chains onto A's and retarget B onto A. The two
// pieces share no triangle edge, so the merged loop is mixed.
static bool mergeFillets(RegionSet& rs, int ia, int ib) {
    if (ia == ib) return false;
    if (ib < ia) std::swap(ia, ib);
    Region& A = rs.regions[(size_t)ia];
    Region& B = rs.regions[(size_t)ib];
    const int idA = A.id, idB = B.id;
    A.tris.insert(A.tris.end(), B.tris.begin(), B.tris.end());
    A.dVolPredicted += B.dVolPredicted;
    A.maxVertexDev = std::max(A.maxVertexDev, B.maxVertexDev);
    A.vMin = std::min(A.vMin, B.vMin);
    A.vMax = std::max(A.vMax, B.vMax);
    int la = -1, lb = -1;
    for (size_t i = 0; i < A.loops.size(); ++i)
        if (A.loops[i].role == LoopRole::Outer) la = (int)i;
    for (size_t i = 0; i < B.loops.size(); ++i)
        if (B.loops[i].role == LoopRole::Outer) lb = (int)i;
    if (la < 0 && !A.loops.empty()) la = 0;
    if (lb < 0 && !B.loops.empty()) lb = 0;
    if (la < 0 || lb < 0) return false;
    Loop& dst = A.loops[(size_t)la];
    const Loop& src = B.loops[(size_t)lb];
    dst.chainIdx.insert(dst.chainIdx.end(), src.chainIdx.begin(), src.chainIdx.end());
    dst.reversed.insert(dst.reversed.end(), src.reversed.begin(), src.reversed.end());
    for (size_t i = 0; i < B.loops.size(); ++i) {
        if ((int)i == lb) continue;
        A.loops.push_back(B.loops[i]);
    }
    for (auto& ch : rs.chains) {
        if (ch.regA == idB) ch.regA = idA;
        if (ch.regB == idB) ch.regB = idA;
    }
    for (int& tr : rs.triRegion)
        if (tr == idB) tr = idA;
    rs.regions.erase(rs.regions.begin() + ib);
    return true;
}

static int filletIndex(const RegionSet& rs, int which) {
    int seen = 0;
    for (size_t i = 0; i < rs.regions.size(); ++i) {
        const Region& r = rs.regions[i];
        if (r.origin != Origin::FilletStrip || r.type != SurfType::Cylinder || r.closed360)
            continue;
        if (seen == which) return (int)i;
        ++seen;
    }
    return -1;
}

struct Piece {
    std::vector<int> tris;
    double meanZ = 0;
};

// Edge-connected pieces of one fillet. The synthetic bar is two axial pieces
// that share the fitted cylinder and no mesh edge.
static std::vector<Piece> piecesOf(const MeshView& mv, const Region& r) {
    std::unordered_map<int, int> parent;
    for (int t : r.tris)
        if (t >= 0) parent[t] = t;
    auto find = [&](int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root) {
            const int nxt = parent[x];
            parent[x] = root;
            x = nxt;
        }
        return root;
    };
    auto unite = [&](int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[b] = a;
    };
    std::unordered_map<int, int> edgeOwner;
    for (int t : r.tris) {
        if (t < 0 || (size_t)t >= mv.nTri || !parent.count(t) || !mv.triEdges) continue;
        for (int s = 0; s < 3; ++s) {
            const int eid = mv.triEdges[t][s];
            if (eid < 0) continue;
            auto it = edgeOwner.find(eid);
            if (it == edgeOwner.end()) edgeOwner.emplace(eid, t);
            else unite(it->second, t);
        }
    }
    std::unordered_map<int, Piece> byRoot;
    for (int t : r.tris) {
        if (!parent.count(t)) continue;
        Piece& pc = byRoot[find(t)];
        pc.tris.push_back(t);
        if (t < 0 || (size_t)t >= mv.nTri || !mv.compTris || !mv.tris || !mv.pts) continue;
        const int gt = mv.compTris[t];
        if (gt < 0) continue;
        for (int k = 0; k < 3; ++k) pc.meanZ += mv.pts[mv.tris[gt][k]].Z();
    }
    std::vector<Piece> out;
    for (auto& kv : byRoot) {
        if (!kv.second.tris.empty()) kv.second.meanZ /= (double)(kv.second.tris.size() * 3);
        out.push_back(std::move(kv.second));
    }
    std::sort(out.begin(), out.end(),
              [](const Piece& a, const Piece& b) { return a.meanZ < b.meanZ; });
    return out;
}

// One loop whose chains cover every piece. That is mixed=1: the builder must
// split by component instead of exploding the strip.
static bool forceOneLoop(Region& r) {
    Loop one;
    one.role = LoopRole::Outer;
    for (const Loop& lp : r.loops) {
        one.chainIdx.insert(one.chainIdx.end(), lp.chainIdx.begin(), lp.chainIdx.end());
        one.reversed.insert(one.reversed.end(), lp.reversed.begin(), lp.reversed.end());
    }
    if (one.chainIdx.empty() || one.chainIdx.size() != one.reversed.size()) return false;
    r.loops.clear();
    r.loops.push_back(std::move(one));
    return true;
}

static bool buildCensus(const MeshView& mv, RegionSet& rs, FaceCensus& c, bool& ok) {
    ::setenv("STL2STEP_DIAG_REL", "1", 1);
    std::vector<TopoDS_Vertex> verts(mv.nVtx);
    BRep_Builder bb;
    for (size_t i = 0; i < mv.nVtx; ++i) {
        const gp_XYZ p = mv.pts[mv.compVtx[i]];
        bb.MakeVertex(verts[i], gp_Pnt(p), 1e-7);
    }
    std::vector<TopoDS_Face> faces;
    ok = buildFaces(mv, rs, verts, faces, nullptr);
    ::unsetenv("STL2STEP_DIAG_REL");
    c = censusOf(faces);
    std::fprintf(stderr, "build ok=%d faces=%zu cyl=%d plane=%d facet3=%d other=%d\n", (int)ok,
                 faces.size(), c.cyl, c.plane, c.facet3, c.other);
    dumpRegions(rs);
    return ok;
}

// Move vertices that belong only to `localTris` off the cylinder by 0.2 mm.
static int offsetExclusive(HarnessMesh& mesh, const MeshView& mv, const Region& fillet,
                           const std::vector<int>& localTris) {
    std::set<int> mine(localTris.begin(), localTris.end());
    std::set<int> inPiece, inOther;
    auto take = [&](int t, std::set<int>& dst) {
        if (t < 0 || (size_t)t >= mv.nTri || !mv.compTris || !mv.tris) return;
        const int gt = mv.compTris[t];
        if (gt < 0) return;
        for (int k = 0; k < 3; ++k) dst.insert(mv.tris[gt][k]);
    };
    for (int t : localTris) take(t, inPiece);
    for (size_t t = 0; t < mv.nTri; ++t)
        if (!mine.count((int)t)) take((int)t, inOther);
    int n = 0;
    for (int gv : inPiece) {
        if (inOther.count(gv) || gv < 0 || (size_t)gv >= mesh.pts.size()) continue;
        gp_Pnt p(mesh.pts[(size_t)gv]);
        gp_Vec d(fillet.ax.Location(), p);
        const gp_Vec ax(fillet.ax.Direction());
        gp_Vec radial = d - ax * d.Dot(ax);
        if (radial.Magnitude() < 1e-9) continue;
        radial.Normalize();
        const gp_XYZ np = p.XYZ() + radial.XYZ() * 0.2;
        mesh.pts[(size_t)gv] = np;
        ++n;
    }
    return n;
}

int main(int argc, char** argv) {
    const std::string exe = (argc > 0 && argv[0]) ? argv[0] : "fillet_strip_split_unit";
    const auto slash = exe.find_last_of('/');
    const std::string dir = slash == std::string::npos ? std::string(".") : exe.substr(0, slash);
    const std::string path = dir + "/fillet_strip_split_mesh.stl";
    check(writeMesh(path), "write_mesh");

    HarnessMesh mesh;
    MeshView mv{};
    RegionSet rs;
    std::string err;
    const bool seg = prepare(path, mesh, mv, rs, err);
    check(seg, "segment_closed");
    if (!seg) {
        std::fprintf(stderr, "segment err=%s\n", err.c_str());
        return 1;
    }
    int i0 = filletIndex(rs, 0);
    const int i1 = filletIndex(rs, 1);
    std::fprintf(stderr, "fillet idx %d %d q=%.6g tau=%.6g\n", i0, i1, mv.quantFloor,
                 mv.quantFloor > 0 ? 2 * mv.quantFloor : 0);
    check(i0 >= 0, "has_fillet_strip");
    if (i0 < 0) return 1;
    if (i1 >= 0) {
        check(std::fabs(rs.regions[(size_t)i0].radius - rs.regions[(size_t)i1].radius) < 1e-6,
              "same_radius");
        check(mergeFillets(rs, i0, i1), "merge_two_regions");
        i0 = filletIndex(rs, 0);
    }
    check(i0 >= 0 && filletIndex(rs, 1) < 0, "one_strip");
    if (i0 < 0) return 1;

    const std::vector<Piece> parts = piecesOf(mv, rs.regions[(size_t)i0]);
    std::fprintf(stderr, "pieces=%zu", parts.size());
    for (const Piece& pc : parts)
        std::fprintf(stderr, " nTri=%zu z=%.3f", pc.tris.size(), pc.meanZ);
    std::fprintf(stderr, " loops=%zu\n", rs.regions[(size_t)i0].loops.size());
    check(parts.size() == 2, "two_mesh_pieces");
    check(parts.size() == 2 && parts[0].tris.size() >= 3 && parts[1].tris.size() >= 3,
          "pieces_determined");
    if (parts.size() != 2) return 1;

    RegionSet natural = rs;
    FaceCensus natC;
    bool natOk = false;
    check(buildCensus(mv, natural, natC, natOk) && natOk, "natural_build");
    check(natC.cyl == 2, "natural_two_cyl");

    if (!prepare(path, mesh, mv, rs, err)) return 1;
    i0 = filletIndex(rs, 0);
    if (filletIndex(rs, 1) >= 0) {
        if (!mergeFillets(rs, i0, filletIndex(rs, 1))) return 1;
        i0 = filletIndex(rs, 0);
    }
    check(forceOneLoop(rs.regions[(size_t)i0]), "force_mixed_loop");
    FaceCensus bothC;
    bool bothOk = false;
    check(buildCensus(mv, rs, bothC, bothOk) && bothOk, "mixed_build");
    const BuiltAs bothAs = rs.regions[(size_t)i0].builtAs;
    std::fprintf(stderr, "mixed cyl=%d facet3=%d builtAs=%d\n", bothC.cyl, bothC.facet3,
                 (int)bothAs);
    check(bothC.cyl == 2, "mixed_two_faces");
    check(bothC.facet3 == natC.facet3, "mixed_same_facets");
    check(bothAs == BuiltAs::Single, "mixed_single");

    if (!prepare(path, mesh, mv, rs, err)) return 1;
    i0 = filletIndex(rs, 0);
    if (filletIndex(rs, 1) >= 0) {
        if (!mergeFillets(rs, i0, filletIndex(rs, 1))) return 1;
        i0 = filletIndex(rs, 0);
    }
    const std::vector<Piece> again = piecesOf(mv, rs.regions[(size_t)i0]);
    if (again.size() != 2) return 1;
    const std::vector<int> offTris = again[1].tris;
    const Region filletCopy = rs.regions[(size_t)i0];
    check(forceOneLoop(rs.regions[(size_t)i0]), "force_mixed_offset");
    const int nOff = offsetExclusive(mesh, mv, filletCopy, offTris);
    std::fprintf(stderr, "offset verts=%d pieceTris=%zu\n", nOff, offTris.size());
    check(nOff >= 3, "offset_interior");
    FaceCensus offC;
    bool offOk = false;
    check(buildCensus(mv, rs, offC, offOk) && offOk, "offset_build");
    int offAs = -1;
    for (const Region& r : rs.regions)
        if (r.origin == Origin::FilletStrip) offAs = (int)r.builtAs;
    std::fprintf(stderr, "offset cyl=%d facet3=%d builtAs=%d pieceTris=%zu\n", offC.cyl,
                 offC.facet3, offAs, offTris.size());
    check(offC.cyl == 1, "offset_one_face");
    // One analytic lid may leave the 3-edge bucket when its boundary becomes a
    // mesh polyline. The offset piece's triangles are 3-edge planes; nothing
    // else joins them.
    check(offC.facet3 >= (int)offTris.size() &&
              offC.facet3 <= (int)offTris.size() + natC.facet3,
          "offset_piece_facets");
    check(offAs == (int)BuiltAs::Single, "offset_single");

    // D-gid-7.4. The split cut the strip's tris down to the shipped piece; its
    // loops must follow, or a later rebuild (a site-B cascade exploding a
    // neighbour) builds one face on the old two-piece boundary.
    {
        const Region* fr = nullptr;
        for (const Region& r : rs.regions)
            if (r.origin == Origin::FilletStrip) fr = &r;
        bool follow = fr != nullptr && !fr->loops.empty();
        if (fr) {
            std::set<int> own;
            for (int t : fr->tris)
                for (int s = 0; s < 3; ++s) own.insert(mv.triEdges[t][s]);
            for (const Loop& lp : fr->loops)
                for (int ci : lp.chainIdx)
                    for (int eid : rs.chains[(size_t)ci].meshEdges)
                        if (!own.count(eid)) follow = false;
        }
        check(follow, "loops_follow_tris");
    }
    FaceCensus reC;
    bool reOk = false;
    check(buildCensus(mv, rs, reC, reOk) && reOk, "rebuild_build");
    int reAs = -1;
    for (const Region& r : rs.regions)
        if (r.origin == Origin::FilletStrip) reAs = (int)r.builtAs;
    std::fprintf(stderr, "rebuild cyl=%d facet3=%d builtAs=%d\n", reC.cyl, reC.facet3, reAs);
    check(reC.cyl == 1, "rebuild_one_face");
    check(reAs == (int)BuiltAs::Single, "rebuild_single");

    std::fprintf(stderr, "SUMMARY fail=%d\n", gFail);
    return gFail ? 1 : 0;
}
