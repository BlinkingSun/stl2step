// Synthetic cases for the zero-area fin / out-and-back slit heal (D-fc-5.11 F3):
// the fin's bigon topology, the cut from the host's loops, the B-Rep coincidence
// of the two edge runs, and a cylinder face whose wire carries a slit (zero area,
// cut, valid whole face). No mesh.
// SPDX-License-Identifier: MIT

#include "refit_slit_heal.hpp"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <ElCLib.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <GProp_GProps.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

static BoundaryChain chain(int regA, int regB, std::vector<int> verts) {
    BoundaryChain ch;
    ch.regA = regA;
    ch.regB = regB;
    ch.meshVerts = std::move(verts);
    for (size_t i = 0; i + 1 < ch.meshVerts.size(); ++i) ch.meshEdges.push_back((int)i);
    return ch;
}

static Loop loop(LoopRole role, std::vector<int> cis, std::vector<uint8_t> rev = {}) {
    Loop lp;
    lp.role = role;
    lp.chainIdx = std::move(cis);
    lp.reversed = rev.empty() ? std::vector<uint8_t>(lp.chainIdx.size(), 0) : std::move(rev);
    return lp;
}

static Region region(int id, std::vector<Loop> loops) {
    Region r;
    r.id = id;
    r.type = SurfType::Plane;
    r.loops = std::move(loops);
    return r;
}

static TopoDS_Edge arcEdge(const gp_Circ& c, double u0, double u1) {
    Handle(Geom_TrimmedCurve) t = GC_MakeArcOfCircle(c, u0, u1, Standard_True).Value();
    return BRepBuilderAPI_MakeEdge(t).Edge();
}

static TopoDS_Edge reversedArcEdge(const gp_Circ& c, double u0, double u1) {
    // The same arc built as its own edge and traversed b -> a.
    return TopoDS::Edge(arcEdge(c, u0, u1).Reversed());
}

int main() {
    // A host H (id 0) whose outer loop runs 10 -> 11 -> [slit 11 -> 12 -> 11] -> 13 ->
    // 10, and whose inner loop is the bigon 20 <-> 21 of a second fin.
    RegionSet rs;
    rs.chains.push_back(chain(0, 7, {10, 11}));      // 0  host | wall
    rs.chains.push_back(chain(0, 1, {11, 12}));      // 1  host | fin 1 (direct)
    rs.chains.push_back(chain(1, 0, {11, 14, 12}));  // 2  fin 1 | host (around)
    rs.chains.push_back(chain(0, 7, {11, 13}));      // 3  host | wall
    rs.chains.push_back(chain(0, 8, {13, 10}));      // 4  host | other wall
    rs.chains.push_back(chain(0, 2, {20, 21}));      // 5  host | fin 2
    rs.chains.push_back(chain(2, 0, {21, 22, 20}));  // 6  fin 2 | host
    rs.chains.push_back(chain(3, 0, {30, 31}));      // 7  not-a-fin | host
    rs.chains.push_back(chain(3, 8, {31, 30}));      // 8  not-a-fin | other
    rs.chains.push_back(chain(4, 0, {40, 41}));      // 9  open pair, different ends
    rs.chains.push_back(chain(4, 0, {41, 42}));      // 10
    rs.chains.push_back(chain(5, -1, {50, 51}));     // 11 island neighbour
    rs.chains.push_back(chain(5, -1, {51, 50}));     // 12
    BoundaryChain ring = chain(6, 0, {60, 61, 62});  // 13 closed ring
    ring.closedLoop = true;
    rs.chains.push_back(ring);
    rs.chains.push_back(chain(6, 0, {62, 60}));      // 14

    const Region host = region(0, {loop(LoopRole::Outer, {0, 1, 2, 3, 4}, {0, 0, 1, 0, 0}),
                                   loop(LoopRole::Inner, {5, 6}, {0, 1})});
    const Region fin1 = region(1, {loop(LoopRole::Outer, {1, 2})});
    const Region fin2 = region(2, {loop(LoopRole::Outer, {5, 6})});

    // finBigonTopology: the bigon of one host, nothing else.
    {
        FinBigon fb;
        expect(finBigonTopology(fin1, rs, fb) && fb.fin == 1 && fb.host == 0 && fb.c1 == 1 &&
                   fb.c2 == 2,
               "topology: slit fin is a bigon of the host");
        FinBigon fb2;
        expect(finBigonTopology(fin2, rs, fb2) && fb2.host == 0 && fb2.c1 == 5 && fb2.c2 == 6,
               "topology: inner-loop fin is a bigon of the host");
        FinBigon x;
        expect(!finBigonTopology(region(3, {loop(LoopRole::Outer, {7, 8})}), rs, x),
               "topology: two chains to two regions is not a fin");
        expect(!finBigonTopology(region(4, {loop(LoopRole::Outer, {9, 10})}), rs, x),
               "topology: chains with different end pairs are not a bigon");
        expect(!finBigonTopology(region(5, {loop(LoopRole::Outer, {11, 12})}), rs, x),
               "topology: an island neighbour is never a host");
        expect(!finBigonTopology(region(6, {loop(LoopRole::Outer, {13, 14})}), rs, x),
               "topology: a closed chain is not a bigon side");
        expect(!finBigonTopology(region(1, {loop(LoopRole::Outer, {1, 2, 3})}), rs, x),
               "topology: three chains are not a bigon");
        expect(!finBigonTopology(region(1, {loop(LoopRole::Outer, {1, 2}),
                                            loop(LoopRole::Inner, {5, 6})}),
                                 rs, x),
               "topology: a region with two loops is not a fin");
        expect(!finBigonTopology(region(1, {loop(LoopRole::Inner, {1, 2})}), rs, x),
               "topology: the bigon must be the fin's outer loop");
        expect(!finBigonTopology(region(9, {loop(LoopRole::Outer, {1, 2})}), rs, x),
               "topology: chains that do not border the region");
    }

    // cutFinFromLoops: the slit leaves the outer loop, the bigon inner loop goes.
    {
        std::vector<Loop> L = host.loops;
        const bool ok = cutFinFromLoops(L, 1, 2);
        expect(ok && L.size() == 2 && L[0].chainIdx == std::vector<int>({0, 3, 4}) &&
                   L[0].reversed == std::vector<uint8_t>({0, 0, 0}) && L[0].role == LoopRole::Outer,
               "cut: consecutive out-and-back pair leaves the outer loop");
        const bool ok2 = cutFinFromLoops(L, 5, 6);
        expect(ok2 && L.size() == 1 && L[0].chainIdx == std::vector<int>({0, 3, 4}),
               "cut: an inner loop of exactly the pair is dropped");
        std::vector<Loop> W = {loop(LoopRole::Outer, {2, 0, 3, 4, 1}, {1, 0, 0, 0, 0})};
        expect(cutFinFromLoops(W, 1, 2) && W[0].chainIdx == std::vector<int>({0, 3, 4}),
               "cut: a pair across the loop's wrap is consecutive");
        std::vector<Loop> N = {loop(LoopRole::Outer, {1, 0, 2, 3, 4})};
        const std::vector<Loop> N0 = N;
        expect(!cutFinFromLoops(N, 1, 2) && N[0].chainIdx == N0[0].chainIdx,
               "cut: a non-consecutive pair is refused, loops unchanged");
        std::vector<Loop> O = {loop(LoopRole::Outer, {1, 2})};
        expect(!cutFinFromLoops(O, 1, 2) && O.size() == 1,
               "cut: an outer loop is never dropped whole");
        std::vector<Loop> S = {loop(LoopRole::Outer, {0, 1, 3, 4}), loop(LoopRole::Inner, {2, 5})};
        expect(!cutFinFromLoops(S, 1, 2) && S.size() == 2,
               "cut: a pair split across loops is refused");
        std::vector<Loop> M = {loop(LoopRole::Outer, {0, 3, 4})};
        expect(!cutFinFromLoops(M, 1, 2), "cut: a pair the host does not hold is refused");
        std::vector<Loop> D = {loop(LoopRole::Outer, {0, 1, 2, 1, 3})};
        expect(!cutFinFromLoops(D, 1, 2), "cut: a chain held twice is refused");
        std::vector<Loop> E = host.loops;
        expect(!cutFinFromLoops(E, 1, 1) && !cutFinFromLoops(E, -1, 2),
               "cut: a degenerate pair is refused");
    }

    // edgeRunsCoincide: the two runs of a bigon on one circle.
    {
        const double R = 3.0, d = 0.5 * kPi / 180.0;
        const gp_Circ c(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), R);
        const TopoDS_Edge out = arcEdge(c, 0.0, d);
        const TopoDS_Edge back = reversedArcEdge(c, 0.0, d);
        expect(edgeRunsCoincide({out}, {back}), "coincide: one arc out and back");
        const TopoDS_Edge half0 = arcEdge(c, 0.0, 0.5 * d), half1 = arcEdge(c, 0.5 * d, d);
        expect(edgeRunsCoincide({half0, half1}, {back}), "coincide: one arc split in two edges");
        const TopoDS_Edge other = arcEdge(c, d, 2.0 * kPi);
        expect(!edgeRunsCoincide({out}, {other}), "coincide: the complementary arc encloses the disc");
        const TopoDS_Edge chord =
            BRepBuilderAPI_MakeEdge(ElCLib::Value(0.0, c), ElCLib::Value(d, c)).Edge();
        expect(!edgeRunsCoincide({out}, {chord}),
               "coincide: a chord off its arc by more than the edge tolerance");
        const TopoDS_Edge l1 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(0, 1, 0)).Edge();
        const TopoDS_Edge l2 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 1, 0), gp_Pnt(0, 0, 0)).Edge();
        expect(edgeRunsCoincide({l1}, {l2}), "coincide: one segment out and back");
        const TopoDS_Edge k1 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(0.3, 0.5, 0)).Edge();
        const TopoDS_Edge k2 = BRepBuilderAPI_MakeEdge(gp_Pnt(0.3, 0.5, 0), gp_Pnt(0, 1, 0)).Edge();
        expect(!edgeRunsCoincide({l1}, {k1, k2}), "coincide: a dog-leg between the same ends encloses area");
        expect(!edgeRunsCoincide({}, {l2}) && !edgeRunsCoincide({l1}, {}),
               "coincide: an empty run is never zero area");
    }

    // A cylinder face whose outer wire makes an out-and-back slit along a circle into
    // the face. The slit carries no area; with the pair cut the face is valid and whole.
    // (Built on one exact pcurve the pair only overlaps and OCCT's analyzer passes it;
    // built by a component's chain collapse the same pair can read SelfIntersectingWire.)
    {
        const double R = 3.0, H = 10.0, span = 0.25 * kPi, d = 0.5 * kPi / 180.0, vs = 4.0;
        Handle(Geom_CylindricalSurface) cyl =
            new Geom_CylindricalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), R);
        auto uvEdge = [&](gp_Pnt2d a, gp_Pnt2d b) {
            const gp_Vec2d v(a, b);
            Handle(Geom2d_Line) ln = new Geom2d_Line(a, gp_Dir2d(v));
            return BRepBuilderAPI_MakeEdge(ln, cyl, 0.0, v.Magnitude()).Edge();
        };
        // u = 0 generator split at v = vs, slit along the circle v = vs to u = d.
        std::vector<TopoDS_Edge> E;
        E.push_back(uvEdge(gp_Pnt2d(0, 0), gp_Pnt2d(span, 0)));      // 0 low circle
        E.push_back(uvEdge(gp_Pnt2d(span, 0), gp_Pnt2d(span, H)));   // 1 far generator
        E.push_back(uvEdge(gp_Pnt2d(span, H), gp_Pnt2d(0, H)));      // 2 high circle
        E.push_back(uvEdge(gp_Pnt2d(0, H), gp_Pnt2d(0, vs)));        // 3 generator, upper
        E.push_back(uvEdge(gp_Pnt2d(0, vs), gp_Pnt2d(d, vs)));       // 4 slit out
        E.push_back(uvEdge(gp_Pnt2d(d, vs), gp_Pnt2d(0, vs)));       // 5 slit back
        E.push_back(uvEdge(gp_Pnt2d(0, vs), gp_Pnt2d(0, 0)));        // 6 generator, lower
        for (TopoDS_Edge& e : E) BRepLib::BuildCurves3d(e);
        // MakeWire shares the vertices where consecutive edges meet.
        BRepBuilderAPI_MakeWire mw;
        for (const TopoDS_Edge& e : E) mw.Add(e);
        expect(mw.IsDone(), "slit face: the wire with the slit closes");
        E.clear();
        for (TopExp_Explorer ex(mw.Wire(), TopAbs_EDGE); ex.More(); ex.Next())
            E.push_back(TopoDS::Edge(ex.Current()));
        auto faceOf = [&](const std::vector<int>& order) -> TopoDS_Face {
            BRep_Builder B;
            TopoDS_Wire w;
            B.MakeWire(w);
            for (int k : order) B.Add(w, E[(size_t)k]);
            w.Closed(Standard_True);
            BRepBuilderAPI_MakeFace mf(cyl, w, Standard_True);
            return mf.IsDone() ? mf.Face() : TopoDS_Face();
        };
        const TopoDS_Face slit = faceOf({0, 1, 2, 3, 4, 5, 6});
        expect(!slit.IsNull(), "slit face: the face with the slit is made");
        expect(edgeRunsCoincide({E[4]}, {E[5]}), "slit face: the slit's two runs coincide");
        // The host's loop as chains: 0..6 one edge each, the fin owns 4 and 5.
        std::vector<Loop> L = {loop(LoopRole::Outer, {0, 1, 2, 3, 4, 5, 6})};
        expect(cutFinFromLoops(L, 4, 5), "slit face: the pair is cut from the loop");
        const TopoDS_Face healed = faceOf(L[0].chainIdx);
        expect(!healed.IsNull() && BRepCheck_Analyzer(healed).IsValid(),
               "slit face: the healed face is valid");
        GProp_GProps ga, gb;
        BRepGProp::SurfaceProperties(slit, ga);
        BRepGProp::SurfaceProperties(healed, gb);
        const double full = R * span * H;
        expect(E.size() == 7 && std::fabs(gb.Mass() - full) <= 1e-9 * full &&
                   std::fabs(ga.Mass() - gb.Mass()) <= 1e-9 * full,
               "slit face: the slit carries no area; the healed face is the whole patch");
        int nSlit = 0;
        for (TopExp_Explorer ex(healed, TopAbs_EDGE); ex.More(); ex.Next())
            if (ex.Current().IsSame(E[4]) || ex.Current().IsSame(E[5])) ++nSlit;
        expect(nSlit == 0, "slit face: neither slit edge is left in the healed wire");
    }

    if (gFail) std::fprintf(stderr, "%d FAIL\n", gFail);
    else std::fprintf(stderr, "ALL PASS\n");
    return gFail ? 1 : 0;
}
