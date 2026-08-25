// STUB: tests/gates/regiondump/stub_segment.cpp
//
// Compiled ONLY when src/refit_segment.cpp is absent (CMake EXISTS guard on
// that path). This is not the P1 segmenter. It fills a small synthetic but
// schema-valid RegionSet so the dump serializer can be gated in wave 2.
// When src/refit_segment.cpp lands, this TU is not compiled — no file edit.
//
// SPDX-License-Identifier: MIT

#if defined(_MSC_VER)
#pragma message("stl2step_regiondump: STUB refit::segment() — src/refit_segment.cpp is absent; this TU is compiled only until the real one lands")
#else
#warning "stl2step_regiondump: STUB refit::segment() — src/refit_segment.cpp is absent; this TU is compiled only until the real one lands"
#endif

#include "refit.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <utility>

namespace stl2step {
namespace refit {
namespace {

static gp_Ax3 axZ() {
    return gp_Ax3(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0), gp_Dir(1.0, 0.0, 0.0));
}

static Loop makeLoop(int chainIdx, LoopRole role) {
    Loop lp;
    lp.chainIdx = {chainIdx};
    lp.reversed = {0};
    lp.role = role;
    return lp;
}

static BoundaryChain makeChain(int regA, int islandB, bool closed,
                               std::initializer_list<int> edges,
                               std::initializer_list<int> verts) {
    BoundaryChain ch;
    ch.regA = regA;
    ch.regB = -1;
    ch.islandA = -1;
    ch.islandB = islandB;
    ch.tangent = false;
    ch.closedLoop = closed;
    ch.meshEdges.assign(edges);
    ch.meshVerts.assign(verts);
    return ch;
}

static Region baseRegion(int id, SurfType type, Origin origin, bool closed360) {
    Region r;
    r.id = id;
    r.type = type;
    r.origin = origin;
    r.ax = axZ();
    r.radius = 0.0;
    r.uMin = 0.0;
    r.uMax = 0.0;
    r.vMin = 0.0;
    r.vMax = 0.0;
    r.closed360 = closed360;
    r.outwardNormal = true;
    r.reject = Reject::None;
    r.builtAs = BuiltAs::NotBuilt;
    r.filletNbrA = -1;
    r.filletNbrB = -1;
    return r;
}

}  // namespace

bool segment(const MeshView& mv, const SegmentParams& /*p*/, RegionSet& out, WarnFn /*warn*/) {
    out = RegionSet{};

    const int nTri = (int)mv.nTri;
    out.triRegion.assign((size_t)nTri, -1);
    out.triIsland.assign((size_t)nTri, -1);

    // Occupancy: first 3 tris -> plane region 0; next 5 -> cylinder region 1;
    // remainder -> island 0. I1-shaped (region XOR island) when nTri > 0.
    int cursor = 0;
    auto take = [&](int k) {
        const int a = cursor;
        cursor = std::min(nTri, cursor + k);
        return std::pair<int, int>{a, cursor};
    };
    const auto planeSpan = take(3);
    const auto cylSpan = take(5);
    const auto islandSpan = take(nTri);  // rest

    auto fillRange = [](std::vector<int>& dst, int lo, int hi, int value) {
        for (int t = lo; t < hi; t++) dst[(size_t)t] = value;
    };
    fillRange(out.triRegion, planeSpan.first, planeSpan.second, 0);
    fillRange(out.triRegion, cylSpan.first, cylSpan.second, 1);
    fillRange(out.triIsland, islandSpan.first, islandSpan.second, 0);

    auto trisOf = [](int lo, int hi) {
        std::vector<int> t;
        t.reserve((size_t)std::max(0, hi - lo));
        for (int i = lo; i < hi; i++) t.push_back(i);
        return t;
    };

    // closed360 == false: exactly one Outer + one Inner (I7).
    Region plane = baseRegion(0, SurfType::Plane, Origin::PlaneGrow, false);
    plane.tris = trisOf(planeSpan.first, planeSpan.second);
    plane.loops = {makeLoop(0, LoopRole::Outer), makeLoop(1, LoopRole::Inner)};
    plane.maxVertexDev = 0.0;
    plane.rmsVertexDev = 0.0;
    plane.dVolPredicted = 0.0;

    // closed360 == true: CapLow + CapHigh, no Outer (I7b).
    Region cyl = baseRegion(1, SurfType::Cylinder, Origin::CylGrow, true);
    cyl.radius = 5.0;
    cyl.uMin = 0.0;
    cyl.uMax = 2.0 * std::acos(-1.0);
    cyl.vMin = 0.0;
    cyl.vMax = 10.0;
    cyl.outwardNormal = false;
    cyl.tris = trisOf(cylSpan.first, cylSpan.second);
    cyl.loops = {makeLoop(2, LoopRole::CapLow), makeLoop(3, LoopRole::CapHigh)};
    cyl.chordSagitta = 0.190429;
    cyl.nSides = 8;
    cyl.dVolPredicted = -78.292;

    // One rejected region (diagnostics only; occupancy stays on the island).
    Region rej = baseRegion(2, SurfType::Plane, Origin::PlaneGrow, false);
    rej.tris = trisOf(islandSpan.first, islandSpan.second);
    rej.loops = {makeLoop(4, LoopRole::Outer)};
    rej.reject = Reject::VertexResidual;
    rej.maxVertexDev = 0.25;

    out.regions.push_back(std::move(plane));
    out.regions.push_back(std::move(cyl));
    out.rejected.push_back(std::move(rej));

    out.chains.push_back(makeChain(0, -1, true, {0, 1, 2, 3}, {0, 1, 2, 3}));
    out.chains.push_back(makeChain(0, 0, true, {4, 5, 6, 7}, {4, 5, 6, 7}));
    out.chains.push_back(makeChain(1, -1, true, {8, 9, 10, 11}, {8, 9, 10, 11}));
    out.chains.push_back(makeChain(1, -1, true, {12, 13, 14, 15}, {12, 13, 14, 15}));
    out.chains.push_back(makeChain(-1, 0, true, {16, 17}, {16, 17}));

    out.nIslands = 1;
    out.stats.planes = 1;
    out.stats.cylinders = 1;
    out.stats.fillets = 0;
    out.stats.rejected = 1;
    out.stats.facetIslands = 1;
    out.stats.facetTriangles = std::max(0, islandSpan.second - islandSpan.first);
    out.stats.distinctRadii = 1;
    out.stats.maxVertexDev = 0.25;
    out.stats.maxEdgeTol = 0.01;
    out.stats.dVolPredicted = -78.292;
    return true;
}

}  // namespace refit
}  // namespace stl2step
