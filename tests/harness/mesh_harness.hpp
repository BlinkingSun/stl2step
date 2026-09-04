// Test-side reimplementation of engine stages 2–3 (weld + manifold split).
// Owns the storage that a later refit::MeshView will view; does not depend on
// src/refit.hpp. Loops are ported in source order from src/stl2step.cpp.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_MESH_HARNESS_HPP
#define STL2STEP_MESH_HARNESS_HPP

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <gp_XYZ.hxx>

namespace stl2step {
namespace harness {

// One manifold component; owns its storage. Local vertex/edge/triangle ids
// match the engine's first-touch assignment (invariant I5).
struct HarnessComponent {
    int  root = -1;                // engine's UnionFind root id (component key)
    std::vector<int>                    compTris;   // local tri  -> global tri id
    std::vector<int>                    compVtx;    // local vtx  -> global point id
    std::vector<std::pair<int,int>>     compEdges;  // local edge -> (local vLo, vHi)
    std::vector<std::array<int,3>>      triEdges;   // local tri  -> 3 local edge ids
    std::vector<uint8_t>                triDirs;    // local tri  -> direction bits
    int    open = 0, conflict = 0, nonManifold = 0; // engine's health counters
    double vol = 0;                                 // engine's signed component volume
    bool   clean() const { return open == 0 && conflict == 0 && nonManifold == 0; }
};

// One whole STL; owns the global tables. comps is ordered by ascending root id
// (a gate-tool contract; the engine's size-then-id order is build scheduling).
struct HarnessMesh {
    std::vector<gp_XYZ>            pts;    // welded, scaled
    std::vector<std::array<int,3>> tris;   // global triangles, post-remap, degenerates dropped
    double diag = 0, weldTol = 0, sewTol = 0;
    double quantFloor = 0;                 // D-130-12, mm (post-scale); see src/stl_quant.hpp
    std::vector<HarnessComponent> comps;
};

// Reads an STL (RWStl, binary + ASCII autodetected), welds, splits, and fills
// `out`. Returns false with `err` set on an unreadable/empty file. Never throws.
bool loadMesh(const std::string& stlPath, double scale, double weldTolArg,
              double sewTolArg, HarnessMesh& out, std::string& err);

}  // namespace harness
}  // namespace stl2step

#endif
