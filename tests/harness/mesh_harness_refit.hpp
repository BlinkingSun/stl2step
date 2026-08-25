// Non-owning adapter from HarnessMesh/HarnessComponent to refit::MeshView.
// Compiles to nothing until src/refit.hpp exists (p-header lane).
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_MESH_HARNESS_REFIT_HPP
#define STL2STEP_MESH_HARNESS_REFIT_HPP

#include "mesh_harness.hpp"

#if __has_include("../../src/refit.hpp")
#include "../../src/refit.hpp"
#define STL2STEP_HARNESS_HAVE_REFIT 1

namespace stl2step {
namespace harness {

// Fills a non-owning refit::MeshView over `m`/`c`. The MeshView is valid only as
// long as both outlive it; nothing is copied.
inline void toMeshView(const HarnessMesh& m, const HarnessComponent& c,
                       stl2step::refit::MeshView& mv) {
    static_assert(sizeof(std::array<int, 3>) == sizeof(int[3]),
                  "std::array<int,3> must be layout-compatible with int[3]");
    mv.pts       = m.pts.data();
    mv.tris      = reinterpret_cast<const int(*)[3]>(m.tris.data());
    mv.compTris  = c.compTris.data();
    mv.compVtx   = c.compVtx.data();
    mv.compEdges = c.compEdges.data();
    mv.triEdges  = c.triEdges.data();
    mv.triDirs   = c.triDirs.data();
    mv.nTri      = c.compTris.size();
    mv.nVtx      = c.compVtx.size();
    mv.nEdge     = c.compEdges.size();
    mv.diag      = m.diag;
    mv.weldTol   = m.weldTol;
    mv.sewTol    = m.sewTol;
}

}  // namespace harness
}  // namespace stl2step

#endif  // __has_include("../../src/refit.hpp")

#endif
