// stl2step -- written-shape face census for the RESULT counters (D-gid-7.6).
// SPDX-License-Identifier: MIT
//
// smoothBuiltPlanes, smoothBuiltFillets, smoothPlanes, smoothCylinders and
// smoothFacetFaces used to be copied from the segment plan (RefitStats). The
// plan is not the file: a healed or exploded region stays in it, and the
// coplanar unify merges faces after it. These helpers count the shape that is
// written instead, from where each written face came from.
//
// Provenance. buildFaces names, for every face it returns, the region it was
// built for and whether it is one mesh triangle (refit::RegionSet::shipFaceRid,
// shipFaceFacet). The caller turns that into one FaceSrc bit and one "unit" per
// built face, seeds a map written-face -> built-face indices, and carries the
// map through each unify pass with the unify's own BRepTools_History. A merged
// face inherits every source of the faces merged into it.
//
// Definitions, over the faces of the written shape of every component whose
// analytic build shipped (smoothBuiltComponents; reverted and skipped
// components are not counted, as before):
//   smoothBuiltPlanes  = planar faces (GeomAbs_Plane) with at least one source
//                        that is not a mesh facet (an analytic region face or a
//                        stage-P profile face).
//   smoothFacetFaces   = faces whose every source is a mesh facet: one mesh
//                        triangle, or several coplanar ones the unify merged.
//   smoothBuiltFillets = faces with a FilletStrip-region source.
//   smoothPlanes       = distinct plane units with a face in the written shape
//                        (a unit is a Plane region of origin PlaneGrow, CylGrow
//                        or NgonWall -- the regions the plan counted -- or one
//                        stage-P planar face).
//   smoothCylinders    = the same for Cylinder regions (FilletStrip excluded,
//                        as before) and stage-P cylindrical faces.
// smoothBuiltCylinders / Cones / Tori stay surface-type counts of the written
// faces (D-130-10(1)); this census does not change them. smoothFillets is not
// in D-gid-7.6's list and stays the plan's FilletStrip region count.

#ifndef STL2STEP_WRITTEN_CENSUS_HPP
#define STL2STEP_WRITTEN_CENSUS_HPP

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools_History.hxx>
#include <NCollection_DataMap.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <cstdint>
#include <set>
#include <vector>

namespace stl2step {
namespace refit {
struct RegionSet;
}
namespace detail {

// One bit per kind of source; a merged face ORs the bits of its sources.
enum FaceSrc : uint8_t {
    kFaceSrcFacet = 1,     // one mesh triangle (island, exploded/E'-filled region, R2)
    kFaceSrcPlane = 2,     // Plane region (PlaneGrow/CylGrow/NgonWall) or stage-P plane
    kFaceSrcCylinder = 4,  // Cylinder region (PlaneGrow/CylGrow/NgonWall) or stage-P cylinder
    kFaceSrcFillet = 8,    // FilletStrip region
    kFaceSrcOther = 16,    // any other analytic face (cone, torus, sphere corner, ruled quad)
};

// written face (IsSame: TShape + Location, orientation ignored) -> indices of
// the built faces it came from.
using FaceSources = NCollection_DataMap<TopoDS_Shape, std::vector<int>, TopTools_ShapeMapHasher>;

inline void seedFaceSources(const std::vector<TopoDS_Face>& built, FaceSources& src) {
    src.Clear();
    for (int i = 0; i < (int)built.size(); ++i) {
        if (built[(size_t)i].IsNull()) continue;
        std::vector<int>* v = src.ChangeSeek(built[(size_t)i]);
        if (!v) v = src.Bound(built[(size_t)i], std::vector<int>());
        v->push_back(i);
    }
}

// Carry the map across one modification step. A null history is the identity
// (the step kept its input). A face the history modifies hands its sources to
// every image; a face it removes without an image drops them (the unify guard
// rebuilds such a removal, D-gid-2.1); any other face is unchanged.
inline void remapFaceSources(const Handle(BRepTools_History)& hist, FaceSources& src) {
    if (hist.IsNull() || src.IsEmpty()) return;
    FaceSources next;
    for (FaceSources::Iterator it(src); it.More(); it.Next()) {
        const TopoDS_Shape& f = it.Key();
        const std::vector<int>& from = it.Value();
        auto add = [&](const TopoDS_Shape& img) {
            std::vector<int>* v = next.ChangeSeek(img);
            if (!v) v = next.Bound(img, std::vector<int>());
            v->insert(v->end(), from.begin(), from.end());
        };
        const TopTools_ListOfShape& mod = hist->Modified(f);
        if (!mod.IsEmpty()) {
            for (TopTools_ListOfShape::Iterator m(mod); m.More(); m.Next())
                if (m.Value().ShapeType() == TopAbs_FACE) add(m.Value());
            continue;
        }
        if (hist->IsRemoved(f)) continue;
        add(f);
    }
    src = std::move(next);
}

struct WrittenFaceCounts {
    int analyticPlanes = 0;  // smoothBuiltPlanes
    int facetFaces = 0;      // smoothFacetFaces
    int filletFaces = 0;     // smoothBuiltFillets
    int planeUnits = 0;      // smoothPlanes
    int cylinderUnits = 0;   // smoothCylinders
    int faces = 0;           // diagnostics: written faces visited
    int unknownFaces = 0;    // diagnostics: written faces with no source
    int mixedFaces = 0;      // diagnostics: a facet merged with an analytic face
};

// kind[i] / unit[i] describe built face i. unit: the region id, a unique
// negative id below -1 for a face that is its own unit (stage P), -1 for none.
inline WrittenFaceCounts countWrittenFaces(const TopoDS_Shape& written, const FaceSources& src,
                                           const std::vector<uint8_t>& kind,
                                           const std::vector<int>& unit) {
    WrittenFaceCounts c;
    std::set<int> planeU, cylU;
    for (TopExp_Explorer ex(written, TopAbs_FACE); ex.More(); ex.Next()) {
        ++c.faces;
        const std::vector<int>* from = src.Seek(ex.Current());
        uint8_t bits = 0;
        if (from) {
            for (int i : *from) {
                if (i < 0 || (size_t)i >= kind.size()) continue;
                bits |= kind[(size_t)i];
                const int u = (size_t)i < unit.size() ? unit[(size_t)i] : -1;
                if (u == -1) continue;
                if (kind[(size_t)i] & kFaceSrcPlane) planeU.insert(u);
                if (kind[(size_t)i] & kFaceSrcCylinder) cylU.insert(u);
            }
        }
        if (bits == 0) {
            ++c.unknownFaces;
            continue;
        }
        if (bits == kFaceSrcFacet) {
            ++c.facetFaces;
        } else {
            if (bits & kFaceSrcFacet) ++c.mixedFaces;
            const GeomAbs_SurfaceType t =
                BRepAdaptor_Surface(TopoDS::Face(ex.Current()), Standard_False).GetType();
            if (t == GeomAbs_Plane) ++c.analyticPlanes;
        }
        if (bits & kFaceSrcFillet) ++c.filletFaces;
    }
    c.planeUnits = (int)planeU.size();
    c.cylinderUnits = (int)cylU.size();
    return c;
}

// kind[i] / unit[i] for each face buildFaces shipped (stl2step.cpp), from
// rs.shipFaceRid / rs.shipFaceFacet and the region's origin and type. A face
// with no recorded provenance gets kind 0 and is reported as unknownFaces.
void describeShippedFaces(const refit::RegionSet& rs, const std::vector<TopoDS_Face>& built,
                          std::vector<uint8_t>& kind, std::vector<int>& unit);

}  // namespace detail
}  // namespace stl2step

#endif  // STL2STEP_WRITTEN_CENSUS_HPP
