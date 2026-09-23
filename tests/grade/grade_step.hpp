#ifndef GRADE_STEP_HPP
#define GRADE_STEP_HPP

#include "grade_mesh.hpp"
#include "grade_oracle.hpp"

#include <mutex>
#include <string>
#include <vector>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>

namespace grade {

struct StepFace {
    int entity = 0;  // 1-based index in the face map
    SurfClass cls = SurfClass::Other;
    SurfParams S{};
    double area = 0;
    std::vector<Vec3> wireVerts;
    std::vector<int> meshVerts;  // matched mesh vertex indices, -1 if unmatched
    bool facet = false;
    TopoDS_Face face;
};

struct StepModel {
    TopoDS_Shape shape;
    std::vector<StepFace> faces;
    int nFacet = 0, nPlane = 0, nCyl = 0, nCone = 0, nSphere = 0, nTorus = 0, nBSpline = 0,
        nOther = 0;
    bool valid = false;
    bool watertight = false;
    double volume = 0;
    std::vector<int> cover;  // per mesh triangle: facet face entity or -1
};

bool loadStep(const std::string& path, const Mesh& mesh, StepModel& out, std::string& err,
              bool computeVolume);

// Serializes OCCT that is not thread-safe or must not overlap a STEP
// reader: loadStep, BRep_Tool / BRepTopAdaptor_FClass2d /
// GeomAPI_ProjectPointOnSurf / BRepAdaptor_Curve / TopExp::MapShapes on
// STEP faces. gp_* and IntAna_QuadQuadGeo stay outside (stack math).
// STEPControl_*, Interface_Static, Message::DefaultMessenger, and
// OSD::SetSignal are process-global; a reader overlapping BRep/Geom on
// another thread deadlocks on OCCT 8 / MSVC. One mutex, one lock order.
std::mutex& occtMutex();

// Drops the cout printer and calls OSD::SetSignal once on the first
// thread, then OSD::SetThreadLocalSignal on this thread (does not replace
// the process filter). Call from the main thread before a worker pool,
// and from each worker before its first OCCT call.
void prepareOcctThread();

}  // namespace grade

#endif
