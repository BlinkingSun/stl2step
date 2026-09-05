#ifndef GRADE_STEP_HPP
#define GRADE_STEP_HPP

#include "grade_mesh.hpp"
#include "grade_oracle.hpp"

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

}  // namespace grade

#endif
