// Independent STEP census. Reads a STEP with STEPControl_Reader and describes
// what is actually inside it. Does not include or call the stl2step engine.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_STEP_CENSUS_HPP
#define STL2STEP_STEP_CENSUS_HPP

#include <string>
#include <vector>

namespace step_census {

struct Vec3 {
    double x = 0, y = 0, z = 0;
};

struct Axis {
    Vec3 location;
    Vec3 direction;
};

struct CylinderFace {
    double radius = 0;
    Axis   axis;
    double uMin = 0, uMax = 0, vMin = 0, vMax = 0;
    double uSpan = 0;
    bool   periodicU = false;
    bool   closedU = false;
    bool   fullU = false;   // |uSpan| covers a full 2π — Seamed360 discriminator
};

// Cylindrical faces grouped by (radius, axis line). `builtAs` is the G2.5
// headline: "seamed360" (one 2π face) vs "twoHalves" (two ~π faces).
struct CylinderGroup {
    double      radius = 0;
    Axis        axis;
    int         nFaces = 0;
    std::vector<double> uSpans;
    std::string builtAs;    // seamed360 | twoHalves | partial | other
};

struct ShellInfo {
    bool closed = false;
};

struct Result {
    bool        ok = false;
    std::string error;
    std::string input;

    bool   valid = false;
    int    solids = 0;
    int    shells = 0;
    int    faces = 0;
    int    edges = 0;
    int    vertices = 0;
    int    openShells = 0;
    bool   closed = false;          // every shell is closed
    double volume = 0;

    std::vector<ShellInfo> shellList;

    int plane = 0, cylinder = 0, cone = 0, sphere = 0, torus = 0, bspline = 0, otherSurf = 0;

    std::vector<CylinderFace>  cylinders;
    std::vector<CylinderGroup> cylinderGroups;

    int line = 0, circle = 0, ellipse = 0, bsplineCurve = 0, otherCurve = 0;

    int pcurveSlots = 0;     // (edge, adjacent-face) pairs
    int pcurvePresent = 0;   // BRep_Tool::CurveOnSurface non-null
    int pcurveMissing = 0;   // J3 fails when this is > 0 on a closed body

    int contC0 = 0, contG1 = 0, contC1 = 0, contG2 = 0, contC2 = 0, contC3 = 0, contCN = 0;
    int contBoundary = 0;    // edges with < 2 faces (no two-face regularity)
    int contNonManifold = 0;

    double maxFaceTol = 0, maxEdgeTol = 0, maxVertexTol = 0;
};

// Read `path` independently of the stl2step engine. Never throws.
Result censusFile(const std::string& path);

// One JSON object, stable key order, no wall-clock. Byte-identical across
// two runs of the same file. A trailing newline is NOT included.
std::string toJson(const Result& r);

}  // namespace step_census

#endif
