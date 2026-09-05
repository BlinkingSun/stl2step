#ifndef GRADE_COMPARE_HPP
#define GRADE_COMPARE_HPP

#include "grade_mesh.hpp"
#include "grade_oracle.hpp"
#include "grade_step.hpp"

#include <string>
#include <vector>

namespace grade {

enum class Status { Recovered, Split, Sliver, Faceted, FacetedPartial, Spanned, Missing };

inline const char* statusName(Status s) {
    switch (s) {
        case Status::Recovered: return "recovered";
        case Status::Split: return "split";
        case Status::Sliver: return "sliver";
        case Status::Faceted: return "faceted";
        case Status::FacetedPartial: return "faceted-partial";
        case Status::Spanned: return "spanned";
        case Status::Missing: return "missing";
        default: return "missing";
    }
}

struct StepFaceRef {
    int entity = 0;
    double areaMM2 = 0;
};

struct Feature {
    Oracle oracle;
    Status status = Status::Missing;
    int splitN = 0;
    double credit = 0;
    double coverage = 0;
    double areaFraction = 0;
    std::vector<StepFaceRef> stepFaces;
};

struct Intersection {
    std::string a, b;
    int expectedTier = 2;
    std::string shipped;
    std::string verdict;
    int edgeEntity = 0;
};

struct UnmatchedFace {
    int entity = 0;
    std::string type;
    double areaMM2 = 0;
};

struct GradeConfig {
    bool reverseSeeds = false;
    int seedOrder = 0;  // 0 = file order; n != 0 = deterministic permutation n
    bool quiet = false;
    bool skipVolume = false;
    std::string engineBin;
};

struct GradeDocument {
    std::string stlPath, stepPath;
    std::string stlSha, stepSha;
    uint64_t stlBytes = 0, stepBytes = 0;
    Mesh mesh;
    OracleSet oracle;
    StepModel step;
    std::vector<Feature> features;
    std::vector<Intersection> intersections;
    std::vector<UnmatchedFace> unmatched;
    std::vector<std::string> hardZero;
    bool watertight = false;
    bool valid = false;
    std::string determinism;  // "skipped" | "pass" | "fail"
    double Vstep = 0, Vmesh = 0, volumeDelta = 0, chordBudget = 0, volumeQ = 0;
    bool volumeWithin = true;
    // grades: NaN means null
    double gradePlane = 0, gradeCyl = 0, gradeCone = 0, gradeSphere = 0, gradeTorus = 0,
           gradeOverall = 0;
    bool hasPlane = false, hasCyl = false, hasCone = false, hasSphere = false, hasTorus = false,
         hasOverall = false;
    double recPlane = 0, recCyl = 0, recCone = 0, recSphere = 0, recTorus = 0;
    double areaPlane = 0, areaCyl = 0, areaCone = 0, areaSphere = 0, areaTorus = 0;
    int nPlane = 0, nCyl = 0, nCone = 0, nSphere = 0, nTorus = 0;
    int cntPlane[7] = {}, cntCyl[7] = {}, cntCone[7] = {}, cntSphere[7] = {}, cntTorus[7] = {};
    double residueFacetedFraction = 0;
};

bool gradeFiles(const std::string& stl, const std::string& step, const GradeConfig& cfg,
                GradeDocument& doc, std::string& err);

}  // namespace grade

#endif
