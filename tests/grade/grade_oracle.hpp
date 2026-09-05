#ifndef GRADE_ORACLE_HPP
#define GRADE_ORACLE_HPP

#include "grade_mesh.hpp"

#include <string>
#include <vector>

namespace grade {

struct SurfParams {
    SurfClass cls = SurfClass::Plane;
    Vec3 n{};      // plane normal / axis dir
    Vec3 p0{};     // plane origin (foot) / axis loc / sphere centre / torus centre
    Vec3 apex{};   // cone
    double R = 0;  // cylinder/sphere radius, torus major
    double r = 0;  // torus minor
    double alpha = 0;
};

struct Oracle {
    int id = 0;
    SurfClass cls = SurfClass::Plane;
    SurfParams S{};
    std::vector<int> tris;
    std::vector<int> verts;  // distinct welded vertex indices
    double w = 0;
    double maxResid = 0;
    double maxNormalDev = 0;
    int oraclePieces = 1;
    int oracleMinSeparation = 0;
    int minVertIndex = 0;
    Vec3 bboxMin{}, bboxMax{};
    std::string featureId;
};

struct ResidueComp {
    std::vector<int> tris;
    double area = 0;
    bool ruled = false;
    Vec3 direction{};
    Vec3 bboxMin{}, bboxMax{};
};

struct OracleSet {
    std::vector<Oracle> oracles;
    std::vector<int> owner;  // per triangle, -1 residue
    std::vector<ResidueComp> residue;
    int unprovableSingletons = 0;
    double residueArea = 0;
    int residueTris = 0;
};

bool fitClass(const Mesh& m, const std::vector<int>& region, SurfClass c, SurfParams& S);
bool certifies(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S,
               double* maxResidOut = nullptr, double* maxNDevOut = nullptr);
bool admits(const Mesh& m, const std::vector<int>& region, SurfClass c, const SurfParams& S);
double distToSurf(const Vec3& v, const SurfParams& S);
Vec3 normalAt(const SurfParams& S, const Vec3& p);

void buildOracle(const Mesh& m, OracleSet& out, bool reverseSeeds);

void assignFeatureIds(const Mesh& m, OracleSet& set);

}  // namespace grade

#endif
