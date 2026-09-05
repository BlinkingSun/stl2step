#ifndef GRADE_MESH_HPP
#define GRADE_MESH_HPP

#include "stl_quant.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace grade {

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(double s, const Vec3& a) { return a * s; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm2(const Vec3& a) { return dot(a, a); }
inline double norm(const Vec3& a) { return std::sqrt(norm2(a)); }
inline Vec3 normalized(const Vec3& a) {
    const double n = norm(a);
    return n > 0.0 ? a * (1.0 / n) : Vec3{};
}
inline double dist(const Vec3& a, const Vec3& b) { return norm(a - b); }

inline double clamp1(double x) {
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}
inline double angleUnit(const Vec3& a, const Vec3& b) { return std::acos(clamp1(dot(a, b))); }

inline void canonDir(Vec3& d) {
    d = normalized(d);
    if (d.x < 0.0 || (d.x == 0.0 && d.y < 0.0) || (d.x == 0.0 && d.y == 0.0 && d.z < 0.0))
        d = -d;
}

inline Vec3 footFromOrigin(const Vec3& loc, const Vec3& dir) {
    const Vec3 a = normalized(dir);
    return loc - a * dot(loc, a);
}

enum class SurfClass { Plane = 0, Cylinder = 1, Cone = 2, Sphere = 3, Torus = 4, Other = 5 };

inline const char* className(SurfClass c) {
    switch (c) {
        case SurfClass::Plane: return "plane";
        case SurfClass::Cylinder: return "cylinder";
        case SurfClass::Cone: return "cone";
        case SurfClass::Sphere: return "sphere";
        case SurfClass::Torus: return "torus";
        default: return "other";
    }
}

inline int paramCount(SurfClass c) {
    switch (c) {
        case SurfClass::Plane: return 3;
        case SurfClass::Sphere: return 4;
        case SurfClass::Cylinder: return 5;
        case SurfClass::Cone: return 6;
        case SurfClass::Torus: return 7;
        default: return 0;
    }
}

struct Tri {
    int v[3] = {0, 0, 0};
    double area = 0, hmin = 0, perim = 0, circumdiam = 0, thetaQ = 0;
    Vec3 n{}, centroid{};
};

struct Mesh {
    std::vector<Vec3> verts;
    std::vector<Tri> tris;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<int>> vertTris;
    int openEdges = 0;
    int nonManifoldEdges = 0;
    int degenerateDropped = 0;
    double surfaceArea = 0;
    double volume = 0;
    double meanCircumdiam = 0;
    double q = 0;
    double tau = 0;
    stl2step::StlQuantFloor qf{};
    std::string path;

    // Exact-bit lookup for binary (float32 bits) and a q-ball spatial hash.
    struct Bit3 {
        uint32_t a = 0, b = 0, c = 0;
        bool operator==(const Bit3& o) const { return a == o.a && b == o.b && c == o.c; }
    };
    struct Bit3Hash {
        size_t operator()(const Bit3& k) const {
            size_t h = k.a;
            h ^= static_cast<size_t>(k.b) + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= static_cast<size_t>(k.c) + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<Bit3, int, Bit3Hash> bitIndex;
    std::unordered_map<int64_t, std::vector<int>> grid;
    double gridCell = 1.0;

    int nearestVert(const Vec3& p, double maxDist) const;
};

bool loadStl(const std::string& path, Mesh& mesh, std::string& err);

inline double areaQ(const Mesh& m, const std::vector<int>& region) {
    double s = 0;
    for (int t : region) s += m.tris[static_cast<size_t>(t)].perim;
    return (m.q * 0.5) * s;
}

inline double regionArea(const Mesh& m, const std::vector<int>& region) {
    double s = 0;
    for (int t : region) s += m.tris[static_cast<size_t>(t)].area;
    return s;
}

void eigen3(const double S[3][3], double eval[3], double evec[3][3]);
bool solveN(int n, double* A, double* b, double* x);  // row-major n*n

inline double roundToQ(double x, double q) {
    if (!(q > 0.0)) return x;
    return std::round(x / q) * q;
}

inline int printDigits(double q) {
    if (!(q > 0.0)) return 9;
    const int k = static_cast<int>(std::ceil(-std::log10(q))) + 1;
    return k < 0 ? 0 : k;
}

std::string relPath(const std::string& stlPath, const std::string& other);
std::string fmtG(double v);

}  // namespace grade

#endif
