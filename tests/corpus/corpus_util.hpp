// Shared mesh I/O, measurement, and JSON helpers for the P0 corpus lane.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_CORPUS_UTIL_HPP
#define STL2STEP_CORPUS_UTIL_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace corpus {

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    explicit Vec3(const gp_Pnt& p) : x(p.X()), y(p.Y()), z(p.Z()) {}
    explicit Vec3(const gp_Dir& d) : x(d.X()), y(d.Y()), z(d.Z()) {}
};

struct AxisSpec {
    Vec3 loc;
    Vec3 dir;
};

struct Recoverable {
    std::string type;   // "plane" | "cylinder"
    double radius = 0;
    AxisSpec axis;
    int count = 1;
    int nSides = 0;
    bool closed360 = false;
    Vec3 normal;        // planes only
    // Optional axial window (axis parameter) so fillet nSides is measured on
    // the strip, not the adjacent spherical blends. NaN => no window.
    double vMin = std::numeric_limits<double>::quiet_NaN();
    double vMax = std::numeric_limits<double>::quiet_NaN();
};

struct FacetedRegion {
    std::string type;   // "sphere" | "torus" | "freeform" | ...
    int count = 1;
    std::string note;
};

struct Sidecar {
    std::string id;
    std::string description;
    double deflection = 0;
    int triangleCount = 0;
    double meshVolume = 0;
    std::vector<Recoverable> recoverable;
    std::vector<FacetedRegion> mustRemainFaceted;
    std::vector<std::string> expectedRejects;
    bool expectFillet = false;
    std::string rLadderRole;  // S16 sub-fixtures only
    // Authored-mesh manifold (after coincident-vertex weld). S15 is dirty on
    // purpose; S02/S09 inherit OCCT fillet/loft tessellation tears.
    int openEdges = 0;
    int nonManifoldEdges = 0;
    bool watertight = true;
    // Engine convert (smooth off) contract — not the same as mesh watertight.
    int expectedExit = 0;
    int expectedSolids = 1;
    int expectedOpenShells = 0;
};

struct MeshData {
    std::vector<Vec3> verts;
    std::vector<std::array<int, 3>> tris;
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

inline Vec3 normalize(const Vec3& v) {
    const double n = norm(v);
    if (n < 1e-30) return {0, 0, 1};
    return {v.x / n, v.y / n, v.z / n};
}

inline Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3 add(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3 scale(const Vec3& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

inline double triArea(const Vec3& a, const Vec3& b, const Vec3& c) {
    return 0.5 * norm(cross(sub(b, a), sub(c, a)));
}

inline Vec3 triNormal(const Vec3& a, const Vec3& b, const Vec3& c) {
    return normalize(cross(sub(b, a), sub(c, a)));
}

inline Vec3 triCentroid(const Vec3& a, const Vec3& b, const Vec3& c) {
    return {(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0};
}

inline float roundFloat(double v) {
    return static_cast<float>(std::rint(v * 1e6) / 1e6);
}

inline MeshData extractMesh(const TopoDS_Shape& shape) {
    MeshData out;
    TopExp_Explorer exp(shape, TopAbs_FACE);
    for (; exp.More(); exp.Next()) {
        const TopoDS_Face face = TopoDS::Face(exp.Current());
        TopLoc_Location loc;
        Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
        if (tri.IsNull()) continue;
        const gp_Trsf tr = loc.Transformation();
        const int base = static_cast<int>(out.verts.size());
        for (int i = 1; i <= tri->NbNodes(); ++i) {
            gp_Pnt p = tri->Node(i);
            p.Transform(tr);
            out.verts.emplace_back(p);
        }
        const bool reversed = (face.Orientation() == TopAbs_REVERSED);
        for (int i = 1; i <= tri->NbTriangles(); ++i) {
            int n1, n2, n3;
            tri->Triangle(i).Get(n1, n2, n3);
            if (reversed) std::swap(n1, n2);
            out.tris.push_back({base + n1 - 1, base + n2 - 1, base + n3 - 1});
        }
    }
    return out;
}

inline void tessellate(TopoDS_Shape& shape, double linearDeflection, double angularDeflection = 0.5) {
    IMeshTools_Parameters params;
    params.Deflection = linearDeflection;
    params.Angle = angularDeflection;
    params.Relative = Standard_False;
    params.InParallel = Standard_False;
    params.AllowQualityDecrease = Standard_False;
    BRepMesh_IncrementalMesh mesher;
    mesher.SetShape(shape);
    mesher.ChangeParameters() = params;
    mesher.Perform();
}

inline double brepVolume(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();
}

inline double meshVolume(const MeshData& mesh) {
    double vol = 0;
    for (const auto& t : mesh.tris) {
        const Vec3& a = mesh.verts[t[0]];
        const Vec3& b = mesh.verts[t[1]];
        const Vec3& c = mesh.verts[t[2]];
        vol += dot(a, cross(b, c)) / 6.0;
    }
    return vol;
}

inline MeshData quantizeMesh(const MeshData& mesh) {
    MeshData out = mesh;
    for (auto& v : out.verts) {
        v.x = roundFloat(v.x);
        v.y = roundFloat(v.y);
        v.z = roundFloat(v.z);
    }
    return out;
}

// Weld vertices that share a quantized key and drop index-degenerate triangles.
inline MeshData weldQuantized(const MeshData& mesh) {
    MeshData out;
    std::map<std::array<float, 3>, int> index;
    auto getIdx = [&](const Vec3& p) {
        const std::array<float, 3> key = {roundFloat(p.x), roundFloat(p.y), roundFloat(p.z)};
        auto it = index.find(key);
        if (it != index.end()) return it->second;
        const int id = static_cast<int>(out.verts.size());
        out.verts.push_back({key[0], key[1], key[2]});
        index.emplace(key, id);
        return id;
    };
    for (const auto& t : mesh.tris) {
        const int a = getIdx(mesh.verts[t[0]]);
        const int b = getIdx(mesh.verts[t[1]]);
        const int c = getIdx(mesh.verts[t[2]]);
        if (a == b || b == c || a == c) continue;
        out.tris.push_back({a, b, c});
    }
    return out;
}

struct EdgeStats {
    int edges = 0;
    int open = 0;
    int nonManifold = 0;
};

inline EdgeStats meshEdgeStats(const MeshData& mesh) {
    auto ek = [](int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
    std::map<std::pair<int, int>, int> use;
    for (const auto& t : mesh.tris) {
        use[ek(t[0], t[1])]++;
        use[ek(t[1], t[2])]++;
        use[ek(t[2], t[0])]++;
    }
    EdgeStats s;
    s.edges = static_cast<int>(use.size());
    for (const auto& kv : use) {
        if (kv.second == 1) s.open++;
        else if (kv.second > 2) s.nonManifold++;
    }
    return s;
}

inline bool writeBinaryStl(const std::string& path, const MeshData& mesh, const char* label) {
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;
    char header[80] = {};
    std::snprintf(header, sizeof(header), "%.79s", label);
    out.write(header, 80);
    const uint32_t nTri = static_cast<uint32_t>(mesh.tris.size());
    out.write(reinterpret_cast<const char*>(&nTri), 4);
    for (const auto& t : mesh.tris) {
        const Vec3& a = mesh.verts[t[0]];
        const Vec3& b = mesh.verts[t[1]];
        const Vec3& c = mesh.verts[t[2]];
        const Vec3 n = triNormal(a, b, c);
        const float rec[12] = {
            roundFloat(n.x), roundFloat(n.y), roundFloat(n.z),
            roundFloat(a.x), roundFloat(a.y), roundFloat(a.z),
            roundFloat(b.x), roundFloat(b.y), roundFloat(b.z),
            roundFloat(c.x), roundFloat(c.y), roundFloat(c.z)};
        out.write(reinterpret_cast<const char*>(rec), 48);
        const uint16_t attr = 0;
        out.write(reinterpret_cast<const char*>(&attr), 2);
    }
    return static_cast<bool>(out);
}

inline std::string jsonNum(double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

inline std::string jsonVec3(const Vec3& v) {
  return "[" + jsonNum(v.x) + "," + jsonNum(v.y) + "," + jsonNum(v.z) + "]";
}

inline std::string writeSidecarJson(const Sidecar& s) {
    std::ostringstream os;
    os << std::setprecision(17);
    os << "{\n";
    os << "  \"id\": \"" << s.id << "\",\n";
    os << "  \"description\": \"" << s.description << "\",\n";
    os << "  \"deflection\": " << s.deflection << ",\n";
    os << "  \"triangleCount\": " << s.triangleCount << ",\n";
    os << "  \"meshVolume\": " << s.meshVolume << ",\n";
    if (!s.rLadderRole.empty())
        os << "  \"rLadderRole\": \"" << s.rLadderRole << "\",\n";
    os << "  \"expectFillet\": " << (s.expectFillet ? "true" : "false") << ",\n";
    os << "  \"watertight\": " << (s.watertight ? "true" : "false") << ",\n";
    os << "  \"openEdges\": " << s.openEdges << ",\n";
    os << "  \"nonManifoldEdges\": " << s.nonManifoldEdges << ",\n";
    os << "  \"expectedExit\": " << s.expectedExit << ",\n";
    os << "  \"expectedSolids\": " << s.expectedSolids << ",\n";
    os << "  \"expectedOpenShells\": " << s.expectedOpenShells << ",\n";
    os << "  \"recoverable\": [\n";
    for (size_t i = 0; i < s.recoverable.size(); ++i) {
        const auto& r = s.recoverable[i];
        os << "    {\n";
        os << "      \"type\": \"" << r.type << "\",\n";
        if (r.type == "cylinder") {
            os << "      \"radius\": " << r.radius << ",\n";
            os << "      \"axis\": {\"loc\": " << jsonVec3(r.axis.loc)
               << ", \"dir\": " << jsonVec3(r.axis.dir) << "},\n";
            os << "      \"count\": " << r.count << ",\n";
            os << "      \"nSides\": " << r.nSides << ",\n";
            os << "      \"closed360\": " << (r.closed360 ? "true" : "false");
            if (std::isfinite(r.vMin) && std::isfinite(r.vMax)) {
                os << ",\n      \"vMin\": " << r.vMin << ",\n      \"vMax\": " << r.vMax << "\n";
            } else {
                os << "\n";
            }
        } else {
            os << "      \"count\": " << r.count << ",\n";
            os << "      \"normal\": " << jsonVec3(r.normal) << ",\n";
            os << "      \"nSides\": 0\n";
        }
        os << "    }" << (i + 1 < s.recoverable.size() ? "," : "") << "\n";
    }
    os << "  ],\n";
    os << "  \"mustRemainFaceted\": [\n";
    for (size_t i = 0; i < s.mustRemainFaceted.size(); ++i) {
        const auto& f = s.mustRemainFaceted[i];
        os << "    {\"type\": \"" << f.type << "\", \"count\": " << f.count;
        if (!f.note.empty()) os << ", \"note\": \"" << f.note << "\"";
        os << "}" << (i + 1 < s.mustRemainFaceted.size() ? "," : "") << "\n";
    }
    os << "  ],\n";
    os << "  \"expectedRejects\": [";
    for (size_t i = 0; i < s.expectedRejects.size(); ++i) {
        if (i) os << ", ";
        os << "\"" << s.expectedRejects[i] << "\"";
    }
    os << "]\n";
    os << "}\n";
    return os.str();
}

inline bool writeTextFile(const std::string& path, const std::string& text) {
    std::ofstream out(path);
    if (!out) return false;
    out << text;
    return static_cast<bool>(out);
}

// D2 nSides from facet-normal bands + vertex span (DECISION-p1-math §2.2).
inline int measureCylinderNSides(const MeshData& mesh, const AxisSpec& axis, double radius,
                                 double tolFrac = 0.08,
                                 double vMin = std::numeric_limits<double>::quiet_NaN(),
                                 double vMax = std::numeric_limits<double>::quiet_NaN()) {
    const Vec3 a = normalize(axis.dir);
  const Vec3 w = (std::fabs(a.x) <= std::fabs(a.y) && std::fabs(a.x) <= std::fabs(a.z))
                     ? Vec3(1, 0, 0)
                 : (std::fabs(a.y) <= std::fabs(a.z) ? Vec3(0, 1, 0) : Vec3(0, 0, 1));
    const Vec3 u0 = normalize(cross(a, w));
    const Vec3 v0 = cross(a, u0);

    auto radialDist = [&](const Vec3& p) {
        const Vec3 d = sub(p, axis.loc);
        const Vec3 perp = sub(d, scale(a, dot(d, a)));
        return norm(perp);
    };

    const double bandTol = std::max(0.05 * radius, radius * tolFrac);
    const bool windowed = std::isfinite(vMin) && std::isfinite(vMax) && vMax > vMin;
    std::vector<int> cylTris;
    std::vector<Vec3> cylVerts;
    std::vector<Vec3> cylNormals;
    for (size_t ti = 0; ti < mesh.tris.size(); ++ti) {
        const auto& t = mesh.tris[ti];
        const Vec3& p0 = mesh.verts[t[0]];
        const Vec3& p1 = mesh.verts[t[1]];
        const Vec3& p2 = mesh.verts[t[2]];
        const Vec3 n = triNormal(p0, p1, p2);
        if (std::fabs(dot(n, a)) > 0.35) continue;
        const Vec3 cent = triCentroid(p0, p1, p2);
        if (windowed) {
            const double v = dot(sub(cent, axis.loc), a);
            if (v < vMin || v > vMax) continue;
        }
        const double r0 = radialDist(p0);
        const double r1 = radialDist(p1);
        const double r2 = radialDist(p2);
        const double rMean = (r0 + r1 + r2) / 3.0;
        if (std::fabs(rMean - radius) > bandTol) continue;
        cylTris.push_back(static_cast<int>(ti));
        cylNormals.push_back(n);
        cylVerts.push_back(p0);
        cylVerts.push_back(p1);
        cylVerts.push_back(p2);
    }
    if (cylTris.empty()) return 0;

    Vec3 centroid{0, 0, 0};
    double areaSum = 0;
    for (int ti : cylTris) {
        const auto& t = mesh.tris[ti];
        const double ar = triArea(mesh.verts[t[0]], mesh.verts[t[1]], mesh.verts[t[2]]);
        const Vec3 c = triCentroid(mesh.verts[t[0]], mesh.verts[t[1]], mesh.verts[t[2]]);
        centroid = add(centroid, scale(c, ar));
        areaSum += ar;
    }
    if (areaSum > 0) centroid = scale(centroid, 1.0 / areaSum);
    const Vec3 loc = axis.loc;

    std::vector<double> psi;
    psi.reserve(cylNormals.size());
    for (const Vec3& n : cylNormals) {
        psi.push_back(std::atan2(dot(n, v0), dot(n, u0)));
    }
    std::sort(psi.begin(), psi.end());
    std::vector<double> gaps;
    for (size_t i = 1; i < psi.size(); ++i) gaps.push_back(psi[i] - psi[i - 1]);
    if (!psi.empty()) gaps.push_back((psi.front() + 2.0 * M_PI) - psi.back());
    // D2.2: drop the single largest circular gap BEFORE the P75. The largest
    // gap is the uncovered span (or the wrap) and must not inflate the bin.
    std::vector<double> sortedGaps = gaps;
    std::sort(sortedGaps.begin(), sortedGaps.end());
    if (sortedGaps.size() > 1) sortedGaps.pop_back();
    double p75 = 0;
    if (!sortedGaps.empty()) {
        const size_t idx = static_cast<size_t>(0.75 * (sortedGaps.size() - 1));
        p75 = sortedGaps[idx];
    }
    const double thetaBin = std::max(M_PI / 720.0, 0.5 * p75);
    int nBands = 1;
    for (size_t i = 1; i < psi.size(); ++i) {
        if (psi[i] - psi[i - 1] > thetaBin) ++nBands;
    }

    std::vector<double> chi;
    chi.reserve(cylVerts.size());
    for (const Vec3& p : cylVerts) {
        const Vec3 d = sub(p, loc);
        chi.push_back(std::atan2(dot(d, v0), dot(d, u0)));
    }
    std::sort(chi.begin(), chi.end());
    double maxGap = 0;
    for (size_t i = 1; i < chi.size(); ++i) maxGap = std::max(maxGap, chi[i] - chi[i - 1]);
    if (!chi.empty()) maxGap = std::max(maxGap, (chi.front() + 2.0 * M_PI) - chi.back());
    const bool closed360 = (nBands >= 3) && (maxGap <= 1.5 * (2.0 * M_PI / nBands));
    const double span = closed360 ? (2.0 * M_PI) : (2.0 * M_PI - maxGap);
    if (span <= 0 || nBands < 1) return 0;
    return static_cast<int>(std::llround(2.0 * M_PI * nBands / span));
}

inline bool axesParallel(const Vec3& d1, const Vec3& d2, double angTolDeg = 3.0) {
    const double c = std::fabs(dot(normalize(d1), normalize(d2)));
    return c >= std::cos(angTolDeg * M_PI / 180.0);
}

inline bool pointOnAxis(const Vec3& p, const AxisSpec& ax, double tol) {
    const Vec3 d = sub(p, ax.loc);
    return norm(cross(d, ax.dir)) <= tol;
}

inline int measureRecoverableNSides(const MeshData& mesh, const Recoverable& r) {
    if (r.type != "cylinder") return 0;
    return measureCylinderNSides(mesh, r.axis, r.radius, 0.08, r.vMin, r.vMax);
}

}  // namespace corpus

#endif
