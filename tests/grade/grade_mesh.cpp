#include "grade_mesh.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <sstream>

namespace grade {
namespace {

uint32_t f32bits(float f) {
    uint32_t u;
    std::memcpy(&u, &f, 4);
    return u;
}

float bitsF32(uint32_t u) {
    float f;
    std::memcpy(&f, &u, 4);
    return f;
}

struct D3 {
    double x, y, z;
    bool operator==(const D3& o) const { return x == o.x && y == o.y && z == o.z; }
};
struct D3Hash {
    size_t operator()(const D3& k) const {
        uint64_t ax, ay, az;
        std::memcpy(&ax, &k.x, 8);
        std::memcpy(&ay, &k.y, 8);
        std::memcpy(&az, &k.z, 8);
        size_t h = static_cast<size_t>(ax);
        h ^= static_cast<size_t>(ay) + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= static_cast<size_t>(az) + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};

void fillMetrics(Mesh& m) {
    m.surfaceArea = 0;
    m.volume = 0;
    m.meanCircumdiam = 0;
    m.adj.assign(m.tris.size(), {});
    m.vertTris.assign(m.verts.size(), {});
    double circumSum = 0;
    int circumN = 0;
    for (size_t ti = 0; ti < m.tris.size(); ++ti) {
        Tri& t = m.tris[ti];
        const Vec3& a = m.verts[static_cast<size_t>(t.v[0])];
        const Vec3& b = m.verts[static_cast<size_t>(t.v[1])];
        const Vec3& c = m.verts[static_cast<size_t>(t.v[2])];
        const Vec3 ab = b - a, ac = c - a, bc = c - b;
        const Vec3 cr = cross(ab, ac);
        const double twice = norm(cr);
        t.area = 0.5 * twice;
        t.n = twice > 0.0 ? cr * (1.0 / twice) : Vec3{0, 0, 1};
        t.centroid = (a + b + c) * (1.0 / 3.0);
        const double la = norm(bc), lb = norm(ac), lc = norm(ab);
        t.perim = la + lb + lc;
        const double longest = std::max(la, std::max(lb, lc));
        t.hmin = longest > 0.0 ? (2.0 * t.area) / longest : 0.0;
        t.thetaQ = std::asin(std::min(1.0, (t.hmin > 0.0) ? (2.0 * m.q / t.hmin) : 1.0));
        t.circumdiam = (t.area > 0.0) ? (la * lb * lc) / (2.0 * t.area) : 0.0;
        m.surfaceArea += t.area;
        m.volume += dot(a, cross(b, c));
        if (t.circumdiam > 0.0) {
            circumSum += t.circumdiam;
            ++circumN;
        }
        m.vertTris[static_cast<size_t>(t.v[0])].push_back(static_cast<int>(ti));
        m.vertTris[static_cast<size_t>(t.v[1])].push_back(static_cast<int>(ti));
        m.vertTris[static_cast<size_t>(t.v[2])].push_back(static_cast<int>(ti));
    }
    m.volume /= 6.0;
    m.meanCircumdiam = circumN ? circumSum / circumN : 1.0;

    struct EH {
        size_t operator()(uint64_t k) const { return static_cast<size_t>(k ^ (k >> 32)); }
    };
    std::unordered_map<uint64_t, std::vector<int>, EH> e2t;
    e2t.reserve(m.tris.size() * 2);
    auto ek = [](int i, int j) -> uint64_t {
        if (i > j) std::swap(i, j);
        return (uint64_t(uint32_t(i)) << 32) | uint32_t(j);
    };
    for (size_t ti = 0; ti < m.tris.size(); ++ti) {
        const Tri& t = m.tris[ti];
        for (int s = 0; s < 3; ++s) {
            e2t[ek(t.v[s], t.v[(s + 1) % 3])].push_back(static_cast<int>(ti));
        }
    }
    m.openEdges = 0;
    m.nonManifoldEdges = 0;
    for (auto& kv : e2t) {
        auto& ts = kv.second;
        if (ts.size() == 1) {
            ++m.openEdges;
        } else if (ts.size() == 2) {
            m.adj[static_cast<size_t>(ts[0])].push_back(ts[1]);
            m.adj[static_cast<size_t>(ts[1])].push_back(ts[0]);
        } else {
            ++m.nonManifoldEdges;
        }
    }

    // q-scaled cells (D-140-1 / cycle-2 F12): nearestVert is a q-ball query.
    // Mean circumdiameter buckets scan hundreds of verts per face vertex.
    m.gridCell = (m.q > 0.0) ? m.q : (m.meanCircumdiam > 0.0 ? m.meanCircumdiam : 1.0);
    m.grid.clear();
    auto celli = [&](double x) -> int64_t {
        return static_cast<int64_t>(std::floor(x / m.gridCell));
    };
    auto key = [](int64_t ix, int64_t iy, int64_t iz) -> int64_t {
        return (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791);
    };
    for (int i = 0; i < static_cast<int>(m.verts.size()); ++i) {
        const Vec3& p = m.verts[static_cast<size_t>(i)];
        m.grid[key(celli(p.x), celli(p.y), celli(p.z))].push_back(i);
    }
    for (auto& kv : m.grid) {
        std::sort(kv.second.begin(), kv.second.end());
    }

    m.bitIndex.clear();
    for (int i = 0; i < static_cast<int>(m.verts.size()); ++i) {
        const Vec3& p = m.verts[static_cast<size_t>(i)];
        Mesh::Bit3 b{f32bits(static_cast<float>(p.x)), f32bits(static_cast<float>(p.y)),
                     f32bits(static_cast<float>(p.z))};
        m.bitIndex.emplace(b, i);
    }
}

bool loadBinary(std::FILE* f, uint32_t nTri, Mesh& m, std::string& err) {
    std::unordered_map<Mesh::Bit3, int, Mesh::Bit3Hash> weld;
    weld.reserve(nTri);
    m.verts.clear();
    m.tris.clear();
    m.degenerateDropped = 0;
    std::vector<unsigned char> buf(50);
    for (uint32_t i = 0; i < nTri; ++i) {
        if (std::fread(buf.data(), 1, 50, f) != 50) {
            err = "truncated binary STL";
            return false;
        }
        int idx[3];
        for (int k = 0; k < 3; ++k) {
            uint32_t bx, by, bz;
            std::memcpy(&bx, buf.data() + 12 + 12 * k, 4);
            std::memcpy(&by, buf.data() + 16 + 12 * k, 4);
            std::memcpy(&bz, buf.data() + 20 + 12 * k, 4);
            Mesh::Bit3 key{bx, by, bz};
            auto it = weld.find(key);
            if (it == weld.end()) {
                const int id = static_cast<int>(m.verts.size());
                weld.emplace(key, id);
                m.verts.push_back(Vec3{static_cast<double>(bitsF32(bx)),
                                       static_cast<double>(bitsF32(by)),
                                       static_cast<double>(bitsF32(bz))});
                idx[k] = id;
            } else {
                idx[k] = it->second;
            }
        }
        if (idx[0] == idx[1] || idx[1] == idx[2] || idx[2] == idx[0]) {
            ++m.degenerateDropped;
            continue;
        }
        Tri t;
        t.v[0] = idx[0];
        t.v[1] = idx[1];
        t.v[2] = idx[2];
        m.tris.push_back(t);
    }
    return true;
}

bool loadAscii(std::FILE* f, Mesh& m, std::string& err) {
    std::unordered_map<D3, int, D3Hash> weld;
    m.verts.clear();
    m.tris.clear();
    m.degenerateDropped = 0;
    std::rewind(f);
    char line[1024];
    int pending[3];
    int got = 0;
    while (std::fgets(line, sizeof line, f)) {
        const char* p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r') ++p;
        if (std::strncmp(p, "vertex", 6) != 0) continue;
        p += 6;
        char* end = nullptr;
        const double x = std::strtod(p, &end);
        if (end == p) continue;
        p = end;
        const double y = std::strtod(p, &end);
        if (end == p) continue;
        p = end;
        const double z = std::strtod(p, &end);
        if (end == p) continue;
        D3 key{x, y, z};
        auto it = weld.find(key);
        int id;
        if (it == weld.end()) {
            id = static_cast<int>(m.verts.size());
            weld.emplace(key, id);
            m.verts.push_back(Vec3{x, y, z});
        } else {
            id = it->second;
        }
        pending[got++] = id;
        if (got == 3) {
            if (pending[0] == pending[1] || pending[1] == pending[2] || pending[2] == pending[0]) {
                ++m.degenerateDropped;
            } else {
                Tri t;
                t.v[0] = pending[0];
                t.v[1] = pending[1];
                t.v[2] = pending[2];
                m.tris.push_back(t);
            }
            got = 0;
        }
    }
    (void)err;
    return !m.tris.empty();
}

}  // namespace

int Mesh::nearestVert(const Vec3& p, double maxDist) const {
    Mesh::Bit3 b{f32bits(static_cast<float>(p.x)), f32bits(static_cast<float>(p.y)),
                 f32bits(static_cast<float>(p.z))};
    auto it = bitIndex.find(b);
    if (it != bitIndex.end()) return it->second;
    if (!(gridCell > 0.0) || grid.empty()) {
        int best = -1;
        double bestD = maxDist;
        for (int i = 0; i < static_cast<int>(verts.size()); ++i) {
            const double d = dist(verts[static_cast<size_t>(i)], p);
            if (d <= bestD) {
                bestD = d;
                best = i;
            }
        }
        return best;
    }
    auto celli = [&](double x) -> int64_t { return static_cast<int64_t>(std::floor(x / gridCell)); };
    auto key = [](int64_t ix, int64_t iy, int64_t iz) -> int64_t {
        return (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791);
    };
    const int64_t ix = celli(p.x), iy = celli(p.y), iz = celli(p.z);
    int best = -1;
    double bestD = maxDist;
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
                auto git = grid.find(key(ix + dx, iy + dy, iz + dz));
                if (git == grid.end()) continue;
                for (int vi : git->second) {
                    const double d = dist(verts[static_cast<size_t>(vi)], p);
                    if (d <= bestD) {
                        bestD = d;
                        best = vi;
                    }
                }
            }
    return best;
}

bool loadStl(const std::string& path, Mesh& mesh, std::string& err) {
    mesh = Mesh{};
    mesh.path = path;
    mesh.qf = stl2step::stlQuantFloor(path);
    if (!mesh.qf.ok || !(mesh.qf.q > 0.0)) {
        err = "q is 0 (unreadable STL or no coordinates)";
        return false;
    }
    mesh.q = mesh.qf.q;
    mesh.tau = 2.0 * mesh.q;

    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) {
        err = "cannot open STL";
        return false;
    }
    if (std::fseek(f, 0, SEEK_END) != 0) {
        std::fclose(f);
        err = "STL seek failed";
        return false;
    }
    const long long size = static_cast<long long>(std::ftell(f));
    std::rewind(f);
    bool binary = false;
    uint32_t nTri = 0;
    if (size >= 84) {
        unsigned char head[84];
        if (std::fread(head, 1, 84, f) == 84) {
            nTri = static_cast<uint32_t>(head[80]) | (static_cast<uint32_t>(head[81]) << 8) |
                   (static_cast<uint32_t>(head[82]) << 16) | (static_cast<uint32_t>(head[83]) << 24);
            if (size == 84LL + 50LL * static_cast<long long>(nTri)) binary = true;
        }
        std::rewind(f);
    }
    bool ok;
    if (binary) {
        if (std::fseek(f, 84, SEEK_SET) != 0) {
            std::fclose(f);
            err = "binary STL seek failed";
            return false;
        }
        ok = loadBinary(f, nTri, mesh, err);
    } else {
        ok = loadAscii(f, mesh, err);
        if (!ok) err = "ASCII STL produced no triangles";
    }
    std::fclose(f);
    if (!ok) return false;
    fillMetrics(mesh);
    return true;
}

void eigen3(const double S[3][3], double eval[3], double evec[3][3]) {
    double A[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) A[i][j] = S[i][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) evec[i][j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 64; ++sweep) {
        int p = 0, q = 1;
        double maxv = std::fabs(A[0][1]);
        if (std::fabs(A[0][2]) > maxv) {
            maxv = std::fabs(A[0][2]);
            p = 0;
            q = 2;
        }
        if (std::fabs(A[1][2]) > maxv) {
            maxv = std::fabs(A[1][2]);
            p = 1;
            q = 2;
        }
        if (maxv == 0.0) break;
        const double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = A[q][p] = 0;
        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q) continue;
            const double akp = A[k][p], akq = A[k][q];
            A[k][p] = A[p][k] = c * akp - s * akq;
            A[k][q] = A[q][k] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
            const double eip = evec[k][p], eiq = evec[k][q];
            evec[k][p] = c * eip - s * eiq;
            evec[k][q] = s * eip + c * eiq;
        }
    }
    eval[0] = A[0][0];
    eval[1] = A[1][1];
    eval[2] = A[2][2];
}

bool solveN(int n, double* A, double* b, double* x) {
    const int N = n;
    std::vector<int> piv(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i) piv[static_cast<size_t>(i)] = i;
    for (int k = 0; k < N; ++k) {
        int best = k;
        double bestv = std::fabs(A[k * N + k]);
        for (int i = k + 1; i < N; ++i) {
            const double v = std::fabs(A[i * N + k]);
            if (v > bestv) {
                bestv = v;
                best = i;
            }
        }
        if (bestv == 0.0) return false;
        if (best != k) {
            for (int j = 0; j < N; ++j) std::swap(A[k * N + j], A[best * N + j]);
            std::swap(b[k], b[best]);
        }
        for (int i = k + 1; i < N; ++i) {
            const double f = A[i * N + k] / A[k * N + k];
            b[i] -= f * b[k];
            for (int j = k; j < N; ++j) A[i * N + j] -= f * A[k * N + j];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < N; ++j) s -= A[i * N + j] * x[j];
        x[i] = s / A[i * N + i];
    }
    return true;
}

std::string relPath(const std::string& stlPath, const std::string& other) {
    std::string dir = stlPath;
    const auto slash = dir.find_last_of("/\\");
    if (slash == std::string::npos) dir.clear();
    else dir.resize(slash + 1);
    std::string o = other;
    for (char& c : o)
        if (c == '\\') c = '/';
    std::string d = dir;
    for (char& c : d)
        if (c == '\\') c = '/';
    if (!d.empty() && o.size() >= d.size() && o.compare(0, d.size(), d) == 0) {
        o = o.substr(d.size());
    }
    return o;
}

std::string fmtG(double v) {
    if (v == 0.0) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.9g", v);
    return buf;
}

}  // namespace grade
