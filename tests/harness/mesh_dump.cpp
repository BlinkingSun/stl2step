// stl2step_meshharness — dump welded/split mesh tables as one JSON object.
//
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using stl2step::harness::HarnessComponent;
using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;

static uint64_t fnv1a64(const uint8_t* p, size_t n) {
    uint64_t h = 14695981039346656037ull;
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ull;
    }
    return h;
}

static std::string hex64(uint64_t h) {
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", (unsigned long long)h);
    return buf;
}

static std::string checksumTriEdges(const HarnessComponent& c) {
    static_assert(sizeof(std::array<int, 3>) == sizeof(int[3]),
                  "std::array<int,3> must be layout-compatible with int[3]");
    const uint8_t* p = reinterpret_cast<const uint8_t*>(c.triEdges.data());
    return hex64(fnv1a64(p, c.triEdges.size() * sizeof(std::array<int, 3>)));
}

static std::string checksumTriDirs(const HarnessComponent& c) {
    return hex64(fnv1a64(c.triDirs.data(), c.triDirs.size()));
}

static std::string jsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if ((unsigned char)c >= 0x20) o += c;
    }
    return o;
}

static std::string jsonNumber(double v) {
    if (!std::isfinite(v)) return "null";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

static void usage(const char* argv0) {
    std::fprintf(stderr,
                 "usage: %s <file.stl> [--scale S] [--weld-tol T] [--sew-tol T]\n",
                 argv0);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }
    std::string path;
    double scale = 1.0, weldTol = 0.0, sewTol = 0.0;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--scale") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 1; }
            scale = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--weld-tol") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 1; }
            weldTol = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--sew-tol") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 1; }
            sewTol = std::atof(argv[++i]);
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            return 1;
        } else if (path.empty()) {
            path = argv[i];
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (path.empty()) {
        usage(argv[0]);
        return 1;
    }

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(path, scale, weldTol, sewTol, mesh, err)) {
        std::fprintf(stderr, "stl2step_meshharness: %s\n", err.c_str());
        return 1;
    }

    std::printf("{\"file\":\"%s\",\"triangles\":%zu,\"vertices\":%zu,"
                "\"diag\":%s,\"weldTol\":%s,\"sewTol\":%s,\"components\":[",
                jsonEscape(path).c_str(),
                mesh.tris.size(), mesh.pts.size(),
                jsonNumber(mesh.diag).c_str(),
                jsonNumber(mesh.weldTol).c_str(),
                jsonNumber(mesh.sewTol).c_str());
    for (size_t i = 0; i < mesh.comps.size(); i++) {
        const HarnessComponent& c = mesh.comps[i];
        if (i) std::printf(",");
        std::printf("{\"root\":%d,\"tris\":%zu,\"vtx\":%zu,\"edges\":%zu,"
                    "\"open\":%d,\"conflict\":%d,\"nonManifold\":%d,"
                    "\"clean\":%s,\"volume\":%s,"
                    "\"triEdgesChecksum\":\"%s\",\"triDirsChecksum\":\"%s\"}",
                    c.root,
                    c.compTris.size(), c.compVtx.size(), c.compEdges.size(),
                    c.open, c.conflict, c.nonManifold,
                    c.clean() ? "true" : "false",
                    jsonNumber(c.vol).c_str(),
                    checksumTriEdges(c).c_str(),
                    checksumTriDirs(c).c_str());
    }
    std::printf("]}\n");
    return 0;
}
