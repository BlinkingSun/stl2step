// D-142-1: threads=1 constructs zero std::thread; threads=4 spawns and
// produces identical canonical STEP DATA (HEADER discarded).
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using stl2step::MeshOptions;
using stl2step::Options;
using stl2step::detail::resetThreadsSpawnedForTest;
using stl2step::detail::threadsSpawnedForTest;

static int gFail = 0;

static void check(bool cond, const char* what) {
    if (cond) return;
    std::fprintf(stderr, "FAIL %s\n", what);
    ++gFail;
}

static std::string readAll(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

static std::string dataSection(const std::string& step) {
    // Geometry starts at the B-Rep representation. PRODUCT/AP context before
    // that carries the output stem and OCCT's per-process " 2" suffix.
    const auto pos = step.find("ADVANCED_BREP_SHAPE_REPRESENTATION");
    return pos == std::string::npos ? step : step.substr(pos);
}

static bool convertSmooth(const fs::path& stl, const fs::path& out, int threads) {
    Options opt;
    opt.input = stl.string();
    opt.output = out.string();
    opt.smooth = true;
    opt.threads = threads;
    opt.verify = false;
    opt.productName = "part";  // output stem is not geometry; keep PRODUCT identical
    resetThreadsSpawnedForTest();
    auto r = stl2step::convert(opt, nullptr);
    if (!r.ok) {
        std::fprintf(stderr, "convert failed %s: %s\n", stl.c_str(), r.error.c_str());
        return false;
    }
    return true;
}

static int runNoSpawn(const fs::path& corpus, const fs::path& meshStep, const fs::path& scratch) {
    fs::create_directories(scratch);
    const char* names[] = {"shelf_bracket.stl", "Body244.stl", "handle-pickup.stl"};
    for (const char* name : names) {
        const fs::path stl = corpus / name;
        if (!fs::exists(stl)) {
            std::fprintf(stderr, "FAIL missing %s\n", stl.c_str());
            return 1;
        }
        const fs::path out = scratch / (std::string(name) + ".t1.step");
        if (!convertSmooth(stl, out, 1)) return 1;
        const auto spawned = threadsSpawnedForTest();
        std::fprintf(stderr, "no_spawn convert %s spawned=%llu\n", name,
                     (unsigned long long)spawned);
        check(spawned == 0, "convert threads=1 spawned==0");
    }

    if (!fs::exists(meshStep)) {
        std::fprintf(stderr, "FAIL missing mesh fixture %s\n", meshStep.c_str());
        return 1;
    }
    MeshOptions mopt;
    mopt.input = meshStep.string();
    mopt.output = (scratch / "mesh.t1.stl").string();
    mopt.threads = 1;
    resetThreadsSpawnedForTest();
    auto mr = stl2step::meshFromStep(mopt, nullptr);
    if (!mr.ok) {
        std::fprintf(stderr, "meshFromStep failed: %s\n", mr.error.c_str());
        return 1;
    }
    const auto spawned = threadsSpawnedForTest();
    std::fprintf(stderr, "no_spawn mesh spawned=%llu\n", (unsigned long long)spawned);
    check(spawned == 0, "mesh threads=1 spawned==0");
    return gFail ? 1 : 0;
}

static int runDefaultSpawns(const fs::path& corpus, const fs::path& /*meshStep*/,
                            const fs::path& scratch) {
    fs::create_directories(scratch);
    const fs::path stl = corpus / "handle-pickup.stl";
    if (!fs::exists(stl)) {
        std::fprintf(stderr, "FAIL missing %s\n", stl.c_str());
        return 1;
    }
    // Same output stem both runs so OCCT PRODUCT/FILE_NAME in DATA match.
    const fs::path dir4 = scratch / "t4";
    const fs::path dir1 = scratch / "t1";
    fs::create_directories(dir4);
    fs::create_directories(dir1);
    const fs::path t4 = dir4 / "handle-pickup.step";
    const fs::path t1 = dir1 / "handle-pickup.step";
    if (!convertSmooth(stl, t4, 4)) return 1;
    const auto spawned4 = threadsSpawnedForTest();
    std::fprintf(stderr, "default_spawns threads=4 spawned=%llu\n",
                 (unsigned long long)spawned4);
    check(spawned4 > 0, "convert threads=4 spawned>0");
    if (!convertSmooth(stl, t1, 1)) return 1;
    check(threadsSpawnedForTest() == 0, "convert threads=1 spawned==0 (identity run)");
    const std::string a = dataSection(readAll(t4));
    const std::string b = dataSection(readAll(t1));
    check(a == b, "canonical STEP DATA threads=4 vs threads=1");
    if (a != b)
        std::fprintf(stderr, "DATA length t4=%zu t1=%zu\n", a.size(), b.size());
    return gFail ? 1 : 0;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: %s no_spawn|default_spawns <corpus> <mesh.step> <scratch>\n",
                     argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    const fs::path corpus = argv[2];
    const fs::path meshStep = argv[3];
    const fs::path scratch = argv[4];
    if (mode == "no_spawn") return runNoSpawn(corpus, meshStep, scratch);
    if (mode == "default_spawns") return runDefaultSpawns(corpus, meshStep, scratch);
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
