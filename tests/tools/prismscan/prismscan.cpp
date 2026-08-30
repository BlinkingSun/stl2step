// Corpus routing census: detectPrismatic over every body/component.
// Emits _team/reports/ac3/prismscan.json (SPEC §6).
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "posix_compat.hpp"
#include "refit.hpp"
#include "refit_prism.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dirent.h>
#endif

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::MeshView;
using stl2step::refit::PrismLevels;
using stl2step::refit::PrismTols;
using stl2step::refit::Region;
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;
using stl2step::refit::detectPrismatic;

namespace stl2step {
namespace refit {
void derivePrismTols(const MeshView& mv, const RegionSet& rs, PrismTols& t);
}
}

namespace {

struct CompRow {
    int index = 0;
    bool ok = false;
    int failedCond = 1;
    int nCyl = 0, nCap = 0, nLat = 0, nOblique = 0;
    std::vector<double> levels;
    std::vector<double> spacings;
    double tauAx = 0.0;
    double tauLvl = 0.0;
};

struct BodyRow {
    std::string file;
    int components = 0;
    std::string route;
    std::vector<CompRow> comps;
};

std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (char c : s) {
        if (c == '\\' || c == '"') {
            o.push_back('\\');
            o.push_back(c);
        } else {
            o.push_back(c);
        }
    }
    return o;
}

std::string relPath(const std::string& path, const std::string& repo) {
    if (path.size() >= repo.size() && path.compare(0, repo.size(), repo) == 0) {
        size_t i = repo.size();
        while (i < path.size() && (path[i] == '/' || path[i] == '\\')) ++i;
        return path.substr(i);
    }
    return path;
}

std::vector<std::string> listStl(const std::string& dir) {
    std::vector<std::string> out;
#if defined(_WIN32)
    WIN32_FIND_DATAA fd;
    const std::string pat = dir + "\\*.stl";
    HANDLE h = FindFirstFileA(pat.c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return out;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
            out.push_back(dir + "/" + fd.cFileName);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR* d = opendir(dir.c_str());
    if (!d) return out;
    while (dirent* e = readdir(d)) {
        const char* n = e->d_name;
        const size_t len = std::strlen(n);
        if (len < 4) continue;
        if (n[len - 4] == '.' && (n[len - 3] == 's' || n[len - 3] == 'S') &&
            (n[len - 2] == 't' || n[len - 2] == 'T') &&
            (n[len - 1] == 'l' || n[len - 1] == 'L'))
            out.push_back(dir + "/" + n);
    }
    closedir(d);
#endif
    std::sort(out.begin(), out.end());
    return out;
}

std::string gitSha(const std::string& repo) {
    const std::string cmd = "git -C \"" + repo + "\" rev-parse HEAD 2>/dev/null";
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) return "unknown";
    char buf[128] = {};
    const char* got = std::fgets(buf, sizeof buf, p);
    pclose(p);
    if (!got) return "unknown";
    std::string s(buf);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s.empty() ? "unknown" : s;
}

void classifyPlanes(const RegionSet& rs, const PrismLevels& lv, const PrismTols& t,
                    int failedCond, int& nCap, int& nLat, int& nOblique) {
    nCap = nLat = nOblique = 0;
    // Axis is only established after condition 2; do not invent a split before that.
    if (failedCond == 1 || failedCond == 2) return;
    const gp_XYZ a = lv.axis.XYZ();
    const double tau = t.tauAx > 0.0 ? t.tauAx : 1e-6;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Plane) continue;
        const double nd = std::abs(r.ax.Direction().XYZ().Dot(a));
        if (nd > 1.0 - tau) ++nCap;
        else if (nd < tau) ++nLat;
        else ++nOblique;
    }
}

CompRow scanComp(const HarnessMesh& mesh, size_t ci) {
    CompRow row;
    row.index = static_cast<int>(ci);
    const auto& c = mesh.comps[ci];
    MeshView mv{};
    toMeshView(mesh, c, mv);
    RegionSet rs;
    rs.compRoot = static_cast<int>(ci);
    if (!c.clean()) {
        PrismTols t;
        stl2step::refit::derivePrismTols(mv, rs, t);
        row.tauAx = t.tauAx;
        row.tauLvl = t.tauLvl;
        row.failedCond = 1;
        return row;
    }
    SegmentParams params;
    (void)stl2step::refit::segment(mv, params, rs, nullptr);
    rs.compRoot = static_cast<int>(ci);
    PrismTols t;
    stl2step::refit::derivePrismTols(mv, rs, t);
    row.tauAx = t.tauAx;
    row.tauLvl = t.tauLvl;
    for (const Region& r : rs.regions)
        if (r.type == SurfType::Cylinder) ++row.nCyl;
    const PrismLevels lv = detectPrismatic(mv, rs, t);
    row.ok = lv.ok;
    row.failedCond = lv.failedCond;
    classifyPlanes(rs, lv, t, row.failedCond, row.nCap, row.nLat, row.nOblique);
    row.levels = lv.y;
    for (size_t i = 1; i < lv.y.size(); ++i)
        row.spacings.push_back(lv.y[i] - lv.y[i - 1]);
    return row;
}

static std::mutex gOcctMu;

BodyRow scanBody(const std::string& path, const std::string& repo) {
    BodyRow body;
    body.file = relPath(path, repo);
    HarnessMesh mesh;
    std::string err;
    std::lock_guard<std::mutex> occt(gOcctMu);
    if (!loadMesh(path, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        body.components = 0;
        body.route = "NON-PRISMATIC";
        CompRow c;
        c.failedCond = 1;
        body.comps.push_back(c);
        return body;
    }
    body.components = static_cast<int>(mesh.comps.size());
    body.comps.resize(mesh.comps.size());
    // segment() is OCCT-backed; keep it on this thread (caller holds gOcctMu).
    for (size_t i = 0; i < mesh.comps.size(); ++i)
        body.comps[i] = scanComp(mesh, i);
    int nOk = 0;
    for (const CompRow& c : body.comps)
        if (c.ok) ++nOk;
    if (nOk == body.components && nOk > 0) body.route = "PRISMATIC";
    else if (nOk > 0) body.route = "PARTIALLY-PRISMATIC";
    else body.route = "NON-PRISMATIC";
    return body;
}

}  // namespace

int main(int argc, char** argv) {
    std::string corpus = "tests/corpus";
    std::string outPath = "_team/reports/ac3/prismscan.json";
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--corpus") == 0 && i + 1 < argc)
            corpus = argv[++i];
        else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc)
            outPath = argv[++i];
    }

    std::string repo = corpus;
    const auto pos = repo.find("tests/corpus");
    if (pos != std::string::npos) repo = repo.substr(0, pos);
    while (!repo.empty() && (repo.back() == '/' || repo.back() == '\\')) repo.pop_back();
    if (repo.empty()) repo = ".";

    const std::vector<std::string> files = listStl(corpus);
    std::vector<BodyRow> bodies(files.size());
    std::atomic<size_t> next{0};
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    const unsigned w = static_cast<unsigned>(std::min<size_t>(hw, files.empty() ? 1 : files.size()));
    std::vector<std::thread> pool;
    pool.reserve(w);
    for (unsigned k = 0; k < w; ++k) {
        pool.emplace_back([&]() {
            for (;;) {
                const size_t i = next.fetch_add(1);
                if (i >= files.size()) break;
                bodies[i] = scanBody(files[i], repo);
            }
        });
    }
    for (auto& th : pool) th.join();

    FILE* f = std::fopen(outPath.c_str(), "wb");
    if (!f) {
        std::fprintf(stderr, "prismscan: cannot write %s\n", outPath.c_str());
        return 1;
    }
    const std::string sha = gitSha(repo);
    std::fprintf(f, "{\n  \"ok\": true,\n  \"commit\": \"%s\",\n  \"bodies\": [\n",
                 jsonEscape(sha).c_str());
    for (size_t b = 0; b < bodies.size(); ++b) {
        const BodyRow& body = bodies[b];
        std::fprintf(f, "    { \"file\": \"%s\", \"components\": %d,\n",
                     jsonEscape(body.file).c_str(), body.components);
        std::fprintf(f, "      \"route\": \"%s\",\n", body.route.c_str());
        std::fprintf(f, "      \"comps\": [");
        for (size_t i = 0; i < body.comps.size(); ++i) {
            const CompRow& c = body.comps[i];
            if (i) std::fprintf(f, ",");
            std::fprintf(f, "\n        { \"index\": %d, \"ok\": %s, \"failedCond\": %d,\n",
                         c.index, c.ok ? "true" : "false", c.failedCond);
            std::fprintf(f, "          \"nCyl\": %d, \"nCap\": %d, \"nLat\": %d, \"nOblique\": %d,\n",
                         c.nCyl, c.nCap, c.nLat, c.nOblique);
            std::fprintf(f, "          \"levels\": [");
            for (size_t k = 0; k < c.levels.size(); ++k) {
                if (k) std::fprintf(f, ", ");
                std::fprintf(f, "%.6f", c.levels[k]);
            }
            std::fprintf(f, "],\n          \"spacings\": [");
            for (size_t k = 0; k < c.spacings.size(); ++k) {
                if (k) std::fprintf(f, ", ");
                std::fprintf(f, "%.6f", c.spacings[k]);
            }
            std::fprintf(f, "],\n          \"tauAx\": %.3e, \"tauLvl\": %.3e }",
                         c.tauAx, c.tauLvl);
        }
        std::fprintf(f, " ] }");
        if (b + 1 < bodies.size()) std::fprintf(f, ",");
        std::fprintf(f, "\n");
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    std::fprintf(stderr, "prismscan: wrote %zu bodies to %s\n", bodies.size(), outPath.c_str());
    return 0;
}
