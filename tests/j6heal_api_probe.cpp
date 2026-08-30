// SPEC AC-J6HEAL — exercise stl2step::convert() (the public C++ API).
// Not registered with ctest (AC requires 27/27 unchanged).
//
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace {

int parseFreeEdges(const std::vector<std::string>& warnings) {
    int last = -1;
    for (const std::string& w : warnings) {
        const char* p = std::strstr(w.c_str(), "freeEdges=");
        if (!p) continue;
        last = std::atoi(p + 10);
    }
    return last;
}

struct Case {
    const char* name;
    std::string stl;
    double volCeil;
    int freeCeil;  // -1 = must be closed (no J6)
};

bool fileExists(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return static_cast<bool>(in);
}

std::string joinCorpus(const char* root, const char* leaf) {
    std::string out(root);
    if (!out.empty() && out.back() != '/' && out.back() != '\\') out += '/';
    out += leaf;
    return out;
}

int runOne(const Case& c, const std::string& outDir, std::string& line) {
    stl2step::Options opt;
    opt.input = c.stl;
    opt.output = outDir + "/" + c.name + ".step";
    opt.smooth = true;
    opt.verify = true;
    std::mutex logMu;
    std::vector<std::string> extra;
    auto log = [&](stl2step::Severity sev, const std::string& msg) {
        if (sev == stl2step::Severity::Warning || sev == stl2step::Severity::Error) {
            std::lock_guard<std::mutex> g(logMu);
            extra.push_back(msg);
        }
    };
    stl2step::Result r = stl2step::convert(opt, log);
    const int fe = parseFreeEdges(r.warnings);
    char buf[512];
    std::snprintf(buf, sizeof(buf),
                  "%s ok=%d exit=%d volΔ=%.6f wt=%d builtCyl=%d rev=%d freeEdges=%d "
                  "warn=%zu err=%s",
                  c.name, (int)r.ok, r.exitCode, r.volumeDeltaPct, (int)r.watertight,
                  r.smoothBuiltCylinders, r.smoothRevertedComponents, fe, r.warnings.size(),
                  r.error.c_str());
    line = buf;
    if (!r.ok) return 1;
    if (c.volCeil >= 0.0 && r.volumeDeltaPct > c.volCeil + 1e-9) return 2;
    if (c.freeCeil < 0) {
        if (fe >= 0) return 3;
        if (!r.warnings.empty()) return 4;
    } else if (fe >= 0 && fe >= c.freeCeil) {
        return 5;
    }
    return 0;
}

}  // namespace

int main() {
    const char* corpus = std::getenv("STL2STEP_PRIVATE_CORPUS");
    if (!corpus || !corpus[0]) {
        std::fprintf(stderr, "SKIP: STL2STEP_PRIVATE_CORPUS unset\n");
        return 77;
    }

    const char* outDir = std::getenv("J6HEAL_OUT");
    if (!outDir || !outDir[0]) outDir = "/tmp/j6heal-r2-api";
    std::string mkdir = std::string("mkdir -p '") + outDir + "'";
    if (std::system(mkdir.c_str()) != 0) return 2;

    const Case wanted[] = {
        {"handle-lock", "tests/diag/handle-lock/handle-lock.stl", 0.0, -1},
        {"Body9", joinCorpus(corpus, "Body9.stl"), 0.030937, 34},
        {"Body12", joinCorpus(corpus, "Body12.stl"), 0.0, 81},
        {"Body18", joinCorpus(corpus, "Body18.stl"), 0.0, 18},
        {"Body20", joinCorpus(corpus, "Body20.stl"), 0.0, 32},
    };

    std::vector<Case> cases;
    cases.reserve(sizeof(wanted) / sizeof(wanted[0]));
    for (const Case& c : wanted) {
        if (!fileExists(c.stl)) {
            std::fprintf(stderr, "SKIP %s: missing %s\n", c.name, c.stl.c_str());
            continue;
        }
        cases.push_back(c);
    }
    if (cases.empty()) {
        std::fprintf(stderr, "SKIP: no probe STLs under STL2STEP_PRIVATE_CORPUS\n");
        return 77;
    }

    std::vector<std::string> lines(cases.size());
    std::vector<int> codes(cases.size(), 1);
    std::vector<std::thread> pool;
    pool.reserve(cases.size());
    for (size_t i = 0; i < cases.size(); i++) {
        pool.emplace_back([&, i]() { codes[i] = runOne(cases[i], outDir, lines[i]); });
    }
    for (auto& t : pool) t.join();

    int fail = 0;
    for (size_t i = 0; i < cases.size(); i++) {
        std::printf("%s rc=%d\n", lines[i].c_str(), codes[i]);
        if (codes[i] != 0) fail++;
    }
    return fail == 0 ? 0 : 1;
}
