#include "grade_compare.hpp"
#include "grade_report.hpp"

#include <clocale>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

std::string dirOf(const std::string& p) {
    const auto sl = p.find_last_of("/\\");
    if (sl == std::string::npos) return std::string();
    return p.substr(0, sl);
}
std::string stemOf(const std::string& p) {
    std::string s = p;
    const auto sl = s.find_last_of("/\\");
    if (sl != std::string::npos) s = s.substr(sl + 1);
    const auto dot = s.find_last_of('.');
    if (dot != std::string::npos) s = s.substr(0, dot);
    return s;
}
std::string join(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    const char c = a.back();
    if (c == '/' || c == '\\') return a + b;
    return a + "/" + b;
}

}  // namespace

int main(int argc, char** argv) {
    std::setlocale(LC_ALL, "C");
    std::string stl, step, outDir, jsonPath, mdPath, engine;
    bool quiet = false;
    int seedOrder = 0;          // §8 case 9 third permutation (test-only probe)
    bool reverseSeeds = false;  // §8 case 9 growth-order probe, exposed for the
                                // training set (the selftest can only reach the
                                // corpus fixtures). Diagnostic only: it changes
                                // the seed walk, never a tolerance.
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (std::strcmp(a, "-o") == 0 && i + 1 < argc) {
            outDir = argv[++i];
        } else if (std::strcmp(a, "--json") == 0 && i + 1 < argc) {
            jsonPath = argv[++i];
        } else if (std::strcmp(a, "--md") == 0 && i + 1 < argc) {
            mdPath = argv[++i];
        } else if (std::strcmp(a, "--engine-bin") == 0 && i + 1 < argc) {
            engine = argv[++i];
        } else if (std::strcmp(a, "--quiet") == 0) {
            quiet = true;
        } else if (std::strcmp(a, "--reverse-seeds") == 0) {
            reverseSeeds = true;
        } else if (std::strcmp(a, "--seed-order") == 0 && i + 1 < argc) {
            seedOrder = std::atoi(argv[++i]);
        } else if (a[0] != '-') {
            if (stl.empty()) stl = a;
            else if (step.empty()) step = a;
        } else {
            std::fprintf(stderr, "unknown flag %s\n", a);
            return 1;
        }
    }
    if (stl.empty() || step.empty()) {
        std::fprintf(stderr,
                     "usage: stl2step_grade <stl> <step> [-o dir] [--json p] [--md p] "
                     "[--engine-bin p] [--quiet] [--reverse-seeds] [--seed-order n]\n");
        return 1;
    }
    const std::string stem = stemOf(step);
    if (jsonPath.empty() || mdPath.empty()) {
        std::string base = outDir.empty() ? dirOf(step) : outDir;
        if (base.empty()) base = ".";
        if (jsonPath.empty()) jsonPath = join(base, stem + ".grade.json");
        if (mdPath.empty()) mdPath = join(base, stem + ".grade.md");
    }

    grade::GradeConfig cfg;
    cfg.quiet = quiet;
    cfg.engineBin = engine;
    cfg.reverseSeeds = reverseSeeds;
    cfg.seedOrder = seedOrder;
    grade::GradeDocument doc;
    std::string err;
    if (!grade::gradeFiles(stl, step, cfg, doc, err)) {
        std::fprintf(stderr, "grade failed: %s\n", err.c_str());
        return 1;
    }
    const std::string js = grade::writeJson(doc);
    const std::string md = grade::writeMd(doc);
    if (!grade::writeFile(jsonPath, js) || !grade::writeFile(mdPath, md)) {
        std::fprintf(stderr, "failed to write grade artefacts\n");
        return 1;
    }
    if (!quiet) std::fputs(md.c_str(), stdout);
    return 0;
}
