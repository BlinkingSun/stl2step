// SPEC AC-CASCADE — exercise stl2step::convert() (the public C++ API).
// Not a ctest (lane pins 27/27 + hl_census_ratchet + cascade_no_waiver).
//
// Modes via CASCADE_EXPECT (default unset):
//   unset — C1: built >= 10, watertight, volΔ <= 6.013, warnings empty
//   11    — C3(ii): FAIL_RID=11 already in the environment
//   0     — C3(iii): FAIL_RID=0; must not collapse to the old blanket
//   bad   — malformed FAIL_RID; convert still runs, one warning
//
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

namespace {

const double kVolBudget = 6.013;

int runConvert(const std::string& stl, const std::string& step, stl2step::Result& r,
               std::vector<std::string>& extra) {
    stl2step::Options opt;
    opt.input = stl;
    opt.output = step;
    opt.smooth = true;
    opt.verify = true;
    std::mutex logMu;
    auto log = [&](stl2step::Severity sev, const std::string& msg) {
        if (sev == stl2step::Severity::Warning || sev == stl2step::Severity::Error) {
            std::lock_guard<std::mutex> g(logMu);
            extra.push_back(msg);
        }
    };
    r = stl2step::convert(opt, log);
    return r.ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    const char* expect = std::getenv("CASCADE_EXPECT");
    if (!expect || !expect[0]) expect = "unset";

    const char* stl = std::getenv("CASCADE_STL");
    if (!stl || !stl[0]) {
        if (argc >= 2) stl = argv[1];
        else
            stl = "tests/corpus/handle-lock.stl";
    }
    const char* outDir = std::getenv("CASCADE_OUT");
    if (!outDir || !outDir[0]) outDir = "/tmp/cascade-r4-api";
    std::string mkdir = std::string("mkdir -p '") + outDir + "'";
    if (std::system(mkdir.c_str()) != 0) return 2;

    std::string step = std::string(outDir) + "/hl-" + expect + ".step";
    stl2step::Result r;
    std::vector<std::string> extra;
    const int rc = runConvert(stl, step, r, extra);

    const int built = r.smoothBuiltCylinders;
    const double vol = r.volumeDeltaPct;
    const int nWarn = (int)r.warnings.size();
    std::printf(
        "CASCADE_API mode=%s ok=%d exit=%d builtCyl=%d volΔ=%.6f wt=%d open=%d "
        "warn=%d err=%s\n",
        expect, (int)r.ok, r.exitCode, built, vol, (int)r.watertight, r.openShells, nWarn,
        r.error.c_str());
    for (const auto& w : r.warnings)
        std::printf("  WARN %s\n", w.c_str());
    for (const auto& w : extra)
        std::printf("  LOG %s\n", w.c_str());

    if (std::strcmp(expect, "bad") == 0) {
        bool saw = false;
        for (const auto& w : r.warnings)
            if (w.find("STL2STEP_FAIL_RID") != std::string::npos) saw = true;
        if (!saw) {
            std::fprintf(stderr, "CASCADE_API FAIL: malformed knob produced no warning\n");
            return 1;
        }
        if (!r.ok) {
            std::fprintf(stderr, "CASCADE_API FAIL: malformed knob aborted convert\n");
            return 1;
        }
        std::printf("CASCADE_API PASS malformed-disabled\n");
        return 0;
    }

    if (rc != 0 || !r.ok) {
        std::fprintf(stderr, "CASCADE_API FAIL: convert failed\n");
        return 1;
    }
    if (!r.watertight || r.openShells != 0) {
        std::fprintf(stderr, "CASCADE_API FAIL: not watertight\n");
        return 1;
    }
    if (vol < 0.0 || vol > kVolBudget + 1e-9) {
        std::fprintf(stderr, "CASCADE_API FAIL: volumeDeltaPct=%.6f > %.3f\n", vol, kVolBudget);
        return 1;
    }
    if (std::strcmp(expect, "0") == 0) {
        // Hub-plane control: a blanket-equivalent census of 1 is AC failure.
        if (built <= 1) {
            std::fprintf(stderr, "CASCADE_API FAIL: FAIL_RID=0 collapsed to blanket (built=%d)\n",
                         built);
            return 1;
        }
        std::printf("CASCADE_API PASS hub-plane built=%d (not blanket)\n", built);
        return 0;
    }
    if (std::strcmp(expect, "11") == 0) {
        if (built < 10) {
            std::fprintf(stderr, "CASCADE_API FAIL: FAIL_RID=11 built=%d < 10\n", built);
            return 1;
        }
        std::printf("CASCADE_API PASS inject-11 built=%d\n", built);
        return 0;
    }
    // unset / C1
    if (built < 10) {
        std::fprintf(stderr, "CASCADE_API FAIL: C1 built=%d < 10\n", built);
        return 1;
    }
    if (!r.warnings.empty()) {
        std::fprintf(stderr, "CASCADE_API FAIL: C1 warnings not empty\n");
        return 1;
    }
    std::printf("CASCADE_API PASS C1 built=%d volΔ=%.6f\n", built, vol);
    return 0;
}
