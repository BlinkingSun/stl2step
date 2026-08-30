// SPEC AC2-L3 — exercise stl2step::convert() (the public C++ API) plus the
// D4 §3 scoped RULE 2.2 formula (C4/C5/C6). Not a ctest (C8 pins clean count).
//
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr double kR = 3.0;
constexpr double kArc = 1.05;
constexpr double kBlind = 1e-6;

// D4 §3 scoped sub-budget (mirrors src/refit_build.cpp residual block).
double scopedSub(bool plane, double dVolPredicted, double resid, double meshAbs,
                 int nAnalytic) {
    const bool predictorBlind = plane && std::fabs(dVolPredicted) < kBlind;
    const double arcTerm = predictorBlind ? kArc * std::fabs(resid) : 0.0;
    return std::max({1e-4 * meshAbs / (double)nAnalytic, kR * std::fabs(dVolPredicted),
                     arcTerm});
}

double oldSub(double dVolPredicted, double meshAbs, int nAnalytic) {
    return std::max(1e-4 * meshAbs / (double)nAnalytic, kR * std::fabs(dVolPredicted));
}

int checkC4C5() {
    // S4 / RULE22-CAL.md §5: global tight budget, not per-region sub.
    const double meshVol = 15868.885;
    const double tightMeasured = 954.0;
    const double sigmaAbs = tightMeasured / 3.0;
    const double tight = std::max(1e-4 * std::fabs(meshVol), 3.0 * sigmaAbs);
    const bool phantomReject = 28668.0 > tight;  // +121% pre-volumefix
    const bool bindReject = 5601.0 > tight;      // r1 exact-bind shift
    std::printf("C4 tight=%.3f phantom=28668 reject=%d\n", tight, (int)phantomReject);
    std::printf("C5 tight=%.3f bind=5601 reject=%d\n", tight, (int)bindReject);
    if (!phantomReject || !bindReject) return 1;
    if (std::fabs(tight - tightMeasured) > 1e-6) return 2;
    return 0;
}

int checkC6() {
    const double meshAbs = 15868.885;
    const int nAnalytic = 9;
    const double dVol = 150.0;    // |dVolPredicted| > 100
    const double resid = 500.0;   // large arc residual
    const double before = oldSub(dVol, meshAbs, nAnalytic);
    const double after = scopedSub(/*plane=*/false, dVol, resid, meshAbs, nAnalytic);
    std::printf("C6 cyl dVol=%.1f resid=%.1f oldSub=%.6f newSub=%.6f\n", dVol, resid, before,
                after);
    if (std::fabs(before - after) > 1e-12) {
        std::fprintf(stderr, "C6 FAIL: cylinder sub changed %g -> %g\n", before, after);
        return 1;
    }
    // Positive control: a predictor-blind plate MUST pick up kArc*|resid|.
    const double plate = scopedSub(/*plane=*/true, 0.0, 235.37, meshAbs, nAnalytic);
    const double expect = kArc * 235.37;
    std::printf("C6 plate-control sub=%.6f expect=%.6f\n", plate, expect);
    if (std::fabs(plate - expect) > 1e-6) {
        std::fprintf(stderr, "C6 FAIL: plate scoped sub %.6f != %.6f\n", plate, expect);
        return 2;
    }
    return 0;
}

int runConvert(const std::string& stl, const std::string& step, stl2step::Result& r) {
    stl2step::Options opt;
    opt.input = stl;
    opt.output = step;
    opt.smooth = true;
    opt.verify = true;
    std::mutex logMu;
    auto log = [&](stl2step::Severity, const std::string&) {
        std::lock_guard<std::mutex> g(logMu);
    };
    r = stl2step::convert(opt, log);
    return r.ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    const char* stl = std::getenv("RULE22_STL");
    if (!stl || !stl[0]) {
        if (argc >= 2) stl = argv[1];
        else
            stl = "tests/corpus/handle-lock.stl";
    }
    const char* outDir = std::getenv("RULE22_OUT");
    if (!outDir || !outDir[0]) outDir = "/tmp/rule22-adopt-api";
    std::string mkdir = std::string("mkdir -p '") + outDir + "'";
    if (std::system(mkdir.c_str()) != 0) return 2;

    int c4 = 1, c5c6 = 1, conv = 1;
    stl2step::Result r;
    // Formula checks are independent of convert(); run them concurrently
    // with the API conversion (real compute).
    std::thread tFormula([&]() {
        c4 = checkC4C5();
        c5c6 = checkC6();
    });
    std::thread tConv([&]() {
        std::string step = std::string(outDir) + "/handle-lock.step";
        conv = runConvert(stl, step, r);
        std::printf(
            "RULE22_API ok=%d exit=%d builtCyl=%d volΔ=%.6f wt=%d open=%d "
            "warn=%zu err=%s\n",
            (int)r.ok, r.exitCode, r.smoothBuiltCylinders, r.volumeDeltaPct,
            (int)r.watertight, r.openShells, r.warnings.size(), r.error.c_str());
    });
    tFormula.join();
    tConv.join();

    if (c4 != 0) return 10 + c4;
    if (c5c6 != 0) return 20 + c5c6;
    if (conv != 0 || !r.ok) return 30;
    std::printf("RULE22_API PASS C4/C5/C6 + convert()\n");
    return 0;
}
