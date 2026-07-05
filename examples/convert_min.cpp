// Minimal example of embedding the stl2step engine in another C++ program.
//
//   c++ -std=c++17 convert_min.cpp -lstl2step ... -o convert_min   (see README)
//   ./convert_min part.stl part.step
//
// SPDX-License-Identifier: MIT

#include <stl2step/stl2step.hpp>

#include <cstdio>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <input.stl> [output.step]\n", argv[0]);
        return 2;
    }

    stl2step::Options opt;
    opt.input  = argv[1];
    if (argc > 2) opt.output = argv[2];
    opt.verify = false;   // the host is about to import the STEP, so skip the re-read

    // Optional: forward engine progress to your own logger.
    auto log = [](stl2step::Severity sev, const std::string& msg) {
        const char* tag = sev == stl2step::Severity::Warning ? "[warn] "
                        : sev == stl2step::Severity::Error   ? "[err ] "
                                                             : "       ";
        std::fprintf(stderr, "%s%s\n", tag, msg.c_str());
    };

    stl2step::Result r = stl2step::convert(opt, log);

    if (!r.ok) {
        std::fprintf(stderr, "conversion failed: %s\n", r.error.c_str());
        return 1;
    }

    std::printf("wrote %s: %d solid(s), %d face(s), %.3f mm^3, %.2fs\n",
                r.output.c_str(), r.solids, r.facesAfterUnify, r.meshVolumeMM3, r.seconds);
    if (!r.warnings.empty())
        std::printf("  (%zu warning(s); e.g. \"%s\")\n",
                    r.warnings.size(), r.warnings.front().c_str());

    // r.toJson() would give you the same machine-readable line the CLI prints.
    return r.exitCode;   // 0 clean, 2 written-with-warnings
}
