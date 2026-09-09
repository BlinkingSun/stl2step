// D-I-4 / GitHub #7: Result::toJson must not use two-pass snprintf(nullptr,0)
// + string((size_t)n+1). On MSVC, n == -1 wraps the allocation to 0 and the
// write overflows the 16-byte SSO. The exact bad size is that wrap.
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>

using stl2step::Result;

static int gFail = 0;

static void check(bool cond, const char* what) {
    if (cond) return;
    std::fprintf(stderr, "FAIL %s\n", what);
    ++gFail;
}

int main() {
    // Exact bad index/size of the old writer.
    const int n = -1;
    const std::size_t alloc = static_cast<std::size_t>(n) + 1;
    check(alloc == 0, "(size_t)(-1)+1 wraps to 0 (MSVC SSO overflow size)");
    check(sizeof(std::string) >= 16, "std::string is at least MSVC SSO (16)");

    Result off;
    off.ok = true;
    off.input.assign(17, 'A');   // one past MSVC SSO
    off.output.assign(17, 'B');
    off.watertight = true;
    off.seconds = 1.25;
    off.warnings.emplace_back(std::string(17, 'W'));
    const std::string j = off.toJson();
    check(!j.empty() && j.front() == '{' && j.back() == '}', "off-path toJson is a JSON object");
    check(j.size() == std::strlen(j.c_str()), "toJson has no embedded NUL (old snprintf+resize)");
    check(j.find(std::string(17, 'A')) != std::string::npos, "17-byte input survives (past SSO 16)");
    check(j.find(std::string(17, 'B')) != std::string::npos, "17-byte output survives");
    check(j.find("smoothPlanes") == std::string::npos, "off-path omits smooth* keys");

    Result sm = off;
    sm.facesAfterSmooth = 4;
    sm.smoothPlanes = 1;
    sm.radiusDriftN = 3;
    sm.radiusDriftMaxAbs = 0.001;
    sm.radiusDriftMaxRel = 0.002;
    const std::string js = sm.toJson();
    check(!js.empty() && js.front() == '{' && js.back() == '}', "smooth toJson is a JSON object");
    check(js.size() == std::strlen(js.c_str()), "smooth toJson has no embedded NUL");
    check(js.find("\"smoothPlanes\":1") != std::string::npos, "smoothPlanes present");
    check(js.find("\"radiusDrift\"") != std::string::npos, "radiusDrift present");

    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "result_tojson_unit ok\n");
    return 0;
}
