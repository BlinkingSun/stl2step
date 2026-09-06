// SPEC-volid §4.3 / D-140-10 — volume-attribution identity on S01 / S02 / S03.
// The test is an instrument: it is never edited to pass. Assertions are the
// §2 identities (partition + bracket) and the generator-exact 10 mm cube
// volume on S01; no engine-output goldens.
//
// SPDX-License-Identifier: MIT

#include <cstdio>
#include <cstdlib>

#include <cmath>
#include <cstring>
#include <string>

namespace {

#ifdef _WIN32
std::string quoteArg(const std::string& p) {
    // cmd.exe quoting: wrap in double quotes; double any embedded quotes.
    std::string o = "\"";
    for (char c : p) {
        if (c == '"') o += "\"\"";
        else o += c;
    }
    o += "\"";
    return o;
}
#else
std::string quoteArg(const std::string& p) {
    std::string o = "'";
    for (char c : p) {
        if (c == '\'') o += "'\\''";
        else o += c;
    }
    o += "'";
    return o;
}
#endif

bool setVolidEnv() {
#ifdef _WIN32
    return _putenv_s("STL2STEP_VOLID", "1") == 0;
#else
    return setenv("STL2STEP_VOLID", "1", 1) == 0;
#endif
}

const char* tempDir() {
#ifdef _WIN32
    const char* tmp = std::getenv("TEMP");
    if (!tmp || !tmp[0]) tmp = std::getenv("TMP");
    if (!tmp || !tmp[0]) tmp = ".";
    return tmp;
#else
    const char* tmp = std::getenv("TMPDIR");
    if (!tmp || !tmp[0]) tmp = "/tmp";
    return tmp;
#endif
}

FILE* launch(const char* cmd) {
#ifdef _WIN32
    return _popen(cmd, "rb");
#else
    return popen(cmd, "r");
#endif
}

int finish(FILE* fp) {
#ifdef _WIN32
    return _pclose(fp);
#else
    return pclose(fp);
#endif
}

bool getD(const std::string& line, const char* key, double& out) {
    const std::string pat = std::string(key) + "=";
    const auto pos = line.find(pat);
    if (pos == std::string::npos) return false;
    return std::sscanf(line.c_str() + pos + pat.size(), "%lf", &out) == 1;
}

bool getI(const std::string& line, const char* key, int& out) {
    const std::string pat = std::string(key) + "=";
    const auto pos = line.find(pat);
    if (pos == std::string::npos) return false;
    return std::sscanf(line.c_str() + pos + pat.size(), "%d", &out) == 1;
}

int gPass = 0;
int gFail = 0;

void check(bool ok, const char* name) {
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

int runOne(const char* bin, const char* stl, const char* tag) {
    const std::string out = std::string(tempDir()) + "/volid-" + tag + ".step";
    // No env-var / shell prefix. 2>&1 is valid in cmd.exe and POSIX sh so
    // DIAG_VOLID_SUM on stderr is read through the pipe. On Windows, wrap the
    // whole command in an extra pair of quotes so cmd.exe /c (used by _popen)
    // does not strip the argv0 quotes (KB 156212).
    std::string cmd = quoteArg(bin) + " --engine trueform --threads 1 --no-verify --quiet " +
                       quoteArg(stl) + " -o " + quoteArg(out) + " 2>&1";
#ifdef _WIN32
    cmd = "\"" + cmd + "\"";
#endif
    if (!setVolidEnv()) {
        std::fprintf(stderr, "FAIL %s: setenv STL2STEP_VOLID\n", tag);
        ++gFail;
        return 1;
    }
    FILE* fp = launch(cmd.c_str());
    if (!fp) {
        std::fprintf(stderr, "FAIL %s: popen\n", tag);
        ++gFail;
        return 1;
    }
    std::string ship, last;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), fp)) {
        if (std::strncmp(buf, "DIAG_VOLID_SUM ", 15) == 0) {
            last = buf;
            if (std::strstr(buf, "site=ship")) ship = buf;
        }
    }
    const int rc = finish(fp);
    const std::string& line = !ship.empty() ? ship : last;
    if (line.empty()) {
        std::fprintf(stderr, "FAIL %s: no DIAG_VOLID_SUM (pclose=%d)\n", tag, rc);
        ++gFail;
        return 1;
    }
    std::fprintf(stderr, "%s", line.c_str());

    int uncovered = -1, dbl = -1, nNull = -1, nUnsup = -1;
    int closesA = 0, closesB = 0;
    double chordResid = 0, ulpSum = 0, residPcurve = 0, quadSup = 0, shellVol = 0;
    const bool parsed =
        getI(line, "uncoveredTri", uncovered) && getI(line, "doubleTri", dbl) &&
        getI(line, "nNull", nNull) && getI(line, "nUnsup", nUnsup) &&
        getD(line, "chordResid", chordResid) && getD(line, "ulpSum", ulpSum) &&
        getI(line, "closesA", closesA) && getI(line, "closesB", closesB) &&
        getD(line, "residPcurve", residPcurve) && getD(line, "quadSup", quadSup) &&
        getD(line, "shellVol", shellVol);
    check(parsed, (std::string(tag) + " parsed").c_str());
    if (!parsed) return 1;

    const std::string pfx = std::string(tag) + " ";
    check(uncovered == 0 && dbl == 0 && nNull == 0 && nUnsup == 0,
          (pfx + "uncoveredTri=doubleTri=nNull=nUnsup=0").c_str());
    check(std::fabs(chordResid) <= ulpSum, (pfx + "|chordResid|<=ulpSum").c_str());
    check(closesA == 1, (pfx + "closesA=1").c_str());
    check(closesB == 1, (pfx + "closesB=1").c_str());
    if (std::strcmp(tag, "S01") == 0) {
        // Generator-exact 10 mm cube; S01.exact.step is the same oracle.
        check(std::fabs(shellVol - 1000.0) <= ulpSum, (pfx + "shellVol=1000 cube").c_str());
    }
    if (std::strcmp(tag, "S02") == 0) {
        check(std::fabs(residPcurve) <= quadSup, (pfx + "|residPcurve|<=quadSup").c_str());
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: volid_identity_test <stl2step> <S01.stl> <S02.stl> <S03.stl>\n");
        return 2;
    }
    runOne(argv[1], argv[2], "S01");
    runOne(argv[1], argv[3], "S02");
    runOne(argv[1], argv[4], "S03");
    std::fprintf(stderr, "volid_identity %d pass %d fail\n", gPass, gFail);
    return gFail ? 1 : 0;
}
