// dxf_check — round-trip reader + API exercise for writeProfileDxf.
//
//   dxf_check <files...> [--expect-closed] [--radius-tol 1e-9]
//   dxf_check --self-test <dir>     (constructs profiles, writes, re-reads)
//
// SPDX-License-Identifier: MIT

#include "dxf_read.hpp"
#include "dxf_export.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <gp_Dir.hxx>
#include <gp_Pnt2d.hxx>

using stl2step::refit::PrismLevels;
using stl2step::refit::ProfLoop;
using stl2step::refit::ProfSeg;
using stl2step::refit::Profile;
using stl2step::refit::emitProfilesDxf;
using stl2step::refit::makeProfileDxfName;
using stl2step::refit::writeProfileDxf;

namespace {

static double kPi() { return std::acos(-1.0); }

static ProfSeg makeLine(double x0, double y0, double x1, double y1, bool declined = false) {
    ProfSeg s;
    s.isArc = false;
    s.a = gp_Pnt2d(x0, y0);
    s.b = gp_Pnt2d(x1, y1);
    s.declinedAmbiguous = declined;
    return s;
}

static ProfSeg makeArc(double cx, double cy, double R, double a0, double a1, bool ccw) {
    ProfSeg s;
    s.isArc = true;
    s.center = gp_Pnt2d(cx, cy);
    s.R = R;
    s.ccw = ccw;
    double sweep = a1 - a0;
    if (!ccw) sweep = a0 - a1;
    if (sweep < 0.0) sweep += 2.0 * kPi();
    s.phi = sweep;
    s.a = gp_Pnt2d(cx + R * std::cos(a0), cy + R * std::sin(a0));
    s.b = gp_Pnt2d(cx + R * std::cos(a1), cy + R * std::sin(a1));
    return s;
}

static ProfSeg makeCircle(double cx, double cy, double R) {
    ProfSeg s;
    s.isArc = true;
    s.center = gp_Pnt2d(cx, cy);
    s.R = R;
    s.phi = 2.0 * kPi();
    s.ccw = true;
    s.a = s.b = gp_Pnt2d(cx + R, cy);
    return s;
}

static ProfLoop closedRect(double x0, double y0, double x1, double y1, bool outer) {
    ProfLoop lp;
    lp.outer = outer;
    lp.area = std::fabs((x1 - x0) * (y1 - y0));
    lp.segs.push_back(makeLine(x0, y0, x1, y0));
    lp.segs.push_back(makeLine(x1, y0, x1, y1));
    lp.segs.push_back(makeLine(x1, y1, x0, y1));
    lp.segs.push_back(makeLine(x0, y1, x0, y0));
    return lp;
}

static ProfLoop circleHole(double cx, double cy, double R) {
    ProfLoop lp;
    lp.outer = false;
    lp.area = kPi() * R * R;
    lp.segs.push_back(makeCircle(cx, cy, R));
    return lp;
}

// Closed hole: one ARC plus the chord LINE so the loop closes.
static ProfLoop arcHole(double cx, double cy, double R, double a0, double a1) {
    ProfLoop lp;
    lp.outer = false;
    ProfSeg a = makeArc(cx, cy, R, a0, a1, true);
    lp.segs.push_back(a);
    lp.segs.push_back(makeLine(a.b.X(), a.b.Y(), a.a.X(), a.a.Y()));
    return lp;
}

// Handle-lock-shaped synthetic profiles (test data, not src/ thresholds).
// Levels and radii are the STEP-TRUTH falsifiers this lane measures against.
static void buildHandleLockProfiles(PrismLevels& lv, std::vector<Profile>& profs) {
    lv.axis = gp_Dir(0.0, 1.0, 0.0);
    lv.y = {-288.825758, -280.450758, -277.625758};
    lv.capRegion = {2, 4, 21};
    lv.ok = true;
    lv.failedCond = 0;

    const double tenDeg = 10.0 * kPi() / 180.0;

    Profile p0;
    p0.slab = 0;
    p0.loops.push_back(closedRect(-40, -20, 40, 20, true));
    p0.loops.push_back(circleHole(20, 0, 5.75));          // f12 full360
    p0.loops.push_back(circleHole(-12, -6, 5.0));          // f5 through
    p0.loops.push_back(circleHole(-12, 8, 10.0));          // f13 through
    p0.loops.push_back(arcHole(-30, 0, 0.5, 0.0, kPi()));
    p0.loops.push_back(arcHole(-24, 10, 3.7872, -0.4, 0.4));
    p0.loops.push_back(arcHole(0, 14, 6.0, 0.2, 1.7));
    p0.loops.push_back(arcHole(28, 8, 15.0, 0.0, 1.1214662292));
    p0.loops.push_back(arcHole(8, -12, 20.0, -0.5 * tenDeg, 0.5 * tenDeg));  // f7 one ARC
    p0.loops.push_back(arcHole(32, -8, 30.0, 3.14159265359, 3.60050560803));
    p0.loops.push_back(arcHole(16, 12, 10.0, -0.785398163397, 0.785398163397));
    profs.push_back(p0);

    Profile p1;
    p1.slab = 1;
    p1.loops.push_back(closedRect(-18, -14, 18, 14, true));
    p1.loops.push_back(circleHole(-12, -6, 5.0));          // f5 through (RULE 5.2a)
    p1.loops.push_back(circleHole(-12, 8, 10.0));          // f13 through
    p1.loops.push_back(arcHole(6, 0, 16.0, 0.410611824223, 1.95453174413));  // f8 R16
    p1.loops.push_back(arcHole(10, -4, 9.0, 3.14159265359, 5.09674827659));  // f20 R9
    profs.push_back(p1);
}

static Profile buildDeclinedProfile() {
    Profile p;
    p.slab = 0;
    ProfLoop lp;
    lp.outer = true;
    lp.segs.push_back(makeLine(0, 0, 4, 0, false));
    lp.segs.push_back(makeLine(4, 0, 4, 1, true));
    lp.segs.push_back(makeLine(4, 1, 3, 2, true));
    lp.segs.push_back(makeLine(3, 2, 0, 2, false));
    lp.segs.push_back(makeLine(0, 2, 0, 0, false));
    p.loops.push_back(lp);
    return p;
}

static int countDeclined(const Profile& p) {
    int n = 0;
    for (const auto& lp : p.loops)
        for (const auto& s : lp.segs)
            if (s.declinedAmbiguous) ++n;
    return n;
}

static bool radiusMatch(double a, double b, double relTol) {
    const double scale = std::max({std::fabs(a), std::fabs(b), 1.0});
    return std::fabs(a - b) <= relTol * scale;
}

static bool compareToProfile(const stl2step::dxfcheck::DxfFile& dxf,
                             const Profile& src, double relTol, double closeTol,
                             std::string& err) {
    if (static_cast<int>(dxf.loops.size()) != static_cast<int>(src.loops.size())) {
        err = dxf.path + ": loop count " + std::to_string(dxf.loops.size()) +
              " != " + std::to_string(src.loops.size());
        return false;
    }
    for (size_t i = 0; i < src.loops.size(); ++i) {
        if (dxf.loops[i].segs.size() != src.loops[i].segs.size()) {
            err = dxf.path + ": loop " + std::to_string(i) + " seg count " +
                  std::to_string(dxf.loops[i].segs.size()) + " != " +
                  std::to_string(src.loops[i].segs.size());
            return false;
        }
        if (!dxf.loops[i].closed) {
            err = dxf.path + ": loop " + std::to_string(i) + " not closed";
            return false;
        }
        for (size_t k = 0; k < src.loops[i].segs.size(); ++k) {
            const auto& ss = src.loops[i].segs[k];
            if (!ss.isArc || ss.declinedAmbiguous) continue;
            const auto& de = dxf.loops[i].segs[k];
            if (de.type == stl2step::dxfcheck::EntType::Line) {
                err = dxf.path + ": arc emitted as line";
                return false;
            }
            if (!radiusMatch(de.R, ss.R, relTol)) {
                err = dxf.path + ": radius drift";
                return false;
            }
        }
    }
    (void)closeTol;
    return true;
}

struct CheckResult {
    std::string path;
    bool ok = false;
    std::string err;
    stl2step::dxfcheck::DxfFile dxf;
};

static CheckResult checkOne(const std::string& path, bool expectClosed,
                            double closeTol) {
    CheckResult r;
    r.path = path;
    std::string err;
    if (!stl2step::dxfcheck::readDxf(path, r.dxf, err)) {
        r.err = err;
        return r;
    }
    if (r.dxf.nBanned > 0) {
        r.err = path + ": banned entity present";
        return r;
    }
    stl2step::dxfcheck::reconstructLoops(r.dxf, closeTol);
    if (expectClosed && !stl2step::dxfcheck::loopsClosed(r.dxf, closeTol, err)) {
        r.err = err;
        return r;
    }
    r.ok = true;
    return r;
}

static int runCheck(const std::vector<std::string>& files, bool expectClosed,
                    double radiusTol) {
    const double closeTol = std::max(1e-6, radiusTol);
    std::vector<CheckResult> results(files.size());
    std::vector<std::thread> workers;
    workers.reserve(files.size());
    for (size_t i = 0; i < files.size(); ++i) {
        workers.emplace_back([&, i]() {
            results[i] = checkOne(files[i], expectClosed, closeTol);
        });
    }
    for (auto& t : workers) t.join();

    int nLine = 0, nArc = 0, nCircle = 0, nDec = 0;
    std::set<double> radii;
    bool allOk = true;
    for (const auto& r : results) {
        if (!r.ok) {
            std::fprintf(stderr, "FAIL %s: %s\n", r.path.c_str(), r.err.c_str());
            allOk = false;
            continue;
        }
        nLine += r.dxf.nLine;
        nArc += r.dxf.nArc;
        nCircle += r.dxf.nCircle;
        nDec += r.dxf.nDeclined;
        for (double R : r.dxf.radii) radii.insert(R);
        std::printf("%s: LINE=%d ARC=%d CIRCLE=%d declined=%d loops=%zu closed=%s\n",
                    r.path.c_str(), r.dxf.nLine, r.dxf.nArc, r.dxf.nCircle,
                    r.dxf.nDeclined, r.dxf.loops.size(),
                    expectClosed ? "ok" : "n/a");
    }
    std::printf("TOTAL LINE=%d ARC=%d CIRCLE=%d declined=%d distinctRadii=%zu radiusTol=%g\n",
                nLine, nArc, nCircle, nDec, radii.size(), radiusTol);
    return allOk ? 0 : 1;
}

static int runSelfTest(const std::string& dir) {
    PrismLevels lv;
    std::vector<Profile> profs;
    buildHandleLockProfiles(lv, profs);

    const int n = emitProfilesDxf(profs, lv, dir);
    if (n != 2) {
        std::fprintf(stderr, "self-test: expected 2 files, emit returned %d\n", n);
        return 1;
    }

    const double relTol = 1e-9;
    const double closeTol = 1e-6;
    bool ok = true;
    int nLine = 0, nArc = 0, nCircle = 0;

    for (const auto& p : profs) {
        const std::string path =
            (std::filesystem::path(dir) / makeProfileDxfName(p, lv)).string();
        stl2step::dxfcheck::DxfFile dxf;
        std::string err;
        if (!stl2step::dxfcheck::readDxf(path, dxf, err)) {
            std::fprintf(stderr, "self-test read: %s\n", err.c_str());
            return 1;
        }
        if (dxf.nBanned) {
            std::fprintf(stderr, "self-test: banned entities in %s\n", path.c_str());
            return 1;
        }
        stl2step::dxfcheck::reconstructLoops(dxf, closeTol);
        if (!compareToProfile(dxf, p, relTol, closeTol, err)) {
            std::fprintf(stderr, "self-test round-trip: %s\n", err.c_str());
            ok = false;
        }
        nLine += dxf.nLine;
        nArc += dxf.nArc;
        nCircle += dxf.nCircle;
        std::printf("%s: LINE=%d ARC=%d CIRCLE=%d declined=%d loops=%zu\n",
                    path.c_str(), dxf.nLine, dxf.nArc, dxf.nCircle,
                    dxf.nDeclined, dxf.loops.size());
    }

    // Determinism: write again to a sibling dir and cmp via byte read.
    const std::string dir2 = dir + "_b";
    if (emitProfilesDxf(profs, lv, dir2) != 2) {
        std::fprintf(stderr, "self-test: second emit failed\n");
        return 1;
    }
    for (const auto& p : profs) {
        const std::string a = (std::filesystem::path(dir) / makeProfileDxfName(p, lv)).string();
        const std::string b = (std::filesystem::path(dir2) / makeProfileDxfName(p, lv)).string();
        std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
        const std::string sa((std::istreambuf_iterator<char>(fa)), {});
        const std::string sb((std::istreambuf_iterator<char>(fb)), {});
        if (sa != sb) {
            std::fprintf(stderr, "self-test: byte drift %s\n", a.c_str());
            ok = false;
        }
    }

    // Declined layer
    Profile dec = buildDeclinedProfile();
    const std::string decPath = (std::filesystem::path(dir) / "declined.dxf").string();
    if (!writeProfileDxf(dec, lv, decPath)) {
        std::fprintf(stderr, "self-test: declined write failed\n");
        return 1;
    }
    stl2step::dxfcheck::DxfFile ddxf;
    std::string err;
    if (!stl2step::dxfcheck::readDxf(decPath, ddxf, err) || ddxf.nDeclined != 2) {
        std::fprintf(stderr, "self-test: declined count %d (want 2)\n", ddxf.nDeclined);
        ok = false;
    }
    std::printf("declined.dxf: DECLINED layer ents=%d header=%d src=%d\n",
                ddxf.nDeclined, ddxf.declinedHeader, countDeclined(dec));

    // f12 must be CIRCLE
    bool sawCircle575 = false;
    for (const auto& p : profs) {
        const std::string path =
            (std::filesystem::path(dir) / makeProfileDxfName(p, lv)).string();
        stl2step::dxfcheck::DxfFile dxf;
        if (!stl2step::dxfcheck::readDxf(path, dxf, err)) continue;
        for (const auto& e : dxf.ents) {
            if (e.type == stl2step::dxfcheck::EntType::Circle &&
                radiusMatch(e.R, 5.75, relTol))
                sawCircle575 = true;
        }
    }
    if (!sawCircle575) {
        std::fprintf(stderr, "self-test: f12 R=5.75 was not a CIRCLE\n");
        ok = false;
    }

    std::printf("SELF-TEST files=%d LINE=%d ARC=%d CIRCLE=%d handle-lock-declined=0\n",
                n, nLine, nArc, nCircle);
    return ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    bool expectClosed = false;
    bool selfTest = false;
    double radiusTol = 1e-9;
    std::string selfDir;
    std::vector<std::string> files;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--expect-closed") expectClosed = true;
        else if (a == "--radius-tol") {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: --radius-tol needs a value\n");
                return 1;
            }
            radiusTol = std::atof(argv[++i]);
        } else if (a == "--self-test") {
            selfTest = true;
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: --self-test needs a directory\n");
                return 1;
            }
            selfDir = argv[++i];
        } else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
            return 1;
        } else {
            files.push_back(a);
        }
    }

    if (selfTest) return runSelfTest(selfDir);
    if (files.empty()) {
        std::fprintf(stderr, "usage: dxf_check <files...> [--expect-closed] [--radius-tol 1e-9]\n"
                             "       dxf_check --self-test <dir>\n");
        return 1;
    }
    return runCheck(files, expectClosed, radiusTol);
}
