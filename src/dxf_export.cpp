// writeProfileDxf — DXF R12 serializer for a fitted 2D profile.
// Write-only. No read path. No timestamps. No geometry decisions.
//
// SPDX-License-Identifier: MIT

#include "dxf_export.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace stl2step {
namespace refit {
namespace {

static double twoPi() { return 2.0 * std::acos(-1.0); }

// DXF group 50/51 are degrees CCW from +X. Unit conversion only — not a
// prismaticity or fitting threshold (RULE 4.2a).
static double radToDxfDeg(double rad) {
    return rad * (180.0 / std::acos(-1.0));
}

// DXF full-turn period is 360 degrees. Wrap only; not a decision boundary.
static double wrapDxfDeg(double deg) {
    const double turn = 360.0;
    deg = std::fmod(deg, turn);
    if (deg < 0.0) deg += turn;
    return deg;
}

// Numeric-closure outer bound on |phi| vs 2*pi. IEEE remainder on a full
// turn, not a geometric fit gate.
static bool isFullTurn(double phi) {
    const double t = twoPi();
    const double bound = 1e-12 * t;
    const double ap = std::fabs(phi);
    return std::fabs(ap - t) <= bound || ap >= t - bound;
}

// Locale-independent float writer: classic() imbue, no thread-locale mutation.
static std::string fmtClassic(double v, std::ios_base::fmtflags field, int prec,
                              bool showPos = false) {
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
    oss.setf(field, std::ios::floatfield);
    oss << std::setprecision(prec);
    if (showPos) oss << std::showpos;
    oss << v;
    return oss.str();
}

static std::string fmtG16(double v) {
    return fmtClassic(v, std::ios::fmtflags(0), 16);
}

static void appendLine(std::string& s, const char* t) {
    s += t;
    s += '\n';
}

static void appendInt(std::string& s, int v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%d", v);
    appendLine(s, buf);
}

static void appendNum(std::string& s, double v) {
    if (v == 0.0) {
        appendLine(s, "0");
        return;
    }
    const std::string t = fmtG16(v);
    if (t == "-0" || t == "-0.0") {
        appendLine(s, "0");
        return;
    }
    appendLine(s, t.c_str());
}

static void pairS(std::string& s, int code, const char* val) {
    appendInt(s, code);
    appendLine(s, val);
}

static void pairI(std::string& s, int code, int val) {
    appendInt(s, code);
    appendInt(s, val);
}

static void pairD(std::string& s, int code, double val) {
    appendInt(s, code);
    appendNum(s, val);
}

static const char* layerFor(const ProfLoop& loop, const ProfSeg& seg) {
    if (seg.declinedAmbiguous) return "DECLINED";
    if (loop.outer) return "OUTER";
    return "HOLES";
}

static void emitLine(std::string& s, const char* layer,
                     double x0, double y0, double x1, double y1) {
    pairS(s, 0, "LINE");
    pairS(s, 8, layer);
    pairD(s, 10, x0);
    pairD(s, 20, y0);
    pairD(s, 30, 0.0);
    pairD(s, 11, x1);
    pairD(s, 21, y1);
    pairD(s, 31, 0.0);
}

static void emitCircle(std::string& s, const char* layer,
                       double cx, double cy, double R) {
    pairS(s, 0, "CIRCLE");
    pairS(s, 8, layer);
    pairD(s, 10, cx);
    pairD(s, 20, cy);
    pairD(s, 30, 0.0);
    pairD(s, 40, R);
}

static void emitArc(std::string& s, const char* layer, const ProfSeg& seg) {
    const double cx = seg.center.X();
    const double cy = seg.center.Y();
    const double angA = std::atan2(seg.a.Y() - cy, seg.a.X() - cx);
    double start = radToDxfDeg(angA);
    double sweep = radToDxfDeg(std::fabs(seg.phi));
    if (!seg.ccw) sweep = -sweep;
    double end = start + sweep;
    // DXF ARC is always CCW from 50 to 51. A CW segment is the same curve
    // as CCW from b to a, so swap the group-50/51 ends.
    if (!seg.ccw) {
        const double tmp = start;
        start = end;
        end = tmp;
    }
    pairS(s, 0, "ARC");
    pairS(s, 8, layer);
    pairD(s, 10, cx);
    pairD(s, 20, cy);
    pairD(s, 30, 0.0);
    pairD(s, 40, seg.R);
    pairD(s, 50, wrapDxfDeg(start));
    pairD(s, 51, wrapDxfDeg(end));
}

static void emitLayer(std::string& s, const char* name, int color) {
    pairS(s, 0, "LAYER");
    pairS(s, 2, name);
    pairI(s, 70, 0);
    pairI(s, 62, color);
    pairS(s, 6, "CONTINUOUS");
}

static int countDeclined(const Profile& p) {
    int n = 0;
    for (const auto& loop : p.loops) {
        for (const auto& seg : loop.segs) {
            if (seg.declinedAmbiguous) ++n;
        }
    }
    return n;
}

static std::string buildDxf(const Profile& p, const PrismLevels& lv) {
    double y0 = 0.0, y1 = 0.0;
    if (p.slab >= 0 && p.slab + 1 < static_cast<int>(lv.y.size())) {
        y0 = lv.y[static_cast<size_t>(p.slab)];
        y1 = lv.y[static_cast<size_t>(p.slab) + 1];
    }
    const int nDec = countDeclined(p);

    std::string s;
    s.reserve(4096);

    pairS(s, 0, "SECTION");
    pairS(s, 2, "HEADER");
    pairS(s, 9, "$ACADVER");
    pairS(s, 1, "AC1009");
    pairS(s, 9, "$INSUNITS");
    pairI(s, 70, 4);  // millimetres
    pairS(s, 9, "$PROJECTNAME");
    {
        std::string proj = "slab=";
        char ib[16];
        std::snprintf(ib, sizeof(ib), "%d", p.slab);
        proj += ib;
        proj += " y0=";
        proj += fmtG16(y0);
        proj += " y1=";
        proj += fmtG16(y1);
        pairS(s, 1, proj.c_str());
    }
    pairS(s, 0, "ENDSEC");

    pairS(s, 999, "stl2step profile dxf");
    {
        char buf[80];
        std::snprintf(buf, sizeof(buf), "declined=%d", nDec);
        pairS(s, 999, buf);
    }
    {
        std::string axis = "axis=";
        axis += fmtG16(lv.axis.X());
        axis += ',';
        axis += fmtG16(lv.axis.Y());
        axis += ',';
        axis += fmtG16(lv.axis.Z());
        pairS(s, 999, axis.c_str());
    }

    pairS(s, 0, "SECTION");
    pairS(s, 2, "TABLES");
    pairS(s, 0, "TABLE");
    pairS(s, 2, "LAYER");
    pairI(s, 70, 3);
    emitLayer(s, "OUTER", 7);
    emitLayer(s, "HOLES", 1);
    emitLayer(s, "DECLINED", 2);
    pairS(s, 0, "ENDTAB");
    pairS(s, 0, "ENDSEC");

    pairS(s, 0, "SECTION");
    pairS(s, 2, "ENTITIES");
    for (const auto& loop : p.loops) {
        for (const auto& seg : loop.segs) {
            const char* layer = layerFor(loop, seg);
            if (seg.declinedAmbiguous || !seg.isArc) {
                emitLine(s, layer, seg.a.X(), seg.a.Y(), seg.b.X(), seg.b.Y());
                continue;
            }
            if (isFullTurn(seg.phi)) {
                emitCircle(s, layer, seg.center.X(), seg.center.Y(), seg.R);
                continue;
            }
            emitArc(s, layer, seg);
        }
    }
    pairS(s, 0, "ENDSEC");
    pairS(s, 0, "EOF");
    return s;
}

static std::string fmtNameCoord(double v) {
    return fmtClassic(v, std::ios::fixed, 9, true);
}

}  // namespace

bool writeProfileDxf(const Profile& p, const PrismLevels& lv, const std::string& path) {
    if (path.empty()) return false;
    const std::string body = buildDxf(p, lv);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(body.data(), static_cast<std::streamsize>(body.size()));
    return static_cast<bool>(out);
}

std::string makeProfileDxfName(const Profile& p, const PrismLevels& lv) {
    double y0 = 0.0, y1 = 0.0;
    if (p.slab >= 0 && p.slab + 1 < static_cast<int>(lv.y.size())) {
        y0 = lv.y[static_cast<size_t>(p.slab)];
        y1 = lv.y[static_cast<size_t>(p.slab) + 1];
    }
    std::string name = "slab";
    char idx[8];
    std::snprintf(idx, sizeof(idx), "%02d", p.slab);
    name += idx;
    name += "_y";
    name += fmtNameCoord(y0);
    name += "_y";
    name += fmtNameCoord(y1);
    name += ".dxf";
    return name;
}

int emitProfilesDxf(const std::vector<Profile>& profs, const PrismLevels& lv,
                    const std::string& dir) {
    if (dir.empty()) return -1;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) return -1;

    const size_t n = profs.size();
    if (n == 0) return 0;

    std::vector<int> ok(n, 0);
    std::vector<std::thread> workers;
    workers.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        workers.emplace_back([&, i]() {
            const std::string path =
                (std::filesystem::path(dir) / makeProfileDxfName(profs[i], lv)).string();
            ok[i] = writeProfileDxf(profs[i], lv, path) ? 1 : 0;
        });
    }
    for (auto& t : workers) t.join();

    int written = 0;
    for (int v : ok) {
        if (!v) return -1;
        written += v;
    }
    return written;
}

bool consumeEmitDxfFlag(const std::string& flag, const std::string& val,
                        std::string& dirOut) {
    if (flag != "--emit-dxf") return false;
    if (val.empty()) return false;
    dirOut = val;
    return true;
}

}  // namespace refit
}  // namespace stl2step
