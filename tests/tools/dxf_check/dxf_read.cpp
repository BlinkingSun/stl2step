// Minimal DXF R12 ENTITIES reader. Test tool — never linked into the engine.
// SPDX-License-Identifier: MIT

#include "dxf_read.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace stl2step {
namespace dxfcheck {
namespace {

static double kPi() { return std::acos(-1.0); }

static std::string trim(std::string s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' '))
        s.pop_back();
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    return s.substr(i);
}

static double degToRad(double deg) { return deg * (kPi() / 180.0); }

static void endpointsFromArc(DxfEnt& e) {
    const double a0 = degToRad(e.a0deg);
    const double a1 = degToRad(e.a1deg);
    e.x0 = e.cx + e.R * std::cos(a0);
    e.y0 = e.cy + e.R * std::sin(a0);
    e.x1 = e.cx + e.R * std::cos(a1);
    e.y1 = e.cy + e.R * std::sin(a1);
}

static void endpointsFromCircle(DxfEnt& e) {
    e.x0 = e.cx + e.R;
    e.y0 = e.cy;
    e.x1 = e.x0;
    e.y1 = e.y0;
}

static double dist2(double ax, double ay, double bx, double by) {
    const double dx = ax - bx, dy = ay - by;
    return dx * dx + dy * dy;
}

}  // namespace

bool readDxf(const std::string& path, DxfFile& out, std::string& err) {
    out = DxfFile{};
    out.path = path;
    std::ifstream in(path);
    if (!in) {
        err = "cannot open " + path;
        return false;
    }

    std::vector<std::pair<int, std::string>> pairs;
    std::string line;
    while (std::getline(in, line)) {
        const std::string codeS = trim(line);
        if (codeS.empty() && in.eof()) break;
        if (!std::getline(in, line)) break;
        const std::string val = trim(line);
        char* end = nullptr;
        const long code = std::strtol(codeS.c_str(), &end, 10);
        if (end == codeS.c_str()) continue;
        pairs.emplace_back(static_cast<int>(code), val);
    }

    bool inEntities = false;
    DxfEnt cur;
    bool have = false;
    std::string curType;

    auto flush = [&]() {
        if (!have) return;
        if (curType == "LINE") {
            cur.type = EntType::Line;
            out.nLine++;
        } else if (curType == "ARC") {
            cur.type = EntType::Arc;
            endpointsFromArc(cur);
            out.nArc++;
            out.radii.push_back(cur.R);
        } else if (curType == "CIRCLE") {
            cur.type = EntType::Circle;
            endpointsFromCircle(cur);
            out.nCircle++;
            out.radii.push_back(cur.R);
        } else {
            have = false;
            return;
        }
        if (cur.layer == "DECLINED") out.nDeclined++;
        out.ents.push_back(cur);
        have = false;
    };

    for (size_t i = 0; i < pairs.size(); ++i) {
        const int code = pairs[i].first;
        const std::string& val = pairs[i].second;

        if (code == 999) {
            if (val.rfind("declined=", 0) == 0)
                out.declinedHeader = std::atoi(val.c_str() + 9);
            continue;
        }
        if (code == 9 && val == "$PROJECTNAME") {
            if (i + 1 < pairs.size() && pairs[i + 1].first == 1)
                out.projectName = pairs[i + 1].second;
            continue;
        }
        if (code == 0 && val == "SECTION") {
            if (i + 1 < pairs.size() && pairs[i + 1].first == 2)
                inEntities = (pairs[i + 1].second == "ENTITIES");
            continue;
        }
        if (code == 0 && (val == "ENDSEC" || val == "EOF")) {
            if (inEntities) flush();
            inEntities = false;
            continue;
        }
        if (!inEntities) continue;

        if (code == 0) {
            flush();
            cur = DxfEnt{};
            curType = val;
            have = (val == "LINE" || val == "ARC" || val == "CIRCLE");
            if (val == "SPLINE" || val == "LWPOLYLINE" || val == "ELLIPSE" ||
                val == "INSERT" || val == "BLOCK") {
                out.nBanned++;
            }
            continue;
        }
        if (!have) continue;
        const double d = std::atof(val.c_str());
        if (code == 8) cur.layer = val;
        else if (code == 10) {
            if (curType == "LINE") cur.x0 = d;
            else cur.cx = d;
        } else if (code == 20) {
            if (curType == "LINE") cur.y0 = d;
            else cur.cy = d;
        } else if (code == 11) cur.x1 = d;
        else if (code == 21) cur.y1 = d;
        else if (code == 40) cur.R = d;
        else if (code == 50) cur.a0deg = d;
        else if (code == 51) cur.a1deg = d;
    }
    flush();
    return true;
}

void reconstructLoops(DxfFile& f, double closeTol) {
    f.loops.clear();
    const double tol2 = closeTol * closeTol;
    std::vector<char> used(f.ents.size(), 0);

    for (size_t i = 0; i < f.ents.size(); ++i) {
        if (used[i]) continue;
        if (f.ents[i].type == EntType::Circle) {
            DxfLoop lp;
            lp.segs.push_back(f.ents[i]);
            lp.closed = true;
            lp.outer = (f.ents[i].layer == "OUTER");
            f.loops.push_back(lp);
            used[i] = 1;
            continue;
        }

        DxfLoop lp;
        lp.segs.push_back(f.ents[i]);
        lp.outer = (f.ents[i].layer == "OUTER");
        used[i] = 1;
        double hx = f.ents[i].x0, hy = f.ents[i].y0;
        double tx = f.ents[i].x1, ty = f.ents[i].y1;

        bool grew = true;
        while (grew) {
            grew = false;
            for (size_t j = 0; j < f.ents.size(); ++j) {
                if (used[j] || f.ents[j].type == EntType::Circle) continue;
                const DxfEnt& e = f.ents[j];
                if (dist2(tx, ty, e.x0, e.y0) <= tol2) {
                    lp.segs.push_back(e);
                    tx = e.x1;
                    ty = e.y1;
                    used[j] = 1;
                    grew = true;
                    break;
                }
                if (dist2(tx, ty, e.x1, e.y1) <= tol2) {
                    DxfEnt rev = e;
                    std::swap(rev.x0, rev.x1);
                    std::swap(rev.y0, rev.y1);
                    lp.segs.push_back(rev);
                    tx = e.x0;
                    ty = e.y0;
                    used[j] = 1;
                    grew = true;
                    break;
                }
            }
        }
        lp.closed = dist2(hx, hy, tx, ty) <= tol2;
        f.loops.push_back(lp);
    }
}

bool loopsClosed(const DxfFile& f, double closeTol, std::string& err) {
    (void)closeTol;
    for (size_t i = 0; i < f.loops.size(); ++i) {
        if (!f.loops[i].closed) {
            err = f.path + ": loop " + std::to_string(i) + " not closed";
            return false;
        }
    }
    return true;
}

}  // namespace dxfcheck
}  // namespace stl2step
