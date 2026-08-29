// Minimal band-anatomy.json reader (cylinder bands + inverse rows).
#pragma once

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct AnatomyBand {
    int face = 0;
    int N = 0;
    int closed360 = 0;
    double radius = 0;
    double extentDeg = 0;
    double axis[3] = {0, 1, 0};
    double origin[3] = {0, 0, 0};
    std::vector<int> triIds;
    std::vector<double> thetaDeg;
    std::vector<double> wMm;
    std::vector<double> dihedralDeg;
    double RfromTheta = 0;
    double RfromDihedral = 0;
};

inline bool readAll(const std::string& path, std::string& out) {
    std::ifstream f(path.c_str());
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return !out.empty();
}

inline const char* skipWs(const char* p) {
    while (p && *p && std::isspace(static_cast<unsigned char>(*p))) ++p;
    return p;
}

inline bool parseIntArray(const char* p, std::vector<int>& v) {
    p = skipWs(p);
    if (!p || *p != '[') return false;
    ++p;
    v.clear();
    while (p && *p) {
        p = skipWs(p);
        if (*p == ']') return true;
        char* end = nullptr;
        const long x = std::strtol(p, &end, 10);
        if (end == p) return false;
        v.push_back(static_cast<int>(x));
        p = skipWs(end);
        if (*p == ',') {
            ++p;
            continue;
        }
        if (*p == ']') return true;
    }
    return false;
}

inline bool parseDblArray(const char* p, std::vector<double>& v) {
    p = skipWs(p);
    if (!p || *p != '[') return false;
    ++p;
    v.clear();
    while (p && *p) {
        p = skipWs(p);
        if (*p == ']') return true;
        char* end = nullptr;
        const double x = std::strtod(p, &end);
        if (end == p) return false;
        v.push_back(x);
        p = skipWs(end);
        if (*p == ',') {
            ++p;
            continue;
        }
        if (*p == ']') return true;
    }
    return false;
}

inline bool loadAnatomy(const std::string& path, std::vector<AnatomyBand>& bands) {
    std::string txt;
    if (!readAll(path, txt)) return false;
    bands.clear();
    const char* bandsKey = std::strstr(txt.c_str(), "\"bands\"");
    const char* invKey = std::strstr(txt.c_str(), "\"inverse\"");
    if (!bandsKey) return false;
    const char* stop = invKey ? invKey : (txt.c_str() + txt.size());
    const char* p = bandsKey;
    while (p < stop && (p = std::strstr(p, "\"true_face_id\"")) != nullptr && p < stop) {
        const char* end = std::strstr(p + 14, "\"true_face_id\"");
        if (!end || end > stop) end = stop;
        std::string chunk(p, static_cast<size_t>(end - p));
        AnatomyBand b;
        b.face = std::atoi(std::strchr(p + 14, ':') + 1);
        const char* t = std::strstr(chunk.c_str(), "\"type\"");
        if (!t || !std::strstr(t, "cylinder")) {
            p = end;
            continue;
        }
        t = std::strstr(chunk.c_str(), "\"N\"");
        if (t) b.N = std::atoi(std::strchr(t + 3, ':') + 1);
        if (b.N <= 0) {
            p = end;
            continue;
        }
        t = std::strstr(chunk.c_str(), "\"radius\"");
        if (t) b.radius = std::strtod(std::strchr(t + 8, ':') + 1, nullptr);
        t = std::strstr(chunk.c_str(), "\"closed360\"");
        if (t) b.closed360 = std::atoi(std::strchr(t + 11, ':') + 1);
        t = std::strstr(chunk.c_str(), "\"extent_deg\"");
        if (t) b.extentDeg = std::strtod(std::strchr(t + 12, ':') + 1, nullptr);
        t = std::strstr(chunk.c_str(), "\"axis\"");
        if (t) {
            const char* br = std::strchr(t, '[');
            if (br) {
                std::vector<double> ax;
                parseDblArray(br, ax);
                if (ax.size() >= 3) {
                    b.axis[0] = ax[0];
                    b.axis[1] = ax[1];
                    b.axis[2] = ax[2];
                }
            }
        }
        t = std::strstr(chunk.c_str(), "\"origin\"");
        if (t) {
            const char* br = std::strchr(t, '[');
            if (br) {
                std::vector<double> o;
                parseDblArray(br, o);
                if (o.size() >= 3) {
                    b.origin[0] = o[0];
                    b.origin[1] = o[1];
                    b.origin[2] = o[2];
                }
            }
        }
        t = std::strstr(chunk.c_str(), "\"tri_ids\"");
        if (t) parseIntArray(std::strchr(t, '['), b.triIds);
        t = std::strstr(chunk.c_str(), "\"theta_deg\"");
        if (t) parseDblArray(std::strchr(t, '['), b.thetaDeg);
        t = std::strstr(chunk.c_str(), "\"w_mm\"");
        if (t) parseDblArray(std::strchr(t, '['), b.wMm);
        t = std::strstr(chunk.c_str(), "\"dihedral_deg\"");
        if (t) parseDblArray(std::strchr(t, '['), b.dihedralDeg);
        bands.push_back(b);
        p = end;
    }

    // inverse block — R_from_w_theta / R_from_w_dihedral
    const char* inv = std::strstr(txt.c_str(), "\"inverse\"");
    if (inv) {
        for (AnatomyBand& b : bands) {
            char key[64];
            std::snprintf(key, sizeof key, "\"true_face_id\": %d", b.face);
            const char* hit = std::strstr(inv, key);
            if (!hit) {
                std::snprintf(key, sizeof key, "\"true_face_id\":%d", b.face);
                hit = std::strstr(inv, key);
            }
            if (!hit) continue;
            const char* th = std::strstr(hit, "\"R_from_w_theta\"");
            const char* dh = std::strstr(hit, "\"R_from_w_dihedral\"");
            if (th && th < hit + 800)
                b.RfromTheta = std::strtod(std::strchr(th, ':') + 1, nullptr);
            if (dh && dh < hit + 800)
                b.RfromDihedral = std::strtod(std::strchr(dh, ':') + 1, nullptr);
        }
    }
    return !bands.empty();
}

struct FaceLabel {
    int tri = -1;
    int face = -1;
};

inline bool loadTriLabels(const std::string& path, std::vector<FaceLabel>& labs) {
    std::string txt;
    if (!readAll(path, txt)) return false;
    labs.clear();
    const char* p = txt.c_str();
    while ((p = std::strstr(p, "\"tri_id\"")) != nullptr) {
        FaceLabel L;
        L.tri = std::atoi(std::strchr(p + 8, ':') + 1);
        const char* f = std::strstr(p, "\"true_face_id\"");
        if (f) L.face = std::atoi(std::strchr(f + 14, ':') + 1);
        labs.push_back(L);
        p += 8;
    }
    return !labs.empty();
}
