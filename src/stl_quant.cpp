// D-130-12 — the mesh's precision is its number format's resolution at its own
// coordinate magnitude.
//
// BINARY.  Every coordinate in the file is an IEEE-754 binary32. A value the CAD
// kernel held in double was rounded to the nearest float32, so each axis carries
// an error of at most half an ulp AT ITS OWN MAGNITUDE, and ulp grows in binades.
// The bound that has to cover every vertex of the mesh is therefore the ulp at
// the largest |coordinate| the file contains; the three axes round independently,
// so the RADIAL displacement of a vertex is bounded by
//     sqrt(3) * ulp/2  ==  ulp * sqrt(3)/2.
// That is q: the radius of the ball a vertex may sit anywhere inside, purely from
// having been written to this file. No surface fitted through those vertices can
// be certified tighter, and nothing about the part is chosen.
//
// ASCII.  The same statement about a decimal grid. The writer emitted P
// significant digits (its format's precision); at magnitude maxAbs the printed
// grid therefore has spacing 10^(floor(log10 maxAbs) - P + 1), and P is read off
// the file as the largest significant-digit count any coordinate token shows.
// Counting per token and taking the MAXIMUM is what makes it robust to
// trailing-zero stripping ("12.5" out of a %g does not mean one decimal of
// precision; some other token in the file spends the format's full width).
//
// SPDX-License-Identifier: MIT

#include "stl_quant.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace stl2step {
namespace {

// sqrt(3)/2 — three independent half-ulp axis roundings, as a radius.
const double kRadialFromAxis = std::sqrt(3.0) / 2.0;

// The spacing of binary32 at |x|: the distance to the next representable value.
// Written with nextafterf so it is the format's own answer, not a reconstructed
// exponent (and so subnormals and 0 come out right).
double ulpFloat32(double x) {
    const float f = static_cast<float>(std::fabs(x));
    if (!std::isfinite(f)) return 0.0;
    const float up = std::nextafterf(f, std::numeric_limits<float>::infinity());
    const double u = static_cast<double>(up) - static_cast<double>(f);
    return u > 0.0 ? u : 0.0;
}

float readF32LE(const unsigned char* p) {
    uint32_t bits = static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) |
                    (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}

// Significant decimal digits of one printed token: from the first non-zero digit
// through the last digit of the mantissa, exponent excluded. "0.00123" -> 3,
// "1.234500e+01" -> 7, "0" / "0.000" -> 0 (carries no precision).
int sigDigitsOfToken(const char* s, const char* end) {
    int seen = 0;
    bool started = false;
    for (const char* c = s; c != end; ++c) {
        if (*c == 'e' || *c == 'E') break;
        if (*c == '+' || *c == '-' || *c == '.') continue;
        if (*c < '0' || *c > '9') break;
        if (!started) {
            if (*c == '0') continue;   // leading zeros are placement, not precision
            started = true;
        }
        ++seen;
    }
    return seen;
}

}  // namespace

StlQuantFloor stlQuantFloor(const std::string& path) {
    StlQuantFloor out;
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return out;

    if (std::fseek(f, 0, SEEK_END) != 0) { std::fclose(f); return out; }
    const long long size = static_cast<long long>(std::ftell(f));
    if (size < 0) { std::fclose(f); return out; }
    std::rewind(f);

    // Format: the binary layout is self-describing (an 80-byte header, a uint32
    // count, then exactly 50 bytes per triangle). If the file size matches that
    // for its own stored count it IS that file; anything else is parsed as text.
    // Same discriminator RWStl applies, and it does not trust the "solid" word,
    // which binary writers also emit.
    bool binary = false;
    long long nTri = 0;
    if (size >= 84) {
        unsigned char head[84];
        if (std::fread(head, 1, 84, f) == 84) {
            const uint32_t n = static_cast<uint32_t>(head[80]) |
                               (static_cast<uint32_t>(head[81]) << 8) |
                               (static_cast<uint32_t>(head[82]) << 16) |
                               (static_cast<uint32_t>(head[83]) << 24);
            if (size == 84LL + 50LL * static_cast<long long>(n)) {
                binary = true;
                nTri = static_cast<long long>(n);
            }
        }
        std::rewind(f);
    }

    double maxAbs = 0.0;
    long long tris = 0;
    int sigDigits = 0;

    if (binary) {
        if (std::fseek(f, 84, SEEK_SET) != 0) { std::fclose(f); return out; }
        std::vector<unsigned char> buf(50);
        for (long long i = 0; i < nTri; ++i) {
            if (std::fread(buf.data(), 1, 50, f) != 50) break;
            // Bytes 0..11 are the facet normal — a direction, not a coordinate.
            for (int k = 0; k < 9; ++k) {
                const double v = std::fabs(static_cast<double>(readF32LE(buf.data() + 12 + 4 * k)));
                if (std::isfinite(v) && v > maxAbs) maxAbs = v;
            }
            ++tris;
        }
    } else {
        std::string line;
        int c;
        auto flushLine = [&]() {
            // "vertex <x> <y> <z>" is the only line carrying coordinates.
            const char* p = line.c_str();
            while (*p == ' ' || *p == '\t' || *p == '\r') ++p;
            if (std::strncmp(p, "vertex", 6) != 0) {
                if (std::strncmp(p, "facet", 5) == 0) ++tris;
                return;
            }
            p += 6;
            for (int k = 0; k < 3 && *p; ++k) {
                while (*p == ' ' || *p == '\t') ++p;
                const char* tok = p;
                while (*p && *p != ' ' && *p != '\t' && *p != '\r') ++p;
                if (tok == p) break;
                const std::string t(tok, static_cast<size_t>(p - tok));
                const double v = std::fabs(std::strtod(t.c_str(), nullptr));
                if (std::isfinite(v) && v > maxAbs) maxAbs = v;
                const int sd = sigDigitsOfToken(tok, p);
                if (sd > sigDigits) sigDigits = sd;
            }
        };
        while ((c = std::fgetc(f)) != EOF) {
            if (c == '\n') { flushLine(); line.clear(); }
            else if (line.size() < 4096) line.push_back(static_cast<char>(c));
        }
        if (!line.empty()) flushLine();
    }
    std::fclose(f);

    out.ascii = !binary;
    out.nTri = tris;
    out.maxAbs = maxAbs;
    out.sigDigits = binary ? 0 : sigDigits;
    if (!(maxAbs > 0.0)) return out;   // ok stays false: nothing was measured

    if (binary) {
        out.res = ulpFloat32(maxAbs);
    } else {
        if (sigDigits <= 0) return out;
        const int decade = static_cast<int>(std::floor(std::log10(maxAbs)));
        out.res = std::pow(10.0, static_cast<double>(decade - sigDigits + 1));
    }
    if (!(out.res > 0.0) || !std::isfinite(out.res)) return out;
    out.q = out.res * kRadialFromAxis;
    out.ok = true;
    return out;
}

}  // namespace stl2step
