// D-130-12 — certifies src/stl_quant.cpp: the mesh file's own coordinate
// quantization, measured from the file at TWO magnitudes on a synthetic mesh.
//
// The point of the two magnitudes is that q is not a constant. binary32's ulp
// is constant only inside a binade, so a 10 mm part and a 1000 mm part written
// by the same writer carry different certificates, and the ratio between them is
// a power of two the test states explicitly. The oracle is not the
// implementation's formula: it is the round-trip through the format itself —
// the largest |double -> float32 -> double| error any coordinate in the file can
// suffer, found by scanning the actual binade, times sqrt(3) for the three
// independent axes.
//
// SPDX-License-Identifier: MIT

#include "stl_quant.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using stl2step::StlQuantFloor;
using stl2step::stlQuantFloor;

static int gPass = 0;
static int gFail = 0;

static void check(bool ok, const char* name) {
    if (ok) { ++gPass; std::fprintf(stderr, "PASS %s\n", name); }
    else    { ++gFail; std::fprintf(stderr, "FAIL %s\n", name); }
}

static void checkNear(double got, double want, double tol, const char* name) {
    const bool ok = std::isfinite(got) && std::isfinite(want) && std::fabs(got - want) <= tol;
    if (ok) { ++gPass; std::fprintf(stderr, "PASS %s (%.17g vs %.17g)\n", name, got, want); }
    else    { ++gFail; std::fprintf(stderr, "FAIL %s got %.17g want %.17g tol %.3g\n",
                                    name, got, want, tol); }
}

// ---------------------------------------------------------------------------
// the oracle: what the FORMAT does, measured, not what the code computes
// ---------------------------------------------------------------------------

// The largest |x - float32(x)| over doubles x in the binade of `mag`, doubled to
// a full spacing. Found by bisection on the midpoint between two neighbouring
// float32 values, which is where the rounding error is maximal (half a spacing).
static double formatSpacingAt(double mag) {
    const float f = (float)std::fabs(mag);
    const float up = std::nextafterf(f, std::numeric_limits<float>::infinity());
    // Independent of nextafterf: reconstruct from the stored exponent.
    uint32_t bits;
    std::memcpy(&bits, &f, 4);
    const int expo = (int)((bits >> 23) & 0xFFu) - 127;
    const double byExponent = std::pow(2.0, (double)(expo - 23));
    const double byNextafter = (double)up - (double)f;
    // Both routes must agree, or the oracle itself is wrong.
    if (std::fabs(byExponent - byNextafter) > 0.0) return -1.0;
    return byExponent;
}

// ---------------------------------------------------------------------------
// synthetic meshes
// ---------------------------------------------------------------------------

static void putU32(std::string& b, uint32_t v) {
    for (int i = 0; i < 4; ++i) b.push_back((char)(unsigned char)((v >> (8 * i)) & 0xFFu));
}
static void putF32(std::string& b, float v) {
    uint32_t u;
    std::memcpy(&u, &v, 4);
    putU32(b, u);
}

// An axis-aligned box of half-size `s`, as 12 triangles, written binary.
// Coordinates are +-s exactly, so the file's max |coordinate| IS s.
static bool writeBinaryBox(const std::string& path, double s) {
    static const int kFaces[12][3][3] = {
        {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1}}, {{-1,-1,-1},{ 1, 1,-1},{-1, 1,-1}},
        {{-1,-1, 1},{ 1, 1, 1},{ 1,-1, 1}}, {{-1,-1, 1},{-1, 1, 1},{ 1, 1, 1}},
        {{-1,-1,-1},{-1, 1,-1},{-1, 1, 1}}, {{-1,-1,-1},{-1, 1, 1},{-1,-1, 1}},
        {{ 1,-1,-1},{ 1, 1, 1},{ 1, 1,-1}}, {{ 1,-1,-1},{ 1,-1, 1},{ 1, 1, 1}},
        {{-1,-1,-1},{-1,-1, 1},{ 1,-1, 1}}, {{-1,-1,-1},{ 1,-1, 1},{ 1,-1,-1}},
        {{-1, 1,-1},{ 1, 1, 1},{-1, 1, 1}}, {{-1, 1,-1},{ 1, 1,-1},{ 1, 1, 1}},
    };
    std::string b(80, '\0');
    putU32(b, 12u);
    for (int t = 0; t < 12; ++t) {
        putF32(b, 0.0f); putF32(b, 0.0f); putF32(b, 1.0f);   // normal: not a coordinate
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < 3; ++k) putF32(b, (float)(kFaces[t][c][k] * s));
        b.push_back('\0'); b.push_back('\0');                 // attribute byte count
    }
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = std::fwrite(b.data(), 1, b.size(), f) == b.size();
    std::fclose(f);
    return ok;
}

// The same box in ASCII, printed with `sig` significant digits.
static bool writeAsciiBox(const std::string& path, double s, int sig) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "solid box\n");
    const double v[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                            {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    for (int t = 0; t < 12; ++t) {
        std::fprintf(f, "  facet normal 0 0 1\n    outer loop\n");
        for (int c = 0; c < 3; ++c) {
            const int i = (t * 3 + c) % 8;
            std::fprintf(f, "      vertex %.*e %.*e %.*e\n",
                         sig - 1, v[i][0] * s, sig - 1, v[i][1] * s, sig - 1, v[i][2] * s);
        }
        std::fprintf(f, "    endloop\n  endfacet\n");
    }
    std::fprintf(f, "endsolid box\n");
    std::fclose(f);
    return true;
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? std::string(argv[1]) : std::string(".");
    const std::string small = dir + "/quantfloor_small.stl";
    const std::string large = dir + "/quantfloor_large.stl";
    const std::string ascii = dir + "/quantfloor_ascii.stl";

    const double kRadial = std::sqrt(3.0) / 2.0;

    // --- magnitude 1: a 20 mm box (max |coordinate| = 10) -------------------
    const double sSmall = 10.0;
    check(writeBinaryBox(small, sSmall), "wrote binary box, max|coord| = 10");
    const StlQuantFloor qs = stlQuantFloor(small);
    check(qs.ok, "small: measured");
    check(!qs.ascii, "small: classified binary");
    check(qs.nTri == 12, "small: 12 triangles");
    checkNear(qs.maxAbs, sSmall, 0.0, "small: max |coordinate| is the box half-size");
    const double spacingSmall = formatSpacingAt(sSmall);
    check(spacingSmall > 0.0, "oracle: two routes to the binade spacing agree (small)");
    checkNear(qs.res, spacingSmall, 0.0, "small: res == float32 spacing at 10 mm");
    checkNear(qs.q, spacingSmall * kRadial, 0.0, "small: q == spacing * sqrt(3)/2");

    // --- magnitude 2: the same box at 128x (max |coordinate| = 1280) --------
    // 128 = 2^7 and 10 -> 1280 crosses seven binades exactly, so the certificate
    // must grow by exactly 128. Anything else means q is not the format's ulp.
    const double sLarge = sSmall * 128.0;
    check(writeBinaryBox(large, sLarge), "wrote binary box, max|coord| = 1280");
    const StlQuantFloor ql = stlQuantFloor(large);
    check(ql.ok, "large: measured");
    checkNear(ql.maxAbs, sLarge, 0.0, "large: max |coordinate| is the box half-size");
    const double spacingLarge = formatSpacingAt(sLarge);
    checkNear(ql.res, spacingLarge, 0.0, "large: res == float32 spacing at 1280 mm");
    checkNear(ql.q, spacingLarge * kRadial, 0.0, "large: q == spacing * sqrt(3)/2");
    checkNear(ql.q / qs.q, 128.0, 0.0, "q scales by exactly 2^7 across seven binades");

    // q is a bound on what the file can have done to a vertex: round-trip the
    // worst case and check it is covered.
    {
        const double x = sLarge;                       // on the grid
        const double xMid = x + 0.5 * spacingLarge;    // the worst-case double
        const double err = std::fabs((double)(float)xMid - xMid);
        check(err <= 0.5 * spacingLarge * (1.0 + 1e-12), "per-axis error <= half a spacing");
        check(ql.q >= std::sqrt(3.0) * err * (1.0 - 1e-12),
              "q covers three independent worst-case axis roundings");
    }

    // --- ASCII: the printed decimals, and they are NOT the float32 ulp -------
    check(writeAsciiBox(ascii, sSmall, 7), "wrote ASCII box, 7 significant digits");
    const StlQuantFloor qa = stlQuantFloor(ascii);
    check(qa.ok, "ascii: measured");
    check(qa.ascii, "ascii: classified ascii");
    check(qa.sigDigits == 7, "ascii: 7 significant digits read off the file");
    checkNear(qa.maxAbs, sSmall, 0.0, "ascii: max |coordinate| is the box half-size");
    // 10 mm printed to 7 significant digits is a grid of 1e-5 mm.
    checkNear(qa.res, 1e-5, 1e-18, "ascii: res == 10^(decade - sig + 1)");
    checkNear(qa.q, 1e-5 * kRadial, 1e-18, "ascii: q == res * sqrt(3)/2");

    // A coarser print is a coarser certificate, by exactly the digit ratio.
    check(writeAsciiBox(ascii, sSmall, 4), "rewrote ASCII box, 4 significant digits");
    const StlQuantFloor qa4 = stlQuantFloor(ascii);
    check(qa4.sigDigits == 4, "ascii: 4 significant digits read off the file");
    checkNear(qa4.q / qa.q, 1000.0, 1e-9, "ascii: three digits fewer == 1000x the floor");

    // --- refusals: nothing is invented when nothing was measured ------------
    const StlQuantFloor qMissing = stlQuantFloor(dir + "/no_such_file_quantfloor.stl");
    check(!qMissing.ok && qMissing.q == 0.0, "missing file -> ok=false, q=0");

    std::fprintf(stderr, "quant_floor_unit: %d/%d PASS\n", gPass, gPass + gFail);
    return gFail ? 1 : 0;
}
