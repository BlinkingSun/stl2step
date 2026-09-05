// 140-UNION-CENSUS (SPEC-union PART B, D-140-8 U-R5..U-R8, U-R13) -- certifies
// `claimChordSagitta` and `unionCensus` in src/refit_union_census.cpp against
// independent ground truth on synthetic rings whose MeshView is built here,
// by hand, with the engine's own welded-edge convention (stl2step.cpp:447).
//
// Two oracles for sigma, neither sharing a line with the implementation:
//   * the closed form of a regular N-gon ring, R (1 - cos(pi / N)) -- every
//     circumferential chord and every quad diagonal spans 2 pi / N;
//   * a numerical supremum along every welded edge: R - dist(axis, P(s)) is
//     sampled on a dense grid of s in [0, 1] and the best bracket is refined
//     by golden-section search to machine precision. It never assumes where
//     the maximum sits.
// Three radii x three sampling densities, in the world frame and under a
// rigid placement; assertions to Precision::Confusion() (SPEC-union B.2).
//
// The topology counts (edgePieces / punctures / domainFaces / pinchVertices)
// are asserted on rings whose interruptions are constructed and whose expected
// counts are derived by hand in the comments (U-R6, U-R7, U-R8).
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

#include <Precision.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

using stl2step::refit::MeshView;
using stl2step::refit::UnionCensus;
using stl2step::refit::claimChordSagitta;
using stl2step::refit::unionCensus;

static const double kPi = 3.14159265358979323846264338327950288;

static int gPass = 0;
static int gFail = 0;

static void check(bool ok, const char* name) {
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

// ---------------------------------------------------------------------------
// A hand-built MeshView: one component, local == global, welded edges keyed on
// (lo, hi) exactly as stl2step.cpp builds CompStat (bit s of triDirs set when
// side s runs global-lo -> hi). Owns its storage; the view points into it.
struct SynthMesh {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;  // 3 per triangle
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;

    void addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
    }
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; t++) {
            compTris[t] = static_cast<int>(t);
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2], triFlat[3 * t]};
            for (int s = 0; s < 3; s++) {
                const int lo = std::min(v[s], v[s + 1]), hi = std::max(v[s], v[s + 1]);
                auto it = ids.find({lo, hi});
                int id;
                if (it == ids.end()) {
                    id = static_cast<int>(edges.size());
                    ids.emplace(std::make_pair(lo, hi), id);
                    edges.emplace_back(lo, hi);
                } else {
                    id = it->second;
                }
                triEdges[t][static_cast<std::size_t>(s)] = id;
                if (v[s] < v[s + 1]) triDirs[t] |= static_cast<uint8_t>(1 << s);
            }
        }
        compVtx.resize(pts.size());
        for (std::size_t i = 0; i < pts.size(); i++) compVtx[i] = static_cast<int>(i);
    }
    MeshView view() const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int(*)[3]>(triFlat.data());
        mv.compTris = compTris.data();
        mv.compVtx = compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = triEdges.data();
        mv.triDirs = triDirs.data();
        mv.nTri = triFlat.size() / 3;
        mv.nVtx = pts.size();
        mv.nEdge = edges.size();
        mv.diag = 1.0;
        mv.weldTol = 0.0;
        mv.sewTol = 0.0;
        mv.quantFloor = 0.0;
        return mv;
    }
};

// A frame for placing the ring: origin O, orthonormal X, Y, Z (Z = axis).
struct Frame {
    gp_XYZ O, X, Y, Z;
    gp_XYZ at(double rho, double u, double z) const {
        return O + X * (rho * std::cos(u)) + Y * (rho * std::sin(u)) + Z * z;
    }
    gp_Ax1 axis() const { return gp_Ax1(gp_Pnt(O.X(), O.Y(), O.Z()), gp_Dir(Z.X(), Z.Y(), Z.Z())); }
};

static Frame worldFrame() { return Frame{gp_XYZ(0, 0, 0), gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0), gp_XYZ(0, 0, 1)}; }

static Frame placedFrame() {
    gp_XYZ z(1.0, 2.0, 3.0);
    z.Normalize();
    gp_XYZ x = gp_XYZ(0.3, -0.7, 0.2).Crossed(z);
    x.Normalize();
    gp_XYZ y = z.Crossed(x);
    return Frame{gp_XYZ(5.0, -7.0, 3.0), x, y, z};
}

// Two-row ring (one strip of N quads) on radius R, height H; vertex i of the
// bottom row is b_i = i, of the top row t_i = N + i; azimuth u_i = 2 pi i / N
// (+ jitter_i when given). Quad i is split T1_i = (b_i, b_i+1, t_i+1) and
// T2_i = (b_i, t_i+1, t_i), so triangle ids are 2i and 2i+1.
static SynthMesh makeRing(const Frame& f, int N, double R, double H,
                          const std::vector<double>* jitter = nullptr,
                          const std::vector<int>* dropTris = nullptr) {
    SynthMesh m;
    for (int row = 0; row < 2; row++) {
        for (int i = 0; i < N; i++) {
            const double u = 2.0 * kPi * i / N + (jitter ? (*jitter)[static_cast<std::size_t>(i)] : 0.0);
            m.pts.push_back(f.at(R, u, row ? H : 0.0));
        }
    }
    for (int i = 0; i < N; i++) {
        const int b0 = i, b1 = (i + 1) % N, t0 = N + i, t1 = N + (i + 1) % N;
        const int id1 = 2 * i, id2 = 2 * i + 1;
        const bool d1 = dropTris && std::find(dropTris->begin(), dropTris->end(), id1) != dropTris->end();
        const bool d2 = dropTris && std::find(dropTris->begin(), dropTris->end(), id2) != dropTris->end();
        // a dropped triangle is still emitted (the mesh keeps it) -- it is simply
        // not part of the CLAIM; the claim is built by the caller from ids.
        (void)d1;
        (void)d2;
        m.addTri(b0, b1, t1);
        m.addTri(b0, t1, t0);
    }
    m.finish();
    return m;
}

static std::vector<int> allTris(const SynthMesh& m) {
    std::vector<int> v(m.triFlat.size() / 3);
    for (std::size_t i = 0; i < v.size(); i++) v[i] = static_cast<int>(i);
    return v;
}

static std::vector<int> allBut(const SynthMesh& m, const std::vector<int>& drop) {
    std::vector<int> v;
    for (int t : allTris(m))
        if (std::find(drop.begin(), drop.end(), t) == drop.end()) v.push_back(t);
    return v;
}

// ---------------------------------------------------------------------------
// Oracle 2: numerical supremum of R - dist(axis, P(s)) along every welded edge
// of the claim. Dense grid, then golden-section refinement on the best bracket.
static double axisDistO(const gp_Ax1& ax, const gp_XYZ& p) {
    const gp_XYZ o = ax.Location().XYZ();
    const gp_XYZ a = ax.Direction().XYZ();
    const gp_XYZ d = p - o;
    return (d - a * d.Dot(a)).Modulus();
}

static double numericChordSup(const SynthMesh& m, const std::vector<int>& claim, const gp_Ax1& ax,
                              double R) {
    std::vector<char> seen(m.edges.size(), 0);
    double best = -1e300;
    bool any = false;
    const int kGrid = 2001;
    for (int t : claim) {
        for (int s = 0; s < 3; s++) {
            const int e = m.triEdges[static_cast<std::size_t>(t)][static_cast<std::size_t>(s)];
            if (seen[static_cast<std::size_t>(e)]) continue;
            seen[static_cast<std::size_t>(e)] = 1;
            const gp_XYZ A = m.pts[static_cast<std::size_t>(m.edges[static_cast<std::size_t>(e)].first)];
            const gp_XYZ B = m.pts[static_cast<std::size_t>(m.edges[static_cast<std::size_t>(e)].second)];
            auto f = [&](double u) { return R - axisDistO(ax, A * (1.0 - u) + B * u); };
            double bi = 0, bv = f(0.0);
            for (int k = 1; k < kGrid; k++) {
                const double u = static_cast<double>(k) / (kGrid - 1);
                const double v = f(u);
                if (v > bv) {
                    bv = v;
                    bi = u;
                }
            }
            double lo = std::max(0.0, bi - 1.0 / (kGrid - 1));
            double hi = std::min(1.0, bi + 1.0 / (kGrid - 1));
            const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
            double c = hi - gr * (hi - lo), d = lo + gr * (hi - lo);
            double fc = f(c), fd = f(d);
            for (int it = 0; it < 200; it++) {
                if (fc > fd) {
                    hi = d;
                    d = c;
                    fd = fc;
                    c = hi - gr * (hi - lo);
                    fc = f(c);
                } else {
                    lo = c;
                    c = d;
                    fc = fd;
                    d = lo + gr * (hi - lo);
                    fd = f(d);
                }
            }
            const double v = std::max(std::max(fc, fd), bv);
            if (!any || v > best) {
                best = v;
                any = true;
            }
        }
    }
    return any ? best : 0.0;
}

int main() {
    const double eps = Precision::Confusion();
    char name[256];

    // -- sigma: three radii x three densities, world frame and placed frame ----
    const double radii[3] = {2.5, 10.0, 100.0};
    const int dens[3] = {12, 48, 180};
    double worstClosed = 0.0, worstNum = 0.0, worstPlace = 0.0;
    for (double R : radii) {
        for (int N : dens) {
            const double closed = R * (1.0 - std::cos(kPi / N));
            const Frame fw = worldFrame();
            const SynthMesh mw = makeRing(fw, N, R, 0.5 * R);
            const std::vector<int> cw = allTris(mw);
            const double sw = claimChordSagitta(mw.view(), cw, fw.axis(), R);
            const double nw = numericChordSup(mw, cw, fw.axis(), R);
            worstClosed = std::max(worstClosed, std::abs(sw - closed));
            worstNum = std::max(worstNum, std::abs(sw - nw));
            std::snprintf(name, sizeof name, "sigma R=%g N=%d world: sigma=%.12g closed=%.12g numeric=%.12g",
                          R, N, sw, closed, nw);
            check(std::abs(sw - closed) <= eps && std::abs(sw - nw) <= eps, name);

            const Frame fp = placedFrame();
            const SynthMesh mp = makeRing(fp, N, R, 0.5 * R);
            const std::vector<int> cp = allTris(mp);
            const double sp = claimChordSagitta(mp.view(), cp, fp.axis(), R);
            const double np = numericChordSup(mp, cp, fp.axis(), R);
            worstPlace = std::max(worstPlace, std::abs(sp - sw));
            std::snprintf(name, sizeof name, "sigma R=%g N=%d placed: sigma=%.12g numeric=%.12g |placed-world|=%.3g",
                          R, N, sp, np, std::abs(sp - sw));
            check(std::abs(sp - closed) <= eps && std::abs(sp - np) <= eps && std::abs(sp - sw) <= eps,
                  name);
        }
    }
    std::fprintf(stderr, "sigma worst |closed| %.3g  worst |numeric| %.3g  worst |placement| %.3g  (eps %.3g)\n",
                 worstClosed, worstNum, worstPlace, eps);

    // -- sigma on a jittered (non-regular) ring: only the numeric oracle applies
    {
        const int N = 36;
        const double R = 10.0;
        std::vector<double> jit(static_cast<std::size_t>(N));
        uint32_t st = 0x9E3779B9u;
        for (int i = 0; i < N; i++) {
            st = st * 1664525u + 1013904223u;
            jit[static_cast<std::size_t>(i)] = ((st >> 8) / 16777216.0 - 0.5) * (0.6 * 2.0 * kPi / N);
        }
        const Frame fp = placedFrame();
        const SynthMesh m = makeRing(fp, N, R, 4.0, &jit);
        const std::vector<int> c = allTris(m);
        const double s = claimChordSagitta(m.view(), c, fp.axis(), R);
        const double n = numericChordSup(m, c, fp.axis(), R);
        const double closedRegular = R * (1.0 - std::cos(kPi / N));
        std::snprintf(name, sizeof name, "sigma jittered N=36 R=10: sigma=%.12g numeric=%.12g (regular would be %.12g)",
                      s, n, closedRegular);
        check(std::abs(s - n) <= eps && s > closedRegular, name);
        // a partial claim (half the ring) reads the supremum of ITS edges only
        std::vector<int> half;
        for (int t = 0; t < N; t++) half.push_back(t);
        const double sh = claimChordSagitta(m.view(), half, fp.axis(), R);
        const double nh = numericChordSup(m, half, fp.axis(), R);
        std::snprintf(name, sizeof name, "sigma half claim: sigma=%.12g numeric=%.12g <= whole %.12g", sh, nh, s);
        check(std::abs(sh - nh) <= eps && sh <= s + eps, name);
    }

    // -- sigma of a pure-generator claim is exactly 0 (axial edges only) ------
    {
        // two axial edges and nothing circumferential: a degenerate 'ring' of
        // one quad whose circumferential chords are zero-length is not a mesh;
        // instead assert the contribution of generator edges via a ring whose
        // claim is empty -> 0, and via the identity sigma(claim) >= 0 on rings.
        const Frame fw = worldFrame();
        const SynthMesh m = makeRing(fw, 24, 10.0, 5.0);
        const std::vector<int> none;
        check(claimChordSagitta(m.view(), none, fw.axis(), 10.0) == 0.0, "sigma empty claim == 0");
    }

    // -- topology: N=16, R=10, two rows. Ids: T1_i = 2i, T2_i = 2i+1. ---------
    // star(b_i) = { T1_{i-1}, T2_i, T1_i } = { 2(i-1), 2i+1, 2i }.
    {
        const int N = 16;
        const double R = 10.0;
        const Frame fp = placedFrame();
        const double sigma = R * (1.0 - std::cos(kPi / N));

        // (a) the intact ring
        {
            const SynthMesh m = makeRing(fp, N, R, 5.0);
            const UnionCensus u = unionCensus(m.view(), allTris(m), fp.axis(), R);
            std::snprintf(name, sizeof name,
                          "intact ring: nTri=%zu edgePieces=%zu punctures=%zu domainFaces=%zu pinch=%zu domainTris=%zu sigma=%.9g",
                          u.nTri, u.edgePieces, u.punctures, u.domainFaces, u.pinchVertices, u.domainTris, u.sigma);
            check(u.nTri == 2u * N && u.edgePieces == 1 && u.punctures == 0 && u.domainFaces == 1 &&
                      u.pinchVertices == 0 && u.domainTris == 2u * N && std::abs(u.sigma - sigma) <= eps &&
                      u.sizes.size() == 1 && u.sizes[0] == 2u * N && u.pieces.size() == 1 && u.pieces[0] == 1,
                  name);
        }

        // (b)/(c) two opposite stars refused (i = 3 and i = 11): the claim is two
        // open strips of N-3 triangles sharing no vertex (edgePieces = 2). The
        // two centre vertices are then the only non-claim vertices adjacent to
        // the claim. Pulled inward by sigma/2 they are PUNCTURES (U-R6): both
        // stars rejoin the domain, one face, no pinch (b). Pulled inward by
        // 2 sigma they are EXTERIOR: two faces, still no pinch (c).
        const int i1 = 3, i2 = 11;
        auto star = [&](int i) {
            return std::vector<int>{2 * ((i - 1 + N) % N), 2 * i + 1, 2 * i};
        };
        std::vector<int> drop = star(i1);
        for (int t : star(i2)) drop.push_back(t);
        for (int variant = 0; variant < 2; variant++) {
            SynthMesh m = makeRing(fp, N, R, 5.0);
            const double pull = variant == 0 ? 0.5 * sigma : 2.0 * sigma;
            for (int i : {i1, i2}) {
                const double u = 2.0 * kPi * i / N;
                m.pts[static_cast<std::size_t>(i)] = fp.at(R - pull, u, 0.0);
            }
            const std::vector<int> claim = allBut(m, drop);
            const UnionCensus u = unionCensus(m.view(), claim, fp.axis(), R);
            std::snprintf(name, sizeof name,
                          "two stars refused, centres at %s sigma: nTri=%zu edgePieces=%zu punctures=%zu domainFaces=%zu pinch=%zu domainTris=%zu",
                          variant == 0 ? "0.5" : "2.0", u.nTri, u.edgePieces, u.punctures, u.domainFaces,
                          u.pinchVertices, u.domainTris);
            if (variant == 0) {
                check(u.nTri == 2u * N - 6 && u.edgePieces == 2 && u.punctures == 2 && u.domainFaces == 1 &&
                          u.pinchVertices == 0 && u.domainTris == 2u * N && u.sizes.size() == 1 &&
                          u.sizes[0] == 2u * N - 6 && u.pieces.size() == 1 && u.pieces[0] == 2 &&
                          std::abs(u.sigma - sigma) <= eps,
                      name);
            } else {
                check(u.nTri == 2u * N - 6 && u.edgePieces == 2 && u.punctures == 0 && u.domainFaces == 2 &&
                          u.pinchVertices == 0 && u.domainTris == 2u * N - 6 && u.sizes.size() == 2 &&
                          u.sizes[0] == static_cast<std::size_t>(N - 3) && u.sizes[1] == static_cast<std::size_t>(N - 3) &&
                          u.pieces.size() == 2 && u.pieces[0] == 1 && u.pieces[1] == 1,
                      name);
            }
        }

        // (d) one triangle refused, T2_5 = (b_5, t_6, t_5): every vertex stays a
        // claim vertex, so there is no puncture candidate. At b_5 the remaining
        // claim triangles T1_4 = (b_4, b_5, t_5) and T1_5 = (b_5, b_6, t_6) share
        // no edge through b_5 -> two fans -> ONE pinch vertex (U-R8). At t_5 and
        // t_6 the refused triangle was the END of the fan -> one fan each.
        {
            const SynthMesh m = makeRing(fp, N, R, 5.0);
            const std::vector<int> claim = allBut(m, {2 * 5 + 1});
            const UnionCensus u = unionCensus(m.view(), claim, fp.axis(), R);
            std::snprintf(name, sizeof name,
                          "one triangle refused: nTri=%zu edgePieces=%zu punctures=%zu domainFaces=%zu pinch=%zu domainTris=%zu",
                          u.nTri, u.edgePieces, u.punctures, u.domainFaces, u.pinchVertices, u.domainTris);
            check(u.nTri == 2u * N - 1 && u.edgePieces == 1 && u.punctures == 0 && u.domainFaces == 1 &&
                      u.pinchVertices == 1 && u.domainTris == 2u * N - 1,
                  name);
        }

        // (e) two stacked rings on ONE cylinder, claimed together: they share no
        // vertex and no puncture can bridge them (the gap vertices ARE claim
        // vertices) -> two domain faces. U-R7 counts, U-R3 is the caller's.
        {
            SynthMesh m = makeRing(fp, N, R, 5.0);
            const int base = static_cast<int>(m.pts.size());
            for (int row = 0; row < 2; row++)
                for (int i = 0; i < N; i++)
                    m.pts.push_back(fp.at(R, 2.0 * kPi * i / N, 20.0 + (row ? 5.0 : 0.0)));
            for (int i = 0; i < N; i++) {
                const int b0 = base + i, b1 = base + (i + 1) % N, t0 = base + N + i, t1 = base + N + (i + 1) % N;
                m.addTri(b0, b1, t1);
                m.addTri(b0, t1, t0);
            }
            m.edges.clear();
            m.finish();
            const UnionCensus u = unionCensus(m.view(), allTris(m), fp.axis(), R);
            std::snprintf(name, sizeof name,
                          "two stacked rings: nTri=%zu edgePieces=%zu punctures=%zu domainFaces=%zu pinch=%zu",
                          u.nTri, u.edgePieces, u.punctures, u.domainFaces, u.pinchVertices);
            check(u.nTri == 4u * N && u.edgePieces == 2 && u.punctures == 0 && u.domainFaces == 2 &&
                      u.pinchVertices == 0,
                  name);
        }
    }

    std::fprintf(stderr, "union_census_unit: %d passed, %d failed\n", gPass, gFail);
    return gFail == 0 ? 0 : 1;
}
