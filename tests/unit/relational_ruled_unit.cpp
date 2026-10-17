// RULED QUAD relation (D-fc-5.9): one case per certificate step (rails, caps,
// corners, line first, interior evidence, vertex test, same surface, dVol),
// the stage on a synthetic twisted square prism, and the prism end to end
// through convert(). No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_relational_ruled.hpp"
#include "stl2step/stl2step.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kQ = 1e-6;
static constexpr double kTau = 2 * kQ;
static constexpr int kIters = 8;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;

    int add(const gp_XYZ& p) {
        pts.push_back(p);
        return (int)pts.size() - 1;
    }
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
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2],
                              triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]);
                const int hi = std::max(v[s], v[s + 1]);
                auto it = ids.find({lo, hi});
                int id;
                if (it == ids.end()) {
                    id = (int)edges.size();
                    ids.emplace(std::make_pair(lo, hi), id);
                    edges.emplace_back(lo, hi);
                } else {
                    id = it->second;
                }
                triEdges[t][(std::size_t)s] = id;
                if (v[s] < v[s + 1]) triDirs[t] = (uint8_t)(triDirs[t] | (uint8_t)(1u << s));
            }
        }
        compVtx.resize(pts.size());
        for (std::size_t i = 0; i < pts.size(); ++i) compVtx[i] = (int)i;
    }
    MeshView view() const {
        MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int (*)[3]>(triFlat.data());
        mv.compTris = compTris.data();
        mv.compVtx = compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = triEdges.data();
        mv.triDirs = triDirs.data();
        mv.nTri = triFlat.size() / 3;
        mv.nVtx = pts.size();
        mv.nEdge = edges.size();
        mv.diag = 100.0;
        mv.quantFloor = kQ;
        return mv;
    }
};

// A square prism of circumradius R and height H whose top is turned by twist
// radians: four doubly ruled walls (nu x nv grids, or boundary-only walls when
// noInterior), two planar caps fanned from their centres. Outward winding.
struct Prism {
    std::array<gp_XYZ, 4> bot, top;
    // wall k: P00 = bot[k], P10 = bot[k+1], P01 = top[k], P11 = top[k+1]
    RuledPoles wall(int k) const {
        return {bot[(std::size_t)k], bot[(std::size_t)((k + 1) % 4)], top[(std::size_t)k],
                top[(std::size_t)((k + 1) % 4)]};
    }
    std::vector<int> wallInterior[4];  // interior vertex ids per wall
    std::vector<int> railMid;          // a mid vertex of rail k (shared by walls k-1, k)
    int botCentre = -1;
};

static Prism buildPrism(Synth& s, double R, double H, double twist, int nu, int nv,
                        bool noInterior = false) {
    Prism P;
    for (int k = 0; k < 4; ++k) {
        const double a = 0.25 * kPi + 0.5 * kPi * k;
        P.bot[(std::size_t)k] = gp_XYZ(R * std::cos(a), R * std::sin(a), 0.0);
        P.top[(std::size_t)k] = gp_XYZ(R * std::cos(a + twist), R * std::sin(a + twist), H);
    }
    // id[k][i][j]; column i = nu of wall k is column 0 of wall k+1.
    std::vector<std::vector<std::vector<int>>> id(
        4, std::vector<std::vector<int>>((std::size_t)nu + 1, std::vector<int>((std::size_t)nv + 1, -1)));
    for (int k = 0; k < 4; ++k) {
        const RuledPoles W = P.wall(k);
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j <= nv; ++j) {
                const bool interior = i > 0 && j > 0 && j < nv;
                if (noInterior && interior && !(i == 0 || j == 0 || j == nv)) continue;
                id[(std::size_t)k][(std::size_t)i][(std::size_t)j] =
                    s.add(ruledPoint(W, (double)i / nu, (double)j / nv));
                if (interior) P.wallInterior[k].push_back(id[(std::size_t)k][(std::size_t)i][(std::size_t)j]);
            }
    }
    for (int k = 0; k < 4; ++k)
        for (int j = 0; j <= nv; ++j)
            id[(std::size_t)k][(std::size_t)nu][(std::size_t)j] =
                id[(std::size_t)((k + 1) % 4)][0][(std::size_t)j];
    for (int k = 0; k < 4; ++k) P.railMid.push_back(id[(std::size_t)k][0][(std::size_t)(nv / 2)]);
    for (int k = 0; k < 4; ++k) {
        auto at = [&](int i, int j) { return id[(std::size_t)k][(std::size_t)i][(std::size_t)j]; };
        if (noInterior) {
            // nu = nv = 2 with the centre removed: boundary vertices only.
            s.addTri(at(0, 0), at(1, 0), at(0, 1));
            s.addTri(at(1, 0), at(2, 0), at(2, 1));
            s.addTri(at(2, 1), at(2, 2), at(1, 2));
            s.addTri(at(1, 2), at(0, 2), at(0, 1));
            s.addTri(at(1, 0), at(2, 1), at(0, 1));
            s.addTri(at(2, 1), at(1, 2), at(0, 1));
            continue;
        }
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j < nv; ++j) {
                s.addTri(at(i, j), at(i + 1, j), at(i + 1, j + 1));
                s.addTri(at(i, j), at(i + 1, j + 1), at(i, j + 1));
            }
    }
    std::vector<int> ringB, ringT;
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < nu; ++i) {
            ringB.push_back(id[(std::size_t)k][(std::size_t)i][0]);
            ringT.push_back(id[(std::size_t)k][(std::size_t)i][(std::size_t)nv]);
        }
    P.botCentre = s.add(gp_XYZ(0, 0, 0));
    const int ct = s.add(gp_XYZ(0, 0, H));
    const int n = (int)ringB.size();
    for (int m = 0; m < n; ++m) {
        s.addTri(P.botCentre, ringB[(std::size_t)((m + 1) % n)], ringB[(std::size_t)m]);
        s.addTri(ct, ringT[(std::size_t)m], ringT[(std::size_t)((m + 1) % n)]);
    }
    s.finish();
    return P;
}

static void stageTols(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = kIters;
    tol = DerivedTols{};
    tol.epsMesh = 0.01;
    tol.epsPlane = 0.02;
    tol.thetaPlane = 2.0 * kPi / 180.0;
    tol.thetaSharp = 30.0 * kPi / 180.0;
    tol.thetaCylLo = 5.0 * kPi / 180.0;
    tol.thetaCylHi = 60.0 * kPi / 180.0;
    tol.thetaBin = 0.25 * kPi / 180.0;
}

static std::vector<const Region*> runStage(const Synth& s, SegmentWork& work) {
    SegmentParams p;
    DerivedTols tol;
    stageTols(p, tol);
    const MeshView mv = s.view();
    std::vector<const Region*> out;
    if (!chartsA1(mv, p, tol, work) || !growProvisionalA2(mv, p, tol, work) ||
        !claimRuledQuads(mv, p, tol, work))
        return out;
    for (const Region& r : work.accepted)
        if (r.type == SurfType::BSpline && r.origin == Origin::RuledQuad) out.push_back(&r);
    return out;
}

static bool samePoint(const gp_XYZ& a, const gp_XYZ& b, double tol) {
    return (a - b).Modulus() <= tol;
}

// Exact surface-minus-chord volume of one parameter cell [u0,u1] x [v0,v1] of S
// against its two chord triangles (a,b,c), (a,c,d) -- the cell's iso-lines are
// the triangles' outer edges, so the region between them is closed and the
// divergence theorem gives it: (1/3) int S . (S_u x S_v) du dv (a polynomial of
// degree <= 3 in each of u, v: 3-point Gauss-Legendre is exact) minus the two
// triangles' cone volumes.
static double dVolCellExact(const RuledPoles& P, double u0, double u1, double v0, double v1) {
    const double g[3] = {-std::sqrt(0.6), 0.0, std::sqrt(0.6)}, w[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    double surf = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const double u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * g[i];
            const double v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * g[j];
            surf += w[i] * w[j] * ruledPoint(P, u, v).Dot(ruledNormal(P, u, v));
        }
    surf *= 0.25 * (u1 - u0) * (v1 - v0) / 3.0;
    const gp_XYZ a = ruledPoint(P, u0, v0), b = ruledPoint(P, u1, v0), c = ruledPoint(P, u1, v1),
                 d = ruledPoint(P, u0, v1);
    return surf - (a.Dot(b.Crossed(c)) + a.Dot(c.Crossed(d))) / 6.0;
}

static void predicateCases() {
    // Rails / caps: line-first certificate of a chain.
    {
        const gp_XYZ o(1, 2, 3), d = gp_XYZ(1, 2, 2) / 3.0, e(2.0 / 3, -2.0 / 3, 1.0 / 3);
        std::vector<gp_XYZ> pts;
        for (int i = 0; i < 6; ++i) pts.push_back(o + d * (2.0 * i) + e * ((i % 2 ? 0.4 : -0.4) * kTau));
        ChainLineCert C = certifyChainLine(pts, kTau);
        expect(C.ok && C.stations == 6 && C.maxDev <= kTau, "line.straight-chain-certifies");
        expect(std::fabs(std::fabs(C.d.Dot(d)) - 1.0) < 1e-12, "line.direction-recovered");
        pts[3] = pts[3] + e * (3 * kTau);
        C = certifyChainLine(pts, kTau);
        expect(!C.ok && std::strcmp(C.why, "vertex-off") == 0, "line.bent-chain-vertex-off");
        std::vector<gp_XYZ> two{o, o + d * 5.0};
        C = certifyChainLine(two, kTau);
        expect(!C.ok && std::strcmp(C.why, "underdetermined") == 0, "line.two-stations-underdetermined");
        std::vector<gp_XYZ> merged{o, o + d * (0.5 * kTau), o + d * 5.0};
        C = certifyChainLine(merged, kTau);
        expect(!C.ok && C.stations == 2, "line.stations-merge-at-tau");
    }
    // Caps: plane anchor.
    expect(ruledCapAnchor(true, true, 0.5 * kTau, kTau), "cap.certified-anchor");
    expect(!ruledCapAnchor(true, true, 2 * kTau, kTau), "cap.uncertified-plane-refused");
    expect(!ruledCapAnchor(false, true, 0.0, kTau), "cap.claimed-provisional-refused");
    expect(!ruledCapAnchor(true, false, 0.0, kTau), "cap.shatter-class-refused");
    // Corners: rail ∩ cap plane.
    {
        gp_XYZ c;
        const bool ok = railCapCorner(gp_XYZ(1, 1, 0), gp_XYZ(0.6, 0, 0.8), gp_XYZ(0, 0, 8),
                                      gp_XYZ(0, 0, 1), c);
        expect(ok && samePoint(c, gp_XYZ(7, 1, 8), 1e-12), "corner.rail-meets-cap");
        expect(!railCapCorner(gp_XYZ(1, 1, 0), gp_XYZ(1, 0, 0), gp_XYZ(0, 0, 8), gp_XYZ(0, 0, 1), c),
               "corner.parallel-rail-refused");
    }
    // Line first on the surface.
    const RuledPoles flat{gp_XYZ(0, 0, 0), gp_XYZ(10, 0, 0), gp_XYZ(0, 0, 10), gp_XYZ(10, 0, 10)};
    const RuledPoles twisted{gp_XYZ(0, 0, 0), gp_XYZ(10, 0, 0), gp_XYZ(0, 1, 10),
                             gp_XYZ(10, -1, 10)};
    // z = h x y over [-1, 1]^2: the corners sit +-h off their LS plane z = 0.
    const double h = 0.7;
    const RuledPoles saddle{gp_XYZ(-1, -1, h), gp_XYZ(1, -1, -h), gp_XYZ(-1, 1, -h),
                            gp_XYZ(1, 1, h)};
    expect(ruledCornerPlanarity(flat) <= 1e-12, "line-first.planar-corners");
    expect(std::fabs(ruledCornerPlanarity(saddle) - h) < 1e-12, "line-first.twisted-corners");
    // Vertex test (max, on the face) and same surface (unbounded S).
    {
        std::vector<gp_XYZ> on;
        for (int i = 0; i <= 8; ++i)
            for (int j = 0; j <= 8; ++j) on.push_back(ruledPoint(twisted, i / 8.0, j / 8.0));
        RuledVertexCert V = certifyRuledVertices(twisted, on, kTau, kIters);
        expect(V.ok && V.maxDev <= 1e-9, "vertex.on-surface-certifies");
        gp_XYZ nrm = ruledNormal(twisted, 0.3, 0.6);
        nrm /= nrm.Modulus();
        on.push_back(ruledPoint(twisted, 0.3, 0.6) + nrm * (0.6 * kTau));
        V = certifyRuledVertices(twisted, on, kTau, kIters);
        expect(V.ok && std::fabs(V.maxDev - 0.6 * kTau) < 1e-3 * kTau, "vertex.offset-within-tau");
        on.back() = ruledPoint(twisted, 0.3, 0.6) + nrm * (3 * kTau);
        V = certifyRuledVertices(twisted, on, kTau, kIters);
        expect(!V.ok, "vertex.one-vertex-off-refuses");
        std::vector<gp_XYZ> ext;
        for (int i = 0; i <= 4; ++i)
            for (int j = 0; j <= 4; ++j) ext.push_back(ruledPoint(twisted, 1 + i / 4.0, j / 4.0));
        expect(ruledSameSurface(twisted, ext, kTau, kIters), "same-surface.extension-joins");
        expect(!certifyRuledVertices(twisted, ext, kTau, kIters).ok,
               "vertex.outside-face-domain-refuses");
        gp_XYZ ne = ruledNormal(twisted, 1.5, 0.25);
        ext[7] = ext[7] + ne / ne.Modulus() * (3 * kTau);
        expect(!ruledSameSurface(twisted, ext, kTau, kIters), "same-surface.off-chart-stays");
    }
    // dVol: the 7-point rule over a cell's two chord triangles against the
    // exact cell volume; a planar patch gives 0.
    {
        const RuledPoles W{gp_XYZ(0, 0, 0), gp_XYZ(10, 0, 0), gp_XYZ(0, 3, 10), gp_XYZ(10, -3, 10)};
        const double u0 = 0.25, u1 = 0.375, v0 = 0.5, v1 = 0.625;
        const gp_XYZ a = ruledPoint(W, u0, v0), b = ruledPoint(W, u1, v0),
                     c = ruledPoint(W, u1, v1), d = ruledPoint(W, u0, v1);
        const double d7 =
            dVolRuledTriangle(W, a, b, c, kIters) + dVolRuledTriangle(W, a, c, d, kIters);
        const double dx = dVolCellExact(W, u0, u1, v0, v1);
        expect(std::fabs(dx) > 0 && std::fabs(d7 - dx) <= 1e-3 * std::fabs(dx),
               "dvol.seven-point-matches-exact-cell");
        expect(std::fabs(dVolRuledTriangle(flat, ruledPoint(flat, 0, 0), ruledPoint(flat, 1, 0),
                                           ruledPoint(flat, 1, 1), kIters)) <= 1e-12,
               "dvol.planar-patch-zero");
    }
}

static void stageCases() {
    const double R = 10, H = 20, tw = 30 * kPi / 180;
    {
        Synth s;
        const Prism P = buildPrism(s, R, H, tw, 8, 10);
        SegmentWork work;
        const auto got = runStage(s, work);
        expect(got.size() == 4, "stage.twisted-prism-four-quads");
        bool poles = got.size() == 4;
        double worst = 0, dvErr = 0;
        for (const Region* r : got) {
            int matched = -1;
            for (int k = 0; k < 4 && matched < 0; ++k) {
                const RuledPoles W = P.wall(k);
                int hit = 0;
                for (const gp_XYZ& p : r->ruledPoles)
                    for (const gp_XYZ& q : W)
                        if (samePoint(p, q, 1e-9)) ++hit;
                if (hit == 4) matched = k;
            }
            // P00 / P11 are a diagonal: S is the same patch.
            const gp_XYZ mid = ruledPoint(r->ruledPoles, 0.5, 0.5);
            if (matched < 0 || !samePoint(mid, ruledPoint(P.wall(matched), 0.5, 0.5), 1e-9))
                poles = false;
            worst = std::max(worst, r->maxVertexDev);
            double exact = 0;
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 10; ++j)
                    exact += dVolCellExact(P.wall(matched < 0 ? 0 : matched), i / 8.0,
                                           (i + 1) / 8.0, j / 10.0, (j + 1) / 10.0);
            dvErr = std::max(dvErr, std::fabs(r->dVolPredicted - exact) / std::fabs(exact));
            if (!r->outwardNormal) poles = false;
        }
        expect(poles, "stage.poles-are-rail-cap-corners");
        expect(worst <= kTau, "stage.vertex-certificate-within-tau");
        expect(dvErr <= 1e-3, "stage.dvol-surface-minus-chord");
        int consumed = 0;
        for (const Provisional& pr : work.provisionals)
            if (pr.claim == ProvClaim::ConsumedRuled) ++consumed;
        expect(consumed > 0, "stage.wall-provisionals-consumed");
    }
    {
        Synth s;
        buildPrism(s, R, H, 0.0, 8, 10);
        SegmentWork work;
        expect(runStage(s, work).empty(), "stage.untwisted-walls-are-planes-line-first");
    }
    {
        Synth s;
        const Prism P = buildPrism(s, R, H, tw, 8, 10);
        const int v = P.wallInterior[0][P.wallInterior[0].size() / 2];
        gp_XYZ n = ruledNormal(P.wall(0), 0.5, 0.5);
        s.pts[(std::size_t)v] += n / n.Modulus() * (3 * kTau);
        SegmentWork work;
        expect(runStage(s, work).size() == 3, "stage.vertex-off-wall-refused-whole");
    }
    {
        Synth s;
        const Prism P = buildPrism(s, R, H, tw, 8, 10);
        s.pts[(std::size_t)P.railMid[1]] += gp_XYZ(0, 0, 0) + gp_XYZ(3 * kTau, -3 * kTau, 0);
        SegmentWork work;
        expect(runStage(s, work).size() == 2, "stage.bent-rail-refuses-both-walls");
    }
    {
        Synth s;
        const Prism P = buildPrism(s, R, H, tw, 8, 10);
        s.pts[(std::size_t)P.botCentre] += gp_XYZ(0, 0, 3 * kTau);
        SegmentWork work;
        expect(runStage(s, work).empty(), "stage.uncertified-cap-no-anchor");
    }
    {
        Synth s;
        buildPrism(s, R, H, tw, 8, 1);
        SegmentWork work;
        expect(runStage(s, work).empty(), "stage.two-station-rails-underdetermined");
    }
    {
        Synth s;
        buildPrism(s, R, H, 10 * kPi / 180, 2, 2, true);
        SegmentWork work;
        expect(runStage(s, work).empty(), "stage.no-interior-vertex-underdetermined");
    }
}

// The twisted prism as a float32 binary STL, end to end: four
// B_SPLINE_SURFACE_WITH_KNOTS faces, a shell that validates (no revert), and a
// B-Rep volume equal to the mesh plus the predicted surface-minus-chord volume.
static void endToEnd(const std::string& dir) {
    Synth s;
    // Coarse on purpose: a twisted wall breaks into many A2 plane anchors and
    // the extrude stage (H0) evaluates every anchor pair.
    buildPrism(s, 10, 20, 20 * kPi / 180, 6, 8);
    const std::string stl = dir + "/ruled_prism.stl", step = dir + "/ruled_prism.step";
    {
        std::ofstream f(stl, std::ios::binary);
        char hdr[80] = {};
        f.write(hdr, 80);
        const uint32_t n = (uint32_t)(s.triFlat.size() / 3);
        f.write(reinterpret_cast<const char*>(&n), 4);
        for (uint32_t t = 0; t < n; ++t) {
            const gp_XYZ& a = s.pts[(std::size_t)s.triFlat[3 * t]];
            const gp_XYZ& b = s.pts[(std::size_t)s.triFlat[3 * t + 1]];
            const gp_XYZ& c = s.pts[(std::size_t)s.triFlat[3 * t + 2]];
            gp_XYZ nn = (b - a).Crossed(c - a);
            nn /= nn.Modulus();
            const float rec[12] = {(float)nn.X(), (float)nn.Y(), (float)nn.Z(),
                                   (float)a.X(), (float)a.Y(), (float)a.Z(),
                                   (float)b.X(), (float)b.Y(), (float)b.Z(),
                                   (float)c.X(), (float)c.Y(), (float)c.Z()};
            f.write(reinterpret_cast<const char*>(rec), sizeof rec);
            const uint16_t attr = 0;
            f.write(reinterpret_cast<const char*>(&attr), 2);
        }
    }
    stl2step::Options opt;
    opt.input = stl;
    opt.output = step;
    opt.smooth = true;
    opt.verify = true;
    opt.threads = 1;
    const stl2step::Result r = stl2step::convert(opt);
    std::fprintf(stderr,
                 "e2e ok=%d solids=%d reverted=%d mesh=%.9f step=%.9f pred=%.9f dPct=%.6g "
                 "faces=%d\n",
                 r.ok ? 1 : 0, r.solids, r.smoothRevertedComponents, r.meshVolumeMM3,
                 r.stepVolumeMM3, r.smoothVolPredictedMM3, r.volumeDeltaPct, r.facesAfterUnify);
    expect(r.ok && r.solids == 1 && r.openShells == 0, "e2e.one-solid");
    expect(r.smoothRevertedComponents == 0, "e2e.shell-validates-no-revert");
    std::ifstream in(step);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    int nb = 0;
    for (std::size_t at = text.find("B_SPLINE_SURFACE_WITH_KNOTS"); at != std::string::npos;
         at = text.find("B_SPLINE_SURFACE_WITH_KNOTS", at + 1))
        ++nb;
    expect(nb == 4, "e2e.four-bspline-walls");
    expect(r.facesAfterUnify == 6, "e2e.six-faces");
    const double gap = std::fabs(r.stepVolumeMM3 - (r.meshVolumeMM3 + r.smoothVolPredictedMM3));
    expect(r.smoothVolPredictedMM3 != 0.0 && gap <= 1e-3 * std::fabs(r.smoothVolPredictedMM3),
           "e2e.volume-is-mesh-plus-predicted");
}

int main(int argc, char** argv) {
    predicateCases();
    stageCases();
    if (argc > 1) endToEnd(argv[1]);
    else expect(false, "e2e.work-directory-argument");
    std::fprintf(stderr, "relational_ruled_unit: %s (%d failures)\n", gFail ? "FAIL" : "PASS",
                 gFail);
    return gFail ? 1 : 0;
}
