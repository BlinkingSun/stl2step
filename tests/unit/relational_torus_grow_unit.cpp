// Synthetic cases for the same-surface growth of a certified torus (D-gid-7.12):
// (a) unclaimed tau-triangles edge-connected to the band join it; an off-surface
// vertex, a one-station cap chord and a withheld triangle stay out; (b) an
// uncertified B1 chip lying wholly on the torus becomes one face with the
// unclaimed tau-triangles joining it, a certified or off-surface chip stays;
// the heal-not-revert helpers name a grown band's growth first; (c) the fill
// enters an uncertified junction chip triangle by triangle when its off-torus
// triangle lies on a certified plane across an edge, the chip dissolves and that
// triangle is released; an uncertified neighbour, no neighbour or a withheld
// triangle keeps the chip whole (D-gid-8.4). Generic quarter-round geometry
// only; no training mesh.
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit_relational.hpp"
#include "refit_relational_torus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

#include <gp_Ax3.hxx>

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;

struct Synth {
    std::vector<gp_XYZ> pts;
    std::vector<int> triFlat;
    std::vector<int> compTris, compVtx;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<uint8_t> triDirs;
    double q = 1e-6;

    int add(const gp_XYZ& p) {
        pts.push_back(p);
        return (int)pts.size() - 1;
    }
    int addTri(int a, int b, int c) {
        triFlat.push_back(a);
        triFlat.push_back(b);
        triFlat.push_back(c);
        return (int)(triFlat.size() / 3) - 1;
    }
    void finish() {
        std::map<std::pair<int, int>, int> ids;
        const std::size_t nT = triFlat.size() / 3;
        compTris.resize(nT);
        triEdges.resize(nT);
        triDirs.assign(nT, 0);
        for (std::size_t t = 0; t < nT; ++t) {
            compTris[t] = (int)t;
            const int v[4] = {triFlat[3 * t], triFlat[3 * t + 1], triFlat[3 * t + 2], triFlat[3 * t]};
            for (int s = 0; s < 3; ++s) {
                const int lo = std::min(v[s], v[s + 1]), hi = std::max(v[s], v[s + 1]);
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
        mv.quantFloor = q;
        return mv;
    }
};

// A quarter-round torus (major Rmaj about +Z through (0, 0, z0), minor r) as
// kRows rows x kCols sectors of two triangles, closed in azimuth. ring[j][i] is
// the vertex of station j at column i.
constexpr int kRows = 6, kCols = 24;
constexpr double kRmaj = 12.0, kRmin = 2.0, kZ0 = 8.0;

struct Part {
    Synth mesh;
    std::vector<std::vector<int>> ring;
    // tri[j][i][0..1]: the two triangles of row j, sector i.
    std::vector<std::vector<std::array<int, 2>>> tri;
    SegmentWork work;
};

static gp_XYZ onTorus(int j, int i, double off) {
    const double phi = kPi + 0.5 * kPi * (double)j / (double)kRows;
    const double th = 2.0 * kPi * (double)i / (double)kCols;
    const double rho = kRmaj + (kRmin + off) * std::cos(phi);
    const double z = kZ0 + (kRmin + off) * std::sin(phi);
    return gp_XYZ(rho * std::cos(th), rho * std::sin(th), z);
}

// `offRow`/`offCol` names one vertex moved `off` mm off the surface along its normal.
static Part quarterRound(int offRow = -1, int offCol = -1, double off = 0.0) {
    Part p;
    for (int j = 0; j <= kRows; ++j) {
        p.ring.emplace_back();
        for (int i = 0; i < kCols; ++i)
            p.ring.back().push_back(p.mesh.add(onTorus(j, i, (j == offRow && i == offCol) ? off : 0.0)));
    }
    for (int j = 0; j < kRows; ++j) {
        p.tri.emplace_back();
        for (int i = 0; i < kCols; ++i) {
            const int i1 = (i + 1) % kCols;
            const auto& a = p.ring[(std::size_t)j];
            const auto& b = p.ring[(std::size_t)j + 1];
            const int t0 = p.mesh.addTri(a[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i]);
            const int t1 = p.mesh.addTri(b[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i1]);
            p.tri.back().push_back({t0, t1});
        }
    }
    return p;
}

static std::vector<int> sectors(const Part& p, int lo, int hi) {
    std::vector<int> out;
    for (int j = 0; j < kRows; ++j)
        for (int i = lo; i <= hi; ++i)
            for (int t : p.tri[(std::size_t)j][(std::size_t)i]) out.push_back(t);
    std::sort(out.begin(), out.end());
    return out;
}

static Region band(const std::vector<int>& tris) {
    Region t;
    t.type = SurfType::Torus;
    t.origin = Origin::TorusBlend;
    t.ax = gp_Ax3(gp_Pnt(0, 0, kZ0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    t.radius = kRmaj;
    t.radius2 = kRmin;
    t.closed360 = false;
    t.nSides = kCols;
    t.tris = tris;
    t.maxVertexDev = 0.0;
    return t;
}

static Region chip(const std::vector<int>& tris, double maxDev) {
    Region a;
    a.type = SurfType::Cylinder;
    a.origin = Origin::CylGrow;
    a.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0));
    a.radius = kRmin;
    a.closed360 = false;
    a.tris = tris;
    a.maxVertexDev = maxDev;
    return a;
}

static Region plane(const std::vector<int>& tris) {
    Region r;
    r.type = SurfType::Plane;
    r.origin = Origin::PlaneGrow;
    r.tris = tris;
    return r;
}

static bool has(const std::vector<int>& v, int t) { return std::find(v.begin(), v.end(), t) != v.end(); }

static std::vector<const Region*> tori(const SegmentWork& w) {
    std::vector<const Region*> out;
    for (const Region& r : w.accepted)
        if (r.origin == Origin::TorusBlend) out.push_back(&r);
    return out;
}

// (a) A band on sectors 0..7; sectors 8..23 unclaimed. One vertex of station 3
// at column 16 sits 10 tau off the surface: its 6 triangles stay out. A cap
// chord on the top station (three vertices on one ring) stays out. Every other
// unclaimed triangle joins the band, which keeps its surface.
static void caseFill() {
    const double q = 1e-6, tau = 2.0 * q;
    Part p = quarterRound(3, 16, 10.0 * tau);
    const int cap = p.mesh.addTri(p.ring[kRows][1], p.ring[kRows][0], p.ring[kRows][12]);
    p.mesh.finish();
    p.work.accepted.push_back(band(sectors(p, 0, 7)));
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    expect(claimTorusSameSurface(mv, sp, p.work), "fill: stage returns true");
    const auto ts = tori(p.work);
    expect(ts.size() == 1, "fill: one torus");
    if (ts.size() != 1) return;
    const Region& T = *ts[0];
    std::vector<int> off;
    for (std::size_t t = 0; t < mv.nTri; ++t)
        for (int k = 0; k < 3; ++k)
            if (p.mesh.triFlat[3 * t + (std::size_t)k] == p.ring[3][16]) off.push_back((int)t);
    expect(off.size() == 6u, "fill: the off vertex has 6 triangles");
    expect(T.tris.size() == (std::size_t)(kRows * kCols * 2) - off.size(),
           "fill: every on-surface triangle joins");
    bool anyOff = false;
    for (int t : off) anyOff = anyOff || has(T.tris, t);
    expect(!anyOff, "fill: a vertex 10 tau off keeps its triangles out");
    expect(!has(T.tris, cap), "fill: a one-station cap chord stays out (F10)");
    expect(T.torusGrowTris.size() == T.tris.size() - 8u * 2u * kRows, "fill: growth recorded");
    expect(T.maxVertexDev <= tau, "fill: union within tau");
    expect(T.radius == kRmaj && T.radius2 == kRmin, "fill: surface unchanged (no refit)");
    expect(!T.closed360, "fill: closure unchanged");
}

// (a) under U-3: a withheld triangle is never taken.
static void caseWithheld() {
    const double tau = 2e-6;
    Part p = quarterRound();
    p.mesh.finish();
    p.work.accepted.push_back(band(sectors(p, 0, 7)));
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    const std::vector<int> held = {p.tri[0][8][0], p.tri[0][8][1]};
    sp.torusWithheldTris = held;
    claimTorusSameSurface(mv, sp, p.work);
    const auto ts = tori(p.work);
    expect(ts.size() == 1 && !has(ts[0]->tris, held[0]) && !has(ts[0]->tris, held[1]),
           "withheld: held triangles stay out");
    expect(ts.size() == 1 && ts[0]->tris.size() == (std::size_t)(kRows * kCols * 2) - 2u,
           "withheld: the rest still joins around them");
    expect(ts.size() == 1 && ts[0]->maxVertexDev <= tau, "withheld: within tau");
}

// (b) Band on sectors 0..7; planes own sectors 8..9 and 15..16 (they cut the
// surface). An uncertified chip (100 tau) on sectors 10..13 and the unclaimed
// sector 14 beside it become one torus face of their own; the chip dissolves.
// A certified chip (0.5 tau) on sectors 17..18 is not robbed; the unclaimed
// sectors 19..23 join the band across the azimuth wrap.
static void caseChipPiece() {
    const double tau = 2e-6;
    Part p = quarterRound();
    p.mesh.finish();
    p.work.accepted.push_back(band(sectors(p, 0, 7)));
    p.work.accepted.push_back(plane(sectors(p, 8, 9)));
    p.work.accepted.push_back(chip(sectors(p, 10, 13), 100.0 * tau));
    p.work.accepted.push_back(plane(sectors(p, 15, 16)));
    p.work.accepted.push_back(chip(sectors(p, 17, 18), 0.5 * tau));
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    claimTorusSameSurface(mv, sp, p.work);
    const auto ts = tori(p.work);
    expect(ts.size() == 2, "chip: band + one new piece");
    if (ts.size() != 2) return;
    const Region* B = has(ts[0]->tris, p.tri[0][0][0]) ? ts[0] : ts[1];
    const Region* N = (B == ts[0]) ? ts[1] : ts[0];
    std::vector<int> wantB = sectors(p, 0, 7), wantN = sectors(p, 10, 14);
    const std::vector<int> wrap = sectors(p, 19, 23);
    wantB.insert(wantB.end(), wrap.begin(), wrap.end());
    std::sort(wantB.begin(), wantB.end());
    std::vector<int> gotB = B->tris, gotN = N->tris;
    std::sort(gotB.begin(), gotB.end());
    std::sort(gotN.begin(), gotN.end());
    expect(gotB == wantB, "chip: band grows over the unclaimed wrap sectors only");
    expect(gotN == wantN, "chip: chip + unclaimed neighbour are one face");
    expect(N->torusGrowTris == wantN, "chip: the piece is all growth");
    expect(N->radius == kRmaj && N->radius2 == kRmin && N->maxVertexDev <= tau,
           "chip: the piece is on the band's certified surface");
    int cyl = 0, cylCert = 0;
    for (const Region& r : p.work.accepted)
        if (r.origin == Origin::CylGrow) {
            ++cyl;
            if (r.maxVertexDev <= tau) ++cylCert;
        }
    expect(cyl == 1 && cylCert == 1, "chip: uncertified chip dissolved, certified chip kept");
}

// (b) An uncertified chip with one vertex 10 tau off the torus is not robbed.
static void caseChipOff() {
    const double tau = 2e-6;
    Part p = quarterRound(3, 12, 10.0 * tau);
    p.mesh.finish();
    p.work.accepted.push_back(band(sectors(p, 0, 7)));
    p.work.accepted.push_back(plane(sectors(p, 8, 9)));
    p.work.accepted.push_back(chip(sectors(p, 10, 13), 100.0 * tau));
    p.work.accepted.push_back(plane(sectors(p, 14, 23)));
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    claimTorusSameSurface(mv, sp, p.work);
    int cyl = 0;
    for (const Region& r : p.work.accepted)
        if (r.origin == Origin::CylGrow) ++cyl;
    expect(tori(p.work).size() == 1 && cyl == 1, "chip off: a chip with a vertex off tau stays");
}

// (c) D-gid-8.4 fixtures. The quarter round ends at station kRows tangent to the
// plane z = kZ0 - kRmin; a skirt ring at rho = kRmaj + kSkirt on that plane
// gives two triangles per sector that lie on the plane and are off the torus.
constexpr double kSkirt = 0.5;

static std::vector<std::array<int, 2>> addSkirt(Part& p) {
    std::vector<int> P;
    for (int i = 0; i < kCols; ++i) {
        const double th = 2.0 * kPi * (double)i / (double)kCols;
        const double rho = kRmaj + kSkirt;
        P.push_back(p.mesh.add(gp_XYZ(rho * std::cos(th), rho * std::sin(th), kZ0 - kRmin)));
    }
    std::vector<std::array<int, 2>> sk;
    const auto& b = p.ring[(std::size_t)kRows];
    for (int i = 0; i < kCols; ++i) {
        const int i1 = (i + 1) % kCols;
        const int s0 = p.mesh.addTri(b[(std::size_t)i], P[(std::size_t)i], b[(std::size_t)i1]);
        const int s1 = p.mesh.addTri(b[(std::size_t)i1], P[(std::size_t)i], P[(std::size_t)i1]);
        sk.push_back({s0, s1});
    }
    return sk;
}

static Region skirtPlane(const std::vector<int>& tris, double maxDev) {
    Region r = plane(tris);
    r.ax = gp_Ax3(gp_Pnt(0, 0, kZ0 - kRmin), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.maxVertexDev = maxDev;
    return r;
}

// Band on sectors 0..7; an uncertified chip (100 tau) = torus sectors 8..13 plus
// the skirt triangle s0 of sector 12 (on the plane, off the torus); a plane owns
// every other skirt triangle. Torus sectors 14..23 are unclaimed.
struct Junction {
    Part p;
    int jtri = -1;
};
static Junction junction(double planeDev) {
    Junction j;
    j.p = quarterRound();
    const auto sk = addSkirt(j.p);
    j.p.mesh.finish();
    j.jtri = sk[12][0];
    std::vector<int> chipTris = sectors(j.p, 8, 13), planeTris;
    chipTris.push_back(j.jtri);
    std::sort(chipTris.begin(), chipTris.end());
    for (const auto& pr : sk)
        for (int t : pr)
            if (t != j.jtri) planeTris.push_back(t);
    std::sort(planeTris.begin(), planeTris.end());
    j.p.work.accepted.push_back(band(sectors(j.p, 0, 7)));
    j.p.work.accepted.push_back(chip(chipTris, 100.0 * 2e-6));
    j.p.work.accepted.push_back(skirtPlane(planeTris, planeDev));
    return j;
}

static int countOrigin(const SegmentWork& w, Origin o) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.origin == o) ++n;
    return n;
}

// The fill enters the junction chip, takes its torus triangles; the skirt
// triangle is within tau of the certified plane across its edge and is released.
static void caseChipJunction() {
    const double tau = 2e-6;
    Junction j = junction(0.0);
    MeshView mv = j.p.mesh.view();
    SegmentParams sp{};
    expect(claimTorusSameSurface(mv, sp, j.p.work), "junction: stage returns true");
    const auto ts = tori(j.p.work);
    expect(ts.size() == 1, "junction: one torus, no new piece");
    if (ts.size() != 1) return;
    std::vector<int> got = ts[0]->tris;
    std::sort(got.begin(), got.end());
    expect(got == sectors(j.p, 0, kCols - 1), "junction: the band takes every torus triangle");
    expect(ts[0]->maxVertexDev <= tau && ts[0]->radius == kRmaj && ts[0]->radius2 == kRmin,
           "junction: band within tau, surface unchanged");
    expect(countOrigin(j.p.work, Origin::CylGrow) == 0, "junction: the chip dissolves");
    bool owned = false;
    for (const Region& r : j.p.work.accepted) owned = owned || has(r.tris, j.jtri);
    expect(!owned, "junction: the off-torus triangle is released to no region");
    expect(countOrigin(j.p.work, Origin::PlaneGrow) == 1, "junction: the plane is untouched");
}

// Same, but the plane does not certify (3 tau): no certified surface explains the
// skirt triangle, so the chip stays whole and the band cannot pass it.
static void caseChipJunctionUncertified() {
    Junction j = junction(3.0 * 2e-6);
    MeshView mv = j.p.mesh.view();
    SegmentParams sp{};
    claimTorusSameSurface(mv, sp, j.p.work);
    const auto ts = tori(j.p.work);
    int cylTris = -1;
    for (const Region& r : j.p.work.accepted)
        if (r.origin == Origin::CylGrow) cylTris = (int)r.tris.size();
    expect(cylTris == 6 * 2 * kRows + 1, "junction uncertified: the chip stays whole");
    bool chipTaken = false;
    for (const Region* t : ts) chipTaken = chipTaken || has(t->tris, j.p.tri[0][8][0]);
    expect(!chipTaken, "junction uncertified: no chip triangle joins the torus");
}

// A chip whose off-torus triangles have no certified neighbour (one vertex 10 tau
// off, inside the chip) stays whole.
static void caseChipOffNoNeighbour() {
    const double tau = 2e-6;
    Part p = quarterRound(3, 12, 10.0 * tau);
    p.mesh.finish();
    p.work.accepted.push_back(band(sectors(p, 0, 7)));
    p.work.accepted.push_back(chip(sectors(p, 8, 13), 100.0 * tau));
    p.work.accepted.push_back(plane(sectors(p, 14, 23)));
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    claimTorusSameSurface(mv, sp, p.work);
    const auto ts = tori(p.work);
    expect(ts.size() == 1 && ts[0]->tris.size() == 8u * 2u * kRows &&
               countOrigin(p.work, Origin::CylGrow) == 1,
           "off no neighbour: the chip stays whole, the band does not grow");
}

// (c) under U-3: the junction chip with a withheld triangle is not entered.
static void caseChipJunctionWithheld() {
    Junction j = junction(0.0);
    MeshView mv = j.p.mesh.view();
    SegmentParams sp{};
    sp.torusWithheldTris = {j.p.tri[5][13][1]};
    claimTorusSameSurface(mv, sp, j.p.work);
    const auto ts = tori(j.p.work);
    int cylTris = -1;
    for (const Region& r : j.p.work.accepted)
        if (r.origin == Origin::CylGrow) cylTris = (int)r.tris.size();
    expect(cylTris == 6 * 2 * kRows + 1, "junction withheld: the chip stays whole");
    bool chipTaken = false;
    for (const Region* t : ts) chipTaken = chipTaken || has(t->tris, j.p.tri[0][8][0]);
    expect(!chipTaken, "junction withheld: the band does not enter the chip");
}

// The stage does nothing without a certified torus (an uncertified torus counts as none).
static void caseNoCert() {
    const double tau = 2e-6;
    Part p = quarterRound();
    p.mesh.finish();
    Region b = band(sectors(p, 0, 7));
    b.maxVertexDev = 3.0 * tau;
    p.work.accepted.push_back(b);
    MeshView mv = p.mesh.view();
    SegmentParams sp{};
    claimTorusSameSurface(mv, sp, p.work);
    expect(tori(p.work).size() == 1 && tori(p.work)[0]->tris.size() == 8u * 2u * kRows,
           "no cert: an uncertified torus does not grow");
}

// Heal-not-revert naming for the U-3 retry.
static void caseFailedNames() {
    Region t;
    t.origin = Origin::TorusBlend;
    t.type = SurfType::Torus;
    t.tris = {1, 2, 3, 4, 5};
    expect(torusFailedTrisOf(t) == t.tris, "names: an ungrown torus names all");
    t.torusGrowTris = {4, 5};
    expect(torusFailedTrisOf(t) == std::vector<int>({4, 5}), "names: a grown band names its growth");
    expect(!torusBandHeld(t, {4, 5}), "held: growth withheld keeps the band");
    expect(torusBandHeld(t, {2}), "held: a band triangle withheld holds it");
    t.torusGrowTris = t.tris;
    expect(torusFailedTrisOf(t) == t.tris, "names: an all-growth piece names all");
    expect(torusBandHeld(t, {5}), "held: a withheld piece triangle holds the piece");
    Region u;
    u.tris = {7, 8};
    expect(torusBandHeld(u, {8}) && !torusBandHeld(u, {9}), "held: ungrown = any triangle withheld");
}

int main() {
    caseFill();
    caseWithheld();
    caseChipPiece();
    caseChipOff();
    caseChipJunction();
    caseChipJunctionUncertified();
    caseChipOffNoNeighbour();
    caseChipJunctionWithheld();
    caseNoCert();
    caseFailedNames();
    if (gFail) std::fprintf(stderr, "relational_torus_grow_unit: %d failure(s)\n", gFail);
    return gFail ? 1 : 0;
}
