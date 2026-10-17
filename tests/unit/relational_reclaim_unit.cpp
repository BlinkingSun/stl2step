// Synthetic relational-reclaim cases: Line mode (D-fc-4.4, D-fc-4.5, the ring donor of
// spike-fc3-lowfold F2, the closed-by-construction wall of D-fc-5.7 O-1) and the
// Torus-mode cap-chord peel. No training mesh.
// SPDX-License-Identifier: MIT

#include "posix_compat.hpp"
#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#define rc_dup _dup
#define rc_dup2 _dup2
#define rc_fileno _fileno
#define rc_close _close
#else
// dup/dup2/close are libc. Declared here so this file does not include <unistd.h>.
// fileno comes from <cstdio>.
extern "C" {
int dup(int);
int dup2(int, int);
int close(int);
}
#define rc_dup dup
#define rc_dup2 dup2
#define rc_fileno fileno
#define rc_close close
#endif

using namespace stl2step::refit;

static int gFail = 0;

static void expect(bool ok, const char* name) {
    std::fprintf(stderr, "%s %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++gFail;
}

static constexpr double kPi = 3.14159265358979323846264338327950288;
static constexpr double kDeg = kPi / 180.0;

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
        mv.weldTol = 0.0;
        mv.sewTol = 0.0;
        mv.quantFloor = q;
        return mv;
    }
};

// A cylinder (axis +Z through the origin) tessellated as nSect sectors of `step`
// from th0, one quad (2 triangles) per sector per level gap. full => the last
// sector closes onto the first column. Returns the triangles of each sector.
struct Grid {
    std::vector<std::vector<int>> col;  // col[j][i]: level j, column i
    std::vector<std::vector<int>> sect; // sect[i]: triangles of sector i
};

static Grid cylGrid(Synth& s, double R, double th0, double step, int nSect, bool full,
                    const std::vector<double>& zs) {
    Grid g;
    const int nCol = full ? nSect : nSect + 1;
    for (double z : zs) {
        std::vector<int> row;
        for (int i = 0; i < nCol; ++i) {
            const double th = th0 + step * i;
            row.push_back(s.add(gp_XYZ(R * std::cos(th), R * std::sin(th), z)));
        }
        g.col.push_back(row);
    }
    g.sect.resize((std::size_t)nSect);
    for (int i = 0; i < nSect; ++i) {
        const int i1 = (i + 1) % nCol;
        for (std::size_t j = 0; j + 1 < zs.size(); ++j) {
            const int a = g.col[j][(std::size_t)i], b = g.col[j][(std::size_t)i1];
            const int c = g.col[j + 1][(std::size_t)i], d = g.col[j + 1][(std::size_t)i1];
            g.sect[(std::size_t)i].push_back(s.addTri(a, b, d));
            g.sect[(std::size_t)i].push_back(s.addTri(a, d, c));
        }
    }
    return g;
}

static std::vector<int> sectors(const Grid& g, int lo, int hi) {
    std::vector<int> out;
    for (int i = lo; i <= hi; ++i)
        out.insert(out.end(), g.sect[(std::size_t)i].begin(), g.sect[(std::size_t)i].end());
    std::sort(out.begin(), out.end());
    return out;
}

static void pushRadial(Synth& s, int v, double dr) {
    gp_XYZ& p = s.pts[(std::size_t)v];
    const double m = std::hypot(p.X(), p.Y());
    p.SetX(p.X() * (m + dr) / m);
    p.SetY(p.Y() * (m + dr) / m);
}

static Region cylRegion(const std::vector<int>& tris, double R, double maxDev) {
    Region r;
    r.id = 0;
    r.type = SurfType::Cylinder;
    r.origin = Origin::CylGrow;
    r.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    r.radius = R;
    r.closed360 = false;
    r.tris = tris;
    r.maxVertexDev = maxDev;
    return r;
}

static Region planeRegion(const std::vector<int>& tris) {
    Region r;
    r.id = 0;
    r.type = SurfType::Plane;
    r.origin = Origin::PlaneGrow;
    r.tris = tris;
    return r;
}

static void addProv(SegmentWork& w, const std::vector<int>& tris, ProvClaim c) {
    Provisional pr;
    pr.claim = c;
    pr.tris = tris;
    w.provisionals.push_back(pr);
}

static void params(SegmentParams& p, DerivedTols& tol) {
    p = SegmentParams{};
    p.maxRefineIters = 8;
    tol = DerivedTols{};
    tol.thetaBin = 0.25 * kDeg;
    tol.thetaCylLo = p.thetaCylLoDeg * kDeg;  // as deriveTols sets it
}

struct Grab {
    int saved = -1;
    bool on = false;
    Grab() {
        ::setenv("STL2STEP_DIAG_REL", "1", 1);
        saved = rc_dup(rc_fileno(stderr));
        FILE* f = std::fopen("rel-reclaim-diag.txt", "w");
        if (!f) return;
        rc_dup2(rc_fileno(f), rc_fileno(stderr));
        std::fclose(f);
        on = true;
    }
    std::string stop() {
        std::fflush(stderr);
        if (on && saved >= 0) rc_dup2(saved, rc_fileno(stderr));
        if (saved >= 0) rc_close(saved);
        saved = -1;
        ::unsetenv("STL2STEP_DIAG_REL");
        std::ifstream in("rel-reclaim-diag.txt");
        std::stringstream ss;
        ss << in.rdbuf();
        in.close();
        std::remove("rel-reclaim-diag.txt");
        return ss.str();
    }
};

static bool has(const std::string& diag, const std::string& key) {
    return diag.find(key) != std::string::npos;
}

static int nCyl(const SegmentWork& w) {
    int n = 0;
    for (const Region& r : w.accepted)
        if (r.type == SurfType::Cylinder) ++n;
    return n;
}

static const Region* firstCyl(const SegmentWork& w) {
    for (const Region& r : w.accepted)
        if (r.type == SurfType::Cylinder) return &r;
    return nullptr;
}

static std::string run(Synth& s, SegmentWork& w, RelRevolveMode mode = RelRevolveMode::Line) {
    const MeshView mv = s.view();
    SegmentParams p;
    DerivedTols tol;
    params(p, tol);
    Grab g;
    reclaimRelational(mv, p, tol, w, mode);
    const std::string out = g.stop();
    if (std::getenv("RECLAIM_UNIT_ECHO")) std::fprintf(stderr, "%s", out.c_str());
    return out;
}

// ---- D-fc-4.4 same-surface closure ------------------------------------------------

// Two certified partials of one R3 surface with an unclaimed gap between them; one
// gap vertex sits 50 tau inside the surface (an intersection-curve vertex), so the
// triangles on it are chord-bounded cover, the rest of the gap is I.
struct Closure {
    Synth s;
    SegmentWork w;
    Grid g;
};

static Closure closurePart(bool gapClaimed, bool secondUncertified) {
    Closure c;
    const double tau = 2 * c.s.q;
    c.g = cylGrid(c.s, 3.0, -90 * kDeg, 10 * kDeg, 18, false, {0.0, 5.0, 10.0});
    pushRadial(c.s, c.g.col[1][9], -50 * tau);
    c.s.finish();
    c.w.accepted.push_back(cylRegion(sectors(c.g, 0, 6), 3.0, 0.0));
    c.w.accepted.push_back(cylRegion(sectors(c.g, 11, 17), 3.0, secondUncertified ? 2 * tau : 0.0));
    if (gapClaimed) c.w.accepted.push_back(planeRegion(sectors(c.g, 7, 10)));
    else addProv(c.w, sectors(c.g, 7, 10), ProvClaim::Unclaimed);
    return c;
}

static void caseClosureAccept() {
    Closure c = closurePart(false, false);
    const double tau = 2 * c.s.q;
    const std::string d = run(c.s, c.w);
    const Region* r = firstCyl(c.w);
    expect(c.w.accepted.size() == 1 && r && r->tris.size() == 72, "closure_one_region_all_triangles");
    expect(has(d, "stage=reclaim mode=closure") && has(d, "why=accept"), "closure_diag_accept");
    expect(r && std::abs(r->radius - 3.0) <= tau, "closure_radius_within_tau");
    // The cover triangles carry the 50 tau vertex; the certified set stays at tau.
    expect(r && r->maxVertexDev > tau && r->maxVertexDev <= 51 * tau, "closure_cover_dev_measured");
    expect(r && !r->closed360 && r->dVolPredicted != 0.0, "closure_fields_partial_with_dvol");
    expect(c.w.provisionals[0].claim == ProvClaim::ConsumedCylinder, "closure_gap_peeled");
}

static void caseClosureNotJoined() {
    Closure c = closurePart(true, false);
    const std::string d = run(c.s, c.w);
    expect(c.w.accepted.size() == 3 && nCyl(c.w) == 2, "closure_not_joined_unchanged");
    expect(has(d, "why=not-joined"), "closure_diag_not_joined");
}

static void caseClosureUncertified() {
    Closure c = closurePart(false, true);
    const std::string d = run(c.s, c.w);
    expect(c.w.accepted.size() == 2 && c.w.provisionals[0].claim == ProvClaim::Unclaimed,
           "closure_uncertified_member_ignored");
    expect(!has(d, "mode=closure"), "closure_uncertified_no_group");
}

// Two halves of one full ring merge to a band with two boundary loops: deferred.
static void caseClosureNotDisk() {
    Synth s;
    SegmentWork w;
    const Grid g = cylGrid(s, 3.0, 0.0, 10 * kDeg, 36, true, {0.0, 5.0});
    s.finish();
    w.accepted.push_back(cylRegion(sectors(g, 0, 17), 3.0, 0.0));
    w.accepted.push_back(cylRegion(sectors(g, 18, 35), 3.0, 0.0));
    const std::string d = run(s, w);
    expect(w.accepted.size() == 2 && w.accepted[0].tris.size() == 36, "closure_not_disk_unchanged");
    expect(has(d, "why=domain-not-disk"), "closure_diag_domain_not_disk");
}

// Stacked coaxial walls (no strict axial overlap) and a different radius never join.
static void caseClosureStackedAndRadius() {
    {
        Synth s;
        SegmentWork w;
        const Grid a = cylGrid(s, 3.0, 0.0, 10 * kDeg, 9, false, {0.0, 5.0});
        const Grid b = cylGrid(s, 3.0, 0.0, 10 * kDeg, 9, false, {5.0, 10.0});
        s.finish();
        w.accepted.push_back(cylRegion(sectors(a, 0, 8), 3.0, 0.0));
        w.accepted.push_back(cylRegion(sectors(b, 0, 8), 3.0, 0.0));
        const std::string d = run(s, w);
        expect(w.accepted.size() == 2 && !has(d, "mode=closure"), "closure_stacked_apart");
    }
    {
        Synth s;
        SegmentWork w;
        const double tau = 2 * s.q;
        const Grid a = cylGrid(s, 3.0, 0.0, 10 * kDeg, 9, false, {0.0, 5.0});
        const Grid b = cylGrid(s, 3.0 + 10 * tau, 180 * kDeg, 10 * kDeg, 9, false, {0.0, 5.0});
        s.finish();
        w.accepted.push_back(cylRegion(sectors(a, 0, 8), 3.0, 0.0));
        w.accepted.push_back(cylRegion(sectors(b, 0, 8), 3.0 + 10 * tau, 0.0));
        const std::string d = run(s, w);
        expect(w.accepted.size() == 2 && !has(d, "mode=closure"), "closure_radius_apart");
    }
}

// ---- D-fc-4.5 narrow donor heal ---------------------------------------------------

// R3.5 half-ring donor (sectors 0..17) with one interior top vertex `off` inside the
// surface (its triangles are the donor's cover); free sectors 18..23 on the surface;
// sectors 24..29 a certified partner reached only through the gap (or a plane);
// sectors 30..35 belong to a plane region.
struct Heal {
    Synth s;
    SegmentWork w;
    Grid g;
};

static Heal healPart(double off, bool partner) {
    Heal h;
    h.g = cylGrid(h.s, 3.5, 0.0, 10 * kDeg, 36, true, {0.0, 5.0});
    if (off != 0.0) pushRadial(h.s, h.g.col[1][9], -off);
    h.s.finish();
    h.w.accepted.push_back(cylRegion(sectors(h.g, 0, 17), 3.5, off));
    if (partner) h.w.accepted.push_back(cylRegion(sectors(h.g, 24, 29), 3.5, 0.0));
    else h.w.accepted.push_back(planeRegion(sectors(h.g, 24, 29)));
    h.w.accepted.push_back(planeRegion(sectors(h.g, 30, 35)));
    addProv(h.w, sectors(h.g, 18, 23), ProvClaim::Unclaimed);
    return h;
}

static void caseHealAccept() {
    const double tau = 2e-6;
    Heal h = healPart(3 * tau, true);
    const std::string d = run(h.s, h.w);
    const Region* r = firstCyl(h.w);
    expect(nCyl(h.w) == 1 && h.w.accepted.size() == 2, "heal_partner_merged");
    expect(r && r->tris.size() == 60, "heal_grows_over_gap");
    expect(r && std::abs(r->radius - 3.5) <= tau, "heal_radius_within_tau");
    expect(r && r->maxVertexDev > tau, "heal_cover_dev_kept");
    expect(r && !r->closed360 && has(d, "annulus=0"), "heal_partial_not_closed");
    expect(has(d, "mode=heal") && has(d, "why=accept"), "heal_diag_accept");
    expect(h.w.provisionals[0].claim == ProvClaim::ConsumedCylinder, "heal_gap_peeled");
}

static void caseHealNoPartner() {
    const double tau = 2e-6;
    Heal h = healPart(3 * tau, false);
    run(h.s, h.w);
    const Region* r = firstCyl(h.w);
    expect(h.w.accepted.size() == 3 && r && r->tris.size() == 48, "heal_grows_without_partner");
}

static void caseHealNarrowTrigger() {
    Heal h = healPart(0.0, false);
    const std::string d = run(h.s, h.w);
    expect(h.w.accepted.size() == 3 && h.w.accepted[0].tris.size() == 36 &&
               h.w.provisionals[0].claim == ProvClaim::Unclaimed,
           "heal_certified_region_untouched");
    expect(!has(d, "mode=heal"), "heal_certified_no_band");
}

static void caseHealNotCover() {
    // 2x the 10 degree chord sagitta: no chord bound admits it.
    const double sag = 3.5 * (1 - std::cos(5 * kDeg));
    Heal h = healPart(2 * sag, false);
    const std::string d = run(h.s, h.w);
    expect(h.w.accepted.size() == 3 && h.w.accepted[0].tris.size() == 36 &&
               h.w.provisionals[0].claim == ProvClaim::Unclaimed,
           "heal_off_surface_refused");
    expect(has(d, "mode=heal") && !has(d, "why=accept"), "heal_off_surface_diag_refusal");
}

// ---- D-fc-5.7 O-1: a healed wall that is an annulus is closed by construction ------

static void caseWindingAndAnnulus() {
    std::vector<double> ring;
    for (int i = 0; i < 8; ++i) ring.push_back(-kPi + 2 * kPi * i / 8.0 + 0.1);
    expect(windingNumber(ring) == 1, "winding_ring_once");
    std::vector<double> back(ring.rbegin(), ring.rend());
    expect(windingNumber(back) == -1, "winding_ring_reversed");
    expect(windingNumber({0.1, 0.4, 0.3, 0.2}) == 0, "winding_window_not_about_axis");
    expect(windingNumber({3.0, -3.0, -1.0, 1.0}) == 1, "winding_across_branch_cut");
    expect(windingNumber({0.0, 1.0}) == 0, "winding_two_points_none");
    expect(annulusAboutAxis(1, {1, -1}) && annulusAboutAxis(1, {-1, -1}), "annulus_two_cycles_once");
    expect(!annulusAboutAxis(2, {1, -1}), "annulus_two_components_refused");
    expect(!annulusAboutAxis(1, {1}), "annulus_one_cycle_is_a_disk");
    expect(!annulusAboutAxis(1, {1, 0}), "annulus_window_cycle_refused");
    expect(!annulusAboutAxis(1, {2, -1}), "annulus_double_wind_refused");
    expect(!annulusAboutAxis(1, {1, -1, 0}), "annulus_three_cycles_refused");
}

// An R4 ring (z 0..6) tessellated unevenly: one 20 degree sector, then 5 degree
// sectors round to 360 (69 sectors). D2's band rule reads it open (largest vertex gap
// 20 degrees > 1.5 x 360/69). The donor (sectors 0..34) carries one top vertex 3 tau
// inside the surface (its cover); sectors 35..lastFree are free, the rest a plane.
struct Uneven {
    Synth s;
    SegmentWork w;
    std::vector<std::vector<int>> sect;
};

static Uneven unevenRing(int lastFree) {
    Uneven u;
    const double R = 4.0, tau = 2 * u.s.q;
    std::vector<double> az{0.0};
    for (double a = 20.0; a < 359.0; a += 5.0) az.push_back(a * kDeg);
    const int n = (int)az.size();
    std::vector<int> bot, top;
    for (double a : az) bot.push_back(u.s.add(gp_XYZ(R * std::cos(a), R * std::sin(a), 0.0)));
    for (double a : az) top.push_back(u.s.add(gp_XYZ(R * std::cos(a), R * std::sin(a), 6.0)));
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        u.sect.push_back({u.s.addTri(bot[i], bot[j], top[j]), u.s.addTri(bot[i], top[j], top[i])});
    }
    pushRadial(u.s, top[17], -3 * tau);
    u.s.finish();
    auto span = [&](int lo, int hi) {
        std::vector<int> out;
        for (int i = lo; i <= hi; ++i) out.insert(out.end(), u.sect[i].begin(), u.sect[i].end());
        std::sort(out.begin(), out.end());
        return out;
    };
    u.w.accepted.push_back(cylRegion(span(0, 34), R, 3 * tau));
    addProv(u.w, span(35, lastFree), ProvClaim::Unclaimed);
    if (lastFree + 1 < n) u.w.accepted.push_back(planeRegion(span(lastFree + 1, n - 1)));
    return u;
}

static void caseHealClosedByConstruction() {
    Uneven u = unevenRing(68);
    const std::string d = run(u.s, u.w);
    const Region* r = firstCyl(u.w);
    expect(nCyl(u.w) == 1 && r && r->tris.size() == 138, "heal_ring_whole_wall");
    expect(r && r->closed360 && r->uMin == 0.0 && r->uMax > 6.28, "heal_ring_closed_by_construction");
    expect(has(d, "closed360=1 annulus=1"), "heal_ring_diag_annulus");
}

static void caseHealOpenStaysOpen() {
    Uneven u = unevenRing(67);
    const std::string d = run(u.s, u.w);
    const Region* r = firstCyl(u.w);
    expect(r && r->tris.size() == 136 && !r->closed360, "heal_disk_stays_open");
    expect(has(d, "closed360=0 annulus=0"), "heal_disk_diag_not_annulus");
}

// ---- Ring donor (spike-fc3-lowfold F2) --------------------------------------------

// A plate at z = 0 whose hole is a circular loop of nSect vertices at R (one Unclaimed
// chart, outer rim at R + 5) over a peg wall at the levels zs (down to z = -20 by
// default), one Unclaimed chart per wall sector. The wall hangs under the first nWall
// sectors only (0: all of them, a closed wall). capFan closes the peg's end with a fan
// on its rim vertices only (every vertex on the cylinder, one axial level).
struct Ring {
    Synth s;
    SegmentWork w;
    Grid peg;
    std::vector<int> plate, cap;
};

static Ring ringPart(double R, int nSect, bool capFan, int nWall = 0,
                     const std::vector<double>& zs = {0.0, -20.0}) {
    Ring r;
    const double step = 2 * kPi / nSect;
    const bool closed = nWall <= 0 || nWall >= nSect;
    const int nW = closed ? nSect : nWall;
    r.peg = cylGrid(r.s, R, 0.0, step, nW, closed, zs);
    std::vector<int> hole = r.peg.col[0];
    for (int i = (int)hole.size(); i < nSect; ++i)
        hole.push_back(r.s.add(gp_XYZ(R * std::cos(step * i), R * std::sin(step * i), 0.0)));
    std::vector<int> outer;
    for (int i = 0; i < nSect; ++i)
        outer.push_back(r.s.add(gp_XYZ((R + 5) * std::cos(step * i), (R + 5) * std::sin(step * i), 0.0)));
    for (int i = 0; i < nSect; ++i) {
        const int i1 = (i + 1) % nSect;
        const int h0 = hole[(std::size_t)i], h1 = hole[(std::size_t)i1];
        r.plate.push_back(r.s.addTri(h0, outer[(std::size_t)i1], h1));
        r.plate.push_back(r.s.addTri(h0, outer[(std::size_t)i], outer[(std::size_t)i1]));
    }
    if (capFan) {
        const std::vector<int>& b = r.peg.col.back();
        for (int i = 1; i + 1 < nSect; ++i)
            r.cap.push_back(r.s.addTri(b[0], b[(std::size_t)i + 1], b[(std::size_t)i]));
    }
    r.s.finish();
    addProv(r.w, r.plate, ProvClaim::Unclaimed);
    for (int i = 0; i < nW; ++i) addProv(r.w, sectors(r.peg, i, i), ProvClaim::Unclaimed);
    if (capFan) addProv(r.w, r.cap, ProvClaim::Unclaimed);
    return r;
}

// The why= of every DIAG_REL_RING line on a loop of n vertices, comma-joined in order.
static std::string ringWhys(const std::string& diag, int n) {
    std::string out;
    std::istringstream in(diag);
    const std::string key = " loop=" + std::to_string(n) + " ";
    for (std::string ln; std::getline(in, ln);) {
        if (ln.rfind("DIAG_REL_RING", 0) != 0 || ln.find(key) == std::string::npos) continue;
        const std::size_t w = ln.find("why=");
        if (w == std::string::npos) continue;
        if (!out.empty()) out += ",";
        out += ln.substr(w + 4);
    }
    return out;
}

static void caseRingAccept() {
    const double R = 7.5;
    Ring r = ringPart(R, 120, false);
    const double tau = 2 * r.s.q;
    const std::string d = run(r.s, r.w);
    const Region* c = firstCyl(r.w);
    expect(nCyl(r.w) == 1 && c && c->tris.size() == 240, "ring_peg_one_closed_cylinder");
    expect(c && c->closed360 && c->origin == Origin::CylGrow, "ring_peg_closed360_cylgrow");
    expect(c && std::abs(c->radius - R) <= tau && c->maxVertexDev <= tau, "ring_peg_certified_tau");
    expect(has(d, "mode=ring") && has(d, "why=accept"), "ring_diag_accept");
    expect(r.w.provisionals[0].claim == ProvClaim::Unclaimed && r.w.provisionals[0].tris.size() == 240 &&
               r.w.provisionals[1].claim == ProvClaim::ConsumedCylinder,
           "ring_anchor_kept_wall_charts_consumed");
}

// A square hole: its corners are concyclic, its walls fold 90 degrees (a prism).
static void caseRingPrismRefused() {
    Ring r = ringPart(5.0, 4, false);
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 0 && has(d, "why=not-blend"), "ring_prism_refused_by_fold_scope");
}

// A fan on the rim vertices: every vertex on the cylinder, one axial level.
static void caseRingCapChordExcluded() {
    Ring r = ringPart(7.5, 120, true);
    run(r.s, r.w);
    const Region* c = firstCyl(r.w);
    expect(nCyl(r.w) == 1 && c && c->tris.size() == 240, "ring_cap_chords_not_in_wall");
    expect(r.w.provisionals.back().claim == ProvClaim::Unclaimed &&
               r.w.provisionals.back().tris.size() == r.cap.size(),
           "ring_cap_chart_untouched");
}

// A rectangle's corners are concyclic too; triangles hanging from each side, each on
// the circumscribed cylinder, touch the loop only one by one.
static void caseRingNotJoined() {
    Synth s;
    SegmentWork w;
    const gp_XYZ c[4] = {gp_XYZ(15, 5, 0), gp_XYZ(-15, 5, 0), gp_XYZ(-15, -5, 0), gp_XYZ(15, -5, 0)};
    int top[4], bot[4];
    for (int i = 0; i < 4; ++i) top[i] = s.add(c[i]);
    for (int i = 0; i < 4; ++i) bot[i] = s.add(c[i] + gp_XYZ(0, 0, -5));
    std::vector<int> plate{s.addTri(top[0], top[1], top[2]), s.addTri(top[0], top[2], top[3])};
    std::vector<int> hang;
    for (int i = 0; i < 4; ++i) hang.push_back(s.addTri(top[i], bot[i], top[(i + 1) % 4]));
    s.finish();
    addProv(w, plate, ProvClaim::Unclaimed);
    for (int t : hang) addProv(w, {t}, ProvClaim::Unclaimed);
    const std::string d = run(s, w);
    expect(nCyl(w) == 0 && has(d, "why=not-joined"), "ring_disconnected_band_refused");
}

// The loop's chart is already another stage's: no anchor, nothing claimed.
static void caseRingClaimedAnchor() {
    Ring r = ringPart(7.5, 120, false);
    r.w.provisionals[0].claim = ProvClaim::ConsumedCylinder;
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 0 && has(d, "why=no-anchor"), "ring_claimed_chart_is_no_anchor");
}

// The loop is only a hypothesis source. One hole vertex 3 tau off the circle, or off
// the plane: the loop names no axis, no band is grown and no ring line is printed.
static void caseRingLoopNotCircle() {
    for (int k = 0; k < 2; ++k) {
        Ring r = ringPart(7.5, 120, false);
        const double tau = 2 * r.s.q;
        const int v = r.peg.col[0][0];
        if (k == 0) pushRadial(r.s, v, 3 * tau);
        else r.s.pts[(std::size_t)v].SetZ(3 * tau);
        const std::string d = run(r.s, r.w);
        expect(nCyl(r.w) == 0 && ringWhys(d, 120).empty(),
               k == 0 ? "ring_loop_not_circular_no_hypothesis" : "ring_loop_not_planar_no_hypothesis");
    }
}

// Every wall vertex is within tau of the loop's circle, but no cylinder fits them all
// at tau: the far rim stands 0.95 tau out and one of its vertices 0.95 tau in.
static void caseRingVertexOff() {
    Ring r = ringPart(7.5, 120, false);
    const double tau = 2 * r.s.q;
    for (std::size_t i = 0; i < r.peg.col[1].size(); ++i)
        pushRadial(r.s, r.peg.col[1][i], (i == 0 ? -0.95 : 0.95) * tau);
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 0 && ringWhys(d, 120) == "vertex-off", "ring_band_refit_off_tau_refused");
}

// A third of the wall is already a certified cylinder: the band reaches it on its own
// surface, so the wall is that region's surface (one face per surface).
static void caseRingSameSurface() {
    Ring r = ringPart(7.5, 120, false);
    r.w.accepted.push_back(cylRegion(sectors(r.peg, 0, 39), 7.5, 0.0));
    for (std::size_t i = 1; i <= 40; ++i) r.w.provisionals[i].claim = ProvClaim::ConsumedCylinder;
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 1 && r.w.accepted[0].tris.size() == 80 && ringWhys(d, 120) == "same-surface",
           "ring_band_on_claimed_surface_refused");
    expect(r.w.provisionals[41].claim == ProvClaim::Unclaimed, "ring_same_surface_claims_nothing");
}

// A wall under two sectors of the hole: three stations cannot over-determine the
// circle (minStations(kParamsCircle)).
static void caseRingUnderdetermined() {
    Ring r = ringPart(7.5, 120, false, 2);
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 0 && ringWhys(d, 120) == "underdetermined", "ring_short_wall_underdetermined");
}

// The wall levels of the dent cases: a dent on the middle level is an interior hole.
static const std::vector<double> kDentLevels{0.0, -5.0, -10.0, -15.0, -20.0};

// A dent 5 tau into a closed wall: the triangles on the dented vertex leave the band,
// whose hole is a third boundary loop.
static void caseRingNotAnnulus() {
    Ring r = ringPart(7.5, 120, false, 0, kDentLevels);
    const double tau = 2 * r.s.q;
    pushRadial(r.s, r.peg.col[2][60], -5 * tau);
    const std::string d = run(r.s, r.w);
    expect(nCyl(r.w) == 0 && ringWhys(d, 120) == "domain-not-annulus", "ring_band_with_hole_refused");
}

// A C-shaped wall (100 of 120 sectors) is a disk: one boundary loop. With a dent it
// has the two loops of an annulus but does not close round the axis.
static void caseRingOpenWall() {
    {
        Ring r = ringPart(7.5, 120, false, 100, kDentLevels);
        const std::string d = run(r.s, r.w);
        expect(nCyl(r.w) == 0 && ringWhys(d, 120) == "domain-not-annulus", "ring_open_wall_is_a_disk");
    }
    {
        Ring r = ringPart(7.5, 120, false, 100, kDentLevels);
        const double tau = 2 * r.s.q;
        pushRadial(r.s, r.peg.col[2][50], -5 * tau);
        const std::string d = run(r.s, r.w);
        expect(nCyl(r.w) == 0 && ringWhys(d, 120) == "degenerate", "ring_open_annulus_not_closed360");
    }
}

// ---- Torus mode: cap chords leave the torus (F10) ---------------------------------

// A closed quarter torus (v from -90 to 0 degrees) and a fan cap on its v = -90 degree
// ring. The TorusBlend region and its charts hold the band and / or the cap.
struct TorusCap {
    Synth s;
    SegmentWork w;
    std::vector<int> band, cap;
};

static TorusCap torusCapPart(bool withBand, bool withCap) {
    TorusCap tc;
    Synth& s = tc.s;
    SegmentWork& w = tc.w;
    std::vector<int>& band = tc.band;
    std::vector<int>& cap = tc.cap;
    const double Rm = 5.5, rm = 1.5;
    const int nU = 72;
    const double vs[4] = {-90 * kDeg, -60 * kDeg, -30 * kDeg, 0.0};
    std::vector<std::vector<int>> ring(4);
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < nU; ++i) {
            const double u = 2 * kPi * i / nU, rho = Rm + rm * std::cos(vs[j]);
            ring[(std::size_t)j].push_back(
                s.add(gp_XYZ(rho * std::cos(u), rho * std::sin(u), rm * std::sin(vs[j]))));
        }
    for (int j = 0; j + 1 < 4; ++j)
        for (int i = 0; i < nU; ++i) {
            const int i1 = (i + 1) % nU;
            const auto& a = ring[(std::size_t)j];
            const auto& b = ring[(std::size_t)j + 1];
            band.push_back(s.addTri(a[(std::size_t)i], a[(std::size_t)i1], b[(std::size_t)i1]));
            band.push_back(s.addTri(a[(std::size_t)i], b[(std::size_t)i1], b[(std::size_t)i]));
        }
    for (int i = 1; i + 1 < nU; ++i)
        cap.push_back(s.addTri(ring[0][0], ring[0][(std::size_t)i + 1], ring[0][(std::size_t)i]));
    s.finish();
    Region T;
    T.type = SurfType::Torus;
    T.origin = Origin::TorusBlend;
    T.ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    T.radius = Rm;
    T.radius2 = rm;
    T.closed360 = true;
    T.uMin = 0;
    T.uMax = 2 * kPi;
    T.vMin = -90 * kDeg;
    T.vMax = 0.0;
    T.nSides = nU;
    T.outwardNormal = true;
    if (withBand) T.tris = band;
    if (withCap) T.tris.insert(T.tris.end(), cap.begin(), cap.end());
    w.accepted.push_back(T);
    if (withBand) addProv(w, band, ProvClaim::ConsumedTorus);
    if (withCap) addProv(w, cap, ProvClaim::ConsumedTorus);
    return tc;
}

static void caseTorusCapPeel() {
    TorusCap tc = torusCapPart(true, true);
    SegmentWork& w = tc.w;
    const std::string d = run(tc.s, w, RelRevolveMode::Torus);
    expect(w.accepted.size() == 1 && w.accepted[0].tris.size() == tc.band.size(), "torus_cap_chords_peeled");
    expect(w.provisionals[1].claim == ProvClaim::Unclaimed && w.provisionals[1].tris.size() == tc.cap.size(),
           "torus_cap_chart_unclaimed");
    expect(w.provisionals[0].claim == ProvClaim::ConsumedTorus, "torus_band_chart_kept");
    expect(has(d, "DIAG_REL_TORUS_CAP"), "torus_cap_diag");
}

// F10 peels cap chords only off a torus that keeps a band: a region of cap chords
// alone, or of band triangles alone, stays as it stood.
static void caseTorusCapGuards() {
    for (int k = 0; k < 2; ++k) {
        TorusCap tc = torusCapPart(k == 1, k == 0);
        const std::size_t n = (k == 0 ? tc.cap : tc.band).size();
        const std::string d = run(tc.s, tc.w, RelRevolveMode::Torus);
        expect(tc.w.accepted.size() == 1 && tc.w.accepted[0].tris.size() == n &&
                   tc.w.provisionals[0].claim == ProvClaim::ConsumedTorus && !has(d, "DIAG_REL_TORUS_CAP"),
               k == 0 ? "torus_all_cap_region_kept" : "torus_no_cap_region_kept");
    }
}

// ---- gates ------------------------------------------------------------------------

static void caseNoQuantum() {
    Closure c = closurePart(false, false);
    c.s.q = 0.0;
    const std::string d = run(c.s, c.w);
    expect(c.w.accepted.size() == 2 && !has(d, "stage=reclaim"), "q_zero_admits_nothing");
}

static void caseTorusModeIsNotLine() {
    Closure c = closurePart(false, false);
    const std::string d = run(c.s, c.w, RelRevolveMode::Torus);
    expect(c.w.accepted.size() == 2 && !has(d, "mode=closure"), "torus_mode_runs_no_line_pass");
}

int main() {
    caseClosureAccept();
    caseClosureNotJoined();
    caseClosureUncertified();
    caseClosureNotDisk();
    caseClosureStackedAndRadius();
    caseHealAccept();
    caseHealNoPartner();
    caseHealNarrowTrigger();
    caseHealNotCover();
    caseWindingAndAnnulus();
    caseHealClosedByConstruction();
    caseHealOpenStaysOpen();
    caseRingAccept();
    caseRingPrismRefused();
    caseRingCapChordExcluded();
    caseRingNotJoined();
    caseRingClaimedAnchor();
    caseRingLoopNotCircle();
    caseRingVertexOff();
    caseRingSameSurface();
    caseRingUnderdetermined();
    caseRingNotAnnulus();
    caseRingOpenWall();
    caseTorusCapPeel();
    caseTorusCapGuards();
    caseNoQuantum();
    caseTorusModeIsNotLine();
    if (gFail) std::fprintf(stderr, "relational_reclaim_unit: %d FAIL\n", gFail);
    return gFail ? 1 : 0;
}
