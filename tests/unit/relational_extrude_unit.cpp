// Pins the shared closed forms (D-fc-1.4, D-fc-1.6) and the extrude certificate
// on synthetic in-memory meshes. No STL paths.
// SPDX-License-Identifier: MIT

#include "refit_relational.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

static int gFail = 0;

static constexpr double kPi = 3.14159265358979323846264338327950288;

static void expect(bool ok, const char* name) {
    if (ok) {
        std::fprintf(stderr, "PASS %s\n", name);
        return;
    }
    std::fprintf(stderr, "FAIL %s\n", name);
    ++gFail;
}

static std::string captureStderr(const std::function<void()>& fn) {
    std::fflush(stderr);
    const int saved = dup(fileno(stderr));
    char path[] = "/tmp/rel-extrude-unit-XXXXXX";
    const int fd = mkstemp(path);
    dup2(fd, fileno(stderr));
    fn();
    std::fflush(stderr);
    dup2(saved, fileno(stderr));
    close(saved);
    std::string out;
    FILE* f = fdopen(fd, "r");
    if (f) {
        std::rewind(f);
        char buf[1024];
        while (std::fgets(buf, sizeof buf, f)) out += buf;
        std::fclose(f);
    } else {
        close(fd);
    }
    std::remove(path);
    return out;
}

static std::string lineWithWhy(const std::string& diag, const char* why) {
    const std::string key = std::string("why=") + why;
    std::istringstream in(diag);
    std::string line, last;
    while (std::getline(in, line)) {
        if (line.find("DIAG_REL_BAND") != std::string::npos && line.find(key) != std::string::npos)
            last = line;
    }
    return last;
}

static bool hasWhy(const std::string& diag, const char* why) {
    return !lineWithWhy(diag, why).empty();
}

struct Synth {
    struct Tri3 {
        int v[3];
    };
    std::vector<gp_XYZ> pts;
    std::vector<Tri3> tris;
    std::vector<std::array<int, 3>> triEdges;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> compTris;
    std::vector<int> compVtx;
    std::map<std::uint64_t, int> edgeOf;

    int addV(const gp_XYZ& p) {
        pts.push_back(p);
        return static_cast<int>(pts.size()) - 1;
    }

    int edgeId(int a, int b) {
        const int lo = std::min(a, b);
        const int hi = std::max(a, b);
        const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32) |
                                  static_cast<std::uint32_t>(hi);
        const auto it = edgeOf.find(key);
        if (it != edgeOf.end()) return it->second;
        const int id = static_cast<int>(edges.size());
        edges.push_back({lo, hi});
        edgeOf.emplace(key, id);
        return id;
    }

    int tri(int a, int b, int c) {
        const int t = static_cast<int>(tris.size());
        tris.push_back(Tri3{{a, b, c}});
        triEdges.push_back(
            std::array<int, 3>{{edgeId(a, b), edgeId(b, c), edgeId(c, a)}});
        return t;
    }

    void addQuad(int a, int b, int c, int d, std::vector<int>& out) {
        // a-c is one generator, b-d the other; a-b and c-d are the ends.
        out.push_back(tri(a, c, d));
        out.push_back(tri(a, d, b));
    }

    void seal() {
        compTris.resize(tris.size());
        compVtx.resize(pts.size());
        for (std::size_t i = 0; i < tris.size(); ++i) compTris[i] = static_cast<int>(i);
        for (std::size_t i = 0; i < pts.size(); ++i) compVtx[i] = static_cast<int>(i);
    }

    stl2step::refit::MeshView view(double q) const {
        stl2step::refit::MeshView mv{};
        mv.pts = pts.data();
        mv.tris = reinterpret_cast<const int(*)[3]>(tris.data());
        mv.compTris = compTris.empty() ? nullptr : compTris.data();
        mv.compVtx = compVtx.empty() ? nullptr : compVtx.data();
        mv.compEdges = edges.data();
        mv.triEdges = triEdges.data();
        mv.nTri = tris.size();
        mv.nVtx = pts.size();
        mv.nEdge = edges.size();
        mv.diag = 100.0;
        mv.quantFloor = q;
        return mv;
    }
};

static void addStrip(stl2step::refit::SegmentWork& w, const std::vector<int>& tris) {
    using namespace stl2step::refit;
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = tris;
    pr.area = 1.0;
    pr.maxVertexDev = 0.0;
    pr.plane = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    w.provisionals.push_back(std::move(pr));
}

static void addAnchor(stl2step::refit::SegmentWork& w, const std::vector<int>& tris,
                      const gp_XYZ& n, const gp_XYZ& p0, double maxDev) {
    using namespace stl2step::refit;
    Provisional pr;
    pr.claim = ProvClaim::Unclaimed;
    pr.tris = tris;
    pr.area = 1000.0;
    pr.maxVertexDev = maxDev;
    pr.plane = gp_Ax3(gp_Pnt(p0.X(), p0.Y(), p0.Z()), gp_Dir(n.X(), n.Y(), n.Z()));
    w.provisionals.push_back(std::move(pr));
}

struct Run {
    stl2step::refit::SegmentWork work;
    std::string diag;
    double tau = 0;
};

static Run runClaim(Synth& s, stl2step::refit::SegmentWork work, double q, int threads) {
    using namespace stl2step::refit;
    s.seal();
    SegmentParams p;
    p.requestedThreads = threads;
    p.maxRefineIters = 8;
    DerivedTols tol;
    tol.epsPlane = 1.0;
    tol.thetaBin = 0.25 * kPi / 180.0;
    const auto mv = s.view(q);
    Run r;
    r.tau = (q > 0.0) ? 2.0 * q : 0.0;
    r.diag = captureStderr([&] { claimRelExtrude(mv, p, tol, work); });
    r.work = std::move(work);
    return r;
}

static void dumpFail(const char* name, const Run& r) {
    std::fprintf(stderr, "---- %s accepted=%d ----\n%s", name,
                 static_cast<int>(r.work.accepted.size()), r.diag.c_str());
}

// Concave quarter fillet along +Y. Centre (R, y, R), planes x=0 and z=0.
struct FilletMesh {
    Synth synth;
    stl2step::refit::SegmentWork work;
    std::vector<int> bandTris;
    int wallA = -1;
    int wallB = -1;
};

static FilletMesh makeFillet(double R, double y0, double y1, double ang0, double ang1, int nStrips,
                             double anchorDev) {
    FilletMesh m;
    const int nAng = nStrips + 1;
    std::vector<int> lo((std::size_t)nAng), hi((std::size_t)nAng);
    for (int i = 0; i < nAng; ++i) {
        const double th = ang0 + (ang1 - ang0) * (double)i / (double)nStrips;
        const double c = std::cos(th);
        const double s = std::sin(th);
        lo[(std::size_t)i] = m.synth.addV(gp_XYZ(R * (1.0 - c), y0, R * (1.0 - s)));
        hi[(std::size_t)i] = m.synth.addV(gp_XYZ(R * (1.0 - c), y1, R * (1.0 - s)));
    }
    for (int i = 0; i < nStrips; ++i) {
        std::vector<int> tris;
        m.synth.addQuad(lo[(std::size_t)i], hi[(std::size_t)i], lo[(std::size_t)i + 1],
                        hi[(std::size_t)i + 1], tris);
        for (int t : tris) m.bandTris.push_back(t);
        addStrip(m.work, tris);
    }
    const gp_XYZ midA(0.0, 0.5 * (y0 + y1), R + 30.0);
    const gp_XYZ midB(R + 30.0, 0.5 * (y0 + y1), 0.0);
    m.wallA = m.synth.tri(lo.front(), hi.front(), m.synth.addV(midA));
    m.wallB = m.synth.tri(lo.back(), hi.back(), m.synth.addV(midB));
    addAnchor(m.work, {m.wallA}, gp_XYZ(-1, 0, 0), gp_XYZ(0, 0, 0), anchorDev);
    addAnchor(m.work, {m.wallB}, gp_XYZ(0, 0, -1), gp_XYZ(0, 0, 0), anchorDev);
    return m;
}

static void pins() {
    using namespace stl2step::refit;
    MeshView qHalf{};
    qHalf.quantFloor = 0.5;
    MeshView qZero{};
    qZero.quantFloor = 0.0;
    const double tauHalf = tauOf(qHalf);
    const double tauZero = tauOf(qZero);
    std::fprintf(stderr, "tauOf(0.5)=%.17g tauOf(0)=%.17g\n", tauHalf, tauZero);
    expect(tauHalf == 1.0, "tauOf(0.5) = 1");
    expect(tauZero == 0.0, "tauOf(0) = 0");

    const double a1 = anchorTol(1.0, 3.0);
    const double a2 = anchorTol(4.0, 1.0);
    std::fprintf(stderr, "anchorTol(1,3)=%.17g anchorTol(4,1)=%.17g\n", a1, a2);
    expect(a1 == 3.0, "anchorTol(1, 3) = 3");
    expect(a2 == 4.0, "anchorTol(4, 1) = 4");

    const double s1 = spanTol(2.0, 8.0, 1.0);
    const double s2 = spanTol(0.5, 0.5, 1.0);
    std::fprintf(stderr, "spanTol(2,8,1)=%.17g spanTol(0.5,0.5,1)=%.17g\n", s1, s2);
    expect(s1 == 6.0, "spanTol(2, 8, 1) = 6");
    expect(s2 == 2.0, "spanTol(0.5, 0.5, 1) = 2");

    MeshView mv{};
    mv.diag = 8.0;
    DerivedTols tol;
    tol.epsPlane = 0.25;
    Provisional two;
    two.claim = ProvClaim::Unclaimed;
    two.tris = {0, 1};
    two.area = 1.0;
    Provisional three = two;
    three.tris = {0, 1, 2};
    const bool c2 = a3WouldCommit(mv, tol, two, false);
    const bool c3 = a3WouldCommit(mv, tol, three, false);
    std::fprintf(stderr, "a3WouldCommit(n=2,area=1)=%d a3WouldCommit(n=3,area=1)=%d\n",
                 c2 ? 1 : 0, c3 ? 1 : 0);
    expect(!c2, "a3WouldCommit 2-triangle tiny = false");
    expect(c3, "a3WouldCommit 3-triangle = true");
}

static bool nearRadius(const stl2step::refit::Region& r, double want, double tau) {
    return std::abs(r.radius - want) <= tau && r.maxVertexDev <= tau;
}

static void case1() {
    using namespace stl2step::refit;
    auto built = makeFillet(3.0, 0.0, 8.0, 0.0, 0.5 * kPi, 24, 0.0);
    const double q = 1e-4;
    Run r = runClaim(built.synth, std::move(built.work), q, 1);
    const bool one = r.work.accepted.size() == 1;
    bool ok = one;
    if (one) {
        const Region& a = r.work.accepted[0];
        const double span = a.uMax - a.uMin;
        ok = a.origin == Origin::FilletStrip && nearRadius(a, 3.0, r.tau) &&
             std::abs(span - 0.5 * kPi) <= 1e-3 && !a.closed360 &&
             a.tris.size() == built.bandTris.size();
        std::fprintf(stderr,
                     "case1 r=%.9g dev=%.3g span=%.6f tris=%d origin=%d\n", a.radius,
                     a.maxVertexDev, span, (int)a.tris.size(), (int)a.origin);
    }
    const std::string acc = lineWithWhy(r.diag, "accept");
    ok = ok && acc.find("fit=anchored") != std::string::npos &&
         acc.find("params=1") != std::string::npos;
    if (!ok) dumpFail("case1", r);
    expect(ok, "case1 concave R3 24-strip fillet anchored span 90");
}

static void case2() {
    using namespace stl2step::refit;
    Synth s;
    SegmentWork w;
    const double L = 8.0;
    const int a0 = s.addV(gp_XYZ(0, 0, 1));
    const int a1 = s.addV(gp_XYZ(0, L, 1));
    const int b0 = s.addV(gp_XYZ(1, 0, 0));
    const int b1 = s.addV(gp_XYZ(1, L, 0));
    std::vector<int> chamfer;
    s.addQuad(a0, a1, b0, b1, chamfer);
    addStrip(w, chamfer);
    const int wa = s.tri(a0, a1, s.addV(gp_XYZ(0, 0.5 * L, 31)));
    const int wb = s.tri(b0, b1, s.addV(gp_XYZ(31, 0.5 * L, 0)));
    addAnchor(w, {wa}, gp_XYZ(-1, 0, 0), gp_XYZ(0, 0, 0), 0.0);
    addAnchor(w, {wb}, gp_XYZ(0, 0, -1), gp_XYZ(0, 0, 0), 0.0);
    Run r = runClaim(s, std::move(w), 1e-4, 1);
    const bool ok = r.work.accepted.empty() && hasWhy(r.diag, "line");
    if (!ok) dumpFail("case2", r);
    expect(ok, "case2 chamfer facet why=line");
}

static void case3() {
    using namespace stl2step::refit;
    auto built = makeFillet(3.0, 0.0, 8.0, 0.0, 0.5 * kPi, 2, 0.0);
    Run r = runClaim(built.synth, std::move(built.work), 1e-4, 1);
    const std::string acc = lineWithWhy(r.diag, "accept");
    bool ok = r.work.accepted.size() == 1 && !acc.empty() &&
              acc.find("fit=anchored") != std::string::npos &&
              acc.find("params=1") != std::string::npos;
    if (ok) {
        const Region& a = r.work.accepted[0];
        ok = a.origin == Origin::FilletStrip && nearRadius(a, 3.0, r.tau);
        std::fprintf(stderr, "case3 r=%.9g stations-line=%s\n", a.radius, acc.c_str());
    }
    if (!ok) dumpFail("case3", r);
    expect(ok, "case3 3-station coarse fillet anchored params=1");
}

static void case4() {
    using namespace stl2step::refit;
    Synth s;
    SegmentWork w;
    const double L = 8.0;
    const int nStrips = 12;
    const int nAng = nStrips + 1;
    std::vector<int> lo((std::size_t)nAng), hi((std::size_t)nAng);
    for (int i = 0; i < nAng; ++i) {
        const double phi = kPi - kPi * (double)i / (double)nStrips;
        lo[(std::size_t)i] = s.addV(gp_XYZ(0.0, 1.0 + std::cos(phi), std::sin(phi)));
        hi[(std::size_t)i] = s.addV(gp_XYZ(L, 1.0 + std::cos(phi), std::sin(phi)));
    }
    for (int i = 0; i < nStrips; ++i) {
        std::vector<int> tris;
        s.addQuad(lo[(std::size_t)i], hi[(std::size_t)i], lo[(std::size_t)i + 1],
                  hi[(std::size_t)i + 1], tris);
        addStrip(w, tris);
    }
    const int wa = s.tri(lo.front(), hi.front(), s.addV(gp_XYZ(0.5 * L, 0.0, -30.0)));
    const int wb = s.tri(lo.back(), hi.back(), s.addV(gp_XYZ(0.5 * L, 2.0, -30.0)));
    addAnchor(w, {wa}, gp_XYZ(0, -1, 0), gp_XYZ(0, 0, 0), 0.0);
    addAnchor(w, {wb}, gp_XYZ(0, 1, 0), gp_XYZ(0, 2, 0), 0.0);
    Run r = runClaim(s, std::move(w), 1e-4, 1);
    bool ok = r.work.accepted.size() == 1;
    if (ok) {
        const Region& a = r.work.accepted[0];
        ok = nearRadius(a, 1.0, r.tau) && a.origin == Origin::FilletStrip;
        std::fprintf(stderr, "case4 r=%.9g dev=%.3g origin=%d span=%.6f\n", a.radius,
                     a.maxVertexDev, (int)a.origin, a.uMax - a.uMin);
    }
    if (!ok) dumpFail("case4", r);
    expect(ok, "case4 antiparallel hem r = gap/2");
}

static void case5() {
    using namespace stl2step::refit;
    Synth s;
    SegmentWork w;
    const double L = 20.0;
    const int nY = 4;
    const int nAng = 3;
    std::vector<std::vector<int>> id((std::size_t)nY + 1, std::vector<int>((std::size_t)nAng));
    for (int j = 0; j <= nY; ++j) {
        const double y = L * (double)j / (double)nY;
        for (int i = 0; i < nAng; ++i) {
            const double th = 0.5 * kPi * (double)i / (double)(nAng - 1);
            // Constant along the ruling, not one circle: the middle generator
            // sits at radius 8 while the walls are the R3 contact lines.
            const double rad = (i == 0 || i == nAng - 1) ? 3.0 : 8.0;
            id[(std::size_t)j][(std::size_t)i] =
                s.addV(gp_XYZ(rad * (1.0 - std::cos(th)), y, rad * (1.0 - std::sin(th))));
        }
    }
    std::vector<int> all;
    for (int j = 0; j < nY; ++j) {
        for (int i = 0; i < nAng - 1; ++i) {
            s.addQuad(id[(std::size_t)j][(std::size_t)i], id[(std::size_t)j + 1][(std::size_t)i],
                      id[(std::size_t)j][(std::size_t)i + 1],
                      id[(std::size_t)j + 1][(std::size_t)i + 1], all);
        }
    }
    addStrip(w, all);
    // Share a real ruling edge of the blend. A chord from the first vertex to
    // the last is not a mesh edge, so it would not make the blend a seed.
    const int wa = s.tri(id[0][0], id[1][0], s.addV(gp_XYZ(0, 0.5 * L / nY, 40)));
    const int wb = s.tri(id[0][nAng - 1], id[1][nAng - 1], s.addV(gp_XYZ(40, 0.5 * L / nY, 0)));
    addAnchor(w, {wa}, gp_XYZ(-1, 0, 0), gp_XYZ(0, 0, 0), 0.0);
    addAnchor(w, {wb}, gp_XYZ(0, 0, -1), gp_XYZ(0, 0, 0), 0.0);
    Run r = runClaim(s, std::move(w), 1e-4, 1);
    const bool ok = r.work.accepted.empty() && hasWhy(r.diag, "vertex-off");
    if (!ok) dumpFail("case5", r);
    expect(ok, "case5 variable-radius blend refused");
}

static void case6() {
    using namespace stl2step::refit;
    const double ang = 4.0 * kPi / 180.0;
    auto built = makeFillet(3.0, 0.0, 8.0, 0.0, ang, 4, 0.0);
    const double q = 0.01;
    Run r = runClaim(built.synth, std::move(built.work), q, 1);
    const bool ok = r.work.accepted.empty() && hasWhy(r.diag, "line");
    if (!ok) dumpFail("case6", r);
    expect(ok, "case6 sagitta within tau why=line");
}

static void case7() {
    using namespace stl2step::refit;
    auto alongY = makeFillet(3.0, 4.0, 12.0, 0.0, 0.5 * kPi, 8, 0.0);
    Synth s = std::move(alongY.synth);
    SegmentWork w = std::move(alongY.work);
    std::vector<int> yTris = alongY.bandTris;
    const int sphereA = yTris.empty() ? -1 : yTris.front();
    // Second fillet along +X between y=0 and z=0. Its own anchors.
    const double R = 3.0;
    const int nStrips = 8;
    std::vector<int> lo((std::size_t)nStrips + 1), hi((std::size_t)nStrips + 1);
    for (int i = 0; i <= nStrips; ++i) {
        const double th = 0.5 * kPi * (double)i / (double)nStrips;
        lo[(std::size_t)i] = s.addV(gp_XYZ(4.0, R * (1.0 - std::cos(th)), R * (1.0 - std::sin(th))));
        hi[(std::size_t)i] = s.addV(gp_XYZ(12.0, R * (1.0 - std::cos(th)), R * (1.0 - std::sin(th))));
    }
    std::vector<int> xTris;
    for (int i = 0; i < nStrips; ++i) {
        std::vector<int> tris;
        s.addQuad(lo[(std::size_t)i], hi[(std::size_t)i], lo[(std::size_t)i + 1],
                  hi[(std::size_t)i + 1], tris);
        for (int t : tris) xTris.push_back(t);
        addStrip(w, tris);
    }
    const int wallY = s.tri(lo.front(), hi.front(), s.addV(gp_XYZ(8.0, 0.0, R + 30.0)));
    const int wallZ = s.tri(lo.back(), hi.back(), s.addV(gp_XYZ(8.0, R + 30.0, 0.0)));
    addAnchor(w, {wallY}, gp_XYZ(0, -1, 0), gp_XYZ(0, 0, 0), 0.0);
    addAnchor(w, {wallZ}, gp_XYZ(0, 0, -1), gp_XYZ(0, 0, 0), 0.0);
    // Sphere-corner triangles hung on the Y-fillet's low end (edge perpendicular
    // to the ruling). The free vertex is off the cylinder.
    const int end0 = -1;
    (void)end0;
    (void)sphereA;
    // Recover the Y-fillet's y=4 arc vertices: they are the first 9 lo-points,
    // created before anything else. Indices 0,2,4,... are y0 of makeFillet.
    std::vector<int> sphereTris;
    const int off = s.addV(gp_XYZ(1.0, 1.0, 1.0));
    for (int i = 0; i < 8; ++i) {
        const int a = i * 2;
        const int b = (i + 1) * 2;
        sphereTris.push_back(s.tri(a, b, off));
    }
    addStrip(w, sphereTris);
    Run r = runClaim(s, std::move(w), 1e-4, 1);
    int fillets = 0;
    bool sphereIn = false;
    for (const Region& a : r.work.accepted) {
        if (a.origin == Origin::FilletStrip) ++fillets;
        for (int t : a.tris) {
            if (std::find(sphereTris.begin(), sphereTris.end(), t) != sphereTris.end())
                sphereIn = true;
        }
    }
    const bool ok = fillets == 2 && !sphereIn;
    std::fprintf(stderr, "case7 fillets=%d sphereIn=%d accepted=%d\n", fillets, sphereIn ? 1 : 0,
                 (int)r.work.accepted.size());
    if (!ok) dumpFail("case7", r);
    expect(ok, "case7 edge fillets claimed, sphere corner faceted");
}

static void case8() {
    using namespace stl2step::refit;
    Synth s;
    SegmentWork w;
    const double R = 4.0;
    const int nStrips = 16;
    const int nY = 2;
    std::vector<std::vector<int>> id((std::size_t)nY + 1,
                                     std::vector<int>((std::size_t)nStrips + 1));
    for (int j = 0; j <= nY; ++j) {
        const double y = 10.0 * (double)j;
        for (int i = 0; i <= nStrips; ++i) {
            const double th = 2.0 * kPi * (double)i / (double)nStrips;
            id[(std::size_t)j][(std::size_t)i] =
                s.addV(gp_XYZ(R * std::cos(th), y, R * std::sin(th)));
        }
    }
    std::vector<int> consumed, open;
    for (int i = 0; i < nStrips; ++i) {
        std::vector<int>* dest = (i < nStrips / 2) ? &consumed : &open;
        for (int j = 0; j < nY; ++j) {
            std::vector<int> tris;
            s.addQuad(id[(std::size_t)j][(std::size_t)i], id[(std::size_t)j + 1][(std::size_t)i],
                      id[(std::size_t)j][(std::size_t)i + 1],
                      id[(std::size_t)j + 1][(std::size_t)i + 1], tris);
            for (int t : tris) dest->push_back(t);
            Provisional pr;
            pr.tris = tris;
            pr.area = 1.0;
            pr.maxVertexDev = 0.0;
            pr.plane = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
            if (i < nStrips / 2) {
                pr.claim = ProvClaim::ConsumedCylinder;
                w.provisionals.push_back(pr);
            } else {
                pr.claim = ProvClaim::Unclaimed;
                w.provisionals.push_back(pr);
            }
        }
    }
    Region acc;
    acc.id = 0;
    acc.type = SurfType::Cylinder;
    acc.origin = Origin::CylGrow;
    acc.radius = R;
    acc.tris = consumed;
    acc.ax = gp_Ax3(gp_Pnt(0, 10, 0), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
    w.accepted.push_back(acc);
    Run r = runClaim(s, std::move(w), 1e-4, 1);
    // The pre-seeded cylinder stays. Nothing new is claimed.
    const bool ok = r.work.accepted.size() == 1 && hasWhy(r.diag, "same-surface");
    if (!ok) dumpFail("case8", r);
    expect(ok, "case8 coaxial accepted cylinder why=same-surface");
}

static void case9() {
    using namespace stl2step::refit;
    const double q = 1e-4;
    const double tau = 2.0 * q;
    auto built = makeFillet(3.0, 0.0, 8.0, 0.0, 0.5 * kPi, 2, 10.0 * tau);
    Run r = runClaim(built.synth, std::move(built.work), q, 1);
    const bool ok = r.work.accepted.empty() && hasWhy(r.diag, "underdetermined");
    if (!ok) dumpFail("case9", r);
    expect(ok, "case9 free stations < params+1 why=underdetermined");
}

static void case10() {
    using namespace stl2step::refit;
    Synth s;
    SegmentWork w;
    const double R = 5.0;
    const double sweep = 356.0 * kPi / 180.0;
    const int nStrips = 40;
    const int nY = 2;
    std::vector<std::vector<int>> id((std::size_t)nY + 1,
                                     std::vector<int>((std::size_t)nStrips + 1));
    for (int j = 0; j <= nY; ++j) {
        const double y = 6.0 * (double)j;
        for (int i = 0; i <= nStrips; ++i) {
            const double th = sweep * (double)i / (double)nStrips;
            id[(std::size_t)j][(std::size_t)i] =
                s.addV(gp_XYZ(R * std::cos(th), y, R * std::sin(th)));
        }
    }
    std::vector<int> band;
    for (int i = 0; i < nStrips; ++i) {
        std::vector<int> tris;
        for (int j = 0; j < nY; ++j) {
            s.addQuad(id[(std::size_t)j][(std::size_t)i], id[(std::size_t)j + 1][(std::size_t)i],
                      id[(std::size_t)j][(std::size_t)i + 1],
                      id[(std::size_t)j + 1][(std::size_t)i + 1], tris);
        }
        for (int t : tris) band.push_back(t);
        addStrip(w, tris);
    }
    const double q = 1e-4;
    const double dev = 10.0 * 2.0 * q;
    const gp_XYZ n0(1, 0, 0);
    const gp_XYZ p0(R, 0, 0);
    const double th1 = sweep;
    const gp_XYZ n1(std::cos(th1), 0, std::sin(th1));
    const gp_XYZ p1(R * n1.X(), 0, R * n1.Z());
    const int w0 = s.tri(id[0][0], id[(std::size_t)nY][0],
                         s.addV(gp_XYZ(p0.X() + 30.0 * n0.X(), 6.0, p0.Z() + 30.0 * n0.Z())));
    const int w1 = s.tri(id[0][(std::size_t)nStrips], id[(std::size_t)nY][(std::size_t)nStrips],
                         s.addV(p1 + n1 * 30.0 + gp_XYZ(0, 6, 0)));
    addAnchor(w, {w0}, n0, p0, dev);
    addAnchor(w, {w1}, n1, p1, dev);
    Run r = runClaim(s, std::move(w), q, 1);
    bool ok = r.work.accepted.size() == 1;
    if (ok) {
        const Region& a = r.work.accepted[0];
        ok = a.origin == Origin::CylGrow && a.closed360 && nearRadius(a, R, r.tau) &&
             a.origin != Origin::FilletStrip;
        std::fprintf(stderr, "case10 r=%.9g dev=%.3g closed=%d origin=%d tris=%d\n", a.radius,
                     a.maxVertexDev, a.closed360 ? 1 : 0, (int)a.origin, (int)a.tris.size());
    }
    if (!ok) dumpFail("case10", r);
    expect(ok, "case10 356 deg cylinder is CylGrow");
}

static void case11() {
    using namespace stl2step::refit;
    auto a = makeFillet(3.0, 0.0, 8.0, 0.0, 0.5 * kPi, 24, 0.0);
    auto b = makeFillet(3.0, 0.0, 8.0, 0.0, 0.5 * kPi, 24, 0.0);
    Run ra = runClaim(a.synth, std::move(a.work), 1e-4, 1);
    Run rb = runClaim(b.synth, std::move(b.work), 1e-4, 4);
    bool ok = ra.work.accepted.size() == 1 && rb.work.accepted.size() == 1;
    if (ok) {
        const Region& A = ra.work.accepted[0];
        const Region& B = rb.work.accepted[0];
        const gp_XYZ dA = A.ax.Direction().XYZ();
        const gp_XYZ dB = B.ax.Direction().XYZ();
        ok = A.tris == B.tris && std::abs(A.radius - B.radius) <= 0.0 &&
             A.origin == B.origin &&
             (A.ax.Location().XYZ() - B.ax.Location().XYZ()).Modulus() <= 1e-12 &&
             (dA - dB).Modulus() <= 1e-12;
        std::fprintf(stderr, "case11 r=%.9g tris=%d\n", A.radius, (int)A.tris.size());
    }
    if (!ok) {
        dumpFail("case11-a", ra);
        dumpFail("case11-b", rb);
    }
    expect(ok, "case11 threads 1 and 4 identical");
}

int main() {
    setenv("STL2STEP_DIAG_REL", "1", 1);
    pins();
    case1();
    case2();
    case3();
    case4();
    case5();
    case6();
    case7();
    case8();
    case9();
    case10();
    case11();
    if (gFail) {
        std::fprintf(stderr, "%d check(s) failed\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "relational_extrude_unit ok\n");
    return 0;
}
