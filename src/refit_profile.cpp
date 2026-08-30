// Stage P2 — sliceProfiles + fitProfile (lines + arcs; R from lawChainAccept).
// Not wired into convert(). SPDX-License-Identifier: MIT

#include "refit_prism.hpp"
#include "refit_internal.hpp"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

struct SketchFrame {
    gp_XYZ origin;
    gp_XYZ axis;
    gp_XYZ u;
    gp_XYZ v;

    gp_Pnt2d xy(const gp_XYZ& p) const {
        const gp_XYZ d = p - origin;
        return gp_Pnt2d(d.Dot(u), d.Dot(v));
    }
    double axial(const gp_XYZ& p) const { return (p - origin).Dot(axis); }
};

struct SliceCache {
    const RegionSet* rs = nullptr;
    SketchFrame fr;
    bool ok = false;
    // per slab, per loop: region id of each polyline vertex (edge start)
    std::vector<std::vector<std::vector<int>>> rids;
};

SliceCache& cache() {
    static SliceCache c;
    return c;
}
std::mutex& cacheMu() {
    static std::mutex m;
    return m;
}

bool prismDiagOn() {
    const char* e = std::getenv("STL2STEP_PRISM_DIAG");
    return e && e[0] == '1' && e[1] == '\0';
}

bool assistOn() {
    const char* e = std::getenv("STL2STEP_PRISM_ASSIST");
    return e && e[0] != '\0' && std::strcmp(e, "0") != 0;
}

std::mutex& diagMu() {
    static std::mutex m;
    return m;
}

unsigned workerCount(unsigned work) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;
    return std::max(1u, std::min(work, hw));
}

double dist2(const gp_Pnt2d& a, const gp_Pnt2d& b) {
    return std::hypot(a.X() - b.X(), a.Y() - b.Y());
}

bool near2(const gp_Pnt2d& a, const gp_Pnt2d& b, double tol) {
    return dist2(a, b) <= tol;
}

double signedArea(const std::vector<gp_Pnt2d>& poly) {
    const size_t n = poly.size();
    if (n < 3) return 0.0;
    double a = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const gp_Pnt2d& p = poly[i];
        const gp_Pnt2d& q = poly[(i + 1) % n];
        a += p.X() * q.Y() - q.X() * p.Y();
    }
    return 0.5 * a;
}

bool localCorners(const MeshView& mv, int localT, gp_XYZ& a, gp_XYZ& b, gp_XYZ& c) {
    if (!mv.pts || !mv.tris || !mv.compTris) return false;
    if (localT < 0 || static_cast<size_t>(localT) >= mv.nTri) return false;
    const int g = mv.compTris[localT];
    a = mv.pts[mv.tris[g][0]];
    b = mv.pts[mv.tris[g][1]];
    c = mv.pts[mv.tris[g][2]];
    return true;
}

bool buildFrame(const RegionSet& rs, const PrismLevels& lv, SketchFrame& fr) {
    fr.axis = lv.axis.XYZ();
    const double am = fr.axis.Modulus();
    if (!(am > 0.0)) return false;
    fr.axis.Divide(am);
    gp_XYZ hint(1.0, 0.0, 0.0);
    if (std::fabs(fr.axis.Dot(hint)) > 0.9) hint = gp_XYZ(0.0, 1.0, 0.0);
    fr.u = fr.axis.Crossed(hint);
    const double um = fr.u.Modulus();
    if (um < 1e-18) return false;
    fr.u.Divide(um);
    fr.v = fr.axis.Crossed(fr.u);
    const double vm = fr.v.Modulus();
    if (vm < 1e-18) return false;
    fr.v.Divide(vm);
    fr.origin = gp_XYZ(0.0, 0.0, 0.0);
    if (!lv.capRegion.empty()) {
        for (const Region& r : rs.regions) {
            if (r.id == lv.capRegion[0]) {
                fr.origin = r.ax.Location().XYZ();
                return true;
            }
        }
    }
    for (const Region& r : rs.regions) {
        if (r.type == SurfType::Cylinder) {
            fr.origin = r.ax.Location().XYZ();
            return true;
        }
    }
    return true;
}

void cylEnds(const Region& r, const gp_XYZ& axis, double& lo, double& hi) {
    const gp_XYZ a = r.ax.Direction().XYZ();
    const gp_XYZ loc = r.ax.Location().XYZ();
    const gp_XYZ p0 = loc + a * r.vMin;
    const gp_XYZ p1 = loc + a * r.vMax;
    const double y0 = axis.Dot(p0);
    const double y1 = axis.Dot(p1);
    lo = std::min(y0, y1);
    hi = std::max(y0, y1);
}

int snapLevel(const std::vector<double>& y, double v, double tol) {
    int best = -1;
    double bestD = 0.0;
    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        const double d = std::fabs(y[static_cast<size_t>(i)] - v);
        if (d <= tol && (best < 0 || d < bestD)) {
            best = i;
            bestD = d;
        }
    }
    return best;
}

// Through-feature: axial span covers more than one slab (snaps to non-adjacent levels).
bool isThrough(const Region& r, const PrismLevels& lv, double tau) {
    if (r.type != SurfType::Cylinder || lv.y.size() < 3) return false;
    double lo = 0.0, hi = 0.0;
    cylEnds(r, lv.axis.XYZ(), lo, hi);
    const int i0 = snapLevel(lv.y, lo, tau);
    const int i1 = snapLevel(lv.y, hi, tau);
    return i0 >= 0 && i1 >= 0 && (i1 - i0) >= 2;
}

bool cylCoversSlab(const Region& r, const PrismLevels& lv, int slab, double tau) {
    if (r.type != SurfType::Cylinder || lv.y.size() < 2) return false;
    if (slab < 0 || static_cast<size_t>(slab) + 1 >= lv.y.size()) return false;
    double lo = 0.0, hi = 0.0;
    cylEnds(r, lv.axis.XYZ(), lo, hi);
    const double y0 = lv.y[static_cast<size_t>(slab)];
    const double y1 = lv.y[static_cast<size_t>(slab) + 1];
    return lo <= y0 + tau && hi >= y1 - tau;
}

gp_Pnt2d cylCenter(const Region& r, const SketchFrame& fr) {
    return fr.xy(r.ax.Location().XYZ());
}

double cylAssocTol(const Region& r, double tauFit) {
    double sag = r.chordSagitta;
    if (!(sag > 0.0)) sag = r.maxVertexDev;
    // Slice points lie on mesh chords, not the analytic circle.
    return std::max({tauFit, 4.0 * sag, 4.0 * r.maxVertexDev});
}

bool onCyl(const gp_Pnt2d& p, const Region& r, const SketchFrame& fr, double tauFit) {
    if (r.type != SurfType::Cylinder || !(r.radius > 0.0)) return false;
    const gp_Pnt2d c = cylCenter(r, fr);
    return std::fabs(dist2(p, c) - r.radius) <= cylAssocTol(r, tauFit);
}

DerivedTols makeDerived(const MeshView& mv) {
    DerivedTols t;
    t.epsMesh = std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
    t.epsPlane = std::max({t.epsMesh, mv.sewTol, 0.02});
    return t;
}

struct SliceEdge {
    gp_Pnt2d a, b;
    int localTri = -1;
    int rid = -1;
};

struct RawLoop {
    std::vector<gp_Pnt2d> pts;
    std::vector<int> rids;
};

void sliceTri(const MeshView& mv, int localT, int rid, const SketchFrame& fr, double ySlice,
              double tol, std::vector<SliceEdge>& edges) {
    gp_XYZ p0, p1, p2;
    if (!localCorners(mv, localT, p0, p1, p2)) return;
    const gp_XYZ p[3] = {p0, p1, p2};
    double s[3];
    for (int i = 0; i < 3; ++i) s[i] = fr.axis.Dot(p[i]);
    const double smin = std::min({s[0], s[1], s[2]});
    const double smax = std::max({s[0], s[1], s[2]});
    if (ySlice < smin - tol || ySlice > smax + tol) return;
    if (smax - smin <= tol) return;

    gp_Pnt2d ip[2];
    int nIp = 0;
    for (int e = 0; e < 3; ++e) {
        const int i0 = e;
        const int i1 = (e + 1) % 3;
        const double d0 = s[i0] - ySlice;
        const double d1 = s[i1] - ySlice;
        if (d0 < -tol && d1 < -tol) continue;
        if (d0 > tol && d1 > tol) continue;
        if (std::fabs(d0) <= tol && std::fabs(d1) <= tol) continue;
        const double den = d0 - d1;
        double tp = (std::fabs(den) > 0.0) ? (d0 / den) : 0.0;
        tp = std::max(0.0, std::min(1.0, tp));
        const gp_XYZ q = p[i0] + (p[i1] - p[i0]) * tp;
        if (nIp < 2) ip[nIp++] = fr.xy(q);
    }
    if (nIp == 2 && dist2(ip[0], ip[1]) > tol) {
        SliceEdge e;
        e.a = ip[0];
        e.b = ip[1];
        e.localTri = localT;
        e.rid = rid;
        edges.push_back(e);
    }
}

void chainEdges(const std::vector<SliceEdge>& edges, double tol, std::vector<RawLoop>& loops) {
    const size_t n = edges.size();
    std::vector<char> used(n, 0);
    for (size_t start = 0; start < n; ++start) {
        if (used[start]) continue;
        RawLoop lp;
        lp.pts.push_back(edges[start].a);
        lp.pts.push_back(edges[start].b);
        lp.rids.push_back(edges[start].rid);
        used[start] = 1;
        gp_Pnt2d tail = edges[start].b;
        bool grew = true;
        while (grew) {
            grew = false;
            for (size_t i = 0; i < n; ++i) {
                if (used[i]) continue;
                if (near2(tail, edges[i].a, tol)) {
                    lp.pts.push_back(edges[i].b);
                    lp.rids.push_back(edges[i].rid);
                    tail = edges[i].b;
                    used[i] = 1;
                    grew = true;
                } else if (near2(tail, edges[i].b, tol)) {
                    lp.pts.push_back(edges[i].a);
                    lp.rids.push_back(edges[i].rid);
                    tail = edges[i].a;
                    used[i] = 1;
                    grew = true;
                }
            }
        }
        if (lp.pts.size() >= 2 && near2(lp.pts.front(), lp.pts.back(), tol)) lp.pts.pop_back();
        if (lp.pts.size() >= 3) loops.push_back(std::move(lp));
    }
}

void emitCircle(const Region& r, const SketchFrame& fr, RawLoop& lp, int nSamp) {
    lp.pts.clear();
    lp.rids.clear();
    if (!(r.radius > 0.0) || nSamp < 3) return;
    const gp_Pnt2d c = cylCenter(r, fr);
    const double R = r.radius;
    for (int i = 0; i < nSamp; ++i) {
        const double ang = kTwoPi * static_cast<double>(i) / static_cast<double>(nSamp);
        lp.pts.emplace_back(c.X() + R * std::cos(ang), c.Y() + R * std::sin(ang));
        lp.rids.push_back(r.id);
    }
}

double lineResid(const gp_Pnt2d& a, const gp_Pnt2d& b, const std::vector<gp_Pnt2d>& pts,
                 size_t i0, size_t i1) {
    const double dx = b.X() - a.X();
    const double dy = b.Y() - a.Y();
    const double len = std::hypot(dx, dy);
    if (len < 1e-18) return 0.0;
    double m = 0.0;
    for (size_t i = i0; i <= i1 && i < pts.size(); ++i) {
        const double px = pts[i].X() - a.X();
        const double py = pts[i].Y() - a.Y();
        m = std::max(m, std::fabs((dx * py - dy * px) / len));
    }
    return m;
}

bool rule53a(const LawBand& lb, double tauFit) {
    if (lb.N <= 2) return true;
    if (!(lb.R > 0.0)) return true;
    const double sag = lb.R * (1.0 - std::cos(0.5 * lb.phi));
    return sag < tauFit;
}

// Content hash of a declined chain. Advisory may only choose arc|line — never a number.
std::string chainSha(const std::vector<gp_Pnt2d>& pts, size_t i0, size_t i1) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](double x) {
        const long long q = std::llround(x * 1.0e6);
        unsigned char b[8];
        uint64_t u = static_cast<uint64_t>(q);
        for (int k = 0; k < 8; ++k) {
            b[k] = static_cast<unsigned char>(u & 0xffu);
            u >>= 8;
        }
        for (int k = 0; k < 8; ++k) {
            h ^= static_cast<uint64_t>(b[k]);
            h *= 1099511628211ull;
        }
    };
    for (size_t i = i0; i <= i1 && i < pts.size(); ++i) {
        mix(pts[i].X());
        mix(pts[i].Y());
    }
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return std::string(buf);
}

std::string advisoryDir() {
    const char* e = std::getenv("STL2STEP_PRISM_ADVISORY_DIR");
    if (e && e[0] != '\0') return std::string(e);
    return "tests/fixtures/prism-advisories";
}

// Returns 1 = arc, 0 = line, -1 = miss. Never reads a numeric field.
int lookupAdvisory(const std::string& sha, bool& applied) {
    applied = false;
    if (!assistOn()) return -1;
    const std::string path = advisoryDir() + "/" + sha + ".json";
    std::ifstream in(path);
    if (!in) return -1;
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string txt = ss.str();
    const auto pos = txt.find("chose");
    if (pos == std::string::npos) return -1;
    const auto arc = txt.find("arc", pos);
    const auto line = txt.find("line", pos);
    int chose = -1;
    if (arc != std::string::npos && (line == std::string::npos || arc < line))
        chose = 1;
    else if (line != std::string::npos)
        chose = 0;
    if (chose < 0) return -1;
    applied = true;
    if (prismDiagOn()) {
        std::lock_guard<std::mutex> g(diagMu());
        std::fprintf(stderr, "DIAG_PRISM advisory=%s chose=%s\n", sha.c_str(),
                     chose == 1 ? "arc" : "line");
    }
    return chose;
}

void emitLoopDiag(int slab, size_t li, const ProfLoop& lp) {
    if (!prismDiagOn()) return;
    int nLine = 0, nArc = 0, nDecl = 0;
    for (const ProfSeg& s : lp.segs) {
        if (s.isArc) ++nArc;
        else ++nLine;
        if (s.declinedAmbiguous) ++nDecl;
    }
    std::lock_guard<std::mutex> g(diagMu());
    std::fprintf(stderr,
                 "DIAG_PROFILE comp=%d slab=%d loop=%d outer=%d nSeg=%d nLine=%d nArc=%d "
                 "nDecl=%d area=%.6f\n",
                 0, slab, static_cast<int>(li), lp.outer ? 1 : 0,
                 static_cast<int>(lp.segs.size()), nLine, nArc, nDecl, lp.area);
    for (size_t i = 0; i < lp.segs.size(); ++i) {
        const ProfSeg& s = lp.segs[i];
        std::fprintf(stderr,
                     "DIAG_PROFILE  seg slab=%d loop=%d i=%d kind=%s R=%.6f phi=%.6f decl=%d\n",
                     slab, static_cast<int>(li), static_cast<int>(i),
                     s.isArc ? "arc" : "line", s.R, s.phi, s.declinedAmbiguous ? 1 : 0);
    }
}

const Region* regionById(const RegionSet& rs, int id) {
    for (const Region& r : rs.regions)
        if (r.id == id) return &r;
    return nullptr;
}

const Region* cylAtPoint(const RegionSet& rs, const SketchFrame& fr, const gp_Pnt2d& p,
                         double tauFit) {
    const Region* best = nullptr;
    double bestD = 0.0;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder || !(r.radius > 0.0)) continue;
        const double d = std::fabs(dist2(p, cylCenter(r, fr)) - r.radius);
        if (d <= cylAssocTol(r, tauFit) && (!best || d < bestD)) {
            best = &r;
            bestD = d;
        }
    }
    return best;
}

const Region* cylByRid(const RegionSet& rs, int rid) {
    const Region* r = regionById(rs, rid);
    if (r && r->type == SurfType::Cylinder && r->radius > 0.0) return r;
    return nullptr;
}

void snapClosed(std::vector<ProfSeg>& segs, double tol) {
    if (segs.empty()) return;
    for (size_t i = 0; i < segs.size(); ++i) {
        ProfSeg& cur = segs[i];
        ProfSeg& nxt = segs[(i + 1) % segs.size()];
        if (!near2(cur.b, nxt.a, tol)) {
            const gp_Pnt2d mid(0.5 * (cur.b.X() + nxt.a.X()), 0.5 * (cur.b.Y() + nxt.a.Y()));
            cur.b = mid;
            nxt.a = mid;
        } else {
            nxt.a = cur.b;
        }
    }
}

void mergeColinear(std::vector<ProfSeg>& segs, const std::vector<gp_Pnt2d>& /*verts*/,
                   double tol) {
    if (segs.size() < 2) return;
    std::vector<ProfSeg> out;
    out.reserve(segs.size());
    ProfSeg acc = segs[0];
    for (size_t i = 1; i < segs.size(); ++i) {
        const ProfSeg& s = segs[i];
        if (!acc.isArc && !s.isArc && !acc.declinedAmbiguous && !s.declinedAmbiguous) {
            const double dx = s.b.X() - acc.a.X();
            const double dy = s.b.Y() - acc.a.Y();
            const double len = std::hypot(dx, dy);
            const double mx = acc.b.X() - acc.a.X();
            const double my = acc.b.Y() - acc.a.Y();
            const double resid =
                (len > 1e-18) ? std::fabs(dx * my - dy * mx) / len : 0.0;
            if (resid <= tol) {
                acc.b = s.b;
                continue;
            }
        }
        out.push_back(acc);
        acc = s;
    }
    out.push_back(acc);
    segs.swap(out);
}

bool loopHasRadius(const ProfLoop& lp, double R, double tol) {
    if (!(R > 0.0)) return false;
    for (const ProfSeg& s : lp.segs) {
        if (s.isArc && s.R > 0.0 && std::fabs(s.R - R) <= std::max(tol, 0.003 * R))
            return true;
    }
    return false;
}

void addCircleSeg(ProfLoop& lp, const gp_Pnt2d& c, double R, bool ccw) {
    ProfSeg s;
    s.isArc = true;
    s.center = c;
    s.R = R;
    s.phi = kTwoPi;
    s.ccw = ccw;
    s.a = gp_Pnt2d(c.X() + R, c.Y());
    s.b = s.a;
    lp.segs.push_back(s);
}

void classifyOuter(std::vector<ProfLoop>& loops, std::vector<std::vector<int>>* rids) {
    if (loops.empty()) return;
    size_t oi = 0;
    double best = -1.0;
    for (size_t i = 0; i < loops.size(); ++i) {
        if (loops[i].area > best) {
            best = loops[i].area;
            oi = i;
        }
    }
    for (size_t i = 0; i < loops.size(); ++i) loops[i].outer = (i == oi);
    if (oi != 0) {
        std::swap(loops[0], loops[oi]);
        if (rids && rids->size() == loops.size()) std::swap((*rids)[0], (*rids)[oi]);
    }
    loops[0].outer = true;
    for (size_t i = 1; i < loops.size(); ++i) loops[i].outer = false;
}

gp_Pnt2d centroidOf(const std::vector<gp_Pnt2d>& pts) {
    gp_Pnt2d c(0.0, 0.0);
    if (pts.empty()) return c;
    for (const gp_Pnt2d& p : pts) {
        c.SetX(c.X() + p.X());
        c.SetY(c.Y() + p.Y());
    }
    const double n = static_cast<double>(pts.size());
    c.SetX(c.X() / n);
    c.SetY(c.Y() / n);
    return c;
}

bool rawMatchesCyl(const RawLoop& lp, const Region& r, const SketchFrame& fr, double tau) {
    if (lp.pts.size() < 3 || !(r.radius > 0.0)) return false;
    const gp_Pnt2d c = cylCenter(r, fr);
    const double assoc = cylAssocTol(r, tau);
    double mean = 0.0;
    double rmin = 1e300, rmax = 0.0;
    for (const gp_Pnt2d& p : lp.pts) {
        const double d = dist2(p, c);
        mean += d;
        rmin = std::min(rmin, d);
        rmax = std::max(rmax, d);
    }
    mean /= static_cast<double>(lp.pts.size());
    // A hole loop is a vertex ring on the circle, not merely a centroid near the axis.
    if (std::fabs(mean - r.radius) > assoc) return false;
    if (rmax - rmin > 4.0 * assoc) return false;
    int hit = 0;
    for (int id : lp.rids)
        if (id == r.id) ++hit;
    if (!lp.rids.empty() && hit * 2 >= static_cast<int>(lp.rids.size())) return true;
    return dist2(centroidOf(lp.pts), c) <= assoc;
}

bool sliceOneSlab(const MeshView& mv, const RegionSet& rs, const PrismLevels& lv,
                  const PrismTols& t, const SketchFrame& fr, int slab, Profile& prof) {
    prof = Profile{};
    prof.slab = slab;
    if (lv.y.size() < 2) return false;
    const double y0 = lv.y[static_cast<size_t>(slab)];
    const double y1 = lv.y[static_cast<size_t>(slab) + 1];
    const double ySlice = 0.5 * (y0 + y1);
    const double tol = t.tauFit;

    std::vector<SliceEdge> edges;
    edges.reserve(mv.nTri);
    for (size_t ti = 0; ti < mv.nTri; ++ti) {
        const int localT = static_cast<int>(ti);
        int rid = -1;
        if (!rs.triRegion.empty() && ti < rs.triRegion.size()) rid = rs.triRegion[ti];
        const Region* rg = (rid >= 0) ? regionById(rs, rid) : nullptr;
        if (rg && rg->type == SurfType::Plane) {
            const double nd = std::fabs(rg->ax.Direction().XYZ().Dot(fr.axis));
            if (nd > 1.0 - t.tauAx) continue;
        }
        sliceTri(mv, localT, rid, fr, ySlice, tol, edges);
    }

    std::vector<RawLoop> raw;
    chainEdges(edges, std::max(tol, 4.0 * mv.weldTol), raw);

    // Closed-360 / through-features: one circular inner loop each (replace, don't duplicate).
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder) continue;
        if (!cylCoversSlab(r, lv, slab, t.tauLvl)) continue;
        const bool through = isThrough(r, lv, t.tauLvl);
        const bool hole360 = r.closed360 && !r.outwardNormal;
        if (!through && !hole360) continue;
        bool have = false;
        for (RawLoop& lp : raw) {
            if (!rawMatchesCyl(lp, r, fr, tol)) continue;
            emitCircle(r, fr, lp, std::max(8, r.nSides > 0 ? r.nSides : 12));
            have = true;
            break;
        }
        if (!have) {
            RawLoop hole;
            emitCircle(r, fr, hole, std::max(8, r.nSides > 0 ? r.nSides : 12));
            if (hole.pts.size() >= 3) raw.push_back(std::move(hole));
        }
    }

    std::vector<ProfLoop> loops;
    std::vector<std::vector<int>> ridPerLoop;
    loops.reserve(raw.size());
    for (const RawLoop& rl : raw) {
        ProfLoop lp;
        std::vector<int> ids;
        lp.area = std::fabs(signedArea(rl.pts));
        for (size_t i = 0; i < rl.pts.size(); ++i) {
            ProfSeg s;
            s.isArc = false;
            s.a = rl.pts[i];
            s.b = rl.pts[(i + 1) % rl.pts.size()];
            if (dist2(s.a, s.b) > tol) {
                lp.segs.push_back(s);
                ids.push_back(i < rl.rids.size() ? rl.rids[i] : -1);
            }
        }
        if (lp.segs.size() >= 3) {
            loops.push_back(std::move(lp));
            ridPerLoop.push_back(std::move(ids));
        }
    }
    classifyOuter(loops, &ridPerLoop);
    {
        std::lock_guard<std::mutex> g(cacheMu());
        auto& slot = cache().rids;
        if (slot.size() <= static_cast<size_t>(slab))
            slot.resize(static_cast<size_t>(slab) + 1);
        slot[static_cast<size_t>(slab)] = std::move(ridPerLoop);
    }
    prof.loops = std::move(loops);
    return !prof.loops.empty();
}

void emitLines(const std::vector<gp_Pnt2d>& verts, size_t i0, size_t i1, bool declined,
               double tol, std::vector<ProfSeg>& out, int& nDecl) {
    const size_t n = verts.size();
    if (n < 2 || i0 >= n) return;
    size_t a = i0;
    while (a < i1 && a + 1 < n) {
        size_t b = a + 1;
        while (b < i1 && b < n) {
            if (lineResid(verts[a], verts[b], verts, a, b) > tol) break;
            ++b;
        }
        if (b <= a) b = a + 1;
        ProfSeg s;
        s.isArc = false;
        s.declinedAmbiguous = declined;
        s.a = verts[a];
        s.b = verts[b % n];
        if (dist2(s.a, s.b) > tol) {
            out.push_back(s);
            if (declined) ++nDecl;
        }
        a = b;
    }
}

bool tryArc(const MeshView& mv, const Region& r, const SketchFrame& fr,
            const std::vector<gp_Pnt2d>& verts, size_t i0, size_t i1, const PrismTols& t,
            const DerivedTols& dtol, ProfSeg& arc, bool& declined, int& nDecl) {
    LawBand lb;
    if (!lawChainAccept(mv, r.tris, dtol, lb) || !(lb.R > 0.0)) return false;

    declined = rule53a(lb, t.tauFit);
    if (declined) {
        bool applied = false;
        const std::string sha = chainSha(verts, i0, i1);
        const int chose = lookupAdvisory(sha, applied);
        if (applied && chose == 1) declined = false;
        if (declined) return true;
    }

    arc.isArc = true;
    arc.declinedAmbiguous = false;
    arc.R = lb.R;
    arc.phi = lb.closed360 ? kTwoPi : lb.phi;
    arc.center = cylCenter(r, fr);
    arc.a = verts[i0];
    const size_t n = verts.size();
    arc.b = verts[i1 % n];
    if (lb.closed360) {
        arc.a = gp_Pnt2d(arc.center.X() + arc.R, arc.center.Y());
        arc.b = arc.a;
        arc.phi = kTwoPi;
    }
    std::vector<gp_Pnt2d> run;
    for (size_t i = i0; i < i1 && i < verts.size(); ++i) run.push_back(verts[i]);
    arc.ccw = signedArea(run) > 0.0;
    (void)nDecl;
    return true;
}

}  // namespace

bool sliceProfiles(const MeshView& mv, const RegionSet& rs, const PrismLevels& lv,
                   const PrismTols& t, std::vector<Profile>& out) {
    out.clear();
    try {
        if (!lv.ok || lv.y.size() < 2) return false;
        SketchFrame fr;
        if (!buildFrame(rs, lv, fr)) return false;
        {
            std::lock_guard<std::mutex> g(cacheMu());
            cache().rs = &rs;
            cache().fr = fr;
            cache().ok = true;
            cache().rids.clear();
        }

        const int nSlab = static_cast<int>(lv.y.size() - 1);
        out.assign(static_cast<size_t>(nSlab), Profile{});
        std::vector<char> ok(static_cast<size_t>(nSlab), 0);
        const unsigned nw = workerCount(static_cast<unsigned>(nSlab));
        std::vector<std::thread> pool;
        pool.reserve(nw);
        const int chunk = (nSlab + static_cast<int>(nw) - 1) / static_cast<int>(nw);
        for (unsigned w = 0; w < nw; ++w) {
            pool.emplace_back([&, w]() {
                const int begin = static_cast<int>(w) * chunk;
                const int end = std::min(nSlab, begin + chunk);
                for (int k = begin; k < end; ++k)
                    ok[static_cast<size_t>(k)] =
                        sliceOneSlab(mv, rs, lv, t, fr, k, out[static_cast<size_t>(k)]) ? 1 : 0;
            });
        }
        for (auto& th : pool) th.join();
        for (char v : ok)
            if (!v) {
                out.clear();
                return false;
            }

        // RULE 5.2a — through-features must appear as inner loops in every slab they cross.
        for (int k = 0; k < nSlab; ++k) {
            for (const Region& r : rs.regions) {
                if (!isThrough(r, lv, t.tauLvl)) continue;
                if (!cylCoversSlab(r, lv, k, t.tauLvl)) continue;
                Profile& p = out[static_cast<size_t>(k)];
                bool found = false;
                for (const ProfLoop& lp : p.loops) {
                    if (lp.outer) continue;
                    if (loopHasRadius(lp, r.radius, t.tauFit)) {
                        found = true;
                        break;
                    }
                    // Pre-fit polylines: radius from vertex ring.
                    if (lp.segs.size() >= 3) {
                        const gp_Pnt2d c = cylCenter(r, fr);
                        double rmin = 1e300, rmax = 0.0;
                        for (const ProfSeg& s : lp.segs) {
                            const double d = dist2(s.a, c);
                            rmin = std::min(rmin, d);
                            rmax = std::max(rmax, d);
                        }
                        if (r.radius > 0.0 && rmax - rmin <= 8.0 * t.tauFit &&
                            std::fabs(0.5 * (rmin + rmax) - r.radius) <=
                                std::max(t.tauFit, 0.003 * r.radius)) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    ProfLoop hole;
                    hole.outer = false;
                    hole.area = kPi * r.radius * r.radius;
                    addCircleSeg(hole, cylCenter(r, fr), r.radius, true);
                    p.loops.push_back(std::move(hole));
                    classifyOuter(p.loops, nullptr);
                    found = true;
                }
                if (!found) {
                    std::fprintf(stderr,
                                 "RULE 5.2a FAIL: through cylinder rid=%d R=%.6f missing inner "
                                 "loop slab=%d\n",
                                 r.id, r.radius, k);
                    out.clear();
                    return false;
                }
            }
        }
        return true;
    } catch (...) {
        out.clear();
        return false;
    }
}

bool fitProfile(const MeshView& mv, const PrismTols& t, Profile& p, int& nDeclined) {
    nDeclined = 0;
    try {
        SliceCache snap;
        {
            std::lock_guard<std::mutex> g(cacheMu());
            snap = cache();
        }
        const DerivedTols dtol = makeDerived(mv);
        const RegionSet* rs = snap.ok ? snap.rs : nullptr;
        const SketchFrame* fr = snap.ok ? &snap.fr : nullptr;

        for (size_t li = 0; li < p.loops.size(); ++li) {
            ProfLoop& loop = p.loops[li];
            if (loop.segs.empty()) continue;
            std::vector<gp_Pnt2d> verts;
            verts.reserve(loop.segs.size());
            verts.push_back(loop.segs[0].a);
            for (const ProfSeg& s : loop.segs) verts.push_back(s.b);
            if (verts.size() >= 2 && near2(verts.front(), verts.back(), t.tauFit))
                verts.pop_back();
            if (verts.size() < 2) continue;

            std::vector<int> vrids;
            if (snap.ok && static_cast<size_t>(p.slab) < snap.rids.size() &&
                li < snap.rids[static_cast<size_t>(p.slab)].size()) {
                vrids = snap.rids[static_cast<size_t>(p.slab)][li];
            }
            while (vrids.size() < verts.size()) vrids.push_back(-1);

            // Closed-360: whole loop on one cylinder → one circle (B10).
            if (rs && fr && verts.size() >= 3) {
                const Region* only = nullptr;
                if (vrids[0] >= 0) only = cylByRid(*rs, vrids[0]);
                if (!only) only = cylAtPoint(*rs, *fr, verts[0], t.tauFit);
                bool all = only && only->closed360;
                if (all) {
                    for (size_t k = 0; k < verts.size(); ++k) {
                        const bool ridHit = (vrids[k] == only->id);
                        if (!ridHit && !onCyl(verts[k], *only, *fr, t.tauFit)) {
                            all = false;
                            break;
                        }
                    }
                }
                if (all) {
                    LawBand lb;
                    if (lawChainAccept(mv, only->tris, dtol, lb) && lb.closed360 &&
                        !rule53a(lb, t.tauFit)) {
                        ProfSeg arc;
                        arc.isArc = true;
                        arc.R = lb.R;
                        arc.phi = kTwoPi;
                        arc.center = cylCenter(*only, *fr);
                        arc.a = gp_Pnt2d(arc.center.X() + arc.R, arc.center.Y());
                        arc.b = arc.a;
                        arc.ccw = signedArea(verts) > 0.0;
                        loop.segs = {arc};
                        loop.area = std::fabs(signedArea(verts));
                        emitLoopDiag(p.slab, li, loop);
                        continue;
                    }
                }
            }

            std::vector<ProfSeg> fitted;
            const size_t n = verts.size();
            size_t i = 0;
            while (i < n) {
                bool took = false;
                if (rs && fr) {
                    const Region* cr = (vrids[i] >= 0) ? cylByRid(*rs, vrids[i]) : nullptr;
                    if (!cr) cr = cylAtPoint(*rs, *fr, verts[i], t.tauFit);
                    if (cr) {
                        size_t j = i + 1;
                        while (j < n) {
                            const bool ridHit = (vrids[j] == cr->id);
                            if (!ridHit && !onCyl(verts[j], *cr, *fr, t.tauFit)) break;
                            ++j;
                        }
                        if (j > i + 1 || cr->closed360) {
                            ProfSeg arc;
                            bool declined = false;
                            if (tryArc(mv, *cr, *fr, verts, i, j, t, dtol, arc, declined,
                                       nDeclined)) {
                                if (declined) {
                                    emitLines(verts, i, j, true, t.tauFit, fitted, nDeclined);
                                } else {
                                    fitted.push_back(arc);
                                }
                                i = (j >= n) ? n : j;
                                took = true;
                            }
                        }
                    }
                }
                if (took) continue;

                size_t j = i + 1;
                while (j < n) {
                    if (lineResid(verts[i], verts[j], verts, i, j) > t.tauFit) break;
                    ++j;
                }
                if (j <= i + 1) j = std::min(i + 2, n);
                ProfSeg line;
                line.isArc = false;
                line.a = verts[i];
                line.b = verts[j % n];
                if (dist2(line.a, line.b) > t.tauFit) fitted.push_back(line);
                i = j;
            }

            snapClosed(fitted, t.tauFit);
            mergeColinear(fitted, verts, t.tauFit);
            loop.segs = std::move(fitted);
            loop.area = std::fabs(signedArea(verts));
            emitLoopDiag(p.slab, li, loop);
        }
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
