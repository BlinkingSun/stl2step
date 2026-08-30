// Stage P — buildPrismSolid: one prism per slab, holes in the profile, union.
// Caps are outputs (RULE 5.2c). Deterministic fuse order. SPDX-License-Identifier: MIT

#include "refit_prism.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepLib.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <ElCLib.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GProp_GProps.hxx>
#include <Precision.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {

void prismBindSketchOrigin(const gp_XYZ& o);
void prismResetStageFlags();
bool prismStagePUsed();
int prismPlatePathHits();
int prismLastReverted();
void prismNoteStageP(bool used, bool reverted, bool plate);

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

thread_local gp_XYZ tOrigin(0.0, 0.0, 0.0);
thread_local bool tHaveOrigin = false;
thread_local const char* tFail = "";
thread_local std::string tDxfDir;
thread_local double tUnifyV0 = 0.0, tUnifyV1 = 0.0;
thread_local int tUnifyP0 = 0, tUnifyC0 = 0, tUnifyP1 = 0, tUnifyC1 = 0;
thread_local int tUnifyG2 = 1, tUnifyG4 = 1;

void failAt(const char* why) {
    tFail = why;
    const char* e = std::getenv("STL2STEP_PRISM_DIAG");
    if (e && e[0] && e[0] != '0')
        std::fprintf(stderr, "DIAG_PRISMBUILD fail=%s\n", why);
}

std::atomic<int> gStageUsed{0};
std::atomic<int> gPlateHits{0};
std::atomic<int> gLastReverted{0};

bool prismDiagOn() {
    const char* e = std::getenv("STL2STEP_PRISM_DIAG");
    return e && e[0] && e[0] != '0';
}

template <class Fn>
void parallelFor(size_t n, Fn fn) {
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    if (n < 2 || hw < 2) {
        for (size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    const unsigned w = static_cast<unsigned>(std::min<size_t>(hw, n));
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(w);
    for (unsigned k = 0; k < w; ++k) {
        pool.emplace_back([&]() {
            for (;;) {
                const size_t i = next.fetch_add(1);
                if (i >= n) break;
                fn(i);
            }
        });
    }
    for (auto& th : pool) th.join();
}

struct Frame {
    gp_XYZ origin;
    gp_XYZ axis;
    gp_XYZ u;
    gp_XYZ v;
    double yRef = 0.0;
};

bool makeFrame(const PrismLevels& lv, Frame& fr) {
    fr.axis = lv.axis.XYZ();
    const double am = fr.axis.Modulus();
    if (!(am > 0.0) || !std::isfinite(am)) return false;
    fr.axis.Divide(am);
    // Hint pick: |n·e1| > this floor means the axis is nearly the X basis, so
    // use e2. Outer bound for a degenerate cross product, not a fit gate.
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
    fr.origin = tHaveOrigin ? tOrigin : gp_XYZ(0.0, 0.0, 0.0);
    fr.yRef = fr.axis.Dot(fr.origin);
    return true;
}

gp_Pnt to3(const Frame& fr, const gp_Pnt2d& p, double yLevel) {
    const gp_XYZ q = fr.origin + fr.u * p.X() + fr.v * p.Y() + fr.axis * (yLevel - fr.yRef);
    return gp_Pnt(q);
}

double dist2(const gp_Pnt2d& a, const gp_Pnt2d& b) {
    return std::hypot(a.X() - b.X(), a.Y() - b.Y());
}

double loopSignedArea(const ProfLoop& lp) {
    double a = 0.0;
    for (const ProfSeg& s : lp.segs)
        a += s.a.X() * s.b.Y() - s.b.X() * s.a.Y();
    return 0.5 * a;
}

bool isFullCircle(const ProfSeg& s) {
    if (!s.isArc || !(s.R > 0.0)) return false;
    if (s.phi >= kTwoPi * 0.5 && dist2(s.a, s.b) <= std::max(1e-9 * s.R, Precision::Confusion()))
        return true;
    return s.phi >= kTwoPi - 1e-9;
}

bool makeLineEdge(const gp_Pnt& a, const gp_Pnt& b, TopoDS_Edge& e) {
    if (a.Distance(b) <= Precision::Confusion()) return false;
    try {
        BRepBuilderAPI_MakeEdge me(a, b);
        if (!me.IsDone()) return false;
        e = me.Edge();
        return !e.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

gp_Pnt2d onCirc2(const ProfSeg& s, const gp_Pnt2d& p) {
    if (!s.isArc || !(s.R > 0.0)) return p;
    const double ang = std::atan2(p.Y() - s.center.Y(), p.X() - s.center.X());
    return gp_Pnt2d(s.center.X() + s.R * std::cos(ang), s.center.Y() + s.R * std::sin(ang));
}

bool makeArcEdge(const Frame& fr, const ProfSeg& s, const gp_Pnt2d& a2, const gp_Pnt2d& b2,
                 double yLevel, TopoDS_Edge& e) {
    if (!(s.R > 0.0)) return false;
    try {
        const gp_Pnt c = to3(fr, s.center, yLevel);
        const gp_Ax2 ax(c, gp_Dir(fr.axis), gp_Dir(fr.u));
        const gp_Circ circ(ax, s.R);
        Handle(Geom_Circle) gc = new Geom_Circle(circ);
        if (isFullCircle(s) || dist2(a2, b2) <= Precision::Confusion()) {
            BRepBuilderAPI_MakeEdge me(gc);
            if (!me.IsDone()) return false;
            e = me.Edge();
            return !e.IsNull();
        }
        const double u0 = std::atan2(a2.Y() - s.center.Y(), a2.X() - s.center.X());
        double sweep = s.phi;
        if (!(sweep > 0.0)) {
            double u1guess = std::atan2(b2.Y() - s.center.Y(), b2.X() - s.center.X());
            sweep = u1guess - u0;
            if (s.ccw) {
                while (sweep <= 0.0) sweep += kTwoPi;
            } else {
                while (sweep >= 0.0) sweep -= kTwoPi;
                sweep = -sweep;
            }
        }
        if (!(sweep > 0.0)) return false;
        const double u1 = s.ccw ? (u0 + sweep) : (u0 - sweep);
        const gp_Pnt pA = ElCLib::Value(u0, circ);
        const gp_Pnt pM = ElCLib::Value(0.5 * (u0 + u1), circ);
        const gp_Pnt pB = ElCLib::Value(u1, circ);
        GC_MakeArcOfCircle mk(pA, pM, pB);
        if (!mk.IsDone()) {
            BRepBuilderAPI_MakeEdge me(gc, std::min(u0, u1), std::max(u0, u1));
            if (!me.IsDone()) return false;
            e = me.Edge();
            return !e.IsNull();
        }
        BRepBuilderAPI_MakeEdge me(mk.Value());
        if (!me.IsDone()) return false;
        e = me.Edge();
        return !e.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

void mergeSameCircleArcs(ProfLoop& lp) {
    if (lp.segs.size() < 2) return;
    std::vector<ProfSeg> out;
    out.reserve(lp.segs.size());
    ProfSeg acc = lp.segs[0];
    auto sameCirc = [](const ProfSeg& a, const ProfSeg& b) {
        if (!a.isArc || !b.isArc || !(a.R > 0.0) || !(b.R > 0.0)) return false;
        if (std::fabs(a.R - b.R) > std::max(1e-6 * a.R, Precision::Confusion()))
            return false;
        return dist2(a.center, b.center) <= std::max(1e-6 * a.R, Precision::Confusion());
    };
    for (size_t i = 1; i < lp.segs.size(); ++i) {
        const ProfSeg& s = lp.segs[i];
        if (sameCirc(acc, s) && acc.ccw == s.ccw) {
            acc.b = s.b;
            acc.phi += s.phi;
            continue;
        }
        out.push_back(acc);
        acc = s;
    }
    if (!out.empty() && sameCirc(acc, out.front()) && acc.ccw == out.front().ccw) {
        out.front().a = acc.a;
        out.front().phi += acc.phi;
    } else {
        out.push_back(acc);
    }
    lp.segs = std::move(out);
}

void snapLineEndsToArcs(ProfLoop& lp) {
    const size_t n = lp.segs.size();
    if (n < 2) return;
    for (size_t i = 0; i < n; ++i) {
        ProfSeg& s = lp.segs[i];
        if (s.isArc && s.R > 0.0) {
            s.a = onCirc2(s, s.a);
            s.b = onCirc2(s, s.b);
        }
    }
    // Snap only LINE ends onto neighboring arcs. Do not copy an arc's
    // endpoint onto a different circle — that opens a 3D gap and a bridge.
    for (size_t i = 0; i < n; ++i) {
        ProfSeg& s = lp.segs[i];
        if (s.isArc) continue;
        ProfSeg& nxt = lp.segs[(i + 1) % n];
        ProfSeg& prev = lp.segs[(i + n - 1) % n];
        if (prev.isArc && prev.R > 0.0) {
            s.a = onCirc2(prev, s.a);
            prev.b = s.a;
        }
        if (nxt.isArc && nxt.R > 0.0) {
            s.b = onCirc2(nxt, s.b);
            nxt.a = s.b;
        }
    }
}

bool buildLoopWire(const Frame& fr, const ProfLoop& lpIn, double yLevel, bool wantCcw,
                   TopoDS_Wire& w) {
    ProfLoop lp = lpIn;
    mergeSameCircleArcs(lp);
    snapLineEndsToArcs(lp);
    if (lp.segs.empty()) return false;
    try {
        if (lp.segs.size() == 1 && isFullCircle(lp.segs[0])) {
            ProfSeg circ = lp.segs[0];
            circ.a = gp_Pnt2d(circ.center.X() + circ.R, circ.center.Y());
            circ.b = circ.a;
            circ.phi = kTwoPi;
            TopoDS_Edge e;
            if (!makeArcEdge(fr, circ, circ.a, circ.b, yLevel, e)) return false;
            BRepBuilderAPI_MakeWire mw(e);
            if (!mw.IsDone()) return false;
            w = mw.Wire();
            // Geom_Circle is CCW about the axis; a hole wants the opposite sense.
            if (!wantCcw) w.Reverse();
            w.Closed(Standard_True);
            return true;
        }
        if (lp.segs.size() == 1 && lp.segs[0].isArc && !isFullCircle(lp.segs[0]))
            return false;
        const size_t n = lp.segs.size();
        BRepBuilderAPI_MakeWire mw;
        gp_Pnt lastEnd;
        bool haveLast = false;
        gp_Pnt firstStart;
        int nBridge = 0;
        for (size_t i = 0; i < n; ++i) {
            const ProfSeg& s = lp.segs[i];
            gp_Pnt2d a2 = s.a;
            gp_Pnt2d b2 = s.b;
            if (s.isArc && s.R > 0.0) {
                a2 = onCirc2(s, s.a);
                b2 = onCirc2(s, s.b);
            }
            TopoDS_Edge e;
            bool ok = false;
            if (s.isArc && s.R > 0.0 && !s.declinedAmbiguous)
                ok = makeArcEdge(fr, s, a2, b2, yLevel, e);
            if (!ok) {
                gp_Pnt pA = to3(fr, a2, yLevel);
                gp_Pnt pB = to3(fr, b2, yLevel);
                if (haveLast) pA = lastEnd;
                ok = makeLineEdge(pA, pB, e);
            }
            if (!ok) return false;
            TopoDS_Vertex va, vb;
            TopExp::Vertices(e, va, vb, Standard_True);
            const gp_Pnt e0 = BRep_Tool::Pnt(va);
            const gp_Pnt e1 = BRep_Tool::Pnt(vb);
            if (haveLast && lastEnd.Distance(e0) > Precision::Confusion()) {
                TopoDS_Edge bridge;
                if (!makeLineEdge(lastEnd, e0, bridge)) return false;
                mw.Add(bridge);
                if (!mw.IsDone()) return false;
                ++nBridge;
                if (prismDiagOn())
                    std::fprintf(stderr, "DIAG_PRISMBUILD gap=%.6e i=%zu\n",
                                 lastEnd.Distance(e0), i);
            }
            mw.Add(e);
            if (!mw.IsDone()) return false;
            if (!haveLast) firstStart = e0;
            lastEnd = e1;
            haveLast = true;
        }
        if (haveLast && lastEnd.Distance(firstStart) > Precision::Confusion()) {
            TopoDS_Edge bridge;
            if (!makeLineEdge(lastEnd, firstStart, bridge)) return false;
            mw.Add(bridge);
            if (!mw.IsDone()) return false;
            ++nBridge;
            if (prismDiagOn())
                std::fprintf(stderr, "DIAG_PRISMBUILD close-gap=%.6e\n",
                             lastEnd.Distance(firstStart));
        }
        if (prismDiagOn() && nBridge)
            std::fprintf(stderr, "DIAG_PRISMBUILD nBridge=%d nSeg=%zu\n", nBridge, n);
        if (!mw.IsDone()) return false;
        w = mw.Wire();
        const double sa = loopSignedArea(lp);
        const bool ccw = sa > 0.0;
        if (ccw != wantCcw) w.Reverse();
        w.Closed(Standard_True);
        return !w.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool buildSlabFace(const Frame& fr, const Profile& p, double yLevel, TopoDS_Face& face) {
    if (p.loops.empty() || !p.loops[0].outer) {
        failAt("no-outer");
        return false;
    }
    try {
        const gp_Pnt o = to3(fr, gp_Pnt2d(0.0, 0.0), yLevel);
        const gp_Pln pln(o, gp_Dir(fr.axis));
        TopoDS_Wire outer;
        if (!buildLoopWire(fr, p.loops[0], yLevel, true, outer)) {
            failAt("outer-wire");
            return false;
        }
        Handle(Geom_Plane) surf = new Geom_Plane(pln);
        BRep_Builder B;
        TopoDS_Face f;
        B.MakeFace(f, surf, Precision::Confusion());
        B.Add(f, outer);
        for (size_t i = 1; i < p.loops.size(); ++i) {
            if (p.loops[i].outer) {
                failAt("extra-outer");
                return false;
            }
            // 5.2a may inject a full circle for a cylinder that is already a
            // scallop on the outer. Cutting it again is an invalid hole.
            bool onOuter = false;
            for (const ProfSeg& hs : p.loops[i].segs) {
                if (!hs.isArc || !(hs.R > 0.0)) continue;
                for (const ProfSeg& os : p.loops[0].segs) {
                    if (os.isArc && os.R > 0.0 &&
                        std::fabs(os.R - hs.R) <= std::max(1e-6 * hs.R, Precision::Confusion()) &&
                        dist2(os.center, hs.center) <=
                            std::max(1e-6 * hs.R, Precision::Confusion())) {
                        onOuter = true;
                        break;
                    }
                }
                if (onOuter) break;
            }
            if (onOuter) continue;
            TopoDS_Wire inner;
            if (!buildLoopWire(fr, p.loops[i], yLevel, false, inner)) {
                // Partial-arc leftovers of an outer scallop are skipped.
                continue;
            }
            B.Add(f, inner);
        }
        try {
            BRepLib::BuildCurves3d(f);
        } catch (const Standard_Failure&) {
        }
        face = f;
        return !face.IsNull();
    } catch (const Standard_Failure&) {
        failAt("slab-face");
        return false;
    }
}

bool buildOnePrism(const Frame& fr, const Profile& p, const PrismLevels& lv,
                   TopoDS_Shape& solid) {
    if (p.slab < 0 || static_cast<size_t>(p.slab) + 1 >= lv.y.size()) return false;
    const double y0 = lv.y[static_cast<size_t>(p.slab)];
    const double y1 = lv.y[static_cast<size_t>(p.slab) + 1];
    const double h = y1 - y0;
    if (!(h > 0.0) || !std::isfinite(h)) return false;
    TopoDS_Face face;
    if (!buildSlabFace(fr, p, y0, face)) return false;
    try {
        const gp_Vec vec = gp_Vec(gp_Dir(fr.axis)) * h;
        BRepPrimAPI_MakePrism mk(face, vec, Standard_False);
        if (!mk.IsDone()) return false;
        solid = mk.Shape();
        return !solid.IsNull();
    } catch (const Standard_Failure&) {
        return false;
    }
}

double shapeVolOf(const TopoDS_Shape& s) {
    if (s.IsNull()) return 0.0;
    try {
        GProp_GProps gp;
        BRepGProp::VolumeProperties(s, gp);
        return gp.Mass();
    } catch (const Standard_Failure&) {
        return 0.0;
    }
}

void countSurf(const TopoDS_Shape& s, int& nP, int& nC) {
    nP = 0;
    nC = 0;
    if (s.IsNull()) return;
    for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) {
        try {
            BRepAdaptor_Surface sa(TopoDS::Face(fx.Current()), Standard_False);
            if (sa.GetType() == GeomAbs_Plane) ++nP;
            else if (sa.GetType() == GeomAbs_Cylinder) ++nC;
        } catch (const Standard_Failure&) {
        }
    }
}

bool shapeClosedOf(const TopoDS_Shape& s) {
    if (s.IsNull()) return false;
    try {
        if (BRep_Tool::IsClosed(s)) return true;
        for (TopExp_Explorer sx(s, TopAbs_SHELL); sx.More(); sx.Next()) {
            if (BRep_Tool::IsClosed(TopoDS::Shell(sx.Current()))) return true;
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

bool shapeValidOf(const TopoDS_Shape& s) {
    if (s.IsNull()) return false;
    try {
        BRepCheck_Analyzer an(s, Standard_True);
        if (an.IsValid() == Standard_True) return true;
        for (TopExp_Explorer sx(s, TopAbs_SHELL); sx.More(); sx.Next()) {
            BRepCheck_Analyzer as(sx.Current(), Standard_True);
            if (as.IsValid() == Standard_True) return true;
        }
    } catch (const Standard_Failure&) {
    }
    return false;
}

double meshVolOf(const MeshView& mv) {
    double v = 0.0;
    if (!mv.pts || !mv.tris || !mv.compTris) return 0.0;
    for (size_t k = 0; k < mv.nTri; ++k) {
        const int gt = mv.compTris[k];
        if (gt < 0) continue;
        const gp_XYZ a = mv.pts[mv.tris[gt][0]];
        const gp_XYZ b = mv.pts[mv.tris[gt][1]];
        const gp_XYZ c = mv.pts[mv.tris[gt][2]];
        v += a.Dot(b.Crossed(c)) / 6.0;
    }
    return v;
}

int nSidesForR(double R, const RegionSet& rs) {
    int best = 0;
    double bestD = 1e300;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder || !(r.radius > 0.0) || r.nSides < 3)
            continue;
        const double assoc =
            std::max({Precision::Confusion(), 4.0 * r.maxVertexDev, 4.0 * r.chordSagitta});
        const double d = std::fabs(r.radius - R);
        if (d <= assoc && d < bestD) {
            bestD = d;
            best = r.nSides;
        }
    }
    return best;
}

double segSweep(const ProfSeg& s) {
    if (!s.isArc || !(s.R > 0.0)) return 0.0;
    if (s.phi > 0.0) return s.phi;
    const double u0 = std::atan2(s.a.Y() - s.center.Y(), s.a.X() - s.center.X());
    const double u1 = std::atan2(s.b.Y() - s.center.Y(), s.b.X() - s.center.X());
    double sw = u1 - u0;
    if (s.ccw) {
        while (sw <= 0.0) sw += kTwoPi;
    } else {
        while (sw >= 0.0) sw -= kTwoPi;
        sw = -sw;
    }
    return sw;
}

// D9 §1: defect(loop) = Σ_strips (R²/2)(θ − sin θ)·h_slab from fitted profiles.
double loopArcDefect(const ProfLoop& lp, double h, const RegionSet& rs) {
    if (!(h > 0.0) || !std::isfinite(h)) return 0.0;
    double d = 0.0;
    for (const ProfSeg& s : lp.segs) {
        if (!s.isArc || !(s.R > 0.0) || s.declinedAmbiguous) continue;
        const double theta = segSweep(s);
        if (!(theta > 0.0) || !std::isfinite(theta)) continue;
        const int n = nSidesForR(s.R, rs);
        if (n >= 3) {
            const double g = kTwoPi / static_cast<double>(n);
            const double nStrips = std::max(1.0, theta / g);
            d += nStrips * 0.5 * s.R * s.R * (g - std::sin(g)) * h;
        } else if (theta < kTwoPi - 1e-9) {
            d += 0.5 * s.R * s.R * (theta - std::sin(theta)) * h;
        }
    }
    return d;
}

void profileArcDefects(const std::vector<Profile>& profs, const PrismLevels& lv,
                       const RegionSet& rs, double& dSigned, double& dAbs) {
    dSigned = 0.0;
    dAbs = 0.0;
    for (const Profile& p : profs) {
        if (p.slab < 0 || static_cast<size_t>(p.slab) + 1 >= lv.y.size()) continue;
        const double h = lv.y[static_cast<size_t>(p.slab) + 1] - lv.y[static_cast<size_t>(p.slab)];
        for (const ProfLoop& lp : p.loops) {
            const double def = loopArcDefect(lp, h, rs);
            // σ = +1 inner (hole: mesh holds MORE), −1 outer (boss: mesh holds LESS).
            const double sig = lp.outer ? -1.0 : 1.0;
            dSigned += sig * def;
            dAbs += std::fabs(def);
        }
    }
}

// D9 §2: UnifySameDomain once, G1 ceiling = Precision::Confusion(). No ShapeFix.
bool unifySameOnce(TopoDS_Shape& s) {
    tUnifyV0 = 0.0;
    tUnifyV1 = 0.0;
    tUnifyP0 = 0;
    tUnifyC0 = 0;
    tUnifyP1 = 0;
    tUnifyC1 = 0;
    tUnifyG2 = 1;
    tUnifyG4 = 1;
    if (s.IsNull()) return false;
    tUnifyV0 = shapeVolOf(s);
    countSurf(s, tUnifyP0, tUnifyC0);
    tUnifyV1 = tUnifyV0;
    tUnifyP1 = tUnifyP0;
    tUnifyC1 = tUnifyC0;
    try {
        ShapeUpgrade_UnifySameDomain usd(s, Standard_True, Standard_True, Standard_False);
        usd.SetLinearTolerance(Precision::Confusion());
        usd.SetAngularTolerance(Precision::Confusion());
        usd.SetSafeInputMode(Standard_True);
        usd.Build();
        const TopoDS_Shape u = usd.Shape();
        if (u.IsNull()) return true;
        const double v1 = shapeVolOf(u);
        int nP = 0, nC = 0;
        countSurf(u, nP, nC);
        const double vFloor = 1e-6 * std::fabs(tUnifyV0);
        if (std::fabs(v1 - tUnifyV0) > vFloor) {
            tUnifyG2 = 0;
            return true;
        }
        if (nC != tUnifyC0) {
            tUnifyG4 = 0;
            return true;
        }
        s = u;
        tUnifyV1 = v1;
        tUnifyP1 = nP;
        tUnifyC1 = nC;
        return true;
    } catch (const Standard_Failure&) {
        return true;
    }
}

}  // namespace

void prismBindSketchOrigin(const gp_XYZ& o) {
    tOrigin = o;
    tHaveOrigin = true;
}

void prismResetStageFlags() {
    gStageUsed.store(0);
    gPlateHits.store(0);
    gLastReverted.store(0);
}

bool prismStagePUsed() { return gStageUsed.load() != 0; }

int prismPlatePathHits() { return gPlateHits.load(); }

int prismLastReverted() { return gLastReverted.load(); }

void prismNoteStageP(bool used, bool reverted, bool plate) {
    if (used) gStageUsed.store(1);
    gLastReverted.store(reverted ? 1 : 0);
    if (plate) gPlateHits.fetch_add(1);
    else if (used && !reverted) gPlateHits.store(0);
}

void prismBindDxfDir(const std::string& d) { tDxfDir = d; }

bool buildPrismSolid(const std::vector<Profile>& profs, const PrismLevels& lv,
                     TopoDS_Shape& out) {
    out = TopoDS_Shape();
    try {
        if (!lv.ok || lv.y.size() < 2) {
            failAt("levels");
            return false;
        }
        const int nSlab = static_cast<int>(lv.y.size() - 1);
        if (nSlab < 1) {
            failAt("nslab");
            return false;
        }

        std::vector<const Profile*> bySlab(static_cast<size_t>(nSlab), nullptr);
        for (const Profile& p : profs) {
            if (p.slab < 0 || p.slab >= nSlab) {
                failAt("slab-id");
                return false;
            }
            if (bySlab[static_cast<size_t>(p.slab)]) {
                failAt("dup-slab");
                return false;
            }
            if (p.loops.empty() || !p.loops[0].outer) {
                failAt("profile-outer");
                return false;
            }
            bySlab[static_cast<size_t>(p.slab)] = &p;
        }
        for (const Profile* p : bySlab)
            if (!p) {
                failAt("missing-slab");
                return false;
            }

        Frame fr;
        if (!makeFrame(lv, fr)) {
            failAt("frame");
            return false;
        }

        std::vector<TopoDS_Shape> slabs(static_cast<size_t>(nSlab));
        // OCCT BRep construction is not re-entrant; slabs are built in slab order.
        // Parallel work is the 2D prep already done by P2 and the later census.
        for (int k = 0; k < nSlab; ++k) {
            try {
                if (!buildOnePrism(fr, *bySlab[static_cast<size_t>(k)], lv,
                                   slabs[static_cast<size_t>(k)])) {
                    failAt("slab-prism");
                    return false;
                }
                if (prismDiagOn()) {
                    int nf = 0;
                    for (TopExp_Explorer fx(slabs[static_cast<size_t>(k)], TopAbs_FACE);
                         fx.More(); fx.Next())
                        ++nf;
                    double sv = 0.0;
                    int okv = 0;
                    try {
                        GProp_GProps gp;
                        BRepGProp::VolumeProperties(slabs[static_cast<size_t>(k)], gp);
                        sv = gp.Mass();
                    } catch (...) {
                    }
                    try {
                        BRepCheck_Analyzer an(slabs[static_cast<size_t>(k)], Standard_True);
                        okv = an.IsValid() ? 1 : 0;
                    } catch (...) {
                    }
                    std::fprintf(stderr, "DIAG_PRISMBUILD slab=%d faces=%d vol=%.6f valid=%d\n",
                                 k, nf, sv, okv);
                }
            } catch (...) {
                failAt("slab-throw");
                return false;
            }
        }

        double hMin = 0.0;
        for (int k = 0; k < nSlab; ++k) {
            const double h = lv.y[static_cast<size_t>(k) + 1] - lv.y[static_cast<size_t>(k)];
            if (h > 0.0 && (hMin <= 0.0 || h < hMin)) hMin = h;
        }
        const double linTol =
            std::max(Precision::Confusion(), (hMin > 0.0 ? 1e-6 * hMin : 0.0));

        auto countF = [](const TopoDS_Shape& s) {
            int n = 0;
            if (s.IsNull()) return 0;
            for (TopExp_Explorer fx(s, TopAbs_FACE); fx.More(); fx.Next()) ++n;
            return n;
        };
        TopoDS_Shape acc = slabs[0];
        for (int k = 1; k < nSlab; ++k) {
            BRepAlgoAPI_Fuse fuse(acc, slabs[static_cast<size_t>(k)]);
            fuse.SetRunParallel(Standard_False);
            fuse.SetFuzzyValue(std::max(linTol, Precision::Confusion()));
            fuse.SetGlue(BOPAlgo_GlueShift);
            fuse.Build();
            if (!fuse.IsDone() || fuse.HasErrors()) {
                // GlueShift can empty a slightly-overlapping pair; retry plain fuse.
                BRepAlgoAPI_Fuse fuse2(acc, slabs[static_cast<size_t>(k)]);
                fuse2.SetRunParallel(Standard_False);
                fuse2.SetFuzzyValue(std::max(linTol, Precision::Confusion()));
                fuse2.Build();
                if (!fuse2.IsDone()) {
                    failAt("fuse");
                    return false;
                }
                acc = fuse2.Shape();
            } else {
                acc = fuse.Shape();
            }
            if (acc.IsNull() || countF(acc) < 1) {
                failAt("fuse-null");
                return false;
            }
            if (prismDiagOn())
                std::fprintf(stderr, "DIAG_PRISMBUILD fuse-k=%d faces=%d\n", k, countF(acc));
        }
        if (!unifySameOnce(acc)) {
            failAt("unify");
            return false;
        }
        if (prismDiagOn())
            std::fprintf(stderr,
                         "DIAG_PRISMBUILD unify faces=%d G1=1 G2=%d G4=%d "
                         "V_before=%.6f V_after=%.6f planes=%d->%d cyls=%d->%d\n",
                         countF(acc), tUnifyG2, tUnifyG4, tUnifyV0, tUnifyV1,
                         tUnifyP0, tUnifyP1, tUnifyC0, tUnifyC1);
        int nF = 0;
        for (TopExp_Explorer fx(acc, TopAbs_FACE); fx.More(); fx.Next()) ++nF;
        if (nF < 1) {
            failAt("no-faces");
            out = TopoDS_Shape();
            return false;
        }
        out = acc;
        return !out.IsNull();
    } catch (...) {
        out = TopoDS_Shape();
        return false;
    }
}

bool tryStageP(const MeshView& mv, RegionSet& rs, std::vector<TopoDS_Face>& out) {
    int nFaces = 0, nPlanes = 0, nCyls = 0, nSlabs = 0;
    double vol = 0.0;
    int wt = 0, valid = 0, reverted = 1;
    auto emitBuild = [&]() {
        if (!prismDiagOn()) return;
        std::fprintf(stderr,
                     "DIAG_PRISMBUILD comp=%d slabs=%d faces=%d planes=%d cyls=%d "
                     "vol=%.6f watertight=%d valid=%d reverted=%d\n",
                     rs.compRoot >= 0 ? rs.compRoot : 0, nSlabs, nFaces, nPlanes, nCyls,
                     vol, wt, valid, reverted);
    };
    try {
        PrismTols pt;
        const PrismLevels lv = detectPrismatic(mv, rs, pt);
        if (!lv.ok || lv.y.size() < 2) return false;
        nSlabs = static_cast<int>(lv.y.size() - 1);

        gp_XYZ origin(0.0, 0.0, 0.0);
        if (!lv.capRegion.empty()) {
            for (const Region& r : rs.regions) {
                if (r.id == lv.capRegion[0]) {
                    origin = r.ax.Location().XYZ();
                    break;
                }
            }
        } else {
            for (const Region& r : rs.regions) {
                if (r.type == SurfType::Cylinder) {
                    origin = r.ax.Location().XYZ();
                    break;
                }
            }
        }
        prismBindSketchOrigin(origin);

        std::vector<Profile> profs;
        bool ready = sliceProfiles(mv, rs, lv, pt, profs);
        int nDecl = 0;
        if (ready) {
            for (Profile& p : profs) {
                int d = 0;
                if (!fitProfile(mv, pt, p, d)) {
                    ready = false;
                    break;
                }
                nDecl += d;
            }
        }
        (void)nDecl;
        // D8 §3.2: guarded, return-discarded; inert when dxfDir unset.
        struct {
            std::string dxfDir;
        } opt;
        opt.dxfDir = tDxfDir;
        if (!opt.dxfDir.empty()) {
            const int comp = rs.compRoot >= 0 ? rs.compRoot : 0;
            for (const Profile& prof : profs) {
                const std::string path = opt.dxfDir + "/profile-comp" +
                                         std::to_string(comp) + "-slab" +
                                         std::to_string(prof.slab) + ".dxf";
                (void)writeProfileDxf(prof, lv, path);
            }
        }
        const char* inj = std::getenv("STL2STEP_PRISM_INJECT_BAD");
        if (ready && inj && inj[0] && inj[0] != '0') {
            if (!profs.empty() && profs[0].loops.size() > 1)
                profs[0].loops.resize(1);
        }
        auto snapLvl = [&](double yv) -> int {
            int best = -1;
            double bestD = 0.0;
            for (size_t k = 0; k < lv.y.size(); ++k) {
                const double d = std::fabs(yv - lv.y[k]);
                if (d <= pt.tauLvl && (best < 0 || d < bestD)) {
                    best = static_cast<int>(k);
                    bestD = d;
                }
            }
            return best;
        };
        auto cylEnds = [&](const Region& r, double& lo, double& hi) {
            const gp_XYZ a = r.ax.Direction().XYZ();
            const gp_XYZ loc = r.ax.Location().XYZ();
            const gp_XYZ p0 = loc + a * r.vMin;
            const gp_XYZ p1 = loc + a * r.vMax;
            const gp_XYZ ah = lv.axis.XYZ();
            const double y0 = ah.Dot(p0);
            const double y1 = ah.Dot(p1);
            lo = std::min(y0, y1);
            hi = std::max(y0, y1);
        };
        if (ready) {
            for (const Region& r : rs.regions) {
                if (r.type != SurfType::Cylinder || !(r.radius > 0.0)) continue;
                double lo = 0.0, hi = 0.0;
                cylEnds(r, lo, hi);
                const int i0 = snapLvl(lo);
                const int i1 = snapLvl(hi);
                if (i0 < 0 || i1 < 0 || (i1 - i0) < 2) continue;
                for (int k = i0; k < i1; ++k) {
                    bool found = false;
                    if (k < 0 || static_cast<size_t>(k) >= profs.size()) {
                        ready = false;
                        break;
                    }
                    const Profile& p = profs[static_cast<size_t>(k)];
                    const double assoc =
                        std::max({pt.tauFit, 4.0 * r.maxVertexDev, 4.0 * r.chordSagitta});
                    for (const ProfLoop& lp : p.loops) {
                        if (lp.outer) continue;
                        for (const ProfSeg& s : lp.segs) {
                            if (s.isArc && s.R > 0.0 &&
                                std::fabs(s.R - r.radius) <= assoc) {
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                        if (lp.segs.size() >= 3) {
                            double cx = 0.0, cy = 0.0;
                            for (const ProfSeg& s : lp.segs) {
                                cx += s.a.X();
                                cy += s.a.Y();
                            }
                            cx /= static_cast<double>(lp.segs.size());
                            cy /= static_cast<double>(lp.segs.size());
                            double rmin = 1e300, rmax = 0.0;
                            for (const ProfSeg& s : lp.segs) {
                                const double d = std::hypot(s.a.X() - cx, s.a.Y() - cy);
                                rmin = std::min(rmin, d);
                                rmax = std::max(rmax, d);
                            }
                            if (rmax - rmin <= 8.0 * assoc &&
                                std::fabs(0.5 * (rmin + rmax) - r.radius) <= assoc)
                                found = true;
                        }
                    }
                    if (!found) {
                        ready = false;
                        break;
                    }
                }
                if (!ready) break;
            }
        }
        TopoDS_Shape solid;
        if (ready && buildPrismSolid(profs, lv, solid) && !solid.IsNull()) {
            std::vector<TopoDS_Face> faces;
            for (TopExp_Explorer ex(solid, TopAbs_FACE); ex.More(); ex.Next())
                faces.push_back(TopoDS::Face(ex.Current()));
            wt = shapeClosedOf(solid) ? 1 : 0;
            valid = shapeValidOf(solid) ? 1 : 0;
            vol = shapeVolOf(solid);
            nFaces = static_cast<int>(faces.size());
            countSurf(solid, nPlanes, nCyls);
            double dVolAbs = 0.0;
            for (const Region& r : rs.regions) dVolAbs += std::fabs(r.dVolPredicted);
            const double meshVol = std::fabs(meshVolOf(mv));
            const double budget = std::max(1e-4 * meshVol, 3.0 * dVolAbs);
            double dSigned = 0.0, dAbs = 0.0;
            profileArcDefects(profs, lv, rs, dSigned, dAbs);
            const double vRef = meshVol - dSigned;
            // D4 §3 K_arc reused as D9 K_prism: 5% margin. Not load-bearing.
            const double kPrism = 1.0 + 5.0 / 100.0;
            const double envelope = kPrism * dAbs;
            const bool cond1 = std::fabs(vol - vRef) <= budget;
            const bool cond2 = std::fabs(vol - meshVol) <= envelope;
            const bool volOk = cond1 && cond2;
            if (prismDiagOn()) {
                std::fprintf(stderr,
                             "DIAG_PRISMBUILD vol D_signed=%.6f D_abs=%.6f V_ref=%.6f "
                             "budget=%.6f envelope=%.6f |s-Vref|=%.6f |s-mesh|=%.6f "
                             "cond1=%d cond2=%d\n",
                             dSigned, dAbs, vRef, budget, envelope,
                             std::fabs(vol - vRef), std::fabs(vol - meshVol),
                             cond1 ? 1 : 0, cond2 ? 1 : 0);
                std::fprintf(stderr,
                             "DIAG_PRISMBUILD G1=1 G2=%d G3_planes=%d G4=%d "
                             "cyls=%d->%d V_before=%.6f V_after=%.6f\n",
                             tUnifyG2, nPlanes, tUnifyG4, tUnifyC0, tUnifyC1,
                             tUnifyV0, tUnifyV1);
            }
            if (wt && valid && volOk && !faces.empty()) {
                reverted = 0;
                rs.stats.planes = nPlanes;
                rs.stats.cylinders = nCyls;
                rs.stats.facetIslands = 0;
                rs.stats.facetTriangles = 0;
                prismNoteStageP(true, false, false);
                emitBuild();
                out = std::move(faces);
                return true;
            }
        }
        prismNoteStageP(false, true, false);
        emitBuild();
        return false;
    } catch (...) {
        prismNoteStageP(false, true, false);
        emitBuild();
        return false;
    }
}

}  // namespace refit
}  // namespace stl2step
