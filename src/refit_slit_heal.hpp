// Zero-area fins / out-and-back slits in a certified face's boundary (D-fc-5.11 F3,
// the builder extension of the D-fc-1.8 heal). A fin is a region whose whole
// boundary is one BIGON -- two open chains joining the same two vertices a and b,
// both shared with one host region -- and whose two built edge runs coincide, so the
// face it would ship encloses no area. The host carries the same pair either as a
// two-edge inner loop or as a consecutive a -> b -> a excursion of a loop: two
// coincident edges, which BRepCheck reads as SelfIntersectingWire. The builder
// absorbs the fin into its host: the host's loops lose the pair, the fin ships no
// face, and every edge of the pair leaves the shell with both of its uses.
// No threshold: coincidence is read at the tolerance the edges themselves carry.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_SLIT_HEAL_HPP
#define STL2STEP_REFIT_SLIT_HEAL_HPP

#include "refit.hpp"

#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace stl2step {
namespace refit {

struct FinBigon {
    int fin = -1, host = -1;
    int c1 = -1, c2 = -1;  // the bigon's two chain indices, in the fin's loop order
};

// The fin's topology, no geometry read: F has one loop of exactly two open chains;
// each chain is shared with the same other region (never an island, never F itself)
// and both join the same two distinct terminal vertices.
inline bool finBigonTopology(const Region& F, const RegionSet& rs, FinBigon& out) {
    if (F.loops.size() != 1) return false;
    const Loop& lp = F.loops[0];
    if (lp.role != LoopRole::Outer || lp.chainIdx.size() != 2 || lp.reversed.size() != 2)
        return false;
    int host = -1;
    int ends[2][2] = {{-1, -1}, {-1, -1}};
    for (int k = 0; k < 2; ++k) {
        const int ci = lp.chainIdx[(size_t)k];
        if (ci < 0 || (size_t)ci >= rs.chains.size()) return false;
        const BoundaryChain& ch = rs.chains[(size_t)ci];
        if (ch.closedLoop || ch.meshVerts.size() < 2) return false;
        int other = -1;
        if (ch.regA == F.id) other = ch.regB;
        else if (ch.regB == F.id) other = ch.regA;
        else return false;
        if (other < 0 || other == F.id) return false;
        if (host >= 0 && other != host) return false;
        host = other;
        ends[k][0] = ch.meshVerts.front();
        ends[k][1] = ch.meshVerts.back();
        if (ends[k][0] == ends[k][1]) return false;
    }
    const bool same = (ends[0][0] == ends[1][0] && ends[0][1] == ends[1][1]) ||
                      (ends[0][0] == ends[1][1] && ends[0][1] == ends[1][0]);
    if (!same) return false;
    out.fin = F.id;
    out.host = host;
    out.c1 = lp.chainIdx[0];
    out.c2 = lp.chainIdx[1];
    return true;
}

// The host's loops with the bigon's chains c1, c2 removed. An inner loop made of
// exactly the pair is dropped. A pair that is cyclically consecutive in a longer loop
// is cut: both chains join a and b, so the loop runs a -> b -> a there and its
// neighbours still meet at a. False, loops unchanged, when either chain is missing,
// the pair is split across loops or not consecutive, or a two-chain loop that is not
// an inner loop would be dropped.
inline bool cutFinFromLoops(std::vector<Loop>& loops, int c1, int c2) {
    if (c1 < 0 || c2 < 0 || c1 == c2) return false;
    const size_t npos = (size_t)-1;
    std::vector<Loop> out;
    out.reserve(loops.size());
    bool cut = false;
    for (const Loop& lp : loops) {
        const size_t n = lp.chainIdx.size();
        if (lp.reversed.size() != n) return false;
        size_t p1 = npos, p2 = npos;
        for (size_t k = 0; k < n; ++k) {
            if (lp.chainIdx[k] == c1) {
                if (p1 != npos) return false;
                p1 = k;
            } else if (lp.chainIdx[k] == c2) {
                if (p2 != npos) return false;
                p2 = k;
            }
        }
        if (p1 == npos && p2 == npos) {
            out.push_back(lp);
            continue;
        }
        if (p1 == npos || p2 == npos || cut) return false;
        if (n == 2) {
            if (lp.role != LoopRole::Inner) return false;
            cut = true;
            continue;
        }
        if ((p1 + 1) % n != p2 && (p2 + 1) % n != p1) return false;
        Loop nl;
        nl.role = lp.role;
        nl.chainIdx.reserve(n - 2);
        nl.reversed.reserve(n - 2);
        for (size_t k = 0; k < n; ++k) {
            if (k == p1 || k == p2) continue;
            nl.chainIdx.push_back(lp.chainIdx[k]);
            nl.reversed.push_back(lp.reversed[k]);
        }
        out.push_back(std::move(nl));
        cut = true;
    }
    if (!cut) return false;
    loops = std::move(out);
    return true;
}

// p lies on one edge of `run` within the larger of tolP and that edge's own
// tolerance (its interior foot point or either end).
inline bool pointOnEdgeRun(const gp_Pnt& p, double tolP, const std::vector<TopoDS_Edge>& run) {
    for (const TopoDS_Edge& e : run) {
        if (e.IsNull()) continue;
        Standard_Real f = 0, l = 0;
        const Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull()) continue;
        const double tol = std::max(tolP, BRep_Tool::Tolerance(e));
        double d = std::min(p.Distance(c->Value(f)), p.Distance(c->Value(l)));
        try {
            GeomAPI_ProjectPointOnCurve pr(p, c, f, l);
            if (pr.NbPoints() > 0) d = std::min(d, (double)pr.LowerDistance());
        } catch (const Standard_Failure&) {
        }
        if (d <= tol) return true;
    }
    return false;
}

// Zero area in the B-Rep, for two runs that join the same two vertices (the
// bigon): the parameter midpoint of every edge of each run lies on the other run
// within the tolerance the two edges carry. The runs are surface-intersection
// curves of one region pair (line, circle or ellipse arcs between common ends):
// one interior point tells the same arc from the other arc of a closed curve.
inline bool edgeRunsCoincide(const std::vector<TopoDS_Edge>& a, const std::vector<TopoDS_Edge>& b) {
    if (a.empty() || b.empty()) return false;
    auto covered = [](const std::vector<TopoDS_Edge>& from, const std::vector<TopoDS_Edge>& onto) {
        for (const TopoDS_Edge& e : from) {
            if (e.IsNull()) return false;
            Standard_Real f = 0, l = 0;
            const Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
            if (c.IsNull()) return false;
            if (!pointOnEdgeRun(c->Value(0.5 * (f + l)), BRep_Tool::Tolerance(e), onto))
                return false;
        }
        return true;
    };
    return covered(a, b) && covered(b, a);
}

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_REFIT_SLIT_HEAL_HPP
