// src/refit_union_census.cpp -- D-140-8 U-R5..U-R8, U-R13: the same-surface
// union census (declared in refit_internal.hpp). Measurement only: every number is
// read off the welded mesh against the cylinder the caller certified; nothing
// is written back and no Region is touched.
//
// SPDX-License-Identifier: MIT

#include "refit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <gp_XYZ.hxx>

namespace stl2step {
namespace refit {
namespace {

bool haveTables(const MeshView& mv) {
    return mv.pts && mv.compVtx && mv.compEdges && mv.triEdges && mv.nTri > 0 && mv.nVtx > 0;
}

bool triInRange(const MeshView& mv, int t) {
    return t >= 0 && static_cast<std::size_t>(t) < mv.nTri;
}

bool vtxInRange(const MeshView& mv, int v) {
    return v >= 0 && static_cast<std::size_t>(v) < mv.nVtx;
}

gp_XYZ localPt(const MeshView& mv, int lv) { return mv.pts[mv.compVtx[lv]]; }

// distance from a point to the axis LINE
double axisDist(const gp_XYZ& o, const gp_XYZ& a, const gp_XYZ& p) {
    const gp_XYZ d = p - o;
    return (d - a * d.Dot(a)).Modulus();
}

// The three local vertices of a local triangle, read off its welded edge
// table: the union of its three edges' endpoints. Corner order is not needed
// by anything here, so triDirs is never consulted.
std::array<int, 3> triVerts(const MeshView& mv, int t) {
    int v[6];
    for (int s = 0; s < 3; s++) {
        const std::pair<int, int>& e = mv.compEdges[mv.triEdges[t][s]];
        v[2 * s] = e.first;
        v[2 * s + 1] = e.second;
    }
    std::sort(v, v + 6);
    std::array<int, 3> out{-1, -1, -1};
    int k = 0;
    for (int i = 0; i < 6 && k < 3; i++) {
        if (i > 0 && v[i] == v[i - 1]) continue;
        out[static_cast<std::size_t>(k++)] = v[i];
    }
    return out;
}

// local edge -> the (at most two) local triangles on it; -1 when absent
std::vector<std::array<int, 2>> buildEdgeTris(const MeshView& mv) {
    std::vector<std::array<int, 2>> et(mv.nEdge, std::array<int, 2>{-1, -1});
    for (int t = 0; t < static_cast<int>(mv.nTri); t++) {
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || static_cast<std::size_t>(e) >= mv.nEdge) continue;
            if (et[static_cast<std::size_t>(e)][0] < 0) et[static_cast<std::size_t>(e)][0] = t;
            else if (et[static_cast<std::size_t>(e)][1] < 0) et[static_cast<std::size_t>(e)][1] = t;
        }
    }
    return et;
}

// Edge-connected pieces of `ts` where adjacency crosses a shared mesh edge and
// both sides carry `mask`. Pieces are counted, not returned.
std::size_t edgePiecesOf(const MeshView& mv, const std::vector<std::array<int, 2>>& et,
                         const std::vector<int>& ts, const std::vector<char>& mask) {
    std::vector<char> seen(mv.nTri, 0);
    std::size_t n = 0;
    for (int t0 : ts) {
        if (!triInRange(mv, t0) || !mask[static_cast<std::size_t>(t0)] ||
            seen[static_cast<std::size_t>(t0)])
            continue;
        ++n;
        std::vector<int> stk{t0};
        seen[static_cast<std::size_t>(t0)] = 1;
        while (!stk.empty()) {
            const int x = stk.back();
            stk.pop_back();
            for (int s = 0; s < 3; s++) {
                const int e = mv.triEdges[x][s];
                if (e < 0 || static_cast<std::size_t>(e) >= mv.nEdge) continue;
                const std::array<int, 2>& pr = et[static_cast<std::size_t>(e)];
                const int u = (pr[0] == x) ? pr[1] : pr[0];
                if (u < 0 || !mask[static_cast<std::size_t>(u)] || seen[static_cast<std::size_t>(u)])
                    continue;
                seen[static_cast<std::size_t>(u)] = 1;
                stk.push_back(u);
            }
        }
    }
    return n;
}

}  // namespace

double claimChordSagitta(const MeshView& mv, const std::vector<int>& claim, const gp_Ax1& axis,
                         double R) {
    if (!haveTables(mv) || claim.empty()) return 0.0;
    const gp_XYZ o = axis.Location().XYZ();
    const gp_XYZ a = axis.Direction().XYZ();
    std::vector<char> seenE(mv.nEdge, 0);
    double sigma = 0.0;
    bool any = false;
    for (int t : claim) {
        if (!triInRange(mv, t)) continue;
        for (int s = 0; s < 3; s++) {
            const int e = mv.triEdges[t][s];
            if (e < 0 || static_cast<std::size_t>(e) >= mv.nEdge || seenE[static_cast<std::size_t>(e)])
                continue;
            seenE[static_cast<std::size_t>(e)] = 1;
            const std::pair<int, int>& pr = mv.compEdges[e];
            if (!vtxInRange(mv, pr.first) || !vtxInRange(mv, pr.second)) continue;
            const gp_XYZ mid = (localPt(mv, pr.first) + localPt(mv, pr.second)) * 0.5;
            const double dev = R - axisDist(o, a, mid);
            if (!any || dev > sigma) {
                sigma = dev;
                any = true;
            }
        }
    }
    return any ? sigma : 0.0;
}

UnionCensus unionCensus(const MeshView& mv, const std::vector<int>& claim, const gp_Ax1& axis,
                        double R) {
    UnionCensus uc;
    if (!haveTables(mv)) return uc;
    const std::size_t nT = mv.nTri;
    const std::size_t nV = mv.nVtx;
    const gp_XYZ o = axis.Location().XYZ();
    const gp_XYZ a = axis.Direction().XYZ();

    // membership, deduplicated and in range
    std::vector<char> inClaim(nT, 0);
    for (int t : claim)
        if (triInRange(mv, t)) inClaim[static_cast<std::size_t>(t)] = 1;
    std::vector<int> cl;
    for (std::size_t t = 0; t < nT; t++)
        if (inClaim[t]) cl.push_back(static_cast<int>(t));
    uc.nTri = cl.size();
    if (cl.empty()) return uc;

    uc.sigma = claimChordSagitta(mv, cl, axis, R);

    // mesh incidence
    std::vector<std::array<int, 3>> tv(nT);
    std::vector<std::vector<int>> triAtVtx(nV);
    for (std::size_t t = 0; t < nT; t++) {
        tv[t] = triVerts(mv, static_cast<int>(t));
        for (int v : tv[t])
            if (vtxInRange(mv, v)) triAtVtx[static_cast<std::size_t>(v)].push_back(static_cast<int>(t));
    }
    const std::vector<std::array<int, 2>> et = buildEdgeTris(mv);

    std::vector<char> claimVtx(nV, 0);
    for (int t : cl)
        for (int v : tv[static_cast<std::size_t>(t)])
            if (vtxInRange(mv, v)) claimVtx[static_cast<std::size_t>(v)] = 1;

    uc.edgePieces = edgePiecesOf(mv, et, cl, inClaim);

    // U-R6 -- punctures, to a fixed point. A vertex incident to the domain's
    // neighbourhood (a vertex of a non-domain triangle that touches a domain or
    // puncture vertex) is a puncture iff dist(v, S) <= sigma; an exterior vertex
    // is classified once and never revisited. The domain is the claim plus the
    // part of every puncture's mesh star that lies ON S: the star triangles
    // whose every vertex is a claim vertex or a puncture (U1's 24 on the plate).
    // A star triangle that reaches an exterior vertex leaves the surface and
    // stays out; so does a refused triangle in no puncture's star.
    std::vector<char> inDom = inClaim;
    std::vector<int> dom = cl;
    std::vector<signed char> vClass(nV, 0);  // 0 unseen, +1 puncture, -1 exterior
    std::vector<int> punct;
    auto onS = [&](int w) {
        return vtxInRange(mv, w) &&
               (claimVtx[static_cast<std::size_t>(w)] || vClass[static_cast<std::size_t>(w)] > 0);
    };
    for (bool changed = true; changed;) {
        changed = false;
        std::vector<char> domVtx(nV, 0);
        for (int t : dom)
            for (int v : tv[static_cast<std::size_t>(t)])
                if (vtxInRange(mv, v)) domVtx[static_cast<std::size_t>(v)] = 1;
        for (int v : punct) domVtx[static_cast<std::size_t>(v)] = 1;
        for (std::size_t v = 0; v < nV; v++) {
            if (!domVtx[v]) continue;
            for (int t : triAtVtx[v]) {
                if (inDom[static_cast<std::size_t>(t)]) continue;
                for (int w : tv[static_cast<std::size_t>(t)]) {
                    if (!vtxInRange(mv, w)) continue;
                    const std::size_t wi = static_cast<std::size_t>(w);
                    if (claimVtx[wi] || vClass[wi] != 0) continue;
                    const double dev = std::abs(axisDist(o, a, localPt(mv, w)) - R);
                    if (dev <= uc.sigma) {
                        vClass[wi] = 1;
                        punct.push_back(w);
                    } else {
                        vClass[wi] = -1;
                    }
                    changed = true;
                }
            }
        }
        for (int v : punct) {
            for (int s : triAtVtx[static_cast<std::size_t>(v)]) {
                if (inDom[static_cast<std::size_t>(s)]) continue;
                const std::array<int, 3>& w3 = tv[static_cast<std::size_t>(s)];
                if (!onS(w3[0]) || !onS(w3[1]) || !onS(w3[2])) continue;
                inDom[static_cast<std::size_t>(s)] = 1;
                dom.push_back(s);
                changed = true;
            }
        }
    }
    const std::size_t nPunct = punct.size();
    uc.punctures = nPunct;
    std::sort(dom.begin(), dom.end());
    uc.domainTris = dom.size();

    // U-R7 -- components of the domain graph under VERTEX adjacency, ascending
    // minimum triangle id.
    std::vector<int> compId(nT, -1);
    std::vector<std::vector<int>> comps;
    for (int t0 : dom) {
        if (compId[static_cast<std::size_t>(t0)] >= 0) continue;
        const int k = static_cast<int>(comps.size());
        std::vector<int> comp, stk{t0};
        compId[static_cast<std::size_t>(t0)] = k;
        while (!stk.empty()) {
            const int x = stk.back();
            stk.pop_back();
            comp.push_back(x);
            for (int v : tv[static_cast<std::size_t>(x)]) {
                if (!vtxInRange(mv, v)) continue;
                for (int u : triAtVtx[static_cast<std::size_t>(v)]) {
                    if (!inDom[static_cast<std::size_t>(u)] || compId[static_cast<std::size_t>(u)] >= 0)
                        continue;
                    compId[static_cast<std::size_t>(u)] = k;
                    stk.push_back(u);
                }
            }
        }
        std::sort(comp.begin(), comp.end());
        comps.push_back(std::move(comp));
    }
    uc.domainFaces = comps.size();
    for (const std::vector<int>& comp : comps) {
        std::vector<int> ownClaim;
        for (int t : comp)
            if (inClaim[static_cast<std::size_t>(t)]) ownClaim.push_back(t);
        uc.sizes.push_back(ownClaim.size());
        uc.pieces.push_back(edgePiecesOf(mv, et, ownClaim, inClaim));
    }

    // U-R8 -- pinch vertices: the domain triangles incident to v form two or
    // more edge-connected fans (fan adjacency crosses an edge that contains v).
    // Every triangle incident to v lies in one component by construction.
    std::size_t nPinch = 0;
    for (std::size_t v = 0; v < nV; v++) {
        std::vector<int> inc;
        for (int t : triAtVtx[v])
            if (inDom[static_cast<std::size_t>(t)]) inc.push_back(t);
        if (inc.size() < 2) continue;
        std::vector<int> parent(inc.size());
        for (std::size_t i = 0; i < inc.size(); i++) parent[i] = static_cast<int>(i);
        auto find = [&](int x) {
            while (parent[static_cast<std::size_t>(x)] != x) {
                parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
                x = parent[static_cast<std::size_t>(x)];
            }
            return x;
        };
        auto indexOf = [&](int t) {
            const auto it = std::find(inc.begin(), inc.end(), t);
            return it == inc.end() ? -1 : static_cast<int>(it - inc.begin());
        };
        for (std::size_t i = 0; i < inc.size(); i++) {
            const int t = inc[i];
            for (int s = 0; s < 3; s++) {
                const int e = mv.triEdges[t][s];
                if (e < 0 || static_cast<std::size_t>(e) >= mv.nEdge) continue;
                const std::pair<int, int>& pr = mv.compEdges[e];
                if (pr.first != static_cast<int>(v) && pr.second != static_cast<int>(v)) continue;
                const std::array<int, 2>& tt = et[static_cast<std::size_t>(e)];
                const int u = (tt[0] == t) ? tt[1] : tt[0];
                if (u < 0 || !inDom[static_cast<std::size_t>(u)]) continue;
                const int j = indexOf(u);
                if (j < 0) continue;
                const int ra = find(static_cast<int>(i)), rb = find(j);
                if (ra != rb) parent[static_cast<std::size_t>(ra)] = rb;
            }
        }
        std::size_t fans = 0;
        for (std::size_t i = 0; i < inc.size(); i++)
            if (find(static_cast<int>(i)) == static_cast<int>(i)) ++fans;
        if (fans >= 2) ++nPinch;
    }
    uc.pinchVertices = nPinch;
    return uc;
}

}  // namespace refit
}  // namespace stl2step
