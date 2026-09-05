#!/usr/bin/env python3
"""RegionSet I1-I9 invariant checker (P0).

stdlib only — jsonschema is not imported. Structural validation is a hand-rolled
subset of tests/gates/regionset.schema.json (required keys, types, enums).
I7 and I7b are two separate functions keyed on closed360; never one function
with a branch. P1 consumes this checker and cannot argue with it.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# Enums / constants matching the frozen schema + refit.hpp LoopRole order
# ---------------------------------------------------------------------------

SURF_TYPES = ("plane", "cylinder", "cone", "sphere", "torus")
ORIGINS = ("planeGrow", "cylGrow", "filletStrip", "ngonWall", "chamferCone")
LOOP_ROLES = ("outer", "inner", "capLow", "capHigh")
# C++ enum class LoopRole : uint8_t { Outer, Inner, CapLow, CapHigh };
ROLE_INT = {"outer": 0, "inner": 1, "capLow": 2, "capHigh": 3}
REJECTS = (
    "none",
    "gaussPlanarity",
    "vertexResidual",
    "chordConsistency",
    "radiusSanity",
    "span",
    "filletConsensus",
    "neighborNotAnalytic",
    "stripWidth",
    "torusNYI",
    "coneNYI",
    "sphereNYI",
    "dirtyComponent",
    "faceBuildFailed",
    "chainUnstable",
)
BUILT_AS = ("notBuilt", "single", "seamed360", "twoHalves", "explodedToFacets")

REGION_REQUIRED = (
    "id",
    "type",
    "origin",
    "ax",
    "radius",
    "uMin",
    "uMax",
    "vMin",
    "vMax",
    "closed360",
    "outwardNormal",
    "tris",
    "loops",
    "maxVertexDev",
    "rmsVertexDev",
    "chordSagitta",
    "nSides",
    "dVolPredicted",
    "maxVertexSnap",
    "reject",
    "builtAs",
    "filletNbrA",
    "filletNbrB",
)
CHAIN_REQUIRED = (
    "regA",
    "regB",
    "islandA",
    "islandB",
    "tangent",
    "closedLoop",
    "meshEdges",
    "meshVerts",
)
STATS_REQUIRED = (
    "planes",
    "cylinders",
    "fillets",
    "rejected",
    "facetIslands",
    "facetTriangles",
    "distinctRadii",
    "maxVertexDev",
    "maxEdgeTol",
    "dVolPredicted",
)
DUMP_REQUIRED = (
    "compRoot",
    "regions",
    "rejected",
    "chains",
    "triRegion",
    "triIsland",
    "nIslands",
    "stats",
)

# Default rule set. I7 and I7b are distinct ids. G5 is first-class so a silent
# decline wrongly filed as rejected[] has its own id. SIDECAR is added when
# --sidecar is present.
DEFAULT_RULES = (
    "I1",
    "I2",
    "I3",
    "I4",
    "I5",
    "I6",
    "I7",
    "I7b",
    "I8",
    "I9",
    "G5",
)
ALL_RULES = DEFAULT_RULES + ("SIDECAR", "SIDECAR_RECOVERABLE")

EPS = 1e-12


# ---------------------------------------------------------------------------
# Tiny result helper
# ---------------------------------------------------------------------------

class Findings:
    """offender-id list plus messages; empty => pass."""

    def __init__(self):
        self.offenders = []  # stable unique ids
        self.notes = []      # human lines, parallel-ish

    def add(self, oid, msg):
        oid = str(oid)
        if oid not in self.offenders:
            self.offenders.append(oid)
        self.notes.append("%s: %s" % (oid, msg))

    def ok(self):
        return not self.offenders


# ---------------------------------------------------------------------------
# Dump accessors
# ---------------------------------------------------------------------------

def _as_list(x):
    return x if isinstance(x, list) else []


def _as_int(x, default=-1):
    if isinstance(x, bool) or x is None:
        return default
    if isinstance(x, int):
        return x
    if isinstance(x, float) and x == int(x):
        return int(x)
    return default


def _as_float(x, default=0.0):
    if isinstance(x, bool) or x is None:
        return default
    if isinstance(x, (int, float)):
        return float(x)
    return default


def _as_bool(x, default=False):
    if isinstance(x, bool):
        return x
    return default


def _min_or_none(xs):
    return min(xs) if xs else None


def regions_by_id(dump):
    out = {}
    for r in _as_list(dump.get("regions")):
        if isinstance(r, dict) and "id" in r:
            out[_as_int(r.get("id"))] = r
    return out


def chains_of(dump):
    return [c for c in _as_list(dump.get("chains")) if isinstance(c, dict)]


def accepted_of(dump):
    return [r for r in _as_list(dump.get("regions")) if isinstance(r, dict)]


def rejected_of(dump):
    return [r for r in _as_list(dump.get("rejected")) if isinstance(r, dict)]


def chain_sides(chain):
    """Return the two (kind, id) sides of a chain.

    kind is 'r' (region), 'i' (island), or 'ext' (both ids -1: open/exterior).
    """
    def one(reg, island):
        reg = _as_int(reg, -1)
        island = _as_int(island, -1)
        if reg >= 0:
            return ("r", reg)
        if island >= 0:
            return ("i", island)
        return ("ext", -1)

    return (
        one(chain.get("regA"), chain.get("islandA")),
        one(chain.get("regB"), chain.get("islandB")),
    )


def chain_touches_region(chain, rid):
    return _as_int(chain.get("regA"), -1) == rid or _as_int(chain.get("regB"), -1) == rid


def loop_min_chain(loop):
    idxs = [_as_int(x, None) for x in _as_list(loop.get("chainIdx"))]
    idxs = [x for x in idxs if x is not None]
    return _min_or_none(idxs)


def reversed_flag(val):
    if isinstance(val, bool):
        return 1 if val else 0
    return 1 if _as_int(val, 0) else 0


def chain_walk_ends(chain, rev):
    verts = _as_list(chain.get("meshVerts"))
    if not verts:
        return None, None
    if rev:
        return verts[-1], verts[0]
    return verts[0], verts[-1]


def touching_chain_indices(dump, rid):
    out = []
    for i, c in enumerate(chains_of(dump)):
        if chain_touches_region(c, rid):
            out.append(i)
    return out


def loop_chain_multiset(region):
    """chain index -> count of appearances across this region's loops."""
    counts = defaultdict(int)
    for loop in _as_list(region.get("loops")):
        if not isinstance(loop, dict):
            continue
        for ci in _as_list(loop.get("chainIdx")):
            counts[_as_int(ci)] += 1
    return counts


def coverage_offenders(dump, region, findings, prefix):
    """Every chain touching the region appears in exactly one of its loops."""
    rid = _as_int(region.get("id"))
    n_chains = len(chains_of(dump))
    touching = set(touching_chain_indices(dump, rid))
    counts = loop_chain_multiset(region)
    for ci, n in sorted(counts.items()):
        if ci < 0 or ci >= n_chains:
            findings.add("%s.loopChain%d" % (prefix, ci),
                         "loop chainIdx %d out of range" % ci)
            continue
        if ci not in touching:
            findings.add("%s.loopChain%d" % (prefix, ci),
                         "loop cites chain %d which does not touch region %d" % (ci, rid))
        elif n != 1:
            findings.add("%s.chain%d" % (prefix, ci),
                         "chain %d appears in %d loops of region %d (need exactly 1)"
                         % (ci, n, rid))
    for ci in sorted(touching):
        if counts.get(ci, 0) == 0:
            findings.add("%s.chain%d" % (prefix, ci),
                         "chain %d touches region %d but is in none of its loops"
                         % (ci, rid))
    for li, loop in enumerate(_as_list(region.get("loops"))):
        if not isinstance(loop, dict):
            findings.add("%s.loop%d" % (prefix, li), "loop is not an object")
            continue
        a = _as_list(loop.get("chainIdx"))
        b = _as_list(loop.get("reversed"))
        if len(a) != len(b):
            findings.add("%s.loop%d" % (prefix, li),
                         "chainIdx len %d != reversed len %d" % (len(a), len(b)))


# ---------------------------------------------------------------------------
# SCHEMA — hand-rolled subset (not jsonschema)
# ---------------------------------------------------------------------------

def check_SCHEMA(dump, ctx):
    f = Findings()
    if not isinstance(dump, dict):
        f.add("dump", "top level is not an object")
        return f
    for k in DUMP_REQUIRED:
        if k not in dump:
            f.add("dump", "missing required key %s" % k)
    for k in ("regions", "rejected", "chains", "triRegion", "triIsland"):
        if k in dump and not isinstance(dump[k], list):
            f.add(k, "%s is not an array" % k)
    if "nIslands" in dump and not isinstance(dump["nIslands"], int):
        f.add("nIslands", "nIslands is not an integer")
    if "compRoot" in dump and not isinstance(dump["compRoot"], int):
        f.add("compRoot", "compRoot is not an integer")

    enums = ctx.get("enums") or {}
    surf = enums.get("type") or SURF_TYPES
    orig = enums.get("origin") or ORIGINS
    roles = enums.get("role") or LOOP_ROLES
    rejs = enums.get("reject") or REJECTS
    built = enums.get("builtAs") or BUILT_AS

    def check_region(r, label):
        if not isinstance(r, dict):
            f.add(label, "not an object")
            return
        for k in REGION_REQUIRED:
            if k not in r:
                f.add(label, "missing %s" % k)
        t = r.get("type")
        if t is not None and t not in surf:
            f.add(label, "type %r not in enum" % t)
        o = r.get("origin")
        if o is not None and o not in orig:
            f.add(label, "origin %r not in enum" % o)
        rj = r.get("reject")
        if rj is not None and rj not in rejs:
            f.add(label, "reject %r not in enum" % rj)
        ba = r.get("builtAs")
        if ba is not None and ba not in built:
            f.add(label, "builtAs %r not in enum" % ba)
        if "closed360" in r and not isinstance(r["closed360"], bool):
            f.add(label, "closed360 is not a boolean")
        if "outwardNormal" in r and not isinstance(r["outwardNormal"], bool):
            f.add(label, "outwardNormal is not a boolean")
        if "tris" in r and not isinstance(r["tris"], list):
            f.add(label, "tris is not an array")
        if "loops" in r and not isinstance(r["loops"], list):
            f.add(label, "loops is not an array")
        ax = r.get("ax")
        if ax is not None:
            if not isinstance(ax, dict):
                f.add(label, "ax is not an object")
            else:
                for kk in ("loc", "dir", "xdir"):
                    v = ax.get(kk)
                    if not (isinstance(v, list) and len(v) == 3):
                        f.add(label, "ax.%s is not a length-3 array" % kk)
        for li, loop in enumerate(_as_list(r.get("loops"))):
            lp = "%s.loop%d" % (label, li)
            if not isinstance(loop, dict):
                f.add(lp, "not an object")
                continue
            for k in ("chainIdx", "reversed", "role"):
                if k not in loop:
                    f.add(lp, "missing %s" % k)
            role = loop.get("role")
            if role is not None and role not in roles:
                f.add(lp, "role %r not in enum" % role)

    for i, r in enumerate(_as_list(dump.get("regions"))):
        check_region(r, "region%d" % i)
    for i, r in enumerate(_as_list(dump.get("rejected"))):
        check_region(r, "rejected%d" % i)
    for i, c in enumerate(_as_list(dump.get("chains"))):
        label = "chain%d" % i
        if not isinstance(c, dict):
            f.add(label, "not an object")
            continue
        for k in CHAIN_REQUIRED:
            if k not in c:
                f.add(label, "missing %s" % k)
        if "meshEdges" in c and not isinstance(c["meshEdges"], list):
            f.add(label, "meshEdges is not an array")
        if "meshVerts" in c and not isinstance(c["meshVerts"], list):
            f.add(label, "meshVerts is not an array")
        if "closedLoop" in c and not isinstance(c["closedLoop"], bool):
            f.add(label, "closedLoop is not a boolean")
    stats = dump.get("stats")
    if stats is not None:
        if not isinstance(stats, dict):
            f.add("stats", "stats is not an object")
        else:
            for k in STATS_REQUIRED:
                if k not in stats:
                    f.add("stats", "missing %s" % k)
    tri_r = dump.get("triRegion")
    tri_i = dump.get("triIsland")
    if isinstance(tri_r, list) and isinstance(tri_i, list):
        for i, v in enumerate(tri_r):
            if not isinstance(v, int):
                f.add("triRegion[%d]" % i, "not an integer")
        for i, v in enumerate(tri_i):
            if not isinstance(v, int):
                f.add("triIsland[%d]" % i, "not an integer")
    return f


# ---------------------------------------------------------------------------
# I1 — total XOR partition
# ---------------------------------------------------------------------------

def check_I1(dump, ctx):
    f = Findings()
    tri_r = _as_list(dump.get("triRegion"))
    tri_i = _as_list(dump.get("triIsland"))
    if len(tri_r) != len(tri_i):
        f.add("triRegion/triIsland",
              "length mismatch %d vs %d" % (len(tri_r), len(tri_i)))
    n = min(len(tri_r), len(tri_i))
    by_id = regions_by_id(dump)
    n_islands = _as_int(dump.get("nIslands"), 0)
    claimed = defaultdict(list)
    island_claimed = defaultdict(list)
    for t in range(n):
        r = _as_int(tri_r[t], -1)
        i = _as_int(tri_i[t], -1)
        in_r = r >= 0
        in_i = i >= 0
        if in_r and in_i:
            f.add("t%d" % t, "triangle %d in both triRegion=%d and triIsland=%d" % (t, r, i))
        elif (not in_r) and (not in_i):
            f.add("t%d" % t, "triangle %d in neither triRegion nor triIsland" % t)
        if in_r:
            if r not in by_id:
                f.add("t%d" % t, "triRegion[%d]=%d is not an accepted region id" % (t, r))
            claimed[r].append(t)
        if in_i:
            if i < 0 or i >= n_islands:
                f.add("t%d" % t, "triIsland[%d]=%d outside [0, nIslands=%d)" % (t, i, n_islands))
            island_claimed[i].append(t)
    for rid, r in by_id.items():
        listed = [_as_int(t) for t in _as_list(r.get("tris"))]
        listed_set = set(listed)
        mapped = set(claimed.get(rid, []))
        if listed_set != mapped:
            f.add("region%d" % rid,
                  "region %d tris %s != triRegion membership %s"
                  % (rid, sorted(listed_set), sorted(mapped)))
    if n_islands < 0:
        f.add("nIslands", "nIslands is negative")
    used_islands = set(island_claimed.keys())
    expected = set(range(n_islands)) if n_islands > 0 else set()
    if used_islands != expected:
        # holes or extras in the island id space
        extra = used_islands - expected
        missing = expected - used_islands
        if extra:
            f.add("nIslands", "triIsland uses ids %s but nIslands=%d" % (sorted(extra), n_islands))
        if missing:
            f.add("nIslands", "nIslands=%d but island ids %s never appear" % (n_islands, sorted(missing)))
    return f


# ---------------------------------------------------------------------------
# I2 — each mesh edge in at most one chain
# ---------------------------------------------------------------------------

def check_I2(dump, ctx):
    """Dump-level I2: a mesh edge appears in at most one chain.

    The 'exactly one if the two triangles are in different regions/islands'
    and 'none if interior' directions need mesh adjacency (triEdges), which
    the frozen dump schema does not carry. See the lane report.
    """
    f = Findings()
    owners = defaultdict(list)
    for ci, c in enumerate(chains_of(dump)):
        seen_here = set()
        for e in _as_list(c.get("meshEdges")):
            e = _as_int(e)
            if e in seen_here:
                f.add("e%d" % e, "edge %d repeated inside chain %d" % (e, ci))
            seen_here.add(e)
            owners[e].append(ci)
    for e, cs in sorted(owners.items()):
        if len(cs) > 1:
            f.add("e%d" % e, "edge %d appears in chains %s (need <= 1)" % (e, cs))
    return f


# ---------------------------------------------------------------------------
# I3 — chain length relation, endpoint walk, -1 island convention
# ---------------------------------------------------------------------------

def check_I3(dump, ctx):
    """Dump-level I3 (geometry-free).

    CAN verify from JSON:
      - closedLoop: len(meshVerts)==len(meshEdges); open: == len(meshEdges)+1
      - consecutive chains in a loop share the connecting vertex (reversed-aware)
      - per side: (reg>=0) and (island>=0) is forbidden (the -1 island convention);
        both -1 is allowed (open/exterior)

    CANNOT verify from JSON:
      - walking meshVerts keeps regA on the left (needs vertex coordinates and
        triangle winding)
    """
    f = Findings()
    chains = chains_of(dump)
    for ci, c in enumerate(chains):
        edges = _as_list(c.get("meshEdges"))
        verts = _as_list(c.get("meshVerts"))
        closed = _as_bool(c.get("closedLoop"), False)
        expect = len(edges) if closed else len(edges) + 1
        if len(verts) != expect:
            f.add("chain%d" % ci,
                  "meshVerts len %d != %s (closedLoop=%s, meshEdges=%d)"
                  % (len(verts), expect, closed, len(edges)))
        if closed and len(edges) == 0:
            f.add("chain%d" % ci, "closed chain with empty meshEdges")
        # -1 island convention per side
        for side, (reg_k, isl_k) in (("A", ("regA", "islandA")), ("B", ("regB", "islandB"))):
            reg = _as_int(c.get(reg_k), -1)
            isl = _as_int(c.get(isl_k), -1)
            if reg >= 0 and isl >= 0:
                f.add("chain%d.%s" % (ci, side),
                      "reg%s=%d and island%s=%d both set (need XOR / both -1 for exterior)"
                      % (side, reg, side, isl))
            if reg < -1:
                f.add("chain%d.%s" % (ci, side), "reg%s=%d is < -1" % (side, reg))
            if isl < -1:
                f.add("chain%d.%s" % (ci, side), "island%s=%d is < -1" % (side, isl))

    # Endpoint chaining inside each accepted region's loops.
    for r in accepted_of(dump):
        rid = _as_int(r.get("id"))
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                continue
            idxs = [_as_int(x) for x in _as_list(loop.get("chainIdx"))]
            revs = [reversed_flag(x) for x in _as_list(loop.get("reversed"))]
            if len(idxs) != len(revs) or not idxs:
                continue
            n = len(idxs)
            # A loop is a closed walk of chains. One closed chain is a ring.
            if n == 1:
                # Length/convention already checked per chain. A single-chain
                # loop has no inter-chain join to verify. Closure of that
                # chain is I7/I7b ("loops are complete and closed").
                continue
            for k in range(n):
                ci = idxs[k]
                cj = idxs[(k + 1) % n]
                if not (0 <= ci < len(chains) and 0 <= cj < len(chains)):
                    continue
                _s0, e0 = chain_walk_ends(chains[ci], revs[k])
                s1, _e1 = chain_walk_ends(chains[cj], revs[(k + 1) % n])
                if e0 is None or s1 is None:
                    continue
                if e0 != s1:
                    f.add("region%d.loop%d" % (rid, li),
                          "chain %d end vert %s does not meet chain %d start vert %s"
                          % (ci, e0, cj, s1))
    return f


# ---------------------------------------------------------------------------
# I4 — rms <= max, max <= tolerance (sidecar / stats)
# ---------------------------------------------------------------------------

def check_I4(dump, ctx):
    """Dump-level I4.

    CAN verify: rmsVertexDev <= maxVertexDev, both finite and >= 0, and
    maxVertexDev <= stats.maxVertexDev. If the sidecar carries epsPlane /
    smoothTolMM / tolerance, also maxVertexDev <= that eps.

    CANNOT verify from JSON: per-vertex dev <= maxVertexDev (no vertex list
    against the fitted surface in the dump).
    """
    f = Findings()
    sidecar = ctx.get("sidecar") or {}
    eps = None
    for key in ("epsPlane", "smoothTolMM", "tolerance", "eps"):
        if key in sidecar and isinstance(sidecar[key], (int, float)):
            eps = float(sidecar[key])
            break
    stats = dump.get("stats") if isinstance(dump.get("stats"), dict) else {}
    stats_max = stats.get("maxVertexDev")
    stats_max = _as_float(stats_max, None) if stats_max is not None else None

    for r in accepted_of(dump):
        rid = _as_int(r.get("id"))
        mx = _as_float(r.get("maxVertexDev"), 0.0)
        rms = _as_float(r.get("rmsVertexDev"), 0.0)
        if mx < -EPS:
            f.add("region%d" % rid, "maxVertexDev %g is negative" % mx)
        if rms < -EPS:
            f.add("region%d" % rid, "rmsVertexDev %g is negative" % rms)
        if rms > mx + EPS:
            f.add("region%d" % rid, "rmsVertexDev %g > maxVertexDev %g" % (rms, mx))
        if stats_max is not None and mx > stats_max + EPS:
            f.add("region%d" % rid,
                  "maxVertexDev %g > stats.maxVertexDev %g" % (mx, stats_max))
        if eps is not None and mx > eps + EPS:
            f.add("region%d" % rid,
                  "maxVertexDev %g > sidecar tolerance %g" % (mx, eps))
    return f


# ---------------------------------------------------------------------------
# I5 — §5.6 sort keys, no ties
# ---------------------------------------------------------------------------

def check_I5(dump, ctx):
    """Determinism checked structurally via DECISION-p1-math §5.6 sort keys.

    Fully checkable from JSON:
      - chains[] by minLocalMeshEdgeId, strictly increasing (a tie is a violation)
      - Region::loops by ((int)role, minChainIdx), strictly increasing
      - Region::tris strictly ascending local tri id (a duplicate is a tie)
      - accepted regions[i].id == i (dense ids assigned after the sort)
      - island ids used are exactly {0 .. nIslands-1}
      - unique minLocalTriId across regions; unique dense id across rejected
        (D5.6-A1 / I5-A1); every rejected entry non-empty

    NOT checkable from JSON (frozen schema has no `area` field):
      - regions[] / rejected[] / islands ordered by (-area, minLocalTriId, id).
        If a region object carries an extra `area` property we honour the full
        key; otherwise we do not invent triangle-count as area.
    """
    f = Findings()

    def strict_sorted(keys, label, oid):
        for i in range(len(keys) - 1):
            a, b = keys[i], keys[i + 1]
            if a is None or b is None:
                continue
            if a == b:
                f.add(oid, "%s tie at positions %d/%d key=%s" % (label, i, i + 1, a))
            elif a > b:
                f.add(oid, "%s mis-ordered at positions %d/%d %s > %s" % (label, i, i + 1, a, b))

    # regions[] dense ids + unique minLocalTriId + optional area order
    mins = []
    area_keys = []
    have_area = True
    for i, r in enumerate(accepted_of(dump)):
        rid = _as_int(r.get("id"))
        if rid != i:
            f.add("regions", "regions[%d].id=%d (need dense id==index after sort)" % (i, rid))
        tris = [_as_int(t) for t in _as_list(r.get("tris"))]
        strict_sorted(tris, "region%d.tris" % rid, "region%d.tris" % rid)
        mn = _min_or_none(tris)
        mins.append((mn, rid))
        if "area" in r and isinstance(r["area"], (int, float)):
            area_keys.append(((-float(r["area"]), mn), rid))
        else:
            have_area = False
        # loops
        loop_keys = []
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                continue
            role = loop.get("role")
            ri = ROLE_INT.get(role)
            mc = loop_min_chain(loop)
            if ri is None:
                f.add("region%d.loop%d" % (rid, li), "unknown role %r" % role)
                continue
            loop_keys.append((ri, mc if mc is not None else -1))
        strict_sorted(loop_keys, "region%d.loops" % rid, "region%d.loops" % rid)

    seen_min = {}
    for mn, rid in mins:
        if mn is None:
            continue
        if mn in seen_min:
            f.add("regions",
                  "minLocalTriId %d ties region %d and region %d" % (mn, seen_min[mn], rid))
        else:
            seen_min[mn] = rid
    if have_area and area_keys:
        strict_sorted([k for k, _ in area_keys], "regions[] (-area, minLocalTriId)", "regions")

    # rejected[] — D1.3-A6b / D5.6-A1 / I5-A1: non-empty tris, dense unique id,
    # minLocalTriId ties are legal; optional (-area, minLocalTriId, id) order
    r_area_keys = []
    r_have_area = True
    r_ids = []
    for i, r in enumerate(rejected_of(dump)):
        rid = _as_int(r.get("id"))
        r_ids.append(rid)
        tris = [_as_int(t) for t in _as_list(r.get("tris"))]
        if not tris:
            f.add("rejected%d" % i,
                  "rejected[] entry carries no triangles — D1.3-A6b requires "
                  "the claim's full pre-peel tri set")
        strict_sorted(tris, "rejected%d.tris" % i, "rejected%d.tris" % i)
        mn = _min_or_none(tris)
        if "area" in r and isinstance(r["area"], (int, float)):
            r_area_keys.append(((-float(r["area"]), mn, rid), i))
        else:
            r_have_area = False
        loop_keys = []
        for loop in _as_list(r.get("loops")):
            if not isinstance(loop, dict):
                continue
            role = loop.get("role")
            ri = ROLE_INT.get(role)
            mc = loop_min_chain(loop)
            if ri is None:
                continue
            loop_keys.append((ri, mc if mc is not None else -1))
        strict_sorted(loop_keys, "rejected[%d].loops" % i, "rejected%d.loops" % i)
    if r_have_area and r_area_keys:
        strict_sorted([k for k, _ in r_area_keys],
                      "rejected[] (-area, minLocalTriId, id)", "rejected")
    if sorted(r_ids) != list(range(len(r_ids))):
        f.add("rejected",
              "rejected[].id must be dense 0..n-1 (have %s)" % r_ids)

    # S12-b: vertexResidual must cover all hex-side provisionals (growx §5.7)
    sidecar = ctx.get("sidecar") or {}
    sid = _sidecar_id(sidecar, ctx.get("sidecar_path") or "").lower().replace("_", "-")
    if "s12-b" in sid or "s12b" in sid:
        side_tri_count = None
        for key in ("sideTriCount", "sideTri", "hexSideTriCount", "hexSideTris"):
            v = sidecar.get(key)
            if isinstance(v, int) and not isinstance(v, bool):
                side_tri_count = v
                break
        if side_tri_count is not None:
            vr = [r for r in rejected_of(dump) if r.get("reject") == "vertexResidual"]
            if len(vr) != 1:
                f.add("rejected",
                      "S12-b tri-coverage needs exactly one vertexResidual (have %d)"
                      % len(vr))
            else:
                n_tris = len(_as_list(vr[0].get("tris")))
                if n_tris != side_tri_count:
                    f.add("rejected",
                          "S12-b vertexResidual tris len %d != sidecar side-tri count %d"
                          % (n_tris, side_tri_count))

    # chains[] by minLocalMeshEdgeId
    chain_keys = []
    for ci, c in enumerate(chains_of(dump)):
        edges = [_as_int(e) for e in _as_list(c.get("meshEdges"))]
        chain_keys.append(_min_or_none(edges))
        # within-chain start vertex is the lowest local vertex id among
        # terminals (open) or among all verts (closed) — DECISION §5.6.
        verts = _as_list(c.get("meshVerts"))
        if verts:
            closed = _as_bool(c.get("closedLoop"), False)
            if closed:
                lowest = min(verts)
                if verts[0] != lowest:
                    f.add("chain%d" % ci,
                          "closed chain start vertex %s is not the lowest local vertex id %s"
                          % (verts[0], lowest))
            else:
                if len(verts) >= 1:
                    terminals = (verts[0], verts[-1])
                    lowest = min(terminals)
                    if verts[0] != lowest:
                        f.add("chain%d" % ci,
                              "open chain start vertex %s is not the lowest terminal id %s"
                              % (verts[0], lowest))
    strict_sorted(chain_keys, "chains[] minLocalMeshEdgeId", "chains")
    return f


# ---------------------------------------------------------------------------
# I6 — never called on a dirty component
# ---------------------------------------------------------------------------

def check_I6(dump, ctx):
    """Dump-level I6.

    The frozen schema has no `clean`/`dirty` field. The flag the dump *does*
    carry is Region.reject: accepted regions must be reject=='none', and
    'dirtyComponent' must not appear in regions[] or rejected[] (P1 is never
    called on a dirty component, DECISION-insertion §4.1 / I6). Optional
    non-schema keys `clean` / `dirty` / `dirtyComponent` are honoured if present.
    """
    f = Findings()
    if dump.get("clean") is False:
        f.add("clean", "dump.clean is false")
    if dump.get("dirty") is True:
        f.add("dirty", "dump.dirty is true")
    if dump.get("dirtyComponent") is True:
        f.add("dirtyComponent", "dump.dirtyComponent is true")
    for r in accepted_of(dump):
        rid = _as_int(r.get("id"))
        rej = r.get("reject")
        if rej != "none":
            f.add("region%d" % rid,
                  "accepted region %d has reject=%r (need none)" % (rid, rej))
        if rej == "dirtyComponent":
            f.add("region%d" % rid, "accepted region carries dirtyComponent")
    for i, r in enumerate(rejected_of(dump)):
        if r.get("reject") == "dirtyComponent":
            f.add("rejected%d" % i,
                  "rejected[] files dirtyComponent — P1 is never called on a dirty component")
    return f


# ---------------------------------------------------------------------------
# I7 — closed360 == false  (SEPARATE function, SEPARATE rule id)
# ---------------------------------------------------------------------------

def check_I7(dump, ctx):
    """I7 (amended): closed360 == false regions only.

    Exactly one role=='outer', plus zero or more 'inner'. No capLow/capHigh.
    Every chain touching the region appears in exactly one of its loops.
    """
    f = Findings()
    for r in accepted_of(dump):
        if _as_bool(r.get("closed360"), False) is True:
            continue  # keyed: I7 does not judge closed360 regions
        rid = _as_int(r.get("id"))
        prefix = "region%d" % rid
        roles = []
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                f.add("%s.loop%d" % (prefix, li), "loop is not an object")
                continue
            roles.append(loop.get("role"))
        n_outer = sum(1 for x in roles if x == "outer")
        n_cap_lo = sum(1 for x in roles if x == "capLow")
        n_cap_hi = sum(1 for x in roles if x == "capHigh")
        n_inner = sum(1 for x in roles if x == "inner")
        n_unknown = sum(1 for x in roles if x not in LOOP_ROLES)
        if n_outer != 1:
            f.add(prefix, "closed360=false needs exactly one outer, have %d" % n_outer)
        if n_cap_lo != 0 or n_cap_hi != 0:
            f.add(prefix,
                  "closed360=false must not carry capLow/capHigh (have %d/%d)"
                  % (n_cap_lo, n_cap_hi))
        if n_unknown:
            f.add(prefix, "closed360=false has %d loop(s) with unknown role" % n_unknown)
        _ = n_inner
        chains = chains_of(dump)
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                continue
            idxs = [_as_int(x) for x in _as_list(loop.get("chainIdx"))]
            if len(idxs) == 1 and 0 <= idxs[0] < len(chains):
                if not _as_bool(chains[idxs[0]].get("closedLoop"), False):
                    f.add("%s.loop%d" % (prefix, li),
                          "I7 loop is a single open chain %d; loops are complete and closed"
                          % idxs[0])
        coverage_offenders(dump, r, f, prefix)
    return f


# ---------------------------------------------------------------------------
# I7b — closed360 == true  (SEPARATE function, SEPARATE rule id)
# ---------------------------------------------------------------------------

def check_I7b(dump, ctx):
    """I7b (new): closed360 == true regions only.

    Exactly one capLow, exactly one capHigh, zero or more inner, NO outer.
    Every chain touching the region appears in exactly one loop.
    No loop corresponds to the seam (P1 emits no seam chain).
    """
    f = Findings()
    chains = chains_of(dump)
    for r in accepted_of(dump):
        if _as_bool(r.get("closed360"), False) is not True:
            continue  # keyed: I7b does not judge open regions
        rid = _as_int(r.get("id"))
        prefix = "region%d" % rid
        roles = []
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                f.add("%s.loop%d" % (prefix, li), "loop is not an object")
                continue
            role = loop.get("role")
            roles.append(role)
            if role == "outer":
                f.add("%s.loop%d" % (prefix, li),
                      "closed360=true must not carry role=outer")
        n_outer = sum(1 for x in roles if x == "outer")
        n_cap_lo = sum(1 for x in roles if x == "capLow")
        n_cap_hi = sum(1 for x in roles if x == "capHigh")
        n_unknown = sum(1 for x in roles if x not in LOOP_ROLES)
        if n_cap_lo != 1:
            f.add(prefix, "closed360=true needs exactly one capLow, have %d" % n_cap_lo)
        if n_cap_hi != 1:
            f.add(prefix, "closed360=true needs exactly one capHigh, have %d" % n_cap_hi)
        if n_outer != 0:
            f.add(prefix, "closed360=true needs zero outer, have %d" % n_outer)
        if n_unknown:
            f.add(prefix, "closed360=true has %d loop(s) with unknown role" % n_unknown)
        for li, loop in enumerate(_as_list(r.get("loops"))):
            if not isinstance(loop, dict):
                continue
            idxs = [_as_int(x) for x in _as_list(loop.get("chainIdx"))]
            if len(idxs) == 1 and 0 <= idxs[0] < len(chains):
                if not _as_bool(chains[idxs[0]].get("closedLoop"), False):
                    f.add("%s.loop%d" % (prefix, li),
                          "I7b loop is a single open chain %d; loops are complete and closed"
                          % idxs[0])

        # No seam chain: a chain with both sides equal to this region is the
        # u=0 generator, which P1 must not emit (DECISION-p1-math §3.3).
        for ci in touching_chain_indices(dump, rid):
            c = chains[ci]
            a = _as_int(c.get("regA"), -1)
            b = _as_int(c.get("regB"), -1)
            if a == rid and b == rid:
                f.add("%s.chain%d" % (prefix, ci),
                      "chain %d has regA=regB=%d — that is the seam; I7b forbids it"
                      % (ci, rid))

        coverage_offenders(dump, r, f, prefix)
    return f


# ---------------------------------------------------------------------------
# I8 — split at >=3-region vertices (collapse-safety linchpin)
# ---------------------------------------------------------------------------

def check_I8(dump, ctx):
    """Dump-level I8.

    CAN verify: an interior vertex of a chain (open: not a terminal; closed:
    every vertex) must not be incident to >=3 distinct regions/islands, and
    must not appear in any other chain. That is the unsplit 3-region vertex.

    CANNOT verify from JSON: splits at vertices that touch a non-manifold or
    open mesh edge, because the dump does not carry edge valence / openness.
    """
    f = Findings()
    chains = chains_of(dump)

    # vertex -> list of (chain_index, is_interior)
    at = defaultdict(list)
    sides_at = defaultdict(set)
    for ci, c in enumerate(chains):
        verts = _as_list(c.get("meshVerts"))
        closed = _as_bool(c.get("closedLoop"), False)
        n = len(verts)
        s0, s1 = chain_sides(c)
        for vi, v in enumerate(verts):
            if closed:
                interior = n > 0
            else:
                interior = (vi != 0 and vi != n - 1)
            at[v].append((ci, interior))
            sides_at[v].add(s0)
            sides_at[v].add(s1)

    for v, hits in at.items():
        n_chains = len({ci for ci, _ in hits})
        interior_in = [ci for ci, interior in hits if interior]
        kinds = sides_at[v]
        # exterior does not count as a region/island
        distinct = {s for s in kinds if s[0] != "ext"}
        if interior_in and n_chains > 1:
            f.add("v%s" % v,
                  "vertex %s is interior to chain(s) %s but also appears in %d chains — unsplit junction"
                  % (v, interior_in, n_chains))
        if interior_in and len(distinct) >= 3:
            f.add("v%s" % v,
                  "vertex %s is interior to chain(s) %s and incident to %d regions/islands %s"
                  % (v, interior_in, len(distinct), sorted(distinct)))
    return f


# ---------------------------------------------------------------------------
# I9 — signed dVolPredicted, sigma = outwardNormal ? +1 : -1
# ---------------------------------------------------------------------------

def check_I9(dump, ctx):
    """I9 / D4.4.

    sigma = outwardNormal ? +1 : -1. For cylinder / fillet (non-plane):
    sign(dVolPredicted) must equal sigma, and dVolPredicted==0 only when the
    D4 guard fires (nSides<3 or R<=0 or empty tris / zero area).

    Planes use D4.3's own formula and may be ~0 even with outwardNormal true
    (the frozen min example is a plane with dVolPredicted=0). Their sign is
    not required to match sigma; a strictly-positive dVol on outwardNormal
    false is still a violation.
    """
    f = Findings()
    for r in accepted_of(dump):
        rid = _as_int(r.get("id"))
        outward = _as_bool(r.get("outwardNormal"), True)
        sigma = 1.0 if outward else -1.0
        dvol = _as_float(r.get("dVolPredicted"), 0.0)
        typ = r.get("type")
        origin = r.get("origin")
        n_sides = _as_int(r.get("nSides"), 0)
        radius = _as_float(r.get("radius"), 0.0)
        tris = _as_list(r.get("tris"))
        zero_area = len(tris) == 0
        guard = (n_sides < 3 or radius <= EPS)
        is_plane = typ == "plane"

        if is_plane:
            # D4.3: flattening an outward bulge removes volume. A positive
            # dVolPredicted with outwardNormal false is still inconsistent.
            if (not outward) and dvol > EPS:
                f.add("region%d" % rid,
                      "plane dVolPredicted=%g > 0 with outwardNormal=false (sigma=-1)"
                      % dvol)
            continue

        # Cylinder / cone / sphere / torus / fillet strip
        if abs(dvol) <= EPS:
            if zero_area or guard:
                continue
            # non-zero area analytic region should have a signed term unless
            # the D4 guard applied. FilletStrip uses the same formula.
            if origin == "filletStrip" or typ in ("cylinder", "cone", "sphere", "torus"):
                f.add("region%d" % rid,
                      "dVolPredicted=0 on a non-zero-area %s (nSides=%d R=%g); zero only for zero area / D4 guard"
                      % (typ, n_sides, radius))
            continue
        sign = 1.0 if dvol > 0 else -1.0
        if sign != sigma:
            f.add("region%d" % rid,
                  "dVolPredicted=%g has sign %+g but outwardNormal=%s so sigma=%+g"
                  % (dvol, sign, outward, sigma))
    return f


# ---------------------------------------------------------------------------
# G5 — declines produce NO rejected[] entry (D5.5 / D1.3)
# ---------------------------------------------------------------------------

def check_G5(dump, ctx):
    """G5 span/extent declines are silent (DECISION-p1-math D1.3 / D5.5).

    Dump-level: rejected[].reject must never be 'span' (that enumerator is the
    G5 span gate, which must not be filed). Sidecar-level: a fixture that
    declares g5Declines must have no extra rejected[] entry for that decline.
    """
    f = Findings()
    for i, r in enumerate(rejected_of(dump)):
        if r.get("reject") == "span":
            f.add("rejected%d" % i,
                  "rejected[] files reject=span — G5 declines produce NO rejected[] entry")
    sidecar = ctx.get("sidecar") or {}
    g5 = sidecar.get("g5Declines", sidecar.get("g5Decline", sidecar.get("G5Declines")))
    if g5:
        # any rejected[] entry is wrong if the sidecar says this fixture is a
        # G5-decline (and lists no expected G1–G4 codes).
        expected = _sidecar_rejected_codes(sidecar)
        if not expected:
            for i, r in enumerate(rejected_of(dump)):
                f.add("rejected%d" % i,
                      "sidecar declares a G5 decline but rejected[] has an entry (reject=%r)"
                      % r.get("reject"))
        else:
            for i, r in enumerate(rejected_of(dump)):
                if r.get("reject") not in expected:
                    f.add("rejected%d" % i,
                          "G5 decline filed as rejected[] reject=%r (not in sidecar expected %s)"
                          % (r.get("reject"), expected))
    return f


# ---------------------------------------------------------------------------
# SIDECAR — fixture expected.json assertions
# ---------------------------------------------------------------------------

def _sidecar_id(sidecar, path):
    for k in ("id", "fixture", "name"):
        v = sidecar.get(k)
        if isinstance(v, str) and v:
            return v
    base = os.path.basename(path or "")
    return base


def _sidecar_rejected_codes(sidecar):
    raw = sidecar.get("rejected", sidecar.get("rejectedCodes", sidecar.get("reject")))
    if raw is None:
        return None
    if isinstance(raw, str):
        return [raw]
    if isinstance(raw, list):
        out = []
        for x in raw:
            if isinstance(x, str):
                out.append(x)
            elif isinstance(x, dict) and "reject" in x:
                out.append(x["reject"])
        return out
    return None


def _infer_sidecar_codes(sid):
    """Fixture-name heuristics from SPEC-P0 / DECISION-p1-math D6.3."""
    s = sid.lower().replace("_", "-")
    codes = None
    exactly = False
    forbid_origin = []
    if "s12-b" in s or "s12b" in s:
        codes = ["vertexResidual"]
        exactly = True
    elif "s11-b" in s or "s11b" in s:
        codes = []
    elif "s11" in s:
        codes = ["filletConsensus"]
        forbid_origin = ["filletStrip"]
    # D-130-15(2): the "cone" name heuristic asserted a premise D-130-3 retired
    # (chamfer cones are built since b50328c); boss_cone_chamfer's row was red
    # for having SUCCEEDED. Retired by that ruling. Sphere stays (SphereNYI).
    elif "sphere" in s:
        codes = ["sphereNYI"]
    return codes, exactly, forbid_origin


def check_SIDECAR(dump, ctx):
    f = Findings()
    sidecar = ctx.get("sidecar")
    if not sidecar:
        return f
    path = ctx.get("sidecar_path") or ""
    sid = _sidecar_id(sidecar, path)
    inferred, inferred_exactly, inferred_forbid = _infer_sidecar_codes(sid)

    codes = _sidecar_rejected_codes(sidecar)
    exactly = bool(sidecar.get("rejectedExactly", sidecar.get("exactly", False)))
    if codes is None:
        codes = inferred
        exactly = exactly or inferred_exactly
    if codes is None:
        codes = []

    have = [r.get("reject") for r in rejected_of(dump) if r.get("reject") and r.get("reject") != "none"]
    # expected codes present
    have_bag = list(have)
    for code in codes:
        if code not in have_bag:
            f.add("rejected", "sidecar expects reject=%r but it is absent (have %s)" % (code, have))
        else:
            have_bag.remove(code)
    if exactly:
        extra = list((Counter(have) - Counter(codes)).elements())
        if extra:
            f.add("rejected",
                  "sidecar wants exactly %s but dump has extra %s" % (codes, extra))
        if "s12-b" in sid.lower().replace("_", "-") or "s12b" in sid.lower():
            n_vr = sum(1 for c in have if c == "vertexResidual")
            if n_vr != 1:
                f.add("rejected", "S12-b must carry exactly one vertexResidual (have %d)" % n_vr)

    forbid = list(sidecar.get("forbidOrigin", sidecar.get("originForbidden", inferred_forbid or [])))
    if "s11" in sid.lower().replace("_", "-") and "s11-b" not in sid.lower().replace("_", "-") and "s11b" not in sid.lower():
        if "filletStrip" not in forbid:
            forbid.append("filletStrip")
        # S11 carries filletConsensus
        if "filletConsensus" not in codes:
            if "filletConsensus" not in have:
                f.add("rejected", "S11 sidecar: rejected[] must carry filletConsensus")
    for r in accepted_of(dump) + rejected_of(dump):
        if r.get("origin") in forbid:
            f.add("region%s" % r.get("id"),
                  "origin=%r is forbidden by sidecar (%s)" % (r.get("origin"), sid))

    # recoverable primitives -> matching accepted region (D-130-22(3))
    rec = sidecar.get("recoverable") or sidecar.get("mustRecover") or []
    if isinstance(rec, list) and rec and not ctx.get("skip_recoverable"):
        comp_idx = ctx.get("component_index")
        if comp_idx is None or _component_has_analytic_census(sidecar, comp_idx):
            reg_map = regions_by_id(dump)
            pool = [(r, reg_map) for r in accepted_of(dump)]
            _check_recoverable_entries(rec, pool, f)

    # mustRemainFaceted -> triIsland, not triRegion
    faceted = sidecar.get("mustRemainFaceted") or sidecar.get("mustFallback") or []
    tri_r = _as_list(dump.get("triRegion"))
    tri_i = _as_list(dump.get("triIsland"))
    if isinstance(faceted, list):
        for i, entry in enumerate(faceted):
            tris = []
            typ = None
            if isinstance(entry, dict):
                tris = _as_list(entry.get("tris"))
                typ = entry.get("type")
            if tris:
                for t in tris:
                    t = _as_int(t)
                    if t < 0 or t >= len(tri_r) or t >= len(tri_i):
                        f.add("mustRemainFaceted%d" % i, "triangle %d out of range" % t)
                        continue
                    if _as_int(tri_r[t], -1) >= 0:
                        f.add("t%d" % t,
                              "mustRemainFaceted triangle %d maps to triRegion=%d, need triIsland"
                              % (t, tri_r[t]))
                    if _as_int(tri_i[t], -1) < 0:
                        f.add("t%d" % t,
                              "mustRemainFaceted triangle %d is in neither (need triIsland)" % t)
            elif typ:
                for r in accepted_of(dump):
                    if r.get("type") == typ:
                        f.add("region%d" % _as_int(r.get("id")),
                              "mustRemainFaceted type %s is an accepted region" % typ)
    return f


# ---------------------------------------------------------------------------
# Schema enum extraction (optional --schema)
# ---------------------------------------------------------------------------

def _vec3(v):
    if not (isinstance(v, list) and len(v) == 3):
        return None
    try:
        return tuple(float(x) for x in v)
    except (TypeError, ValueError):
        return None


def _norm(v):
    if v is None:
        return None
    m = math.sqrt(sum(x * x for x in v))
    if m <= 0.0:
        return None
    return tuple(x / m for x in v)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _geo_tol(ref):
    return max(1e-6, 0.005 * abs(_as_float(ref, 0.0)))


def _dirs_agree(a, b, tol=0.01):
    na, nb = _norm(a), _norm(b)
    if na is None or nb is None:
        return False
    d = abs(_dot(na, nb))
    return d >= 1.0 - tol


def _plane_normal_from_region(r, regions_by_id):
    """Normal of the surface that ships for this region (D-130-22(3))."""
    if r.get("type") == "plane":
        ax = r.get("ax") or {}
        return _norm(_vec3(ax.get("dir")))
    if r.get("origin") != "filletStrip":
        return None
    na = nb = None
    for key in ("filletNbrA", "filletNbrB"):
        nid = _as_int(r.get(key), -1)
        if nid < 0:
            continue
        nbr = regions_by_id.get(nid)
        if not nbr or nbr.get("type") != "plane":
            continue
        n = _norm(_vec3((nbr.get("ax") or {}).get("dir")))
        if n is None:
            continue
        if key == "filletNbrA":
            na = n
        else:
            nb = n
    if na is None or nb is None:
        return None
    return _norm(tuple(na[i] + nb[i] for i in range(3)))


def _shipped_surface_type(r, regions_by_id):
    """Detector label vs shipped surface (D-130-22(3))."""
    origin = r.get("origin")
    if origin == "chamferCone" or r.get("type") == "cone":
        return "cone"
    if r.get("type") == "plane":
        return "plane"
    if origin == "filletStrip":
        return "cylinder"
    return r.get("type")


def _cone_radii(r):
    r0 = _as_float(r.get("radius"), 0.0)
    r1 = _as_float(r.get("vMax"), 0.0)
    r2 = _as_float(r.get("vMin"), 0.0)
    vals = [x for x in (r0, r1, r2) if x > 0.0]
    if len(vals) < 2:
        return (r0, r1) if r0 > 0.0 else (None, None)
    vals = sorted(set(vals))
    if len(vals) == 1:
        return vals[0], vals[0]
    return vals[0], vals[-1]


def _coaxial(loc, axis_loc, axis_dir, tol):
    if loc is None or axis_loc is None or axis_dir is None:
        return True
    ad = _norm(axis_dir)
    if ad is None:
        return False
    d = tuple(loc[i] - axis_loc[i] for i in range(3))
    along = _dot(d, ad)
    perp = math.sqrt(sum((d[i] - along * ad[i]) ** 2 for i in range(3)))
    return perp <= tol


def _axis_of(entry):
    ax = entry.get("axis") if isinstance(entry, dict) else None
    if not isinstance(ax, dict):
        return None, None
    return _vec3(ax.get("loc")), _vec3(ax.get("dir"))


def _cone_entry_radii(entry):
    if "radiusLo" in entry or "radiusHi" in entry:
        return _as_float(entry.get("radiusHi"), 0.0), _as_float(entry.get("radiusLo"), 0.0)
    e0 = _as_float(entry.get("radius"), 0.0)
    e1 = _as_float(entry.get("radius2"), e0)
    return e0, e1


def _live_row_for_component(sidecar, component_index):
    if component_index is None:
        return None
    want = _as_int(component_index)
    for row in sidecar.get("live") or []:
        if isinstance(row, dict) and _as_int(row.get("component")) == want:
            return row
    return None


def _component_has_analytic_census(sidecar, component_index):
    """True when live[] says this component ships plane/cylinder/cone analytics."""
    row = _live_row_for_component(sidecar, component_index)
    if row is None:
        return True
    census = row.get("surfaceCensus") or {}
    for key in ("plane", "cylinder", "cone"):
        if _as_int(census.get(key)) > 0:
            return True
    return False


def _recoverable_pool_from_paths(paths):
    pool = []
    for path in paths:
        with open(path, "r", encoding="utf-8") as fh:
            dump = json.load(fh)
        reg_map = regions_by_id(dump)
        for r in accepted_of(dump):
            pool.append((r, reg_map))
    return pool


def _check_recoverable_entries(rec, pool, findings=None):
    """Match recoverable[] against a pool of (region, reg_map) pairs."""
    f = findings or Findings()
    if not isinstance(rec, list):
        return f
    unused = list(pool)
    for i, entry in enumerate(rec):
        hit_idx = None
        for j, (r, reg_map) in enumerate(unused):
            if _match_recoverable_entry(entry, r, reg_map):
                hit_idx = j
                break
        if hit_idx is None:
            f.add("recoverable%d" % i,
                  "sidecar recoverable %s has no matching accepted region" % entry)
        else:
            unused.pop(hit_idx)
    return f


def check_SIDECAR_RECOVERABLE(dump, ctx):
    """Fixture-wide recoverable[] vs union of all components' shipped surfaces."""
    sidecar = ctx.get("sidecar")
    if not sidecar:
        return Findings()
    pool = ctx.get("recoverable_pool")
    if not pool:
        return Findings()
    rec = sidecar.get("recoverable") or sidecar.get("mustRecover") or []
    return _check_recoverable_entries(rec, pool)


def _match_recoverable_entry(entry, r, regions_by_id):
    """Match a sidecar recoverable entry to a region by shipped geometry."""
    if not isinstance(entry, dict):
        return _shipped_surface_type(r, regions_by_id) == entry
    want = entry.get("type")
    if "id" in entry and _as_int(r.get("id")) != _as_int(entry.get("id")):
        return False
    if want == "plane":
        want_n = _norm(_vec3(entry.get("normal")))
        if want_n is None:
            return False
        if r.get("type") == "plane":
            got_n = _plane_normal_from_region(r, regions_by_id)
        elif r.get("origin") == "filletStrip":
            got_n = _plane_normal_from_region(r, regions_by_id)
        else:
            return False
        if got_n is None:
            return False
        return _dirs_agree(want_n, got_n) or _dirs_agree(want_n, tuple(-x for x in got_n))
    shipped = _shipped_surface_type(r, regions_by_id)
    if want and shipped != want:
        return False
    if want == "cylinder":
        if "radius" in entry:
            tol = _geo_tol(entry.get("radius"))
            if abs(_as_float(r.get("radius")) - _as_float(entry.get("radius"))) > tol:
                return False
        if "closed360" in entry and _as_bool(r.get("closed360")) != bool(entry["closed360"]):
            return False
        eloc, edir = _axis_of(entry)
        rloc, rdir = _vec3((r.get("ax") or {}).get("loc")), _vec3((r.get("ax") or {}).get("dir"))
        if edir is not None and rdir is not None and not _dirs_agree(edir, rdir):
            return False
        if not _coaxial(eloc, rloc, rdir, max(_geo_tol(entry.get("radius", 1.0)), 0.05)):
            return False
        return True
    if want == "cone":
        r_lo, r_hi = _cone_radii(r)
        if r_lo is None:
            return False
        e0, e1 = _cone_entry_radii(entry)
        tol = max(_geo_tol(e0), _geo_tol(e1))
        got = sorted((r_lo, r_hi))
        want_r = sorted((e0, e1))
        if abs(got[0] - want_r[0]) > tol or abs(got[1] - want_r[1]) > tol:
            return False
        eloc, edir = _axis_of(entry)
        rloc, rdir = _vec3((r.get("ax") or {}).get("loc")), _vec3((r.get("ax") or {}).get("dir"))
        if edir is not None and rdir is not None and not _dirs_agree(edir, rdir):
            return False
        span = abs(want_r[1] - want_r[0])
        if not _coaxial(
            eloc, rloc, rdir,
            max(span + 0.05, _geo_tol(max(abs(x) for x in eloc)) if eloc else 0.05),
        ):
            return False
        if "closed360" in entry and _as_bool(r.get("closed360")) != bool(entry["closed360"]):
            return False
        return True
    if "radius" in entry:
        tol = _geo_tol(entry.get("radius"))
        if abs(_as_float(r.get("radius")) - _as_float(entry.get("radius"))) > tol:
            return False
    if "closed360" in entry and _as_bool(r.get("closed360")) != bool(entry["closed360"]):
        return False
    return True


def load_schema_enums(path):
    with open(path, "r", encoding="utf-8") as fh:
        schema = json.load(fh)
    defs = schema.get("$defs") or schema.get("definitions") or {}
    enums = {}
    region = defs.get("Region") or {}
    props = region.get("properties") or {}
    for key in ("type", "origin", "reject", "builtAs"):
        e = (props.get(key) or {}).get("enum")
        if e:
            enums[key] = tuple(e)
    loop = defs.get("Loop") or {}
    role = ((loop.get("properties") or {}).get("role") or {}).get("enum")
    if role:
        enums["role"] = tuple(role)
    # D-130-22(2): dump_regionset emits chamferCone / ngonWall even if the
    # frozen schema predates them.
    if "origin" in enums:
        enums["origin"] = tuple(dict.fromkeys((*enums["origin"], *ORIGINS)))
    return enums


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

RULE_FUNCS = {
    "I1": check_I1,
    "I2": check_I2,
    "I3": check_I3,
    "I4": check_I4,
    "I5": check_I5,
    "I6": check_I6,
    "I7": check_I7,
    "I7b": check_I7b,
    "I8": check_I8,
    "I9": check_I9,
    "G5": check_G5,
    "SIDECAR": check_SIDECAR,
    "SIDECAR_RECOVERABLE": check_SIDECAR_RECOVERABLE,
}


def parse_rules(s):
    if not s:
        return None
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if part not in RULE_FUNCS:
            raise SystemExit("unknown --rule id %r (known: %s)" % (part, ", ".join(ALL_RULES)))
        if part not in out:
            out.append(part)
    return out


def default_schema_path():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "regionset.schema.json")


def run_rules(dump, rules, ctx, verbose=False):
    report = {}
    # SCHEMA is a precondition, not a numbered I-rule. Always run it; surface
    # it in the report under SCHEMA so a structurally broken dump is not a
    # silent crash. --rule filtering does not skip SCHEMA unless SCHEMA is
    # the only thing the caller asked for... actually always run it.
    schema_f = check_SCHEMA(dump, ctx)
    report["SCHEMA"] = {
        "pass": schema_f.ok(),
        "offenders": list(schema_f.offenders),
        "notes": list(schema_f.notes),
    }
    if verbose and not schema_f.ok():
        sys.stderr.write("SCHEMA FAIL\n")
        for n in schema_f.notes:
            sys.stderr.write("  %s\n" % n)

    for rid in rules:
        fn = RULE_FUNCS[rid]
        findings = fn(dump, ctx)
        report[rid] = {
            "pass": findings.ok(),
            "offenders": list(findings.offenders),
            "notes": list(findings.notes),
        }
        if verbose:
            flag = "PASS" if findings.ok() else "FAIL"
            sys.stderr.write("%s %s\n" % (rid, flag))
            for n in findings.notes:
                sys.stderr.write("  %s\n" % n)
    return report


def report_ok(report, rules):
    if not report.get("SCHEMA", {}).get("pass", False):
        return False
    for rid in rules:
        if not report.get(rid, {}).get("pass", False):
            return False
    return True


def emit_json(dump_path, report, rules, ok):
    # Stable: sort_keys. Machine contract for run_gates.py:
    # rule id -> pass/fail -> offending ids.
    slim = {}
    # include SCHEMA plus selected rules, in a fixed order
    order = ["SCHEMA"] + [r for r in ALL_RULES if r in rules]
    seen = set()
    for rid in order:
        if rid in report and rid not in seen:
            slim[rid] = {
                "pass": bool(report[rid]["pass"]),
                "offenders": list(report[rid]["offenders"]),
            }
            seen.add(rid)
    payload = {
        "ok": bool(ok),
        "dump": dump_path,
        "rules": slim,
    }
    json.dump(payload, sys.stdout, sort_keys=True, indent=2)
    sys.stdout.write("\n")


def emit_human(report, rules, ok):
    order = ["SCHEMA"] + [r for r in ALL_RULES if r in rules]
    seen = set()
    for rid in order:
        if rid not in report or rid in seen:
            continue
        seen.add(rid)
        rec = report[rid]
        flag = "PASS" if rec["pass"] else "FAIL"
        extra = ""
        if rec["offenders"]:
            extra = " " + ",".join(rec["offenders"])
        sys.stdout.write("%s %s%s\n" % (rid, flag, extra))
        for n in rec.get("notes") or []:
            sys.stdout.write("  %s\n" % n)
    sys.stdout.write("OK\n" if ok else "FAIL\n")


def build_parser():
    p = argparse.ArgumentParser(
        prog="check_regionset.py",
        description=(
            "I1-I9 RegionSet invariant checker. I7 and I7b are two separate "
            "rules keyed on closed360. python3 stdlib only (no jsonschema)."
        ),
    )
    p.add_argument("dump", metavar="dump.json", help="RegionSet JSON dump")
    p.add_argument("--sidecar", metavar="fixture.expected.json",
                   help="fixture expected.json (rejected codes, G5, recoverable)")
    p.add_argument("--component", type=int, metavar="N",
                   help="component index for live[] surfaceCensus (D-130-22(3) addendum)")
    p.add_argument("--skip-recoverable", action="store_true",
                   help="SIDECAR: omit recoverable[] (fixture-wide union pass)")
    p.add_argument("--recoverable-union", metavar="PATH[,PATH...]",
                   help="comma-separated bare dumps for SIDECAR_RECOVERABLE union pool")
    p.add_argument("--schema", metavar="SCHEMA",
                   default=default_schema_path(),
                   help="regionset.schema.json (enum source / documentation; default: alongside this script)")
    p.add_argument("--rule", dest="rule", metavar="I1,I7,I7b",
                   help="comma-separated rule ids to run (default: I1-I9,I7b,G5)")
    p.add_argument("--json", action="store_true",
                   help="machine report on stdout (rule id -> pass/fail -> offending ids)")
    p.add_argument("-v", action="store_true", help="verbose (rule notes on stderr)")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    dump_path = args.dump
    try:
        with open(dump_path, "r", encoding="utf-8") as fh:
            dump = json.load(fh)
    except OSError as e:
        sys.stderr.write("cannot read dump: %s\n" % e)
        return 1
    except json.JSONDecodeError as e:
        sys.stderr.write("invalid JSON in dump: %s\n" % e)
        return 1

    ctx = {
        "enums": {},
        "sidecar": None,
        "sidecar_path": None,
        "component_index": args.component,
        "skip_recoverable": bool(args.skip_recoverable),
        "recoverable_pool": None,
    }
    if args.recoverable_union:
        union_paths = [p.strip() for p in args.recoverable_union.split(",") if p.strip()]
        missing = [p for p in union_paths if not os.path.isfile(p)]
        if missing:
            sys.stderr.write("recoverable-union path(s) not found: %s\n" % ", ".join(missing))
            return 1
        ctx["recoverable_pool"] = _recoverable_pool_from_paths(union_paths)
    schema_path = args.schema
    if schema_path:
        if not os.path.isfile(schema_path):
            sys.stderr.write("schema not found: %s\n" % schema_path)
            return 1
        try:
            ctx["enums"] = load_schema_enums(schema_path)
        except (OSError, json.JSONDecodeError) as e:
            sys.stderr.write("cannot read schema: %s\n" % e)
            return 1

    if args.sidecar:
        try:
            with open(args.sidecar, "r", encoding="utf-8") as fh:
                ctx["sidecar"] = json.load(fh)
            ctx["sidecar_path"] = args.sidecar
        except OSError as e:
            sys.stderr.write("cannot read sidecar: %s\n" % e)
            return 1
        except json.JSONDecodeError as e:
            sys.stderr.write("invalid JSON in sidecar: %s\n" % e)
            return 1

    selected = parse_rules(args.rule)
    if selected is None:
        selected = list(DEFAULT_RULES)
        if ctx["sidecar"] is not None:
            selected.append("SIDECAR")
    elif "SIDECAR" not in selected and ctx["sidecar"] is not None and args.rule is None:
        selected.append("SIDECAR")

    report = run_rules(dump, selected, ctx, verbose=args.v)
    ok = report_ok(report, selected)
    if args.json:
        emit_json(dump_path, report, selected, ok)
    else:
        emit_human(report, selected, ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
