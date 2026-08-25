#!/usr/bin/env python3
"""Verify I1, I2, I5, I7/I7b, I8, I9 across build_fixtures RegionSet + mesh."""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple

from mesh_topology import Component, count_manifold_components, read_binary_stl, vdot, vsub, vunit

HERE = os.path.dirname(os.path.abspath(__file__))

ROLE_ORDER = {"outer": 0, "inner": 1, "capLow": 2, "capHigh": 3}


def owner(tr: List[int], ti: List[int], t: int) -> str:
    if tr[t] >= 0:
        return f"r{tr[t]}"
    if ti[t] >= 0:
        return f"i{ti[t]}"
    return "x"


def tri_area(pts, tri):
    a, b, c = (pts[tri[0]], pts[tri[1]], pts[tri[2]])
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    return 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)


def check_i1(tr, ti, n_tri, errors):
    for t in range(n_tri):
        has_r = tr[t] >= 0
        has_i = ti[t] >= 0
        if has_r == has_i:
            errors.append(f"I1 tri {t}: triRegion/triIsland XOR violated ({tr[t]}, {ti[t]})")
    unassigned = sum(1 for t in range(n_tri) if tr[t] < 0 and ti[t] < 0)
    if unassigned:
        errors.append(f"I1: {unassigned} triangles unassigned")


def check_i2(tr, ti, comp: Component, chains, errors):
    tris = comp.tris
    edges = comp.comp_edges
    tri_edges = comp.tri_edges
    n_tri = len(tris)
    boundary: Set[int] = set()
    edge_info: Dict[int, List[str]] = defaultdict(list)
    for t in range(n_tri):
        o = owner(tr, ti, t)
        for s in range(3):
            eid = tri_edges[t][s]
            edge_info[eid].append(o)
    for eid, uses in edge_info.items():
        if len(uses) == 1:
            boundary.add(eid)
        elif len(uses) == 2 and uses[0] != uses[1]:
            boundary.add(eid)
    chain_edges: Dict[int, int] = {}
    for ci, ch in enumerate(chains):
        for e in ch["meshEdges"]:
            if e in chain_edges:
                errors.append(f"I2 edge {e} in multiple chains ({chain_edges[e]}, {ci})")
            chain_edges[e] = ci
    for eid in boundary:
        if eid not in chain_edges:
            errors.append(f"I2 boundary edge {eid} missing from chains")
    for eid, ci in chain_edges.items():
        if eid not in boundary:
            errors.append(f"I2 interior edge {eid} must not be in chain {ci}")


def check_i5(regions, chains, pts, tris, errors):
    def area(tris_idx):
        return sum(tri_area(pts, tris[t]) for t in tris_idx)

    keys = [(-area(r["tris"]), min(r["tris"])) for r in regions]
    if keys != sorted(keys):
        errors.append("I5: regions[] not sorted by (-area, minLocalTriId)")
    ckeys = [min(c["meshEdges"]) if c["meshEdges"] else 0 for c in chains]
    if ckeys != sorted(ckeys):
        errors.append("I5: chains[] not sorted by minLocalMeshEdgeId")
    for ri, r in enumerate(regions):
        lkeys = [(ROLE_ORDER[lp["role"]], min(lp["chainIdx"]) if lp["chainIdx"] else 0) for lp in r["loops"]]
        if lkeys != sorted(lkeys):
            errors.append(f"I5: region {ri} loops not sorted by ((int)role, minChainIdx)")


def check_i7(regions, chains, pts, errors):
    """I7/I7b role counts PLUS complete loop coverage and distinct CapLow/CapHigh."""
    for r in regions:
        roles = [lp["role"] for lp in r["loops"]]
        used: List[int] = []
        for lp in r["loops"]:
            if len(lp["chainIdx"]) != len(lp["reversed"]):
                errors.append(f"I7 region {r['id']}: chainIdx/reversed length mismatch")
            used.extend(lp["chainIdx"])
        if len(used) != len(set(used)):
            errors.append(f"I7 region {r['id']}: a chain appears in more than one loop ({used})")
        touching = [i for i, c in enumerate(chains) if r["id"] in (c["regA"], c["regB"])]
        if sorted(used) != sorted(touching):
            errors.append(
                f"I7 region {r['id']}: loop coverage {sorted(used)} != touching chains {sorted(touching)}"
            )
        for ci in used:
            if ci < 0 or ci >= len(chains):
                errors.append(f"I7 region {r['id']}: chainIdx {ci} out of range")
            elif r["id"] not in (chains[ci]["regA"], chains[ci]["regB"]):
                errors.append(f"I7 region {r['id']}: loop cites chain {ci} which does not touch the region")
        if r["closed360"]:
            if roles.count("capLow") != 1 or roles.count("capHigh") != 1 or roles.count("outer") != 0:
                errors.append(f"I7b region {r['id']}: invalid loop roles {roles}")
            lo = next((lp for lp in r["loops"] if lp["role"] == "capLow"), None)
            hi = next((lp for lp in r["loops"] if lp["role"] == "capHigh"), None)
            if lo and hi:
                if set(lo["chainIdx"]) & set(hi["chainIdx"]):
                    errors.append(
                        f"I7b region {r['id']}: CapLow and CapHigh share chain(s) "
                        f"{sorted(set(lo['chainIdx']) & set(hi['chainIdx']))}"
                    )
                if not lo["chainIdx"] or not hi["chainIdx"]:
                    errors.append(f"I7b region {r['id']}: empty cap loop")
                loc = tuple(r["ax"]["loc"])
                axis = vunit(tuple(r["ax"]["dir"]))
                vmin, vmax = r["vMin"], r["vMax"]
                mid = 0.5 * (vmin + vmax)

                def mean_v(idxs):
                    vs = []
                    for ci in idxs:
                        for v in chains[ci]["meshVerts"]:
                            vs.append(vdot(vsub(pts[v], loc), axis))
                    return sum(vs) / max(len(vs), 1)

                if lo["chainIdx"] and hi["chainIdx"]:
                    mv_lo, mv_hi = mean_v(lo["chainIdx"]), mean_v(hi["chainIdx"])
                    if mv_lo > mid + 1e-6:
                        errors.append(
                            f"I7b region {r['id']}: CapLow mean v={mv_lo:.6g} is not at vMin={vmin}"
                        )
                    if mv_hi < mid - 1e-6:
                        errors.append(
                            f"I7b region {r['id']}: CapHigh mean v={mv_hi:.6g} is not at vMax={vmax}"
                        )
        else:
            if roles.count("outer") != 1 or roles.count("capLow") or roles.count("capHigh"):
                errors.append(f"I7 region {r['id']}: invalid loop roles {roles}")


def check_id_consistency(regions, tri_region, tri_island, n_islands, errors):
    """triRegion[t] agrees with regions[id].tris (and islands) everywhere."""
    by_id = {}
    for r in regions:
        if r["id"] in by_id:
            errors.append(f"id: duplicate region id {r['id']}")
        by_id[r["id"]] = r
    for i, r in enumerate(regions):
        if r["id"] != i:
            errors.append(f"id: regions[{i}].id is {r['id']}, expected {i}")
        for t in r["tris"]:
            if t < 0 or t >= len(tri_region):
                errors.append(f"id: region {r['id']} tri {t} out of triRegion range")
            elif tri_region[t] != r["id"]:
                errors.append(
                    f"id: region {r['id']} lists tri {t} but triRegion[{t}]={tri_region[t]}"
                )
    for t, rid in enumerate(tri_region):
        if rid < 0:
            continue
        if rid not in by_id:
            errors.append(f"id: triRegion[{t}]={rid} has no regions[id]")
        elif t not in by_id[rid]["tris"]:
            errors.append(f"id: triRegion[{t}]={rid} but region.tris does not contain {t}")
    island_tris: Dict[int, List[int]] = defaultdict(list)
    for t, iid in enumerate(tri_island):
        if iid >= 0:
            island_tris[iid].append(t)
    if n_islands == 0 and island_tris:
        errors.append(f"id: nIslands=0 but triIsland assigns {sorted(island_tris)}")
    if n_islands > 0:
        expected = set(range(n_islands))
        got = set(island_tris)
        if got != expected:
            errors.append(f"id: island ids {sorted(got)} != range({n_islands})")


def check_island_census(rs, exp, errors):
    n_facet = sum(1 for x in rs.get("triIsland", []) if x >= 0)
    if n_facet <= 0:
        return
    census = (exp or {}).get("surfaceCensus") or {}
    fc = (exp or {}).get("faceCount")
    n_reg = len(rs.get("regions") or [])
    if "planar_facets" not in census:
        errors.append("expected.surfaceCensus missing planar_facets for island triangles")
    if isinstance(fc, int) and fc < n_reg + n_facet:
        errors.append(
            f"expected.faceCount={fc} under-counts {n_facet} island tris ({n_reg} regions)"
        )


def check_i8(tr, ti, comp: Component, chains, errors):
    tris = comp.tris
    edges = comp.comp_edges
    tri_edges = comp.tri_edges
    edge_info: Dict[int, List[str]] = defaultdict(list)
    edge_cnt: Dict[int, int] = defaultdict(int)
    for t in range(len(tris)):
        o = owner(tr, ti, t)
        for s in range(3):
            eid = tri_edges[t][s]
            edge_info[eid].append(o)
            edge_cnt[eid] += 1
    boundary_verts: Dict[int, Set[str]] = defaultdict(set)
    open_or_nm_verts: Set[int] = set()
    pair_adj: Dict[frozenset, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for eid, uses in edge_info.items():
        a, b = edges[eid]
        if edge_cnt[eid] != 2:
            open_or_nm_verts.add(a)
            open_or_nm_verts.add(b)
        owners = set(uses)
        if edge_cnt[eid] == 1:
            owners.add("x")
        if len(owners) < 2 and edge_cnt[eid] != 1:
            continue
        if len(owners) >= 2 or edge_cnt[eid] == 1:
            boundary_verts[a].update(owners)
            boundary_verts[b].update(owners)
            pair = frozenset(owners)
            pair_adj[pair][a] += 1
            pair_adj[pair][b] += 1
    split = {v for v, ow in boundary_verts.items() if len(ow) >= 3}
    # I8: also split at vertices whose owner-pair degree != 2, and at nm vertices.
    deg_split: Set[int] = set()
    for adj in pair_adj.values():
        for v, deg in adj.items():
            if deg != 2:
                deg_split.add(v)
    must_split = split | {v for v in open_or_nm_verts if v in deg_split}
    for ci, ch in enumerate(chains):
        verts = ch["meshVerts"]
        n = len(verts)
        if n == 0:
            continue
        interior = range(1, n) if ch["closedLoop"] else range(1, n - 1)
        for i in interior:
            v = verts[i]
            if v in must_split:
                errors.append(f"I8 chain {ci} not split at vertex {v} (>=3 owners or junction)")


def check_i4_surface(rs, pts, tris, errors, eps: float = 1e-3):
    """Region vertices must lie on the declared surface; vMin/vMax match the mesh (E1)."""
    for r in rs["regions"]:
        loc = tuple(r["ax"]["loc"])
        direction = vunit(tuple(r["ax"]["dir"]))
        verts = set()
        for t in r["tris"]:
            if 0 <= t < len(tris):
                verts.update(tris[t])
        if not verts:
            errors.append(f"I4 region {r['id']}: empty tris")
            continue
        if r["type"] == "plane":
            for v in verts:
                h = abs(vdot(vsub(pts[v], loc), direction))
                if h > eps:
                    errors.append(
                        f"I4 region {r['id']} plane: vertex {v} is {h:.4g} mm off the plane"
                    )
        elif r["type"] == "cylinder":
            rad = r["radius"]
            vs = []
            for v in verts:
                w = vsub(pts[v], loc)
                vs.append(vdot(w, direction))
                radial = vsub(w, (direction[0] * vs[-1], direction[1] * vs[-1], direction[2] * vs[-1]))
                rho = (radial[0] ** 2 + radial[1] ** 2 + radial[2] ** 2) ** 0.5
                if abs(rho - rad) > eps:
                    errors.append(
                        f"I4 region {r['id']} cylinder: vertex {v} |ρ−R|={abs(rho - rad):.4g} mm (R={rad})"
                    )
            if vs:
                if abs(min(vs) - r["vMin"]) > eps:
                    errors.append(
                        f"I4 region {r['id']}: vMin={r['vMin']} != mesh min axial {min(vs)}"
                    )
                if abs(max(vs) - r["vMax"]) > eps:
                    errors.append(
                        f"I4 region {r['id']}: vMax={r['vMax']} != mesh max axial {max(vs)}"
                    )


def check_i9(regions, errors):
    for r in regions:
        dv = r["dVolPredicted"]
        if r["outwardNormal"] and dv < -1e-9:
            errors.append(f"I9 region {r['id']}: outwardNormal=true but dVolPredicted={dv} < 0")
        if not r["outwardNormal"] and dv > 1e-9:
            errors.append(f"I9 region {r['id']}: outwardNormal=false but dVolPredicted={dv} > 0")


def check_tangent_fillet(rs, errors):
    """G1: every FilletStrip region has a tangent chain on a fillet|neighbour-plane interface.

    Dropped B's `len(meshEdges) >= 2` (wrong on a genuine 2-tri strip). Keyed off
    origin==filletStrip, never the filename. Requires regA,regB >= 0 and the other
    side in {filletNbrA, filletNbrB}.
    """
    chains = rs["chains"]
    for r in rs["regions"]:
        if r.get("origin") != "filletStrip":
            continue
        nbrs = {r["filletNbrA"], r["filletNbrB"]} - {-1}
        found = []
        for i, c in enumerate(chains):
            if not c.get("tangent"):
                continue
            if c["regA"] < 0 or c["regB"] < 0:
                continue
            sides = {c["regA"], c["regB"]}
            if r["id"] not in sides:
                continue
            if not (sides & nbrs):
                continue
            found.append(i)
        if not found:
            errors.append(
                f"G1 fillet region {r['id']}: no tangent fillet|neighbour-plane chain "
                f"(filletNbrA={r['filletNbrA']}, filletNbrB={r['filletNbrB']})"
            )


def check_r8(rs, exp, errors):
    """DECISION §5 / audit R8: explode un-collapses a collapsed gp_Circ cap to N chords
    and rebuilds the planar inner from those chords."""
    block = exp.get("r8")
    if not block:
        errors.append("R8: expected.json missing 'r8' un-collapse block")
        return
    cyls = [r for r in rs["regions"] if r.get("closed360")]
    planes = [r for r in rs["regions"] if r.get("type") == "plane"]
    if len(cyls) != 1 or len(planes) != 1:
        errors.append(f"R8: expected 1 closed360 + 1 plane, got cyls={len(cyls)} planes={len(planes)}")
        return
    cyl, plane = cyls[0], planes[0]
    chains = rs["chains"]
    lo = next((lp for lp in cyl["loops"] if lp["role"] == "capLow"), None)
    hi = next((lp for lp in cyl["loops"] if lp["role"] == "capHigh"), None)
    inner = next((lp for lp in plane["loops"] if lp["role"] == "inner"), None)
    if not lo or not hi or not inner:
        errors.append("R8: missing CapLow/CapHigh/plane Inner")
        return
    if lo["chainIdx"] != block.get("capLowChainIdx"):
        errors.append(f"R8: capLowChainIdx {lo['chainIdx']} != expected {block.get('capLowChainIdx')}")
    if hi["chainIdx"] != block.get("capHighChainIdx"):
        errors.append(f"R8: capHighChainIdx {hi['chainIdx']} != expected {block.get('capHighChainIdx')}")
    if inner["chainIdx"] != hi["chainIdx"]:
        errors.append(
            f"R8: plane inner {inner['chainIdx']} must be the cylinder CapHigh {hi['chainIdx']} "
            "(neighbour inner-wire rebuild)"
        )
    if not block.get("planeInnerUsesCapHigh"):
        errors.append("R8: planeInnerUsesCapHigh must be true")
    n = cyl["nSides"]
    if block.get("nSides") != n:
        errors.append(f"R8: nSides {block.get('nSides')} != region nSides {n}")
    if block.get("explodeUncollapsesGpCircToNChords") != n:
        errors.append(
            f"R8: explodeUncollapsesGpCircToNChords={block.get('explodeUncollapsesGpCircToNChords')} "
            f"!= nSides={n}"
        )
    for role, lp in (("capLow", lo), ("capHigh", hi)):
        n_chords = 0
        for ci in lp["chainIdx"]:
            ch = chains[ci]
            if not ch.get("closedLoop"):
                errors.append(f"R8: {role} chain {ci} is not closedLoop (P1 collapsed gp_Circ ring)")
            n_chords += len(ch["meshEdges"])
        if n_chords != n:
            errors.append(f"R8: {role} has {n_chords} chords, need nSides={n} (un-collapse target)")


def d45_budget(vol: float, dvols: List[float]) -> Tuple[float, float, float]:
    signed = sum(dvols)
    abs_sum = sum(abs(x) for x in dvols)
    return signed, abs_sum, max(1e-4 * abs(vol), 3.0 * abs_sum)


def check_expected(rs, exp, comp: Component, name: str, errors):
    """Sidecar *.expected.json is a gate, not a sticker."""
    if not exp:
        errors.append("expected.json: missing sidecar")
        return
    if exp.get("closed") != (comp.open == 0):
        errors.append(f"expected.closed={exp.get('closed')} but mesh open={comp.open}")
    if comp.conflict != 0:
        errors.append(f"mesh conflict={comp.conflict} (mis-winding); must be 0 (open>0 is ok)")
    if comp.non_manifold != 0:
        errors.append(f"mesh nonManifold={comp.non_manifold}; must be 0")
    for key in ("faceCount", "surfaceCensus", "rung", "rounds", "meshTriangles", "volumeMM3"):
        if key not in exp:
            errors.append(f"expected.json missing '{key}'")
    if exp.get("meshTriangles") != len(comp.tris):
        errors.append(f"expected.meshTriangles={exp.get('meshTriangles')} != mesh {len(comp.tris)}")
    fc = exp.get("faceCount")
    if not isinstance(fc, int) or fc < len(rs["regions"]):
        errors.append(f"expected.faceCount={fc} is not >= region count {len(rs['regions'])}")
    census = exp.get("surfaceCensus")
    if not isinstance(census, dict) or not census:
        errors.append("expected.surfaceCensus missing or empty")
    exp_regs = exp.get("regions")
    if not isinstance(exp_regs, list) or len(exp_regs) != len(rs["regions"]):
        errors.append("expected.regions length != RegionSet regions")
    else:
        for i, (er, r) in enumerate(zip(exp_regs, rs["regions"])):
            if er.get("id") != r["id"]:
                errors.append(f"expected.regions[{i}].id={er.get('id')} != {r['id']}")
            if "builtAs" not in er:
                errors.append(f"expected.regions[{i}] missing builtAs")
    signed, abs_sum, budget = d45_budget(comp.vol, [r["dVolPredicted"] for r in rs["regions"]])
    if "dVolPredicted" in exp and abs(exp["dVolPredicted"] - signed) > 1e-6 * max(1.0, abs(signed)):
        errors.append(
            f"D4.5: expected.dVolPredicted={exp['dVolPredicted']} != signed sum {signed}"
        )
    if "dVolPredAbs" in exp and abs(exp["dVolPredAbs"] - abs_sum) > 1e-6 * max(1.0, abs_sum):
        errors.append(f"D4.5: expected.dVolPredAbs={exp['dVolPredAbs']} != Σ|dVol| {abs_sum}")
    if "volumeBudgetMM3" not in exp:
        errors.append("D4.5: expected.volumeBudgetMM3 missing")
    elif abs(exp["volumeBudgetMM3"] - budget) > 1e-6 * max(1.0, budget):
        errors.append(
            f"D4.5: expected.volumeBudgetMM3={exp['volumeBudgetMM3']} != max(1e-4·|vol|, 3·Σ|dVol|)={budget}"
        )
    if name == "r8_360_explode_caps" or "r8" in exp:
        check_r8(rs, exp, errors)
    if name == "r2_forced":
        cyls = [r for r in rs["regions"] if r.get("closed360")]
        if not cyls or any(r.get("nSides", 0) != 3 for r in cyls):
            errors.append("r2_forced: closed360 cylinders must be nSides=3 (sagitta/R=0.5)")
        if len(cyls) < 3:
            errors.append("r2_forced: need >=3 closed360 cylinders (2-round R1 cap / sequential explode)")
        sag = 1.0 - math.cos(math.pi / 3.0)
        if abs((exp or {}).get("sagittaOverR", 0) - sag) > 1e-6:
            errors.append("r2_forced: expected.sagittaOverR must be 1-cos(pi/3)=0.5 (vs hole 0.293)")
    check_island_census(rs, exp, errors)


def check_mesh_ids(rs, comp: Component, errors):
    pts, tris = comp.pts, comp.tris
    n_tri, n_edge = len(tris), comp.n_edge
    for r in rs["regions"]:
        for t in r["tris"]:
            if t < 0 or t >= n_tri:
                errors.append(f"mesh: region {r['id']} tri {t} out of range")
    for ci, c in enumerate(rs["chains"]):
        for e in c["meshEdges"]:
            if e < 0 or e >= n_edge:
                errors.append(f"mesh: chain {ci} edge {e} out of range")
        for v in c["meshVerts"]:
            if v < 0 or v >= len(pts):
                errors.append(f"mesh: chain {ci} vert {v} out of range")
    if len(rs["triRegion"]) != n_tri or len(rs["triIsland"]) != n_tri:
        errors.append("mesh: triRegion/triIsland length != mesh triangle count")


def validate_file(rs_path: str) -> Tuple[bool, List[str]]:
    errors: List[str] = []
    with open(rs_path, encoding="utf-8") as f:
        rs = json.load(f)
    stl_path = rs_path.replace(".regionset.json", ".stl")
    if not os.path.isfile(stl_path):
        errors.append(f"missing mesh {stl_path}")
        return False, errors
    pts, tris = read_binary_stl(stl_path)
    comp = Component(pts=list(pts), tris=list(tris))
    comp.build_topology()
    tr, ti = rs["triRegion"], rs["triIsland"]
    check_mesh_ids(rs, comp, errors)
    ncomp = count_manifold_components(comp)
    if ncomp != 1:
        errors.append(f"mesh: {ncomp} manifold components; comps[0] cannot consume these ids")
    check_i1(tr, ti, len(tris), errors)
    check_i2(tr, ti, comp, rs["chains"], errors)
    check_i5(rs["regions"], rs["chains"], pts, tris, errors)
    check_id_consistency(rs["regions"], tr, ti, rs.get("nIslands", 0), errors)
    check_i7(rs["regions"], rs["chains"], pts, errors)
    check_i8(tr, ti, comp, rs["chains"], errors)
    check_i9(rs["regions"], errors)
    check_i4_surface(rs, pts, tris, errors)
    check_tangent_fillet(rs, errors)
    name = os.path.basename(os.path.dirname(os.path.abspath(rs_path)))
    exp_path = rs_path.replace(".regionset.json", ".expected.json")
    exp = None
    if not os.path.isfile(exp_path):
        errors.append(f"missing sidecar {exp_path}")
    else:
        with open(exp_path, encoding="utf-8") as f:
            exp = json.load(f)
        check_expected(rs, exp, comp, name, errors)
    return len(errors) == 0, errors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*", default=[HERE])
    args = ap.parse_args()
    files: List[str] = []
    for root in args.paths:
        if os.path.isdir(root):
            files.extend(sorted(glob.glob(os.path.join(root, "**", "*.regionset.json"), recursive=True)))
        else:
            files.append(root)
    failed = 0
    for fp in files:
        ok, errors = validate_file(fp)
        if ok:
            print(f"PASS {fp}")
        else:
            failed += 1
            print(f"FAIL {fp}")
            for e in errors:
                print(f"  {e}")
    print(f"\n{len(files) - failed}/{len(files)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
