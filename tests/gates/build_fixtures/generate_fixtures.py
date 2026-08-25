#!/usr/bin/env python3
"""Generate the eleven blocking P2 build_fixtures (mesh + RegionSet + expected).

Each mesh is one manifold component (harness loadMesh comps[0] consumes the JSON
ids as-is). Loops honour I7/I7b in full: every touching chain appears in exactly
one loop; closed360 CapLow/CapHigh are distinct chains at vMin/vMax.
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

from mesh_topology import (
    Component,
    MeshBuilder,
    Tri,
    Vec3,
    chain_g1,
    count_manifold_components,
    cyl_dvol,
    earclip_xy,
    extract_chains,
    fillet_dvol,
    lowest_id_xdir,
    plane_dvol,
    poly_area2,
    read_binary_stl,
    region_vertices,
    tri_area,
    vdot,
    vnorm,
    vscale,
    vsub,
    vunit,
    write_binary_stl,
)

HERE = os.path.dirname(os.path.abspath(__file__))
TWOPI = 2.0 * math.pi
ROLE_ORDER = {"outer": 0, "inner": 1, "capLow": 2, "capHigh": 3}


def ax3(loc: Vec3, direction: Vec3, xdir: Vec3) -> Dict[str, Any]:
    return {"loc": list(loc), "dir": list(direction), "xdir": list(xdir)}


def loop(chain_idx: Sequence[int], reversed_: Sequence[int], role: str) -> Dict[str, Any]:
    return {"chainIdx": list(chain_idx), "reversed": list(reversed_), "role": role}


def base_region(
    rid: int,
    *,
    stype: str,
    origin: str,
    tris: Sequence[int],
    loops: List[Dict[str, Any]],
    closed360: bool,
    outward: bool,
    radius: float = 0.0,
    umin: float = 0.0,
    umax: float = 0.0,
    vmin: float = 0.0,
    vmax: float = 0.0,
    axis: Optional[Dict[str, Any]] = None,
    n_sides: int = 0,
    chord: float = 0.0,
    dvol: float = 0.0,
    fil_a: int = -1,
    fil_b: int = -1,
    max_dev: float = 0.0,
    exp: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    rec = {
        "id": rid,
        "type": stype,
        "origin": origin,
        "ax": axis or ax3((0, 0, 0), (0, 0, 1), (1, 0, 0)),
        "radius": radius,
        "uMin": umin,
        "uMax": umax,
        "vMin": vmin,
        "vMax": vmax,
        "closed360": closed360,
        "outwardNormal": outward,
        "tris": sorted(tris),
        "loops": loops,
        "maxVertexDev": max_dev,
        "rmsVertexDev": max_dev * 0.5,
        "chordSagitta": chord,
        "nSides": n_sides,
        "dVolPredicted": dvol,
        "maxVertexSnap": 0.0,
        "reject": "none",
        "builtAs": "notBuilt",
        "filletNbrA": fil_a,
        "filletNbrB": fil_b,
    }
    if exp:
        rec["_exp"] = exp
    return rec


def ring(m: MeshBuilder, r: float, z: float, n: int, phase: float = 0.0) -> List[int]:
    return [
        m.v(r * math.cos(phase + TWOPI * i / n), r * math.sin(phase + TWOPI * i / n), z)
        for i in range(n)
    ]


def stitch_cyl(m: MeshBuilder, bot: Sequence[int], top: Sequence[int], outward: bool = True) -> List[int]:
    n = len(bot)
    ids: List[int] = []
    for i in range(n):
        j = (i + 1) % n
        if outward:
            ids.extend(m.quad(bot[i], bot[j], top[j], top[i]))
        else:
            ids.extend(m.quad(bot[i], top[i], top[j], bot[j]))
    return ids


def stitch_annulus(m: MeshBuilder, outer: Sequence[int], inner: Sequence[int], normal_up: bool) -> List[int]:
    n = len(outer)
    ids: List[int] = []
    for i in range(n):
        j = (i + 1) % n
        if normal_up:
            ids.extend(m.quad(outer[i], inner[i], inner[j], outer[j]))
        else:
            ids.extend(m.quad(outer[i], outer[j], inner[j], inner[i]))
    return ids


def cap_fan(m: MeshBuilder, loop_ids: Sequence[int], reverse: bool) -> List[int]:
    if len(loop_ids) == 3:
        a, b, c = loop_ids
        return [m.tri(a, c, b) if reverse else m.tri(a, b, c)]
    ids: List[int] = []
    origin = loop_ids[0]
    for i in range(1, len(loop_ids) - 1):
        if reverse:
            ids.append(m.tri(origin, loop_ids[i + 1], loop_ids[i]))
        else:
            ids.append(m.tri(origin, loop_ids[i], loop_ids[i + 1]))
    return ids


def chord(r: float, n_sides: int) -> float:
    if n_sides < 2:
        return 0.0
    return r * (1.0 - math.cos(math.pi / n_sides))


def touching(chains: Sequence[dict], rid: int) -> List[int]:
    return [i for i, c in enumerate(chains) if rid in (c["regA"], c["regB"])]


def rev_flags(chain_idxs: Sequence[int], rid: int, chains: Sequence[dict]) -> List[int]:
    out: List[int] = []
    for ci in chain_idxs:
        c = chains[ci]
        out.append(0 if c["regA"] == rid else 1)
    return out


def chain_mean_v(ch: dict, pts: Sequence[Vec3], loc: Vec3, axis: Vec3) -> float:
    ax = vunit(axis)
    vs = [vdot(vsub(pts[v], loc), ax) for v in ch["meshVerts"]]
    return sum(vs) / max(len(vs), 1)


def group_connected(chain_idxs: Sequence[int], chains: Sequence[dict]) -> List[List[int]]:
    parent = {i: i for i in chain_idxs}

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def unite(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    verts: Dict[int, List[int]] = defaultdict(list)
    for i in chain_idxs:
        for v in chains[i]["meshVerts"]:
            verts[v].append(i)
    for ids in verts.values():
        for k in range(1, len(ids)):
            unite(ids[0], ids[k])
    groups: Dict[int, List[int]] = defaultdict(list)
    for i in chain_idxs:
        groups[find(i)].append(i)
    return [sorted(g) for g in groups.values()]


def loop_shoelace(chain_idxs: Sequence[int], chains: Sequence[dict], pts: Sequence[Vec3],
                  loc: Vec3, n: Vec3, xdir: Vec3) -> float:
    n = vunit(n)
    x = vunit(xdir)
    y = (
        n[1] * x[2] - n[2] * x[1],
        n[2] * x[0] - n[0] * x[2],
        n[0] * x[1] - n[1] * x[0],
    )
    seq: List[int] = []
    remaining = list(chain_idxs)
    if not remaining:
        return 0.0
    cur = remaining.pop(0)
    seq.extend(chains[cur]["meshVerts"])
    guard = 0
    while remaining and guard < 64:
        guard += 1
        tail = seq[-1]
        head = seq[0]
        found = None
        for i, ci in enumerate(remaining):
            vs = chains[ci]["meshVerts"]
            if vs[0] == tail:
                seq.extend(vs[1:])
                found = i
                break
            if vs[-1] == tail:
                seq.extend(list(reversed(vs[:-1])))
                found = i
                break
            if vs[0] == head:
                seq = list(reversed(vs[1:])) + seq
                found = i
                break
            if vs[-1] == head:
                seq = vs[:-1] + seq
                found = i
                break
        if found is None:
            break
        remaining.pop(found)
    if seq and seq[0] == seq[-1]:
        seq = seq[:-1]
    poly = []
    for v in seq:
        d = vsub(pts[v], loc)
        poly.append((vdot(d, x), vdot(d, y)))
    return poly_area2(poly)


def assign_loops(region: dict, chains: Sequence[dict], pts: Sequence[Vec3]) -> None:
    rid = region["id"]
    idxs = touching(chains, rid)
    if not idxs:
        raise RuntimeError(f"region {rid} has no touching chains")
    if region["closed360"]:
        loc = tuple(region["ax"]["loc"])
        axis = tuple(region["ax"]["dir"])
        vmin, vmax = region["vMin"], region["vMax"]
        mid = 0.5 * (vmin + vmax)
        cap_lo: List[int] = []
        cap_hi: List[int] = []
        inner: List[int] = []
        for i in idxs:
            vs = [vdot(vsub(pts[v], loc), vunit(axis)) for v in chains[i]["meshVerts"]]
            lo, hi = min(vs), max(vs)
            if hi <= mid + 1e-6:
                cap_lo.append(i)
            elif lo >= mid - 1e-6:
                cap_hi.append(i)
            else:
                inner.append(i)
        if not cap_lo or not cap_hi:
            raise RuntimeError(
                f"region {rid} closed360 needs distinct caps, got lo={cap_lo} hi={cap_hi} inner={inner}"
            )
        if set(cap_lo) & set(cap_hi):
            raise RuntimeError(f"region {rid} CapLow and CapHigh share a chain")
        cap_lo.sort()
        cap_hi.sort()
        loops = [
            loop(cap_lo, rev_flags(cap_lo, rid, chains), "capLow"),
            loop(cap_hi, rev_flags(cap_hi, rid, chains), "capHigh"),
        ]
        for g in group_connected(inner, chains):
            loops.append(loop(g, rev_flags(g, rid, chains), "inner"))
        region["loops"] = loops
        return

    loc = tuple(region["ax"]["loc"])
    n = tuple(region["ax"]["dir"])
    xdir = tuple(region["ax"]["xdir"])
    groups = group_connected(idxs, chains)
    scored = []
    for g in groups:
        scored.append((abs(loop_shoelace(g, chains, pts, loc, n, xdir)), min(g), g))
    scored.sort(reverse=True)
    outer = scored[0][2]
    loops = [loop(outer, rev_flags(outer, rid, chains), "outer")]
    for _a, _m, g in scored[1:]:
        loops.append(loop(g, rev_flags(g, rid, chains), "inner"))
    region["loops"] = loops


def sort_and_remap(
    regions: List[dict],
    chains: List[dict],
    tri_region: List[int],
    pts: Sequence[Vec3],
    tris: Sequence[Tri],
) -> Tuple[List[dict], List[dict], List[int]]:
    def key(r: dict) -> Tuple[float, int]:
        area = sum(tri_area(pts, tris[t]) for t in r["tris"])
        return (-area, min(r["tris"]))

    order = sorted(range(len(regions)), key=lambda i: key(regions[i]))
    old_to_new = {old: new for new, old in enumerate(order)}
    new_regions = []
    for new, old in enumerate(order):
        r = regions[old]
        r["id"] = new
        if r["filletNbrA"] >= 0:
            r["filletNbrA"] = old_to_new[r["filletNbrA"]]
        if r["filletNbrB"] >= 0:
            r["filletNbrB"] = old_to_new[r["filletNbrB"]]
        new_regions.append(r)
    tri_region = [old_to_new[t] if t >= 0 else -1 for t in tri_region]
    for c in chains:
        if c["regA"] >= 0:
            c["regA"] = old_to_new[c["regA"]]
        if c["regB"] >= 0:
            c["regB"] = old_to_new[c["regB"]]
    chains = sorted(chains, key=lambda c: min(c["meshEdges"]) if c["meshEdges"] else 0)
    return new_regions, chains, tri_region


def finalize_regionset(
    comp: Component,
    regions: List[dict],
    chains: List[dict],
    *,
    n_islands: int,
    tri_region: List[int],
    tri_island: List[int],
) -> Dict[str, Any]:
    for r in regions:
        r["loops"].sort(key=lambda lp: (ROLE_ORDER[lp["role"]], min(lp["chainIdx"]) if lp["chainIdx"] else 0))
    planes = sum(1 for r in regions if r["type"] == "plane")
    cylinders = sum(1 for r in regions if r["type"] == "cylinder")
    fillets = sum(1 for r in regions if r["origin"] == "filletStrip")
    radii = sorted({round(r["radius"], 6) for r in regions if r["type"] == "cylinder" and r["radius"] > 0})
    dvol = sum(r["dVolPredicted"] for r in regions)
    return {
        "compRoot": 0,
        "regions": regions,
        "rejected": [],
        "chains": chains,
        "triRegion": tri_region,
        "triIsland": tri_island,
        "nIslands": n_islands,
        "stats": {
            "planes": planes,
            "cylinders": cylinders,
            "fillets": fillets,
            "rejected": 0,
            "facetIslands": n_islands,
            "facetTriangles": sum(1 for x in tri_island if x >= 0),
            "distinctRadii": len(radii),
            "maxVertexDev": max((r["maxVertexDev"] for r in regions), default=0.0),
            "maxEdgeTol": 0.0,
            "dVolPredicted": dvol,
        },
    }


def budget(vol: float, dvols: Sequence[float]) -> float:
    return max(1e-4 * abs(vol), 3.0 * sum(abs(x) for x in dvols))


def snap_region_to_mesh(r: dict, comp: Component) -> None:
    """Re-derive loc/vMin/vMax/radius/dVol from the region's own triangles (E1)."""
    tris = r["tris"]
    if not tris:
        return
    verts = region_vertices(comp, tris)
    pts = [comp.pts[v] for v in verts]
    if r["type"] == "plane":
        n = vunit(tuple(r["ax"]["dir"]))
        loc = pts[0]
        r["ax"]["loc"] = list(loc)
        r["vMin"] = 0.0
        r["vMax"] = 0.0
        r["dVolPredicted"] = plane_dvol(comp.pts, [comp.tris[t] for t in tris], n, loc)
        devs = [abs(vdot(vsub(p, loc), n)) for p in pts]
        r["maxVertexDev"] = max(devs) if devs else 0.0
        r["rmsVertexDev"] = (sum(d * d for d in devs) / len(devs)) ** 0.5 if devs else 0.0
        return
    if r["type"] == "cylinder":
        ax = vunit(tuple(r["ax"]["dir"]))
        loc0 = tuple(r["ax"]["loc"])
        # Keep axis direction; place loc at the axis point nearest the vertex centroid.
        c = (
            sum(p[0] for p in pts) / len(pts),
            sum(p[1] for p in pts) / len(pts),
            sum(p[2] for p in pts) / len(pts),
        )
        loc = (
            loc0[0] + ax[0] * vdot(vsub(c, loc0), ax),
            loc0[1] + ax[1] * vdot(vsub(c, loc0), ax),
            loc0[2] + ax[2] * vdot(vsub(c, loc0), ax),
        )
        vs = [vdot(vsub(p, loc), ax) for p in pts]
        r["vMin"] = min(vs)
        r["vMax"] = max(vs)
        rhos = []
        for p in pts:
            w = vsub(p, loc)
            radv = vsub(w, vscale(ax, vdot(w, ax)))
            rhos.append(vnorm(radv))
        # Keep the authored radius; only snap frame/extent to the mesh.
        rad = r["radius"]
        r["ax"]["loc"] = list(loc)
        r["ax"]["xdir"] = list(lowest_id_xdir(comp.pts, verts, loc, ax))
        n_sides = r["nSides"] or 3
        r["chordSagitta"] = chord(rad, n_sides)
        axial = r["vMax"] - r["vMin"]
        n_bands = n_sides if r["closed360"] else max(1, round(n_sides * (r["uMax"] - r["uMin"]) / TWOPI))
        r["dVolPredicted"] = cyl_dvol(rad, n_sides, axial, r["outwardNormal"], n_bands=n_bands)
        if r["origin"] == "filletStrip":
            r["dVolPredicted"] = fillet_dvol(rad, axial, r["outwardNormal"], r["uMax"] - r["uMin"] or math.pi / 2)
        devs = [abs(rho - rad) for rho in rhos]
        r["maxVertexDev"] = max(devs) if devs else 0.0
        r["rmsVertexDev"] = (sum(d * d for d in devs) / len(devs)) ** 0.5 if devs else 0.0


def encode_r8(rs: Dict[str, Any]) -> Dict[str, Any]:
    """P1 emits one closed cap ring of N mesh chords; P2 collapses that to gp_Circ
    and on explode must un-collapse to those N chords and rebuild the planar inner."""
    cyl = next(r for r in rs["regions"] if r["closed360"])
    plane = next(r for r in rs["regions"] if r["type"] == "plane")
    lo = next(lp for lp in cyl["loops"] if lp["role"] == "capLow")
    hi = next(lp for lp in cyl["loops"] if lp["role"] == "capHigh")
    inner = next(lp for lp in plane["loops"] if lp["role"] == "inner")
    chains = rs["chains"]
    hi_ci = hi["chainIdx"][0]
    n_chords = len(chains[hi_ci]["meshEdges"])
    return {
        "nSides": cyl["nSides"],
        "capLowChainIdx": lo["chainIdx"],
        "capHighChainIdx": hi["chainIdx"],
        "capLowClosedLoop": all(chains[i]["closedLoop"] for i in lo["chainIdx"]),
        "capHighClosedLoop": all(chains[i]["closedLoop"] for i in hi["chainIdx"]),
        "capLowChordCount": sum(len(chains[i]["meshEdges"]) for i in lo["chainIdx"]),
        "capHighChordCount": n_chords,
        "planeInnerChainIdx": inner["chainIdx"],
        "planeInnerUsesCapHigh": inner["chainIdx"] == hi["chainIdx"],
        "explodeUncollapsesGpCircToNChords": n_chords,
    }


def emit(
    name: str,
    m: MeshBuilder,
    tri_region: List[int],
    tri_island: List[int],
    regions: List[dict],
    expected: Dict[str, Any],
    n_islands: int,
) -> None:
    d = os.path.join(HERE, name)
    os.makedirs(d, exist_ok=True)
    stl_path = os.path.join(d, f"{name}.stl")
    m.orient_manifold()
    write_binary_stl(stl_path, m.pts, m.tris)
    pts, tris = read_binary_stl(stl_path)
    if len(tris) != len(m.tris):
        raise RuntimeError(f"{name}: STL reload dropped triangles {len(m.tris)} -> {len(tris)}")
    comp = Component(pts=list(pts), tris=list(tris))
    comp.build_topology()
    ncomp = count_manifold_components(comp)
    if ncomp != 1:
        raise RuntimeError(f"{name}: harness would split into {ncomp} components (need 1)")
    if comp.conflict != 0:
        raise RuntimeError(f"{name}: conflict={comp.conflict} after orient (mis-winding)")
    if len(tri_region) != len(tris) or len(tri_island) != len(tris):
        raise RuntimeError(f"{name}: occupancy length mismatch")
    chains = extract_chains(comp, tri_region, tri_island)
    regions, chains, tri_region = sort_and_remap(regions, chains, tri_region, comp.pts, comp.tris)
    for r in regions:
        snap_region_to_mesh(r, comp)
        assign_loops(r, chains, comp.pts)
    for c in chains:
        c["tangent"] = chain_g1(c, regions)
    exp_regions = []
    for r in regions:
        extra = r.pop("_exp", {})
        rec = {"id": r["id"], "builtAs": extra.get("builtAs", "single")}
        rec.update({k: v for k, v in extra.items() if k != "builtAs"})
        exp_regions.append(rec)
    rs = finalize_regionset(
        comp, regions, chains, n_islands=n_islands, tri_region=tri_region, tri_island=tri_island
    )
    dvols = [r["dVolPredicted"] for r in regions]
    signed = sum(dvols)
    abs_sum = sum(abs(x) for x in dvols)
    expected.setdefault("fixture", name)
    expected.setdefault("meshTriangles", len(tris))
    expected.setdefault("closed", comp.open == 0)
    expected.setdefault("volumeMM3", comp.vol)
    # D4.5: signed sum and magnitude sum are distinct; budget uses magnitudes.
    expected["dVolPredicted"] = signed
    expected["dVolPredAbs"] = abs_sum
    expected["volumeBudgetMM3"] = budget(comp.vol, dvols)
    expected["regions"] = exp_regions
    if name == "r8_360_explode_caps":
        expected["r8"] = encode_r8(rs)
    with open(os.path.join(d, f"{name}.regionset.json"), "w", encoding="utf-8") as f:
        json.dump(rs, f, indent=2)
        f.write("\n")
    with open(os.path.join(d, f"{name}.expected.json"), "w", encoding="utf-8") as f:
        json.dump(expected, f, indent=2)
        f.write("\n")
    print(
        f"  {name}: tris={len(tris)} open={comp.open} conflict={comp.conflict} "
        f"nm={comp.non_manifold} comps=1 chains={len(chains)} regs={len(regions)} "
        f"closed={expected['closed']}"
    )


# ---------------------------------------------------------------------------
# F01 full_360_hole — closed square washer, n=4 inner cylinder
# ---------------------------------------------------------------------------
def fixture_full_360_hole() -> None:
    n, r_in, r_out, h = 4, 4.0, 8.0, 4.0
    m = MeshBuilder()
    ib, it = ring(m, r_in, 0.0, n), ring(m, r_in, h, n)
    ob, ot = ring(m, r_out, 0.0, n), ring(m, r_out, h, n)
    cyl = stitch_cyl(m, ib, it, outward=True)
    top = stitch_annulus(m, ot, it, normal_up=True)
    bot = stitch_annulus(m, ob, ib, normal_up=False)
    outer = stitch_cyl(m, ob, ot, outward=False)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyl:
        tr[t] = 0
    for t in top:
        tr[t] = 1
    for t in bot:
        tr[t] = 2
    for t in outer:
        ti[t] = 0
    dvol = cyl_dvol(r_in, n, h, False)
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=cyl, loops=[], closed360=True, outward=False,
            radius=r_in, umin=0.0, umax=TWOPI, vmin=0.0, vmax=h,
            axis=ax3(loc, axis, xdir), n_sides=n, chord=chord(r_in, n), dvol=dvol,
            exp={"builtAs": "seamed360"},
        ),
        base_region(
            1, stype="plane", origin="planeGrow", tris=top, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, h), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in top], (0.0, 0.0, 1.0), (0.0, 0.0, h)),
            exp={"builtAs": "single"},
        ),
        base_region(
            2, stype="plane", origin="planeGrow", tris=bot, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in bot], (0.0, 0.0, -1.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("full_360_hole", m, tr, ti, regs, {
        "faceCount": 3 + len(outer),
        "surfaceCensus": {"cylinder": 1, "plane": 2, "planar_facets": len(outer)},
        "rung": 1, "rounds": 0,
        "note": "Closed washer; inner closed360 cylinder CapLow≠CapHigh, no Outer, no seam.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F02 plate_half_hole — one plate with a ~180° notch, extruded (single solid)
# ---------------------------------------------------------------------------
def fixture_plate_half_hole() -> None:
    w, depth, rad, thick = 16.0, 10.0, 4.0, 2.0
    n_arc = 4
    cx, cy = w / 2.0, depth
    poly: List[Tuple[float, float]] = [(0.0, 0.0), (w, 0.0), (w, depth)]
    for i in range(n_arc + 1):
        ang = 0.0 - math.pi * i / n_arc  # 0 -> -pi (into the plate)
        poly.append((cx + rad * math.cos(ang), cy + rad * math.sin(ang)))
    poly.append((0.0, depth))
    m = MeshBuilder()
    bot = [m.v(x, y, 0.0) for x, y in poly]
    top = [m.v(x, y, thick) for x, y in poly]
    ears = earclip_xy(poly)
    if len(ears) < 3:
        raise RuntimeError(f"plate_half_hole earclip failed ({len(ears)} tris)")
    top_ids = [m.tri(top[a], top[b], top[c]) for a, b, c in ears]
    bot_ids = [m.tri(bot[a], bot[c], bot[b]) for a, b, c in ears]
    # poly: 0=(0,0) 1=(w,0) 2=(w,depth) 3..3+n_arc = arc (n_arc+1 verts on the cylinder) then last=(0,depth)
    # Cylinder walls are the n_arc edges of the semicircle, not the (w,depth)→arc-start chord.
    arc_start = 3
    walls_cyl: List[int] = []
    walls_other: List[int] = []
    nv = len(poly)
    for i in range(nv):
        j = (i + 1) % nv
        quads = m.quad(bot[i], bot[j], top[j], top[i])
        # arc walls: edges from vertex arc_start .. arc_start+n_arc-1
        if arc_start <= i < arc_start + n_arc:
            walls_cyl.extend(quads)
        else:
            walls_other.extend(quads)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in walls_cyl:
        tr[t] = 0
    for t in top_ids:
        tr[t] = 1
    for t in bot_ids + walls_other:
        ti[t] = 0
    n_sides = max(4, round(TWOPI * n_arc / math.pi))
    loc = (cx, cy, 0.0)
    axis = (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[ti][k] for ti in walls_cyl for k in range(3)], loc, axis)
    dvol = cyl_dvol(rad, n_sides, thick, False, n_bands=n_arc)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=walls_cyl, loops=[], closed360=False, outward=False,
            radius=rad, umin=-math.pi, umax=0.0, vmin=0.0, vmax=thick,
            axis=ax3(loc, axis, xdir), n_sides=n_sides, chord=chord(rad, n_sides), dvol=dvol,
            exp={"builtAs": "single"},
        ),
        base_region(
            1, stype="plane", origin="planeGrow", tris=top_ids, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, thick), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[ti] for ti in top_ids], (0.0, 0.0, 1.0), (0.0, 0.0, thick)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("plate_half_hole", m, tr, ti, regs, {
        "faceCount": 2 + len(bot_ids) + len(walls_other),
        "surfaceCensus": {"cylinder": 1, "plane": 1, "planar_facets": len(bot_ids) + len(walls_other)},
        "rung": 1, "rounds": 0,
        "note": "Isolated ~180° partial depression on a plate; not S05 slot. One manifold solid.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F03 slotted_stadium — two 180° half-pipes + tangent flats, shared generators
# ---------------------------------------------------------------------------
def fixture_slotted_stadium() -> None:
    r, h, gap, n = 3.0, 8.0, 10.0, 4
    m = MeshBuilder()

    def half(center_x: float, a0: float, a1: float) -> Tuple[List[int], List[int]]:
        bot, top = [], []
        for i in range(n + 1):
            ang = a0 + (a1 - a0) * i / n
            bot.append(m.v(center_x + r * math.cos(ang), r * math.sin(ang), 0.0))
            top.append(m.v(center_x + r * math.cos(ang), r * math.sin(ang), h))
        return bot, top

    # left: π/2 → 3π/2 through -x (outer left). right: -π/2 → π/2 through +x.
    lb, lt = half(0.0, math.pi / 2, 3 * math.pi / 2)
    rb, rt = half(gap, -math.pi / 2, math.pi / 2)
    cyl_l: List[int] = []
    cyl_r: List[int] = []
    for i in range(n):
        cyl_l.extend(m.quad(lb[i], lb[i + 1], lt[i + 1], lt[i]))
        cyl_r.extend(m.quad(rb[i], rb[i + 1], rt[i + 1], rt[i]))
    # flats: y=+r connects left first to right last; y=-r connects left last to right first
    flat_top = m.quad(lb[0], rb[-1], rt[-1], lt[0])
    flat_bot = m.quad(lb[-1], lt[-1], rt[0], rb[0])
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyl_l:
        tr[t] = 0
    for t in cyl_r:
        tr[t] = 1
    for t in flat_bot:
        tr[t] = 2
    for t in flat_top:
        tr[t] = 3
    n_sides = round(TWOPI * n / math.pi)
    dvl = cyl_dvol(r, n_sides, h, False, n_bands=n)
    loc_l, loc_r = (0.0, 0.0, 0.0), (gap, 0.0, 0.0)
    axis = (0.0, 0.0, 1.0)
    xl = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl_l for k in range(3)], loc_l, axis)
    xr = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl_r for k in range(3)], loc_r, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=cyl_l, loops=[], closed360=False, outward=False,
            radius=r, umin=math.pi / 2, umax=3 * math.pi / 2, vmin=0.0, vmax=h,
            axis=ax3(loc_l, axis, xl), n_sides=n_sides, chord=chord(r, n_sides), dvol=dvl,
            exp={"builtAs": "single"},
        ),
        base_region(
            1, stype="cylinder", origin="cylGrow", tris=cyl_r, loops=[], closed360=False, outward=False,
            radius=r, umin=-math.pi / 2, umax=math.pi / 2, vmin=0.0, vmax=h,
            axis=ax3(loc_r, axis, xr), n_sides=n_sides, chord=chord(r, n_sides), dvol=dvl,
            exp={"builtAs": "single"},
        ),
        base_region(
            2, stype="plane", origin="planeGrow", tris=flat_bot, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, -r, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in flat_bot], (0.0, -1.0, 0.0), (0.0, -r, 0.0)),
            exp={"builtAs": "single"},
        ),
        base_region(
            3, stype="plane", origin="planeGrow", tris=flat_top, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, r, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in flat_top], (0.0, 1.0, 0.0), (0.0, r, 0.0)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("slotted_stadium", m, tr, ti, regs, {
        "faceCount": 4,
        "surfaceCensus": {"cylinder": 2, "plane": 2},
        "rung": 1, "rounds": 0,
        "note": "Two ~180° partials + connecting flats sharing generators; G1 tangent on cyl|flat.",
    }, n_islands=0)


# ---------------------------------------------------------------------------
# F04 quarter_round_1tri — D7 R=2, vertices on the cylinder, one G1 generator
# ---------------------------------------------------------------------------
def fixture_quarter_round_1tri() -> None:
    R, L, W = 2.0, 6.0, 6.0
    m = MeshBuilder()
    # Tangent walls occupy y>=R (plane x=0) and x>=R (plane y=0). Axis at (R,R).
    a0R, a0W = m.v(0, R, 0), m.v(0, R + W, 0)
    a0Rz, a0Wz = m.v(0, R, L), m.v(0, R + W, L)
    bR0, bW0 = m.v(R, 0, 0), m.v(R + W, 0, 0)
    bR0z, bW0z = m.v(R, 0, L), m.v(R + W, 0, L)
    plane_x = m.quad(a0R, a0W, a0Wz, a0Rz)
    plane_y = m.quad(bR0, bR0z, bW0z, bW0)
    # 1-tri fillet: generator on plane x=0. Island chord at z=0 ties in plane y=0.
    fil = [m.tri(a0R, a0Rz, bR0)]
    island = [m.tri(a0R, bR0, bW0)]
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in plane_x:
        tr[t] = 0
    for t in plane_y:
        tr[t] = 1
    for t in fil:
        tr[t] = 2
    for t in island:
        ti[t] = 0
    n_sides = 4
    dvol_f = fillet_dvol(R, L, False, math.pi / 2)
    loc = (R, R, 0.0)
    axis = (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[t][k] for t in fil for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="plane", origin="planeGrow", tris=plane_x, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in plane_x], (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
        base_region(
            1, stype="plane", origin="planeGrow", tris=plane_y, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in plane_y], (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
        base_region(
            2, stype="cylinder", origin="filletStrip", tris=fil, loops=[], closed360=False, outward=False,
            radius=R, umin=math.pi, umax=1.5 * math.pi, vmin=0.0, vmax=L,
            axis=ax3(loc, axis, xdir), n_sides=n_sides, chord=chord(R, n_sides), dvol=dvol_f,
            fil_a=0, fil_b=1, max_dev=0.0, exp={"builtAs": "single", "filletRadiusMM": R},
        ),
    ]
    emit("quarter_round_1tri", m, tr, ti, regs, {
        "faceCount": 3 + len(island),
        "surfaceCensus": {"plane": 2, "cylinder": 1, "planar_facets": len(island)},
        "rung": 1, "rounds": 0,
        "note": "D7 ground-truth R=2.0; all fillet verts on the cylinder. Withdrawn median rule → R/2.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F05 fillet_strip_2tri — one-band quad (2 tris), both long sides G1
# ---------------------------------------------------------------------------
def fixture_fillet_strip_2tri() -> None:
    R, L, W = 3.0, 8.0, 6.0
    m = MeshBuilder()
    a0R, a0W = m.v(0, R, 0), m.v(0, R + W, 0)
    a0Rz, a0Wz = m.v(0, R, L), m.v(0, R + W, L)
    bR0, bW0 = m.v(R, 0, 0), m.v(R + W, 0, 0)
    bR0z, bW0z = m.v(R, 0, L), m.v(R + W, 0, L)
    plane_x = m.quad(a0R, a0W, a0Wz, a0Rz)
    plane_y = m.quad(bR0, bR0z, bW0z, bW0)
    strip = m.quad(a0R, a0Rz, bR0z, bR0)  # 2 tris, both generators
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in plane_x:
        tr[t] = 0
    for t in plane_y:
        tr[t] = 1
    for t in strip:
        tr[t] = 2
    n_sides = 4
    dvol_f = fillet_dvol(R, L, False, math.pi / 2)
    loc = (R, R, 0.0)
    axis = (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[t][k] for t in strip for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="plane", origin="planeGrow", tris=plane_x, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in plane_x], (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
        base_region(
            1, stype="plane", origin="planeGrow", tris=plane_y, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in plane_y], (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
        base_region(
            2, stype="cylinder", origin="filletStrip", tris=strip, loops=[], closed360=False, outward=False,
            radius=R, umin=math.pi, umax=1.5 * math.pi, vmin=0.0, vmax=L,
            axis=ax3(loc, axis, xdir), n_sides=n_sides, chord=chord(R, n_sides), dvol=dvol_f,
            fil_a=0, fil_b=1, exp={"builtAs": "single"},
        ),
    ]
    emit("fillet_strip_2tri", m, tr, ti, regs, {
        "faceCount": 3,
        "surfaceCensus": {"plane": 2, "cylinder": 1},
        "rung": 1, "rounds": 0,
        "note": "2-tri (one-band) strip; both long chains tangent:true (constructed generators).",
    }, n_islands=0)


# ---------------------------------------------------------------------------
# F06 counterbore — two Seamed360 cylinders + annular plane with inner wire
# ---------------------------------------------------------------------------
def fixture_counterbore() -> None:
    n, r_o, r_i, h1, h2 = 6, 6.0, 3.0, 4.0, 6.0
    m = MeshBuilder()
    ob, ot = ring(m, r_o, 0.0, n), ring(m, r_o, h1, n)
    ib, it = ring(m, r_i, h1, n), ring(m, r_i, h1 + h2, n)
    # share ot with annulus outer, ib with annulus inner
    outer = stitch_cyl(m, ob, ot, outward=True)
    inner = stitch_cyl(m, ib, it, outward=True)
    cap = stitch_annulus(m, ot, ib, normal_up=True)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in outer:
        tr[t] = 0
    for t in inner:
        tr[t] = 1
    for t in cap:
        tr[t] = 2
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    xo = lowest_id_xdir(m.pts, [m.tris[t][k] for t in outer for k in range(3)], loc, axis)
    xi = lowest_id_xdir(m.pts, [m.tris[t][k] for t in inner for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=outer, loops=[], closed360=True, outward=False,
            radius=r_o, umin=0.0, umax=TWOPI, vmin=0.0, vmax=h1,
            axis=ax3(loc, axis, xo), n_sides=n, chord=chord(r_o, n),
            dvol=cyl_dvol(r_o, n, h1, False),
            exp={"builtAs": "seamed360"},
        ),
        base_region(
            1, stype="cylinder", origin="cylGrow", tris=inner, loops=[], closed360=True, outward=False,
            radius=r_i, umin=0.0, umax=TWOPI, vmin=h1, vmax=h1 + h2,
            axis=ax3(loc, axis, xi), n_sides=n, chord=chord(r_i, n),
            dvol=cyl_dvol(r_i, n, h2, False),
            exp={"builtAs": "seamed360"},
        ),
        base_region(
            2, stype="plane", origin="planeGrow", tris=cap, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, h1), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in cap], (0.0, 0.0, 1.0), (0.0, 0.0, h1)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("counterbore", m, tr, ti, regs, {
        "faceCount": 3,
        "surfaceCensus": {"cylinder": 2, "plane": 1},
        "rung": 1, "rounds": 0,
        "note": "Annular cap has outer (outer cyl CapHigh) and inner (inner cyl CapLow) wires.",
    }, n_islands=0)


# ---------------------------------------------------------------------------
# F07 S09 — analytic plane sharing an edge with a faceted tet island
# ---------------------------------------------------------------------------
def fixture_s09_mixed() -> None:
    m = MeshBuilder()
    a, b, c, d = m.v(0, 0, 0), m.v(10, 0, 0), m.v(10, 10, 0), m.v(0, 10, 0)
    plane = [m.tri(a, b, c), m.tri(a, c, d)]
    # Open tet-corner: three faces, missing base (the base IS the plane edge b-c).
    e, f = m.v(13, 4, 3), m.v(12, 7, 6)
    island = [m.tri(b, c, e), m.tri(b, e, f), m.tri(c, f, e)]
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in plane:
        tr[t] = 0
    for t in island:
        ti[t] = 0
    regs = [
        base_region(
            0, stype="plane", origin="planeGrow", tris=plane, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in plane], (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("s09_mixed", m, tr, ti, regs, {
        "faceCount": 1 + len(island),
        "surfaceCensus": {"plane": 1, "planar_facets": len(island)},
        "maxEdgeTolMM": 0.05,
        "rung": 1, "rounds": 0,
        "note": "Mixed analytic|faceted chain is a shared mesh edge; chords kept verbatim.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F08 R1-success — closed triangular prism; 3-sided closed360 explodes, shell stays closed
# ---------------------------------------------------------------------------
def fixture_r1_success() -> None:
    n, r, h = 3, 4.0, 5.0
    m = MeshBuilder()
    bot, top = ring(m, r, 0.0, n), ring(m, r, h, n)
    cyl = stitch_cyl(m, bot, top, outward=True)
    cap_t = cap_fan(m, top, reverse=False)
    cap_b = cap_fan(m, bot, reverse=True)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyl:
        tr[t] = 0
    for t in cap_t + cap_b:
        ti[t] = 0
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=cyl, loops=[], closed360=True, outward=False,
            radius=r, umin=0.0, umax=TWOPI, vmin=0.0, vmax=h,
            axis=ax3(loc, axis, xdir), n_sides=n, chord=chord(r, n),
            dvol=cyl_dvol(r, n, h, False),
            exp={"builtAs": "explodedToFacets", "rung": 3},
        ),
    ]
    emit("r1_success", m, tr, ti, regs, {
        "faceCount": len(cyl) + len(cap_t) + len(cap_b),
        "surfaceCensus": {"planar_facets": len(cyl) + len(cap_t) + len(cap_b)},
        "rung": 3, "rounds": 1,
        "note": "Closed triangular prism; 3-sided closed360 forces MakeFace/seam failure; R1 explode; shell stays closed.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F09 R1 round-2 cascade — two stacked n=3 cylinders + annulus, closed
# ---------------------------------------------------------------------------
def fixture_r1_round2_cascade() -> None:
    n, r0, r1, h = 3, 5.0, 2.5, 4.0
    m = MeshBuilder()
    b0, t0 = ring(m, r0, 0.0, n), ring(m, r0, h, n)
    b1, t1 = ring(m, r1, h, n), ring(m, r1, 2 * h, n)
    cyl0 = stitch_cyl(m, b0, t0, outward=True)
    cyl1 = stitch_cyl(m, b1, t1, outward=True)
    annulus = stitch_annulus(m, t0, b1, normal_up=True)
    cap_b = cap_fan(m, b0, reverse=True)
    cap_t = cap_fan(m, t1, reverse=False)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyl0:
        tr[t] = 0
    for t in cyl1:
        tr[t] = 1
    for t in annulus:
        tr[t] = 2
    for t in cap_b + cap_t:
        ti[t] = 0
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    x0 = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl0 for k in range(3)], loc, axis)
    x1 = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl1 for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=cyl0, loops=[], closed360=True, outward=False,
            radius=r0, umin=0.0, umax=TWOPI, vmin=0.0, vmax=h,
            axis=ax3(loc, axis, x0), n_sides=n, chord=chord(r0, n),
            dvol=cyl_dvol(r0, n, h, False),
            exp={"builtAs": "explodedToFacets"},
        ),
        base_region(
            1, stype="cylinder", origin="cylGrow", tris=cyl1, loops=[], closed360=True, outward=False,
            radius=r1, umin=0.0, umax=TWOPI, vmin=h, vmax=2 * h,
            axis=ax3(loc, axis, x1), n_sides=n, chord=chord(r1, n),
            dvol=cyl_dvol(r1, n, h, False),
            exp={"builtAs": "explodedToFacets"},
        ),
        base_region(
            2, stype="plane", origin="planeGrow", tris=annulus, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, h), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in annulus], (0.0, 0.0, 1.0), (0.0, 0.0, h)),
            exp={"builtAs": "single"},
        ),
    ]
    emit("r1_round2_cascade", m, tr, ti, regs, {
        "faceCount": len(cyl0) + len(cyl1) + 1 + len(cap_b) + len(cap_t),
        "surfaceCensus": {"planar_facets": len(cyl0) + len(cyl1) + len(cap_b) + len(cap_t), "plane": 1},
        "rung": 3, "rounds": 2,
        "note": "First explode forces annular cap rebuild; second explode in round 2. Distinct caps on both cylinders.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F10 R2-forced — closed square prism (NOT a clone of R1); J6/ChainUnstable
# ---------------------------------------------------------------------------
def fixture_r2_forced() -> None:
    # Distinct from full_360_hole (nSides=4, sagitta/R=1-cos(pi/4)≈0.293, Seamed360).
    # Three nSides=3 closed360 cylinders (sagitta/R=0.5) stacked with two annuli:
    # F3 residual gate fails rung 1 on every cylinder; each explode un-collapses a
    # neighbour cap, so a 2-round R1 cap cannot fixpoint → ChainUnstable / J6.
    n = 3
    radii = (5.0, 3.5, 2.0)
    h = 3.0
    m = MeshBuilder()
    rings_bot = [ring(m, radii[i], i * h, n) for i in range(3)]
    rings_top = [ring(m, radii[i], (i + 1) * h, n) for i in range(3)]
    cyls = [stitch_cyl(m, rings_bot[i], rings_top[i], outward=True) for i in range(3)]
    an0 = stitch_annulus(m, rings_top[0], rings_bot[1], normal_up=True)
    an1 = stitch_annulus(m, rings_top[1], rings_bot[2], normal_up=True)
    cap_b = cap_fan(m, rings_bot[0], reverse=True)
    cap_t = cap_fan(m, rings_top[2], reverse=False)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyls[0]:
        tr[t] = 0
    for t in cyls[1]:
        tr[t] = 1
    for t in cyls[2]:
        tr[t] = 2
    for t in an0:
        tr[t] = 3
    for t in an1:
        tr[t] = 4
    for t in cap_b + cap_t:
        ti[t] = 0
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    xs = [
        lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyls[i] for k in range(3)], loc, axis)
        for i in range(3)
    ]
    regs = []
    for i, (rad, tix, xdir) in enumerate(zip(radii, cyls, xs)):
        regs.append(base_region(
            i, stype="cylinder", origin="cylGrow", tris=tix, loops=[], closed360=True, outward=False,
            radius=rad, umin=0.0, umax=TWOPI, vmin=i * h, vmax=(i + 1) * h,
            axis=ax3(loc, axis, xdir), n_sides=n, chord=chord(rad, n),
            dvol=cyl_dvol(rad, n, h, False),
            exp={"builtAs": "explodedToFacets"},
        ))
    regs.append(base_region(
        3, stype="plane", origin="planeGrow", tris=an0, loops=[], closed360=False, outward=True,
        axis=ax3((0.0, 0.0, h), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        dvol=plane_dvol(m.pts, [m.tris[t] for t in an0], (0.0, 0.0, 1.0), (0.0, 0.0, h)),
        exp={"builtAs": "single"},
    ))
    regs.append(base_region(
        4, stype="plane", origin="planeGrow", tris=an1, loops=[], closed360=False, outward=True,
        axis=ax3((0.0, 0.0, 2 * h), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        dvol=plane_dvol(m.pts, [m.tris[t] for t in an1], (0.0, 0.0, 1.0), (0.0, 0.0, 2 * h)),
        exp={"builtAs": "single"},
    ))
    n_exploded = sum(len(c) for c in cyls) + len(cap_b) + len(cap_t)
    emit("r2_forced", m, tr, ti, regs, {
        "faceCount": n_exploded + 2,
        "surfaceCensus": {"planar_facets": n_exploded, "plane": 2},
        "rung": 3, "rounds": 1, "abandonComponent": False,
        "sagittaOverR": 1.0 - math.cos(math.pi / n),
        "nSides": n,
        "note": "Three stacked nSides=3 closed360 cylinders (sagitta/R=0.5) + 2 annuli. Distinct from full_360_hole (nSides=4, sagitta/R≈0.293, Seamed360). F3 residual fails every cylinder; 2-round R1 cap cannot fixpoint → ChainUnstable/J6.",
    }, n_islands=1)


# ---------------------------------------------------------------------------
# F11 R8 — n=3 closed washer; explode un-collapses cap rings + planar inner wires
# ---------------------------------------------------------------------------
def fixture_r8_360_explode_caps() -> None:
    n, r_in, r_out, h = 3, 5.0, 8.0, 6.0
    m = MeshBuilder()
    ib, it = ring(m, r_in, 0.0, n), ring(m, r_in, h, n)
    ob, ot = ring(m, r_out, 0.0, n), ring(m, r_out, h, n)
    cyl = stitch_cyl(m, ib, it, outward=True)
    top = stitch_annulus(m, ot, it, normal_up=True)
    bot = stitch_annulus(m, ob, ib, normal_up=False)
    outer = stitch_cyl(m, ob, ot, outward=False)
    n_tri = len(m.tris)
    tr = [-1] * n_tri
    ti = [-1] * n_tri
    for t in cyl:
        tr[t] = 0
    for t in top:
        tr[t] = 1
    for t in bot + outer:
        ti[t] = 0
    loc, axis = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    xdir = lowest_id_xdir(m.pts, [m.tris[t][k] for t in cyl for k in range(3)], loc, axis)
    regs = [
        base_region(
            0, stype="cylinder", origin="cylGrow", tris=cyl, loops=[], closed360=True, outward=False,
            radius=r_in, umin=0.0, umax=TWOPI, vmin=0.0, vmax=h,
            axis=ax3(loc, axis, xdir), n_sides=n, chord=chord(r_in, n),
            dvol=cyl_dvol(r_in, n, h, False),
            exp={"builtAs": "explodedToFacets", "capRingsUncollapsed": True},
        ),
        base_region(
            1, stype="plane", origin="planeGrow", tris=top, loops=[], closed360=False, outward=True,
            axis=ax3((0.0, 0.0, h), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            dvol=plane_dvol(m.pts, [m.tris[t] for t in top], (0.0, 0.0, 1.0), (0.0, 0.0, h)),
            exp={"builtAs": "single", "innerWireRebuiltFromChords": True},
        ),
    ]
    emit("r8_360_explode_caps", m, tr, ti, regs, {
        "faceCount": len(cyl) + 1 + len(bot) + len(outer),
        "surfaceCensus": {"planar_facets": len(cyl) + len(bot) + len(outer), "plane": 1},
        "rung": 3, "rounds": 1,
        "note": "P1 cap rings are one closedLoop chain of nSides mesh chords (the gp_Circ collapse). Explode un-collapses each ring to those N chords and rebuilds the planar inner from CapHigh.",
    }, n_islands=1)


def main() -> None:
    print("Generating 11 fixtures under", HERE)
    fixture_full_360_hole()
    fixture_plate_half_hole()
    fixture_slotted_stadium()
    fixture_quarter_round_1tri()
    fixture_fillet_strip_2tri()
    fixture_counterbore()
    fixture_s09_mixed()
    fixture_r1_success()
    fixture_r1_round2_cascade()
    fixture_r2_forced()
    fixture_r8_360_explode_caps()
    print("Generated 11 fixtures")


if __name__ == "__main__":
    main()
