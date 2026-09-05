#!/usr/bin/env python3
"""gate_130 — 1.3.0 battery (B1–B6). B6 is a recorded deferral (D-130-23).

For every corpus sidecar with `"battery": "130"`: convert with `--smooth
--no-verify`, assert ok / solids=1 / openShells=0 / reverted=0, distinct
cylindrical surfaces (D-130-20) == GT wall count, distinct conical surfaces
== GT cone count (DIAG_130_CENSUS ChamferCone is diagnostic only), volume
delta ≤ 0.01% (census B-Rep volume vs mesh), (exactVolume where the sidecar
records it, D-130-15(1)), and when the sidecar lists intersections, each is
represented by its D-130-2 tier:
tier 1 (`cylplane`, `coneplane`, `conecyl-coaxial`, `cylcyl-coaxial`) ships
a CIRCLE / ELLIPSE / LINE; tier 2 (`cylcyl`, the general skew quartic, and
`conecyl`) ships the MESH POLYLINE as ONE edge shared by the two analytic
faces, counted by the 130-BIND census under `edgeClasses.polylineTier2`,
with no `unhandled`, `overTol` or `overCap` on any analytic|analytic edge.

Volume cell (D-130-15(1)): where the sidecar carries `exactVolume` (the
generator's exact analytic volume, B2–B6) the STEP B-Rep volume is compared
against it; otherwise against the mesh volume as before. Threshold unchanged.

D-130-23: B6 (`cyl_meets_chamfer`) is listed in
`tests/gates/baseline/expected-red.json`. A listed FAIL is XFAIL (exit 0);
a listed row that turns green prints ``XPASS <fixture> — shrink
expected-red.json`` and does not fail; an unlisted FAIL exits 1. CMake
does not invert this test.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = REPO / "tests" / "corpus"
EXPECTED_RED_FILE = REPO / "tests" / "gates" / "baseline" / "expected-red.json"
ENV_STRICT = "GATE_130_STRICT"
MARKER_UNEXPECTED = "GATE_130_UNEXPECTED"
MARKER_FLOORS_MET = "GATE_130_FLOORS_MET"

BATTERY_ORDER = (
    "linkage_bores_chamfer",
    "cross_bores",
    "chamfer_straight_ring",
    "boss_cone_chamfer",
    "counterbore_chamfer",
    "cyl_meets_chamfer",
)
BATTERY_LABEL = {
    "linkage_bores_chamfer": "B1",
    "cross_bores": "B2",
    "chamfer_straight_ring": "B3",
    "boss_cone_chamfer": "B4",
    "counterbore_chamfer": "B5",
    "cyl_meets_chamfer": "B6",
}

DIAG_130_RE = re.compile(
    r"DIAG_130_CENSUS\b[^\n]*\bcone=(\d+)",
    re.MULTILINE,
)
POLYLINE_RE = re.compile(
    r"\bPOLYLINE\s*\(\s*(?:'[^']*'|\"[^\"]*\"|\$)?\s*,?\s*\(([^)]*)\)",
    re.IGNORECASE,
)
COMPOSITE_SEGS_RE = re.compile(
    r"\bCOMPOSITE_CURVE\s*\([^;]*?\(([^)]*)\)",
    re.IGNORECASE | re.DOTALL,
)


def strict_enabled(env: Optional[Mapping[str, str]] = None) -> bool:
    environ: Mapping[str, str] = os.environ if env is None else env
    raw = environ.get(ENV_STRICT)
    if raw is None or str(raw).strip() == "":
        return True
    return str(raw).strip().lower() not in ("0", "false", "off", "no")


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise RuntimeError("no RESULT line")


def gt_counts(sc: Dict[str, Any]) -> Tuple[int, int]:
    rec = sc.get("recoverable") or []
    cyl = sum(int(p.get("count") or 0) for p in rec if p.get("type") == "cylinder")
    cone = sum(int(p.get("count") or 0) for p in rec if p.get("type") == "cone")
    live = next((row for row in (sc.get("live") or []) if isinstance(row, dict)), {})
    census = live.get("surfaceCensus") or {}
    if cyl == 0:
        cyl = int(census.get("cylinder") or live.get("builtCylindersFloor") or 0)
    if cone == 0:
        cone = int(census.get("cone") or live.get("builtConesFloor") or 0)
    return cyl, cone


def discover_battery(corpus: Path) -> List[Tuple[str, Path, Path, Dict[str, Any]]]:
    found: Dict[str, Tuple[Path, Path, Dict[str, Any]]] = {}
    for sidecar in sorted(corpus.glob("*.expected.json")):
        try:
            doc = json.loads(sidecar.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if str(doc.get("battery") or "") != "130":
            continue
        fid = sidecar.name[: -len(".expected.json")]
        stl = corpus / f"{fid}.stl"
        found[fid] = (stl, sidecar, doc)
    ordered: List[Tuple[str, Path, Path, Dict[str, Any]]] = []
    seen = set()
    for fid in BATTERY_ORDER:
        if fid in found:
            stl, scp, doc = found[fid]
            ordered.append((fid, stl, scp, doc))
            seen.add(fid)
    for fid in sorted(found):
        if fid in seen:
            continue
        stl, scp, doc = found[fid]
        ordered.append((fid, stl, scp, doc))
    return ordered


TIER1_KINDS = ("cylplane", "coneplane", "conecyl-coaxial", "cylcyl-coaxial")
TIER2_KINDS = ("cylcyl", "conecyl")

ENTITY_RE = re.compile(r"^#(\d+)\s*=\s*(.*?);\s*$", re.MULTILINE | re.DOTALL)
REF_RE = re.compile(r"#(\d+)")
NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")
# D-130-20: axis line and placement within 1e-5 mm; radius within 1e-5 mm on
# the fitted value (adjacent face pieces of one wall can differ by ~1.1e-5 on
# B1 slot ends — still one surface).
SURFACE_Q_MM = 1e-5
SURFACE_RAD_TOL_MM = 1.2e-5


def _step_entities(step_text: str) -> Dict[int, str]:
    """#id -> entity body (name + args), one entry per STEP instance."""
    out: Dict[int, str] = {}
    for m in ENTITY_RE.finditer(step_text):
        out[int(m.group(1))] = " ".join(m.group(2).split())
    return out


def _entity_name(body: str) -> str:
    return body.split("(", 1)[0].strip().upper()


def _split_step_args(args: str) -> List[str]:
    out: List[str] = []
    depth = 0
    in_str = False
    start = 0
    for i, c in enumerate(args):
        if in_str:
            if c == "'" and i + 1 < len(args) and args[i + 1] == "'":
                continue
            if c == "'":
                in_str = False
            continue
        if c == "'":
            in_str = True
            continue
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "," and depth == 0:
            out.append(args[start:i].strip())
            start = i + 1
    tail = args[start:].strip()
    if tail:
        out.append(tail)
    return out


def _parse_cartesian_point(body: str) -> Optional[Tuple[float, float, float]]:
    m = re.search(r"CARTESIAN_POINT\s*\([^,]*,\(([^)]*)\)", body)
    if not m:
        return None
    vals = [float(x.strip()) for x in m.group(1).split(",") if x.strip()]
    if len(vals) != 3:
        return None
    return vals[0], vals[1], vals[2]


def _parse_direction(body: str) -> Optional[Tuple[float, float, float]]:
    m = re.search(r"DIRECTION\s*\([^,]*,\(([^)]*)\)", body)
    if not m:
        return None
    vals = [float(x.strip()) for x in m.group(1).split(",") if x.strip()]
    if len(vals) != 3:
        return None
    return vals[0], vals[1], vals[2]


def _parse_axis2_placement(
    ents: Dict[int, str], ref: int
) -> Tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    body = ents.get(ref, "")
    refs = [int(x) for x in REF_RE.findall(body)]
    if len(refs) < 2:
        return None, None
    loc = _parse_cartesian_point(ents.get(refs[0], ""))
    direction = _parse_direction(ents.get(refs[1], ""))
    return loc, direction


def _canonical_direction(d: Tuple[float, float, float]) -> Tuple[float, float, float]:
    x, y, z = d
    n = math.sqrt(x * x + y * y + z * z)
    if n > 0:
        x, y, z = x / n, y / n, z / n
    ax, ay, az = abs(x), abs(y), abs(z)
    if az >= ax and az >= ay:
        if z < 0:
            x, y, z = -x, -y, -z
    elif ay >= ax:
        if y < 0:
            x, y, z = -x, -y, -z
    else:
        if x < 0:
            x, y, z = -x, -y, -z
    return x, y, z


def _axis_line_distance(
    a_loc: Tuple[float, float, float],
    a_dir: Tuple[float, float, float],
    b_loc: Tuple[float, float, float],
    b_dir: Tuple[float, float, float],
) -> float:
    dx = a_loc[0] - b_loc[0]
    dy = a_loc[1] - b_loc[1]
    dz = a_loc[2] - b_loc[2]
    cx = a_dir[1] * b_dir[2] - a_dir[2] * b_dir[1]
    cy = a_dir[2] * b_dir[0] - a_dir[0] * b_dir[2]
    cz = a_dir[0] * b_dir[1] - a_dir[1] * b_dir[0]
    cn = math.sqrt(cx * cx + cy * cy + cz * cz)
    if cn < 1e-12:
        fx = b_loc[0] - a_loc[0]
        fy = b_loc[1] - a_loc[1]
        fz = b_loc[2] - a_loc[2]
        px = fy * a_dir[2] - fz * a_dir[1]
        py = fz * a_dir[0] - fx * a_dir[2]
        pz = fx * a_dir[1] - fy * a_dir[0]
        return math.sqrt(px * px + py * py + pz * pz)
    return abs(dx * cx + dy * cy + dz * cz) / cn


def _dirs_parallel(
    a_dir: Tuple[float, float, float], b_dir: Tuple[float, float, float]
) -> bool:
    dot = a_dir[0] * b_dir[0] + a_dir[1] * b_dir[1] + a_dir[2] * b_dir[2]
    return abs(abs(dot) - 1.0) < 1e-8


def _along_axis_separation(
    a_loc: Tuple[float, float, float],
    b_loc: Tuple[float, float, float],
    direction: Tuple[float, float, float],
) -> float:
    vx = b_loc[0] - a_loc[0]
    vy = b_loc[1] - a_loc[1]
    vz = b_loc[2] - a_loc[2]
    t = vx * direction[0] + vy * direction[1] + vz * direction[2]
    return abs(t)


def _quantize_mm(v: float, q: float = SURFACE_Q_MM) -> float:
    return round(v / q) * q


def _radii_equal(r1: float, r2: float) -> bool:
    return abs(r1 - r2) <= SURFACE_RAD_TOL_MM


def _parse_surface_entity(
    ents: Dict[int, str], eid: int, surface_type: str
) -> Optional[Tuple[float, Tuple[float, float, float], Tuple[float, float, float], float]]:
    body = ents.get(eid, "")
    if _entity_name(body) != surface_type:
        return None
    open_paren = body.find("(")
    if open_paren < 0:
        return None
    args = _split_step_args(body[open_paren + 1 : body.rfind(")")])
    if len(args) < 3:
        return None
    pref = REF_RE.search(args[1])
    if not pref:
        return None
    loc, direction = _parse_axis2_placement(ents, int(pref.group(1)))
    if loc is None or direction is None:
        return None
    radius = float(NUM_RE.search(args[2]).group(0)) if NUM_RE.search(args[2]) else None
    if radius is None:
        return None
    semi = 0.0
    if surface_type == "CONICAL_SURFACE":
        if len(args) < 4:
            return None
        semi_m = NUM_RE.search(args[3])
        if not semi_m:
            return None
        semi = float(semi_m.group(0))
    return radius, loc, direction, semi


def count_distinct_cylinder_surfaces(step_text: str) -> Tuple[int, int]:
    """Return (distinct cylindrical surfaces, CYLINDRICAL_SURFACE entity count).

    D-130-20 / D-130-17: group CYLINDRICAL_SURFACE entities that share an axis
    line (parallel directions, axis-line distance <= q) and radius (within q) and
    are close along the axis (<= max(r)) — multi-face pieces of one wall. Do
    not merge separate walls that happen to be coaxial (B2 R5 halves, two-plate
    holes).
    """
    ents = _step_entities(step_text)
    surfaces: List[Tuple[float, Tuple[float, float, float], Tuple[float, float, float]]] = []
    for eid, body in ents.items():
        if _entity_name(body) != "CYLINDRICAL_SURFACE":
            continue
        parsed = _parse_surface_entity(ents, eid, "CYLINDRICAL_SURFACE")
        if parsed is None:
            continue
        radius, loc, direction, _semi = parsed
        surfaces.append((radius, loc, _canonical_direction(direction)))

    if not surfaces:
        return 0, 0

    parent = list(range(len(surfaces)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for i in range(len(surfaces)):
        r1, loc1, dir1 = surfaces[i]
        for j in range(i + 1, len(surfaces)):
            r2, loc2, dir2 = surfaces[j]
            if not _radii_equal(r1, r2):
                continue
            if not _dirs_parallel(dir1, dir2):
                continue
            if _axis_line_distance(loc1, dir1, loc2, dir2) > SURFACE_Q_MM:
                continue
            if _along_axis_separation(loc1, loc2, dir1) > max(r1, r2):
                continue
            union(i, j)

    groups = {find(i) for i in range(len(surfaces))}
    return len(groups), len(surfaces)


def count_distinct_cone_surfaces(step_text: str) -> Tuple[int, int]:
    """Return (distinct conical surfaces, CONICAL_SURFACE entity count).

    Group by axis placement (position + direction, direction sign-insensitive),
    reference radius, and semi-angle — each within q. Unlike cylinders, cones on
    the same axis line but at different placements (B1 chamfer frustums) stay
    distinct.
    """
    ents = _step_entities(step_text)
    keys: set = set()
    entities = 0
    for eid, body in ents.items():
        if _entity_name(body) != "CONICAL_SURFACE":
            continue
        parsed = _parse_surface_entity(ents, eid, "CONICAL_SURFACE")
        if parsed is None:
            continue
        entities += 1
        radius, loc, direction, semi = parsed
        cd = _canonical_direction(direction)
        key = (
            _quantize_mm(loc[0]),
            _quantize_mm(loc[1]),
            _quantize_mm(loc[2]),
            _quantize_mm(cd[0]),
            _quantize_mm(cd[1]),
            _quantize_mm(cd[2]),
            _quantize_mm(radius),
            _quantize_mm(semi),
        )
        keys.add(key)
    return len(keys), entities


def shared_tier2_edges(step_text: str) -> Tuple[List[str], List[str]]:
    """D-130-2 tier 2, structurally on the STEP file.

    Returns (shared, unshared): the EDGE_CURVEs whose underlying 3D curve is
    a degree-1 B_SPLINE_CURVE_WITH_KNOTS (the mesh polyline) and that are
    referenced by the loops of exactly two faces whose surfaces are both
    CYLINDRICAL_SURFACE (or CONICAL_SURFACE), versus such polyline edges that
    only ONE analytic face references. D-130-19: sharing is the EDGE_CURVE
    entity id on both faces' loops (EDGE_LOOP → ORIENTED_EDGE → EDGE_CURVE),
    whether that EDGE_CURVE's curve is the B-spline itself or a SURFACE_CURVE
    / SEAM_CURVE / INTERSECTION_CURVE pcurve pair whose 3D curve is the
    B-spline. Degree 1 is required of the underlying 3D curve, not of the
    wrapper. A tier-2 seam is honest only when it is one TShape on both
    faces; a polyline each face carries alone is the shell opening.
    """
    ents = _step_entities(step_text)
    name = {k: _entity_name(v) for k, v in ents.items()}

    # SURFACE_CURVE('', #3d, (#pcurve, ...), .PCURVE_S1.) — first ref is 3D.
    curve_wrappers = ("SURFACE_CURVE", "SEAM_CURVE", "INTERSECTION_CURVE")

    def underlying_3d(curve_id: int) -> int:
        seen: set = set()
        cid = curve_id
        while cid not in seen:
            seen.add(cid)
            if name.get(cid, "") not in curve_wrappers:
                return cid
            refs = [int(x) for x in REF_RE.findall(ents.get(cid, ""))]
            if not refs:
                return cid
            cid = refs[0]
        return cid

    # face -> surface entity, face -> set of edge_curve ids
    face_surface: Dict[int, int] = {}
    face_edges: Dict[int, set] = {}
    for fid, body in ents.items():
        if name[fid] != "ADVANCED_FACE":
            continue
        refs = [int(x) for x in REF_RE.findall(body)]
        if not refs:
            continue
        face_surface[fid] = refs[-1]  # ADVANCED_FACE('', (bounds), surface, sense)
        edges = set()
        for b in refs[:-1]:
            bb = ents.get(b, "")
            if name.get(b, "") not in ("FACE_BOUND", "FACE_OUTER_BOUND"):
                continue
            for loop in (int(x) for x in REF_RE.findall(bb)):
                lb = ents.get(loop, "")
                if name.get(loop, "") != "EDGE_LOOP":
                    continue
                for oe in (int(x) for x in REF_RE.findall(lb)):
                    ob = ents.get(oe, "")
                    if name.get(oe, "") != "ORIENTED_EDGE":
                        continue
                    ers = [int(x) for x in REF_RE.findall(ob)]
                    if ers:
                        edges.add(ers[-1])
        face_edges[fid] = edges
    analytic_surf = ("CYLINDRICAL_SURFACE", "CONICAL_SURFACE")
    edge_faces: Dict[int, List[int]] = {}
    for fid, edges in face_edges.items():
        if name.get(face_surface.get(fid, -1), "") not in analytic_surf:
            continue
        for e in edges:
            edge_faces.setdefault(e, []).append(fid)
    shared: List[str] = []
    unshared: List[str] = []
    for e, faces in edge_faces.items():
        eb = ents.get(e, "")
        if name.get(e, "") != "EDGE_CURVE":
            continue
        crefs = [int(x) for x in REF_RE.findall(eb)]
        if len(crefs) < 3:
            continue
        curve = underlying_3d(crefs[2])
        cb = ents.get(curve, "")
        if name.get(curve, "") != "B_SPLINE_CURVE_WITH_KNOTS":
            continue
        m = re.search(r"B_SPLINE_CURVE_WITH_KNOTS\s*\(\s*'[^']*'\s*,\s*(\d+)", cb)
        if not m or int(m.group(1)) != 1:
            continue
        nseg = len(REF_RE.findall(cb.split(",", 2)[2].split(")")[0])) - 1
        tag = f"#{e} polyline segments={nseg} faces={sorted(faces)}"
        if len(faces) == 2:
            shared.append(tag)
        else:
            unshared.append(tag)
    return shared, unshared


def tier1_curves_present(census: Optional[Dict[str, Any]]) -> bool:
    curves = (census or {}).get("curves") or {}
    return int(curves.get("circle", 0)) + int(curves.get("ellipse", 0)) > 0


def run_convert(binary: Path, stl: Path, step: Path) -> Tuple[int, Dict[str, Any], str, str]:
    env = os.environ.copy()
    env["STL2STEP_DIAG_130"] = "1"
    # Spec fallback: RESULT has no cone field. DIAG_130_CENSUS is the ChamferCone
    # census the engine already prints; P2_DIAG is set as the documented backup.
    env["STL2STEP_P2_DIAG"] = "1"
    proc = subprocess.run(
        [str(binary), str(stl), str(step), "--smooth", "--no-verify", "--quiet"],
        capture_output=True,
        text=True,
        env=env,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    try:
        result = parse_result_line(proc.stdout or "")
    except (RuntimeError, json.JSONDecodeError) as exc:
        result = {"ok": False, "error": f"RESULT parse: {exc}"}
    return proc.returncode, result, proc.stdout or "", combined


def run_census(census_bin: Path, step: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [str(census_bin), str(step)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout or "census failed")
    return json.loads(proc.stdout)


def evaluate_one(
    fid: str,
    stl: Path,
    sc: Dict[str, Any],
    binary: Path,
    census_bin: Optional[Path],
    work: Path,
) -> Dict[str, Any]:
    gt_cyl, gt_cone = gt_counts(sc)
    row: Dict[str, Any] = {
        "id": fid,
        "label": BATTERY_LABEL.get(fid, fid),
        "gt_cyl": gt_cyl,
        "gt_cone": gt_cone,
        "ok": False,
        "solids": -1,
        "openShells": -1,
        "reverted": -1,
        "builtCyl": -1,
        "cylSurfaces": -1,
        "cylFaces": -1,
        "coneSurfaces": -1,
        "diagCone": -1,
        "stepCone": -1,
        "volPct": -1.0,
        "volRef": "none",
        "curves": "",
        "fails": [],
        "infra": [],
    }
    if not stl.is_file():
        row["infra"].append(f"missing STL {stl}")
        return row

    step = work / f"{fid}.step"
    rc, result, _stdout, combined = run_convert(binary, stl, step)
    row["ok"] = bool(result.get("ok"))
    row["solids"] = int(result.get("solids") or 0)
    row["openShells"] = int(result.get("openShells") or 0)
    row["reverted"] = int(result.get("smoothRevertedComponents") or 0)
    row["builtCyl"] = int(result.get("smoothBuiltCylinders") or 0)

    m = DIAG_130_RE.search(combined)
    if m:
        row["diagCone"] = int(m.group(1))
    row["coneSource"] = "DIAG_130_CENSUS Origin::ChamferCone (RESULT has no cone field)"

    census: Optional[Dict[str, Any]] = None
    if census_bin is not None and step.is_file():
        try:
            census = run_census(census_bin, step)
        except (RuntimeError, json.JSONDecodeError) as exc:
            row["infra"].append(f"census: {exc}")

    if census:
        row["stepCone"] = int((census.get("surfaces") or {}).get("cone") or 0)
        row["cylFaces"] = int((census.get("surfaces") or {}).get("cylinder") or 0)
        curves = census.get("curves") or {}
        row["curves"] = (
            f"LINE={curves.get('line', 0)} CIRCLE={curves.get('circle', 0)} "
            f"ELLIPSE={curves.get('ellipse', 0)} B_SPLINE={curves.get('bspline', 0)} "
            f"other={curves.get('other', 0)}"
        )
        mesh_vol = float(result.get("meshVolumeMM3") or sc.get("meshVolume") or 0.0)
        step_vol = float(census.get("volume") or 0.0)
        # D-130-15(1): where the sidecar records the generator's exact analytic
        # volume, the cell compares against THAT. An inscribed N-gon mesh is short
        # of the analytic solid by the sagitta volume by construction, so a
        # mesh-referenced cell is unreachable by a correct recovery (B3: 0.029 %,
        # B4: 0.044 %, B5: 0.015 % measured) and passes only a faceted one.
        exact_vol = float(sc.get("exactVolume") or 0.0)
        if exact_vol > 0.0 and step_vol != 0.0:
            row["volPct"] = abs(step_vol - exact_vol) / abs(exact_vol) * 100.0
            row["volRef"] = "exact"
        elif mesh_vol > 0.0 and step_vol != 0.0:
            row["volPct"] = abs(step_vol - mesh_vol) / abs(mesh_vol) * 100.0
            row["volRef"] = "mesh"
        elif float(result.get("volumeDeltaPct") or -1) >= 0:
            row["volPct"] = float(result["volumeDeltaPct"])
            row["volRef"] = "mesh"
    else:
        vd = float(result.get("volumeDeltaPct") or -1)
        row["volPct"] = vd

    step_text = ""
    if step.is_file():
        try:
            step_text = step.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            row["infra"].append(f"STEP read: {exc}")
    if step_text:
        cyl_distinct, cyl_entities = count_distinct_cylinder_surfaces(step_text)
        row["cylSurfaces"] = cyl_distinct
        if row["cylFaces"] < 0:
            row["cylFaces"] = cyl_entities
        cone_distinct, _cone_entities = count_distinct_cone_surfaces(step_text)
        row["coneSurfaces"] = cone_distinct

    cone_got = row["coneSurfaces"] if row["coneSurfaces"] >= 0 else (
        row["diagCone"] if row["diagCone"] >= 0 else row["stepCone"]
    )

    if not row["ok"]:
        row["fails"].append(f"ok=false {result.get('error', '')}".strip())
    if row["solids"] != 1:
        row["fails"].append(f"solids={row['solids']} != 1")
    if row["openShells"] != 0:
        row["fails"].append(f"openShells={row['openShells']} != 0")
    if row["reverted"] != 0:
        row["fails"].append(f"smoothRevertedComponents={row['reverted']} != 0")
    if row["cylSurfaces"] < 0:
        row["fails"].append("cylSurfaces unmeasured (STEP read/parse failed)")
    elif row["cylSurfaces"] != gt_cyl:
        faces = row["cylFaces"] if row["cylFaces"] >= 0 else row["builtCyl"]
        row["fails"].append(
            f"cylSurfaces={row['cylSurfaces']}/{gt_cyl} faces={faces} "
            f"(smoothBuiltCylinders={row['builtCyl']})"
        )
    if cone_got < 0:
        row["fails"].append("cones unmeasured (STEP read/parse failed)")
    elif cone_got != gt_cone:
        row["fails"].append(
            f"cones={cone_got}/{gt_cone} ({row.get('coneSource', 'STEP CONICAL_SURFACE')})"
        )
    if row["volPct"] < 0 or row["volPct"] > 0.01 + 1e-12:
        row["fails"].append(f"volumeDelta={row['volPct']}% > 0.01% (or unmeasured)")

    intersections = sc.get("intersections") or []
    if intersections and step_text:
        text = step_text
        if text:
            # D-130-2 tier rule (SPEC-130-cylcyl addendum 2026-09-03 23:40).
            kinds = [str(x.get("kind") or "") for x in intersections]
            shared, unshared = shared_tier2_edges(text)
            ec = result.get("edgeClasses") if isinstance(result.get("edgeClasses"), dict) else None
            if any(k in TIER2_KINDS for k in kinds):
                if not shared:
                    row["fails"].append(
                        "tier-2 intersection: no polyline edge shared by two analytic faces"
                        + (f" (unshared: {', '.join(unshared)})" if unshared else "")
                    )
                if unshared:
                    row["fails"].append(
                        "tier-2 polyline carried by ONE analytic face only: "
                        + ", ".join(unshared)
                    )
                if ec is None:
                    row["fails"].append(
                        "edgeClasses census absent from RESULT (130-BIND, D-130-2)"
                    )
                else:
                    if int(ec.get("polylineTier2") or 0) < 1:
                        row["fails"].append(
                            f"edgeClasses.polylineTier2={ec.get('polylineTier2')} < 1"
                        )
                    for red in ("unhandled", "overTol", "overCap"):
                        if int(ec.get(red) or 0) != 0:
                            row["fails"].append(f"edgeClasses.{red}={ec.get(red)} != 0")
            if any(k in TIER1_KINDS for k in kinds) and not tier1_curves_present(census):
                row["fails"].append(
                    "tier-1 intersection but no CIRCLE/ELLIPSE in the STEP census"
                )
            if not row["curves"]:
                row["curves"] = "STEP curves not censused"
            kinds_s = ", ".join(
                f"{x.get('a')}∩{x.get('b')}={x.get('kind')}" for x in intersections
            )
            tier2 = f" tier2 shared={len(shared)} unshared={len(unshared)}"
            row["curves"] = f"{row['curves']}  intersections[{kinds_s}]{tier2}"

    return row


def format_table(rows: List[Dict[str, Any]]) -> str:
    headers = (
        "id",
        "ok",
        "solids",
        "open",
        "reverted",
        "cylSurfaces/GT faces",
        "cones got/GT",
        "vol%",
        "curves / fails",
    )
    lines = ["gate_130 battery (D-130-23 recorded deferrals are XFAIL):", "  " + "  ".join(headers)]
    for r in rows:
        cone_got = (
            r["coneSurfaces"]
            if r["coneSurfaces"] >= 0
            else (r["diagCone"] if r["diagCone"] >= 0 else r["stepCone"])
        )
        cyl_faces = r["cylFaces"] if r["cylFaces"] >= 0 else r["builtCyl"]
        cyl_s = r["cylSurfaces"] if r["cylSurfaces"] >= 0 else r["builtCyl"]
        status = "FAIL" if (r["fails"] or r["infra"]) else "PASS"
        tail_bits = []
        if r["curves"]:
            tail_bits.append(r["curves"])
        if r["infra"] or r["fails"]:
            tail_bits.append("; ".join(r["infra"] + r["fails"]))
        tail = " | ".join(tail_bits) or "ok"
        lines.append(
            f"  {r['label']:<3} {r['id']:<24} {status:<4} "
            f"ok={int(r['ok'])} sol={r['solids']} open={r['openShells']} "
            f"rev={r['reverted']} "
            f"cylSurfaces={cyl_s}/{r['gt_cyl']} faces={cyl_faces} "
            f"cones={cone_got}/{r['gt_cone']} vol={r['volPct']:.4g}({r['volRef']})  {tail}"
        )
    return "\n".join(lines)


def current_blob(rows: List[Dict[str, Any]]) -> str:
    parts = []
    by_id = {r["id"]: r for r in rows}
    for fid, label in BATTERY_LABEL.items():
        r = by_id.get(fid)
        if not r:
            parts.append(f"{label}=?/?/?")
            continue
        cone_got = (
            r["coneSurfaces"]
            if r["coneSurfaces"] >= 0
            else (r["diagCone"] if r["diagCone"] >= 0 else r["stepCone"])
        )
        cyl_s = r["cylSurfaces"] if r["cylSurfaces"] >= 0 else r["builtCyl"]
        parts.append(f"{label}={cyl_s}/{cone_got}/{r['reverted']}")
    return " ".join(parts)


def load_expected_red_section(gate: str, path: Path = EXPECTED_RED_FILE) -> Dict[str, str]:
    """D-130-23 recorded-deferral list for one gate."""
    if not path.is_file():
        raise RuntimeError(f"required file missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(doc, dict):
        raise RuntimeError(f"{path}: expected a JSON object")
    section = doc.get(gate) or {}
    if not isinstance(section, dict):
        raise RuntimeError(f"{path}: {gate} must be an object")
    return {str(k): str(v) for k, v in section.items()}


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    p.add_argument("--census", type=Path, help="stl2step_census (volume + STEP surfaces/curves)")
    args = p.parse_args(argv)

    is_strict = strict_enabled()
    binary = args.binary
    if not binary.is_file():
        print(f"{MARKER_UNEXPECTED}  missing binary {binary}", file=sys.stderr)
        return 1

    fixtures = discover_battery(args.corpus)
    if not fixtures:
        print(
            f"{MARKER_UNEXPECTED}  no battery=130 sidecars under {args.corpus}",
            file=sys.stderr,
        )
        return 1

    census_bin = args.census if args.census and args.census.is_file() else None
    rows: List[Dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="gate_130_") as tmp:
        work = Path(tmp)
        for fid, stl, _scp, sc in fixtures:
            rows.append(evaluate_one(fid, stl, sc, binary, census_bin, work))

    print(
        f"gate_130  corpus={args.corpus}  binary={binary}  "
        f"census={census_bin or 'none'}  strict={is_strict}  "
        f"cone census=DIAG_130_CENSUS (RESULT has no cone field; "
        f"STL2STEP_P2_DIAG=1 also set)",
        flush=True,
    )
    print(format_table(rows), flush=True)
    print(f"current: {current_blob(rows)}", flush=True)

    infra = [r for r in rows if r["infra"]]
    battery_fail = [r for r in rows if r["fails"] and not r["infra"]]
    if infra:
        named = ", ".join(r["id"] for r in infra)
        print(f"{MARKER_UNEXPECTED}  infra fail(s): {named}", file=sys.stderr, flush=True)
        return 1

    try:
        listed = load_expected_red_section("gate_130")
    except (RuntimeError, json.JSONDecodeError, OSError) as exc:
        print(f"{MARKER_UNEXPECTED}  expected-red.json: {exc}", file=sys.stderr, flush=True)
        return 1

    unlisted_fail = [r for r in battery_fail if r["id"] not in listed]
    listed_fail = [r for r in battery_fail if r["id"] in listed]
    by_id = {r["id"]: r for r in rows}
    if unlisted_fail:
        named = ", ".join(f"{r['label']}:{r['id']}" for r in unlisted_fail)
        print(f"{MARKER_UNEXPECTED}  unmet battery: {named}", file=sys.stderr, flush=True)
        print(f"gate_130 FAIL  unlisted red: {named}", file=sys.stderr, flush=True)
        return 1
    for r in listed_fail:
        print(f"XFAIL {r['id']} — {listed[r['id']]}", flush=True)
    for fid in listed:
        r = by_id.get(fid)
        if r is not None and not r["fails"] and not r["infra"]:
            print(f"XPASS {fid} — shrink expected-red.json", flush=True)
    if not battery_fail:
        print(f"{MARKER_FLOORS_MET}  battery B1–B6 MET", flush=True)
        print("gate_130 PASS", flush=True)
        return 0
    print("gate_130 PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
