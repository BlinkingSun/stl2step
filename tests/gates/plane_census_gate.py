#!/usr/bin/env python3
"""plane_census_gate — per-GT-plane STEP face-count ratchet (D-FL-6).

Shaped on hl_census_ratchet.py. Floor authority is
tests/gates/baseline/plane-census-ratchet.json. Never expected-red.

The gate converts tests/corpus/handle-pickup.stl with --engine trueform,
lists every shipped STEP PLANE face (equation + bound-polygon area +
centroid + loop count), and reports per GT plane: nFaces (each shipped
face maps to exactly one GT face — D-S3-143), total area, hasLiveRegion.

Assignment (D-S3-142/146): equation match first; among co-equation GT
faces, nearest-match by gtFaceIdentity (nLoops, then area/centroid). A
tie is a named failure, never a coin flip. Co-equation ids are never
double-counted.

Cells:
  (a) per-plane nFaces <= ratcheted ceiling
  (b) phantom-plane count (unmatched PLANE, area > phantomAreaFloorMM2)
      <= ratcheted ceiling
  (c) unique faces matching live-region GT equations <= sumNFacesCeiling

Anti-confusion (mandatory): print RESULT smoothBuiltPlanes (regions) beside
the per-plane face counts. Regions and faces move in opposite directions
under consolidation.

Usage:
  plane_census_gate.py --self-test
  plane_census_gate.py --binary ./build/stl2step
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
DEFAULT_STL = REPO / "tests" / "corpus" / "handle-pickup.stl"
DEFAULT_GT = REPO / "tests" / "diag" / "handle-pickup" / "ground-truth.json"
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "plane-census-ratchet.json"
RATCHET_NAME = "plane-census-ratchet.json"
ANTI_CONFUSION = (
    "regions and faces move in opposite directions under consolidation; "
    "do not read a falling region count as regression"
)

REF_RE = re.compile(r"#(\d+)\b")
NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")


class GateError(Exception):
    """Hard gate failure."""


def load_json(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise GateError("no RESULT line in stl2step stdout")


def vdot(a: Sequence[float], b: Sequence[float]) -> float:
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def vnorm(u: Sequence[float]) -> Tuple[float, float, float]:
    L = math.sqrt(vdot(u, u)) or 1.0
    return (float(u[0]) / L, float(u[1]) / L, float(u[2]) / L)


def split_args(args: str) -> List[str]:
    out: List[str] = []
    depth = 0
    in_str = False
    start = 0
    i = 0
    while i < len(args):
        c = args[i]
        if in_str:
            if c == "'" and i + 1 < len(args) and args[i + 1] == "'":
                i += 2
                continue
            if c == "'":
                in_str = False
            i += 1
            continue
        if c == "'":
            in_str = True
            i += 1
            continue
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "," and depth == 0:
            out.append(args[start:i].strip())
            start = i + 1
        i += 1
    tail = args[start:].strip()
    if tail:
        out.append(tail)
    return out


def first_ref(arg: Optional[str]) -> Optional[int]:
    if not arg:
        return None
    m = REF_RE.search(arg)
    return int(m.group(1)) if m else None


def all_refs(arg: Optional[str]) -> List[int]:
    if not arg:
        return []
    return [int(x) for x in REF_RE.findall(arg)]


def parse_entities(text: str) -> Dict[int, Tuple[str, str]]:
    entities: Dict[int, Tuple[str, str]] = {}
    i, n = 0, len(text)
    while i < n:
        if text[i] != "#":
            i += 1
            continue
        j = i + 1
        while j < n and text[j].isdigit():
            j += 1
        if j == i + 1:
            i += 1
            continue
        eid = int(text[i + 1 : j])
        while j < n and text[j] in " \t\r\n":
            j += 1
        if j >= n or text[j] != "=":
            i = j
            continue
        j += 1
        while j < n and text[j] in " \t\r\n":
            j += 1
        start = j
        depth = 0
        in_str = False
        while j < n:
            c = text[j]
            if in_str:
                if c == "'" and j + 1 < n and text[j + 1] == "'":
                    j += 2
                    continue
                if c == "'":
                    in_str = False
                j += 1
                continue
            if c == "'":
                in_str = True
                j += 1
                continue
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            elif c == ";" and depth == 0:
                break
            j += 1
        body = text[start:j].strip()
        p = body.find("(")
        if p < 0:
            etype, args = body, ""
        else:
            etype = body[:p].strip()
            close = len(body) - 1
            while close > p and body[close] in " \t":
                close -= 1
            args = body[p + 1 : close] if body[close] == ")" else body[p + 1 :]
        entities[eid] = (etype, args)
        i = j + 1
    return entities


def cart(entities: Dict[int, Tuple[str, str]], ref: Optional[int]) -> Optional[Tuple[float, float, float]]:
    if ref is None:
        return None
    _et, args = entities.get(ref, ("", ""))
    nums = [float(x) for x in NUM_RE.findall(args)]
    if len(nums) >= 3:
        return (nums[0], nums[1], nums[2])
    return None


def shoelace(pts: Sequence[Sequence[float]], n: Sequence[float]) -> float:
    clean = [tuple(p) for p in pts if p is not None]
    if len(clean) < 3:
        return 0.0
    ax, ay, az = n
    if abs(ax) < 0.9:
        rx, ry, rz = 1.0, 0.0, 0.0
    else:
        rx, ry, rz = 0.0, 1.0, 0.0
    ux, uy, uz = ay * rz - az * ry, az * rx - ax * rz, ax * ry - ay * rx
    lu = math.sqrt(ux * ux + uy * uy + uz * uz) or 1.0
    ux, uy, uz = ux / lu, uy / lu, uz / lu
    vx, vy, vz = n[1] * uz - n[2] * uy, n[2] * ux - n[0] * uz, n[0] * uy - n[1] * ux
    uv = [(p[0] * ux + p[1] * uy + p[2] * uz, p[0] * vx + p[1] * vy + p[2] * vz) for p in clean]
    cl: List[Tuple[float, float]] = []
    for p in uv:
        if not cl or abs(cl[-1][0] - p[0]) + abs(cl[-1][1] - p[1]) > 1e-12:
            cl.append(p)
    if len(cl) >= 2 and abs(cl[0][0] - cl[-1][0]) + abs(cl[0][1] - cl[-1][1]) < 1e-12:
        cl = cl[:-1]
    if len(cl) < 3:
        return 0.0
    cl.append(cl[0])
    a = 0.0
    for i in range(len(cl) - 1):
        a += cl[i][0] * cl[i + 1][1] - cl[i + 1][0] * cl[i][1]
    return abs(a) * 0.5


def _edge_pts(entities: Dict[int, Tuple[str, str]], eref: int) -> List[Tuple[float, float, float]]:
    et, args = entities.get(eref, ("", ""))
    if et == "ORIENTED_EDGE":
        curve = None
        for r in all_refs(args):
            if entities.get(r, ("", ""))[0] == "EDGE_CURVE":
                curve = r
                break
        if curve is None:
            refs = all_refs(args)
            curve = refs[-1] if refs else None
        if curve is None:
            return []
        return _edge_pts(entities, curve)
    if et != "EDGE_CURVE":
        return []
    parts = split_args(args)
    out: List[Tuple[float, float, float]] = []
    for r in (first_ref(parts[1]) if len(parts) > 1 else None,
              first_ref(parts[2]) if len(parts) > 2 else None):
        vet, vargs = entities.get(r or -1, ("", ""))
        if vet == "VERTEX_POINT":
            p = cart(entities, first_ref(vargs))
            if p:
                out.append(p)
    return out


def _face_geom(
    entities: Dict[int, Tuple[str, str]], bound_arg: str, n: Sequence[float]
) -> Tuple[float, Optional[Tuple[float, float, float]], int]:
    """Bound-polygon area, area-weighted centroid, loop count (FACE_*BOUND)."""
    bound_list = all_refs(bound_arg)
    kinds = [entities.get(b, ("", ""))[0] for b in bound_list]
    has_outer = "FACE_OUTER_BOUND" in kinds
    area = 0.0
    outer_used = False
    n_loops = 0
    c_acc = [0.0, 0.0, 0.0]
    w_acc = 0.0
    for bref in bound_list:
        bet, bargs = entities.get(bref, ("", ""))
        if bet not in ("FACE_OUTER_BOUND", "FACE_BOUND"):
            continue
        n_loops += 1
        bparts = split_args(bargs)
        lref = first_ref(bparts[1] if len(bparts) > 1 else (bparts[0] if bparts else None))
        if lref is None:
            continue
        pts: List[Tuple[float, float, float]] = []
        for eref in all_refs(entities.get(lref, ("", ""))[1]):
            for p in _edge_pts(entities, eref):
                if not pts or pts[-1] != p:
                    pts.append(p)
        a, c = _loop_area_centroid(pts, n)
        is_outer = bet == "FACE_OUTER_BOUND" or (not has_outer and not outer_used)
        if is_outer:
            area += a
            outer_used = True
            w = a
        else:
            area -= a
            w = a
        if c is not None and w > 1e-18:
            c_acc[0] += c[0] * w
            c_acc[1] += c[1] * w
            c_acc[2] += c[2] * w
            w_acc += w
    area = max(0.0, area)
    centroid: Optional[Tuple[float, float, float]] = None
    if w_acc > 1e-18:
        centroid = (c_acc[0] / w_acc, c_acc[1] / w_acc, c_acc[2] / w_acc)
    return area, centroid, n_loops


def _face_area(entities: Dict[int, Tuple[str, str]], bound_arg: str, n: Sequence[float]) -> float:
    return _face_geom(entities, bound_arg, n)[0]


@dataclass
class PlaneFace:
    id: int
    normal: Tuple[float, float, float]
    offset: float
    area: float
    centroid: Optional[Tuple[float, float, float]] = None
    nLoops: int = 1


def _uv_basis(n: Sequence[float]) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    ax, ay, az = n
    if abs(ax) < 0.9:
        rx, ry, rz = 1.0, 0.0, 0.0
    else:
        rx, ry, rz = 0.0, 1.0, 0.0
    ux, uy, uz = ay * rz - az * ry, az * rx - ax * rz, ax * ry - ay * rx
    lu = math.sqrt(ux * ux + uy * uy + uz * uz) or 1.0
    u = (ux / lu, uy / lu, uz / lu)
    v = (
        n[1] * u[2] - n[2] * u[1],
        n[2] * u[0] - n[0] * u[2],
        n[0] * u[1] - n[1] * u[0],
    )
    return u, v


def _loop_area_centroid(
    pts: Sequence[Sequence[float]], n: Sequence[float]
) -> Tuple[float, Optional[Tuple[float, float, float]]]:
    """Signed-magnitude polygon area and 3D centroid in the face plane (vertex shoelace)."""
    area = shoelace(pts, n)
    clean = [tuple(p) for p in pts if p is not None]
    if len(clean) < 3 or area <= 1e-18:
        if clean:
            m = len(clean)
            return area, (sum(p[0] for p in clean) / m, sum(p[1] for p in clean) / m, sum(p[2] for p in clean) / m)
        return area, None
    u, v = _uv_basis(n)
    uv = [(p[0] * u[0] + p[1] * u[1] + p[2] * u[2], p[0] * v[0] + p[1] * v[1] + p[2] * v[2]) for p in clean]
    cl: List[Tuple[float, float]] = []
    for p in uv:
        if not cl or abs(cl[-1][0] - p[0]) + abs(cl[-1][1] - p[1]) > 1e-12:
            cl.append(p)
    if len(cl) >= 2 and abs(cl[0][0] - cl[-1][0]) + abs(cl[0][1] - cl[-1][1]) < 1e-12:
        cl = cl[:-1]
    if len(cl) < 3:
        m = len(clean)
        return area, (sum(p[0] for p in clean) / m, sum(p[1] for p in clean) / m, sum(p[2] for p in clean) / m)
    cl.append(cl[0])
    a2 = cx = cy = 0.0
    for i in range(len(cl) - 1):
        cross = cl[i][0] * cl[i + 1][1] - cl[i + 1][0] * cl[i][1]
        a2 += cross
        cx += (cl[i][0] + cl[i + 1][0]) * cross
        cy += (cl[i][1] + cl[i + 1][1]) * cross
    if abs(a2) > 1e-18:
        cx /= 3.0 * a2
        cy /= 3.0 * a2
    else:
        body = cl[:-1]
        cx = sum(p[0] for p in body) / len(body)
        cy = sum(p[1] for p in body) / len(body)
    c3 = (cx * u[0] + cy * v[0], cx * u[1] + cy * v[1], cx * u[2] + cy * v[2])
    return area, c3


def list_step_planes(step: Path) -> List[PlaneFace]:
    """Every ADVANCED_FACE whose surface is PLANE, with equation and area."""
    text = step.read_text(encoding="utf-8", errors="replace")
    entities = parse_entities(text)
    out: List[PlaneFace] = []
    for eid, (et, args) in entities.items():
        if et != "ADVANCED_FACE":
            continue
        parts = split_args(args)
        sref = first_ref(parts[2]) if len(parts) >= 3 else None
        st, sargs = entities.get(sref or -1, ("", ""))
        if st != "PLANE":
            continue
        sp = split_args(sargs)
        ax = first_ref(sp[1] if len(sp) > 1 else (sp[0] if sp else None))
        _aet, aargs = entities.get(ax or -1, ("", ""))
        ap = split_args(aargs)
        loc = cart(entities, first_ref(ap[1]) if len(ap) > 1 else None)
        direc = cart(entities, first_ref(ap[2]) if len(ap) > 2 else None)
        if loc is None or direc is None:
            continue
        n = vnorm(direc)
        d = vdot(loc, n)
        area, centroid, n_loops = _face_geom(entities, parts[1] if len(parts) > 1 else "", n)
        if centroid is not None:
            t = d - vdot(centroid, n)
            centroid = (centroid[0] + t * n[0], centroid[1] + t * n[1], centroid[2] + t * n[2])
        out.append(PlaneFace(id=eid, normal=n, offset=d, area=area, centroid=centroid, nLoops=max(1, n_loops)))
    return out


def match_gt_plane(
    face: PlaneFace,
    gt_n: Sequence[float],
    gt_d: float,
    *,
    eps_plane: float,
    normal_agree: float,
) -> bool:
    gn = vnorm(gt_n)
    dp = vdot(face.normal, gn)
    if abs(dp) < normal_agree:
        return False
    gd = gt_d if dp >= 0.0 else -gt_d
    return abs(face.offset - gd) <= eps_plane + 1e-12


def load_gt_planes(gt: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    planes = gt.get("planes")
    match = gt.get("planeMatch") or {}
    if not isinstance(planes, list) or len(planes) != 18:
        raise GateError(f"ground-truth.json must list 18 GT planes (got {0 if not isinstance(planes, list) else len(planes)})")
    if "epsPlaneMM" not in match:
        raise GateError("ground-truth.json planeMatch.epsPlaneMM missing")
    return list(planes), dict(match)


def assert_ratchet_authority(ratchet: Dict[str, Any], path: Path) -> None:
    if path.name != RATCHET_NAME:
        raise GateError(f"floor authority must be {RATCHET_NAME}, opened {path}")
    if str(ratchet.get("authority") or "") != RATCHET_NAME:
        raise GateError(
            f"{path}: authority must be {RATCHET_NAME!r} (got {ratchet.get('authority')!r})"
        )
    if "perPlaneNFacesCeiling" not in ratchet:
        raise GateError(f"{path}: missing perPlaneNFacesCeiling")


IDENTITY_LOOP_WEIGHT = 1.0e6
IDENTITY_TIE_EPS = 1e-12


def _identity_score(face: PlaneFace, gid: int, identity: Dict[str, Any]) -> Optional[float]:
    """Distance from a shipped face to one co-equation GT. None = identity missing."""
    gi = identity.get(str(gid))
    if not isinstance(gi, dict) or not gi:
        return None
    score = 0.0
    used = False
    if "nLoops" in gi:
        used = True
        score += abs(int(gi["nLoops"]) - int(face.nLoops)) * IDENTITY_LOOP_WEIGHT
    if "areaMM2" in gi:
        used = True
        ga = float(gi["areaMM2"])
        score += abs(face.area - ga) / max(abs(ga), 1.0)
    if "centroid" in gi and face.centroid is not None:
        used = True
        c = gi["centroid"]
        score += math.sqrt(sum((float(a) - float(b)) ** 2 for a, b in zip(face.centroid, c))) / 100.0
    return score if used else None


def evaluate(
    faces: Sequence[PlaneFace],
    gt_planes: Sequence[Dict[str, Any]],
    match: Dict[str, Any],
    ratchet: Dict[str, Any],
    *,
    smooth_built_planes: int,
) -> Tuple[List[str], Dict[str, Any]]:
    eps = float(match["epsPlaneMM"])
    nmin = float(match.get("normalAgreeMin") or 0.99)
    area_floor = float(match["phantomAreaFloorMM2"])
    ceilings = {int(k): int(v) for k, v in (ratchet.get("perPlaneNFacesCeiling") or {}).items()}
    sum_ceil = int(ratchet["sumNFacesCeiling"])
    ph_ceil = int(ratchet["phantomPlaneCeiling"])
    identity = {str(k): v for k, v in (ratchet.get("gtFaceIdentity") or {}).items()}

    failures: List[str] = []
    assignment: Dict[int, int] = {}
    for f in faces:
        cands = [
            int(gp["id"])
            for gp in gt_planes
            if match_gt_plane(f, gp["normal"], float(gp["offset"]), eps_plane=eps, normal_agree=nmin)
        ]
        if not cands:
            continue
        if len(cands) == 1:
            assignment[f.id] = cands[0]
            continue
        scored: List[Tuple[float, int]] = []
        missing = False
        for gid in cands:
            sc = _identity_score(f, gid, identity)
            if sc is None:
                missing = True
                break
            scored.append((sc, gid))
        if missing or not scored:
            failures.append(
                f"STEP#{f.id} co-equation GTs {cands} have no gtFaceIdentity (D-S3-142)"
            )
            continue
        scored.sort(key=lambda t: (t[0], t[1]))
        if len(scored) >= 2 and abs(scored[0][0] - scored[1][0]) <= IDENTITY_TIE_EPS:
            failures.append(
                f"STEP#{f.id} identity tied between GT#{scored[0][1]} and GT#{scored[1][1]} "
                f"(score={scored[0][0]:.6g}; D-S3-142 named failure, not a coin flip)"
            )
            continue
        assignment[f.id] = scored[0][1]

    per: List[Dict[str, Any]] = []
    assigned: set[int] = set(assignment.keys())
    live_unique: set[int] = set()
    live_ids = {int(gp["id"]) for gp in gt_planes if gp.get("hasLiveRegion")}
    for gp in gt_planes:
        gid = int(gp["id"])
        hits = [f for f in faces if assignment.get(f.id) == gid]
        for f in hits:
            if gid in live_ids:
                live_unique.add(f.id)
        row = {
            "id": gid,
            "nFaces": len(hits),
            "area": sum(f.area for f in hits),
            "hasLiveRegionExpected": bool(gp.get("hasLiveRegion")),
            "hasLiveRegion": len(hits) > 0,
            "ceiling": int(ceilings.get(gid, 1)),
            "nLoops": [f.nLoops for f in hits],
        }
        per.append(row)

    phantoms = [f for f in faces if f.id not in assigned and f.area > area_floor + 1e-12]
    sum_live = len(live_unique)

    for row in per:
        if row["nFaces"] > row["ceiling"]:
            failures.append(
                f"GT#{row['id']} nFaces={row['nFaces']} > ceiling {row['ceiling']}"
            )
    if len(phantoms) > ph_ceil:
        failures.append(
            f"phantom-plane count={len(phantoms)} > ceiling {ph_ceil} "
            f"(area floor {area_floor} mm^2)"
        )
    if sum_live > sum_ceil:
        failures.append(
            f"sum nFaces (unique, live-region GT eqs)={sum_live} > ceiling {sum_ceil}"
        )

    details: Dict[str, Any] = {
        "epsPlaneMM": eps,
        "normalAgreeMin": nmin,
        "phantomAreaFloorMM2": area_floor,
        "smoothBuiltPlanes": int(smooth_built_planes),
        "nStepPlanes": len(faces),
        "perPlane": per,
        "sumNFacesLiveUnique": sum_live,
        "sumNFacesCeiling": sum_ceil,
        "phantomPlaneCount": len(phantoms),
        "phantomPlaneCeiling": ph_ceil,
        "phantomAreas": [round(f.area, 3) for f in phantoms[:12]],
        "antiConfusion": ANTI_CONFUSION,
        "assignment": {str(k): v for k, v in assignment.items()},
    }
    return failures, details


def run_trueform(binary: Path, stl: Path, step: Path, threads: int) -> Dict[str, Any]:
    cmd = [
        str(binary),
        str(stl),
        "-o",
        str(step),
        "--engine",
        "trueform",
        "--quiet",
        "--threads",
        str(threads),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        result = parse_result_line(proc.stdout)
    except GateError as exc:
        raise GateError(
            f"{exc}\nstdout={proc.stdout!r}\nstderr={proc.stderr!r}"
        ) from exc
    result["_exit"] = proc.returncode
    return result


def format_report(details: Dict[str, Any]) -> str:
    lines = [
        f"plane_census_gate  STEP PLANE faces={details['nStepPlanes']}  "
        f"smoothBuiltPlanes(regions)={details['smoothBuiltPlanes']}",
        f"  {ANTI_CONFUSION}",
        f"  epsPlaneMM={details['epsPlaneMM']} (planeMatch.epsPlaneMM)  "
        f"normalAgreeMin={details['normalAgreeMin']}  "
        f"phantomAreaFloorMM2={details['phantomAreaFloorMM2']}",
        "  per GT plane:",
    ]
    for row in details["perPlane"]:
        lines.append(
            f"    GT#{row['id']}: nFaces={row['nFaces']} <= {row['ceiling']}  "
            f"area={row['area']:.3f}  hasLiveRegion={row['hasLiveRegion']} "
            f"(expected {row['hasLiveRegionExpected']})"
        )
    lines.append(
        f"  sum unique live-region GT faces={details['sumNFacesLiveUnique']} "
        f"<= {details['sumNFacesCeiling']}"
    )
    lines.append(
        f"  phantom-plane (area>{details['phantomAreaFloorMM2']})="
        f"{details['phantomPlaneCount']} <= {details['phantomPlaneCeiling']}  "
        f"areas={details['phantomAreas']}"
    )
    return "\n".join(lines)


def run_gate(args: argparse.Namespace) -> int:
    gt_path = args.ground_truth.resolve()
    ratchet_path = args.ratchet.resolve()
    stl = args.stl.resolve()
    binary = args.binary.resolve() if args.binary else None
    gt = load_json(gt_path)
    ratchet = load_json(ratchet_path)
    gt_planes, match = load_gt_planes(gt)
    assert_ratchet_authority(ratchet, ratchet_path)
    if not stl.is_file():
        raise GateError(f"handle-pickup STL missing: {stl}")
    if binary is None or not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")

    with tempfile.TemporaryDirectory(prefix="plane_census_gate_") as td:
        step = Path(td) / "handle-pickup.step"
        result = run_trueform(binary, stl, step, args.threads)
        if not step.is_file():
            raise GateError(f"TrueForm wrote no STEP: {result.get('error') or result}")
        faces = list_step_planes(step)

    built = int(result.get("smoothBuiltPlanes") or 0)
    failures, details = evaluate(
        faces, gt_planes, match, ratchet, smooth_built_planes=built
    )
    details["ok"] = result.get("ok")
    print(format_report(details), flush=True)
    if not result.get("ok"):
        failures.append(f"conversion failed: {result.get('error')}")
    if failures:
        print("plane_census_gate FAIL", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print("plane_census_gate PASS")
    return 0


def _synthetic_step_two_coplanar() -> str:
    """Two 10x10 coplanar squares on z=0 (FACE_BOUND, OCCT style)."""
    return """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('plane_census_gate synthetic'),'2;1');
FILE_NAME('synth.step','2026-09-01',('gate'),('gate'),'x','x','x');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = DIRECTION('',(0.,0.,1.));
#3 = DIRECTION('',(1.,0.,0.));
#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5 = PLANE('',#4);
#10 = CARTESIAN_POINT('',(0.,0.,0.));
#11 = CARTESIAN_POINT('',(10.,0.,0.));
#12 = CARTESIAN_POINT('',(10.,10.,0.));
#13 = CARTESIAN_POINT('',(0.,10.,0.));
#14 = VERTEX_POINT('',#10);
#15 = VERTEX_POINT('',#11);
#16 = VERTEX_POINT('',#12);
#17 = VERTEX_POINT('',#13);
#18 = DIRECTION('',(1.,0.,0.));
#19 = VECTOR('',#18,1.);
#20 = LINE('',#10,#19);
#21 = DIRECTION('',(0.,1.,0.));
#22 = VECTOR('',#21,1.);
#23 = LINE('',#11,#22);
#24 = DIRECTION('',(-1.,0.,0.));
#25 = VECTOR('',#24,1.);
#26 = LINE('',#12,#25);
#27 = DIRECTION('',(0.,-1.,0.));
#28 = VECTOR('',#27,1.);
#29 = LINE('',#13,#28);
#30 = EDGE_CURVE('',#14,#15,#20,.T.);
#31 = EDGE_CURVE('',#15,#16,#23,.T.);
#32 = EDGE_CURVE('',#16,#17,#26,.T.);
#33 = EDGE_CURVE('',#17,#14,#29,.T.);
#34 = ORIENTED_EDGE('',*,*,#30,.T.);
#35 = ORIENTED_EDGE('',*,*,#31,.T.);
#36 = ORIENTED_EDGE('',*,*,#32,.T.);
#37 = ORIENTED_EDGE('',*,*,#33,.T.);
#38 = EDGE_LOOP('',(#34,#35,#36,#37));
#39 = FACE_BOUND('',#38,.T.);
#40 = ADVANCED_FACE('',(#39),#5,.T.);
#50 = CARTESIAN_POINT('',(20.,0.,0.));
#51 = CARTESIAN_POINT('',(30.,0.,0.));
#52 = CARTESIAN_POINT('',(30.,10.,0.));
#53 = CARTESIAN_POINT('',(20.,10.,0.));
#54 = VERTEX_POINT('',#50);
#55 = VERTEX_POINT('',#51);
#56 = VERTEX_POINT('',#52);
#57 = VERTEX_POINT('',#53);
#58 = LINE('',#50,#19);
#59 = LINE('',#51,#22);
#60 = LINE('',#52,#25);
#61 = LINE('',#53,#28);
#62 = EDGE_CURVE('',#54,#55,#58,.T.);
#63 = EDGE_CURVE('',#55,#56,#59,.T.);
#64 = EDGE_CURVE('',#56,#57,#60,.T.);
#65 = EDGE_CURVE('',#57,#54,#61,.T.);
#66 = ORIENTED_EDGE('',*,*,#62,.T.);
#67 = ORIENTED_EDGE('',*,*,#63,.T.);
#68 = ORIENTED_EDGE('',*,*,#64,.T.);
#69 = ORIENTED_EDGE('',*,*,#65,.T.);
#70 = EDGE_LOOP('',(#66,#67,#68,#69));
#71 = FACE_BOUND('',#70,.T.);
#72 = ADVANCED_FACE('',(#71),#5,.T.);
ENDSEC;
END-ISO-10303-21;
"""


def _self_test() -> int:
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    gt = load_json(DEFAULT_GT)
    ratchet = load_json(DEFAULT_RATCHET)
    planes, match = load_gt_planes(gt)
    assert_ratchet_authority(ratchet, DEFAULT_RATCHET)
    check(len(planes) == 18, "GT lists 18 planes")
    ids = [int(p["id"]) for p in planes]
    check(ids[0] == 1466 and ids[-1] == 1538, "GT plane ids bookend 1466..1538")
    check(1530 in ids and 1475 in ids and 1467 in ids, "named flats 1530/1475/1467 present")
    check(float(match["epsPlaneMM"]) > 0, "epsPlaneMM named in ground-truth.json")
    check(match.get("epsPlaneFormula") == "max(epsMesh, sewTol, 0.02)", "epsPlane formula is engine DerivedTols")
    check(float(match["phantomAreaFloorMM2"]) == 192.6, "phantom floor is GT#1505 192.6 mm^2")
    check(float(match.get("normalAgreeMin") or 0) == 0.99, "normalAgreeMin 0.99 from census §1")
    live = [p for p in planes if p.get("hasLiveRegion")]
    check(len(live) == 6, f"census hasLiveRegion true on 6 planes (got {len(live)})")
    check(any(int(p["id"]) == 1530 and p.get("hasLiveRegion") for p in planes), "1530 hasLiveRegion")
    check(any(int(p["id"]) == 1475 and p.get("hasLiveRegion") for p in planes), "1475 hasLiveRegion")

    ceil = ratchet["perPlaneNFacesCeiling"]
    check(int(ceil["1475"]) == 28, "1475 ceiling ratcheted at 28 (hub-2 exploded facets; co-eq 1484 split off)")
    check(int(ceil["1484"]) == 1, "1484 ceiling 1 (rid 3 two-loop face; D-S3-142 identity)")
    check(int(ceil["1467"]) == 32, "1467 ceiling ratcheted at 32 (epsPlane-measured; census §1 listed 33)")
    check(int(ceil["1530"]) == 1, "1530 ceiling 1")
    ident = ratchet.get("gtFaceIdentity") or {}
    check(int(ident["1475"]["nLoops"]) == 1, "gtFaceIdentity 1475 nLoops=1 (CAD no hole)")
    check(int(ident["1484"]["nLoops"]) == 2, "gtFaceIdentity 1484 nLoops=2 (CAD has hole)")
    check(int(ratchet["phantomPlaneCeiling"]) == 1, "phantom ceiling 1 (today non-zero)")
    check(int(ratchet["sumNFacesCeiling"]) == 77, "sum unique live-region faces 77")

    z_n = (0.0, 0.0, 1.0)
    gt_one = [{"id": 1, "normal": list(z_n), "offset": 0.0, "hasLiveRegion": True}]
    match_s = {
        "epsPlaneMM": 0.02,
        "normalAgreeMin": 0.99,
        "phantomAreaFloorMM2": 50.0,
    }
    two = [
        PlaneFace(id=10, normal=z_n, offset=0.0, area=100.0),
        PlaneFace(id=11, normal=z_n, offset=0.0, area=100.0),
    ]
    split_r = {
        "authority": RATCHET_NAME,
        "perPlaneNFacesCeiling": {1: 1},
        "sumNFacesCeiling": 1,
        "phantomPlaneCeiling": 0,
    }
    fails_split, det_split = evaluate(two, gt_one, match_s, split_r, smooth_built_planes=1)
    check(any("nFaces=2 > ceiling 1" in f for f in fails_split), "split plane trips nFaces ceiling")
    check(det_split["smoothBuiltPlanes"] == 1, "smoothBuiltPlanes printed in details")

    ok_r = dict(split_r)
    ok_r["perPlaneNFacesCeiling"] = {1: 2}
    ok_r["sumNFacesCeiling"] = 2
    fails_ok, _ = evaluate(two, gt_one, match_s, ok_r, smooth_built_planes=1)
    check(fails_ok == [], "matching nFaces=2 <= ceiling 2 PASSES")

    phantom = [
        PlaneFace(id=10, normal=z_n, offset=0.0, area=100.0),
        PlaneFace(id=99, normal=(1.0, 0.0, 0.0), offset=0.0, area=80.0),
    ]
    fails_ph, det_ph = evaluate(phantom, gt_one, match_s, ok_r, smooth_built_planes=1)
    check(det_ph["phantomPlaneCount"] == 1, "large unmatched face is a phantom")
    check(any("phantom-plane" in f for f in fails_ph), "phantom above ceiling 0 FAILS")
    tiny = [
        PlaneFace(id=10, normal=z_n, offset=0.0, area=100.0),
        PlaneFace(id=99, normal=(1.0, 0.0, 0.0), offset=0.0, area=1.0),
    ]
    fails_tiny, det_tiny = evaluate(tiny, gt_one, match_s, ok_r, smooth_built_planes=1)
    check(det_tiny["phantomPlaneCount"] == 0, "tiny unmatched face is not a phantom")
    check(fails_tiny == [], "tiny leftover does not fail phantom cell")

    with tempfile.TemporaryDirectory(prefix="plane_census_self_") as td:
        sp = Path(td) / "two.step"
        sp.write_text(_synthetic_step_two_coplanar(), encoding="utf-8")
        parsed = list_step_planes(sp)
        check(len(parsed) == 2, f"synthetic STEP has 2 PLANE faces (got {len(parsed)})")
        if parsed:
            check(abs(parsed[0].normal[2]) > 0.99, "synthetic normal is +Z")
            check(abs(parsed[0].offset) <= 0.02, "synthetic offset ~0")
            check(abs(parsed[0].area - 100.0) < 1.0, f"synthetic square area ~100 (got {parsed[0].area})")
        fails_syn, _ = evaluate(parsed, gt_one, match_s, split_r, smooth_built_planes=2)
        check(any("nFaces=2 > ceiling 1" in f for f in fails_syn), "synthetic STEP split trips ceiling")
        fails_syn_ok, _ = evaluate(parsed, gt_one, match_s, ok_r, smooth_built_planes=2)
        check(fails_syn_ok == [], "synthetic STEP with ceiling 2 PASSES")
        if parsed:
            check(all(f.nLoops == 1 for f in parsed), "synthetic squares have one loop each")

    gt_pair = [
        {"id": 10, "normal": list(z_n), "offset": 0.0, "hasLiveRegion": True},
        {"id": 11, "normal": list(z_n), "offset": 0.0, "hasLiveRegion": True},
    ]
    pair_r = {
        "authority": RATCHET_NAME,
        "perPlaneNFacesCeiling": {10: 1, 11: 1},
        "sumNFacesCeiling": 2,
        "phantomPlaneCeiling": 0,
        "gtFaceIdentity": {"10": {"nLoops": 1}, "11": {"nLoops": 2}},
    }
    two_id = [
        PlaneFace(id=100, normal=z_n, offset=0.0, area=100.0, nLoops=1),
        PlaneFace(id=101, normal=z_n, offset=0.0, area=200.0, nLoops=2),
    ]
    fails_id, det_id = evaluate(two_id, gt_pair, match_s, pair_r, smooth_built_planes=2)
    by_id = {r["id"]: r["nFaces"] for r in det_id["perPlane"]}
    check(fails_id == [], "co-equation unique nLoops assignment PASSES")
    check(by_id.get(10) == 1 and by_id.get(11) == 1, "co-equation nFaces are 1 and 1 (never 2/2)")
    check(det_id["assignment"].get("100") == 10 and det_id["assignment"].get("101") == 11,
          "1-loop face → GT#10, 2-loop face → GT#11")

    pair_no = {k: v for k, v in pair_r.items() if k != "gtFaceIdentity"}
    fails_no, _ = evaluate(two_id, gt_pair, match_s, pair_no, smooth_built_planes=2)
    check(any("gtFaceIdentity" in f or "tied" in f for f in fails_no),
          "co-equation without identity is a named failure")

    two_same = [
        PlaneFace(id=100, normal=z_n, offset=0.0, area=100.0, nLoops=1),
        PlaneFace(id=101, normal=z_n, offset=0.0, area=100.0, nLoops=1),
    ]
    pair_tie = dict(pair_r)
    pair_tie["gtFaceIdentity"] = {"10": {"nLoops": 1}, "11": {"nLoops": 1}}
    fails_tie, _ = evaluate(two_same, gt_pair, match_s, pair_tie, smooth_built_planes=2)
    check(any("tied" in f for f in fails_tie), "equal identity scores are a named tie, not a coin flip")

    check(ANTI_CONFUSION in format_report(det_split), "report includes anti-confusion sentence")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, help="stl2step CLI")
    p.add_argument("--stl", type=Path, default=DEFAULT_STL)
    p.add_argument("--ground-truth", type=Path, default=DEFAULT_GT)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--threads", type=int, default=0)
    p.add_argument("--self-test", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        return run_gate(args)
    except GateError as exc:
        print(f"plane_census_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
