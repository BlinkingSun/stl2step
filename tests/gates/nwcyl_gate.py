#!/usr/bin/env python3
"""nwcyl_gate — NW-CYL never-worse (D-S3-15, v4 D-S3-34 multiplicity-capped).

Rationale (D-S3-34): demoting 6-tri blend patches beyond the top-m is free;
killing a true round drops a large face out of the top-m and fires; losing 8 of
S02's 12 still fires; consolidation preserves coverage under a surviving top-m
face. Option (b) — exempting engine predicates from the gate — is refused.

Per part with a GT cylinder multiset in its sidecar, assert shipped analytic
cylinders do not regress vs tests/gates/nwcyl_baseline.json:

  nTriUnderGTMatchedCylinders(now) >= baseline  (multiplicity-capped coverage)
  nPhantom(now) <= baseline.nPhantom

GT-matched coverage (v4): per GT radius class g with multiplicity m(g), let
C(g) be shipped analytic cylinder faces matching g @ 0.3%, sorted by covered
triangle count descending; coverage = sum of triangles over the first m(g)
faces of C(g), summed over all g. Triangles come from regiondump cylinder
regions attributed to each census face (axis + UV overlap). uncappedCoverage
and faces are informational only (uncapped = sum all GT-matched face triangles
without the top-m cap).

Non-recoverable GT (meshRecoverable:false, e.g. HP hex-boss #1533/#1535) is
excluded from gtTotal / the missing listing but still matches as real GT —
never a phantom (D-S3-23 / hexnote).

GT matching reuses partial_recovery_gate.py (match_gt_radius @ 0.3 %,
split_cylinder_radii_all — do not fork rules).

RATChet procedure (landing lane updates baseline when a part improves):
  1. Run gate; note per-part capped tri / uncapped / faces / phantoms in PASS.
  2. If capped coverage rose or nPhantom fell vs baseline, edit
     nwcyl_baseline.json (bump engineRef to the landing commit).
  3. Commit baseline + gate in the same lane; never lower capped coverage or
     raise nPhantom in baseline without a measured regression fix elsewhere.

Usage:
  nwcyl_gate.py --self-test
  nwcyl_gate.py --measure --binary ./build/stl2step
  nwcyl_gate.py --binary ./build/stl2step [--baseline PATH] [--jobs N]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
CORPUS = REPO / "tests" / "corpus"
DEFAULT_BASELINE = REPO / "tests" / "gates" / "nwcyl_baseline.json"
BASELINE_NAME = "nwcyl_baseline.json"
HP_GT = REPO / "tests" / "diag" / "handle-pickup" / "ground-truth.json"

# Parts named in D-S3-15 / DECISION follow-up; only those with GT radii are gated.
CANDIDATE_PARTS: Tuple[str, ...] = (
    "handle-pickup",
    "handle-lock",
    "S02",
    "S04",
    "S05",
    "S09",
    "S11-b",
    "S16-R1-explode-success",
    "Body11",
)


class GateError(Exception):
    """Hard gate failure."""


def _load_prg():
    prg_path = REPO / "tests" / "gates" / "partial_recovery_gate.py"
    if not prg_path.is_file():
        raise GateError(f"partial_recovery_gate.py missing: {prg_path}")
    spec = importlib.util.spec_from_file_location("partial_recovery_gate", prg_path)
    if spec is None or spec.loader is None:
        raise GateError(f"cannot load {prg_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["partial_recovery_gate"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_json(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, doc: Dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


def sibling_bin(binary: Path, name: str) -> Path:
    return binary.parent / name


def load_gt_radii_all(part_id: str, sidecar: Mapping[str, Any]) -> Tuple[List[Any], str]:
    """Raw GT cylinder multiset (floats and/or annotated objects)."""
    if part_id == "handle-pickup" and HP_GT.is_file():
        doc = json.loads(HP_GT.read_text(encoding="utf-8"))
        raw_all = doc.get("cylinder_radii_all")
        if isinstance(raw_all, list) and raw_all:
            return list(raw_all), "diag/handle-pickup/ground-truth.json"

    raw = sidecar.get("cylinder_radii_all")
    if isinstance(raw, list) and raw:
        return list(raw), "sidecar.cylinder_radii_all"

    for row in sidecar.get("live") or []:
        if not isinstance(row, dict):
            continue
        row_all = row.get("cylinder_radii_all")
        if isinstance(row_all, list) and row_all:
            return list(row_all), "sidecar.live.cylinder_radii_all"
        row_r = row.get("cylinder_radii")
        if isinstance(row_r, list) and row_r:
            return list(row_r), "sidecar.live.cylinder_radii"

    out: List[Any] = []
    for item in sidecar.get("recoverable") or []:
        if not isinstance(item, dict) or item.get("type") != "cylinder":
            continue
        try:
            radius = float(item["radius"])
            count = int(item.get("count") or 1)
        except (KeyError, TypeError, ValueError):
            continue
        out.extend([radius] * max(1, count))
    if out:
        return out, "sidecar.recoverable"

    return [], "none"


def discover_gated_parts(corpus: Path) -> List[str]:
    gated: List[str] = []
    for pid in CANDIDATE_PARTS:
        sc_path = corpus / f"{pid}.stl"
        sc_json = corpus / f"{pid}.expected.json"
        if not sc_path.is_file() or not sc_json.is_file():
            continue
        sidecar = load_json(sc_json)
        gt_all, _ = load_gt_radii_all(pid, sidecar)
        if gt_all:
            gated.append(pid)
    return gated


def radii_near(
    a: float,
    b: float,
    *,
    match_gt_radius,
) -> bool:
    """True when two radii match within the GT gate tolerance (0.3%)."""
    return match_gt_radius(float(a), [float(b)]) is not None


def region_shipped_as_cylinder(
    region_radius: float,
    census_face_radii: Sequence[float],
    *,
    match_gt_radius,
) -> bool:
    return any(
        radii_near(region_radius, cr, match_gt_radius=match_gt_radius)
        for cr in census_face_radii
    )


def occt_census_cylinder_faces(cen: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Full cylinder face rows from stl2step_census (radius + axis + UV bounds)."""
    out: List[Dict[str, Any]] = []
    for row in cen.get("cylinders") or []:
        if isinstance(row, dict) and "radius" in row:
            out.append(dict(row))
    if out:
        return out
    for r in cen.get("cylinder_radii") or []:
        try:
            out.append({"radius": float(r)})
        except (TypeError, ValueError):
            continue
    return out


def occt_census_radii(cen: Dict[str, Any]) -> List[float]:
    return [float(f["radius"]) for f in occt_census_cylinder_faces(cen)]


def run_occt_census(census_bin: Path, step: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [str(census_bin), str(step)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise GateError(
            f"stl2step_census failed rc={proc.returncode}: {proc.stderr or proc.stdout}"
        )
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise GateError(f"stl2step_census output is not JSON: {exc}") from exc


def collect_cylinder_regions_from_dump(doc: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Flatten regiondump envelope or bare RegionSet into cylinder rows."""
    rows: List[Dict[str, Any]] = []

    def add_region(reg: Mapping[str, Any]) -> None:
        if str(reg.get("type") or "") != "cylinder":
            return
        tris = reg.get("tris") or []
        try:
            radius = float(reg["radius"])
        except (KeyError, TypeError, ValueError):
            return
        rows.append(
            {
                "radius": radius,
                "nTri": len(tris),
                "id": reg.get("id"),
                "ax": reg.get("ax"),
                "uMin": reg.get("uMin"),
                "uMax": reg.get("uMax"),
                "vMin": reg.get("vMin"),
                "vMax": reg.get("vMax"),
            }
        )

    if isinstance(doc.get("regions"), list):
        for reg in doc["regions"]:
            if isinstance(reg, dict):
                add_region(reg)
        return rows

    for comp in doc.get("comps") or []:
        if not isinstance(comp, dict):
            continue
        rs = comp.get("regionSet")
        if not isinstance(rs, dict):
            continue
        for reg in rs.get("regions") or []:
            if isinstance(reg, dict):
                add_region(reg)
    return rows


def run_regiondump(dump_bin: Path, stl: Path, out_path: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [str(dump_bin), str(stl), "--out", str(out_path)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not out_path.is_file():
        raise GateError(
            f"regiondump failed rc={proc.returncode}: {proc.stderr or proc.stdout}"
        )
    return load_json(out_path)


def _vec3(raw: Any) -> Optional[Tuple[float, float, float]]:
    if isinstance(raw, (list, tuple)) and len(raw) >= 3:
        return (float(raw[0]), float(raw[1]), float(raw[2]))
    if isinstance(raw, dict):
        try:
            return (float(raw["x"]), float(raw["y"]), float(raw["z"]))
        except (KeyError, TypeError, ValueError):
            return None
    return None


def _census_axis(face: Mapping[str, Any]) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    axis = face.get("axis")
    if not isinstance(axis, dict):
        return None
    loc = axis.get("location")
    if loc is None:
        loc = axis.get("loc")
    direction = axis.get("direction")
    if direction is None:
        direction = axis.get("dir")
    o = _vec3(loc)
    d = _vec3(direction)
    if o is None or d is None:
        return None
    return o, d


def _region_axis(reg: Mapping[str, Any]) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    ax = reg.get("ax")
    if not isinstance(ax, dict):
        return None
    o = _vec3(ax.get("loc"))
    d = _vec3(ax.get("dir"))
    if o is None or d is None:
        return None
    return o, d


def _norm_dir(d: Tuple[float, float, float]) -> Tuple[float, float, float]:
    n = (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) ** 0.5
    if n <= 0:
        return (0.0, 0.0, 0.0)
    return (d[0] / n, d[1] / n, d[2] / n)


def _canonical_axis(
    loc: Tuple[float, float, float], direction: Tuple[float, float, float]
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    d = _norm_dir(direction)
    ax, ay, az = abs(d[0]), abs(d[1]), abs(d[2])
    flip = False
    if az >= ax and az >= ay:
        flip = d[2] < 0
    elif ay >= ax:
        flip = d[1] < 0
    else:
        flip = d[0] < 0
    if flip:
        d = (-d[0], -d[1], -d[2])
    t = loc[0] * d[0] + loc[1] * d[1] + loc[2] * d[2]
    o = (loc[0] - t * d[0], loc[1] - t * d[1], loc[2] - t * d[2])
    return o, d


def _dirs_parallel(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> bool:
    dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    return abs(abs(dot) - 1.0) < 1e-8


def _axis_distance(
    a_loc: Tuple[float, float, float],
    a_dir: Tuple[float, float, float],
    b_loc: Tuple[float, float, float],
    b_dir: Tuple[float, float, float],
) -> float:
    dx = b_loc[0] - a_loc[0]
    dy = b_loc[1] - a_loc[1]
    dz = b_loc[2] - a_loc[2]
    cx = a_dir[1] * b_dir[2] - a_dir[2] * b_dir[1]
    cy = a_dir[2] * b_dir[0] - a_dir[0] * b_dir[2]
    cz = a_dir[0] * b_dir[1] - a_dir[1] * b_dir[0]
    cn = (cx * cx + cy * cy + cz * cz) ** 0.5
    if cn < 1e-12:
        px = dy * a_dir[2] - dz * a_dir[1]
        py = dz * a_dir[0] - dx * a_dir[2]
        pz = dx * a_dir[1] - dy * a_dir[0]
        return (px * px + py * py + pz * pz) ** 0.5
    return abs(dx * cx + dy * cy + dz * cz) / cn


def _interval_overlap(a_min: float, a_max: float, b_min: float, b_max: float) -> float:
    if not all(map(math.isfinite, (a_min, a_max, b_min, b_max))):
        return 0.0
    lo = max(min(a_min, a_max), min(b_min, b_max))
    hi = min(max(a_min, a_max), max(b_min, b_max))
    return max(0.0, hi - lo)


def normalize_census_faces(census: Sequence[Any]) -> List[Dict[str, Any]]:
    faces: List[Dict[str, Any]] = []
    for i, item in enumerate(census):
        if isinstance(item, (int, float)):
            faces.append({"radius": float(item), "_idx": i})
        elif isinstance(item, dict):
            faces.append(dict(item))
        else:
            continue
    return faces


def region_matches_census_face(
    reg: Mapping[str, Any],
    face: Mapping[str, Any],
    *,
    match_gt_radius,
    rad_tol: float = 1e-6,
    axis_tol: float = 1e-6,
) -> bool:
    try:
        r_rad = float(reg["radius"])
        f_rad = float(face["radius"])
    except (KeyError, TypeError, ValueError):
        return False
    if not radii_near(r_rad, f_rad, match_gt_radius=match_gt_radius):
        return False
    reg_ax = _region_axis(reg)
    face_ax = _census_axis(face)
    if reg_ax is None or face_ax is None:
        return True
    r_loc, r_dir = _canonical_axis(*reg_ax)
    f_loc, f_dir = _canonical_axis(*face_ax)
    if not _dirs_parallel(r_dir, f_dir):
        return False
    return _axis_distance(r_loc, r_dir, f_loc, f_dir) <= axis_tol


def census_face_match_score(reg: Mapping[str, Any], face: Mapping[str, Any]) -> float:
    """Higher = better geometric overlap for picking one face among same-radius hits."""
    try:
        ru0, ru1 = float(reg.get("uMin", 0)), float(reg.get("uMax", 0))
        rv0, rv1 = float(reg.get("vMin", 0)), float(reg.get("vMax", 0))
        fu0, fu1 = float(face.get("uMin", 0)), float(face.get("uMax", 0))
        fv0, fv1 = float(face.get("vMin", 0)), float(face.get("vMax", 0))
    except (TypeError, ValueError):
        return 0.0
    u_ov = _interval_overlap(ru0, ru1, fu0, fu1)
    v_ov = _interval_overlap(rv0, rv1, fv0, fv1)
    return u_ov * v_ov + u_ov + v_ov


def assign_region_triangles_to_faces(
    census_faces: Sequence[Mapping[str, Any]],
    cylinder_regions: Sequence[Mapping[str, Any]],
    *,
    match_gt_radius,
) -> List[int]:
    """Per census face index: triangle count from attributed mesh regions."""
    face_tri = [0] * len(census_faces)
    for reg in cylinder_regions:
        try:
            n_tri = int(reg.get("nTri") or 0)
        except (TypeError, ValueError):
            continue
        if n_tri <= 0:
            continue
        hits: List[int] = []
        scores: List[float] = []
        for i, face in enumerate(census_faces):
            if region_matches_census_face(
                reg, face, match_gt_radius=match_gt_radius
            ):
                hits.append(i)
                scores.append(census_face_match_score(reg, face))
        if not hits:
            continue
        best = hits[0]
        best_score = scores[0]
        for j, sc in zip(hits[1:], scores[1:]):
            if sc > best_score:
                best, best_score = j, sc
        if best_score <= 0.0:
            best = min(hits, key=lambda i: (face_tri[i], i))
        face_tri[best] += n_tri
    return face_tri


def compute_nwcyl(
    census_input: Sequence[Any],
    cylinder_regions: Sequence[Mapping[str, Any]],
    gt_all: Sequence[float],
    *,
    match_gt_radius,
    extra_match_radii: Sequence[float] = (),
) -> Tuple[int, int, int, int, List[Dict[str, Any]]]:
    """Return (nTriCapped, nTriUncapped, nPhantom, nGTMatchedFaces, per-radius rows).

    gt_all is the mesh-recoverable denominator (gtTotal / missing listing).
    extra_match_radii (non-recoverable GT) still match as real — never phantoms.
    """
    census_faces = normalize_census_faces(census_input)
    gt_counts: Counter = Counter(float(x) for x in gt_all)
    match_keys = list(gt_counts.keys()) + [
        float(x) for x in extra_match_radii if float(x) not in gt_counts
    ]

    n_phantom = 0
    built_match_faces: Counter = Counter()
    for face in census_faces:
        try:
            r = float(face["radius"])
        except (KeyError, TypeError, ValueError):
            n_phantom += 1
            continue
        g = match_gt_radius(r, match_keys)
        if g is None:
            n_phantom += 1
        else:
            built_match_faces[g] += 1

    n_faces_matched = sum(built_match_faces.values())
    face_tri = assign_region_triangles_to_faces(
        census_faces, cylinder_regions, match_gt_radius=match_gt_radius
    )

    n_tri_uncapped = 0
    n_tri_capped = 0
    tri_by_gt_uncapped: Counter = Counter()
    tri_by_gt_capped: Counter = Counter()

    for g, mult in gt_counts.items():
        face_tri_pairs: List[Tuple[int, int]] = []
        for i, face in enumerate(census_faces):
            try:
                r = float(face["radius"])
            except (KeyError, TypeError, ValueError):
                continue
            if match_gt_radius(r, [g]) is None:
                continue
            tri = face_tri[i]
            if tri > 0:
                face_tri_pairs.append((tri, i))
        face_tri_pairs.sort(key=lambda x: (-x[0], x[1]))
        uncapped = sum(t for t, _ in face_tri_pairs)
        capped = sum(t for t, _ in face_tri_pairs[:mult])
        tri_by_gt_uncapped[g] = uncapped
        tri_by_gt_capped[g] = capped
        n_tri_uncapped += uncapped
        n_tri_capped += capped

    per_radius: List[Dict[str, Any]] = []
    for g, mult in sorted(gt_counts.items()):
        per_radius.append(
            {
                "radius": g,
                "built": int(built_match_faces.get(g, 0)),
                "triangles": int(tri_by_gt_capped.get(g, 0)),
                "uncappedTriangles": int(tri_by_gt_uncapped.get(g, 0)),
                "gtMultiplicity": int(mult),
            }
        )
    return n_tri_capped, n_tri_uncapped, n_phantom, n_faces_matched, per_radius


def format_radius_table(per_radius: Sequence[Mapping[str, Any]]) -> str:
    bits = [
        f"R={row['radius']} tri={row['triangles']}/{row.get('uncappedTriangles', row['triangles'])} "
        f"face={row['built']}/{row['gtMultiplicity']}"
        for row in per_radius
    ]
    return " ".join(bits)


def convert_part(
    prg: Any,
    *,
    binary: Path,
    stl: Path,
    work: Path,
    verify: bool,
) -> Dict[str, Any]:
    work.mkdir(parents=True, exist_ok=True)
    step = work / f"{stl.stem}.step"
    art = prg.run_trueform(binary, stl, step, threads=0, verify=verify)
    if art.step is None or not art.step.is_file():
        raise GateError(
            f"{stl.stem}: TrueForm wrote no STEP: {art.result.get('error') or art.result}"
        )
    return {
        "step": art.step,
        "result": art.result,
        "exit_code": art.exit_code,
    }


def measure_part(
    prg: Any,
    *,
    part_id: str,
    corpus: Path,
    binary: Path,
    census_bin: Path,
    dump_bin: Path,
    work: Path,
) -> Dict[str, Any]:
    stl = corpus / f"{part_id}.stl"
    sidecar_path = corpus / f"{part_id}.expected.json"
    if not stl.is_file():
        raise GateError(f"{part_id}: STL missing: {stl}")
    if not sidecar_path.is_file():
        raise GateError(f"{part_id}: sidecar missing: {sidecar_path}")
    sidecar = load_json(sidecar_path)
    raw, gt_src = load_gt_radii_all(part_id, sidecar)
    recoverable, nonrec, all_r = prg.split_cylinder_radii_all(raw)
    if not all_r and not recoverable:
        raise GateError(f"{part_id}: no GT cylinder radii in sidecar")
    extra = [float(e.get("radius")) for e in nonrec]
    verify = part_id not in ("handle-lock", "Body11")
    conv = convert_part(prg, binary=binary, stl=stl, work=work, verify=verify)
    cen = run_occt_census(census_bin, conv["step"])
    dump_path = work / f"{part_id}.regions.json"
    dump_doc = run_regiondump(dump_bin, stl, dump_path)
    census_faces = occt_census_cylinder_faces(cen)
    cyl_regions = collect_cylinder_regions_from_dump(dump_doc)
    n_tri, n_uncapped, n_phantom, n_faces, per_radius = compute_nwcyl(
        census_faces,
        cyl_regions,
        recoverable,
        match_gt_radius=prg.match_gt_radius,
        extra_match_radii=extra,
    )
    return {
        "nTriUnderGTMatchedCylinders": n_tri,
        "uncappedCoverage": n_uncapped,
        "nPhantom": n_phantom,
        "faces": n_faces,
        "nGtTotal": len(recoverable),
        "nNotMeshRecoverable": len(nonrec),
        "gtSource": gt_src,
        "perRadius": per_radius,
        "nBuilt": len(census_faces),
    }


def load_baseline(path: Path) -> Dict[str, Any]:
    if path.name != BASELINE_NAME:
        raise GateError(f"baseline authority must be {BASELINE_NAME}, opened {path}")
    doc = load_json(path)
    if str(doc.get("authority") or "") != BASELINE_NAME:
        raise GateError(f"{path}: authority must be {BASELINE_NAME!r}")
    if "parts" not in doc or not isinstance(doc["parts"], dict):
        raise GateError(f"{path}: missing parts map")
    return doc


def baseline_tri_floor(row: Mapping[str, Any]) -> int:
    if "nTriUnderGTMatchedCylinders" in row:
        return int(row["nTriUnderGTMatchedCylinders"])
    # Legacy face-count baselines before D-S3-28 amendment.
    return int(row.get("nGTMatched") or 0)


def run_measure(args: argparse.Namespace) -> int:
    prg = _load_prg()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    census_bin = args.census.resolve()
    dump_bin = args.dump.resolve()
    if not census_bin.is_file():
        raise GateError(f"stl2step_census missing: {census_bin}")
    if not dump_bin.is_file():
        raise GateError(f"stl2step_regiondump missing: {dump_bin}")
    corpus = args.corpus.resolve()
    parts = discover_gated_parts(corpus)
    if not parts:
        raise GateError(f"no gated parts under {corpus}")
    doc: Dict[str, Any] = {
        "id": "nwcyl-baseline",
        "authority": BASELINE_NAME,
        "engineRef": args.engine_ref,
        "notes": (
            "NW-CYL v4 multiplicity-capped coverage (D-S3-34). Ratchet up when landing "
            "improves nTriUnderGTMatchedCylinders or lowers nPhantom — see nwcyl_gate.py "
            "header. uncappedCoverage and faces are informational. nGtTotal is "
            "mesh-recoverable GT only."
        ),
        "parts": {},
    }
    with tempfile.TemporaryDirectory(prefix="nwcyl_measure_") as td:
        work = Path(td)
        for pid in parts:
            row = measure_part(
                prg,
                part_id=pid,
                corpus=corpus,
                binary=binary,
                census_bin=census_bin,
                dump_bin=dump_bin,
                work=work / pid,
            )
            doc["parts"][pid] = {
                "nTriUnderGTMatchedCylinders": row["nTriUnderGTMatchedCylinders"],
                "uncappedCoverage": row["uncappedCoverage"],
                "nPhantom": row["nPhantom"],
                "faces": row["faces"],
                "nGtTotal": row["nGtTotal"],
                "gtSource": row["gtSource"],
            }
            print(
                f"  {pid}: tri={row['nTriUnderGTMatchedCylinders']} "
                f"uncapped={row['uncappedCoverage']} "
                f"face={row['faces']} phantom={row['nPhantom']} "
                f"gtTotal={row['nGtTotal']}  {format_radius_table(row['perRadius'])}"
            )
            if int(row.get("nNotMeshRecoverable") or 0):
                print(f"    GT not mesh-recoverable (n={row['nNotMeshRecoverable']})")
    out = args.baseline.resolve()
    write_json(out, doc)
    print(f"nwcyl_gate --measure wrote {out} ({len(parts)} parts)")
    return 0


def run_gate(args: argparse.Namespace) -> int:
    prg = _load_prg()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    census_bin = args.census.resolve()
    dump_bin = args.dump.resolve()
    if not census_bin.is_file():
        raise GateError(f"stl2step_census missing: {census_bin}")
    if not dump_bin.is_file():
        raise GateError(f"stl2step_regiondump missing: {dump_bin}")
    baseline = load_baseline(args.baseline.resolve())
    corpus = args.corpus.resolve()
    part_ids = sorted(baseline["parts"].keys())
    failures: List[str] = []

    with tempfile.TemporaryDirectory(prefix="nwcyl_gate_") as td:
        work = Path(td)
        workers = max(1, min(int(args.jobs) or 1, len(part_ids)))
        measured: Dict[str, Dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    measure_part,
                    prg,
                    part_id=pid,
                    corpus=corpus,
                    binary=binary,
                    census_bin=census_bin,
                    dump_bin=dump_bin,
                    work=work / pid,
                ): pid
                for pid in part_ids
            }
            for fut in as_completed(futs):
                pid = futs[fut]
                measured[pid] = fut.result()

    for pid in part_ids:
        base = baseline["parts"][pid]
        now = measured[pid]
        b_tri = baseline_tri_floor(base)
        b_phantom = int(base["nPhantom"])
        n_tri = int(now["nTriUnderGTMatchedCylinders"])
        n_uncapped = int(now.get("uncappedCoverage") or 0)
        n_phantom = int(now["nPhantom"])
        n_faces = int(now["faces"])
        gt_total = int(now["nGtTotal"])
        ok_tri = n_tri >= b_tri
        ok_phantom = n_phantom <= b_phantom
        status = "PASS" if ok_tri and ok_phantom else "FAIL"
        line = (
            f"  {pid}: tri={n_tri}/{b_tri} uncapped={n_uncapped} "
            f"face={n_faces} phantom={n_phantom}/{b_phantom} "
            f"gtTotal={gt_total} [{status}]"
        )
        print(line)
        print(f"    radii: {format_radius_table(now['perRadius'])}")
        n_non = int(now.get("nNotMeshRecoverable") or 0)
        if n_non:
            print(f"    GT not mesh-recoverable (n={n_non})")
        if not ok_tri:
            failures.append(
                f"{pid}: nTriUnderGTMatchedCylinders {n_tri} < baseline {b_tri}"
            )
        if not ok_phantom:
            failures.append(
                f"{pid}: nPhantom {n_phantom} > baseline {b_phantom}"
            )

    if failures:
        print("nwcyl_gate FAIL", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print(
        f"nwcyl_gate PASS  baseline={args.baseline} engineRef={baseline.get('engineRef')}"
    )
    return 0


def _self_test() -> int:
    prg = _load_prg()
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    gt_all = [2.0] * 12 + [5.0, 5.0]
    regions_good = [{"radius": 2.0, "nTri": 80}]
    census_good = [2.0] * 10 + [7.0, 7.0, 7.0]
    tri, unc, n_phantom, n_faces, per = compute_nwcyl(
        census_good,
        regions_good,
        gt_all,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri == 80, f"synthetic: triUnderGT=80 (got {tri})")
    check(unc == 80, f"synthetic: uncapped=80 (got {unc})")
    check(n_phantom == 3, f"synthetic: nPhantom=3 (got {n_phantom})")
    check(n_faces == 10, f"synthetic: faces=10 (got {n_faces})")
    check(per[0]["triangles"] == 80 and per[0]["built"] == 10, "R=2 tri 80 face 10/12")

    # Merge two GT-matched faces of the same radius: coverage unchanged, faces down.
    regions_merged = [{"radius": 2.0, "nTri": 80}]
    census_merged = [2.0] * 8 + [7.0]
    tri_m, unc_m, ph_m, face_m, _ = compute_nwcyl(
        census_merged,
        regions_merged,
        gt_all,
        match_gt_radius=prg.match_gt_radius,
    )
    baseline_row = {"nTriUnderGTMatchedCylinders": 80, "nPhantom": 3, "faces": 10}
    check(tri_m == 80 and face_m == 8, "merge same-GT faces: tri unchanged, faces down")
    check(
        tri_m >= baseline_row["nTriUnderGTMatchedCylinders"]
        and ph_m <= baseline_row["nPhantom"],
        "merge satisfies NW-CYL (coverage unchanged)",
    )

    # Demote a GT round to facets: region still in dump but no shipped cylinder face.
    regions_demote = [{"radius": 2.0, "nTri": 80}]
    census_demote = [7.0]
    tri_d, _, ph_d, face_d, _ = compute_nwcyl(
        census_demote,
        regions_demote,
        gt_all,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri_d == 0, f"demote GT round: triUnderGT=0 (got {tri_d})")
    check(
        not (
            tri_d >= baseline_row["nTriUnderGTMatchedCylinders"]
            and ph_d <= baseline_row["nPhantom"]
        ),
        "demote fails NW-CYL inequality",
    )

    # Phantom-only improvement.
    regions_ph = [{"radius": 2.0, "nTri": 80}]
    census_ph = [2.0] * 10
    tri_p, _, ph_p, face_p, _ = compute_nwcyl(
        census_ph,
        regions_ph,
        gt_all,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri_p == 80 and ph_p == 0, "phantom drop 3→0 with same coverage passes")
    check(
        tri_p >= baseline_row["nTriUnderGTMatchedCylinders"]
        and ph_p <= baseline_row["nPhantom"],
        "phantom-only improvement satisfies NW-CYL",
    )

    # Matched-face decrease without merge (lost GT face, region not shipped).
    regions_regress = [{"radius": 2.0, "nTri": 50}]
    census_regress = [2.0] * 8 + [7.0]
    tri_r, _, ph_r, _, _ = compute_nwcyl(
        census_regress,
        regions_regress,
        gt_all,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri_r == 50, "coverage regression 80→50")
    check(
        not (
            tri_r >= baseline_row["nTriUnderGTMatchedCylinders"]
            and ph_r <= baseline_row["nPhantom"]
        ),
        "coverage decrease fails NW-CYL inequality",
    )

    # v4: demote 50 tiny R=3 patches beyond top-21 — capped coverage unchanged.
    gt_hp_r3 = [3.0] * 21
    regions_hp = [{"radius": 3.0, "nTri": 100, "id": i} for i in range(21)] + [
        {"radius": 3.0, "nTri": 6, "id": 100 + i} for i in range(50)
    ]
    census_hp_full = [3.0] * 71
    tri_full, unc_full, _, _, _ = compute_nwcyl(
        census_hp_full,
        regions_hp,
        gt_hp_r3,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri_full == 2100, f"HP-like capped tri=2100 (got {tri_full})")
    check(unc_full == 2400, f"HP-like uncapped tri=2400 (got {unc_full})")
    census_hp_demoted = [3.0] * 21
    regions_hp_kept = regions_hp[:21]
    tri_dem, _, ph_dem, _, _ = compute_nwcyl(
        census_hp_demoted,
        regions_hp_kept,
        gt_hp_r3,
        match_gt_radius=prg.match_gt_radius,
    )
    check(
        tri_dem >= tri_full and ph_dem <= 27,
        "demote 50 tiny R=3 patches beyond top-21 must PASS",
    )

    # v4: kill one top-21 R=3 face — capped coverage drops.
    census_hp_kill = [3.0] * 20
    regions_hp_kill = regions_hp_kept[:20]
    tri_kill, _, _, _, _ = compute_nwcyl(
        census_hp_kill,
        regions_hp_kill,
        gt_hp_r3,
        match_gt_radius=prg.match_gt_radius,
    )
    check(tri_kill == 2000, f"kill one top-21 face tri=2000 (got {tri_kill})")
    check(
        not (tri_kill >= tri_full),
        "kill one top-21 R=3 face must FAIL NW-CYL",
    )

    sc = load_json(CORPUS / "handle-lock.expected.json")
    gt_hl, src_hl = load_gt_radii_all("handle-lock", sc)
    check(len(gt_hl) == 15, f"handle-lock GT multiset len=15 (got {len(gt_hl)})")
    check("cylinder_radii" in src_hl, "handle-lock GT from live cylinder_radii")

    sc_hp = load_json(CORPUS / "handle-pickup.expected.json")
    gt_hp, src_hp = load_gt_radii_all("handle-pickup", sc_hp)
    rec_hp, non_hp, all_hp = prg.split_cylinder_radii_all(gt_hp)
    check(len(all_hp) == 54, f"handle-pickup GT census still 54 (got {len(all_hp)})")
    check(len(rec_hp) == 52, f"handle-pickup recoverable gtTotal=52 (got {len(rec_hp)})")
    check(len(non_hp) == 2, f"handle-pickup hex note n=2 (got {len(non_hp)})")
    check("ground-truth" in src_hp, "HP multiplicity from diag ground-truth")

    extra = [float(e["radius"]) for e in non_hp]
    n_hex, unc_hex, p_hex, face_hex, per_hex = compute_nwcyl(
        [3.0] * 10 + [11.544154],
        [],
        rec_hp,
        match_gt_radius=prg.match_gt_radius,
        extra_match_radii=extra,
    )
    check(p_hex == 0, "R=11.544 shipped face is not a phantom")
    check(
        not any(abs(float(row["radius"]) - 11.544154) < 1e-6 for row in per_hex),
        "hex radius excluded from missing listing",
    )
    check(face_hex == 11, "hex face counts as a real GT match")

    try:
        load_baseline(Path("/tmp/not-nwcyl.json"))
        check(False, "wrong baseline name must FAIL")
    except GateError:
        check(True, "wrong baseline name FAIL")

    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path, help="stl2step CLI")
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    p.add_argument("--corpus", type=Path, default=CORPUS)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument(
        "--census",
        type=Path,
        help="stl2step_census witness (default: beside --binary)",
    )
    p.add_argument(
        "--dump",
        type=Path,
        help="stl2step_regiondump (default: beside --binary)",
    )
    p.add_argument(
        "--engine-ref",
        default="6948ec156b75aefaa5c3713bd951bc480cef5b52",
        help="commit recorded in baseline (measure mode)",
    )
    p.add_argument(
        "--measure",
        action="store_true",
        help="write baseline JSON from a fresh build (orchestrator only)",
    )
    p.add_argument("--self-test", action="store_true", help="API tests (no engine)")
    args = p.parse_args(argv)
    if args.binary is not None:
        binary = args.binary.resolve()
        if args.census is None:
            args.census = sibling_bin(binary, "stl2step_census")
        if args.dump is None:
            args.dump = sibling_bin(binary, "stl2step_regiondump")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.measure:
            if args.binary is None:
                raise GateError("--binary is required for --measure")
            return run_measure(args)
        if args.binary is None:
            raise GateError("--binary is required")
        return run_gate(args)
    except GateError as exc:
        print(f"nwcyl_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
