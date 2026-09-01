#!/usr/bin/env python3
"""fillet_exclusion_e0 — D3-9 excluded set E, test-side over regiondump JSON.

E is not in product code. This gate evaluates the D3-9 predicate over
`stl2step_regiondump` JSON (stderr `--diag` is the human summary; the JSON
is the datum). Exclusion is computed BEFORE any prismaticity condition.

Hard-green on handle-lock: seed (torus/sphere fit + island seed) is empty,
so |E|==0 by construction, DIAG_PRISM ok=1 failedCond=0, 15/15 cylinders.
That fence must go red if an equal-θ detector starts excluding the 15.

Handle pickup: GT target is 21 faces (expected-red at HEAD). Live |E| is
counted in region ids (over-segmented blends) and ratcheted. Spike S-1
measured 66; this tree ships whatever it measures.

Usage:
  fillet_exclusion_e0.py --self-test
  fillet_exclusion_e0.py --synthetic-pass
  fillet_exclusion_e0.py --binary ./build/stl2step --dump ./build/stl2step_regiondump
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

REPO = Path(__file__).resolve().parents[2]
PICKUP_STL = REPO / "tests" / "corpus" / "handle-pickup.stl"
HL_STL = REPO / "tests" / "corpus" / "handle-lock.stl"
S13_STL = REPO / "tests" / "corpus" / "S13.stl"
S14_STL = REPO / "tests" / "corpus" / "S14.stl"
PICKUP_GT = REPO / "tests" / "diag" / "handle-pickup" / "ground-truth.json"
PICKUP_LABELS = REPO / "tests" / "gates" / "labels" / "handle-pickup.tri-labels.json"
PICKUP_ANATOMY = REPO / "tests" / "gates" / "labels" / "handle-pickup.band-anatomy.json"
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "fillet-e0-ratchet.json"
RATCHET_NAME = "fillet-e0-ratchet.json"

ENV_STRICT = "FILLET_EXCLUSION_E0_STRICT"
MARKER_EXPECTED_RED = "FILLET_EXCLUSION_E0_EXPECTED_RED"
MARKER_UNEXPECTED = "FILLET_EXCLUSION_E0_UNEXPECTED"
MARKER_FLOORS_MET = "FILLET_EXCLUSION_E0_FLOORS_MET"

# 6 tori + 9 spheres + six non-axis rounds. fixture-leftover may land a list;
# this inline set is the documented fallback (SPEC §6.3).
FALLBACK_GT_E_FACES: Tuple[int, ...] = (
    # six non-axis rounds (SPEC inline list)
    1469,
    1479,
    1490,
    1495,
    1527,
    1536,
)
FALLBACK_NON_AXIS: Tuple[Tuple[int, float], ...] = (
    (1469, 3.0),
    (1479, 3.0),
    (1490, 3.0),
    (1495, 3.0),
    (1527, 3.0),
    (1536, 2.55),
)
GT_TARGET_FACES = 21
R_EXCLUDE_MAX = 3.5
EPS_PLANE_FLOOR = 0.05
HL_CYL_COUNT = 15
S13_FILLET_N = 1
S13_FILLET_R = 2.0
S13_FILLET_RTOL = 0.25
S14_FILLET_N = 1
S14_FILLET_R = 2.5
S14_FILLET_RTOL = 0.25

EXPECTED_RED_NAMES = frozenset({"pickup.gtFaces"})
SEED_TYPES = frozenset({"torus", "sphere"})
SEED_REJECTS = frozenset({"torusNYI", "sphereNYI"})


class GateError(Exception):
    """Hard gate failure."""


@dataclass
class Check:
    name: str
    ok: bool
    message: str
    group: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Exclusion:
    seed_ids: Set[int]
    e_ids: Set[int]
    seed_empty_structural: bool
    reason: str
    evaluated_before_prism: bool = True


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
    raise GateError("no RESULT line in stl2step stdout")


def parse_diag_kv(line: str) -> Optional[Tuple[str, Dict[str, str]]]:
    s = line.strip()
    if not s.startswith("DIAG_"):
        return None
    parts = s.split()
    if not parts:
        return None
    kv: Dict[str, str] = {}
    for tok in parts[1:]:
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        kv[k] = v
    return parts[0], kv


def vec3(v: Sequence[Any]) -> Tuple[float, float, float]:
    return float(v[0]), float(v[1]), float(v[2])


def vdot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vcross_norm(a: Sequence[float], b: Sequence[float]) -> float:
    x = a[1] * b[2] - a[2] * b[1]
    y = a[2] * b[0] - a[0] * b[2]
    z = a[0] * b[1] - a[1] * b[0]
    return math.sqrt(x * x + y * y + z * z)


def vnorm(a: Sequence[float]) -> Tuple[float, float, float]:
    m = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    if m <= 0.0:
        return (0.0, 0.0, 0.0)
    return (a[0] / m, a[1] / m, a[2] / m)


def tau_surf(weld: float, diag: float) -> float:
    return max(5e-5, max(4.0 * max(weld, 0.0), 1e-6 * max(diag, 0.0)))


def tau_ax(tau_s: float, regions: Sequence[Mapping[str, Any]]) -> float:
    h_min = 0.0
    for r in regions:
        if r.get("type") != "cylinder":
            continue
        h = abs(float(r.get("vMax") or 0.0) - float(r.get("vMin") or 0.0))
        if h > 0.0 and math.isfinite(h) and (h_min <= 0.0 or h < h_min):
            h_min = h
    if h_min > 0.0:
        return max(1e-6, 2.0 * tau_s / h_min)
    return 1e-6


def eps_plane(weld: float, sew: float, diag: float) -> float:
    eps_mesh = max(max(weld, 1e-4 * max(diag, 1.0)), 1e-3)
    return max(max(eps_mesh, sew), 0.02)


def epsilon_e(tau_s: float, eps_p: float) -> float:
    return max(tau_s, eps_p, EPS_PLANE_FLOOR)


def dominant_axis(regions: Sequence[Mapping[str, Any]]) -> Optional[Tuple[float, float, float]]:
    best = None
    best_h = -1.0
    for r in regions:
        if r.get("type") != "cylinder":
            continue
        h = abs(float(r.get("vMax") or 0.0) - float(r.get("vMin") or 0.0))
        if h > best_h:
            ax = (r.get("ax") or {}).get("dir")
            if isinstance(ax, list) and len(ax) == 3:
                best = vnorm(vec3(ax))
                best_h = h
    return best


def axis_parallel(
    a: Sequence[float], b: Sequence[float], tau: float
) -> bool:
    """Parallel within tauAx: |sin(theta)| < tauAx (prism tilt gate)."""
    return vcross_norm(vnorm(a), vnorm(b)) < tau


def all_regions(rs: Mapping[str, Any]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for r in rs.get("regions") or []:
        if isinstance(r, dict):
            out.append(r)
    for r in rs.get("rejected") or []:
        if isinstance(r, dict):
            out.append(r)
    return out


def is_seed_region(r: Mapping[str, Any]) -> bool:
    if str(r.get("type") or "") in SEED_TYPES:
        return True
    if str(r.get("reject") or "") in SEED_REJECTS:
        return True
    return False


def g1_neighbors(chains: Sequence[Mapping[str, Any]], rid: int) -> Set[int]:
    out: Set[int] = set()
    for ch in chains:
        if not ch.get("tangent"):
            continue
        a, b = ch.get("regA"), ch.get("regB")
        if a == rid and isinstance(b, int) and b >= 0:
            out.add(b)
        if b == rid and isinstance(a, int) and a >= 0:
            out.add(a)
    return out


def point_to_seed_dist(px: Sequence[float], seed: Mapping[str, Any]) -> float:
    """Distance of a point to a torus/sphere fit (sphere: |p-c|-R; else inf)."""
    ax = seed.get("ax") or {}
    loc = ax.get("loc")
    if not isinstance(loc, list) or len(loc) != 3:
        return float("inf")
    c = vec3(loc)
    r = float(seed.get("radius") or 0.0)
    kind = str(seed.get("type") or "")
    d = math.sqrt((px[0] - c[0]) ** 2 + (px[1] - c[1]) ** 2 + (px[2] - c[2]) ** 2)
    if kind == "sphere" or str(seed.get("reject") or "") == "sphereNYI":
        return abs(d - r)
    # Torus: major radius from Region::radius, no minor radius in the dump.
    # Use |radial-in-plane - R| as a lower bound; without Rminor this is the
    # tube-centre residual, not a full torus distance.
    direc = ax.get("dir")
    if not isinstance(direc, list) or len(direc) != 3:
        return abs(d - r)
    n = vnorm(vec3(direc))
    rel = (px[0] - c[0], px[1] - c[1], px[2] - c[2])
    axial = vdot(rel, n)
    rad = math.sqrt(max(0.0, vdot(rel, rel) - axial * axial))
    return abs(rad - r)


def compute_E(
    rs: Mapping[str, Any],
    *,
    weld: float = 0.0,
    sew: float = 0.0,
    diag: float = 0.0,
    island_points: Optional[Mapping[int, Sequence[Tuple[float, float, float]]]] = None,
) -> Exclusion:
    """D3-9 predicate. Evaluated before any prismaticity condition.

    Seed = regions fitted as torus/sphere (or rejected TorusNYI/SphereNYI).
    A component with no torus/sphere has an empty E by construction.
    Island membership uses optional per-island sample points (mesh coords);
    without them, island proximity is reported empty-because-no-seed or
    not-applied (JSON has island ids but no coordinates).
    """
    regions = [r for r in (rs.get("regions") or []) if isinstance(r, dict)]
    rejected = [r for r in (rs.get("rejected") or []) if isinstance(r, dict)]
    chains = [c for c in (rs.get("chains") or []) if isinstance(c, dict)]
    pool = regions + rejected
    by_id = {int(r["id"]): r for r in pool if "id" in r}

    seed_ids = {int(r["id"]) for r in pool if is_seed_region(r) and "id" in r}
    tau_s = tau_surf(weld, diag)
    eps_p = eps_plane(weld, sew, diag)
    eps = epsilon_e(tau_s, eps_p)
    t_ax = tau_ax(tau_s, regions)
    dom = dominant_axis(regions)

    # Facet islands within ε of a seed fit. No seed => no island can be near one.
    island_seed: Set[int] = set()
    if seed_ids and island_points:
        seeds = [by_id[i] for i in seed_ids if i in by_id]
        for iid, pts in island_points.items():
            hit = False
            for p in pts:
                for s in seeds:
                    if point_to_seed_dist(p, s) <= eps:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                island_seed.add(int(iid))

    e_ids: Set[int] = set(seed_ids)
    # Grow: cylinder R<=3.5 G1-tangent to a member of E (or another excluded
    # round) whose axis is NOT parallel, within tauAx, to the dominant axis.
    changed = True
    while changed:
        changed = False
        for r in regions:
            rid = int(r.get("id", -1))
            if rid in e_ids:
                continue
            if r.get("type") != "cylinder":
                continue
            radius = float(r.get("radius") or 0.0)
            if radius > R_EXCLUDE_MAX + 1e-12:
                continue
            nbrs = g1_neighbors(chains, rid)
            if not (nbrs & e_ids):
                # "or to another excluded round" — G1 to a cylinder already in E
                # is covered by nbrs & e_ids. Also allow G1 to another R<=3.5
                # cylinder that is itself a grow candidate once seeded.
                continue
            ax = (r.get("ax") or {}).get("dir")
            if dom is not None and isinstance(ax, list) and len(ax) == 3:
                if axis_parallel(vec3(ax), dom, t_ax):
                    continue
            e_ids.add(rid)
            changed = True

    structural_empty = len(seed_ids) == 0 and len(island_seed) == 0
    if structural_empty:
        reason = (
            "seed empty (no torus/sphere fit, no island seed) => |E|=0 by construction"
        )
        e_ids = set()
    else:
        reason = (
            f"seed={sorted(seed_ids)} islands={sorted(island_seed)} "
            f"grown={sorted(e_ids)}"
        )
    return Exclusion(
        seed_ids=seed_ids,
        e_ids=e_ids,
        seed_empty_structural=structural_empty,
        reason=reason,
        evaluated_before_prism=True,
    )


def load_gt_e_faces(
    gt_path: Optional[Path] = None,
    labels_path: Optional[Path] = None,
) -> Tuple[Set[int], str, bool]:
    """Return (gt face ids in E, source, used_fallback).

    fixture-leftover may land an explicit list on ground-truth.json
    (`excludedFaces` / `eFaces` / `d39ExcludedFaces`). Else tri-labels
    torus+sphere plus the SPEC inline non-axis six.
    """
    search: List[Path] = []
    if gt_path is not None:
        search.append(Path(gt_path))
    search.append(PICKUP_GT)
    for p in search:
        if not p.is_file():
            continue
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for key in ("excludedFaces", "eFaces", "d39ExcludedFaces", "gtExcludedFaces"):
            raw = doc.get(key)
            if isinstance(raw, list) and raw and all(isinstance(x, (int, float)) for x in raw):
                return {int(x) for x in raw}, f"{p}:{key}", False

    faces: Set[int] = set()
    src_bits = []
    labels = Path(labels_path) if labels_path is not None else PICKUP_LABELS
    if labels.is_file():
        try:
            lab = json.loads(labels.read_text(encoding="utf-8"))
            for row in lab.get("triangles") or []:
                if not isinstance(row, dict):
                    continue
                ft = str(row.get("face_type") or "")
                fid = row.get("true_face_id")
                if ft in SEED_TYPES and isinstance(fid, int):
                    faces.add(fid)
            src_bits.append(f"{labels}:torus+sphere")
        except (OSError, json.JSONDecodeError):
            pass
    if PICKUP_ANATOMY.is_file():
        try:
            anat = json.loads(PICKUP_ANATOMY.read_text(encoding="utf-8"))
            n_before = len(faces)
            for b in anat.get("bands") or []:
                if not isinstance(b, dict):
                    continue
                if str(b.get("type") or "") in SEED_TYPES and isinstance(b.get("true_face_id"), int):
                    faces.add(int(b["true_face_id"]))
            if len(faces) > n_before:
                src_bits.append(f"{PICKUP_ANATOMY}:torus+sphere bands")
        except (OSError, json.JSONDecodeError):
            pass
    faces.update(FALLBACK_GT_E_FACES)
    src_bits.append("SPEC inline non-axis 1469/1479/1490/1495/1527/1536")
    used_fallback = True
    return faces, " + ".join(src_bits), used_fallback


def live_e_from_labels(
    rs: Mapping[str, Any],
    gt_faces: Set[int],
    labels: Optional[Mapping[str, Any]],
) -> Set[int]:
    """Live region ids whose triangles map onto a GT E face."""
    if not labels:
        return set()
    tri_to_face: Dict[int, int] = {}
    for row in labels.get("triangles") or []:
        if not isinstance(row, dict):
            continue
        tid, fid = row.get("tri_id"), row.get("true_face_id")
        if isinstance(tid, int) and isinstance(fid, int):
            tri_to_face[tid] = fid
    out: Set[int] = set()
    for r in rs.get("regions") or []:
        if not isinstance(r, dict) or "id" not in r:
            continue
        for t in r.get("tris") or []:
            if int(t) in tri_to_face and tri_to_face[int(t)] in gt_faces:
                out.add(int(r["id"]))
                break
    return out


def gt_faces_hit_by_regions(
    rs: Mapping[str, Any],
    live_ids: Set[int],
    labels: Optional[Mapping[str, Any]],
) -> Set[int]:
    if not labels:
        return set()
    tri_to_face: Dict[int, int] = {}
    for row in labels.get("triangles") or []:
        if not isinstance(row, dict):
            continue
        tid, fid = row.get("tri_id"), row.get("true_face_id")
        if isinstance(tid, int) and isinstance(fid, int):
            tri_to_face[tid] = fid
    by_id = {int(r["id"]): r for r in (rs.get("regions") or []) if "id" in r}
    hit: Set[int] = set()
    for rid in live_ids:
        r = by_id.get(rid)
        if not r:
            continue
        for t in r.get("tris") or []:
            fid = tri_to_face.get(int(t))
            if fid is not None:
                hit.add(fid)
    return hit


def regionset_from_dump(doc: Mapping[str, Any]) -> Dict[str, Any]:
    comps = doc.get("comps") or []
    if not comps:
        return {"regions": [], "rejected": [], "chains": [], "stats": {}}
    # First clean component, else first with a regionSet.
    for c in comps:
        if c.get("clean") and isinstance(c.get("regionSet"), dict):
            return c["regionSet"]
    for c in comps:
        if isinstance(c.get("regionSet"), dict):
            return c["regionSet"]
    return {"regions": [], "rejected": [], "chains": [], "stats": {}}


def fillet_census(rs: Mapping[str, Any]) -> Dict[str, Any]:
    fils = [
        r
        for r in (rs.get("regions") or [])
        if isinstance(r, dict) and r.get("origin") == "filletStrip"
    ]
    radii = [float(r.get("radius") or 0.0) for r in fils]
    payload = {
        "n": len(fils),
        "radii": radii,
        "statsFillets": int((rs.get("stats") or {}).get("fillets") or 0),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["sha256"] = hashlib.sha256(blob.encode("utf-8")).hexdigest()
    return payload


def parse_diag_prism(stderr: str) -> Optional[Dict[str, str]]:
    for line in (stderr or "").splitlines():
        parsed = parse_diag_kv(line)
        if parsed and parsed[0] == "DIAG_PRISM":
            return parsed[1]
    return None


def load_ratchet(path: Path) -> Dict[str, Any]:
    if path.name != RATCHET_NAME:
        raise GateError(f"authority must be {RATCHET_NAME}, opened {path}")
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("id") != "fillet-e0-ratchet":
        raise GateError(f"{path} is not fillet-e0-ratchet (id={doc.get('id')!r})")
    if "liveECeiling" not in doc:
        raise GateError(f"{path}: missing liveECeiling")
    return doc


def evaluate_handle_lock(
    excl: Exclusion,
    prism: Optional[Mapping[str, str]],
    n_cyl: int,
) -> List[Check]:
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(name=name, ok=ok, message=message, group="handle-lock", hard=True, details=details)
        )

    add(
        "seedEmpty",
        excl.seed_empty_structural,
        (
            excl.reason
            if excl.seed_empty_structural
            else f"handle-lock seed NOT empty: {excl.reason} "
            "(equal-θ-style detector would kill the 15/15 fence)"
        ),
        seed=sorted(excl.seed_ids),
        e=sorted(excl.e_ids),
    )
    add(
        "E0",
        len(excl.e_ids) == 0 and excl.seed_empty_structural,
        (
            f"handle-lock |E|={len(excl.e_ids)} structural-empty={excl.seed_empty_structural}"
        ),
        nE=len(excl.e_ids),
    )
    add(
        "beforePrism",
        excl.evaluated_before_prism,
        "exclusion evaluated before any prismaticity condition",
    )
    if prism is None:
        add("DIAG_PRISM", False, "handle-lock DIAG_PRISM missing (not measured)")
    else:
        missing = [k for k in ("ok", "failedCond", "nCyl") if k not in prism]
        if missing:
            add(
                "DIAG_PRISM",
                False,
                f"handle-lock DIAG_PRISM missing key(s) {missing} (not measured)",
            )
        else:
            ok = int(float(prism["ok"])) == 1 and int(float(prism["failedCond"])) == 0
            add(
                "DIAG_PRISM",
                ok,
                (
                    f"DIAG_PRISM ok={prism['ok']} failedCond={prism['failedCond']}"
                    + ("" if ok else " (want ok=1 failedCond=0)")
                ),
            )
    add(
        "cyl15",
        n_cyl == HL_CYL_COUNT,
        (
            f"handle-lock cylinders={n_cyl} != {HL_CYL_COUNT}/15"
            if n_cyl != HL_CYL_COUNT
            else f"handle-lock cylinders={n_cyl}/15"
        ),
        nCyl=n_cyl,
    )
    return checks


def evaluate_pickup(
    excl: Exclusion,
    live_ids: Set[int],
    gt_faces: Set[int],
    gt_hit: Set[int],
    ratchet: Mapping[str, Any],
    *,
    used_fallback: bool,
    gt_src: str,
    strict: bool,
) -> List[Check]:
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, hard: bool = True, **details: Any) -> None:
        checks.append(
            Check(
                name=name,
                ok=ok,
                message=message,
                group="pickup",
                hard=hard,
                details=details,
            )
        )

    add(
        "beforePrism",
        excl.evaluated_before_prism,
        "pickup exclusion evaluated before any prismaticity condition",
    )
    n_live = len(live_ids)
    ceiling = int(ratchet["liveECeiling"])
    add(
        "liveE",
        n_live <= ceiling,
        (
            f"live |E| region-ids={n_live} > ratchet liveECeiling={ceiling}"
            if n_live > ceiling
            else f"live |E| region-ids={n_live} <= ratchet liveECeiling={ceiling}"
        ),
        nLive=n_live,
        liveECeiling=ceiling,
        structuralE=sorted(excl.e_ids),
        liveIds=sorted(live_ids)[:40],
    )
    if n_live < ceiling:
        add(
            "ratchetTighten",
            True,
            f"RATCHET-TIGHTEN: measured live |E|={n_live} < committed {ceiling}; "
            f"set {RATCHET_NAME} liveECeiling to {n_live} (do not auto-rewrite)",
        )
    n_gt = len(gt_faces)
    n_hit = len(gt_hit & gt_faces) if gt_hit else 0
    # The target is 21 GT faces, not live region ids. HEAD over-segments and
    # types no torus/sphere, so the 21-face target is expected-red.
    gt_ok = n_hit == GT_TARGET_FACES and n_gt == GT_TARGET_FACES
    msg = (
        f"pickup GT |E| faces hit={n_hit}/{n_gt} target={GT_TARGET_FACES} "
        f"(live region-ids={n_live}; S-1 measured 66; src={gt_src})"
    )
    if used_fallback:
        msg += " [WARNING: GT face list from SPEC inline fallback — fixture-leftover list not landed]"
    if not gt_ok:
        msg += (
            f" -- expected-red: GT target is {GT_TARGET_FACES} faces "
            f"(6 tori + 9 spheres + 6 non-axis rounds), live engine types no "
            f"torus/sphere and over-segments the blends"
        )
    add("pickup.gtFaces", gt_ok, msg, hard=strict, nHit=n_hit, nGt=n_gt, target=GT_TARGET_FACES)
    return checks


def evaluate_s13_s14(
    s13: Optional[Mapping[str, Any]],
    s14: Optional[Mapping[str, Any]],
    ratchet: Mapping[str, Any],
) -> List[Check]:
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(name=name, ok=ok, message=message, group="s13s14", hard=True, details=details)
        )

    if s13 is None:
        add("S13", False, "S13.stl missing — fillet census not measured")
    else:
        n = int(s13["n"])
        radii = list(s13["radii"])
        r_ok = n == S13_FILLET_N and any(abs(r - S13_FILLET_R) <= S13_FILLET_RTOL for r in radii)
        add(
            "S13",
            r_ok,
            (
                f"S13 fillet census n={n} radii={radii} "
                f"(want n={S13_FILLET_N} R~{S13_FILLET_R})"
            ),
            census=s13,
        )
        want_sha = ratchet.get("s13Sha256")
        if want_sha:
            add(
                "S13.sha",
                s13.get("sha256") == want_sha,
                (
                    f"S13 census sha {s13.get('sha256')} != ratchet {want_sha}"
                    if s13.get("sha256") != want_sha
                    else f"S13 census sha byte-stable {want_sha}"
                ),
            )
    if s14 is None:
        add("S14", False, "S14.stl missing — fillet census not measured")
    else:
        n = int(s14["n"])
        radii = list(s14["radii"])
        r_ok = n == S14_FILLET_N and any(abs(r - S14_FILLET_R) <= S14_FILLET_RTOL for r in radii)
        add(
            "S14",
            r_ok,
            (
                f"S14 fillet census n={n} radii={radii} "
                f"(want n={S14_FILLET_N} R~{S14_FILLET_R})"
            ),
            census=s14,
        )
        want_sha = ratchet.get("s14Sha256")
        if want_sha:
            add(
                "S14.sha",
                s14.get("sha256") == want_sha,
                (
                    f"S14 census sha {s14.get('sha256')} != ratchet {want_sha}"
                    if s14.get("sha256") != want_sha
                    else f"S14 census sha byte-stable {want_sha}"
                ),
            )
    return checks


def format_checks(checks: Sequence[Check]) -> str:
    lines = []
    for c in checks:
        if c.ok:
            tag = "PASS"
        elif c.hard:
            tag = "FAIL"
        else:
            tag = "EXPECTED-RED"
        lines.append(f"  [{c.group}/{tag}] {c.message}")
    return "\n".join(lines)


def failing_names(checks: Sequence[Check], *, hard_only: bool = True) -> List[str]:
    return [c.name for c in checks if (not c.ok) and (c.hard if hard_only else True)]


def split_hard_fails(checks: Sequence[Check]) -> Tuple[List[Check], List[Check]]:
    hard = [c for c in checks if (not c.ok) and c.hard]
    expected = [c for c in hard if c.name in EXPECTED_RED_NAMES]
    unexpected = [c for c in hard if c not in expected]
    return expected, unexpected


def run_dump(dump: Path, stl: Path, out: Path) -> Dict[str, Any]:
    cmd = [str(dump), str(stl), "--diag", "--out", str(out)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise GateError(
            f"regiondump failed rc={proc.returncode} {stl}: {proc.stderr or proc.stdout}"
        )
    doc = json.loads(out.read_text(encoding="utf-8"))
    doc["_stderr"] = proc.stderr
    return doc


def run_trueform_prism(binary: Path, stl: Path, step: Path) -> Tuple[Dict[str, Any], str]:
    cmd = [
        str(binary),
        str(stl),
        "-o",
        str(step),
        "--engine",
        "trueform",
        "--no-verify",
        "--quiet",
        "--threads",
        "0",
    ]
    environ = os.environ.copy()
    environ["STL2STEP_PRISM_DIAG"] = "1"
    proc = subprocess.run(cmd, capture_output=True, text=True, env=environ)
    result = parse_result_line(proc.stdout)
    return result, proc.stderr


def run_live(
    binary: Path,
    dump: Path,
    *,
    pickup_stl: Path = PICKUP_STL,
    hl_stl: Path = HL_STL,
    s13_stl: Path = S13_STL,
    s14_stl: Path = S14_STL,
    ratchet_path: Path = DEFAULT_RATCHET,
    jobs: int = 4,
    strict: Optional[bool] = None,
) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    if not dump.is_file():
        raise GateError(f"stl2step_regiondump missing: {dump}")
    for p, name in ((pickup_stl, "handle-pickup"), (hl_stl, "handle-lock")):
        if not p.is_file():
            raise GateError(f"{name} STL missing: {p}")
    is_strict = strict_enabled() if strict is None else bool(strict)
    ratchet = load_ratchet(ratchet_path)
    gt_faces, gt_src, used_fallback = load_gt_e_faces()
    labels = None
    if PICKUP_LABELS.is_file():
        labels = json.loads(PICKUP_LABELS.read_text(encoding="utf-8"))

    dump_jobs: Dict[str, Path] = {"pickup": pickup_stl, "handle-lock": hl_stl}
    if s13_stl.is_file():
        dump_jobs["S13"] = s13_stl
    if s14_stl.is_file():
        dump_jobs["S14"] = s14_stl

    dumps: Dict[str, Dict[str, Any]] = {}
    prism_art: Optional[Tuple[Dict[str, Any], str]] = None
    with tempfile.TemporaryDirectory(prefix="fillet_exclusion_e0_") as td:
        work = Path(td)
        workers = max(1, min(int(jobs) or 1, len(dump_jobs) + 1))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(run_dump, dump, stl, work / f"{name}.json"): ("dump", name)
                for name, stl in dump_jobs.items()
            }
            futs[
                pool.submit(run_trueform_prism, binary, hl_stl, work / "handle-lock.step")
            ] = ("prism", "handle-lock")
            for fut in as_completed(futs):
                kind, name = futs[fut]
                if kind == "dump":
                    dumps[name] = fut.result()
                else:
                    prism_art = fut.result()

    # Exclusion BEFORE prismaticity: compute E from dumps first, then read DIAG_PRISM.
    hl_rs = regionset_from_dump(dumps["handle-lock"])
    pk_rs = regionset_from_dump(dumps["pickup"])
    hl_excl = compute_E(hl_rs)
    pk_excl = compute_E(pk_rs)
    live_ids = live_e_from_labels(pk_rs, gt_faces, labels)
    if not live_ids:
        # Structural E is empty at HEAD (no torus/sphere typed). S-1's 66 is
        # GT-overlap region ids; without labels, ratchet against structural |E|.
        live_ids = set(pk_excl.e_ids)
    gt_hit = gt_faces_hit_by_regions(pk_rs, live_ids, labels)

    n_cyl = sum(1 for r in (hl_rs.get("regions") or []) if r.get("type") == "cylinder")
    prism = parse_diag_prism(prism_art[1] if prism_art else "")
    checks: List[Check] = []
    checks.extend(evaluate_handle_lock(hl_excl, prism, n_cyl))
    checks.extend(
        evaluate_pickup(
            pk_excl,
            live_ids,
            gt_faces,
            gt_hit,
            ratchet,
            used_fallback=used_fallback,
            gt_src=gt_src,
            strict=is_strict,
        )
    )
    s13c = fillet_census(regionset_from_dump(dumps["S13"])) if "S13" in dumps else None
    s14c = fillet_census(regionset_from_dump(dumps["S14"])) if "S14" in dumps else None
    checks.extend(evaluate_s13_s14(s13c, s14c, ratchet))

    expected_red, unexpected = split_hard_fails(checks)
    # pickup.gtFaces is expected-red even when strict (CMake invert). When
    # strict, it is hard; split_hard_fails then puts it in expected.
    if is_strict:
        expected_red, unexpected = split_hard_fails(checks)
    else:
        expected_red = [c for c in checks if (not c.ok) and c.name in EXPECTED_RED_NAMES]
        unexpected = [c for c in checks if (not c.ok) and c.hard and c.name not in EXPECTED_RED_NAMES]

    print(
        f"fillet_exclusion_e0  gt={gt_src} fallback={used_fallback}  "
        f"strict={is_strict}  liveE={len(live_ids)}  "
        f"gtHit={len(gt_hit & gt_faces)}/{len(gt_faces)}",
        flush=True,
    )
    print(format_checks(checks), flush=True)
    if unexpected:
        named = ", ".join(c.name for c in unexpected)
        print(f"{MARKER_UNEXPECTED}  unexpected fail(s): {named}", file=sys.stderr, flush=True)
        print(f"fillet_exclusion_e0 FAIL  unexpected: {named}", file=sys.stderr, flush=True)
        return 1
    if expected_red or any((not c.ok) and c.name in EXPECTED_RED_NAMES for c in checks):
        named = ", ".join(
            c.name for c in checks if (not c.ok) and c.name in EXPECTED_RED_NAMES
        )
        print(f"{MARKER_EXPECTED_RED}  unmet floors: {named}", flush=True)
        print(
            f"fillet_exclusion_e0 FAIL  failing assertion(s): {named}",
            file=sys.stderr,
            flush=True,
        )
        return 1
    print(f"{MARKER_FLOORS_MET}  D3-9 E-set MET (including 21 GT faces)", flush=True)
    print(
        "FLIP PROTOCOL: remove PASS_REGULAR_EXPRESSION from "
        "tests/gates/CMakeLists.txt (fillet_exclusion_e0). "
        "Keep FILLET_EXCLUSION_E0_STRICT=1. Do not redefine the 21-face target.",
        flush=True,
    )
    print("fillet_exclusion_e0 PASS", flush=True)
    return 0


def synthetic_handle_lock_rs() -> Dict[str, Any]:
    cyls = []
    for i in range(HL_CYL_COUNT):
        cyls.append(
            {
                "id": i,
                "type": "cylinder",
                "origin": "cylGrow",
                "radius": 5.0 + i,
                "ax": {"loc": [0, 0, 0], "dir": [0, 1, 0], "xdir": [1, 0, 0]},
                "vMin": 0.0,
                "vMax": 10.0,
                "tris": [],
                "reject": "none",
            }
        )
    return {"regions": cyls, "rejected": [], "chains": [], "stats": {"cylinders": 15, "fillets": 0}}


def synthetic_pickup_rs_head() -> Dict[str, Any]:
    """HEAD: no torus/sphere typed; blends over-segmented as small cylinders."""
    regs = []
    for i in range(10):
        regs.append(
            {
                "id": i,
                "type": "cylinder",
                "origin": "cylGrow",
                "radius": 3.0,
                "ax": {"loc": [0, 0, 0], "dir": [1, 0, 0], "xdir": [0, 1, 0]},
                "vMin": 0.0,
                "vMax": 2.0,
                "tris": [1000 + i],
                "reject": "none",
            }
        )
    return {"regions": regs, "rejected": [], "chains": [], "stats": {"cylinders": 10}}


def synthetic_pickup_rs_pass() -> Dict[str, Any]:
    regs = []
    for i, typ in enumerate(["torus"] * 6 + ["sphere"] * 9):
        regs.append(
            {
                "id": i,
                "type": typ,
                "origin": "cylGrow",
                "radius": 3.0,
                "ax": {"loc": [0, 0, 0], "dir": [0, 1, 0], "xdir": [1, 0, 0]},
                "vMin": 0.0,
                "vMax": 1.0,
                "tris": [i],
                "reject": "none",
            }
        )
    for j, fid in enumerate(FALLBACK_GT_E_FACES):
        regs.append(
            {
                "id": 100 + j,
                "type": "cylinder",
                "origin": "cylGrow",
                "radius": 3.0 if fid != 1536 else 2.55,
                "ax": {"loc": [0, 0, 0], "dir": [1, 0, 0], "xdir": [0, 1, 0]},
                "vMin": 0.0,
                "vMax": 1.0,
                "tris": [200 + j],
                "reject": "none",
            }
        )
    return {"regions": regs, "rejected": [], "chains": [], "stats": {}}


def synthetic_fillet_rs(n: int, radius: float) -> Dict[str, Any]:
    regs = [
        {
            "id": 0,
            "type": "cylinder",
            "origin": "filletStrip",
            "radius": radius,
            "ax": {"loc": [0, 0, 0], "dir": [1, 0, 0], "xdir": [0, 1, 0]},
            "vMin": 0,
            "vMax": 1,
            "tris": [],
            "reject": "none",
        }
        for _ in range(n)
    ]
    return {"regions": regs, "rejected": [], "chains": [], "stats": {"fillets": n}}


def synthetic_ratchet(**overrides: Any) -> Dict[str, Any]:
    s13 = fillet_census(synthetic_fillet_rs(1, 2.0))
    s14 = fillet_census(synthetic_fillet_rs(1, 2.5))
    doc: Dict[str, Any] = {
        "id": "fillet-e0-ratchet",
        "authority": RATCHET_NAME,
        "liveECeiling": 66,
        "s13Sha256": s13["sha256"],
        "s14Sha256": s14["sha256"],
        "notes": "synthetic",
    }
    doc.update(overrides)
    return doc


def run_synthetic_pass() -> int:
    """Green path: handle-lock E=0 structural + pickup 21 GT faces hit."""
    ratchet = synthetic_ratchet()
    hl = compute_E(synthetic_handle_lock_rs())
    prism = {"ok": "1", "failedCond": "0", "nCyl": "15"}
    checks = evaluate_handle_lock(hl, prism, 15)

    # Labels that map synthetic tris onto the 21 GT faces.
    gt_faces, _, _ = load_gt_e_faces()
    # Force the 21: 6 torus + 9 sphere synthetic ids 0..14 plus 6 non-axis.
    labels = {"triangles": []}
    rs = synthetic_pickup_rs_pass()
    # Map region tris to GT faces: first 15 regions -> 15 torus/sphere faces
    # from labels-or-fallback, last 6 -> FALLBACK_GT_E_FACES.
    torus_sphere = sorted(gt_faces - set(FALLBACK_GT_E_FACES))
    if len(torus_sphere) < 15:
        # Fallback path has only the six non-axis; synthesize 15 dummy GT ids
        # so the 21-count can be demonstrated, then also include the six.
        torus_sphere = list(range(9000, 9015))
        gt_faces = set(torus_sphere) | set(FALLBACK_GT_E_FACES)
    for i, fid in enumerate(torus_sphere[:15]):
        labels["triangles"].append({"tri_id": i, "true_face_id": fid, "face_type": "torus"})
    for j, fid in enumerate(FALLBACK_GT_E_FACES):
        labels["triangles"].append({"tri_id": 200 + j, "true_face_id": fid, "face_type": "cylinder"})
    live = live_e_from_labels(rs, gt_faces, labels)
    hit = gt_faces_hit_by_regions(rs, live, labels)
    excl = compute_E(rs)
    checks.extend(
        evaluate_pickup(
            excl,
            live,
            gt_faces,
            hit,
            {"liveECeiling": max(len(live), 66)},
            used_fallback=False,
            gt_src="synthetic",
            strict=True,
        )
    )
    checks.extend(
        evaluate_s13_s14(
            fillet_census(synthetic_fillet_rs(1, 2.0)),
            fillet_census(synthetic_fillet_rs(1, 2.5)),
            ratchet,
        )
    )
    print("fillet_exclusion_e0 --synthetic-pass")
    print(format_checks(checks))
    named = failing_names(checks)
    if named:
        print(
            f"fillet_exclusion_e0 FAIL  synthetic should pass; failed: {named}",
            file=sys.stderr,
        )
        return 1
    print("fillet_exclusion_e0 PASS  synthetic floors met")
    return 0


def _self_test() -> int:
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    hl = compute_E(synthetic_handle_lock_rs())
    check(hl.seed_empty_structural, "handle-lock seed empty is structural")
    check(len(hl.e_ids) == 0, "handle-lock |E|=0 by construction")
    check(hl.evaluated_before_prism, "E computed before prism")

    # Equal-θ detector simulation: type all 15 as torus → seed non-empty, fence dies.
    poisoned = synthetic_handle_lock_rs()
    for r in poisoned["regions"]:
        r["type"] = "torus"
    poisoned_e = compute_E(poisoned)
    check(not poisoned_e.seed_empty_structural, "torus-typed 15 makes seed non-empty")
    check(len(poisoned_e.e_ids) == 15, "equal-θ detector excludes all 15")
    prism_ok = {"ok": "1", "failedCond": "0", "nCyl": "15"}
    hl_checks = evaluate_handle_lock(hl, prism_ok, 15)
    check(failing_names(hl_checks) == [], "handle-lock healthy PASSES")
    poison_checks = evaluate_handle_lock(poisoned_e, prism_ok, 0)
    check("seedEmpty" in failing_names(poison_checks), "poisoned seed names seedEmpty")
    check("E0" in failing_names(poison_checks), "poisoned |E| names E0")
    check(
        any("equal-θ" in c.message or "seed NOT empty" in c.message for c in poison_checks),
        "poisoned fence message is actionable",
    )

    prism_bad = evaluate_handle_lock(hl, {"ok": "0", "failedCond": "3", "nCyl": "15"}, 15)
    check("DIAG_PRISM" in failing_names(prism_bad), "DIAG_PRISM ok=0 FAILS")
    prism_miss = evaluate_handle_lock(hl, None, 15)
    check("DIAG_PRISM" in failing_names(prism_miss), "missing DIAG_PRISM is not measured")
    cyl_bad = evaluate_handle_lock(hl, prism_ok, 14)
    check("cyl15" in failing_names(cyl_bad), "14/15 FAILS")

    # Grow: seed torus + G1-tangent R=3 cylinder off-axis.
    grow_rs = {
        "regions": [
            {
                "id": 1,
                "type": "torus",
                "radius": 4.0,
                "ax": {"loc": [0, 0, 0], "dir": [0, 1, 0], "xdir": [1, 0, 0]},
                "vMin": 0,
                "vMax": 1,
                "tris": [],
                "reject": "none",
            },
            {
                "id": 2,
                "type": "cylinder",
                "radius": 3.0,
                "ax": {"loc": [0, 0, 0], "dir": [1, 0, 0], "xdir": [0, 1, 0]},
                "vMin": 0,
                "vMax": 5,
                "tris": [],
                "reject": "none",
            },
            {
                "id": 3,
                "type": "cylinder",
                "radius": 20.0,
                "ax": {"loc": [0, 0, 0], "dir": [0, 1, 0], "xdir": [1, 0, 0]},
                "vMin": 0,
                "vMax": 50,
                "tris": [],
                "reject": "none",
            },
        ],
        "rejected": [],
        "chains": [{"regA": 1, "regB": 2, "tangent": True}],
        "stats": {},
    }
    ge = compute_E(grow_rs)
    check(1 in ge.seed_ids, "torus is seed")
    check(2 in ge.e_ids, "R=3 G1-tangent off-axis cylinder joins E")
    check(3 not in ge.e_ids, "R=20 dominant-axis cylinder stays out of E")

    # Head-like pickup: GT 21 expected-red.
    head_rs = synthetic_pickup_rs_head()
    head_e = compute_E(head_rs)
    check(head_e.seed_empty_structural, "HEAD pickup seed empty (no torus/sphere typed)")
    gt_faces, src, used_fb = load_gt_e_faces()
    check(used_fb or len(gt_faces) >= 6, "GT face list loaded (fallback or fixture)")
    check(set(FALLBACK_GT_E_FACES) <= gt_faces or len(gt_faces) == 21,
          "inline non-axis six present or fixture list is 21")
    pk_head = evaluate_pickup(
        head_e,
        set(),
        gt_faces if len(gt_faces) == 21 else (gt_faces | set(FALLBACK_GT_E_FACES)),
        set(),
        {"liveECeiling": 66},
        used_fallback=True,
        gt_src=src,
        strict=True,
    )
    check("pickup.gtFaces" in failing_names(pk_head), "HEAD-red names pickup.gtFaces")
    check(
        any("21" in c.message and "expected-red" in c.message for c in pk_head),
        "HEAD-red message cites the 21-face target and expected-red",
    )
    exp, unexp = split_hard_fails(pk_head)
    check(any(c.name == "pickup.gtFaces" for c in exp), "21-GT is expected-red")
    check(not any(c.name == "pickup.gtFaces" for c in unexp), "21-GT is not unexpected")

    live_hi = evaluate_pickup(
        head_e, set(range(70)), set(FALLBACK_GT_E_FACES), set(),
        {"liveECeiling": 66}, used_fallback=True, gt_src="x", strict=True,
    )
    check("liveE" in failing_names(live_hi), "live |E| > ratchet FAILS (hard)")

    s13 = fillet_census(synthetic_fillet_rs(1, 2.0))
    s14 = fillet_census(synthetic_fillet_rs(1, 2.5))
    rat = synthetic_ratchet()
    s_ok = evaluate_s13_s14(s13, s14, rat)
    check(failing_names(s_ok) == [], "S13/S14 census PASSES")
    s_bad = evaluate_s13_s14(fillet_census(synthetic_fillet_rs(0, 2.0)), s14, rat)
    check("S13" in failing_names(s_bad), "S13 n=0 FAILS")
    s_miss = evaluate_s13_s14(None, None, rat)
    check("S13" in failing_names(s_miss) and "S14" in failing_names(s_miss),
          "missing S13/S14 is not measured")

    check(strict_enabled({}) is True, "STRICT default ON")
    check(strict_enabled({ENV_STRICT: "0"}) is False, "STRICT=0 is OFF")
    rc = run_synthetic_pass()
    check(rc == 0, "--synthetic-pass returns 0")
    check(PICKUP_STL.is_file(), "handle-pickup.stl present")
    check(HL_STL.is_file(), "handle-lock.stl present")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path)
    p.add_argument("--dump", type=Path, help="stl2step_regiondump")
    p.add_argument("--pickup-stl", type=Path, default=PICKUP_STL)
    p.add_argument("--handle-lock-stl", type=Path, default=HL_STL)
    p.add_argument("--s13-stl", type=Path, default=S13_STL)
    p.add_argument("--s14-stl", type=Path, default=S14_STL)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--synthetic-pass", action="store_true")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--strict", dest="strict", action="store_true", default=None)
    g.add_argument("--no-strict", dest="strict", action="store_false")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.synthetic_pass:
            return run_synthetic_pass()
        if not args.binary or not args.dump:
            print("fillet_exclusion_e0: --binary and --dump are required", file=sys.stderr)
            return 1
        return run_live(
            Path(args.binary),
            Path(args.dump),
            pickup_stl=Path(args.pickup_stl),
            hl_stl=Path(args.handle_lock_stl),
            s13_stl=Path(args.s13_stl),
            s14_stl=Path(args.s14_stl),
            ratchet_path=Path(args.ratchet),
            jobs=args.jobs,
            strict=args.strict,
        )
    except GateError as exc:
        print(f"fillet_exclusion_e0 FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
