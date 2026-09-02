#!/usr/bin/env python3
"""partial_recovery_gate — Handle-pickup DECISION-boundary floors + measured PRG ratchet.

ctest name: partial_recovery_gate (owns the G0.1 *lane*; G0.1 itself stays in
run_gates.py and is not rewritten here).

Pickup floors are DECISION-boundary A1/A2 verbatim and stay floors:
smoothRevertedComponents==0, smoothBuiltComponents==1, planes>=18,
cylinders>=30, combined>=48, openShells==0, census.plane/cylinder same
floors, census.sphere/torus==0. At 2971d31 those floors PASS (analytic HP).

The three PRG cells are a MEASURED never-get-worse ratchet
(tests/gates/baseline/prg-ratchet.json), never a pass of the R4 targets:
  tightBudget       |V_step-V_mesh| <= ceiling (R4 formula budget still owed)
  census.radii      phantom count <= ceiling (R4 owes 0 phantoms)
  d3f7 R=3 built    <= ceiling (GT multiplicity 21 still owed; other
                    classes stay hard multiset-containment)

PARTIAL_RECOVERY_STRICT defaults ON. A cell that gets worse FAILS and
prints the current values. EXPECTED-RED invert is retired (D-S3-62).

FLIPPED 2026-09-02 at 2971d31: HP ships analytic; PRG cells ratcheted
from the first landing, R4 owes their reduction. CMake no longer uses
PASS_REGULAR_EXPRESSION. Keep ENVIRONMENT PARTIAL_RECOVERY_STRICT=1.
Do not lower FLOOR_*.

# D3F-4/6 (DECISION-floors-p3): FLOOR_CYLINDERS 50→30, FLOOR_COMBINED 100→48
# (charter combined>=100 STRUCK). Do not raise these here.

Usage:
  partial_recovery_gate.py --self-test
  partial_recovery_gate.py --synthetic-pass
  partial_recovery_gate.py --binary ./build/stl2step [--jobs 4]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
CENSUS_PY = REPO / "tests" / "tools" / "step_census.py"
HL_STL = REPO / "tests" / "corpus" / "handle-lock.stl"
BODY11_STL = REPO / "tests" / "corpus" / "Body11.stl"
PICKUP_CORPUS = REPO / "tests" / "corpus" / "handle-pickup.stl"
PICKUP_DIAG = REPO / "tests" / "diag" / "handle-pickup" / "handle-pickup.stl"
PICKUP_GT = REPO / "tests" / "diag" / "handle-pickup" / "ground-truth.json"
PICKUP_SIDECAR = REPO / "tests" / "corpus" / "handle-pickup.expected.json"
# Soft last fallbacks — skipped silently when the worktree has no _team/.
# Never hard-require repo-external _team/inputs (fresh-clone CI has none).
TEAM_PICKUP_STL = REPO / "_team" / "inputs" / "Handle pickup.stl"
TEAM_PICKUP_GT = REPO / "_team" / "inputs" / "ground-truth.json"

ENV_STRICT = "PARTIAL_RECOVERY_STRICT"
ENV_PRIVATE = "STL2STEP_PRIVATE_CORPUS"

# D3F-4/5/6 — Handle-pickup floors. combined = 18+30; charter >=100 STRUCK.
FLOOR_PLANES = 18
FLOOR_CYLINDERS = 30
FLOOR_COMBINED = 48
FLOOR_REVERTED = 0
FLOOR_BUILT_COMPONENTS = 1
FLOOR_OPEN_SHELLS = 0
ALLOWED_EXITS = (0, 2)
GT_REL_TOL = 0.003  # 0.3%
D3F7_R3 = 3.0

# Measured PRG ratchet (D-S3-62 flip at 2971d31). Authority is the JSON.
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "prg-ratchet.json"
RATCHET_NAME = "prg-ratchet.json"
RATCHET_ID = "prg-ratchet"
FLIP_REASON = (
    "re-derived 2026-09-02 on landing-#3 shell 9b693c8 (D-S3-133): "
    "phantom ceiling 18 attributed rid 164 (R=3.0594 BAD_B blend register)"
)

# Historical 18/30/48 family: still named on a revert (no CMake invert).
# census.sphere / census.torus are NOT in this set.
EXPECTED_RED_FLOOR_NAMES = frozenset(
    {
        "smoothRevertedComponents",
        "smoothBuiltComponents",
        "smoothBuiltPlanes",
        "smoothBuiltCylinders",
        "combined",
        "census.cylinder",
    }
)
MARKER_EXPECTED_RED = "PARTIAL_RECOVERY_EXPECTED_RED"
MARKER_UNEXPECTED = "PARTIAL_RECOVERY_UNEXPECTED"
MARKER_FLOORS_MET = "PARTIAL_RECOVERY_FLOORS_MET"
MARKER_RATCHET = "PRG_RATCHET"

# GROUND-TRUTH.md 9 spheres + 6 tori — mustRemainFaceted sidecar fallback
# (S02/S04 shape). Empty sidecar [] is treated as unpopulated (FINDINGS-0
# GAP4); fixture-pickup owns the file, this gate fills from GT.
FALLBACK_MUST_REMAIN_FACETED: Tuple[Dict[str, Any], ...] = (
    {
        "type": "sphere",
        "count": 9,
        "note": "corner balls; island/reject, not SphereNYI",
    },
    {
        "type": "torus",
        "count": 6,
        "note": "toroidal fillets; same",
    },
)

# GROUND-TRUTH.md unique cylinder radii (mm) — inline fallback.
FALLBACK_GT_RADII: Tuple[float, ...] = (
    0.5,
    2.0,
    2.55,
    2.825018,
    3.0,
    4.0,
    5.0,
    6.0,
    9.0,
    10.0,
    11.544154,
    25.0,
    72.5,
    87.5,
    92.5,
    100.0,
)

HEXNOTE_REASON = (
    "circumscribed cylinder of a hex boss; CAD faces #1533/#1535; "
    "mesh has 6 flats — see eprime-v2-step1.md §4"
)
HEXNOTE_RADIUS = 11.544154
HEXNOTE_ENTRY: Dict[str, Any] = {
    "radius": HEXNOTE_RADIUS,
    "meshRecoverable": False,
    "reason": HEXNOTE_REASON,
    "cadFaces": [1533, 1535],
}

# tests/diag/handle-pickup/ground-truth.json cylinder_radii_all (54, multiplicity).
# Two R=11.544154 entries are CAD circumscribed hex-boss cylinders (D-S3-23);
# annotated meshRecoverable:false — not deleted. Single source is the JSON;
# this fallback mirrors it when the file is absent.
FALLBACK_GT_RADII_ALL: Tuple[Any, ...] = (
    0.5,
    2.0,
    2.55,
    2.825018,
    2.825018,
    *([3.0] * 21),
    4.0,
    *([5.0] * 11),
    6.0,
    9.0,
    10.0,
    dict(HEXNOTE_ENTRY),
    dict(HEXNOTE_ENTRY),
    25.0,
    25.0,
    72.5,
    72.5,
    87.5,
    92.5,
    92.5,
    92.5,
    100.0,
    100.0,
    100.0,
)

# D3F-8 named gap — reported, not a floor. 12 deferred incl. three R=92.5
# perimeter walls, 2 over-budget, 5 no-bound, 6 rounds.
NAMED_GAP = (
    "D3F-8 named gap: 12 deferred (1468,1471,1473,1478,1481,1486,1489,1492,"
    "1494,1497,1501,1526 incl. three R=92.5 perimeter walls) / "
    "2 over-budget (1469,1477) / 5 no-bound (1532,1533,1534,1535,1539) / "
    "6 rounds (1469/1479/1490/1495/1527 R=3 + 1536 R=2.55). "
    "Charter combined>=100 STRUCK (GT max 72)."
)

HL_CENSUS_CYL_FLOOR = 15
BODY11_VOL_MAX = 0.023832
BODY11_SOLIDS = 2
BODY11_BUILT_FLOOR = 0  # never pin the stale 127 figure
BODY11_STALE_127 = 127  # named only so self-test can prove it is NOT a pin
BODY11_VERIFY = True  # D-S3-7: measure volumeDeltaPct; never --no-verify
BODY18_FREEEDGES_MAX = 18
BODY18_STL_NAME = "Body18.stl"

TIGHT_MESH_K = 1e-4
TIGHT_DVOL_K = 3.0

FREE_EDGES_RE = re.compile(r"freeEdges=(\d+)")
J6_MARK = "J6:"
REVERT_MARK = "analytic rebuild reverted on one component"


class GateError(Exception):
    """Hard gate failure (binary missing, no RESULT, census crash)."""


@dataclass
class Check:
    name: str
    ok: bool
    message: str
    group: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ConvertArt:
    id: str
    exit_code: int
    result: Dict[str, Any]
    stdout: str
    stderr: str
    step: Optional[Path] = None
    census: Optional[Dict[str, Any]] = None
    skipped: str = ""


def strict_enabled(env: Optional[Mapping[str, str]] = None) -> bool:
    """PARTIAL_RECOVERY_STRICT defaults ON (unset / empty => strict)."""
    environ: Mapping[str, str] = os.environ if env is None else env
    raw = environ.get(ENV_STRICT)
    if raw is None or str(raw).strip() == "":
        return True
    return str(raw).strip().lower() not in ("0", "false", "off", "no")


def load_prg_ratchet(path: Path = DEFAULT_RATCHET) -> Dict[str, Any]:
    """Floor authority for the three PRG cells. Missing/wrong file is a hard fail."""
    if path.name != RATCHET_NAME:
        raise GateError(f"floor authority must be {RATCHET_NAME}, opened {path}")
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("id") != RATCHET_ID:
        raise GateError(f"{path} is not the prg-ratchet (id={doc.get('id')!r})")
    if str(doc.get("authority") or "") != RATCHET_NAME:
        raise GateError(f"{path}: authority must be {RATCHET_NAME!r}")
    required = (
        "tightBudgetAbsDeltaMM3Ceiling",
        "censusRadiiPhantomCeiling",
        "d3f7R3BuiltCeiling",
        "reason",
    )
    missing = [k for k in required if k not in doc]
    if missing:
        raise GateError(f"{path}: missing {missing}")
    return doc


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise GateError("no RESULT line in stl2step stdout")


def parse_j6_free_edges(
    warnings: Optional[Sequence[str]] = None,
    stderr: str = "",
) -> int:
    vals: List[int] = []
    for w in warnings or []:
        if J6_MARK in w:
            m = FREE_EDGES_RE.search(w)
            if m:
                vals.append(int(m.group(1)))
    if not vals:
        for line in (stderr or "").splitlines():
            if J6_MARK in line:
                m = FREE_EDGES_RE.search(line)
                if m:
                    vals.append(int(m.group(1)))
    return vals[-1] if vals else 0


def pickup_stl_candidate_paths(extra: Optional[Path] = None) -> List[Path]:
    """In-repo copies first; worktree _team/inputs last (skipped if absent)."""
    out: List[Path] = []
    if extra is not None:
        out.append(Path(extra))
    out.extend([PICKUP_CORPUS, PICKUP_DIAG, TEAM_PICKUP_STL])
    return out


def resolve_pickup_stl(
    extra: Optional[Path] = None,
    *,
    candidates: Optional[Sequence[Path]] = None,
) -> Path:
    search = (
        [Path(p) for p in candidates]
        if candidates is not None
        else pickup_stl_candidate_paths(extra)
    )
    for p in search:
        if p.is_file():
            return p
    raise GateError(
        "handle-pickup STL missing; looked for "
        + ", ".join(str(p) for p in search)
    )


def pickup_gt_candidate_paths(gt_path: Optional[Path] = None) -> List[Path]:
    """Explicit path, in-repo diag GT, then _team/inputs (soft last fallback)."""
    out: List[Path] = []
    if gt_path is not None:
        out.append(Path(gt_path))
    out.append(PICKUP_GT)
    out.append(TEAM_PICKUP_GT)
    return out


def load_pickup_gt_radii(
    gt_path: Optional[Path] = None,
    *,
    candidates: Optional[Sequence[Path]] = None,
) -> Tuple[List[float], str]:
    """Read GT radii from diag ground-truth.json; else GROUND-TRUTH.md inline.

    Runtime search: explicit path, in-repo diag GT, then _team/inputs as a
    soft last fallback (skipped silently when absent). Pass `candidates` to
    isolate the search — required for the missing-GT self-test, which must
    not depend on repo or repo-external _team state.
    """
    search = (
        [Path(p) for p in candidates]
        if candidates is not None
        else pickup_gt_candidate_paths(gt_path)
    )
    for p in search:
        if not p.is_file():
            continue
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        raw = doc.get("cylinder_radii_unique") or doc.get("cylinder_radii")
        if not raw:
            continue
        return [float(x) for x in raw], str(p)
    return list(FALLBACK_GT_RADII), "GROUND-TRUTH.md inline fallback"


def load_pickup_gt_radii_all(
    gt_path: Optional[Path] = None,
    *,
    candidates: Optional[Sequence[Path]] = None,
) -> Tuple[List[float], str]:
    """54-entry cylinder_radii_all multiset (mixed floats / annotated objects)."""
    search = (
        [Path(p) for p in candidates]
        if candidates is not None
        else pickup_gt_candidate_paths(gt_path)
    )
    for p in search:
        if not p.is_file():
            continue
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        raw = doc.get("cylinder_radii_all")
        if not isinstance(raw, list) or not raw:
            continue
        return list(raw), str(p)
    return list(FALLBACK_GT_RADII_ALL), "GROUND-TRUTH.md inline fallback (radii_all)"


def load_must_remain_faceted(
    sidecar_path: Optional[Path] = None,
    *,
    candidates: Optional[Sequence[Path]] = None,
) -> Tuple[List[Dict[str, Any]], str]:
    """Sidecar mustRemainFaceted, or GROUND-TRUTH.md 9 spheres + 6 tori.

    Empty list is unpopulated (FINDINGS-0 GAP4), not "nothing stays faceted".
    Pass `candidates` to isolate the search from the in-repo sidecar.
    """
    search = (
        [Path(p) for p in candidates]
        if candidates is not None
        else (
            ([Path(sidecar_path)] if sidecar_path is not None else [])
            + [PICKUP_SIDECAR]
        )
    )
    for p in search:
        if not p.is_file():
            continue
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        raw = doc.get("mustRemainFaceted")
        if not isinstance(raw, list) or not raw:
            continue
        out: List[Dict[str, Any]] = []
        for item in raw:
            if isinstance(item, dict) and item.get("type"):
                out.append(dict(item))
        if out:
            return out, str(p)
    return [dict(x) for x in FALLBACK_MUST_REMAIN_FACETED], (
        "GROUND-TRUTH.md inline fallback"
    )


def match_gt_radius(radius: float, gt_radii: Sequence[float]) -> Optional[float]:
    best: Optional[float] = None
    best_err: Optional[float] = None
    for g in gt_radii:
        denom = abs(g)
        if denom == 0.0:
            if abs(radius) == 0.0:
                return 0.0
            continue
        err = abs(radius - g) / denom
        if err <= GT_REL_TOL and (best_err is None or err < best_err):
            best, best_err = g, err
    return best


def census_surface_count(cen: Dict[str, Any], kind: str) -> int:
    surfaces = cen.get("surfaces")
    if isinstance(surfaces, dict) and kind in surfaces:
        return int(surfaces.get(kind) or 0)
    return int(cen.get(kind) or 0)


def census_radii(cen: Dict[str, Any]) -> List[float]:
    out: List[float] = []
    for r in cen.get("cylinder_radii") or []:
        try:
            out.append(float(r))
        except (TypeError, ValueError):
            continue
    if out:
        return out
    for grp in cen.get("cylinderGroups") or []:
        try:
            out.append(float(grp.get("radius")))
        except (TypeError, ValueError):
            continue
    return out


def split_cylinder_radii_all(
    raw: Sequence[Any],
) -> Tuple[List[float], List[Dict[str, Any]], List[float]]:
    """Split cylinder_radii_all into recoverable / non-recoverable / all floats.

    Annotated objects with meshRecoverable:false stay in the JSON (not deleted)
    but are excluded from the d3f7 denominator and NW-CYL gtTotal / missing
    listing. Their radii still match as real GT (never phantoms).
    """
    recoverable: List[float] = []
    non_recoverable: List[Dict[str, Any]] = []
    all_radii: List[float] = []
    for item in raw:
        if isinstance(item, dict):
            try:
                radius = float(item["radius"])
            except (KeyError, TypeError, ValueError):
                continue
            rec = item.get("meshRecoverable", True)
            all_radii.append(radius)
            if rec is False:
                non_recoverable.append(dict(item))
            else:
                recoverable.append(radius)
        else:
            try:
                radius = float(item)
            except (TypeError, ValueError):
                continue
            recoverable.append(radius)
            all_radii.append(radius)
    return recoverable, non_recoverable, all_radii


def _d3f7_class_ceiling(
    g: float,
    n_gt: int,
    *,
    r3_built_ceiling: int,
    extra: Mapping[Any, Any],
) -> int:
    if match_gt_radius(float(g), (D3F7_R3,)) is not None:
        return int(r3_built_ceiling)
    for key, cap in (extra or {}).items():
        try:
            ek = float(key)
            ev = int(cap)
        except (TypeError, ValueError):
            continue
        if match_gt_radius(float(g), (ek,)) is not None:
            return ev
    return int(n_gt)


def radii_multiset_containment(
    built: Sequence[float],
    gt_all: Sequence[float],
    *,
    r3_built_ceiling: int,
    extra_class_ceilings: Optional[Mapping[Any, Any]] = None,
) -> Tuple[bool, str, Dict[str, Any]]:
    """D3F-7: per-GT-radius multiset-containment (multiplicity-aware, 0.3%).

    R=3 is a never-get-worse ratchet (built <= r3_built_ceiling). Other
    classes that over-segmented at the first landing may carry a measured
    extra ceiling; unlisted classes stay hard at GT multiplicity.
    Phantoms stay on census.radii.
    """
    extra = extra_class_ceilings or {}
    gt_counts: Counter = Counter(float(x) for x in gt_all)
    keys = list(gt_counts.keys())
    built_counts: Counter = Counter()
    for r in built:
        g = match_gt_radius(float(r), keys)
        if g is None:
            continue
        built_counts[g] += 1
    r3_key = match_gt_radius(D3F7_R3, keys)
    r3_built = int(built_counts.get(r3_key, 0)) if r3_key is not None else 0
    r3_gt_n = int(gt_counts.get(r3_key, 0)) if r3_key is not None else 0
    over: List[Dict[str, Any]] = []
    for g, n_gt in sorted(gt_counts.items()):
        n_b = int(built_counts.get(g, 0))
        cap = _d3f7_class_ceiling(
            g, n_gt, r3_built_ceiling=r3_built_ceiling, extra=extra
        )
        if n_b > cap:
            over.append(
                {
                    "gt": g,
                    "built": n_b,
                    "gtMultiplicity": n_gt,
                    "ceiling": cap,
                }
            )
    details = {
        "nBuilt": len(built),
        "nGtAll": len(gt_all),
        "over": over[:12],
        "r3Built": r3_built,
        "r3Ceiling": int(r3_built_ceiling),
        "r3GtMultiplicity": r3_gt_n,
    }
    if over:
        bits = ", ".join(
            f"R={row['gt']} built={row['built']}>ceil={row['ceiling']}"
            f"(gt={row['gtMultiplicity']})"
            for row in over[:6]
        )
        return (
            False,
            f"d3f7.radiiMultiset: over-segmentation (R=3 ratchet "
            f"{r3_built_ceiling}; other classes hard or landing-ratchet "
            f"@ 0.3%): {bits}  current R=3 built={r3_built}",
            details,
        )
    return (
        True,
        f"d3f7.radiiMultiset: {len(built)} built; R=3 built={r3_built} "
        f"<= ratchet {r3_built_ceiling} (gt multiplicity {r3_gt_n}; "
        f"R4 owes reduction)",
        details,
    )


def import_step_census():
    tools = str(CENSUS_PY.parent)
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import step_census  # type: ignore

    return step_census


def run_step_census(step: Path) -> Dict[str, Any]:
    """Census a STEP via tests/tools/step_census.py (same module SPEC names)."""
    if not CENSUS_PY.is_file():
        raise GateError(f"step_census.py missing: {CENSUS_PY}")
    mod = import_step_census()
    doc = mod.census_path(step)
    if not doc.get("ok"):
        raise GateError(doc.get("error") or f"census failed for {step}")
    return doc


def tight_budget(
    result: Mapping[str, Any],
    *,
    abs_ceiling: float,
) -> Tuple[bool, str, Dict[str, Any]]:
    """PRG ratchet: |V_step-V_mesh| vs landed ceiling (R4 formula still owed)."""
    vd = float(result.get("volumeDeltaPct", -1))
    mesh = float(result.get("meshVolumeMM3") or 0.0)
    step = float(result.get("stepVolumeMM3") or 0.0)
    pred = float(result.get("smoothVolPredictedMM3") or 0.0)
    details = {
        "volumeDeltaPct": vd,
        "meshVolumeMM3": mesh,
        "stepVolumeMM3": step,
        "smoothVolPredictedMM3": pred,
        "absCeilingMM3": float(abs_ceiling),
    }
    if vd < 0:
        return False, "tightBudget: volume not measured (need verify on)", details
    mesh_abs = abs(mesh)
    budget = max(TIGHT_MESH_K * mesh_abs, TIGHT_DVOL_K * abs(pred))
    abs_delta = abs(step - mesh)
    details["absDeltaMM3"] = abs_delta
    details["budgetMM3"] = budget
    ceil = float(abs_ceiling)
    if abs_delta > ceil + 1e-9:
        return (
            False,
            f"tightBudget: |V_step-V_mesh|={abs_delta:.6g} mm^3 > ratchet "
            f"{ceil:.6g} (formula {budget:.6g}; R4 owes reduction)",
            details,
        )
    return (
        True,
        f"tightBudget: |V_step-V_mesh|={abs_delta:.6g} <= ratchet {ceil:.6g} "
        f"(formula {budget:.6g}; R4 owes reduction)",
        details,
    )


def evaluate_pickup(
    result: Mapping[str, Any],
    census: Mapping[str, Any],
    gt_radii: Sequence[float],
    exit_code: int,
    *,
    strict: bool = True,
    gt_radii_all: Optional[Sequence[float]] = None,
    ratchet: Optional[Mapping[str, Any]] = None,
) -> List[Check]:
    """DECISION-boundary pickup floors + measured PRG ratchet. Each FAIL names the cell."""
    rch = dict(ratchet) if ratchet is not None else load_prg_ratchet()
    tb_ceil = float(rch["tightBudgetAbsDeltaMM3Ceiling"])
    ph_ceil = int(rch["censusRadiiPhantomCeiling"])
    r3_ceil = int(rch["d3f7R3BuiltCeiling"])
    extra_ceil = rch.get("d3f7ExtraClassCeilings") or {}
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(
                name=name,
                ok=ok,
                message=message,
                group="pickup",
                hard=strict,
                details=details,
            )
        )

    if not result.get("ok"):
        add(
            "conversion",
            False,
            f"conversion: pickup ok=false error={result.get('error')}",
            error=result.get("error"),
        )
        return checks

    exit_ok = int(exit_code) in ALLOWED_EXITS
    add(
        "exit",
        exit_ok,
        f"exit={exit_code} (want {ALLOWED_EXITS})"
        if not exit_ok
        else f"exit={exit_code} in {ALLOWED_EXITS}",
        exit=exit_code,
    )

    rev = int(result.get("smoothRevertedComponents") or 0)
    add(
        "smoothRevertedComponents",
        rev == FLOOR_REVERTED,
        f"smoothRevertedComponents={rev} != {FLOOR_REVERTED}"
        if rev != FLOOR_REVERTED
        else f"smoothRevertedComponents={rev}",
        got=rev,
        floor=FLOOR_REVERTED,
    )

    built_co = int(result.get("smoothBuiltComponents") or 0)
    add(
        "smoothBuiltComponents",
        built_co == FLOOR_BUILT_COMPONENTS,
        f"smoothBuiltComponents={built_co} != {FLOOR_BUILT_COMPONENTS}"
        if built_co != FLOOR_BUILT_COMPONENTS
        else f"smoothBuiltComponents={built_co}",
        got=built_co,
        floor=FLOOR_BUILT_COMPONENTS,
    )

    planes = int(result.get("smoothBuiltPlanes") or 0)
    add(
        "smoothBuiltPlanes",
        planes >= FLOOR_PLANES,
        f"smoothBuiltPlanes={planes} < floor {FLOOR_PLANES}"
        if planes < FLOOR_PLANES
        else f"smoothBuiltPlanes={planes} >= {FLOOR_PLANES}",
        got=planes,
        floor=FLOOR_PLANES,
    )

    cyls = int(result.get("smoothBuiltCylinders") or 0)
    add(
        "smoothBuiltCylinders",
        cyls >= FLOOR_CYLINDERS,
        f"smoothBuiltCylinders={cyls} < floor {FLOOR_CYLINDERS}"
        if cyls < FLOOR_CYLINDERS
        else f"smoothBuiltCylinders={cyls} >= {FLOOR_CYLINDERS}",
        got=cyls,
        floor=FLOOR_CYLINDERS,
    )

    combined = planes + cyls
    add(
        "combined",
        combined >= FLOOR_COMBINED,
        f"combined={combined} (planes {planes}+cylinders {cyls}) < floor {FLOOR_COMBINED}"
        if combined < FLOOR_COMBINED
        else f"combined={combined} >= {FLOOR_COMBINED}",
        got=combined,
        floor=FLOOR_COMBINED,
        planes=planes,
        cylinders=cyls,
    )

    open_s = int(result.get("openShells") or 0)
    add(
        "openShells",
        open_s == FLOOR_OPEN_SHELLS,
        f"openShells={open_s} != {FLOOR_OPEN_SHELLS}"
        if open_s != FLOOR_OPEN_SHELLS
        else f"openShells={open_s}",
        got=open_s,
        floor=FLOOR_OPEN_SHELLS,
    )

    tb_ok, tb_msg, tb_det = tight_budget(result, abs_ceiling=tb_ceil)
    add("tightBudget", tb_ok, tb_msg, **tb_det)

    cen_plane = census_surface_count(dict(census), "plane")
    add(
        "census.plane",
        cen_plane >= FLOOR_PLANES,
        f"census.plane={cen_plane} < floor {FLOOR_PLANES}"
        if cen_plane < FLOOR_PLANES
        else f"census.plane={cen_plane} >= {FLOOR_PLANES}",
        got=cen_plane,
        floor=FLOOR_PLANES,
    )

    cen_cyl = census_surface_count(dict(census), "cylinder")
    add(
        "census.cylinder",
        cen_cyl >= FLOOR_CYLINDERS,
        f"census.cylinder={cen_cyl} < floor {FLOOR_CYLINDERS}"
        if cen_cyl < FLOOR_CYLINDERS
        else f"census.cylinder={cen_cyl} >= {FLOOR_CYLINDERS}",
        got=cen_cyl,
        floor=FLOOR_CYLINDERS,
    )

    cen_sphere = census_surface_count(dict(census), "sphere")
    add(
        "census.sphere",
        cen_sphere == 0,
        (
            f"census.sphere={cen_sphere} != 0 (DECISION A1: spheres remain faceted)"
            if cen_sphere != 0
            else f"census.sphere={cen_sphere} (must remain faceted)"
        ),
        got=cen_sphere,
        floor=0,
    )

    cen_torus = census_surface_count(dict(census), "torus")
    add(
        "census.torus",
        cen_torus == 0,
        (
            f"census.torus={cen_torus} != 0 (DECISION A1: tori remain faceted)"
            if cen_torus != 0
            else f"census.torus={cen_torus} (must remain faceted)"
        ),
        got=cen_torus,
        floor=0,
    )

    radii = census_radii(dict(census))
    phantoms = [r for r in radii if match_gt_radius(r, gt_radii) is None]
    n_ph = len(phantoms)
    ph_ok = n_ph <= ph_ceil
    add(
        "census.radii",
        ph_ok,
        (
            f"census.radii: {n_ph} phantom(s) > ratchet {ph_ceil}: {phantoms[:8]}"
            if not ph_ok
            else (
                f"census.radii: {n_ph} phantom(s) <= ratchet {ph_ceil} "
                f"(nRadii={len(radii)}; R4 owes 0)"
            )
        ),
        phantoms=phantoms[:12],
        nRadii=len(radii),
        nPhantoms=n_ph,
        ceiling=ph_ceil,
    )

    gt_raw = (
        list(gt_radii_all) if gt_radii_all is not None else list(FALLBACK_GT_RADII_ALL)
    )
    recoverable, nonrec, _all_r = split_cylinder_radii_all(gt_raw)
    if nonrec:
        add(
            "gt.notMeshRecoverable",
            True,
            f"GT not mesh-recoverable (n={len(nonrec)})",
            n=len(nonrec),
            radii=[float(e.get("radius")) for e in nonrec],
        )
    ms_ok, ms_msg, ms_det = radii_multiset_containment(
        radii,
        recoverable,
        r3_built_ceiling=r3_ceil,
        extra_class_ceilings=extra_ceil,
    )
    # D3F-7 R=3 is the measured ratchet; other-class over-seg stays hard.
    checks.append(
        Check(
            name="d3f7.radiiMultiset",
            ok=ms_ok,
            message=ms_msg,
            group="pickup",
            hard=True,
            details=ms_det,
        )
    )
    return checks


def evaluate_handle_lock(
    result: Mapping[str, Any],
    census: Mapping[str, Any],
    exit_code: int,
) -> List[Check]:
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(
                name=name,
                ok=ok,
                message=message,
                group="handle-lock",
                hard=True,
                details=details,
            )
        )

    if not result.get("ok"):
        add(
            "conversion",
            False,
            f"handle-lock conversion failed: {result.get('error')}",
        )
        return checks
    if int(exit_code) not in ALLOWED_EXITS:
        add("exit", False, f"handle-lock exit={exit_code} not in {ALLOWED_EXITS}")
    cen_cyl = census_surface_count(dict(census), "cylinder")
    add(
        "census.cylinder",
        cen_cyl >= HL_CENSUS_CYL_FLOOR,
        f"handle-lock census.cylinder={cen_cyl} < floor {HL_CENSUS_CYL_FLOOR}"
        if cen_cyl < HL_CENSUS_CYL_FLOOR
        else f"handle-lock census.cylinder={cen_cyl} >= {HL_CENSUS_CYL_FLOOR}",
        got=cen_cyl,
        floor=HL_CENSUS_CYL_FLOOR,
    )
    return checks


def evaluate_body11(
    result: Mapping[str, Any],
    exit_code: int,
) -> List[Check]:
    """NEVER-GET-WORSE. built floor 0; do NOT pin the stale 127 figure."""
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(
                name=name,
                ok=ok,
                message=message,
                group="Body11",
                hard=True,
                details=details,
            )
        )

    if not result.get("ok"):
        add("conversion", False, f"Body11 conversion failed: {result.get('error')}")
        return checks
    # Never-get-worse: a clean exit 0 is better than today's 2; exit 1 is worse.
    if int(exit_code) not in ALLOWED_EXITS:
        add("exit", False, f"Body11 exit={exit_code} not in {ALLOWED_EXITS} (live HEAD is 2)")
    else:
        add("exit", True, f"Body11 exit={exit_code} (never-get-worse; live HEAD is 2)")

    vd = float(result.get("volumeDeltaPct", -1))
    built = int(result.get("smoothBuiltCylinders") or 0)
    if vd < 0 and built <= BODY11_BUILT_FLOOR:
        # Faceted ship is volume-exact. Live Body11 now converts with
        # verify=True (D-S3-7) so built>0 always has a measured vd.
        add(
            "volumeDeltaPct",
            True,
            f"Body11 volumeDeltaPct unmeasured; built={built} faceted path "
            f"(pin {BODY11_VOL_MAX} applies once anything is built)",
            got=vd,
            floor=BODY11_VOL_MAX,
            built=built,
        )
    else:
        vol_ok = vd >= 0 and vd <= BODY11_VOL_MAX + 1e-9
        add(
            "volumeDeltaPct",
            vol_ok,
            (
                f"Body11 volumeDeltaPct={vd} not measured while built={built}"
                if vd < 0
                else (
                    f"Body11 volumeDeltaPct={vd} > {BODY11_VOL_MAX}"
                    if not vol_ok
                    else f"Body11 volumeDeltaPct={vd} <= {BODY11_VOL_MAX}"
                )
            ),
            got=vd,
            floor=BODY11_VOL_MAX,
            built=built,
        )

    solids = int(result.get("solids") or 0)
    add(
        "solids",
        solids == BODY11_SOLIDS,
        f"Body11 solids={solids} != {BODY11_SOLIDS}"
        if solids != BODY11_SOLIDS
        else f"Body11 solids={solids}",
        got=solids,
        floor=BODY11_SOLIDS,
    )

    built = int(result.get("smoothBuiltCylinders") or 0)
    add(
        "smoothBuiltCylinders",
        built >= BODY11_BUILT_FLOOR,
        f"Body11 smoothBuiltCylinders={built} < floor {BODY11_BUILT_FLOOR}"
        if built < BODY11_BUILT_FLOOR
        else f"Body11 smoothBuiltCylinders={built} >= floor {BODY11_BUILT_FLOOR} "
        f"(127 is stale; not pinned)",
        got=built,
        floor=BODY11_BUILT_FLOOR,
        stale127=BODY11_STALE_127,
    )
    return checks


def evaluate_body18(free_edges: int) -> List[Check]:
    ok = free_edges <= BODY18_FREEEDGES_MAX
    return [
        Check(
            name="freeEdges",
            ok=ok,
            message=(
                f"Body18 freeEdges={free_edges} > {BODY18_FREEEDGES_MAX}"
                if not ok
                else f"Body18 freeEdges={free_edges} <= {BODY18_FREEEDGES_MAX}"
            ),
            group="Body18",
            hard=True,
            details={"got": free_edges, "floor": BODY18_FREEEDGES_MAX},
        )
    ]


def synthetic_pass_result() -> Dict[str, Any]:
    return {
        "ok": True,
        "solids": 1,
        "openShells": 0,
        "watertight": True,
        "meshVolumeMM3": 100000.0,
        "stepVolumeMM3": 100000.0,
        "volumeDeltaPct": 0.0,
        "smoothBuiltPlanes": FLOOR_PLANES,
        "smoothBuiltCylinders": FLOOR_COMBINED - FLOOR_PLANES,
        "smoothBuiltComponents": FLOOR_BUILT_COMPONENTS,
        "smoothRevertedComponents": FLOOR_REVERTED,
        "smoothVolPredictedMM3": 0.0,
        "warnings": [],
    }


def synthetic_pass_census() -> Dict[str, Any]:
    # First FLOOR_CYLINDERS recoverable GT radii — within multiplicity.
    recoverable, _non, _all = split_cylinder_radii_all(FALLBACK_GT_RADII_ALL)
    radii = list(recoverable[:FLOOR_CYLINDERS])
    return {
        "ok": True,
        "surfaces": {
            "plane": FLOOR_PLANES,
            "cylinder": FLOOR_CYLINDERS,
            "sphere": 0,
            "torus": 0,
        },
        "cylinder_radii": radii,
    }


def synthetic_head_fail_result() -> Dict[str, Any]:
    """HEAD @ 7cf77a2 shape: total collapse, component reverted."""
    return {
        "ok": True,
        "solids": 1,
        "openShells": 0,
        "watertight": True,
        "meshVolumeMM3": 100000.0,
        "stepVolumeMM3": 100000.0,
        "volumeDeltaPct": 0.0,
        "smoothBuiltPlanes": 0,
        "smoothBuiltCylinders": 0,
        "smoothBuiltComponents": 0,
        "smoothRevertedComponents": 1,
        "smoothVolPredictedMM3": 0.0,
        "warnings": ["smooth: analytic rebuild reverted on one component -- kept faceted"],
    }


def synthetic_head_fail_census() -> Dict[str, Any]:
    return {
        "ok": True,
        "surfaces": {"plane": 0, "cylinder": 0, "sphere": 0, "torus": 0},
        "cylinder_radii": [],
    }


def failing_names(checks: Sequence[Check], *, hard_only: bool = True) -> List[str]:
    out: List[str] = []
    for c in checks:
        if c.ok:
            continue
        if hard_only and not c.hard:
            continue
        out.append(c.name)
    return out


def split_hard_fails(
    checks: Sequence[Check],
) -> Tuple[List[Check], List[Check]]:
    """Split hard FAILs into expected-red floors vs unexpected (must stay loud)."""
    hard = [c for c in checks if (not c.ok) and c.hard]
    expected = [
        c
        for c in hard
        if c.group == "pickup" and c.name in EXPECTED_RED_FLOOR_NAMES
    ]
    unexpected = [c for c in hard if c not in expected]
    return expected, unexpected


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


def prg_ratchet_line(checks: Sequence[Check], ratchet: Mapping[str, Any]) -> str:
    """Always print current PRG cell values (PASS and FAIL)."""
    by_name = {c.name: c for c in checks if c.group == "pickup"}
    tb = by_name.get("tightBudget")
    cr = by_name.get("census.radii")
    d3 = by_name.get("d3f7.radiiMultiset")
    tb_d = tb.details if tb else {}
    cr_d = cr.details if cr else {}
    d3_d = d3.details if d3 else {}
    reason = str(ratchet.get("reason") or FLIP_REASON)
    return (
        f"{MARKER_RATCHET}  {reason}  "
        f"tightBudget |dV|={tb_d.get('absDeltaMM3', 'n/a')} "
        f"(ceil {ratchet.get('tightBudgetAbsDeltaMM3Ceiling')}, "
        f"formula {tb_d.get('budgetMM3', 'n/a')}); "
        f"census.radii phantoms={cr_d.get('nPhantoms', 'n/a')} "
        f"(ceil {ratchet.get('censusRadiiPhantomCeiling')}); "
        f"d3f7 R=3 built={d3_d.get('r3Built', 'n/a')} "
        f"(ceil {ratchet.get('d3f7R3BuiltCeiling')}, "
        f"gt {d3_d.get('r3GtMultiplicity', 'n/a')})"
    )


def resolve_body18_stl(
    private_corpus: Optional[Path] = None,
    env: Optional[Mapping[str, str]] = None,
) -> Optional[Path]:
    if private_corpus is not None:
        corpus = Path(private_corpus).expanduser()
    else:
        environ = os.environ if env is None else env
        raw = (environ.get(ENV_PRIVATE) or "").strip()
        if not raw:
            return None
        corpus = Path(raw).expanduser()
    if not corpus.is_dir():
        return None
    p = corpus / BODY18_STL_NAME
    return p if p.is_file() else None


def run_trueform(
    binary: Path,
    stl: Path,
    step: Path,
    *,
    threads: int,
    verify: bool,
) -> ConvertArt:
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
    if not verify:
        cmd.append("--no-verify")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        result = parse_result_line(proc.stdout)
    except GateError as exc:
        raise GateError(
            f"{exc}\nstdout={proc.stdout!r}\nstderr={proc.stderr!r}"
        ) from exc
    return ConvertArt(
        id=stl.stem,
        exit_code=proc.returncode,
        result=result,
        stdout=proc.stdout,
        stderr=proc.stderr,
        step=step if step.is_file() else None,
    )


def _convert_job(
    *,
    job_id: str,
    binary: Path,
    stl: Path,
    work: Path,
    threads: int,
    verify: bool,
    do_census: bool,
) -> ConvertArt:
    step = work / f"{job_id}.step"
    art = run_trueform(binary, stl, step, threads=threads, verify=verify)
    art.id = job_id
    if do_census and art.step is not None:
        art.census = run_step_census(art.step)
    return art


def run_live(
    binary: Path,
    *,
    pickup_stl: Optional[Path] = None,
    hl_stl: Path = HL_STL,
    body11_stl: Path = BODY11_STL,
    private_corpus: Optional[Path] = None,
    jobs: int = 4,
    gt_path: Optional[Path] = None,
    strict: Optional[bool] = None,
    ratchet_path: Optional[Path] = None,
) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    pickup = resolve_pickup_stl(pickup_stl)
    gt_radii, gt_src = load_pickup_gt_radii(gt_path)
    gt_radii_all, gt_all_src = load_pickup_gt_radii_all(gt_path)
    faceted, faceted_src = load_must_remain_faceted()
    is_strict = strict_enabled() if strict is None else bool(strict)
    rch = load_prg_ratchet(Path(ratchet_path) if ratchet_path else DEFAULT_RATCHET)

    tasks: List[Dict[str, Any]] = [
        {
            "job_id": "pickup",
            "stl": pickup,
            "verify": True,
            "do_census": True,
        },
        {
            "job_id": "handle-lock",
            "stl": hl_stl,
            "verify": False,
            "do_census": True,
        },
        {
            "job_id": "Body11",
            "stl": body11_stl,
            # D-S3-7: measure volumeDeltaPct (BODY11_VOL_MAX=0.023832).
            # --no-verify left RESULT vd=-1, then the built>0 branch
            # demanded a measured value — tests-only defect.
            "verify": BODY11_VERIFY,
            "do_census": False,
        },
    ]
    body18 = resolve_body18_stl(private_corpus)
    body18_skip = ""
    if body18 is not None:
        tasks.append(
            {
                "job_id": "Body18",
                "stl": body18,
                "verify": False,
                "do_census": False,
            }
        )
    else:
        body18_skip = (
            f"Body18 skipped: {ENV_PRIVATE} unset or {BODY18_STL_NAME} absent"
        )

    for t in tasks:
        stl = Path(t["stl"])
        if not stl.is_file():
            raise GateError(f"{t['job_id']} STL missing: {stl}")

    workers = max(1, min(int(jobs) or 1, len(tasks)))
    # 0 = all cores inside each CLI (last job is not starved at cpu//jobs).
    per_file = 0

    arts: Dict[str, ConvertArt] = {}
    with tempfile.TemporaryDirectory(prefix="partial_recovery_gate_") as td:
        work = Path(td)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    _convert_job,
                    job_id=t["job_id"],
                    binary=binary,
                    stl=Path(t["stl"]),
                    work=work,
                    threads=per_file,
                    verify=bool(t["verify"]),
                    do_census=bool(t["do_census"]),
                ): t["job_id"]
                for t in tasks
            }
            for fut in as_completed(futs):
                job_id = futs[fut]
                arts[job_id] = fut.result()

        checks: List[Check] = []
        pk = arts["pickup"]
        checks.extend(
            evaluate_pickup(
                pk.result,
                pk.census or {},
                gt_radii,
                pk.exit_code,
                strict=is_strict,
                gt_radii_all=gt_radii_all,
                ratchet=rch,
            )
        )
        hl = arts["handle-lock"]
        checks.extend(
            evaluate_handle_lock(hl.result, hl.census or {}, hl.exit_code)
        )
        b11 = arts["Body11"]
        checks.extend(evaluate_body11(b11.result, b11.exit_code))
        if "Body18" in arts:
            b18 = arts["Body18"]
            fe = parse_j6_free_edges(
                b18.result.get("warnings") or [], b18.stderr
            )
            checks.extend(evaluate_body18(fe))

        expected_red, unexpected = split_hard_fails(checks)
        soft_fails = [c for c in checks if (not c.ok) and not c.hard]
        n_sph = sum(int(e.get("count") or 0) for e in faceted if e.get("type") == "sphere")
        n_tor = sum(int(e.get("count") or 0) for e in faceted if e.get("type") == "torus")
        print(
            f"partial_recovery_gate  pickup={pickup}  gt={gt_src}  "
            f"gtAll={gt_all_src} n={len(split_cylinder_radii_all(gt_radii_all)[0])}  "
            f"faceted={faceted_src} (sphere={n_sph} torus={n_tor})  "
            f"strict={is_strict}  jobs={workers} threads/file={per_file}",
            flush=True,
        )
        print(f"  {NAMED_GAP}", flush=True)
        print(format_checks(checks), flush=True)
        print(prg_ratchet_line(checks, rch), flush=True)
        if body18_skip:
            print(f"  [Body18/SKIP] {body18_skip}", flush=True)
        if unexpected:
            named = ", ".join(c.name for c in unexpected)
            print(
                f"{MARKER_UNEXPECTED}  unexpected fail(s): {named}",
                file=sys.stderr,
                flush=True,
            )
            print(
                f"partial_recovery_gate FAIL  unexpected: {named}",
                file=sys.stderr,
                flush=True,
            )
            return 1
        if expected_red:
            named = ", ".join(c.name for c in expected_red)
            print(f"{MARKER_EXPECTED_RED}  unmet floors: {named}", flush=True)
            print(
                f"partial_recovery_gate FAIL  failing floor(s): {named}",
                file=sys.stderr,
                flush=True,
            )
            return 1
        if soft_fails:
            named = ", ".join(c.name for c in soft_fails)
            print(
                f"partial_recovery_gate EXPECTED-RED (PARTIAL_RECOVERY_STRICT=0)  "
                f"pickup floor(s) still open: {named}  "
                f"— contain race turns this green at integration",
                flush=True,
            )
            return 0
        print(f"{MARKER_FLOORS_MET}  pickup floors MET", flush=True)
        print("partial_recovery_gate PASS", flush=True)
        return 0


def run_synthetic_pass(strict: bool = True) -> int:
    gt_radii, gt_src = load_pickup_gt_radii()
    checks = evaluate_pickup(
        synthetic_pass_result(),
        synthetic_pass_census(),
        gt_radii,
        0,
        strict=strict,
    )
    fails = [c for c in checks if not c.ok]
    print(f"partial_recovery_gate --synthetic-pass  gt={gt_src}")
    print(format_checks(checks))
    if fails:
        named = ", ".join(c.name for c in fails)
        print(
            f"partial_recovery_gate FAIL  synthetic floors should pass; failed: {named}",
            file=sys.stderr,
        )
        return 1
    print("partial_recovery_gate PASS  synthetic floors met")
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

    check(FLOOR_PLANES == 18, f"FLOOR_PLANES is DECISION 18 (got {FLOOR_PLANES})")
    check(FLOOR_CYLINDERS == 30, f"FLOOR_CYLINDERS is DECISION 30 (got {FLOOR_CYLINDERS})")
    check(FLOOR_COMBINED == 48, f"FLOOR_COMBINED is DECISION 48 (got {FLOOR_COMBINED})")
    check(
        FLOOR_COMBINED == FLOOR_PLANES + FLOOR_CYLINDERS,
        "combined floor is 18+30=48 (charter combined>=100 STRUCK)",
    )
    check(len(FALLBACK_GT_RADII_ALL) == 54, f"GT radii_all multiplicity 54 (got {len(FALLBACK_GT_RADII_ALL)})")
    fb_rec, fb_non, fb_all = split_cylinder_radii_all(FALLBACK_GT_RADII_ALL)
    check(len(fb_all) == 54, f"split all-radii still 54 (got {len(fb_all)})")
    check(len(fb_rec) == 52, f"d3f7 denominator excludes hex (got {len(fb_rec)})")
    check(len(fb_non) == 2, f"GT not mesh-recoverable n=2 (got {len(fb_non)})")
    check(
        all(float(e.get("radius")) == HEXNOTE_RADIUS for e in fb_non),
        "non-recoverable entries are R=11.544154",
    )
    check("12 deferred" in NAMED_GAP and "R=92.5" in NAMED_GAP, "named gap cites 12 deferred + R=92.5")
    check("2 over-budget" in NAMED_GAP, "named gap cites 2 over-budget")
    check("5 no-bound" in NAMED_GAP, "named gap cites 5 no-bound")
    check("6 rounds" in NAMED_GAP, "named gap cites 6 rounds")
    check(BODY11_BUILT_FLOOR == 0, "Body11 built floor is 0 (not 127)")
    check(BODY11_VERIFY is True, "Body11 converts with verify so volumeDeltaPct is measured")
    check(BODY11_VOL_MAX == 0.023832, "Body11 vol pin is 0.023832 (not widened)")
    check(BODY11_STALE_127 == 127, "stale 127 constant exists only as a non-pin")
    check(
        FALLBACK_GT_RADII[-1] == 100.0 and FALLBACK_GT_RADII[0] == 0.5,
        "GROUND-TRUTH.md fallback radii bookends",
    )

    rch = load_prg_ratchet()
    check(rch["authority"] == RATCHET_NAME, "prg-ratchet authority is the JSON filename")
    check(str(rch.get("engineRef") or "") == "9b693c8", "prg-ratchet engineRef is 9b693c8")
    check(
        abs(float(rch["tightBudgetAbsDeltaMM3Ceiling"]) - 10.81525) < 1e-9,
        "tightBudget ceiling is superseded-shell 10.81525 mm3 (OP2e 10.3094 not ratcheted)",
    )
    check(int(rch["censusRadiiPhantomCeiling"]) == 18, "census.radii phantom ceiling 18 (rid 164)")
    check(int(rch["d3f7R3BuiltCeiling"]) == 75, "d3f7 R=3 built ceiling 75")
    extra = rch.get("d3f7ExtraClassCeilings") or {}
    check(int(extra.get("2.55", 0)) == 2, "d3f7 R=2.55 landing extra ceiling 2")
    check("rid 164" in str(rch.get("reason") or ""), "ratchet reason names rid 164 (D-S3-133)")
    try:
        load_prg_ratchet(Path("/tmp/not-prg-ratchet.json"))
        check(False, "wrong ratchet name must FAIL")
    except GateError:
        check(True, "wrong ratchet name FAIL")

    runtime_gt = pickup_gt_candidate_paths()
    check(runtime_gt[-1] == TEAM_PICKUP_GT, "_team/inputs GT is last runtime candidate")
    runtime_stl = pickup_stl_candidate_paths()
    check(runtime_stl[-1] == TEAM_PICKUP_STL, "_team/inputs STL is last runtime candidate")

    # Isolated missing-GT scenario: tempdir with no GT candidates reachable.
    # Must not consult in-repo diag GT or repo-external _team/ state.
    with tempfile.TemporaryDirectory() as td:
        isolated = Path(td)
        missing_gt = isolated / "no-such-gt.json"
        isolated_team_gt = isolated / "_team" / "inputs" / "ground-truth.json"
        gt_fb, src_fb = load_pickup_gt_radii(
            missing_gt,
            candidates=[missing_gt, isolated_team_gt],
        )
        check("inline" in src_fb, f"missing GT file uses inline fallback (src={src_fb})")
        check(gt_fb == list(FALLBACK_GT_RADII), "fallback radii match GROUND-TRUTH.md")
        check(not missing_gt.exists() and not isolated_team_gt.exists(),
              "isolated GT scenario created no GT files")
        missing_sc = isolated / "no-such-sidecar.json"
        faceted_iso, src_iso = load_must_remain_faceted(
            missing_sc,
            candidates=[missing_sc],
        )
        check("inline" in src_iso, f"isolated missing sidecar uses GT fallback (src={src_iso})")
        iso_types = {str(e.get("type")): int(e.get("count") or 0) for e in faceted_iso}
        check(iso_types.get("sphere") == 9 and iso_types.get("torus") == 6,
              "isolated missing sidecar yields 9 spheres + 6 tori")

    check(match_gt_radius(5.0, FALLBACK_GT_RADII) == 5.0, "GT R5 exact")
    # 0.3% of 100 = 0.3; 100.3 ok, 100.4 not
    check(match_gt_radius(100.3, FALLBACK_GT_RADII) == 100.0, "R100 at +0.3% matches")
    check(match_gt_radius(100.4, FALLBACK_GT_RADII) is None, "R100 at +0.4% is phantom")
    check(match_gt_radius(7.0, FALLBACK_GT_RADII) is None, "R7 is not a pickup GT radius")

    with tempfile.TemporaryDirectory() as td:
        empty_sc = Path(td) / "empty.expected.json"
        empty_sc.write_text('{"mustRemainFaceted": []}\n', encoding="utf-8")
        empty_ents, empty_src = load_must_remain_faceted(
            empty_sc, candidates=[empty_sc]
        )
        check("inline" in empty_src, "empty mustRemainFaceted uses GT fallback")
        empty_types = {str(e.get("type")): int(e.get("count") or 0) for e in empty_ents}
        check(empty_types.get("sphere") == 9 and empty_types.get("torus") == 6,
              "empty sidecar still yields 9 spheres + 6 tori")
        pop_sc = Path(td) / "populated.expected.json"
        pop_sc.write_text(
            json.dumps(
                {
                    "mustRemainFaceted": [
                        {"type": "sphere", "count": 9, "note": "from sidecar"},
                        {"type": "torus", "count": 6, "note": "from sidecar"},
                    ]
                }
            ),
            encoding="utf-8",
        )
        pop_ents, pop_src = load_must_remain_faceted(
            pop_sc, candidates=[pop_sc]
        )
        check(str(pop_sc) == pop_src, "populated sidecar is used when present")
        check(pop_ents[0].get("note") == "from sidecar", "populated sidecar entries win")

    gt_radii, _ = load_pickup_gt_radii()
    pass_checks = evaluate_pickup(
        synthetic_pass_result(),
        synthetic_pass_census(),
        gt_radii,
        0,
        strict=True,
    )
    check(all(c.ok for c in pass_checks), "synthetic RESULT/census satisfying floors PASSES")
    check(failing_names(pass_checks) == [], "synthetic-pass has no failing floors")
    check(
        any(c.name == "census.sphere" and c.ok for c in pass_checks)
        and any(c.name == "census.torus" and c.ok for c in pass_checks),
        "synthetic-pass asserts census.sphere==0 and census.torus==0",
    )

    head_checks = evaluate_pickup(
        synthetic_head_fail_result(),
        synthetic_head_fail_census(),
        gt_radii,
        0,
        strict=True,
    )
    named = failing_names(head_checks)
    check("smoothRevertedComponents" in named, "HEAD-red names smoothRevertedComponents")
    check("smoothBuiltComponents" in named, "HEAD-red names smoothBuiltComponents")
    check("smoothBuiltPlanes" in named, "HEAD-red names smoothBuiltPlanes")
    check("smoothBuiltCylinders" in named, "HEAD-red names smoothBuiltCylinders")
    check("combined" in named, "HEAD-red names combined")
    check("census.plane" in named, "HEAD-red names census.plane")
    check("census.cylinder" in named, "HEAD-red names census.cylinder")
    check("census.sphere" not in named, "HEAD-red does not fail census.sphere (0 at HEAD)")
    check("census.torus" not in named, "HEAD-red does not fail census.torus (0 at HEAD)")
    msgs = " ".join(c.message for c in head_checks if not c.ok)
    check("smoothBuiltPlanes=0 < floor 18" in msgs, "HEAD-red message cites plane floor 18")
    check("smoothBuiltCylinders=0 < floor 30" in msgs, "HEAD-red message cites cylinder floor 30")
    check(any(c.name == "d3f7.radiiMultiset" and c.ok for c in head_checks),
          "HEAD-red D3F-7 vacuously PASSES (0 built <= multiplicity)")

    head_soft = evaluate_pickup(
        synthetic_head_fail_result(),
        synthetic_head_fail_census(),
        gt_radii,
        0,
        strict=False,
    )
    check(
        all((not c.ok) and (not c.hard) for c in head_soft if c.group == "pickup" and not c.ok),
        "STRICT=0 marks pickup floor FAILs as non-hard EXPECTED-RED",
    )
    check(failing_names(head_soft, hard_only=True) == [], "STRICT=0 has no hard pickup fails")

    combo_low = synthetic_pass_result()
    combo_low["smoothBuiltPlanes"] = FLOOR_PLANES - 1
    combo_checks = evaluate_pickup(
        combo_low, synthetic_pass_census(), gt_radii, 0, strict=True
    )
    check("combined" in failing_names(combo_checks), "combined=47 fails when planes=17 (18+30=48)")
    check("smoothBuiltPlanes" in failing_names(combo_checks), "planes=17 also fails plane floor 18")

    phantom_cen = synthetic_pass_census()
    phantom_cen["cylinder_radii"] = list(phantom_cen["cylinder_radii"]) + [7.0]
    ph = evaluate_pickup(
        synthetic_pass_result(), phantom_cen, gt_radii, 0, strict=True
    )
    check("census.radii" not in failing_names(ph), "1 phantom within ratchet 18 PASSES")

    ph17_cen = synthetic_pass_census()
    ph17_cen["cylinder_radii"] = list(ph17_cen["cylinder_radii"]) + [7.0 + 0.01 * i for i in range(17)]
    ph17 = evaluate_pickup(
        synthetic_pass_result(), ph17_cen, gt_radii, 0, strict=True
    )
    check("census.radii" not in failing_names(ph17), "17 phantoms < ratchet 18 PASSES")

    ph18_cen = synthetic_pass_census()
    ph18_cen["cylinder_radii"] = list(ph18_cen["cylinder_radii"]) + [7.0 + 0.01 * i for i in range(18)]
    ph18 = evaluate_pickup(
        synthetic_pass_result(), ph18_cen, gt_radii, 0, strict=True
    )
    check("census.radii" not in failing_names(ph18), "18 phantoms == ratchet 18 PASSES")

    ph19_cen = synthetic_pass_census()
    ph19_cen["cylinder_radii"] = list(ph19_cen["cylinder_radii"]) + [7.0 + 0.01 * i for i in range(19)]
    ph19 = evaluate_pickup(
        synthetic_pass_result(), ph19_cen, gt_radii, 0, strict=True
    )
    check("census.radii" in failing_names(ph19), "19 phantoms > ratchet 18 FAILS")
    check(
        any("19 phantom" in c.message and "18" in c.message for c in ph19 if not c.ok),
        "census.radii FAIL prints current phantom count vs ceiling",
    )

    file_raw, file_src = load_pickup_gt_radii_all()
    file_rec, file_non, file_all = split_cylinder_radii_all(file_raw)
    check(len(file_all) == 54, f"JSON cylinder_radii_all still 54 (got {len(file_all)})")
    check(len(file_rec) == 52, f"JSON recoverable denominator 52 (got {len(file_rec)})")
    check(len(file_non) == 2, f"JSON hex note n=2 (got {len(file_non)})")
    check(
        all(e.get("meshRecoverable") is False for e in file_non)
        and all("hex boss" in str(e.get("reason") or "") for e in file_non),
        "JSON hex entries carry meshRecoverable:false + reason",
    )
    check(
        any(c.name == "gt.notMeshRecoverable" and "n=2" in c.message for c in pass_checks),
        "synthetic-pass prints GT not mesh-recoverable (n=2)",
    )

    hex_cen = synthetic_pass_census()
    hex_cen["cylinder_radii"] = list(hex_cen["cylinder_radii"]) + [HEXNOTE_RADIUS]
    hex_ch = evaluate_pickup(
        synthetic_pass_result(), hex_cen, gt_radii, 0, strict=True
    )
    check("census.radii" not in failing_names(hex_ch), "R=11.544 face is a real match, not a phantom")
    check(
        "d3f7.radiiMultiset" not in failing_names(hex_ch),
        "R=11.544 is excluded from d3f7 denominator (not over-seg)",
    )

    over_cen = synthetic_pass_census()
    over_cen["cylinder_radii"] = list(over_cen["cylinder_radii"]) + [0.5]
    # GT multiplicity of R=0.5 is 1; synthetic-pass already includes one 0.5.
    over = evaluate_pickup(
        synthetic_pass_result(), over_cen, gt_radii, 0, strict=True
    )
    check("d3f7.radiiMultiset" in failing_names(over), "D3F-7 over-segmentation FAILS")
    check(
        any("over-segmentation" in c.message and "R=3 ratchet" in c.message for c in over if not c.ok),
        "D3F-7 message names over-segmentation / R=3 ratchet",
    )

    r3_ok_cen = synthetic_pass_census()
    r3_ok_cen["cylinder_radii"] = [3.0] * 75
    r3_ok_cen["surfaces"] = dict(r3_ok_cen["surfaces"])
    r3_ok_cen["surfaces"]["cylinder"] = 75
    r3_ok = evaluate_pickup(
        synthetic_pass_result(), r3_ok_cen, gt_radii, 0, strict=True
    )
    check("d3f7.radiiMultiset" not in failing_names(r3_ok), "R=3 built=75 == ratchet PASSES")

    r3_bad_cen = synthetic_pass_census()
    r3_bad_cen["cylinder_radii"] = [3.0] * 76
    r3_bad_cen["surfaces"] = dict(r3_bad_cen["surfaces"])
    r3_bad_cen["surfaces"]["cylinder"] = 76
    r3_bad = evaluate_pickup(
        synthetic_pass_result(), r3_bad_cen, gt_radii, 0, strict=True
    )
    check("d3f7.radiiMultiset" in failing_names(r3_bad), "R=3 built=76 > ratchet 75 FAILS")
    check(
        any("built=76" in c.message and "75" in c.message for c in r3_bad if not c.ok),
        "d3f7 FAIL prints current R=3 built vs ceiling",
    )

    r255_ok_cen = synthetic_pass_census()
    r255_ok_cen["cylinder_radii"] = [2.55] * 2 + [3.0] * 28
    r255_ok_cen["surfaces"] = dict(r255_ok_cen["surfaces"])
    r255_ok_cen["surfaces"]["cylinder"] = 30
    r255_ok = evaluate_pickup(
        synthetic_pass_result(), r255_ok_cen, gt_radii, 0, strict=True
    )
    check(
        "d3f7.radiiMultiset" not in failing_names(r255_ok),
        "R=2.55 built=2 == landing extra ceiling PASSES",
    )
    r255_bad_cen = synthetic_pass_census()
    r255_bad_cen["cylinder_radii"] = [2.55] * 3 + [3.0] * 27
    r255_bad_cen["surfaces"] = dict(r255_bad_cen["surfaces"])
    r255_bad_cen["surfaces"]["cylinder"] = 30
    r255_bad = evaluate_pickup(
        synthetic_pass_result(), r255_bad_cen, gt_radii, 0, strict=True
    )
    check(
        "d3f7.radiiMultiset" in failing_names(r255_bad),
        "R=2.55 built=3 > landing extra ceiling 2 FAILS",
    )
    over_exp, over_unexp = split_hard_fails(over)
    check(
        any(c.name == "d3f7.radiiMultiset" for c in over_unexp)
        and not any(c.name == "d3f7.radiiMultiset" for c in over_exp),
        "D3F-7 over-segmentation is UNEXPECTED (not PASS_REGEX inverted)",
    )

    sphere_cen = synthetic_pass_census()
    sphere_cen["surfaces"] = dict(sphere_cen["surfaces"])
    sphere_cen["surfaces"]["sphere"] = 1
    sph = evaluate_pickup(
        synthetic_pass_result(), sphere_cen, gt_radii, 0, strict=True
    )
    check("census.sphere" in failing_names(sph), "analytic sphere fails census.sphere")
    sph_exp, sph_unexp = split_hard_fails(sph)
    check(
        any(c.name == "census.sphere" for c in sph_unexp)
        and not any(c.name == "census.sphere" for c in sph_exp),
        "census.sphere leak is UNEXPECTED (not PASS_REGEX inverted)",
    )

    torus_cen = synthetic_pass_census()
    torus_cen["surfaces"] = dict(torus_cen["surfaces"])
    torus_cen["surfaces"]["torus"] = 1
    tor = evaluate_pickup(
        synthetic_pass_result(), torus_cen, gt_radii, 0, strict=True
    )
    check("census.torus" in failing_names(tor), "analytic torus fails census.torus")

    head_exp, head_unexp = split_hard_fails(head_checks)
    check(
        EXPECTED_RED_FLOOR_NAMES <= {c.name for c in head_exp},
        "HEAD-like fail classifies 18/30/48 family as expected-red",
    )
    check(
        {c.name for c in head_unexp} <= {"census.plane"},
        "synthetic HEAD-red unexpected is only census.plane "
        "(live HEAD is vacuous PASS at 1905; self-test uses plane=0 to prove the floor)",
    )

    pass_exp, pass_unexp = split_hard_fails(pass_checks)
    check(not pass_exp and not pass_unexp, "synthetic-pass has no hard fails to invert")

    tb_bad = synthetic_pass_result()
    tb_bad["stepVolumeMM3"] = 200000.0
    tb_bad["volumeDeltaPct"] = 50.0
    tb = evaluate_pickup(tb_bad, synthetic_pass_census(), gt_radii, 0, strict=True)
    check("tightBudget" in failing_names(tb), "tightBudget floor is named on volume miss")

    tb_at = synthetic_pass_result()
    tb_at["meshVolumeMM3"] = 100000.0
    tb_at["stepVolumeMM3"] = 100000.0 + 10.81525
    tb_at["volumeDeltaPct"] = 0.01
    tb_at_ch = evaluate_pickup(tb_at, synthetic_pass_census(), gt_radii, 0, strict=True)
    check("tightBudget" not in failing_names(tb_at_ch), "tightBudget |dV|=10.81525 == ratchet PASSES")

    tb_over = synthetic_pass_result()
    tb_over["meshVolumeMM3"] = 100000.0
    tb_over["stepVolumeMM3"] = 100000.0 + 10.816
    tb_over["volumeDeltaPct"] = 0.01
    tb_over_ch = evaluate_pickup(tb_over, synthetic_pass_census(), gt_radii, 0, strict=True)
    check("tightBudget" in failing_names(tb_over_ch), "tightBudget |dV|=10.816 > ratchet FAILS")
    check(
        any("10.816" in c.message and "ratchet" in c.message for c in tb_over_ch if not c.ok),
        "tightBudget FAIL prints current |dV| vs ceiling",
    )

    tb_unmeasured = synthetic_pass_result()
    tb_unmeasured["volumeDeltaPct"] = -1
    tb2 = evaluate_pickup(
        tb_unmeasured, synthetic_pass_census(), gt_radii, 0, strict=True
    )
    check("tightBudget" in failing_names(tb2), "tightBudget fails when volume not measured")

    hl_ok = evaluate_handle_lock(
        {"ok": True},
        {"surfaces": {"cylinder": 15, "plane": 13}},
        0,
    )
    check(all(c.ok for c in hl_ok), "handle-lock census.cylinder>=15 PASSES")
    hl_bad = evaluate_handle_lock(
        {"ok": True},
        {"surfaces": {"cylinder": 14, "plane": 13}},
        0,
    )
    check("census.cylinder" in failing_names(hl_bad), "handle-lock 14 cylinders FAILS")

    b11_zero = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": 0.01,
            "solids": 2,
            "smoothBuiltCylinders": 0,
        },
        2,
    )
    check(all(c.ok for c in b11_zero), "Body11 built=0 (live HEAD) PASSES never-get-worse")
    b11_127 = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": 0.01,
            "solids": 2,
            "smoothBuiltCylinders": BODY11_STALE_127,
        },
        2,
    )
    check(all(c.ok for c in b11_127), "Body11 built=127 also PASSES (floor is 0, 127 not pinned)")
    b11_vol = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": 0.03,
            "solids": 2,
            "smoothBuiltCylinders": 0,
        },
        2,
    )
    check("volumeDeltaPct" in failing_names(b11_vol), "Body11 volΔ>0.023832 FAILS")
    b11_solids = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": 0.0,
            "solids": 1,
            "smoothBuiltCylinders": 0,
        },
        2,
    )
    check("solids" in failing_names(b11_solids), "Body11 solids!=2 FAILS")
    b11_skip_vol = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": -1,
            "solids": 2,
            "smoothBuiltCylinders": 0,
        },
        2,
    )
    check(
        all(c.ok for c in b11_skip_vol),
        "Body11 built=0 + unmeasured volume PASSES (faceted path)",
    )
    b11_need_vol = evaluate_body11(
        {
            "ok": True,
            "volumeDeltaPct": -1,
            "solids": 2,
            "smoothBuiltCylinders": 5,
        },
        2,
    )
    check(
        "volumeDeltaPct" in failing_names(b11_need_vol),
        "Body11 built>0 requires a measured volumeDeltaPct",
    )
    # Equality-to-127 is NOT encoded: a built count of 0 must not mention a miss vs 127.
    zero_msgs = " ".join(c.message for c in b11_zero)
    check("!= 127" not in zero_msgs and "< 127" not in zero_msgs, "Body11 does not pin 127")

    check(evaluate_body18(18)[0].ok, "Body18 freeEdges=18 PASSES")
    check(not evaluate_body18(19)[0].ok, "Body18 freeEdges=19 FAILS")
    check("freeEdges" in failing_names(evaluate_body18(19)), "Body18 names freeEdges")

    check(strict_enabled({}) is True, "STRICT default ON when env empty")
    check(strict_enabled({ENV_STRICT: "1"}) is True, "STRICT=1 is ON")
    check(strict_enabled({ENV_STRICT: "0"}) is False, "STRICT=0 is OFF")
    check(strict_enabled({ENV_STRICT: "off"}) is False, "STRICT=off is OFF")

    check(
        parse_result_line('progress\nRESULT {"ok":true,"openShells":0}\n')["ok"] is True,
        "parse_result_line",
    )
    check(
        parse_j6_free_edges(
            ["J6: shell not closed freeEdges=36 faces=957 recover=0"]
        )
        == 36,
        "parse_j6_free_edges",
    )

    try:
        stl = resolve_pickup_stl()
        check(stl.is_file(), f"pickup STL fallback resolves ({stl})")
    except GateError as exc:
        check(False, f"pickup STL fallback: {exc}")

    check(CENSUS_PY.is_file(), f"step_census.py present ({CENSUS_PY})")
    check(HL_STL.is_file(), "handle-lock.stl present")
    check(BODY11_STL.is_file(), "Body11.stl present")

    rc = run_synthetic_pass()
    check(rc == 0, "--synthetic-pass returns 0")

    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path, help="stl2step CLI (required unless --self-test/--synthetic-pass)")
    p.add_argument("--pickup-stl", type=Path, help="Override handle-pickup STL")
    p.add_argument("--handle-lock-stl", type=Path, default=HL_STL)
    p.add_argument("--body11-stl", type=Path, default=BODY11_STL)
    p.add_argument("--ground-truth", type=Path, help="diag ground-truth.json (optional)")
    p.add_argument(
        "--private-corpus",
        type=Path,
        help=f"Private STL dir for Body18 (else ${ENV_PRIVATE}; missing SKIPs Body18 only)",
    )
    p.add_argument("--jobs", type=int, default=4, help="Parallel conversions (default 4)")
    p.add_argument(
        "--ratchet",
        type=Path,
        default=DEFAULT_RATCHET,
        help="PRG ratchet JSON (tests/gates/baseline/prg-ratchet.json)",
    )
    p.add_argument("--self-test", action="store_true", help="Exercise gate API (no engine)")
    p.add_argument(
        "--synthetic-pass",
        action="store_true",
        help="Evaluate a RESULT/census fixture that satisfies the floors",
    )
    g = p.add_mutually_exclusive_group()
    g.add_argument(
        "--strict",
        dest="strict",
        action="store_true",
        default=None,
        help="Pickup floors are hard (default; also PARTIAL_RECOVERY_STRICT=1)",
    )
    g.add_argument(
        "--no-strict",
        dest="strict",
        action="store_false",
        help="Pickup floors are EXPECTED-RED (PARTIAL_RECOVERY_STRICT=0)",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.synthetic_pass:
            return run_synthetic_pass()
        if not args.binary:
            print("partial_recovery_gate: --binary is required", file=sys.stderr)
            return 1
        return run_live(
            Path(args.binary),
            pickup_stl=args.pickup_stl,
            hl_stl=Path(args.handle_lock_stl),
            body11_stl=Path(args.body11_stl),
            private_corpus=args.private_corpus,
            jobs=args.jobs,
            gt_path=args.ground_truth,
            strict=args.strict,
            ratchet_path=Path(args.ratchet) if args.ratchet else DEFAULT_RATCHET,
        )
    except GateError as exc:
        print(f"partial_recovery_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
