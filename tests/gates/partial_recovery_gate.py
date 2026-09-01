#!/usr/bin/env python3
"""partial_recovery_gate — Handle-pickup DECISION-boundary floors + non-regression.

ctest name: partial_recovery_gate (owns the G0.1 *lane*; G0.1 itself stays in
run_gates.py and is not rewritten here).

Pickup floors are DECISION-boundary A1/A2 verbatim. At current HEAD they FAIL
(EXPECTED-RED until a later phase ships analytic faces). PARTIAL_RECOVERY_STRICT
defaults ON so the Python process is red-at-HEAD with a message that names the
failing floor. PARTIAL_RECOVERY_STRICT=0 keeps pickup floors documented but
non-fatal (CLI / debug only).

CMake EXPECTED-RED (FINDINGS-0 GAP2 / ci-local-gate.sh):
  grep of scripts/ci-local-gate.sh + tests/gates found no WILL_FAIL property
  (run_gates.py XFAIL is "phase not landed" inside gates_full; j6 / stress_sweep
  use SKIP_RETURN_CODE 77; census/check_error.cmake expects exit 1 via a
  wrapper; harness_dump_cube uses PASS_REGULAR_EXPRESSION).

  WILL_FAIL TRUE alone would let `ctest` pass while Python still exits 1, but
  it also inverts unexpected fails (handle-lock, census.sphere leak) into
  greens — measured. PASS_REGULAR_EXPRESSION "PARTIAL_RECOVERY_EXPECTED_RED"
  is the CTest invert that keeps ci-local-gate.sh green on unmet 18/50/100
  without swallowing GAP4 / non-regression.

FLIP PROTOCOL (the day 18/50/100 are met — FAIL LOUDLY if unmarked):
  1. Live Python prints PARTIAL_RECOVERY_FLOORS_MET and exits 0.
  2. ctest FAILS: Required regular expression not found
     Regex=[PARTIAL_RECOVERY_EXPECTED_RED].
  3. Remove PASS_REGULAR_EXPRESSION from tests/gates/CMakeLists.txt
     (partial_recovery_gate PROPERTIES).
  4. Keep ENVIRONMENT PARTIAL_RECOVERY_STRICT=1. Do not lower FLOOR_*.

# TODO(PHASE-CLOSE): tighten FLOOR_PLANES / FLOOR_CYLINDERS / FLOOR_COMBINED
# from the landed contain-race winner. Knob: the three constants below in this
# file. Do not raise them in the gates-partial lane.

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
MAIN_PICKUP_STL = Path(
    "/Users/jroberts/Desktop/Internal Development/Tools/stl2step"
    "/_team/inputs/Handle pickup.stl"
)

ENV_STRICT = "PARTIAL_RECOVERY_STRICT"
ENV_PRIVATE = "STL2STEP_PRIVATE_CORPUS"

# DECISION-boundary A2 — provisional floors. PHASE-CLOSE knob is these three.
FLOOR_PLANES = 18
FLOOR_CYLINDERS = 50
FLOOR_COMBINED = 100
FLOOR_REVERTED = 0
FLOOR_BUILT_COMPONENTS = 1
FLOOR_OPEN_SHELLS = 0
ALLOWED_EXITS = (0, 2)
GT_REL_TOL = 0.003  # 0.3%

# Pickup floors honestly unmet at HEAD 7cf77a2 (FINDINGS-0 GAP2). CTest
# PASS_REGULAR_EXPRESSION inverts only these (see module docstring).
# census.sphere / census.torus are NOT in this set — they pass at HEAD
# (analytic count 0) and must stay hard so a sphere/torus leak is not
# inverted into a green.
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

HL_CENSUS_CYL_FLOOR = 15
BODY11_VOL_MAX = 0.023832
BODY11_SOLIDS = 2
BODY11_BUILT_FLOOR = 0  # never pin the stale 127 figure
BODY11_STALE_127 = 127  # named only so self-test can prove it is NOT a pin
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


def resolve_pickup_stl(extra: Optional[Path] = None) -> Path:
    candidates: List[Path] = []
    if extra is not None:
        candidates.append(Path(extra))
    candidates.extend([PICKUP_CORPUS, PICKUP_DIAG, MAIN_PICKUP_STL])
    for p in candidates:
        if p.is_file():
            return p
    raise GateError(
        "handle-pickup STL missing; looked for "
        + ", ".join(str(p) for p in candidates)
    )


def load_pickup_gt_radii(gt_path: Optional[Path] = None) -> Tuple[List[float], str]:
    """Read GT radii from diag ground-truth.json; else GROUND-TRUTH.md inline."""
    candidates: List[Path] = []
    if gt_path is not None:
        candidates.append(Path(gt_path))
    candidates.append(PICKUP_GT)
    for p in candidates:
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


def load_must_remain_faceted(
    sidecar_path: Optional[Path] = None,
) -> Tuple[List[Dict[str, Any]], str]:
    """Sidecar mustRemainFaceted, or GROUND-TRUTH.md 9 spheres + 6 tori.

    Empty list is unpopulated (FINDINGS-0 GAP4), not "nothing stays faceted".
    """
    candidates: List[Path] = []
    if sidecar_path is not None:
        candidates.append(Path(sidecar_path))
    candidates.append(PICKUP_SIDECAR)
    for p in candidates:
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
) -> Tuple[bool, str, Dict[str, Any]]:
    """S5 / host D4.5 proxy: |V_step-V_mesh| vs max(1e-4*|V_mesh|, 3*|dVol|)."""
    vd = float(result.get("volumeDeltaPct", -1))
    mesh = float(result.get("meshVolumeMM3") or 0.0)
    step = float(result.get("stepVolumeMM3") or 0.0)
    pred = float(result.get("smoothVolPredictedMM3") or 0.0)
    details = {
        "volumeDeltaPct": vd,
        "meshVolumeMM3": mesh,
        "stepVolumeMM3": step,
        "smoothVolPredictedMM3": pred,
    }
    if vd < 0:
        return False, "tightBudget: volume not measured (need verify on)", details
    mesh_abs = abs(mesh)
    budget = max(TIGHT_MESH_K * mesh_abs, TIGHT_DVOL_K * abs(pred))
    abs_delta = abs(step - mesh)
    details["absDeltaMM3"] = abs_delta
    details["budgetMM3"] = budget
    if abs_delta > budget + 1e-9:
        return (
            False,
            f"tightBudget: |V_step-V_mesh|={abs_delta:.6g} mm^3 > budget {budget:.6g} "
            f"(max({TIGHT_MESH_K}*|V_mesh|, {TIGHT_DVOL_K}*|dVolPredicted|))",
            details,
        )
    return True, f"tightBudget: |V_step-V_mesh|={abs_delta:.6g} <= {budget:.6g}", details


def evaluate_pickup(
    result: Mapping[str, Any],
    census: Mapping[str, Any],
    gt_radii: Sequence[float],
    exit_code: int,
    *,
    strict: bool = True,
) -> List[Check]:
    """DECISION-boundary pickup floors. Each FAIL message names the floor."""
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

    tb_ok, tb_msg, tb_det = tight_budget(result)
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
    add(
        "census.radii",
        not phantoms,
        (
            f"census.radii: {len(phantoms)} radius(es) outside GT 0.3%: {phantoms[:8]}"
            if phantoms
            else f"census.radii: {len(radii)} all within 0.3% of GT"
        ),
        phantoms=phantoms[:12],
        nRadii=len(radii),
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
        # Faceted ship is volume-exact; skip the OCCT re-read so the gate
        # stays inside the 2-minute add. Re-enable when built > 0.
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
    radii: List[float] = []
    gt = list(FALLBACK_GT_RADII)
    while len(radii) < FLOOR_CYLINDERS:
        radii.append(gt[len(radii) % len(gt)])
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
) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    pickup = resolve_pickup_stl(pickup_stl)
    gt_radii, gt_src = load_pickup_gt_radii(gt_path)
    faceted, faceted_src = load_must_remain_faceted()
    is_strict = strict_enabled() if strict is None else bool(strict)

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
            "verify": False,
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
            f"faceted={faceted_src} (sphere={n_sph} torus={n_tor})  "
            f"strict={is_strict}  jobs={workers} threads/file={per_file}",
            flush=True,
        )
        print(format_checks(checks), flush=True)
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
        print(
            "FLIP PROTOCOL: remove PASS_REGULAR_EXPRESSION from "
            "tests/gates/CMakeLists.txt (partial_recovery_gate). "
            "Keep PARTIAL_RECOVERY_STRICT=1. Do not lower FLOOR_*.",
            flush=True,
        )
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
    check(FLOOR_CYLINDERS == 50, f"FLOOR_CYLINDERS is DECISION 50 (got {FLOOR_CYLINDERS})")
    check(FLOOR_COMBINED == 100, f"FLOOR_COMBINED is DECISION 100 (got {FLOOR_COMBINED})")
    check(BODY11_BUILT_FLOOR == 0, "Body11 built floor is 0 (not 127)")
    check(BODY11_STALE_127 == 127, "stale 127 constant exists only as a non-pin")
    check(
        FALLBACK_GT_RADII[-1] == 100.0 and FALLBACK_GT_RADII[0] == 0.5,
        "GROUND-TRUTH.md fallback radii bookends",
    )

    gt_fb, src_fb = load_pickup_gt_radii(Path("/no/such/handle-pickup-gt.json"))
    check("inline" in src_fb, f"missing GT file uses inline fallback (src={src_fb})")
    check(gt_fb == list(FALLBACK_GT_RADII), "fallback radii match GROUND-TRUTH.md")

    check(match_gt_radius(5.0, FALLBACK_GT_RADII) == 5.0, "GT R5 exact")
    # 0.3% of 100 = 0.3; 100.3 ok, 100.4 not
    check(match_gt_radius(100.3, FALLBACK_GT_RADII) == 100.0, "R100 at +0.3% matches")
    check(match_gt_radius(100.4, FALLBACK_GT_RADII) is None, "R100 at +0.4% is phantom")
    check(match_gt_radius(7.0, FALLBACK_GT_RADII) is None, "R7 is not a pickup GT radius")

    faceted_fb, src_faceted = load_must_remain_faceted(Path("/no/such-sidecar.json"))
    check("inline" in src_faceted, f"missing sidecar uses GT fallback (src={src_faceted})")
    by_type = {str(e.get("type")): int(e.get("count") or 0) for e in faceted_fb}
    check(by_type.get("sphere") == 9, f"GT fallback 9 spheres (got {by_type})")
    check(by_type.get("torus") == 6, f"GT fallback 6 tori (got {by_type})")

    with tempfile.TemporaryDirectory() as td:
        empty_sc = Path(td) / "empty.expected.json"
        empty_sc.write_text('{"mustRemainFaceted": []}\n', encoding="utf-8")
        empty_ents, empty_src = load_must_remain_faceted(empty_sc)
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
        pop_ents, pop_src = load_must_remain_faceted(pop_sc)
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
    check("smoothBuiltCylinders=0 < floor 50" in msgs, "HEAD-red message cites cylinder floor 50")

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
    combo_low["smoothBuiltCylinders"] = FLOOR_COMBINED - FLOOR_PLANES - 1
    combo_checks = evaluate_pickup(
        combo_low, synthetic_pass_census(), gt_radii, 0, strict=True
    )
    check("combined" in failing_names(combo_checks), "combined floor is independent of 18/50")

    phantom_cen = synthetic_pass_census()
    phantom_cen["cylinder_radii"] = list(phantom_cen["cylinder_radii"]) + [7.0]
    ph = evaluate_pickup(
        synthetic_pass_result(), phantom_cen, gt_radii, 0, strict=True
    )
    check("census.radii" in failing_names(ph), "phantom radius fails census.radii")

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
        "HEAD-like fail classifies 18/50/100 family as expected-red",
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
        )
    except GateError as exc:
        print(f"partial_recovery_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
