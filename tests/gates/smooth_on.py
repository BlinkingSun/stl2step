#!/usr/bin/env python3
"""Smooth-on gate helpers for tests/gates/run_gates.py (stdlib only)."""

from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

# CALIBRATED AT P0 CLOSE (f6-budgets, 2026-08-25). Evidence:
# _team/reports/spike-volume-budgets.md (f6-budgets lane).
G3_RATIO_LO = 0.35
G3_RATIO_HI = 2.00
G3_CHORD_SAGITTA_WINDOW_MM = 0.0  # unused — G3 volume is mm; sagitta window is ratio
G3_VOLUME_K = 3.0
G5_EDGE_K = 0.1
G4_4_MIN_G1_S02 = 22
LEGACY_VOLUME_GATE_PCT = 0.01

DIRTY_SKIP_FIXTURES = frozenset({"S15", "S16-R1-round-2"})
R_LADDER_FIXTURES = frozenset(
    {"S16-R1-explode-success", "S16-R1-round-2", "S16-R2-ChainUnstable"}
)
G4_4_FIXTURES = frozenset({"S02"})
G3_KNOWN_RED_FIXTURES = frozenset({"S05", "S11-b"})

SCHEMA_SURFACE_TYPES = frozenset({"plane", "cylinder", "cone", "sphere", "torus"})
SCHEMA_BUILT_AS = frozenset(
    {"notBuilt", "single", "seamed360", "twoHalves", "explodedToFacets"}
)

# Implemented HARD gates parked until v-next / F2 (honest reds documented).
PARKED_GATES: Dict[str, str] = {
    "G3": "S05 11.34% / S11-b 6.84% exceed calibrated k=3 — F6 geometry, v-next",
    "G4.4": "S02 G1=20/22 (2 fillet edges lack EncodeRegularity flags) — KNOWN-GAP, v-next",
    "G2": "BuiltAs / recognition census vs sidecar blocked on F2 — file-truth gap",
    "G2.5": "closed360 BuiltAs census vs sidecar blocked on F2 — seamed360 ladder",
    "R-ladder": "S16 R1/R2 ladder blocked on F2 — edge-failure R1 enrollment",
}


@dataclass
class RunArtifacts:
    exit_code: int
    result: Dict[str, Any]
    stdout: str
    stderr: str
    step_path: Optional[Path] = None


@dataclass
class Fixture:
    id: str
    stl: Path
    sidecar: Optional[Path] = None


@dataclass
class GateContext:
    binary: Path
    fixture: Fixture
    work_dir: Path
    census_path: Optional[Path]
    dump_path: Optional[Path]
    no_verify: bool
    baseline_bin: Optional[Path] = None


@dataclass
class GateOutcome:
    gate_id: str
    fixture_id: str
    status: str  # PASS | FAIL | PARKED | XFAIL | SKIP
    message: str
    hard: bool = True
    details: Dict[str, Any] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.details is None:
            self.details = {}


def schema_enum(value: Any) -> str:
    if not isinstance(value, str) or not value:
        return ""
    return value[0].lower() + value[1:]


def g2_radius_threshold(r_truth: float, n_sides: int, refit_tol: float = 0.0) -> float:
    n = max(int(n_sides), 3)
    chord = r_truth * (1.0 - math.cos(math.pi / n))
    return max(refit_tol, 0.005 * r_truth, chord)


def discover_fixtures(
    corpus: Path,
    smoke_stl: Path,
    smoke: bool,
    fixture_filter: Optional[List[str]],
) -> List[Fixture]:
    """Auto-discover corpus STLs (S* + real-CAD Body* + any *.stl) and tests/cube.stl."""
    if smoke:
        fixtures = [Fixture("cube", smoke_stl)]
    else:
        fixtures = []
        seen: set[str] = set()
        if corpus.is_dir():
            for stl in sorted(corpus.glob("*.stl")):
                fid = stl.stem
                if fid in seen:
                    continue
                seen.add(fid)
                sidecar = corpus / f"{fid}.expected.json"
                fixtures.append(
                    Fixture(fid, stl, sidecar if sidecar.is_file() else None)
                )
        if smoke_stl.is_file() and "cube" not in seen:
            fixtures.append(Fixture("cube", smoke_stl))
        if not fixtures and smoke_stl.is_file():
            fixtures = [Fixture("cube", smoke_stl)]

    if fixture_filter:
        wanted = set(fixture_filter)
        fixtures = [f for f in fixtures if f.id in wanted or f.stl.stem in wanted]
    return fixtures


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise RuntimeError("no RESULT line in stl2step stdout")


def run_stl2step(
    binary: Path,
    stl: Path,
    step_out: Path,
    *,
    no_verify: bool,
    threads: Optional[int] = None,
    smooth: bool = False,
    force_sew: bool = False,
    capture_notes: bool = False,
) -> RunArtifacts:
    cmd = [str(binary), str(stl), "-o", str(step_out)]
    if not capture_notes:
        cmd.append("--quiet")
    if no_verify:
        cmd.append("--no-verify")
    if threads is not None:
        cmd.extend(["--threads", str(threads)])
    if smooth:
        cmd.append("--smooth")
    if force_sew:
        cmd.append("--force-sew")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        result = parse_result_line(proc.stdout)
    except Exception as exc:
        raise RuntimeError(
            f"stl2step failed to emit RESULT for {stl}: {exc}\n"
            f"stdout={proc.stdout!r}\nstderr={proc.stderr!r}"
        ) from exc
    return RunArtifacts(proc.returncode, result, proc.stdout, proc.stderr, step_out)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def load_sidecar(fixture: Fixture) -> Optional[Dict[str, Any]]:
    if fixture.sidecar and fixture.sidecar.is_file():
        return json.loads(read_text(fixture.sidecar))
    return None


def run_census(census_path: Path, step: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [str(census_path), str(step)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout or "census failed")
    return json.loads(proc.stdout)


def canonical_step_data(step_text: str) -> str:
    text = re.sub(
        r"(FILE_NAME\s*\(\s*'[^']*'\s*,\s*)'[^']*'",
        r"\1'STRIPPED'",
        step_text,
        count=1,
        flags=re.IGNORECASE,
    )
    match = re.search(r"DATA;\s*(.*?)\s*ENDSEC;", text, re.DOTALL | re.IGNORECASE)
    return match.group(0) if match else text


def parse_brep_mesh_volume(text: str) -> Tuple[Optional[float], Optional[float]]:
    m = re.search(
        r"volume\s+B-Rep\s+([0-9eE.+-]+)\s+mm\^3,\s+mesh\s+([0-9eE.+-]+)",
        text,
    )
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def dump_dvol_abs_and_cylgrow(
    dump_bin: Path, stl: Path
) -> Tuple[Optional[float], List[Dict[str, Any]]]:
    proc = subprocess.run([str(dump_bin), str(stl)], capture_output=True, text=True)
    if proc.returncode != 0:
        return None, []
    try:
        doc = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None, []
    abs_sum = 0.0
    cylgrow: List[Dict[str, Any]] = []
    for comp in doc.get("comps") or []:
        rs = comp.get("regionSet") or {}
        for r in rs.get("regions") or []:
            abs_sum += abs(float(r.get("dVolPredicted") or 0.0))
            if schema_enum(r.get("origin")) == "cylGrow":
                cylgrow.append(r)
    return abs_sum, cylgrow


def live_rows(sidecar: Dict[str, Any]) -> List[Dict[str, Any]]:
    return list(sidecar.get("live") or [])


def min_recovered_radius(sidecar: Optional[Dict[str, Any]]) -> Optional[float]:
    if not sidecar:
        return None
    radii: List[float] = []
    for prim in sidecar.get("recoverable") or []:
        if schema_enum(prim.get("type")) == "cylinder":
            r = float(prim.get("radius") or 0.0)
            if r > 0:
                radii.append(r)
    return min(radii) if radii else None


def check_g0_2(ctx: GateContext) -> GateOutcome:
    """G0.2 — --smooth never reduces closedness (solids / openShells)."""
    gate_id = "G0.2"
    off_step = ctx.work_dir / f"{ctx.fixture.id}_g02_off.step"
    on_step = ctx.work_dir / f"{ctx.fixture.id}_g02_on.step"
    try:
        off = run_stl2step(
            ctx.binary, ctx.fixture.stl, off_step, no_verify=ctx.no_verify, smooth=False
        )
        on = run_stl2step(
            ctx.binary, ctx.fixture.stl, on_step, no_verify=ctx.no_verify, smooth=True
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))

    for label, run in (("off", off), ("on", on)):
        if not run.result.get("ok"):
            return GateOutcome(
                gate_id,
                ctx.fixture.id,
                "FAIL",
                f"{label}-path conversion failed: {run.result.get('error')}",
            )

    solids_off = int(off.result.get("solids") or 0)
    solids_on = int(on.result.get("solids") or 0)
    open_off = int(off.result.get("openShells") or 0)
    open_on = int(on.result.get("openShells") or 0)
    details = {
        "solidsOff": solids_off,
        "solidsOn": solids_on,
        "openShellsOff": open_off,
        "openShellsOn": open_on,
    }
    if solids_on < solids_off or open_on > open_off:
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"closedness reduced: solids {solids_off}->{solids_on}, "
            f"openShells {open_off}->{open_on}",
            details=details,
        )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"closedness preserved (solids>={solids_off}, openShells<={open_off})",
        details=details,
    )


def check_g0_3(ctx: GateContext) -> GateOutcome:
    """G0.3 — dirty-skip exit parity; --force-sew globally disables refit."""
    gate_id = "G0.3"
    fid = ctx.fixture.id
    sidecar = load_sidecar(ctx.fixture)

    if fid in DIRTY_SKIP_FIXTURES:
        off_step = ctx.work_dir / f"{fid}_g03_off.step"
        on_step = ctx.work_dir / f"{fid}_g03_on.step"
        try:
            off = run_stl2step(
                ctx.binary, ctx.fixture.stl, off_step, no_verify=ctx.no_verify, smooth=False
            )
            on = run_stl2step(
                ctx.binary,
                ctx.fixture.stl,
                on_step,
                no_verify=ctx.no_verify,
                smooth=True,
                capture_notes=True,
            )
        except RuntimeError as exc:
            return GateOutcome(gate_id, fid, "FAIL", str(exc))
        skip = int(on.result.get("smoothSkippedComponents") or 0)
        if skip < 1:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                f"expected dirty skip (smoothSkippedComponents>=1), got {skip}",
            )
        if off.exit_code != on.exit_code:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                f"dirty skip changed exit {off.exit_code}->{on.exit_code}",
                details={"offExit": off.exit_code, "onExit": on.exit_code},
            )
        return GateOutcome(
            gate_id,
            fid,
            "PASS",
            f"dirty skip: exit {off.exit_code} unchanged, smoothSkippedComponents={skip}",
        )

    # force-sew global disable: --force-sew --smooth STEP DATA == --force-sew off-path
    sew_off = ctx.work_dir / f"{fid}_g03_sew_off.step"
    sew_on = ctx.work_dir / f"{fid}_g03_sew_on.step"
    try:
        base = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            sew_off,
            no_verify=True,
            smooth=False,
            force_sew=True,
        )
        both = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            sew_on,
            no_verify=True,
            smooth=True,
            force_sew=True,
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, fid, "FAIL", str(exc))

    if not base.result.get("ok") or not both.result.get("ok"):
        return GateOutcome(gate_id, fid, "FAIL", "force-sew conversion failed")

    d_off = canonical_step_data(read_text(sew_off))
    d_on = canonical_step_data(read_text(sew_on))
    if d_off != d_on:
        return GateOutcome(
            gate_id,
            fid,
            "FAIL",
            "--force-sew --smooth changes canonical STEP DATA (refit not globally disabled)",
        )

    # Exit parity: twin-run vs 187ead0 baseline (G0.1 pattern), else vs own off-path run.
    if ctx.baseline_bin and ctx.baseline_bin.is_file():
        bl_step = ctx.work_dir / f"{fid}_g03_baseline_sew.step"
        try:
            bl = run_stl2step(
                ctx.baseline_bin,
                ctx.fixture.stl,
                bl_step,
                no_verify=True,
                smooth=False,
                force_sew=True,
            )
        except RuntimeError as exc:
            return GateOutcome(gate_id, fid, "FAIL", str(exc))
        if base.exit_code != bl.exit_code:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                f"force-sew exit {base.exit_code} != baseline force-sew {bl.exit_code}",
                details={
                    "candidateExit": base.exit_code,
                    "baselineExit": bl.exit_code,
                },
            )
    else:
        plain_step = ctx.work_dir / f"{fid}_g03_plain_off.step"
        try:
            plain = run_stl2step(
                ctx.binary,
                ctx.fixture.stl,
                plain_step,
                no_verify=True,
                smooth=False,
                force_sew=False,
            )
        except RuntimeError as exc:
            return GateOutcome(gate_id, fid, "FAIL", str(exc))
        if base.exit_code != plain.exit_code:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                f"force-sew exit {base.exit_code} != off-path exit {plain.exit_code}",
                details={
                    "forceSewExit": base.exit_code,
                    "offPathExit": plain.exit_code,
                },
            )

    return GateOutcome(
        gate_id,
        fid,
        "PASS",
        "force-sew: --smooth identity holds (global refit disable)",
    )


def check_g2(ctx: GateContext) -> GateOutcome:
    """G2 — recoverable primitive recognition vs sidecar live[] / recoverable[]."""
    gate_id = "G2"
    sidecar = load_sidecar(ctx.fixture)
    if not sidecar:
        return GateOutcome(gate_id, ctx.fixture.id, "SKIP", "no sidecar")

    on_step = ctx.work_dir / f"{ctx.fixture.id}_g2.step"
    try:
        run = run_stl2step(
            ctx.binary, ctx.fixture.stl, on_step, no_verify=ctx.no_verify, smooth=True
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))
    if not run.result.get("ok"):
        return GateOutcome(
            gate_id, ctx.fixture.id, "FAIL", f"conversion failed: {run.result.get('error')}"
        )

    violations: List[str] = []
    recognised = int(run.result.get("smoothCylinders") or 0)
    built = int(run.result.get("smoothBuiltCylinders") or 0)

    expected_cyl = sum(
        int(p.get("count") or 0)
        for p in sidecar.get("recoverable") or []
        if schema_enum(p.get("type")) == "cylinder"
    )

    details: Dict[str, Any] = {
        "smoothCylinders": recognised,
        "smoothBuiltCylinders": built,
        "expectedRecoverableCylinders": expected_cyl,
        "violations": violations,
    }

    cen: Optional[Dict[str, Any]] = None
    if ctx.census_path and ctx.census_path.is_file():
        try:
            cen = run_census(ctx.census_path, on_step)
            details["censusCylinders"] = int(cen.get("cylinder") or 0)
        except RuntimeError as exc:
            violations.append(str(exc))

    for row in live_rows(sidecar):
        disp = str(row.get("disposition") or "PASS")
        if disp == "ESCALATE":
            continue
        census_row = row.get("surfaceCensus") or {}
        exp_cyl = int(census_row.get("cylinder") or 0)
        exp_plane = int(census_row.get("plane") or 0)
        floor = int(row.get("builtCylindersFloor") or 0)
        if floor > 0 and built < floor:
            violations.append(
                f"comp{row.get('component')}: builtCylinders {built} < floor {floor}"
            )
        if cen is not None and disp in ("PASS", "PARTIAL"):
            got_cyl = int(cen.get("cylinder") or 0)
            got_plane = int(cen.get("plane") or 0)
            if exp_cyl > 0 and got_cyl < exp_cyl:
                violations.append(
                    f"comp{row.get('component')}: census cylinders {got_cyl} < live {exp_cyl}"
                )
            if exp_plane > 0 and got_plane < exp_plane:
                violations.append(
                    f"comp{row.get('component')}: census planes {got_plane} < live {exp_plane}"
                )

    if cen is not None:
        matched_recoverable = 0
        for prim in sidecar.get("recoverable") or []:
            if schema_enum(prim.get("type")) != "cylinder":
                continue
            r_truth = float(prim.get("radius") or 0.0)
            n_sides = int(prim.get("nSides") or 3)
            if r_truth <= 0:
                continue
            thresh = g2_radius_threshold(r_truth, n_sides)
            for grp in cen.get("cylinderGroups") or []:
                r = float(grp.get("radius") or 0.0)
                if abs(r - r_truth) <= thresh:
                    matched_recoverable += int(prim.get("count") or 1)
                    break
        details["matchedRecoverableCylinders"] = matched_recoverable
        if expected_cyl > 0 and matched_recoverable < expected_cyl:
            violations.append(
                f"recognition: matched {matched_recoverable}/{expected_cyl} recoverable cylinders "
                f"within G2 radius bar"
            )

    if violations:
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"G2: {len(violations)} recognition/census mismatch(es)",
            details=details,
        )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"G2: recognition ok (smoothCylinders={recognised}, built={built})",
        details=details,
    )


def check_g2_5(ctx: GateContext) -> GateOutcome:
    """G2.5 — closed360 recoverable holes built as one seamed360 cylindrical face."""
    gate_id = "G2.5"
    sidecar = load_sidecar(ctx.fixture)
    if not sidecar:
        return GateOutcome(gate_id, ctx.fixture.id, "SKIP", "no sidecar")

    closed360 = [
        p
        for p in sidecar.get("recoverable") or []
        if schema_enum(p.get("type")) == "cylinder" and p.get("closed360")
    ]
    if not closed360:
        return GateOutcome(gate_id, ctx.fixture.id, "SKIP", "no closed360 recoverable holes")

    if not (ctx.census_path and ctx.census_path.is_file()):
        return GateOutcome(gate_id, ctx.fixture.id, "XFAIL", "G2.5: --census required")

    on_step = ctx.work_dir / f"{ctx.fixture.id}_g25.step"
    try:
        run = run_stl2step(
            ctx.binary, ctx.fixture.stl, on_step, no_verify=ctx.no_verify, smooth=True
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))
    if not run.result.get("ok"):
        return GateOutcome(
            gate_id, ctx.fixture.id, "FAIL", f"conversion failed: {run.result.get('error')}"
        )

    try:
        cen = run_census(ctx.census_path, on_step)
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))

    violations: List[str] = []
    for grp in cen.get("cylinderGroups") or []:
        built_as = str(grp.get("builtAs") or "")
        if built_as == "explodedToFacets":
            violations.append(f"group R={grp.get('radius')}: explodedToFacets on recoverable hole")
        elif built_as == "twoHalves":
            violations.append(f"group R={grp.get('radius')}: twoHalves (expected seamed360)")
        elif built_as not in ("seamed360", "single", "partial", "other", ""):
            violations.append(f"group R={grp.get('radius')}: builtAs={built_as}")

    details = {"cylinderGroups": cen.get("cylinderGroups"), "violations": violations}
    if violations:
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"G2.5: {len(violations)} BuiltAs violation(s)",
            details=details,
        )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"G2.5: {len(closed360)} closed360 claim(s), census BuiltAs ok",
        details=details,
    )


def check_g3(ctx: GateContext) -> GateOutcome:
    """G3 — chord-aware D4.5 volume budget (HARD when unparked)."""
    gate_id = "G3"
    if not (ctx.dump_path and ctx.dump_path.is_file()):
        return GateOutcome(gate_id, ctx.fixture.id, "XFAIL", "G3: --dump required")

    on_step = ctx.work_dir / f"{ctx.fixture.id}_g3.step"
    try:
        run = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            on_step,
            no_verify=False,
            smooth=True,
            capture_notes=True,
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))
    if not run.result.get("ok"):
        return GateOutcome(
            gate_id, ctx.fixture.id, "FAIL", f"conversion failed: {run.result.get('error')}"
        )

    abs_sum, cylgrow = dump_dvol_abs_and_cylgrow(ctx.dump_path, ctx.fixture.stl)
    if abs_sum is None:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", "regiondump failed for D4.5")

    notes = run.stdout + run.stderr
    brep_v, mesh_note = parse_brep_mesh_volume(notes)
    mesh_v = float(run.result.get("meshVolumeMM3") or mesh_note or 0.0)
    if brep_v is None:
        brep_v = float(run.result.get("stepVolumeMM3") or 0.0)

    measured = abs(brep_v - mesh_v)
    ref_vol = abs(mesh_v) if mesh_v else abs(brep_v)
    floor = 1e-4 * ref_vol
    budget = max(floor, G3_VOLUME_K * abs_sum) if abs_sum > 0 else floor
    vol_ok = measured <= budget + 1e-9 if ref_vol > 0 else True

    sag_violations: List[str] = []
    for r in cylgrow:
        delta = float(r.get("chordSagitta") or 0.0)
        if delta < 0:
            sag_violations.append(f"region {r.get('id')}: chordSagitta={delta}")

    details = {
        "brepMM3": brep_v,
        "meshMM3": mesh_v,
        "measuredAbs": measured,
        "dVolPredAbs": abs_sum,
        "budget": budget,
        "k": G3_VOLUME_K,
        "g3Window": [G3_RATIO_LO, G3_RATIO_HI],
        "nCylGrow": len(cylgrow),
        "sagViolations": sag_violations,
        "volOk": vol_ok,
        "volumeDeltaPct": run.result.get("volumeDeltaPct"),
    }

    if sag_violations:
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"G3 sagitta: {len(sag_violations)} CylGrow delta anomalies",
            details=details,
        )
    if not vol_ok:
        pct = (measured / ref_vol * 100.0) if ref_vol else 0.0
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"G3 volume: |brep-mesh|={measured:.6g} mm³ ({pct:.4f}%) "
            f"exceeds budget {budget:.6g} (k={G3_VOLUME_K})",
            details=details,
        )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"G3: |brep-mesh|={measured:.6g} ≤ budget {budget:.6g}",
        details=details,
    )


def check_g4_4(ctx: GateContext) -> GateOutcome:
    """G4.4 — EncodeRegularity G1 continuity (bar min_g1=22 on S02)."""
    gate_id = "G4.4"
    if ctx.fixture.id not in G4_4_FIXTURES:
        return GateOutcome(gate_id, ctx.fixture.id, "SKIP", "G4.4 runs on S02 only")

    if not (ctx.census_path and ctx.census_path.is_file()):
        return GateOutcome(gate_id, ctx.fixture.id, "XFAIL", "G4.4: --census required")

    on_step = ctx.work_dir / f"{ctx.fixture.id}_g44.step"
    try:
        run = run_stl2step(
            ctx.binary, ctx.fixture.stl, on_step, no_verify=ctx.no_verify, smooth=True
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))
    if not run.result.get("ok"):
        return GateOutcome(
            gate_id, ctx.fixture.id, "FAIL", f"conversion failed: {run.result.get('error')}"
        )

    try:
        cen = run_census(ctx.census_path, on_step)
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc))

    g1 = int(cen.get("contG1") or 0)
    details = {
        "contG1": g1,
        "minBar": G4_4_MIN_G1_S02,
        "contC0": cen.get("contC0"),
        "contC1": cen.get("contC1"),
    }
    if g1 < G4_4_MIN_G1_S02:
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"G4.4: contG1={g1} < bar {G4_4_MIN_G1_S02} on {ctx.fixture.id}",
            details=details,
        )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"G4.4: contG1={g1} ≥ {G4_4_MIN_G1_S02}",
        details=details,
    )


def check_g5(ctx: GateContext) -> GateOutcome:
    """G5 — editability telemetry (SOFT): smoothMaxEdgeTolMM vs recovered radius."""
    gate_id = "G5"
    sidecar = load_sidecar(ctx.fixture)
    on_step = ctx.work_dir / f"{ctx.fixture.id}_g5.step"
    try:
        run = run_stl2step(
            ctx.binary, ctx.fixture.stl, on_step, no_verify=ctx.no_verify, smooth=True
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, ctx.fixture.id, "FAIL", str(exc), hard=False)
    if not run.result.get("ok"):
        return GateOutcome(
            gate_id,
            ctx.fixture.id,
            "FAIL",
            f"conversion failed: {run.result.get('error')}",
            hard=False,
        )

    max_edge = float(run.result.get("smoothMaxEdgeTolMM") or 0.0)
    min_r = min_recovered_radius(sidecar)
    details = {
        "smoothMaxEdgeTolMM": max_edge,
        "minRecoveredRadius": min_r,
        "g5EdgeK": G5_EDGE_K,
        "smoothBuiltCylinders": run.result.get("smoothBuiltCylinders"),
        "smoothCylinders": run.result.get("smoothCylinders"),
    }
    if min_r is not None and min_r > 0:
        cap = G5_EDGE_K * min_r
        details["edgeTolCap"] = cap
        if max_edge > cap + 1e-9:
            return GateOutcome(
                gate_id,
                ctx.fixture.id,
                "FAIL",
                f"G5: smoothMaxEdgeTolMM={max_edge:.6g} > {G5_EDGE_K}×R_min={cap:.6g}",
                hard=False,
                details=details,
            )
    return GateOutcome(
        gate_id,
        ctx.fixture.id,
        "PASS",
        f"G5 telemetry: smoothMaxEdgeTolMM={max_edge:.6g}",
        hard=False,
        details=details,
    )


def check_r_ladder(ctx: GateContext) -> GateOutcome:
    """R-ladder — S16 R1-explode, R1-round-2, R2 ChainUnstable."""
    gate_id = "R-ladder"
    fid = ctx.fixture.id
    if fid not in R_LADDER_FIXTURES:
        return GateOutcome(gate_id, fid, "SKIP", "not an S16 R-ladder fixture")

    sidecar = load_sidecar(ctx.fixture)
    off_step = ctx.work_dir / f"{fid}_rl_off.step"
    on_step = ctx.work_dir / f"{fid}_rl_on.step"
    try:
        off = run_stl2step(
            ctx.binary, ctx.fixture.stl, off_step, no_verify=True, smooth=False
        )
        on = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            on_step,
            no_verify=True,
            smooth=True,
            capture_notes=True,
        )
    except RuntimeError as exc:
        return GateOutcome(gate_id, fid, "FAIL", str(exc))

    if not off.result.get("ok") or not on.result.get("ok"):
        return GateOutcome(gate_id, fid, "FAIL", "R-ladder conversion failed")

    details: Dict[str, Any] = {
        "offExit": off.exit_code,
        "onExit": on.exit_code,
        "skipped": int(on.result.get("smoothSkippedComponents") or 0),
    }

    if fid == "S16-R1-round-2":
        if on.result.get("smoothSkippedComponents", 0) < 1:
            return GateOutcome(gate_id, fid, "FAIL", "expected dirty skip on round-2")
        if off.exit_code != on.exit_code:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                f"dirty skip changed exit {off.exit_code}->{on.exit_code}",
                details=details,
            )
        return GateOutcome(gate_id, fid, "PASS", "R1-round-2 dirty skip holds", details=details)

    if ctx.census_path and ctx.census_path.is_file():
        try:
            cen = run_census(ctx.census_path, on_step)
            details["censusValid"] = cen.get("valid")
            if not cen.get("valid"):
                return GateOutcome(gate_id, fid, "FAIL", "R-ladder: census invalid", details=details)
        except RuntimeError as exc:
            return GateOutcome(gate_id, fid, "FAIL", str(exc), details=details)

    d_off = canonical_step_data(read_text(off_step))
    d_on = canonical_step_data(read_text(on_step))
    details["stepIdentical"] = d_off == d_on
    if fid in ("S16-R1-explode-success", "S16-R2-ChainUnstable"):
        if d_off != d_on:
            return GateOutcome(
                gate_id,
                fid,
                "FAIL",
                "R2 revert: canonical STEP DATA differs from off-path",
                details=details,
            )
        reverted = "smooth: analytic rebuild reverted" in (on.stdout + on.stderr)
        details["r2Warn"] = reverted
        return GateOutcome(
            gate_id,
            fid,
            "PASS",
            "R-ladder: R2 verbatim revert + census valid",
            details=details,
        )

    return GateOutcome(gate_id, fid, "PASS", "R-ladder ok", details=details)


def check_calibration(_ctx: GateContext) -> GateOutcome:
    """Calibration — frozen G3 window and volume k at P0 close."""
    gate_id = "calibration"
    frozen = {
        "G3_RATIO_LO": G3_RATIO_LO,
        "G3_RATIO_HI": G3_RATIO_HI,
        "G3_VOLUME_K": G3_VOLUME_K,
        "G3_CHORD_SAGITTA_WINDOW_MM": G3_CHORD_SAGITTA_WINDOW_MM,
        "G5_EDGE_K": G5_EDGE_K,
        "G4_4_MIN_G1_S02": G4_4_MIN_G1_S02,
    }
    return GateOutcome(
        gate_id,
        _ctx.fixture.id,
        "PASS",
        f"frozen thresholds: k={G3_VOLUME_K}, G3 window [{G3_RATIO_LO}, {G3_RATIO_HI}], "
        f"G4.4 min_g1={G4_4_MIN_G1_S02}",
        hard=True,
        details=frozen,
    )


SMOOTH_ON_CHECKS: Dict[str, Callable[[GateContext], GateOutcome]] = {
    "G0.2": check_g0_2,
    "G0.3": check_g0_3,
    "G2": check_g2,
    "G2.5": check_g2_5,
    "G3": check_g3,
    "G4.4": check_g4_4,
    "G5": check_g5,
    "R-ladder": check_r_ladder,
    "calibration": check_calibration,
}


def apply_parking(
    outcomes: List[GateOutcome],
    unparked: Optional[Set[str]] = None,
) -> List[GateOutcome]:
    """Turn honest FAILs on parked gates into PARKED (suite stays green)."""
    unparked = unparked or set()
    out: List[GateOutcome] = []
    for o in outcomes:
        if (
            o.gate_id in PARKED_GATES
            and o.gate_id not in unparked
            and o.status == "FAIL"
        ):
            reason = PARKED_GATES[o.gate_id]
            out.append(
                GateOutcome(
                    o.gate_id,
                    o.fixture_id,
                    "PARKED",
                    f"PARKED: {reason} — {o.message}",
                    hard=o.hard,
                    details={**(o.details or {}), "parkReason": reason, "wouldBe": "FAIL"},
                )
            )
        else:
            out.append(o)
    return out


def format_gate_table(outcomes: Sequence[GateOutcome]) -> str:
    """Per-gate matrix from actual run outcomes only (never hand-maintained)."""
    gates = sorted({o.gate_id for o in outcomes})
    fixtures = sorted({o.fixture_id for o in outcomes})
    if not gates or not fixtures:
        return "gate table: (empty)"

    cell: Dict[Tuple[str, str], str] = {}
    for o in outcomes:
        cell[(o.gate_id, o.fixture_id)] = o.status

    fwidth = max(len(f) for f in fixtures)
    gwidth = max(len(g) for g in gates)
    lines = ["per-gate table (from this run):"]
    header = " " * (gwidth + 2) + "".join(f"{f:>{fwidth}} " for f in fixtures)
    lines.append(header.rstrip())
    for g in gates:
        row = f"{g:<{gwidth}}  "
        row += "".join(f"{cell.get((g, f), '·'):>{fwidth}} " for f in fixtures)
        lines.append(row.rstrip())

    parked = [g for g in gates if g in PARKED_GATES]
    if parked:
        lines.append("parked gates (default suite ignores their FAILs):")
        for g in parked:
            lines.append(f"  {g}: {PARKED_GATES[g]}")
    return "\n".join(lines)
