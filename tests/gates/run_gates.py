#!/usr/bin/env python3
"""SPEC-P0 acceptance gate runner for stl2step (stdlib only)."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set

# CALIBRATED AT P0 CLOSE
G3_CHORD_SAGITTA_WINDOW_MM = 0.0  # placeholder — freeze at P0 close
G3_VOLUME_K = 0.0  # placeholder — freeze at P0 close
G5_EDGE_K = 0.1  # start 0.1 × min recovered radius (SOFT)

LEGACY_VOLUME_GATE_PCT = 0.01  # applies with smooth off only

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = REPO_ROOT / "tests" / "corpus"
DEFAULT_BASELINE = REPO_ROOT / "tests" / "gates" / "baseline"
SMOKE_STL = REPO_ROOT / "tests" / "cube.stl"

ALL_GATE_IDS = (
    "G0.1",
    "G0.2",
    "G0.3",
    "G1",
    "G2",
    "G2.5",
    "G3",
    "G4",
    "G4.4",
    "G5",
    "I-checker",
    "include-allowlist",
    "calibration",
    "R-ladder",
)

# Gates whose implementation phase has landed in wave 2.
LIVE_GATES: Set[str] = {"G0.1", "G1", "G4", "include-allowlist"}

HARD_GATES: Set[str] = {
    "G0.1",
    "G0.2",
    "G0.3",
    "G1",
    "G2",
    "G2.5",
    "G4",
    "G4.4",
    "I-checker",
    "include-allowlist",
    "R-ladder",
}

SOFT_GATES: Set[str] = {"G3", "G5"}

INCLUDE_ALLOWLIST_FILES = (
    REPO_ROOT / "src" / "refit_segment.cpp",
    REPO_ROOT / "src" / "refit_math.cpp",
    REPO_ROOT / "src" / "refit_grow.cpp",
    REPO_ROOT / "src" / "refit_fillet.cpp",
    REPO_ROOT / "src" / "refit_chains.cpp",
    REPO_ROOT / "src" / "refit_internal.hpp",
)

# DECISION-p1-math D5.3 — exact allowlist for P1 TUs + refit_internal.hpp
INCLUDE_ALLOWED_STDLIB = frozenset(
    {
        "algorithm",
        "array",
        "cmath",
        "cstdint",
        "functional",
        "limits",
        "numeric",
        "string",
        "utility",
        "vector",
    }
)
INCLUDE_ALLOWED_PROJECT = frozenset({"refit.hpp", "refit_internal.hpp"})
INCLUDE_BANNED_PREFIXES = (
    "GProp_",
    "GeomConvert_",
    "ShapeAnalysis_",
    "TopoDS",
    "BRep",
    "Geom",
    "Poly_",
)


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
    baseline_dir: Optional[Path]
    census_path: Optional[Path]
    no_verify: bool
    smooth: bool = False


@dataclass
class GateOutcome:
    gate_id: str
    fixture_id: str
    status: str  # PASS | FAIL | XFAIL
    message: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="SPEC-P0 acceptance gate runner for stl2step.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Shard examples:
              run_gates.py --smoke --binary ./build/stl2step
              run_gates.py --fixture S01,S06,S09 --gate G0.1,G4 --binary ./build/stl2step
            """
        ),
    )
    p.add_argument("--fixture", help="Comma-separated fixture ids (e.g. S01,S06,S09)")
    p.add_argument("--gate", help="Comma-separated gate ids (e.g. G0.1,G4)")
    p.add_argument("--smoke", action="store_true", help="Smoke run on tests/cube.stl")
    p.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS, help="Corpus directory")
    p.add_argument("--binary", type=Path, required=True, help="Path to stl2step CLI")
    p.add_argument("--census", type=Path, help="Path to stl2step_census binary")
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE, help="G0.1 baseline dir")
    p.add_argument("--no-verify", action="store_true", help="Pass --no-verify to stl2step")
    p.add_argument("--jobs", type=int, default=1, help="Parallel fixture workers")
    p.add_argument(
        "--report",
        type=Path,
        default=Path("gates-report.json"),
        help="Machine-readable report path",
    )
    return p.parse_args(argv)


def split_csv(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    return [part.strip() for part in value.split(",") if part.strip()]


def discover_fixtures(corpus: Path, smoke: bool, fixture_filter: Optional[List[str]]) -> List[Fixture]:
    if smoke:
        fixtures = [Fixture("cube", SMOKE_STL)]
    else:
        fixtures = []
        if corpus.is_dir():
            for stl in sorted(corpus.glob("S*.stl")):
                fid = stl.stem
                sidecar = corpus / f"{fid}.expected.json"
                fixtures.append(
                    Fixture(fid, stl, sidecar if sidecar.is_file() else None)
                )
        if not fixtures and SMOKE_STL.is_file():
            fixtures = [Fixture("cube", SMOKE_STL)]

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


def strip_file_name_timestamp(step_text: str) -> str:
    return re.sub(
        r"(FILE_NAME\s*\(\s*'[^']*'\s*,\s*)'[^']*'",
        r"\1'STRIPPED'",
        step_text,
        count=1,
        flags=re.IGNORECASE,
    )


def extract_data_section(step_text: str) -> str:
    match = re.search(r"DATA;\s*(.*?)\s*ENDSEC;", step_text, re.DOTALL | re.IGNORECASE)
    if not match:
        return step_text
    return match.group(0)


def canonical_step_data(step_text: str) -> str:
    return extract_data_section(strip_file_name_timestamp(step_text))


def canonicalize_result(result: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, value in result.items():
        if key == "seconds":
            continue
        if key in {"input", "output"} and isinstance(value, str):
            out[key] = Path(value).name
        else:
            out[key] = value
    return out


def canonicalize_result_no_paths(result: Dict[str, Any]) -> Dict[str, Any]:
    skip = {"seconds", "input", "output"}
    return {key: value for key, value in result.items() if key not in skip}


def result_has_smooth_keys(result: Dict[str, Any]) -> List[str]:
    return sorted(k for k in result if k.startswith("smooth"))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def g2_radius_threshold(r_truth: float, n_sides: int, refit_tol: float) -> float:
    n = max(int(n_sides), 3)
    chord = r_truth * (1.0 - math.cos(math.pi / n))
    return max(refit_tol, 0.005 * r_truth, chord)


def include_is_allowed(header: str) -> bool:
    """Return True when header is on the D5.3 allowlist."""
    base = Path(header).name
    if base in INCLUDE_ALLOWED_PROJECT:
        return True
    if base in INCLUDE_ALLOWED_STDLIB:
        return True
    # Angle-bracket stdlib generally: extensionless lowercase names.
    if not base.endswith(".hxx") and re.fullmatch(r"[a-z][a-z0-9_]*", base):
        return True
    if base == "Precision.hxx":
        return True
    if base.startswith("Standard_") and base.endswith(".hxx"):
        return True
    if base.startswith("gp_") and base.endswith(".hxx"):
        return True
    if base.startswith("math_") and base.endswith(".hxx"):
        return True
    return False


def include_is_banned(header: str) -> bool:
    """Return True for explicitly banned or other OCCT topology/algorithm headers."""
    base = Path(header).name
    if any(base.startswith(prefix) for prefix in INCLUDE_BANNED_PREFIXES):
        return True
    # Any other OCCT .hxx not on the D5.3 allowlist is banned topology/algorithm.
    if base.endswith(".hxx") and not include_is_allowed(header):
        return True
    return False


def xfail_not_landed(gate_id: str, fixture_id: str, phase: str) -> GateOutcome:
    return GateOutcome(
        gate_id,
        fixture_id,
        "XFAIL",
        f"{gate_id}: phase {phase} not landed — expected fail until downstream lane ships",
        hard=gate_id in HARD_GATES,
    )


def go(
    ctx: GateContext,
    gate_id: str,
    status: str,
    message: str,
    *,
    hard: bool = True,
    details: Optional[Dict[str, Any]] = None,
) -> GateOutcome:
    return GateOutcome(gate_id, ctx.fixture.id, status, message, hard=hard, details=details or {})


def gate_g0_1_off_path_identity(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.1 — off-path identity, canonicalized vs 187ead0 baseline."""
    gate_id = "G0.1"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P0-baseline")

    baseline_dir = ctx.baseline_dir
    if baseline_dir is None or not baseline_dir.is_dir():
        return go(
            ctx,
            gate_id,
            "XFAIL",
            "baseline directory absent (tests/gates/baseline/ from p0-golden not landed yet)",
        )

    fixture = ctx.fixture
    baseline_step = baseline_dir / "step" / f"{fixture.id}.step"
    baseline_result = baseline_dir / "result" / f"{fixture.id}.json"
    if not baseline_step.is_file() or not baseline_result.is_file():
        return go(
            ctx,
            gate_id,
            "XFAIL",
            f"baseline artifacts missing for fixture {fixture.id}",
        )

    out_step = ctx.work_dir / f"{fixture.id}.step"
    try:
        run = run_stl2step(
            ctx.binary,
            fixture.stl,
            out_step,
            no_verify=ctx.no_verify,
            smooth=False,
        )
    except RuntimeError as exc:
        return go(ctx, gate_id, "FAIL", str(exc))

    if not run.result.get("ok"):
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"off-path conversion failed: {run.result.get('error')}",
        )

    smooth_keys = result_has_smooth_keys(run.result)
    if smooth_keys:
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"smooth* keys present with --smooth absent: {smooth_keys}",
        )

    cur_data = canonical_step_data(read_text(out_step))
    base_data = canonical_step_data(read_text(baseline_step))
    if cur_data != base_data:
        diff = "\n".join(
            difflib.unified_diff(
                base_data.splitlines(),
                cur_data.splitlines(),
                fromfile="baseline",
                tofile="current",
                lineterm="",
            )
        )
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonical STEP DATA section differs from 187ead0 baseline",
            details={"diff_head": diff.splitlines()[:40]},
        )

    cur_res = canonicalize_result(run.result)
    base_res = canonicalize_result(json.loads(read_text(baseline_result)))
    if cur_res != base_res:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonical RESULT differs from 187ead0 baseline",
            details={"current": cur_res, "baseline": base_res},
        )

    return go(ctx, gate_id, "PASS", "off-path STEP DATA + RESULT match baseline")


def gate_g0_2_refit_closedness(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.2 — refit never reduces closedness (solids/openShells)."""
    return xfail_not_landed("G0.2", ctx.fixture.id, "P3-engine/smooth")


def gate_g0_3_dirty_skip(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.3 — dirty skip on S15; force-sew global disable."""
    return xfail_not_landed("G0.3", ctx.fixture.id, "P3-engine/smooth")


def gate_g1_validity(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G1 — BRepCheck valid; ShapeFix rewrite path did not fire."""
    gate_id = "G1"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P0-census")

    out_step = ctx.work_dir / f"{ctx.fixture.id}_g1.step"
    try:
        run = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            out_step,
            no_verify=ctx.no_verify,
            smooth=False,
            capture_notes=True,
        )
    except RuntimeError as exc:
        return go(ctx, gate_id, "FAIL", str(exc))

    if not run.result.get("ok"):
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"conversion failed: {run.result.get('error')}",
        )

    combined = run.stdout + run.stderr
    if re.search(r"\bfix\b.*ShapeFix|ShapeFix.*\bfix\b", combined, re.IGNORECASE):
        return go(
            ctx,
            gate_id,
            "FAIL",
            "ShapeFix invalid-rewrite path fired (fix note in CLI output)",
        )

    for warning in run.result.get("warnings", []):
        if "ShapeFix" in warning or re.search(r"\bfix\b", warning, re.IGNORECASE):
            return go(ctx, gate_id, "FAIL", f"ShapeFix-related warning: {warning}")

    if ctx.census_path and ctx.census_path.is_file():
        proc = subprocess.run(
            [str(ctx.census_path), str(out_step)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            return go(
                ctx,
                gate_id,
                "FAIL",
                f"census tool failed: {proc.stderr or proc.stdout}",
            )
        try:
            census = json.loads(proc.stdout)
        except json.JSONDecodeError:
            return go(ctx, gate_id, "FAIL", "census output is not JSON")
        if not census.get("valid", False):
            return go(ctx, gate_id, "FAIL", f"BRepCheck invalid per census: {census}")
        return go(ctx, gate_id, "PASS", "BRepCheck valid (census)")

    if run.result.get("watertight") is False or run.result.get("openShells", 0) > 0:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "mesh not watertight / open shells present (proxy validity without census)",
        )

    return go(
        ctx,
        gate_id,
        "PASS",
        "no ShapeFix rewrite; conversion ok (census unavailable — stderr/RESULT proxy)",
    )


def gate_g2_recognition(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G2 — recoverable primitive recognition vs sidecar."""
    if ctx.fixture.sidecar and ctx.fixture.sidecar.is_file():
        sidecar = json.loads(read_text(ctx.fixture.sidecar))
        for prim in sidecar.get("recoverable", []):
            if prim.get("type") != "Cylinder":
                continue
            n_sides = prim.get("nSides")
            if n_sides is None:
                continue
            g2_radius_threshold(float(prim.get("radius", 0.0)), int(n_sides), 0.0)
    return xfail_not_landed("G2", ctx.fixture.id, "P2-build + p0-census")


def gate_g2_5_built_as(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G2.5 — closed360 holes built as one cylindrical face."""
    return xfail_not_landed("G2.5", ctx.fixture.id, "P2-build + p0-census")


def gate_g3_sagitta_volume(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G3 — sagitta sanity + volume budget (SOFT until P3)."""
    # Off-path mesh-volume gate uses LEGACY_VOLUME_GATE_PCT; retired on smooth runs.
    _off_path_volume_pct = LEGACY_VOLUME_GATE_PCT
    return xfail_not_landed("G3", ctx.fixture.id, "P3-engine/smooth")


def gate_g4_determinism(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G4 — determinism: --threads 1 vs --threads 8 canonical STEP DATA."""
    gate_id = "G4"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P0-gates")

    step1 = ctx.work_dir / f"{ctx.fixture.id}_t1.step"
    step8 = ctx.work_dir / f"{ctx.fixture.id}_t8.step"
    try:
        run1 = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            step1,
            no_verify=ctx.no_verify,
            threads=1,
            smooth=False,
        )
        run8 = run_stl2step(
            ctx.binary,
            ctx.fixture.stl,
            step8,
            no_verify=ctx.no_verify,
            threads=8,
            smooth=False,
        )
    except RuntimeError as exc:
        return go(ctx, gate_id, "FAIL", str(exc))

    for label, run in (("threads=1", run1), ("threads=8", run8)):
        if not run.result.get("ok"):
            return go(
                ctx,
                gate_id,
                "FAIL",
                f"{label} conversion failed: {run.result.get('error')}",
            )

    data1 = canonical_step_data(read_text(step1))
    data8 = canonical_step_data(read_text(step8))
    if data1 != data8:
        diff = "\n".join(
            difflib.unified_diff(
                data1.splitlines(),
                data8.splitlines(),
                fromfile="threads1",
                tofile="threads8",
                lineterm="",
            )
        )
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonical STEP DATA differs between --threads 1 and --threads 8",
            details={
                "sha256_t1": sha256_text(data1),
                "sha256_t8": sha256_text(data8),
                "diff_head": diff.splitlines()[:40],
            },
        )

    res1 = canonicalize_result_no_paths(run1.result)
    res8 = canonicalize_result_no_paths(run8.result)
    if res1 != res8:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonical RESULT differs between thread counts",
            details={"threads1": res1, "threads8": res8},
        )

    return go(
        ctx,
        gate_id,
        "PASS",
        "canonical STEP DATA + RESULT identical for --threads 1 vs 8",
        details={"data_sha256": sha256_text(data1)},
    )


def gate_g4_4_encode_regularity(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G4.4 — EncodeRegularity G1 continuity on S02 fillet edges."""
    return xfail_not_landed("G4.4", ctx.fixture.id, "P3-engine/release")


def gate_g5_editability(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G5 — editability censuses; smoothMaxEdgeTolMM threshold (SOFT)."""
    return xfail_not_landed("G5", ctx.fixture.id, "P3-engine/smooth")


def gate_i_checker(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 I-checker — RegionSet dump validation (check_regionset.py)."""
    return xfail_not_landed("I-checker", ctx.fixture.id, "p0-ichecker + P1 dump")


def gate_include_allowlist(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 include-allowlist grep (DECISION D5.3) on all P1 TUs + refit_internal.hpp."""
    gate_id = "include-allowlist"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P1-segment")

    existing = [p for p in INCLUDE_ALLOWLIST_FILES if p.is_file()]
    if not existing:
        return go(
            ctx,
            gate_id,
            "PASS",
            "P1 TUs + refit_internal.hpp absent — vacuous pass",
        )

    violations: List[str] = []
    for path in existing:
        for line_no, line in enumerate(read_text(path).splitlines(), 1):
            m = re.match(r'^\s*#\s*include\s*[<"]([^">]+)[">]', line)
            if not m:
                continue
            header = m.group(1)
            if include_is_banned(header):
                violations.append(f"{path}:{line_no}: banned include {header}")
            elif not include_is_allowed(header):
                violations.append(f"{path}:{line_no}: non-allowlisted include {header}")

    if violations:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "banned or non-allowlisted #include in P1 sources",
            details={"violations": violations},
        )

    return go(ctx, gate_id, "PASS", "all #includes within D5.3 allowlist")


def gate_calibration(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 calibration — freeze G3 window and volume k at P0 close."""
    return xfail_not_landed("calibration", ctx.fixture.id, "P0-corpus-close")


def gate_r_ladder(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 R-ladder — S16 R1-explode, R1-round-2, R2 revert fixtures."""
    return xfail_not_landed("R-ladder", ctx.fixture.id, "P0-corpus S16 + P2-build")


GATE_REGISTRY: Dict[str, Callable[[GateContext], GateOutcome]] = {
    "G0.1": gate_g0_1_off_path_identity,
    "G0.2": gate_g0_2_refit_closedness,
    "G0.3": gate_g0_3_dirty_skip,
    "G1": gate_g1_validity,
    "G2": gate_g2_recognition,
    "G2.5": gate_g2_5_built_as,
    "G3": gate_g3_sagitta_volume,
    "G4": gate_g4_determinism,
    "G4.4": gate_g4_4_encode_regularity,
    "G5": gate_g5_editability,
    "I-checker": gate_i_checker,
    "include-allowlist": gate_include_allowlist,
    "calibration": gate_calibration,
    "R-ladder": gate_r_ladder,
}


def selected_gates(gate_filter: Optional[List[str]]) -> List[str]:
    if gate_filter:
        unknown = [g for g in gate_filter if g not in GATE_REGISTRY]
        if unknown:
            raise SystemExit(f"unknown gate ids: {', '.join(unknown)}")
        return gate_filter
    return list(ALL_GATE_IDS)


def run_fixture_gates(
    fixture: Fixture,
    gate_ids: Sequence[str],
    binary: Path,
    baseline_dir: Optional[Path],
    census_path: Optional[Path],
    no_verify: bool,
) -> List[GateOutcome]:
    outcomes: List[GateOutcome] = []
    with tempfile.TemporaryDirectory(prefix=f"gates_{fixture.id}_") as tmp:
        work_dir = Path(tmp)
        ctx = GateContext(
            binary=binary,
            fixture=fixture,
            work_dir=work_dir,
            baseline_dir=baseline_dir,
            census_path=census_path,
            no_verify=no_verify,
        )
        for gate_id in gate_ids:
            fn = GATE_REGISTRY[gate_id]
            outcomes.append(fn(ctx))
    return outcomes


def summarize(outcomes: List[GateOutcome]) -> str:
    lines = ["gates summary:"]
    width = max(len(o.gate_id) for o in outcomes) if outcomes else 10
    for o in outcomes:
        lines.append(f"  {o.gate_id:<{width}}  {o.status:<6}  {o.message}")
    counts = {"PASS": 0, "FAIL": 0, "XFAIL": 0}
    for o in outcomes:
        counts[o.status] = counts.get(o.status, 0) + 1
    lines.append(
        f"totals: PASS={counts.get('PASS', 0)} "
        f"FAIL={counts.get('FAIL', 0)} XFAIL={counts.get('XFAIL', 0)}"
    )
    return "\n".join(lines)


def build_report(
    fixtures: List[Fixture],
    gate_ids: List[str],
    outcomes: List[GateOutcome],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    return {
        "binary": str(args.binary),
        "smoke": args.smoke,
        "fixtures": [f.id for f in fixtures],
        "gates": gate_ids,
        "results": [
            {
                "gate": o.gate_id,
                "fixture": o.fixture_id,
                "status": o.status,
                "message": o.message,
                "hard": o.hard,
                "details": o.details,
            }
            for o in outcomes
        ],
        "summary": {
            "PASS": sum(1 for o in outcomes if o.status == "PASS"),
            "FAIL": sum(1 for o in outcomes if o.status == "FAIL"),
            "XFAIL": sum(1 for o in outcomes if o.status == "XFAIL"),
        },
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    binary = args.binary.resolve()
    if not binary.is_file():
        print(f"error: binary not found: {binary}", file=sys.stderr)
        return 1

    fixture_filter = split_csv(args.fixture)
    gate_filter = split_csv(args.gate)
    gate_ids = selected_gates(gate_filter)
    fixtures = discover_fixtures(args.corpus, args.smoke, fixture_filter)
    if not fixtures:
        print("error: no fixtures selected", file=sys.stderr)
        return 1

    baseline_dir = args.baseline.resolve() if args.baseline else None
    census_path = args.census.resolve() if args.census else None

    all_outcomes: List[GateOutcome] = []
    if args.jobs > 1 and len(fixtures) > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(
                    run_fixture_gates,
                    fixture,
                    gate_ids,
                    binary,
                    baseline_dir,
                    census_path,
                    args.no_verify,
                ): fixture
                for fixture in fixtures
            }
            for fut in as_completed(futures):
                all_outcomes.extend(fut.result())
    else:
        for fixture in fixtures:
            all_outcomes.extend(
                run_fixture_gates(
                    fixture,
                    gate_ids,
                    binary,
                    baseline_dir,
                    census_path,
                    args.no_verify,
                )
            )

    report = build_report(fixtures, gate_ids, all_outcomes, args)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(summarize(all_outcomes))
    print(f"report: {args.report}")

    hard_failures = [
        o
        for o in all_outcomes
        if o.status == "FAIL" and o.hard and o.gate_id not in SOFT_GATES
    ]
    return 1 if hard_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
