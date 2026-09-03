#!/usr/bin/env python3
"""SPEC-P0 acceptance gate runner for stl2step (stdlib only)."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set

import smooth_on as so

# Re-export frozen thresholds (single source: smooth_on.py).
G3_RATIO_LO = so.G3_RATIO_LO
G3_RATIO_HI = so.G3_RATIO_HI
G3_CHORD_SAGITTA_WINDOW_MM = so.G3_CHORD_SAGITTA_WINDOW_MM
G3_VOLUME_K = so.G3_VOLUME_K
G5_EDGE_K = so.G5_EDGE_K
G4_4_MIN_G1_S02 = so.G4_4_MIN_G1_S02
LEGACY_VOLUME_GATE_PCT = so.LEGACY_VOLUME_GATE_PCT
PARKED_GATES = so.PARKED_GATES

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
    "G-LAW",
)

# Gates whose implementation has landed (off-path + smooth-on).
LIVE_GATES: Set[str] = {
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
    "R-ladder",
    "calibration",
    "G-LAW",
}

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
    "G-LAW",
}

SOFT_GATES: Set[str] = {"G5"}

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
    dump_path: Optional[Path]
    ichecker_path: Optional[Path]
    no_verify: bool
    smooth: bool = False
    baseline_bin: Optional[Path] = None
    baseline_error: Optional[str] = None


@dataclass
class GateOutcome:
    gate_id: str
    fixture_id: str
    status: str  # PASS | FAIL | XFAIL | SKIP
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
              run_gates.py --unpark G3,G4.4 --binary ./build/stl2step  # surface parked reds
            """
        ),
    )
    p.add_argument("--fixture", help="Comma-separated fixture ids (e.g. S01,S06,S09)")
    p.add_argument("--gate", help="Comma-separated gate ids (e.g. G0.1,G4)")
    p.add_argument("--smoke", action="store_true", help="Smoke run on tests/cube.stl")
    p.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS, help="Corpus directory")
    p.add_argument("--binary", type=Path, required=True, help="Path to stl2step CLI")
    p.add_argument("--census", type=Path, help="Path to stl2step_census binary")
    p.add_argument("--dump", type=Path,
                   help="Path to stl2step_regiondump binary (I-checker gate)")
    p.add_argument("--ichecker", type=Path,
                   help="Path to tests/gates/check_regionset.py (I-checker gate)")
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE, help="G0.1 baseline dir")
    p.add_argument("--no-verify", action="store_true", help="Pass --no-verify to stl2step")
    p.add_argument("--jobs", type=int, default=1, help="Parallel fixture workers")
    p.add_argument(
        "--report",
        type=Path,
        default=Path("gates-report.json"),
        help="Machine-readable report path",
    )
    p.add_argument(
        "--unpark",
        help="Comma-separated parked gate ids to surface honest FAILs (e.g. G3,G4.4)",
    )
    return p.parse_args(argv)


def split_csv(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    return [part.strip() for part in value.split(",") if part.strip()]


def discover_fixtures(corpus: Path, smoke: bool, fixture_filter: Optional[List[str]]) -> List[Fixture]:
    raw = so.discover_fixtures(corpus, SMOKE_STL, smoke, fixture_filter)
    return [Fixture(f.id, f.stl, f.sidecar) for f in raw]


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
    # G0.1 / canonicalize.py: smooth* glob plus the one sibling that does not match it.
    return sorted(k for k in result if k.startswith("smooth") or k == "facesAfterSmooth")


# Frozen tests/gates/regionset.schema.json vocabulary: lowerCamelCase of the
# C++ enumerator. Producers have emitted PascalCase ("Cylinder"); compare only
# after this mapping so a future producer cannot reintroduce the mismatch.
SCHEMA_SURFACE_TYPES = frozenset({"plane", "cylinder", "cone", "sphere", "torus"})
SCHEMA_ORIGINS = frozenset({"planeGrow", "cylGrow", "filletStrip"})
SCHEMA_ROLES = frozenset({"outer", "inner", "capLow", "capHigh"})
SCHEMA_BUILT_AS = frozenset(
    {"notBuilt", "single", "seamed360", "twoHalves", "explodedToFacets"}
)


def schema_enum(value: Any) -> str:
    """Map PascalCase / lowerCamelCase enumerator spellings to the frozen schema.

    Cylinder → cylinder, FilletStrip → filletStrip, seamed360 stays.
    First-character lower is the C++ enumerator → schema mapping.
    """
    if not isinstance(value, str) or not value:
        return ""
    return value[0].lower() + value[1:]


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


_CONFIG_BUILD_DIRS = frozenset({"Release", "Debug", "RelWithDebInfo", "MinSizeRel"})


def current_build_from_binary(binary: Path) -> Path:
    """CMake build dir for the candidate CLI (--current-build for the baseline)."""
    parent = binary.resolve().parent
    if parent.name in _CONFIG_BUILD_DIRS:
        return parent.parent
    return parent


def _baseline_bin_exists(path: Path) -> bool:
    """Accept Unix, .exe, and multi-config layout paths."""
    if path.is_file():
        return True
    if path.suffix.lower() != ".exe" and Path(str(path) + ".exe").is_file():
        return True
    return False


def _resolve_baseline_bin(path: Path) -> Optional[Path]:
    if path.is_file():
        return path
    exe = Path(str(path) + ".exe")
    if path.suffix.lower() != ".exe" and exe.is_file():
        return exe
    return None


_WIN_GIT_BASH = (
    r"C:\Program Files\Git\bin\bash.exe",
    r"C:\Program Files\Git\usr\bin\bash.exe",
    r"C:\Program Files (x86)\Git\bin\bash.exe",
)


def _under_windows_system_dir(path: str) -> bool:
    """True if *path* lives under %SystemRoot% (the WSL stub is System32\\bash.exe)."""
    root = os.environ.get("SystemRoot", r"C:\Windows")
    abs_path = os.path.normcase(os.path.abspath(path))
    abs_root = os.path.normcase(os.path.abspath(root))
    prefix = abs_root if abs_root.endswith(os.sep) else abs_root + os.sep
    return abs_path == abs_root or abs_path.startswith(prefix)


def _bash_tried_locations() -> List[str]:
    tried: List[str] = []
    env = os.environ.get("STL2STEP_BASH")
    tried.append(f"STL2STEP_BASH={env}" if env else "STL2STEP_BASH (unset)")
    tried.extend(_WIN_GIT_BASH)
    found = shutil.which("bash")
    if found:
        note = " (Windows system directory, skipped)" if _under_windows_system_dir(found) else ""
        tried.append(f"shutil.which(bash)={found}{note}")
    else:
        tried.append("shutil.which(bash)=(none)")
    return tried


def _bash_executable() -> Optional[str]:
    """Resolve bash for build_baseline.sh.

    On Windows, CreateProcess searches System32 before PATH, so a bare
    ``bash`` name always hits the WSL stub even when Git Bash is on PATH.
    Prefer an explicit Git for Windows path (or STL2STEP_BASH).
    """
    if sys.platform != "win32":
        return shutil.which("bash") or "bash"

    env = os.environ.get("STL2STEP_BASH")
    for cand in ((env,) if env else ()) + _WIN_GIT_BASH:
        if cand and os.path.isfile(cand):
            return cand
    found = shutil.which("bash")
    if found and not _under_windows_system_dir(found):
        return found
    return None


def ensure_baseline(
    baseline_dir: Optional[Path],
    current_build: Optional[Path] = None,
) -> "tuple[Optional[Path], Optional[str]]":
    """Build/locate the 187ead0 CLI once per run via p0-golden's build_baseline.sh.

    XFAIL is legal only when the baseline genuinely cannot be built. Never
    XFAIL on missing canned STEP goldens — those files are gitignored on
    purpose.
    """
    if baseline_dir is None or not baseline_dir.is_dir():
        return None, (
            "baseline cannot be built: tests/gates/baseline/ directory is absent"
        )
    script = baseline_dir / "build_baseline.sh"
    if not script.is_file():
        return None, f"baseline cannot be built: {script} is absent"

    print(
        "building 187ead0 baseline via tests/gates/baseline/build_baseline.sh ...",
        file=sys.stderr,
        flush=True,
    )
    bash_exe = _bash_executable()
    if bash_exe is None:
        tried = "; ".join(_bash_tried_locations())
        return None, (
            "baseline cannot be built: no Git Bash found; tried: " + tried
        )
    cmd = [bash_exe, str(script)]
    if current_build is not None:
        cmd.extend(["--current-build", current_build.resolve().as_posix()])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    combined = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        tail = "\n".join(combined.strip().splitlines()[-30:])
        return None, (
            f"baseline cannot be built: build_baseline.sh exited {proc.returncode}\n"
            f"{tail}"
        )

    bin_path: Optional[Path] = None
    for line in (proc.stdout or "").splitlines():
        if line.startswith("BASELINE_BIN="):
            candidate = line.split("=", 1)[1].strip()
            if candidate:
                bin_path = _resolve_baseline_bin(Path(candidate))
            break
    if bin_path is None:
        build = baseline_dir / ".build"
        for fallback in (
            build / "stl2step",
            build / "stl2step.exe",
            build / "Release" / "stl2step",
            build / "Release" / "stl2step.exe",
            build / "Debug" / "stl2step",
            build / "Debug" / "stl2step.exe",
            build / "RelWithDebInfo" / "stl2step",
            build / "RelWithDebInfo" / "stl2step.exe",
            build / "MinSizeRel" / "stl2step",
            build / "MinSizeRel" / "stl2step.exe",
        ):
            hit = _resolve_baseline_bin(fallback)
            if hit is not None:
                bin_path = hit
                break
    if bin_path is None or not _baseline_bin_exists(bin_path):
        tail = combined.strip()[-2000:]
        return None, (
            "baseline cannot be built: build_baseline.sh exited 0 but "
            f"BASELINE_BIN is missing\n{tail}"
        )
    print(f"baseline ready: {bin_path}", file=sys.stderr, flush=True)
    return bin_path, None


def run_canonicalize(
    tool: Path, mode: str, path_a: Path, path_b: Path
) -> subprocess.CompletedProcess:
    """Call the golden lane's canonicalize.py; never reimplement its strip set."""
    return subprocess.run(
        [sys.executable, str(tool), mode, str(path_a), str(path_b)],
        capture_output=True,
        text=True,
    )


def gate_g0_1_off_path_identity(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.1 — live twin-run vs 187ead0, DATA-compared through canonicalize.py."""
    gate_id = "G0.1"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P0-baseline")

    if ctx.baseline_error:
        return go(ctx, gate_id, "XFAIL", ctx.baseline_error)

    if ctx.baseline_bin is None or not ctx.baseline_bin.is_file():
        return go(
            ctx,
            gate_id,
            "XFAIL",
            "baseline cannot be built: 187ead0 binary not available "
            "(build_baseline.sh did not yield BASELINE_BIN)",
        )

    baseline_dir = ctx.baseline_dir or DEFAULT_BASELINE
    canon = baseline_dir / "canonicalize.py"
    if not canon.is_file():
        return go(
            ctx,
            gate_id,
            "XFAIL",
            f"baseline cannot be built: canonicalize.py missing at {canon}",
        )

    fixture = ctx.fixture
    cur_step = ctx.work_dir / f"{fixture.id}.cur.step"
    bas_step = ctx.work_dir / f"{fixture.id}.bas.step"
    cur_txt = ctx.work_dir / f"{fixture.id}.cur.txt"
    bas_txt = ctx.work_dir / f"{fixture.id}.bas.txt"

    try:
        cur = run_stl2step(
            ctx.binary,
            fixture.stl,
            cur_step,
            no_verify=ctx.no_verify,
            smooth=False,
        )
        bas = run_stl2step(
            ctx.baseline_bin,
            fixture.stl,
            bas_step,
            no_verify=ctx.no_verify,
            smooth=False,
        )
    except RuntimeError as exc:
        return go(ctx, gate_id, "FAIL", str(exc))

    cur_txt.write_text(cur.stdout, encoding="utf-8")
    bas_txt.write_text(bas.stdout, encoding="utf-8")

    if not cur.result.get("ok"):
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"off-path conversion failed: {cur.result.get('error')}",
        )
    if not bas.result.get("ok"):
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"187ead0 baseline conversion failed: {bas.result.get('error')}",
        )

    smooth_keys = result_has_smooth_keys(cur.result)
    if smooth_keys:
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"smooth* / facesAfterSmooth keys present with --smooth absent: {smooth_keys}",
        )

    step_cmp = run_canonicalize(canon, "step", cur_step, bas_step)
    if step_cmp.returncode == 2:
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"canonicalize.py step IO/usage error: {step_cmp.stderr or step_cmp.stdout}",
        )
    if step_cmp.returncode != 0:
        diff = (step_cmp.stdout or "") + (step_cmp.stderr or "")
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonicalize.py step: DATA section differs from 187ead0",
            details={"canonicalize": diff[:4000]},
        )

    res_cmp = run_canonicalize(canon, "result", cur_txt, bas_txt)
    if res_cmp.returncode == 2:
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"canonicalize.py result IO/usage error: {res_cmp.stderr or res_cmp.stdout}",
        )
    if res_cmp.returncode != 0:
        diff = (res_cmp.stdout or "") + (res_cmp.stderr or "")
        return go(
            ctx,
            gate_id,
            "FAIL",
            "canonicalize.py result: RESULT differs from 187ead0",
            details={"canonicalize": diff[:4000]},
        )

    return go(
        ctx,
        gate_id,
        "PASS",
        "off-path twin-run vs 187ead0: canonicalize.py step+result IDENTICAL",
    )


def smooth_ctx(ctx: GateContext) -> so.GateContext:
    return so.GateContext(
        binary=ctx.binary,
        fixture=so.Fixture(ctx.fixture.id, ctx.fixture.stl, ctx.fixture.sidecar),
        work_dir=ctx.work_dir,
        census_path=ctx.census_path,
        dump_path=ctx.dump_path,
        no_verify=ctx.no_verify,
        baseline_bin=ctx.baseline_bin,
    )


def outcome_from_so(o: so.GateOutcome) -> GateOutcome:
    return GateOutcome(
        o.gate_id,
        o.fixture_id,
        o.status,
        o.message,
        hard=o.hard,
        details=o.details or {},
    )


def gate_g0_2_refit_closedness(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.2 — refit never reduces closedness (solids/openShells)."""
    if "G0.2" not in LIVE_GATES:
        return xfail_not_landed("G0.2", ctx.fixture.id, "P3-engine/smooth")
    return outcome_from_so(so.check_g0_2(smooth_ctx(ctx)))


def gate_g0_3_dirty_skip(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G0.3 — dirty skip on S15; force-sew global disable."""
    if "G0.3" not in LIVE_GATES:
        return xfail_not_landed("G0.3", ctx.fixture.id, "P3-engine/smooth")
    return outcome_from_so(so.check_g0_3(smooth_ctx(ctx)))


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
    shapefix_note = bool(
        re.search(r"\bfix\b.*ShapeFix|ShapeFix.*\bfix\b", combined, re.IGNORECASE)
    )
    shapefix_warn = [
        w
        for w in run.result.get("warnings", [])
        if "ShapeFix" in w or re.search(r"\bfix\b", w, re.IGNORECASE)
    ]

    # Census is the G1 authority (CTest default path passes --census).
    # 1.0.0 itself rewrites S09 via ShapeFix; G0.1 already covers off-path
    # identity. When the witness says the written file is BRepCheck-valid,
    # that is G1 PASS. The fix-note is the no-census detector only.
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
            extra = ""
            if shapefix_note or shapefix_warn:
                extra = " (ShapeFix rewrite noted)"
            return go(
                ctx,
                gate_id,
                "FAIL",
                f"BRepCheck invalid per census{extra}: {census}",
            )
        return go(ctx, gate_id, "PASS", "BRepCheck valid (census)")

    # Census genuinely unavailable: do not proxy via watertight/openShells
    # (that hard-fails open-but-valid B-Reps). SKIP the BRepCheck sub-check.
    if shapefix_note:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "ShapeFix invalid-rewrite path fired (fix note in CLI output; census unavailable)",
        )
    if shapefix_warn:
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"ShapeFix-related warning (census unavailable): {shapefix_warn[0]}",
        )
    return go(
        ctx,
        gate_id,
        "SKIP",
        "G1 BRepCheck SKIPPED: census binary not provided (no watertight proxy)",
        hard=False,
    )


def gate_g2_recognition(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G2 — recoverable primitive recognition vs sidecar."""
    if "G2" not in LIVE_GATES:
        return xfail_not_landed("G2", ctx.fixture.id, "P2-build + p0-census")
    return outcome_from_so(so.check_g2(smooth_ctx(ctx)))


def gate_g2_5_built_as(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G2.5 — closed360 holes built as one cylindrical face."""
    if "G2.5" not in LIVE_GATES:
        return xfail_not_landed("G2.5", ctx.fixture.id, "P2-build + p0-census")
    return outcome_from_so(so.check_g2_5(smooth_ctx(ctx)))


def gate_g3_sagitta_volume(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G3 — chord-aware D4.5 volume budget (HARD)."""
    if "G3" not in LIVE_GATES:
        return xfail_not_landed("G3", ctx.fixture.id, "P3-engine/smooth")
    return outcome_from_so(so.check_g3(smooth_ctx(ctx)))


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
    if "G4.4" not in LIVE_GATES:
        return xfail_not_landed("G4.4", ctx.fixture.id, "P3-engine/release")
    return outcome_from_so(so.check_g4_4(smooth_ctx(ctx)))


def gate_g5_editability(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 G5 — editability censuses; smoothMaxEdgeTolMM threshold (SOFT)."""
    if "G5" not in LIVE_GATES:
        return xfail_not_landed("G5", ctx.fixture.id, "P3-engine/smooth")
    return outcome_from_so(so.check_g5(smooth_ctx(ctx)))


def gate_i_checker(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 I-checker — RegionSet dump validation (check_regionset.py).

    LIVE since the wave-4 staging assembly: p1-compose-fix landed the
    envelope->bare unwrap (`stl2step_regiondump --component N --bare --out`),
    so the runner drives the dump per CLEAN component and validates each bare
    RegionSet with check_regionset.py, including the fixture sidecar where one
    exists. A fixture with no clean components (S15, S16-R1-round-2) is SKIP,
    not PASS — I6 says P1 is never called on a dirty component.
    """
    gate_id = "I-checker"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "p0-ichecker + P1 dump")
    if not (ctx.dump_path and ctx.dump_path.is_file()):
        return go(ctx, gate_id, "XFAIL",
                  "I-checker: --dump <stl2step_regiondump> not supplied or missing; "
                  "build the target and pass --dump", hard=True)
    if not (ctx.ichecker_path and ctx.ichecker_path.is_file()):
        return go(ctx, gate_id, "XFAIL",
                  "I-checker: --ichecker <check_regionset.py> not supplied or missing",
                  hard=True)

    stl = ctx.fixture.stl
    envelope = subprocess.run([str(ctx.dump_path), str(stl)],
                              check=False, capture_output=True, text=True)
    if envelope.returncode != 0:
        return go(ctx, gate_id, "FAIL",
                  f"regiondump failed on {stl.name}: {envelope.stderr.strip()[:400]}")
    try:
        doc = json.loads(envelope.stdout)
    except json.JSONDecodeError as exc:
        return go(ctx, gate_id, "FAIL", f"regiondump emitted invalid JSON on {stl.name}: {exc}")

    comps = [int(c["index"]) for c in doc.get("comps", []) if c.get("clean")]
    if not comps:
        return go(ctx, gate_id, "SKIP",
                  f"{stl.name}: no clean components (I6 — P1 is never called on a dirty component)",
                  hard=True)

    sidecar = ctx.fixture.sidecar if getattr(ctx.fixture, "sidecar", None) else None
    failures: List[str] = []
    for comp in comps:
        bare = ctx.work_dir / f"{ctx.fixture.id}-comp{comp}.regionset.json"
        made = subprocess.run(
            [str(ctx.dump_path), str(stl), "--component", str(comp), "--bare",
             "--out", str(bare)],
            check=False, capture_output=True, text=True)
        if made.returncode != 0 or not bare.is_file():
            failures.append(f"comp{comp}: --bare dump failed: {made.stderr.strip()[:200]}")
            continue
        cmd = [sys.executable, str(ctx.ichecker_path), str(bare)]
        if sidecar and Path(sidecar).is_file():
            cmd += ["--sidecar", str(sidecar)]
        chk = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if chk.returncode != 0:
            detail = (chk.stdout.strip() or chk.stderr.strip())[:400]
            failures.append(f"comp{comp}: {detail}")

    if failures:
        return go(ctx, gate_id, "FAIL",
                  f"I-checker: {len(failures)}/{len(comps)} component(s) failed — "
                  + " | ".join(failures),
                  details={"components": comps, "failures": failures})
    return go(ctx, gate_id, "PASS",
              f"I-checker: I1-I9 (I7/I7b separately) hold on {len(comps)} clean "
              f"component(s) of {stl.name}",
              details={"components": comps})


def gate_include_allowlist(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 include-allowlist grep (DECISION D5.3) on all P1 TUs + refit_internal.hpp."""
    gate_id = "include-allowlist"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "P1-segment")

    missing = [str(p) for p in INCLUDE_ALLOWLIST_FILES if not p.is_file()]
    existing = [p for p in INCLUDE_ALLOWLIST_FILES if p.is_file()]
    if missing:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "include-allowlist files missing — refusing vacuous pass",
            details={"missing": missing, "checked": [str(p) for p in existing]},
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

    checked = [str(p) for p in existing]
    if violations:
        return go(
            ctx,
            gate_id,
            "FAIL",
            "banned or non-allowlisted #include in P1 sources",
            details={"violations": violations, "checked": checked},
        )

    return go(
        ctx,
        gate_id,
        "PASS",
        "all #includes within D5.3 allowlist: " + ", ".join(p.name for p in existing),
        details={"checked": checked},
    )


def gate_calibration(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 calibration — freeze G3 window and volume k at P0 close."""
    if "calibration" not in LIVE_GATES:
        return xfail_not_landed("calibration", ctx.fixture.id, "P0-corpus-close")
    return outcome_from_so(so.check_calibration(smooth_ctx(ctx)))


def gate_r_ladder(ctx: GateContext) -> GateOutcome:
    """SPEC-P0 R-ladder — S16 R1-explode, R1-round-2, R2 revert fixtures."""
    if "R-ladder" not in LIVE_GATES:
        return xfail_not_landed("R-ladder", ctx.fixture.id, "P0-corpus S16 + P2-build")
    return outcome_from_so(so.check_r_ladder(smooth_ctx(ctx)))


def gate_g_law(ctx: GateContext) -> GateOutcome:
    """G-LAW — handle-lock supervised recognition (regiondump vs tri-labels)."""
    gate_id = "G-LAW"
    if gate_id not in LIVE_GATES:
        return xfail_not_landed(gate_id, ctx.fixture.id, "ac2-l4")
    if ctx.fixture.id != "handle-lock":
        return go(
            ctx,
            gate_id,
            "SKIP",
            "G-LAW: handle-lock only (single-fixture supervised gate)",
            hard=True,
        )
    if not (ctx.dump_path and ctx.dump_path.is_file()):
        return go(
            ctx,
            gate_id,
            "XFAIL",
            "G-LAW: --dump <stl2step_regiondump> not supplied or missing",
            hard=True,
        )
    script = REPO_ROOT / "tests" / "gates" / "law_recognition.py"
    if not script.is_file():
        return go(ctx, gate_id, "FAIL", f"G-LAW driver missing: {script}")

    bare = ctx.work_dir / f"{ctx.fixture.id}.regionset.json"
    made = subprocess.run(
        [
            str(ctx.dump_path),
            str(ctx.fixture.stl),
            "--component",
            "0",
            "--bare",
            "--out",
            str(bare),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if made.returncode != 0 or not bare.is_file():
        return go(
            ctx,
            gate_id,
            "FAIL",
            f"G-LAW regiondump failed: {(made.stderr or made.stdout).strip()[:400]}",
        )

    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--dump",
            str(ctx.dump_path),
            "--regiondump",
            str(bare),
            "--stl",
            str(ctx.fixture.stl),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        detail = (proc.stdout + proc.stderr).strip()[-2000:]
        return go(ctx, gate_id, "FAIL", f"G-LAW: {detail}", hard=True)
    return go(
        ctx,
        gate_id,
        "PASS",
        "G-LAW: band recall / purity / radius ratchet + RULE 4.2a PASS",
        hard=True,
    )


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
    "G-LAW": gate_g_law,
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
    dump_path: Optional[Path],
    ichecker_path: Optional[Path],
    no_verify: bool,
    baseline_bin: Optional[Path] = None,
    baseline_error: Optional[str] = None,
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
            dump_path=dump_path,
            ichecker_path=ichecker_path,
            no_verify=no_verify,
            baseline_bin=baseline_bin,
            baseline_error=baseline_error,
        )
        for gate_id in gate_ids:
            fn = GATE_REGISTRY[gate_id]
            outcomes.append(fn(ctx))
    return outcomes


def is_hard_failure(o: GateOutcome) -> bool:
    """True when an outcome should fail the process (excludes SOFT gates and PARKED)."""
    if o.status != "FAIL":
        return False
    if o.gate_id in SOFT_GATES or not o.hard:
        return False
    return True


def summarize(outcomes: List[GateOutcome]) -> str:
    lines = ["gates summary:"]
    width = max(len(o.gate_id) for o in outcomes) if outcomes else 10
    for o in outcomes:
        lines.append(f"  {o.gate_id:<{width}}  {o.status:<7}  {o.message}")
    counts: Dict[str, int] = {}
    for o in outcomes:
        counts[o.status] = counts.get(o.status, 0) + 1
    soft_fail = sum(
        1
        for o in outcomes
        if o.status == "FAIL" and (o.gate_id in SOFT_GATES or not o.hard)
    )
    hard_fail = sum(1 for o in outcomes if is_hard_failure(o))
    lines.append(
        f"totals: PASS={counts.get('PASS', 0)} "
        f"FAIL(hard)={hard_fail} FAIL(soft)={soft_fail} "
        f"PARKED={counts.get('PARKED', 0)} "
        f"XFAIL={counts.get('XFAIL', 0)} SKIP={counts.get('SKIP', 0)}"
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
            "FAIL_hard": sum(1 for o in outcomes if is_hard_failure(o)),
            "FAIL_soft": sum(
                1
                for o in outcomes
                if o.status == "FAIL" and (o.gate_id in SOFT_GATES or not o.hard)
            ),
            "PARKED": sum(1 for o in outcomes if o.status == "PARKED"),
            "XFAIL": sum(1 for o in outcomes if o.status == "XFAIL"),
            "SKIP": sum(1 for o in outcomes if o.status == "SKIP"),
        },
        "parkedGates": dict(PARKED_GATES),
        "parkedGateFixtures": {
            f"{g}@{f}": reason for (g, f), reason in so.PARKED_GATE_FIXTURES.items()
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
    unpark = set(split_csv(args.unpark) or [])
    gate_ids = selected_gates(gate_filter)
    fixtures = discover_fixtures(args.corpus, args.smoke, fixture_filter)
    if not fixtures:
        print("error: no fixtures selected", file=sys.stderr)
        return 1

    baseline_dir = args.baseline.resolve() if args.baseline else None
    census_path = args.census.resolve() if args.census else None
    dump_path = args.dump.resolve() if args.dump else None
    ichecker_path = (args.ichecker.resolve() if args.ichecker
                     else (REPO_ROOT / "tests/gates/check_regionset.py"))

    baseline_bin: Optional[Path] = None
    baseline_error: Optional[str] = None
    if "G0.1" in gate_ids:
        baseline_bin, baseline_error = ensure_baseline(
            baseline_dir, current_build_from_binary(binary)
        )

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
                    dump_path,
                    ichecker_path,
                    args.no_verify,
                    baseline_bin,
                    baseline_error,
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
                    dump_path,
                    ichecker_path,
                    args.no_verify,
                    baseline_bin,
                    baseline_error,
                )
            )

    all_outcomes = so.apply_parking(all_outcomes, unpark)

    report = build_report(fixtures, gate_ids, all_outcomes, args)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(summarize(all_outcomes))
    print(so.format_gate_table(all_outcomes))
    print(f"report: {args.report}")

    hard_failures = [o for o in all_outcomes if is_hard_failure(o)]
    return 1 if hard_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
