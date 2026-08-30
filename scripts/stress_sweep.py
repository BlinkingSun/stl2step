#!/usr/bin/env python3
"""TrueForm 11-file stress-sweep harness (stdlib only).

Replaces hand-rolled STRESS-SWEEP tables. One stl2step process per file;
files run in parallel up to ``-j``. Private STLs are skip-if-missing via ``STL2STEP_PRIVATE_CORPUS``.
Unset/empty env or a missing corpus directory → exit 77.

Exit codes:
  0  — table written (every present file produced a RESULT)
  1  — hard failure (binary missing, handle-lock missing, no RESULT)
  77 — private corpus absent (ctest SKIP_RETURN_CODE)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
PRIVATE_ENV = "STL2STEP_PRIVATE_CORPUS"
PRIVATE_BODY_IDS: Tuple[int, ...] = (1, 2, 9, 10, 11, 12, 13, 18, 20, 28)
SKIP_RC = 77

# volΔ percentage-point tolerance vs hand measurements (AC3).
VOL_TOL_PP = 0.05

# STRESS-SWEEP-2.md (hl/sweep2 @ 369173a) — drift source, not main truth.
SWEEP2_HAND: Dict[str, Dict[str, Any]] = {
    "Body12.stl": {
        "tris": 7918,
        "exit": 2,
        "watertight": True,
        "built": 0,
        "reverted": 1,
        "segCyl": 201,
        "builtCyl": 0,
        "censusCyl": 0,
        "volumeDeltaPct": 0.0,
        "j6FreeEdges": 14,
    },
    "handle-lock.stl": {
        "tris": 908,
        "exit": 0,
        "watertight": True,
        "built": 1,
        "reverted": 0,
        "segCyl": 16,
        "builtCyl": 1,
        "censusCyl": 1,
        "volumeDeltaPct": 0.0,
        "j6FreeEdges": 0,
    },
    "Body9.stl": {
        "tris": 4668,
        "exit": 2,
        "watertight": True,
        "built": 0,
        "reverted": 1,
        "segCyl": 183,
        "builtCyl": 0,
        "censusCyl": 0,
        "volumeDeltaPct": 0.031,
        "j6FreeEdges": 103,
    },
}

# spike-diag.md main @ fd4954d — current hand truth for AC3.
SPIKE_HAND: Dict[str, Dict[str, Any]] = {
    "Body12.stl": {
        "tris": 7918,
        "exit": 2,
        "watertight": True,
        "built": 0,
        "reverted": 1,
        "segCyl": 201,
        "builtCyl": 0,
        "censusCyl": 0,
        "volumeDeltaPct": 0.0,
        "j6FreeEdges": 81,
    },
    "handle-lock.stl": {
        "tris": 908,
        "exit": 0,
        "watertight": True,
        "built": 1,
        "reverted": 0,
        "segCyl": 16,
        "builtCyl": 1,
        "censusCyl": 1,
        "volumeDeltaPct": 0.0,
        "j6FreeEdges": 0,
    },
    "Body9.stl": {
        "tris": 4668,
        "exit": 2,
        "watertight": True,
        "built": 0,
        "reverted": 1,
        "segCyl": 183,
        "builtCyl": 0,
        "censusCyl": 0,
        "volumeDeltaPct": 0.030937,
        "j6FreeEdges": 34,
    },
}

_MAKEEDGE_RE = re.compile(r"analytic MakeEdge failed")
_INTANA_CYL_RE = re.compile(r"IntAna cyl\|cyl empty/same")
_INTANA_PLN_RE = re.compile(r"IntAna plane\|cyl empty/same")
_SEAMED_RE = re.compile(r"seamed360:")
_J6_RE = re.compile(r"J6: shell not closed freeEdges=(\d+)")
_REVERT_RE = re.compile(r"analytic rebuild reverted")
_RESULT_RE = re.compile(r"^RESULT\s+(\{.*\})\s*$")


@dataclass
class CorpusFile:
    name: str
    path: Path
    private: bool


@dataclass
class WarningDigest:
    make_edge: int = 0
    intana_cyl: int = 0
    intana_plane: int = 0
    seamed360: int = 0
    j6_hits: int = 0
    j6_free_edges: int = 0
    explicit_revert: int = 0

    def format(self) -> str:
        parts: List[str] = []
        if self.make_edge:
            parts.append(f"MakeEdge×{self.make_edge}")
        if self.intana_cyl:
            parts.append(f"IntAna cyl|cyl×{self.intana_cyl}")
        if self.intana_plane:
            parts.append(f"plane|cyl×{self.intana_plane}")
        if self.seamed360:
            parts.append(f"seamed360×{self.seamed360}")
        if self.j6_hits:
            parts.append(f"J6 freeEdges={self.j6_free_edges}")
        if self.explicit_revert:
            parts.append("explicit revert")
        return "; ".join(parts) if parts else "(none)"


@dataclass
class SweepRow:
    file: str
    path: str
    skipped: bool = False
    skip_reason: str = ""
    tris: int = 0
    exit: int = -1
    watertight: bool = False
    built: int = 0
    reverted: int = 0
    segCyl: int = 0
    builtCyl: int = 0
    censusCyl: int = 0
    volumeDeltaPct: float = -1.0
    warningsDigest: str = ""
    warningCounts: Dict[str, int] = field(default_factory=dict)
    seconds: float = 0.0
    error: str = ""

    def to_json(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SweepReport:
    ok: bool
    jobs: int
    seconds: float
    privateCorpus: str
    handleLock: str
    skipped: List[str] = field(default_factory=list)
    rows: List[SweepRow] = field(default_factory=list)
    error: str = ""

    def to_json(self) -> Dict[str, Any]:
        return {
            "ok": self.ok,
            "jobs": self.jobs,
            "seconds": round(self.seconds, 3),
            "privateCorpus": self.privateCorpus,
            "handleLock": self.handleLock,
            "skipped": list(self.skipped),
            "error": self.error,
            "rows": [r.to_json() for r in self.rows],
        }


def resolve_private_corpus(
    env: Optional[Mapping[str, str]] = None,
    override: Optional[Path] = None,
) -> Path:
    """Resolve the private STL directory. Unset/empty env → empty path (SKIP)."""
    if override is not None:
        return Path(override).expanduser()
    src = os.environ if env is None else env
    raw = (src.get(PRIVATE_ENV) or "").strip()
    if not raw:
        return Path()
    return Path(raw).expanduser()


def resolve_handle_lock(repo: Path = REPO_ROOT, override: Optional[Path] = None) -> Path:
    if override is not None:
        return Path(override)
    corpus = repo / "tests" / "corpus" / "handle-lock.stl"
    if corpus.is_file():
        return corpus
    return repo / "tests" / "diag" / "handle-lock" / "handle-lock.stl"


def list_corpus_files(
    private_root: Path,
    handle_lock: Path,
) -> List[CorpusFile]:
    """Deterministic 11-file order: Body1, Body2, Body9…Body28, then handle-lock."""
    files = [
        CorpusFile(f"Body{n}.stl", private_root / f"Body{n}.stl", True)
        for n in PRIVATE_BODY_IDS
    ]
    files.append(CorpusFile("handle-lock.stl", handle_lock, False))
    return files


def private_corpus_present(private_root: Path) -> bool:
    return any((private_root / f"Body{n}.stl").is_file() for n in PRIVATE_BODY_IDS)


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
        m = _RESULT_RE.match(line)
        if m:
            return json.loads(m.group(1))
    raise RuntimeError("no RESULT line in stl2step stdout")


def digest_warnings(stderr: str, result_warnings: Optional[Sequence[str]] = None) -> WarningDigest:
    # Hand tables count stderr lines (emit(warn) already mirrors into RESULT).
    blob = stderr or ""
    if not blob and result_warnings:
        blob = "\n".join(str(w) for w in result_warnings)
    j6_vals = [int(m.group(1)) for m in _J6_RE.finditer(blob)]
    return WarningDigest(
        make_edge=len(_MAKEEDGE_RE.findall(blob)),
        intana_cyl=len(_INTANA_CYL_RE.findall(blob)),
        intana_plane=len(_INTANA_PLN_RE.findall(blob)),
        seamed360=len(_SEAMED_RE.findall(blob)),
        j6_hits=len(j6_vals),
        j6_free_edges=j6_vals[-1] if j6_vals else 0,
        explicit_revert=len(_REVERT_RE.findall(blob)),
    )


def _import_census():
    tools = REPO_ROOT / "tests" / "tools"
    if str(tools) not in sys.path:
        sys.path.insert(0, str(tools))
    import step_census  # type: ignore

    return step_census


def census_cylinders(step_path: Path) -> int:
    census = _import_census()
    doc = census.census_path(step_path)
    if not doc.get("ok"):
        raise RuntimeError(doc.get("error") or f"census failed for {step_path}")
    return int((doc.get("surfaces") or {}).get("cylinder") or 0)


def run_one(
    item: CorpusFile,
    binary: Path,
    out_dir: Path,
    *,
    threads: int = 1,
) -> SweepRow:
    if not item.path.is_file():
        return SweepRow(
            file=item.name,
            path=str(item.path),
            skipped=True,
            skip_reason="missing",
        )

    step_path = out_dir / (item.path.stem + ".step")
    cmd = [
        str(binary),
        str(item.path),
        "-o",
        str(step_path),
        "--engine",
        "trueform",
        "--quiet",
        "--threads",
        str(threads),
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    try:
        result = parse_result_line(proc.stdout)
    except Exception as exc:
        return SweepRow(
            file=item.name,
            path=str(item.path),
            exit=proc.returncode,
            seconds=elapsed,
            error=f"no RESULT: {exc}",
        )

    digest = digest_warnings(proc.stderr or "", result.get("warnings") or [])
    census_cyl = 0
    if step_path.is_file():
        try:
            census_cyl = census_cylinders(step_path)
        except Exception as exc:
            return SweepRow(
                file=item.name,
                path=str(item.path),
                tris=int(result.get("triangles") or 0),
                exit=proc.returncode,
                seconds=elapsed,
                error=f"census failed: {exc}",
            )

    vol = result.get("volumeDeltaPct", -1)
    try:
        vol_f = float(vol)
    except (TypeError, ValueError):
        vol_f = -1.0

    return SweepRow(
        file=item.name,
        path=str(item.path),
        tris=int(result.get("triangles") or 0),
        exit=int(proc.returncode),
        watertight=bool(result.get("watertight", False)),
        built=int(result.get("smoothBuiltComponents") or 0),
        reverted=int(result.get("smoothRevertedComponents") or 0),
        segCyl=int(result.get("smoothCylinders") or 0),
        builtCyl=int(result.get("smoothBuiltCylinders") or 0),
        censusCyl=census_cyl,
        volumeDeltaPct=vol_f,
        warningsDigest=digest.format(),
        warningCounts={
            "makeEdge": digest.make_edge,
            "intanaCyl": digest.intana_cyl,
            "intanaPlane": digest.intana_plane,
            "seamed360": digest.seamed360,
            "j6Hits": digest.j6_hits,
            "j6FreeEdges": digest.j6_free_edges,
            "explicitRevert": digest.explicit_revert,
        },
        seconds=elapsed,
    )


def run_sweep(
    files: Sequence[CorpusFile],
    binary: Path,
    out_dir: Path,
    *,
    jobs: int = 4,
    threads_per_file: int = 1,
) -> List[SweepRow]:
    """Convert files in parallel (max ``jobs`` processes). Output order is input order."""
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = max(1, int(jobs))
    indexed = list(enumerate(files))
    rows: List[Optional[SweepRow]] = [None] * len(files)

    if workers == 1 or len(files) <= 1:
        for i, item in indexed:
            rows[i] = run_one(item, binary, out_dir, threads=threads_per_file)
        return [r for r in rows if r is not None]

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(run_one, item, binary, out_dir, threads=threads_per_file): i
            for i, item in indexed
        }
        for fut in as_completed(futs):
            rows[futs[fut]] = fut.result()
    return [r if r is not None else SweepRow(file="?", path="", error="missing row") for r in rows]


def format_markdown(report: SweepReport) -> str:
    lines = [
        "| File | Tris | Exit | Watertight | Built / Reverted | "
        "Seg cyl → built cyl → census cyl | volΔ% | Warnings digest |",
        "|------|-----:|-----:|:----------:|:----------------:|"
        "----------------------------------:|------:|-----------------|",
    ]
    for r in report.rows:
        if r.skipped:
            lines.append(
                f"| {r.file} | — | — | — | — | — | — | skipped ({r.skip_reason}) |"
            )
            continue
        wt = "yes" if r.watertight else "no"
        vol = f"{r.volumeDeltaPct:.3f}" if r.volumeDeltaPct >= 0 else "—"
        lines.append(
            f"| {r.file} | {r.tris} | {r.exit} | {wt} | "
            f"{r.built} / {r.reverted} | "
            f"{r.segCyl} → {r.builtCyl} → {r.censusCyl} | "
            f"{vol} | {r.warningsDigest or '(none)'} |"
        )
    lines.append("")
    lines.append(f"jobs={report.jobs}  wall_s={report.seconds:.2f}  ok={report.ok}")
    if report.skipped:
        lines.append("skipped: " + ", ".join(report.skipped))
    return "\n".join(lines) + "\n"


def _row_metrics(row: SweepRow) -> Dict[str, Any]:
    return {
        "tris": row.tris,
        "exit": row.exit,
        "watertight": row.watertight,
        "built": row.built,
        "reverted": row.reverted,
        "segCyl": row.segCyl,
        "builtCyl": row.builtCyl,
        "censusCyl": row.censusCyl,
        "volumeDeltaPct": row.volumeDeltaPct,
        "j6FreeEdges": int((row.warningCounts or {}).get("j6FreeEdges") or 0),
    }


def compare_row(
    row: SweepRow,
    expected: Mapping[str, Any],
    *,
    vol_tol: float = VOL_TOL_PP,
) -> List[str]:
    """Return mismatch strings. Counters must match exactly; volΔ within vol_tol pp."""
    if row.skipped:
        return [f"{row.file}: skipped ({row.skip_reason})"]
    if row.error:
        return [f"{row.file}: {row.error}"]
    got = _row_metrics(row)
    drifts: List[str] = []
    for key in (
        "tris",
        "exit",
        "watertight",
        "built",
        "reverted",
        "segCyl",
        "builtCyl",
        "censusCyl",
    ):
        if key in expected and got[key] != expected[key]:
            drifts.append(f"{row.file} {key}: got {got[key]} want {expected[key]}")
    if "volumeDeltaPct" in expected:
        want = float(expected["volumeDeltaPct"])
        if abs(got["volumeDeltaPct"] - want) > vol_tol:
            drifts.append(
                f"{row.file} volumeDeltaPct: got {got['volumeDeltaPct']} "
                f"want {want} ±{vol_tol}"
            )
    if "j6FreeEdges" in expected and got["j6FreeEdges"] != expected["j6FreeEdges"]:
        drifts.append(
            f"{row.file} j6FreeEdges: got {got['j6FreeEdges']} "
            f"want {expected['j6FreeEdges']}"
        )
    return drifts


def compare_report(
    report: SweepReport,
    expected_by_file: Mapping[str, Mapping[str, Any]],
    *,
    vol_tol: float = VOL_TOL_PP,
) -> List[str]:
    by_name = {r.file: r for r in report.rows}
    drifts: List[str] = []
    for name, exp in expected_by_file.items():
        row = by_name.get(name)
        if row is None:
            drifts.append(f"{name}: not in sweep")
            continue
        drifts.extend(compare_row(row, exp, vol_tol=vol_tol))
    return drifts


def find_binary(explicit: Optional[Path] = None, repo: Path = REPO_ROOT) -> Path:
    if explicit is not None:
        return Path(explicit)
    for cand in (
        repo / "build" / "stl2step",
        repo / "build" / "stl2step.exe",
        repo / "build" / "Release" / "stl2step.exe",
    ):
        if cand.is_file():
            return cand
    return Path("stl2step")


def _self_test() -> int:
    """Exercise the harness API without the engine or private STLs."""
    errors: List[str] = []

    result = parse_result_line(
        'progress\nRESULT {"ok":true,"triangles":908,"watertight":true,'
        '"smoothCylinders":16,"smoothBuiltCylinders":1,'
        '"smoothBuiltComponents":1,"smoothRevertedComponents":0,'
        '"volumeDeltaPct":0.0,"warnings":[]}\n'
    )
    if result.get("smoothBuiltCylinders") != 1 or result.get("triangles") != 908:
        errors.append(f"parse_result_line: {result}")

    stderr = (
        "smooth: analytic MakeEdge failed — keeping mesh polyline\n"
        "smooth: analytic MakeEdge failed — keeping mesh polyline\n"
        "smooth: IntAna cyl|cyl empty/same — keeping mesh polyline\n"
        "smooth: IntAna plane|cyl empty/same — keeping mesh polyline\n"
        "seamed360: BRepCheck invalid on seamed face\n"
        "J6: shell not closed freeEdges=103 faces=2995 recover=0\n"
        "J6: shell not closed freeEdges=81 faces=3799 recover=1\n"
        "smooth: analytic rebuild reverted on one component -- kept faceted\n"
    )
    d = digest_warnings(stderr)
    if d.make_edge != 2 or d.intana_cyl != 1 or d.intana_plane != 1:
        errors.append(f"digest counts: {d}")
    if d.seamed360 != 1 or d.j6_hits != 2 or d.j6_free_edges != 81:
        errors.append(f"digest j6/seamed: {d}")
    if d.explicit_revert != 1:
        errors.append(f"digest revert: {d}")
    want_txt = (
        "MakeEdge×2; IntAna cyl|cyl×1; plane|cyl×1; seamed360×1; "
        "J6 freeEdges=81; explicit revert"
    )
    if d.format() != want_txt:
        errors.append(f"digest format: {d.format()!r}")

    empty = Path("/no/such/stl2step-private-corpus")
    if private_corpus_present(empty):
        errors.append("empty corpus reported present")
    unset = resolve_private_corpus(env={}, override=None)
    if unset != Path() and str(unset) != ".":
        errors.append("unset env did not yield empty path")
    forced = resolve_private_corpus(env={PRIVATE_ENV: str(empty)})
    if forced != empty:
        errors.append(f"env override ignored: {forced}")

    hl = Path("/tmp/handle-lock.stl")
    names = [f.name for f in list_corpus_files(empty, hl)]
    expect_names = [f"Body{n}.stl" for n in PRIVATE_BODY_IDS] + ["handle-lock.stl"]
    if names != expect_names:
        errors.append(f"corpus order: {names}")

    fake_hl = SweepRow(
        file="handle-lock.stl",
        path="x",
        tris=908,
        exit=0,
        watertight=True,
        built=1,
        reverted=0,
        segCyl=16,
        builtCyl=1,
        censusCyl=1,
        volumeDeltaPct=0.0,
        warningCounts={"j6FreeEdges": 0},
        warningsDigest="(none)",
    )
    fake_b12 = SweepRow(
        file="Body12.stl",
        path="y",
        tris=7918,
        exit=2,
        watertight=True,
        built=0,
        reverted=1,
        segCyl=201,
        builtCyl=0,
        censusCyl=0,
        volumeDeltaPct=0.0,
        warningCounts={"j6FreeEdges": 81},
        warningsDigest="J6 freeEdges=81",
    )
    fake_b9 = SweepRow(
        file="Body9.stl",
        path="z",
        tris=4668,
        exit=2,
        watertight=True,
        built=0,
        reverted=1,
        segCyl=183,
        builtCyl=0,
        censusCyl=0,
        volumeDeltaPct=0.030937,
        warningCounts={"j6FreeEdges": 34},
        warningsDigest="J6 freeEdges=34",
    )
    fake = SweepReport(
        ok=True,
        jobs=4,
        seconds=1.25,
        privateCorpus=str(empty),
        handleLock=str(hl),
        rows=[fake_b9, fake_b12, fake_hl],
    )
    md = format_markdown(fake)
    for col in (
        "Tris",
        "Exit",
        "Watertight",
        "Built / Reverted",
        "Seg cyl → built cyl → census cyl",
        "volΔ%",
        "Warnings digest",
    ):
        if col not in md:
            errors.append(f"markdown missing column {col!r}")
    if "16 → 1 → 1" not in md or "201 → 0 → 0" not in md:
        errors.append(f"markdown missing cylinder chain:\n{md}")

    spike_drifts = compare_report(fake, SPIKE_HAND)
    if spike_drifts:
        errors.append("compare vs spike-diag should be clean: " + "; ".join(spike_drifts))
    sweep2_drifts = compare_report(fake, SWEEP2_HAND)
    if not any("j6FreeEdges" in x and "Body12" in x for x in sweep2_drifts):
        errors.append(f"expected Body12 J6 drift vs sweep2, got {sweep2_drifts}")

    mutated = SweepRow(**{**fake_hl.to_json(), "volumeDeltaPct": 0.2})
    vol_drifts = compare_row(mutated, SPIKE_HAND["handle-lock.stl"])
    if not vol_drifts:
        errors.append("volume tolerance failed to fire at +0.2 pp")

    def _square(n: int) -> int:
        return n * n

    with ThreadPoolExecutor(max_workers=4) as pool:
        got = list(pool.map(_square, range(8)))
    if got != [i * i for i in range(8)]:
        errors.append(f"thread pool: {got}")

    if errors:
        print("SELF-TEST FAIL", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        return 1
    print("SELF-TEST PASS")
    return 0


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, help="stl2step executable")
    p.add_argument("-j", "--jobs", type=int, default=4, help="parallel files (default 4)")
    p.add_argument("--threads", type=int, default=1, help="engine threads per file")
    p.add_argument("--private-corpus", type=Path, help=f"override ${PRIVATE_ENV}")
    p.add_argument("--handle-lock", type=Path, help="override handle-lock.stl path")
    p.add_argument("--out-dir", type=Path, help="directory for STEP + table outputs")
    p.add_argument("--json", type=Path, help="write JSON table")
    p.add_argument("--markdown", type=Path, help="write markdown table")
    p.add_argument(
        "--compare",
        action="store_true",
        help="print drift vs STRESS-SWEEP-2 and spike-diag hand measurements",
    )
    p.add_argument("--self-test", action="store_true", help="exercise API; no engine")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    if args.self_test:
        return _self_test()

    private = resolve_private_corpus(override=args.private_corpus)
    handle_lock = resolve_handle_lock(override=args.handle_lock)
    files = list_corpus_files(private, handle_lock)

    if not handle_lock.is_file():
        print(f"error: handle-lock missing: {handle_lock}", file=sys.stderr)
        return 1
    if not private_corpus_present(private):
        print(
            f"SKIP: private corpus absent ({PRIVATE_ENV}={private})",
            file=sys.stderr,
        )
        return SKIP_RC

    binary = find_binary(args.binary)
    if not Path(binary).is_file():
        print(f"error: stl2step binary not found: {binary}", file=sys.stderr)
        return 1

    out_dir = args.out_dir or Path(tempfile.mkdtemp(prefix="stress-sweep-"))
    t0 = time.perf_counter()
    rows = run_sweep(
        files,
        Path(binary),
        Path(out_dir),
        jobs=max(1, args.jobs),
        threads_per_file=max(1, args.threads),
    )
    elapsed = time.perf_counter() - t0

    skipped = [r.file for r in rows if r.skipped]
    hard = [r for r in rows if (not r.skipped) and r.error]
    report = SweepReport(
        ok=not hard,
        jobs=max(1, args.jobs),
        seconds=elapsed,
        privateCorpus=str(private),
        handleLock=str(handle_lock),
        skipped=skipped,
        rows=rows,
        error="; ".join(f"{r.file}: {r.error}" for r in hard),
    )

    md = format_markdown(report)
    sys.stdout.write(md)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(md, encoding="utf-8")
    payload = json.dumps(report.to_json(), indent=2, sort_keys=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload + "\n", encoding="utf-8")

    if args.compare:
        print("\n## compare vs spike-diag (main @ fd4954d, AC3 authority)", file=sys.stderr)
        spike = compare_report(report, SPIKE_HAND)
        if spike:
            for line in spike:
                print(f"  DRIFT {line}", file=sys.stderr)
        else:
            print("  handle-lock + Body12 (+Body9) match within tolerance", file=sys.stderr)
        print("\n## compare vs STRESS-SWEEP-2.md (stacked branch; expect J6 drift)", file=sys.stderr)
        sweep2 = compare_report(report, SWEEP2_HAND)
        if sweep2:
            for line in sweep2:
                print(f"  DRIFT {line}", file=sys.stderr)
        else:
            print("  no drift vs STRESS-SWEEP-2", file=sys.stderr)

    return 0 if report.ok else 1


if __name__ == "__main__":
    sys.exit(main())
