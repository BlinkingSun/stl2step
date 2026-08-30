#!/usr/bin/env python3
"""J6 freeEdges watch — private-corpus ctest gate (lane gates-j6watch).

Red-line ceilings from `_team/reports/spike-diag.md` (main @ fd4954d):
Body9=34, Body12=81, Body18=18, Body20=32 (never worse). Lane-I stretch
targets 14/26/14/8 are PARKED unless `--unpark j6-lane-i` (same register
convention as run_gates.py).

Private meshes resolve via env `STL2STEP_PRIVATE_CORPUS` (or
`--private-corpus`). Unset/empty env, missing corpus, or missing STLs →
exit 77 (ctest SKIP_RETURN_CODE), never FAIL. Paths are never required
inputs in committed fixtures.

Stdlib only. Files convert in parallel (ThreadPoolExecutor); each CLI
gets a slice of cores.
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
from typing import Any, Dict, List, Optional, Sequence, Set

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASELINE = Path(__file__).resolve().parent / "j6-watch.json"
ENV_PRIVATE_CORPUS = "STL2STEP_PRIVATE_CORPUS"
SKIP_EXIT = 77
FREE_EDGES_RE = re.compile(r"freeEdges=(\d+)")
J6_MARK = "J6:"

# Spike-diag main @ fd4954d — duplicated here so --self-test can pin file truth
# without silently drifting if the JSON is rewritten incorrectly.
SPIKE_REDLINES = {"Body9": 34, "Body12": 81, "Body18": 18, "Body20": 32}
LANE_I_TARGETS = {"Body9": 14, "Body12": 26, "Body18": 14, "Body20": 8}


@dataclass
class BodySpec:
    id: str
    stl_name: str
    redline: int
    parked_lane_i: int


@dataclass
class Baseline:
    path: Path
    parked_gate_id: str
    park_reason: str
    bodies: List[BodySpec]

    def body_ids(self) -> List[str]:
        return [b.id for b in self.bodies]


@dataclass
class Outcome:
    gate_id: str
    fixture_id: str
    status: str  # PASS | FAIL | PARKED | SKIP
    message: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


def resolve_corpus(
    cli_path: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
) -> Path:
    """Resolve the private STL directory.

    Precedence: --private-corpus, then STL2STEP_PRIVATE_CORPUS.
    Unset or empty env is treated as absent (caller SKIPs). No home-path default.
    """
    if cli_path is not None:
        return Path(cli_path).expanduser()
    environ = os.environ if env is None else env
    raw = (environ.get(ENV_PRIVATE_CORPUS) or "").strip()
    if not raw:
        return Path()
    return Path(raw).expanduser()


def corpus_is_absent(corpus: Path) -> bool:
    if not corpus or str(corpus) == ".":
        return True
    return not corpus.is_dir()


def load_baseline(path: Path = DEFAULT_BASELINE) -> Baseline:
    """Load red-line JSON. Asserts this file, not a corpus expected.json."""
    path = Path(path)
    if path.name.endswith("expected.json"):
        raise RuntimeError("j6-watch must not read a corpus expected.json")
    doc = json.loads(path.read_text())
    if doc.get("id") != "j6-watch":
        raise RuntimeError(f"{path} is not the j6-watch baseline (id={doc.get('id')!r})")
    bodies: List[BodySpec] = []
    for bid, row in (doc.get("bodies") or {}).items():
        bodies.append(
            BodySpec(
                id=str(bid),
                stl_name=str(row["stl"]),
                redline=int(row["redlineFreeEdges"]),
                parked_lane_i=int(row["parkedLaneI"]),
            )
        )
    if not bodies:
        raise RuntimeError(f"{path} has no bodies")
    return Baseline(
        path=path,
        parked_gate_id=str(doc.get("parkedGateId") or "j6-lane-i"),
        park_reason=str(doc.get("parkReason") or "Lane I stretch parked"),
        bodies=bodies,
    )


def stl_for(corpus: Path, spec: BodySpec) -> Optional[Path]:
    p = corpus / spec.stl_name
    return p if p.is_file() else None


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise RuntimeError("no RESULT line in stl2step stdout")


def parse_j6_free_edges(
    warnings: Optional[Sequence[str]] = None,
    stderr: str = "",
) -> int:
    """Last J6 freeEdges count (0 if the analytic shell closed / no J6 line).

    Body12 on main emits two J6 lines (103 then 81); the red-line is the
    final recover-pass value from spike-diag.
    """
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


def split_csv(value: Optional[str]) -> Set[str]:
    if not value:
        return set()
    return {part.strip() for part in value.split(",") if part.strip()}


def evaluate_body(
    spec: BodySpec,
    free_edges: int,
    *,
    parked_gate_id: str,
    park_reason: str,
    unparked: Optional[Set[str]] = None,
) -> List[Outcome]:
    """Red-line is LIVE. Lane-I is PARKED unless --unpark names the gate id."""
    unparked = unparked or set()
    out: List[Outcome] = []
    if free_edges > spec.redline:
        out.append(
            Outcome(
                "j6-redline",
                spec.id,
                "FAIL",
                f"{spec.id} freeEdges={free_edges} > red-line {spec.redline}",
                details={"freeEdges": free_edges, "redline": spec.redline},
            )
        )
    else:
        out.append(
            Outcome(
                "j6-redline",
                spec.id,
                "PASS",
                f"{spec.id} freeEdges={free_edges} <= red-line {spec.redline}",
                details={"freeEdges": free_edges, "redline": spec.redline},
            )
        )

    lane_fail = free_edges > spec.parked_lane_i
    lane_msg = (
        f"{spec.id} freeEdges={free_edges} "
        f"{'>' if lane_fail else '<='} Lane-I {spec.parked_lane_i}"
    )
    if lane_fail and parked_gate_id not in unparked:
        out.append(
            Outcome(
                parked_gate_id,
                spec.id,
                "PARKED",
                f"PARKED: {park_reason} — {lane_msg}",
                details={
                    "freeEdges": free_edges,
                    "laneI": spec.parked_lane_i,
                    "parkReason": park_reason,
                    "wouldBe": "FAIL",
                },
            )
        )
    elif lane_fail:
        out.append(
            Outcome(
                parked_gate_id,
                spec.id,
                "FAIL",
                lane_msg,
                details={"freeEdges": free_edges, "laneI": spec.parked_lane_i},
            )
        )
    else:
        out.append(
            Outcome(
                parked_gate_id,
                spec.id,
                "PASS",
                lane_msg,
                details={"freeEdges": free_edges, "laneI": spec.parked_lane_i},
            )
        )
    return out


def run_trueform(
    binary: Path,
    stl: Path,
    step_out: Path,
    *,
    threads: int,
) -> Dict[str, Any]:
    cmd = [
        str(binary),
        str(stl),
        "-o",
        str(step_out),
        "--quiet",
        "--smooth",
        "--no-verify",
        "--threads",
        str(threads),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    result = parse_result_line(proc.stdout)
    warnings = result.get("warnings") or []
    free_edges = parse_j6_free_edges(warnings, proc.stderr)
    return {
        "exit_code": proc.returncode,
        "result": result,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "free_edges": free_edges,
        "step": str(step_out),
    }


def convert_one(
    binary: Path,
    spec: BodySpec,
    stl: Path,
    work_dir: Path,
    threads: int,
) -> Dict[str, Any]:
    step_out = work_dir / f"{spec.id}.step"
    art = run_trueform(binary, stl, step_out, threads=threads)
    art["body"] = spec.id
    art["stl"] = str(stl)
    return art


def present_bodies(
    corpus: Path, baseline: Baseline
) -> List[tuple[BodySpec, Path]]:
    found: List[tuple[BodySpec, Path]] = []
    for spec in baseline.bodies:
        stl = stl_for(corpus, spec)
        if stl is not None:
            found.append((spec, stl))
    return found


def format_table(outcomes: Sequence[Outcome]) -> str:
    lines = [
        f"{'gate':<12} {'body':<10} {'status':<8} message",
        "-" * 72,
    ]
    for o in outcomes:
        lines.append(f"{o.gate_id:<12} {o.fixture_id:<10} {o.status:<8} {o.message}")
    return "\n".join(lines)


def self_test() -> int:
    """Exercise the gate API without a private corpus or engine binary."""
    errors: List[str] = []

    def check(cond: bool, msg: str) -> None:
        if not cond:
            errors.append(msg)

    baseline = load_baseline(DEFAULT_BASELINE)
    check(baseline.path == DEFAULT_BASELINE, "baseline path is j6-watch.json")
    check("expected.json" not in str(baseline.path), "must not read fixture expected.json")
    check(baseline.parked_gate_id == "j6-lane-i", f"parked id {baseline.parked_gate_id}")
    got_red = {b.id: b.redline for b in baseline.bodies}
    got_li = {b.id: b.parked_lane_i for b in baseline.bodies}
    check(got_red == SPIKE_REDLINES, f"red-lines {got_red} != spike {SPIKE_REDLINES}")
    check(got_li == LANE_I_TARGETS, f"Lane-I {got_li} != {LANE_I_TARGETS}")

    sample_warns = [
        "J6: shell not closed freeEdges=103 faces=2995 recover=0",
        "J6: shell not closed freeEdges=81 faces=3799 recover=1",
    ]
    check(parse_j6_free_edges(sample_warns) == 81, "final freeEdges from two J6 lines")
    check(parse_j6_free_edges([]) == 0, "no J6 => 0")
    check(
        parse_j6_free_edges(
            None, "note: other\nJ6: shell not closed freeEdges=34 faces=1 recover=0\n"
        )
        == 34,
        "stderr fallback",
    )
    check(
        parse_result_line('progress\nRESULT {"ok":true,"warnings":[]}\n')["ok"] is True,
        "parse_result_line",
    )

    b9 = next(b for b in baseline.bodies if b.id == "Body9")
    red_pass = evaluate_body(
        b9, 34, parked_gate_id=baseline.parked_gate_id, park_reason=baseline.park_reason
    )
    check(red_pass[0].status == "PASS" and red_pass[0].gate_id == "j6-redline", "red-line ==")
    check(red_pass[1].status == "PARKED", "Lane-I 34>14 parked")
    red_fail = evaluate_body(
        b9, 35, parked_gate_id=baseline.parked_gate_id, park_reason=baseline.park_reason
    )
    check(red_fail[0].status == "FAIL", "red-line 35>34 FAIL")
    unparked = evaluate_body(
        b9,
        34,
        parked_gate_id=baseline.parked_gate_id,
        park_reason=baseline.park_reason,
        unparked={baseline.parked_gate_id},
    )
    check(unparked[1].status == "FAIL", "--unpark surfaces Lane-I FAIL")
    lane_ok = evaluate_body(
        b9, 10, parked_gate_id=baseline.parked_gate_id, park_reason=baseline.park_reason
    )
    check(lane_ok[1].status == "PASS", "Lane-I met is PASS even while parked")

    missing = resolve_corpus(env={ENV_PRIVATE_CORPUS: ""})
    check(corpus_is_absent(missing), "empty env => absent corpus")
    missing_dir = resolve_corpus(env={ENV_PRIVATE_CORPUS: "/no/such/stl2step-private-corpus"})
    check(corpus_is_absent(missing_dir), "missing dir => absent")
    check(SKIP_EXIT == 77, "SKIP exit is 77")

    defaulted = resolve_corpus(env={})
    check(corpus_is_absent(defaulted), "unset env => absent corpus")

    if errors:
        print("SELF-TEST FAIL:")
        for e in errors:
            print(f"  - {e}")
        return 1
    print("SELF-TEST PASS")
    print(f"  baseline={baseline.path}")
    print(f"  red-lines={got_red}")
    print(f"  lane-i={got_li} parked as {baseline.parked_gate_id}")
    print(f"  skip_exit={SKIP_EXIT}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path, help="Path to stl2step CLI (required unless --self-test)")
    p.add_argument(
        "--private-corpus",
        type=Path,
        help=f"Private STL directory (else ${ENV_PRIVATE_CORPUS}; unset SKIPs)",
    )
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE, help="j6-watch JSON")
    p.add_argument(
        "--unpark",
        help="Comma-separated parked gate ids to surface honest FAILs (e.g. j6-lane-i)",
    )
    p.add_argument("--jobs", type=int, default=4, help="Parallel file conversions")
    p.add_argument("--self-test", action="store_true", help="Exercise the gate API and exit")
    p.add_argument(
        "--report",
        type=Path,
        help="Optional JSON report path",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.self_test:
        return self_test()

    if not args.binary:
        print("j6_freeedges_watch: --binary is required", file=sys.stderr)
        return 1
    binary = Path(args.binary)
    if not binary.is_file():
        print(f"j6_freeedges_watch: binary not found: {binary}", file=sys.stderr)
        return 1

    baseline = load_baseline(args.baseline)
    corpus = resolve_corpus(args.private_corpus)
    if corpus_is_absent(corpus):
        shown = str(corpus) if corpus and str(corpus) != "." else "<empty>"
        print(f"SKIP: private corpus absent ({ENV_PRIVATE_CORPUS}={shown})")
        return SKIP_EXIT

    found = present_bodies(corpus, baseline)
    missing = [b.id for b in baseline.bodies if all(b.id != s.id for s, _ in found)]
    if not found:
        print(
            f"SKIP: no watch STLs in {corpus} "
            f"(looked for {', '.join(b.stl_name for b in baseline.bodies)})"
        )
        return SKIP_EXIT
    for mid in missing:
        print(f"SKIP {mid}: mesh not present under {corpus}")

    unparked = split_csv(args.unpark)
    jobs = max(1, int(args.jobs) or 1)
    workers = min(jobs, len(found))
    cpu = os.cpu_count() or 4
    per_file_threads = max(1, cpu // workers)

    outcomes: List[Outcome] = []
    arts: Dict[str, Dict[str, Any]] = {}
    with tempfile.TemporaryDirectory(prefix="j6watch-") as tmp:
        work = Path(tmp)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    convert_one, binary, spec, stl, work, per_file_threads
                ): spec
                for spec, stl in found
            }
            for fut in as_completed(futs):
                spec = futs[fut]
                try:
                    art = fut.result()
                except Exception as exc:
                    outcomes.append(
                        Outcome(
                            "j6-redline",
                            spec.id,
                            "FAIL",
                            f"{spec.id} convert error: {exc}",
                        )
                    )
                    continue
                arts[spec.id] = art
                outcomes.extend(
                    evaluate_body(
                        spec,
                        int(art["free_edges"]),
                        parked_gate_id=baseline.parked_gate_id,
                        park_reason=baseline.park_reason,
                        unparked=unparked,
                    )
                )

    # Deterministic print order (baseline body order).
    order = {b.id: i for i, b in enumerate(baseline.bodies)}
    outcomes.sort(key=lambda o: (order.get(o.fixture_id, 99), o.gate_id))
    print(format_table(outcomes))
    n_fail = sum(1 for o in outcomes if o.status == "FAIL")
    n_park = sum(1 for o in outcomes if o.status == "PARKED")
    n_pass = sum(1 for o in outcomes if o.status == "PASS")
    print(
        f"SUMMARY PASS={n_pass} PARKED={n_park} FAIL={n_fail} "
        f"bodies={len(found)}/{len(baseline.bodies)} "
        f"jobs={workers} threads/file={per_file_threads}"
    )
    if args.report:
        payload = {
            "ok": n_fail == 0,
            "baseline": str(baseline.path),
            "corpus": str(corpus),
            "unparked": sorted(unparked),
            "outcomes": [
                {
                    "gate": o.gate_id,
                    "body": o.fixture_id,
                    "status": o.status,
                    "message": o.message,
                    "details": o.details,
                }
                for o in outcomes
            ],
            "runs": {
                bid: {
                    "freeEdges": art["free_edges"],
                    "exit": art["exit_code"],
                    "warnings": (art["result"].get("warnings") or []),
                }
                for bid, art in arts.items()
            },
        }
        args.report.write_text(json.dumps(payload, indent=2) + "\n")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
