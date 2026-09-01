#!/usr/bin/env python3
"""collapse_fail_ratchet — DIAG_COLLAPSE fail= is a down-only ratchet.

Green at dc49e45 once the JSON datum is seeded from a full build of THIS
tree (G-3: never copy a number from another tree). The brief's 56 is a
hint only.

Asserts:
  * DIAG_COLLAPSE present with required keys (missing line/key = not measured)
  * fail <= ratchet.failCeiling
  * a lower live fail prints the tightening instruction (does not rewrite JSON)
  * chainEdgeFail class set is a subset of the recorded baseline
  * handle-lock fail stays at the recorded baseline

Usage:
  collapse_fail_ratchet.py --self-test
  collapse_fail_ratchet.py --synthetic-pass
  collapse_fail_ratchet.py --measure --binary ./build/stl2step
  collapse_fail_ratchet.py --binary ./build/stl2step
"""

from __future__ import annotations

import argparse
import json
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
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "collapse-ratchet.json"
RATCHET_NAME = "collapse-ratchet.json"

COLLAPSE_REQUIRED = ("mix", "none", "fail", "ok", "total", "recover", "rounds")

# Warning text that enrolls chainEdgeFail (src/refit_build.cpp rebuildCollapsed).
# cyl|cyl IntAna none is IA_CYLCYL_NOGEOM and is NEVER a chainEdgeFail class.
CLASS_PATTERNS: Tuple[Tuple[str, str], ...] = (
    ("intAna-plane-plane", "IntAna plane|plane"),
    ("intAna-plane-cyl", "IntAna plane|cyl"),
    ("makeEdgeFail", "analytic MakeEdge failed"),
)


class GateError(Exception):
    """Hard gate failure."""


@dataclass
class CollapseLine:
    mix: int
    none: int
    fail: int
    ok: int
    total: int
    recover: int
    rounds: int
    extra: Dict[str, str] = field(default_factory=dict)


@dataclass
class Check:
    name: str
    ok: bool
    message: str
    group: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


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
    kind = parts[0]
    kv: Dict[str, str] = {}
    for tok in parts[1:]:
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        kv[k] = v
    return kind, kv


def parse_collapse_lines(stderr: str) -> Tuple[List[CollapseLine], Optional[str]]:
    """Return parsed DIAG_COLLAPSE rows, or a not-measured reason."""
    rows: List[CollapseLine] = []
    for line in (stderr or "").splitlines():
        parsed = parse_diag_kv(line)
        if parsed is None:
            continue
        kind, kv = parsed
        if kind != "DIAG_COLLAPSE":
            continue
        missing = [k for k in COLLAPSE_REQUIRED if k not in kv]
        if missing:
            return [], f"DIAG_COLLAPSE missing key(s) {missing} (not measured)"
        extra = {k: v for k, v in kv.items() if k not in COLLAPSE_REQUIRED}
        rows.append(
            CollapseLine(
                mix=int(float(kv["mix"])),
                none=int(float(kv["none"])),
                fail=int(float(kv["fail"])),
                ok=int(float(kv["ok"])),
                total=int(float(kv["total"])),
                recover=int(float(kv["recover"])),
                rounds=int(float(kv["rounds"])),
                extra=extra,
            )
        )
    if not rows:
        return [], "zero DIAG_COLLAPSE lines (not measured)"
    return rows, None


def first_pass_fail(rows: Sequence[CollapseLine]) -> int:
    """fail= of the first recover=0 rounds=0 line (chainEdgeFail enrollment).

    rebuildCollapsed may re-run at recover=0; summing would triple-count. The
    first DIAG_COLLAPSE line is the birth pass.
    """
    for r in rows:
        if r.recover == 0 and r.rounds == 0:
            return r.fail
    return rows[0].fail


def chain_edge_fail_classes(stderr: str, rows: Sequence[CollapseLine]) -> Set[str]:
    """Classes actually seen. fail>0 without a warning is still makeEdgeFail."""
    found: Set[str] = set()
    text = stderr or ""
    for name, needle in CLASS_PATTERNS:
        if needle in text:
            found.add(name)
    if first_pass_fail(rows) > 0 and "makeEdgeFail" not in found:
        # DIAG_COLLAPSE fail= counts MakeEdge failures; the warning text is the
        # class name when present, but fail>0 is itself the class.
        found.add("makeEdgeFail")
    return found


def load_ratchet(path: Path) -> Dict[str, Any]:
    if path.name != RATCHET_NAME:
        raise GateError(f"floor authority must be {RATCHET_NAME}, opened {path}")
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("id") != "collapse-ratchet":
        raise GateError(f"{path} is not the collapse-ratchet (id={doc.get('id')!r})")
    if str(doc.get("authority") or "") != RATCHET_NAME:
        raise GateError(f"{path}: authority must be {RATCHET_NAME!r}")
    if "failCeiling" not in doc:
        raise GateError(f"{path}: missing failCeiling")
    if "handleLockFail" not in doc:
        raise GateError(f"{path}: missing handleLockFail")
    if "chainEdgeFailClasses" not in doc:
        raise GateError(f"{path}: missing chainEdgeFailClasses")
    return doc


def evaluate_pickup(
    rows: Sequence[CollapseLine],
    classes: Set[str],
    ratchet: Mapping[str, Any],
    *,
    measure_err: Optional[str] = None,
) -> List[Check]:
    checks: List[Check] = []

    def add(name: str, ok: bool, message: str, **details: Any) -> None:
        checks.append(
            Check(name=name, ok=ok, message=message, group="pickup", details=details)
        )

    if measure_err:
        add("measured", False, f"handle-pickup DIAG_COLLAPSE: {measure_err}")
        return checks
    ceiling = int(ratchet["failCeiling"])
    fail = first_pass_fail(rows)
    add(
        "failCeiling",
        fail <= ceiling,
        (
            f"DIAG_COLLAPSE fail={fail} > ratchet failCeiling={ceiling}"
            if fail > ceiling
            else f"DIAG_COLLAPSE fail={fail} <= ratchet failCeiling={ceiling}"
        ),
        fail=fail,
        failCeiling=ceiling,
    )
    if fail < ceiling:
        add(
            "ratchetTighten",
            True,
            f"RATCHET-TIGHTEN: measured fail={fail} < committed {ceiling}; "
            f"set tests/gates/baseline/{RATCHET_NAME} failCeiling to {fail} "
            f"(do not auto-rewrite inside ctest)",
            fail=fail,
            failCeiling=ceiling,
        )
    else:
        add(
            "ratchetTighten",
            True,
            f"fail={fail} equals ceiling {ceiling} (no tighten)",
        )
    baseline = {str(x) for x in (ratchet.get("chainEdgeFailClasses") or [])}
    extra = sorted(classes - baseline)
    add(
        "chainEdgeFailClasses",
        not extra,
        (
            f"new chainEdgeFail class(es) {extra} not in baseline {sorted(baseline)}"
            if extra
            else f"chainEdgeFail classes {sorted(classes)}subseteq {sorted(baseline)}"
        ),
        seen=sorted(classes),
        baseline=sorted(baseline),
        extra=extra,
    )
    return checks


def evaluate_handle_lock(
    rows: Sequence[CollapseLine],
    ratchet: Mapping[str, Any],
    *,
    measure_err: Optional[str] = None,
) -> List[Check]:
    """handle-lock at this tree takes Stage-P and emits zero DIAG_COLLAPSE.

    Recorded baseline: fail class 0. Absence of the line is the measurement
    (prism path, collapse not entered) — not a silent pass. A later
    DIAG_COLLAPSE fail>0 is a regression.
    """
    checks: List[Check] = []
    want = int(ratchet["handleLockFail"])
    expect_lines = int(ratchet.get("handleLockCollapseLines", 0))
    if measure_err:
        if expect_lines == 0 and want == 0:
            checks.append(
                Check(
                    "handle-lock.fail",
                    True,
                    "handle-lock DIAG_COLLAPSE absent (prism path); "
                    "fail class 0 (recorded baseline)",
                    "handle-lock",
                    details={"got": 0, "want": want, "lines": 0},
                )
            )
            return checks
        checks.append(
            Check(
                "handle-lock.fail",
                False,
                f"handle-lock DIAG_COLLAPSE: {measure_err}",
                "handle-lock",
            )
        )
        return checks
    got = first_pass_fail(rows)
    checks.append(
        Check(
            "handle-lock.fail",
            got == want,
            (
                f"handle-lock DIAG_COLLAPSE fail={got} != baseline {want}"
                if got != want
                else f"handle-lock DIAG_COLLAPSE fail={got} (baseline {want})"
            ),
            "handle-lock",
            details={"got": got, "want": want, "lines": len(rows)},
        )
    )
    return checks


def format_checks(checks: Sequence[Check]) -> str:
    lines = []
    for c in checks:
        tag = "PASS" if c.ok else "FAIL"
        lines.append(f"  [{c.group}/{tag}] {c.message}")
    return "\n".join(lines)


def failing_names(checks: Sequence[Check]) -> List[str]:
    return [c.name for c in checks if not c.ok]


def run_trueform(
    binary: Path,
    stl: Path,
    step: Path,
    *,
    threads: int,
) -> Tuple[Dict[str, Any], str, int]:
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
        str(threads),
    ]
    environ = os.environ.copy()
    environ["STL2STEP_COLLAPSE_DIAG"] = "1"
    proc = subprocess.run(cmd, capture_output=True, text=True, env=environ)
    result = parse_result_line(proc.stdout)
    return result, proc.stderr, proc.returncode


def convert_pair(
    binary: Path,
    pickup_stl: Path,
    hl_stl: Path,
    jobs: int,
) -> Dict[str, Tuple[Dict[str, Any], str, int]]:
    workers = max(1, min(int(jobs) or 1, 2))
    arts: Dict[str, Tuple[Dict[str, Any], str, int]] = {}
    with tempfile.TemporaryDirectory(prefix="collapse_fail_ratchet_") as td:
        work = Path(td)
        tasks = {"pickup": pickup_stl, "handle-lock": hl_stl}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    run_trueform, binary, stl, work / f"{name}.step", threads=0
                ): name
                for name, stl in tasks.items()
            }
            for fut in as_completed(futs):
                arts[futs[fut]] = fut.result()
    return arts


def measure_payload(
    arts: Mapping[str, Tuple[Dict[str, Any], str, int]],
) -> Dict[str, Any]:
    pk_res, pk_err, _ = arts["pickup"]
    hl_res, hl_err, _ = arts["handle-lock"]
    pk_rows, pk_err_m = parse_collapse_lines(pk_err)
    hl_rows, hl_err_m = parse_collapse_lines(hl_err)
    pk_fail = first_pass_fail(pk_rows) if pk_rows else None
    hl_fail = first_pass_fail(hl_rows) if hl_rows else None
    classes = chain_edge_fail_classes(pk_err, pk_rows) if pk_rows else set()
    return {
        "pickupOk": bool(pk_res.get("ok")),
        "handleLockOk": bool(hl_res.get("ok")),
        "pickupFail": pk_fail,
        "handleLockFail": hl_fail,
        "pickupMeasureError": pk_err_m,
        "handleLockMeasureError": hl_err_m,
        "chainEdgeFailClasses": sorted(classes),
        "pickupRows": [r.__dict__ for r in pk_rows],
        "handleLockRows": [r.__dict__ for r in hl_rows],
        "pickupStderrTail": "\n".join(
            ln for ln in (pk_err or "").splitlines() if ln.startswith("DIAG_COLLAPSE")
        ),
        "handleLockStderrTail": "\n".join(
            ln for ln in (hl_err or "").splitlines() if ln.startswith("DIAG_COLLAPSE")
        ),
    }


def run_measure(binary: Path, pickup_stl: Path, hl_stl: Path, jobs: int) -> int:
    arts = convert_pair(binary, pickup_stl, hl_stl, jobs)
    payload = measure_payload(arts)
    print("collapse_fail_ratchet --measure")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if payload["pickupMeasureError"]:
        print("collapse_fail_ratchet FAIL: pickup not measured", file=sys.stderr)
        return 1
    # handle-lock may emit zero DIAG_COLLAPSE lines (Stage-P). That is the
    # measurement: fail class 0. Record it; do not treat as a silent pass.
    print(
        f"SEED failCeiling={payload['pickupFail']} "
        f"handleLockFail={payload['handleLockFail']} "
        f"classes={payload['chainEdgeFailClasses']}"
    )
    return 0


def run_live(
    binary: Path,
    *,
    pickup_stl: Path = PICKUP_STL,
    hl_stl: Path = HL_STL,
    ratchet_path: Path = DEFAULT_RATCHET,
    jobs: int = 2,
) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    if not pickup_stl.is_file():
        raise GateError(f"handle-pickup STL missing: {pickup_stl}")
    if not hl_stl.is_file():
        raise GateError(f"handle-lock STL missing: {hl_stl}")
    ratchet = load_ratchet(ratchet_path)
    arts = convert_pair(binary, pickup_stl, hl_stl, jobs)
    pk_res, pk_err, pk_rc = arts["pickup"]
    hl_res, hl_err, hl_rc = arts["handle-lock"]
    checks: List[Check] = []
    if not pk_res.get("ok"):
        checks.append(
            Check(
                "conversion",
                False,
                f"pickup conversion failed: {pk_res.get('error')}",
                "pickup",
            )
        )
    pk_rows, pk_m = parse_collapse_lines(pk_err)
    pk_classes = chain_edge_fail_classes(pk_err, pk_rows) if pk_rows else set()
    checks.extend(evaluate_pickup(pk_rows, pk_classes, ratchet, measure_err=pk_m))
    if not hl_res.get("ok"):
        checks.append(
            Check(
                "conversion",
                False,
                f"handle-lock conversion failed: {hl_res.get('error')}",
                "handle-lock",
            )
        )
    hl_rows, hl_m = parse_collapse_lines(hl_err)
    checks.extend(evaluate_handle_lock(hl_rows, ratchet, measure_err=hl_m))

    print(
        f"collapse_fail_ratchet  ratchet={ratchet_path}  "
        f"ceiling={ratchet.get('failCeiling')}  pickup_rc={pk_rc} lock_rc={hl_rc}",
        flush=True,
    )
    print(format_checks(checks), flush=True)
    named = failing_names(checks)
    if named:
        print(
            f"collapse_fail_ratchet FAIL  failing assertion(s): {', '.join(named)}",
            file=sys.stderr,
            flush=True,
        )
        return 1
    print("collapse_fail_ratchet PASS", flush=True)
    return 0


def synthetic_pass_stderr(fail: int = 4, classes: bool = True) -> str:
    warn = "warning: smooth: analytic MakeEdge failed — keeping mesh polyline\n" if classes else ""
    return (
        warn
        + f"DIAG_COLLAPSE mix=1 none=2 fail={fail} ok=10 total=13 recover=0 rounds=0 extra=1\n"
        + "DIAG_COLLAPSE mix=0 none=0 fail=0 ok=13 total=13 recover=1 rounds=1\n"
    )


def synthetic_lock_stderr(fail: int = 0) -> str:
    return f"DIAG_COLLAPSE mix=0 none=0 fail={fail} ok=68 total=68 recover=0 rounds=0\n"


def synthetic_ratchet(**overrides: Any) -> Dict[str, Any]:
    doc: Dict[str, Any] = {
        "id": "collapse-ratchet",
        "authority": RATCHET_NAME,
        "fixtureId": "handle-pickup",
        "failCeiling": 4,
        "handleLockFail": 0,
        "handleLockCollapseLines": 0,
        "chainEdgeFailClasses": ["makeEdgeFail", "intAna-plane-plane", "intAna-plane-cyl"],
        "notes": "synthetic",
    }
    doc.update(overrides)
    return doc


def run_synthetic_pass() -> int:
    ratchet = synthetic_ratchet()
    pk_rows, pk_m = parse_collapse_lines(synthetic_pass_stderr())
    hl_rows, hl_m = parse_collapse_lines(synthetic_lock_stderr())
    classes = chain_edge_fail_classes(synthetic_pass_stderr(), pk_rows)
    checks = evaluate_pickup(pk_rows, classes, ratchet, measure_err=pk_m)
    checks.extend(evaluate_handle_lock(hl_rows, ratchet, measure_err=hl_m))
    print("collapse_fail_ratchet --synthetic-pass")
    print(format_checks(checks))
    named = failing_names(checks)
    if named:
        print(
            f"collapse_fail_ratchet FAIL  synthetic should pass; failed: {named}",
            file=sys.stderr,
        )
        return 1
    print("collapse_fail_ratchet PASS  synthetic floors met")
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

    rows, err = parse_collapse_lines(synthetic_pass_stderr(fail=4))
    check(err is None, "synthetic DIAG_COLLAPSE parses")
    check(len(rows) == 2, "multiple DIAG_COLLAPSE lines kept")
    check(first_pass_fail(rows) == 4, "first-pass fail uses recover=0 rounds=0")
    check(rows[0].extra.get("extra") == "1", "unknown extra keys tolerated")

    empty_rows, empty_err = parse_collapse_lines("no diag here\n")
    check(not empty_rows and empty_err is not None and "not measured" in empty_err,
          "zero lines is not measured")
    miss_rows, miss_err = parse_collapse_lines("DIAG_COLLAPSE mix=1 none=2 fail=3\n")
    check(miss_err is not None and "missing key" in miss_err, "missing key is not measured")

    ratchet = synthetic_ratchet()
    classes = chain_edge_fail_classes(synthetic_pass_stderr(), rows)
    check("makeEdgeFail" in classes, "MakeEdge warning classifies makeEdgeFail")
    pk = evaluate_pickup(rows, classes, ratchet)
    check(failing_names(pk) == [], "fail==ceiling PASSES")

    high_rows, _ = parse_collapse_lines(synthetic_pass_stderr(fail=5))
    high = evaluate_pickup(high_rows, classes, ratchet)
    check("failCeiling" in failing_names(high), "fail=5 > ceiling=4 FAILS and names failCeiling")
    check(any("fail=5" in c.message and "failCeiling=4" in c.message for c in high),
          "red path names the failing assertion with numbers")

    low_rows, _ = parse_collapse_lines(synthetic_pass_stderr(fail=2))
    low = evaluate_pickup(low_rows, classes, ratchet)
    check(failing_names(low) == [], "fail < ceiling still PASSES")
    check(any("RATCHET-TIGHTEN" in c.message and "fail=2" in c.message for c in low),
          "lower fail prints ratchet-tightening instruction")

    new_cls = evaluate_pickup(rows, classes | {"brand-new-class"}, ratchet)
    check("chainEdgeFailClasses" in failing_names(new_cls), "new chainEdgeFail class FAILS")
    check(any("brand-new-class" in c.message for c in new_cls), "new class is named")

    hl_ok = evaluate_handle_lock(
        parse_collapse_lines(synthetic_lock_stderr(0))[0], ratchet
    )
    check(failing_names(hl_ok) == [], "handle-lock fail=0 PASSES")
    hl_bad = evaluate_handle_lock(
        parse_collapse_lines(synthetic_lock_stderr(1))[0], ratchet
    )
    check("handle-lock.fail" in failing_names(hl_bad), "handle-lock fail!=baseline FAILS")
    hl_abs = evaluate_handle_lock([], ratchet, measure_err="zero DIAG_COLLAPSE lines (not measured)")
    check(failing_names(hl_abs) == [], "handle-lock absent DIAG_COLLAPSE is fail class 0")
    check(any("prism path" in c.message for c in hl_abs), "absent line names prism path")

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "not-the-ratchet.json"
        p.write_text("{}", encoding="utf-8")
        try:
            load_ratchet(p)
            check(False, "mis-named ratchet file must FAIL")
        except GateError:
            check(True, "mis-named ratchet file FAIL")

    check(
        parse_result_line('RESULT {"ok":true}\n')["ok"] is True,
        "parse_result_line",
    )
    rc = run_synthetic_pass()
    check(rc == 0, "--synthetic-pass returns 0")
    check(PICKUP_STL.is_file(), "handle-pickup.stl present")
    check(HL_STL.is_file(), "handle-lock.stl present")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path)
    p.add_argument("--pickup-stl", type=Path, default=PICKUP_STL)
    p.add_argument("--handle-lock-stl", type=Path, default=HL_STL)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--synthetic-pass", action="store_true")
    p.add_argument(
        "--measure",
        action="store_true",
        help="print live DIAG_COLLAPSE fail= / classes (does not write JSON)",
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
            print("collapse_fail_ratchet: --binary is required", file=sys.stderr)
            return 1
        binary = Path(args.binary)
        if args.measure:
            return run_measure(binary, Path(args.pickup_stl), Path(args.handle_lock_stl), args.jobs)
        return run_live(
            binary,
            pickup_stl=Path(args.pickup_stl),
            hl_stl=Path(args.handle_lock_stl),
            ratchet_path=Path(args.ratchet),
            jobs=args.jobs,
        )
    except GateError as exc:
        print(f"collapse_fail_ratchet FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
