#!/usr/bin/env python3
"""edge_class_gate — D-130-2's red lines on the shipped edge census.

The engine classifies every edge two analytic faces share (RESULT
``edgeClasses``, one ``DIAG_EDGECLASS`` line per edge under
``STL2STEP_P2_DIAG=1``).  This gate holds the three red lines DECISION-130
D-130-2 states, over the whole corpus plus the P2 build fixtures:

  * ``unhandled`` on an analytic|analytic edge  == 0
  * ``overTol``  (deviation above the recorded tolerance)   == 0
  * ``overCap``  (recorded tolerance above ``meshTolCap``)  == 0

and ratchets the tier-2 mesh-polyline count per fixture in
``tests/gates/baseline/edge-class-ratchet.json``: an increase FAILS (a
polyline that was analytic yesterday was silently converted), a decrease
PASSES and prints the re-baselining instruction.  The gate never rewrites
the ratchet itself.

A fixture whose red line is already violated on the branch that seeded the
ratchet is listed in ``expectedRed`` with the lines it violates.  That entry
is an admission, not a waiver: it is loud on every run, a fixture that heals
prints the instruction to delete its entry, and a NEW violation — a fixture
not listed, or a listed fixture violating a line its entry does not name —
is a hard FAIL.  Widening ``expectedRed`` to make a run green is a defect.

Usage:
  edge_class_gate.py --self-test
  edge_class_gate.py --synthetic-pass
  edge_class_gate.py --measure --binary ./build/stl2step      # seeding view
  edge_class_gate.py --binary ./build/stl2step
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
CORPUS = REPO / "tests" / "corpus"
BUILD_FIXTURES = REPO / "tests" / "gates" / "build_fixtures"
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "edge-class-ratchet.json"
RATCHET_NAME = "edge-class-ratchet.json"

RED_LINES = ("unhandled", "overTol", "overCap")


class GateError(Exception):
    """Hard gate failure."""


@dataclass
class Census:
    analytic: int
    polylineTier2: int
    unhandled: int
    overTol: int
    overCap: int

    @staticmethod
    def from_result(result: Mapping[str, Any]) -> "Census":
        ec = result.get("edgeClasses")
        if not isinstance(ec, dict):
            raise GateError("RESULT has no edgeClasses object (census not measured)")
        missing = [k for k in ("analytic", "polylineTier2") + RED_LINES if k not in ec]
        if missing:
            raise GateError(f"RESULT edgeClasses missing key(s) {missing} (not measured)")
        return Census(
            analytic=int(ec["analytic"]),
            polylineTier2=int(ec["polylineTier2"]),
            unhandled=int(ec["unhandled"]),
            overTol=int(ec["overTol"]),
            overCap=int(ec["overCap"]),
        )

    def red(self) -> Dict[str, int]:
        return {k: getattr(self, k) for k in RED_LINES if getattr(self, k) > 0}

    def asdict(self) -> Dict[str, int]:
        return {
            "analytic": self.analytic,
            "polylineTier2": self.polylineTier2,
            "unhandled": self.unhandled,
            "overTol": self.overTol,
            "overCap": self.overCap,
        }


@dataclass
class Check:
    name: str
    ok: bool
    message: str
    group: str
    details: Dict[str, Any] = field(default_factory=dict)


def fixtures() -> List[Tuple[str, Path]]:
    """Every STL the gate covers: the corpus, then the P2 build fixtures."""
    out: List[Tuple[str, Path]] = []
    for p in sorted(CORPUS.glob("*.stl")):
        out.append((p.stem, p))
    for p in sorted(BUILD_FIXTURES.glob("*/*.stl")):
        out.append((f"build_fixtures/{p.stem}", p))
    return out


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise GateError("no RESULT line in stl2step stdout")


def run_one(binary: Path, stl: Path, work: Path) -> Dict[str, Any]:
    step = work / (stl.stem + ".step")
    cmd = [str(binary), str(stl), "-o", str(step), "--smooth", "--no-verify", "--quiet"]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=os.environ.copy())
    try:
        result = parse_result_line(proc.stdout)
    except GateError as exc:
        return {"ok": False, "error": str(exc)}
    step.unlink(missing_ok=True)
    return result


def measure(binary: Path, jobs: int) -> Dict[str, Any]:
    """name -> RESULT dict (or an {'ok': False} stub)."""
    fx = fixtures()
    out: Dict[str, Any] = {}
    with tempfile.TemporaryDirectory(prefix="edge_class_gate_") as td:
        work = Path(td)
        with ThreadPoolExecutor(max_workers=max(1, int(jobs) or 1)) as pool:
            futs = {pool.submit(run_one, binary, stl, work): name for name, stl in fx}
            for fut, name in futs.items():
                out[name] = fut.result()
    return out


def load_ratchet(path: Path) -> Dict[str, Any]:
    if path.name != RATCHET_NAME:
        raise GateError(f"ratchet authority must be {RATCHET_NAME}, opened {path}")
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("id") != "edge-class-ratchet":
        raise GateError(f"{path} is not the edge-class ratchet (id={doc.get('id')!r})")
    if str(doc.get("authority") or "") != RATCHET_NAME:
        raise GateError(f"{path}: authority must be {RATCHET_NAME!r}")
    if not isinstance(doc.get("tier2Ceiling"), dict):
        raise GateError(f"{path}: missing tier2Ceiling map")
    if not isinstance(doc.get("expectedRed"), dict):
        raise GateError(f"{path}: missing expectedRed map")
    for name, lines in doc["expectedRed"].items():
        bad = [x for x in (lines or []) if x not in RED_LINES]
        if bad:
            raise GateError(f"{path}: expectedRed[{name}] names non-red-line(s) {bad}")
    return doc


def evaluate(
    measured: Mapping[str, Any], ratchet: Mapping[str, Any]
) -> Tuple[List[Check], Dict[str, Census]]:
    checks: List[Check] = []
    censuses: Dict[str, Census] = {}
    ceilings: Mapping[str, Any] = ratchet.get("tier2Ceiling") or {}
    expected_red: Mapping[str, Any] = ratchet.get("expectedRed") or {}

    def add(name: str, ok: bool, message: str, group: str, **details: Any) -> None:
        checks.append(Check(name=name, ok=ok, message=message, group=group, details=details))

    for name in sorted(measured):
        result = measured[name]
        if not result.get("ok"):
            add(
                f"{name}.conversion",
                False,
                f"{name}: conversion failed ({result.get('error')})",
                name,
            )
            continue
        try:
            c = Census.from_result(result)
        except GateError as exc:
            add(f"{name}.census", False, f"{name}: {exc}", name)
            continue
        censuses[name] = c

        violated = c.red()
        allowed = set(expected_red.get(name) or [])
        unexpected = sorted(set(violated) - allowed)
        if unexpected:
            add(
                f"{name}.redline",
                False,
                f"{name}: red line(s) {', '.join(f'{k}={violated[k]}' for k in unexpected)} "
                f"— not in expectedRed{sorted(allowed) if allowed else ''}",
                name,
                violated=violated,
                allowed=sorted(allowed),
            )
        elif violated:
            add(
                f"{name}.redline",
                True,
                f"{name}: EXPECTED-RED {', '.join(f'{k}={violated[k]}' for k in sorted(violated))} "
                f"(recorded in {RATCHET_NAME}; do not widen)",
                name,
                violated=violated,
                allowed=sorted(allowed),
            )
        else:
            add(f"{name}.redline", True, f"{name}: red lines clear", name)
            if allowed:
                add(
                    f"{name}.redlineCleared",
                    True,
                    f"RATCHET-CLEAR: {name} no longer violates {sorted(allowed)}; "
                    f"delete its expectedRed entry in tests/gates/baseline/{RATCHET_NAME}",
                    name,
                )

        if name not in ceilings:
            add(
                f"{name}.tier2",
                False,
                f"{name}: no tier2Ceiling in {RATCHET_NAME} (a fixture the ratchet "
                f"never saw is not measured); seed it with --measure",
                name,
                polylineTier2=c.polylineTier2,
            )
            continue
        ceil = int(ceilings[name])
        if c.polylineTier2 > ceil:
            add(
                f"{name}.tier2",
                False,
                f"{name}: polylineTier2={c.polylineTier2} > ceiling={ceil} — an edge that "
                f"was analytic now ships as a mesh polyline",
                name,
                polylineTier2=c.polylineTier2,
                ceiling=ceil,
            )
        elif c.polylineTier2 < ceil:
            add(
                f"{name}.tier2",
                True,
                f"RATCHET-TIGHTEN: {name} polylineTier2={c.polylineTier2} < committed {ceil}; "
                f"lower it in tests/gates/baseline/{RATCHET_NAME} with a note",
                name,
                polylineTier2=c.polylineTier2,
                ceiling=ceil,
            )
        else:
            add(
                f"{name}.tier2",
                True,
                f"{name}: polylineTier2={c.polylineTier2} == ceiling",
                name,
                polylineTier2=c.polylineTier2,
                ceiling=ceil,
            )

    stale = sorted(set(expected_red) - set(censuses))
    if stale:
        add(
            "expectedRedFixtures",
            False,
            f"expectedRed names fixture(s) the gate did not measure: {stale}",
            "ratchet",
        )
    return checks, censuses


def format_checks(checks: Sequence[Check]) -> str:
    return "\n".join(
        f"  [{c.group}/{'PASS' if c.ok else 'FAIL'}] {c.message}" for c in checks
    )


def failing_names(checks: Sequence[Check]) -> List[str]:
    return [c.name for c in checks if not c.ok]


def totals(censuses: Mapping[str, Census]) -> Dict[str, int]:
    t = {k: 0 for k in ("analytic", "polylineTier2") + RED_LINES}
    for c in censuses.values():
        for k in t:
            t[k] += getattr(c, k)
    return t


def run_live(binary: Path, ratchet_path: Path, jobs: int) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    ratchet = load_ratchet(ratchet_path)
    measured = measure(binary, jobs)
    checks, censuses = evaluate(measured, ratchet)
    print(f"edge_class_gate  ratchet={ratchet_path}  fixtures={len(measured)}", flush=True)
    print(format_checks(checks), flush=True)
    t = totals(censuses)
    print(
        "edge_class_gate TOTALS "
        + " ".join(f"{k}={v}" for k, v in t.items()),
        flush=True,
    )
    named = failing_names(checks)
    if named:
        print(
            f"edge_class_gate FAIL  failing assertion(s): {', '.join(named)}",
            file=sys.stderr,
            flush=True,
        )
        return 1
    print("edge_class_gate PASS", flush=True)
    return 0


def run_measure(binary: Path, jobs: int) -> int:
    measured = measure(binary, jobs)
    seed_ceiling: Dict[str, int] = {}
    seed_red: Dict[str, List[str]] = {}
    rows: Dict[str, Any] = {}
    bad = []
    for name in sorted(measured):
        result = measured[name]
        if not result.get("ok"):
            bad.append(name)
            rows[name] = {"error": result.get("error")}
            continue
        try:
            c = Census.from_result(result)
        except GateError as exc:
            bad.append(name)
            rows[name] = {"error": str(exc)}
            continue
        rows[name] = c.asdict()
        seed_ceiling[name] = c.polylineTier2
        red = sorted(c.red())
        if red:
            seed_red[name] = red
    print("edge_class_gate --measure")
    print(json.dumps({"perFixture": rows}, indent=2, sort_keys=True))
    print("SEED tier2Ceiling = " + json.dumps(seed_ceiling, sort_keys=True))
    print("SEED expectedRed  = " + json.dumps(seed_red, sort_keys=True))
    if bad:
        print(f"edge_class_gate FAIL: not measured on {bad}", file=sys.stderr)
        return 1
    return 0


# --------------------------------------------------------------- self-test


def synthetic_result(analytic=10, tier2=3, unhandled=0, over_tol=0, over_cap=0):
    return {
        "ok": True,
        "edgeClasses": {
            "analytic": analytic,
            "polylineTier2": tier2,
            "unhandled": unhandled,
            "overTol": over_tol,
            "overCap": over_cap,
        },
    }


def synthetic_ratchet(**overrides: Any) -> Dict[str, Any]:
    doc: Dict[str, Any] = {
        "id": "edge-class-ratchet",
        "authority": RATCHET_NAME,
        "tier2Ceiling": {"alpha": 3, "beta": 0},
        "expectedRed": {"beta": ["unhandled"]},
        "notes": "synthetic",
    }
    doc.update(overrides)
    return doc


def run_synthetic_pass() -> int:
    ratchet = synthetic_ratchet()
    measured = {
        "alpha": synthetic_result(),
        "beta": synthetic_result(analytic=4, tier2=0, unhandled=2),
    }
    checks, censuses = evaluate(measured, ratchet)
    print("edge_class_gate --synthetic-pass")
    print(format_checks(checks))
    named = failing_names(checks)
    if named:
        print(f"edge_class_gate FAIL  synthetic should pass; failed: {named}", file=sys.stderr)
        return 1
    print("edge_class_gate PASS  synthetic red lines met")
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

    # `r` lists `beta` in expectedRed, so a measurement that does not include
    # `beta` legitimately trips the staleness check. Cases that are about ONE
    # fixture use `r1`, which lists nothing; the staleness check has its own case.
    r = synthetic_ratchet()
    r1 = synthetic_ratchet(expectedRed={})
    ok, _ = evaluate({"alpha": synthetic_result()}, r1)
    check(failing_names(ok) == [], "clean fixture at its ceiling PASSES")

    grew, _ = evaluate({"alpha": synthetic_result(tier2=4)}, r1)
    check("alpha.tier2" in failing_names(grew), "tier-2 increase FAILS")
    check(
        any("polylineTier2=4" in c.message and "ceiling=3" in c.message for c in grew),
        "tier-2 increase names both numbers",
    )

    shrank, _ = evaluate({"alpha": synthetic_result(tier2=1)}, r1)
    check(failing_names(shrank) == [], "tier-2 decrease PASSES")
    check(
        any("RATCHET-TIGHTEN" in c.message for c in shrank),
        "tier-2 decrease prints the re-baselining instruction",
    )

    for line in RED_LINES:
        kw = {"unhandled": {"unhandled": 1}, "overTol": {"over_tol": 1}, "overCap": {"over_cap": 1}}[line]
        red, _ = evaluate({"alpha": synthetic_result(**kw)}, r1)
        check(f"alpha.redline" in failing_names(red), f"{line}>0 on an unlisted fixture FAILS")
        check(any(line in c.message for c in red), f"{line} is named in the failure")

    both = {"alpha": synthetic_result(), "beta": synthetic_result(analytic=4, tier2=0, unhandled=5)}
    xr, _ = evaluate(both, r)
    check(failing_names(xr) == [], "listed expected-red fixture PASSES on its listed line")
    check(any("EXPECTED-RED" in c.message for c in xr), "expected-red stays loud")

    xr2, _ = evaluate({"alpha": synthetic_result(),
                       "beta": synthetic_result(analytic=4, tier2=0, over_tol=1)}, r)
    check(
        "beta.redline" in failing_names(xr2),
        "a listed fixture violating an UNLISTED line still FAILS",
    )

    healed, _ = evaluate({"alpha": synthetic_result(),
                          "beta": synthetic_result(analytic=4, tier2=0)}, r)
    check(failing_names(healed) == [], "a healed expected-red fixture PASSES")
    check(
        any("RATCHET-CLEAR" in c.message for c in healed),
        "a healed expected-red prints the delete instruction",
    )

    unseen, _ = evaluate({"gamma": synthetic_result()}, r1)
    check("gamma.tier2" in failing_names(unseen), "a fixture with no ceiling FAILS (not measured)")

    stale, _ = evaluate({"alpha": synthetic_result()}, r)
    check(
        "expectedRedFixtures" in failing_names(stale),
        "expectedRed naming an unmeasured fixture FAILS",
    )

    nocensus, _ = evaluate({"alpha": {"ok": True}}, r)
    check("alpha.census" in failing_names(nocensus), "RESULT without edgeClasses FAILS")
    partial, _ = evaluate({"alpha": {"ok": True, "edgeClasses": {"analytic": 1}}}, r)
    check("alpha.census" in failing_names(partial), "edgeClasses missing keys FAILS")

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "not-the-ratchet.json"
        p.write_text("{}", encoding="utf-8")
        try:
            load_ratchet(p)
            check(False, "mis-named ratchet file must FAIL")
        except GateError:
            check(True, "mis-named ratchet file FAIL")
        q = Path(td) / RATCHET_NAME
        q.write_text(
            json.dumps(synthetic_ratchet(expectedRed={"beta": ["notALine"]})), encoding="utf-8"
        )
        try:
            load_ratchet(q)
            check(False, "expectedRed with a bogus line name must FAIL")
        except GateError:
            check(True, "expectedRed with a bogus line name FAIL")

    check(
        parse_result_line('RESULT {"ok":true}\n')["ok"] is True, "parse_result_line"
    )
    check(run_synthetic_pass() == 0, "--synthetic-pass returns 0")
    fx = dict(fixtures())
    check(len(fx) >= 32, f"fixture list covers the corpus (+ build fixtures): {len(fx)}")
    check(DEFAULT_RATCHET.is_file(), f"{RATCHET_NAME} present")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--synthetic-pass", action="store_true")
    p.add_argument(
        "--measure",
        action="store_true",
        help="print the live per-fixture census and the seed values (writes nothing)",
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
            print("edge_class_gate: --binary is required", file=sys.stderr)
            return 1
        if args.measure:
            return run_measure(Path(args.binary), args.jobs)
        return run_live(Path(args.binary), Path(args.ratchet), args.jobs)
    except GateError as exc:
        print(f"edge_class_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
