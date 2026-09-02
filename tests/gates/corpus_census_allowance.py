#!/usr/bin/env python3
"""corpus_census_allowance — anti-silent-waiver for censusValidExpected:false.

The set of corpus parts with any live[] row censusValidExpected == false
must be exactly {Body11}. Each such part must have a non-empty
censusValidReason and a non-zero builtCylindersFloor.

Usage:
  corpus_census_allowance.py
  corpus_census_allowance.py --corpus DIR
  corpus_census_allowance.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = REPO / "tests" / "corpus"
ALLOWED_FALSE: Set[str] = {"Body11"}


class GateError(Exception):
    """Hard gate failure."""


def load_sidecars(corpus: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for p in sorted(corpus.glob("*.expected.json")):
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise GateError(f"sidecar unreadable {p}: {exc}") from exc
        out[p.name[: -len(".expected.json")]] = doc
    if not out:
        raise GateError(f"no *.expected.json under {corpus}")
    return out


def false_parts(sidecars: Dict[str, Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    found: Dict[str, List[Dict[str, Any]]] = {}
    for pid, doc in sidecars.items():
        rows = [
            row
            for row in (doc.get("live") or [])
            if isinstance(row, dict) and row.get("censusValidExpected") is False
        ]
        if rows:
            found[pid] = rows
    return found


def evaluate(sidecars: Dict[str, Dict[str, Any]]) -> List[str]:
    fails: List[str] = []
    found = false_parts(sidecars)
    got = set(found)
    if got != ALLOWED_FALSE:
        extra = sorted(got - ALLOWED_FALSE)
        missing = sorted(ALLOWED_FALSE - got)
        if extra:
            fails.append(
                f"censusValidExpected:false parts {extra} not in allowed {sorted(ALLOWED_FALSE)}"
            )
        if missing:
            fails.append(f"allowed parts missing censusValidExpected:false: {missing}")
    for pid, rows in sorted(found.items()):
        reasons = [str(r.get("censusValidReason") or "").strip() for r in rows]
        if not any(reasons):
            fails.append(f"{pid}: censusValidReason empty")
        floor = max((int(r.get("builtCylindersFloor") or 0) for r in rows), default=0)
        if floor <= 0:
            fails.append(f"{pid}: builtCylindersFloor={floor} must be non-zero")
    return fails


def run_gate(corpus: Path) -> int:
    sidecars = load_sidecars(corpus)
    fails = evaluate(sidecars)
    found = false_parts(sidecars)
    print(
        f"corpus_census_allowance  corpus={corpus}  "
        f"false_parts={sorted(found)}  allowed={sorted(ALLOWED_FALSE)}"
    )
    if fails:
        print("corpus_census_allowance FAIL", file=sys.stderr)
        for f in fails:
            print(f"  {f}", file=sys.stderr)
        return 1
    for pid, rows in sorted(found.items()):
        floor = max(int(r.get("builtCylindersFloor") or 0) for r in rows)
        reason = next(
            (str(r.get("censusValidReason") or "") for r in rows if r.get("censusValidReason")),
            "",
        )
        print(f"  {pid}: floor={floor} reason={reason!r}")
    print("corpus_census_allowance PASS")
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

    ok_doc = {
        "live": [
            {
                "component": 0,
                "censusValidExpected": False,
                "censusValidReason": "R2 validity stop: 159 analytic cylinders kept, census BRepCheck invalid pre-R2",
                "builtCylindersFloor": 159,
            }
        ]
    }
    check(evaluate({"Body11": ok_doc}) == [], "Body11-only with reason+floor PASSES")
    extra = evaluate({"Body11": ok_doc, "S09": ok_doc})
    check(any("S09" in f for f in extra), "second false part FAILS")
    missing = evaluate({})
    check(any("missing" in f for f in missing), "missing Body11 FAILS")
    empty_reason = {
        "live": [
            {
                "censusValidExpected": False,
                "censusValidReason": "",
                "builtCylindersFloor": 159,
            }
        ]
    }
    check(
        any("censusValidReason empty" in f for f in evaluate({"Body11": empty_reason})),
        "empty reason FAILS",
    )
    zero_floor = {
        "live": [
            {
                "censusValidExpected": False,
                "censusValidReason": "R2 validity stop",
                "builtCylindersFloor": 0,
            }
        ]
    }
    check(
        any("non-zero" in f for f in evaluate({"Body11": zero_floor})),
        "zero floor FAILS",
    )
    absent = evaluate({"S01": {"live": [{"builtCylindersFloor": 0}]}})
    check(any("missing" in f for f in absent), "absent key is not a false part")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    p.add_argument("--self-test", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        return run_gate(Path(args.corpus))
    except GateError as exc:
        print(f"corpus_census_allowance FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
