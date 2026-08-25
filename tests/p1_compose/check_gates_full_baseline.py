#!/usr/bin/env python3
"""AC2-A1 guard: ctest -R '^gates_full$' failure list must match the inherited baseline.

Baseline (audit-p1-compose-fix / faceted 1.0.0 engine in src/stl2step.cpp):
  G1 FAIL on S09 (ShapeFix rewrite), S15 (open), S16-R1-round-2 (open).
Any G1 FAIL outside that set is this lane's and is blocking.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys

BASELINE_G1_FAIL = frozenset({"S09", "S15", "S16-R1-round-2"})


def parse_g1_fails(text: str) -> set[str]:
    current = None
    fails: set[str] = set()
    fx = re.compile(r"\b(S\d+(?:-[A-Za-z0-9]+)*)\b")
    for line in text.splitlines():
        for m in fx.finditer(line):
            current = m.group(1)
        if re.search(r"\bG1\s+FAIL\b", line) and current:
            fails.add(current)
    return fails


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True)
    args = ap.parse_args()
    proc = subprocess.run(
        ["ctest", "--test-dir", args.build_dir, "-R", "^gates_full$",
         "--output-on-failure"],
        capture_output=True, text=True,
    )
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    got = parse_g1_fails(text)
    print("gates_full G1 FAIL fixtures:", sorted(got))
    print("inherited baseline:         ", sorted(BASELINE_G1_FAIL))
    extra = got - BASELINE_G1_FAIL
    missing = BASELINE_G1_FAIL - got
    if extra:
        print("NEW G1 FAIL outside baseline (blocking):", sorted(extra),
              file=sys.stderr)
        return 1
    if missing:
        print("note: baseline fixtures not observed as G1 FAIL:", sorted(missing))
    print("AC2-A1 guard PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
