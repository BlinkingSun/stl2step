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


def nested_ran(text: str) -> bool:
    return bool(re.search(r"Test\s+#\d+:\s+gates_full", text)) and bool(
        re.search(r"[0-9]+% tests passed", text)
    )


def evaluate(text: str, cfg: str) -> int:
    """Guard body. `cfg` is the nested ctest -C value, empty when unset."""
    if not nested_ran(text):
        tail = text.splitlines()[-20:]
        if tail:
            print("\n".join(tail), file=sys.stderr)
        print(
            f"nested ctest ran no gates_full (config={cfg}) — AC2-A1 did not execute",
            file=sys.stderr,
        )
        return 1
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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True)
    ap.add_argument("--config", default="")
    args = ap.parse_args(argv)
    cmd = [
        "ctest", "--test-dir", args.build_dir, "-R", "^gates_full$",
        "--output-on-failure",
    ] + (["-C", args.config] if args.config else [])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return evaluate(text, args.config)


if __name__ == "__main__":
    sys.exit(main())
