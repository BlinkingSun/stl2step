#!/usr/bin/env python3
"""cascade_no_waiver — fail the build if forbidden acceptance predicates reappear.

Greps src/ for keepValidPartials and for nPartial-conditioned acceptance
(D1 §5.1). API self-test needs no engine.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
SRC_SUFFIXES = {".cpp", ".hpp", ".h", ".cc", ".cxx"}
NPARTIAL_ACCEPT = re.compile(r"\bnPartial\b\s*[<>=!]")


def scan(src: Path) -> list[str]:
    fails: list[str] = []
    if not src.is_dir():
        return [f"src missing: {src}"]
    for p in sorted(src.rglob("*")):
        if not p.is_file() or p.suffix not in SRC_SUFFIXES:
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        rel = p.relative_to(src.parent)
        if "keepValidPartials" in text:
            fails.append(f"{rel}: keepValidPartials")
        for i, line in enumerate(text.splitlines(), 1):
            if NPARTIAL_ACCEPT.search(line):
                fails.append(f"{rel}:{i}: nPartial-conditioned acceptance: {line.strip()}")
    return fails


def _self_test() -> int:
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    hits = scan(SRC)
    check(hits == [], f"src/ is clean (hits={hits})")
    check(NPARTIAL_ACCEPT.search("if (nPartial >= 10)") is not None, "detects nPartial >= 10")
    check(NPARTIAL_ACCEPT.search("int nPartialFaces = 0;") is None, "ignores nPartialFaces name")
    check("keepValidPartials" in "bool keepValidPartials = true;", "detects keepValidPartials")
    return 1 if fails else 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--src", type=Path, default=SRC)
    args = p.parse_args(argv)
    if args.self_test:
        return _self_test()
    hits = scan(args.src)
    if hits:
        print("cascade_no_waiver FAIL", file=sys.stderr)
        for h in hits:
            print(f"  {h}", file=sys.stderr)
        return 1
    print("cascade_no_waiver PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
