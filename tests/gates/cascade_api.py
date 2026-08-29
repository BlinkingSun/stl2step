#!/usr/bin/env python3
"""cascade_api — SPEC C3 three-state FAIL_RID via the public C++ convert() API.

Spawns isolated processes (FAIL_RID is parsed-once/cached) and runs them
concurrently. Each child calls stl2step::convert().
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_STL = REPO / "tests" / "corpus" / "handle-lock.stl"


def run_probe(probe: Path, stl: Path, out_dir: Path, expect: str, fail_rid: str | None) -> int:
    env = os.environ.copy()
    env["CASCADE_EXPECT"] = expect
    env["CASCADE_STL"] = str(stl)
    env["CASCADE_OUT"] = str(out_dir)
    if fail_rid is None:
        env.pop("STL2STEP_FAIL_RID", None)
    else:
        env["STL2STEP_FAIL_RID"] = fail_rid
    proc = subprocess.run([str(probe)], capture_output=True, text=True, env=env)
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return proc.returncode


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--probe", type=Path, required=True, help="cascade_api_probe binary")
    p.add_argument("--stl", type=Path, default=DEFAULT_STL)
    args = p.parse_args(argv)
    probe = args.probe.resolve()
    stl = args.stl.resolve()
    if not probe.is_file():
        print(f"cascade_api FAIL: probe missing {probe}", file=sys.stderr)
        return 1
    if not stl.is_file():
        print(f"cascade_api FAIL: STL missing {stl}", file=sys.stderr)
        return 1

    cases = [
        ("unset", None),
        ("11", "11"),
        ("0", "0"),
        ("bad", "not-an-int"),
    ]
    with tempfile.TemporaryDirectory(prefix="cascade_api_") as td:
        root = Path(td)
        with ThreadPoolExecutor(max_workers=len(cases)) as pool:
            futs = [
                pool.submit(run_probe, probe, stl, root / expect, expect, rid)
                for expect, rid in cases
            ]
            codes = [f.result() for f in futs]
    if any(c != 0 for c in codes):
        print(f"cascade_api FAIL codes={codes}", file=sys.stderr)
        return 1
    print("cascade_api PASS C3 unset/11/0 + malformed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
