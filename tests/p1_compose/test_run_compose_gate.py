#!/usr/bin/env python3
"""p1_compose_gate_selftest — the expected-red XFAIL scope (ruling 17:20, G3).

A fixture listed under expected-red.json p1_compose_gate may XFAIL only its
recoverable[] misses. This drives run_compose_gate.main() end to end with the
real check_regionset.py; only the region dump is stubbed (a committed
RegionSet stands in for stl2step_regiondump), so no engine binary is needed.

Forced cases: a listed fixture whose only failure is a recoverable miss is
XFAIL; a listed fixture with an unlisted failure class (an absent expected
reject, a recognition failure) still FAILS; an unlisted fixture with a
recoverable miss FAILS; a listed fixture that passes prints XPASS.
"""
from __future__ import annotations

import io
import json
import shutil
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import run_compose_gate as gate  # noqa: E402

ICHECKER = ROOT / "tests" / "gates" / "check_regionset.py"
DUMP = ROOT / "tests" / "gates" / "build_fixtures" / "counterbore" / "counterbore.regionset.json"
FID = "S90_forced"
MISS = {"type": "cylinder", "radius": 50.0,
        "axis": {"loc": [0.0, 0.0, 0.0], "dir": [0, 0, 1]}, "count": 1}


def check(cond: bool, msg: str) -> None:
    if not cond:
        print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
        raise SystemExit(1)
    print(f"SELFTEST PASS: {msg}")


def run_gate(sidecar: dict | None, listed: bool, stl_name: str = FID) -> tuple[int, str]:
    with tempfile.TemporaryDirectory(prefix="p1_compose_selftest_") as td:
        tmp = Path(td)
        corpus = tmp / "corpus"
        corpus.mkdir()
        (corpus / f"{stl_name}.stl").write_bytes(b"")
        if sidecar is not None:
            (corpus / f"{stl_name}.expected.json").write_text(json.dumps(sidecar))
        cube = tmp / "cube.stl"
        cube.write_bytes(b"")
        baseline = tmp / "repo" / "tests" / "gates" / "baseline"
        baseline.mkdir(parents=True)
        red = {"p1_compose_gate": {Path(stl_name).stem: "forced row"} if listed else {}}
        (baseline / "expected-red.json").write_text(json.dumps(red))

        orig = (gate.envelope_clean_indices, gate.dump_bare)
        gate.envelope_clean_indices = lambda dump_bin, stl: [0]
        gate.dump_bare = lambda dump_bin, stl, comp, threads, out: shutil.copyfile(DUMP, out)
        argv = sys.argv
        sys.argv = ["run_compose_gate.py", "--repo", str(tmp / "repo"),
                    "--build-dir", str(tmp), "--dump", str(tmp / "no-dump"),
                    "--ichecker", str(ICHECKER), "--corpus", str(corpus),
                    "--cube", str(cube)]
        out, err = io.StringIO(), io.StringIO()
        try:
            with redirect_stdout(out), redirect_stderr(err):
                rc = gate.main()
        finally:
            gate.envelope_clean_indices, gate.dump_bare = orig
            sys.argv = argv
        return rc, out.getvalue() + err.getvalue()


def main() -> int:
    # Classifier on the real I-checker's report lines.
    check(gate.recoverable_misses_only("SIDECAR FAIL recoverable0,recoverable2\n"),
          "recoverable-only SIDECAR line is XFAIL-able")
    check(gate.recoverable_misses_only("SIDECAR_RECOVERABLE FAIL recoverable1\n"),
          "recoverable-only SIDECAR_RECOVERABLE line is XFAIL-able")
    check(not gate.recoverable_misses_only("SIDECAR FAIL rejected,recoverable0\n"),
          "SIDECAR rejected offender is hard")
    check(not gate.recoverable_misses_only("SIDECAR FAIL recoverable0\nI3 FAIL region4\n"),
          "I-rule failure beside a recoverable miss is hard")
    check(not gate.recoverable_misses_only("SCHEMA FAIL\nSIDECAR_RECOVERABLE FAIL recoverable0\n"),
          "SCHEMA failure is hard")
    check(not gate.recoverable_misses_only("Traceback (most recent call last):\n"),
          "non-zero exit with no FAIL line is hard")

    rc, text = run_gate({"recoverable": [MISS]}, listed=True)
    check(rc == 0 and f"XFAIL {FID}.stl comp0: forced row" in text,
          "listed fixture, recoverable miss only -> XFAIL, exit 0")

    rc, text = run_gate({"recoverable": [MISS], "rejected": ["sphereNYI"]}, listed=True)
    check(rc == 1 and f"FAIL {FID}.stl comp0: ichecker" in text and "XFAIL" not in text,
          "listed fixture, absent expected reject (unlisted class) -> FAIL, exit 1")

    rc, text = run_gate({"recoverable": [MISS]}, listed=True, stl_name="S13")
    check(rc == 1 and "FAIL S13.stl comp0: recognition" in text and "XFAIL" not in text,
          "listed fixture, recoverable miss + recognition failure -> FAIL, exit 1")

    rc, text = run_gate({"recoverable": [MISS]}, listed=False)
    check(rc == 1 and f"FAIL {FID}.stl comp0: ichecker" in text,
          "unlisted fixture, recoverable miss -> FAIL, exit 1")

    rc, text = run_gate({}, listed=True)
    check(rc == 0 and f"XPASS {FID}" in text,
          "listed fixture that passes -> XPASS, exit 0")

    print("p1_compose_gate_selftest PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
