#!/usr/bin/env python3
"""p1_ac2_guard_selftest — the AC2-A1 guard against generator-exact fixtures.

Hosted ctest does not echo a passing test's stdout, so run 34171646888's
job log has the gates_full / p1_ac2 timing lines and the suite summary,
and zero "G1 FAIL" lines (measured). Those lines are not invented into
linux-real-nest.log. The parser contract for the three baseline names is
checked on an in-memory transcript appended to that verbatim excerpt.
"""
from __future__ import annotations

import io
import re
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_gates_full_baseline as guard  # noqa: E402

FIX = HERE / "fixtures" / "ac2"

# D-train-push-4, copied: the nest is real iff b >= a/1.5.
_TIME = re.compile(
    r"Test #\d+:\s+(gates_full|p1_ac2_gates_baseline)\b.*?"
    r"(Passed|\*\*\*Timeout|\*\*\*Skipped|\*\*\*Failed|\*\*\*Not Run)"
    r"\s+([0-9]+(?:\.[0-9]+)?)\s+sec"
)

# Prefix-stripped verbatim lines from run 34343987075 job 102441119099.
_WIN_FALSE = """\
11/87 Test #11: gates_full ...........................   Passed  1952.55 sec
57/87 Test #57: p1_ac2_gates_baseline ................   Passed    0.33 sec
100% tests passed out of 87
"""


def check(cond: bool, msg: str) -> None:
    if not cond:
        print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
        raise SystemExit(1)


def nest_ratio_ok(a: float, b: float) -> bool:
    """D-train-push-4 floor. The only recorded factor, read as a minimum."""
    return b >= a / 1.5


def scrape_times(text: str) -> dict[str, tuple[str, str]]:
    found: dict[str, tuple[str, str]] = {}
    for line in text.splitlines():
        m = _TIME.search(line)
        if m:
            found[m.group(1)] = (m.group(2), m.group(3))
    return found


def call_main(text: str, config: str | None) -> tuple[int, str, str, list]:
    seen: list = []

    def fake_run(cmd, **kwargs):
        seen.append(list(cmd))

        class P:
            returncode = 0
            stdout = text
            stderr = ""

        return P()

    argv = ["--build-dir", "/nonexistent"]
    if config is not None:
        argv += ["--config", config]
    orig = guard.subprocess.run
    guard.subprocess.run = fake_run
    out, err = io.StringIO(), io.StringIO()
    try:
        with redirect_stdout(out), redirect_stderr(err):
            rc = guard.main(argv)
    finally:
        guard.subprocess.run = orig
    return rc, out.getvalue(), err.getvalue(), seen


def main() -> int:
    linux = (FIX / "linux-real-nest.log").read_text(encoding="utf-8")
    vacuous = (FIX / "vacuous-no-tests.log").read_text(encoding="utf-8")
    skipped = (FIX / "skipped-gates-full.log").read_text(encoding="utf-8")
    new_fail = (FIX / "new-g1-fail.log").read_text(encoding="utf-8")

    check(guard.nested_ran(linux) is True, "linux excerpt: nested_ran")
    check("G1 FAIL" not in linux, "linux excerpt must stay free of invented G1 FAIL lines")
    check(guard.parse_g1_fails(linux) == set(),
          f"linux excerpt parse got {guard.parse_g1_fails(linux)!r}, want empty")
    rc, out, err, cmd = call_main(linux, None)
    check(rc == 0, f"linux verdict rc={rc} err={err}")
    check("AC2-A1 guard PASS" in out, "linux verdict did not PASS")
    check("-C" not in cmd[0], f"empty config added a -C flag: {cmd[0]}")

    transcript = linux + "S09 G1 FAIL\nS15 G1 FAIL\nS16-R1-round-2 G1 FAIL\n"
    check(
        guard.parse_g1_fails(transcript) == {"S09", "S15", "S16-R1-round-2"},
        f"parser contract got {guard.parse_g1_fails(transcript)!r}",
    )
    rc, out, err, cmd = call_main(transcript, "")
    check(rc == 0 and "AC2-A1 guard PASS" in out and "NEW G1 FAIL" not in err,
          f"baseline names must not be extra rc={rc} out={out} err={err}")
    check("-C" not in cmd[0], f"explicit empty --config must add no -C: {cmd[0]}")

    rc, out, err, cmd = call_main(vacuous, None)
    check(guard.nested_ran(vacuous) is False, "vacuous nested_ran should be false")
    check(rc == 1, f"vacuous main rc={rc}")
    check("nested ctest ran no gates_full (config=)" in err, f"vacuous message missing: {err}")
    check("No tests were found!!!" in err, "vacuous tail was not printed")

    check(guard.nested_ran(skipped) is True, "skipped gates_full must still count as nested_ran")
    rc, out, err, _ = call_main(skipped, None)
    check(rc == 0 and "AC2-A1 guard PASS" in out,
          f"skipped path became a failure rc={rc} err={err}")
    check("nested ctest ran no gates_full" not in err, "skip was treated as vacuous")

    rc, out, err, _ = call_main(new_fail, None)
    check(rc == 1, f"new G1 FAIL rc={rc}")
    check("S11" in err, f"extra fixture S11 not named: {err}")
    check("NEW G1 FAIL outside baseline" in err, f"extra message missing: {err}")

    rc, _, _, cmd = call_main(linux, "Release")
    check(cmd[0][-2:] == ["-C", "Release"], f"--config Release argv={cmd[0]}")
    check(rc == 0, f"Release config changed the verdict rc={rc}")

    live = scrape_times(linux)
    check(live.get("gates_full", (None, None))[0] == "Passed", f"linux gates_full {live}")
    check(live.get("p1_ac2_gates_baseline", (None, None))[0] == "Passed", f"linux p1_ac2 {live}")
    a = float(live["gates_full"][1])
    b = float(live["p1_ac2_gates_baseline"][1])
    check(nest_ratio_ok(a, b), f"real nest {b} >= {a}/1.5 failed")

    false = scrape_times(_WIN_FALSE)
    fa = float(false["gates_full"][1])
    fb = float(false["p1_ac2_gates_baseline"][1])
    check(false["p1_ac2_gates_baseline"][0] == "Passed", "windows false status")
    check(not nest_ratio_ok(fa, fb), f"false green {fb} vs {fa}/1.5 must fail the floor")

    print("p1_ac2_guard_selftest PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
