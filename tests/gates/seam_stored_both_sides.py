#!/usr/bin/env python3
"""seam_stored_both_sides — D3-2: every seam stored on both adjacent surfaces.

Expected-red at dc49e45: DIAG_SEAM / DIAG_G7 / DIAG_CYLG7 are emitted by
h1-seams behind STL2STEP_P2_DIAG=1 and do not exist at this tree. Zero
DIAG_SEAM lines is **not measured**, never a silent pass.

CMake EXPECTED-RED (copy of partial_recovery_gate protocol):
  SEAM_STORED_BOTH_SIDES_STRICT defaults ON (unset/empty => strict).
  PASS_REGULAR_EXPRESSION inverts only the named expected-red set.
  Anything outside that set prints SEAM_STORED_BOTH_SIDES_UNEXPECTED.

FLIP PROTOCOL (the day DIAG_SEAM lands with storedOwner=storedConsumer=1):
  1. Live Python prints SEAM_STORED_BOTH_SIDES_FLOORS_MET and exits 0.
  2. ctest FAILS: Required regular expression not found.
  3. Remove PASS_REGULAR_EXPRESSION from tests/gates/CMakeLists.txt.
  4. Keep ENVIRONMENT SEAM_STORED_BOTH_SIDES_STRICT=1. Do not weaken asserts.

Usage:
  seam_stored_both_sides.py --self-test
  seam_stored_both_sides.py --synthetic-pass
  seam_stored_both_sides.py --binary ./build/stl2step [--jobs 2]
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
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
PICKUP_STL = REPO / "tests" / "corpus" / "handle-pickup.stl"
HL_STL = REPO / "tests" / "corpus" / "handle-lock.stl"

ENV_STRICT = "SEAM_STORED_BOTH_SIDES_STRICT"
MARKER_EXPECTED_RED = "SEAM_STORED_BOTH_SIDES_EXPECTED_RED"
MARKER_UNEXPECTED = "SEAM_STORED_BOTH_SIDES_UNEXPECTED"
MARKER_FLOORS_MET = "SEAM_STORED_BOTH_SIDES_FLOORS_MET"

# HEAD @ dc49e45: instrument not yet emitted. These names are the invert set.
EXPECTED_RED_NAMES = frozenset(
    {
        "nSeam",
        "storedOwner",
        "storedConsumer",
        "nCreated",
        "covered",
        "g7",
        "g7Coverage",
        "handle-lock.storedNotMeasured",
    }
)

SEAM_REQUIRED = (
    "ci",
    "owner",
    "consumer",
    "ownerType",
    "consumerType",
    "kind",
    "f",
    "l",
    "storedOwner",
    "storedConsumer",
)
G7_REQUIRED = (
    "rid",
    "nEdgesIter",
    "nEdgesExp",
    "nBound",
    "nSkippedHasPc",
    "nSkippedNoCurve",
    "covered",
    "g7",
)
CYLG7_REQUIRED = G7_REQUIRED + ("nCreated",)
OWNER_TYPES = frozenset({"plane", "cyl", "facet"})
KIND_TYPES = frozenset({"circ", "lin", "elips", "poly"})


class GateError(Exception):
    """Hard gate failure (binary missing, no RESULT)."""


@dataclass
class Check:
    name: str
    ok: bool
    message: str
    group: str
    hard: bool = True
    details: Dict[str, Any] = field(default_factory=dict)


def strict_enabled(env: Optional[Mapping[str, str]] = None) -> bool:
    """SEAM_STORED_BOTH_SIDES_STRICT defaults ON (unset / empty => strict)."""
    environ: Mapping[str, str] = os.environ if env is None else env
    raw = environ.get(ENV_STRICT)
    if raw is None or str(raw).strip() == "":
        return True
    return str(raw).strip().lower() not in ("0", "false", "off", "no")


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise GateError("no RESULT line in stl2step stdout")


def parse_diag_kv(line: str) -> Optional[Tuple[str, Dict[str, str]]]:
    """Parse `DIAG_* key=value ...`. Extra keys tolerated; no fallbacks."""
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


def require_keys(kv: Mapping[str, str], required: Sequence[str], kind: str) -> Optional[str]:
    missing = [k for k in required if k not in kv]
    if missing:
        return f"{kind}: missing key(s) {missing} (not measured)"
    return None


def as_int(kv: Mapping[str, str], key: str) -> int:
    return int(float(kv[key]))


def collect_diag(stderr: str, kind: str) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    for line in (stderr or "").splitlines():
        parsed = parse_diag_kv(line)
        if parsed is None:
            continue
        k, kv = parsed
        if k == kind:
            out.append(kv)
    return out


def evaluate_pickup(stderr: str, *, strict: bool = True) -> List[Check]:
    checks: List[Check] = []

    def add(
        name: str,
        ok: bool,
        message: str,
        hard: Optional[bool] = None,
        **details: Any,
    ) -> None:
        checks.append(
            Check(
                name=name,
                ok=ok,
                message=message,
                group="pickup",
                hard=strict if hard is None else hard,
                details=details,
            )
        )

    seams = collect_diag(stderr, "DIAG_SEAM")
    g7s = collect_diag(stderr, "DIAG_G7")
    cyls = collect_diag(stderr, "DIAG_CYLG7")

    n_seam = len(seams)
    if n_seam == 0:
        add(
            "nSeam",
            False,
            "nSeam=0 DIAG_SEAM lines absent -- not measured",
            nSeam=0,
        )
    else:
        add("nSeam", True, f"nSeam={n_seam}", nSeam=n_seam)

    stored_owner_bad: List[int] = []
    stored_consumer_bad: List[int] = []
    missing_seam = 0
    for kv in seams:
        miss = require_keys(kv, SEAM_REQUIRED, "DIAG_SEAM")
        if miss:
            missing_seam += 1
            add("nSeam", False, miss)
            continue
        if as_int(kv, "storedOwner") != 1:
            stored_owner_bad.append(as_int(kv, "ci"))
        if as_int(kv, "storedConsumer") != 1:
            stored_consumer_bad.append(as_int(kv, "ci"))
        ot, ct = kv["ownerType"], kv["consumerType"]
        if ot not in OWNER_TYPES or ct not in OWNER_TYPES:
            add(
                "seamGrammar",
                False,
                f"DIAG_SEAM ci={kv.get('ci')} ownerType/consumerType "
                f"not in {sorted(OWNER_TYPES)} (got {ot}/{ct})",
                hard=True,
            )
        if kv["kind"] not in KIND_TYPES:
            add(
                "seamGrammar",
                False,
                f"DIAG_SEAM ci={kv.get('ci')} kind={kv['kind']!r} "
                f"not in {sorted(KIND_TYPES)}",
                hard=True,
            )

    if n_seam > 0:
        add(
            "storedOwner",
            not stored_owner_bad and missing_seam == 0,
            (
                f"storedOwner!=1 on ci={stored_owner_bad[:12]}"
                if stored_owner_bad
                else "every DIAG_SEAM storedOwner=1"
            ),
            bad=stored_owner_bad[:12],
        )
        add(
            "storedConsumer",
            not stored_consumer_bad and missing_seam == 0,
            (
                f"storedConsumer!=1 on ci={stored_consumer_bad[:12]}"
                if stored_consumer_bad
                else "every DIAG_SEAM storedConsumer=1"
            ),
            bad=stored_consumer_bad[:12],
        )

    created_bad: List[int] = []
    if not cyls and n_seam == 0:
        add(
            "nCreated",
            False,
            "nCreated: zero DIAG_CYLG7 lines -- not measured",
        )
    for kv in cyls:
        miss = require_keys(kv, CYLG7_REQUIRED, "DIAG_CYLG7")
        if miss:
            add("nCreated", False, miss)
            continue
        if as_int(kv, "nCreated") != 0:
            created_bad.append(as_int(kv, "rid"))
    if cyls:
        add(
            "nCreated",
            not created_bad,
            (
                f"DIAG_CYLG7 nCreated!=0 on rid={created_bad[:12]} "
                "(D3-2: consumers verify, never create)"
                if created_bad
                else "every DIAG_CYLG7 nCreated=0"
            ),
            bad=created_bad[:12],
        )

    cover_bad: List[str] = []
    g7_bad: List[str] = []
    covsum_bad: List[str] = []
    g7_lines = [("DIAG_G7", g7s, G7_REQUIRED), ("DIAG_CYLG7", cyls, CYLG7_REQUIRED)]
    any_g7 = bool(g7s or cyls)
    if not any_g7:
        add("covered", False, "covered: zero DIAG_G7 / DIAG_CYLG7 lines -- not measured")
        add("g7", False, "g7: zero DIAG_G7 / DIAG_CYLG7 lines -- not measured")
        add(
            "g7Coverage",
            False,
            "g7Coverage: zero DIAG_G7 / DIAG_CYLG7 lines -- not measured "
            "(nBound+nSkippedHasPc+nSkippedNoCurve == nEdgesIter)",
        )
    for kind, rows, req in g7_lines:
        for kv in rows:
            miss = require_keys(kv, req, kind)
            if miss:
                cover_bad.append(miss)
                continue
            n_iter = as_int(kv, "nEdgesIter")
            covered = as_int(kv, "covered")
            g7 = as_int(kv, "g7")
            n_bound = as_int(kv, "nBound")
            n_pc = as_int(kv, "nSkippedHasPc")
            n_nc = as_int(kv, "nSkippedNoCurve")
            rid = kv.get("rid", "?")
            if covered != n_iter:
                cover_bad.append(f"{kind} rid={rid} covered={covered} != nEdgesIter={n_iter}")
            if g7 != 1:
                g7_bad.append(f"{kind} rid={rid} g7={g7} != 1")
            if n_bound + n_pc + n_nc != n_iter:
                covsum_bad.append(
                    f"{kind} rid={rid} nBound+nSkippedHasPc+nSkippedNoCurve="
                    f"{n_bound + n_pc + n_nc} != nEdgesIter={n_iter}"
                )
    if any_g7:
        add(
            "covered",
            not cover_bad,
            (
                f"G-7 covered!=nEdgesIter: {cover_bad[:6]}"
                if cover_bad
                else "every DIAG_G7/DIAG_CYLG7 covered==nEdgesIter"
            ),
        )
        add(
            "g7",
            not g7_bad,
            f"G-7 g7!=1: {g7_bad[:6]}" if g7_bad else "every DIAG_G7/DIAG_CYLG7 g7=1",
        )
        add(
            "g7Coverage",
            not covsum_bad,
            (
                f"G-7 nBound+nSkippedHasPc+nSkippedNoCurve != nEdgesIter: {covsum_bad[:6]}"
                if covsum_bad
                else "every wire nBound+nSkippedHasPc+nSkippedNoCurve==nEdgesIter"
            ),
        )
    return checks


def evaluate_handle_lock(stderr: str, *, strict: bool = True) -> List[Check]:
    """handle-lock must not regress: no DIAG_SEAM line with a 0 stored flag."""
    checks: List[Check] = []
    seams = collect_diag(stderr, "DIAG_SEAM")
    if not seams:
        # Instrument absent: expected-red at HEAD. A later 0-stored flag is a
        # distinct unexpected fail (handle-lock.stored).
        checks.append(
            Check(
                name="handle-lock.storedNotMeasured",
                ok=False,
                message="handle-lock: zero DIAG_SEAM lines -- not measured",
                group="handle-lock",
                hard=strict,
            )
        )
        return checks
    bad: List[int] = []
    for kv in seams:
        miss = require_keys(kv, SEAM_REQUIRED, "DIAG_SEAM")
        if miss:
            checks.append(
                Check(
                    name="handle-lock.stored",
                    ok=False,
                    message=f"handle-lock {miss}",
                    group="handle-lock",
                    hard=True,
                )
            )
            continue
        if as_int(kv, "storedOwner") != 1 or as_int(kv, "storedConsumer") != 1:
            bad.append(as_int(kv, "ci"))
    checks.append(
        Check(
            name="handle-lock.stored",
            ok=not bad,
            message=(
                f"handle-lock DIAG_SEAM stored flag 0 on ci={bad[:12]}"
                if bad
                else f"handle-lock nSeam={len(seams)} storedOwner=storedConsumer=1"
            ),
            group="handle-lock",
            hard=True,
            details={"bad": bad[:12], "nSeam": len(seams)},
        )
    )
    return checks


def synthetic_pass_stderr() -> str:
    """A fixture that satisfies D3-2 / G-7 on both bodies."""
    return "\n".join(
        [
            "DIAG_SEAM ci=0 owner=1 consumer=2 ownerType=plane consumerType=cyl "
            "kind=circ f=0 l=6.283185307 storedOwner=1 storedConsumer=1",
            "DIAG_SEAM ci=1 owner=3 consumer=4 ownerType=cyl consumerType=plane "
            "kind=lin f=0 l=1 storedOwner=1 storedConsumer=1 extra=ignored",
            "DIAG_G7 rid=1 nEdgesIter=4 nEdgesExp=4 nBound=4 nSkippedHasPc=0 "
            "nSkippedNoCurve=0 covered=4 g7=1",
            "DIAG_CYLG7 rid=2 nEdgesIter=4 nEdgesExp=4 nBound=2 nSkippedHasPc=1 "
            "nSkippedNoCurve=1 covered=4 g7=1 nCreated=0",
            "",
        ]
    )


def synthetic_head_stderr() -> str:
    """HEAD @ dc49e45: P2 diag knob is live, seam records are not."""
    return "note: other progress\nRESULT skipped here\n"


def failing_names(checks: Sequence[Check], *, hard_only: bool = True) -> List[str]:
    out: List[str] = []
    for c in checks:
        if c.ok:
            continue
        if hard_only and not c.hard:
            continue
        out.append(c.name)
    return out


def split_hard_fails(
    checks: Sequence[Check],
) -> Tuple[List[Check], List[Check]]:
    hard = [c for c in checks if (not c.ok) and c.hard]
    expected = [c for c in hard if c.name in EXPECTED_RED_NAMES]
    unexpected = [c for c in hard if c not in expected]
    return expected, unexpected


def format_checks(checks: Sequence[Check]) -> str:
    lines = []
    for c in checks:
        if c.ok:
            tag = "PASS"
        elif c.hard:
            tag = "FAIL"
        else:
            tag = "EXPECTED-RED"
        lines.append(f"  [{c.group}/{tag}] {c.message}")
    return "\n".join(lines)


def run_trueform(
    binary: Path,
    stl: Path,
    step: Path,
    *,
    threads: int,
    env: Optional[Mapping[str, str]] = None,
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
    if env:
        environ.update(env)
    environ["STL2STEP_P2_DIAG"] = "1"
    proc = subprocess.run(cmd, capture_output=True, text=True, env=environ)
    result = parse_result_line(proc.stdout)
    return result, proc.stderr, proc.returncode


def run_live(
    binary: Path,
    *,
    pickup_stl: Path = PICKUP_STL,
    hl_stl: Path = HL_STL,
    jobs: int = 2,
    strict: Optional[bool] = None,
) -> int:
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    if not pickup_stl.is_file():
        raise GateError(f"handle-pickup STL missing: {pickup_stl}")
    if not hl_stl.is_file():
        raise GateError(f"handle-lock STL missing: {hl_stl}")
    is_strict = strict_enabled() if strict is None else bool(strict)
    workers = max(1, min(int(jobs) or 1, 2))
    per_file = 0
    arts: Dict[str, Tuple[Dict[str, Any], str, int]] = {}
    with tempfile.TemporaryDirectory(prefix="seam_stored_both_sides_") as td:
        work = Path(td)
        tasks = {
            "pickup": pickup_stl,
            "handle-lock": hl_stl,
        }
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    run_trueform,
                    binary,
                    stl,
                    work / f"{name}.step",
                    threads=per_file,
                ): name
                for name, stl in tasks.items()
            }
            for fut in as_completed(futs):
                arts[futs[fut]] = fut.result()

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
                hard=True,
            )
        )
    else:
        checks.extend(evaluate_pickup(pk_err, strict=is_strict))
    if not hl_res.get("ok"):
        checks.append(
            Check(
                "conversion",
                False,
                f"handle-lock conversion failed: {hl_res.get('error')}",
                "handle-lock",
                hard=True,
            )
        )
    else:
        checks.extend(evaluate_handle_lock(hl_err, strict=is_strict))

    expected_red, unexpected = split_hard_fails(checks)
    print(
        f"seam_stored_both_sides  pickup={pickup_stl}  handle-lock={hl_stl}  "
        f"strict={is_strict}  jobs={workers}  pickup_rc={pk_rc} lock_rc={hl_rc}",
        flush=True,
    )
    print(format_checks(checks), flush=True)
    if unexpected:
        named = ", ".join(c.name for c in unexpected)
        print(f"{MARKER_UNEXPECTED}  unexpected fail(s): {named}", file=sys.stderr, flush=True)
        print(f"seam_stored_both_sides FAIL  unexpected: {named}", file=sys.stderr, flush=True)
        return 1
    if expected_red:
        named = ", ".join(c.name for c in expected_red)
        print(f"{MARKER_EXPECTED_RED}  unmet floors: {named}", flush=True)
        print(
            f"seam_stored_both_sides FAIL  failing assertion(s): {named}",
            file=sys.stderr,
            flush=True,
        )
        return 1
    print(f"{MARKER_FLOORS_MET}  D3-2 stored-both-sides MET", flush=True)
    print(
        "FLIP PROTOCOL: remove PASS_REGULAR_EXPRESSION from "
        "tests/gates/CMakeLists.txt (seam_stored_both_sides). "
        "Keep SEAM_STORED_BOTH_SIDES_STRICT=1.",
        flush=True,
    )
    print("seam_stored_both_sides PASS", flush=True)
    return 0


def run_synthetic_pass() -> int:
    pk = evaluate_pickup(synthetic_pass_stderr(), strict=True)
    hl = evaluate_handle_lock(synthetic_pass_stderr(), strict=True)
    checks = pk + hl
    fails = [c for c in checks if not c.ok]
    print("seam_stored_both_sides --synthetic-pass")
    print(format_checks(checks))
    if fails:
        named = ", ".join(c.name for c in fails)
        print(
            f"seam_stored_both_sides FAIL  synthetic should pass; failed: {named}",
            file=sys.stderr,
        )
        return 1
    print("seam_stored_both_sides PASS  synthetic floors met")
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

    parsed = parse_diag_kv(
        "DIAG_SEAM ci=3 owner=1 consumer=2 ownerType=cyl consumerType=plane "
        "kind=circ f=0.1 l=0.2 storedOwner=1 storedConsumer=0 extra=1"
    )
    check(parsed is not None, "parse_diag_kv DIAG_SEAM")
    assert parsed is not None
    kind, kv = parsed
    check(kind == "DIAG_SEAM", "kind is DIAG_SEAM")
    check(kv.get("extra") == "1", "unknown extra keys tolerated")
    check("storedOwner" in kv and "storedConsumer" in kv, "required stored flags present")
    miss = require_keys({"ci": "1"}, SEAM_REQUIRED, "DIAG_SEAM")
    check(miss is not None and "missing key" in miss, "missing key is not measured")

    head = evaluate_pickup(synthetic_head_stderr(), strict=True)
    named = failing_names(head)
    check("nSeam" in named, "HEAD-red names nSeam (not measured)")
    check(
        any("not measured" in c.message for c in head if c.name == "nSeam"),
        "HEAD-red message says exactly that DIAG_SEAM is not measured",
    )
    check("nCreated" in named, "HEAD-red names nCreated (not measured)")
    check("covered" in named, "HEAD-red names covered (not measured)")
    check("g7" in named, "HEAD-red names g7 (not measured)")
    exp, unexp = split_hard_fails(head)
    check(
        {c.name for c in exp} <= EXPECTED_RED_NAMES,
        "HEAD-red classifies instrument-absent as expected-red",
    )
    check(not unexp, "HEAD-red has no unexpected fails")

    syn = evaluate_pickup(synthetic_pass_stderr(), strict=True) + evaluate_handle_lock(
        synthetic_pass_stderr(), strict=True
    )
    check(all(c.ok for c in syn), "synthetic-pass satisfies every D3-2 / G-7 assert")
    check(failing_names(syn) == [], "synthetic-pass has no failing names")

    bad_store = synthetic_pass_stderr().replace("storedConsumer=1", "storedConsumer=0", 1)
    bs = evaluate_pickup(bad_store, strict=True)
    check("storedConsumer" in failing_names(bs), "storedConsumer=0 is named FAIL")

    created = synthetic_pass_stderr().replace("nCreated=0", "nCreated=2")
    cr = evaluate_pickup(created, strict=True)
    check("nCreated" in failing_names(cr), "nCreated!=0 is named FAIL")

    g7 = synthetic_pass_stderr().replace("g7=1", "g7=0", 1)
    g7c = evaluate_pickup(g7, strict=True)
    check("g7" in failing_names(g7c), "g7=0 is named FAIL")

    cov = synthetic_pass_stderr().replace("covered=4", "covered=3", 1)
    cc = evaluate_pickup(cov, strict=True)
    check("covered" in failing_names(cc), "covered!=nEdgesIter is named FAIL")

    hl_bad = evaluate_handle_lock(
        "DIAG_SEAM ci=9 owner=1 consumer=2 ownerType=plane consumerType=cyl "
        "kind=circ f=0 l=1 storedOwner=1 storedConsumer=0\n",
        strict=True,
    )
    check("handle-lock.stored" in failing_names(hl_bad), "handle-lock 0 stored flag FAILS")
    hl_exp, hl_unexp = split_hard_fails(hl_bad)
    check(
        any(c.name == "handle-lock.stored" for c in hl_unexp)
        and not any(c.name == "handle-lock.stored" for c in hl_exp),
        "handle-lock 0 stored flag is UNEXPECTED (not inverted)",
    )
    bad_type = synthetic_pass_stderr().replace("ownerType=plane", "ownerType=torus", 1)
    bt = evaluate_pickup(bad_type, strict=True)
    check("seamGrammar" in failing_names(bt), "illegal ownerType names seamGrammar")
    check("storedOwner" not in failing_names(bt), "grammar fail is not named storedOwner")
    bt_exp, bt_unexp = split_hard_fails(bt)
    check(
        any(c.name == "seamGrammar" for c in bt_unexp)
        and not any(c.name == "seamGrammar" for c in bt_exp),
        "seamGrammar is UNEXPECTED (not PASS_REGEX inverted)",
    )
    bad_kind = synthetic_pass_stderr().replace("kind=circ", "kind=nurbs", 1)
    bk = evaluate_pickup(bad_kind, strict=True)
    check("seamGrammar" in failing_names(bk), "illegal kind names seamGrammar")
    hl_miss = evaluate_handle_lock("", strict=True)
    check(
        "handle-lock.storedNotMeasured" in failing_names(hl_miss),
        "handle-lock zero DIAG_SEAM is not measured",
    )

    check(strict_enabled({}) is True, "STRICT default ON")
    check(strict_enabled({ENV_STRICT: "0"}) is False, "STRICT=0 is OFF")
    check(
        parse_result_line('progress\nRESULT {"ok":true}\n')["ok"] is True,
        "parse_result_line",
    )
    rc = run_synthetic_pass()
    check(rc == 0, "--synthetic-pass returns 0")
    check(PICKUP_STL.is_file(), f"handle-pickup.stl present ({PICKUP_STL})")
    check(HL_STL.is_file(), "handle-lock.stl present")
    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path, help="stl2step CLI")
    p.add_argument("--pickup-stl", type=Path, default=PICKUP_STL)
    p.add_argument("--handle-lock-stl", type=Path, default=HL_STL)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--synthetic-pass", action="store_true")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--strict", dest="strict", action="store_true", default=None)
    g.add_argument("--no-strict", dest="strict", action="store_false")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.synthetic_pass:
            return run_synthetic_pass()
        if not args.binary:
            print("seam_stored_both_sides: --binary is required", file=sys.stderr)
            return 1
        return run_live(
            Path(args.binary),
            pickup_stl=Path(args.pickup_stl),
            hl_stl=Path(args.handle_lock_stl),
            jobs=args.jobs,
            strict=args.strict,
        )
    except GateError as exc:
        print(f"seam_stored_both_sides FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
