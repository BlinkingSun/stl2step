#!/usr/bin/env python3
"""hl_census_ratchet — handle-lock census-growth + GT-radii machine gate.

Floor authority is tests/gates/baseline/hl-ratchet.json (start floor=1).
tests/corpus/handle-lock.expected.json live[].builtCylindersFloor stays 0.

Usage:
  hl_census_ratchet.py --binary ./build/stl2step --census ./build/stl2step_census
  hl_census_ratchet.py --ratchet-up          # orchestrator: floor += 1
  hl_census_ratchet.py --ratchet-up 2        # orchestrator: set floor
  hl_census_ratchet.py --self-test           # API unit tests (no engine)
  hl_census_ratchet.py --waiver-audit --binary ... --census ... --dump ...
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
DEFAULT_STL = REPO / "tests" / "corpus" / "handle-lock.stl"
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "hl-ratchet.json"
DEFAULT_FIXTURE = REPO / "tests" / "corpus" / "handle-lock.expected.json"
DEFAULT_GT = REPO / "tests" / "diag" / "handle-lock" / "ground-truth.json"
DEFAULT_WAIVER = REPO / "tests" / "diag" / "handle-lock" / "waiver_audit.py"
RATCHET_NAME = "hl-ratchet.json"
GT_REL_TOL = 0.003  # 0.3%
NAMED_RADII: Tuple[Tuple[str, float], ...] = (("R20", 20.0), ("R30", 30.0))


class GateError(Exception):
    """Hard gate failure."""


def load_json(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, doc: Dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


def fixture_built_floor(fixture: Dict[str, Any]) -> int:
    live = fixture.get("live") or []
    if not live:
        return 0
    return max(int(row.get("builtCylindersFloor") or 0) for row in live)


def fixture_volume_budget(fixture: Dict[str, Any]) -> float:
    return float(fixture.get("smoothVolumeDeltaPctMax", 6.013))


def assert_floor_authority(
    ratchet: Dict[str, Any], fixture: Dict[str, Any], ratchet_path: Path
) -> int:
    """Prove the moving floor came from hl-ratchet.json, not the sidecar."""
    if ratchet_path.name != RATCHET_NAME:
        raise GateError(
            f"floor authority must be {RATCHET_NAME}, opened {ratchet_path}"
        )
    if "floor" not in ratchet:
        raise GateError(f"{ratchet_path}: missing 'floor'")
    if str(ratchet.get("authority") or "") != RATCHET_NAME:
        raise GateError(
            f"{ratchet_path}: authority must be {RATCHET_NAME!r} "
            f"(got {ratchet.get('authority')!r}) — refusing fixture floor"
        )
    floor = int(ratchet["floor"])
    fx_floor = fixture_built_floor(fixture)
    if fx_floor != 0:
        raise GateError(
            f"fixture builtCylindersFloor={fx_floor} must stay frozen at 0; "
            f"authority is {ratchet_path} floor={floor}"
        )
    if floor == fx_floor:
        raise GateError(
            f"ratchet floor {floor} equals fixture builtCylindersFloor "
            f"{fx_floor} — cannot prove ratchet (not fixture) was read"
        )
    return floor


def ratchet_up(ratchet_path: Path, target: Optional[int] = None) -> int:
    """Orchestrator-only: increment floor, or set it to target."""
    doc = load_json(ratchet_path)
    old = int(doc["floor"])
    doc["floor"] = old + 1 if target is None else int(target)
    write_json(ratchet_path, doc)
    print(f"ratchet-up: floor {old} -> {doc['floor']} ({ratchet_path})")
    return int(doc["floor"])


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise GateError("no RESULT line in stl2step stdout")


def census_cylinder_count(cen: Dict[str, Any]) -> int:
    surfaces = cen.get("surfaces")
    if isinstance(surfaces, dict) and "cylinder" in surfaces:
        return int(surfaces.get("cylinder") or 0)
    return int(cen.get("cylinder") or 0)


def census_radii(cen: Dict[str, Any]) -> List[float]:
    out: List[float] = []
    for face in cen.get("cylinders") or []:
        try:
            out.append(float(face.get("radius")))
        except (TypeError, ValueError):
            continue
    if not out:
        for grp in cen.get("cylinderGroups") or []:
            try:
                out.append(float(grp.get("radius")))
            except (TypeError, ValueError):
                continue
    return out


def match_gt_radius(radius: float, gt_radii: Sequence[float]) -> Optional[float]:
    """Return the nearest GT radius within 0.3%, or None."""
    best: Optional[float] = None
    best_err: Optional[float] = None
    for g in gt_radii:
        denom = abs(g)
        if denom == 0.0:
            if abs(radius) == 0.0:
                return 0.0
            continue
        err = abs(radius - g) / denom
        if err <= GT_REL_TOL and (best_err is None or err < best_err):
            best, best_err = g, err
    return best


def named_radius_report(
    census_rs: Sequence[float], gt_radii: Sequence[float]
) -> Dict[str, Dict[str, Any]]:
    """R20/R30 named explicitly — present only when a census face matches GT."""
    report: Dict[str, Dict[str, Any]] = {}
    for name, nominal in NAMED_RADII:
        gt = match_gt_radius(nominal, gt_radii)
        want = gt if gt is not None else nominal
        hits = [r for r in census_rs if match_gt_radius(r, [want]) is not None]
        report[name] = {
            "nominal": nominal,
            "gt": want,
            "present": bool(hits),
            "censusRadii": hits,
            "tol": GT_REL_TOL,
        }
    return report


def check_gt_radii(
    census_rs: Sequence[float], gt_radii: Sequence[float]
) -> Tuple[List[str], Dict[str, Dict[str, Any]]]:
    """Every census cylinder radius must match a GT entry within 0.3%."""
    violations: List[str] = []
    for r in census_rs:
        hit = match_gt_radius(r, gt_radii)
        if hit is None:
            violations.append(
                f"census R={r} matches no ground-truth.json radius "
                f"within {GT_REL_TOL * 100:.1f}% (anti-phantom)"
            )
    named = named_radius_report(census_rs, gt_radii)
    return violations, named


def check_file_truth(result_built: int, census_cyl: int) -> Optional[str]:
    if result_built != census_cyl:
        return (
            f"file-truth: census cylinders={census_cyl} != "
            f"RESULT smoothBuiltCylinders={result_built}"
        )
    return None


def run_trueform(
    binary: Path, stl: Path, step: Path, threads: int
) -> Dict[str, Any]:
    import subprocess

    cmd = [
        str(binary),
        str(stl),
        "-o",
        str(step),
        "--engine",
        "trueform",
        "--quiet",
        "--threads",
        str(threads),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        result = parse_result_line(proc.stdout)
    except GateError as exc:
        raise GateError(
            f"{exc}\nstdout={proc.stdout!r}\nstderr={proc.stderr!r}"
        ) from exc
    result["_exit"] = proc.returncode
    result["_stderr"] = proc.stderr
    return result


def run_census(census_bin: Path, step: Path) -> Dict[str, Any]:
    import subprocess

    proc = subprocess.run(
        [str(census_bin), str(step)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise GateError(
            f"census failed rc={proc.returncode}: {proc.stderr or proc.stdout}"
        )
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise GateError(f"census output is not JSON: {exc}") from exc


def run_gate(args: argparse.Namespace) -> int:
    ratchet_path = args.ratchet.resolve()
    fixture_path = args.fixture.resolve()
    gt_path = args.ground_truth.resolve()
    stl = args.stl.resolve()
    binary = args.binary.resolve() if args.binary else None
    census_bin = args.census.resolve() if args.census else None

    ratchet = load_json(ratchet_path)
    fixture = load_json(fixture_path)
    gt = load_json(gt_path)
    floor = assert_floor_authority(ratchet, fixture, ratchet_path)
    vol_budget = fixture_volume_budget(fixture)
    gt_radii = [float(x) for x in (gt.get("cylinder_radii") or [])]

    if not stl.is_file():
        raise GateError(f"handle-lock STL missing: {stl}")
    if binary is None or not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    if census_bin is None or not census_bin.is_file():
        raise GateError(f"stl2step_census binary missing: {census_bin}")

    with tempfile.TemporaryDirectory(prefix="hl_census_ratchet_") as td:
        step = Path(td) / "handle-lock.step"
        result = run_trueform(binary, stl, step, args.threads)
        if not step.is_file():
            raise GateError(
                f"TrueForm wrote no STEP: {result.get('error') or result}"
            )
        # Real compute after convert: census is independent I/O; keep the
        # worker pool so a future second witness can join without serialising.
        with ThreadPoolExecutor(max_workers=2) as pool:
            cen_f = pool.submit(run_census, census_bin, step)
            cen = cen_f.result()

    failures: List[str] = []
    details: Dict[str, Any] = {
        "floorSource": str(ratchet_path),
        "floor": floor,
        "fixtureFloor": fixture_built_floor(fixture),
        "usedFixtureFloor": False,
        "authority": ratchet.get("authority"),
    }

    if not result.get("ok"):
        failures.append(f"conversion failed: {result.get('error')}")

    built = int(result.get("smoothBuiltCylinders") or 0)
    details["smoothBuiltCylinders"] = built
    if built < floor:
        failures.append(f"smoothBuiltCylinders={built} < ratchet floor={floor}")

    census_cyl = census_cylinder_count(cen)
    details["censusCylinders"] = census_cyl
    ft = check_file_truth(built, census_cyl)
    if ft:
        failures.append(ft)

    vd = float(result.get("volumeDeltaPct", -1))
    details["volumeDeltaPct"] = vd
    details["volumeBudgetPct"] = vol_budget
    if vd < 0 or vd > vol_budget + 1e-9:
        failures.append(
            f"volumeDeltaPct={vd} outside fixture budget {vol_budget}"
        )

    if not result.get("watertight", False):
        failures.append("RESULT watertight=false")
    if int(result.get("openShells") or 0) != 0:
        failures.append(f"RESULT openShells={result.get('openShells')}")

    if not cen.get("valid", False):
        failures.append("census BRepCheck valid=false")
    if not cen.get("closed", False):
        failures.append("census closed=false")
    if int(cen.get("openShells") or 0) != 0:
        failures.append(f"census openShells={cen.get('openShells')}")

    warnings = list(result.get("warnings") or [])
    details["warnings"] = warnings
    if warnings:
        failures.append(f"warnings not empty: {warnings}")

    radii = census_radii(cen)
    details["censusRadii"] = radii
    gt_violations, named = check_gt_radii(radii, gt_radii)
    details["namedRadii"] = named
    failures.extend(gt_violations)

    named_bits = " ".join(
        f"{n}={'present' if named[n]['present'] else 'absent'}"
        for n, _ in NAMED_RADII
    )

    if failures:
        print("hl_census_ratchet FAIL", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        print(json.dumps(details, indent=2, sort_keys=True))
        return 1

    print(
        f"hl_census_ratchet PASS  floor={floor} from {ratchet_path} "
        f"(not fixture builtCylindersFloor={details['fixtureFloor']})  "
        f"built={built} census={census_cyl} volΔ={vd}  {named_bits}"
    )
    print(json.dumps(details, indent=2, sort_keys=True))
    return 0


def run_waiver_audit(args: argparse.Namespace) -> int:
    waiver = Path(args.waiver).resolve() if args.waiver else DEFAULT_WAIVER
    if not waiver.is_file():
        print(f"waiver_audit missing: {waiver}", file=sys.stderr)
        return 1
    spec = importlib.util.spec_from_file_location("waiver_audit", waiver)
    if spec is None or spec.loader is None:
        print(f"cannot load {waiver}", file=sys.stderr)
        return 1
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if args.binary:
        mod.STL2STEP = Path(args.binary).resolve()
    if args.census:
        mod.CENSUS = Path(args.census).resolve()
    if args.dump:
        mod.REGIONDUMP = Path(args.dump).resolve()
    return int(mod.main())


def _self_test() -> int:
    """Exercise the gate API without the engine."""
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    ratchet = load_json(DEFAULT_RATCHET)
    fixture = load_json(DEFAULT_FIXTURE)
    gt = load_json(DEFAULT_GT)
    floor = assert_floor_authority(ratchet, fixture, DEFAULT_RATCHET)
    check(floor == 1, f"start floor is 1 (got {floor})")
    check(
        DEFAULT_RATCHET.name == RATCHET_NAME,
        "floor source basename is hl-ratchet.json",
    )
    check(
        fixture_built_floor(fixture) == 0,
        "fixture builtCylindersFloor frozen at 0",
    )
    check(floor != fixture_built_floor(fixture), "ratchet floor != fixture floor")

    gt_radii = [float(x) for x in gt["cylinder_radii"]]
    check(match_gt_radius(5.75, gt_radii) == 5.75, "R5.75 matches GT")
    check(match_gt_radius(20.0, gt_radii) == 20.0, "R20 exact matches GT")
    check(match_gt_radius(30.0, gt_radii) == 30.0, "R30 exact matches GT")
    # 0.3% of 20 = 0.06; 20.06 ok, 20.07 not
    check(match_gt_radius(20.06, gt_radii) == 20.0, "R20 at +0.3% still matches")
    check(match_gt_radius(20.07, gt_radii) is None, "R20 at +0.35% is phantom")
    check(match_gt_radius(7.0, gt_radii) is None, "R7 is not a GT radius")

    viol, named = check_gt_radii([5.75], gt_radii)
    check(viol == [], "today's built R5.75 is GT-legal")
    check(named["R20"]["present"] is False, "R20 named absent at floor=1")
    check(named["R30"]["present"] is False, "R30 named absent at floor=1")

    viol20, named20 = check_gt_radii([5.75, 20.0, 30.0], gt_radii)
    check(viol20 == [], "R20+R30 exact faces are GT-legal")
    check(named20["R20"]["present"] is True, "R20 named present when built")
    check(named20["R30"]["present"] is True, "R30 named present when built")

    phantom, _ = check_gt_radii([5.75, 7.0], gt_radii)
    check(len(phantom) == 1, "phantom radius is a GT-radii FAIL")

    check(check_file_truth(1, 1) is None, "file-truth 1==1")
    check(check_file_truth(16, 1) is not None, "file-truth 16!=1 is FAIL")

    # Authority refuses a mis-named floor file.
    try:
        assert_floor_authority(ratchet, fixture, Path("/tmp/not-the-ratchet.json"))
        check(False, "mis-named floor file must FAIL")
    except GateError:
        check(True, "mis-named floor file FAIL")

    spoof = dict(ratchet)
    spoof["authority"] = "handle-lock.expected.json"
    try:
        assert_floor_authority(spoof, fixture, DEFAULT_RATCHET)
        check(False, "fixture-authority spoof must FAIL")
    except GateError:
        check(True, "fixture-authority spoof FAIL")

    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, help="stl2step CLI")
    p.add_argument("--census", type=Path, help="stl2step_census witness")
    p.add_argument("--dump", type=Path, help="stl2step_regiondump (waiver_audit)")
    p.add_argument("--stl", type=Path, default=DEFAULT_STL)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    p.add_argument("--ground-truth", type=Path, default=DEFAULT_GT)
    p.add_argument("--waiver", type=Path, default=DEFAULT_WAIVER)
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="stl2step --threads (0 = all cores)",
    )
    p.add_argument(
        "--ratchet-up",
        nargs="?",
        const="__inc__",
        metavar="N",
        help="bump floor by 1, or set to N (orchestrator only)",
    )
    p.add_argument(
        "--self-test",
        action="store_true",
        help="exercise gate API (no engine)",
    )
    p.add_argument(
        "--waiver-audit",
        action="store_true",
        help="run tests/diag/handle-lock/waiver_audit.py with injected bins",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.ratchet_up is not None:
            target = None if args.ratchet_up == "__inc__" else int(args.ratchet_up)
            ratchet_up(args.ratchet.resolve(), target)
            return 0
        if args.waiver_audit:
            return run_waiver_audit(args)
        return run_gate(args)
    except GateError as exc:
        print(f"hl_census_ratchet FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
