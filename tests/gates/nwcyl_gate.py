#!/usr/bin/env python3
"""nwcyl_gate — NW-CYL never-worse (D-S3-15).

Per part with a GT cylinder multiset in its sidecar, assert shipped analytic
cylinders do not regress vs tests/gates/nwcyl_baseline.json:

  nGTMatched(now) >= baseline.nGTMatched
  nPhantom(now) <= baseline.nPhantom

GT matching reuses partial_recovery_gate.py (match_gt_radius @ 0.3 %,
census_radii / step_census — do not fork rules).

RATChet procedure (landing lane updates baseline when a part improves):
  1. Run gate; note per-part nGTMatched / nPhantom in the PASS line.
  2. If nGTMatched rose or nPhantom fell vs baseline, edit nwcyl_baseline.json:
       - set parts[<id>].nGTMatched = measured nGTMatched
       - set parts[<id>].nPhantom = measured nPhantom
       - bump engineRef to the landing commit
  3. Commit baseline + gate in the same lane; never lower nGTMatched or raise
     nPhantom in baseline without a measured regression fix elsewhere.

Usage:
  nwcyl_gate.py --self-test
  nwcyl_gate.py --measure --binary ./build/stl2step
  nwcyl_gate.py --binary ./build/stl2step [--baseline PATH] [--jobs N]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
CORPUS = REPO / "tests" / "corpus"
DEFAULT_BASELINE = REPO / "tests" / "gates" / "nwcyl_baseline.json"
BASELINE_NAME = "nwcyl_baseline.json"
HP_GT = REPO / "tests" / "diag" / "handle-pickup" / "ground-truth.json"

# Parts named in D-S3-15 / DECISION follow-up; only those with GT radii are gated.
CANDIDATE_PARTS: Tuple[str, ...] = (
    "handle-pickup",
    "handle-lock",
    "S02",
    "S04",
    "S05",
    "S09",
    "S11-b",
    "S16-R1-explode-success",
    "Body11",
)


class GateError(Exception):
    """Hard gate failure."""


def _load_prg():
    prg_path = REPO / "tests" / "gates" / "partial_recovery_gate.py"
    if not prg_path.is_file():
        raise GateError(f"partial_recovery_gate.py missing: {prg_path}")
    spec = importlib.util.spec_from_file_location("partial_recovery_gate", prg_path)
    if spec is None or spec.loader is None:
        raise GateError(f"cannot load {prg_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["partial_recovery_gate"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_json(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, doc: Dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


def load_gt_radii_all(part_id: str, sidecar: Mapping[str, Any]) -> Tuple[List[float], str]:
    """GT cylinder multiset from sidecar (or HP diag ground-truth for multiplicity)."""
    if part_id == "handle-pickup" and HP_GT.is_file():
        doc = json.loads(HP_GT.read_text(encoding="utf-8"))
        raw_all = doc.get("cylinder_radii_all")
        if isinstance(raw_all, list) and raw_all:
            return [float(x) for x in raw_all], "diag/handle-pickup/ground-truth.json"

    raw = sidecar.get("cylinder_radii_all")
    if isinstance(raw, list) and raw:
        return [float(x) for x in raw], "sidecar.cylinder_radii_all"

    for row in sidecar.get("live") or []:
        if not isinstance(row, dict):
            continue
        row_all = row.get("cylinder_radii_all")
        if isinstance(row_all, list) and row_all:
            return [float(x) for x in row_all], "sidecar.live.cylinder_radii_all"
        row_r = row.get("cylinder_radii")
        if isinstance(row_r, list) and row_r:
            return [float(x) for x in row_r], "sidecar.live.cylinder_radii"

    out: List[float] = []
    for item in sidecar.get("recoverable") or []:
        if not isinstance(item, dict) or item.get("type") != "cylinder":
            continue
        try:
            radius = float(item["radius"])
            count = int(item.get("count") or 1)
        except (KeyError, TypeError, ValueError):
            continue
        out.extend([radius] * max(1, count))
    if out:
        return out, "sidecar.recoverable"

    return [], "none"


def discover_gated_parts(corpus: Path) -> List[str]:
    gated: List[str] = []
    for pid in CANDIDATE_PARTS:
        sc_path = corpus / f"{pid}.expected.json"
        if not sc_path.is_file():
            continue
        sidecar = load_json(sc_path)
        gt_all, _ = load_gt_radii_all(pid, sidecar)
        if gt_all:
            gated.append(pid)
    return gated


def compute_nwcyl(
    built: Sequence[float],
    gt_all: Sequence[float],
    *,
    match_gt_radius,
) -> Tuple[int, int, List[Dict[str, Any]]]:
    """Return (nGTMatched, nPhantom, per-radius rows)."""
    gt_counts: Counter = Counter(float(x) for x in gt_all)
    gt_keys = list(gt_counts.keys())
    built_match: Counter = Counter()
    n_phantom = 0
    for r in built:
        g = match_gt_radius(float(r), gt_keys)
        if g is None:
            n_phantom += 1
        else:
            built_match[g] += 1
    n_matched = sum(built_match.values())
    per_radius: List[Dict[str, Any]] = []
    for g, mult in sorted(gt_counts.items()):
        per_radius.append(
            {
                "radius": g,
                "built": int(built_match.get(g, 0)),
                "gtMultiplicity": int(mult),
            }
        )
    return n_matched, n_phantom, per_radius


def format_radius_table(per_radius: Sequence[Mapping[str, Any]]) -> str:
    bits = [
        f"R={row['radius']} {row['built']}/{row['gtMultiplicity']}"
        for row in per_radius
    ]
    return " ".join(bits)


def convert_part(
    prg: Any,
    *,
    binary: Path,
    stl: Path,
    work: Path,
    verify: bool,
) -> Dict[str, Any]:
    work.mkdir(parents=True, exist_ok=True)
    step = work / f"{stl.stem}.step"
    art = prg.run_trueform(binary, stl, step, threads=0, verify=verify)
    if art.step is None or not art.step.is_file():
        raise GateError(
            f"{stl.stem}: TrueForm wrote no STEP: {art.result.get('error') or art.result}"
        )
    census = prg.run_step_census(art.step)
    built = prg.census_radii(census)
    return {
        "built": built,
        "census": census,
        "result": art.result,
        "exit_code": art.exit_code,
    }


def measure_part(
    prg: Any,
    *,
    part_id: str,
    corpus: Path,
    binary: Path,
    work: Path,
) -> Dict[str, Any]:
    stl = corpus / f"{part_id}.stl"
    sidecar_path = corpus / f"{part_id}.expected.json"
    if not stl.is_file():
        raise GateError(f"{part_id}: STL missing: {stl}")
    if not sidecar_path.is_file():
        raise GateError(f"{part_id}: sidecar missing: {sidecar_path}")
    sidecar = load_json(sidecar_path)
    gt_all, gt_src = load_gt_radii_all(part_id, sidecar)
    if not gt_all:
        raise GateError(f"{part_id}: no GT cylinder radii in sidecar")
    verify = part_id not in ("handle-lock", "Body11")
    conv = convert_part(prg, binary=binary, stl=stl, work=work, verify=verify)
    n_matched, n_phantom, per_radius = compute_nwcyl(
        conv["built"], gt_all, match_gt_radius=prg.match_gt_radius
    )
    return {
        "nGTMatched": n_matched,
        "nPhantom": n_phantom,
        "nGtTotal": len(gt_all),
        "gtSource": gt_src,
        "perRadius": per_radius,
        "nBuilt": len(conv["built"]),
    }


def load_baseline(path: Path) -> Dict[str, Any]:
    if path.name != BASELINE_NAME:
        raise GateError(f"baseline authority must be {BASELINE_NAME}, opened {path}")
    doc = load_json(path)
    if str(doc.get("authority") or "") != BASELINE_NAME:
        raise GateError(f"{path}: authority must be {BASELINE_NAME!r}")
    if "parts" not in doc or not isinstance(doc["parts"], dict):
        raise GateError(f"{path}: missing parts map")
    return doc


def run_measure(args: argparse.Namespace) -> int:
    prg = _load_prg()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    corpus = args.corpus.resolve()
    parts = discover_gated_parts(corpus)
    if not parts:
        raise GateError(f"no gated parts under {corpus}")
    doc: Dict[str, Any] = {
        "id": "nwcyl-baseline",
        "authority": BASELINE_NAME,
        "engineRef": args.engine_ref,
        "notes": (
            "NW-CYL per-part baselines (D-S3-15). Ratchet up when landing improves "
            "nGTMatched or lowers nPhantom — see nwcyl_gate.py header."
        ),
        "parts": {},
    }
    with tempfile.TemporaryDirectory(prefix="nwcyl_measure_") as td:
        work = Path(td)
        for pid in parts:
            row = measure_part(
                prg, part_id=pid, corpus=corpus, binary=binary, work=work
            )
            doc["parts"][pid] = {
                "nGTMatched": row["nGTMatched"],
                "nPhantom": row["nPhantom"],
                "nGtTotal": row["nGtTotal"],
                "gtSource": row["gtSource"],
            }
            row["_pid"] = pid
            print(
                f"  {pid}: nGTMatched={row['nGTMatched']} nPhantom={row['nPhantom']} "
                f"gtTotal={row['nGtTotal']}  {format_radius_table(row['perRadius'])}"
            )
    out = args.baseline.resolve()
    write_json(out, doc)
    print(f"nwcyl_gate --measure wrote {out} ({len(parts)} parts)")
    return 0


def run_gate(args: argparse.Namespace) -> int:
    prg = _load_prg()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise GateError(f"stl2step binary missing: {binary}")
    baseline = load_baseline(args.baseline.resolve())
    corpus = args.corpus.resolve()
    part_ids = sorted(baseline["parts"].keys())
    failures: List[str] = []
    lines: List[str] = []

    with tempfile.TemporaryDirectory(prefix="nwcyl_gate_") as td:
        work = Path(td)
        workers = max(1, min(int(args.jobs) or 1, len(part_ids)))
        measured: Dict[str, Dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    measure_part,
                    prg,
                    part_id=pid,
                    corpus=corpus,
                    binary=binary,
                    work=work / pid,
                ): pid
                for pid in part_ids
            }
            for fut in as_completed(futs):
                pid = futs[fut]
                measured[pid] = fut.result()

    for pid in part_ids:
        base = baseline["parts"][pid]
        now = measured[pid]
        b_match = int(base["nGTMatched"])
        b_phantom = int(base["nPhantom"])
        n_match = int(now["nGTMatched"])
        n_phantom = int(now["nPhantom"])
        gt_total = int(now["nGtTotal"])
        ok_match = n_match >= b_match
        ok_phantom = n_phantom <= b_phantom
        status = "PASS" if ok_match and ok_phantom else "FAIL"
        line = (
            f"  {pid}: nGTMatched={n_match}/{b_match} nPhantom={n_phantom}/{b_phantom} "
            f"gtTotal={gt_total} [{status}]"
        )
        lines.append(line)
        print(line)
        print(f"    radii: {format_radius_table(now['perRadius'])}")
        if not ok_match:
            failures.append(
                f"{pid}: nGTMatched {n_match} < baseline {b_match}"
            )
        if not ok_phantom:
            failures.append(
                f"{pid}: nPhantom {n_phantom} > baseline {b_phantom}"
            )

    if failures:
        print("nwcyl_gate FAIL", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print(
        f"nwcyl_gate PASS  baseline={args.baseline} engineRef={baseline.get('engineRef')}"
    )
    return 0


def _self_test() -> int:
    prg = _load_prg()
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    gt_all = [2.0] * 12 + [5.0, 5.0]
    built_good = [2.0] * 10 + [7.0, 7.0, 7.0]
    n_match, n_phantom, per = compute_nwcyl(
        built_good, gt_all, match_gt_radius=prg.match_gt_radius
    )
    check(n_match == 10, f"synthetic: nGTMatched=10 (got {n_match})")
    check(n_phantom == 3, f"synthetic: nPhantom=3 (got {n_phantom})")
    check(per[0]["built"] == 10 and per[0]["gtMultiplicity"] == 12, "R=2 built 10/12")

    # GT-matched decrease must trip RED vs baseline 10/3.
    baseline_row = {"nGTMatched": 10, "nPhantom": 3}
    built_regress = [2.0] * 8 + [7.0]
    n_reg, p_reg, _ = compute_nwcyl(
        built_regress, gt_all, match_gt_radius=prg.match_gt_radius
    )
    check(n_reg < baseline_row["nGTMatched"], "matched decrease is a regression")
    check(n_phantom <= baseline_row["nPhantom"], "phantom decrease is allowed")

    built_phantom_down = [2.0] * 10
    n2, p2, _ = compute_nwcyl(
        built_phantom_down, gt_all, match_gt_radius=prg.match_gt_radius
    )
    check(n2 == 10 and p2 == 0, "phantom drop 3→0 with same matched passes NW-CYL")
    check(
        n2 >= baseline_row["nGTMatched"] and p2 <= baseline_row["nPhantom"],
        "phantom-only improvement satisfies NW-CYL",
    )

    built_regress2 = [2.0] * 8 + [7.0]
    n3, p3, _ = compute_nwcyl(
        built_regress2, gt_all, match_gt_radius=prg.match_gt_radius
    )
    check(n3 == 8, "matched regression 10→8")
    check(
        not (n3 >= baseline_row["nGTMatched"] and p3 <= baseline_row["nPhantom"]),
        "matched decrease fails NW-CYL inequality",
    )

    sc = load_json(CORPUS / "handle-lock.expected.json")
    gt_hl, src_hl = load_gt_radii_all("handle-lock", sc)
    check(len(gt_hl) == 15, f"handle-lock GT multiset len=15 (got {len(gt_hl)})")
    check("cylinder_radii" in src_hl, "handle-lock GT from live cylinder_radii")

    sc_hp = load_json(CORPUS / "handle-pickup.expected.json")
    gt_hp, src_hp = load_gt_radii_all("handle-pickup", sc_hp)
    check(len(gt_hp) == 54, f"handle-pickup GT multiset len=54 (got {len(gt_hp)})")
    check("ground-truth" in src_hp, "HP multiplicity from diag ground-truth")

    try:
        load_baseline(Path("/tmp/not-nwcyl.json"))
        check(False, "wrong baseline name must FAIL")
    except GateError:
        check(True, "wrong baseline name FAIL")

    return 1 if fails else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--binary", type=Path, help="stl2step CLI")
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    p.add_argument("--corpus", type=Path, default=CORPUS)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument(
        "--engine-ref",
        default="ab6a7ee667d78d95e62022fb474f6191741c3540",
        help="commit recorded in baseline (measure mode)",
    )
    p.add_argument(
        "--measure",
        action="store_true",
        help="write baseline JSON from a fresh build (orchestrator only)",
    )
    p.add_argument("--self-test", action="store_true", help="API tests (no engine)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.measure:
            if args.binary is None:
                raise GateError("--binary is required for --measure")
            return run_measure(args)
        if args.binary is None:
            raise GateError("--binary is required")
        return run_gate(args)
    except GateError as exc:
        print(f"nwcyl_gate FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
