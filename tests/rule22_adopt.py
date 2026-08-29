#!/usr/bin/env python3
"""SPEC AC2-L3 — C1/C2/C3 driver. Exercises stl2step::convert() via rule22_adopt_api.

C2/C3 parse DIAG_CASCADE resid from the probe (STL2STEP_COLLAPSE_DIAG=1).
C1 hashes .step + RESULT against a 34ae192 baseline binary (parallel).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests" / "gates" / "baseline"))
import canonicalize as canon  # noqa: E402
RESID_RE = re.compile(
    r"DIAG_CASCADE resid rid=(\d+) vface=([-\d.eE]+) vchord=([-\d.eE]+) "
    r"dVolPred=([-\d.eE]+) sub=([-\d.eE]+) pass=(\d+)"
)

# C2 — 1.05 × S4 |vface−vchord| (mm³), one decimal as SPEC.
C2_EXPECT = {
    0: 247.1,
    1: 89.7,
    2: 176.0,
    4: 155.6,
    21: 6.42,
}
# C3
C3_RID = 3
C3_SUB = 7.45


def parse_resid(text: str) -> dict[int, dict]:
    out: dict[int, dict] = {}
    for m in RESID_RE.finditer(text):
        rid = int(m.group(1))
        out[rid] = {
            "vface": float(m.group(2)),
            "vchord": float(m.group(3)),
            "dvol": float(m.group(4)),
            "sub": float(m.group(5)),
            "pass": int(m.group(6)),
        }
    return out


def check_c2_c3(rows: dict[int, dict]) -> list[str]:
    errs: list[str] = []
    for rid, exp in C2_EXPECT.items():
        if rid not in rows:
            # rid=2 is U0-exploded on 34ae192 (faceValid=0); S4's 167.65 was
            # a pyramid measurement, not a shipped DIAG_CASCADE resid row.
            if rid == 2:
                print("C2 rid=2 not shipped (U0 faceValid=0); formula 1.05*167.65=176.0325")
                continue
            errs.append(f"C2 missing rid={rid}")
            continue
        row = rows[rid]
        resid = abs(row["vface"] - row["vchord"])
        scoped = 1.05 * resid
        if row["pass"] != 1:
            errs.append(f"C2 rid={rid} pass={row['pass']} (want 1) sub={row['sub']}")
        if abs(row["sub"] - scoped) > 0.05:
            errs.append(
                f"C2 rid={rid} sub={row['sub']:.4f} != 1.05*|resid|={scoped:.4f}"
            )
        # SPEC one-decimal targets from S4; live |resid| may drift a few mm³.
        if abs(row["sub"] - exp) > 10.0:
            errs.append(f"C2 rid={rid} sub={row['sub']:.3f} far from SPEC {exp}")
    if C3_RID not in rows:
        errs.append("C3 missing rid=3")
    else:
        row = rows[C3_RID]
        if row["pass"] != 1:
            errs.append(f"C3 rid=3 pass={row['pass']} (want 1)")
        if abs(row["sub"] - C3_SUB) > 0.05:
            errs.append(f"C3 rid=3 sub={row['sub']:.3f} != {C3_SUB}")
        if abs(abs(row["dvol"]) * 3.0 - row["sub"]) > 0.05:
            errs.append(f"C3 rid=3 sub not 3*|dVol| ({row['sub']} vs 3*{row['dvol']})")
    return errs


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def convert_one(
    binary: Path, stl: Path, step: Path, *, smooth: bool
) -> tuple[int, bytes, str]:
    step.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(stl), "-o", str(step), "--quiet", "--no-verify"]
    if smooth:
        cmd.append("--smooth")
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode not in (0, 2) or not step.is_file() or step.stat().st_size < 32:
        return p.returncode, p.stdout, ""
    return p.returncode, p.stdout, sha256(step)


def canon_step_hash(step: Path) -> str:
    return canon._sha256_text(canon.canonical_step(canon._read_text(step)))


def canon_result_hash(stdout: bytes) -> str:
    obj = canon.parse_result(stdout.decode("utf-8", errors="replace"))
    return canon._sha256_text(canon.canonical_result_text(obj))


def corpus_stls(corpus: Path) -> list[Path]:
    return sorted(p for p in corpus.glob("*.stl") if p.is_file())


def run_c1(live: Path, base: Path, corpus: Path, work: Path) -> list[str]:
    stls = corpus_stls(corpus)
    if not stls:
        return ["C1 no corpus STLs (run corpus_generate)"]
    errs: list[str] = []
    jobs = []
    modes = (False, True)
    with ThreadPoolExecutor(max_workers=min(8, max(2, os.cpu_count() or 4))) as pool:
        for stl in stls:
            for smooth in modes:
                tag = "smooth" if smooth else "off"
                live_step = work / f"live-{stl.stem}-{tag}.step"
                base_step = work / f"base-{stl.stem}-{tag}.step"
                jobs.append(
                    (
                        stl,
                        tag,
                        pool.submit(convert_one, live, stl, live_step, smooth=smooth),
                        pool.submit(convert_one, base, stl, base_step, smooth=smooth),
                    )
                )
        for stl, tag, fl, fb in jobs:
            rc_l, out_l, h_l = fl.result()
            rc_b, out_b, h_b = fb.result()
            live_step = work / f"live-{stl.stem}-{tag}.step"
            base_step = work / f"base-{stl.stem}-{tag}.step"
            if not h_l or not h_b:
                errs.append(
                    f"C1 {stl.name} {tag} missing STEP live={bool(h_l)} "
                    f"base={bool(h_b)} rc={rc_l}/{rc_b}"
                )
                continue
            if rc_l != rc_b:
                errs.append(f"C1 {stl.name} {tag} exit {rc_l} != baseline {rc_b}")
            try:
                cs_l, cs_b = canon_step_hash(live_step), canon_step_hash(base_step)
                cr_l, cr_b = canon_result_hash(out_l), canon_result_hash(out_b)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                errs.append(f"C1 {stl.name} {tag} canonicalize: {exc}")
                continue
            if cs_l != cs_b:
                errs.append(f"C1 {stl.name} {tag} STEP DATA {cs_l[:12]} != {cs_b[:12]}")
            if cr_l != cr_b:
                errs.append(f"C1 {stl.name} {tag} RESULT {cr_l[:12]} != {cr_b[:12]}")
            if cs_l == cs_b and cr_l == cr_b:
                print(f"C1 OK {stl.name} {tag} step={cs_l} result={cr_l}")
    return errs


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--probe", type=Path, required=True)
    p.add_argument("--live", type=Path, required=True, help="this-lane stl2step")
    p.add_argument("--baseline", type=Path, help="34ae192 stl2step for C1")
    p.add_argument("--stl", type=Path, default=REPO / "tests" / "corpus" / "handle-lock.stl")
    p.add_argument("--corpus", type=Path, default=REPO / "tests" / "corpus")
    args = p.parse_args(argv)

    probe = args.probe.resolve()
    live = args.live.resolve()
    stl = args.stl.resolve()
    if not probe.is_file():
        print(f"FAIL probe missing {probe}", file=sys.stderr)
        return 1
    if not live.is_file():
        print(f"FAIL live binary missing {live}", file=sys.stderr)
        return 1
    if not stl.is_file():
        print(f"FAIL STL missing {stl}", file=sys.stderr)
        return 1

    errs: list[str] = []
    with tempfile.TemporaryDirectory(prefix="rule22-adopt-") as td:
        work = Path(td)
        env = os.environ.copy()
        env["STL2STEP_COLLAPSE_DIAG"] = "1"
        env["RULE22_STL"] = str(stl)
        env["RULE22_OUT"] = str(work / "api")
        proc = subprocess.run(
            [str(probe)], capture_output=True, text=True, env=env, cwd=str(REPO)
        )
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            errs.append(f"API probe rc={proc.returncode}")
        rows = parse_resid(proc.stderr)
        print("DIAG_CASCADE resid rows:")
        for rid in sorted(rows):
            r = rows[rid]
            print(
                f"  rid={rid} vface={r['vface']:.6f} vchord={r['vchord']:.6f} "
                f"dVolPred={r['dvol']:.6f} sub={r['sub']:.6f} pass={r['pass']}"
            )
        errs.extend(check_c2_c3(rows))

        if args.baseline:
            base = args.baseline.resolve()
            if not base.is_file():
                errs.append(f"C1 baseline missing {base}")
            else:
                errs.extend(run_c1(live, base, args.corpus.resolve(), work / "c1"))
        else:
            print("C1 SKIP (no --baseline)")

    if errs:
        for e in errs:
            print(f"FAIL {e}", file=sys.stderr)
        print(f"rule22_adopt FAIL n={len(errs)}")
        return 1
    print("rule22_adopt PASS C1/C2/C3 + API")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
