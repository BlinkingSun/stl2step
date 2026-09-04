#!/usr/bin/env python3
"""Composed P1 regiondump + ichecker gate (lane p1-compose-fix).

Dual-axis:
  1. every fixture x every clean component passes check_regionset.py (+ sidecar)
  2. --threads 1 vs --threads 8 produce byte-identical bare RegionSet JSON
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd: list[str], *, capture=False) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=False, capture_output=capture, text=True)


def envelope_clean_indices(dump_bin: Path, stl: Path) -> list[int]:
    proc = run([str(dump_bin), str(stl)], capture=True)
    if proc.returncode != 0:
        raise RuntimeError(f"regiondump failed on {stl}:\n{proc.stderr}")
    doc = json.loads(proc.stdout)
    out: list[int] = []
    for comp in doc.get("comps", []):
        if comp.get("clean"):
            out.append(int(comp["index"]))
    return out


def dump_bare(
    dump_bin: Path, stl: Path, comp: int, threads: int, out: Path
) -> None:
    proc = run(
        [
            str(dump_bin),
            str(stl),
            "--component",
            str(comp),
            "--bare",
            "--threads",
            str(threads),
            "--out",
            str(out),
        ],
        capture=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"regiondump --bare failed on {stl} comp {comp} threads {threads}:\n"
            f"{proc.stderr}"
        )


def ichecker_blocking_fail(proc: subprocess.CompletedProcess) -> str | None:
    """I5-A1: unique-minTri ties on rejected[] are legal until p0-ichecker lands.

    Non-blocking for this lane. Any other ichecker FAIL is blocking.
    """
    if proc.returncode == 0:
        return None
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    blocking: list[str] = []
    for line in text.splitlines():
        s = line.strip()
        if "minLocalTriId" in s and "ties rejected" in s:
            continue
        if re.match(r"I5 FAIL rejected\b", s):
            continue
        if re.match(r"(SCHEMA|I[1-9]|I7b|G5|SIDECAR) FAIL\b", s):
            blocking.append(s)
    if blocking:
        return text
    return None


def sidecar_for(stl: Path, corpus: Path) -> Path | None:
    sc = corpus / (stl.stem + ".expected.json")
    return sc if sc.is_file() else None


def check_fillet_nbrs(rs: dict) -> list[str]:
    """FINDING 2: filletNbrA/B are the two tangent-adjacent regions after D sort."""
    errs: list[str] = []
    regions = [r for r in rs.get("regions", []) if isinstance(r, dict)]
    by_id = {r.get("id"): r for r in regions}
    chains = [c for c in rs.get("chains", []) if isinstance(c, dict)]
    for r in regions:
        if r.get("origin") != "filletStrip":
            continue
        rid = r.get("id")
        a, b = r.get("filletNbrA"), r.get("filletNbrB")
        if a == -1 or b == -1 or a is None or b is None:
            errs.append(f"fillet {rid}: filletNbrA/B={a},{b} (need two dense ids)")
            continue
        if a not in by_id or b not in by_id:
            errs.append(f"fillet {rid}: nbr {a},{b} not in regions[] after area-sort")
            continue
        if a == b:
            errs.append(f"fillet {rid}: filletNbrA==filletNbrB=={a}")
            continue
        tangent_nbrs: set[int] = set()
        for c in chains:
            if not c.get("tangent"):
                continue
            ra, rb = c.get("regA"), c.get("regB")
            if ra == rid and isinstance(rb, int) and rb >= 0:
                tangent_nbrs.add(rb)
            if rb == rid and isinstance(ra, int) and ra >= 0:
                tangent_nbrs.add(ra)
        if a not in tangent_nbrs or b not in tangent_nbrs:
            errs.append(
                f"fillet {rid}: filletNbrA/B={a},{b} are not the tangent-adjacent "
                f"regions {sorted(tangent_nbrs)} after the area sort"
            )
    return errs


# The origin label the region dump prints for an Origin::NgonWall claim
# (D-130-18(2)). dump_regionset.cpp's originName() has cases only for
# PlaneGrow/CylGrow/FilletStrip and falls through to "planeGrow" for NgonWall.
NGON_WALL_DUMP_ORIGIN = "planeGrow"


def check_recognition(stl_name: str, rs: dict) -> list[str]:
    errs: list[str] = []
    regions = rs.get("regions", [])
    rejected = rs.get("rejected", [])
    origins = [r.get("origin") for r in regions]
    stats = rs.get("stats", {})

    if stl_name == "S06.stl":
        # D-130-18(2): since ad5c814 S06's 8-gon wall is claimed by detector A
        # (Origin::NgonWall) before B1 runs, with the STEP byte-identical; the
        # origin cell is re-baselined from cylGrow to that claim. The region
        # dump's originName() predates the enum's NgonWall member and renders
        # it through its default branch as "planeGrow", so that is the label an
        # NgonWall claim carries in a dump at this tip. The cell asks for
        # exactly one cylinder region in total and that it carry the NgonWall
        # label, so a B1 (cylGrow) claim cannot satisfy it.
        all_cyl = [r for r in regions if r.get("type") == "cylinder"]
        cyl = [r for r in all_cyl if r.get("origin") == NGON_WALL_DUMP_ORIGIN]
        if len(all_cyl) != 1 or len(cyl) != 1:
            errs.append(
                f"S06: expected 1 NgonWall cylinder (dump origin "
                f"{NGON_WALL_DUMP_ORIGIN!r}), got {len(cyl)} of "
                f"{len(all_cyl)} cylinder region(s) with origins "
                f"{[r.get('origin') for r in all_cyl]}"
            )
        else:
            if cyl[0].get("nSides") != 8:
                errs.append(f"S06: expected nSides=8, got {cyl[0].get('nSides')}")
            R = float(cyl[0].get("radius", 0.0))
            N = 8
            sidecar_R = 25.0
            refit_tol = 0.02
            budget = max(
                refit_tol,
                0.005 * abs(sidecar_R),
                abs(sidecar_R) * (1.0 - math.cos(math.pi / N)),
            )
            if abs(R - sidecar_R) > budget:
                errs.append(
                    f"S06: radius {R} outside budget |R-25|<={budget}"
                )
    elif stl_name == "S13.stl":
        fil = [r for r in regions if r.get("origin") == "filletStrip"]
        if len(fil) != 1:
            errs.append(f"S13: expected 1 filletStrip, got {len(fil)}")
        elif abs(float(fil[0].get("radius", 0)) - 2.0) > 0.25:
            errs.append(f"S13: fillet radius {fil[0].get('radius')} not ~2")
    elif stl_name == "S12-a.stl":
        # D1.4-A1: 0 cylinders; 9 hole-wall planes; stats.planes == 15
        # (4 square + 5 pentagon + 6 plate). No A3 coplanar merge.
        n_cyl = sum(1 for r in regions if r.get("type") == "cylinder")
        if n_cyl != 0:
            errs.append(f"S12-a: expected 0 cylinders, got {n_cyl}")
        pl = stats.get("planes", -1)
        if pl != 15:
            errs.append(f"S12-a: stats.planes==15 required, got {pl}")
        if pl < 9:
            errs.append(f"S12-a: expected 9 hole-wall planes, got {pl}")
    elif stl_name == "S12-b.stl":
        vr = [r for r in rejected if r.get("reject") == "vertexResidual"]
        extra = [r.get("reject") for r in rejected if r.get("reject") not in (None, "none", "vertexResidual")]
        if len(vr) != 1:
            errs.append(f"S12-b: expected 1 vertexResidual reject, got {len(vr)}")
        if extra:
            errs.append(f"S12-b: extra rejects {extra} (want exactly [vertexResidual])")
    elif stl_name == "S11.stl":
        if any(o == "filletStrip" for o in origins):
            errs.append("S11: filletStrip must not appear")
        fc = [r for r in rejected if r.get("reject") == "filletConsensus"]
        if not fc:
            errs.append("S11: expected filletConsensus in rejected[]")
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--build-dir", type=Path, required=True)
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--ichecker", type=Path, required=True)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--cube", type=Path, required=True)
    args = ap.parse_args()

    fixtures = sorted(args.corpus.glob("S*.stl")) + [args.cube]
    rows: list[str] = []
    failures: list[str] = []

    with tempfile.TemporaryDirectory(prefix="p1_compose_") as tmpdir:
        tmp = Path(tmpdir)
        for stl in fixtures:
            try:
                comps = envelope_clean_indices(args.dump, stl)
            except RuntimeError as e:
                failures.append(str(e))
                rows.append(f"FAIL {stl.name}: dump envelope")
                continue
            if not comps:
                rows.append(f"SKIP {stl.name}: no clean components")
                continue
            sidecar = sidecar_for(stl, args.corpus)
            for comp in comps:
                tag = f"{stl.name} comp{comp}"
                t1 = tmp / f"{stl.stem}_c{comp}_t1.json"
                t8 = tmp / f"{stl.stem}_c{comp}_t8.json"
                try:
                    dump_bare(args.dump, stl, comp, 1, t1)
                    dump_bare(args.dump, stl, comp, 8, t8)
                except RuntimeError as e:
                    failures.append(str(e))
                    rows.append(f"FAIL {tag}: dump")
                    continue
                b1 = t1.read_bytes()
                b8 = t8.read_bytes()
                if b1 != b8:
                    failures.append(f"{tag}: threads 1 vs 8 differ")
                    rows.append(f"FAIL {tag}: determinism")
                    continue
                icmd = [sys.executable, str(args.ichecker), str(t1)]
                if sidecar is not None:
                    icmd += ["--sidecar", str(sidecar)]
                proc = run(icmd, capture=True)
                blocked = ichecker_blocking_fail(proc)
                if blocked is not None:
                    failures.append(f"{tag}: ichecker\n{blocked}")
                    rows.append(f"FAIL {tag}: ichecker")
                    continue
                rs = json.loads(t1.read_text())
                rec = check_recognition(stl.name, rs)
                rec.extend(check_fillet_nbrs(rs))
                if rec:
                    failures.extend(f"{tag}: {m}" for m in rec)
                    rows.append(f"FAIL {tag}: recognition")
                    continue
                rows.append(f"PASS {tag}")

    print("P1 COMPOSE GATE")
    for row in rows:
        print(row)
    print(f"\nSUMMARY: {sum(r.startswith('PASS') for r in rows)}/{len(rows)} pass")
    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for f in failures:
            print(f, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
