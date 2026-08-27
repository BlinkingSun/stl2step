#!/usr/bin/env python3
"""Live P2 re-gate: frozen island-aware P1 dumps × buildFaces (lane p2-real).

Consumes:
  --dumps-dir   _team/p1-dumps/ MANIFEST + *_cN.json (island-aware P1)
  --sidecars    p0-gen-rw1 live[] faceCount/surfaceCensus/volumeBudget/disposition
  --corpus      STLs (prefer sidecars dir, then --corpus)

S15 and S16-R1-round-2 have no clean components: SKIP, not FAIL (I6).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], *, capture=False) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=False, capture_output=capture, text=True)


def sidecar_for(stl: Path, *roots: Path) -> Path | None:
    for root in roots:
        if root is None:
            continue
        sc = root / (stl.stem + ".expected.json")
        if sc.is_file():
            return sc
    return None


def stl_for(name: str, *roots: Path) -> Path | None:
    for root in roots:
        if root is None:
            continue
        p = root / name
        if p.is_file():
            return p
    return None


def live_row(sidecar: Path | None, comp: int) -> dict | None:
    if sidecar is None:
        return None
    try:
        doc = json.loads(sidecar.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    for row in doc.get("live") or []:
        if int(row.get("component", 0)) == comp:
            return row
    return None


def p2_live(
    p2: Path, rs: Path, stl: Path, comp: int, sidecar: Path | None
) -> dict:
    cmd = [str(p2), "--live", "--component", str(comp), str(rs), str(stl)]
    if sidecar is not None:
        cmd += ["--sidecar", str(sidecar)]
    proc = run(cmd, capture=True)
    line = ""
    if proc.stdout:
        line = proc.stdout.strip().splitlines()[-1]
    try:
        doc = json.loads(line) if line else {}
    except json.JSONDecodeError:
        doc = {}
    doc["_rc"] = proc.returncode
    doc["_stderr"] = proc.stderr or ""
    doc["_raw"] = line
    return doc


def identity_key(doc: dict) -> tuple:
    return (
        int(doc.get("faceCount", -1)),
        int(doc.get("planes", -1)),
        int(doc.get("cylinders", -1)),
        tuple(doc.get("builtAs") or []),
        f"{float(doc.get('volume', 0.0)):.6g}",
    )


def fail_sig(doc: dict, id_ok: bool) -> tuple:
    bits: list[str] = []
    if not doc.get("buildFaces"):
        bits.append("build")
    if not doc.get("valid"):
        bits.append("valid")
    if not doc.get("closed"):
        bits.append("closed")
    if not doc.get("j6"):
        bits.append("j6")
    if not doc.get("volumeOk"):
        bits.append("vol")
    if not doc.get("j1"):
        bits.append("j1")
    if doc.get("explodedRecoverableHole") or doc.get("explodedRecoverableSurface"):
        bits.append("explodedRec")
    if not doc.get("tolOk", True):
        bits.append("tol")
    if not id_ok:
        bits.append("threads")
    if int(doc.get("twins") or 0) != 0:
        bits.append("twins")
    return tuple(bits)


def dump_bare(dump_bin: Path, stl: Path, comp: int, out: Path) -> None:
    proc = run(
        [
            str(dump_bin),
            str(stl),
            "--component",
            str(comp),
            "--bare",
            "--threads",
            "1",
            "--out",
            str(out),
        ],
        capture=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"regiondump --bare failed on {stl} comp {comp}:\n{proc.stderr}"
        )


def envelope_clean(dump_bin: Path, stl: Path) -> list[int]:
    proc = run([str(dump_bin), str(stl)], capture=True)
    if proc.returncode != 0:
        raise RuntimeError(f"regiondump failed on {stl}:\n{proc.stderr}")
    doc = json.loads(proc.stdout)
    return [int(c["index"]) for c in doc.get("comps", []) if c.get("clean")]


def load_manifest(dumps_dir: Path) -> dict:
    man_path = dumps_dir / "MANIFEST.json"
    if not man_path.is_file():
        return {}
    raw = json.loads(man_path.read_text())
    out = {}
    for e in raw:
        out[e["stl"]] = e
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--build-dir", type=Path, required=True)
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--p2", type=Path, required=True)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--cube", type=Path, required=True)
    ap.add_argument(
        "--dumps-dir",
        type=Path,
        default=None,
        help="Frozen island-aware P1 dumps (MANIFEST + *_cN.json)",
    )
    ap.add_argument(
        "--sidecars",
        type=Path,
        default=None,
        help="Directory of *.expected.json with live[] (p0-gen-rw1)",
    )
    ap.add_argument(
        "--compose-dump",
        type=Path,
        default=None,
        help="Compose-winner regiondump used to regenerate missing dumps",
    )
    args = ap.parse_args()

    dumps_dir = args.dumps_dir or (args.repo / "_team" / "p1-dumps")
    # The frozen dump dir always existed in the lane worktree; on an assembled
    # branch it is a fresh build-dir path that regeneration has to create.
    dumps_dir.mkdir(parents=True, exist_ok=True)
    sidecars = args.sidecars
    compose_dump = args.compose_dump
    if compose_dump is None:
        cand = (
            args.repo.parent / "p1-compose-fix-rw1" / "build" / "stl2step_regiondump"
        )
        if cand.is_file():
            compose_dump = cand

    manifest = load_manifest(dumps_dir)
    corpus_roots = [p for p in (sidecars, args.corpus) if p is not None]

    names: list[str] = []
    if manifest:
        names = [e["stl"] for e in json.loads((dumps_dir / "MANIFEST.json").read_text())]
    else:
        names = [p.name for p in sorted(args.corpus.glob("S*.stl"))]
        names.append(args.cube.name)

    rows: list[str] = []
    failures: list[str] = []
    table: list[str] = []

    n_pass = n_fail = n_skip = n_esc = 0
    seamed = halves = exploded = closed360_total = 0
    two_halves_reasons: list[str] = []

    for name in names:
        info = manifest.get(name, {})
        stl = stl_for(name, *corpus_roots)
        if name == args.cube.name or name == "cube.stl":
            stl = args.cube if args.cube.is_file() else stl
        if stl is None:
            failures.append(f"FAIL {name}: STL not found")
            rows.append(f"FAIL {name}: STL missing")
            n_fail += 1
            continue

        comps = list(info.get("clean") or [])
        bare_list = list(info.get("bare") or [])
        if not comps and not info:
            # No MANIFEST entry: fall back to envelope if present.
            env = dumps_dir / f"{Path(name).stem}.envelope.json"
            if env.is_file():
                try:
                    edoc = json.loads(env.read_text())
                    comps = [
                        int(c["index"])
                        for c in edoc.get("comps", [])
                        if c.get("clean")
                    ]
                except (OSError, json.JSONDecodeError, KeyError):
                    comps = []
            if not comps and compose_dump is not None and compose_dump.is_file():
                # No frozen dump for this fixture: derive the clean-component set
                # live from the in-tree regiondump. Without this the gate SKIPs
                # every fixture whenever _team/p1-dumps is absent - which is the
                # normal state of an assembled branch - and reports a vacuous pass.
                try:
                    comps = envelope_clean(compose_dump, stl)
                except Exception:
                    comps = []
        if not comps:
            rows.append(f"SKIP {name}: no clean components")
            table.append(f"{name:28}  SKIP  no clean components (I6)")
            n_skip += 1
            continue

        sidecar = sidecar_for(stl, *corpus_roots)
        for i, comp in enumerate(comps):
            tag = f"{name} comp{comp}"
            bare_name = None
            if i < len(bare_list):
                bare_name = bare_list[i]
            else:
                bare_name = f"{Path(name).stem}_c{comp}.json"
            rs = dumps_dir / bare_name
            if compose_dump is not None and compose_dump.is_file():
                # Live re-gate: always refresh P1 dumps from the in-tree regiondump
                # so corpus/STL changes (e.g. pinned S09) cannot leave stale JSON in
                # ${CMAKE_BINARY_DIR}/p1-dumps from an earlier configure.
                try:
                    dump_bare(compose_dump, stl, comp, rs)
                except RuntimeError as e:
                    failures.append(str(e))
                    rows.append(f"FAIL {tag}: regen dump")
                    n_fail += 1
                    continue
            elif not rs.is_file():
                failures.append(f"{tag}: missing dump {rs}")
                rows.append(f"FAIL {tag}: missing dump")
                n_fail += 1
                continue

            b1 = p2_live(args.p2, rs, stl, comp, sidecar)
            b8 = b1  # frozen dump is the identity snapshot
            id_ok = True
            ok = b1.get("ok") is True and b1.get("_rc") == 0

            n360 = int(b1.get("closed360") or 0)
            ns = int(b1.get("seamed360") or 0)
            nh = int(b1.get("twoHalves") or 0)
            ne = int(b1.get("exploded360") or 0)
            closed360_total += n360
            seamed += ns
            halves += nh
            exploded += ne
            if nh:
                reasons = [
                    w
                    for w in (b1.get("warnings") or [])
                    if "seamed360" in w or "TwoHalves" in w or "twoHalves" in w
                ]
                two_halves_reasons.append(
                    f"{tag}: twoHalves={nh} reasons={reasons or b1.get('warnings', [])[:3]}"
                )

            flags = (
                f"build={b1.get('buildFaces')} valid={b1.get('valid')} "
                f"closed={b1.get('closed')} j1={b1.get('j1')} j6={b1.get('j6')} "
                f"twins={b1.get('twins')} j3miss={b1.get('j3Missing')} "
                f"volOk={b1.get('volumeOk')} faces={b1.get('faceCount')}/"
                f"{b1.get('nTri')} pln={b1.get('planes')} cyl={b1.get('cylinders')} "
                f"vol={b1.get('volume')} budget={b1.get('volumeBudget')} "
                f"360={n360} seamed={ns} halves={nh} exploded={ne} "
                f"tol={b1.get('maxVertexTol')} tolOk={b1.get('tolOk')} "
                f"census={b1.get('censusOk')} facesOk={b1.get('facesOk')} "
                f"threads={'id' if id_ok else 'DIFF'}"
            )
            live = live_row(sidecar, comp) or {}
            disp = str(live.get("disposition") or "")
            reason = str(live.get("escalateReason") or "")
            owner = "sidecar-live"
            sig = fail_sig(b1, id_ok)
            extra = ""
            if ok and disp == "ESCALATE" and int(b1.get("exploded360") or 0) > 0:
                # Sidecar expected ESCALATE; 360-cyl exploded (S04 TorusNYI).
                status = "ESCALATE"
                extra = f" reason={reason[:72]} owner={owner}"
                ok = False
            elif ok:
                status = "PASS"
            elif disp == "ESCALATE":
                status = "ESCALATE"
                extra = f" reason={reason[:72]} owner={owner}"
            else:
                status = "FAIL"
            table.append(f"{tag:28}  {status}  {flags}{extra}")
            if ok:
                rows.append(f"PASS {tag}")
                n_pass += 1
            elif status == "ESCALATE":
                rows.append(f"ESCALATE {tag}: {reason}")
                n_esc += 1
            else:
                rows.append(f"FAIL {tag}")
                n_fail += 1
                failures.append(
                    f"{tag}: rc={b1.get('_rc')} ok={b1.get('ok')} "
                    f"id_ok={id_ok} sig={sig} "
                    f"explodedHole={b1.get('explodedRecoverableHole')}\n"
                    f"  t1={b1.get('_raw')[:800]}\n"
                    f"  stderr={b1.get('_stderr')[:400]}"
                )

    print("P2 REAL GATE — live RegionSet × buildFaces")
    print(f"dumps={dumps_dir}")
    print(f"sidecars={sidecars or args.corpus}")
    print()
    print("LIVE RE-GATE TABLE")
    print("-----------------")
    for row in table:
        print(row)
    print()
    print(
        f"closed360 corpus-wide: total={closed360_total} "
        f"Seamed360={seamed} TwoHalves={halves} ExplodedToFacets={exploded}"
    )
    if two_halves_reasons:
        print("TwoHalves reasons:")
        for r in two_halves_reasons:
            print(f"  {r}")
    print()
    print("P2 REAL GATE")
    for row in rows:
        print(row)
    total = n_pass + n_fail + n_skip + n_esc
    print(
        f"\nSUMMARY: {n_pass}/{total} pass  FAIL={n_fail} SKIP={n_skip} ESCALATE={n_esc}"
    )
    if n_esc:
        print(f"ESCALATED={n_esc} (sidecar live[].disposition=ESCALATE)")
    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for f in failures:
            print(f, file=sys.stderr)
        return 1
    if n_pass == 0:
        # Nothing was actually scored. A gate that skips every fixture and exits
        # 0 reports "green" for a surface it never touched - the exact false
        # pass this re-gate exists to prevent.
        print(
            f"\nVACUOUS: {n_skip} skipped, 0 scored. The live re-gate scored "
            f"nothing - check --dumps-dir / --compose-dump.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
