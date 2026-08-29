#!/usr/bin/env python3
"""Adversarial audit: keepValidPartials + probe waiver must not rescue broken shells.

Run from repo root after building:
  python3 tests/diag/handle-lock/waiver_audit.py

Positive control: intact handle-lock reverts to faceted when analytic shell is
invalid or volume-busts D4.5 budget (no waiver rescue).
Negative: delete triangles from partial-cylinder region rid=16 (15 mm bore) — revert.
"""
from __future__ import annotations

import json
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
STL = Path(__file__).resolve().parent / "handle-lock.stl"
STL2STEP = REPO / "build" / "stl2step"
REGIONDUMP = REPO / "build" / "stl2step_regiondump"
CENSUS = REPO / "build" / "stl2step_census"


def read_binary_stl(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    data = path.read_bytes()
    if len(data) < 84:
        raise ValueError("STL too short")
    n = struct.unpack_from("<I", data, 80)[0]
    off = 84
    verts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    vmap: dict[tuple[float, float, float], int] = {}

    def vid(p: tuple[float, float, float]) -> int:
        if p not in vmap:
            vmap[p] = len(verts)
            verts.append(p)
        return vmap[p]

    for _ in range(n):
        tri = struct.unpack_from("<12fH", data, off)
        off += 50
        i0 = vid((tri[3], tri[4], tri[5]))
        i1 = vid((tri[6], tri[7], tri[8]))
        i2 = vid((tri[9], tri[10], tri[11]))
        tris.append((i0, i1, i2))
    return verts, tris


def write_binary_stl(path: Path, verts: list[tuple[float, float, float]], tris: list[tuple[int, int, int]]) -> None:
    with path.open("wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tris)))
        for i0, i1, i2 in tris:
            p0, p1, p2 = verts[i0], verts[i1], verts[i2]
            nx = (p1[1] - p0[1]) * (p2[2] - p0[2]) - (p1[2] - p0[2]) * (p2[1] - p0[1])
            ny = (p1[2] - p0[2]) * (p2[0] - p0[0]) - (p1[0] - p0[0]) * (p2[2] - p0[2])
            nz = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])
            ln = (nx * nx + ny * ny + nz * nz) ** 0.5 or 1.0
            nx, ny, nz = nx / ln, ny / ln, nz / ln
            f.write(
                struct.pack(
                    "<12fH",
                    nx,
                    ny,
                    nz,
                    *p0,
                    *p1,
                    *p2,
                    0,
                )
            )


def region_tri_indices() -> list[int]:
    """Mesh triangle indices of the 15 mm bore (rid was 16 at 34ae192)."""
    proc = subprocess.run(
        [str(REGIONDUMP), str(STL)],
        check=True,
        capture_output=True,
        text=True,
    )
    doc = json.loads(proc.stdout)
    for r in doc["comps"][0]["regionSet"]["regions"]:
        if str(r.get("type", "")).lower() != "cylinder":
            continue
        rad = float(r.get("radius", 0.0))
        if rad < 14.95 or rad > 15.05:
            continue
        tris = r.get("tris") or []
        if len(tris) >= 20:
            return list(tris)
    return []


def run_trueform(stl: Path, step: Path, *, verify: bool = True) -> dict:
    cmd = [str(STL2STEP), str(stl), "-o", str(step), "--engine", "trueform", "--quiet"]
    if not verify:
        cmd.append("--no-verify")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    line = proc.stdout.strip().splitlines()[-1]
    if not line.startswith("RESULT "):
        raise RuntimeError(f"no RESULT line: {proc.stdout!r} {proc.stderr!r}")
    return json.loads(line[len("RESULT ") :])


def per_face_brep_valid(step: Path) -> tuple[bool, int, int]:
  """Whole-shape + per-face BRepCheck via stl2step_census internals (shape-level)."""
  proc = subprocess.run([str(CENSUS), str(step)], check=True, capture_output=True, text=True)
  r = json.loads(proc.stdout)
  return bool(r.get("valid")), int(r.get("faces", 0)), int(r.get("openShells", 0))


def main() -> int:
    if not STL2STEP.is_file():
        print("build stl2step first", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        good_step = td_path / "good.step"
        good = run_trueform(STL, good_step)
        valid, faces, open_shells = per_face_brep_valid(good_step)

        print("=== positive control (intact handle-lock) ===")
        print(
            f"smoothBuiltCylinders={good.get('smoothBuiltCylinders')} "
            f"smoothRevertedComponents={good.get('smoothRevertedComponents')} "
            f"watertight={good.get('watertight')} openShells={good.get('openShells')} "
            f"volumeDeltaPct={good.get('volumeDeltaPct')}"
        )
        print(f"census valid={valid} faces={faces} openShells={open_shells}")

        vol_pct = float(good.get("volumeDeltaPct", 999))
        reverted = good.get("smoothRevertedComponents", 0) > 0
        vol_ok = vol_pct >= 0 and vol_pct <= 6.02
        if not good.get("watertight") or good.get("openShells", 1) != 0:
            print("FAIL: positive control not watertight")
            return 1
        if not reverted and not vol_ok:
            print("FAIL: positive control kept invalid analytic build")
            return 1
        if reverted and not vol_ok:
            print("FAIL: reverted but volume still busts budget")
            return 1
        print("PASS: positive control volume-safe (reverted or valid analytic)")

        # Negative: strip all tris from rid=16 partial — breaks bore patch + opens shell.
        rid16 = region_tri_indices()
        if len(rid16) < 4:
            print(f"FAIL: expected rid=16 tris, got {len(rid16)}")
            return 1
        verts, tris = read_binary_stl(STL)
        drop = set(rid16)
        broken_tris = [t for i, t in enumerate(tris) if i not in drop]
        broken_stl = td_path / "broken-partial.stl"
        write_binary_stl(broken_stl, verts, broken_tris)
        broken_step = td_path / "broken.step"
        broken = run_trueform(broken_stl, broken_step, verify=False)

        print("=== negative (rid=16 partial deleted, nTri drop=%d) ===" % len(rid16))
        print(
            f"smoothBuiltCylinders={broken.get('smoothBuiltCylinders')} "
            f"smoothRevertedComponents={broken.get('smoothRevertedComponents')} "
            f"watertight={broken.get('watertight')} openShells={broken.get('openShells')}"
        )
        if broken.get("warnings"):
            print("warnings:", broken["warnings"][:3])

        # Waiver must NOT rescue: no analytic census on broken mesh.
        rescued = broken.get("smoothBuiltCylinders", 0) >= 10 or (
            broken.get("smoothBuiltComponents", 0) > 0
            and broken.get("smoothRevertedComponents", 0) == 0
            and broken.get("watertight")
        )
        if rescued:
            print("FAIL: waiver rescued broken partial")
            return 1

        print("PASS: broken partial not rescued by waiver")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
