#!/usr/bin/env python3
"""ctest result_validity_check — RESULT reports a skipped validity check (#24).

BRepCheck runs only up to 50,000 faces. Above that the check is skipped, and
RESULT must say so with "validityCheck":"skipped" without turning the skip
into a warning (the exit code is unchanged). At or under the cap the key is
absent, so RESULT stays character-identical to earlier releases.

Large case: a closed sphere whose rings are offset by half a step, so no two
neighbouring triangles are coplanar and the coplanar merge keeps every face.
Stdlib only.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

CAP = 50000


def write_sphere_stl(path: Path, rings: int, segs: int, radius: float = 50.0) -> int:
    """Binary STL of a closed, outward-wound sphere. Returns the triangle count."""

    def ring_pt(i: int, j: int) -> tuple[float, float, float]:
        theta = math.pi * i / rings
        # Wrap j so the seam reuses the exact same coordinates (welds cleanly).
        phi = 2.0 * math.pi * (j % segs + 0.5 * i) / segs
        s = math.sin(theta)
        return (radius * s * math.cos(phi), radius * s * math.sin(phi), radius * math.cos(theta))

    north = (0.0, 0.0, radius)
    south = (0.0, 0.0, -radius)
    tris: list[tuple[tuple[float, float, float], ...]] = []
    for j in range(segs):
        tris.append((north, ring_pt(1, j), ring_pt(1, j + 1)))
    for i in range(1, rings - 1):
        for j in range(segs):
            a0, a1 = ring_pt(i, j), ring_pt(i, j + 1)
            b0, b1 = ring_pt(i + 1, j), ring_pt(i + 1, j + 1)
            tris.append((a0, b0, a1))
            tris.append((a1, b0, b1))
    for j in range(segs):
        tris.append((south, ring_pt(rings - 1, j + 1), ring_pt(rings - 1, j)))

    with path.open("wb") as fh:
        fh.write(b"\0" * 80)
        fh.write(struct.pack("<I", len(tris)))
        for t in tris:
            fh.write(struct.pack("<3f", 0.0, 0.0, 0.0))
            for p in t:
                fh.write(struct.pack("<3f", *p))
            fh.write(b"\0\0")
    return len(tris)


def convert(binary: str, stl: Path, step: Path) -> tuple[int, dict, str]:
    p = subprocess.run(
        [binary, str(stl), "-o", str(step), "--no-verify", "--engine", "verbatim"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    lines = [ln for ln in p.stdout.replace("\r\n", "\n").split("\n") if ln.startswith("RESULT ")]
    if not lines:
        raise SystemExit(f"FAIL {stl.name}: no RESULT line (exit={p.returncode})\n{p.stderr}")
    return p.returncode, json.loads(lines[-1][len("RESULT "):]), p.stdout + p.stderr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--small-stl", required=True)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory(prefix="result-validity-") as td:
        td_path = Path(td)

        code, res, _ = convert(args.binary, Path(args.small_stl), td_path / "small.step")
        if code not in (0, 2) or not res.get("ok"):
            raise SystemExit(f"FAIL small: exit={code} result={res}")
        if res["facesAfterUnify"] > CAP:
            raise SystemExit(f"FAIL small: fixture has {res['facesAfterUnify']} faces, above the cap")
        if "validityCheck" in res:
            raise SystemExit(f"FAIL small: validityCheck emitted under the cap: {res}")

        big = td_path / "sphere.stl"
        ntri = write_sphere_stl(big, rings=120, segs=220)
        code, res, log = convert(args.binary, big, td_path / "sphere.step")
        if code not in (0, 2) or not res.get("ok"):
            raise SystemExit(f"FAIL large: exit={code} result={res}")
        faces = res["facesAfterUnify"]
        if not res["watertight"] or res["solids"] != 1:
            raise SystemExit(f"FAIL large: generated sphere is not one closed solid: {res}")
        if faces <= CAP:
            raise SystemExit(f"FAIL large: only {faces} faces from {ntri} triangles; need > {CAP}")
        if res.get("validityCheck") != "skipped":
            raise SystemExit(f"FAIL large: expected validityCheck=skipped at {faces} faces: {res}")
        if "check     skipped" not in log:
            raise SystemExit("FAIL large: progress log no longer reports the skip")
        if any("validity" in w.lower() for w in res["warnings"]):
            raise SystemExit(f"FAIL large: the skip became a warning: {res['warnings']}")
        if code != (2 if res["warnings"] else 0):
            raise SystemExit(f"FAIL large: exit {code} does not follow warnings {res['warnings']}")

    print(f"result_validity_check: PASS (small: no key; sphere {faces} faces: skipped, exit {code})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
