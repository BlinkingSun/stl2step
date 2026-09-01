#!/usr/bin/env python3
"""ctest drivers for stl2step --mesh mode (U3a-2 fixtures lane).

Exercises mesh_step_smoke, mesh_edges_seam, and mesh_fail_closed via --test.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


def last_stdout_line(proc: subprocess.CompletedProcess[str]) -> str:
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return lines[-1] if lines else ""


def parse_mesh_result(line: str) -> dict:
    if not line.startswith("MESH_RESULT "):
        raise AssertionError(f"expected MESH_RESULT line, got {line!r}")
    return json.loads(line[len("MESH_RESULT ") :])


def run(binary: str, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [binary, *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


def ensure_step(binary: str, dest: Path, step_src: Path | None, stl_src: Path | None) -> None:
    if step_src is not None and step_src.is_file():
        shutil.copy2(step_src, dest)
        return
    if stl_src is None or not stl_src.is_file():
        raise SystemExit("FAIL: need --step (existing) or --make-step-from <stl>")
    stl_copy = dest.with_suffix(".src.stl")
    shutil.copy2(stl_src, stl_copy)
    proc = run(binary, [str(stl_copy), "-o", str(dest), "--quiet", "--no-verify"])
    if proc.returncode not in (0, 2) or not dest.is_file():
        raise SystemExit(
            "FAIL: could not make a STEP for mesh_step_smoke\n"
            f"exit={proc.returncode}\nstdout={proc.stdout}\nstderr={proc.stderr}"
        )


def assert_binary_stl(path: Path, triangles: int) -> None:
    data = path.read_bytes()
    if data[:5].lower() == b"solid":
        raise AssertionError("STL is ASCII (starts with 'solid')")
    expected = 84 + 50 * triangles
    if len(data) != expected:
        raise AssertionError(f"STL size {len(data)} != 84 + 50*{triangles} = {expected}")


def facet_edge_upper_bound(triangles: int) -> int:
    return triangles * 3


def count_edge_segments(edges_path: Path) -> int:
    size = edges_path.stat().st_size
    if size % 24 != 0:
        raise AssertionError(f"edges file size {size} not divisible by 24")
    return size // 24


def mesh_step_smoke(binary: str, step_src: Path | None, stl_src: Path | None) -> int:
    with tempfile.TemporaryDirectory(prefix="mesh-step-smoke-") as td:
        td_path = Path(td)
        step = td_path / "input.step"
        stl_out = td_path / "out.stl"
        edges_out = td_path / "out.edges"

        ensure_step(binary, step, step_src, stl_src)

        proc = run(
            binary,
            ["--mesh", str(step), "-o", str(stl_out), "--edges", str(edges_out), "--quiet", "--threads", "4"],
        )
        if proc.returncode != 0:
            print("FAIL smoke: exit != 0")
            print("stdout:", proc.stdout)
            print("stderr:", proc.stderr)
            return 1

        out_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
        if len(out_lines) != 1:
            print(f"FAIL smoke: expected exactly 1 stdout line, got {len(out_lines)}")
            return 1
        line = out_lines[0]
        if not line.startswith("MESH_RESULT "):
            print(f"FAIL smoke: line does not start with MESH_RESULT: {line!r}")
            return 1

        r = parse_mesh_result(line)
        if r.get("ok") is not True:
            print("FAIL smoke: ok is not true", r)
            return 1
        if r.get("faces", 0) <= 0 or r.get("triangles", 0) <= 0:
            print("FAIL smoke: faces/triangles must be > 0", r)
            return 1
        if not stl_out.is_file():
            print("FAIL smoke: STL not written")
            return 1
        if not edges_out.is_file():
            print("FAIL smoke: edges file not written")
            return 1

        triangles = int(r["triangles"])
        try:
            assert_binary_stl(stl_out, triangles)
        except AssertionError as e:
            print(f"FAIL smoke: {e}")
            return 1
        if edges_out.stat().st_size % 24 != 0:
            print("FAIL smoke: edges size not multiple of 24")
            return 1

        print("PASS mesh_step_smoke")
        print(f"  faces={r['faces']} triangles={triangles} edges={r.get('edges')}")
        return 0


def mesh_edges_seam(binary: str, cylinder_step: Path, *, expect_edges: int = 3) -> int:
    with tempfile.TemporaryDirectory(prefix="mesh-edges-seam-") as td:
        td_path = Path(td)
        stl_out = td_path / "cyl.stl"
        edges_out = td_path / "cyl.edges"

        proc = run(
            binary,
            [
                "--mesh",
                str(cylinder_step),
                "-o",
                str(stl_out),
                "--edges",
                str(edges_out),
                "--quiet",
                "--threads",
                "4",
            ],
        )
        if proc.returncode != 0:
            print("FAIL seam: exit != 0")
            print(proc.stdout, proc.stderr)
            return 1

        line = last_stdout_line(proc)
        r = parse_mesh_result(line)
        if r.get("ok") is not True:
            print("FAIL seam: ok is not true", r)
            return 1
        if r.get("edges") != expect_edges:
            print(f"FAIL seam: edges={r.get('edges')} expected {expect_edges}")
            return 1

        segments = count_edge_segments(edges_out)
        triangles = int(r.get("triangles", 0))
        facet_edges = facet_edge_upper_bound(triangles)
        if segments >= facet_edges // 2:
            print(
                f"FAIL seam: segments={segments} not << facet_edges={facet_edges} "
                "(buffer looks like facet boundaries)"
            )
            return 1

        # Classify: 2 tessellated circles + 1 straight seam (parallel checks).
        data = edges_out.read_bytes()
        seam = cap0 = cap_h = other = 0
        for off in range(0, len(data), 24):
            x0, y0, z0, x1, y1, z1 = struct.unpack("<6f", data[off : off + 24])
            za, zb = z0, z1
            ra = (x0 * x0 + y0 * y0) ** 0.5
            rb = (x1 * x1 + y1 * y1) ** 0.5
            if abs(x0 - x1) < 1e-2 and abs(y0 - y1) < 1e-2 and abs(abs(za - zb) - 30.0) < 0.6:
                seam += 1
            elif abs(za) < 0.6 and abs(zb) < 0.6 and abs(ra - 10.0) < 0.6 and abs(rb - 10.0) < 0.6:
                cap0 += 1
            elif abs(za - 30.0) < 0.6 and abs(zb - 30.0) < 0.6 and abs(ra - 10.0) < 0.6 and abs(rb - 10.0) < 0.6:
                cap_h += 1
            else:
                other += 1

        if seam != 1 or cap0 < 1 or cap_h < 1 or other != 0:
            print(
                f"FAIL seam: classify seam={seam} cap0={cap0} cap_h={cap_h} other={other} "
                f"(expected 1 seam + 2 caps, no facet edges)"
            )
            return 1

        print("PASS mesh_edges_seam")
        print(f"  edges={expect_edges} segments={segments} facet_upper={facet_edges}")
        print(f"  classify seam={seam} cap0={cap0} cap_h={cap_h}")
        return 0


def mesh_fail_closed(binary: str, step_src: Path | None, stl_src: Path | None) -> int:
    with tempfile.TemporaryDirectory(prefix="mesh-fail-closed-") as td:
        step = Path(td) / "input.step"
        ensure_step(binary, step, step_src, stl_src)

        proc = run(binary, ["--mesh", str(step), "--engine", "trueform", "--quiet"])
        if proc.returncode != 1:
            print(f"FAIL fail_closed: expected exit 1, got {proc.returncode}")
            return 1
        for ln in proc.stdout.splitlines():
            if ln.startswith("RESULT ") or ln.startswith("MESH_RESULT "):
                print(f"FAIL fail_closed: unexpected stdout contract line: {ln!r}")
                return 1
        if "not valid with --mesh" not in proc.stderr:
            print("FAIL fail_closed: expected rejection on stderr")
            print(proc.stderr)
            return 1

        print("PASS mesh_fail_closed")
        return 0


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--binary", required=True)
    p.add_argument("--test", required=True, choices=["smoke", "edges_seam", "fail_closed"])
    p.add_argument("--step", type=Path, default=None)
    p.add_argument("--make-step-from", type=Path, default=None)
    p.add_argument("--cylinder-step", type=Path, default=None)
    p.add_argument("--expect-edges", type=int, default=3)
    args = p.parse_args()

    binary = args.binary
    if not os.path.isfile(binary):
        raise SystemExit(f"FAIL: binary not found: {binary}")

    if args.test == "smoke":
        return mesh_step_smoke(binary, args.step, args.make_step_from)
    if args.test == "edges_seam":
        if args.cylinder_step is None or not args.cylinder_step.is_file():
            raise SystemExit("FAIL: --cylinder-step required for edges_seam")
        return mesh_edges_seam(binary, args.cylinder_step, expect_edges=args.expect_edges)
    return mesh_fail_closed(binary, args.step, args.make_step_from)


if __name__ == "__main__":
    sys.exit(main())
