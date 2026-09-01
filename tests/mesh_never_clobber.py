#!/usr/bin/env python3
"""Never-clobber default for `stl2step --mesh` without `-o`.

(a) Next to a pre-existing <stem>.stl, `--mesh <stem>.step --quiet` must leave
    that STL byte-identical and write <stem>.mesh.stl.
(b) A second run (default path now exists) must exit 1 with MESH_RESULT
    ok=false and error "output exists: <path> — pass -o".

All work is in a temp directory. Corpus fixture bytes are never written.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def last_stdout_line(proc: subprocess.CompletedProcess[str]) -> str:
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return lines[-1] if lines else ""


def parse_mesh_result(line: str) -> dict:
    if not line.startswith("MESH_RESULT "):
        raise SystemExit(f"FAIL: expected MESH_RESULT line, got {line!r}")
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
    # Convert a copy of the STL in the temp dir so corpus bytes stay untouched.
    stl_copy = dest.with_suffix(".src.stl")
    shutil.copy2(stl_src, stl_copy)
    proc = run(binary, [str(stl_copy), "-o", str(dest), "--quiet", "--no-verify"])
    if proc.returncode not in (0, 2) or not dest.is_file():
        raise SystemExit(
            "FAIL: could not make a STEP for the never-clobber test\n"
            f"exit={proc.returncode}\nstdout={proc.stdout}\nstderr={proc.stderr}"
        )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--binary", required=True)
    p.add_argument("--step", type=Path, default=None)
    p.add_argument("--make-step-from", type=Path, default=None)
    args = p.parse_args()

    binary = args.binary
    if not os.path.isfile(binary):
        raise SystemExit(f"FAIL: binary not found: {binary}")

    with tempfile.TemporaryDirectory(prefix="mesh-never-clobber-") as td:
        td_path = Path(td)
        step = td_path / "handle-lock.step"
        sibling_stl = td_path / "handle-lock.stl"
        default_out = td_path / "handle-lock.mesh.stl"

        ensure_step(binary, step, args.step, args.make_step_from)

        marker = b"SOURCE-MESH-MUST-NOT-BE-CLOBBERED\n" * 32
        sibling_stl.write_bytes(marker)
        before = sha256(sibling_stl)

        # (a) first run: default output is <stem>.mesh.stl
        a = run(binary, ["--mesh", str(step), "--quiet"])
        if a.returncode != 0:
            print("FAIL a: expected exit 0")
            print("stdout:", a.stdout)
            print("stderr:", a.stderr)
            return 1
        if sha256(sibling_stl) != before:
            print("FAIL a: sibling <stem>.stl was modified")
            return 1
        if sibling_stl.read_bytes() != marker:
            print("FAIL a: sibling STL bytes changed")
            return 1
        if not default_out.is_file():
            print(f"FAIL a: output did not land at {default_out}")
            print("stdout:", a.stdout)
            return 1
        line = last_stdout_line(a)
        r = parse_mesh_result(line)
        if r.get("ok") is not True:
            print("FAIL a: ok is not true", r)
            return 1
        out_name = Path(r.get("output", "")).name
        if out_name != "handle-lock.mesh.stl":
            print(f"FAIL a: MESH_RESULT output name {out_name!r}")
            return 1
        if "edgesFile" in r:
            print("FAIL a: edgesFile present without --edges", r)
            return 1

        # (b) second run: default path exists → fail closed, no rewrite of source
        size_before = default_out.stat().st_size
        b = run(binary, ["--mesh", str(step), "--quiet"])
        if b.returncode != 1:
            print(f"FAIL b: expected exit 1, got {b.returncode}")
            print("stdout:", b.stdout)
            print("stderr:", b.stderr)
            return 1
        line_b = last_stdout_line(b)
        r2 = parse_mesh_result(line_b)
        if r2.get("ok") is not False:
            print("FAIL b: ok is not false", r2)
            return 1
        err = r2.get("error", "")
        expect_prefix = "output exists: "
        expect_suffix = " — pass -o"
        if not err.startswith(expect_prefix) or not err.endswith(expect_suffix):
            print(f"FAIL b: error {err!r}")
            return 1
        if "handle-lock.mesh.stl" not in err:
            print(f"FAIL b: error path missing default name: {err!r}")
            return 1
        if sha256(sibling_stl) != before:
            print("FAIL b: sibling <stem>.stl was modified on second run")
            return 1
        if default_out.stat().st_size != size_before:
            print("FAIL b: existing default output was rewritten")
            return 1
        # No RESULT line (convert contract) on a mesh-mode fail-closed run.
        if any(ln.startswith("RESULT ") for ln in b.stdout.splitlines()):
            print("FAIL b: unexpected RESULT line on stdout")
            return 1

        print("PASS mesh_never_clobber")
        print(f"  sibling_stl sha256={before} unchanged")
        print(f"  default_out={default_out} bytes={size_before}")
        print(f"  b.error={err}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
