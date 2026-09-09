#!/usr/bin/env python3
"""ctest result_line_redirected — RESULT / MESH_RESULT survive a pipe and a file (#7).

Runs the built CLI with stdout redirected to a FILE and with stdout=PIPE.
Asserts the contract line is present and is the last line (convert: RESULT;
--mesh: MESH_RESULT). Stdlib only.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def last_line(text: str) -> str:
    body = text.replace("\r\n", "\n").replace("\r", "\n")
    if body.endswith("\n"):
        body = body[:-1]
    if not body:
        return ""
    return body.split("\n")[-1]


def assert_contract(stdout: str, prefix: str, where: str) -> None:
    last = last_line(stdout)
    if not last.startswith(prefix + " "):
        raise SystemExit(
            f"FAIL {where}: last line is not {prefix} …\n"
            f"  last={last!r}\n  stdout={stdout!r}"
        )
    blob = last[len(prefix) + 1 :]
    try:
        obj = json.loads(blob)
    except json.JSONDecodeError as e:
        raise SystemExit(f"FAIL {where}: {prefix} JSON: {e}\n  blob={blob!r}") from e
    if not isinstance(obj, dict):
        raise SystemExit(f"FAIL {where}: {prefix} JSON is not an object: {obj!r}")


def run_to_file(cmd: list[str], path: Path) -> subprocess.CompletedProcess[str]:
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        return subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )


def run_to_pipe(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--stl", required=True)
    args = ap.parse_args()
    binary = args.binary
    stl = Path(args.stl)
    if not stl.is_file():
        raise SystemExit(f"FAIL: missing STL {stl}")

    with tempfile.TemporaryDirectory(prefix="result-line-redir-") as td:
        td_path = Path(td)
        step = td_path / "out.step"
        stl_out = td_path / "out.mesh.stl"
        file_out = td_path / "stdout.file.txt"

        convert = [
            binary,
            str(stl),
            "-o",
            str(step),
            "--quiet",
            "--no-verify",
            "--engine",
            "verbatim",
        ]

        p_file = run_to_file(convert, file_out)
        if p_file.returncode not in (0, 2):
            raise SystemExit(
                f"FAIL convert>file: exit={p_file.returncode} stderr={p_file.stderr}"
            )
        assert_contract(file_out.read_text(encoding="utf-8"), "RESULT", "convert>file")
        if not step.is_file():
            raise SystemExit("FAIL convert>file: no STEP written")

        p_pipe = run_to_pipe(convert)
        if p_pipe.returncode not in (0, 2):
            raise SystemExit(
                f"FAIL convert|PIPE: exit={p_pipe.returncode} stderr={p_pipe.stderr}"
            )
        assert_contract(p_pipe.stdout, "RESULT", "convert|PIPE")

        mesh = [
            binary,
            "--mesh",
            str(step),
            "-o",
            str(stl_out),
            "--quiet",
        ]
        mesh_file = td_path / "mesh.stdout.file.txt"
        p_mfile = run_to_file(mesh, mesh_file)
        if p_mfile.returncode != 0:
            raise SystemExit(
                f"FAIL mesh>file: exit={p_mfile.returncode} stderr={p_mfile.stderr}"
            )
        assert_contract(mesh_file.read_text(encoding="utf-8"), "MESH_RESULT", "mesh>file")

        p_mpipe = run_to_pipe(mesh)
        if p_mpipe.returncode != 0:
            raise SystemExit(
                f"FAIL mesh|PIPE: exit={p_mpipe.returncode} stderr={p_mpipe.stderr}"
            )
        assert_contract(p_mpipe.stdout, "MESH_RESULT", "mesh|PIPE")

    print("result_line_redirected: PASS (file + pipe, RESULT + MESH_RESULT)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
