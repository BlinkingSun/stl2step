#!/usr/bin/env python3
"""Run Body11-class repro STLs on the tree's stl2step. Exit 0 if each
mesh still emits its class warn (today's engine). Prints command+output."""
from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BIN = ROOT / "build" / "stl2step"
HERE = Path(__file__).resolve().parent

CASES = [
    ("me_plnpln_transversal.stl", "smooth: analytic MakeEdge failed"),
    ("me_cylpln_ellipse.stl", "smooth: analytic MakeEdge failed"),
    ("ia_cylcyl_nogeom.stl", "smooth: IntAna cyl|cyl empty/same"),
    ("seamed360_capseam.stl", "seamed360: BRepCheck invalid on seamed face"),
]


def main() -> int:
    if not BIN.is_file():
        print(f"missing {BIN}", file=sys.stderr)
        return 2
    rc = 0
    for name, needle in CASES:
        stl = HERE / name
        step = Path("/tmp") / f"f2diag-{name}.step"
        cmd = [str(BIN), str(stl), "-o", str(step), "--smooth", "--threads", "1",
               "--no-verify", "--quiet"]
        print("CMD", " ".join(cmd))
        p = subprocess.run(cmd, capture_output=True, text=True)
        out = (p.stdout or "") + (p.stderr or "")
        print("EXIT", p.returncode)
        print(out.strip()[-2000:])
        line = ""
        for L in (p.stdout or "").splitlines():
            if L.startswith("RESULT "):
                line = L
        if not line:
            print("FAIL no RESULT")
            rc = 1
            continue
        r = json.loads(line[len("RESULT "):])
        c = Counter(r.get("warnings") or [])
        hit = any(needle in w for w in c)
        print("CLASS", needle, "HIT" if hit else "MISS", dict(c))
        if not hit:
            rc = 1
        print()
    return rc


if __name__ == "__main__":
    sys.exit(main())
