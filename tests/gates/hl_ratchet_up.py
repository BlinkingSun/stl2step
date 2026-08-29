#!/usr/bin/env python3
"""One-line orchestrator bumper: increment (or set) hl-ratchet.json floor."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hl_census_ratchet import main  # noqa: E402

args = ["--ratchet-up", *sys.argv[1:]] if "--ratchet-up" not in sys.argv[1:] else sys.argv[1:]
raise SystemExit(main(args))
