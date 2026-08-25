#!/usr/bin/env python3
"""G0.1 identity comparator for stl2step (python3, stdlib only).

Compare two STEP files or two RESULT payloads the way SPEC-P0 G0.1
requires: DATA-section STEP (HEADER discarded; FILE_NAME timestamp is
in the HEADER so it is not compared), and RESULT key-set + values with
`seconds` and resolved absolute paths excluded. Key order is part of
the contract.

Exit 0 = canonical forms identical, 1 = different, 2 = usage/IO error.

This tool is the off-path instrument. It FAILS if any `smooth*` key
(or `facesAfterSmooth`) appears in either RESULT: G0.1 is `--smooth`
absent. Lane p0-gates should call this via `--baseline DIR`:

    python3 $BASELINE/canonicalize.py step  a.step b.step
    python3 $BASELINE/canonicalize.py result a.txt b.txt
    python3 $BASELINE/canonicalize.py step  a.step b.step --report json

Do not commit golden .step bytes; live-twin the 187ead0 binary this
directory builds against the current CLI.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


# ISO 10303-21 FILE_NAME (name, timestamp, author, organisation,
# preprocessor_version, originating_system, authorization). Field 2 is
# a wall-clock timestamp. Applied to the whole file before DATA extract
# so a future HEADER-including compare cannot accidentally reintroduce
# it; DATA itself does not contain FILE_NAME.
_FILE_NAME_TS = re.compile(
    r"(FILE_NAME\s*\(\s*'(?:[^']|'')*'\s*,\s*)'(?:[^']|'')*'",
    re.IGNORECASE | re.DOTALL,
)

# RESULT JSON keys whose VALUES are resolved paths (or whatever the
# caller passed). Keys stay in the ordered key list; values are dropped.
_PATH_KEYS = ("input", "output")

# Wall-clock; value dropped, key kept.
_TIME_KEYS = ("seconds",)

# G0.1: --smooth absent ⇒ none of these keys may appear.
# `smooth*` covers the DECISION §7 RESULT block; facesAfterSmooth is
# the one sibling that does not match the glob but is equally new.
_SMOOTH_EXACT = {"facesAfterSmooth"}


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _read_text(path: Path) -> str:
    # STEP is ASCII; RESULT is UTF-8 JSON. Replace keeps the tool from
    # crashing on a truncated write; a replacement character is itself
    # a difference.
    return path.read_text(encoding="utf-8", errors="replace")


def _normalize_newlines(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def strip_file_name_timestamp(text: str) -> str:
    """Replace FILE_NAME field 2 (timestamp) with an empty string.

    Minimum HEADER strip that two runs of the same binary require.
    Other HEADER fields are not rewritten here: they are discarded
    wholesale by extract_data_section. See README / the lane report
    for the per-field reasons.
    """
    return _FILE_NAME_TS.sub(r"\1''", text, count=1)


def extract_data_section(text: str) -> str:
    """Return `DATA; ... ENDSEC;` (inclusive). Raises ValueError if missing."""
    upper = text.upper()
    start: Optional[int] = None
    idx = 0
    while True:
        pos = upper.find("DATA;", idx)
        if pos < 0:
            break
        if pos == 0 or text[pos - 1] in "\n\r\t ":
            # Token start: not a suffix of e.g. FOODATA;
            if pos == 0 or text[pos - 1] in "\n\r":
                start = pos
                break
            # Allow leading whitespace on the DATA; line.
            line_start = text.rfind("\n", 0, pos) + 1
            if text[line_start:pos].strip() == "":
                start = pos
                break
        idx = pos + 5
    if start is None:
        raise ValueError("STEP file has no DATA; section")
    end = upper.find("ENDSEC;", start + 5)
    if end < 0:
        raise ValueError("STEP DATA section is not terminated by ENDSEC;")
    end += len("ENDSEC;")
    return text[start:end]


def canonical_step(text: str) -> str:
    """Canonical STEP form: DATA section of a newline-normalized file.

    FILE_NAME timestamp is stripped first (HEADER). Nothing in DATA is
    rewritten: entity `#N` ids are compared verbatim. A canonicalizer
    that renumbers can mask a topology change; we refuse to do that
    unless two runs of the same binary prove the writer's numbering
    unstable (see the lane report).
    """
    norm = _normalize_newlines(text)
    stripped = strip_file_name_timestamp(norm)
    data = extract_data_section(stripped)
    if not data.endswith("\n"):
        # Trailing newline is not significant inside DATA; keep the
        # extracted bytes as written so a real missing ENDSEC char
        # still diffs. Do not add/remove whitespace.
        pass
    return data


def _extract_json_blob(text: str) -> str:
    """Pull the RESULT JSON out of CLI stdout or a bare object."""
    text = _normalize_newlines(text)
    lines = text.split("\n")
    # Last non-empty line that is a RESULT payload wins (the contract:
    # the last stdout line is always `RESULT {json}`).
    for raw in reversed(lines):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("RESULT "):
            return line[len("RESULT ") :].strip()
        if line.startswith("{"):
            return line
    stripped = text.strip()
    if stripped.startswith("{"):
        return stripped
    raise ValueError("no RESULT JSON object found (expected 'RESULT {...}' or '{...}')")


def parse_result(text: str) -> Dict[str, Any]:
    blob = _extract_json_blob(text)
    obj = json.loads(blob)
    if not isinstance(obj, dict):
        raise ValueError("RESULT JSON is not an object")
    return obj


def is_smooth_key(key: str) -> bool:
    return key.startswith("smooth") or key in _SMOOTH_EXACT


def result_smooth_keys(obj: Dict[str, Any]) -> List[str]:
    return [k for k in obj.keys() if is_smooth_key(k)]


def canonical_result_text(obj: Dict[str, Any]) -> str:
    """Stable text form used for hashing and unified diffs.

    KEYS lists every key in emission order, including `seconds` / path
    keys (order is the contract). VALUES omits `seconds` and path
    values. Python 3.7+ json.loads preserves object key order.
    """
    keys = list(obj.keys())
    values: Dict[str, Any] = {}
    for k, v in obj.items():
        if k in _TIME_KEYS or k in _PATH_KEYS:
            continue
        values[k] = v
    payload = {"keys": keys, "values": values}
    return json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=False) + "\n"


def _unified(a: str, b: str, fa: str, fb: str) -> str:
    return "\n".join(
        difflib.unified_diff(
            a.splitlines(),
            b.splitlines(),
            fromfile=fa,
            tofile=fb,
            lineterm="",
        )
    )


def _print_report(report: Dict[str, Any], as_json: bool, human: str) -> None:
    if as_json:
        json.dump(report, sys.stdout, indent=2, sort_keys=False)
        sys.stdout.write("\n")
        if human:
            sys.stderr.write(human)
            if not human.endswith("\n"):
                sys.stderr.write("\n")
    else:
        sys.stdout.write(human)
        if not human.endswith("\n"):
            sys.stdout.write("\n")


def compare_step(path_a: Path, path_b: Path, as_json: bool) -> int:
    try:
        raw_a = _read_text(path_a)
        raw_b = _read_text(path_b)
        can_a = canonical_step(raw_a)
        can_b = canonical_step(raw_b)
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2

    ha, hb = _sha256_text(can_a), _sha256_text(can_b)
    identical = can_a == can_b
    diff = "" if identical else _unified(can_a, can_b, str(path_a), str(path_b))
    report = {
        "mode": "step",
        "identical": identical,
        "a": {"path": str(path_a), "sha256": ha, "bytes": len(can_a.encode("utf-8"))},
        "b": {"path": str(path_b), "sha256": hb, "bytes": len(can_b.encode("utf-8"))},
        "diff": diff or None,
        "notes": {
            "compared": "STEP DATA section only (HEADER discarded)",
            "stripped": [
                "FILE_NAME field 2 (wall-clock timestamp) — HEADER, run-variant",
            ],
            "entity_ids": "verbatim; no renumbering",
        },
    }
    if identical:
        human = f"IDENTICAL step sha256={ha}\n"
        _print_report(report, as_json, human)
        return 0
    human = (
        f"DIFFER step\n"
        f"  a sha256={ha}\n"
        f"  b sha256={hb}\n"
        f"{diff}\n"
    )
    _print_report(report, as_json, human)
    return 1


def compare_result(path_a: Path, path_b: Path, as_json: bool) -> int:
    try:
        obj_a = parse_result(_read_text(path_a))
        obj_b = parse_result(_read_text(path_b))
        can_a = canonical_result_text(obj_a)
        can_b = canonical_result_text(obj_b)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2

    keys_a, keys_b = list(obj_a.keys()), list(obj_b.keys())
    smooth_a, smooth_b = result_smooth_keys(obj_a), result_smooth_keys(obj_b)
    ha, hb = _sha256_text(can_a), _sha256_text(can_b)

    reasons: List[str] = []
    if smooth_a or smooth_b:
        reasons.append(
            "smooth* keys present with --smooth absent: "
            f"a={smooth_a} b={smooth_b}"
        )
    if keys_a != keys_b:
        if set(keys_a) != set(keys_b):
            reasons.append(
                "RESULT key set differs: "
                f"only_a={ [k for k in keys_a if k not in obj_b] } "
                f"only_b={ [k for k in keys_b if k not in obj_a] }"
            )
        else:
            reasons.append(f"RESULT key ORDER differs: a={keys_a} b={keys_b}")
    if can_a != can_b and "RESULT key" not in " ".join(reasons):
        reasons.append("RESULT values differ (seconds/paths excluded)")

    identical = not reasons and can_a == can_b
    diff = "" if can_a == can_b else _unified(can_a, can_b, str(path_a), str(path_b))
    report = {
        "mode": "result",
        "identical": identical,
        "a": {"path": str(path_a), "sha256": ha, "keys": keys_a, "smooth_keys": smooth_a},
        "b": {"path": str(path_b), "sha256": hb, "keys": keys_b, "smooth_keys": smooth_b},
        "reasons": reasons,
        "diff": diff or None,
        "excluded_values": list(_TIME_KEYS) + list(_PATH_KEYS),
    }
    if identical:
        human = f"IDENTICAL result sha256={ha} keys={','.join(keys_a)}\n"
        _print_report(report, as_json, human)
        return 0
    human_lines = ["DIFFER result"] + [f"  {r}" for r in reasons]
    human_lines.append(f"  a sha256={ha}")
    human_lines.append(f"  b sha256={hb}")
    if diff:
        human_lines.append(diff)
    human = "\n".join(human_lines) + "\n"
    _print_report(report, as_json, human)
    return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description=(
            "G0.1 canonical compare: STEP DATA sections, or RESULT key-set "
            "+ values (seconds and absolute paths excluded)."
        )
    )
    p.add_argument("mode", choices=("step", "result"))
    p.add_argument("a", type=Path, help="left file")
    p.add_argument("b", type=Path, help="right file")
    p.add_argument(
        "--report",
        choices=("json",),
        help="write a JSON report to stdout (human text on stderr if different)",
    )
    args = p.parse_args(argv)
    as_json = args.report == "json"
    if args.mode == "step":
        return compare_step(args.a, args.b, as_json)
    return compare_result(args.a, args.b, as_json)


if __name__ == "__main__":
    sys.exit(main())
