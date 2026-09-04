#!/usr/bin/env python3
"""gate_130 — 1.3.0 battery (B1–B6). Expected-red until Wave 2 detectors land.

For every corpus sidecar with `"battery": "130"`: convert with `--smooth
--no-verify`, assert ok / solids=1 / openShells=0 / reverted=0, built
cylinders == GT, cones == GT (DIAG_130_CENSUS ChamferCone; RESULT has no
cone field), volume delta ≤ 0.01% (census B-Rep volume vs mesh), and when
the sidecar lists intersections, each is represented by its D-130-2 tier:
tier 1 (`cylplane`, `coneplane`, `conecyl-coaxial`, `cylcyl-coaxial`) ships
a CIRCLE / ELLIPSE / LINE; tier 2 (`cylcyl`, the general skew quartic, and
`conecyl`) ships the MESH POLYLINE as ONE edge shared by the two analytic
faces, counted by the 130-BIND census under `edgeClasses.polylineTier2`,
with no `unhandled`, `overTol` or `overCap` on any analytic|analytic edge.

CTest invert: LABELS gates;expected-red + PASS_REGULAR_EXPRESSION
GATE_130_EXPECTED_RED. Flip protocol: remove PASS_REGULAR_EXPRESSION when
the battery is green; keep GATE_130_STRICT=1.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = REPO / "tests" / "corpus"
ENV_STRICT = "GATE_130_STRICT"
MARKER_EXPECTED_RED = "GATE_130_EXPECTED_RED"
MARKER_UNEXPECTED = "GATE_130_UNEXPECTED"
MARKER_FLOORS_MET = "GATE_130_FLOORS_MET"

BATTERY_ORDER = (
    "linkage_bores_chamfer",
    "cross_bores",
    "chamfer_straight_ring",
    "boss_cone_chamfer",
    "counterbore_chamfer",
    "cyl_meets_chamfer",
)
BATTERY_LABEL = {
    "linkage_bores_chamfer": "B1",
    "cross_bores": "B2",
    "chamfer_straight_ring": "B3",
    "boss_cone_chamfer": "B4",
    "counterbore_chamfer": "B5",
    "cyl_meets_chamfer": "B6",
}

DIAG_130_RE = re.compile(
    r"DIAG_130_CENSUS\b[^\n]*\bcone=(\d+)",
    re.MULTILINE,
)
POLYLINE_RE = re.compile(
    r"\bPOLYLINE\s*\(\s*(?:'[^']*'|\"[^\"]*\"|\$)?\s*,?\s*\(([^)]*)\)",
    re.IGNORECASE,
)
COMPOSITE_SEGS_RE = re.compile(
    r"\bCOMPOSITE_CURVE\s*\([^;]*?\(([^)]*)\)",
    re.IGNORECASE | re.DOTALL,
)


def strict_enabled(env: Optional[Mapping[str, str]] = None) -> bool:
    environ: Mapping[str, str] = os.environ if env is None else env
    raw = environ.get(ENV_STRICT)
    if raw is None or str(raw).strip() == "":
        return True
    return str(raw).strip().lower() not in ("0", "false", "off", "no")


def parse_result_line(stdout: str) -> Dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT ") :])
        if line.startswith("RESULT{"):
            return json.loads(line[len("RESULT") :])
    raise RuntimeError("no RESULT line")


def gt_counts(sc: Dict[str, Any]) -> Tuple[int, int]:
    rec = sc.get("recoverable") or []
    cyl = sum(int(p.get("count") or 0) for p in rec if p.get("type") == "cylinder")
    cone = sum(int(p.get("count") or 0) for p in rec if p.get("type") == "cone")
    live = next((row for row in (sc.get("live") or []) if isinstance(row, dict)), {})
    census = live.get("surfaceCensus") or {}
    if cyl == 0:
        cyl = int(census.get("cylinder") or live.get("builtCylindersFloor") or 0)
    if cone == 0:
        cone = int(census.get("cone") or live.get("builtConesFloor") or 0)
    return cyl, cone


def discover_battery(corpus: Path) -> List[Tuple[str, Path, Path, Dict[str, Any]]]:
    found: Dict[str, Tuple[Path, Path, Dict[str, Any]]] = {}
    for sidecar in sorted(corpus.glob("*.expected.json")):
        try:
            doc = json.loads(sidecar.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if str(doc.get("battery") or "") != "130":
            continue
        fid = sidecar.name[: -len(".expected.json")]
        stl = corpus / f"{fid}.stl"
        found[fid] = (stl, sidecar, doc)
    ordered: List[Tuple[str, Path, Path, Dict[str, Any]]] = []
    seen = set()
    for fid in BATTERY_ORDER:
        if fid in found:
            stl, scp, doc = found[fid]
            ordered.append((fid, stl, scp, doc))
            seen.add(fid)
    for fid in sorted(found):
        if fid in seen:
            continue
        stl, scp, doc = found[fid]
        ordered.append((fid, stl, scp, doc))
    return ordered


TIER1_KINDS = ("cylplane", "coneplane", "conecyl-coaxial", "cylcyl-coaxial")
TIER2_KINDS = ("cylcyl", "conecyl")

ENTITY_RE = re.compile(r"^#(\d+)\s*=\s*(.*?);\s*$", re.MULTILINE | re.DOTALL)
REF_RE = re.compile(r"#(\d+)")


def _step_entities(step_text: str) -> Dict[int, str]:
    """#id -> entity body (name + args), one entry per STEP instance."""
    out: Dict[int, str] = {}
    for m in ENTITY_RE.finditer(step_text):
        out[int(m.group(1))] = " ".join(m.group(2).split())
    return out


def _entity_name(body: str) -> str:
    return body.split("(", 1)[0].strip().upper()


def shared_tier2_edges(step_text: str) -> Tuple[List[str], List[str]]:
    """D-130-2 tier 2, structurally on the STEP file.

    Returns (shared, unshared): the EDGE_CURVEs whose 3D curve is a degree-1
    B_SPLINE_CURVE_WITH_KNOTS (the mesh polyline) and that are referenced by
    the loops of exactly two faces whose surfaces are both CYLINDRICAL_SURFACE
    (or CONICAL_SURFACE), versus such polyline edges that only ONE analytic
    face references. A tier-2 seam is honest only when it is one TShape on
    both faces; a polyline each face carries alone is the shell opening.
    """
    ents = _step_entities(step_text)
    name = {k: _entity_name(v) for k, v in ents.items()}
    # face -> surface entity, face -> set of edge_curve ids
    face_surface: Dict[int, int] = {}
    face_edges: Dict[int, set] = {}
    for fid, body in ents.items():
        if name[fid] != "ADVANCED_FACE":
            continue
        refs = [int(x) for x in REF_RE.findall(body)]
        if not refs:
            continue
        face_surface[fid] = refs[-1]  # ADVANCED_FACE('', (bounds), surface, sense)
        edges = set()
        for b in refs[:-1]:
            bb = ents.get(b, "")
            if name.get(b, "") not in ("FACE_BOUND", "FACE_OUTER_BOUND"):
                continue
            for loop in (int(x) for x in REF_RE.findall(bb)):
                lb = ents.get(loop, "")
                if name.get(loop, "") != "EDGE_LOOP":
                    continue
                for oe in (int(x) for x in REF_RE.findall(lb)):
                    ob = ents.get(oe, "")
                    if name.get(oe, "") != "ORIENTED_EDGE":
                        continue
                    ers = [int(x) for x in REF_RE.findall(ob)]
                    if ers:
                        edges.add(ers[-1])
        face_edges[fid] = edges
    analytic_surf = ("CYLINDRICAL_SURFACE", "CONICAL_SURFACE")
    edge_faces: Dict[int, List[int]] = {}
    for fid, edges in face_edges.items():
        if name.get(face_surface.get(fid, -1), "") not in analytic_surf:
            continue
        for e in edges:
            edge_faces.setdefault(e, []).append(fid)
    shared: List[str] = []
    unshared: List[str] = []
    for e, faces in edge_faces.items():
        eb = ents.get(e, "")
        if name.get(e, "") != "EDGE_CURVE":
            continue
        crefs = [int(x) for x in REF_RE.findall(eb)]
        if len(crefs) < 3:
            continue
        curve = crefs[2]
        cb = ents.get(curve, "")
        if name.get(curve, "") != "B_SPLINE_CURVE_WITH_KNOTS":
            continue
        m = re.search(r"B_SPLINE_CURVE_WITH_KNOTS\s*\(\s*'[^']*'\s*,\s*(\d+)", cb)
        if not m or int(m.group(1)) != 1:
            continue
        nseg = len(REF_RE.findall(cb.split(",", 2)[2].split(")")[0])) - 1
        tag = f"#{e} polyline segments={nseg} faces={sorted(faces)}"
        if len(faces) == 2:
            shared.append(tag)
        else:
            unshared.append(tag)
    return shared, unshared


def tier1_curves_present(census: Optional[Dict[str, Any]]) -> bool:
    curves = (census or {}).get("curves") or {}
    return int(curves.get("circle", 0)) + int(curves.get("ellipse", 0)) > 0


def run_convert(binary: Path, stl: Path, step: Path) -> Tuple[int, Dict[str, Any], str, str]:
    env = os.environ.copy()
    env["STL2STEP_DIAG_130"] = "1"
    # Spec fallback: RESULT has no cone field. DIAG_130_CENSUS is the ChamferCone
    # census the engine already prints; P2_DIAG is set as the documented backup.
    env["STL2STEP_P2_DIAG"] = "1"
    proc = subprocess.run(
        [str(binary), str(stl), str(step), "--smooth", "--no-verify", "--quiet"],
        capture_output=True,
        text=True,
        env=env,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    try:
        result = parse_result_line(proc.stdout or "")
    except (RuntimeError, json.JSONDecodeError) as exc:
        result = {"ok": False, "error": f"RESULT parse: {exc}"}
    return proc.returncode, result, proc.stdout or "", combined


def run_census(census_bin: Path, step: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [str(census_bin), str(step)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout or "census failed")
    return json.loads(proc.stdout)


def evaluate_one(
    fid: str,
    stl: Path,
    sc: Dict[str, Any],
    binary: Path,
    census_bin: Optional[Path],
    work: Path,
) -> Dict[str, Any]:
    gt_cyl, gt_cone = gt_counts(sc)
    row: Dict[str, Any] = {
        "id": fid,
        "label": BATTERY_LABEL.get(fid, fid),
        "gt_cyl": gt_cyl,
        "gt_cone": gt_cone,
        "ok": False,
        "solids": -1,
        "openShells": -1,
        "reverted": -1,
        "builtCyl": -1,
        "diagCone": -1,
        "stepCone": -1,
        "volPct": -1.0,
        "curves": "",
        "fails": [],
        "infra": [],
    }
    if not stl.is_file():
        row["infra"].append(f"missing STL {stl}")
        return row

    step = work / f"{fid}.step"
    rc, result, _stdout, combined = run_convert(binary, stl, step)
    row["ok"] = bool(result.get("ok"))
    row["solids"] = int(result.get("solids") or 0)
    row["openShells"] = int(result.get("openShells") or 0)
    row["reverted"] = int(result.get("smoothRevertedComponents") or 0)
    row["builtCyl"] = int(result.get("smoothBuiltCylinders") or 0)

    m = DIAG_130_RE.search(combined)
    if m:
        row["diagCone"] = int(m.group(1))
    row["coneSource"] = "DIAG_130_CENSUS Origin::ChamferCone (RESULT has no cone field)"

    census: Optional[Dict[str, Any]] = None
    if census_bin is not None and step.is_file():
        try:
            census = run_census(census_bin, step)
        except (RuntimeError, json.JSONDecodeError) as exc:
            row["infra"].append(f"census: {exc}")

    if census:
        row["stepCone"] = int((census.get("surfaces") or {}).get("cone") or 0)
        curves = census.get("curves") or {}
        row["curves"] = (
            f"LINE={curves.get('line', 0)} CIRCLE={curves.get('circle', 0)} "
            f"ELLIPSE={curves.get('ellipse', 0)} B_SPLINE={curves.get('bspline', 0)} "
            f"other={curves.get('other', 0)}"
        )
        mesh_vol = float(result.get("meshVolumeMM3") or sc.get("meshVolume") or 0.0)
        step_vol = float(census.get("volume") or 0.0)
        if mesh_vol > 0.0 and step_vol != 0.0:
            row["volPct"] = abs(step_vol - mesh_vol) / abs(mesh_vol) * 100.0
        elif float(result.get("volumeDeltaPct") or -1) >= 0:
            row["volPct"] = float(result["volumeDeltaPct"])
    else:
        vd = float(result.get("volumeDeltaPct") or -1)
        row["volPct"] = vd

    cone_got = row["diagCone"] if row["diagCone"] >= 0 else row["stepCone"]

    if not row["ok"]:
        row["fails"].append(f"ok=false {result.get('error', '')}".strip())
    if row["solids"] != 1:
        row["fails"].append(f"solids={row['solids']} != 1")
    if row["openShells"] != 0:
        row["fails"].append(f"openShells={row['openShells']} != 0")
    if row["reverted"] != 0:
        row["fails"].append(f"smoothRevertedComponents={row['reverted']} != 0")
    if row["builtCyl"] != gt_cyl:
        row["fails"].append(f"smoothBuiltCylinders={row['builtCyl']} != GT {gt_cyl}")
    if cone_got != gt_cone:
        row["fails"].append(f"cones={cone_got} != GT {gt_cone} ({row['coneSource']})")
    if row["volPct"] < 0 or row["volPct"] > 0.01 + 1e-12:
        row["fails"].append(f"volumeDelta={row['volPct']}% > 0.01% (or unmeasured)")

    intersections = sc.get("intersections") or []
    if intersections and step.is_file():
        try:
            text = step.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            row["infra"].append(f"STEP read: {exc}")
            text = ""
        if text:
            # D-130-2 tier rule (SPEC-130-cylcyl addendum 2026-09-03 23:40).
            kinds = [str(x.get("kind") or "") for x in intersections]
            shared, unshared = shared_tier2_edges(text)
            ec = result.get("edgeClasses") if isinstance(result.get("edgeClasses"), dict) else None
            if any(k in TIER2_KINDS for k in kinds):
                if not shared:
                    row["fails"].append(
                        "tier-2 intersection: no polyline edge shared by two analytic faces"
                        + (f" (unshared: {', '.join(unshared)})" if unshared else "")
                    )
                if unshared:
                    row["fails"].append(
                        "tier-2 polyline carried by ONE analytic face only: "
                        + ", ".join(unshared)
                    )
                if ec is None:
                    row["fails"].append(
                        "edgeClasses census absent from RESULT (130-BIND, D-130-2)"
                    )
                else:
                    if int(ec.get("polylineTier2") or 0) < 1:
                        row["fails"].append(
                            f"edgeClasses.polylineTier2={ec.get('polylineTier2')} < 1"
                        )
                    for red in ("unhandled", "overTol", "overCap"):
                        if int(ec.get(red) or 0) != 0:
                            row["fails"].append(f"edgeClasses.{red}={ec.get(red)} != 0")
            if any(k in TIER1_KINDS for k in kinds) and not tier1_curves_present(census):
                row["fails"].append(
                    "tier-1 intersection but no CIRCLE/ELLIPSE in the STEP census"
                )
            if not row["curves"]:
                row["curves"] = "STEP curves not censused"
            kinds_s = ", ".join(
                f"{x.get('a')}∩{x.get('b')}={x.get('kind')}" for x in intersections
            )
            tier2 = f" tier2 shared={len(shared)} unshared={len(unshared)}"
            row["curves"] = f"{row['curves']}  intersections[{kinds_s}]{tier2}"

    return row


def format_table(rows: List[Dict[str, Any]]) -> str:
    headers = (
        "id",
        "ok",
        "solids",
        "open",
        "reverted",
        "cyl got/GT",
        "cone got/GT",
        "vol%",
        "curves / fails",
    )
    lines = ["gate_130 battery (expected-red until Wave 2):", "  " + "  ".join(headers)]
    for r in rows:
        cone_got = r["diagCone"] if r["diagCone"] >= 0 else r["stepCone"]
        status = "FAIL" if (r["fails"] or r["infra"]) else "PASS"
        tail_bits = []
        if r["curves"]:
            tail_bits.append(r["curves"])
        if r["infra"] or r["fails"]:
            tail_bits.append("; ".join(r["infra"] + r["fails"]))
        tail = " | ".join(tail_bits) or "ok"
        lines.append(
            f"  {r['label']:<3} {r['id']:<24} {status:<4} "
            f"ok={int(r['ok'])} sol={r['solids']} open={r['openShells']} "
            f"rev={r['reverted']} cyl={r['builtCyl']}/{r['gt_cyl']} "
            f"cone={cone_got}/{r['gt_cone']} vol={r['volPct']:.4g}  {tail}"
        )
    return "\n".join(lines)


def current_blob(rows: List[Dict[str, Any]]) -> str:
    parts = []
    by_id = {r["id"]: r for r in rows}
    for fid, label in BATTERY_LABEL.items():
        r = by_id.get(fid)
        if not r:
            parts.append(f"{label}=?/?/?")
            continue
        cone_got = r["diagCone"] if r["diagCone"] >= 0 else r["stepCone"]
        parts.append(f"{label}={r['builtCyl']}/{cone_got}/{r['reverted']}")
    return " ".join(parts)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    p.add_argument("--census", type=Path, help="stl2step_census (volume + STEP surfaces/curves)")
    args = p.parse_args(argv)

    is_strict = strict_enabled()
    binary = args.binary
    if not binary.is_file():
        print(f"{MARKER_UNEXPECTED}  missing binary {binary}", file=sys.stderr)
        return 1

    fixtures = discover_battery(args.corpus)
    if not fixtures:
        print(
            f"{MARKER_UNEXPECTED}  no battery=130 sidecars under {args.corpus}",
            file=sys.stderr,
        )
        return 1

    census_bin = args.census if args.census and args.census.is_file() else None
    rows: List[Dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="gate_130_") as tmp:
        work = Path(tmp)
        for fid, stl, _scp, sc in fixtures:
            rows.append(evaluate_one(fid, stl, sc, binary, census_bin, work))

    print(
        f"gate_130  corpus={args.corpus}  binary={binary}  "
        f"census={census_bin or 'none'}  strict={is_strict}  "
        f"cone census=DIAG_130_CENSUS (RESULT has no cone field; "
        f"STL2STEP_P2_DIAG=1 also set)",
        flush=True,
    )
    print(format_table(rows), flush=True)
    print(f"current: {current_blob(rows)}", flush=True)

    infra = [r for r in rows if r["infra"]]
    battery_fail = [r for r in rows if r["fails"] and not r["infra"]]
    if infra:
        named = ", ".join(r["id"] for r in infra)
        print(f"{MARKER_UNEXPECTED}  infra fail(s): {named}", file=sys.stderr, flush=True)
        return 1
    if battery_fail:
        named = ", ".join(f"{r['label']}:{r['id']}" for r in battery_fail)
        print(f"{MARKER_EXPECTED_RED}  unmet battery: {named}", flush=True)
        print(f"gate_130 FAIL  expected-red: {named}", file=sys.stderr, flush=True)
        return 1
    print(f"{MARKER_FLOORS_MET}  battery B1–B6 MET", flush=True)
    print(
        "FLIP PROTOCOL: remove PASS_REGULAR_EXPRESSION from "
        "tests/gates/CMakeLists.txt (gate_130). Keep GATE_130_STRICT=1.",
        flush=True,
    )
    print("gate_130 PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
