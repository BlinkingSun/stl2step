#!/usr/bin/env python3
"""G-LAW — supervised recognition gate for handle-lock (regiondump vs tri-labels).

G-LAW is a **single-fixture supervised recognition gate** (handle-lock only — it is the only mesh
with labels). It is **not** a generalization gate: generalization is guarded by D4 RULE 4.2e
(per-body calibration intervals, the zero-accept negative control, J6 non-regression, ctest, and
decline-path byte-identity). G-LAW **never** gates shipped census.

Scores three ratcheted metrics from a regiondump against tri-labels:
  band recall — true cylinders with exactly one cylinder region at >=95% face coverage
  purity — min majority-label fraction over cylinder-only cylinder regions
  radius exactness — max |R_ours - R_true|/R_true over primary recognized-band regions

RULE 4.2a is a hard sub-check (not ratcheted): banned recognition-threshold literals in src/.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parents[2]
DEFAULT_STL = REPO / "tests" / "corpus" / "handle-lock.stl"
DEFAULT_LABELS = REPO / "tests" / "gates" / "labels" / "handle-lock.tri-labels.json"
DEFAULT_RATCHET = REPO / "tests" / "gates" / "baseline" / "law-ratchet.json"
SOURCE_LABELS = REPO / "_team" / "reports" / "ac2" / "tri-labels.json"
RATCHET_NAME = "law-ratchet.json"
SRC = REPO / "src"
SRC_SUFFIXES = {".cpp", ".hpp", ".h", ".cc", ".cxx"}

RECALL_PURITY_FRAC = 0.95
POST_L2_RECALL_NUM = 15
POST_L2_RECALL_DEN = 15
POST_L2_MIN_PURITY = 1.0
POST_L2_MAX_RADIUS_REL_ERR = 1e-4

# RULE 4.2a — recognition-threshold literals (DECISION tessmath D4 / L1 A16).
BAN_0125 = re.compile(r"\b0\.0125\b")
BAN_261799 = re.compile(r"\b0\.261799\b")
# Bare 15 mm / 15° recognition constants — not 1e-15 epsilons or 0x…C15 hashes.
BAN_BARE_15 = re.compile(r"(?<![0-9.eE+-])15(?:\.0)?(?![0-9eE-])")
LAWBAND_SOURCES = frozenset({"refit_lawband.cpp"})
JUSTIFIED_BOUND = re.compile(
    r"(outer\s+bound|upper\s+bound|lower\s+bound|max\s+span|cap\s+at|clamp)",
    re.IGNORECASE,
)


class GateError(Exception):
    """Hard gate failure."""


@dataclass
class TriLabel:
    tri_id: int
    true_face_id: int
    face_type: str
    radius_mm: float


@dataclass
class Metrics:
    band_recall_num: int
    band_recall_den: int
    min_purity: float
    max_radius_rel_err: float
    per_face: List[Dict[str, Any]]
    per_region: List[Dict[str, Any]]


@dataclass
class ScoreOutcome:
    ok: bool
    metrics: Metrics
    failures: List[str]
    rule42a_hits: List[str]


def load_json(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required file missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, doc: Dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


def load_labels(path: Path) -> Tuple[Dict[int, TriLabel], Dict[int, List[int]]]:
    doc = load_json(path)
    tri_info: Dict[int, TriLabel] = {}
    cyl_faces: Dict[int, List[int]] = defaultdict(list)
    for row in doc.get("triangles") or []:
        tl = TriLabel(
            int(row["tri_id"]),
            int(row["true_face_id"]),
            str(row["face_type"]),
            float(row.get("radius_mm") or 0.0),
        )
        tri_info[tl.tri_id] = tl
        if tl.face_type == "cylinder":
            cyl_faces[tl.true_face_id].append(tl.tri_id)
    return tri_info, dict(cyl_faces)


def region_maps(dump: Dict[str, Any]) -> Tuple[Dict[int, str], Dict[int, float], List[int]]:
    region_type: Dict[int, str] = {}
    region_radius: Dict[int, float] = {}
    for reg in dump.get("regions") or []:
        rid = int(reg["id"])
        region_type[rid] = str(reg.get("type") or "")
        if region_type[rid] == "cylinder":
            region_radius[rid] = float(reg.get("radius") or 0.0)
    tri_region = [int(x) for x in dump.get("triRegion") or []]
    return region_type, region_radius, tri_region


def build_region_tris(tri_region: Sequence[int]) -> Dict[int, List[int]]:
    out: Dict[int, List[int]] = defaultdict(list)
    for tid, rid in enumerate(tri_region):
        if rid >= 0:
            out[rid].append(tid)
    return dict(out)


def face_primary_region(
    cyl_faces: Dict[int, List[int]], tri_region: Sequence[int]
) -> Dict[int, int]:
    primary: Dict[int, int] = {}
    for fid, tris in cyl_faces.items():
        rc = Counter(tri_region[t] for t in tris if 0 <= tri_region[t])
        if rc:
            primary[fid] = int(rc.most_common(1)[0][0])
    return primary


def score_regiondump(
    dump: Dict[str, Any],
    tri_info: Dict[int, TriLabel],
    cyl_faces: Dict[int, List[int]],
) -> Metrics:
    region_type, region_radius, tri_region = region_maps(dump)
    region_tris = build_region_tris(tri_region)
    face_primary = face_primary_region(cyl_faces, tri_region)

    per_face: List[Dict[str, Any]] = []
    recall_num = 0
    for fid in sorted(cyl_faces):
        tris = cyl_faces[fid]
        rc = Counter(tri_region[t] for t in tris if 0 <= tri_region[t])
        n = len(tris)
        winners = [
            rid
            for rid, c in rc.items()
            if n and (c / n) >= RECALL_PURITY_FRAC and region_type.get(rid) == "cylinder"
        ]
        recalled = len(winners) == 1
        if recalled:
            recall_num += 1
        meta = tri_info[tris[0]]
        per_face.append(
            {
                "trueFaceId": fid,
                "radiusMm": meta.radius_mm,
                "nTris": n,
                "regionCounts": {str(k): v for k, v in sorted(rc.items())},
                "winningRegion": winners[0] if recalled else None,
                "recalled": recalled,
            }
        )

    per_region: List[Dict[str, Any]] = []
    min_purity = 1.0
    for rid in sorted(region_tris):
        if region_type.get(rid) != "cylinder":
            continue
        tris = region_tris[rid]
        types = {tri_info[t].face_type for t in tris}
        if types != {"cylinder"}:
            continue
        hist = Counter(tri_info[t].true_face_id for t in tris)
        maj_f, maj_n = hist.most_common(1)[0]
        purity = maj_n / len(tris)
        min_purity = min(min_purity, purity)
        per_region.append(
            {
                "regionId": rid,
                "nTris": len(tris),
                "majorityFaceId": maj_f,
                "purity": purity,
                "radiusMm": region_radius.get(rid, 0.0),
            }
        )

    max_radius_rel_err = 0.0
    for rid, r_ours in region_radius.items():
        tris = region_tris.get(rid, [])
        if not tris:
            continue
        hist = Counter(tri_info[t].true_face_id for t in tris)
        geom_primary = int(hist.most_common(1)[0][0])
        if face_primary.get(geom_primary) != rid:
            continue
        r_true = tri_info[next(t for t in tris if tri_info[t].true_face_id == geom_primary)].radius_mm
        if r_true <= 0:
            continue
        err = abs(r_ours - r_true) / r_true
        max_radius_rel_err = max(max_radius_rel_err, err)

    if not per_region:
        min_purity = 0.0

    return Metrics(
        band_recall_num=recall_num,
        band_recall_den=len(cyl_faces),
        min_purity=min_purity,
        max_radius_rel_err=max_radius_rel_err,
        per_face=per_face,
        per_region=per_region,
    )


def format_face_table(per_face: Sequence[Dict[str, Any]]) -> str:
    lines = [
        "per-true-cylinder-face (tri-labels step face id):",
        f"{'face':>6} {'R_mm':>8} {'nTri':>5} {'recall':>6}  regionCounts",
    ]
    for row in per_face:
        rc = row.get("regionCounts") or {}
        rc_s = " ".join(f"rid{k}:{v}" for k, v in rc.items()) or "-"
        lines.append(
            f"{row['trueFaceId']:>6} {row['radiusMm']:8.4g} {row['nTris']:5d} "
            f"{'PASS' if row['recalled'] else 'FAIL':>6}  {rc_s}"
        )
    return "\n".join(lines)


def scan_rule_42a(src: Path) -> List[str]:
    hits: List[str] = []
    if not src.is_dir():
        return [f"src missing: {src}"]
    for path in sorted(src.rglob("*")):
        if not path.is_file() or path.suffix not in SRC_SUFFIXES:
            continue
        rel = path.relative_to(src.parent)
        lawband_file = path.name in LAWBAND_SOURCES
        for i, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            stripped = line.split("//")[0]
            if not stripped.strip():
                continue
            if BAN_0125.search(stripped):
                hits.append(f"{rel}:{i}: RULE 4.2a banned literal 0.0125: {line.strip()}")
            if BAN_261799.search(stripped):
                hits.append(f"{rel}:{i}: RULE 4.2a banned literal 0.261799: {line.strip()}")
            if lawband_file and BAN_BARE_15.search(stripped):
                if JUSTIFIED_BOUND.search(line):
                    continue
                hits.append(f"{rel}:{i}: RULE 4.2a banned bare 15/15.0: {line.strip()}")
    return hits


def assert_ratchet_authority(ratchet: Dict[str, Any], ratchet_path: Path) -> None:
    if ratchet_path.name != RATCHET_NAME:
        raise GateError(f"floor authority must be {RATCHET_NAME}, opened {ratchet_path}")
    if str(ratchet.get("authority") or "") != RATCHET_NAME:
        raise GateError(f"{ratchet_path}: authority must be {RATCHET_NAME!r}")


def check_against_ratchet(metrics: Metrics, ratchet: Dict[str, Any]) -> List[str]:
    fails: List[str] = []
    floor_num = int(ratchet["bandRecallNum"])
    floor_den = int(ratchet["bandRecallDen"])
    if metrics.band_recall_num < floor_num:
        fails.append(
            f"band recall {metrics.band_recall_num}/{metrics.band_recall_den} "
            f"< floor {floor_num}/{floor_den}"
        )
    if metrics.min_purity + 1e-12 < float(ratchet["minPurity"]):
        fails.append(
            f"min purity {metrics.min_purity:.6f} < floor {float(ratchet['minPurity']):.6f}"
        )
    if metrics.max_radius_rel_err > float(ratchet["maxRadiusRelErr"]) + max(1e-9, float(ratchet["maxRadiusRelErr"]) * 1e-6):
        fails.append(
            f"max radius rel err {metrics.max_radius_rel_err:.6f} "
            f"> floor {float(ratchet['maxRadiusRelErr']):.6f}"
        )
    return fails


def run_regiondump(dump_bin: Path, stl: Path, out_path: Path) -> Dict[str, Any]:
    proc = subprocess.run(
        [
            str(dump_bin),
            str(stl),
            "--component",
            "0",
            "--bare",
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not out_path.is_file():
        raise GateError(
            f"regiondump failed rc={proc.returncode}: {proc.stderr or proc.stdout}"
        )
    return load_json(out_path)


def evaluate(
    dump: Dict[str, Any],
    labels_path: Path,
    ratchet_path: Path,
    *,
    src: Path = SRC,
) -> ScoreOutcome:
    tri_info, cyl_faces = load_labels(labels_path)
    ratchet = load_json(ratchet_path)
    assert_ratchet_authority(ratchet, ratchet_path)
    metrics = score_regiondump(dump, tri_info, cyl_faces)
    failures = check_against_ratchet(metrics, ratchet)
    rule_hits = scan_rule_42a(src)
    if rule_hits:
        failures.extend(f"RULE 4.2a: {h}" for h in rule_hits)
    return ScoreOutcome(not failures, metrics, failures, rule_hits)


def ratchet_up(ratchet_path: Path) -> Dict[str, Any]:
    doc = load_json(ratchet_path)
    old = dict(doc)
    doc["bandRecallNum"] = POST_L2_RECALL_NUM
    doc["bandRecallDen"] = POST_L2_RECALL_DEN
    doc["minPurity"] = POST_L2_MIN_PURITY
    doc["maxRadiusRelErr"] = POST_L2_MAX_RADIUS_REL_ERR
    write_json(ratchet_path, doc)
    print(
        f"ratchet-up: recall {old['bandRecallNum']}/{old['bandRecallDen']} -> "
        f"{doc['bandRecallNum']}/{doc['bandRecallDen']}, "
        f"purity {old['minPurity']} -> {doc['minPurity']}, "
        f"radius {old['maxRadiusRelErr']} -> {doc['maxRadiusRelErr']} ({ratchet_path})"
    )
    return doc


def make_regression_dump(base: Dict[str, Any]) -> Dict[str, Any]:
    """Synthetic worse dump: lower recall, purity, higher radius error."""
    dump = copy.deepcopy(base)
    # Merge two recalled cylinder faces into one polluted region.
    dump["regions"] = copy.deepcopy(base["regions"])
    # pick rid=12 (R5 F4) and rid=15 (R5 F17) — both pure cylinder regions
    r12 = next(r for r in dump["regions"] if r["id"] == 12)
    r15 = next(r for r in dump["regions"] if r["id"] == 15)
    merged_tris = list(r12["tris"]) + list(r15["tris"])
    r12["tris"] = merged_tris
    r12["radius"] = 99.0
    dump["regions"] = [r for r in dump["regions"] if r["id"] != 15]
    tri_region = list(dump["triRegion"])
    for i, rid in enumerate(tri_region):
        if rid == 15:
            tri_region[i] = 12
    dump["triRegion"] = tri_region
    return dump


def make_improvement_dump(base: Dict[str, Any], tri_info: Dict[int, TriLabel], cyl_faces: Dict[int, List[int]]) -> Dict[str, Any]:
    """Synthetic 15/15 dump: one cylinder region per true cylinder face."""
    dump = copy.deepcopy(base)
    tri_region = list(dump["triRegion"])
    new_regions = [r for r in dump["regions"] if r["type"] != "cylinder"]
    next_id = max(int(r["id"]) for r in dump["regions"]) + 1
    for fid in sorted(cyl_faces):
        tris = cyl_faces[fid]
        meta = tri_info[tris[0]]
        rid = next_id
        next_id += 1
        for t in tris:
            tri_region[t] = rid
        new_regions.append(
            {
                "id": rid,
                "type": "cylinder",
                "origin": "cylGrow",
                "ax": {
                    "loc": [0.0, 0.0, 0.0],
                    "dir": [0.0, 1.0, 0.0],
                    "xdir": [1.0, 0.0, 0.0],
                },
                "radius": meta.radius_mm,
                "uMin": 0.0,
                "uMax": 1.0,
                "vMin": 0.0,
                "vMax": 1.0,
                "closed360": False,
                "outwardNormal": True,
                "tris": tris,
                "loops": [],
                "maxVertexDev": 0.0,
                "rmsVertexDev": 0.0,
                "chordSagitta": 0.01,
                "nSides": 8,
                "dVolPredicted": 0.0,
                "maxVertexSnap": 0.0,
                "reject": "none",
                "builtAs": "notBuilt",
                "filletNbrA": -1,
                "filletNbrB": -1,
            }
        )
    dump["regions"] = new_regions
    dump["triRegion"] = tri_region
    return dump


def verify_fixture(labels_path: Path, source_path: Path) -> List[str]:
    fails: List[str] = []
    if not labels_path.is_file():
        return [f"fixture missing: {labels_path}"]
    if not source_path.is_file():
        return [f"source labels missing: {source_path}"]
    a = labels_path.read_bytes()
    b = source_path.read_bytes()
    if a != b:
        fails.append(
            f"fixture not byte-identical to {source_path} "
            f"(fixture {len(a)} bytes, source {len(b)} bytes)"
        )
    doc = json.loads(a.decode("utf-8"))
    if int(doc.get("n_tri") or 0) != 908:
        fails.append(f"n_tri={doc.get('n_tri')} expected 908")
    if int(doc.get("n_ambiguous") or 0) != 0:
        fails.append(f"n_ambiguous={doc.get('n_ambiguous')} expected 0")
    if int(doc.get("n_faces") or 0) != 28:
        fails.append(f"n_faces={doc.get('n_faces')} expected 28")
    if len(doc.get("triangles") or []) != 908:
        fails.append(f"triangles rows={len(doc.get('triangles') or [])} expected 908")
    return fails


def _self_test() -> int:
    fails = 0

    def check(cond: bool, msg: str) -> None:
        nonlocal fails
        if not cond:
            print(f"SELFTEST FAIL: {msg}", file=sys.stderr)
            fails += 1
        else:
            print(f"SELFTEST PASS: {msg}")

    fixture_fails = verify_fixture(DEFAULT_LABELS, SOURCE_LABELS)
    check(fixture_fails == [], f"fixture integrity ({fixture_fails or 'ok'})")

    tri_info, cyl_faces = load_labels(DEFAULT_LABELS)
    check(len(cyl_faces) == 15, f"15 true cylinder faces (got {len(cyl_faces)})")

    dump_bin = REPO / "build" / "stl2step_regiondump"
    if dump_bin.is_file() and DEFAULT_STL.is_file():
        with tempfile.TemporaryDirectory(prefix="glaw_self_") as td:
            base = run_regiondump(dump_bin, DEFAULT_STL, Path(td) / "base.json")
        metrics = score_regiondump(base, tri_info, cyl_faces)
        check(
            metrics.band_recall_num == 10 and metrics.band_recall_den == 15,
            f"baseline recall 10/15 (got {metrics.band_recall_num}/{metrics.band_recall_den})",
        )
        check(
            abs(metrics.min_purity - 0.857) < 0.001,
            f"baseline min purity ~0.857 (got {metrics.min_purity:.3f})",
        )
        check(
            abs(metrics.max_radius_rel_err - 0.349785) < 0.00001,
            f"baseline radius err ~0.349785 (got {metrics.max_radius_rel_err:.6f})",
        )
        outcome = evaluate(base, DEFAULT_LABELS, DEFAULT_RATCHET)
        check(outcome.ok, f"baseline green on ratchet floors ({outcome.failures})")

        reg_dump = make_regression_dump(base)
        reg = evaluate(reg_dump, DEFAULT_LABELS, DEFAULT_RATCHET)
        check(not reg.ok, "synthetic regression FAIL")
        check(bool(reg.failures), "regression has readable failures")
        print("--- injected regression transcript ---")
        print(format_face_table(reg.metrics.per_face))

        imp_dump = make_improvement_dump(base, tri_info, cyl_faces)
        imp_metrics = score_regiondump(imp_dump, tri_info, cyl_faces)
        check(
            imp_metrics.band_recall_num == 15,
            f"improvement recall 15/15 (got {imp_metrics.band_recall_num})",
        )
        check(
            imp_metrics.min_purity >= 1.0 - 1e-12,
            f"improvement purity 1.0 (got {imp_metrics.min_purity})",
        )
        check(
            imp_metrics.max_radius_rel_err < 1e-4,
            f"improvement radius <1e-4 (got {imp_metrics.max_radius_rel_err})",
        )

        with tempfile.TemporaryDirectory(prefix="glaw_ratchet_") as td:
            ratchet_copy = Path(td) / RATCHET_NAME
            ratchet_copy.write_text(DEFAULT_RATCHET.read_text(encoding="utf-8"), encoding="utf-8")
            ratchet_up(ratchet_copy)
            doc = load_json(ratchet_copy)
            check(doc["bandRecallNum"] == 15, "ratchet-up recall num 15")
            check(float(doc["minPurity"]) == 1.0, "ratchet-up purity 1.0")
            check(float(doc["maxRadiusRelErr"]) == POST_L2_MAX_RADIUS_REL_ERR, "ratchet-up radius 1e-4")
            failures_tight = check_against_ratchet(imp_metrics, doc)
            check(not failures_tight, f"improvement PASS on ratchet-up floors ({failures_tight})")
    else:
        print("SELFTEST SKIP: build/stl2step_regiondump or handle-lock.stl missing", file=sys.stderr)

    # RULE 4.2a API
    with tempfile.TemporaryDirectory(prefix="glaw_42a_") as td:
        src = Path(td) / "src"
        src.mkdir()
        clean = src / "ok.cpp"
        clean.write_text(
            "const double span = 15.0; // outer bound for chart span\n",
            encoding="utf-8",
        )
        bad = src / "refit_lawband.cpp"
        bad.write_text("const double d = 0.0125;\nconst double r = 15.0;\n", encoding="utf-8")
        hits = scan_rule_42a(src)
        check(any("refit_lawband.cpp" in h and "0.0125" in h for h in hits), "RULE 4.2a catches planted 0.0125")
        check(any("refit_lawband.cpp" in h and "15.0" in h for h in hits), "RULE 4.2a catches planted 15.0 in lawband TU")
        check(not any("ok.cpp" in h for h in hits), "RULE 4.2a allows justified bounded 15.0")

    hits_live = scan_rule_42a(SRC)
    check(hits_live == [], f"src/ RULE 4.2a clean (hits={hits_live})")

    return 1 if fails else 0


def run_gate(args: argparse.Namespace) -> int:
    labels_path = args.labels.resolve()
    ratchet_path = args.ratchet.resolve()
    stl = args.stl.resolve()
    dump_bin = args.dump.resolve() if args.dump else None
    dump_path = args.regiondump.resolve() if args.regiondump else None

    fixture_fails = verify_fixture(labels_path, SOURCE_LABELS)
    if fixture_fails:
        for f in fixture_fails:
            print(f"G-LAW FAIL: {f}", file=sys.stderr)
        return 1

    if dump_path and dump_path.is_file():
        dump = load_json(dump_path)
    else:
        if dump_bin is None or not dump_bin.is_file():
            raise GateError("stl2step_regiondump missing (pass --dump)")
        if not stl.is_file():
            raise GateError(f"handle-lock STL missing: {stl}")
        with tempfile.TemporaryDirectory(prefix="glaw_gate_") as td:
            out = Path(td) / "regionset.json"
            with ThreadPoolExecutor(max_workers=2) as pool:
                dump_fut = pool.submit(run_regiondump, dump_bin, stl, out)
                rule_fut = pool.submit(scan_rule_42a, args.src.resolve())
                dump = dump_fut.result()
                rule_hits_precheck = rule_fut.result()
            tri_info, cyl_faces = load_labels(labels_path)
            ratchet = load_json(ratchet_path)
            assert_ratchet_authority(ratchet, ratchet_path)
            metrics = score_regiondump(dump, tri_info, cyl_faces)
            failures = check_against_ratchet(metrics, ratchet)
            if rule_hits_precheck:
                failures.extend(f"RULE 4.2a: {h}" for h in rule_hits_precheck)
            outcome = ScoreOutcome(not failures, metrics, failures, rule_hits_precheck)
            print(
                f"G-LAW metrics: recall={outcome.metrics.band_recall_num}/"
                f"{outcome.metrics.band_recall_den} "
                f"purity={outcome.metrics.min_purity:.6f} "
                f"radiusErr={outcome.metrics.max_radius_rel_err:.6f}"
            )
            print(format_face_table(outcome.metrics.per_face))
            if outcome.failures:
                print("G-LAW FAIL", file=sys.stderr)
                for f in outcome.failures:
                    print(f"  {f}", file=sys.stderr)
                return 1
            print("G-LAW PASS")
            return 0

    outcome = evaluate(dump, labels_path, ratchet_path, src=args.src.resolve())

    print(
        f"G-LAW metrics: recall={outcome.metrics.band_recall_num}/"
        f"{outcome.metrics.band_recall_den} "
        f"purity={outcome.metrics.min_purity:.6f} "
        f"radiusErr={outcome.metrics.max_radius_rel_err:.6f}"
    )
    print(format_face_table(outcome.metrics.per_face))

    if outcome.failures:
        print("G-LAW FAIL", file=sys.stderr)
        for f in outcome.failures:
            print(f"  {f}", file=sys.stderr)
        return 1

    print("G-LAW PASS")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dump", type=Path, help="stl2step_regiondump binary")
    p.add_argument("--regiondump", type=Path, help="precomputed bare regionset JSON")
    p.add_argument("--stl", type=Path, default=DEFAULT_STL)
    p.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    p.add_argument("--ratchet", type=Path, default=DEFAULT_RATCHET)
    p.add_argument("--src", type=Path, default=SRC)
    p.add_argument("--threads", type=int, default=2, help="parallel sub-check workers")
    p.add_argument("--ratchet-up", action="store_true", help="tighten floors after green (orchestrator)")
    p.add_argument("--self-test", action="store_true", help="exercise gate API (uses engine when built)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            return _self_test()
        if args.ratchet_up:
            ratchet_up(args.ratchet.resolve())
            return 0
        return run_gate(args)
    except GateError as exc:
        print(f"G-LAW FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
