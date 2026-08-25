#!/usr/bin/env python3
"""Throwaway RegionSet conformance checker (lane p1-dump).

Independent of tests/gates/check_regionset.py (p0-ichecker). Validates
emitted JSON against the frozen regionset.schema.json: required keys,
enum spellings, additionalProperties: false, and the I7 / I7b if/then
rules keyed on closed360.

Stdlib only. Usage:

  python3 check_dump.py [dump.json|-] [--schema PATH] [--examples-dir DIR]
  ./stl2step_regiondump tests/cube.stl | python3 check_dump.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable

HERE = Path(__file__).resolve().parent
GATES = HERE.parent
DEFAULT_SCHEMA = GATES / "regionset.schema.json"
DEFAULT_EXAMPLES = GATES

TYPE_ENUM = {"plane", "cylinder", "cone", "sphere", "torus"}
ORIGIN_ENUM = {"planeGrow", "cylGrow", "filletStrip"}
ROLE_ENUM = {"outer", "inner", "capLow", "capHigh"}
REJECT_ENUM = {
    "none",
    "gaussPlanarity",
    "vertexResidual",
    "chordConsistency",
    "radiusSanity",
    "span",
    "filletConsensus",
    "neighborNotAnalytic",
    "stripWidth",
    "torusNYI",
    "coneNYI",
    "sphereNYI",
    "dirtyComponent",
    "faceBuildFailed",
    "chainUnstable",
}
BUILTAS_ENUM = {"notBuilt", "single", "seamed360", "twoHalves", "explodedToFacets"}

REGIONSET_KEYS = [
    "compRoot",
    "regions",
    "rejected",
    "chains",
    "triRegion",
    "triIsland",
    "nIslands",
    "stats",
]
REGION_KEYS = [
    "id",
    "type",
    "origin",
    "ax",
    "radius",
    "uMin",
    "uMax",
    "vMin",
    "vMax",
    "closed360",
    "outwardNormal",
    "tris",
    "loops",
    "maxVertexDev",
    "rmsVertexDev",
    "chordSagitta",
    "nSides",
    "dVolPredicted",
    "maxVertexSnap",
    "reject",
    "builtAs",
    "filletNbrA",
    "filletNbrB",
]
LOOP_KEYS = ["chainIdx", "reversed", "role"]
CHAIN_KEYS = [
    "regA",
    "regB",
    "islandA",
    "islandB",
    "tangent",
    "closedLoop",
    "meshEdges",
    "meshVerts",
]
STATS_KEYS = [
    "planes",
    "cylinders",
    "fillets",
    "rejected",
    "facetIslands",
    "facetTriangles",
    "distinctRadii",
    "maxVertexDev",
    "maxEdgeTol",
    "dVolPredicted",
]
AX_KEYS = ["loc", "dir", "xdir"]

# Wrapper keys the dump tool may add around RegionSet. "stub" is stub-path only.
WRAPPER_KEYS = {
    "stub",
    "file",
    "triangles",
    "vertices",
    "components",
    "dirtySkipped",
    "dumped",
    "comps",
    "index",
    "root",
    "tris",
    "vtx",
    "edges",
    "clean",
    "segmentOk",
    "compVtx",
    "regionSet",
    "regionSets",
}


class Violations:
    def __init__(self) -> None:
        self.items: list[str] = []

    def add(self, path: str, msg: str) -> None:
        self.items.append(f"{path}: {msg}")

    def __bool__(self) -> bool:
        return bool(self.items)


def is_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def is_num(v: Any) -> bool:
    return (isinstance(v, (int, float)) and not isinstance(v, bool))


def expect_obj(v: Any, path: str, viol: Violations) -> dict[str, Any] | None:
    if not isinstance(v, dict):
        viol.add(path, f"expected object, got {type(v).__name__}")
        return None
    return v


def expect_arr(v: Any, path: str, viol: Violations) -> list[Any] | None:
    if not isinstance(v, list):
        viol.add(path, f"expected array, got {type(v).__name__}")
        return None
    return v


def check_keys(
    obj: dict[str, Any],
    required: list[str],
    allowed: set[str],
    path: str,
    viol: Violations,
) -> None:
    keys = set(obj.keys())
    missing = [k for k in required if k not in keys]
    extra = sorted(keys - allowed)
    if missing:
        viol.add(path, f"missing required keys: {missing}")
    if extra:
        viol.add(path, f"additionalProperties: extra keys {extra}")


def check_int_array(v: Any, path: str, viol: Violations) -> None:
    arr = expect_arr(v, path, viol)
    if arr is None:
        return
    for i, item in enumerate(arr):
        if not is_int(item):
            viol.add(f"{path}[{i}]", f"expected integer, got {item!r}")


def check_vec3(v: Any, path: str, viol: Violations) -> None:
    arr = expect_arr(v, path, viol)
    if arr is None:
        return
    if len(arr) != 3:
        viol.add(path, f"vec3 must have 3 items, got {len(arr)}")
    for i, item in enumerate(arr):
        if not is_num(item):
            viol.add(f"{path}[{i}]", f"expected number, got {item!r}")


def check_ax(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, AX_KEYS, set(AX_KEYS), path, viol)
    for k in AX_KEYS:
        if k in obj:
            check_vec3(obj[k], f"{path}.{k}", viol)


def check_loop(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, LOOP_KEYS, set(LOOP_KEYS), path, viol)
    if "chainIdx" in obj:
        check_int_array(obj["chainIdx"], f"{path}.chainIdx", viol)
    if "reversed" in obj:
        check_int_array(obj["reversed"], f"{path}.reversed", viol)
    if "role" in obj:
        if obj["role"] not in ROLE_ENUM:
            viol.add(f"{path}.role", f"invalid enum {obj['role']!r}")
    if "outer" in obj and not isinstance(obj.get("outer"), str):
        # schema has no outer boolean; a leftover bool is an extra key (caught
        # above) — keep the comment so a future bool cannot sneak in silently.
        viol.add(f"{path}.outer", "there is no outer boolean")


def check_region(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, REGION_KEYS, set(REGION_KEYS), path, viol)
    if "id" in obj and not is_int(obj["id"]):
        viol.add(f"{path}.id", f"expected integer, got {obj['id']!r}")
    if "type" in obj and obj["type"] not in TYPE_ENUM:
        viol.add(f"{path}.type", f"invalid enum {obj['type']!r}")
    if "origin" in obj and obj["origin"] not in ORIGIN_ENUM:
        viol.add(f"{path}.origin", f"invalid enum {obj['origin']!r}")
    if "ax" in obj:
        check_ax(obj["ax"], f"{path}.ax", viol)
    for nk in ("radius", "uMin", "uMax", "vMin", "vMax",
               "maxVertexDev", "rmsVertexDev", "chordSagitta",
               "dVolPredicted", "maxVertexSnap"):
        if nk in obj and not is_num(obj[nk]):
            viol.add(f"{path}.{nk}", f"expected number, got {obj[nk]!r}")
    for bk in ("closed360", "outwardNormal"):
        if bk in obj and not isinstance(obj[bk], bool):
            viol.add(f"{path}.{bk}", f"expected boolean, got {obj[bk]!r}")
    if "tris" in obj:
        check_int_array(obj["tris"], f"{path}.tris", viol)
    loops = obj.get("loops")
    arr = expect_arr(loops, f"{path}.loops", viol) if "loops" in obj else None
    if arr is not None:
        for i, lp in enumerate(arr):
            check_loop(lp, f"{path}.loops[{i}]", viol)
    if "nSides" in obj and not is_int(obj["nSides"]):
        viol.add(f"{path}.nSides", f"expected integer, got {obj['nSides']!r}")
    if "reject" in obj and obj["reject"] not in REJECT_ENUM:
        viol.add(f"{path}.reject", f"invalid enum {obj['reject']!r}")
    if "builtAs" in obj and obj["builtAs"] not in BUILTAS_ENUM:
        viol.add(f"{path}.builtAs", f"invalid enum {obj['builtAs']!r}")
    for nk in ("filletNbrA", "filletNbrB"):
        if nk in obj and not is_int(obj[nk]):
            viol.add(f"{path}.{nk}", f"expected integer, got {obj[nk]!r}")

    # I7 / I7b — two if/then rules keyed on closed360.
    if arr is not None and "closed360" in obj:
        roles = [lp.get("role") if isinstance(lp, dict) else None for lp in arr]
        n_outer = roles.count("outer")
        n_low = roles.count("capLow")
        n_high = roles.count("capHigh")
        if obj["closed360"] is False:
            if n_outer != 1:
                viol.add(path, f"I7: closed360=false needs exactly 1 outer, got {n_outer}")
            if n_low != 0:
                viol.add(path, f"I7: closed360=false forbids capLow, got {n_low}")
            if n_high != 0:
                viol.add(path, f"I7: closed360=false forbids capHigh, got {n_high}")
        elif obj["closed360"] is True:
            if n_low != 1:
                viol.add(path, f"I7b: closed360=true needs exactly 1 capLow, got {n_low}")
            if n_high != 1:
                viol.add(path, f"I7b: closed360=true needs exactly 1 capHigh, got {n_high}")
            if n_outer != 0:
                viol.add(path, f"I7b: closed360=true forbids outer, got {n_outer}")


def check_chain(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, CHAIN_KEYS, set(CHAIN_KEYS), path, viol)
    for nk in ("regA", "regB", "islandA", "islandB"):
        if nk in obj and not is_int(obj[nk]):
            viol.add(f"{path}.{nk}", f"expected integer, got {obj[nk]!r}")
    for bk in ("tangent", "closedLoop"):
        if bk in obj and not isinstance(obj[bk], bool):
            viol.add(f"{path}.{bk}", f"expected boolean, got {obj[bk]!r}")
    if "meshEdges" in obj:
        check_int_array(obj["meshEdges"], f"{path}.meshEdges", viol)
    if "meshVerts" in obj:
        check_int_array(obj["meshVerts"], f"{path}.meshVerts", viol)


def check_stats(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, STATS_KEYS, set(STATS_KEYS), path, viol)
    for nk in ("planes", "cylinders", "fillets", "rejected",
               "facetIslands", "facetTriangles", "distinctRadii"):
        if nk in obj and not is_int(obj[nk]):
            viol.add(f"{path}.{nk}", f"expected integer, got {obj[nk]!r}")
    for nk in ("maxVertexDev", "maxEdgeTol", "dVolPredicted"):
        if nk in obj and not is_num(obj[nk]):
            viol.add(f"{path}.{nk}", f"expected number, got {obj[nk]!r}")


def check_regionset(v: Any, path: str, viol: Violations) -> None:
    obj = expect_obj(v, path, viol)
    if obj is None:
        return
    check_keys(obj, REGIONSET_KEYS, set(REGIONSET_KEYS), path, viol)
    if "compRoot" in obj and not is_int(obj["compRoot"]):
        viol.add(f"{path}.compRoot", f"expected integer, got {obj['compRoot']!r}")
    if "nIslands" in obj and not is_int(obj["nIslands"]):
        viol.add(f"{path}.nIslands", f"expected integer, got {obj['nIslands']!r}")
    for name in ("regions", "rejected"):
        if name not in obj:
            continue
        arr = expect_arr(obj[name], f"{path}.{name}", viol)
        if arr is None:
            continue
        for i, r in enumerate(arr):
            check_region(r, f"{path}.{name}[{i}]", viol)
    if "chains" in obj:
        arr = expect_arr(obj["chains"], f"{path}.chains", viol)
        if arr is not None:
            for i, ch in enumerate(arr):
                check_chain(ch, f"{path}.chains[{i}]", viol)
    if "triRegion" in obj:
        check_int_array(obj["triRegion"], f"{path}.triRegion", viol)
    if "triIsland" in obj:
        check_int_array(obj["triIsland"], f"{path}.triIsland", viol)
    if "stats" in obj:
        check_stats(obj["stats"], f"{path}.stats", viol)


def looks_like_regionset(obj: dict[str, Any]) -> bool:
    return "compRoot" in obj and "regions" in obj and "stats" in obj and "chains" in obj


def extract_regionsets(doc: Any, path: str = "$") -> list[tuple[str, dict[str, Any]]]:
    found: list[tuple[str, dict[str, Any]]] = []
    if isinstance(doc, list):
        for i, item in enumerate(doc):
            found.extend(extract_regionsets(item, f"{path}[{i}]"))
        return found
    if not isinstance(doc, dict):
        return found
    if looks_like_regionset(doc):
        found.append((path, doc))
        return found
    if "regionSet" in doc:
        found.extend(extract_regionsets(doc["regionSet"], f"{path}.regionSet"))
    if "regionSets" in doc:
        found.extend(extract_regionsets(doc["regionSets"], f"{path}.regionSets"))
    if "comps" in doc:
        found.extend(extract_regionsets(doc["comps"], f"{path}.comps"))
    if "components" in doc and isinstance(doc["components"], list):
        found.extend(extract_regionsets(doc["components"], f"{path}.components"))
    return found


def collect_keys(obj: Any, bucket: dict[str, set[str]], kind: str) -> None:
    if not isinstance(obj, dict):
        return
    bucket.setdefault(kind, set()).update(obj.keys())
    if kind == "RegionSet":
        for r in obj.get("regions", []) if isinstance(obj.get("regions"), list) else []:
            collect_keys(r, bucket, "Region")
        for r in obj.get("rejected", []) if isinstance(obj.get("rejected"), list) else []:
            collect_keys(r, bucket, "Region")
        for ch in obj.get("chains", []) if isinstance(obj.get("chains"), list) else []:
            collect_keys(ch, bucket, "BoundaryChain")
        collect_keys(obj.get("stats"), bucket, "RefitStats")
    elif kind == "Region":
        collect_keys(obj.get("ax"), bucket, "ax3")
        for lp in obj.get("loops", []) if isinstance(obj.get("loops"), list) else []:
            collect_keys(lp, bucket, "Loop")


def keyset_account(
    dumps: Iterable[dict[str, Any]],
    examples: dict[str, dict[str, Any]],
    wrapper: dict[str, Any] | None,
) -> list[str]:
    lines: list[str] = []
    dump_keys: dict[str, set[str]] = {}
    for rs in dumps:
        collect_keys(rs, dump_keys, "RegionSet")
    for name, ex in examples.items():
        ex_keys: dict[str, set[str]] = {}
        collect_keys(ex, ex_keys, "RegionSet")
        lines.append(f"--- vs {name} ---")
        for kind in ("RegionSet", "Region", "Loop", "BoundaryChain", "RefitStats", "ax3"):
            a = dump_keys.get(kind, set())
            b = ex_keys.get(kind, set())
            if a == b:
                lines.append(f"  {kind} keys: identical {sorted(a)}")
            else:
                lines.append(f"  {kind} keys: dump-only {sorted(a - b)} example-only {sorted(b - a)}")
                lines.append(f"           dump {sorted(a)}")
                lines.append(f"           example {sorted(b)}")
        # value-shape (not keys) so every structural difference is named
        if dumps:
            d0 = list(dumps)[0] if not isinstance(dumps, list) else dumps[0]
            # dumps is already a list
        d0 = next(iter(dumps)) if dumps else {}
        def nreg(o: dict[str, Any], k: str) -> int:
            v = o.get(k, [])
            return len(v) if isinstance(v, list) else -1
        lines.append(
            f"  value-shape: dump regions={nreg(d0,'regions')} rejected={nreg(d0,'rejected')} "
            f"chains={nreg(d0,'chains')} nIslands={d0.get('nIslands')} "
            f"| example regions={nreg(ex,'regions')} rejected={nreg(ex,'rejected')} "
            f"chains={nreg(ex,'chains')} nIslands={ex.get('nIslands')}"
        )
        dump_types = [r.get("type") for r in d0.get("regions", []) if isinstance(r, dict)]
        dump_closed = [r.get("closed360") for r in d0.get("regions", []) if isinstance(r, dict)]
        ex_types = [r.get("type") for r in ex.get("regions", []) if isinstance(r, dict)]
        ex_closed = [r.get("closed360") for r in ex.get("regions", []) if isinstance(r, dict)]
        lines.append(f"  dump region types={dump_types} closed360={dump_closed}")
        lines.append(f"  example region types={ex_types} closed360={ex_closed}")
    if wrapper is not None:
        extra = sorted(set(wrapper.keys()) - set(REGIONSET_KEYS))
        lines.append("--- wrapper envelope (not a RegionSet) ---")
        lines.append(f"  wrapper keys: {sorted(wrapper.keys())}")
        lines.append(f"  keys outside RegionSet schema: {extra}")
        unknown = sorted(set(extra) - WRAPPER_KEYS)
        if unknown:
            lines.append(f"  UNEXPECTED wrapper keys: {unknown}")
        else:
            lines.append("  all wrapper extras are the dump envelope (file/mesh/stub/comps)")
        if wrapper.get("stub") is True:
            lines.append("  stub marker present (strip with del(.stub); real path omits the key)")
    return lines


def load_json(path: str) -> Any:
    if path == "-":
        return json.load(sys.stdin)
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump", nargs="?", default="-",
                    help="dump JSON file, or - for stdin (default)")
    ap.add_argument("--schema", default=str(DEFAULT_SCHEMA),
                    help="frozen regionset.schema.json")
    ap.add_argument("--examples-dir", default=str(DEFAULT_EXAMPLES))
    ap.add_argument("--allow-empty", action="store_true")
    args = ap.parse_args()

    schema_path = Path(args.schema)
    if not schema_path.is_file():
        print(f"FAIL: schema not found: {schema_path}", file=sys.stderr)
        return 2
    # Touch the schema so we fail loudly if it cannot be parsed, even though
    # this checker encodes the rules rather than interpreting JSON Schema.
    with schema_path.open(encoding="utf-8") as f:
        schema = json.load(f)
    if schema.get("additionalProperties") is not False:
        print("FAIL: schema additionalProperties is not false", file=sys.stderr)
        return 2

    examples_dir = Path(args.examples_dir)
    examples = {
        "regionset.example.min.json": json.loads(
            (examples_dir / "regionset.example.min.json").read_text()
        ),
        "regionset.example.closed360.json": json.loads(
            (examples_dir / "regionset.example.closed360.json").read_text()
        ),
    }

    viol = Violations()
    for name, ex in examples.items():
        check_regionset(ex, name, viol)
    if viol:
        print("FAIL: frozen examples themselves violate the checker:", file=sys.stderr)
        for item in viol.items:
            print(f"  {item}", file=sys.stderr)
        return 1

    try:
        doc = load_json(args.dump)
    except json.JSONDecodeError as e:
        print(f"FAIL: dump is not JSON: {e}", file=sys.stderr)
        return 1

    wrapper = doc if isinstance(doc, dict) and not looks_like_regionset(doc) else None
    sets = extract_regionsets(doc)
    if not sets:
        msg = "no RegionSet objects found in dump"
        if args.allow_empty:
            print(f"OK (empty): {msg}")
            return 0
        print(f"FAIL: {msg}", file=sys.stderr)
        return 1

    dump_viol = Violations()
    for path, rs in sets:
        check_regionset(rs, path, dump_viol)
    if dump_viol:
        print(f"FAIL: {len(dump_viol.items)} violation(s)", file=sys.stderr)
        for item in dump_viol.items:
            print(f"  {item}", file=sys.stderr)
        return 1

    print(f"OK: {len(sets)} RegionSet(s) validated against {schema_path.name} "
          f"(required keys, enums, additionalProperties:false, I7/I7b)")
    for path, rs in sets:
        nreg = len(rs.get("regions", []))
        nrej = len(rs.get("rejected", []))
        print(f"  {path}: compRoot={rs.get('compRoot')} regions={nreg} "
              f"rejected={nrej} nIslands={rs.get('nIslands')}")
    if wrapper is not None:
        print(f"  wrapper: triangles={wrapper.get('triangles')} "
              f"vertices={wrapper.get('vertices')} "
              f"components={wrapper.get('components')} "
              f"stub={wrapper.get('stub', '<absent>')}")

    print()
    print("KEY-SET ACCOUNTING")
    for line in keyset_account([rs for _, rs in sets], examples, wrapper):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
