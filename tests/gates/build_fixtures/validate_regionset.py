#!/usr/bin/env python3
"""Stdlib-only RegionSet schema validator (no jsonschema dependency)."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from typing import Any, Dict, List, Set, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_PATH = os.path.normpath(os.path.join(HERE, "..", "regionset.schema.json"))


def load_frozen_schema() -> Dict[str, Any]:
    with open(SCHEMA_PATH, encoding="utf-8") as f:
        schema = json.load(f)
    # Cross-check that this validator still matches the frozen enum sets.
    defs = schema["$defs"]
    frozen_roles = set(defs["Loop"]["properties"]["role"]["enum"])
    frozen_surf = set(defs["Region"]["properties"]["type"]["enum"])
    frozen_origin = set(defs["Region"]["properties"]["origin"]["enum"])
    frozen_reject = set(defs["Region"]["properties"]["reject"]["enum"])
    frozen_built = set(defs["Region"]["properties"]["builtAs"]["enum"])
    if frozen_roles != ROLE_ENUM:
        raise RuntimeError(f"ROLE_ENUM drift vs {SCHEMA_PATH}: {frozen_roles ^ ROLE_ENUM}")
    if frozen_surf != SURF_ENUM:
        raise RuntimeError(f"SURF_ENUM drift vs {SCHEMA_PATH}: {frozen_surf ^ SURF_ENUM}")
    if frozen_origin != ORIGIN_ENUM:
        raise RuntimeError(f"ORIGIN_ENUM drift vs {SCHEMA_PATH}: {frozen_origin ^ ORIGIN_ENUM}")
    if frozen_reject != REJECT_ENUM:
        raise RuntimeError(f"REJECT_ENUM drift vs {SCHEMA_PATH}: {frozen_reject ^ REJECT_ENUM}")
    if frozen_built != BUILT_ENUM:
        raise RuntimeError(f"BUILT_ENUM drift vs {SCHEMA_PATH}: {frozen_built ^ BUILT_ENUM}")
    return schema

VEC3 = {"type": "array", "minItems": 3, "maxItems": 3, "items": "number"}
ROLE_ENUM = {"outer", "inner", "capLow", "capHigh"}
SURF_ENUM = {"plane", "cylinder", "cone", "sphere", "torus"}
ORIGIN_ENUM = {"planeGrow", "cylGrow", "filletStrip"}
REJECT_ENUM = {
    "none", "gaussPlanarity", "vertexResidual", "chordConsistency", "radiusSanity",
    "span", "filletConsensus", "neighborNotAnalytic", "stripWidth", "torusNYI",
    "coneNYI", "sphereNYI", "dirtyComponent", "faceBuildFailed", "chainUnstable",
}
BUILT_ENUM = {"notBuilt", "single", "seamed360", "twoHalves", "explodedToFacets"}


def err(path: str, msg: str, errors: List[str]) -> None:
    errors.append(f"{path}: {msg}")


def check_type(val: Any, want: str, path: str, errors: List[str]) -> None:
    ok = (
        (want == "object" and isinstance(val, dict))
        or (want == "array" and isinstance(val, list))
        or (want == "string" and isinstance(val, str))
        or (want == "integer" and isinstance(val, int) and not isinstance(val, bool))
        or (want == "number" and isinstance(val, (int, float)) and not isinstance(val, bool))
        or (want == "boolean" and isinstance(val, bool))
    )
    if not ok:
        err(path, f"expected {want}, got {type(val).__name__}", errors)


def check_vec3(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "array", path, errors)
    if not isinstance(val, list) or len(val) != 3:
        err(path, "vec3 must have exactly 3 numbers", errors)
        return
    for i, x in enumerate(val):
        check_type(x, "number", f"{path}[{i}]", errors)


def check_ax3(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "object", path, errors)
    if not isinstance(val, dict):
        return
    allowed = {"loc", "dir", "xdir"}
    extra = set(val) - allowed
    for k in extra:
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in allowed:
        if k not in val:
            err(path, f"missing required key '{k}'", errors)
        else:
            check_vec3(val[k], f"{path}.{k}", errors)


def check_loop(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "object", path, errors)
    if not isinstance(val, dict):
        return
    allowed = {"chainIdx", "reversed", "role"}
    for k in set(val) - allowed:
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in allowed:
        if k not in val:
            err(path, f"missing required key '{k}'", errors)
    if "role" in val and val["role"] not in ROLE_ENUM:
        err(f"{path}.role", f"invalid enum {val['role']!r}", errors)
    for k in ("chainIdx", "reversed"):
        if k in val:
            check_type(val[k], "array", f"{path}.{k}", errors)
            if isinstance(val[k], list):
                for i, x in enumerate(val[k]):
                    check_type(x, "integer", f"{path}.{k}[{i}]", errors)


def count_role(loops: List[Any], role: str) -> int:
    n = 0
    for lp in loops:
        if isinstance(lp, dict) and lp.get("role") == role:
            n += 1
    return n


def check_region(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "object", path, errors)
    if not isinstance(val, dict):
        return
    required = [
        "id", "type", "origin", "ax", "radius", "uMin", "uMax", "vMin", "vMax",
        "closed360", "outwardNormal", "tris", "loops", "maxVertexDev", "rmsVertexDev",
        "chordSagitta", "nSides", "dVolPredicted", "maxVertexSnap", "reject", "builtAs",
        "filletNbrA", "filletNbrB",
    ]
    allowed = set(required)
    for k in set(val) - allowed:
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in required:
        if k not in val:
            err(path, f"missing required key '{k}'", errors)
    if val.get("type") not in SURF_ENUM:
        err(f"{path}.type", f"invalid enum {val.get('type')!r}", errors)
    if val.get("origin") not in ORIGIN_ENUM:
        err(f"{path}.origin", f"invalid enum {val.get('origin')!r}", errors)
    if val.get("reject") not in REJECT_ENUM:
        err(f"{path}.reject", f"invalid enum {val.get('reject')!r}", errors)
    if val.get("builtAs") not in BUILT_ENUM:
        err(f"{path}.builtAs", f"invalid enum {val.get('builtAs')!r}", errors)
    if "ax" in val:
        check_ax3(val["ax"], f"{path}.ax", errors)
    if "loops" in val:
        check_type(val["loops"], "array", f"{path}.loops", errors)
        if isinstance(val["loops"], list):
            for i, lp in enumerate(val["loops"]):
                check_loop(lp, f"{path}.loops[{i}]", errors)
            closed360 = val.get("closed360")
            if closed360 is False:
                if count_role(val["loops"], "outer") != 1:
                    err(path, "I7: exactly one outer loop required when closed360=false", errors)
                if count_role(val["loops"], "capLow") != 0 or count_role(val["loops"], "capHigh") != 0:
                    err(path, "I7: capLow/capHigh forbidden when closed360=false", errors)
            if closed360 is True:
                if count_role(val["loops"], "capLow") != 1 or count_role(val["loops"], "capHigh") != 1:
                    err(path, "I7b: exactly one capLow and one capHigh when closed360=true", errors)
                if count_role(val["loops"], "outer") != 0:
                    err(path, "I7b: outer forbidden when closed360=true", errors)


def check_chain(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "object", path, errors)
    if not isinstance(val, dict):
        return
    required = ["regA", "regB", "islandA", "islandB", "tangent", "closedLoop", "meshEdges", "meshVerts"]
    for k in set(val) - set(required):
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in required:
        if k not in val:
            err(path, f"missing required key '{k}'", errors)
    for k in ("regA", "regB", "islandA", "islandB"):
        if k in val:
            check_type(val[k], "integer", f"{path}.{k}", errors)
    for k in ("tangent", "closedLoop"):
        if k in val:
            check_type(val[k], "boolean", f"{path}.{k}", errors)
    for k in ("meshEdges", "meshVerts"):
        if k in val:
            check_type(val[k], "array", f"{path}.{k}", errors)
            if isinstance(val[k], list):
                for i, x in enumerate(val[k]):
                    check_type(x, "integer", f"{path}.{k}[{i}]", errors)


def check_stats(val: Any, path: str, errors: List[str]) -> None:
    check_type(val, "object", path, errors)
    if not isinstance(val, dict):
        return
    required = [
        "planes", "cylinders", "fillets", "rejected", "facetIslands", "facetTriangles",
        "distinctRadii", "maxVertexDev", "maxEdgeTol", "dVolPredicted",
    ]
    for k in set(val) - set(required):
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in required:
        if k not in val:
            err(path, f"missing required key '{k}'", errors)


def validate_regionset(doc: Any, path: str) -> List[str]:
    errors: List[str] = []
    check_type(doc, "object", path, errors)
    if not isinstance(doc, dict):
        return errors
    required = ["compRoot", "regions", "rejected", "chains", "triRegion", "triIsland", "nIslands", "stats"]
    for k in set(doc) - set(required):
        err(f"{path}.{k}", "additional property not allowed", errors)
    for k in required:
        if k not in doc:
            err(path, f"missing required key '{k}'", errors)
    if "regions" in doc and isinstance(doc["regions"], list):
        for i, r in enumerate(doc["regions"]):
            check_region(r, f"{path}.regions[{i}]", errors)
    if "rejected" in doc and isinstance(doc["rejected"], list):
        for i, r in enumerate(doc["rejected"]):
            check_region(r, f"{path}.rejected[{i}]", errors)
    if "chains" in doc and isinstance(doc["chains"], list):
        for i, c in enumerate(doc["chains"]):
            check_chain(c, f"{path}.chains[{i}]", errors)
    if "stats" in doc:
        check_stats(doc["stats"], f"{path}.stats", errors)
    for k in ("triRegion", "triIsland"):
        if k in doc:
            check_type(doc[k], "array", f"{path}.{k}", errors)
            if isinstance(doc[k], list):
                for i, x in enumerate(doc[k]):
                    check_type(x, "integer", f"{path}.{k}[{i}]", errors)
    if isinstance(doc.get("triRegion"), list) and isinstance(doc.get("triIsland"), list):
        if len(doc["triRegion"]) != len(doc["triIsland"]):
            err(path, "triRegion and triIsland must have the same length", errors)
    if "nIslands" in doc:
        check_type(doc["nIslands"], "integer", f"{path}.nIslands", errors)
    if "compRoot" in doc:
        check_type(doc["compRoot"], "integer", f"{path}.compRoot", errors)
    return errors


def main() -> int:
    ap = argparse.ArgumentParser(description="Validate RegionSet JSON against frozen schema rules")
    ap.add_argument("paths", nargs="*", help="files or directories (default: build_fixtures)")
    args = ap.parse_args()
    load_frozen_schema()
    roots = args.paths or [HERE]
    files: List[str] = []
    for root in roots:
        if os.path.isdir(root):
            files.extend(sorted(glob.glob(os.path.join(root, "**", "*.regionset.json"), recursive=True)))
        elif root.endswith(".regionset.json"):
            files.append(root)
    if not files:
        print("No *.regionset.json files found", file=sys.stderr)
        return 1
    failed = 0
    for fp in files:
        with open(fp, encoding="utf-8") as f:
            doc = json.load(f)
        errors = validate_regionset(doc, fp)
        if errors:
            failed += 1
            print(f"FAIL {fp}")
            for e in errors:
                print(f"  {e}")
        else:
            print(f"PASS {fp}")
    print(f"\n{len(files) - failed}/{len(files)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
