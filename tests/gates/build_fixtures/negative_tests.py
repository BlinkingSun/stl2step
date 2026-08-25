#!/usr/bin/env python3
"""Eight negative tests from the p2-fixtures adjudication. All must FAIL the checker."""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from check_invariants import validate_file  # noqa: E402
from validate_regionset import validate_regionset  # noqa: E402

GOOD = os.path.join(HERE, "counterbore", "counterbore.regionset.json")
STL = os.path.join(HERE, "counterbore", "counterbore.stl")
EXP = os.path.join(HERE, "counterbore", "counterbore.expected.json")


def stage(rs_src: str, stl_src: str, exp_src: str | None = None):
    tmp = tempfile.mkdtemp()
    shutil.copy(stl_src, os.path.join(tmp, "t.stl"))
    if exp_src and os.path.isfile(exp_src):
        shutil.copy(exp_src, os.path.join(tmp, "t.expected.json"))
    p = os.path.join(tmp, "t.regionset.json")
    rs = json.load(open(rs_src, encoding="utf-8"))
    return tmp, p, rs


def must_fail(label: str, ok: bool, errors: list[str], needle: str | None = None) -> bool:
    if ok:
        print(f"FAIL {label}: checker passed a mutation")
        return False
    if needle and not any(needle.lower() in e.lower() for e in errors):
        print(f"FAIL {label}: no error matching {needle!r}")
        for e in errors[:6]:
            print(f"  {e}")
        return False
    print(f"PASS {label} ({len(errors)} errors)")
    return True


def main() -> int:
    caught = 0

    # 1 dup_cap
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if x["closed360"])
    lo = next(lp for lp in r["loops"] if lp["role"] == "capLow")
    hi = next(lp for lp in r["loops"] if lp["role"] == "capHigh")
    hi["chainIdx"] = list(lo["chainIdx"])
    hi["reversed"] = list(lo["reversed"])
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("dup_cap", ok, err, "cap")
    shutil.rmtree(tmp)

    # 2 drop_chain
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if not x["closed360"])
    r["loops"] = [lp for lp in r["loops"] if lp["role"] == "outer"]
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("drop_chain", ok, err, "coverage")
    shutil.rmtree(tmp)

    # 3 scramble_triregion
    tmp, p, rs = stage(GOOD, STL, EXP)
    rs["triRegion"] = [1 if x == 0 else (0 if x == 1 else x) for x in rs["triRegion"]]
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("scramble_triregion", ok, err, "triregion")
    shutil.rmtree(tmp)

    # 4 drop_inner
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if any(lp["role"] == "inner" for lp in x["loops"]))
    r["loops"] = [lp for lp in r["loops"] if lp["role"] != "inner"]
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("drop_inner", ok, err, "coverage")
    shutil.rmtree(tmp)

    # 5 extra_outer
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if not x["closed360"])
    r["loops"].append({"chainIdx": list(r["loops"][0]["chainIdx"]), "reversed": [0], "role": "outer"})
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("extra_outer", ok, err, "outer")
    shutil.rmtree(tmp)

    # 6 role_swap (capLow -> outer on closed360)
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if x["closed360"])
    for lp in r["loops"]:
        if lp["role"] == "capLow":
            lp["role"] = "outer"
            break
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("role_swap", ok, err, "i7b")
    shutil.rmtree(tmp)

    # 7 cap_order (swap CapLow/CapHigh chains — no longer at vMin/vMax)
    tmp, p, rs = stage(GOOD, STL, EXP)
    r = next(x for x in rs["regions"] if x["closed360"])
    lo = next(lp for lp in r["loops"] if lp["role"] == "capLow")
    hi = next(lp for lp in r["loops"] if lp["role"] == "capHigh")
    lo["chainIdx"], hi["chainIdx"] = hi["chainIdx"], lo["chainIdx"]
    lo["reversed"], hi["reversed"] = hi["reversed"], lo["reversed"]
    json.dump(rs, open(p, "w"))
    ok, err = validate_file(p)
    caught += must_fail("cap_order", ok, err, "vmin")
    shutil.rmtree(tmp)

    # 8 frozen-schema role-enum drift (capMid injected into ROLE_ENUM)
    import validate_regionset as vr
    old = set(vr.ROLE_ENUM)
    vr.ROLE_ENUM = old | {"capMid"}
    try:
        vr.load_frozen_schema()
        print("FAIL schema_enum_drift: load_frozen_schema accepted capMid")
    except RuntimeError as e:
        if "ROLE_ENUM drift" in str(e):
            print("PASS schema_enum_drift")
            caught += 1
        else:
            print(f"FAIL schema_enum_drift: unexpected {e}")
    finally:
        vr.ROLE_ENUM = old

    # bonus: document with capMid also fails validate_regionset
    doc = json.load(open(GOOD, encoding="utf-8"))
    doc["regions"][0]["loops"][0]["role"] = "capMid"
    verr = validate_regionset(doc, "mut.regionset.json")
    if verr:
        print(f"PASS document_capMid ({len(verr)} schema errors)")
    else:
        print("FAIL document_capMid: validator accepted capMid role")

    print(f"\n{caught}/8 negatives caught")
    return 0 if caught == 8 else 1


if __name__ == "__main__":
    sys.exit(main())
