#!/usr/bin/env bash
# Convert every corpus fixture and assert the sidecar's expectedExit /
# expectedSolids / expectedOpenShells (smooth off). When smoothExpectedExit is
# set, also run --smooth and check warning/volume/solids never-get-worse bounds.
set -euo pipefail
STL2STEP="$1"
CENSUS="${CENSUS:-}"
CORPUS="${CORPUS:-tests/corpus}"
FAIL=0
for stl in "$CORPUS"/*.stl; do
  [[ -f "$stl" ]] || continue
  id=$(basename "$stl" .stl)
  sidecar="$CORPUS/${id}.expected.json"
  out="/tmp/${id}.step"
  set +e
  "$STL2STEP" "$stl" "$out" --quiet --no-verify >/tmp/"${id}"_result.json 2>/tmp/"${id}"_stderr.txt
  rc=$?
  set -e
  if [[ ! -f "$sidecar" ]]; then
    echo "FAIL $id: missing sidecar"
    FAIL=1
    continue
  fi
  python3 - "$id" "$rc" "$sidecar" /tmp/"${id}"_result.json <<'PY' || { FAIL=1; continue; }
import json, sys
id, rc, sc_path, res_path = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
sc = json.load(open(sc_path))
text = open(res_path).read().strip()
if text.startswith("RESULT "):
    text = text[len("RESULT "):]
r = json.loads(text.splitlines()[-1])
exp_rc = int(sc.get("expectedExit", 0))
exp_sol = int(sc.get("expectedSolids", 1))
exp_open = int(sc.get("expectedOpenShells", 0))
ok = True
if rc != exp_rc:
    print(f"FAIL {id}: exit {rc} != expectedExit {exp_rc}")
    ok = False
if not r.get("ok", False):
    print(f"FAIL {id}: ok!=true")
    ok = False
if int(r.get("solids", -1)) != exp_sol:
    print(f"FAIL {id}: solids {r.get('solids')} != expectedSolids {exp_sol}")
    ok = False
if int(r.get("openShells", -1)) != exp_open:
    print(f"FAIL {id}: openShells {r.get('openShells')} != expectedOpenShells {exp_open}")
    ok = False
if ok:
    print(f"OK {id} exit={rc} solids={exp_sol} openShells={exp_open}")
sys.exit(0 if ok else 1)
PY

  smooth_exit=$(python3 -c "import json; print(json.load(open('$sidecar')).get('smoothExpectedExit', -1))")
  if [[ "$smooth_exit" != "-1" ]]; then
    sout="/tmp/${id}_smooth.step"
    set +e
    "$STL2STEP" "$stl" "$sout" --quiet --smooth >/tmp/"${id}"_smooth_result.json 2>/tmp/"${id}"_smooth_stderr.txt
    src=$?
    set -e
    python3 - "$id" "$src" "$sidecar" /tmp/"${id}"_smooth_result.json <<'PY' || { FAIL=1; continue; }
import json, sys
id, rc, sc_path, res_path = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
sc = json.load(open(sc_path))
text = open(res_path).read().strip()
if text.startswith("RESULT "):
    text = text[len("RESULT "):]
r = json.loads(text.splitlines()[-1])
exp_rc = int(sc["smoothExpectedExit"])
exp_sol = int(sc.get("expectedSolids", 1))
exp_open = int(sc.get("expectedOpenShells", 0))
warn_max = int(sc.get("smoothWarningCountMax", -1))
vol_max = float(sc.get("smoothVolumeDeltaPctMax", -1))
rec_cyl = int(sc.get("smoothRecognisedCylinders", -1))
ok = True
if rc != exp_rc:
    print(f"FAIL {id} smooth: exit {rc} != smoothExpectedExit {exp_rc}")
    ok = False
if not r.get("ok", False):
    print(f"FAIL {id} smooth: ok!=true")
    ok = False
nw = len(r.get("warnings") or [])
if warn_max >= 0 and nw > warn_max:
    print(f"FAIL {id} smooth: warnings {nw} > smoothWarningCountMax {warn_max}")
    ok = False
if vol_max >= 0:
    vd = float(r.get("volumeDeltaPct", -1))
    if vd < 0 or vd > vol_max + 1e-9:
        print(f"FAIL {id} smooth: volumeDeltaPct {vd} > smoothVolumeDeltaPctMax {vol_max}")
        ok = False
if int(r.get("solids", -1)) != exp_sol:
    print(f"FAIL {id} smooth: solids {r.get('solids')} != expectedSolids {exp_sol}")
    ok = False
if int(r.get("openShells", -1)) != exp_open:
    print(f"FAIL {id} smooth: openShells {r.get('openShells')} != expectedOpenShells {exp_open}")
    ok = False
if not r.get("watertight", False):
    print(f"FAIL {id} smooth: watertight!=true")
    ok = False
# File-truth census: built cylinders must meet every live[] floor (today 0).
live = sc.get("live") or []
floor = max((int(row.get("builtCylindersFloor", 0)) for row in live), default=0)
smooth_cyl = int(r.get("smoothCylinders", 0))
if rec_cyl >= 0 and smooth_cyl != rec_cyl:
    print(f"NOTE {id} smooth: smoothCylinders recognised={smooth_cyl} (sidecar {rec_cyl})")
if ok:
    print(f"OK {id} smooth exit={rc} warnings={nw} solids={exp_sol} recognisedCyl={smooth_cyl} floor={floor}")
sys.exit(0 if ok else 1)
PY
    if [[ -n "$CENSUS" && -x "$CENSUS" ]]; then
      python3 - "$id" "$sidecar" "$sout" "$CENSUS" <<'PY' || { FAIL=1; continue; }
import json, subprocess, sys
id, sc_path, step, census = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
sc = json.load(open(sc_path))
out = subprocess.check_output([census, step], text=True)
c = json.loads(out)
built_cyl = int(c["surfaces"]["cylinder"])
live = sc.get("live") or []
floor = max((int(row.get("builtCylindersFloor", 0)) for row in live), default=0)
if built_cyl < floor:
    print(f"FAIL {id} census: built cylinders {built_cyl} < floor {floor}")
    sys.exit(1)
if not c.get("valid") or not c.get("closed"):
    print(f"FAIL {id} census: valid={c.get('valid')} closed={c.get('closed')}")
    sys.exit(1)
print(f"OK {id} census cylinders={built_cyl} planes={c['surfaces']['plane']} valid closed")
PY
    fi
  fi
done
exit $FAIL
