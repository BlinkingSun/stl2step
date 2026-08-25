#!/usr/bin/env bash
# Convert every corpus fixture and assert the sidecar's expectedExit /
# expectedSolids / expectedOpenShells. Capture $? immediately.
set -euo pipefail
STL2STEP="$1"
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
done
exit $FAIL
