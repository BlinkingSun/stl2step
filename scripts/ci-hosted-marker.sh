#!/usr/bin/env bash
# ci-hosted-marker.sh — mint .ci-local/<sha>.linux.green and .windows.green
# from a green hosted run of THIS sha. Predicates are D-train-push-3 (R1–R8)
# and D-train-push-4 (b >= a/1.5). macOS is never minted here.
#
# Usage: scripts/ci-hosted-marker.sh <run-id> [repo]
# Default repo: BlinkingSun/stl2step
# Every refusal: exit 1, reason on stderr, no marker file written.
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "usage: scripts/ci-hosted-marker.sh <run-id> [repo]" >&2
  exit 1
fi
case "${1}" in
  --*) echo "REFUSE: no override flag may exist" >&2; exit 1 ;;
esac
if [ "$#" -eq 2 ]; then
  case "${2}" in
    --*) echo "REFUSE: no override flag may exist" >&2; exit 1 ;;
  esac
fi

RUN_ID="$1"
REPO="${2:-BlinkingSun/stl2step}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

HEAD_SHA="$(git rev-parse HEAD)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

if ! gh api "repos/${REPO}/actions/runs/${RUN_ID}" > "$tmp/api.json" 2>"$tmp/gh.err"; then
  echo "REFUSE: gh api failed for ${REPO} run ${RUN_ID}" >&2
  cat "$tmp/gh.err" >&2
  exit 1
fi
if ! gh run view "$RUN_ID" -R "$REPO" --json status,conclusion,headSha,jobs > "$tmp/run.json" 2>"$tmp/gh.err"; then
  echo "REFUSE: gh run view failed for ${REPO} run ${RUN_ID}" >&2
  cat "$tmp/gh.err" >&2
  exit 1
fi

# R1–R4. Prints three job ids (linux, windows, macos) or exits 1.
python3 - "$tmp/api.json" "$tmp/run.json" "$HEAD_SHA" "$REPO" <<'PY' > "$tmp/jobs.txt"
import json, re, sys
api_path, run_path, head, repo = sys.argv[1:5]
api = json.load(open(api_path))
run = json.load(open(run_path))

def refuse(msg):
    print(f"REFUSE: {msg}", file=sys.stderr)
    sys.exit(1)

if repo != "BlinkingSun/stl2step":
    refuse(f"R3 repo {repo} is not BlinkingSun/stl2step")
status = api.get("status") or run.get("status")
conclusion = api.get("conclusion") or run.get("conclusion")
if status != "completed" or conclusion != "success":
    refuse(f"R1 status={status} conclusion={conclusion} (need completed/success)")
sha = api.get("head_sha") or ""
if not re.fullmatch(r"[0-9a-f]{40}", sha):
    refuse(f"R2 headSha {sha!r} is not 40 lowercase hex")
if sha != head:
    refuse(f"R2 headSha {sha} != git rev-parse HEAD {head}")
if run.get("headSha") not in (None, sha):
    refuse(f"R2 run headSha {run.get('headSha')} != api head_sha {sha}")
path = api.get("path") or ""
wid = api.get("workflow_id")
if path != ".github/workflows/ci.yml" or wid != 307298879:
    refuse(f"R3 workflow path={path!r} id={wid} (need .github/workflows/ci.yml id 307298879)")
jobs = {}
for job in run.get("jobs") or []:
    jobs[job.get("name")] = job
ids = []
for name in ("linux", "windows", "macos"):
    job = jobs.get(name)
    if job is None:
        refuse(f"R4 job {name} missing")
    conc = job.get("conclusion")
    if conc != "success":
        refuse(f"R4 job {name} conclusion={conc} (need success)")
    jid = job.get("databaseId")
    if jid is None:
        refuse(f"R4 job {name} has no databaseId")
    ids.append(str(jid))
print("\n".join(ids))
PY

LINUX_JID="$(sed -n '1p' "$tmp/jobs.txt")"
WIN_JID="$(sed -n '2p' "$tmp/jobs.txt")"
MAC_JID="$(sed -n '3p' "$tmp/jobs.txt")"

# R7 before any marker is written. (R1–R4 already passed.)
if [ -n "$(git status --porcelain)" ]; then
  echo "REFUSE: R7 git status --porcelain is not empty" >&2
  exit 1
fi

gh run view "$RUN_ID" -R "$REPO" --log --job "$LINUX_JID" > "$tmp/linux.log"
gh run view "$RUN_ID" -R "$REPO" --log --job "$WIN_JID" > "$tmp/windows.log"
gh run view "$RUN_ID" -R "$REPO" --log --job "$MAC_JID" > "$tmp/macos.log"

# R5–R6. stdout is JSON for linux and windows only. macOS is checked and dropped.
PAYLOAD="$(python3 - "$tmp/linux.log" "$tmp/windows.log" "$tmp/macos.log" <<'PY'
import json, re, sys

PREFIX = re.compile(r"^[^\t]*\t[^\t]*\t[0-9]{4}-[0-9]{2}-[0-9]{2}T\S+ ")
SUMMARY = re.compile(r"[0-9]+% tests passed")
TIME = re.compile(
    r"Test #\d+:\s+(gates_full|p1_ac2_gates_baseline)\b.*?"
    r"(Passed|\*\*\*Timeout|\*\*\*Skipped|\*\*\*Failed|\*\*\*Not Run)"
    r"\s+([0-9]+(?:\.[0-9]+)?)\s+sec"
)

def refuse(msg):
    print(f"REFUSE: {msg}", file=sys.stderr)
    sys.exit(1)

def scrape(path, name):
    lines = []
    for raw in open(path, encoding="utf-8", errors="replace"):
        lines.append(PREFIX.sub("", raw.rstrip("\n"), count=1))
    summaries = [ln for ln in lines if SUMMARY.search(ln)]
    if not summaries:
        refuse(f"R5 job {name} has no '[0-9]+% tests passed' line")
    times = {}
    for ln in lines:
        m = TIME.search(ln)
        if m:
            times[m.group(1)] = (m.group(2), m.group(3))
    if "gates_full" not in times or "p1_ac2_gates_baseline" not in times:
        refuse(f"R6 job {name} missing gates_full or p1_ac2_gates_baseline timing ({sorted(times)})")
    g_status, a_s = times["gates_full"]
    p_status, b_s = times["p1_ac2_gates_baseline"]
    if p_status != "Passed":
        refuse(f"R6 job {name} p1_ac2 status={p_status} (need Passed)")
    a = float(a_s)
    b = float(b_s)
    if a <= 0 or not (b >= a / 1.5):
        refuse(
            f"R6 job {name} nest not real: gates_full={a_s}s p1_ac2={b_s}s "
            f"(need p1_ac2 >= gates_full/1.5)"
        )
    return {
        "summary": summaries[-1],
        "a": a_s,
        "b": b_s,
        "ratio": format(b / a, ".6f"),
    }

linux = scrape(sys.argv[1], "linux")
windows = scrape(sys.argv[2], "windows")
scrape(sys.argv[3], "macos")  # required corroboration; never returned, never written
json.dump({"linux": linux, "windows": windows}, sys.stdout)
PY
)"

SHA="$HEAD_SHA"
MARKER_DIR="$REPO_ROOT/.ci-local"
mkdir -p "$MARKER_DIR"
TS="$(date +%Y-%m-%dT%H:%M:%S%z)"

field() {
  python3 -c 'import json,sys; print(json.load(sys.stdin)[sys.argv[1]][sys.argv[2]])' "$1" "$2" <<<"$PAYLOAD"
}

write_one() {
  local os="$1"
  local summary a b ratio jid f
  summary="$(field "$os" summary)"
  a="$(field "$os" a)"
  b="$(field "$os" b)"
  ratio="$(field "$os" ratio)"
  if [ "$os" = "linux" ]; then jid="$LINUX_JID"; else jid="$WIN_JID"; fi
  f="${MARKER_DIR}/${SHA}.${os}.green"
  {
    printf '%s\n' "$summary"
    printf '%s\n' "$TS"
    printf 'hosted-run=%s repo=%s workflow=ci.yml job=%s job-id=%s conclusion=success headSha=%s\n' \
      "$RUN_ID" "$REPO" "$os" "$jid" "$SHA"
    printf 'nest-real: gates_full=%ss p1_ac2_gates_baseline=%ss ratio=%s (D-train-push-3/-4)\n' \
      "$a" "$b" "$ratio"
  } > "$f"
}

# Only these two. There is no macos write.
write_one linux
write_one windows

marker_ok_copy() {
  # marker_ok() from scripts/ci-local-gate.sh — the two regexes, copied.
  local f="$1"
  grep -Eq '[0-9]+% tests passed' "$f" || return 1
  grep -Eq '[0-9]{4}-[0-9]{2}-[0-9]{2}T' "$f" || return 1
  return 0
}

for os in linux windows; do
  f="${MARKER_DIR}/${SHA}.${os}.green"
  if ! marker_ok_copy "$f"; then
    rm -f "${MARKER_DIR}/${SHA}.linux.green" "${MARKER_DIR}/${SHA}.windows.green"
    echo "REFUSE: R8 marker_ok failed for ${f}" >&2
    exit 1
  fi
done

echo "minted ${MARKER_DIR}/${SHA}.linux.green and ${MARKER_DIR}/${SHA}.windows.green"
