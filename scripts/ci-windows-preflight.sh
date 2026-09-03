#!/usr/bin/env bash
# ci-windows-preflight.sh — run the Windows CI leg on the shop PC BEFORE
# pushing, so hosted GitHub minutes are never spent discovering failures
# (Windows counts 2x on metered private-repo Actions; the vcpkg OCCT
# build alone is ~3.5h on a cold hosted runner).
#
#   scripts/ci-windows-preflight.sh <user@host>
#   STL2STEP_WIN_HOST=<user@host> scripts/ci-windows-preflight.sh
#
# Host is required ($1 or STL2STEP_WIN_HOST). No baked-in user, hostname,
# or LAN address.
#
# One-time prep on the PC (already done 2026-08-26): VS 2022 Build Tools
# with MSVC + SDK + bundled CMake, and D:\stl2step-ci\install-occt.cmd
# (vcpkg + opencascade x64-windows into D:\vcpkg).
#
# Repo lands wholesale in D:\stl2step-ci\repo each run (no rsync on
# Windows); the build tree D:\stl2step-ci\build persists for incremental
# rebuilds. Mirrors ci.yml (windows job): vcpkg toolchain, VS multi-config
# generator, Release, ctest -C Release.
#
# On FULL SUCCESS, writes .ci-local/<sync-sha>.windows.green (suite summary
# + timestamp) for scripts/ci-local-gate.sh. The SHA is the local HEAD
# recorded at sync time. A dirty working tree vs that HEAD refuses the
# marker — the remote run tested uncommitted / untracked files, not git.
set -uo pipefail

HOST="${1:-${STL2STEP_WIN_HOST:-}}"
if [ -z "$HOST" ]; then
  echo "usage: scripts/ci-windows-preflight.sh <user@host>" >&2
  echo "   or: STL2STEP_WIN_HOST=<user@host> scripts/ci-windows-preflight.sh" >&2
  exit 1
fi
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Record the exact local HEAD the tree is synced from.
SYNC_SHA="$(git -C "$REPO_ROOT" rev-parse HEAD)"
SYNC_DIRTY="$(git -C "$REPO_ROOT" status --porcelain || true)"
SYNC_SHA_SHORT="${SYNC_SHA:0:7}"

# Single-tenant node lock (NODE_LOCK_PROTOCOL=1; scripts/node-lock.sh).
# --windows lock dir D:\team-node-lock. Owner is this preflight (local
# hostname + pid). Equivalent to:
#   scripts/node-lock.sh acquire --node "$HOST" --windows --name "windows-preflight-$SYNC_SHA_SHORT" --wait 1800
LOCK_NAME="windows-preflight-${SYNC_SHA_SHORT}"
LOCK_TTL=7200
LOCK_WAIT=1800
LOCK_HELD=0
LOG=""

_nl_ssh() {
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$HOST" "$@"
}
_nl_parse() {
  printf '%s' "$1" | tr ' ' '\n' | sed -n "s/^${2}=//p" | head -1
}
_nl_owner() {
  printf 'name=%s host=%s pid=%s started=%s ttl=%s' \
    "$LOCK_NAME" "$(hostname 2>/dev/null || echo unknown)" "$$" "$(date +%s)" "$LOCK_TTL"
}
_nl_read() {
  _nl_ssh 'cmd /c "if exist D:\team-node-lock\owner type D:\team-node-lock\owner"' 2>/dev/null | tr -d '\r' | head -1 || true
}
_nl_stale() {
  local line="$1" started ttl age host pid now
  started="$(_nl_parse "$line" started)"
  ttl="$(_nl_parse "$line" ttl)"
  host="$(_nl_parse "$line" host)"
  pid="$(_nl_parse "$line" pid)"
  [ -n "$started" ] || return 0
  [ -n "$ttl" ] || ttl=7200
  now=$(date +%s)
  age=$((now - started))
  if [ "$age" -ge "$ttl" ]; then
    echo "stale-lock-broken name=$(_nl_parse "$line" name) age=${age}" >&2
    return 0
  fi
  if [ "$host" = "$(hostname 2>/dev/null || echo unknown)" ] && [ -n "$pid" ] && [ "$pid" != 0 ] && ! kill -0 "$pid" 2>/dev/null; then
    echo "stale-lock-broken name=$(_nl_parse "$line" name) age=${age}" >&2
    return 0
  fi
  return 1
}
_nl_acquire() {
  local start now elapsed owner retry=0 sleep_for
  start=$(date +%s)
  while :; do
    owner="$(_nl_owner)"
    if _nl_ssh "cmd /c \"mkdir D:\\team-node-lock 2>nul && echo ${owner}>D:\\team-node-lock\\owner\""; then
      echo "acquired ${HOST} ${LOCK_NAME}"
      return 0
    fi
    owner="$(_nl_read)"
    if [ -n "$owner" ] && _nl_stale "$owner" && [ "$retry" -eq 0 ]; then
      _nl_ssh 'cmd /c "rmdir /s /q D:\team-node-lock"' || true
      retry=1
      continue
    fi
    now=$(date +%s)
    elapsed=$((now - start))
    if [ "$LOCK_WAIT" -gt 0 ] && [ "$elapsed" -lt "$LOCK_WAIT" ]; then
      sleep_for=$((LOCK_WAIT - elapsed))
      [ "$sleep_for" -gt 15 ] && sleep_for=15
      sleep "$sleep_for"
      retry=0
      continue
    fi
    echo "busy: ${owner}"
    return 3
  done
}
_nl_release() {
  local owner name
  owner="$(_nl_read)"
  [ -n "$owner" ] || return 0
  name="$(_nl_parse "$owner" name)"
  [ "$name" = "$LOCK_NAME" ] || return 4
  _nl_ssh 'cmd /c "rmdir /s /q D:\team-node-lock"'
}
_nl_status() {
  local owner
  owner="$(_nl_read)"
  if [ -z "$owner" ]; then
    echo free
    return 0
  fi
  echo "$owner"
  return 3
}
_nl_cleanup() {
  if [ "$LOCK_HELD" = 1 ]; then
    _nl_release >/dev/null 2>&1 || true
    LOCK_HELD=0
  fi
  [ -n "${LOG:-}" ] && rm -f "$LOG"
}

if ! _nl_acquire; then
  echo "REFUSE: node busy: $(_nl_status || true)" >&2
  exit 5
fi
LOCK_HELD=1
trap _nl_cleanup EXIT

echo "== preflight: sync repo -> $HOST D:\\stl2step-ci\\repo (HEAD $SYNC_SHA)"
ssh "$HOST" "cmd /c (if exist D:\\stl2step-ci\\repo rmdir /s /q D:\\stl2step-ci\\repo) & mkdir D:\\stl2step-ci\\repo"
# tar on both ends. NO --exclude patterns: bsdtar matches them against
# nested path components, which silently drops tests/gates/build_fixtures
# and build_local.cpp (proven 2026-08-26, twice). Instead, enumerate the
# top-level entries and skip build*/_team/.stl2step-* by name.
entries=()
while IFS= read -r e; do
  case "$e" in build*|_team|.stl2step-*|.ci-local) continue ;; esac
  if [[ "$e" == "tests" ]]; then
    while IFS= read -r rel; do
      entries+=("$rel")
    done < <(cd "$REPO_ROOT" && find tests \
      \( -path 'tests/gates/baseline/.worktree-*' -o -path 'tests/gates/baseline/.worktree-*/*' \
          -o -path 'tests/gates/baseline/.build' -o -path 'tests/gates/baseline/.build/*' \) -prune \
      -o -print)
    continue
  fi
  entries+=("$e")
done < <(ls -A "$REPO_ROOT")
# COPYFILE_DISABLE: stop macOS bsdtar emitting AppleDouble ._* companions —
# they extract as literal files on Windows and corrupt *.expected.json globs.
COPYFILE_DISABLE=1 tar -czf - -C "$REPO_ROOT" "${entries[@]}" \
  | ssh "$HOST" "cmd /c tar -xzf - -C D:\\stl2step-ci\\repo"

echo "== preflight: build + ctest (Release)"
cat > /tmp/stl2step-baseline-cleanup.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
COMMIT=187ead0d8cf3d3694153cbcff9314d65324fec63
for wt in repo/tests/gates/baseline/.worktree-*; do
  [ -e "$wt" ] || continue
  head="$(git -C "$wt" rev-parse HEAD 2>/dev/null || true)"
  if [[ -z "$head" || "$head" != "$COMMIT" ]]; then
    echo "== preflight: removing invalid baseline worktree $wt"
    rm -rf "$wt"
  fi
done
if [[ -f repo/tests/gates/baseline/.build/CMakeCache.txt ]] \
    && grep -q '/Users/' repo/tests/gates/baseline/.build/CMakeCache.txt 2>/dev/null; then
  echo "== preflight: removing Mac-synced baseline build cache"
  rm -rf repo/tests/gates/baseline/.build
fi
EOF
cat > /tmp/stl2step-win-preflight.cmd <<'EOF'
@echo off
set BT=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools
set CMAKE=%BT%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe
set CTEST=%BT%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\ctest.exe
call "%BT%\VC\Auxiliary\Build\vcvars64.bat" >nul
rem Git Bash + unix tools (bash for build_baseline.sh / run_engine_check.sh,
rem diff for the corpus determinism check). Appended, so nothing is shadowed.
set PATH=%PATH%;C:\Program Files\Git\bin;C:\Program Files\Git\usr\bin
rem Windows Python defaults stdout to cp1252; the gate runner prints U+2264.
set "PYTHONUTF8=1"
cd /d D:\stl2step-ci
"%CMAKE%" -S repo -B build ^
  -DCMAKE_TOOLCHAIN_FILE=D:\vcpkg\scripts\buildsystems\vcpkg.cmake ^
  -DSTL2STEP_BUILD_EXAMPLES=ON -DSTL2STEP_BUILD_TESTS=ON || exit /b 1
"%CMAKE%" --build build --config Release -j 14 || exit /b 1
bash baseline-cleanup.sh || exit /b 1
"%CTEST%" --test-dir build -C Release --output-on-failure
EOF
scp -q /tmp/stl2step-baseline-cleanup.sh "$HOST:D:/stl2step-ci/baseline-cleanup.sh"
scp -q /tmp/stl2step-win-preflight.cmd "$HOST:D:/stl2step-ci/preflight.cmd"

LOG="$(mktemp)"
if ! ssh "$HOST" "cmd /c D:\\stl2step-ci\\preflight.cmd" 2>&1 | tee "$LOG"; then
  rm -f "$LOG"
  exit 1
fi

SUMMARY="$(tr -d '\r' < "$LOG" | grep -E '[0-9]+% tests passed' | tail -1 || true)"
rm -f "$LOG"
if [ -z "$SUMMARY" ]; then
  echo "REFUSE marker: no ctest suite summary in windows preflight output" >&2
  exit 1
fi
if [ -n "$SYNC_DIRTY" ]; then
  echo "REFUSE marker: working tree dirty vs HEAD ${SYNC_SHA} at sync time; not writing .ci-local/${SYNC_SHA}.windows.green" >&2
  echo "Dirtiness (git status --porcelain at sync):" >&2
  echo "$SYNC_DIRTY" >&2
  echo "The remote run tested a dirty working tree, not git HEAD. Commit or stash (and disposition untracked files), then re-run." >&2
  exit 1
fi

mkdir -p "$REPO_ROOT/.ci-local"
{
  echo "$SUMMARY"
  date +%Y-%m-%dT%H:%M:%S%z
} > "$REPO_ROOT/.ci-local/${SYNC_SHA}.windows.green"
echo "== preflight: wrote .ci-local/${SYNC_SHA}.windows.green"
echo "$SUMMARY"
