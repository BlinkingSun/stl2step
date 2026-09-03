#!/usr/bin/env bash
# ci-linux-preflight.sh — run the Linux CI leg on a local machine BEFORE
# pushing, so hosted GitHub minutes are never spent discovering failures
# (private-repo Actions are metered; macOS counts 10x, Windows 2x).
#
# Mirrors .github/workflows/ci.yml (linux job): conda-forge occt=7.9.3
# via micromamba + the machine's system g++/cmake.
#
#   scripts/ci-linux-preflight.sh [host]     # default host: cnc (shop NUC)
#
# Remote workspace: ~/stl2step-ci (repo copy, micromamba, build tree).
# Everything runs nice -n 10 with two cores held back so the machine's
# real workload (cnc-node) keeps priority.
#
# On FULL SUCCESS, writes .ci-local/<sync-sha>.linux.green (suite summary
# + timestamp) for scripts/ci-local-gate.sh. The SHA is the local HEAD
# recorded at sync time. A dirty working tree vs that HEAD refuses the
# marker — the remote run tested uncommitted / untracked files, not git.
set -euo pipefail

HOST="${1:-cnc}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WS="stl2step-ci"

# Record the exact local HEAD the tree is synced from.
SYNC_SHA="$(git -C "$REPO_ROOT" rev-parse HEAD)"
SYNC_DIRTY="$(git -C "$REPO_ROOT" status --porcelain || true)"
SYNC_SHA_SHORT="${SYNC_SHA:0:7}"

# Single-tenant node lock (NODE_LOCK_PROTOCOL=1; scripts/node-lock.sh).
# The vendored twin is local-only; the lock directory lives on the node
# (default ~/.team-node-lock) and the owner record is this preflight
# (local hostname + pid) so a fire-and-forget remote pid cannot look
# stale to the next acquire. Equivalent to:
#   scripts/node-lock.sh acquire --node "$HOST" --name "linux-preflight-$SYNC_SHA_SHORT" --wait 1800
LOCK_NAME="linux-preflight-${SYNC_SHA_SHORT}"
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
  _nl_ssh 'if [ -f "$HOME/.team-node-lock/owner" ]; then tr -d "\r" < "$HOME/.team-node-lock/owner" | head -1; fi' 2>/dev/null || true
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
    if _nl_ssh "mkdir \"\$HOME/.team-node-lock\" 2>/dev/null && printf '%s\n' $(printf '%q' "$owner") > \"\$HOME/.team-node-lock/owner\""; then
      echo "acquired ${HOST} ${LOCK_NAME}"
      return 0
    fi
    owner="$(_nl_read)"
    if [ -n "$owner" ] && _nl_stale "$owner" && [ "$retry" -eq 0 ]; then
      _nl_ssh 'rm -rf "$HOME/.team-node-lock"' || true
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
  _nl_ssh 'rm -rf "$HOME/.team-node-lock"'
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

echo "== preflight: sync repo -> $HOST:~/$WS/repo (HEAD $SYNC_SHA)"
ssh "$HOST" mkdir -p "$WS/repo"
# .git rides along: the gate suite's baseline build checks out a pinned
# ancestor commit and needs full history.
# Excludes are ANCHORED (leading /): an unanchored 'build*/' silently drops
# tests/gates/build_fixtures/ and 11 p2buildtest cases vanish — the exact
# historical trap the fixtures dir was renamed to dodge.
rsync -a --delete \
  --exclude '/build*/' --exclude '/_team/' --exclude '/.stl2step-*' \
  --exclude '/.ci-local/' \
  --exclude '/tests/gates/baseline/.worktree-*' \
  --exclude '/tests/gates/baseline/.build/' \
  "$REPO_ROOT/" "$HOST:$WS/repo/"

LOG="$(mktemp)"

ssh "$HOST" bash -s <<'REMOTE' 2>&1 | tee "$LOG"
set -euo pipefail
cd ~/stl2step-ci

# micromamba + pinned OCCT (idempotent; ~25MB binary package, no sudo)
if [ ! -x bin/micromamba ]; then
  mkdir -p bin
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$HOME/stl2step-ci/mamba"
if [ ! -d "$MAMBA_ROOT_PREFIX/envs/occt" ]; then
  bin/micromamba create -y -n occt -c conda-forge occt=7.9.3
fi
export CMAKE_PREFIX_PATH="$MAMBA_ROOT_PREFIX/envs/occt"

# Compiler parity with hosted CI (2026-09-03): ubuntu-latest builds with its
# distro GCC 13 while this node's default g++ is 15, and a handle initialisation
# that GCC 15 accepted broke the hosted build after a green local round. Build
# with the distro's g++-<LINUX_CI_GCC_MAJOR> (same frontend major as the hosted
# runner, system C++ runtime like the runner) and refuse the marker when the
# major differs. conda-forge's GCC was tried first and links against an older
# C++ runtime than the env's OCCT needs, so the distro package is the right one.
# A build tree configured with another compiler is discarded once.
LINUX_CI_GCC_MAJOR="${LINUX_CI_GCC_MAJOR:-13}"
NODE_CXX="/usr/bin/g++-${LINUX_CI_GCC_MAJOR}"
NODE_CC="/usr/bin/gcc-${LINUX_CI_GCC_MAJOR}"
if [ ! -x "$NODE_CXX" ]; then
  if sudo -n true 2>/dev/null; then
    echo "== preflight: installing g++-${LINUX_CI_GCC_MAJOR} (distro package) to match hosted CI"
    sudo -n apt-get install -y -q "g++-${LINUX_CI_GCC_MAJOR}" "gcc-${LINUX_CI_GCC_MAJOR}"
  else
    echo "REFUSE: ${NODE_CXX} is missing and sudo needs a password; run: sudo apt-get install -y g++-${LINUX_CI_GCC_MAJOR} gcc-${LINUX_CI_GCC_MAJOR}" >&2
    exit 6
  fi
fi
GCC_MAJOR="$("$NODE_CXX" -dumpversion 2>/dev/null | cut -d. -f1)"
echo "== preflight: compiler $("$NODE_CXX" --version 2>/dev/null | head -1) (major ${GCC_MAJOR:-?}; hosted CI major ${LINUX_CI_GCC_MAJOR})"
if [ "$GCC_MAJOR" != "$LINUX_CI_GCC_MAJOR" ]; then
  echo "REFUSE: node compiler major ${GCC_MAJOR:-?} != hosted CI major ${LINUX_CI_GCC_MAJOR}; not building" >&2
  exit 6
fi
export CC="$NODE_CC" CXX="$NODE_CXX"
if [ -f build/CMakeCache.txt ] && ! grep -q "^CMAKE_CXX_COMPILER:FILEPATH=${NODE_CXX}$" build/CMakeCache.txt; then
  echo "== preflight: build tree was configured with another compiler; discarding it once"
  rm -rf build
fi

JOBS=$(( $(nproc) - 2 )); [ "$JOBS" -lt 1 ] && JOBS=1
echo "== preflight: configure + build (-j$JOBS, nice 10)"
# Link against the env's C++ runtime (the one the OCCT libraries were built
# against and the one the rpath loads at run time). Ubuntu's non-default g++-13
# would otherwise resolve -lstdc++ to its own GCC 13 development runtime, which
# predates the CXXABI symbols conda-forge's OCCT needs. CMAKE_BUILD_RPATH pins the
# env lib dir into every build-tree executable, so tools that run DURING the
# build (corpus generation) load OCCT without LD_LIBRARY_PATH (hosted CI exports
# LD_LIBRARY_PATH globally instead; the preflight keeps cmake/curl off conda's libs).
ENV_LIB="$MAMBA_ROOT_PREFIX/envs/occt/lib"
nice -n 10 cmake -S repo -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ENV_LIB}" -DCMAKE_SHARED_LINKER_FLAGS="-L${ENV_LIB}" \
  -DCMAKE_BUILD_RPATH="${ENV_LIB}" \
  -DSTL2STEP_BUILD_EXAMPLES=ON -DSTL2STEP_BUILD_TESTS=ON
nice -n 10 cmake --build build -j "$JOBS"

# Drop Mac-synced baseline artifacts; keep healthy node-built trees.
for wt in repo/tests/gates/baseline/.worktree-*; do
  [ -e "$wt" ] || continue
  head="$(git -C "$wt" rev-parse HEAD 2>/dev/null || true)"
  if [[ -z "$head" || "$head" != 187ead0d8cf3d3694153cbcff9314d65324fec63 ]]; then
    echo "== preflight: removing invalid baseline worktree $wt"
    rm -rf "$wt"
  fi
done
build_cache="repo/tests/gates/baseline/.build/CMakeCache.txt"
if [[ -f "$build_cache" ]] && grep -q '/Users/' "$build_cache" 2>/dev/null; then
  echo "== preflight: removing Mac-synced baseline build cache"
  rm -rf repo/tests/gates/baseline/.build
fi

echo "== preflight: ctest"
# LD_LIBRARY_PATH only for the test binaries — exporting it globally makes
# the system cmake/curl load conda's libs (harmless but noisy).
LD_LIBRARY_PATH="$MAMBA_ROOT_PREFIX/envs/occt/lib" \
  nice -n 10 ctest --test-dir build --output-on-failure
REMOTE

SUMMARY="$(grep -E '[0-9]+% tests passed' "$LOG" | tail -1 || true)"
if [ -z "$SUMMARY" ]; then
  echo "REFUSE marker: no ctest suite summary in linux preflight output" >&2
  exit 1
fi
if [ -n "$SYNC_DIRTY" ]; then
  echo "REFUSE marker: working tree dirty vs HEAD ${SYNC_SHA} at sync time; not writing .ci-local/${SYNC_SHA}.linux.green" >&2
  echo "Dirtiness (git status --porcelain at sync):" >&2
  echo "$SYNC_DIRTY" >&2
  echo "The remote run tested a dirty working tree, not git HEAD. Commit or stash (and disposition untracked files), then re-run." >&2
  exit 1
fi

mkdir -p "$REPO_ROOT/.ci-local"
{
  echo "$SUMMARY"
  date +%Y-%m-%dT%H:%M:%S%z
} > "$REPO_ROOT/.ci-local/${SYNC_SHA}.linux.green"
echo "== preflight: wrote .ci-local/${SYNC_SHA}.linux.green"
echo "$SUMMARY"
