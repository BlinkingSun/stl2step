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
# shellcheck disable=SC2064
trap "rm -f '$LOG'" EXIT

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

JOBS=$(( $(nproc) - 2 )); [ "$JOBS" -lt 1 ] && JOBS=1
echo "== preflight: configure + build (-j$JOBS, nice 10)"
nice -n 10 cmake -S repo -B build -DCMAKE_BUILD_TYPE=Release \
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
