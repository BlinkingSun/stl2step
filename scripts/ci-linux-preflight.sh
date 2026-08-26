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
set -euo pipefail

HOST="${1:-cnc}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WS="stl2step-ci"

echo "== preflight: sync repo -> $HOST:~/$WS/repo"
ssh "$HOST" mkdir -p "$WS/repo"
# .git rides along: the gate suite's baseline build checks out a pinned
# ancestor commit and needs full history.
# Excludes are ANCHORED (leading /): an unanchored 'build*/' silently drops
# tests/gates/build_fixtures/ and 11 p2buildtest cases vanish — the exact
# historical trap the fixtures dir was renamed to dodge.
rsync -a --delete \
  --exclude '/build*/' --exclude '/_team/' --exclude '/.stl2step-*' \
  "$REPO_ROOT/" "$HOST:$WS/repo/"

ssh "$HOST" bash -s <<'REMOTE'
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

echo "== preflight: ctest"
# LD_LIBRARY_PATH only for the test binaries — exporting it globally makes
# the system cmake/curl load conda's libs (harmless but noisy).
LD_LIBRARY_PATH="$MAMBA_ROOT_PREFIX/envs/occt/lib" \
  nice -n 10 ctest --test-dir build --output-on-failure
REMOTE
