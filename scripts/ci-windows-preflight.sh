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
SYNC_DIRTY="$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no || true)"

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
# Windows preflight keeps the main build tree as a sibling of repo/ (not
# repo/build). Point the default CURRENT_BUILD at the MSVC tree and wrap
# build_baseline.sh so multi-config outputs land where the gate runner looks.
if [[ -n "${MSYSTEM:-}" || "$(uname -s 2>/dev/null)" == MINGW* ]]; then
  rm -rf repo/build
  ln -sf ../build repo/build
  real=repo/tests/gates/baseline/build_baseline.sh
  if [[ ! -f "${real}.real" ]] || ! grep -q build_baseline.sh.real "$real" 2>/dev/null; then
    cp -f "$real" "${real}.real"
  fi
  cat > "$real" <<'WRAPPER'
#!/usr/bin/env bash
set -euo pipefail
dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CMAKE_TOOLCHAIN_FILE="${CMAKE_TOOLCHAIN_FILE:-D:/vcpkg/scripts/buildsystems/vcpkg.cmake}"
export STL2STEP_CURRENT_BUILD="${STL2STEP_CURRENT_BUILD:-$(cd "$dir/../../.." && pwd)/../build}"
args=("$@")
"${dir}/build_baseline.sh.real" "${args[@]}" || true
for cfg in Release Debug; do
  src="${dir}/.build/${cfg}/stl2step.exe"
  if [[ -f "$src" ]]; then
    cp -f "$src" "${dir}/.build/stl2step.exe"
    cp -f "$src" "${dir}/.build/stl2step"
    ver="$("$src" --version)"
    echo ""
    echo "BASELINE_BIN=${src}"
    echo "BASELINE_COMMIT=187ead0d8cf3d3694153cbcff9314d65324fec63"
    echo "BASELINE_GENERATOR=Visual Studio 17 2022"
    echo "BASELINE_BUILD_TYPE=${cfg}"
    echo "BASELINE_VERSION=${ver}"
    exit 0
  fi
done
echo "error: gate baseline binary missing after build" >&2
exit 1
WRAPPER
  chmod +x "$real"
  echo "== preflight: pre-building gate baseline (Visual Studio)"
  export CMAKE_TOOLCHAIN_FILE="D:/vcpkg/scripts/buildsystems/vcpkg.cmake"
  export STL2STEP_CURRENT_BUILD="$(pwd)/build"
  bash "$real" --generator "Visual Studio 17 2022" --build-type Release --current-build "$(pwd)/build"
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
set STL2STEP_CURRENT_BUILD=D:\stl2step-ci\build
set CMAKE_TOOLCHAIN_FILE=D:\vcpkg\scripts\buildsystems\vcpkg.cmake
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
