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
set -uo pipefail

HOST="${1:-${STL2STEP_WIN_HOST:-}}"
if [ -z "$HOST" ]; then
  echo "usage: scripts/ci-windows-preflight.sh <user@host>" >&2
  echo "   or: STL2STEP_WIN_HOST=<user@host> scripts/ci-windows-preflight.sh" >&2
  exit 1
fi
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "== preflight: sync repo -> $HOST D:\\stl2step-ci\\repo"
ssh "$HOST" "cmd /c (if exist D:\\stl2step-ci\\repo rmdir /s /q D:\\stl2step-ci\\repo) & mkdir D:\\stl2step-ci\\repo"
# tar on both ends. NO --exclude patterns: bsdtar matches them against
# nested path components, which silently drops tests/gates/build_fixtures
# and build_local.cpp (proven 2026-08-26, twice). Instead, enumerate the
# top-level entries and skip build*/_team/.stl2step-* by name.
entries=()
while IFS= read -r e; do
  case "$e" in build*|_team|.stl2step-*) continue ;; esac
  entries+=("$e")
done < <(ls -A "$REPO_ROOT")
# COPYFILE_DISABLE: stop macOS bsdtar emitting AppleDouble ._* companions —
# they extract as literal files on Windows and corrupt *.expected.json globs.
COPYFILE_DISABLE=1 tar -czf - -C "$REPO_ROOT" "${entries[@]}" \
  | ssh "$HOST" "cmd /c tar -xzf - -C D:\\stl2step-ci\\repo"

echo "== preflight: build + ctest (Release)"
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
set PYTHONUTF8=1
cd /d D:\stl2step-ci
"%CMAKE%" -S repo -B build ^
  -DCMAKE_TOOLCHAIN_FILE=D:\vcpkg\scripts\buildsystems\vcpkg.cmake ^
  -DSTL2STEP_BUILD_EXAMPLES=ON -DSTL2STEP_BUILD_TESTS=ON || exit /b 1
"%CMAKE%" --build build --config Release -j 14 || exit /b 1
"%CTEST%" --test-dir build -C Release --output-on-failure
EOF
scp -q /tmp/stl2step-win-preflight.cmd "$HOST:D:/stl2step-ci/preflight.cmd"
ssh "$HOST" "cmd /c D:\\stl2step-ci\\preflight.cmd"
