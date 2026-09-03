#!/usr/bin/env bash
# ci-local-gate.sh — decide whether HEAD may be pushed.
#
# After the 2026-08-30 posix_compat.hpp incident (git add -u dropped a NEW
# header; working-tree / tar-synced preflights passed; every clean Actions
# checkout failed), a push that touches compiled paths must already have
# three-platform local CI green markers. Private-repo Actions minutes are
# the stake: never discover an OS-support failure on the hosted runner.
#
# Usage:
#   scripts/ci-local-gate.sh              # classify + enforce (pre-push)
#   scripts/ci-local-gate.sh --fresh-clone  # produce .ci-local/<sha>.macos.green
#
# Classification of `git diff --name-only <upstream>...HEAD`:
#   BUILD-AFFECTING  src/ include/ tests/ CMakeLists.txt CMakePresets.json cmake/
#                    plus any other non-docs, non-scripts/*.sh path (e.g. examples/)
#   DOCS-ONLY        docs/  *.md  images/  image files  .github/  .gitignore
#   SCRIPT           scripts/*.sh and scripts/hooks/* — gate-relevant, not compiled
#                    (macOS fresh-clone only; no linux/windows markers)
#   Empty range (HEAD already matches upstream) is NOT docs-only — treated as
#   BUILD-AFFECTING so a standalone gate of HEAD still requires markers.
#
# BUILD-AFFECTING → require .ci-local/<HEAD>.{macos,linux,windows}.green
#   (suite summary line + timestamp). Missing any → REFUSE with the exact
#   command to produce it. When all three exist, always run the FRESH-CLONE
#   check (file:// clone + cmake + ctest) and rewrite macos.green.
# DOCS-ONLY → pass (no markers, no fresh-clone).
# SCRIPT → fresh-clone only.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MODE="gate"
if [ "${1:-}" = "--fresh-clone" ]; then
  MODE="fresh-clone"
elif [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  sed -n '2,30p' "$0"
  exit 0
elif [ -n "${1:-}" ]; then
  echo "usage: scripts/ci-local-gate.sh [--fresh-clone]" >&2
  exit 2
fi

SHA="$(git rev-parse HEAD)"
MARKER_DIR="$REPO_ROOT/.ci-local"
mkdir -p "$MARKER_DIR"

upstream=""
if git rev-parse --abbrev-ref '@{upstream}' >/dev/null 2>&1; then
  upstream="$(git rev-parse --abbrev-ref '@{upstream}')"
elif git rev-parse --verify origin/main >/dev/null 2>&1; then
  upstream="origin/main"
fi

# --- classification --------------------------------------------------------

is_build_path() {
  case "$1" in
    src/*|include/*|tests/*|cmake/*|examples/*) return 0 ;;
    CMakeLists.txt|CMakePresets.json) return 0 ;;
    */CMakeLists.txt) return 0 ;;
  esac
  return 1
}

is_docs_path() {
  case "$1" in
    docs/*|.github/*|images/*|.gitignore) return 0 ;;
    *.md|*.png|*.jpg|*.jpeg|*.gif|*.svg|*.webp) return 0 ;;
  esac
  return 1
}

is_script_path() {
  case "$1" in
    scripts/*.sh|scripts/hooks/*) return 0 ;;
  esac
  return 1
}

classify() {
  local build=0 script=0 docs=0 other=0
  local f
  if [ -z "$upstream" ]; then
    echo "BUILD-AFFECTING"
    return
  fi
  local files
  files="$(git diff --name-only "${upstream}...HEAD")"
  if [ -z "$files" ]; then
    # Nothing unpushed: cannot prove docs-only. Gate HEAD as compiled.
    echo "BUILD-AFFECTING"
    return
  fi
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    if is_build_path "$f"; then
      build=1
    elif is_script_path "$f"; then
      script=1
    elif is_docs_path "$f"; then
      docs=1
    else
      other=1
    fi
  done <<< "$files"
  if [ "$build" -eq 1 ] || [ "$other" -eq 1 ]; then
    echo "BUILD-AFFECTING"
  elif [ "$script" -eq 1 ]; then
    echo "SCRIPT"
  else
    echo "DOCS-ONLY"
  fi
}

# --- markers ---------------------------------------------------------------

marker_ok() {
  local f="$1"
  [ -f "$f" ] || return 1
  grep -Eq '[0-9]+% tests passed' "$f" || return 1
  grep -Eq '[0-9]{4}-[0-9]{2}-[0-9]{2}T' "$f" || return 1
  return 0
}

win_cmd() {
  if [ -n "${STL2STEP_WIN_HOST:-}" ]; then
    echo "scripts/ci-windows-preflight.sh ${STL2STEP_WIN_HOST}"
  else
    echo "scripts/ci-windows-preflight.sh <user@host>   # or STL2STEP_WIN_HOST=<user@host> scripts/ci-windows-preflight.sh"
  fi
}

refuse_missing() {
  local missing=0
  echo "REFUSE: build-affecting push of ${SHA} requires three-platform local CI markers." >&2
  echo >&2
  if ! marker_ok "${MARKER_DIR}/${SHA}.macos.green"; then
    echo "  missing .ci-local/${SHA}.macos.green" >&2
    echo "    produce: scripts/ci-local-gate.sh --fresh-clone" >&2
    missing=1
  fi
  if ! marker_ok "${MARKER_DIR}/${SHA}.linux.green"; then
    echo "  missing .ci-local/${SHA}.linux.green" >&2
    echo "    produce: scripts/ci-linux-preflight.sh" >&2
    missing=1
  fi
  if ! marker_ok "${MARKER_DIR}/${SHA}.windows.green"; then
    echo "  missing .ci-local/${SHA}.windows.green" >&2
    echo "    produce: $(win_cmd)" >&2
    missing=1
  fi
  echo >&2
  echo "Working-tree / tar-synced preflights are not enough (untracked files ride along)." >&2
  echo "Fresh-clone (--fresh-clone) tests git's contents — what Actions checks out." >&2
  exit 1
}

# --- fresh-clone (catches the git add -u / untracked-header class) ---------

# Reuse a fresh-clone macOS marker for THIS HEAD when it is younger than
# STL2STEP_FRESH_MARKER_MAX_MIN minutes (default 120): the marker is only ever
# written by a completed fresh clone of the identical commit, so re-cloning
# minutes later proves nothing new and costs ~15 min per push. Older, missing
# or malformed markers fall through to a real fresh clone. `--fresh-clone`
# itself always clones.
fresh_clone_or_reuse() {
  local f="${MARKER_DIR}/${SHA}.macos.green" max="${STL2STEP_FRESH_MARKER_MAX_MIN:-120}" stamp now then age
  if marker_ok "$f"; then
    stamp="$(grep -Eo '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}[+-][0-9]{4}' "$f" | tail -n1)"
    if [ -n "$stamp" ]; then
      then="$(date -j -f '%Y-%m-%dT%H:%M:%S%z' "$stamp" +%s 2>/dev/null || date -d "$stamp" +%s 2>/dev/null || echo 0)"
      now="$(date +%s)"
      if [ "$then" -gt 0 ]; then
        age=$(( (now - then) / 60 ))
        if [ "$age" -ge 0 ] && [ "$age" -lt "$max" ]; then
          echo "== ci-local-gate: reusing fresh-clone marker for ${SHA} (age ${age} min < ${max})"
          return 0
        fi
      fi
    fi
  fi
  fresh_clone
}

fresh_clone() {
  local tmp dest build log summary
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/stl2step-fresh-XXXXXX")"
  dest="$tmp/repo"
  build="$tmp/build"
  log="$tmp/ctest.log"
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp'" EXIT

  echo "== ci-local-gate: FRESH-CLONE file://${REPO_ROOT} → $dest"
  git clone "file://${REPO_ROOT}" "$dest"

  echo "== ci-local-gate: cmake + build + ctest (Release, examples+tests)"
  cmake -S "$dest" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSTL2STEP_BUILD_EXAMPLES=ON -DSTL2STEP_BUILD_TESTS=ON
  cmake --build "$build" -j
  # pipefail: ctest failure aborts; summary line is captured from the log.
  ctest --test-dir "$build" --output-on-failure | tee "$log"
  summary="$(grep -E '[0-9]+% tests passed' "$log" | tail -1 || true)"
  if [ -z "$summary" ]; then
    echo "REFUSE: fresh-clone ctest produced no suite summary line" >&2
    exit 1
  fi

  mkdir -p "$MARKER_DIR"
  {
    echo "$summary"
    date +%Y-%m-%dT%H:%M:%S%z
  } > "${MARKER_DIR}/${SHA}.macos.green"
  echo "== ci-local-gate: wrote .ci-local/${SHA}.macos.green"
  echo "$summary"
}

# --- main ------------------------------------------------------------------

if [ "$MODE" = "fresh-clone" ]; then
  fresh_clone
  exit 0
fi

CLASS="$(classify)"
echo "== ci-local-gate: HEAD=${SHA} vs ${upstream:-<no-upstream>}  class=${CLASS}"

case "$CLASS" in
  DOCS-ONLY)
    echo "== ci-local-gate: PASS (docs-only; no markers required)"
    exit 0
    ;;
  SCRIPT)
    echo "== ci-local-gate: scripts/*.sh only — macOS fresh-clone (no linux/windows markers)"
    fresh_clone_or_reuse
    echo "== ci-local-gate: PASS (script / fresh-clone)"
    exit 0
    ;;
  BUILD-AFFECTING)
    if ! marker_ok "${MARKER_DIR}/${SHA}.macos.green" \
      || ! marker_ok "${MARKER_DIR}/${SHA}.linux.green" \
      || ! marker_ok "${MARKER_DIR}/${SHA}.windows.green"; then
      refuse_missing
    fi
    echo "== ci-local-gate: three-platform markers present; running FRESH-CLONE"
    fresh_clone_or_reuse
    echo "== ci-local-gate: PASS"
    exit 0
    ;;
  *)
    echo "REFUSE: unknown class ${CLASS}" >&2
    exit 1
    ;;
esac
