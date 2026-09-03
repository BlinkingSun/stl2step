#!/usr/bin/env bash
# build_baseline.sh — reproducible stl2step 1.0.0 (commit 187ead0) CLI.
#
# Checks out 187ead0 without touching the current worktree, configures it
# with the same CMake generator / build type as the current tree (or
# CMake's platform default + Release), and builds tests/gates/baseline/.build/stl2step.
#
# Never commit the checkout, the build tree, or the binary.
# Idempotent: a second run reuses the worktree and an incremental build.
set -euo pipefail

COMMIT_SHORT="187ead0"
COMMIT="187ead0d8cf3d3694153cbcff9314d65324fec63"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
SRC_DIR="${STL2STEP_BASELINE_SRC:-$SCRIPT_DIR/.worktree-187ead0}"
BUILD_DIR="${STL2STEP_BASELINE_BUILD:-$SCRIPT_DIR/.build}"
CURRENT_BUILD="${STL2STEP_CURRENT_BUILD:-$REPO_ROOT/build}"
FORCE=0
JOBS=""
GENERATOR=""
BUILD_TYPE=""

usage() {
    cat <<EOF
usage: $0 [--force] [--jobs N] [--generator G] [--build-type T] [--current-build DIR]

Build the released 1.0.0 engine (git commit ${COMMIT_SHORT}) into:
  ${BUILD_DIR}/stl2step

--current-build DIR   read CMAKE_GENERATOR / CMAKE_BUILD_TYPE from this
                      cache (default: <repo>/build). Used so the baseline
                      matches the tree G0.1 will compare against.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --force) FORCE=1; shift ;;
        --jobs) JOBS="${2:?}"; shift 2 ;;
        --generator) GENERATOR="${2:?}"; shift 2 ;;
        --build-type) BUILD_TYPE="${2:?}"; shift 2 ;;
        --current-build) CURRENT_BUILD="${2:?}"; shift 2 ;;
        *) echo "error: unknown argument $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$JOBS" ]]; then
    if command -v sysctl >/dev/null 2>&1; then
        JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
    else
        JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
    fi
fi

cache_get() {
    # cache_get <CMakeCache.txt> <VAR>
    local cache="$1" var="$2"
    [[ -f "$cache" ]] || return 1
    # CMAKE_GENERATOR is INTERNAL; CMAKE_BUILD_TYPE is STRING; compiler FILEPATH.
    sed -n "s/^${var}:[A-Z]*=//p" "$cache" | head -n1
}

is_vs_generator() {
    [[ "${1:-}" == "Visual Studio"* ]]
}

emit_path() {
    local p="$1"
    if command -v cygpath >/dev/null 2>&1; then
        cygpath -m "$p"
    else
        printf '%s\n' "$p"
    fi
}

# Match the current tree's generator / build type when the caller did not
# pick one. If the cache and flag are both empty, omit -G so CMake uses
# $CMAKE_GENERATOR or its platform default.
if [[ -z "$GENERATOR" ]]; then
    GENERATOR="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" CMAKE_GENERATOR || true)"
fi
if [[ -z "$BUILD_TYPE" ]]; then
    BUILD_TYPE="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" CMAKE_BUILD_TYPE || true)"
fi
if [[ -z "$BUILD_TYPE" ]]; then
    BUILD_TYPE="Release"
fi

CXX_COMPILER="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" CMAKE_CXX_COMPILER || true)"
C_COMPILER="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" CMAKE_C_COMPILER || true)"
OCCT_DIR="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" OpenCASCADE_DIR || true)"
TOOLCHAIN="$(cache_get "$CURRENT_BUILD/CMakeCache.txt" CMAKE_TOOLCHAIN_FILE || true)"
if [[ -z "$TOOLCHAIN" && -n "${CMAKE_TOOLCHAIN_FILE:-}" ]]; then
    TOOLCHAIN="$CMAKE_TOOLCHAIN_FILE"
fi
if [[ -z "$OCCT_DIR" && -d /opt/homebrew/lib/cmake/opencascade ]]; then
    OCCT_DIR="/opt/homebrew/lib/cmake/opencascade"
fi

src_commit() {
    if [[ -e "$SRC_DIR/.git" ]]; then
        git -C "$SRC_DIR" rev-parse HEAD 2>/dev/null || true
    elif [[ -f "$SRC_DIR/.stl2step-baseline-commit" ]]; then
        tr -d '[:space:]' < "$SRC_DIR/.stl2step-baseline-commit"
    fi
}

ensure_src() {
    local have
    have="$(src_commit || true)"
    if [[ "$FORCE" -eq 1 && -e "$SRC_DIR" ]]; then
        # Drop a previous checkout. If it is a registered worktree, remove
        # it through git so the worktree list stays accurate.
        if [[ -e "$SRC_DIR/.git" ]]; then
            git -C "$REPO_ROOT" worktree remove --force "$SRC_DIR" 2>/dev/null \
                || rm -rf "$SRC_DIR"
        else
            rm -rf "$SRC_DIR"
        fi
        have=""
    fi
    if [[ -n "$have" && "$have" == "$COMMIT" && -f "$SRC_DIR/CMakeLists.txt" ]]; then
        echo "baseline source already at ${COMMIT_SHORT} ($SRC_DIR)"
        return 0
    fi
    if [[ -e "$SRC_DIR" && -n "$(ls -A "$SRC_DIR" 2>/dev/null || true)" ]]; then
        echo "error: $SRC_DIR exists but is not commit $COMMIT_SHORT" >&2
        echo "       re-run with --force to replace it" >&2
        exit 2
    fi

    echo "checking out ${COMMIT_SHORT} without disturbing the current worktree"
    # Prefer a detached worktree. Nested-worktree adds (this directory
    # lives inside another worktree of the same repo) can be refused; the
    # equivalent that still does not touch HEAD is `git archive`.
    if git -C "$REPO_ROOT" worktree add --detach "$SRC_DIR" "$COMMIT"; then
        echo "baseline source: git worktree (detached $COMMIT_SHORT)"
        return 0
    fi
    echo "worktree add refused (likely nested); falling back to git archive"
    mkdir -p "$SRC_DIR"
    git -C "$REPO_ROOT" archive "$COMMIT" | tar -x -C "$SRC_DIR"
    printf '%s\n' "$COMMIT" > "$SRC_DIR/.stl2step-baseline-commit"
    echo "baseline source: git archive $COMMIT_SHORT -> $SRC_DIR"
}

ensure_src

if [[ "$FORCE" -eq 1 && -d "$BUILD_DIR" ]]; then
    rm -rf "$BUILD_DIR"
fi
mkdir -p "$BUILD_DIR"

need_configure=0
if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
    need_configure=1
else
    have_gen="$(cache_get "$BUILD_DIR/CMakeCache.txt" CMAKE_GENERATOR || true)"
    have_bt="$(cache_get "$BUILD_DIR/CMakeCache.txt" CMAKE_BUILD_TYPE || true)"
    mismatch=0
    if [[ -n "$GENERATOR" && "$have_gen" != "$GENERATOR" ]]; then
        mismatch=1
    fi
    # Multi-config generators leave CMAKE_BUILD_TYPE empty in the cache.
    if [[ -n "$have_bt" && "$have_bt" != "$BUILD_TYPE" ]]; then
        mismatch=1
    fi
    if [[ "$mismatch" -eq 1 ]]; then
        echo "build cache generator/type mismatch ($have_gen/$have_bt vs ${GENERATOR:-<default>}/$BUILD_TYPE); reconfiguring"
        need_configure=1
    fi
fi

CMAKE_ARGS=(
    -S "$SRC_DIR"
    -B "$BUILD_DIR"
)
if [[ -n "$GENERATOR" ]]; then
    CMAKE_ARGS+=(-G "$GENERATOR")
fi
CMAKE_ARGS+=(
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    -DSTL2STEP_BUILD_CLI=ON
    -DSTL2STEP_BUILD_TESTS=OFF
    -DSTL2STEP_BUILD_EXAMPLES=OFF
    -DSTL2STEP_INSTALL=OFF
)
if ! is_vs_generator "$GENERATOR"; then
    if [[ -n "$C_COMPILER" ]]; then
        CMAKE_ARGS+=(-DCMAKE_C_COMPILER="$C_COMPILER")
    fi
    if [[ -n "$CXX_COMPILER" ]]; then
        CMAKE_ARGS+=(-DCMAKE_CXX_COMPILER="$CXX_COMPILER")
    fi
fi
if [[ -n "$TOOLCHAIN" ]]; then
    CMAKE_ARGS+=(-DCMAKE_TOOLCHAIN_FILE="$TOOLCHAIN")
fi
if [[ -n "$OCCT_DIR" ]]; then
    CMAKE_ARGS+=(-DOpenCASCADE_DIR="$OCCT_DIR")
fi

if [[ "$need_configure" -eq 1 ]]; then
    echo "configuring baseline (${GENERATOR:-<cmake default>} / $BUILD_TYPE)"
    cmake "${CMAKE_ARGS[@]}"
else
    echo "reusing existing baseline CMake cache (${GENERATOR:-<cmake default>} / $BUILD_TYPE)"
fi

echo "building baseline stl2step (-j $JOBS)"
cmake --build "$BUILD_DIR" --target stl2step --config "$BUILD_TYPE" -j "$JOBS"

# After a no -G configure, report the generator CMake actually picked.
if [[ -z "$GENERATOR" ]]; then
    GENERATOR="$(cache_get "$BUILD_DIR/CMakeCache.txt" CMAKE_GENERATOR || true)"
fi

BIN=""
for candidate in \
    "$BUILD_DIR/stl2step" \
    "$BUILD_DIR/stl2step.exe" \
    "$BUILD_DIR/$BUILD_TYPE/stl2step" \
    "$BUILD_DIR/$BUILD_TYPE/stl2step.exe"
do
    if [[ -f "$candidate" ]]; then
        BIN="$candidate"
        break
    fi
done
if [[ -z "$BIN" ]]; then
    echo "error: expected executable under $BUILD_DIR (stl2step[.exe] or $BUILD_TYPE/stl2step[.exe])" >&2
    exit 1
fi

BIN_EMIT="$(emit_path "$BIN")"
VERSION_OUT="$("$BIN" --version)"
{
    echo "COMMIT=$COMMIT"
    echo "GENERATOR=${GENERATOR}"
    echo "BUILD_TYPE=$BUILD_TYPE"
    echo "SRC=$SRC_DIR"
    echo "BIN=$BIN_EMIT"
    echo "VERSION=$VERSION_OUT"
} > "$BUILD_DIR/.baseline-meta"

echo
echo "BASELINE_BIN=$BIN_EMIT"
echo "BASELINE_COMMIT=$COMMIT"
echo "BASELINE_GENERATOR=${GENERATOR}"
echo "BASELINE_BUILD_TYPE=$BUILD_TYPE"
echo "BASELINE_VERSION=$VERSION_OUT"
echo "$BIN_EMIT" > "$BUILD_DIR/BIN_PATH"
