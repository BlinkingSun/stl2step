#!/usr/bin/env bash
# install-git-hooks.sh — copy committed hooks into this clone's .git/hooks.
# Git does not honour scripts/hooks/ by itself; each clone must run this once.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "install-git-hooks: not a git work tree: $REPO_ROOT" >&2
  exit 1
fi

HOOKS="$(git rev-parse --git-path hooks)"
mkdir -p "$HOOKS"

src="$REPO_ROOT/scripts/hooks/pre-push"
if [ ! -f "$src" ]; then
  echo "install-git-hooks: missing $src" >&2
  exit 1
fi
cp "$src" "$HOOKS/pre-push"
chmod +x "$HOOKS/pre-push"

echo "install-git-hooks: installed $HOOKS/pre-push"
echo "  -> runs scripts/ci-local-gate.sh on every push"
