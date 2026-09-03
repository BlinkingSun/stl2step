#!/usr/bin/env bash
# vendored from agent-team-v5 templates/node-lock.sh, protocol 1
# keep byte-identical below this line
# Vendorable single-tenant node lock — same protocol as bin/team-node-lock.
# Public repos may copy to scripts/node-lock.sh (NODE_LOCK_PROTOCOL=1).
NODE_LOCK_PROTOCOL=1
set -euo pipefail

usage() {
    cat <<'EOF'
node-lock.sh — single-tenant node lock (NODE_LOCK_PROTOCOL=1)

Usage:
  node-lock.sh acquire --name <run-name> [--lock-dir <path>] [--wait <sec>] [--ttl <sec>]
  node-lock.sh release --name <run-name> [--lock-dir <path>]
  node-lock.sh status  --name <run-name> [--lock-dir <path>]
  node-lock.sh heartbeat --name <run-name> [--lock-dir <path>]

Exit codes: 0 success/free; 3 busy; 4 release refused (wrong owner name)
EOF
}

NAME=""
LOCK_DIR="${HOME:-~}/.team-node-lock"
WAIT_SEC=0
TTL=7200
VERB=""

while [ $# -gt 0 ]; do
    case "$1" in
        acquire|release|status|heartbeat)
            VERB="$1"; shift ;;
        --name) NAME="${2:?}"; shift 2 ;;
        --lock-dir) LOCK_DIR="${2:?}"; shift 2 ;;
        --wait) WAIT_SEC="${2:?}"; shift 2 ;;
        --ttl) TTL="${2:?}"; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) echo "node-lock.sh: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ -z "$VERB" ]; then
    echo "node-lock.sh: verb required" >&2
    exit 2
fi
if [ -z "$NAME" ]; then
    echo "node-lock.sh: --name is required" >&2
    exit 2
fi

_parse_owner_field() {
    local line="$1" key="$2"
    printf '%s' "$line" | tr ' ' '\n' | sed -n "s/^${key}=//p" | head -1
}

_local_hostname() {
    hostname 2>/dev/null || hostname -s 2>/dev/null || echo unknown
}

_pid_alive() {
    local pid="$1"
    [ -n "$pid" ] && [ "$pid" != "0" ] && kill -0 "$pid" 2>/dev/null
}

_is_stale() {
    local owner_line="$1" now started ttl age owner_host owner_pid
    now=$(date +%s)
    started="$(_parse_owner_field "$owner_line" started)"
    ttl="$(_parse_owner_field "$owner_line" ttl)"
    owner_host="$(_parse_owner_field "$owner_line" host)"
    owner_pid="$(_parse_owner_field "$owner_line" pid)"
    [ -n "$started" ] || return 0
    [ -n "$ttl" ] || ttl=7200
    age=$((now - started))
    if [ "$age" -ge "$ttl" ]; then
        STALE_AGE="$age"
        STALE_NAME="$(_parse_owner_field "$owner_line" name)"
        return 0
    fi
    if [ "$owner_host" = "$(_local_hostname)" ] && ! _pid_alive "$owner_pid"; then
        STALE_AGE="$age"
        STALE_NAME="$(_parse_owner_field "$owner_line" name)"
        return 0
    fi
    return 1
}

_break_stale() {
    local owner_line="$1"
    local stale_name="${STALE_NAME:-$(_parse_owner_field "$owner_line" name)}"
    local stale_age="${STALE_AGE:-0}"
    echo "stale-lock-broken name=${stale_name} age=${stale_age}" >&2
    rm -rf "$LOCK_DIR"
}

_read_owner() {
    if [ -f "${LOCK_DIR}/owner" ]; then
        tr -d '\r' < "${LOCK_DIR}/owner" | head -1
    fi
}

_write_owner() {
    printf '%s\n' "$1" > "${LOCK_DIR}/owner"
}

ACQUIRE_LINE=""
_acquire_once() {
    local host pid started owner_line
    if mkdir "$LOCK_DIR" 2>/dev/null; then
        host="$(_local_hostname)"
        pid=$PPID
        started=$(date +%s)
        owner_line="name=${NAME} host=${host} pid=${pid} started=${started} ttl=${TTL}"
        _write_owner "$owner_line"
        ACQUIRE_LINE="$owner_line"
        return 0
    fi
    owner_line="$(_read_owner || true)"
    if [ -z "$owner_line" ]; then
        ACQUIRE_LINE=""
        return 1
    fi
    if _is_stale "$owner_line"; then
        _break_stale "$owner_line"
        ACQUIRE_LINE=""
        return 2
    fi
    ACQUIRE_LINE="$owner_line"
    return 1
}

_acquire() {
    local retry=0 rc
    set +e
    while :; do
        _acquire_once
        rc=$?
        if [ "$rc" -eq 0 ]; then
            set -e
            echo "acquired local ${NAME}"
            return 0
        fi
        if [ "$rc" -eq 2 ] && [ "$retry" -eq 0 ]; then
            retry=1
            continue
        fi
        set -e
        printf '%s' "$ACQUIRE_LINE"
        return 1
    done
}

_release() {
    local owner_line owner_name
    if [ ! -d "$LOCK_DIR" ]; then
        return 0
    fi
    owner_line="$(_read_owner || true)"
    owner_name="$(_parse_owner_field "${owner_line:-}" name)"
    if [ "$owner_name" != "$NAME" ]; then
        return 4
    fi
    rm -rf "$LOCK_DIR"
    return 0
}

_status() {
    local owner_line
    if [ ! -d "$LOCK_DIR" ]; then
        echo free
        return 0
    fi
    owner_line="$(_read_owner || true)"
    if [ -z "$owner_line" ]; then
        echo free
        return 0
    fi
    echo "$owner_line"
    return 3
}

_heartbeat() {
    local owner_line owner_name host pid started ttl
    owner_line="$(_read_owner || true)"
    owner_name="$(_parse_owner_field "${owner_line:-}" name)"
    if [ "$owner_name" != "$NAME" ]; then
        return 4
    fi
    host="$(_parse_owner_field "$owner_line" host)"
    pid="$(_parse_owner_field "$owner_line" pid)"
    started=$(date +%s)
    ttl="$(_parse_owner_field "$owner_line" ttl)"
    [ -n "$ttl" ] || ttl="$TTL"
    _write_owner "name=${NAME} host=${host} pid=${pid} started=${started} ttl=${ttl}"
    return 0
}

_poll_acquire() {
    local out rc start now elapsed
    start=$(date +%s)
    while :; do
        set +e
        out="$(_acquire 2>&1)"
        rc=$?
        set -e
        if [ "$rc" -eq 0 ]; then
            echo "$out"
            return 0
        fi
        now=$(date +%s)
        elapsed=$((now - start))
        if [ "$WAIT_SEC" -gt 0 ] && [ "$elapsed" -lt "$WAIT_SEC" ]; then
            sleep_for=$((WAIT_SEC - elapsed))
            [ "$sleep_for" -gt 15 ] && sleep_for=15
            sleep "$sleep_for"
            continue
        fi
        echo "busy: ${ACQUIRE_LINE:-$out}"
        return 3
    done
}

case "$VERB" in
    acquire)
        if [ "$WAIT_SEC" -gt 0 ]; then
            _poll_acquire
            exit $?
        fi
        set +e
        out="$(_acquire 2>&1)"
        rc=$?
        set -e
        if [ "$rc" -eq 0 ]; then
            echo "$out"
            exit 0
        fi
        echo "busy: ${ACQUIRE_LINE:-$out}"
        exit 3
        ;;
    release)
        set +e
        _release >/dev/null 2>&1
        rc=$?
        set -e
        exit "$rc"
        ;;
    status)
        set +e
        out="$(_status 2>/dev/null)"
        rc=$?
        set -e
        echo "$out"
        exit "$rc"
        ;;
    heartbeat)
        set +e
        _heartbeat >/dev/null 2>&1
        rc=$?
        set -e
        exit "$rc"
        ;;
esac
