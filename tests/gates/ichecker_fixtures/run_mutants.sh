#!/bin/sh
# Mutant runner for check_regionset.py.
#
# tests/gates/CMakeLists.txt is owned by p0-gates and is not in this tree,
# so this lane does not register ctest. Run:
#   sh tests/gates/ichecker_fixtures/run_mutants.sh
set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
CHECKER="$HERE/../check_regionset.py"
fail=0
n=0

expect_fail() {
    dump="$1"
    rule="$2"
    shift 2
    n=$((n + 1))
    json=$(python3 "$CHECKER" "$dump" --json "$@" 2>/dev/null) || true
    got=$(printf '%s\n' "$json" | python3 -c "
import json,sys
r=json.load(sys.stdin)
rule=sys.argv[1]
rec=r.get('rules',{}).get(rule,{})
print('ok' if (not r.get('ok') and rec.get('pass') is False) else 'bad')
print(rule, 'pass=' + str(rec.get('pass')), 'offenders=' + ','.join(rec.get('offenders') or []))
" "$rule")
    status=$(printf '%s\n' "$got" | head -n1)
    detail=$(printf '%s\n' "$got" | tail -n1)
    if [ "$status" = "ok" ]; then
        printf 'PASS  %s  failed %s (%s)\n' "$(basename "$dump")" "$rule" "$detail"
    else
        printf 'FAIL  %s  expected rule %s to fail; %s\n' "$(basename "$dump")" "$rule" "$detail"
        fail=$((fail + 1))
    fi
}

expect_pass() {
    dump="$1"
    shift 1
    n=$((n + 1))
    if python3 "$CHECKER" "$dump" "$@" >/dev/null; then
        printf 'PASS  %s  exit 0\n' "$(basename "$dump")"
    else
        printf 'FAIL  %s  expected exit 0\n' "$(basename "$dump")"
        fail=$((fail + 1))
    fi
}

expect_fail "$HERE/I1_both.json" I1
expect_fail "$HERE/I2_two_chains.json" I2
expect_fail "$HERE/I3_len.json" I3
expect_fail "$HERE/I4_rms.json" I4
expect_fail "$HERE/I5_order.json" I5
expect_fail "$HERE/I5_tie.json" I5
expect_fail "$HERE/I6_dirty.json" I6
expect_fail "$HERE/I7_two_outer.json" I7
expect_fail "$HERE/I7b_outer.json" I7b
expect_fail "$HERE/I8_unsplit.json" I8
expect_fail "$HERE/I9_sign.json" I9
expect_fail "$HERE/G5_filed.json" G5
expect_fail "$HERE/G5_filed.json" G5 --sidecar "$HERE/G5_filed.expected.json"
expect_fail "$HERE/sidecar_s11_filletstrip.json" SIDECAR --sidecar "$HERE/sidecar_s11.expected.json"

expect_pass "$HERE/../regionset.example.min.json"
expect_pass "$HERE/../regionset.example.closed360.json"
expect_pass "$HERE/valid_closed360_inner.json"
expect_pass "$HERE/sidecar_cone.json" --sidecar "$HERE/sidecar_cone.expected.json"
expect_pass "$HERE/sidecar_sphere.json" --sidecar "$HERE/sidecar_sphere.expected.json"
expect_pass "$HERE/sidecar_s12b.json" --sidecar "$HERE/sidecar_s12b.expected.json"
expect_pass "$HERE/sidecar_s11.json" --sidecar "$HERE/sidecar_s11.expected.json"

# --rule filtering
n=$((n + 1))
if python3 "$CHECKER" "$HERE/I7b_outer.json" --rule I7b >/dev/null; then
    printf 'FAIL  --rule I7b on I7b mutant: expected nonzero\n'
    fail=$((fail + 1))
else
    printf 'PASS  --rule I7b on I7b mutant (nonzero)\n'
fi
n=$((n + 1))
if python3 "$CHECKER" "$HERE/I1_both.json" --rule I7b >/dev/null; then
    printf 'PASS  --rule I7b on I1 mutant (zero)\n'
else
    printf 'FAIL  --rule I7b on I1 mutant: expected zero\n'
    fail=$((fail + 1))
fi

printf '\n%d checks, %d failed\n' "$n" "$fail"
[ "$fail" -eq 0 ]
