# Lane J — TrueForm flat consolidation (2026-08-26)

## Root cause

Fusion float32 STL normals jitter beyond the frozen Verbatim coplanar gate
(`unifyAngleDeg = 0.001°`). `ShapeUpgrade_UnifySameDomain` at that angle leaves
handle-lock at 434 planar faces (908 triangles → 434 after unify). TrueForm
revert fallback inherited the same faceted shell — analytic rebuild failed
(J6, 49 free edges) but output stayed at 434 faces.

## Fix

TrueForm-only second `UnifySameDomain` pass after the standard unify step.
Angular tolerance derived from segmentation `maxVertexDev` (same θ = dev/L
model as `refit_grow.cpp`), with `epsMesh` floor when dev is zero, and a
**0.01°** minimum. Analytic faces stay protected via existing `KeepShape`.
Verbatim path unchanged (G0.1).

Implementation: `src/stl2step.cpp` — `smooth-flat` unify pass inside
`if (smooth)` only.

## Proof (macOS arm64, hl/J worktree)

### handle-lock TrueForm revert (unchanged analytic outcome)

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl -o /tmp/o.step \
  --engine trueform --quiet
```

```
facesAfterUnify=230 facesAfterSmooth=230 smoothRevertedComponents=1
volumeDeltaPct=0.000000 meshVolumeMM3=15868.884516 stepVolumeMM3=15868.884516
smoothMaxDevMM=0.000000
```

STEP census: `faces=230 plane=230 cylinder=0 line=683` (was 434/888).

### handle-lock Verbatim (byte-identical off-path)

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl -o /tmp/v.step \
  --quiet --no-verify
```

`facesAfterUnify=434` (unchanged).

### G0.1 S01–S03

```sh
python3 tests/gates/run_gates.py --fixture S01,S02,S03 --gate G0.1 \
  --binary ./build/stl2step --no-verify
```

`G0.1 PASS` on S01, S02, S03 (canonicalize.py step+result IDENTICAL vs 187ead0).

### Gates

`ctest -R 'gates_smoke|p2buildtest'` — 12/12 PASS.
