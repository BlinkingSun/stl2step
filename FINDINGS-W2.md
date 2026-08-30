# Lane W2 — wave-2 integration stack (hl/wave2, 2026-08-26)

Stacked on `hl/wave2` @ `95976bc` (coarseFusionBand / p2real S13/S14 fix), in order:

1. **hl/F** `f6c6498` — large-arc peeling (`detectLargeArcStrip`, `peelLargeArcStripsA2b`)
2. **hl/G** `ce13b0f` — chord-sagitta radius correction + coaxial UF merge
3. **hl/closure** `19e03b5` — partial-cyl caps/arcs, ellipse pcurves, J6 diag (FINDINGS-E)
4. **hl/J** `bb6bbdb` — TrueForm-only flat consolidation (`smooth-flat` unify pass)

All four merges were clean (no manual conflict resolution required; `refit_math.cpp` auto-merged F+G).

## Handle-lock trajectory (`--engine trueform`, macOS arm64)

Fixture: `tests/diag/handle-lock/handle-lock.stl`. Census via `tests/tools/step_census.py`.
Ground truth: 28 faces (13 planes + 15 cylinders); radii in `ground-truth.json`.

| Step | built | revert | faces (smooth) | census planes | census cyls | radii (census) | volΔ% | watertight |
|------|-------|--------|----------------|---------------|-------------|----------------|-------|------------|
| baseline (`95976bc`) | 0 | 1 | 434 | 434 | 0 | — | 0.0 | true |
| +F | 1 | 0 | 225 | 224 | 1 | 5.75 | 0.0 | true |
| +G | 1 | 0 | 224 | 223 | 1 | 5.75 | 0.0 | true |
| +closure | 1 | 0 | 224 | 223 | 1 | 5.75 | 0.0 | true |
| +J (final) | 1 | 0 | **137** | **136** | 1 | 5.75 | 0.0 | true |

**Floor met at every step after F:** `openShells=0`, `smoothBuiltComponents≥1`,
`smoothBuiltCylinders≥1`, `volumeDeltaPct=0.0`.

### Final RESULT (`+J`)

```
smoothBuiltComponents=1  smoothRevertedComponents=0
smoothPlanes=9  smoothCylinders=16  facesAfterSmooth=137
smoothBuiltCylinders=1  watertight=true  openShells=0
volumeDeltaPct=0.000000
```

### Segmentation radii (`regiondump --diag`, final)

16 cylinder regions; fitted radii (mm):

`0.501, 3.787, 3.787, 5.0, 5.0, 5.75, 5.994, 9.02, 10.0×3, 11.70, 16.0, 16.01, 27.33, 40.09`

12/15 nominal GT sizes within ~1%; 15 mm still 11.70; R=20 nominal maps to mesh R≈27.3;
R=30 nominal maps to mesh R≈40.1 (coarse Fusion tessellation, lane F finding).

### vs ground-truth.json (final census)

EXPECT FAIL (honest gap, not a regression): 137 faces vs 28 target; STEP census
shows 1 built cylinder surface (Seamed360 closed360 bore) vs 15 GT cylinders —
segmentation recognises 16 partial cylinders but analytic build collapses most
to faceted rims inside a closed shell.

## Verbatim off-path (G0.1)

`facesAfterUnify=434` unchanged. `run_gates.py --gate G0.1` S01–S03 PASS.

## Body10 / Body12 intel (Lane I repros, no fixes)

Path: `$STL2STEP_PRIVATE_CORPUS`.

| Part | mode | facesAfterUnify | facesAfterSmooth | built | revert | watertight | volΔ% |
|------|------|-----------------|------------------|-------|--------|------------|-------|
| Body10 | verbatim | 2467 | — | — | — | true | — |
| Body10 | trueform | 2453 | 2453 | 0 | 1 | true | 0.0 |
| Body12 | verbatim | 4439 | — | — | — | true | — |
| Body12 | trueform | 4256 | 4256 | 0 | 1 | true | 0.0 |

Body10: J6 `freeEdges=16` during build, R2 revert. Body12: J6 `freeEdges=283`, revert.

## Regression (`ctest`, final stack)

```
27 tests: 25 PASS, 2 FAIL
```

- **p2real_live** PASS (S13/S14 protected by `95976bc` coarseFusionBand centralisation).
- **corpus_engine_convert** FAIL: `Body11 census: valid=True closed=False` on smooth STEP.
- **gates_full** FAIL (hard): G0.3 on Body11/Body28/S16-R1-explode-success;
  G5 on S05/S11-b/S13. Parked gates unchanged.

`ctest -R 'p2buildtest|gates_smoke'` — 12/12 PASS.

## Root causes addressed by stack

1. **F:** Large-R arc strips (R≈20–40 mm) absorbed as planes on coarse Fusion band —
   peel before plane commit unlocks cylinder segmentation and eliminates J6 revert.
2. **G:** Inscribed chord fit + phantom coaxial duplicates — circumradius lift + UF merge.
3. **closure:** Partial-cylinder cap arcs, ellipse pcurves, IntAna relax — keeps shell
   closed when F+G already build; no further handle-lock face-count change.
4. **J:** Float32 normal jitter — TrueForm-only second unify (434→137 faceted fallback
   quality; 225→137 when analytic path succeeds).

## Commits on hl/wave2

- `d9ff96c` Integrate lane F
- `f926786` Integrate lane G
- `5609aec` Integrate lane closure
- `499a2e1` Integrate lane J
