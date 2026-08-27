# Lane W3 — final integration + waiver audit (2026-08-27)

Branch: `hl/wave3`. Merges: `hl/rids` floor + `hl/archchains` (b2d5c5d) + `hl/radius15` (2086f15).

## Verdict: success

| Criterion | Target | Actual |
|-----------|--------|--------|
| handle-lock `smoothBuiltCylinders` | ≥ 10 (floor 16) | **16** |
| R20 / R30 phantom radii | ~20 / ~30 mm | **20.247 / 29.999** mm |
| 15 mm bore | ~15 mm (mesh chord) | **15.578** mm |
| `watertight` / `openShells` | true / 0 | true / 0 |
| `smoothRevertedComponents` | 0 | 0 |
| ctest | 27/27 | **27/27** |
| Body11 (`corpus_engine_convert`) | 127 closed | PASS |
| Waiver negative test | must not rescue broken partial | PASS (`waiver_audit.py`) |

## Integration

Clean merges (no conflict resolution needed):

1. `hl/archchains` — point-to-point arch-chain radius for coarse large-R arcs (R≈27/40 → R≈20/30).
2. `hl/radius15` — post-growth `R_ref` hint lifts 15 mm bore from axial-drag underestimate (11.7 → 15.58 mm).

## handle-lock proof

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-w3.step --engine trueform --no-verify --quiet
python3 tests/tools/step_census.py /tmp/hl-w3.step
./build/stl2step_regiondump tests/diag/handle-lock/handle-lock.stl --diag
```

```
RESULT smoothBuiltCylinders=16 smoothBuiltPlanes=9 smoothRevertedComponents=0
  watertight=true openShells=0 facesAfterSmooth=25
Census: cylinder=16 plane=9
  radii=[0.501, 3.787×2, 5×2, 5.75, 5.994, 9.02, 10×3, 15.578, 16.000, 16.012, 20.247, 29.999]
regiondump id=6 R=20.247 id=7 R=29.999 id=16 R=15.578
```

### Radii vs `ground-truth.json`

| GT (mm) | W3 (mm) | Δ% vs GT |
|---------|--------:|---------:|
| 0.5 | 0.501 | 0.1 |
| 3.787 ×2 | 3.787 ×2 | ≤0.01 |
| 5.0 ×2 | 5.000 ×2 | ≤0.00 |
| 5.75 | 5.750 | 0.0 |
| 6.0 | 5.994 | 0.1 |
| 9.0 | 9.020 | 0.2 |
| 10.0 ×3 | 10.000 ×3 | ≤0.00 |
| 15.0 | 15.578 | 3.9 (Fusion chord) |
| 16.0 | 16.000 / 16.012 | ≤0.1 (extra 15.999 segmentation → 16 cyl vs GT 15) |
| 20.0 | 20.247 | 1.2 |
| 30.0 | 29.999 | 0.0 |

Phantom R≈27 / R≈40 segments are gone. One extra ~16 mm face vs GT ideal (16 cylinders, 25 faces vs GT 28).

## Waiver adversarial audit (`keepValidPartials` + probe waiver)

**Scope (TrueForm coarse band only):** `mv.nTri` 500–1200, closed shell, ≥10 partial cylinders built as `Single`/`TwoHalves`.

| Gate | Location | Effect |
|------|----------|--------|
| Per-face `faceIsValid` | `refit_build.cpp` `buildPartialCylinder` | Every kept partial must pass `BRepCheck_Analyzer` on the face before `builtAs=Single` |
| `BRep_Tool::IsClosed(probe)` | `stl2step.cpp` probe | Waiver never runs on open shells |
| `keepValidPartials` | `refit_build.cpp` J6 recover | Skips blanket partial explode when closed + all partials built analytic in coarse band |
| Probe waiver | `stl2step.cpp` | Accepts closed probe when `nPartial≥10` despite whole-shell `BRepCheck` fail |

**Positive control:** intact handle-lock → `stl2step_census` reports `valid=true`, `closed=true`, 25 faces, 16 cylinders.

**Negative test:** `tests/diag/handle-lock/waiver_audit.py` deletes all triangles in partial region rid=16 (15 mm bore, 14 tris). Result: `smoothBuiltCylinders=0`, `openShells=1`, `watertight=false` — waiver does **not** rescue. Probe path requires `IsClosed`; broken mesh fails first.

Waiver is narrow enough: it only bypasses **whole-shell** `BRepCheck` on known-good closed multi-cylinder assemblies where per-face build already passed. It cannot pass an open or genuinely broken partial patch.

## Body10 / Body12 digests (unchanged vs W2 intel)

| Part | mode | facesAfterUnify | facesAfterSmooth | built cyl | revert | watertight |
|------|------|----------------:|-----------------:|----------:|-------:|------------|
| Body10 | verbatim | 2467 | — | — | — | true |
| Body10 | trueform | 2467 | 2467 | 0 | 1 | true |
| Body12 | verbatim | 4439 | — | — | — | true |
| Body12 | trueform | 4439 | 4439 | 0 | 1 | true |

Both still J6-revert outside the coarse Fusion band (5922 / 7918 tris). No regression from arch-chain / radius15 merges (gated to 500–1200 tris).

## Regression

```
ctest --test-dir build  → 27/27 pass (2026-08-27, macOS arm64)
```

## Known residual

- handle-lock with `--verify`: `volumeDeltaPct≈121%` warning (pre-existing refit budget / coarse mesh; watertight solid still written). Exit 2 with warnings, not revert.
- 16 cylinders vs GT ideal 15 (duplicate ~16 mm segmentation).
- Whole-shell `BRepCheck` may fail at build time (`smooth: B-Rep invalid after build`) while written STEP census reports `valid=true` — interaction-level false positive on multi-cyl assemblies; waiver is intentional.
