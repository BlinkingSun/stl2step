# Lane CYLEDGES — cyl|cyl edge construction (2026-08-27)

Branch: `hl/cyledges`. Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris).

## Root cause

Recognition finds **16 cylinders**; build kept **1** analytic cylinder in STEP because
**7/8** analytic cyl|cyl boundary chains hit `IntAna cyl|cyl empty/same` on coarse
Fusion exports. Fitted axes are a few degrees off true parallel/tangent, so
`IntAna_QuadQuadGeo` returns empty/degenerate and `g1Tangent` (3° + tight `epsPlane`)
never routed chains to `constructedCylCylGenerator`. Failed chains enrolled
`chainEdgeFail` → J6/R1 exploded partial cylinders → census stayed at 1×
`CYLINDRICAL_SURFACE` (R≈5.75 closed360 only).

## Fixes (TrueForm path only — `refit_build.cpp`, `refit_chains.cpp`)

1. **`g1Tangent` cyl|cyl** — widen parallel gate to 8°; use
   `max(epsPlane, maxVertexDev)` for tangent distance (matches plane|cyl).
2. **`intAnaAcceptResidual`** — higher band for cyl|cyl and short plane|cyl chains.
3. **`intersectSurfaces` cyl|cyl** — after IntAna miss: `cylCylTangentContact` →
   `constructedCylCylGenerator`; `cylCylParallelOffset` → `bestCylCylConstructed`
   (both surface branches + mesh-anchored generator, residual-gated).
4. **Line snap cap** — cyl|cyl generator edges use the same relaxed terminal snap
   as plane|plane.
5. **Oblique plane|cyl ellipse** — second `pickIntAna` pass with loose residual
   when `TypeInter==Ellipse` and `nV≤3` (unblocks reg 4|14 chain that was
   cascading R1).

## Proof (macOS arm64, this worktree)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest --test-dir build   # 27/27 pass
```

### handle-lock TrueForm

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-cy.step --engine trueform --quiet
```

**Before (branch tip @ start):** 7× `IntAna cyl|cyl empty/same`, 1× plane|cyl,
`DIAG_COLLAPSE ok=60 none=8`, `smoothBuiltCylinders=1`.

**After:**

```
RESULT … "warnings":[],"smoothBuiltCylinders":1,"smoothBuiltComponents":1,
"smoothRevertedComponents":0,"volumeDeltaPct":0.000000,"watertight":true
DIAG_COLLAPSE mix=0 none=0 fail=0 ok=68 total=68 recover=0 rounds=0
```

- cyl|cyl IntAna warnings: **7 → 0**
- Analytic chain collapse: **60/68 → 68/68**
- Total TrueForm warnings: **8 → 0**
- `volumeDeltaPct` 0.0, watertight, no R2 revert

```sh
python3 tests/tools/step_census.py /tmp/hl-cy.step
# surfaces.cylinder=1 (rid=13 partial build still ensure-invalid → R1 drops partials)
```

`smoothBuiltCylinders` remains 1: all cyl|cyl **edges** now build, but
`buildPartialCylinder` rid=13 (`R≈6 mm`) still `ensure-invalid` — partial-face
wiring is Lane E / closure follow-on, not edge construction.

### Body12.stl (second fixture)

```sh
./build/stl2step "$STL2STEP_PRIVATE_CORPUS/Body12.stl" \
  -o /tmp/b12.step --engine trueform --quiet
# cyl|cyl warnings: 108 (skew/non-parallel mother lode — needs further conic accept)
```

## Remaining gaps

1. **Partial-cylinder face validity** (rid 13 on handle-lock) blocks census growth
   after edge collapse is complete.
2. **Body12** skew cyl|cyl pairs still miss constructed fallbacks (108 hits).
3. **Ground-truth census** (15 cylinders) not reached on handle-lock yet.
