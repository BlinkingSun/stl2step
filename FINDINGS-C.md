# Lane C — handle-lock plane|cyl IntAna / J6 (2026-08-26)

## Offending surface pair (numerical)

| Field | reg 55 (plane) | reg 17 (cylinder) |
|---|---|---|
| **chain** | ci=162, 1 edge, 2 verts | same |
| **R** | — | 9.0195 mm |
| **plane normal · cyl axis** | \|dot\| ≈ 5.6×10⁻¹² (≈90°, side-grazing) | |
| **dist(cyl axis, plane)** | 9.0417 mm | |
| **\|dist − R\|** | **0.0222 mm** | |
| **epsPlane floor** | 0.02 mm | |
| **cyl maxVertexDev** | — | 0.0238 mm |
| **IntAna TypeInter** | 7 (`IntAna_Ellipse`, degenerate near generator) | |
| **pickIntAna** | rejects (wrong branch / residual gate) | |

reg 55 is a 2-triangle sliver plane grazing the R≈9 mm partial cylinder.

## Mechanism

1. **P1 (`g1Tangent`)** marks plane|cyl chains tangent only when
   `|dist(axis,plane) − R| ≤ epsPlane`. Here `0.022 > 0.02` → `tangent=false`
   even though `0.022 ≤ cyl.maxVertexDev (0.024)`.

2. **P2 (`intersectSurfaces`)** calls `IntAna_QuadQuadGeo` plane|cyl. For this
   near-tangent side contact OCCT returns `IntAna_Ellipse`; `pickIntAna` returns
   `None` (degenerate ellipse / 1 mm residual gate on a 2-vertex chain).

3. **Dead fallback**: `constructedGenerator` was gated on
   `ch.tangent && planePerpCylinder` — mutually exclusive geometry (side-graze
   vs cap). Generator never ran.

4. **chainEdgeFail** enrolled reg 55|17 on pass 0 (before fix) → J6 R1 explode.
   Full component still ends at **J6 freeEdges=49 recover=2 → R2 revert**
   (parameter-insensitive; separate mixed-shell closure issue).

## Fix (TrueForm path only)

1. **`refit_chains.cpp`**: `g1Tangent` plane|cyl tolerance
   `max(epsPlane, pl.maxVertexDev, cy.maxVertexDev)`.

2. **`refit_build.cpp`**: `planeCylSideContact` + `constructedGenerator` fallback
   when IntAna/pickIntAna miss; `constructedPlaneCylCap` for cap plane|cyl when
   `planePerpCylinder`. Optional `STL2STEP_P2_DIAG=1` stderr on remaining misses.

## Before / after (macOS arm64, this worktree)

**Before** (`92c2905`):

```
warning: smooth: IntAna plane|cyl empty/same — keeping mesh polyline
warning: J6: shell not closed freeEdges=49 faces=519 recover=2
warning: smooth: analytic rebuild reverted on one component -- kept faceted
chain162 tangent: false
smoothBuiltComponents: 0  smoothRevertedComponents: 1
```

**After**:

```
warning: J6: shell not closed freeEdges=49 faces=519 recover=2
warning: smooth: analytic rebuild reverted on one component -- kept faceted
chain162 tangent: true
smoothBuiltComponents: 0  smoothRevertedComponents: 1
```

Plane|cyl IntAna warning **eliminated**; chain 162 gets analytic generator edge.
Whole-component R2 revert **unchanged** — 49 free edges persist after J6 recover
(mixed analytic + exploded facets; out of scope for this chain-specific fix).

## Gates

```
ctest -R "p2buildtest|gates_smoke"  → 12/12 passed
```
