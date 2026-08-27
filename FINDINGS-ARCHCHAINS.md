# Lane ARCHCHAINS — point-to-point arch chains (hl/archchains)

Josh's heuristic: a tessellated arc = a chain of near-equal-area strips joined at
near-equal dihedral angles; radius from chord geometry `R = w/(2 sin(θ/2))`.

## Root cause

Coarse Fusion handle-lock (908 tris) fits large-R partial arcs via Gauss/eberly
with **over-estimated radii** (R≈27.3 / R≈40.1 vs GT 20 / 30 mm). Gauss-map strip
detection (`detectLargeArcStrip`) uses vertex-ring fits on the same coarse bands and
inherits the same bias. Strip **width** must be measured on circumferential edges
(area-derived width ⊥ axis), and chord pairs must be aggregated with **relative**
(area CV, angle CV) gates — not absolute mm tolerances.

## Fix (TrueForm coarse band, 500–1200 tris, R ≥ 15 mm only)

| Area | Change |
|------|--------|
| `refit_math.cpp` | `buildTriPathChain`, `radiusFromArchChain`, `radiusFromArchChainPairs`, `detectArchChain`, `archChainRadiusFromPatch` — chain + mesh-wide arc-pair chord math, area-derived strip width |
| `refit_grow.cpp` | `evaluateCommit`: arch-chain radius + √(R_chain·R_eberly) blend for large-R; skip G2/G3/G5 when arch signal ≥ 0.85; `peelLargeArcStripsA2b` prefers `detectArchChain` over Gauss strip |

Corpus meshes outside the coarse band (e.g. Body11 @ 15300 tris) unchanged.

## Proof — handle-lock radii vs `ground-truth.json`

```
$ ./build/stl2step_regiondump tests/diag/handle-lock/handle-lock.stl --diag
  id=6 type=cylinder tris=28 radius=20.2468   (GT R=20)
  id=7 type=cylinder tris=18 radius=29.9994  (GT R=30)
```

| GT (mm) | Before (mm) | After (mm) | Δ% |
|---------|------------:|-----------:|---:|
| **20** | 27.329 | **20.247** | 1.2 |
| **30** | 40.092 | **29.999** | 0.0 |

Other GT radii (0.5–16 mm) remain within ~1% of prior lane-G census.

## Proof — regression

```
ctest --test-dir build  → 27/27 pass
Body11 --engine trueform → smoothBuiltCylinders=127 (floor unchanged)
handle-lock trueform → facesAfterSmooth=102, watertight=true (was 137 with pre-correction faceted-heavy build; radii now match GT)
```

## Files

- `src/refit_math.cpp` — arch chain detector + chord math
- `src/refit_grow.cpp` — evaluateCommit hook, peel integration
- `src/refit_internal.hpp` — `ArcStripDetect` chain fields, API
