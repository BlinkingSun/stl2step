# Lane RADIUS15 — post-growth cylinder radius (handle-lock 15 mm bore)

Fixture: `tests/diag/handle-lock/handle-lock.stl` (908 tris). GT cylinder radii include **15.0 mm** and **16.0 mm**.

## Root cause

B1 cylinder growth on the 15 mm bore (seed P=183 Q=220) locks `R_ref` from the first circumferential band (`R≈15.578` at `|S|=3`), then axial bands drag the bulk `eberlyCenterRadius` down (`R≈11.7` on the full patch). `refineCylinderRadius` in `evaluateCommit` could not lift back: chord acceptance was anchored to the dragged radius and the **1.25× lift cap** blocked the ~28% recovery needed. Circumferential facet chords on this patch are mostly **patch boundary** edges (thin axial bands), not mesh-internal edges.

## Fix (TrueForm coarse band only; `refit_math.cpp` + `refit_grow.cpp`)

| Area | Change |
|------|--------|
| `refit_grow.cpp` | Pass growth `R_ref` as `rHint` into `evaluateCommit`; track `R_ref = max(R_ref, radius)` during growth |
| `refit_math.cpp` | `refineCylinderRadius(..., rHint)`: widen chord window/cap toward `rHint`; partial-arc boundary chord pass; large-bore axial-drag lift when `rHint` exceeds bulk fit; pratt center/radius snap after lift |

Corpus / Body11 unchanged — `refineCylinderRadius` still gated to 500–1200 tris.

## Proof — handle-lock radii (`regiondump --diag`)

```
$ ./build/stl2step_regiondump tests/diag/handle-lock/handle-lock.stl --diag
segment comp0 ok regions=25 … cylinders=16
  id=16 type=cylinder … radius=15.5781  (15 mm GT bore)
  cylinder radii: …,15.578,…,16.012,20.247,29.999,…
```

| GT (mm) | Before (mm) | After (mm) | Δ% vs GT |
|---------|------------:|-----------:|---------:|
| **15** | 11.703 | **15.578** | 3.9 |
| 16 | 16.012 | **16.012** | 0.1 |
| 20 | 20.247 | **20.247** | 1.2 |
| 30 | 29.999 | **29.999** | 0.0 |
| 0.5–10 (×12) | unchanged | unchanged | ≤0.2 |

15 mm bore recovered from axial-drag underestimate (22% error → ~4% vs CAD). On this Fusion STL the circumferential band chord geometry reads **15.58 mm** (thin-band eberly matches commit); CAD GT is 15.0 mm — same ~4% overshoot class as coarse tessellation, not the prior 11.7 mm phantom.

## Regression

```
ctest --test-dir build  → 27/27 pass
Body11 --engine trueform → smoothBuiltCylinders=127, openShells=0, watertight=true
```
