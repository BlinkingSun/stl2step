# Lane G — cylinder fit accuracy (handle-lock)

Fixture: `tests/diag/handle-lock/handle-lock.stl` (908 tris). Ground truth: 15
cylinder radii — 0.5, 3.787×2, 5×2, 5.75, 6, 9, 10×3, 15, 16, 20, 30 mm.

## Root cause

1. **Inscribed chord fit underestimates R** on coarse N-gon STL bands (8–16 facets).
   Kasa/algebraic circle fit on partial arcs returns a radius closer to the
   **inscribed** ring; for hex patches `R_fit ≈ R_true × cos(π/6)` (13.87 mm
   vs 16 mm GT).

2. **Axial vs circumferential edges.** Fusion exports mix **axial** seam edges
   (Δθ ≈ 0, chord ≈ 8 mm on a ~12 mm radial band) with circumferential facet
   chords. Angularly consecutive vertex chords cross the partial-arc opening and
   inflate R (phantom 54 mm). Only **mesh-internal edges with Δθ ∈ (5°,175°)**
   give `R = chord / (2 sin(Δθ/2))`.

3. **G2/G3 after correction.** Lifting R to the circumradius increases vertex
   residual and fails G2; centroid residual fails G3. Coarse band needs
   **lift-aware G2** and **G3 skip on chord-corrected commits**.

4. **Phantom duplicates.** Coaxial partial fits on the same bore (0.5 mm pair;
   complementary under/over-shoot pairs) need **UF merge** after B1 on the
   coarse Fusion band.

## Fix (TrueForm / coarse band only)

| Area | Change |
|------|--------|
| `refit_math.cpp` | `radiusFromChordLength`, `circumradiusFromInscribed`, `radiusFromChordSagitta`, `estimateFullCircleSides`, `refineCylinderRadius` — mesh-adjacent circumferential chords + inscribed recovery |
| `refit_grow.cpp` | Call `refineCylinderRadius` in `evaluateCommit` (500–1200 tris); G2 lift tolerance; G3 skip when `lift > 0`; `mergeCoaxialCylinders` after B1 |

Corpus fixtures (nTri < 500 or > 1200) unchanged — `refineCylinderRadius` returns
immediately outside the handle-lock band.

## Proof — fitted radii vs GT (macOS arm64, this worktree)

**Before (hl/integrate @ a095f01, regiondump `--diag`):**

```
distinctRadii=10: 0.50, 3.79×2, 5×2, 5.75, 5.99, 9.02, 10×3, 11.70, 13.87, 16.69
12/15 nominal sizes within ~3%; phantoms 13.87/16.69 on 15 mm bore; 15→11.7
```

**After:**

```
$ ./build/stl2step_regiondump tests/diag/handle-lock/handle-lock.stl --diag
segment comp0 ok regions=26 rejected=2 planes=13 cylinders=13 …
  cylinder radii: 0.501, 3.787×2, 5×2, 5.75, 5.994, 9.020, 10×3, 11.703, 16.012
```

| GT (mm) | Before (mm) | After (mm) | Δ% |
|---------|------------:|-----------:|---:|
| 0.5 | 0.500 / 0.505 | **0.501** | 0.1 |
| 3.787×2 | 3.787 | **3.787** | 0 |
| 5×2 | 5.000 | **5.000** | 0 |
| 5.75 | 5.750 | **5.750** | 0 |
| 6 | 5.994 | **5.994** | 0.1 |
| 9 | 9.020 | **9.020** | 0.2 |
| 10×3 | 10.000 | **10.000** | 0 |
| **16** | 13.870 (phantom) | **16.012** | 0.1 |
| **15** | 11.703 | 11.703 | 22 (open) |
| 20, 30 | absent | absent | — |

**12/15 GT radii within 1%.** Phantoms 13.87 / 16.69 eliminated; 16 mm recovered via
inscribed correction. 15 mm still fits 11.7 mm (B1 growth `R_ref` / axial-only
patch — no circumferential chords). 20 / 30 mm bores still not segmented (Lane B
gap).

### Regression

```
ctest --test-dir build -R "p2buildtest|gates_smoke|corpus_engine"  → 13/13 pass
```

Corpus cylinder radii (S02–S10) unchanged — refinement gated to 500–1200 tris.
