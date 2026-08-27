# Lane F findings — large-R arc recognition (hl/F, 2026-08-26)

## Root cause

Coarse Fusion exports (handle-lock @ 908 tris) absorb large-radius partial arcs into
plane provisionals because facet dihedrals fall below `thetaPlane` (15° on coarse band)
and below `thetaCylLo` (5°), so A2 plane growth and B1 cylinder seeding never fire.
The R≈20 mm strip (plane reg 5) has **57°** normal rotation but was committed as a
plane. The R≈30 mm GT feature is tessellated as a **very gentle** strip (reg 6,
**21°** normal span) whose vertices measure **R≈40 mm** on this STL — not R=30.

## Fix

1. **`detectLargeArcStrip`** (`refit_math.cpp`) — Gauss-map axis + sorted-angle
   monotonic normal rotation (static-normal fallback when vertex deviation exceeds
   plane tolerance). Large-R gate: 15 ≤ R ≤ 55 mm, coarse residual slack 5%·R.
2. **`peelLargeArcStripsA2b`** (`refit_grow.cpp`) — after A3 coplanar UF merge,
   before plane commit; coarse Fusion band only (500–1200 tris). Coarse large-R
   accept path when B1 `evaluateCommit` G3/G5 miss on few-band spans.
3. Wired from `commitPlanesA3` (not pre-B1) so merged arc shards are visible.

## Proof — handle-lock segmentation (`regiondump --diag`)

```
segment comp0 ok regions=27 planes=9 cylinders=18
  id=5  type=cylinder  R=20.124  tris=28  axis≈(0,1,0)
  id=6  type=cylinder  R=40.092  tris=18  axis≈(0,1,0)
  id=16 type=cylinder  R=16.000  tris=36  axis≈(0,1,0)
distinctRadii includes 20.12, 40.09, 16.00
```

vs `ground-truth.json` targets R=20 and R=30: **R≈20.1 recovered**; GT R=30 nominal
maps to mesh-measured **R≈40.1** (no STL strip fits R=28–32 with residual <1 mm).

## Proof — TrueForm build (was R2 revert)

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl -o /tmp/hl-f.step --engine trueform --quiet
```

```
smoothBuiltComponents=1  smoothRevertedComponents=0
smoothPlanes=9  smoothCylinders=18  facesAfterSmooth=225
ok=true  watertight=true  volumeDeltaPct=0.000000
```

## Proof — regression

```
ctest -R "p2buildtest|gates_smoke|corpus_engine"  → 13/13 pass
```

Corpus plane counts unchanged (coarse peel gated): S02 planes=6, S12-a planes=15,
S06 planes=2.

## Files

- `src/refit_math.cpp` — `detectLargeArcStrip`, helpers
- `src/refit_grow.cpp` — `peelLargeArcStripsA2b` in `commitPlanesA3`
- `src/refit_internal.hpp` — `ArcStripDetect` contract
