# Integration findings — hl/integrate (2026-08-26)

Adjudication per `LANE.md`. Branch `hl/integrate` merges lanes D, A, C in full and
lane B segmentation selectively.

## Commits merged (prior to B)

| Lane | Commit   | Scope |
|------|----------|-------|
| D    | f1c6594  | `tests/tools/step_census.py`, handle-lock `ground-truth.json` |
| A    | 9692bd4  | regiondump truth, `--diag`, CMake half-link fix |
| C    | 8d3f78a  | `planeCylSideContact` + mesh-aware `g1Tangent`; cap/generator fallbacks |

## Lane B — merged (segmentation only)

**Taken in full** (scoped to `coarseFusionBand`: 500 ≤ nTri ≤ 1200, handle-lock
band; corpus S02=412 and S13/S14≤20 stay on default gates):

- `adaptCoarseSegmentParams` — `thetaPlaneDeg`≥15°, `thetaCylHiDeg`≥70°
- `evaluateCommit` G2 — chord-sagitta + `epsCylRing` for ≤8-tris patches
- `evaluateCommit` G3 — skip when arc span < 0.35 rad
- `evaluateCommit` G5 — partial span 30° + `g5Micro` (R≥5 mm, span≥12°, nBands≥2)
- `estimateNBandsFromPatch` — dihedral-median D2 fallback
- `commitPlanesA3` — coplanar UF merge (10° on coarse shards)

Files: `src/refit_grow.cpp`, `src/refit_segment.cpp`.

## Lane B — stripped (adjudication)

| Rejected item | Rationale |
|---------------|-----------|
| `stl2step.cpp` open-shell accept when STEP contains cylinders | Violates watertight contract; open shell never ships |
| `refit_build.cpp` J6 “coarse keep” (return built cylinders on J6 fail) | Same — relaxation of J6 closure |
| `planeNearPerpendicular`, `planeCapCylinder`, `capCircleOnCylinder` | Overlap C's principled side-graze + cap fallbacks; **C wins** |
| Looser `pickIntAna` accept (`sewTol*100`) | Overlap C's mesh-aware tolerance path |
| Global `kG5SpanPartialDeg` 40→30 in `refit_internal.hpp` | Replaced by coarse-only 30° in `evaluateCommit` |
| `g1Tangent` revert to `epsPlane`-only | Would undo lane C |
| B's `adaptCoarseSegmentParams` for all nTri≤1200 | Narrowed to 500–1200 band to protect corpus fixtures |

## Proof (a) — build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
→ exit 0, all targets built
```

## Proof (b) — handle-lock TrueForm

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-integrate.step --engine trueform --quiet
```

**RESULT (exit 2 — written with warnings, R2 revert):**

```
smoothBuiltComponents=0  smoothRevertedComponents=1
smoothPlanes=0  smoothCylinders=0  smoothBuiltCylinders=0
openShells=0  watertight=true  facesAfterSmooth=434
warnings: … IntAna plane|cyl … ×5; J6: shell not closed freeEdges=49 faces=479 recover=2;
          analytic rebuild reverted on one component -- kept faceted
```

Segmentation **does** run (regiondump `--diag`):

```
segment comp0 ok regions=27 planes=12 cylinders=15 fillets=0 rejected=0
distinctRadii=10  (0.50, 3.79×2, 5×2, 5.75, 6.0, 9.02, 10×3, 11.70, 13.87, 16.69 mm)
```

J6 still cannot close (49 free edges at plane|partial-cyl chains). Honest R2 revert
fires — output STEP remains faceted. **No closure weakening applied.**

**Census vs ground truth:**

```sh
python3 tests/tools/step_census.py --expect \
  tests/diag/handle-lock/ground-truth.json /tmp/hl-integrate.step
```

```
EXPECT FAIL: faces 434 vs 28; cylinder 0 vs 15; circle 0 vs 30
(post-revert faceted output — expected until J6 closes)
```

Target ground truth: 28 faces (13 planes + 15 cylinders), 30 circles.

## Proof (c) — regression

```sh
ctest --test-dir build
→ 26/27 passed
```

**Failure (pre-existing on main @ 92c2905, not introduced by B merge):**

```
p2real_live: S13.stl comp0 censusOk=False (exploded360=1);
             S14.stl comp0 buildFaces=False
```

Lane-specific re-gates after B+C integration:

```
ctest -R "p2buildtest|gates_smoke|corpus_engine"  → 13/13 pass
ctest -R "p2buildtest|gates_smoke"                 → 12/12 pass
```

## Remaining gaps (ranked)

1. **J6 shell closure on handle-lock** — 49 free edges after 479-face analytic
   attempt; mixed analytic + partial-cylinder polylines need arc-length MakeEdge /
   follow-on P2 work. Segmentation now reaches 12+15 regions (was 0 built).
2. **Large bores missing** — 20 mm / 30 mm radii not segmented (too few facets in
   arc; need wider growth through 90° caps).
3. **p2real_live S13/S14** — corpus gate failures on main baseline; unrelated to
   handle-lock lanes but blocks 27/27 on this box.
4. **IntAna plane|plane / cyl|cyl / remaining plane|cyl misses** — several chains
   still polyline after C fix (reg 18|17 `tangent=false` in regiondump).
