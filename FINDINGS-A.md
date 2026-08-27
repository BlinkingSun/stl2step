# Lane A — segmentation truth (handle-lock)

Date: 2026-08-26. Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris, Fusion STLB ATF 15.8).

## Contradiction (reported)

`stl2step_regiondump` was said to show **0 regions** for comp0 while
`--engine trueform` attempted a 519-face analytic rebuild (J6: 49 free edges)
then reverted. That hid whether segmentation was empty vs buildFaces/J6 failure.

## Root causes

1. **Misread RESULT vs RegionSet.** After R2 revert, RESULT reports
   `smoothCylinders=0` / `smoothPlanes=0` — that is post-revert build stats,
   not `refit::segment()` output. Segmentation did run and produced regions;
   `buildFaces` expanded them to 519 faces before J6 failed.

2. **Half-link in regiondump CMake (FINDINGS-0 D1).** `stl2step_regiondump`
   recompiled all five P1 TUs into the executable while `stl2step::harness`
   already links `stl2step::core` (same `refit::segment()`). Duplicate symbols
   let the linker pick the wrong TU (stub or stale object) vs the engine.

3. **Silent dirty skip.** Dirty components were omitted from `comps[]` entirely
   (`dumped=0`), which reads like “no segmentation” even when `components: N`
   is non-zero. handle-lock comp0 is clean (`open=0 conflict=0 nonManifold=0`).

Mesh preprocessing (weld/split/`MeshView`) matches between harness and engine
(same checksums on triEdges/triDirs). Segment params match defaults:
`epsPlane=0` (auto), `thetaPlaneDeg=2`, `doFillets=true`.

## Fixes (this lane)

- **CMake:** link `refit::segment` from `stl2step::core` only; no duplicate P1
  compile when core + `src/refit_segment.cpp` exist.
- **regiondump:** emit dirty rows with `segmentSkip: "dirtyComponent"`,
  `regionCount` / `totalRegions`, `--include-dirty`, `--diag` stderr summary,
  `--no-smooth-fillets` alias.
- **Engine:** `STL2STEP_SEGMENT_SUMMARY=1` stderr dump after `refit::segment()`.

Recognizer gates unchanged (Lane B).

## Proof (macOS arm64, this worktree)

```text
$ ./build/stl2step_regiondump tests/diag/handle-lock/handle-lock.stl --diag
segment comp0 root=30 ok regions=67 rejected=2 planes=56 cylinders=11 fillets=0 facetIslands=3
  cylinder radii:5.75,9.99999,10.3338,10.0001,5.00001,5.99359,5,9.01954,3.78713,3.78715,0.499774

$ STL2STEP_SEGMENT_SUMMARY=1 ./build/stl2step tests/diag/handle-lock/handle-lock.stl \
    -o /tmp/o.step --engine trueform --quiet 2>&1 | head -3
engine segment root=30 regions=67 rejected=2 planes=56 cylinders=11 fillets=0 facetIslands=3
```

Envelope: `totalRegions=67`, `regionCount=67`, `segmentOk=true`, `clean=true`.

Ground-truth CAD: 28 faces (13 planes + 15 cylinders). Segmentation finds
67 accepted regions (56 micro-planes + 11 cylinders) — oversegmentation, not
empty segmentation. Missing radii vs GT (6, 9, 10×3, 15, 16, 20, 30 mm) are
recognizer/gate issues (Lane B). J6 / IntAna revert is downstream of P2 build.
