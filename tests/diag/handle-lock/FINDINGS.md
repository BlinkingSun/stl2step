# Stress case: handle-lock (2026-08-26, Josh's first live test file)

Fusion export (`STLB ATF 15.8.0.0`), 908 tris / 454 verts — very coarse.
Source: `~/Desktop/Handle lock.stl` (copied here verbatim).

## Behavior (main @ 83733e4 era build, macOS)

`--engine trueform`: converts OK (exit 0, watertight, volumeDeltaPct
0.000000) but the analytic rebuild REVERTS — output is Verbatim-quality:

- `smoothBuiltComponents 0, smoothRevertedComponents 1`
- warnings: `smooth: IntAna plane|cyl empty/same — keeping mesh polyline`;
  `J6: shell not closed freeEdges=49 faces=519 recover=2`;
  `analytic rebuild reverted on one component — kept faceted`

Parameter-insensitive: `--smooth-angle 8|15|25` × `--smooth-tol 0.5|1.0`
all revert identically → structural, not gate tuning.

## The bug-class finding

`stl2step_regiondump` reports **0 regions** for comp0, yet the in-engine
pass clearly attempted a rebuild reaching J6 with 519 faces. Two
segmentation entry points disagree on the same mesh — reconcile first;
whichever is wrong hides the real repair surface. Then: why does the
plane|cyl IntAna edge rebuild degenerate (49 free edges) on a mesh this
coarse — likely near-tangent plane/cylinder pairs at 8–16-segment
resolution where facet normals straddle the cyl gates.

## Status

Not wired into ctest (no ground-truth sidecar yet). Verbatim output is
correct and shippable for this part; TrueForm work item is open.
