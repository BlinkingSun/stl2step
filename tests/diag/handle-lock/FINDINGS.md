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

## Verbatim unify-angle sweep (macOS arm64, main @ hl/D worktree)

Coarse Fusion STL (908 tris); `facesAfterUnify` from RESULT (`--no-verify`):

| `--unify-angle` | `facesAfterUnify` |
|-----------------|-------------------|
| 0.001 (default) | 434 |
| 0.01            | 263 |

Default gate leaves one planar facet per mesh triangle pair that passes the
0.001° coplanar test; float32 STL normal jitter blocks most merges. Raising
to 0.01° collapses many coplanar facet pairs (263 faces) but is **not** a
Verbatim default change (G0.1).

STEP text census of Verbatim output (`tests/tools/step_census.py`):
434 `ADVANCED_FACE`, 434 `PLANE`, 0 cylinders, 888 `LINE` edges.

## Status

`ground-truth.json` sidecar added (28 faces = 13 planes + 15 cylinders per
original CAD). Not wired into ctest yet. Verbatim output is correct and
shippable for this part; TrueForm work item is open.
