# Body11 — known gap after edge-failure J6 deferred R1

Measured on the small component after the **chainEdgeFail ledger** (class-filtered,
deferred to J6 recoverPass 0): **127** analytic cylinders survive in-region census;
J6 **freeEdges** drops **164 → 12**; shell still open.

## Residual 12 free edges

Six **plane|plane** IntAna lines remain as mesh polylines on otherwise-built
plates. Three **copy-fallback** plates fail BRepCheck when built via the
`makeFaceKeep` path:

| Region id | BRepCheck | Wire walk stall |
|-----------|-----------|-----------------|
| 571 | `UnorientableShape` / `SelfIntersectingWire` | 1/7 |
| 707 | `UnorientableShape` / `SelfIntersectingWire` | 1/5 |
| 743 | `UnorientableShape` / `SelfIntersectingWire` | 4/6 |

Each failure is a multi-inner plate where the outer wire is valid under
`makeFaceKeep` but the assembled face is not. The stalled walk counts are
edges visited vs edges in the wire before `BRepTools_WireExplorer` stops.

## Banned copy fallbacks

These paths **copy TShapes** and relocate the J2 defect — they must not ship as
the repair:

- `BRepBuilderAPI_MakeFace(plane, wire)` on refit wires (projection copies edges)
- `ShapeFix_Face` on refit wires (same copy semantics; measured S03/S10 shell opens)

`makeFaceKeep` (Builder `MakeFace` + `Add` wire) is the only legal face constructor
for mixed analytic plates.

## Ruled-but-unlanded repair recipe

| Step | Action |
|------|--------|
| **A1** | Per-wire UV signed-area orientation via `Oriented()` flag flip before `Add` |
| **A2** | `addPcurvesOnFace` before BRepCheck validity gate |
| **A3** | Status-routed at-source wire repair (route on `BRepCheck_Status`, fix wire in place — no copy) |

A0 diagnostic: `STL2STEP_DIAG_PLATES=1` dumps `BRepCheck_Status` after a failing
`makeFaceKeep` (stderr only; zero default-path effect).
