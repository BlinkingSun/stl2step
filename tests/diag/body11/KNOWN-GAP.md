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

## I-checker I7 on large real-part RegionSet dumps (Body28 comp0)

Live `stl2step_regiondump --bare` on Body28 exposes **I7 loop incompleteness** on
two open regions — not a checker false positive:

| Region | Loop | Chain | Symptom |
|--------|------|-------|---------|
| 547 | loop0 | 1124 | single open chain; loops not complete/closed |
| 783 | loop0 | 1031 | single open chain; loops not complete/closed |

Body11 comp0/comp1 pass I-checker today; Body28 comp0 fails. Both real-CAD fixtures
are **I-checker PARKED** in `gates_full` (honest red via `--unpark I-checker`);
all synthetic corpus fixtures remain LIVE.

## handle-lock plate validity → census unlock — PARKED, v-next

**Registered 2026-08-29 by D3 (`_team/reports/DECISION-census.md`) after two unanimous
all-four-fail races (8 lane-attempts) on the same criterion.** Phase `ac1` closed at
**census 1, volume-true**. This is the same defect class as 571/707/743 above, reached from the
**in-place** side rather than the copy side.

handle-lock ships `smoothBuiltCylinders=1` (rid=3, R5.750001, GT-matched), volΔ 0.000%,
closed + BRepCheck-valid + watertight, `warnings=[]`, via the U2 blanket. Recognition finds 16
cylinders; **max GT-distinct is 13** (D2 §0.6 — rid=6 off 1.234%, rid=16 off 2.637%, rid=5 barred
by D1 §M1). The gap 1 → 13 is blocked by planar-plate face validity, not by the ladder.

**Both routes out of `makeFaceKeep` are now measured as defective:**

| route | status | measurement |
|---|---|---|
| copy (`MakeFace(plane,wire)`, `ShapeFix_Face`) | BANNED — **relocates** J2 | S03/S10 shell opens |
| in-place (`GeomAPI::To2d` + `UpdateEdge` + tol absorb) | **PARKED — propagates J2** | J6 free-edges Body9 8→20, **Body12 49→82 (red-line 81)**, Body18 14→18, Body20 18→22 |

**Open chain, in dependency order** (measured, `.stl2step-wt/ac-platefix-r3`):

1. **Shared-TShape propagation (root).** Plates 0/1 are `DIAG_PLATE valid=1` at ~3 µm residual,
   then `DIAG_CASCADE culprit faceValid=0` at the shell walk. Not an ordering bug — the same
   writes regress J6 on four bodies sharing no topology with these plates. `B.SameParameter` /
   `SameRange` restore-flags are insufficient.
2. **rid=2** — `BRepCheck_Analyzer IsValid=0` with an **empty status list**
   (`intersect=0 classify=0 orient=0`, `nBadCtx=0`, all edge d ≤ 0.00013). Needs a subshape walk
   beyond `Status`/`StatusOnShape`.
3. **rid=4** — `st=32`, pcurve UV signed area (−53.21) disagrees in sign with the 3D-plane area;
   wire-child reverse does not reach `BRepCheck`. `SameParameter` on shared edges caused a D4.5
   revert in an earlier iteration — do not reach for it.
4. **Cylinders 7 (R≈29.999, the only R30) and 14** stay `faceValid=0` after plates 0/1/9/10/21 are
   isolated-valid; `nbrs(7)=[0,1,2,11,16,18]`, `nbrs(14)=[1,2,4,18]`. Suspected downstream of (1);
   **never measured as such.** U0 explodes 7 → R30 absent from every shipped census to date.
5. **`RULE 2.2` per-region residual, ANSWERED by AC2-S4, adopted scoped by D4 §3.** |vface − vchord| = 235.4 / 85.4 / 167.7 /
   148.2 / 6.1 mm³ on rids 0/1/2/4/21 against the old `sub=0.176321` (≈643 mm³ ≈ 4.1 % of 15868.885).
   Scoped form: `sub = max(floor, 3·|dVolPredicted|, (plane && |dVol|<1e-6) ? 1.05·|resid| : 0)`.
   Plate `DIAG_CASCADE resid` is now truthful (`pass=1`); cylinder landmine class unchanged by construction.
   Face-validity race stays PARKED (D4 §4).
6. **Fit accuracy (`AC-FIT20`, pre-authorized by D2 §4.4, never dispatched).** rid=6 → 20.000 and
   rid=16 → 15.000 within 0.3 % via growth `R_ref`. Sequential after the above; not measurable
   until the plates ship.

**Standing constraints for any v-next lane** (carried, not re-litigated): `GT_REL_TOL=0.003` is not
loosened; no rid special-casing; no shared-edge tol ≥ the S09 poison threshold; no ladder or gate
edits; `src/stl2step.cpp:614-652` (D4.5 probe) frozen at
`ddf4be5391d662ceabe39c9bac506f9c86e506a345a3762dec31331e96169fd6`; `addPcurvesThenFace` stays
dead; `hl-ratchet.json` floor stays **1** until a fully green unlock.

**Prior art — start here, do not re-derive:** `ac/platefix-r3` (A3 status routing, 5/7 plates at
source, `DO-NOT-MERGE`) and `ac/platefix-r1` (sliver A2/A1; its 10.2 mm gap / ΔV 5601 mm³ exact-bind
is a **dead branch** — r3 achieves the same faces at ~3 µm for free; do not re-attempt it).
Full chain: `DECISION-plates.md` (D2) → `adjudicate-platefix.md` → `DECISION-census.md` (D3).
