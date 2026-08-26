# Body11 failure-class repros (f2-diag)

Minimal meshes that reproduce Body11 `--smooth` failure classes on
`main @ 6e3348d`. Drive with `python3 tests/diag/body11/run_repros.py`.

Canonical Body11 (not copied here):
`/Users/jroberts/Desktop/Internal Development/3D files/STL/Body11.stl`

| Mesh | Class reproduced today | Notes |
|---|---|---|
| `me_plnpln_transversal.stl` | `ME_PLNPLN_*_PROJ` — `smooth: analytic MakeEdge failed` | Shared-lattice 30 mm cube, 3×3 faces, 0.07 mm jitter. Multi-tri plane regions whose LS plane misses terminals by > `meshTolCap` (sewTol). |
| `me_cylpln_ellipse.stl` | `ME_CYLPLN_*` — same MakeEdge warn | 16-gon prism with jittered cap verts. Cap circles / oblique generators fail `MakeEdge`. |
| `ia_cylcyl_nogeom.stl` | `IA_CYLCYL_NOGEOM` — `smooth: IntAna cyl\|cyl empty/same` | Bent pipe, 12-gon, 5° kink. Axes not ∥ within P1's 3°; IntAna returns `NoGeometricSolution`. Also hits seamed360 invalid + J6 + R2. |
| `seamed360_capseam.stl` | `SEAMED360_BREPCHECK_INVALID` | Body11 c0 rid 0 wall (nSides=44, R≈4.9) + synthetic caps. Emits `seamed360: BRepCheck invalid on seamed face`. |

**Cap-seam** (`seamed360: cap wire does not pass seam vertex`) is 1 event on
Body11 c0 rid 0 (nSides=44, §2.4 rung 3 → NotBuilt). No smaller mesh emitted
that exact string; the seamed360 repro is the hole-editability fixture.

What happens today vs what should happen is in each `*.expected.json`.
The race builds to `_team/SPEC-F2.md`, not to “make Body11 work”.
