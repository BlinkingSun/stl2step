# Lane PARTIALFACES — partial-cylinder face construction (2026-08-27)

Branch: `hl/partialfaces`. Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris).

## Root cause

`buildPartialCylinder` failed with `BRepCheck_BadOrientationOfSubshape` / `UnorientableShape`
even when trim vertices were on the fitted cylinder (`maxCylDev≈0`). The failure was **UV
pcurve direction**, not mesh chord sag:

1. **`bindCylPCurves` used 3D curve parameters** (`c3->Value(f/l)`) instead of wire traversal
   order. Loop-reversed cap chains (e.g. rid=13 ch=25) bound UV opposite to the wire → open UV loop.
2. **Generator lines** whose `Geom_Line` param direction opposed chain `ia→ib` (rid=13 ch=3)
   bound bottom→top in UV while the wire walked top→bottom.
3. **`makeArc` complementary branch** called `e.Reverse()`; `appendCollapsed(reversed=1)` reversed
   again on loop-reversed caps → double-reverse / wrong `Range` on partial cap circles.

Moving shared vertex TShapes onto the cylinder (Lane E chord-sag hypothesis) **breaks** arc
`Range`/vertex-parameter wiring and must not be done.

## Fixes (`refit_build.cpp`, TrueForm path only)

1. **`bindCylPCurves`** — `BRepTools_WireExplorer::CurrentVertex()` for UV endpoints; project
   endpoints with `projectPntOnCylinder` before `(u,v)` (mesh-inscribed → fitted circumscribed
   for linear edges only).
2. **`orientEdgeFromTo`** — after collapse, align **line** edges to chain `ia→ib`.
3. **`makeArc` complementary branch** — keep OCCT `e.Reverse()` ladder (plate_half_hole
   volume gate); UV closure for handle-lock rid=13 is fixed by (1–2), not by removing
   `Reverse()` (that regresses `p2buildtest_plate_half_hole`).
4. **`partialFaceTolCap`** — face tolerance from measured `chordSagitta` when `nSides>0`.
5. **UV-trimmed sheet** — `Geom_RectangularTrimmedSurface` + `makeFaceBound` retry with
   `rebindCirclePCurvesOnWire` for all partial cylinders; try reversed outer wire.
6. **`STL2STEP_COLLAPSE_DIAG`** — `DIAG_PARTIAL rid=N maxCylDev cap valid st=…` inventory.

## Proof (macOS arm64, this worktree)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest -R p2buildtest   # 11/11 pass
```

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-pf.step --engine trueform --quiet
```

**Before (CYLEDGES tip):** rid=13 `ensure-invalid`; `smoothBuiltCylinders=1`; 68/68 chains but
partial faces invalid.

**After:**

```
DIAG_COLLAPSE mix=0 none=0 fail=0 ok=68 total=68
DIAG_PARTIAL: rid=14,19,20,23 (rid=13 builds OK — no DIAG_PARTIAL)
RESULT smoothBuiltCylinders=1 smoothBuiltComponents=1 smoothRevertedComponents=0
  watertight=true volumeDeltaPct=0.000000 warnings=[]
```

```sh
python3 tests/tools/step_census.py /tmp/hl-pf.step
# surfaces.cylinder=1 (R=5.75 closed360 only in STEP census)
```

R1 still drops successfully-built partials when rid **6 / 20 / 23** fail and uncollapse shared
chains — census stays at 1 cylinder until those three build.

## Remaining gaps

1. **rid=14 / 19** — UV/orientation on multi-reversed loops (`st=27/32`).
2. **rid=20 / 23** — small-R complementary cap arcs (`st=32`).
3. **rid=6** — large-R partial (was failing intermittently during R1 rounds).
4. **Ground-truth census** (15 cylinders) blocked by R1 cascade from partial build failures.
5. **Body12** skew cyl|cyl — unchanged (CYLEDGES scope).
