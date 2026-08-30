# Lane PARTIALFACES — partial-cylinder face construction (2026-08-27)

Branch: `hl/partialfaces`. Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris).

## Verdict: partial success (target not met)

Lane success criterion was `smoothBuiltCylinders >= 10` with census matching
`tests/diag/handle-lock/ground-truth.json` (15 cylinders). **Actual:
`smoothBuiltCylinders=1`** — only the closed-360 bore (rid=3, R=5.75 mm) survives
R1. That is honest partial progress, not closure.

What **did** land: the proven UV/pcurve core for partial-cylinder faces. rid=13
(R≈6 mm) no longer hits `ensure-invalid`; four other partials still fail face
validity and trigger R1 uncollapse of every shared partial.

## Root cause (confirmed)

`buildPartialCylinder` failed with `BRepCheck_BadOrientationOfSubshape` /
`UnorientableShape` even when trim vertices were on the fitted cylinder
(`maxCylDev≈0`). The failure was **UV pcurve direction**, not mesh chord sag
(Lane E vertex projection onto the cylinder breaks arc Range wiring — do not do it).

1. **`bindCylPCurves`** used 3D curve parameters instead of wire traversal order.
   Loop-reversed cap chains bound UV opposite to the wire → open UV loop.
2. **Generator lines** on partial-cylinder chains: `Geom_Line` param direction can
   oppose loop walk; UV must follow wire order.
3. **`makeArc` complementary branch** keeps OCCT `e.Reverse()` (plate_half_hole
   volume gate); handle-lock rid=13 is fixed by (1–2), not by removing `Reverse()`.

## Fixes (`refit_build.cpp`, TrueForm path only)

1. **`bindCylPCurves`** — `BRepTools_WireExplorer::CurrentVertex()` for UV endpoints;
   `projectPntOnCylinder` on endpoints before `(u,v)` (mesh-inscribed → fitted
   circumscribed for linear edges only).
2. **`orientEdgeFromTo`** — scoped to **partial-cylinder** chains only
   (`Cylinder && !closed360` on either side). Global line orient regressed S09
   comp1 volume (`p2real_live`); partial-only orient preserves rid=13 and corpus gates.
3. **`partialFaceTolCap`** — face tolerance from measured `chordSagitta` when `nSides>0`.
4. **UV-trimmed sheet** — `Geom_RectangularTrimmedSurface` + `makeFaceBound` retry with
   `rebindCirclePCurvesOnWire`; try reversed outer wire.
5. **`STL2STEP_COLLAPSE_DIAG`** — `DIAG_PARTIAL rid=N maxCylDev cap valid st=…` inventory.

## Per-rid partial-cylinder inventory (handle-lock, 16 recognized)

| rid | R (mm) | closed360 | Face build | BRepCheck | Notes |
|-----|--------|-----------|------------|-----------|-------|
| 3 | 5.750 | yes | **OK** | — | Only cylinder in STEP census |
| 5 | 16.012 | no | OK | — | Dropped by R1 when 14/19/20/23 fail |
| 6 | 27.329 | no | OK | — | same |
| 7 | 40.092 | no | OK | — | same |
| 8 | 10.000 | no | OK | — | same |
| 11 | 10.000 | no | OK | — | same |
| 12 | 5.000 | no | OK | — | same |
| **13** | **5.994** | no | **OK (fixed)** | — | Was `ensure-invalid`; UV fix proven |
| **14** | **16.000** | no | **FAIL** | st=27 | BadOrientationOfSubshape |
| 15 | 5.000 | no | OK | — | R1 cascade |
| 16 | 11.703 | no | OK | — | R1 cascade |
| 17 | 10.000 | no | OK | — | R1 cascade |
| 18 | 9.020 | no | OK | — | R1 cascade |
| **19** | **3.787** | no | **FAIL** | st=32 | UnorientableShape |
| **20** | **3.787** | no | **FAIL** | st=32 | UnorientableShape |
| **23** | **0.501** | no | **FAIL** | st=32 | UnorientableShape; smallest bore |

R1 still uncollapses **all** partial cylinders sharing chains with any build failure
among rid **14 / 19 / 20 / 23** — census stays at 1 cylinder until those four build.

## Proof (macOS arm64, this worktree)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest --test-dir build                    # 27/27 pass (2026-08-27)
```

### Handle-lock RESULT + census

```sh
export STL2STEP_COLLAPSE_DIAG=1
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-pf.step --engine trueform --quiet
python3 tests/tools/step_census.py /tmp/hl-pf.step
```

```
DIAG_COLLAPSE mix=0 none=0 fail=0 ok=68 total=68
DIAG_PARTIAL rid=14 … st=27 ensure-invalid
DIAG_PARTIAL rid=19 … st=32 ensure-invalid  (×3 attempts)
DIAG_PARTIAL rid=20 … st=32 ensure-invalid  (×2)
DIAG_PARTIAL rid=23 … st=32 ensure-invalid  (×2)
DIAG_COLLAPSE mix=53 … recover=1            # R1 uncollapse
RESULT smoothBuiltCylinders=1 smoothBuiltComponents=1 smoothRevertedComponents=0
  watertight=true openShells=0 volumeDeltaPct=0.000000 warnings=[]
Census: surfaces.cylinder=1  radii=[5.750001]   (GT: 15 cyl, 13 planes)
```

### Body11 floor

`corpus_engine_convert` (ctest #8) **PASS** — Body11 `smoothBuiltCylinders=127`,
`closed=true`, `openShells=0` intact.

### Body12 digest (unchanged vs CYLEDGES scope)

```sh
./build/stl2step "$STL2STEP_PRIVATE_CORPUS/Body12.stl" \
  -o /tmp/b12-pf.step --engine trueform --quiet
shasum -a 256 /tmp/b12-pf.step
# 70c721d503f2af723d1397b2e6e001dd748da43b6e4dae747e5a25b1c8af68d1
# smoothRevertedComponents=1; IntAna cyl|cyl warnings=108; J6 hits=2
# Skew cyl|cyl — CYLEDGES follow-on, not this lane.
```

## Remaining gaps (next lane)

1. **rid=14** — multi-loop UV/orientation (`st=27`).
2. **rid=19 / 20** — small-R complementary cap arcs (`st=32`).
3. **rid=23** — R≈0.5 mm partial (`st=32`).
4. **R1 cascade policy** — successfully-built partials (incl. rid=13) are discarded
   when unrelated partials fail; census cannot grow until all four build or R1 narrows.
5. **Ground-truth census** (15 cylinders) blocked until (3–4).
