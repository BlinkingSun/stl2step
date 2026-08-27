# Lane E — J6 closure / partial-cylinder arc edges (2026-08-26)

Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris). Branch: `hl/closure`
from `hl/integrate` (a095f01).

## Inventory (49 free edges @ J6 recover=2)

`STL2STEP_J6_DIAG=1` on this worktree:

```
DIAG_J6 summary freeEdges=49 analyticCollapsed=12 analyticPolyline=62
warning: J6: shell not closed freeEdges=49 faces=514 recover=2
```

Pattern (representative):

| Kind | count | detail |
|------|-------|--------|
| Seamed360 hole rim | ~1 | `freeE#0 len=0` degenerate vtx on plane rid=3; `ci=37` matches one rim edge on rid=1 |
| Faceted rim polylines | ~46 | `len≈0.752 mm` mesh chords on plane faces rid=0/1 (`ci=-1`, no collapsed chain) |
| Collapsed chain gap | ~2 | isolated edges where one side collapsed, neighbour still facet |

First-pass collapse (`STL2STEP_COLLAPSE_DIAG=1`):

```
DIAG_COLLAPSE mix=0 none=6 fail=2 ok=66 total=74 recover=0 rounds=0
```

Six chains still `curve=None` (all 2-vertex, `tangent=0`): cyl|cyl `13|8`, `15|4`, `18|17`;
plane|cyl oblique `8|9`, `4|5`, `17|6`. Two `MakeEdge` fails (closed360 cap + oblique arc).

**Root cascade:** partial cylinders **rid 14, 25, 26** fail `buildPartialCylinder` → R1
explode → mixed analytic/facet shell → 49 free edges → honest R2 revert.

`buildPartialCylinder` failure: `BRepCheck` invalid face (`ensure-invalid` / `face-invalid`);
wire vertex deviation up to **2.44 mm** vs cap **1.76 mm** (rid=14, R≈6 mm). IntAna
ellipse trim on coarse plane|cyl boundaries does not land on the face cylinder within
post-fit gate.

## Fixes landed (TrueForm path only)

1. **`partialCylCapArcMid` + `makeArc` for 2-vertex partial-cylinder cap circles**
   (u-span midpoint, not shortest OCCT arc).
2. **`makeFullCircle` for closed360 inner/cap loops** (plane holes around Seamed360).
3. **`intAnaAcceptResidual`** — relax `pickIntAna` gate on short chains using
   `maxVertexDev` (coarse Fusion band).
4. **`cylForIntersect`** in `intersectSurfaces` (u₀-rotated cylinder for partials).
5. **`bindCylPCurves`** — ellipse edges use `ShapeFix_Edge` pcurves, not linear UV chord.
6. **`makeEllipseArc` path** for open ellipse chains (2+ verts).
7. **Seamed360 inner-wire `Reversed()`** in `buildPlanarFace` when sharing full cap edge.
8. **`STL2STEP_J6_DIAG=1`** — `diagJ6FreeEdges` inventory (region/chain/length).

## Proof (macOS arm64, this worktree)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest -R "p2buildtest|gates_smoke"   # 12/12 pass
```

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-e.step --engine trueform --quiet
```

**RESULT (exit 2 — R2 revert unchanged):**

```
smoothBuiltComponents=0  smoothRevertedComponents=1
facesAfterSmooth=434  J6: freeEdges=49 faces=514 recover=2
```

```sh
python3 tests/tools/step_census.py --expect \
  tests/diag/handle-lock/ground-truth.json /tmp/hl-e.step
# EXPECT FAIL: 434 faceted planes (post-revert)
```

Segmentation still **12 planes + 15 cylinders** (regiondump `--diag`); build still
cannot close J6 on this mesh without partial-cylinder face validity.

## Remaining gaps (irreducible this lane)

1. **rid 14 / 25 / 26 `buildPartialCylinder`** — oblique plane|cyl trim wires;
   need exact ellipse-on-rotated-cylinder wiring or alternate partial-face builder.
2. **Six 2-vertex oblique/cyl|cyl chains** — IntAna miss; need constructed trim curves.
3. **Seamed360 hole ci=37** — plane inner wire TShape/orientation vs rid=3 cap (degenerate `len=0` edge).
4. **20 / 30 mm bores** — still absent from segmentation (Lane B growth).

Honest revert fires; no J6 waiver.
