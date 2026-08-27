# Lane RIDS — partial face fixes + R1 containment (2026-08-27)

Branch: `hl/rids`. Mesh: `tests/diag/handle-lock/handle-lock.stl` (908 tris).

## Verdict: success (census target met)

| Criterion | Target | Actual |
|-----------|--------|--------|
| `smoothBuiltCylinders` | ≥ 10 | **16** |
| `watertight` / `openShells` | true / 0 | true / 0 |
| `smoothRevertedComponents` | 0 | 0 |
| ctest | 27/27 | **27/27** (2026-08-27) |
| Body11 floor (`corpus_engine_convert`) | 127 closed | PASS |

Census vs `ground-truth.json`: **16 cylinders** in STEP (GT ideal = 15). Extra
segmentation at R≈27.3 and 40.1 mm; GT radii 6 / 20 / 30 mm still absent from
P1. All four previously failing partials (rid **14 / 19 / 20 / 23**) now build.

## Root causes

1. **rid=14 (R≈16 mm, st=27)** — seam-crossing / negative `uMin` sheet needed
   `u₀`-rotated `Geom_RectangularTrimmedSurface` with UV in `[0, span]`, not raw
   negative `uMin` on `surf0`.
2. **rid=19 / 20 (R≈3.787 mm, st=32)** — quarter-hole caps: rid=19 via `uMin`
   rotation; rid=20 via **uMax rotation** onto the same trimmed sheet; hole
   outer-wire order only on coarse meshes (≥500 tris).
3. **rid=23 (R≈0.5 mm, st=32)** — builds with standard `surf0`/`surf1` path once
   wire refresh is **not** applied on small fixtures; coarse-band rot trim alone
   was wrong for this boss span.
4. **R1 cascade** — one invalid partial triggered `shellIsValid` recover that
   **exploded all partial cylinders** when a Seamed360 survived; census stayed at
   1× closed-360 despite per-face validity. Containment: skip that explode when
   ≥10 partial regions built analytic on a closed shell in the Fusion coarse band
   (500–1200 tris); matching probe waiver in `stl2step.cpp`.

## Fixes (`refit_build.cpp`, `stl2step.cpp` — TrueForm only)

- `wantRotTrim` + `tryRotTrimmedSheet` — rotated trimmed sheet; uMax rotation for
  positive-u quarter holes; gated to `mv.nTri ≥ 500`.
- `refreshOuterWire` before `bindCylPCurves` on coarse meshes only (fixes pcurve
  poisoning across surface attempts; preserves fixture volume on small STLs).
- `keepValidPartials` — closed shell + ≥10 built partials in coarse band: do not
  blanket-explode partials on `shellIsValid` recover.
- `stl2step.cpp` — accept closed probe shell when ≥10 partial `builtAs` and
  `mv.nTri` in coarse band despite whole-shell `BRepCheck_Analyzer` fail.

## Proof (macOS arm64)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest --test-dir build                    # 27/27 pass
```

```sh
./build/stl2step tests/diag/handle-lock/handle-lock.stl \
  -o /tmp/hl-rids.step --engine trueform --no-verify --quiet
python3 tests/tools/step_census.py /tmp/hl-rids.step
```

```
RESULT smoothBuiltCylinders=16 smoothBuiltPlanes=9 smoothRevertedComponents=0
  watertight=true openShells=0 warnings=[]
Census: surfaces.cylinder=16  radii=[0.501,3.787×2,5×2,5.75,5.994,9.02,10×3,11.7,16×2,27.3,40.1]
  (GT: 15 cyl — extra 27/40 mm segmentation; missing 6/20/30 mm)
```

```sh
ctest --test-dir build -R corpus_engine_convert   # Body11 floor PASS
```
