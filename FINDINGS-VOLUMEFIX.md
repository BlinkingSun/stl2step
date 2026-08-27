# FINDINGS-VOLUMEFIX — handle-lock volume blowout (2026-08-27)

## Gate RED (pre-fix, hl/wave3 @ 2e62079)

Registered `tests/corpus/handle-lock.{stl,expected.json}` with
`smoothVolumeDeltaPctMax=6.013` (k=3 budget: max(1e-4·15868.9, 3·318.13) mm³).

Pre-fix TrueForm run (waiver intact):

```
./build/stl2step tests/diag/handle-lock/handle-lock.stl -o /tmp/o.step --engine trueform --quiet
RESULT volumeDeltaPct=121.567807 stepVolumeMM3=44536.658346 meshVolumeMM3=15868.884516
       smoothBuiltCylinders=16 smoothRevertedComponents=0 watertight=true
warnings: "smooth: B-Rep invalid after build", "B-Rep volume differs from mesh volume beyond refit budget"
```

121.57% >> 6.013% → `corpus_engine_convert` FAIL. ctest had no handle-lock fixture.

## Root cause

Two hl/rids waivers suppressed the R2 safety net on the coarse-Fusion band (500–1200 tris):

1. **`keepValidPartials`** (`refit_build.cpp`): when shell closed + ≥10 partial cylinders,
   skipped BRepCheck-invalid recovery instead of exploding bad analytic faces.
2. **Probe waiver** (`stl2step.cpp`): when component probe failed `BRepCheck_Analyzer`,
   re-accepted the build if `nPartial≥10` in the 500–1200 tri band.

Together they shipped a **closed, watertight, BRepCheck-invalid** 25-face shell
(9 plane + 16 cylinder) with phantom volume:

| metric | mesh | pre-fix STEP | post-fix STEP |
|--------|------|--------------|---------------|
| volume mm³ | 15868.9 | 44536.7 (+181%) | 15864.7 (−0.03%) |
| faces | 908 tris | 25 analytic | 102 (9 pln + 1 built cyl + facets) |
| BRepCheck | — | INVALID | valid |
| `smoothBuiltCylinders` | — | 16 | 1 |

Extra ~28 700 mm³ came from **partial-cylinder sheets** whose UV-trimmed faces
extended beyond their boundary wires (and/or inverted orientation on hole walls)
while the shell still reported closed. `stl2step_census` on the bad STEP reported
`valid=true` at shape level but whole-shape volume was 2.8× the mesh.

`regiondump`: 25 regions, 16 cylinders (15 partial), `dVolAbs=318.13`, budget=954.39 mm³.

## Fix

1. **Removed `keepValidPartials`** — BRepCheck-invalid closed shells now enter normal
   J6 recovery (explode invalid analytic faces / partial cylinders).
2. **Removed probe waiver** — component probe must pass `BRepCheck_Analyzer`.
3. **Probe D4.5 volume gate** — after valid+closed probe, reject if
   `|shellVol − meshVol| > max(1e-4·|meshVol|, 3·Σ|dVolPredicted|)`.
4. **G0.3** — `smooth-flat` unify gated on `!forceSew` so handle-lock passes
   force-sew global refit-disable.
5. **Corpus fixture** `handle-lock.expected.json` with live[] census radii from
   `ground-truth.json`, `builtCylindersFloor=0`, volume budget 954.39 mm³.
6. **`waiver_audit.py`** — positive control accepts volume-safe revert OR valid analytic;
   negative control (rid=16 tri delete) still not rescued.

## Proof (post-fix, hl/volumefix)

```
./build/stl2step tests/corpus/handle-lock.stl -o /tmp/hl.step --engine trueform --quiet
RESULT volumeDeltaPct=0.0 stepVolumeMM3=15864.701 meshVolumeMM3=15868.885
       smoothBuiltCylinders=1 smoothRevertedComponents=0 watertight=true warnings=[]
```

- `corpus_engine_convert`: OK handle-lock (verbatim + smooth)
- `waiver_audit.py`: PASS (positive + negative)
- `ctest`: **27/27 GREEN** (incl. gates_full G0.3 on handle-lock)
- Body11: `smoothCylinders=419 built=0 reverted=2 vol=0%` — unchanged

Correctness over count: 1 built cylinder vs ideal 15; volume-true + valid wins.
