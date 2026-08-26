# Lane D — ground-truth compare harness (2026-08-26)

## Deliverable

`tests/tools/step_census.py` — pure-text STEP Part 21 parser (no OCCT). Counts
`ADVANCED_FACE` entities, classifies face surfaces (`PLANE`, `CYLINDRICAL_SURFACE`,
`CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`, `B_SPLINE_*`), classifies
`EDGE_CURVE` 3D curves (`LINE`, `CIRCLE`, `B_SPLINE_CURVE` via `SURFACE_CURVE` /
`SEAM_CURVE` indirection), extracts sorted cylinder radii, emits JSON. Multiple
inputs → JSON array on stdout + human diff table on stderr. `--expect <json>`
asserts against a ground-truth sidecar (exit 1 on mismatch).

## Self-test vs OCCT census (`build/stl2step_census`)

| fixture | py faces | occt | py plane/cyl | occt | py line/circle | occt |
|---------|----------|------|--------------|------|----------------|------|
| S01 verbatim (`/tmp/S01.step`) | 6 | 6 | 6 / 0 | 6 / 0 | 12 / 0 | 12 / 0 |
| S06 `--smooth` (`/tmp/S06-smooth.step`) | 3 | 3 | 2 / 1 | 2 / 1 | 1 / 2 | 1 / 2 |

S06 smooth cylinder radius: py `[25]`, occt `24.999999996806` (6 dp round).

## handle-lock ground truth

`tests/diag/handle-lock/ground-truth.json` — 28 faces (13 planes + 15 cylinders),
30 circles, radii per COMMON.md.

Current stl2step output census (verbatim & trueform revert — identical):

```
faces=434 plane=434 cylinder=0 circle=0 line=888
```

`--expect ground-truth.json` correctly fails (434 vs 28 faces).

Unify-angle sweep added to `tests/diag/handle-lock/FINDINGS.md`:
0.001° → 434 faces; 0.01° → 263 faces.

## Usage

```sh
python3 tests/tools/step_census.py part.step
python3 tests/tools/step_census.py a.step b.step        # + diff table on stderr
python3 tests/tools/step_census.py --expect tests/diag/handle-lock/ground-truth.json out.step
```
