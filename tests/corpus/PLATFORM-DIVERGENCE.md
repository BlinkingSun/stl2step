# Corpus / gate platform divergence register

Precedent: **S09** — `gen_corpus.cpp` `isPinnedCorpusFixture` + committed
`S09.stl` / `S09.expected.json`. OCCT `BRepMesh` of the loft is not
byte-stable across libm/FMA (macOS 7.9.3 vs Windows 8.x). The pinned pair
is the calibration input. The generator must not overwrite it. macOS/Linux
assertions stay exactly as committed; Windows consumes the same bytes.

## Body11 / Body28 file-census cylinders (OCCT 8 / MSVC)

| fixture | assertion | macOS / OCCT 7.9.3 | Windows / OCCT 8.0.1 (measured 2026-08-30) |
|---|---|---|---|
| Body11 `--smooth` | `live[].builtCylindersFloor` | 0 (sidecar) | **0** (`OK Body11 census cylinders=0 planes=10798 valid closed`) |
| Body28 `--smooth` | same | 0 (sidecar) | **0** (`OK Body28 census cylinders=0 planes=8468 valid closed`) |

P1 still recognises cylinders (`smoothCylinders` 419 vs sidecar 583 on
Windows; 340 vs 344 on Body28). File-truth census after OCCT 8 write/read
is **zero analytic cylinders**. That is the documented S09-class numeric
divergence. **Do not raise `builtCylindersFloor` on Windows** and **do not
lower any macOS/Linux gate**. The floor is already 0 on every platform;
`run_engine_check.sh` treats it as never-get-worse. A future Windows-only
expected block would go in the sidecar as `live[].windows` / a waiver
sidecar **only if** a platform started asserting a *higher* floor than
another — not the case today.

## handle-lock `prism_build_unit` (not this class)

API solid on Windows is healthy (35 faces / 17 cyl pre-unify, vol exact).
The 2026-08-30 red was a hardcoded `/tmp/` write path (class a), not
census divergence.
