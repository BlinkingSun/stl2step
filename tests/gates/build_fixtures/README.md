# build_fixtures — P2 adversary RegionSet fixtures

Hand-authored `RegionSet` JSON + reference STL meshes + `*.expected.json` ground truth for the P2 `buildFaces` gate. **Data only** — no C++ target. `p2-harness` loads each `<name>/<name>.stl` via `stl2step::harness::loadMesh()` and pairs it with `<name>.regionset.json`. Every mesh is **one manifold component**; JSON local ids are comps[0] ids.

Regenerate: `python3 generate_fixtures.py` (rewrites all eleven dirs). Validate: `python3 validate_regionset.py`, `python3 check_invariants.py` (opens every `*.expected.json`), `python3 negative_tests.py`.

`check_invariants.py` enforces the full I7/I7b contract (every chain touching a region appears in exactly one loop; `closed360` CapLow and CapHigh are distinct chains at vMin/vMax), `triRegion[t] == regions[id].tris` id consistency, I8 splits at ≥3-owner and owner-pair junctions (including open/non-manifold vertices of degree ≠ 2), single manifold component, `conflict==0` (open>0 is allowed), G1 fillet|neighbour-plane tangent chains (origin==filletStrip; 1-edge long sides legal), D4.5 `volumeBudgetMM3 = max(1e-4·|vol|, 3·Σ|dVol|)` with signed sum reported separately, and the R8 un-collapse block on `r8_360_explode_caps`.

## Fixtures

| Dir | What it breaks | Builder must |
|---|---|---|
| `full_360_hole` | 360° seam ladder rung 1 | One `Seamed360` cylindrical face; caps are `capLow`/`capHigh`, no `outer` |
| `plate_half_hole` | Isolated ~180° partial on a plate (not S05 slot) | `builtAs=single` cylinder + plane; one `outer` loop on partial |
| `slotted_stadium` | Two ~180° partials + connecting flats | Four analytic faces; both cylinders `closed360=false` |
| `quarter_round_1tri` | D7 fillet radius (withdrawn median rule → R/2) | Fillet `radius=2.0` in ground truth; three faces |
| `fillet_strip_2tri` | 2-tri strip, tangent long chains | Constructed generators (`tangent:true`); not intersection |
| `counterbore` | Periodic-U + inner wire on seam | Two `Seamed360` cylinders + annular plane; inner wire on inner cyl |
| `s09_mixed` | Analytic plane abutting faceted island | Mesh chords verbatim on mixed chain; edge tol absorbs sagitta |
| `r1_success` | MakeFace fails, R1 explode, shell stays closed | `explodedToFacets`; facet count = triangle count; 1 round |
| `r1_round2_cascade` | Neighbour rebuild after explode | Converges on round 2; all regions explode |
| `r2_forced` | nSides=3 cascade vs nSides=4 hole | Three stacked nSides=3 cylinders (sagitta/R=0.5) explode; distinct from `full_360_hole` nSides=4 Seamed360 |
| `r8_360_explode_caps` | R8 cap-ring un-collapse + inner wire rebuild | Each cap is one `closedLoop` chain of `nSides` chords (the `gp_Circ` collapse); explode un-collapses to those N chords; plane inner **is** CapHigh |

Each directory contains `<name>.stl`, `<name>.regionset.json`, and `<name>.expected.json` (face census, `builtAs`, D4.5 volume budget + signed/`dVolPredAbs`, R-ladder rung/rounds). `r8_360_explode_caps.expected.json` carries an `r8` block the checker verifies.
