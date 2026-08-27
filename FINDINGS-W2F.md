# Lane W2FIX — wave-2 regression repair (hl/w2fix, 2026-08-26)

## Failures on entry (hl/wave2 + hl/H merge)

`ctest` 25/27:

1. **corpus_engine_convert** — Body11 smooth census `closed=False` (127 built cylinders OK).
2. **gates_full** — G0.3 hard FAIL on Body11, Body28, cube (`--force-sew --smooth` STEP ≠ off-path).

## Bisection (Body11 smooth, macOS arm64)

| commit | census closed | built cyl | G0.3 |
|--------|---------------|-----------|------|
| `95976bc` (pre-F) | True | 127 | — |
| `d9ff96c` (+F) | True | 127 | — |
| `f926786` (+G) | True | 127 | — |
| `5609aec` (+closure) | True | 127 | PASS |
| `499a2e1` (+J) | **False** | 127 | FAIL |
| HEAD (+J, pre-fix) | False | 127 | FAIL |

**Root cause:** lane J’s TrueForm `smooth-flat` second `UnifySameDomain` pass ran on
**every** `--smooth` run. On Body11 (15 300 tris) it over-merged faceted islands,
J6 left `freeEdges=164`, census `openShells=1`, and `--force-sew --smooth` STEP
diverged from off-path (G0.3).

Intended scope (FINDINGS-INT / `coarseFusionBand`): **500–1200 tris** coarse Fusion
exports (handle-lock @ 908 tris), not flagship real-CAD meshes.

## Fixes

1. **Merged hl/H** (`e6f323e`, `4eda8b0`) — pinned S09 corpus + G0.3 S09 exit exemption
   (already satisfied; S09 G0.3 PASS retained).
2. **Gate `smooth-flat`** in `src/stl2step.cpp` to `nTri ∈ [500, 1200]` (sync with
   `coarseFusionBand` in `refit_internal.hpp`).

## Proof (post-fix, macOS arm64)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTL2STEP_BUILD_TESTS=ON
cmake --build build -j 8
ctest --test-dir build --output-on-failure   # 27/27

./build/stl2step tests/corpus/Body11.stl -o /tmp/b11.step --quiet --smooth --no-verify
./build/stl2step_census /tmp/b11.step
# closed=true cylinder=127 openShells=0

./build/stl2step tests/diag/handle-lock/handle-lock.stl -o /tmp/hl.step \
  --engine trueform --quiet --no-verify
# facesAfterUnify=137 watertight=true smoothRevertedComponents=0

python3 tests/gates/run_gates.py --fixture Body11,Body28,cube --gate G0.3 \
  --binary ./build/stl2step --corpus tests/corpus --no-verify
# G0.3 PASS on all three
```
