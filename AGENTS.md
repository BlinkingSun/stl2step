# AGENTS.md — integrating the stl2step engine

This file is written for coding agents (and developers) wiring the **stl2step**
engine into other software. It is self-contained: read this and you can integrate
correctly without reading the source. Keep it in the repo root so agents discover
it.

---

## 1. What this engine does (one paragraph)

stl2step converts a **triangle mesh** (`.stl`, binary or ASCII) into a
**parametric B-Rep solid** (`.step`). It welds vertices, splits the mesh into
manifold bodies, optionally runs **stage 3.5 (TrueForm segment)** when
`Options::smooth` is true (`refit::segment()` on each clean component; dirty or
`--force-sew` components are skipped), then builds each body as an exact B-Rep
in parallel: a component with a refit plan takes the `refit::buildFaces()`
branch (analytic `Geom_Plane` / `Geom_CylindricalSurface` faces); otherwise the
per-triangle **Verbatim** path runs. It repairs dirty meshes (open edges / flipped
facets / non-manifold junctions) via sewing — those components are never refit
— merges coplanar triangles into single planar faces, fits tolerances, writes a
STEP file, and optionally re-reads it to self-verify. It is built on
OpenCASCADE (OCCT). Output is always in **millimetres**. Without `--smooth`
(the 1.x default) curved surfaces remain **faceted** at the mesh resolution,
byte-identical to 1.0.0. With **TrueForm** (`--engine trueform` / `--smooth`), v1 recovers planes, right circular
cylinders (holes/bosses, N ≥ 6), and plane–plane fillet strips (1–3 rows);
cones, spheres, and tori stay faceted.

Use it when your software consumes STEP/B-Rep but a user hands you an STL.

---

## 2. Pick an integration mode

| Your host software is… | Use | Section |
|---|---|---|
| C++ | Link the library, call `stl2step::convert()` | [4](#4-c-in-process-integration) |
| Any other language (Python, Node, Go, Rust, C#, shell) | Spawn the CLI, parse the `RESULT` or `MESH_RESULT` line | [5](#5-subprocess-integration-any-language) |
| A build that already uses CMake | `FetchContent` / `find_package` + link | [4](#4-c-in-process-integration) |

Both modes run the identical engine and expose the identical result data. Prefer
in-process for C++ hosts (no process spawn, structured `Result`); prefer
subprocess everywhere else (zero binding code, process isolation from the CAD
kernel).

---

## 3. The data contract (read this before coding)

### Options (inputs)

| Field / CLI flag | Default | Notes |
|---|---|---|
| `input` / positional | — | **Required.** Path to the `.stl`. |
| `output` / `-o` | `<input>.step` | Output path. |
| `schema` / `--schema` | `AP214` | `AP203` \| `AP214` \| `AP242`. AP214 is the safe default. |
| `inchInput` / `--units in` | false | STL is unitless; set this if it was modelled in inches (×25.4). |
| `scale` / `--scale` | 1.0 | Extra uniform scale. |
| `weldTol` / `--weld` | 0 | 0 = weld exact duplicates only. Raise for per-facet (unshared) vertices. |
| `sewTol` / `--sew-tol` | auto | Repair-pass tolerance; auto-derived from the bounding box. |
| `unifyAngleDeg` / `--unify-angle` | 0.001 | Coplanar-merge angle. Leave it. |
| `unify` / `--no-unify` | true | Merge coplanar facets. Keep ON — this is the main value. |
| `makeSolids` / `--no-solid` | true | Convert closed shells to solids. Keep ON. |
| `verify` / `--no-verify` | true | Re-read + volume-check the output. **Turn OFF if you import the STEP right after** (your import is the check). |
| `forceSew` / `--force-sew` | false | Diagnostics only; slower. |
| `threads` / `--threads` | auto | 0 = all cores. |
| `productName` | output stem | STEP product name. |
| `smooth` / `--engine trueform` / `--smooth` / `--refit` (`--engine verbatim` / `--no-smooth`) | false | **TrueForm** analytic recovery (premium). Default **Verbatim** for the whole 1.x line. Off-path STEP + RESULT are byte-identical to 1.0.0 (gate G0.1, 22/22). |
| `smoothTolMM` / `--smooth-tol <mm>` | 0 (auto) | Surface-fit tolerance in mm. 0 = auto-derived from bbox / weld / sew. Only meaningful when `smooth` is true. |
| `smoothAngleDeg` / `--smooth-angle <deg>` | 2.0 | Near-flat normal gate for segmentation (degrees). Unrelated to `unifyAngleDeg`. |
| `smoothFillets` / `--no-smooth-fillets` | true | Recover plane–plane fillet strips as cylinders. `--no-smooth-fillets` sets this false. Ignored when `smooth` is false. |

### Result (outputs)

Success flag is `ok`. On failure: `ok=false` + `error`.

| Field | Default | Notes |
|---|---|---|
| `ok` | false | True when a STEP file was written (possibly with warnings). |
| `exitCode` | 1 | 0 clean, 2 written-with-warnings, 1 failed. CLI-equivalent. |
| `input`, `output` | "" | Resolved paths. |
| `error` | "" | Human-readable reason when `ok == false`. |
| `triangles` | 0 | Triangles read from the STL. |
| `vertices` | 0 | Vertices after welding. |
| `components` | 0 | Manifold bodies the mesh split into. |
| `solids` | 0 | Solids written. |
| `openShells` | 0 | Components that could not close (written as open shells). |
| `facesBeforeUnify` | 0 | Face count before coplanar merge. |
| `facesAfterUnify` | 0 | Face count after coplanar merge. |
| `meshVolumeMM3` | 0 | Source mesh volume (mm³). |
| `stepVolumeMM3` | 0 | Volume re-read from the STEP (mm³; 0 if verify off). |
| `volumeDeltaPct` | -1 | \|step − brep\| / brep × 100 (−1 if not measured). |
| `watertight` | true | Every component closed with consistent winding. |
| `seconds` | 0 | Wall-clock time of the whole conversion. |
| `warnings` | [] | Every Warning emitted, in order. |

`smooth*` keys — **C++ `Result` members are always present and default to
zero.** The RESULT JSON emits these keys **only when `Options::smooth` was
true**, appended after `warnings`. They are omitted (never emitted as zeros)
when the feature is off, so an off-path RESULT string stays
character-identical to 1.0.0.

| Field | Default | Units / notes |
|---|---|---|
| `smoothPlanes` | 0 | Planar regions recovered as `Geom_Plane`. |
| `smoothCylinders` | 0 | Cylindrical regions recovered. |
| `smoothFillets` | 0 | Fillet strips recovered as cylinders. |
| `smoothDistinctRadii` | 0 | Distinct cylinder radii accepted. |
| `smoothRejected` | 0 | Candidate regions rejected by gates. |
| `smoothFacetFaces` | 0 | Faceted faces left after the smooth pass. |
| `facesAfterSmooth` | 0 | Total face count after the smooth pass. |
| `smoothSkippedComponents` | 0 | Dirty components not refit. |
| `smoothMaxDevMM` | 0 | Max vertex deviation from fit (mm). |
| `smoothMaxEdgeTolMM` | 0 | Max edge tolerance written (mm). |
| `smoothVolPredictedMM3` | 0 | Predicted volume from analytic fits (mm³). |

### Exit codes (CLI) / `Result::exitCode`

- `0` — clean success.
- `2` — **STEP was written, but with warnings** (e.g. an open shell, a volume
  mismatch). The file exists and is usually usable; surface the warnings.
- `1` — **failure, no output written**. `error` explains why.

> Treat `ok == true` (equivalently exit `0` or `2`) as "a file was produced."
> Only exit `1` / `ok == false` means "no file."

---

## 4. C++ in-process integration

The whole API is one header. `convert()` **never throws**.

```cpp
#include <stl2step/stl2step.hpp>

stl2step::Options opt;
opt.input  = inputStlPath;
opt.output = outputStepPath;   // or leave empty -> "<input>.step"
opt.verify = false;            // set false when you import the STEP yourself next
opt.smooth = true;             // opt-in analytic recovery; default false for 1.x

// Optional progress sink (safe to omit / pass nullptr).
auto log = [](stl2step::Severity sev, const std::string& msg) {
    if (sev != stl2step::Severity::Info) myLogger.warn(msg);
};

stl2step::Result r = stl2step::convert(opt, log);
if (!r.ok) {
    reportError(r.error);                 // no file was written
} else {
    if (!r.warnings.empty() || r.openShells > 0 || !r.watertight)
        reportWarning("mesh converted but is not a clean closed solid");
    importStep(r.output);                 // hand r.output to your kernel
}
```

### Wire it into CMake

```cmake
include(FetchContent)
FetchContent_Declare(stl2step
    GIT_REPOSITORY https://github.com/<owner>/stl2step.git
    GIT_TAG        v1.0.0)
FetchContent_MakeAvailable(stl2step)
target_link_libraries(your_target PRIVATE stl2step::core)
```

`stl2step::core` carries its own include path and OCCT link requirements
transitively — you do not repeat them. (Your project still needs OCCT installed;
the same `find_package(OpenCASCADE)` the engine uses will locate it.)

If you only want the geometry engine and not the CLI/tests in your build, set
`-DSTL2STEP_BUILD_CLI=OFF -DSTL2STEP_BUILD_TESTS=OFF` before
`MakeAvailable`/`add_subdirectory`.

### Threading / reentrancy

`convert()` is self-contained: each call owns its state, so independent
conversions may run on separate threads concurrently. The engine already
parallelises internally across cores, so for throughput prefer converting **one
file at a time** on a busy machine and only fan out on idle ones. The
`LogCallback` may be invoked from worker threads, but the engine serialises those
calls for you.

---

## 5. Subprocess integration (any language)

Spawn the CLI and read the **last stdout line**:

- **Convert mode** (default): `RESULT {json}` — exit `0` clean, `2` written-with-warnings, `1` failed.
- **Mesh mode** (`--mesh`): `MESH_RESULT {json}` — exit `0` ok, `1` failed (no exit 2).

Use `--quiet` so the contract line is the *only* stdout line. Errors/warnings
go to stderr; the last line prefix and exit code are the contract.

### Rules

- **Parse only the contract line** (`RESULT ` or `MESH_RESULT `). Everything before
  it is human progress text — do not scrape it. With `--quiet`, stdout is exactly
  one line.
- **Strip the prefix** (`RESULT ` = 7 chars, `MESH_RESULT ` = 13 chars) before
  JSON-parsing, or match the first `{` onward and branch on the prefix.
- **Trust the exit code**: convert `0`/`2` → a STEP exists at `output`; mesh `0` →
  an STL exists at `output`; `1` → no output file.
- **Set `--no-verify`** when you will load the STEP immediately after — your load
  is the verification and it is meaningfully faster on large meshes.
- **Set units.** STL is unitless. If unsure, mm is assumed.
- **`--smooth` is opt-in** (default off). `smooth*` keys appear on the RESULT
  object only when you pass it; use `result.get("smoothPlanes")` / `'smoothPlanes' in r`.
- **Mesh mode rejects convert-only flags** (`--engine`, `--schema`, `--units`, …)
  with exit `1` and **no** stdout contract line — do not reuse convert arg builders.
- **`MESH_RESULT.edges`** is always the drawable B-Rep edge count; `edgesFile` is
  present only when `--edges <path>` was passed.

### Python

```python
import json, subprocess

def stl_to_step(stl_path, step_path, *, inches=False, verify=False, smooth=False):
    cmd = ["stl2step", stl_path, "-o", step_path, "--quiet"]
    if inches:      cmd += ["--units", "in"]
    if not verify:  cmd += ["--no-verify"]
    if smooth:      cmd += ["--smooth"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    line = p.stdout.strip().splitlines()[-1]          # the RESULT line
    result = json.loads(line[len("RESULT "):])
    if not result["ok"]:                              # exit 1
        raise RuntimeError(f"stl2step failed: {result.get('error')}\n{p.stderr}")
    if result["warnings"]:                            # exit 2, file still written
        for w in result["warnings"]:
            log.warning("stl2step: %s", w)
    return step_path, result
```

### Node.js

```js
const { execFile } = require("node:child_process");
function stlToStep(stl, step, { inches = false, verify = false, smooth = false } = {}) {
  const args = [stl, "-o", step, "--quiet"];
  if (inches) args.push("--units", "in");
  if (!verify) args.push("--no-verify");
  if (smooth) args.push("--smooth");
  return new Promise((resolve, reject) => {
    execFile("stl2step", args, (err, stdout) => {
      const line = stdout.trim().split("\n").at(-1);
      const r = JSON.parse(line.replace(/^RESULT /, ""));
      if (!r.ok) return reject(new Error(r.error || "stl2step failed"));
      resolve({ step, result: r });        // r.warnings may be non-empty (exit 2)
    });
  });
}
```

### Shell

```sh
if stl2step "$in" -o "$out" --quiet --no-verify >/tmp/r 2>/dev/null; then
    :   # exit 0 or 2 -> file written; inspect "$(cat /tmp/r)" for warnings
else
    echo "conversion failed: $(sed 's/^RESULT //' /tmp/r)" >&2
fi
```

Because each invocation is independent, batch with `xargs -P` /
`parallel` for many files.

---

## 6. Behavioral contracts an integrator MUST respect

1. **Units.** Output is always mm. STL has none — pass `--units in` / `inchInput`
   or a `scale`, or your solid will be the wrong size. This is the #1 mistake.
2. **Two named modes — Verbatim (default) and TrueForm (`Options::smooth`).**
   **Verbatim** keeps curved surfaces faceted at the mesh resolution — STEP +
   RESULT are byte-identical to 1.0.0. **TrueForm** recovers **planes**, **right
   circular cylinders** (holes/bosses, N ≥ 6), and **plane–plane fillet strips**
   (1–3 rows) as `Geom_Plane` / `Geom_CylindricalSurface` faces with editable
   radii. Components that cannot close analytically fall back to Verbatim per
   component (never corrupt the file). On a real CAD export, TrueForm may build a
   fraction of recognised cylinders (e.g. 127/583 on Body11); see
   `tests/diag/body11/KNOWN-GAP.md`. Cones, spheres, and tori are reported and
   left faceted; there is no freeform/NURBS reconstruction. Refit is skipped on
   any component that needs the sewing repair path. A regular N≥6 prism (e.g. a hex socket) **is** recovered as a
   cylinder, and a symmetric 45° chamfer **is** recovered as a fillet — both
   intended, not bugs. An asymmetric chamfer (`sL/sR ≥ 1.3`) is rejected. Flat
   faces *are* recovered exactly whether or not `--smooth` is on.
3. **Non-watertight input still produces a file.** Open/dirty meshes yield exit
   `2` with `watertight=false` and/or `openShells>0`, and warnings. Decide per
   your use case whether an open shell is acceptable; the file is written either
   way. Check `openShells == 0 && watertight` if you require a true solid.
4. **Multiple bodies → multiple solids.** `components` and `solids` can be > 1
   (e.g. a mesh with several disjoint parts). Handle N solids in the STEP.
5. **`verify=false` is a performance lever, not a correctness risk** when you
   re-read the STEP yourself — your import re-parses the same file.
6. **Determinism.** Same input + options → same output. Component order is
   size-then-id sorted; thread count does not change geometry.
7. **Large meshes are memory-heavy.** >500k triangles emits a warning and a big
   STEP. If you control upstream tessellation, coarser meshes convert faster and
   smaller. There is no built-in decimation.
8. **The `RESULT` / `MESH_RESULT` line is the only stable stdout contract.** Progress
   text is for humans and may change; never parse it. Convert mode ends with
   `RESULT {json}`; mesh mode (`--mesh`) ends with `MESH_RESULT {json}`.

---

## 7. Recommended presets by host scenario

- **CAM / CAD importer that loads the STEP next:** `unify=true` (default),
 `verify=false`, **Verbatim** mode (`smooth=false`), units set. Fastest path; your
 import validates. Verbatim is the default because it optimises for round-trip
 speed into a CAM kernel that is about to re-tessellate anyway.
- **Batch/library conversion where files are archived, not immediately loaded:**
  keep `verify=true` so `volumeDeltaPct` is populated as a quality signal; run
  files concurrently.
- **Diagnosing a bad mesh:** `verify=true`, inspect `warnings`, `watertight`,
  `openShells`, and `facesBeforeUnify` vs `facesAfterUnify` (little/no reduction
  can indicate a noisy or non-planar mesh).

---

## 8. Failure handling checklist

- Exit `1` / `ok=false` → show `error` (and stderr in subprocess mode). No file.
  Common causes: unreadable/empty STL, no usable geometry, write path not
  writable.
- Exit `2` → file written; iterate `warnings` and decide if acceptable
  (open shell? volume mismatch?). `warn()` appends to `Result::warnings` and
  maps to exit 2. Example from the `--smooth` path: `smooth: analytic rebuild
  reverted on one component -- kept faceted` (R2 revert). Dirty-skip
  (`N component(s) skipped (dirty mesh, repaired path)`) and NYI rejects
  (cones, spheres, tori left faceted) are `note()`, **not** warnings — a dirty
  input converted with `--smooth` keeps its exit code versus the same input
  without it.
- Timeout guard (subprocess): conversion time scales with triangle count; set a
  generous timeout for very large meshes rather than a tight one.

---

## 9. Versioning

`stl2step::version()` / `stl2step --version` returns the semantic version
(`1.1.0-pre`). The `Options` field set, the `Result`/`RESULT` field set, and the
exit codes are the compatibility surface — additive changes only within a major
version. Pin a tag (`GIT_TAG v1.0.0`) when embedding.

RESULT addendum: when `Options::smooth` is true, `toJson()` / the CLI `RESULT`
line appends the `smooth*` keys after `warnings`. Those keys are present
**only when `smooth == true`**; they are omitted (never zero-valued) on the
off path so the RESULT string stays character-identical to 1.0.0. The C++
`Result` members are always present and default to zero (the dual contract).

---

## 10. CI discipline for agents

Private-repo Actions minutes are metered (macOS 10×, Windows 2×). After the
2026-08-30 `posix_compat.hpp` miss (`git add -u` never stages **new** files),
this is binding:

- **Never push build-affecting changes** (`src/`, `include/`, `tests/`,
  CMake, `cmake/`) without three-platform green markers:
  `.ci-local/<HEAD>.macos.green` + `.linux.green` + `.windows.green`.
  Produce them with `scripts/ci-local-gate.sh --fresh-clone`,
  `scripts/ci-linux-preflight.sh`, and `scripts/ci-windows-preflight.sh`.
  The pre-push hook (`scripts/install-git-hooks.sh`) refuses otherwise.
- **Batch fix iterations.** One push after local green — not a push per
  compile error. Each intermediate push bills a full matrix and emails
  the owner.
- **New files need an explicit `git add`.** `git add -u` is a trap: it
  stages tracked-file edits only. `git status --porcelain | grep '^??'`
  and disposition every untracked path. Working-tree / tar-synced tests
  cannot see a file that is not in git; the gate's `file://` fresh-clone
  can.
