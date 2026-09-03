# stl2step

**A universal engine that converts triangle meshes (STL) into parametric B-Rep
solids (STEP).** Point it at an `.stl`, get back a clean `.step` solid that CAD
and CAM kernels can consume as real boundary geometry — not triangle soup.

It ships as both an **embeddable C++ library** (one header, one call) and a
**standalone command-line tool**, and builds from the same source on **Linux,
macOS, and Windows**.

```
part.stl  ──►  weld ─► split into solids ─► build B-Rep (parallel) ─► repair
          ─►  merge coplanar facets ─► fit tolerances ─► verify  ──►  part.step
```

- **Watertight → solid.** Closed, consistently wound meshes become oriented STEP
  solids. Flat regions collapse from thousands of triangles into single planar
  faces.
- **Robust to dirty meshes.** Open edges, flipped facets, and non-manifold
  junctions are split into clean bodies and repaired (sewing + shape-fix); what
  still can't close is written as an open shell and flagged, never dropped
  silently.
- **Multiple bodies.** Meshes containing several disjoint solids convert to
  several solids in one file.
- **Fast.** Per-body B-Rep construction is lock-free parallel across all cores,
  and the validity/volume checks overlap the STEP write. Typical small parts
  convert in a fraction of a second.
- **Machine-friendly.** Every run ends with a single-line `RESULT {json}` and a
  meaningful exit code, so it drops cleanly into scripts, pipelines, and agents.

> Integrating this into your own software (including via an AI coding agent)?
> Read **[AGENTS.md](AGENTS.md)** — a focused, copy-paste integration guide.

---

## What it is and isn't

stl2step exposes two named conversion modes:

| Mode | CLI | What you get |
|---|---|---|
| **Verbatim** (default) | `--engine verbatim` / `--no-smooth` | Byte-faithful tessellation. Curved surfaces stay faceted at the STL resolution. STEP + RESULT are **byte-identical to 1.0.0** (gate G0.1, 22/22). |
| **TrueForm** (premium) | `--engine trueform` / `--smooth` / `--refit` | Analytic recovery: tessellation-**law-based** arc recognition (exact radii via the chord/angle inverse), **planes**, **right circular cylinders**, fillet strips — and for **prismatic (2.5D) parts**, full profile-based reconstruction: 2D line/arc profiles per level (optionally emitted as DXF via `--dxf <dir>`), extruded back into a solid whose curved walls are analytic **by construction**. |

![STL mesh to faceted solid to analytic solid](docs/assets/conversion-stages.png)

*The two stages: input mesh → Stage 1 faceted solid (Verbatim) → Stage 2
analytic solid (TrueForm). Stop after Stage 1, or run the full conversion.
The method is described in [`docs/METHOD.md`](docs/METHOD.md).*

**Verbatim is always safe.** TrueForm runs per component; any component that cannot close analytically falls back to Verbatim for that body only — the file is still written, never corrupted.

On a fully prismatic real CAD export (the `handle-lock` fixture: 908 triangles, 28 true faces), TrueForm reconstructs **all 15 cylinders as analytic faces with every radius matching the source CAD within 0.3%**, at 0.000000% volume deviation, watertight and valid. On non-prismatic parts the general path applies; the remaining recognition gaps are documented in [`tests/diag/body11/KNOWN-GAP.md`](tests/diag/body11/KNOWN-GAP.md) — not silent regressions.

Flat faces are recovered exactly in both modes — coplanar triangles merge into single planar faces. This is the main win for downstream CAM even without TrueForm.

TrueForm does **not** recover cones, spheres, or tori (reported, left faceted), does not reconstruct freeform/NURBS, and **skips refit on any component that needs the sewing repair path**.

**Documented limitations, not bugs.** A regular N≥6 prism (e.g. a hex socket) **is** recovered as a cylinder. A symmetric 45° chamfer **is** recovered as a fillet. An asymmetric chamfer (`sL/sR ≥ 1.3`) is rejected.

Output is always written in **millimetres**. STL is unitless, so tell the engine the input units (see `--units` / `Options::inchInput` / `Options::scale`).

## SolidOut: the desktop app (macOS)

**SolidOut** is a desktop app over this engine for people who want the conversion without a terminal: open an STL, look at the mesh, convert with Verbatim or TrueForm, read the result (planes, cylinders, fillets, volume delta, watertight check) and export the STEP. The engine and OpenCASCADE are bundled inside the app, so nothing else needs to be installed and nothing runs in the background: the app launches the bundled engine for each conversion.

Download: [SolidOut 0.1.0 for macOS (Apple silicon)](https://github.com/BlinkingSun/stl2step/releases/download/v1.2.0/SolidOut-0.1.0-macOS-AppleSilicon.dmg) from the [1.2.0 release](https://github.com/BlinkingSun/stl2step/releases/tag/v1.2.0). It is ad-hoc signed, not notarized: on first launch right-click the app and choose Open.

![SolidOut with the Handle pickup mesh loaded (3,338 triangles)](https://github.com/BlinkingSun/stl2step/releases/download/v1.2.0/solidout-import-mesh.png)

![SolidOut after TrueForm conversion: 918 faces, watertight, 0.000% volume delta](https://github.com/BlinkingSun/stl2step/releases/download/v1.2.0/solidout-step-result.png)

*The same part before and after: the STL tessellation on the left screen, the analytic STEP solid (70 planes, 121 cylinders) on the right, converted in under 20 seconds on an M-series Mac.*

---

---

## Dependencies

- A **C++17** compiler (Clang, GCC, or MSVC).
- **[OpenCASCADE Technology (OCCT)]** 7.x — the geometry kernel that does the STL
  reading, B-Rep construction, healing, and STEP writing. (Tested against 7.9;
  7.6+ is expected to work.)
- **CMake** 3.16+.

Install OCCT per platform:

| Platform | Command |
|---|---|
| macOS (Homebrew) | `brew install opencascade` |
| Debian / Ubuntu  | `sudo apt install libocct-foundation-dev libocct-modeling-data-dev libocct-modeling-algorithms-dev libocct-data-exchange-dev` |
| Windows (vcpkg)  | `vcpkg install opencascade` |
| Fedora           | `sudo dnf install opencascade-devel` |
| Arch             | `sudo pacman -S opencascade` |

The build locates OCCT automatically via its own CMake package; on distributions
that don't ship that config, the bundled `cmake/FindOpenCASCADE.cmake` discovers
headers and toolkits manually. Point it at a custom install with
`-DOpenCASCADE_ROOT=/path/to/occt` or `-DCMAKE_PREFIX_PATH=...`.

[OpenCASCADE Technology (OCCT)]: https://dev.opencascade.org/

---

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
# optional:
ctest --test-dir build --output-on-failure     # runs the smoke test
```

Or with the bundled presets:

```sh
cmake --preset dev      # Release + examples + tests
cmake --build --preset dev
ctest --preset dev
```

The CLI lands at `build/stl2step` (`build/Release/stl2step.exe` on Windows).

### Useful CMake options

| Option | Default | Meaning |
|---|---|---|
| `STL2STEP_BUILD_CLI`      | `ON`  | Build the `stl2step` command-line tool |
| `STL2STEP_BUILD_EXAMPLES` | `OFF` | Build the C++ examples in `examples/` |
| `STL2STEP_BUILD_TESTS`    | top-level | Register the CTest smoke test |
| `STL2STEP_INSTALL`        | top-level | Generate install + `find_package` rules |

---

## Command-line usage

```sh
stl2step part.stl                         # writes part.step next to the input (Verbatim)
stl2step part.stl -o out/part.step        # explicit output
stl2step part.stl --engine trueform       # TrueForm analytic recovery
stl2step part.stl --units in              # STL modelled in inches -> scaled to mm
stl2step part.stl --no-verify             # fastest: skip the re-read self-check
stl2step part.stl --smooth                # alias for --engine trueform
stl2step part.stl --schema AP242 --threads 4
stl2step --mesh part.step -o part.stl --edges part.edges   # mesh mode (STEP → STL)
```

Progress prints to stdout (silence with `--quiet`); warnings and errors go to
stderr. In **convert** mode the last stdout line is `RESULT {json}`; in **mesh**
mode it is `MESH_RESULT {json}`.

```
$ stl2step part.stl
stl2step: part.stl -> part.step
  read      12 triangles, 8 vertices ...
  ...
  verify    re-read 1 solid, 6 faces, volume delta 0.000000%  [0.00s]
RESULT {"ok":true,"input":"part.stl","output":"part.step","triangles":12,"solids":1, ... }
```

Run `stl2step --help` for the full option list.

### Mesh mode (`--mesh`)

Tessellate an existing **STEP** solid to a **binary STL**, optionally emitting
**drawable B-Rep edges** (not facet boundaries) for overlay in a viewer.

```sh
stl2step --mesh part.step -o part.stl --edges part.edges --quiet
```

| Flag | Meaning |
|---|---|
| `--mesh <file>` | Input STEP (`.step` / `.stp`). Selects mesh mode. |
| `-o <file>` | Output binary STL (default: input stem + `.mesh.stl`; never clobbers an existing default — pass `-o`) |
| `--edges <file>` | Optional Format A edge buffer (see below) |
| `--quiet` | Suppress progress (`MESH_RESULT` + errors only) |
| `--threads <n>` | Tessellation worker threads (default: all cores) |

**Rejected in mesh mode** (stderr + exit 1, **no** `MESH_RESULT` / `RESULT` line):
`--engine`, `--schema`, `--units`, and all other convert-only flags.

#### `MESH_RESULT` contract

Last stdout line: `MESH_RESULT {json}`. Exit **0** on success, **1** on failure
(no exit 2 tier).

| key | type | meaning |
|---|---|---|
| `ok` | bool | STL was written |
| `faces` | int | B-Rep face count |
| `edges` | int | drawable B-Rep edge count (always present; not zero when `--edges` is omitted) |
| `triangles` | int | triangles in the STL |
| `seconds` | double | wall-clock time |
| `input`, `output` | string | resolved paths |
| `edgesFile` | string | present only when `--edges` was passed |
| `error` | string | present only when `ok` is false |

#### Format A edge buffer

Little-endian **float32** **xyz pairs**: one segment = `x0 y0 z0 x1 y1 z1` =
**24 bytes**. No header, no sentinel, no count prefix. File size must be a
multiple of 24; **0 bytes** means no drawable edges.

### Options (convert mode)

| Flag | Meaning |
|---|---|
| `-o <file>` | Output path (default: input with `.step`) |
| `--schema AP203\|AP214\|AP242` | STEP application protocol (default `AP214`) |
| `--units mm\|in` | Input units; `in` scales ×25.4 to mm |
| `--scale <f>` | Extra uniform scale factor |
| `--weld <tol>` | Weld vertices within `tol` (default: exact duplicates only) |
| `--sew-tol <tol>` | Tolerance for the repair pass (default: auto from bbox) |
| `--unify-angle <deg>` | Max normal deviation treated as coplanar (default `0.001`) |
| `--no-unify` | Keep one face per triangle |
| `--no-solid` | Emit shells only |
| `--force-sew` | Route every body through the repair path |
| `--no-verify` | Skip re-reading the output (see below) |
| `--engine verbatim\|trueform` | Conversion mode (default `verbatim`). `trueform` = analytic recovery |
| `--smooth` | Alias for `--engine trueform` |
| `--refit` | Alias for `--engine trueform` |
| `--no-smooth` | Alias for `--engine verbatim` |
| `--smooth-tol <mm>` | Surface-fit tolerance in mm (default: auto) |
| `--smooth-angle <deg>` | Near-flat normal gate for segmentation (default 2.0) |
| `--no-smooth-fillets` | Skip recovery of fillet strips as cylinders |
| `--threads <n>` | Worker threads (default: all cores) |
| `--quiet` | Suppress progress (RESULT + errors only) |

> **`--no-verify`:** the default re-reads the written STEP and compares volumes as
> a self-check. If your program is about to import the STEP anyway, that import
> *is* the verification — pass `--no-verify` to roughly halve wall time on large
> files.

---

## Embedding the library

The entire API is `include/stl2step/stl2step.hpp`:

```cpp
#include <stl2step/stl2step.hpp>

stl2step::Options opt;
opt.input  = "part.stl";
opt.output = "part.step";
opt.verify = false;                 // we're about to import it ourselves

stl2step::Result r = stl2step::convert(opt, /*log=*/nullptr);
if (r.ok) {
    // r.solids, r.facesAfterUnify, r.meshVolumeMM3, r.warnings, r.toJson(), ...
}
```

`convert()` never throws — failures come back as `r.ok == false` with `r.error`.
Pass a `LogCallback` to receive progress. See `examples/convert_min.cpp`.

### Add it to your CMake project

**FetchContent** (no install needed):

```cmake
include(FetchContent)
FetchContent_Declare(stl2step
    GIT_REPOSITORY https://github.com/<you>/stl2step.git
    GIT_TAG        v1.0.0)
FetchContent_MakeAvailable(stl2step)

target_link_libraries(your_app PRIVATE stl2step::core)
```

**Installed package:**

```sh
cmake --install build --prefix /your/prefix
```
```cmake
find_package(stl2step CONFIG REQUIRED)
target_link_libraries(your_app PRIVATE stl2step::core)
```

**Vendored:** copy `include/stl2step/stl2step.hpp` + `src/stl2step.cpp` into your
tree and link OCCT yourself.

### Calling from another language

No bindings required — spawn the CLI as a subprocess and read the final `RESULT`
line (or pass `--quiet` so it's the only stdout line). Python, Node, Go, Rust,
C#, shell — all work the same way. Recipes are in **[AGENTS.md](AGENTS.md)**.

---

## The RESULT contract

The `RESULT {json}` object (also `Result::toJson()`) has a stable field set:

| Field | Type | Meaning |
|---|---|---|
| `ok` | bool | A STEP file was written (possibly with warnings) |
| `input`, `output` | string | Resolved paths |
| `error` | string | Present only when `ok` is false |
| `triangles`, `vertices`, `components` | int | Mesh + body counts |
| `solids`, `openShells` | int | Solids written; components that couldn't close |
| `facesBeforeUnify`, `facesAfterUnify` | int | Face count around the coplanar merge |
| `meshVolumeMM3`, `stepVolumeMM3` | number | Source vs. re-read volume (mm³) |
| `volumeDeltaPct` | number | Round-trip volume error (`-1` if not measured) |
| `watertight` | bool | Every component closed with consistent winding |
| `seconds` | number | Wall-clock time |
| `warnings` | string[] | Every warning emitted |

When `--smooth` is on, the following `smooth*` keys are **appended after
`warnings`**. They are present **only when `smooth == true`**; they are omitted
(never emitted as zeros) when the flag is off, so an off-path RESULT string
stays character-identical to 1.0.0. The C++ `Result` members are always present
and default to zero.

| Field | Type | Meaning |
|---|---|---|
| `smoothPlanes` | int | Planar regions recovered as `Geom_Plane` |
| `smoothCylinders` | int | Cylindrical regions recovered |
| `smoothFillets` | int | Fillet strips recovered as cylinders |
| `smoothDistinctRadii` | int | Distinct cylinder radii accepted |
| `smoothRejected` | int | Candidate regions rejected by gates |
| `smoothFacetFaces` | int | Faceted faces left after the smooth pass |
| `facesAfterSmooth` | int | Total face count after the smooth pass |
| `smoothSkippedComponents` | int | Dirty components not refit |
| `smoothMaxDevMM` | number | Max vertex deviation from fit (mm) |
| `smoothMaxEdgeTolMM` | number | Max edge tolerance written (mm) |
| `smoothVolPredictedMM3` | number | Predicted volume from analytic fits (mm³) |

**Exit codes:** `0` clean · `2` STEP written but with warnings (open shell,
volume mismatch, …) · `1` failed, no output written.

---

## Performance notes

Vertex weld, component split, and B-Rep construction are parallel across cores;
the validity check and volume integration run on a worker thread *while* the STEP
serialiser writes, so they cost effectively no extra wall time. The remaining
serial floor is OCCT's STEP writer plus the per-body coplanar merge. To go
faster on big jobs, decimate the mesh upstream, pass `--no-verify` when you'll
re-import anyway, or convert many files concurrently (each run is independent —
e.g. `ls *.stl | xargs -P8 -I{} stl2step {}`).

---

## CI

Node preflights (`scripts/ci-linux-preflight.sh` and `scripts/ci-windows-preflight.sh`) are single-tenant and refuse to start while another run holds the node lock (exit 5).

## License

MIT — see [LICENSE](LICENSE). Note that OpenCASCADE, the geometry kernel this
engine links against, is LGPL-2.1-with-exception and is a separate dependency
(not bundled). Distributors of linked binaries should review the OCCT license.
