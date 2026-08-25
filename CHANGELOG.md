# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- CI (Linux): install `libtbb-dev` alongside the OCCT packages — Ubuntu's
  OpenCASCADE CMake config references TBB by absolute path, so linking failed
  with "No rule to make target '/usr/lib/x86_64-linux-gnu/libtbb.so'".
- CI: update `actions/checkout` v4 -> v5 to clear Node 20 deprecation warnings.

## [1.1.0-unreleased]

Optional analytic-surface recovery behind `--smooth`. No tag, no release.

### Added
- `--smooth` opt-in recovery of **planes**, **right circular cylinders**
  (holes/bosses, N ≥ 6), and **plane–plane fillet strips** (1–3 rows) as
  `Geom_Plane` / `Geom_CylindricalSurface` faces with editable radii. Library
  field: `Options::smooth`. Default **OFF** for the whole 1.x line.
- Six CLI flags: `--smooth`, `--refit` (alias), `--no-smooth`,
  `--smooth-tol <mm>`, `--smooth-angle <deg>`, `--no-smooth-fillets`.
  Matching `Options` fields: `smooth`, `smoothTolMM` (0 = auto, mm),
  `smoothAngleDeg` (default 2.0 deg), `smoothFillets` (default true;
  `--no-smooth-fillets` clears it).
- Eleven `Result` fields: `smoothPlanes`, `smoothCylinders`, `smoothFillets`,
  `smoothDistinctRadii`, `smoothRejected`, `smoothFacetFaces`,
  `facesAfterSmooth`, `smoothSkippedComponents`, `smoothMaxDevMM`,
  `smoothMaxEdgeTolMM`, `smoothVolPredictedMM3`. C++ members are always
  present and default to zero; RESULT `smooth*` keys are spliced only when
  `smooth == true`, appended after `warnings` (omitted, never zero-valued,
  on the off path).
- Documented D6.4 limitations: a regular N≥6 prism (e.g. a hex socket) **is**
  recovered as a cylinder; a symmetric 45° chamfer **is** recovered as a
  fillet; an asymmetric chamfer (`sL/sR ≥ 1.3`) is rejected.

### Changed
- Engine version string is `1.1.0-pre` (`STL2STEP_VERSION_*` macros). CMake
  `project(... VERSION …)` is owned separately and is not part of this
  changelog entry.
- **`--smooth` off is byte-identical to 1.0.0** (STEP + RESULT), measured by
  gate G0.1 at 22/22 on every corpus fixture.

## [1.0.0]

First public release as a standalone, cross-platform engine.

### Added
- Reusable engine library (`stl2step::core`) exposing a single
  `stl2step::convert(Options, LogCallback) -> Result` entry point, plus
  `version()`, `schemaName()`, and `parseSchema()`. The whole public interface
  is one header (`include/stl2step/stl2step.hpp`).
- Standalone `stl2step` command-line tool built on the library.
- `Result::toJson()` producing the stable, machine-readable `RESULT {json}`
  payload that the CLI prints as its last stdout line.
- Optional `LogCallback` progress sink (Info / Warning / Error) so host
  applications can capture or silence engine output.
- Portable CMake build with a robust `FindOpenCASCADE.cmake` (prefers OCCT's own
  config package; falls back to manual discovery, handling the 7.8 toolkit
  rename `TKSTEP`/`TKSTL` -> `TKDESTEP`/`TKDESTL`).
- Install + package-config export so downstream projects can
  `find_package(stl2step)`; also usable via `add_subdirectory` / `FetchContent`.
- `AGENTS.md` integration guide, a minimal C++ example, and a CTest smoke test.

### Verified
- Built and smoke-tested on macOS (AppleClang, OpenCASCADE 7.9.3) and Windows
  (MSVC 2022 Build Tools, vcpkg OpenCASCADE 8.0.0). Both produce identical
  conversion results (cube: 1 solid, 6 faces, 1000 mm^3, 0% volume delta). Linux
  uses the same standard CMake + system-OCCT path exercised by CI.

### Changed (from the original single-file tool)
- Refactored the monolithic `main()` into a reentrant library `Converter` plus a
  thin CLI, with no change to the conversion pipeline or numeric behavior.
- Replaced the POSIX `<sys/stat.h>` file-size call with `std::filesystem` and
  `M_PI` with an internal constant, so the engine builds cleanly on MSVC as well
  as Clang/GCC.

### Engine pipeline (unchanged)
- STL read (binary + ASCII) -> exact/tolerance vertex weld -> manifold component
  split -> lock-free parallel per-body B-Rep construction -> sewing + ShapeFix
  repair fallback for dirty meshes -> coplanar facet unify -> planar tolerance
  fit -> validity + volume check overlapped with the STEP write -> optional
  re-read verify. Exit codes: 0 clean, 2 written-with-warnings, 1 failed.
