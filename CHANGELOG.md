# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/).

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
