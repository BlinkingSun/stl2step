# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **N-gon wall cylinders** (`Origin::NgonWall`): closed unclaimed A2 walls (N≥6, generators parallel) are claimed as right circular cylinders after law bands and before B1. Axis from successive plane intersections or the cap normal; radius is the circumradius, accepted when max |ρ−R| ≤ max(sewTol, chordSagitta(R,N)). Turning-axis blends and cones are skipped.
- **Planar CIRCLE loops**: LS-fit CIRCLE edges on already-analytic plane loops (N≥6) where the intersection of the two shipped surfaces yields no curve at all. The same circle handle is reused when a neighbor is already a cylinder or cone; stadium leftovers stay mixed; analytic|analytic is never forced back to polyline. A fitted circle never replaces a curve derived from the shipped surfaces — two distinct planes meet along their line and nowhere else, so an arc put there would lie on neither face; those chains keep the exact line, and the refusals are counted (`DIAG_130_PLCKEEPLIN`).
- **Chamfer cones** (`Origin::ChamferCone`): after A+B1 and before C1, a ring adjacent to one closed-360 cylinder and one cap, sharing an axis, with two radii and constant α (`|n·axis|=sin α`), ships as `Geom_ConicalSurface`. Other cones, tori, and spheres stay NYI.
- Gate fixture `linkage_bores_chamfer` (`linkage_bores_chamfer.stl`, `.expected.json`, `gates/linkage_chamfer_gate.py`) for the Mastercam linkage plate: N-gon bores plus the hole-mouth 45° chamfer.
- **Cone volume prediction**: a `ChamferCone` region predicts its facet-to-surface volume in closed form — the circular-segment gap between each inscribed facet and the frustum, summed as an exact per-triangle density over the triangles the region claims — so a built cone is no longer read as zero-area by the cascade budget and the D4 guard. At equal rim radii the formula collapses to the cylinder's own. Certified against a numeric polygon/circle oracle on the plate's frustum and on the synthetic battery frustums.
- **Planar noise floor measurement** (`DIAG_NOISEFLOOR`, `DIAG_NOISEFLOOR_ACCEPTED`, print-only): the max, p90 and median vertex-to-plane residual over the plane regions a mesh's own segmentation accepts. Measured because a recorded decision proposed that maximum as an admission certificate; on both measured parts it is the plane-growth tolerance showing through on a handful of curved facets (max 1.9e-2 mm against a median of 1.7e-15 mm on the linkage plate) rather than a noise level, so it is reported and nothing consumes it.

### Changed
- **Facet fill planes**: a mesh facet is filled on the plane of its own wire. A shared vertex on the boundary of an analytic region is snapped onto that region's curve at construction, and a plane taken from the unsnapped mesh array then does not contain the edges the face is built from; the face reads `InvalidCurveOnSurface` in its own context once `SameParameter` recomputes the shared edge's tolerance from the stored representations, and becomes unorientable. The refill fires only when a wire vertex has left the mesh plane by more than the tolerance the face will carry, and the mesh normal still decides the facet's side.
- **Cascade containment**: a face that is not an analytic region (an island facet or a fill triangle) explodes its own neighbourhood — the regions sharing its edges — instead of contributing no culprit at all and escalating the component to the blanket rung. The 50 % escalation rule itself is unchanged; only what counts as a culprit is.
- **G2 sagitta in `archChainBand`**: outside coarse, G2 accepts a cylinder commit only when the max vertex residual is within the exact deviation a correct fit of a **counted** N-gon can produce (`N = ev.d2.nSides`, and only when `N >= 6`); the gate is never loosened for a region whose N was not counted (no `tris.size()` substitute, no literal floor of 6). Otherwise `g2Tol = epsCylAccept(R)`. `coarseFusionBand` is not widened (D-130-1).

## [1.2.0] - 2026-09-03

### Changed
- **Plate-wire construction** (TrueForm): each planar region's boundary loop is now ordered by an exact vertex-connectivity walk over its own chain edges and oriented by continuity at construction; the loop's sense (outer counter-clockwise about the outward normal, inner clockwise) is taken from its signed area before orientation, so shipped edge orientations follow the part's geometry rather than an incidental construction order. On the bundled Handle pickup fixture all four large flats now build valid on the first attempt, including the flange, which ships as its two CAD faces.
- **Edge tolerances at construction**: the tolerance written at the pcurve bind site is the exact closed-form maximum of the 3D-to-pcurve deviation over all of an edge's pcurves plus the strict-inequality margin, at one site for planes and cylinders; inherited tolerances equal to a deviation are replaced, fat inherited tolerances are never lowered, and a deviation above the mesh tolerance cap is counted as a construction defect rather than absorbed.
- **Determinism**: plate-wire construction is a pure function of each loop's own chain set (a process-wide static in the collapsed-chain rebuild made single-threaded runs differ from multi-threaded ones); every corpus part now produces identical results and STEP data under any thread count.

### Added
- Gates: the plane census maps each shipped face to one ground-truth face (co-equation faces are no longer double counted); the partial-recovery phantom ratchet is re-derived from the current shell with the one new phantom attributed by region.
- Diagnostics (all off unless `STL2STEP_P2_DIAG` is set): first-keep in-context census, tolerance/supremum prints, wire-order and orientation prints, repair-rung fire counters.

### Fixed
- Builds against OpenCASCADE 8.0 as well as 7.9.
- The four p2buildtest fixtures that had been failing (counterbore, full 360 hole, r8 360 explode caps, s09 mixed) pass: their plate wires close under the connectivity rule.
- Diagnostic prints no longer construct on live shapes (they were altering the build when enabled).

### Known limitations (registered)
- 111 of 638 Handle pickup seam chains ship as mesh polylines rather than analytic lines/arcs (surfaces exact, some edges not).
- 27 small sliver plates reach validity only through the tolerance repair pass (their fitted plane does not contain their own mesh boundary).
- Nine single-closed-edge inner loops decide a traversal reversal that has no carrier on a one-edge loop; no output consequence (reported by an internal assertion; fix scoped).

## [1.1.0] - 2026-08-30

### Added
- **Tessellation-law arc recognition** (TrueForm): exported meshes are
  recognized by the deterministic law their tessellator followed; radii are
  recovered by the closed-form inverse `R = w/(2·sin(θ/2))` with
  self-calibrated, never hard-coded, law parameters. Law-driven segmentation
  claims arc bands before plane absorption, splits mixed-radius chimeras, and
  is scored against labeled triangle↔surface truth in the test suite (G-LAW).
- **Prismatic (2.5D) reconstruction** (TrueForm): prismaticity detection,
  cap-level slicing, 2D line/arc profile fitting, and extrude-and-unite
  rebuild — curved walls become analytic cylinders by construction. On the
  bundled fully-prismatic fixture, all 15 cylinders ship analytic with radii
  within 0.3% of the source CAD at 0.000000% volume deviation.
- **`--dxf <dir>`**: emit one DXF per level of the recognized 2D profiles
  (prismatic parts; inert when unset, output byte-identical).
- Arc-aware volume authority for reconstructed solids (arc-vs-chord defect
  budgets computed from the law); phantom-volume classes remain rejected.
- New gate families: G-LAW supervised recognition, G-PRISM, prism routing,
  census/GT-radii ratchets, per-body free-edge red-lines, cascade waiver
  tripwire, stress-sweep harness.
- `docs/METHOD.md`: an educational description of the two-stage engine and
  the recognition method.

### Changed
- J6 closure improvements on real CAD exports (free edges reduced up to 4×
  on benchmark bodies) with volume strictly preserved.
- Verbatim (Stage 1) remains byte-identical (G0.1) and remains the default;
  the engine is explicitly two-part — stop at Stage 1, or run TrueForm.

### Fixed
- CI (Linux): install `libtbb-dev` alongside the OCCT packages — Ubuntu's
  OpenCASCADE CMake config references TBB by absolute path, so linking failed
  with "No rule to make target '/usr/lib/x86_64-linux-gnu/libtbb.so'".
- CI: update `actions/checkout` v4 -> v5 to clear Node 20 deprecation warnings.

## [1.1.0-unreleased]

Optional analytic-surface recovery behind **TrueForm** mode (`--engine trueform`,
`--smooth`). No tag, no release.

### Added
- **Verbatim / TrueForm** branding: two named conversion modes. **Verbatim**
  (default) is byte-faithful tessellation, identical to 1.0.0. **TrueForm** is
  premium analytic reconstruction (today's `--smooth` path).
- CLI `--engine verbatim|trueform` with `--smooth` / `--refit` / `--no-smooth`
  aliases mapping to TrueForm / Verbatim respectively.
- **TrueForm** opt-in recovery of **planes**, **right circular cylinders**
  (holes/bosses, N ≥ 6), and **plane–plane fillet strips** (1–3 rows) as
  `Geom_Plane` / `Geom_CylindricalSurface` faces with editable radii. Library
  field: `Options::smooth`. Default **Verbatim** for the whole 1.x line.
- Six CLI flags: `--engine`, `--smooth`, `--refit` (TrueForm aliases),
  `--no-smooth` (Verbatim alias), `--smooth-tol <mm>`, `--smooth-angle <deg>`,
  `--no-smooth-fillets`. Matching `Options` fields: `smooth`, `smoothTolMM`
  (0 = auto, mm), `smoothAngleDeg` (default 2.0 deg), `smoothFillets` (default
  true; `--no-smooth-fillets` clears it).
- Eleven `Result` fields: `smoothPlanes`, `smoothCylinders`, `smoothFillets`,
  `smoothDistinctRadii`, `smoothRejected`, `smoothFacetFaces`,
  `facesAfterSmooth`, `smoothSkippedComponents`, `smoothMaxDevMM`,
  `smoothMaxEdgeTolMM`, `smoothVolPredictedMM3`. C++ members are always
  present and default to zero; RESULT `smooth*` keys are spliced only when
  `smooth == true`, appended after `warnings` (omitted, never zero-valued,
  on the off path).
- Five `smoothBuilt*` fields (`smoothBuiltPlanes`, `smoothBuiltCylinders`,
  `smoothBuiltFillets`, `smoothBuiltComponents`, `smoothRevertedComponents`)
  counting post-unify file-truth census on components that delivered cylinders.
- Smooth-on acceptance gates in `tests/gates/run_gates.py` + `smooth_on.py`
  (G0.2, G0.3, G2, G2.5, G3, G4.4, G5, R-ladder, calibration) with a **parked**
  mechanism for documented reds (`--unpark` to surface honest fails).
- Documented D6.4 limitations: a regular N≥6 prism (e.g. a hex socket) **is**
  recovered as a cylinder; a symmetric 45° chamfer **is** recovered as a
  fillet; an asymmetric chamfer (`sL/sR ≥ 1.3`) is rejected.

### Changed
- Engine version string is `1.1.0-pre` (`STL2STEP_VERSION_*` macros). CMake
  `project(... VERSION …)` is owned separately and is not part of this
  changelog entry.
- **Verbatim mode is byte-identical to 1.0.0** (STEP + RESULT), measured by
  gate G0.1 at 22/22 on every corpus fixture.
- README and AGENTS.md present Verbatim/TrueForm by name, including honest
  real-world numbers (127/583 cylinders built on Body11) and a KNOWN-GAP pointer
  to `tests/diag/body11/KNOWN-GAP.md`.

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
