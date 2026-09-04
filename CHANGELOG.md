# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.3.0] - 2026-09-05

### Added
- **Cones as a surface type: chamfer frustums** (`Origin::ChamferCone` → `Geom_ConicalSurface`; D-130-3). Detector C claims a ring adjacent to one accepted closed-360 cylinder and one cap, sharing the axis, with two rim radii and a constant slope, and ships it as a real cone face: the face is built on the same surface handle the bind site uses (`SurfVar::ConeBase`, from the fitted rims via `coneFromRims`), rim and seam pcurves are written explicitly (a seam needs two, at `u0` and `u0 + 2π`), the seam is an explicit generator with both vertices snapped to it, and the apex is excluded by construction because the face's v-range runs between two positive-radius rims (b50328c, 49952bc). The taper sign is carried in a signed height and a signed `SemiAngle`, never by negating the axis, so both frusta of a hole-mouth pair build in one frame (30f67dd). Detector C runs a first pass right after the law bands and detector A — before B1 — so the exact frustum claims the ring before B1 can commit ~60° arcs of it as cylinders; the original post-B1 pass is kept (978ab69; D-130-3, D-130-8). General cones (apex inside the face, oblique sections) are not built (D-130-3).
- **Certified bind classes on cones** — `circle-on-cone` and `line-on-cone` (`src/refit_cone_bind.{hpp,cpp}`, `cone_bind_unit`; 10ce7cd; D-130-3's hard gate). A rim circle binds only when its pcurve is a v-iso line and the circle is coaxial with the cone, in which case `|C3d(t) − S(pc(t))|` is constant in `t` and the endpoint value *is* the supremum (proved, not sampled); a seam binds as an affine generator whose supremum is at an endpoint. Every other cone edge returns an `unhandled-…-cone` class and ships as the mesh polyline — tier 2, counted (D-130-2).
- **Cone math** (`src/refit_cone_math.{hpp,cpp}`, `cone_math_unit`; b3114da, 7f4eb5f, 0d67eba): `circleOnConeMax` (exact for coaxial and perpendicular-plane circles, a provable upper bound otherwise, invariant under rigid placement), the exact `pointConeDist`, tier-1 `IntAna` cone cases (cone∩plane ⟂ axis, coaxial cone∩cylinder, coaxial cone∩cone; everything oblique returns `false`), and the frustum helpers. Kept in its own translation unit so the P1 include allowlist is neither touched nor widened (0d67eba).
- **Cone volume prediction** (19ec800; D-130-11(2)): a `ChamferCone` region predicts its facet-to-surface volume in closed form — the circular-segment gap between each inscribed facet and the frustum, summed as an exact per-triangle density over the triangles it claims — so a built cone is no longer read as zero-area by the cascade budget and the D4 guard (I9 green on the plate). At equal rim radii the formula collapses to the cylinder's own.
- **Mesh quantization floor** (`src/stl_quant.{hpp,cpp}`, `DIAG_QUANTFLOOR`, `quant_floor_unit`; 61fab1d; D-130-12): the resolution of the number format a mesh file was written in, at that mesh's own coordinate magnitude — for a binary STL the binary32 ulp at the largest |coordinate| in the file times √3/2; for an ASCII STL the printed-decimal grid at the same magnitude. Measured from the file; an unreadable file yields 0, which certifies nothing. It is the certificate the absorb below runs on.
- **Law-band extent — the absorb** (4af42c2, 3d25757; D-130-9, D-130-12; supersedes D-130-11(1)): a recognised law band owns every adjacent unclaimed triangle whose vertices lie on the cylinder the band has already certified, admitted only inside `max(the band's own residual, the mesh's quantization floor)` — never a part-scale budget. After each round the band's cylinder is refit by least squares over its members **with the axis direction among the unknowns**, and a round whose refit leaves the certificate is rolled back whole. On the plate this closes the R 8.5 through-bore (6 → 139 triangles) and completes both R 30 profile arcs to exactly the sidecar's 75 and 80 triangles. The planar noise floor D-130-11(1) proposed as the certificate was measured (`DIAG_NOISEFLOOR`, a7f7b99) and refused: on the plate its median is 1.7e-15 mm and its max 0.019 mm — 105× looser than the budget it was meant to tighten.
- **Virtual generators** (e2cf8f5; D-130-8): a law band's chain may be built on the lines where two adjacent facet planes meet rather than on the mesh edges the tessellation happened to draw, so a staggered strip chains at every fold the mesh resolves. The virtual path runs only where the mesh-edge chain refused. On the plate the four R 3 slot ends are recognised for the first time, the R 8.5 bore is claimed as a full 92-generator closed chain, and the cross bore reaches 214 of its 279 ground-truth triangles.
- **N-gon wall cylinders** (`Origin::NgonWall`, detector A; 49952bc, ad5c814): closed unclaimed A2 walls (N ≥ 6, generators parallel) are claimed as right circular cylinders after the law bands and before B1. The accept test is the **vertex residual** against the fitted circle, with no chord-sagitta term (ad5c814): with the sagitta term the non-concyclic hex control `S12-b` was accepted at a 7 %-of-R residual; without it, refused. Detector A can certify a through hole: a link sharing exactly one mesh edge with its neighbour is no longer read as a turning-axis blend, and a single-cap ring (a chamfered through hole) is collected (ad5c814).
- **Planar CIRCLE loops** (detector B; 49952bc, 88b14fd, 74f885c; D-130-2): LS-fit CIRCLE edges on already-analytic plane loops (N ≥ 6) where the two shipped surfaces yield no curve at all. A plane-loop CIRCLE replaces the polyline only when both sides of the chain are real regions that can share it (88b14fd), and it never replaces a curve derived from the shipped surfaces — two distinct planes meet along their line and nowhere else (74f885c, 47bdd57; `DIAG_130_PLCKEEPLIN`).
- **Edge-class census and gate** (`edgeClasses` in RESULT; `tests/gates/edge_class_gate.py` + `tests/gates/baseline/edge-class-ratchet.json`; 9249ab0, 149564f, 082194b; D-130-2, D-130-5, D-130-7(c)). Every edge two analytic faces share is classified on the shipped shell: tier 1 (`analytic`) when the curve is a line, circle or ellipse and the closed-form bind-site supremum exists on both surfaces; tier 2 (`polylineTier2`) when it is a mesh polyline — counted, never absorbed; `unhandled` otherwise. Red lines `unhandled*=0`, `dev>tol=0`, `tol>meshTolCap=0`; the tier-2 count is a never-increase ceiling. **Expected-red on 9 fixtures, recorded not relaxed** (082194b): Body11 (`overCap`), Body28 (`overCap`,`unhandled`), S03 (`overCap`), S04 (`overCap`), S09 (`overCap`), S11-b (`unhandled`), S16-R1-explode-success (`overCap`,`unhandled`), handle-pickup (`overCap`,`overTol`,`unhandled`), linkage_bores_chamfer (`overCap`,`overTol`,`unhandled`).
- **Radius truth census** (`radiusDrift` in RESULT, `DIAG_RADIUS`; 082194b; D-130-8) and **torus-band census** (`DIAG_TORUSBAND`, class `filletStripOnTorusBand`; D-130-7(c)): measurement only. Plate max |R_built−R_lsq| 1.55e-7 mm, torus-band count 0; S11-b's 12 fillets stay off-torus.
- **Fixtures and gates**: the Mastercam linkage plate `linkage_bores_chamfer` (49952bc) with a sidecar measured from the mesh; synthetic battery B2–B6 (`cross_bores`, `chamfer_straight_ring`, `boss_cone_chamfer`, `counterbore_chamfer`, `cyl_meets_chamfer`) via the corpus generator (beb6322; D-130-5). Contract owned by `gate_130` (`tests/gates/gate_130.py`; c49e169, 8a646f5; D-130-18(3)) — there is no `gates/linkage_chamfer_gate.py`. B2–B5 sidecars carry the generator's `exactVolume`; `gate_130`'s volume cell reads it (9ccaf8d; D-130-15(1)). The `coneNYI` fixture-name heuristic and sidecar entries are retired (6952dcc; D-130-15(2)). Unit targets `cone_math_unit`, `cone_bind_unit`, `quant_floor_unit`.

### Changed
- **Body28 builds** (D-130-6; 0361889): main reverted the whole component; on this branch it is a valid closed solid of 499 planes + 263 cylinders (2533 fewer faces than main's faceted answer), enabled by the two plate heals (7a1c822, e7c00a8). Its sidecar cell `smoothVolumeDeltaPctMax` is re-derived from the built shell as the exact measured round-trip value, 0.0 → 0.0033. Its six `analytic MakeEdge failed` chains class as `polylineTier2` (082194b; D-130-6 open).
- **Handle pickup heals, recorded** (D-130-11(1), D-130-13(1)): under the quantization-certified absorb two R≈3 bands admit three triangles that lie inside their own already-certified residual — 918 → 916 faces (PLANE 797 → 795; CYLINDRICAL 121), |V_step − V_mesh| 10.309374 → 10.115732, both radii moving to the least-squares value (4af42c2). Virtual generators then ship two more rims as circles, CIRCLE 95 → 97, |V_step − V_mesh| 10.115732 → 9.732831 (e2cf8f5, b1da5e4). The control is 916 / 795 / 121 / 97, |dV| 9.732831 from b1da5e4 on; PRG, plane-census, nwcyl and hl ratchets untouched and green.
- **Cascade containment** (66608dc; D-130-9(b)): a face that is not an analytic region explodes its own neighbourhood instead of escalating the component to the U2 blanket rung. RULE 1.5's 50 % test is unchanged.
- **RESULT `warnings` order** (b94a094; D-130-10(2)): the list handed back is sorted, so the RESULT line is a property of the mesh rather than of the thread pool; the live log callback still sees each warning at the moment it is raised.
- **RESULT `smoothBuiltCones`** (b94a094; D-130-10(1)): the conical faces a build actually shipped, counted off the written shape the way `smoothBuiltCylinders` is. Reads 0 on every part with no chamfer frustum. RESULT identity is compared modulo this key.
- **I7 loop closure** (b53029b, e523298; D-130-21): a boundary walk that leaves a split vertex and returns to the same vertex is a cycle. If that vertex is incident to ≥ 3 regions/islands it is a chain endpoint (I8), so the cycle is split there. PARKED I-checker rows on Body11, Body28 (regions 538/775) and handle-pickup (region 171) are lifted (038aafb; D-130-21).
- **B2 `cross_bores` pierced bore** (b8b3291, 488cce7; D-130-16, D-130-17, D-130-2): the R8 wall is one closed-360 cylinder face with two tier-2 inner wires (the R5 windows); the R5 surface is two faces (the larger bore severs every generator). `gate_130` counts a tier-2 edge as shared when both faces reference one EDGE_CURVE (8dd2503; D-130-19) and counts distinct cylindrical surfaces, not faces (53c46f6; D-130-20).
- **G2 sagitta in `archChainBand`** (75101e7; D-130-1): outside coarse, G2 accepts a cylinder commit only when the max vertex residual is within the exact deviation a correct fit of a **counted** N-gon can produce (`N = ev.d2.nSides`, only when `N ≥ 6`); never loosened for a region whose N was not counted. `coarseFusionBand` is not widened.
- **Plate fillet strips** (D-130-9(a), D-130-10(3); 1a2343d, 3d25757): once the real R 10 cross bore is a recognised cylinder, 52 two-facet strips C1 had claimed are the chamfer frustum, and 13 two-facet slivers on the toroidal mouth rounds lose fillet status (triangles ship as the same facets). `S11-b`'s 12 real fillets are byte-identical. The 1.4 torus item inherits the 13.
- **Cap plates carry their own measured planar residual** (`plateFaceTolCap = max(meshTolCap, Region::maxVertexDev)`; 7a1c822).
- **S03 builds** (74f885c vs main 3e7636f): main reverted it whole; the branch ships 24 planes + 30 cylinders.

### Fixed
- **A plane-loop CIRCLE never overrides the curve IntAna derived from the shipped surfaces** (74f885c; D-130-2): detector B had replaced 42 exact plane|plane lines on Handle pickup with arcs lying on neither face.
- **Facet fill planes** (0f065e7; D-130-9(b)): a mesh facet is filled on the plane of its own wire, so a shared vertex snapped onto a neighbouring region's curve no longer leaves the fill `UnorientableShape`.
- **One cone frame for both taper signs** (30f67dd): detector C had negated the neighbour cylinder's axis for a frustum whose radius grows against it, which also negates `YDirection` on a right-handed `gp_Ax3`; the sign now lives in a signed height; both plate frusta reach `valid=1`.
- **B2 volume** (9fa4eff; D-130-20): default `BRepGProp::VolumeProperties` under-integrated a cylindrical face whose inner wire is a 68-span polyline (37 mm³ / 0.0598 % vs the generator's exact volume). Adaptive 2D Gauss (`Eps = Precision::Confusion()`) recovers 1.70 mm³ / 0.00272 % — the D-130-2 tier-2 chord, counted not absorbed. Geometry unchanged; `gate_130` B2 volume cell green.
- **Single-closed-edge inner loops keep their sense on the edge, not the wire** (714a6af): closes the 1.2.0 registered limitation.
- **The plate `ShapeFix_Face` rung repairs orientation only, never wire geometry** (e7c00a8).
- **An ellipse arc binds its vertices at the parameters they occupy** (8c1498d): `S11-b` ships 12 cylinders / 12 fillets, one solid (1.2.0 re-read as 0 solids, volume 2.9674 % off).
- **Detector A never claimed a through hole** (ad5c814): `isGeneratorRidge` read a cap-rim link as a generator ridge ⟂ the axis.
- **Sidecar `semiAngleDeg` 45.017° → 45.000004°** on the plate (65d79a1): the recorded value was the facet-normal angle, not the surface's.

### Deferred
- **Torus surfaces → 1.4** (D-130-7, amended by D-130-9(a) and D-130-10(3)): the plate's three toroidal rounds at the cross-bore mouths stay faceted; no lane builds `Geom_ToroidalSurface` in 1.3.0. The 13 two-facet fillet-strip slivers on those bands are the 1.4 item's inherited count.
- **Same-surface islands and the loop-level union → 1.4, with the tori** (D-130-13(2), D-130-14, D-130-16, D-130-17 and the 2026-09-05 addendum). The island claiming landed and was reverted the same lane (228dfbc → 000c29f) because its per-piece emit shipped Handle pickup's single R 4 wall as eleven faces and `partial_recovery_gate` refused it — the instrument was not edited. The union is compiled-in behind `STL2STEP_UNION`, off by default (d3179af, 0477c72). 1.3.0 ships the seed piece(s) per surface: the plate's cross bore as **214 of 279** triangles in 2 faces; each R 3 slot-end top row stays facets (28 of 39 tris per wall in two faces). The plate sidecar's cross-bore entry is the mesh truth (255 on-surface + 24 off-surface, D-130-13(2)) and is **not** edited to match the shipped 214/279; the plate's `gates_full` SIDECAR row is **expected-red** by D-130-17.
- **B6 `cyl_meets_chamfer`** (a cone cut by a second bore: non-coaxial cone∩cyl and a trimmed cone face) stays `expected-red` and may slip to 1.3.1 (D-130-5). `gate_130` remains inverted on `GATE_130_EXPECTED_RED` while B6 is unmet.
- **Edge-class red lines on 9 fixtures** (082194b; D-130-2): Body11, Body28, S03, S04, S09, S11-b, S16-R1-explode-success, handle-pickup, linkage_bores_chamfer — recorded in `edge-class-ratchet.json` `expectedRed`; widening that map is a defect.
- **Intersection edges** (D-130-2): tier 1 analytic only where `IntAna` on the shipped surfaces yields a line, circle or ellipse; general cyl∩cyl and non-coaxial cyl∩cone ship as the mesh polyline, counted and ratcheted, never absorbed; free BSpline approximations and `BRepAlgoAPI_Section` curves are refused as shipped geometry (door D-130-2b).
- **Hand constants named, not moved** (D-130-8): `SegmentParams::thetaCylLoDeg = 5.0` still gates B1's seed and grow; `kTurnDot = 0.7` is unchanged.
- **The tripwire's resolver** (cd42b46, 8c1498d): reports but does not refuse until it is pass-aware.

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
