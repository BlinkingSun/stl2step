// stl2step - a universal mesh-to-B-Rep conversion engine.
//
// Converts a triangle mesh (STL) into a parametric B-Rep solid (STEP), suitable
// for import into CAD/CAM kernels that expect analytic boundary representations
// rather than raw triangle soup. The engine welds vertices, splits the mesh into
// manifold bodies, builds each body as a B-Rep in parallel, repairs dirty meshes
// (open edges / flipped facets), merges coplanar facets into single planar faces,
// fits tolerances to true geometric deviation, and writes a validated STEP file.
//
// This header is the entire public interface. Everything is in namespace
// stl2step. The engine is dependency-light at the API boundary: only the C++
// standard library appears here (OpenCASCADE is an implementation detail).
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_STL2STEP_HPP
#define STL2STEP_STL2STEP_HPP

#include <functional>
#include <string>
#include <vector>

namespace stl2step {

// -------------------------------------------------------------------- version

// Semantic version of the engine (also available as STL2STEP_VERSION macros
// below for preprocessor checks).
#define STL2STEP_VERSION_MAJOR 1
#define STL2STEP_VERSION_MINOR 3
#define STL2STEP_VERSION_PATCH 0
#define STL2STEP_VERSION_STRING "1.3.0"

// Returns the runtime version string ("1.0.0"). Useful for logging which engine
// build a host application linked against.
const char* version();

// ---------------------------------------------------------------- STEP schema

// STEP application protocol used for the written file. AP214 is the default and
// the safest choice for mechanical solids consumed by other CAD/CAM systems.
enum class Schema { AP203, AP214, AP242 };

// "AP203" | "AP214" | "AP242"  <->  enum, case-insensitive parse.
const char* schemaName(Schema s);
bool        parseSchema(const std::string& text, Schema& out);

// ------------------------------------------------------------------- options

// Everything that controls one conversion. Only `input` is required; the rest
// carry sensible defaults that match the standalone CLI.
struct Options {
    // Input STL path. Binary and ASCII STL are both auto-detected. Required.
    std::string input;

    // Output STEP path. If empty, the input path with its extension replaced by
    // ".step" is used (e.g. "part.stl" -> "part.step").
    std::string output;

    // STEP application protocol for the written file.
    Schema schema = Schema::AP214;

    // Extra uniform scale applied to every coordinate. STL carries no units, so
    // use this (or `inchInput`) to convert. Output is always written in mm.
    double scale = 1.0;

    // Convenience: treat the STL as modelled in inches (multiplies scale by 25.4).
    // Equivalent to setting scale = 25.4 yourself.
    bool inchInput = false;

    // Weld vertices whose coordinates fall in the same grid cell of this size.
    // 0 (default) welds only exact-bit duplicates. Raise it for meshes exported
    // with per-facet (unshared) vertices that need snapping.
    double weldTol = 0.0;

    // Tolerance for the sewing repair pass used on dirty (non-manifold / open /
    // flipped-facet) components. 0 (default) auto-derives it from the bounding
    // box diagonal.
    double sewTol = 0.0;

    // Maximum normal deviation (degrees) treated as coplanar when merging facets
    // into single faces. Default 0.001 deg is far above STL float32 normal noise
    // yet far below any real tessellation angle. 0 = OpenCASCADE's strict default.
    double unifyAngleDeg = 0.001;

    // Merge coplanar facets into single planar faces (recommended). When false,
    // every triangle stays its own face.
    bool unify = true;

    // Convert closed shells into solids (recommended). When false, shells are
    // written without solid conversion.
    bool makeSolids = true;

    // Route every component through the sewing repair path, even manifold ones.
    // Slower; useful only for diagnostics.
    bool forceSew = false;

    // Re-read the written STEP and compare its volume against the source as a
    // self-check. Set false when the caller is about to import the STEP anyway
    // (that import is itself the verification) - this can roughly halve wall time
    // on large files.
    bool verify = true;

    // Worker threads for the parallel stages. 0 (default) uses all hardware cores.
    int threads = 0;

    // Overrides the STEP product name written into the file. Empty = the output
    // file's stem (e.g. "part" for "part.step").
    std::string productName;

    // Recognise tessellated radii and near-flat regions and emit true analytic
    // Geom_CylindricalSurface / Geom_Plane faces instead of triangle facets.
    // Default false for the whole 1.x line; with smooth off the engine is
    // bit-identical to 1.0.0.
    bool smooth = false;

    // Surface-fit tolerance in mm. 0 auto-derives from the bounding box, weld
    // and sew tolerances. Only meaningful when smooth is true.
    double smoothTolMM = 0.0;

    // Near-flat normal gate in degrees (segmentation threshold). This is not
    // unifyAngleDeg (0.001 deg, a coplanarity test on already-flat facets);
    // 2.0 deg is a segmentation gate and the two are unrelated.
    double smoothAngleDeg = 2.0;

    // Recover plane-to-plane fillet strips as cylinders. Ignored when smooth
    // is false.
    bool smoothFillets = true;

    // Directory for one DXF per prismatic slab (`--dxf <dir>`). Empty (default)
    // leaves emission off: no directory is created and STEP is unchanged.
    std::string dxfDir;
};

// -------------------------------------------------------------------- logging

// Severity of a progress/diagnostic message reported through the log callback.
enum class Severity { Info, Warning, Error };

// Optional progress sink. The engine calls this for human-readable progress
// (Info), recoverable problems (Warning - also collected into Result::warnings),
// and fatal errors (Error - also placed in Result::error). Messages carry no
// trailing newline. Pass nullptr (the default) for a silent conversion. The
// callback may be invoked from worker threads, but the engine serialises calls,
// so it need not be thread-safe itself.
using LogCallback = std::function<void(Severity severity, const std::string& message)>;

// --------------------------------------------------------------------- result

// Outcome of a conversion. `ok` is the single success flag; `exitCode` mirrors
// the CLI contract (0 = clean, 2 = written with warnings, 1 = failed).
struct Result {
    bool ok = false;       // true when a STEP file was written (possibly with warnings)
    int  exitCode = 1;     // 0 clean, 2 written-with-warnings, 1 failed

    std::string input;     // resolved input path
    std::string output;    // resolved output path (written on success)
    std::string error;     // human-readable reason when ok == false

    int triangles = 0;     // triangles read from the STL
    int vertices  = 0;     // vertices after welding
    int components = 0;    // manifold bodies the mesh split into

    int solids     = 0;    // solids written
    int openShells = 0;    // components that could not close (written as open shells)

    int facesBeforeUnify = 0;  // face count before coplanar merge
    int facesAfterUnify  = 0;  // face count after coplanar merge

    double meshVolumeMM3 = 0.0;   // volume of the source mesh (mm^3)
    double stepVolumeMM3 = 0.0;   // volume re-read from the STEP (mm^3; 0 if verify off)
    double volumeDeltaPct = -1.0; // |step - brep| / brep * 100 (-1 if not measured)
    bool   watertight = true;     // every component closed with consistent winding

    double seconds = 0.0;         // wall-clock time of the whole conversion

    // D-130-10(2): warnings are emitted from parallel component builds, so
    // arrival order is not a property of the mesh. The list is sorted before
    // it is handed back, which makes the RESULT line deterministic; the live
    // log callback still sees each warning at the moment it is raised.
    std::vector<std::string> warnings;  // every Warning emitted, sorted

    // These C++ members are always present and default to zero, so host code can
    // read them unconditionally and ABI does not shift with a flag. The RESULT
    // JSON is the opposite: toJson() emits the smooth* keys only when
    // Options::smooth was true, so an off-path RESULT string stays
    // character-identical to 1.0.0. Keys are omitted when the feature is off,
    // never emitted as zeros.

    int smoothPlanes = 0;              // planar regions recovered as Geom_Plane
    int smoothCylinders = 0;           // cylindrical regions recovered
    int smoothFillets = 0;             // fillet strips recovered as cylinders
    int smoothDistinctRadii = 0;       // distinct cylinder radii accepted
    int smoothRejected = 0;            // candidate regions rejected by gates
    int smoothFacetFaces = 0;          // faceted faces left after smooth pass
    int facesAfterSmooth = 0;          // total face count after smooth pass
    int smoothSkippedComponents = 0;   // dirty components not refit

    double smoothMaxDevMM = 0.0;       // max vertex deviation from fit (mm)
    double smoothMaxEdgeTolMM = 0.0;   // max edge tolerance written (mm)
    double smoothVolPredictedMM3 = 0.0; // predicted volume from analytic fits (mm^3)

    int smoothBuiltPlanes = 0;
    int smoothBuiltCylinders = 0;
    // D-130-10(1): conical faces actually shipped on built components. Reads 0
    // on every part with no chamfer frustum; a new RESULT key, so an audit
    // comparing RESULT identity across this commit compares modulo it.
    int smoothBuiltCones = 0;
    int smoothBuiltFillets = 0;
    int smoothBuiltComponents = 0;
    int smoothRevertedComponents = 0;

    // D-130-2 edge-class census over the shipped shells: every edge two
    // analytic faces share, classified by the closed-form bind-site supremum
    // against the tolerance the edge records and the region's meshTolCap.
    // Emitted as the RESULT object "edgeClasses".
    int edgeClassAnalytic = 0;       // tier 1 (closed-form supremum both sides)
    int edgeClassPolylineTier2 = 0;  // tier 2 (mesh polyline, counted not absorbed)
    int edgeClassUnhandled = 0;      // no class, or an analytic edge no face names
    int edgeClassOverTol = 0;        // deviation above the recorded tolerance
    int edgeClassOverCap = 0;        // recorded tolerance above meshTolCap

    // Radius truth (SPEC-130-bind addendum): the shipped radius of every built
    // cylinder against the least-squares radius of its own claimed vertices.
    // Emitted as the RESULT object "radiusDrift". Measurement, not a gate.
    int    radiusDriftN = 0;
    double radiusDriftMaxAbs = 0.0;
    double radiusDriftMaxRel = 0.0;

    // The machine-readable payload the CLI prints after "RESULT ". Stable field
    // set and ordering; safe to parse. Does not include the "RESULT " prefix.
    std::string toJson() const;
};

// --------------------------------------------------------------------- engine

// Convert one mesh to a STEP file according to `opt`. `log` is an optional
// progress sink (see LogCallback). Never throws: OpenCASCADE failures and
// standard exceptions are caught and reported through Result::error / Result::ok.
//
// The call is self-contained and reentrant - each invocation owns its own state,
// so independent conversions may run on separate threads concurrently (the
// engine parallelises internally, so prefer converting one file at a time on a
// busy machine, or several at once on an idle one).
Result convert(const Options& opt, const LogCallback& log = nullptr);

// --------------------------------------------------------------- mesh-from-STEP

// Sibling of convert(): tessellate a STEP/B-Rep into a binary STL and an
// optional B-Rep edge polyline buffer (Format A). Does not call convert() and
// does not touch Result / Result::toJson().

struct MeshOptions {
    // Input STEP path (.step / .stp). Required.
    std::string input;

    // Output binary STL path. Empty = <input-stem>.mesh.stl beside the input.
    // A derived default never overwrites an existing file (fail closed:
    // error "output exists: <path> — pass -o"). An explicit path may overwrite.
    std::string output;

    // Optional Format A edges path. Empty = do not write an edges file.
    // Never derived: --edges is explicit-only (same never-clobber rule would
    // apply if a default were ever added).
    std::string edgesFile;

    // Worker threads for BRepMesh_IncrementalMesh. 0 (default) = all cores.
    // InParallel is true when the resolved thread count is greater than 1.
    int threads = 0;
};

struct MeshResult {
    bool ok = false;       // true when the binary STL was written
    int  exitCode = 1;     // 0 ok, 1 failed (no warnings tier)

    int faces = 0;         // TopExp::MapShapes(shape, TopAbs_FACE) extent
    int edges = 0;         // drawable (non-degenerate, polygon-bearing) edges
    int triangles = 0;     // triangles written into the binary STL

    double seconds = 0.0;  // wall-clock time of the whole mesh run

    std::string input;     // resolved input path
    std::string output;    // resolved STL path
    std::string edgesFile; // resolved edges path; empty if --edges was not passed
    std::string error;     // human-readable reason when ok == false

    // Machine-readable payload the CLI prints after "MESH_RESULT ".
    // Own writer: never shares Result::toJson().
    std::string toJson() const;
};

// Read a STEP file, tessellate it, write a binary STL, and optionally write a
// Format A edges buffer. Never throws: OCCT Standard_Failure and std::exception
// are caught and reported through MeshResult::error / MeshResult::ok.
MeshResult meshFromStep(const MeshOptions& opt, const LogCallback& log = nullptr);

}  // namespace stl2step

#endif  // STL2STEP_STL2STEP_HPP
