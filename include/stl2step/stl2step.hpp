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
#define STL2STEP_VERSION_MINOR 1
#define STL2STEP_VERSION_PATCH 0
#define STL2STEP_VERSION_STRING "1.1.0-pre"

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

    std::vector<std::string> warnings;  // every Warning emitted, in order

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
    int smoothBuiltFillets = 0;
    int smoothBuiltComponents = 0;
    int smoothRevertedComponents = 0;

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

}  // namespace stl2step

#endif  // STL2STEP_STL2STEP_HPP
