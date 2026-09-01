// stl2step command-line front-end. Parses arguments into stl2step::Options,
// runs the engine, streams progress, and prints the machine-readable RESULT line.
//
// SPDX-License-Identifier: MIT

#include "stl2step/stl2step.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace stl2step;

static void usage() {
    printf(
        "stl2step %s -- convert an STL mesh into a STEP (B-Rep) solid\n\n"
        "usage: stl2step <input.stl> [output.step] [options]\n"
        "       stl2step --mesh <in.step> [-o out.stl]  default <stem>.mesh.stl; never clobbers (pass -o)\n\n"
        "  -o <file>            output path (default: input with .step extension)\n"
        "  --schema <s>         AP203 | AP214 | AP242            (default AP214)\n"
        "  --units <mm|in>      units the STL was modelled in; in -> scaled x25.4 to mm\n"
        "  --scale <f>          extra scale factor applied to all coordinates\n"
        "  --weld <tol>         weld vertices within tol (default: exact duplicates only)\n"
        "  --sew-tol <tol>      tolerance for the sewing repair pass (default: auto)\n"
        "  --unify-angle <deg>  max normal deviation treated as coplanar (default 0.001;\n"
        "                       0 = OCCT strict default)\n"
        "  --no-unify           keep one face per triangle (skip coplanar merging)\n"
        "  --no-solid           emit shells only (skip solid conversion)\n"
        "  --force-sew          route every component through the sewing repair path\n"
        "  --no-verify          skip re-reading the written STEP (use when the caller\n"
        "                       imports the file right away -- that import IS the check)\n"
        "  --engine <mode>      conversion mode: verbatim (default) | trueform\n"
        "                       verbatim = byte-faithful tessellation (same as --no-smooth)\n"
        "                       trueform = premium analytic reconstruction (--smooth)\n"
        "  --smooth             alias for --engine trueform\n"
        "  --refit              alias for --engine trueform\n"
        "  --no-smooth          alias for --engine verbatim\n"
        "  --smooth-tol <mm>    surface-fit tolerance in mm (default: auto)\n"
        "  --smooth-angle <deg> near-flat normal gate for segmentation (default 2.0)\n"
        "  --no-smooth-fillets  skip recovery of fillet strips as cylinders\n"
        "  --dxf <dir>          write one DXF per slab to <dir> (prismatic parts only; off by default)\n"
        "  --threads <n>        worker threads for parallel stages (default: auto)\n"
        "  --quiet              suppress progress output (RESULT line + errors only)\n"
        "  -v, --version        print version and exit\n"
        "  -h, --help           print this help and exit\n\n"
        "STEP files are written in millimetres. Default mode is Verbatim (faceted\n"
        "surfaces at STL resolution). TrueForm (--engine trueform) recovers analytic\n"
        "planes, cylinders, and fillets where the refit ladder succeeds.\n"
        "Exit codes: 0 ok, 2 ok-with-warnings, 1 failed. Last stdout line: RESULT {json}\n",
        version());
}

static bool isMeshRejectedFlag(const std::string& a) {
    return a == "--engine" || a == "--schema" || a == "--units"
        || a == "--scale" || a == "--weld" || a == "--sew-tol"
        || a == "--unify-angle" || a == "--no-unify" || a == "--no-solid"
        || a == "--force-sew" || a == "--no-verify" || a == "--smooth"
        || a == "--refit" || a == "--no-smooth" || a == "--smooth-tol"
        || a == "--smooth-angle" || a == "--no-smooth-fillets" || a == "--dxf";
}

// Mesh mode is selected by the presence of --mesh and is dispatched before
// convert() is ever called. Rejected convert-only flags fail closed with no
// MESH_RESULT / RESULT line (plan-audit R4).
static int runMeshMode(int argc, char** argv) {
    MeshOptions opt;
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto val = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) { fprintf(stderr, "error: %s needs a value\n", flag); exit(1); }
            return argv[++i];
        };
        if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (a == "-v" || a == "--version") { printf("stl2step %s\n", version()); return 0; }
        else if (a == "--mesh")           opt.input = val("--mesh");
        else if (a == "-o")               opt.output = val("-o");
        else if (a == "--edges")          opt.edgesFile = val("--edges");
        else if (a == "--threads")        opt.threads = atoi(val("--threads").c_str());
        else if (a == "--quiet")          quiet = true;
        else if (isMeshRejectedFlag(a)) {
            fprintf(stderr, "error: %s is not valid with --mesh\n", a.c_str());
            return 1;
        }
        else if (a.size() && a[0] == '-') {
            fprintf(stderr, "error: unknown option %s\n", a.c_str());
            return 1;
        }
        else {
            fprintf(stderr, "error: unexpected argument %s\n", a.c_str());
            return 1;
        }
    }
    if (opt.input.empty()) {
        fprintf(stderr, "error: --mesh needs a value\n");
        return 1;
    }

    auto logcb = [&](Severity sev, const std::string& msg) {
        if (sev == Severity::Info) {
            if (!quiet) { fputs(msg.c_str(), stdout); fputc('\n', stdout); fflush(stdout); }
        } else if (sev == Severity::Warning) {
            fprintf(stderr, "warning: %s\n", msg.c_str());
        } else {
            fprintf(stderr, "error: %s\n", msg.c_str());
        }
    };

    MeshResult r = meshFromStep(opt, logcb);
    printf("MESH_RESULT %s\n", r.toJson().c_str());
    return r.exitCode;
}

int main(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--mesh")
            return runMeshMode(argc, argv);
    }

    Options opt;
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto val = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) { fprintf(stderr, "error: %s needs a value\n", flag); exit(1); }
            return argv[++i];
        };
        auto posDouble = [&](const char* flag) -> double {
            std::string s = val(flag);
            char* end = nullptr;
            double v = strtod(s.c_str(), &end);
            if (end == s.c_str() || *end != '\0' || v < 0.0) {
                fprintf(stderr, "error: %s must be a non-negative number\n", flag);
                exit(1);
            }
            return v;
        };
        auto angleDeg = [&](const char* flag) -> double {
            std::string s = val(flag);
            char* end = nullptr;
            double v = strtod(s.c_str(), &end);
            if (end == s.c_str() || *end != '\0' || v <= 0.0 || v > 90.0) {
                fprintf(stderr, "error: %s must be between 0 and 90\n", flag);
                exit(1);
            }
            return v;
        };
        if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (a == "-v" || a == "--version") { printf("stl2step %s\n", version()); return 0; }
        else if (a == "-o")              opt.output = val("-o");
        else if (a == "--schema")      { if (!parseSchema(val("--schema"), opt.schema)) {
                                             fprintf(stderr, "error: --schema must be AP203, AP214 or AP242\n");
                                             return 1; } }
        else if (a == "--units")       { std::string u = val("--units");
                                         if (u == "in" || u == "inch") opt.inchInput = true;
                                         else if (u != "mm") { fprintf(stderr, "error: --units mm|in\n"); return 1; } }
        else if (a == "--scale")         opt.scale *= atof(val("--scale").c_str());
        else if (a == "--weld")          opt.weldTol = atof(val("--weld").c_str());
        else if (a == "--sew-tol")       opt.sewTol = atof(val("--sew-tol").c_str());
        else if (a == "--unify-angle")   opt.unifyAngleDeg = atof(val("--unify-angle").c_str());
        else if (a == "--no-unify")      opt.unify = false;
        else if (a == "--no-solid")      opt.makeSolids = false;
        else if (a == "--force-sew")     opt.forceSew = true;
        else if (a == "--no-verify")     opt.verify = false;
        else if (a == "--engine") {
            std::string mode = val("--engine");
            if (mode == "verbatim" || mode == "faceted") opt.smooth = false;
            else if (mode == "trueform" || mode == "smooth") opt.smooth = true;
            else {
                fprintf(stderr, "error: --engine must be verbatim or trueform\n");
                return 1;
            }
        }
        else if (a == "--smooth" || a == "--refit") opt.smooth = true;
        else if (a == "--no-smooth")     opt.smooth = false;
        else if (a == "--smooth-tol")    opt.smoothTolMM = posDouble("--smooth-tol");
        else if (a == "--smooth-angle")  opt.smoothAngleDeg = angleDeg("--smooth-angle");
        else if (a == "--no-smooth-fillets") opt.smoothFillets = false;
        else if (a == "--dxf")           opt.dxfDir = val("--dxf");
        else if (a == "--threads")       opt.threads = atoi(val("--threads").c_str());
        else if (a == "--quiet")         quiet = true;
        else if (a.size() && a[0] == '-') { fprintf(stderr, "error: unknown option %s\n", a.c_str()); return 1; }
        else if (opt.input.empty())      opt.input = a;
        else if (opt.output.empty())     opt.output = a;
        else { fprintf(stderr, "error: unexpected argument %s\n", a.c_str()); return 1; }
    }
    if (opt.input.empty()) { usage(); return 1; }

    // Progress -> stdout (unless --quiet); warnings/errors -> stderr (always).
    auto logcb = [&](Severity sev, const std::string& msg) {
        if (sev == Severity::Info) {
            if (!quiet) { fputs(msg.c_str(), stdout); fputc('\n', stdout); fflush(stdout); }
        } else if (sev == Severity::Warning) {
            fprintf(stderr, "warning: %s\n", msg.c_str());
        } else {
            fprintf(stderr, "error: %s\n", msg.c_str());
        }
    };

    Result r = convert(opt, logcb);

    // The final stdout line is always the machine-readable contract.
    printf("RESULT %s\n", r.toJson().c_str());
    return r.exitCode;
}
