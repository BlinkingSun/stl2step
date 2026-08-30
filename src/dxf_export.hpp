// DXF emission — P4. Pure serializer + --emit-dxf plumbing.
// Types match the frozen §3 contract so this lane compiles without P1's header.
// If src/refit_prism.hpp is present, that header is authoritative (same types).
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_DXF_EXPORT_HPP
#define STL2STEP_DXF_EXPORT_HPP

#if __has_include("refit_prism.hpp")
#include "refit_prism.hpp"
#else

#include <string>
#include <vector>

#include <gp_Dir.hxx>
#include <gp_Pnt2d.hxx>

#include "refit.hpp"

namespace stl2step {
namespace refit {

struct PrismTols {
    double tauSurf = 0.0;   // mm: max(5e-5, 4*weldTol, 1e-6*diag)
    double tauLvl  = 0.0;   // mm: == tauSurf
    double tauFit  = 0.0;   // mm: == tauSurf
    double tauAx   = 0.0;   // rad: max(1e-6, 2*tauSurf/hMin)
};

struct PrismLevels {
    gp_Dir              axis;         // the common axis
    std::vector<double> y;            // sorted cap offsets n.p, size L >= 2
    std::vector<int>    capRegion;    // region id per level, size L
    bool                ok = false;   // RULE 5.1: all six conditions
    int                 failedCond = 0;  // 1..6, 0 when ok; reported on decline
};

struct ProfSeg {
    bool   isArc = false;
    gp_Pnt2d a, b;                    // endpoints, always set
    gp_Pnt2d center;                  // arc only
    double R = 0.0, phi = 0.0;        // arc only; R from the law inverse (RULE 5.3)
    bool   ccw = false;
    bool   declinedAmbiguous = false; // RULE 5.3a fired; emitted as line
};
struct ProfLoop { std::vector<ProfSeg> segs; bool outer = false; double area = 0.0; };
struct Profile  { int slab = -1; std::vector<ProfLoop> loops; };  // loops[0] outer, rest holes

}  // namespace refit
}  // namespace stl2step

#endif

#include <string>
#include <vector>

namespace stl2step {
namespace refit {

// P4. Pure serializer; never consulted by the engine.
bool writeProfileDxf(const Profile& p, const PrismLevels& lv, const std::string& path);

// Deterministic one-file-per-slab name: slab<KK>_y<yK>_y<yK+1>.dxf
std::string makeProfileDxfName(const Profile& p, const PrismLevels& lv);

// Write every profile into dir (created if needed). Parallel over slabs.
// Returns the number of files written, or -1 on failure.
int emitProfilesDxf(const std::vector<Profile>& profs, const PrismLevels& lv,
                    const std::string& dir);

// CLI plumbing for --emit-dxf <dir>. Returns true when `flag` is that option
// and `val` was stored in dirOut. main.cpp does not call this yet.
bool consumeEmitDxfFlag(const std::string& flag, const std::string& val,
                        std::string& dirOut);

}  // namespace refit
}  // namespace stl2step

#endif  // STL2STEP_DXF_EXPORT_HPP
