#ifndef STL2STEP_REFIT_PRISM_HPP
#define STL2STEP_REFIT_PRISM_HPP

#include "refit.hpp"

#include <gp_Dir.hxx>
#include <gp_Pnt2d.hxx>
#include <TopoDS_Shape.hxx>

#include <string>
#include <vector>

namespace stl2step {
namespace refit {

// src/refit_prism.hpp — Stage P (prism reconstruction), D7 §1-3.
// All tolerances are SELF-COMPUTED from the mesh (RULE 4.2a). No degree constants.
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

// One closed 2D loop in the sketch plane. Alphabet is CLOSED: line or arc.
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

// RULE 5.1 — predicate + level extraction. Never throws; sets ok=false and failedCond.
PrismLevels detectPrismatic(const MeshView& mv, const RegionSet& rs, const PrismTols& t);

// RULE 5.2 — one profile per slab. RULE 5.2a: a through-feature appears in EVERY slab it crosses.
bool sliceProfiles(const MeshView& mv, const RegionSet& rs, const PrismLevels& lv,
                   const PrismTols& t, std::vector<Profile>& out);

// RULE 5.3 — lines + arcs only; arc R from the law inverse. Ambiguous -> line, flagged.
bool fitProfile(const MeshView& mv, const PrismTols& t, Profile& p, int& nDeclined);

// RULE 5.2b — union of axially disjoint prisms. Caps are OUTPUTS (RULE 5.2c).
bool buildPrismSolid(const std::vector<Profile>& profs, const PrismLevels& lv,
                     TopoDS_Shape& out);

// P4. Pure serializer; never consulted by the engine.
bool writeProfileDxf(const Profile& p, const PrismLevels& lv, const std::string& path);

}  // namespace refit
}  // namespace stl2step

#endif
