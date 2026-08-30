// Tool-local detectPrismatic stub (P1 owns the predicate). Weak so P1's
// definition wins on integrate. SPDX-License-Identifier: MIT

#include "refit_prism.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <gp_XYZ.hxx>

#if defined(__GNUC__) || defined(__clang__)
#define STL2STEP_WEAK __attribute__((weak))
#else
#define STL2STEP_WEAK
#endif

namespace stl2step {
namespace refit {

STL2STEP_WEAK PrismTols makePrismTols(const MeshView& mv, const RegionSet& rs) {
    PrismTols t;
    t.tauSurf = std::max({5e-5, 4.0 * mv.weldTol, 1e-6 * mv.diag});
    t.tauLvl = t.tauSurf;
    t.tauFit = t.tauSurf;
    double hMin = 0.0;
    for (const Region& r : rs.regions) {
        if (r.type != SurfType::Cylinder) continue;
        const double h = r.vMax - r.vMin;
        if (h > 0.0 && (hMin <= 0.0 || h < hMin)) hMin = h;
    }
    t.tauAx = (hMin > 0.0) ? std::max(1e-6, 2.0 * t.tauSurf / hMin) : 1e-6;
    return t;
}

STL2STEP_WEAK PrismLevels detectPrismatic(const MeshView& mv, const RegionSet& rs,
                                          const PrismTols& t) {
    (void)mv;
    PrismLevels lv;
    lv.ok = false;
    lv.failedCond = 1;

    PrismTols tols = t;
    if (!(tols.tauSurf > 0.0 && tols.tauAx > 0.0)) tols = makePrismTols(mv, rs);

    std::vector<const Region*> cyls;
    std::vector<const Region*> planes;
    for (const Region& r : rs.regions) {
        if (r.type == SurfType::Cylinder) cyls.push_back(&r);
        else if (r.type == SurfType::Plane) planes.push_back(&r);
    }
    if (cyls.size() < 2) {
        lv.failedCond = 1;
        return lv;
    }

    gp_XYZ sum(0.0, 0.0, 0.0);
    const gp_Dir ref = cyls[0]->ax.Direction();
    for (const Region* r : cyls) {
        gp_Dir d = r->ax.Direction();
        if (d.Dot(ref) < 0.0) d.Reverse();
        sum += d.XYZ();
    }
    if (sum.Modulus() < 1e-18) {
        lv.failedCond = 2;
        return lv;
    }
    sum.Normalize();
    if (sum.Y() < 0.0 && std::fabs(sum.Y()) >= std::fabs(sum.X()) &&
        std::fabs(sum.Y()) >= std::fabs(sum.Z()))
        sum.Reverse();
    lv.axis = gp_Dir(sum);
    const gp_XYZ ah = lv.axis.XYZ();

    struct Cap {
        double y = 0.0;
        int id = -1;
    };
    std::vector<Cap> caps;
    for (const Region* r : planes) {
        const double c = std::fabs(r->ax.Direction().XYZ().Dot(ah));
        if (c > 1.0 - tols.tauAx) {
            Cap h;
            h.y = ah.Dot(r->ax.Location().XYZ());
            h.id = r->id;
            caps.push_back(h);
        } else if (c >= tols.tauAx) {
            lv.failedCond = 3;
            return lv;
        }
    }
    if (caps.size() < 2) {
        lv.failedCond = 4;
        return lv;
    }
    std::sort(caps.begin(), caps.end(), [](const Cap& a, const Cap& b) {
        if (a.y != b.y) return a.y < b.y;
        return a.id < b.id;
    });
    lv.y.push_back(caps[0].y);
    lv.capRegion.push_back(caps[0].id);
    for (size_t i = 1; i < caps.size(); ++i) {
        if (caps[i].y - lv.y.back() > tols.tauLvl) {
            lv.y.push_back(caps[i].y);
            lv.capRegion.push_back(caps[i].id);
        }
    }
    if (lv.y.size() < 2) {
        lv.failedCond = 4;
        return lv;
    }
    lv.ok = true;
    lv.failedCond = 0;
    return lv;
}

}  // namespace refit
}  // namespace stl2step
