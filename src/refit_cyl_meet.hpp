// Partial counterbore walls (D-fc-5.17). Predicates only; no face builder.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_CYL_MEET_HPP
#define STL2STEP_REFIT_CYL_MEET_HPP

#include "refit_relational.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace stl2step {
namespace refit {

// RULE 2.2 volume budget, the two numbers the cascade already used as bare
// literals (shell: max(kTightRelFloor * |mesh|, kCascadeResidFactor * shipped
// |dVol|); per region: the same with |mesh| / nAnalytic). Named, not moved.
constexpr double kTightRelFloor = 1e-4;
constexpr double kCascadeResidFactor = 3.0;

// The RULE 2.2 shell budget: max(kTightRelFloor * |mesh|, kCascadeResidFactor * shipped |dVol|).
inline double shellVolBudget(double meshAbs, double shippedAbs) {
    return std::max(kTightRelFloor * meshAbs, kCascadeResidFactor * shippedAbs);
}

// max(kTightRelFloor * |mesh| / nAnalytic, kCascadeResidFactor * |dVol|).
// nAnalytic < 1 contributes no floor term.
inline double partialSheetVolBudget(double meshAbs, int nAnalytic, double dVolAbs) {
    const double rel = (nAnalytic > 0 && std::isfinite(meshAbs))
                           ? kTightRelFloor * std::fabs(meshAbs) / static_cast<double>(nAnalytic)
                           : 0.0;
    const double pred = std::isfinite(dVolAbs) ? kCascadeResidFactor * std::fabs(dVolAbs) : 0.0;
    return std::max(rel, pred);
}

// |((vFace - vChord) - dVolPredicted)| <= budget. Same residual the cascade scores.
inline bool partialSheetVolumeAgrees(double vFace, double vChord, double dVolPredicted,
                                     double budget) {
    if (!std::isfinite(vFace) || !std::isfinite(vChord) || !std::isfinite(dVolPredicted) ||
        !std::isfinite(budget) || !(budget >= 0.0))
        return false;
    return std::fabs((vFace - vChord) - dVolPredicted) <= budget;
}

// Cylindrical sector the stored radius claims and the certified radius does not.
inline double cylinderSectorGap(double rStored, double rFit, double uSpan, double vSpan) {
    if (!std::isfinite(rStored) || !std::isfinite(rFit) || !std::isfinite(uSpan) ||
        !std::isfinite(vSpan))
        return 0.0;
    if (!(uSpan > 0.0) || !(vSpan != 0.0)) return 0.0;
    return 0.5 * std::fabs(rStored * rStored - rFit * rFit) * uSpan * std::fabs(vSpan);
}

// |a × b| for unit directions. True when the axes separate by more than τ
// across a measured length.
inline bool axesCross(double crossMag, double tau, double length) {
    if (!(tau > 0.0) || !(length > tau) || !std::isfinite(crossMag) || !std::isfinite(length))
        return false;
    return crossMag > tau / length;
}

// Cap plane: its normal is parallel to the cylinder axis within τ across
// `extent` (sine of the tilt times that extent).
inline bool planeNormalParallelAxis(double absDot, double extent, double tau) {
    if (!(tau > 0.0) || !(extent > 0.0) || !std::isfinite(absDot) || !std::isfinite(extent))
        return false;
    const double d = std::min(1.0, std::fabs(absDot));
    return std::sqrt(std::max(0.0, 1.0 - d * d)) * extent <= tau;
}

// Outer radius of a bore-cut wall about a certified axis: the maximum rho.
// Inliers are samples within τ of that radius (the wall vertices). Samples on
// the cut lie inside it and are not inliers. Certified when the inliers meet
// minStations(kParamsCircle) and span two axial levels separated by τ.
struct RhoCertificate {
    double R = 0;
    int inliers = 0;
    int levels = 0;
    bool ok = false;
};

inline bool certifyOuterRho(const double* rho, const double* z, int n, double tau,
                            RhoCertificate& out) {
    out = RhoCertificate{};
    if (!rho || !z || n < minStations(kParamsCircle) || !(tau > 0.0)) return false;
    double R = -1.0;
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(rho[i]) || !std::isfinite(z[i]) || !(rho[i] >= 0.0)) return false;
        if (rho[i] > R) R = rho[i];
    }
    if (!(R > 0.0) || !std::isfinite(R)) return false;
    out.R = R;
    std::vector<double> lev;
    lev.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (std::fabs(rho[i] - R) <= tau) lev.push_back(z[i]);
    }
    out.inliers = static_cast<int>(lev.size());
    std::sort(lev.begin(), lev.end());
    out.levels = 0;
    for (std::size_t i = 0; i < lev.size(); ++i) {
        if (i == 0 || lev[i] - lev[i - 1] > tau) ++out.levels;
    }
    out.ok = out.inliers >= minStations(kParamsCircle) && out.levels >= 2;
    return out.ok;
}

struct CylMeetQuantities {
    bool cylinder = false;
    bool closed360 = false;
    double tau = 0;
    double storedDev = 0;
    double rStored = 0;
    double rFit = 0;
    bool fitOk = false;
    bool crossingMate = false;
    double gap = 0;
    double budget = 0;
};

// Intervene only for an open cylinder whose stored surface misses τ, whose
// certified section is inside that stored radius (the sheet claims area the
// vertices do not), which meets a crossing mate, and whose sector gap exceeds
// the sheet budget.
inline bool cylMeetIntervene(const CylMeetQuantities& q) {
    if (!q.cylinder || q.closed360) return false;
    if (!(q.tau > 0.0) || !std::isfinite(q.tau)) return false;
    if (!std::isfinite(q.storedDev) || !(q.storedDev > q.tau)) return false;
    if (!q.fitOk || !q.crossingMate) return false;
    if (!(q.rStored > q.rFit) || !std::isfinite(q.rStored) || !std::isfinite(q.rFit)) return false;
    if (!std::isfinite(q.gap) || !std::isfinite(q.budget)) return false;
    return q.gap > q.budget;
}

// Mid-azimuth of a partial wall, in (-pi, pi]. The covered arc is the
// complement of the largest gap, including the wrap from the last sample
// to the first. `centre` is that measured midpoint; `span` is the arc.
struct SectorCentre {
    double centre = 0;
    double span = 0;
    bool ok = false;
};

inline bool cylMeetSectorCentre(const double* ang, int n, SectorCentre& out) {
    out = SectorCentre{};
    if (!ang || n < 2) return false;
    const double turn = 6.28318530717958647692;
    const double half = 0.5 * turn;
    std::vector<double> a;
    a.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(ang[i])) return false;
        double t = ang[i];
        while (t > half) t -= turn;
        while (t <= -half) t += turn;
        a.push_back(t);
    }
    std::sort(a.begin(), a.end());
    double maxGap = (a.front() + turn) - a.back();
    int after = -1;
    for (int i = 0; i + 1 < n; ++i) {
        const double g = a[static_cast<std::size_t>(i + 1)] - a[static_cast<std::size_t>(i)];
        if (g > maxGap) {
            maxGap = g;
            after = i;
        }
    }
    const double start = (after < 0) ? a.front() : a[static_cast<std::size_t>(after + 1)];
    const double end = (after < 0) ? a.back() : a[static_cast<std::size_t>(after)] + turn;
    const double span = end - start;
    if (!(span > 0.0) || !(span < turn)) return false;
    double centre = start + 0.5 * span;
    while (centre > half) centre -= turn;
    while (centre <= -half) centre += turn;
    out.centre = centre;
    out.span = span;
    out.ok = true;
    return true;
}

// The segment chord census (refit_chains.cpp chordDevToRegion): distance of
// the chord midpoint (p + q) / 2 from the plane through ax, or from the
// cylinder of radius R about ax.
inline double chordMidDevPlane(const gp_Ax3& ax, const gp_XYZ& p, const gp_XYZ& q) {
    const gp_XYZ m((p.X() + q.X()) * 0.5, (p.Y() + q.Y()) * 0.5, (p.Z() + q.Z()) * 0.5);
    return std::fabs(ax.Direction().XYZ().Dot(m - ax.Location().XYZ()));
}

inline double chordMidDevCylinder(const gp_Ax3& ax, double R, const gp_XYZ& p, const gp_XYZ& q) {
    const gp_XYZ m((p.X() + q.X()) * 0.5, (p.Y() + q.Y()) * 0.5, (p.Z() + q.Z()) * 0.5);
    return std::fabs((m - ax.Location().XYZ()).Crossed(ax.Direction().XYZ()).Modulus() - R);
}

}  // namespace refit
}  // namespace stl2step

#endif
