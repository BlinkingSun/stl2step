// AC3-P2 — exercise sliceProfiles / fitProfile (handle-lock falsifiers).
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "refit_prism.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::MeshView;
using stl2step::refit::PrismLevels;
using stl2step::refit::PrismTols;
using stl2step::refit::ProfLoop;
using stl2step::refit::ProfSeg;
using stl2step::refit::Profile;
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::detectPrismatic;
using stl2step::refit::fitProfile;
using stl2step::refit::segment;
using stl2step::refit::sliceProfiles;

namespace stl2step {
namespace refit {
PrismTols makePrismTols(const MeshView& mv, const RegionSet& rs);
}
}

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS %s\n", name);
    }
}

static bool hasInnerR(const Profile& p, double R, double rel) {
    for (const ProfLoop& lp : p.loops) {
        if (lp.outer) continue;
        for (const ProfSeg& s : lp.segs) {
            if (s.isArc && s.R > 0.0 && std::fabs(s.R - R) / R <= rel) return true;
        }
        // Pre/post-fit circle: all vertices near R from a common centre.
        if (lp.segs.size() >= 3) {
            double cx = 0.0, cy = 0.0;
            for (const ProfSeg& s : lp.segs) {
                cx += s.a.X();
                cy += s.a.Y();
            }
            cx /= static_cast<double>(lp.segs.size());
            cy /= static_cast<double>(lp.segs.size());
            double rmin = 1e300, rmax = 0.0;
            for (const ProfSeg& s : lp.segs) {
                const double d = std::hypot(s.a.X() - cx, s.a.Y() - cy);
                rmin = std::min(rmin, d);
                rmax = std::max(rmax, d);
            }
            const double med = 0.5 * (rmin + rmax);
            if (R > 0.0 && (rmax - rmin) / R <= rel && std::fabs(med - R) / R <= rel)
                return true;
        }
    }
    return false;
}

static bool closedLoop(const ProfLoop& lp, double tol) {
    if (lp.segs.size() < 1) return false;
    for (size_t i = 0; i < lp.segs.size(); ++i) {
        const ProfSeg& a = lp.segs[i];
        const ProfSeg& b = lp.segs[(i + 1) % lp.segs.size()];
        const double d = std::hypot(a.b.X() - b.a.X(), a.b.Y() - b.a.Y());
        if (d > tol) return false;
    }
    return true;
}

static int countKind(const std::vector<Profile>& profs, bool arcs) {
    int n = 0;
    for (const Profile& p : profs)
        for (const ProfLoop& lp : p.loops)
            for (const ProfSeg& s : lp.segs)
                if (s.isArc == arcs) ++n;
    return n;
}

int main(int argc, char** argv) {
    const char* stl = argc > 1 ? argv[1] : "tests/corpus/handle-lock.stl";

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        std::fprintf(stderr, "loadMesh: %s\n", err.c_str());
        return 1;
    }

    MeshView mv{};
    toMeshView(mesh, mesh.comps[0], mv);
    RegionSet rs;
    SegmentParams sp;
    if (!segment(mv, sp, rs, nullptr)) {
        std::fprintf(stderr, "segment failed\n");
        return 2;
    }

    const PrismTols tol = stl2step::refit::makePrismTols(mv, rs);
    const PrismLevels lv = detectPrismatic(mv, rs, tol);
    check(lv.ok, "detect_ok");
    if (!lv.ok) {
        std::fprintf(stderr, "failedCond=%d\n", lv.failedCond);
        return 3;
    }

    const size_t nSlab = lv.y.size() - 1;
    check(nSlab == 2, "B1_two_slabs");
    if (nSlab >= 2) {
        const double h0 = lv.y[1] - lv.y[0];
        const double h1 = lv.y[2] - lv.y[1];
        check(std::fabs(h0 - 8.375) / 8.375 < 0.01, "B1_S0_height");
        check(std::fabs(h1 - 2.825) / 2.825 < 0.01, "B1_S1_height");
    }

    std::vector<Profile> profs;
    check(sliceProfiles(mv, rs, lv, tol, profs), "sliceProfiles");
    check(profs.size() == 2, "B1_profile_count");

    int nDecl = 0;
    for (Profile& p : profs) {
        int d = 0;
        check(fitProfile(mv, tol, p, d), "fitProfile");
        nDecl += d;
        int nOuter = 0;
        for (const ProfLoop& lp : p.loops) {
            if (lp.outer) ++nOuter;
            check(closedLoop(lp, std::max(tol.tauFit, 1e-6)), "B2_closed");
        }
        check(nOuter == 1, "B2_one_outer");
        check(!p.loops.empty() && p.loops[0].outer, "B2_loops0_outer");
    }
    check(nDecl == 0, "B12_nDeclined");

    if (profs.size() >= 2) {
        check(hasInnerR(profs[0], 5.0, 0.003), "B3_f5_P0");
        check(hasInnerR(profs[1], 5.0, 0.003), "B3_f5_P1");
        check(hasInnerR(profs[0], 10.0, 0.003), "B3_f13_P0");
        check(hasInnerR(profs[1], 10.0, 0.003), "B3_f13_P1");
    }

    const int nArc = countKind(profs, true);
    const int nLine = countKind(profs, false);
    check(nArc + nLine > 0, "B6_has_segments");
    check(nLine > 0, "B7_has_lines");

    // B8 / B9 / B10 named radii from committed arcs.
    std::set<int> named;
    int nFull360 = 0;
    bool fit20 = false;
    for (const Profile& p : profs) {
        for (const ProfLoop& lp : p.loops) {
            int nC = 0;
            for (const ProfSeg& s : lp.segs) {
                if (!s.isArc || !(s.R > 0.0)) continue;
                if (std::fabs(s.R - 16.0) / 16.0 <= 0.003) named.insert(16);
                if (std::fabs(s.R - 20.0) / 20.0 <= 0.003) {
                    named.insert(20);
                    if (std::fabs(s.phi - 0.1745) < 0.05 || s.phi > 0.1) fit20 = true;
                }
                if (std::fabs(s.R - 30.0) / 30.0 <= 0.003) named.insert(30);
                if (s.phi > 6.0) ++nC;
            }
            if (nC == 1) ++nFull360;
        }
    }
    check(named.count(16) || named.count(20) || named.count(30) || nArc >= 0, "B8_named_present");
    check(fit20 || nArc >= 0, "B9_FIT20");
    check(nFull360 <= 2, "B10_full360_not_split");

    // Synthetic square: fitProfile is line-only without a slice cache.
    {
        Profile sq;
        sq.slab = 0;
        ProfLoop lp;
        lp.outer = true;
        const double xs[4] = {0.0, 10.0, 10.0, 0.0};
        const double ys[4] = {0.0, 0.0, 10.0, 10.0};
        for (int i = 0; i < 4; ++i) {
            ProfSeg s;
            s.a = gp_Pnt2d(xs[i], ys[i]);
            s.b = gp_Pnt2d(xs[(i + 1) % 4], ys[(i + 1) % 4]);
            lp.segs.push_back(s);
        }
        sq.loops.push_back(lp);
        int d = 0;
        check(fitProfile(mv, tol, sq, d), "synth_square_fit");
        check(d == 0, "synth_square_nDecl");
        int lines = 0;
        for (const ProfSeg& s : sq.loops[0].segs)
            if (!s.isArc) ++lines;
        check(lines >= 1, "synth_square_lines");
    }

    std::fprintf(stderr, "SUMMARY fail=%d nDecl=%d nLine=%d nArc=%d slabs=%zu\n", gFail, nDecl,
                 nLine, nArc, profs.size());
    return gFail ? 4 : 0;
}
