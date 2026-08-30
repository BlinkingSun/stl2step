// Dump Profile as JSON for the AC3-P2 adjudicator.
// SPDX-License-Identifier: MIT

#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "refit_prism.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
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

static double loopR(const ProfLoop& lp) {
    double best = 0.0;
    for (const ProfSeg& s : lp.segs)
        if (s.isArc && s.R > best) best = s.R;
    return best;
}

static void printSeg(const ProfSeg& s) {
    std::printf("{\"kind\":\"%s\",\"R\":%.9f,\"phi\":%.9f,\"decl\":%s,"
                "\"a\":[%.9f,%.9f],\"b\":[%.9f,%.9f]",
                s.isArc ? "arc" : "line", s.R, s.phi, s.declinedAmbiguous ? "true" : "false",
                s.a.X(), s.a.Y(), s.b.X(), s.b.Y());
    if (s.isArc)
        std::printf(",\"center\":[%.9f,%.9f],\"ccw\":%s", s.center.X(), s.center.Y(),
                    s.ccw ? "true" : "false");
    std::printf("}");
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: profile_dump <stl> [--component N] [--json]\n");
        return 1;
    }
    const char* stl = argv[1];
    int comp = 0;
    bool json = false;
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--json") == 0) json = true;
        if (std::strcmp(argv[i], "--component") == 0 && i + 1 < argc)
            comp = std::atoi(argv[++i]);
    }

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        std::fprintf(stderr, "loadMesh: %s\n", err.c_str());
        return 2;
    }
    if (comp < 0 || static_cast<size_t>(comp) >= mesh.comps.size()) {
        std::fprintf(stderr, "bad component %d\n", comp);
        return 3;
    }

    MeshView mv{};
    toMeshView(mesh, mesh.comps[static_cast<size_t>(comp)], mv);
    RegionSet rs;
    SegmentParams sp;
    if (!segment(mv, sp, rs, nullptr)) {
        std::fprintf(stderr, "segment failed\n");
        return 4;
    }

    const PrismTols tol = stl2step::refit::makePrismTols(mv, rs);
    const PrismLevels lv = detectPrismatic(mv, rs, tol);
    if (!lv.ok) {
        std::fprintf(stderr, "not prismatic failedCond=%d\n", lv.failedCond);
        return 5;
    }

    std::vector<Profile> profs;
    if (!sliceProfiles(mv, rs, lv, tol, profs)) {
        std::fprintf(stderr, "sliceProfiles failed\n");
        return 6;
    }

    int nDecl = 0;
    for (Profile& p : profs) {
        int d = 0;
        if (!fitProfile(mv, tol, p, d)) {
            std::fprintf(stderr, "fitProfile failed slab=%d\n", p.slab);
            return 7;
        }
        nDecl += d;
    }

    if (!json) {
        std::printf("profiles=%zu nDeclined=%d levels=%zu\n", profs.size(), nDecl, lv.y.size());
        for (const Profile& p : profs) {
            std::printf("slab=%d loops=%zu\n", p.slab, p.loops.size());
            for (size_t i = 0; i < p.loops.size(); ++i) {
                const ProfLoop& lp = p.loops[i];
                int nLine = 0, nArc = 0;
                for (const ProfSeg& s : lp.segs) {
                    if (s.isArc) ++nArc;
                    else ++nLine;
                }
                std::printf("  loop=%zu outer=%d area=%.6f R=%.4f nLine=%d nArc=%d\n", i,
                            lp.outer ? 1 : 0, lp.area, loopR(lp), nLine, nArc);
            }
        }
        return 0;
    }

    std::printf("{\"ok\":true,\"nDeclined\":%d,\"levels\":[", nDecl);
    for (size_t i = 0; i < lv.y.size(); ++i) {
        if (i) std::printf(",");
        std::printf("%.9f", lv.y[i]);
    }
    std::printf("],\"profiles\":[");
    for (size_t pi = 0; pi < profs.size(); ++pi) {
        if (pi) std::printf(",");
        const Profile& p = profs[pi];
        std::printf("{\"slab\":%d,\"loops\":[", p.slab);
        for (size_t li = 0; li < p.loops.size(); ++li) {
            if (li) std::printf(",");
            const ProfLoop& lp = p.loops[li];
            std::printf("{\"outer\":%s,\"area\":%.9f,\"R\":%.9f,\"segs\":[",
                        lp.outer ? "true" : "false", lp.area, loopR(lp));
            for (size_t si = 0; si < lp.segs.size(); ++si) {
                if (si) std::printf(",");
                printSeg(lp.segs[si]);
            }
            std::printf("]}");
        }
        std::printf("]}");
    }
    std::printf("]}\n");
    return 0;
}
