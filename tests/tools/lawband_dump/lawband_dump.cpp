// Throwaway: load handle-lock + band-anatomy, run lawChainAccept / lawCalibrate.
// SPDX-License-Identifier: MIT

#include "anatomy_io.hpp"
#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "refit.hpp"
#include "refit_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

using stl2step::harness::HarnessMesh;
using stl2step::harness::loadMesh;
using stl2step::harness::toMeshView;
using stl2step::refit::DerivedTols;
using stl2step::refit::LawBand;
using stl2step::refit::MeshView;
using stl2step::refit::SegmentParams;
using stl2step::refit::TessLawInterval;
using stl2step::refit::lawCalibrate;
using stl2step::refit::lawChainAccept;

static DerivedTols makeTols(const MeshView& mv) {
    DerivedTols t;
    t.epsMesh = std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
    t.epsPlane = std::max({t.epsMesh, mv.sewTol, 0.02});
    return t;
}

static std::vector<int> toLocal(const MeshView& mv, const std::vector<int>& global) {
    std::vector<int> loc2g(mv.nTri, -1);
    for (size_t i = 0; i < mv.nTri; ++i)
        if (mv.compTris) loc2g[i] = mv.compTris[i];
    std::vector<int> g2l(4096, -1);
    int mx = 0;
    for (size_t i = 0; i < mv.nTri; ++i) mx = std::max(mx, loc2g[i]);
    if (mx >= 0) g2l.assign(static_cast<size_t>(mx + 1), -1);
    for (size_t i = 0; i < mv.nTri; ++i)
        if (loc2g[i] >= 0 && static_cast<size_t>(loc2g[i]) < g2l.size())
            g2l[static_cast<size_t>(loc2g[i])] = static_cast<int>(i);
    std::vector<int> out;
    out.reserve(global.size());
    for (int g : global) {
        if (g >= 0 && static_cast<size_t>(g) < g2l.size() && g2l[static_cast<size_t>(g)] >= 0)
            out.push_back(g2l[static_cast<size_t>(g)]);
        else
            out.push_back(g);
    }
    return out;
}

int main(int argc, char** argv) {
    const char* stl = argc > 1 ? argv[1] : "tests/corpus/handle-lock.stl";
    const char* anat = argc > 2 ? argv[2]
                                : "tests/gates/labels/handle-lock.band-anatomy.json";
    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        std::fprintf(stderr, "loadMesh failed: %s\n", err.c_str());
        return 1;
    }
    MeshView mv{};
    toMeshView(mesh, mesh.comps[0], mv);
    const DerivedTols tol = makeTols(mv);

    std::vector<AnatomyBand> bands;
    if (!loadAnatomy(anat, bands)) {
        std::fprintf(stderr, "loadAnatomy failed: %s\n", anat);
        return 2;
    }

    std::vector<LawBand> lb(bands.size());
    std::vector<int> acc(bands.size(), 0);
    std::vector<std::thread> pool;
    pool.reserve(bands.size());
    for (size_t i = 0; i < bands.size(); ++i) {
        pool.emplace_back([&, i]() {
            const std::vector<int> loc = toLocal(mv, bands[i].triIds);
            acc[i] = lawChainAccept(mv, loc, tol, lb[i]) ? 1 : 0;
        });
    }
    for (auto& t : pool) t.join();

    std::printf("face N Nhat R Rtrue rel%% phiDeg closed cvT cvR resid accept\n");
    int nAcc = 0;
    double maxRel = 0.0;
    for (size_t i = 0; i < bands.size(); ++i) {
        const double rel =
            (bands[i].radius > 0.0) ? std::abs(lb[i].R - bands[i].radius) / bands[i].radius : 0.0;
        maxRel = std::max(maxRel, rel);
        nAcc += acc[i];
        std::printf("F%d %d %d %.9f %.6f %.6f%% %.4f %d %.3g %.3g %.3g %d\n", bands[i].face,
                    bands[i].N, lb[i].N, lb[i].R, bands[i].radius, rel * 100.0,
                    lb[i].phi * 180.0 / 3.141592653589793, lb[i].closed360 ? 1 : 0, lb[i].cvTheta,
                    lb[i].cvR, lb[i].maxVertResid, acc[i]);
    }
    std::vector<LawBand> accepted;
    for (size_t i = 0; i < lb.size(); ++i)
        if (acc[i]) accepted.push_back(lb[i]);
    const TessLawInterval cal = lawCalibrate(accepted.empty() ? lb : accepted);
    std::printf("CALIB dLo=%.6f dHi=%.6f aLoDeg=%.4f aHiDeg=%.4f nD=%d nA=%d empty=%d\n", cal.dLo,
                cal.dHi, cal.alphaLo * 180.0 / 3.141592653589793,
                cal.alphaHi * 180.0 / 3.141592653589793, cal.nDLimited, cal.nAlphaLimited,
                cal.empty ? 1 : 0);
    std::printf("SUMMARY accept=%d/%zu maxRel=%.3e\n", nAcc, bands.size(), maxRel);
    return nAcc == static_cast<int>(bands.size()) ? 0 : 3;
}
