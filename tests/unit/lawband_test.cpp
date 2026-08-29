// AC2-L1 — exercise lawChainAccept / lawCalibrate / lawNConsistent / lawBandsMergeable.
// SPDX-License-Identifier: MIT

#include "../tools/lawband_dump/anatomy_io.hpp"
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
using stl2step::refit::RegionSet;
using stl2step::refit::SegmentParams;
using stl2step::refit::SurfType;
using stl2step::refit::TessLawInterval;
using stl2step::refit::lawBandsMergeable;
using stl2step::refit::lawCalibrate;
using stl2step::refit::lawChainAccept;
using stl2step::refit::lawNConsistent;
using stl2step::refit::segment;

static int gFail = 0;

static void check(bool ok, const char* name) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        ++gFail;
    } else {
        std::fprintf(stderr, "PASS %s\n", name);
    }
}

static DerivedTols makeTols(const MeshView& mv) {
    DerivedTols t;
    t.epsMesh = std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
    t.epsPlane = std::max({t.epsMesh, mv.sewTol, 0.02});
    return t;
}

static std::vector<int> toLocal(const MeshView& mv, const std::vector<int>& global) {
    int mx = 0;
    for (size_t i = 0; i < mv.nTri; ++i)
        if (mv.compTris) mx = std::max(mx, mv.compTris[i]);
    std::vector<int> g2l(static_cast<size_t>(std::max(0, mx) + 1), -1);
    for (size_t i = 0; i < mv.nTri; ++i)
        if (mv.compTris && mv.compTris[i] >= 0 &&
            static_cast<size_t>(mv.compTris[i]) < g2l.size())
            g2l[static_cast<size_t>(mv.compTris[i])] = static_cast<int>(i);
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

static double inverseR(const std::vector<double>& w, const std::vector<double>& th) {
    std::vector<double> Ri;
    for (size_t i = 0; i < w.size() && i < th.size(); ++i) {
        const double s = std::sin(0.5 * th[i]);
        if (s > 1e-15) Ri.push_back(w[i] / (2.0 * s));
    }
    if (Ri.empty()) return 0.0;
    std::sort(Ri.begin(), Ri.end());
    const size_t n = Ri.size();
    return (n & 1u) ? Ri[n / 2] : 0.5 * (Ri[n / 2 - 1] + Ri[n / 2]);
}

static const AnatomyBand* findFace(const std::vector<AnatomyBand>& b, int face) {
    for (const auto& x : b)
        if (x.face == face) return &x;
    return nullptr;
}

static std::vector<int> byFace(const std::vector<FaceLabel>& labs, int face,
                               const std::vector<int>& pool) {
    int mx = 0;
    for (const auto& L : labs) mx = std::max(mx, L.tri);
    std::vector<int> faceOf(static_cast<size_t>(mx + 1), -1);
    for (const auto& L : labs)
        if (L.tri >= 0 && L.tri <= mx) faceOf[static_cast<size_t>(L.tri)] = L.face;
    std::vector<int> out;
    for (int t : pool) {
        if (t >= 0 && t <= mx && faceOf[static_cast<size_t>(t)] == face) out.push_back(t);
    }
    return out;
}

int main(int argc, char** argv) {
    const char* stl = argc > 1 ? argv[1] : "tests/corpus/handle-lock.stl";
    const char* anatP = argc > 2 ? argv[2] : "_team/reports/ac2/band-anatomy.json";
    const char* labP = argc > 3 ? argv[3] : "_team/reports/ac2/tri-labels.json";
    (void)argc;

    HarnessMesh mesh;
    std::string err;
    if (!loadMesh(stl, 1.0, 0.0, 0.0, mesh, err) || mesh.comps.empty()) {
        std::fprintf(stderr, "loadMesh: %s\n", err.c_str());
        return 2;
    }
    MeshView mv{};
    toMeshView(mesh, mesh.comps[0], mv);
    const DerivedTols tol = makeTols(mv);

    std::vector<AnatomyBand> bands;
    if (!loadAnatomy(anatP, bands) || bands.size() != 15) {
        std::fprintf(stderr, "anatomy: got %zu bands from %s\n", bands.size(), anatP);
        return 3;
    }
    std::vector<FaceLabel> labs;
    (void)loadTriLabels(labP, labs);

    // --- A1 formula: inverse on spike (w,θ) is machine-zero ---
    double maxInv = 0.0;
    for (const auto& b : bands) {
        std::vector<double> th;
        for (double d : b.thetaDeg) th.push_back(d * 3.141592653589793 / 180.0);
        const double R = inverseR(b.wMm, th);
        const double rel = (b.radius > 0.0) ? std::abs(R - b.radius) / b.radius : 0.0;
        maxInv = std::max(maxInv, rel);
    }
    std::fprintf(stderr, "A1a formula max|R-Rtrue|/Rtrue = %.3e\n", maxInv);
    check(maxInv < 1e-9, "A1a inverse formula < 1e-9");

    // --- lawChainAccept on all 15 (threaded) ---
    std::vector<LawBand> lb(bands.size());
    std::vector<int> acc(bands.size(), 0);
    {
        std::vector<std::thread> pool;
        pool.reserve(bands.size());
        for (size_t i = 0; i < bands.size(); ++i) {
            pool.emplace_back([&, i]() {
                acc[i] = lawChainAccept(mv, toLocal(mv, bands[i].triIds), tol, lb[i]) ? 1 : 0;
            });
        }
        for (auto& t : pool) t.join();
    }

    int nAcc = 0;
    double maxRel = 0.0, maxCvt = 0.0, maxResid = 0.0, maxAx = 0.0;
    bool face526 = false;
    for (size_t i = 0; i < bands.size(); ++i) {
        nAcc += acc[i];
        const double rel =
            (bands[i].radius > 0.0) ? std::abs(lb[i].R - bands[i].radius) / bands[i].radius : 1.0;
        maxRel = std::max(maxRel, rel);
        maxCvt = std::max(maxCvt, lb[i].cvTheta);
        maxResid = std::max(maxResid, lb[i].maxVertResid);
        const double axn = std::sqrt(bands[i].axis[0] * bands[i].axis[0] +
                                     bands[i].axis[1] * bands[i].axis[1] +
                                     bands[i].axis[2] * bands[i].axis[2]);
        if (axn > 0.0) {
            const double d = std::abs(lb[i].axis.Direction().X() * bands[i].axis[0] / axn +
                                      lb[i].axis.Direction().Y() * bands[i].axis[1] / axn +
                                      lb[i].axis.Direction().Z() * bands[i].axis[2] / axn);
            const double ang = std::acos(std::min(1.0, d));
            maxAx = std::max(maxAx, ang);
        }
        if (bands[i].face == 526) {
            face526 = lb[i].closed360 && std::abs(lb[i].phi - 2.0 * 3.141592653589793) < 1e-6;
        }
        const double oneStrip = (lb[i].N > 0) ? (lb[i].phi / lb[i].N) : 0.0;
        const double phiTrue = bands[i].extentDeg * 3.141592653589793 / 180.0;
        const bool phiOk = std::abs(lb[i].phi - phiTrue) <= std::max(oneStrip, 1e-3);
        if (!phiOk)
            std::fprintf(stderr, "  phi miss F%d got %.4f deg want %.4f\n", bands[i].face,
                         lb[i].phi * 180.0 / 3.141592653589793, bands[i].extentDeg);
        check(acc[i] == 1, (std::string("A2 accept F") + std::to_string(bands[i].face)).c_str());
        check(lb[i].cvTheta <= 2e-4,
              (std::string("A2 cvT F") + std::to_string(bands[i].face)).c_str());
        check(phiOk, (std::string("A3 phi F") + std::to_string(bands[i].face)).c_str());
        (void)rel;
    }
    check(nAcc == 15, "A2 15/15 accept");
    check(face526, "A3 F526 closed360 360deg");
    check(maxAx < 1e-6, "A4 axis || true < 1e-6 rad");
    // D5.2/D5.3: A1b mesh-recovered < 1e-4; A5 bound is τ_surf from the mesh
    // (handle-lock: max(5e-5, 4·weldTol, 1e-6·diag) = 7.820e-5), not 2.0e-5.
    const double tauSurf = std::max({5e-5, 4.0 * mv.weldTol, 1e-6 * mv.diag});
    check(maxRel < 1e-4, "A1b mesh-recovered |R-Rtrue|/Rtrue < 1e-4");
    check(maxResid < tauSurf, "A5 on-surface residual < tau_surf (mesh-computed)");
    std::fprintf(stderr, "mesh R maxRel=%.3e maxCvt=%.3e maxResid=%.3e maxAx=%.3e tau_surf=%.6e\n",
                 maxRel, maxCvt, maxResid, maxAx, tauSurf);

    std::vector<LawBand> accepted;
    for (size_t i = 0; i < lb.size(); ++i)
        if (acc[i]) accepted.push_back(lb[i]);
    const TessLawInterval cal = lawCalibrate(accepted.empty() ? lb : accepted);
    std::fprintf(stderr, "CAL d=[%.6f,%.6f] aDeg=[%.4f,%.4f) nD=%d nA=%d empty=%d\n", cal.dLo,
                 cal.dHi, cal.alphaLo * 180.0 / 3.141592653589793,
                 cal.alphaHi * 180.0 / 3.141592653589793, cal.nDLimited, cal.nAlphaLimited,
                 cal.empty ? 1 : 0);
    const double dHat = 0.5 * (cal.dLo + cal.dHi);
    const bool a6 = !cal.empty && cal.dLo < 0.0125 && 0.0125 < cal.dHi &&
                    (cal.dHi - cal.dLo) < 0.05 * std::max(dHat, 1e-12);
    check(a6, "A6 d-interval contains 0.0125 width<5%");
    const double aLoDeg = cal.alphaLo * 180.0 / 3.141592653589793;
    const double aHiDeg = cal.alphaHi * 180.0 / 3.141592653589793;
    check(cal.nAlphaLimited == 1 && std::abs(aLoDeg - 12.884) < 0.02 &&
              std::abs(aHiDeg - 15.031) < 0.02 && aHiDeg > aLoDeg + 0.5,
          "A7 alpha interval [12.884, 15.031)");

    // --- A8 RULE 4.2b ---
    {
        LawBand r30;
        r30.R = 30.0;
        r30.phi = 26.294 * 3.141592653589793 / 180.0;
        r30.N = 8;
        r30.theta.assign(8, r30.phi / 8.0);
        r30.w.assign(8, 2.0 * 30.0 * std::sin(0.5 * r30.phi / 8.0));
        // D4 §2 measured proof: ceil(26.294 / 3.2866) = 9 flips a ceil at 1.3% d-error.
        const int nPred = static_cast<int>(std::ceil(26.294 / 3.2866));
        std::fprintf(stderr, "A8 dHat=0.012345 Npred=%d trueN=8 consistent=%d\n", nPred,
                     lawNConsistent(r30, cal) ? 1 : 0);
        check(nPred == 9, "A8 point estimate predicts N=9");
        check(lawNConsistent(r30, cal), "A8 lawNConsistent(N=8) true");
    }

    // --- A11 short-arc F521 ---
    if (const AnatomyBand* f521 = findFace(bands, 521)) {
        LawBand b;
        const bool ok = lawChainAccept(mv, toLocal(mv, f521->triIds), tol, b);
        std::fprintf(stderr, "A11 F521 accept=%d N=%d R=%.9f phiDeg=%.4f w0=%.4f\n", ok ? 1 : 0,
                     b.N, b.R, b.phi * 180.0 / 3.141592653589793, b.w.empty() ? 0.0 : b.w[0]);
        check(ok && b.N == 3, "A11 F521 accept N=3");
        check(std::abs(b.R - 20.0) / 20.0 < 1e-4, "A11b F521 mesh |R-Rtrue|/Rtrue < 1e-4");
        std::vector<double> th;
        for (double d : f521->thetaDeg) th.push_back(d * 3.141592653589793 / 180.0);
        const double Rinv = inverseR(f521->wMm, th);
        std::fprintf(stderr, "A11 anatomy nW=%zu nTh=%zu Rinv=%.12f\n", f521->wMm.size(),
                     th.size(), Rinv);
        check(std::abs(Rinv - 20.0) / 20.0 < 1e-9 ||
                  (f521->RfromTheta > 0.0 && std::abs(f521->RfromTheta - 20.0) / 20.0 < 1e-9),
              "A11a inverse(w,θ) < 1e-9");
    }

    // --- A12 dihedral form NOT used ---
    if (const AnatomyBand* f515 = findFace(bands, 515)) {
        LawBand b;
        (void)lawChainAccept(mv, toLocal(mv, f515->triIds), tol, b);
        std::vector<double> thA, thD;
        for (double d : f515->thetaDeg) thA.push_back(d * 3.141592653589793 / 180.0);
        for (double d : f515->dihedralDeg) thD.push_back(d * 3.141592653589793 / 180.0);
        const double RinvA = inverseR(f515->wMm, thA);
        const double Rdih = inverseR(f515->wMm, thD);
        std::fprintf(stderr, "A12 F515 R_theta=%.9f R_anat=%.12f R_dih=%.4f (want 0.500000 vs 0.4617)\n",
                     b.R, RinvA, Rdih);
        check(std::abs(RinvA - 0.5) / 0.5 < 1e-9, "A12a anatomy inverse < 1e-9");
        check(std::abs(b.R - 0.5) / 0.5 < 1e-4, "A12b mesh |R-Rtrue|/Rtrue < 1e-4");
        check(Rdih < 0.48 && Rdih > 0.45, "A12 dihedral form ~ 0.4617");
        check(std::abs(b.R - Rdih) > 0.02, "A12 code returned θ-form not dihedral");
    }

    // --- A13 empty-intersection decline ---
    {
        LawBand a, b;
        a.R = 5.0;
        a.N = 10;
        a.phi = 80.0 * 3.141592653589793 / 180.0;
        a.theta.assign(10, a.phi / 10.0);
        a.w.assign(10, 2.0 * a.R * std::sin(0.5 * a.phi / 10.0));
        b.R = 5.0;
        b.N = 20;
        b.phi = 80.0 * 3.141592653589793 / 180.0;
        b.theta.assign(20, b.phi / 20.0);
        b.w.assign(20, 2.0 * b.R * std::sin(0.5 * b.phi / 20.0));
        bool threw = false;
        TessLawInterval empty;
        try {
            empty = lawCalibrate({a, b});
        } catch (...) {
            threw = true;
        }
        check(!threw && empty.empty, "A13 mixed-preset empty==true no throw");
    }

    // --- A9 / A10 via segment() region sets + labels ---
    RegionSet rs;
    SegmentParams sp;
    const bool segOk = segment(mv, sp, rs, nullptr);
    check(segOk, "segment() for chimera/merge");
    const stl2step::refit::Region* r5 = nullptr;
    const stl2step::refit::Region* r6 = nullptr;
    const stl2step::refit::Region* r7 = nullptr;
    const stl2step::refit::Region* r16 = nullptr;
    const stl2step::refit::Region* r17 = nullptr;
    for (const auto& r : rs.regions) {
        if (r.id == 5) r5 = &r;
        if (r.id == 6) r6 = &r;
        if (r.id == 7) r7 = &r;
        if (r.id == 16) r16 = &r;
        if (r.id == 17) r17 = &r;
    }
    if (r16 && r5) {
        LawBand chim16, chim5;
        const bool a16 = lawChainAccept(mv, r16->tris, tol, chim16);
        const bool a5 = lawChainAccept(mv, r5->tris, tol, chim5);
        std::fprintf(stderr, "A9 rid16 accept=%d cvR=%.3g N=%d R=%.4f\n", a16 ? 1 : 0, chim16.cvR,
                     chim16.N, chim16.R);
        std::fprintf(stderr, "A9 rid5  accept=%d cvR=%.3g resid=%.3g N=%d\n", a5 ? 1 : 0, chim5.cvR,
                     chim5.maxVertResid, chim5.N);
        check(!a16, "A9 rid=16 chimera rejected");
        check(!a5, "A9 rid=5 chimera rejected");
    } else {
        std::fprintf(stderr, "A9 regions missing r16=%p r5=%p (nReg=%zu)\n", (void*)r16, (void*)r5,
                     rs.regions.size());
        check(false, "A9 region ids present");
    }

    if (r16 && r7 && r17 && r5 && r6 && !labs.empty()) {
        const std::vector<int> r16r30 = byFace(labs, 538, r16->tris);
        const std::vector<int> r16r10 = byFace(labs, 542, r16->tris);
        const std::vector<int> r5f21 = byFace(labs, 536, r5->tris);
        LawBand a, b, c, d, e, f, g;
        (void)lawChainAccept(mv, r16r30, tol, a);
        (void)lawChainAccept(mv, r7->tris, tol, b);
        (void)lawChainAccept(mv, r16r10, tol, c);
        (void)lawChainAccept(mv, r17->tris, tol, d);
        (void)lawChainAccept(mv, r5f21, tol, e);
        (void)lawChainAccept(mv, r6->tris, tol, f);
        (void)lawChainAccept(mv, r7->tris, tol, g);
        const bool m1 = lawBandsMergeable(a, b, tol);
        const bool m2 = lawBandsMergeable(c, d, tol);
        const bool m3 = lawBandsMergeable(e, f, tol);
        const bool m4 = lawBandsMergeable(b, d, tol);
        auto rng = [](const LawBand& x) {
            int lo = 0, hi = -1;
            if (!x.tris.empty()) {
                lo = hi = x.tris[0];
                for (int t : x.tris) {
                    lo = std::min(lo, t);
                    hi = std::max(hi, t);
                }
            }
            std::fprintf(stderr, "    N=%d R=%.6f nTri=%zu ids=[%d,%d] theta=%zu\n", x.N, x.R,
                         x.tris.size(), lo, hi, x.theta.size());
        };
        std::fprintf(stderr, "A10 16R30:\n");
        rng(a);
        std::fprintf(stderr, "A10 rid7:\n");
        rng(b);
        std::fprintf(stderr, "A10 5F21:\n");
        rng(e);
        std::fprintf(stderr, "A10 rid6:\n");
        rng(f);
        std::fprintf(stderr, "A10 merge (16R30,7)=%d (16R10,17)=%d (5F21,6)=%d (7,17)=%d\n",
                     m1 ? 1 : 0, m2 ? 1 : 0, m3 ? 1 : 0, m4 ? 1 : 0);
        check(m1, "A10 merge rid16-R30 + rid7");
        check(m2, "A10 merge rid16-R10 + rid17");
        check(m3, "A10 merge rid5-F21 + rid6");
        check(!m4, "A10 no-merge rid7 + rid17");
    } else {
        check(false, "A10 inputs available");
    }

    if (gFail) {
        std::fprintf(stderr, "FAILED %d checks\n", gFail);
        return 1;
    }
    std::fprintf(stderr, "ALL CHECKS PASSED\n");
    return 0;
}
