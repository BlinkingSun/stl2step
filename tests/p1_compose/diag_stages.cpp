// Stage-by-stage diagnostic for composed P1 pipeline (lane p1-compose-fix).
#include "mesh_harness.hpp"
#include "mesh_harness_refit.hpp"
#include "posix_compat.hpp"
#include "refit_internal.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace stl2step;

static void printWork(const char* tag, const refit::SegmentWork& w) {
    int uncl = 0, seedT = 0, consFil = 0, plane = 0;
    for (const auto& p : w.provisionals) {
        switch (p.claim) {
            case refit::ProvClaim::Unclaimed: ++uncl; break;
            case refit::ProvClaim::InCylinderClaim: break;
            case refit::ProvClaim::ConsumedCylinder: break;
            case refit::ProvClaim::InFilletClaim: break;
            case refit::ProvClaim::ConsumedFillet: ++consFil; break;
            case refit::ProvClaim::CommittedPlane: ++plane; break;
        }
        if (p.seedTried) ++seedT;
    }
    int accCyl = 0, accPln = 0, accFil = 0;
    for (const auto& r : w.accepted) {
        if (r.origin == refit::Origin::CylGrow) ++accCyl;
        else if (r.origin == refit::Origin::FilletStrip) ++accFil;
        else ++accPln;
    }
    std::printf("%s: prov=%zu uncl=%d seedTried=%d acc(cyl/pln/fil)=%d/%d/%d rej=%zu\n",
                tag, w.provisionals.size(), uncl, seedT, accCyl, accPln, accFil,
                w.rejected.size());
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    setenv("STL2STEP_P1_DIAG", "1", 1);
    harness::HarnessMesh mesh;
    std::string err;
    if (!harness::loadMesh(argv[1], 1.0, 0.0, 0.0, mesh, err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    for (size_t ci = 0; ci < mesh.comps.size(); ++ci) {
        if (!mesh.comps[ci].clean()) continue;
        refit::MeshView mv{};
        harness::toMeshView(mesh, mesh.comps[ci], mv);
        refit::SegmentParams p;
        refit::DerivedTols tol = [&]() {
            refit::DerivedTols t;
            t.epsMesh = std::max({mv.weldTol, 1e-4 * mv.diag, 1e-3});
            t.epsPlane = std::max({t.epsMesh, mv.sewTol, 0.02});
            t.thetaPlane = p.thetaPlaneDeg * M_PI / 180.0;
            t.thetaSharp = p.thetaSharpDeg * M_PI / 180.0;
            t.thetaCylLo = p.thetaCylLoDeg * M_PI / 180.0;
            t.thetaCylHi = p.thetaCylHiDeg * M_PI / 180.0;
            t.thetaBin = p.thetaBinDeg * M_PI / 180.0;
            return t;
        }();
        refit::SegmentWork work;
        std::printf("=== %s comp %zu tris=%zu ===\n", argv[1], ci,
                    mesh.comps[ci].compTris.size());
        refit::chartsA1(mv, p, tol, work);
        printWork("A1", work);
        refit::growProvisionalA2(mv, p, tol, work);
        printWork("A2", work);
        refit::claimCylindersB1(mv, p, tol, work);
        printWork("B1", work);
        refit::claimFilletsC1(mv, p, tol, work);
        printWork("C1", work);
        refit::commitPlanesA3(mv, p, tol, work);
        printWork("A3", work);
        refit::RegionSet out;
        refit::buildTopologyD(mv, p, tol, work, out);
        std::printf("D: regions=%zu stats cyl=%d pln=%d fil=%d\n",
                    out.regions.size(), out.stats.cylinders, out.stats.planes,
                    out.stats.fillets);
    }
    return 0;
}
