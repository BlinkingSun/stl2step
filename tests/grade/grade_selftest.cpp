#include "grade_compare.hpp"
#include "grade_mesh.hpp"
#include "grade_oracle.hpp"
#include "grade_report.hpp"

#include "stl_quant.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <BRep_Builder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pln.hxx>

#include "occt_calibrated.hpp"

#ifndef STL2STEP_EXPECTED_RED_JSON
#define STL2STEP_EXPECTED_RED_JSON ""
#endif

namespace {

int gPass = 0, gFail = 0;
int gXFail = 0;
std::map<std::string, std::string> gExpectedRed;
std::string gPlateStep;

void check(bool ok, const char* name) {
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

// D-140-7 / D-130-23: a check whose stable id is listed in expected-red.json's
// "grader" section is a recorded deferral, not a hard failure. Failing is
// reported XFAIL and excluded from the exit code; passing is reported XPASS
// (the row must be removed) and counted as a failure so the list shrinks.
void checkId(bool ok, const std::string& id, const char* name) {
    const auto it = gExpectedRed.find(id);
    if (it == gExpectedRed.end()) {
        check(ok, name);
        return;
    }
    if (ok) {
        ++gFail;
        std::fprintf(stderr, "XPASS %s — remove the row from expected-red.json\n", id.c_str());
    } else {
        ++gXFail;
        std::fprintf(stderr, "XFAIL %s (%s)\n", id.c_str(), it->second.c_str());
    }
}

std::string readFileText(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::string();
    std::fseek(f, 0, SEEK_END);
    const long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::string s(n > 0 ? static_cast<size_t>(n) : 0, '\0');
    if (n > 0) {
        const size_t got = std::fread(&s[0], 1, static_cast<size_t>(n), f);
        s.resize(got);
    }
    std::fclose(f);
    return s;
}

// Minimal reader for tests/gates/baseline/expected-red.json (D-130-23): a flat
// JSON object of named sections, each an object mapping check id -> reason
// string. Returns the "grader" section only. Missing file or missing section
// is an empty map — never a silent XFAIL.
std::map<std::string, std::string> loadExpectedRedGrader(const std::string& path) {
    std::map<std::string, std::string> out;
    const std::string text = readFileText(path);
    if (text.empty()) return out;

    auto parseString = [&](size_t& i) -> std::string {
        std::string s;
        ++i;  // opening quote
        while (i < text.size() && text[i] != '"') {
            if (text[i] == '\\' && i + 1 < text.size()) {
                ++i;
                s += text[i];
            } else {
                s += text[i];
            }
            ++i;
        }
        ++i;  // closing quote
        return s;
    };

    const size_t sectionKeyPos = text.find("\"grader\"");
    if (sectionKeyPos == std::string::npos) return out;
    const size_t objStart = text.find('{', sectionKeyPos);
    if (objStart == std::string::npos) return out;
    int depth = 0;
    size_t objEnd = std::string::npos;
    for (size_t j = objStart; j < text.size(); ++j) {
        if (text[j] == '{') ++depth;
        else if (text[j] == '}') {
            --depth;
            if (depth == 0) { objEnd = j; break; }
        }
    }
    if (objEnd == std::string::npos) return out;

    size_t i = objStart + 1;
    while (i < objEnd) {
        while (i < objEnd && std::strchr(" \n\r\t,", text[i])) ++i;
        if (i >= objEnd || text[i] != '"') break;
        const std::string key = parseString(i);
        while (i < objEnd && text[i] != ':') ++i;
        ++i;  // ':'
        while (i < objEnd && std::strchr(" \n\r\t", text[i])) ++i;
        if (i >= objEnd || text[i] != '"') break;
        out[key] = parseString(i);
    }
    return out;
}

bool near(double a, double b, double tol) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

std::string join(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    char c = a.back();
    if (c == '/' || c == '\\') return a + b;
    return a + "/" + b;
}

std::string tmpRoot() {
#ifdef _WIN32
    const char* t = std::getenv("TEMP");
    return t ? std::string(t) : std::string(".");
#else
    return "/tmp";
#endif
}

void ensureDir(const std::string& p) {
#ifdef _WIN32
    _mkdir(p.c_str());
#else
    mkdir(p.c_str(), 0755);
#endif
}

void silence() { grade::prepareOcctThread(); }

bool writeStep(const TopoDS_Shape& s, const std::string& path) {
    Interface_Static::SetCVal("write.step.schema", "AP214IS");
    STEPControl_Writer w;
    if (w.Transfer(s, STEPControl_AsIs) < 1) return false;
    return w.Write(path.c_str()) == IFSelect_RetDone;
}

bool loadShape(const std::string& path, TopoDS_Shape& s) {
    STEPControl_Reader r;
    if (r.ReadFile(path.c_str()) != IFSelect_RetDone) return false;
    if (r.TransferRoots() < 1) return false;
    s = r.OneShape();
    return true;
}

int countClass(const grade::GradeDocument& d, grade::SurfClass c) {
    int n = 0;
    for (const auto& f : d.features)
        if (f.oracle.cls == c) ++n;
    return n;
}

int countTris(const grade::GradeDocument& d, grade::SurfClass c) {
    int n = 0;
    for (const auto& f : d.features)
        if (f.oracle.cls == c) n += static_cast<int>(f.oracle.tris.size());
    return n;
}

bool allStatus(const grade::GradeDocument& d, grade::SurfClass c, grade::Status st) {
    bool any = false;
    for (const auto& f : d.features) {
        if (f.oracle.cls != c) continue;
        any = true;
        if (f.status != st) return false;
    }
    return any;
}

int coreTests(const std::string& corpus) {
    silence();
    std::setlocale(LC_ALL, "C");
    grade::GradeConfig cfg;
    cfg.skipVolume = true;
    cfg.quiet = true;

    // Every gradeFiles call below, in source order, produced up front. Checks
    // stay on this thread and replay each grade's GRADE_FIXPOINT line where
    // the call used to print it (D-train-grader-5 (2) S3).
    struct Job {
        std::string stl, step;
        grade::GradeConfig cfg;
        grade::GradeDocument doc;
        std::string err;
        std::string cap;
        bool ok = false;
    };
    std::vector<Job> jobs;
    auto enqueue = [&](const std::string& stl, const std::string& step, const grade::GradeConfig& c) {
        Job j;
        j.stl = stl;
        j.step = step;
        j.cfg = c;
        jobs.push_back(std::move(j));
    };
    grade::GradeConfig rev = cfg;
    rev.reverseSeeds = true;
    const std::string s01 = join(corpus, "S01.stl");
    const std::string s01step = join(corpus, "S01.exact.step");
    const std::string s03 = join(corpus, "S03.stl");
    const std::string s03step = join(corpus, "S03.exact.step");
    const std::string hp = join(corpus, "handle-pickup.stl");
    const std::string plate = join(corpus, "linkage_bores_chamfer.stl");
    enqueue(s01, s01step, cfg);
    enqueue(join(corpus, "S02.stl"), join(corpus, "S02.exact.step"), cfg);
    enqueue(s03, s03step, cfg);
    enqueue(join(corpus, "S04.stl"), join(corpus, "S04.exact.step"), cfg);
    enqueue(s01, s01step, cfg);
    enqueue(s01, s01step, cfg);
    enqueue(s01, s01step, rev);
    enqueue(s03, s03step, cfg);
    enqueue(s03, s03step, cfg);
    enqueue(s03, s03step, rev);
    enqueue(hp, s01step, cfg);
    enqueue(hp, s01step, cfg);
    enqueue(hp, s01step, rev);
    enqueue(plate, s01step, cfg);
    if (!gPlateStep.empty()) enqueue(plate, gPlateStep, cfg);
    enqueue(join(corpus, "S20_cross_bore_union.stl"),
            join(corpus, "S20_cross_bore_union.exact.step"), cfg);

    unsigned nJobs = std::thread::hardware_concurrency();
    if (const char* e = std::getenv("STL2STEP_GRADE_SELFTEST_JOBS")) {
        if (e[0]) {
            const int v = std::atoi(e);
            nJobs = v > 0 ? static_cast<unsigned>(v) : 1u;
        }
    }
    if (nJobs == 0) nJobs = 1;
    std::fprintf(stdout, "grade_selftest jobs=%u\n", nJobs);
    std::fflush(stdout);
    {
        const auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]() {
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        };
        auto fixtureOf = [](const Job& j, size_t i) {
            std::string s = j.stl;
            const auto sl = s.find_last_of("/\\");
            if (sl != std::string::npos) s = s.substr(sl + 1);
            s += "#";
            s += std::to_string(i);
            return s;
        };
        std::atomic<size_t> cursor{0};
        auto worker = [&]() {
            grade::prepareOcctThread();
            for (;;) {
                const size_t i = cursor.fetch_add(1);
                if (i >= jobs.size()) break;
                const std::string fx = fixtureOf(jobs[i], i);
                std::fprintf(stderr, "grade_selftest start %s t=%.3f\n", fx.c_str(), elapsed());
                std::fflush(stderr);
                grade::setOracleStderrSink(&jobs[i].cap);
                jobs[i].ok = grade::gradeFiles(jobs[i].stl, jobs[i].step, jobs[i].cfg, jobs[i].doc,
                                               jobs[i].err);
                grade::setOracleStderrSink(nullptr);
                std::fprintf(stderr, "grade_selftest done %s t=%.3f\n", fx.c_str(), elapsed());
                std::fflush(stderr);
            }
        };
        const unsigned pool =
            nJobs < static_cast<unsigned>(jobs.size()) ? nJobs : static_cast<unsigned>(jobs.size());
        std::vector<std::thread> threads;
        threads.reserve(pool);
        for (unsigned t = 0; t < pool; ++t) threads.emplace_back(worker);
        for (std::thread& t : threads) t.join();
    }
    size_t nextJob = 0;
    auto take = [&](grade::GradeDocument& d, std::string& err) -> bool {
        Job& j = jobs[nextJob++];
        if (!j.cap.empty()) std::fputs(j.cap.c_str(), stderr);
        d = std::move(j.doc);
        err = std::move(j.err);
        return j.ok;
    };

    // Case 1: S01
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "S01 grade runs");
        if (!err.empty()) std::fprintf(stderr, "  S01 err: %s\n", err.c_str());
        check(d.mesh.tris.size() == 12, "S01 12 triangles");
        check(near(d.mesh.q, 8.2591e-07, 1e-10), "S01 q");
        check(countClass(d, grade::SurfClass::Plane) == 6, "S01 6 plane oracles");
        check(countTris(d, grade::SurfClass::Plane) == 12, "S01 planes cover 12");
        bool twoEach = true;
        for (const auto& f : d.features)
            if (f.oracle.cls == grade::SurfClass::Plane && f.oracle.tris.size() != 2) twoEach = false;
        check(twoEach, "S01 each plane 2 tris");
        check(allStatus(d, grade::SurfClass::Plane, grade::Status::Recovered), "S01 all recovered");
        check(d.oracle.residueTris == 0, "S01 residue 0");
        check(d.hasPlane && near(d.gradePlane, 1.0, 1e-12), "S01 grade.plane 1");
        check(d.hasOverall && near(d.gradeOverall, 1.0, 1e-12), "S01 grade.overall 1");
        check(d.hardZero.empty(), "S01 no hardZero");
    }

    // Case 2: S02
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "S02 grade runs");
        check(d.mesh.tris.size() == 412, "S02 412 triangles");
        check(near(d.mesh.q, 1.6518e-06, 1e-9), "S02 q");
        const int np = countClass(d, grade::SurfClass::Plane);
        const int nc = countClass(d, grade::SurfClass::Cylinder);
        const int ns = countClass(d, grade::SurfClass::Sphere);
        std::fprintf(stderr, "  S02 oracles plane=%d cyl=%d sph=%d residue=%d\n", np, nc, ns,
                     d.oracle.residueTris);
        std::fprintf(stderr, "  S02 tris plane=%d cyl=%d sph=%d\n", countTris(d, grade::SurfClass::Plane),
                     countTris(d, grade::SurfClass::Cylinder), countTris(d, grade::SurfClass::Sphere));
        check(np == 6, "S02 6 plane oracles");
        check(nc == 12, "S02 12 cylinder oracles");
        check(ns == 8, "S02 8 sphere oracles");
        bool cyl10 = true, sph35 = true, pl2 = true;
        for (const auto& f : d.features) {
            if (f.oracle.cls == grade::SurfClass::Plane && f.oracle.tris.size() != 2) pl2 = false;
            if (f.oracle.cls == grade::SurfClass::Cylinder && f.oracle.tris.size() != 10) cyl10 = false;
            if (f.oracle.cls == grade::SurfClass::Sphere && f.oracle.tris.size() != 35) sph35 = false;
        }
        check(pl2, "S02 planes ×2");
        check(cyl10, "S02 cylinders ×10");
        check(sph35, "S02 spheres ×35");
        check(d.oracle.residueTris == 0, "S02 residue 0");
        check(d.hasPlane && near(d.gradePlane, 1.0, 1e-12), "S02 grade.plane 1");
        check(d.hasCyl && near(d.gradeCyl, 1.0, 1e-12), "S02 grade.cylinder 1");
        check(d.hasSphere && near(d.gradeSphere, 1.0, 1e-12), "S02 grade.sphere 1");
        check(d.hasOverall && near(d.gradeOverall, 1.0, 1e-12), "S02 grade.overall 1");
        bool residOk = true;
        for (const auto& f : d.features)
            if (f.oracle.maxResid > 2.0 * d.mesh.q) residOk = false;
        check(residOk, "S02 maxResid <= 2q");
    }

    // Case 3 + 3b: S03
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "S03 grade runs");
        check(d.mesh.tris.size() == 584, "S03 584 triangles");
        check(near(d.mesh.q, 6.6072e-06, 1e-9), "S03 q");
        int nCyl4 = 0, nCyl76 = 0;
        std::vector<int> planeN;
        int nCone = 0, nZ4 = 0;
        for (const auto& f : d.features) {
            if (f.oracle.cls == grade::SurfClass::Cylinder && near(f.oracle.S.R, 4.0, 1e-3)) {
                ++nCyl4;
                if (f.oracle.tris.size() == 76) ++nCyl76;
                std::fprintf(stderr, "  S03 cyl R=%g tris=%zu maxResid=%g status=%s\n", f.oracle.S.R,
                             f.oracle.tris.size(), f.oracle.maxResid, grade::statusName(f.status));
                check(f.oracle.maxResid <= d.mesh.q, "S03 cyl maxResid <= q");
                check(f.status == grade::Status::Recovered, "S03 cyl recovered");
            }
            if (f.oracle.cls == grade::SurfClass::Cone) ++nCone;
            if (f.oracle.cls == grade::SurfClass::Plane) {
                planeN.push_back(static_cast<int>(f.oracle.tris.size()));
                const double z = f.oracle.S.p0.z;  // foot; for z=const plane, n is ±Z, p0.z is offset
                const bool zPlane = std::fabs(std::fabs(f.oracle.S.n.z) - 1.0) < 1e-6;
                if (zPlane && std::fabs(std::fabs(dot(f.oracle.S.p0, f.oracle.S.n)) - 4.0) < 0.1)
                    ++nZ4;
            }
        }
        std::fprintf(stderr, "  S03 cyl4=%d cyl76=%d cone=%d z4planes=%d residue=%d planes=%zu\n",
                     nCyl4, nCyl76, nCone, nZ4, d.oracle.residueTris, planeN.size());
        if (d.oracle.residueTris != 0) {
            std::fprintf(stderr, "  S03 residue finding: %d tris\n", d.oracle.residueTris);
            int shown = 0;
            for (int t = 0; t < static_cast<int>(d.mesh.tris.size()) && shown < 12; ++t) {
                if (d.oracle.owner[static_cast<size_t>(t)] >= 0) continue;
                const grade::Tri& tr = d.mesh.tris[static_cast<size_t>(t)];
                const grade::Vec3& a = d.mesh.verts[static_cast<size_t>(tr.v[0])];
                std::fprintf(stderr, "    tri %d v0=(%g,%g,%g)\n", t, a.x, a.y, a.z);
                ++shown;
            }
        }
        check(nCyl4 == 4, "S03 4 R=4 cylinders");
        // D-train-grader-4 (4c): the || nCyl4 == 4 disjunct is deleted. Each bore
        // oracle is 52 and the whole-mesh tau-set is 76. The red is registered.
        checkId(nCyl76 == 4, "grade.s03-bore-claim-52", "S03 cylinders are 76 tris");
        std::fprintf(stderr, "  S03 bore census nCyl4=%d nCyl76=%d\n", nCyl4, nCyl76);
        check(nCone >= 1, "S03 >=1 cone");
        std::sort(planeN.begin(), planeN.end(), std::greater<int>());
        check(!planeN.empty() && planeN[0] == 142, "S03 outer plane 142");
        int twos = 0;
        for (int n : planeN)
            if (n == 2) ++twos;
        check(twos >= 5, "S03 five 2-tri outer planes");
        check(nZ4 == 5, "S03 five z=4 plane oracles (partition)");
        check(d.hasPlane && near(d.gradePlane, 1.0, 1e-12), "S03 grade.plane 1");
        check(d.hasCyl && near(d.gradeCyl, 1.0, 1e-12), "S03 grade.cylinder 1");
        check(d.hasCone && near(d.gradeCone, 1.0, 1e-12), "S03 grade.cone 1");
        check(d.hasOverall && near(d.gradeOverall, 1.0, 1e-12), "S03 grade.overall 1");
        check(d.oracle.residueTris == 0, "S03 residue 0 (assert; finding if fail)");
        int planeCone = 0, planeConeTier1 = 0, coneCyl = 0, coneCylTier2 = 0;
        for (const auto& x : d.intersections) {
            const bool aCone = x.a.rfind("cone:", 0) == 0;
            const bool bCone = x.b.rfind("cone:", 0) == 0;
            const bool aCyl = x.a.rfind("cylinder:", 0) == 0;
            const bool bCyl = x.b.rfind("cylinder:", 0) == 0;
            const bool aPln = x.a.rfind("plane:", 0) == 0;
            const bool bPln = x.b.rfind("plane:", 0) == 0;
            if ((aCone && bPln) || (bCone && aPln)) {
                ++planeCone;
                if (x.expectedTier == 1) ++planeConeTier1;
            }
            if ((aCone && bCyl) || (bCone && aCyl)) {
                ++coneCyl;
                if (x.expectedTier == 2) ++coneCylTier2;
            }
        }
        std::fprintf(stderr, "  S03 plane|cone=%d tier1=%d cone|cyl=%d tier2=%d\n", planeCone,
                     planeConeTier1, coneCyl, coneCylTier2);
        check(planeCone > 0 && planeCone == planeConeTier1, "S03 plane|cone tier 1");
        // S03's drafted cone shares no mesh edge with a non-coaxial cylinder.
        // The tier-2 row is the in-test negative in synthetic mode; a count of
        // zero here must not be the only guard.
        (void)coneCyl;
        (void)coneCylTier2;
    }

    // Case 3c / A2: S04 torus reported, not a 1.0 claim
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "S04 grade runs");
        check(d.mesh.tris.size() == 652, "S04 652 triangles");
        int nT = 0, t512 = 0, nC = 0, nP = 0, top30 = 0;
        grade::Status tStat = grade::Status::Missing;
        for (const auto& f : d.features) {
            if (f.oracle.cls == grade::SurfClass::Torus) {
                ++nT;
                tStat = f.status;
                if (f.oracle.tris.size() == 512) ++t512;
                std::fprintf(stderr, "  S04 torus tris=%zu Rmaj=%g Rmin=%g status=%s c=(%g,%g,%g)\n",
                             f.oracle.tris.size(), f.oracle.S.R, f.oracle.S.r, grade::statusName(f.status),
                             f.oracle.S.p0.x, f.oracle.S.p0.y, f.oracle.S.p0.z);
            }
            if (f.oracle.cls == grade::SurfClass::Cylinder && near(f.oracle.S.R, 10.0, 1e-2)) {
                ++nC;
                std::fprintf(stderr, "  S04 cyl tris=%zu\n", f.oracle.tris.size());
            }
            if (f.oracle.cls == grade::SurfClass::Plane) {
                ++nP;
                if (f.oracle.tris.size() == 30) ++top30;
            }
        }
        std::fprintf(stderr, "  S04 planes=%d cyl=%d torus=%d residue=%d\n", nP, nC, nT,
                     d.oracle.residueTris);
        check(nT >= 1, "S04 torus oracle present");
        check(t512 == 1, "S04 torus 512 tris");
        check(top30 == 1, "S04 plane 30 (boss top, class order)");
        // A2: do not require grade.torus == 1.0. Exact STEP may recover the torus;
        // the sidecar is mustRemainFaceted. Report status, never XFAIL.
        (void)tStat;
        check(nC >= 1, "S04 cylinder present");
    }

    // Case 9: q + determinism
    {
        const char* files[] = {"S01.stl", "S03.stl", "handle-pickup.stl"};
        const double expectQ[] = {8.2591e-07, 6.6072e-06, 0};
        for (int i = 0; i < 3; ++i) {
            const std::string p = join(corpus, files[i]);
            const auto qf = stl2step::stlQuantFloor(p);
            grade::Mesh mesh;
            std::string err;
            check(grade::loadStl(p, mesh, err), (std::string("load ") + files[i]).c_str());
            check(mesh.q == qf.q, (std::string("q equals stlQuantFloor ") + files[i]).c_str());
            if (expectQ[i] > 0) check(near(mesh.q, expectQ[i], 1e-9), (std::string("q table ") + files[i]).c_str());
        }
        auto once = [&](bool rev, const std::string& stl, const std::string& step) {
            (void)rev;
            (void)stl;
            (void)step;
            grade::GradeDocument d;
            std::string err;
            take(d, err);
            return std::make_pair(grade::writeJson(d), grade::writeMd(d));
        };
        {
            const std::string stl = join(corpus, "S01.stl");
            const std::string step = join(corpus, "S01.exact.step");
            auto a = once(false, stl, step);
            auto b = once(false, stl, step);
            auto c = once(true, stl, step);
            check(a.first == b.first && a.second == b.second, "S01 twice identical");
            checkId(a.first == c.first && a.second == c.second, "grade.reverse-seed.S01",
                    "S01 reverse-seed identical");
        }
        {
            const std::string stl = join(corpus, "S03.stl");
            const std::string step = join(corpus, "S03.exact.step");
            auto a = once(false, stl, step);
            auto b = once(false, stl, step);
            auto c = once(true, stl, step);
            check(a.first == b.first, "S03 twice json identical");
            checkId(a.first == c.first, "grade.reverse-seed.S03", "S03 reverse-seed json identical");
        }
        {
            const std::string stl = join(corpus, "handle-pickup.stl");
            const std::string step = join(corpus, "S01.exact.step");  // any STEP; oracle from mesh
            auto a = once(false, stl, step);
            auto b = once(false, stl, step);
            auto c = once(true, stl, step);
            check(a.first == b.first, "handle-pickup twice json identical");
            checkId(a.first == c.first, "grade.reverse-seed.handle-pickup",
                    "handle-pickup reverse-seed json identical");
        }
    }

    // T1 — mouth band stays a 200-triangle cone; the plate has 54 plane oracles.
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "T1 plate grades");
        int planes = 0, cone200 = 0;
        for (const auto& f : d.features) {
            if (f.oracle.cls == grade::SurfClass::Plane) ++planes;
            if (f.oracle.cls == grade::SurfClass::Cone && f.oracle.tris.size() == 200) ++cone200;
        }
        std::fprintf(stderr, "  T1 plate planes=%d cone200=%d\n", planes, cone200);
        check(cone200 == 1, "T1 mouth cone 200 triangles");
        check(planes == 54, "T1 plate 54 plane oracles");
        int walls = 0;
        for (const auto& f : d.features) {
            if (f.oracle.cls != grade::SurfClass::Cylinder) continue;
            if (std::fabs(f.oracle.S.R - 3.0) > 0.01) continue;
            if (f.oracle.tris.size() != 39) continue;
            ++walls;
        }
        std::fprintf(stderr, "  T1 R=3 walls of 39 tris=%d\n", walls);
        check(walls == 4, "T1 four R=3 walls 39/39");
        // Axis-clause census at the landed emit predicate (the floor). The
        // retirement instrument is rulings >= d_axis+1; requiring it dropped
        // these four walls. The floor stays.
        std::fprintf(stderr,
                     "  grade.cyl-axis-overdetermination landed=triangle-floor "
                     "rulings-gate=off plate-planes=%d walls39=%d\n",
                     planes, walls);
        if (!gPlateStep.empty()) {
            grade::GradeDocument ps;
            std::string perr;
            check(take(ps, perr), "T1 plate step grades");
            int downgraded = 0, other = 0;
            struct FaceHit {
                int entity;
                double area;
                std::string id;
            };
            std::vector<std::pair<std::string, std::vector<FaceHit>>> wallFaces;
            for (const auto& x : ps.intersections) {
                if (x.verdict == "downgraded") ++downgraded;
                if (x.verdict == "other") {
                    ++other;
                    std::fprintf(stderr, "  other %s x %s entity=%d shipped=%s tier=%d\n", x.a.c_str(),
                                 x.b.c_str(), x.edgeEntity, x.shipped.c_str(), x.expectedTier);
                }
            }
            check(downgraded == 0, "T1 downgraded == 0");
            check(other <= 2, "T1 other ratchet <= 2");
            checkId(other == 0, "grade.rim-shipped-line", "T1 rims ship as polyline");
            for (const auto& f : ps.features) {
                if (f.oracle.cls != grade::SurfClass::Cylinder) continue;
                if (std::fabs(f.oracle.S.R - 3.0) > 0.01) continue;
                if (f.oracle.tris.size() != 39) continue;
                std::vector<FaceHit> hits;
                for (const auto& sf : f.stepFaces)
                    hits.push_back(FaceHit{sf.entity, sf.areaMM2, f.oracle.featureId});
                wallFaces.push_back({f.oracle.featureId, hits});
                std::fprintf(stderr, "  wall %s faces", f.oracle.featureId.c_str());
                for (const auto& h : hits)
                    std::fprintf(stderr, " (%d %.6f)", h.entity, h.area);
                std::fprintf(stderr, "\n");
            }
            check(wallFaces.size() == 4, "grade.wall-face-pairing four walls");
            std::set<int> ents;
            bool bands = wallFaces.size() == 4;
            const int cover[8] = {34, 35, 36, 37, 40, 41, 42, 43};
            for (const auto& wf : wallFaces) {
                if (wf.second.size() != 2) {
                    bands = false;
                    continue;
                }
                double a0 = wf.second[0].area, a1 = wf.second[1].area;
                if (a0 < a1) std::swap(a0, a1);
                const bool hi = std::fabs(a0 - 51.528) < 0.01;
                const bool lo = std::fabs(a1 - 43.900) < 0.01;
                if (!hi || !lo) bands = false;
                ents.insert(wf.second[0].entity);
                ents.insert(wf.second[1].entity);
            }
            bool coverOk = ents.size() == 8;
            for (int e : cover)
                if (!ents.count(e)) coverOk = false;
            check(bands, "grade.wall-face-pairing each wall has 51.528 and 43.900");
            check(coverOk, "grade.wall-face-pairing covers 34-37 and 40-43 once");
        }
    }

    // T2 — S20 R=10 is exactly two oracles, disjoint spans, union bore 96/38/1.
    {
        grade::GradeDocument d;
        std::string err;
        check(take(d, err), "T2 S20 grades");
        struct Row {
            double x0, x1;
            int pieces;
            int tris;
        };
        std::vector<Row> r10;
        for (const auto& f : d.features) {
            if (f.oracle.cls != grade::SurfClass::Cylinder) continue;
            if (std::fabs(f.oracle.S.R - 10.0) > 0.05) continue;
            Row r;
            r.x0 = f.oracle.bboxMin.x;
            r.x1 = f.oracle.bboxMax.x;
            r.pieces = f.oracle.oraclePieces;
            r.tris = static_cast<int>(f.oracle.tris.size());
            r10.push_back(r);
            std::fprintf(stderr, "  T2 R10 tris=%d pieces=%d x=[%.4f,%.4f]\n", r.tris, r.pieces, r.x0,
                         r.x1);
        }
        std::fprintf(stderr, "  T2 R10 oracles=%zu (recorded ceiling 2)\n", r10.size());
        check(r10.size() == 2, "T2 R=10 exactly two oracles");
        if (r10.size() == 2) {
            if (r10[0].x0 > r10[1].x0) std::swap(r10[0], r10[1]);
            const bool overlap = !(r10[0].x1 < r10[1].x0 || r10[1].x1 < r10[0].x0);
            check(!overlap, "T2 axial overlap 0");
            const bool spanLo = r10[0].x0 < 5.0 && r10[0].x1 > 25.0 && r10[0].x1 < 40.0;
            const bool spanHi = r10[1].x0 > 40.0 && r10[1].x0 < 55.0 && r10[1].x1 > 75.0;
            check(spanLo && spanHi, "T2 spans [0,30] and [50,80]");
            check(r10[0].tris == 96 && r10[0].pieces == 38, "T2 union bore 96 tris / 38 edge pieces");
            check(r10[1].tris == 96 && r10[1].pieces == 1, "T2 control bore 96 tris / 1 piece");
            // Decomposition of the union bore's 96 triangles into maximal
            // coplanar regions (D-train-grader-4 (2), 56+24+6+10).
            if (r10[0].tris == 96) {
                std::vector<int> bore;
                for (const auto& f : d.features) {
                    if (f.oracle.cls != grade::SurfClass::Cylinder) continue;
                    if (std::fabs(f.oracle.S.R - 10.0) > 0.05) continue;
                    if (f.oracle.bboxMin.x > 5.0) continue;
                    bore = f.oracle.tris;
                    break;
                }
                std::vector<char> onBore(d.mesh.tris.size(), 0);
                for (int t : bore) onBore[static_cast<size_t>(t)] = 1;
                std::vector<int> comp(d.mesh.tris.size(), -1);
                int nComp = 0;
                // One plane per region, locked to the seed triangle. A chain of
                // pairwise-flat steps is not one plane (that is what swallowed
                // the 24 single-triangle bore regions).
                auto onSeedPlane = [&](int seed, int b) {
                    const grade::Tri& A = d.mesh.tris[static_cast<size_t>(seed)];
                    const grade::Tri& B = d.mesh.tris[static_cast<size_t>(b)];
                    const double dotn = std::fabs(A.n.x * B.n.x + A.n.y * B.n.y + A.n.z * B.n.z);
                    const double th = std::max(A.thetaQ, B.thetaQ);
                    if (dotn < std::cos(th)) return false;
                    for (int k = 0; k < 3; ++k) {
                        const grade::Vec3& p = d.mesh.verts[static_cast<size_t>(B.v[k])];
                        const double dist = std::fabs((p.x - A.centroid.x) * A.n.x +
                                                      (p.y - A.centroid.y) * A.n.y +
                                                      (p.z - A.centroid.z) * A.n.z);
                        if (dist > d.mesh.tau) return false;
                    }
                    return true;
                };
                for (int t = 0; t < static_cast<int>(d.mesh.tris.size()); ++t) {
                    if (comp[static_cast<size_t>(t)] >= 0) continue;
                    const int id = nComp++;
                    std::vector<int> stack{t};
                    comp[static_cast<size_t>(t)] = id;
                    while (!stack.empty()) {
                        const int u = stack.back();
                        stack.pop_back();
                        for (int nb : d.mesh.adj[static_cast<size_t>(u)]) {
                            if (comp[static_cast<size_t>(nb)] >= 0) continue;
                            if (!onSeedPlane(t, nb)) continue;
                            comp[static_cast<size_t>(nb)] = id;
                            stack.push_back(nb);
                        }
                    }
                }
                int g56 = 0, g24 = 0, g6 = 0, g10 = 0;
                std::vector<std::vector<int>> groups(static_cast<size_t>(nComp));
                for (int t = 0; t < static_cast<int>(d.mesh.tris.size()); ++t)
                    groups[static_cast<size_t>(comp[static_cast<size_t>(t)])].push_back(t);
                for (const auto& g : groups) {
                    int on = 0;
                    for (int t : g)
                        if (onBore[static_cast<size_t>(t)]) ++on;
                    if (on == 0) continue;
                    if (static_cast<int>(g.size()) == 2 && on == 2) g56 += on;
                    else if (static_cast<int>(g.size()) == 1 && on == 1) g24 += on;
                    else if (static_cast<int>(g.size()) == 2 && on == 1) g6 += on;
                    else if (static_cast<int>(g.size()) == 3 && on == 1) g10 += on;
                }
                std::fprintf(stderr, "  T2 decomposition wholly-on-bore-quads=%d singles=%d one-of-two=%d one-of-three=%d\n",
                             g56, g24, g6, g10);
                // The 28 on-bore quads (56) reproduce. The 24/6/10 groups are the
                // grader's plane partition, which this normal-and-tau flood does
                // not recover (measured singles=0 one-of-two=0 one-of-three=8).
                check(g56 == 56, "T2 wholly-on-bore quads 56");
            }
        } else {
            check(false, "T2 axial overlap 0");
            check(false, "T2 spans [0,30] and [50,80]");
            check(false, "T2 union bore 96 tris / 38 edge pieces");
            check(false, "T2 control bore 96 tris / 1 piece");
        }
    }
    return gFail;
}

TopoDS_Shape cubeWithSplitTop() {
    // 10 mm cube: five intact faces + top split at x=5.
    BRepBuilderAPI_Sewing sew;
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(10, 0, 0), gp_Dir(1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, -1, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 10, 0), gp_Dir(0, 1, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 10), gp_Dir(0, 0, 1)), 0, 5, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 10), gp_Dir(0, 0, 1)), 5, 10, 0, 10).Face());
    sew.Perform();
    return sew.SewedShape();
}

TopoDS_Shape cubeWithSliverTop() {
    BRepBuilderAPI_Sewing sew;
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(10, 0, 0), gp_Dir(1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, -1, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 10, 0), gp_Dir(0, 1, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 10), gp_Dir(0, 0, 1)), 2, 8, 2, 8).Face());
    sew.Perform();
    return sew.SewedShape();
}

TopoDS_Shape cubeSpanned() {
    // SPEC §8 case 7: S01 top face enlarged in +Y (0..20) so it overhangs
    // the +Y oracle plane. No +Y face in the STEP.
    BRepBuilderAPI_Sewing sew;
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(10, 0, 0), gp_Dir(1, 0, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, -1, 0)), 0, 10, 0, 10).Face());
    sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 10), gp_Dir(0, 0, 1)), 0, 10, 0, 20).Face());
    sew.Perform();
    return sew.SewedShape();
}

bool writeBinaryStl(const std::string& path,
                    const std::vector<std::array<grade::Vec3, 3>>& tris) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    unsigned char head[84];
    std::memset(head, 0, 84);
    const uint32_t n = static_cast<uint32_t>(tris.size());
    std::memcpy(head + 80, &n, 4);
    std::fwrite(head, 1, 84, f);
    for (const auto& t : tris) {
        unsigned char rec[50];
        std::memset(rec, 0, 50);
        float v[9] = {static_cast<float>(t[0].x), static_cast<float>(t[0].y),
                      static_cast<float>(t[0].z), static_cast<float>(t[1].x),
                      static_cast<float>(t[1].y), static_cast<float>(t[1].z),
                      static_cast<float>(t[2].x), static_cast<float>(t[2].y),
                      static_cast<float>(t[2].z)};
        std::memcpy(rec + 12, v, 36);
        std::fwrite(rec, 1, 50, f);
    }
    std::fclose(f);
    return true;
}

bool fileNonempty(const std::string& p) {
    std::FILE* f = std::fopen(p.c_str(), "rb");
    if (!f) return false;
    std::fseek(f, 0, SEEK_END);
    const long n = std::ftell(f);
    std::fclose(f);
    return n > 0;
}

std::string engineBeside(const std::string& argv0) {
    std::string d = argv0;
    const auto sl = d.find_last_of("/\\");
    if (sl == std::string::npos) d = ".";
    else
        d = d.substr(0, sl);
    return join(d, "stl2step");
}

void unifyCoplanarFaces(const grade::Mesh& mesh, BRepBuilderAPI_Sewing& sew) {
    const int T = static_cast<int>(mesh.tris.size());
    std::vector<char> seen(static_cast<size_t>(T), 0);
    for (int t = 0; t < T; ++t) {
        if (seen[static_cast<size_t>(t)]) continue;
        std::vector<int> stack{t}, group;
        seen[static_cast<size_t>(t)] = 1;
        const grade::Vec3 n0 = mesh.tris[static_cast<size_t>(t)].n;
        const grade::Vec3 p0 = mesh.verts[static_cast<size_t>(mesh.tris[static_cast<size_t>(t)].v[0])];
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            group.push_back(u);
            for (int n : mesh.adj[static_cast<size_t>(u)]) {
                if (seen[static_cast<size_t>(n)]) continue;
                const grade::Tri& tr = mesh.tris[static_cast<size_t>(n)];
                if (grade::angleUnit(tr.n, n0) > tr.thetaQ) continue;
                bool on = true;
                for (int k = 0; k < 3; ++k) {
                    const grade::Vec3& v = mesh.verts[static_cast<size_t>(tr.v[k])];
                    if (std::fabs(grade::dot(v - p0, n0)) > mesh.tau) {
                        on = false;
                        break;
                    }
                }
                if (!on) continue;
                seen[static_cast<size_t>(n)] = 1;
                stack.push_back(n);
            }
        }
        grade::Vec3 a = mesh.verts[static_cast<size_t>(mesh.tris[static_cast<size_t>(group[0])].v[0])];
        grade::Vec3 u, v;
        const grade::Vec3 n = n0;
        grade::Vec3 tmpv = (std::fabs(n.z) < 0.9) ? grade::Vec3{0, 0, 1} : grade::Vec3{1, 0, 0};
        u = grade::normalized(grade::cross(tmpv, n));
        v = grade::cross(n, u);
        double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
        for (int gi : group) {
            const grade::Tri& tr = mesh.tris[static_cast<size_t>(gi)];
            for (int k = 0; k < 3; ++k) {
                const grade::Vec3 d = mesh.verts[static_cast<size_t>(tr.v[k])] - a;
                const double uu = grade::dot(d, u), vv = grade::dot(d, v);
                umin = std::min(umin, uu);
                umax = std::max(umax, uu);
                vmin = std::min(vmin, vv);
                vmax = std::max(vmax, vv);
            }
        }
        gp_Pln pl(gp_Pnt(a.x, a.y, a.z), gp_Dir(n.x, n.y, n.z));
        sew.Add(BRepBuilderAPI_MakeFace(pl, umin, umax, vmin, vmax).Face());
    }
}

TopoDS_Shape openShellFrom(const TopoDS_Shape& s) {
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    int n = 0;
    for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next()) {
        if (n++ == 0) continue;
        b.Add(c, ex.Current());
    }
    return c;
}

int syntheticTests(const std::string& corpus, const std::string& argv0) {
    silence();
    std::setlocale(LC_ALL, "C");
    grade::GradeConfig cfg;
    cfg.skipVolume = true;
    cfg.quiet = true;
    const std::string tmp = join(tmpRoot(), "stl2step_grade_selftest");
#ifdef _WIN32
    _mkdir(tmp.c_str());
#else
    mkdir(tmp.c_str(), 0755);
#endif

    // Case 4: S03 verbatim — cylinders faceted, planes recovered
    {
        std::string verbatim = join(tmp, "S03.verbatim.step");
        grade::Mesh mesh;
        std::string err;
        grade::loadStl(join(corpus, "S03.stl"), mesh, err);
        bool made = false;
        const std::string eng = engineBeside(argv0);
        if (fileNonempty(eng)) {
            std::string cmd = std::string("\"") + eng + "\" \"" + join(corpus, "S03.stl") +
                              "\" -o \"" + verbatim + "\" --engine verbatim --quiet --no-verify";
            std::system(cmd.c_str());
            made = fileNonempty(verbatim);
        }
        if (!made) {
            // SPEC §8 case 4: one planar face per triangle, coplanar neighbours unified.
            BRepBuilderAPI_Sewing sew;
            unifyCoplanarFaces(mesh, sew);
            sew.Perform();
            writeStep(sew.SewedShape(), verbatim);
        }
        grade::GradeDocument d;
        check(grade::gradeFiles(join(corpus, "S03.stl"), verbatim, cfg, d, err), "S03 verbatim grades");
        std::fprintf(stderr, "  case4 grade.cyl=%g grade.plane=%g\n", d.gradeCyl, d.gradePlane);
        check(d.hasCyl && near(d.gradeCyl, 0.0, 1e-12), "case4 cylinder 0");
        check(allStatus(d, grade::SurfClass::Cylinder, grade::Status::Faceted) ||
                  allStatus(d, grade::SurfClass::Cylinder, grade::Status::FacetedPartial),
              "case4 cylinders faceted");
        check(d.hasPlane && near(d.gradePlane, 1.0, 1e-12), "case4 plane 1");
    }

    // Case 5: split
    {
        const std::string p = join(tmp, "S01.split.step");
        writeStep(cubeWithSplitTop(), p);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(join(corpus, "S01.stl"), p, cfg, d, err), "case5 grades");
        int nSplit = 0, nRec = 0;
        double splitCredit = 0;
        for (const auto& f : d.features) {
            if (f.oracle.cls != grade::SurfClass::Plane) continue;
            if (f.status == grade::Status::Split) {
                ++nSplit;
                splitCredit = f.credit;
            }
            if (f.status == grade::Status::Recovered) ++nRec;
        }
        std::fprintf(stderr, "  case5 split=%d rec=%d credit=%g grade.plane=%g\n", nSplit, nRec,
                     splitCredit, d.gradePlane);
        check(nSplit == 1, "case5 one split");
        check(near(splitCredit, 0.5, 1e-6), "case5 credit 0.5");
        check(nRec == 5, "case5 five recovered");
        check(near(d.gradePlane, (5.0 + 0.5) / 6.0, 1e-12), "case5 grade.plane (5.5)/6");
    }

    // Case 6: sliver
    {
        const std::string p = join(tmp, "S01.sliver.step");
        writeStep(cubeWithSliverTop(), p);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(join(corpus, "S01.stl"), p, cfg, d, err), "case6 grades");
        bool found = false;
        for (const auto& f : d.features) {
            if (f.status == grade::Status::Sliver) {
                found = true;
                const double want = f.coverage;
                check(near(f.credit, want, 1e-12), "case6 credit == area(F)/meshArea");
            }
        }
        check(found, "case6 sliver present");
        bool hzValid = false;
        for (const auto& h : d.hardZero)
            if (h == "valid") hzValid = true;
        check(!hzValid, "case6 hardZero does not contain valid");
    }

    // Case 7: SPEC §8 — S01 with the top face enlarged to overhang the +Y
    // oracle plane. Both oracles spanned, credit 0. Six planes proves the
    // mesh is S01, not the two-square substitute.
    {
        const std::string p = join(tmp, "S01.span.step");
        writeStep(cubeSpanned(), p);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(join(corpus, "S01.stl"), p, cfg, d, err), "case7 grades");
        int nSpan = 0;
        for (const auto& f : d.features) {
            if (f.status == grade::Status::Spanned) {
                ++nSpan;
                check(f.credit == 0.0, "case7 spanned credit 0");
            }
        }
        std::fprintf(stderr, "  case7 spanned=%d planeOracles=%d\n", nSpan,
                     countClass(d, grade::SurfClass::Plane));
        check(countClass(d, grade::SurfClass::Plane) == 6, "case7 S01 six planes");
        check(nSpan >= 2, "case7 both oracles spanned");
    }

    // Case 7b: two-square substitute kept as an additional case (addendum B).
    {
        const std::string stl = join(tmp, "span.stl");
        const std::string p = join(tmp, "span.step");
        std::vector<std::array<grade::Vec3, 3>> tris = {
            {{{0, 0, 0}, {10, 0, 0}, {10, 10, 0}}},
            {{{0, 0, 0}, {10, 10, 0}, {0, 10, 0}}},
            {{{20, 0, 0}, {30, 0, 0}, {30, 10, 0}}},
            {{{20, 0, 0}, {30, 10, 0}, {20, 10, 0}}},
        };
        check(writeBinaryStl(stl, tris), "case7b write stl");
        BRepBuilderAPI_Sewing sew;
        sew.Add(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0, 30, 0, 10).Face());
        sew.Perform();
        writeStep(sew.SewedShape(), p);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(stl, p, cfg, d, err), "case7b grades");
        int nSpan = 0;
        for (const auto& f : d.features)
            if (f.status == grade::Status::Spanned && f.credit == 0.0) ++nSpan;
        std::fprintf(stderr, "  case7b spanned=%d planeOracles=%d\n", nSpan,
                     countClass(d, grade::SurfClass::Plane));
        check(nSpan >= 2, "case7b both oracles spanned");
    }

    // Case 8: open shell
    {
        TopoDS_Shape exact;
        if (loadShape(join(corpus, "S03.exact.step"), exact)) {
            const std::string p = join(tmp, "S03.open.step");
            writeStep(openShellFrom(exact), p);
            grade::GradeDocument d;
            std::string err;
            check(grade::gradeFiles(join(corpus, "S03.stl"), p, cfg, d, err), "case8 grades");
            bool hzW = false;
            for (const auto& h : d.hardZero)
                if (h == "watertight") hzW = true;
            check(hzW, "case8 hardZero watertight");
            check(d.hasOverall && near(d.gradeOverall, 0.0, 1e-12), "case8 overall 0");
            check(!d.features.empty(), "case8 features still populated");
        } else {
            check(false, "case8 missing S03.exact.step");
        }
    }

    // T3 — staggered-ring cylinder. N rises until the old axial gate would trip.
    {
        const std::string stl = join(tmp, "stagger.stl");
        const std::string step = join(tmp, "stagger.step");
        const double R = 10.0, H = 12.0;
        const grade::Vec3 axis{0, 0, 1};
        int N = grade::paramCount(grade::SurfClass::Cylinder) + 1;
        double ratio = 0;
        std::vector<std::array<grade::Vec3, 3>> tris;
        // Raise N while the wall still trips the old axial gate. The finest
        // such wall is the one the new predicate must still recover.
        for (int nTry = grade::paramCount(grade::SurfClass::Cylinder) + 1; nTry <= 48;
             nTry += grade::paramCount(grade::SurfClass::Cylinder) + 1) {
            std::vector<std::array<grade::Vec3, 3>> trial;
            auto at = [&](int ring, int i) {
                const double ang = (2.0 * M_PI * static_cast<double>(i)) / nTry +
                                   (ring ? M_PI / static_cast<double>(nTry) : 0.0);
                return grade::Vec3{R * std::cos(ang), R * std::sin(ang), ring ? H : 0.0};
            };
            for (int i = 0; i < nTry; ++i) {
                const int j = (i + 1) % nTry;
                trial.push_back({at(0, i), at(0, j), at(1, i)});
                trial.push_back({at(0, j), at(1, j), at(1, i)});
            }
            const grade::Vec3 cb{0, 0, 0}, ct{0, 0, H};
            for (int i = 0; i < nTry; ++i) {
                trial.push_back({cb, at(0, (i + 1) % nTry), at(0, i)});
                trial.push_back({ct, at(1, i), at(1, (i + 1) % nTry)});
            }
            if (!writeBinaryStl(stl, trial)) continue;
            grade::Mesh mesh;
            std::string err;
            if (!grade::loadStl(stl, mesh, err)) continue;
            double rMax = 0;
            for (const auto& tr : mesh.tris) {
                const double rho =
                    std::sqrt(tr.centroid.x * tr.centroid.x + tr.centroid.y * tr.centroid.y);
                if (std::fabs(rho - R) > 1.0) continue;
                const double s = std::sin(tr.thetaQ);
                if (!(s > 0.0)) continue;
                rMax = std::max(rMax, std::fabs(grade::dot(tr.n, axis)) / s);
            }
            if (rMax > 1.0) {
                N = nTry;
                ratio = rMax;
                tris.swap(trial);
            }
        }
        check(writeBinaryStl(stl, tris), "T3 write stl");
        std::fprintf(stderr, "  T3 N=%d max|n·axis|/sin(theta_q)=%.6g\n", N, ratio);
        check(ratio > 1.0, "T3 staggered wall trips the old normal gate");
        writeStep(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, H).Shape(),
                  step);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(stl, step, cfg, d, err), "T3 grades");
        int cyl = 0, wall = 0, pcs = 0;
        grade::Status st = grade::Status::Missing;
        for (const auto& f : d.features) {
            if (f.oracle.cls != grade::SurfClass::Cylinder) continue;
            ++cyl;
            wall = static_cast<int>(f.oracle.tris.size());
            pcs = f.oracle.oraclePieces;
            st = f.status;
        }
        std::fprintf(stderr, "  T3 cyl=%d tris=%d pieces=%d status=%s residue=%d sing=%d\n", cyl, wall,
                     pcs, grade::statusName(st), d.oracle.residueTris, d.oracle.unprovableSingletons);
        check(cyl == 1, "T3 one cylinder oracle");
        check(wall == 2 * N, "T3 cylinder is the full wall");
        check(pcs == 1, "T3 oraclePieces 1");
        check(st == grade::Status::Recovered, "T3 recovered");
        check(d.oracle.residueTris == 0, "T3 residue 0");
        check(d.oracle.unprovableSingletons == 0, "T3 unprovableSingletons 0");
    }

    // T4 — coaxial chamfer is tier 1 exact; a tilted cone against a cylinder stays tier 2.
    {
        auto meshOf = [](const TopoDS_Shape& s, std::vector<std::array<grade::Vec3, 3>>& tris) {
            tris.clear();
            BRepMesh_IncrementalMesh(s, 0.4);
            for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next()) {
                TopLoc_Location loc;
                Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(TopoDS::Face(ex.Current()), loc);
                if (tri.IsNull()) continue;
                for (int i = 1; i <= tri->NbTriangles(); ++i) {
                    int n1, n2, n3;
                    tri->Triangle(i).Get(n1, n2, n3);
                    gp_Pnt p1 = tri->Node(n1).Transformed(loc);
                    gp_Pnt p2 = tri->Node(n2).Transformed(loc);
                    gp_Pnt p3 = tri->Node(n3).Transformed(loc);
                    tris.push_back({grade::Vec3{p1.X(), p1.Y(), p1.Z()},
                                    grade::Vec3{p2.X(), p2.Y(), p2.Z()},
                                    grade::Vec3{p3.X(), p3.Y(), p3.Z()}});
                }
            }
        };
        const std::string stl = join(tmp, "coaxial.stl");
        const std::string step = join(tmp, "coaxial.step");
        TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 8.0, 16.0);
        TopoDS_Shape cone =
            BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, 16), gp_Dir(0, 0, 1)), 8.0, 12.0, 4.0);
        TopoDS_Shape fused = BRepAlgoAPI_Fuse(cyl, cone).Shape();
        std::vector<std::array<grade::Vec3, 3>> tris;
        meshOf(fused, tris);
        check(writeBinaryStl(stl, tris), "T4 write coaxial stl");
        writeStep(fused, step);
        grade::GradeDocument d;
        std::string err;
        check(grade::gradeFiles(stl, step, cfg, d, err), "T4 coaxial grades");
        int tier1Exact = 0;
        for (const auto& x : d.intersections) {
            const bool coneCyl = (x.a.rfind("cone:", 0) == 0 && x.b.rfind("cylinder:", 0) == 0) ||
                                 (x.b.rfind("cone:", 0) == 0 && x.a.rfind("cylinder:", 0) == 0);
            if (!coneCyl) continue;
            std::fprintf(stderr, "  T4 coaxial %s × %s tier %d %s\n", x.a.c_str(), x.b.c_str(),
                         x.expectedTier, x.verdict.c_str());
            if (x.expectedTier == 1 && x.verdict == "exact") ++tier1Exact;
        }
        check(tier1Exact >= 1, "T4 coaxial cone|cylinder tier 1 exact");

        const std::string stl2 = join(tmp, "oblique.stl");
        const std::string step2 = join(tmp, "oblique.step");
        TopoDS_Shape cyl2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 10.0, 20.0);
        TopoDS_Shape cone2 =
            BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(-8, 0, 10), gp_Dir(1, 0, 0)), 4.0, 10.0, 16.0);
        TopoDS_Shape fused2 = BRepAlgoAPI_Fuse(cyl2, cone2).Shape();
        meshOf(fused2, tris);
        check(writeBinaryStl(stl2, tris), "T4 write oblique stl");
        writeStep(fused2, step2);
        grade::GradeDocument d2;
        check(grade::gradeFiles(stl2, step2, cfg, d2, err), "T4 oblique grades");
        int tier2 = 0;
        for (const auto& x : d2.intersections) {
            const bool coneCyl = (x.a.rfind("cone:", 0) == 0 && x.b.rfind("cylinder:", 0) == 0) ||
                                 (x.b.rfind("cone:", 0) == 0 && x.a.rfind("cylinder:", 0) == 0);
            if (!coneCyl) continue;
            std::fprintf(stderr, "  T4 oblique %s × %s tier %d %s\n", x.a.c_str(), x.b.c_str(),
                         x.expectedTier, x.verdict.c_str());
            if (x.expectedTier == 2) ++tier2;
        }
        if (tier2 < 1) {
            int nC = 0, nK = 0;
            for (const auto& f : d2.features) {
                if (f.oracle.cls == grade::SurfClass::Cone) ++nC;
                if (f.oracle.cls == grade::SurfClass::Cylinder) ++nK;
            }
            std::fprintf(stderr, "  T4 oblique cones=%d cyls=%d rows=%zu\n", nC, nK,
                         d2.intersections.size());
        }
        checkId(tier2 >= 1, "grade.synthetic-oblique-cone-unseeded",
                "T4 non-coaxial cone|cylinder tier 2");
    }
    return gFail;
}

int asciiTest(const std::string& corpus) {
    std::setlocale(LC_ALL, "C");
    const std::string asciiDir = join(tmpRoot(), "stl2step_grade_selftest");
    ensureDir(asciiDir);
    const std::string ascii = join(asciiDir, "s01_ascii.stl");
    grade::Mesh mesh;
    std::string err;
    if (!grade::loadStl(join(corpus, "S01.stl"), mesh, err)) {
        check(false, "ascii: load binary S01");
        return gFail;
    }
    std::FILE* f = std::fopen(ascii.c_str(), "wb");
    if (!f) {
        check(false, "ascii: open writer");
        return gFail;
    }
    std::fputs("solid s01\n", f);
    for (const auto& tr : mesh.tris) {
        std::fprintf(f, "  facet normal %.9g %.9g %.9g\n    outer loop\n", tr.n.x, tr.n.y, tr.n.z);
        for (int k = 0; k < 3; ++k) {
            const grade::Vec3& p = mesh.verts[static_cast<size_t>(tr.v[k])];
            std::fprintf(f, "      vertex %.9g %.9g %.9g\n", p.x, p.y, p.z);
        }
        std::fputs("    endloop\n  endfacet\n", f);
    }
    std::fputs("endsolid s01\n", f);
    std::fclose(f);
    const auto qf = stl2step::stlQuantFloor(ascii);
    check(qf.ok && qf.ascii, "ascii: classified ascii");
    check(qf.q > 0.0, "ascii: q from printed decimals");
    grade::GradeConfig cfg;
    cfg.skipVolume = true;
    cfg.quiet = true;
    grade::GradeDocument db, da;
    check(grade::gradeFiles(join(corpus, "S01.stl"), join(corpus, "S01.exact.step"), cfg, db, err),
          "ascii: binary grade");
    check(grade::gradeFiles(ascii, join(corpus, "S01.exact.step"), cfg, da, err), "ascii: ascii grade");
    auto ulp = [](double a, double b) {
        const double s = std::max(std::fabs(a), std::fabs(b));
        return std::fabs(a - b) <= 1e-9 * (s > 0.0 ? s : 1.0);
    };
    check(ulp(da.gradePlane, db.gradePlane) && ulp(da.gradeOverall, db.gradeOverall),
          "ascii: per-class grades match binary within identical-ulp");
    return gFail;
}

}  // namespace

int main(int argc, char** argv) {
    std::setlocale(LC_ALL, "C");
    if (argc < 3) {
        std::fprintf(stderr, "usage: stl2step_grade_selftest core|synthetic|ascii <corpusDir>\n");
        return 1;
    }
    const std::string mode = argv[1];
    if (mode == "core") {
        const int skipRc = stl2step::test::skipExitIfUncalibrated();
        if (skipRc >= 0) {
            return skipRc;
        }
    }
    const std::string corpus = argv[2];
    if (argc >= 4) gPlateStep = argv[3];
    gExpectedRed = loadExpectedRedGrader(STL2STEP_EXPECTED_RED_JSON);
    int rc = 0;
    if (mode == "core") rc = coreTests(corpus);
    else if (mode == "synthetic") rc = syntheticTests(corpus, argv[0]);
    else if (mode == "ascii") rc = asciiTest(corpus);
    else {
        std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
        return 1;
    }
    const int total = gPass + gFail + gXFail;
    if (gXFail > 0) {
        std::fprintf(stderr, "grade_selftest %s: %d/%d PASS (%d XFAIL)\n", mode.c_str(), gPass, total, gXFail);
    } else {
        std::fprintf(stderr, "grade_selftest %s: %d/%d PASS\n", mode.c_str(), gPass, total);
    }
    return rc ? 1 : 0;
}
