#include "grade_compare.hpp"

#include "grade_sha256.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GeomAbs_CurveType.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_ResultType.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>

namespace grade {
namespace {

int statusIndex(Status s) {
    switch (s) {
        case Status::Recovered: return 0;
        case Status::Split: return 1;
        case Status::Sliver: return 2;
        case Status::Faceted: return 3;
        case Status::FacetedPartial: return 4;
        case Status::Spanned: return 5;
        default: return 6;
    }
}

gp_Pln toPln(const SurfParams& S) {
    return gp_Pln(gp_Pnt(S.p0.x, S.p0.y, S.p0.z), gp_Dir(S.n.x, S.n.y, S.n.z));
}
gp_Cylinder toCyl(const SurfParams& S) {
    return gp_Cylinder(gp_Ax3(gp_Pnt(S.p0.x, S.p0.y, S.p0.z), gp_Dir(S.n.x, S.n.y, S.n.z)), S.R);
}
gp_Cone toCone(const SurfParams& S) {
    return gp_Cone(gp_Ax3(gp_Pnt(S.apex.x, S.apex.y, S.apex.z), gp_Dir(S.n.x, S.n.y, S.n.z)),
                   S.alpha, 0.0);
}
gp_Sphere toSph(const SurfParams& S) {
    return gp_Sphere(gp_Ax3(gp_Pnt(S.p0.x, S.p0.y, S.p0.z), gp_Dir(0, 0, 1)), S.R);
}
gp_Torus toTor(const SurfParams& S) {
    return gp_Torus(gp_Ax3(gp_Pnt(S.p0.x, S.p0.y, S.p0.z), gp_Dir(S.n.x, S.n.y, S.n.z)), S.R, S.r);
}

bool sameSurface(const SurfParams& a, const SurfParams& b, const Oracle& O, const Mesh& m) {
    if (a.cls != b.cls) return false;
    for (int vi : O.verts) {
        if (distToSurf(m.verts[static_cast<size_t>(vi)], b) > m.tau) return false;
    }
    return true;
}

double maxDevOn(const Oracle& O, const SurfParams& G, const Mesh& m) {
    double mx = 0;
    for (int vi : O.verts)
        mx = std::max(mx, distToSurf(m.verts[static_cast<size_t>(vi)], G));
    return mx;
}

double sagittaVolume(const Mesh& m, int t, const SurfParams& S) {
    const Tri& tr = m.tris[static_cast<size_t>(t)];
    if (S.cls == SurfClass::Plane) return 0.0;
    const Vec3& a = m.verts[static_cast<size_t>(tr.v[0])];
    const Vec3& b = m.verts[static_cast<size_t>(tr.v[1])];
    const Vec3& c = m.verts[static_cast<size_t>(tr.v[2])];
    auto chordPerp = [&](const Vec3& axis) {
        auto proj = [&](const Vec3& p) {
            const Vec3 w = p - S.p0;
            return w - axis * dot(w, axis);
        };
        const Vec3 pa = proj(a), pb = proj(b), pc = proj(c);
        const double l0 = dist(pb, pc), l1 = dist(pa, pc), l2 = dist(pa, pb);
        return std::max(l0, std::max(l1, l2));
    };
    const double twoThirds = 2.0 / 3.0;
    switch (S.cls) {
        case SurfClass::Cylinder: {
            const double w = chordPerp(S.n);
            const double s = (S.R > 0.0) ? S.R * (1.0 - std::cos(w / (2.0 * S.R))) : 0.0;
            return tr.area * twoThirds * s;
        }
        case SurfClass::Cone: {
            const Vec3 u = tr.centroid - S.apex;
            const double ax = dot(u, S.n);
            const double localR = std::fabs(ax) * std::tan(std::fabs(S.alpha));
            const double w = chordPerp(S.n);
            const double s = (localR > 0.0) ? localR * (1.0 - std::cos(w / (2.0 * localR))) : 0.0;
            return tr.area * twoThirds * s;
        }
        case SurfClass::Torus: {
            const double w = chordPerp(S.n);
            const double s = (S.r > 0.0) ? S.r * (1.0 - std::cos(w / (2.0 * S.r))) : 0.0;
            return tr.area * twoThirds * s;
        }
        case SurfClass::Sphere: {
            const double s = S.R - distToSurf(tr.centroid, S) - 0.0;
            const double sag = S.R - dist(S.p0, tr.centroid);  // R - dist(c, plane-ish)
            const double ss = std::max(0.0, S.R - std::fabs(norm(tr.centroid - S.p0)));
            (void)s;
            (void)sag;
            return tr.area * twoThirds * ss;
        }
        default:
            return 0.0;
    }
}

int expectedTier(const Oracle& A, const Oracle& B, double tau) {
    try {
        IntAna_QuadQuadGeo q;
        auto run = [&]() -> int {
            if (!q.IsDone()) return 2;
            const IntAna_ResultType ty = q.TypeInter();
            if (ty == IntAna_Line || ty == IntAna_Circle || ty == IntAna_Ellipse) return 1;
            return 2;
        };
        if (A.cls == SurfClass::Plane && B.cls == SurfClass::Plane) {
            q.Perform(toPln(A.S), toPln(B.S), tau, tau);
            return run();
        }
        if (A.cls == SurfClass::Plane && B.cls == SurfClass::Cylinder) {
            q.Perform(toPln(A.S), toCyl(B.S), tau, tau);
            return run();
        }
        if (A.cls == SurfClass::Cylinder && B.cls == SurfClass::Plane) {
            q.Perform(toPln(B.S), toCyl(A.S), tau, tau);
            return run();
        }
        if (A.cls == SurfClass::Cylinder && B.cls == SurfClass::Cylinder) {
            q.Perform(toCyl(A.S), toCyl(B.S), tau);
            return run();
        }
        if (A.cls == SurfClass::Plane && B.cls == SurfClass::Cone) {
            q.Perform(toPln(A.S), toCone(B.S), tau, tau);
            return run();
        }
        if (A.cls == SurfClass::Cone && B.cls == SurfClass::Plane) {
            q.Perform(toPln(B.S), toCone(A.S), tau, tau);
            return run();
        }
        if (A.cls == SurfClass::Cylinder && B.cls == SurfClass::Cone) {
            q.Perform(toCyl(A.S), toCone(B.S), tau);
            return run();
        }
        if (A.cls == SurfClass::Cone && B.cls == SurfClass::Cylinder) {
            q.Perform(toCyl(B.S), toCone(A.S), tau);
            return run();
        }
        if (A.cls == SurfClass::Sphere && B.cls == SurfClass::Plane) {
            q.Perform(toPln(B.S), toSph(A.S));
            return run();
        }
        if (A.cls == SurfClass::Plane && B.cls == SurfClass::Sphere) {
            q.Perform(toPln(A.S), toSph(B.S));
            return run();
        }
    } catch (...) {
        return 2;
    }
    return 2;
}

std::string classifyEdge(const TopoDS_Edge& e) {
    BRepAdaptor_Curve ac(e);
    switch (ac.GetType()) {
        case GeomAbs_Line: return "LINE";
        case GeomAbs_Circle: return "CIRCLE";
        case GeomAbs_Ellipse: return "ELLIPSE";
        case GeomAbs_Hyperbola: return "HYPERBOLA";
        case GeomAbs_Parabola: return "PARABOLA";
        case GeomAbs_BSplineCurve: {
            Handle(Geom_BSplineCurve) bs = ac.BSpline();
            if (!bs.IsNull() && bs->Degree() == 1) return "BSPLINE(1)";
            return "BSPLINE";
        }
        default: return "OTHER";
    }
}

bool analyticName(const std::string& s) {
    return s == "LINE" || s == "CIRCLE" || s == "ELLIPSE" || s == "HYPERBOLA" || s == "PARABOLA";
}

bool shareMeshEdge(const Mesh& m, const Oracle& A, const Oracle& B) {
    std::unordered_set<uint64_t> ea;
    auto ek = [](int i, int j) -> uint64_t {
        if (i > j) std::swap(i, j);
        return (uint64_t(uint32_t(i)) << 32) | uint32_t(j);
    };
    for (int t : A.tris) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        for (int s = 0; s < 3; ++s) ea.insert(ek(tr.v[s], tr.v[(s + 1) % 3]));
    }
    for (int t : B.tris) {
        const Tri& tr = m.tris[static_cast<size_t>(t)];
        for (int s = 0; s < 3; ++s)
            if (ea.count(ek(tr.v[s], tr.v[(s + 1) % 3]))) return true;
    }
    return false;
}

bool identUlpFiles(const std::string& a, const std::string& b) {
    std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
    if (!fa || !fb) return false;
    std::string sa((std::istreambuf_iterator<char>(fa)), std::istreambuf_iterator<char>());
    std::string sb((std::istreambuf_iterator<char>(fb)), std::istreambuf_iterator<char>());
    if (sa == sb) return true;
    // Walk both, rounding numeric tokens to 1e-9 relative.
    auto nextNum = [](const std::string& s, size_t& i, bool& isNum, double& v) {
        while (i < s.size() && !(s[i] == '-' || s[i] == '+' || s[i] == '.' ||
                                 (s[i] >= '0' && s[i] <= '9'))) {
            ++i;
        }
        if (i >= s.size()) {
            isNum = false;
            return;
        }
        const size_t start = i;
        if (s[i] == '+' || s[i] == '-') ++i;
        bool seen = false;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            ++i;
            seen = true;
        }
        if (i < s.size() && s[i] == '.') {
            ++i;
            while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
                ++i;
                seen = true;
            }
        }
        if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
            ++i;
            if (i < s.size() && (s[i] == '+' || s[i] == '-')) ++i;
            while (i < s.size() && s[i] >= '0' && s[i] <= '9') ++i;
        }
        if (!seen) {
            isNum = false;
            i = start + 1;
            return;
        }
        isNum = true;
        v = std::strtod(s.c_str() + static_cast<int>(start), nullptr);
    };
    size_t ia = 0, ib = 0;
    while (ia < sa.size() && ib < sb.size()) {
        bool na = false, nb = false;
        double va = 0, vb = 0;
        const size_t ia0 = ia, ib0 = ib;
        nextNum(sa, ia, na, va);
        nextNum(sb, ib, nb, vb);
        if (na != nb) return false;
        if (!na) break;
        const std::string preA = sa.substr(ia0, ia - ia0);
        (void)preA;
        if (na && nb) {
            const double scale = std::max(std::fabs(va), std::fabs(vb));
            const double tol = std::max(1e-9 * scale, 0.0);
            if (std::fabs(va - vb) > tol) return false;
        }
    }
    return true;
}

}  // namespace

bool gradeFiles(const std::string& stl, const std::string& step, const GradeConfig& cfg,
                GradeDocument& doc, std::string& err) {
    doc = GradeDocument{};
    doc.stlPath = stl;
    doc.stepPath = step;
    doc.stlSha = sha256File(stl);
    doc.stepSha = sha256File(step);
    doc.stlBytes = fileSize(stl);
    doc.stepBytes = fileSize(step);
    if (!loadStl(stl, doc.mesh, err)) return false;
    buildOracle(doc.mesh, doc.oracle, cfg.reverseSeeds);
    if (!loadStep(step, doc.mesh, doc.step, err, !cfg.skipVolume)) return false;

    doc.watertight = doc.step.watertight;
    doc.valid = doc.step.valid;
    doc.Vmesh = doc.mesh.volume;
    doc.Vstep = doc.step.volume;
    doc.volumeQ = doc.mesh.q * doc.mesh.surfaceArea;
    doc.chordBudget = doc.volumeQ;
    for (const Oracle& o : doc.oracle.oracles) {
        for (int t : o.tris) doc.chordBudget += sagittaVolume(doc.mesh, t, o.S);
    }
    doc.volumeDelta = doc.Vstep - doc.Vmesh;
    doc.volumeWithin = true;
    if (!cfg.skipVolume) {
        if (doc.volumeDelta < -doc.volumeQ || doc.volumeDelta > doc.chordBudget)
            doc.volumeWithin = false;
    }

    // Group remaining (non-facet) faces by canonical surface.
    std::vector<int> groupOf(doc.step.faces.size(), -1);
    std::vector<std::vector<int>> groups;
    for (int i = 0; i < static_cast<int>(doc.step.faces.size()); ++i) {
        const StepFace& F = doc.step.faces[static_cast<size_t>(i)];
        if (F.facet) continue;
        int found = -1;
        for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
            const StepFace& G0 = doc.step.faces[static_cast<size_t>(groups[static_cast<size_t>(g)][0])];
            if (G0.cls != F.cls) continue;
            // every sampled point of one within tau of the other
            bool ok = true;
            for (const Vec3& p : F.wireVerts) {
                if (distToSurf(p, G0.S) > doc.mesh.tau) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                for (const Vec3& p : G0.wireVerts) {
                    if (distToSurf(p, F.S) > doc.mesh.tau) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                found = g;
                break;
            }
        }
        if (found < 0) {
            found = static_cast<int>(groups.size());
            groups.push_back({});
        }
        groups[static_cast<size_t>(found)].push_back(i);
        groupOf[static_cast<size_t>(i)] = found;
    }

    // Assignment: one-to-one, smallest max-dev, then lower entity id.
    std::vector<int> assign(doc.oracle.oracles.size(), -1);  // group index
    std::vector<char> taken(groups.size(), 0);
    struct Cand {
        int oi, g;
        double dev;
        int entity;
    };
    std::vector<Cand> cands;
    for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
        const Oracle& O = doc.oracle.oracles[static_cast<size_t>(oi)];
        for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
            const StepFace& G0 = doc.step.faces[static_cast<size_t>(groups[static_cast<size_t>(g)][0])];
            if (G0.cls != O.cls) continue;
            if (!sameSurface(O.S, G0.S, O, doc.mesh)) continue;
            int minEnt = G0.entity;
            for (int fi : groups[static_cast<size_t>(g)])
                minEnt = std::min(minEnt, doc.step.faces[static_cast<size_t>(fi)].entity);
            cands.push_back(Cand{oi, g, maxDevOn(O, G0.S, doc.mesh), minEnt});
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
        if (a.dev != b.dev) return a.dev < b.dev;
        if (a.entity != b.entity) return a.entity < b.entity;
        return a.oi < b.oi;
    });
    for (const Cand& c : cands) {
        if (assign[static_cast<size_t>(c.oi)] >= 0) continue;
        if (taken[static_cast<size_t>(c.g)]) continue;
        assign[static_cast<size_t>(c.oi)] = c.g;
        taken[static_cast<size_t>(c.g)] = 1;
    }

    // Spanned: a non-facet face whose untrimmed surface is within tau of two oracles' vertices.
    std::vector<char> spanned(doc.oracle.oracles.size(), 0);
    for (const StepFace& F : doc.step.faces) {
        if (F.facet) continue;
        std::vector<int> hit;
        for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
            const Oracle& O = doc.oracle.oracles[static_cast<size_t>(oi)];
            if (O.cls != F.cls) continue;
            bool all = true;
            for (int vi : O.verts) {
                if (distToSurf(doc.mesh.verts[static_cast<size_t>(vi)], F.S) > doc.mesh.tau) {
                    all = false;
                    break;
                }
            }
            if (all) hit.push_back(oi);
        }
        if (hit.size() >= 2) {
            for (int oi : hit) spanned[static_cast<size_t>(oi)] = 1;
        }
    }

    doc.features.resize(doc.oracle.oracles.size());
    const double totalOracleArea = [&]() {
        double s = 0;
        for (const Oracle& o : doc.oracle.oracles) s += o.w;
        return s;
    }();

    for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
        Feature& f = doc.features[static_cast<size_t>(oi)];
        f.oracle = doc.oracle.oracles[static_cast<size_t>(oi)];
        f.areaFraction = (totalOracleArea > 0.0) ? f.oracle.w / totalOracleArea : 0.0;
        const int g = assign[static_cast<size_t>(oi)];
        if (spanned[static_cast<size_t>(oi)]) {
            f.status = Status::Spanned;
            f.credit = 0;
        } else if (g >= 0) {
            const auto& members = groups[static_cast<size_t>(g)];
            double sumA = 0, amax = 0;
            for (int fi : members) {
                const StepFace& F = doc.step.faces[static_cast<size_t>(fi)];
                f.stepFaces.push_back(StepFaceRef{F.entity, F.area});
                sumA += F.area;
                amax = std::max(amax, F.area);
            }
            std::sort(f.stepFaces.begin(), f.stepFaces.end(),
                      [](const StepFaceRef& a, const StepFaceRef& b) { return a.entity < b.entity; });
            f.coverage = (f.oracle.w > 0.0) ? (sumA / f.oracle.w) : 0.0;
            const double aq = areaQ(doc.mesh, f.oracle.tris);
            if (members.size() > 1) {
                f.status = Status::Split;
                f.splitN = static_cast<int>(members.size());
                f.credit = (sumA > 0.0) ? (amax / sumA) : 0.0;
            } else if (sumA >= f.oracle.w - aq) {
                f.status = Status::Recovered;
                f.credit = 1.0;
            } else {
                f.status = Status::Sliver;
                f.credit = (f.oracle.w > 0.0) ? (sumA / f.oracle.w) : 0.0;
            }
        } else {
            int cov = 0;
            for (int t : f.oracle.tris)
                if (doc.step.cover[static_cast<size_t>(t)] >= 0) ++cov;
            if (cov == static_cast<int>(f.oracle.tris.size()) && cov > 0) {
                f.status = Status::Faceted;
            } else if (cov > 0) {
                f.status = Status::FacetedPartial;
            } else {
                f.status = Status::Missing;
            }
            f.credit = 0;
        }

        auto bump = [&](double& area, double& rec, int& n, int* cnt) {
            area += f.oracle.w;
            rec += f.oracle.w * f.credit;
            ++n;
            ++cnt[statusIndex(f.status)];
        };
        switch (f.oracle.cls) {
            case SurfClass::Plane:
                bump(doc.areaPlane, doc.recPlane, doc.nPlane, doc.cntPlane);
                break;
            case SurfClass::Cylinder:
                bump(doc.areaCyl, doc.recCyl, doc.nCyl, doc.cntCyl);
                break;
            case SurfClass::Cone:
                bump(doc.areaCone, doc.recCone, doc.nCone, doc.cntCone);
                break;
            case SurfClass::Sphere:
                bump(doc.areaSphere, doc.recSphere, doc.nSphere, doc.cntSphere);
                break;
            case SurfClass::Torus:
                bump(doc.areaTorus, doc.recTorus, doc.nTorus, doc.cntTorus);
                break;
            default:
                break;
        }
    }

    auto setG = [](double rec, double area, double& g, bool& has) {
        if (area > 0.0) {
            g = rec / area;
            has = true;
        } else {
            g = 0;
            has = false;
        }
    };
    setG(doc.recPlane, doc.areaPlane, doc.gradePlane, doc.hasPlane);
    setG(doc.recCyl, doc.areaCyl, doc.gradeCyl, doc.hasCyl);
    setG(doc.recCone, doc.areaCone, doc.gradeCone, doc.hasCone);
    setG(doc.recSphere, doc.areaSphere, doc.gradeSphere, doc.hasSphere);
    setG(doc.recTorus, doc.areaTorus, doc.gradeTorus, doc.hasTorus);
    const double recAll =
        doc.recPlane + doc.recCyl + doc.recCone + doc.recSphere + doc.recTorus;
    const double areaAll =
        doc.areaPlane + doc.areaCyl + doc.areaCone + doc.areaSphere + doc.areaTorus;
    setG(recAll, areaAll, doc.gradeOverall, doc.hasOverall);

    // Intersections
    for (int i = 0; i < static_cast<int>(doc.oracle.oracles.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(doc.oracle.oracles.size()); ++j) {
            const Oracle& A = doc.oracle.oracles[static_cast<size_t>(i)];
            const Oracle& B = doc.oracle.oracles[static_cast<size_t>(j)];
            if (!shareMeshEdge(doc.mesh, A, B)) continue;
            Intersection ix;
            ix.a = A.featureId;
            ix.b = B.featureId;
            if (ix.a > ix.b) std::swap(ix.a, ix.b);
            ix.expectedTier = expectedTier(A, B, doc.mesh.tau);
            // shared EDGE of assigned faces
            const Feature& fa = doc.features[static_cast<size_t>(i)];
            const Feature& fb = doc.features[static_cast<size_t>(j)];
            ix.shipped = "absent";
            ix.verdict = "absent";
            if (!fa.stepFaces.empty() && !fb.stepFaces.empty()) {
                std::unordered_map<int, int> edgeCount;
                auto addEdges = [&](int entity) {
                    const StepFace* fp = nullptr;
                    for (const StepFace& F : doc.step.faces)
                        if (F.entity == entity) {
                            fp = &F;
                            break;
                        }
                    if (!fp) return;
                    TopTools_IndexedMapOfShape emap;
                    TopExp::MapShapes(fp->face, TopAbs_EDGE, emap);
                    for (int e = 1; e <= emap.Extent(); ++e) ++edgeCount[e];  // not unique across faces
                };
                (void)addEdges;
                TopTools_IndexedMapOfShape ea, eb;
                const StepFace* Fa = nullptr;
                const StepFace* Fb = nullptr;
                for (const StepFace& F : doc.step.faces) {
                    if (F.entity == fa.stepFaces[0].entity) Fa = &F;
                    if (F.entity == fb.stepFaces[0].entity) Fb = &F;
                }
                if (Fa && Fb) {
                    TopExp::MapShapes(Fa->face, TopAbs_EDGE, ea);
                    TopExp::MapShapes(Fb->face, TopAbs_EDGE, eb);
                    for (int ei = 1; ei <= ea.Extent(); ++ei) {
                        for (int ej = 1; ej <= eb.Extent(); ++ej) {
                            if (ea(ei).IsSame(eb(ej))) {
                                ix.shipped = classifyEdge(TopoDS::Edge(ea(ei)));
                                ix.edgeEntity = ei;
                                goto found_edge;
                            }
                        }
                    }
                found_edge:;
                }
            }
            if (ix.shipped == "absent") {
                ix.verdict = "absent";
            } else if (analyticName(ix.shipped)) {
                ix.verdict = (ix.expectedTier == 1) ? "exact" : "other";
            } else if (ix.shipped == "BSPLINE(1)") {
                ix.verdict = (ix.expectedTier == 1) ? "downgraded" : "polyline";
            } else {
                ix.verdict = "other";
            }
            doc.intersections.push_back(std::move(ix));
        }
    }
    std::sort(doc.intersections.begin(), doc.intersections.end(),
              [](const Intersection& a, const Intersection& b) {
                  if (a.a != b.a) return a.a < b.a;
                  return a.b < b.b;
              });

    // unmatched non-facet faces not in any assigned group
    std::vector<char> usedFace(doc.step.faces.size(), 0);
    for (const Feature& f : doc.features)
        for (const StepFaceRef& r : f.stepFaces)
            for (int i = 0; i < static_cast<int>(doc.step.faces.size()); ++i)
                if (doc.step.faces[static_cast<size_t>(i)].entity == r.entity)
                    usedFace[static_cast<size_t>(i)] = 1;
    for (int i = 0; i < static_cast<int>(doc.step.faces.size()); ++i) {
        const StepFace& F = doc.step.faces[static_cast<size_t>(i)];
        if (F.facet || usedFace[static_cast<size_t>(i)]) continue;
        UnmatchedFace u;
        u.entity = F.entity;
        u.type = (F.cls == SurfClass::Other) ? "BSpline" : className(F.cls);
        u.areaMM2 = F.area;
        doc.unmatched.push_back(u);
    }

    int resCov = 0;
    for (int t = 0; t < static_cast<int>(doc.mesh.tris.size()); ++t) {
        if (doc.oracle.owner[static_cast<size_t>(t)] < 0 &&
            doc.step.cover[static_cast<size_t>(t)] >= 0)
            ++resCov;
    }
    doc.residueFacetedFraction =
        (doc.oracle.residueTris > 0)
            ? (static_cast<double>(resCov) / static_cast<double>(doc.oracle.residueTris))
            : 1.0;

    // features order: area desc, featureId asc
    std::sort(doc.features.begin(), doc.features.end(), [](const Feature& a, const Feature& b) {
        if (a.oracle.w != b.oracle.w) return a.oracle.w > b.oracle.w;
        return a.oracle.featureId < b.oracle.featureId;
    });

    doc.determinism = "skipped";
    if (!cfg.engineBin.empty()) {
        // convert twice; compare under identical-ulp
        auto runOnce = [&](const std::string& outp, const char* extra) -> bool {
            std::ostringstream cmd;
            cmd << '"' << cfg.engineBin << "\" \"" << stl << "\" -o \"" << outp
                << "\" --quiet --no-verify";
            if (extra) cmd << ' ' << extra;
            return std::system(cmd.str().c_str()) == 0 || true;  // exit 2 is ok
        };
        const std::string a = step + ".det_a.step";
        const std::string b = step + ".det_b.step";
        runOnce(a, "--threads 1");
        runOnce(b, nullptr);
        doc.determinism = identUlpFiles(a, b) ? "pass" : "fail";
        std::remove(a.c_str());
        std::remove(b.c_str());
    }

    if (!doc.watertight) doc.hardZero.push_back("watertight");
    if (!doc.valid) doc.hardZero.push_back("valid");
    if (doc.determinism == "fail") doc.hardZero.push_back("determinism");
    if (!cfg.skipVolume && !doc.volumeWithin) doc.hardZero.push_back("volume");
    if (!doc.hardZero.empty() && doc.hasOverall) doc.gradeOverall = 0.0;

    return true;
}

}  // namespace grade
