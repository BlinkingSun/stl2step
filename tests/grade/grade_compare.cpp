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
#include <BRepTopAdaptor_FClass2d.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_ResultType.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt2d.hxx>
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

    // Verbatim cylinder walls unify to |S|>1 planar faces. Those are
    // tessellation facets of a curved oracle, not design planes. Mark them
    // so curved oracles grade as faceted (SPEC §8 case 4).
    for (StepFace& F : doc.step.faces) {
        if (F.facet || F.cls != SurfClass::Plane) continue;
        std::unordered_set<int> wire(F.meshVerts.begin(), F.meshVerts.end());
        wire.erase(-1);
        if (wire.empty()) continue;
        std::vector<int> S;
        for (int t = 0; t < static_cast<int>(doc.mesh.tris.size()); ++t) {
            const Tri& tr = doc.mesh.tris[static_cast<size_t>(t)];
            if (!wire.count(tr.v[0]) || !wire.count(tr.v[1]) || !wire.count(tr.v[2])) continue;
            bool on = true;
            for (int k = 0; k < 3; ++k)
                if (distToSurf(doc.mesh.verts[static_cast<size_t>(tr.v[k])], F.S) > doc.mesh.tau) {
                    on = false;
                    break;
                }
            if (on) S.push_back(t);
        }
        if (S.size() < 2) continue;
        bool curved = !S.empty();
        for (int t : S) {
            const int ow = doc.oracle.owner[static_cast<size_t>(t)];
            if (ow < 0 || doc.oracle.oracles[static_cast<size_t>(ow)].cls == SurfClass::Plane) {
                curved = false;
                break;
            }
        }
        if (!curved) continue;
        F.facet = true;
        ++doc.step.nFacet;
        for (int t : S) doc.step.cover[static_cast<size_t>(t)] = F.entity;
    }
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
        // Two-sided chord budget: bosses make V_step > V_mesh; holes make
        // V_mesh > V_step by the same sagitta sum already in chordBudget.
        if (doc.volumeDelta < -doc.chordBudget || doc.volumeDelta > doc.chordBudget)
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

    // Assignment: a coplanar group of disjoint faces (S03 z=4 discs) must
    // not be handed to a single oracle as split(n). One matching oracle
    // keeps the whole group (case 5 split). Several matching oracles share
    // the faces by point-in-face overlap (SPEC §5.5 vs §6.2).
    std::vector<std::vector<int>> facesOf(doc.oracle.oracles.size());
    auto overlapFace = [&](const Oracle& O, const StepFace& F) -> double {
        std::unordered_set<int> wire(F.meshVerts.begin(), F.meshVerts.end());
        wire.erase(-1);
        double a = 0;
        if (!wire.empty()) {
            for (int t : O.tris) {
                const Tri& tr = doc.mesh.tris[static_cast<size_t>(t)];
                int hit = 0;
                for (int k = 0; k < 3; ++k)
                    if (wire.count(tr.v[k])) ++hit;
                if (hit >= 2) a += tr.area;
            }
        }
        if (a > 0.0 || F.face.IsNull()) return a;
        // Analytic circular discs (exact.step) have a CIRCLE edge and one
        // seam vertex — mesh-vert identity is empty. Classify centroids.
        Handle(Geom_Surface) surf = BRep_Tool::Surface(F.face);
        if (surf.IsNull()) return 0.0;
        BRepTopAdaptor_FClass2d cls(F.face, Precision::PConfusion());
        for (int t : O.tris) {
            const Vec3& c = doc.mesh.tris[static_cast<size_t>(t)].centroid;
            GeomAPI_ProjectPointOnSurf proj(gp_Pnt(c.x, c.y, c.z), surf);
            if (proj.NbPoints() < 1) continue;
            Standard_Real u = 0, v = 0;
            proj.LowerDistanceParameters(u, v);
            const TopAbs_State st = cls.Perform(gp_Pnt2d(u, v));
            if (st == TopAbs_IN || st == TopAbs_ON) a += doc.mesh.tris[static_cast<size_t>(t)].area;
        }
        return a;
    };
    for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
        const StepFace& G0 = doc.step.faces[static_cast<size_t>(groups[static_cast<size_t>(g)][0])];
        std::vector<int> match;
        for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
            const Oracle& O = doc.oracle.oracles[static_cast<size_t>(oi)];
            if (O.cls != G0.cls) continue;
            if (!sameSurface(O.S, G0.S, O, doc.mesh)) continue;
            match.push_back(oi);
        }
        if (match.size() == 1) {
            facesOf[static_cast<size_t>(match[0])].insert(facesOf[static_cast<size_t>(match[0])].end(),
                                                         groups[static_cast<size_t>(g)].begin(),
                                                         groups[static_cast<size_t>(g)].end());
        } else if (match.size() > 1) {
            for (int fi : groups[static_cast<size_t>(g)]) {
                const StepFace& F = doc.step.faces[static_cast<size_t>(fi)];
                int best = -1;
                double bestA = 0;
                std::string bestId;
                for (int oi : match) {
                    const double a = overlapFace(doc.oracle.oracles[static_cast<size_t>(oi)], F);
                    if (a <= 0.0) continue;
                    const std::string& id = doc.oracle.oracles[static_cast<size_t>(oi)].featureId;
                    if (a > bestA || (a == bestA && (best < 0 || id < bestId))) {
                        bestA = a;
                        best = oi;
                        bestId = id;
                    }
                }
                if (best >= 0) facesOf[static_cast<size_t>(best)].push_back(fi);
            }
        }
    }

    // Spanned: one STEP face overlaps two oracle surfaces (mesh-vertex
    // identity, not unordered-wire even-odd — that marked S03's five z=4
    // discs spanned of each other).
    // SPEC §8 case 7: an enlarged face that overhangs a neighbouring
    // oracle (S01 top UV past a shared edge) also spans both — neighbour
    // centroid projects strictly IN (not ON) F. Exact cube edges are ON.
    std::vector<char> spanned(doc.oracle.oracles.size(), 0);
    for (const StepFace& F : doc.step.faces) {
        if (F.facet) continue;
        std::vector<int> hit;
        for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
            const Oracle& O = doc.oracle.oracles[static_cast<size_t>(oi)];
            if (O.cls != F.cls) continue;
            if (!sameSurface(O.S, F.S, O, doc.mesh)) continue;
            if (overlapFace(O, F) <= 0.0) continue;
            hit.push_back(oi);
        }
        if (hit.size() >= 2) {
            for (int oi : hit) spanned[static_cast<size_t>(oi)] = 1;
            continue;
        }
        if (hit.size() != 1 || F.face.IsNull()) continue;
        const int primary = hit[0];
        const Oracle& O0 = doc.oracle.oracles[static_cast<size_t>(primary)];
        // Enlarged face: area(F) exceeds the matched oracle by more than
        // areaQ. Exact recovered faces are within areaQ and must not span.
        if (!(F.area > O0.w + areaQ(doc.mesh, O0.tris))) continue;
        Handle(Geom_Surface) surf = BRep_Tool::Surface(F.face);
        if (surf.IsNull()) continue;
        BRepTopAdaptor_FClass2d cls(F.face, Precision::PConfusion());
        for (int oi = 0; oi < static_cast<int>(doc.oracle.oracles.size()); ++oi) {
            if (oi == primary) continue;
            const Oracle& Op = doc.oracle.oracles[static_cast<size_t>(oi)];
            if (Op.cls != F.cls) continue;
            if (sameSurface(Op.S, F.S, Op, doc.mesh)) continue;
            if (!shareMeshEdge(doc.mesh, O0, Op)) continue;
            // Cube overhang is orthogonal (S01 top vs +Y). Shallow dihedrals
            // project a neighbour centroid into F and must not span.
            double th = 0;
            for (int t : O0.tris) th = std::max(th, doc.mesh.tris[static_cast<size_t>(t)].thetaQ);
            for (int t : Op.tris) th = std::max(th, doc.mesh.tris[static_cast<size_t>(t)].thetaQ);
            if (std::fabs(dot(O0.S.n, Op.S.n)) > std::sin(th)) continue;
            bool inside = false;
            for (int t : Op.tris) {
                const Vec3& c = doc.mesh.tris[static_cast<size_t>(t)].centroid;
                GeomAPI_ProjectPointOnSurf proj(gp_Pnt(c.x, c.y, c.z), surf);
                if (proj.NbPoints() < 1) continue;
                Standard_Real u = 0, v = 0;
                proj.LowerDistanceParameters(u, v);
                if (cls.Perform(gp_Pnt2d(u, v)) == TopAbs_IN) {
                    inside = true;
                    break;
                }
            }
            if (inside) {
                spanned[static_cast<size_t>(primary)] = 1;
                spanned[static_cast<size_t>(oi)] = 1;
            }
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
        const auto& members = facesOf[static_cast<size_t>(oi)];
        if (spanned[static_cast<size_t>(oi)]) {
            f.status = Status::Spanned;
            f.credit = 0;
        } else if (!members.empty()) {
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
                // Hole tessellation: a correct plane face with circular holes
                // has GProp area below the polygonal mesh. If every wire
                // vertex of the face is a vertex of this oracle, the face
                // is the whole design surface, not a sliver (case 6 uses
                // new vertices and stays sliver).
                const StepFace& F0 = doc.step.faces[static_cast<size_t>(members[0])];
                std::unordered_set<int> ov(f.oracle.verts.begin(), f.oracle.verts.end());
                bool covers = !F0.meshVerts.empty();
                for (int mid : F0.meshVerts) {
                    if (mid < 0 || !ov.count(mid)) {
                        covers = false;
                        break;
                    }
                }
                if (covers) {
                    f.status = Status::Recovered;
                    f.credit = 1.0;
                } else {
                    f.status = Status::Sliver;
                    f.credit = (f.oracle.w > 0.0) ? (sumA / f.oracle.w) : 0.0;
                }
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
