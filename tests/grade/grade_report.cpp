#include "grade_report.hpp"

#include <cmath>
#include <cstdio>
#include <sstream>

namespace grade {
namespace {

std::string esc(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') {
            o += '\\';
            o += static_cast<char>(c);
        } else if (c >= 0x20) {
            o += static_cast<char>(c);
        }
    }
    return o;
}

std::string jnum(double v) {
    if (!std::isfinite(v)) return "null";
    return fmtG(v);
}

std::string jbool(bool v) { return v ? "true" : "false"; }

std::string jnull(double v, bool has) { return has ? jnum(v) : "null"; }

std::string jvec(const Vec3& v) {
    return std::string("[") + jnum(v.x) + "," + jnum(v.y) + "," + jnum(v.z) + "]";
}

void countsJson(std::ostringstream& o, const int c[7]) {
    o << "{\"recovered\":" << c[0] << ",\"split\":" << c[1] << ",\"sliver\":" << c[2]
      << ",\"faceted\":" << c[3] << ",\"facetedPartial\":" << c[4] << ",\"spanned\":" << c[5]
      << ",\"missing\":" << c[6] << "}";
}

void byClassOne(std::ostringstream& o, const char* name, int n, double area, double rec,
                const int c[7]) {
    o << "\"" << name << "\":{\"oracles\":" << n << ",\"areaMM2\":" << jnum(area)
      << ",\"recoveredAreaMM2\":" << jnum(rec) << ",\"counts\":";
    countsJson(o, c);
    o << "}";
}

std::string paramsJson(const Oracle& o) {
    std::ostringstream s;
    s << "{";
    switch (o.cls) {
        case SurfClass::Plane:
            s << "\"normal\":" << jvec(o.S.n) << ",\"origin\":" << jvec(o.S.p0);
            break;
        case SurfClass::Cylinder:
            s << "\"radius\":" << jnum(o.S.R) << ",\"axisDir\":" << jvec(o.S.n)
              << ",\"axisLoc\":" << jvec(o.S.p0);
            break;
        case SurfClass::Cone:
            s << "\"alpha\":" << jnum(o.S.alpha) << ",\"axisDir\":" << jvec(o.S.n)
              << ",\"apex\":" << jvec(o.S.apex);
            break;
        case SurfClass::Sphere:
            s << "\"radius\":" << jnum(o.S.R) << ",\"centre\":" << jvec(o.S.p0);
            break;
        case SurfClass::Torus:
            s << "\"Rmaj\":" << jnum(o.S.R) << ",\"Rmin\":" << jnum(o.S.r)
              << ",\"axisDir\":" << jvec(o.S.n) << ",\"centre\":" << jvec(o.S.p0);
            break;
        default:
            break;
    }
    s << "}";
    return s.str();
}

std::string statusPrint(const Feature& f) {
    if (f.status == Status::Split) return std::string("split(") + std::to_string(f.splitN) + ")";
    return statusName(f.status);
}

}  // namespace

std::string writeJson(const GradeDocument& d) {
    std::ostringstream o;
    const std::string stlRel = relPath(d.stlPath, d.stlPath);
    const std::string stepRel = relPath(d.stlPath, d.stepPath);
    const char* fmt = d.mesh.qf.ascii ? "ascii" : "binary";
    o << "{\n";
    o << "\"schema\":\"stl2step.grade/1\",\n";
    o << "\"stl\":{\"path\":\"" << esc(stlRel.empty() ? d.stlPath : stlRel) << "\",\"sha256\":\""
      << d.stlSha << "\",\"bytes\":" << d.stlBytes << "},\n";
    o << "\"step\":{\"path\":\"" << esc(stepRel) << "\",\"sha256\":\"" << d.stepSha
      << "\",\"bytes\":" << d.stepBytes << "},\n";
    o << "\"q\":{\"q\":" << jnum(d.mesh.qf.q) << ",\"res\":" << jnum(d.mesh.qf.res)
      << ",\"maxAbs\":" << jnum(d.mesh.qf.maxAbs) << ",\"format\":\"" << fmt
      << "\",\"sigDigits\":" << d.mesh.qf.sigDigits << "},\n";
    o << "\"mesh\":{\"triangles\":" << d.mesh.tris.size()
      << ",\"weldedVertices\":" << d.mesh.verts.size() << ",\"openEdges\":" << d.mesh.openEdges
      << ",\"nonManifoldEdges\":" << d.mesh.nonManifoldEdges
      << ",\"degenerateDropped\":" << d.mesh.degenerateDropped
      << ",\"surfaceAreaMM2\":" << jnum(d.mesh.surfaceArea)
      << ",\"volumeMM3\":" << jnum(d.mesh.volume) << "},\n";
    o << "\"step_census\":{\"faces\":" << d.step.faces.size() << ",\"facetFaces\":" << d.step.nFacet
      << ",\"planes\":" << d.step.nPlane << ",\"cylinders\":" << d.step.nCyl
      << ",\"cones\":" << d.step.nCone << ",\"spheres\":" << d.step.nSphere
      << ",\"tori\":" << d.step.nTorus << ",\"bsplines\":" << d.step.nBSpline
      << ",\"other\":" << d.step.nOther << "},\n";
    o << "\"hardZero\":[";
    for (size_t i = 0; i < d.hardZero.size(); ++i) {
        if (i) o << ",";
        o << "\"" << d.hardZero[i] << "\"";
    }
    o << "],\n";
    o << "\"checks\":{\"watertight\":" << jbool(d.watertight) << ",\"valid\":" << jbool(d.valid)
      << ",\"determinism\":\"" << d.determinism << "\",\"volume\":{\"stepMM3\":" << jnum(d.Vstep)
      << ",\"meshMM3\":" << jnum(d.Vmesh) << ",\"deltaMM3\":" << jnum(d.volumeDelta)
      << ",\"chordBudgetMM3\":" << jnum(d.chordBudget) << ",\"volumeQMM3\":" << jnum(d.volumeQ)
      << ",\"within\":" << jbool(d.volumeWithin) << "}},\n";
    o << "\"grade\":{\"plane\":" << jnull(d.gradePlane, d.hasPlane)
      << ",\"cylinder\":" << jnull(d.gradeCyl, d.hasCyl) << ",\"cone\":" << jnull(d.gradeCone, d.hasCone)
      << ",\"sphere\":" << jnull(d.gradeSphere, d.hasSphere)
      << ",\"torus\":" << jnull(d.gradeTorus, d.hasTorus)
      << ",\"overall\":" << jnull(d.gradeOverall, d.hasOverall) << "},\n";
    o << "\"byClass\":{";
    byClassOne(o, "plane", d.nPlane, d.areaPlane, d.recPlane, d.cntPlane);
    o << ",";
    byClassOne(o, "cylinder", d.nCyl, d.areaCyl, d.recCyl, d.cntCyl);
    o << ",";
    byClassOne(o, "cone", d.nCone, d.areaCone, d.recCone, d.cntCone);
    o << ",";
    byClassOne(o, "sphere", d.nSphere, d.areaSphere, d.recSphere, d.cntSphere);
    o << ",";
    byClassOne(o, "torus", d.nTorus, d.areaTorus, d.recTorus, d.cntTorus);
    o << "},\n";
    const double frac = (d.mesh.surfaceArea > 0.0) ? (d.oracle.residueArea / d.mesh.surfaceArea) : 0.0;
    o << "\"residue\":{\"triangles\":" << d.oracle.residueTris
      << ",\"areaMM2\":" << jnum(d.oracle.residueArea) << ",\"areaFraction\":" << jnum(frac)
      << ",\"facetedFraction\":" << jnum(d.residueFacetedFraction)
      << ",\"unprovableSingletons\":" << d.oracle.unprovableSingletons << ",\"components\":[";
    for (size_t i = 0; i < d.oracle.residue.size(); ++i) {
        const ResidueComp& c = d.oracle.residue[i];
        if (i) o << ",";
        o << "{\"triangles\":" << c.tris.size() << ",\"areaMM2\":" << jnum(c.area)
          << ",\"ruled\":" << jbool(c.ruled) << ",\"direction\":" << jvec(c.direction)
          << ",\"bbox\":[" << jvec(c.bboxMin) << "," << jvec(c.bboxMax) << "]}";
    }
    o << "]},\n";
    o << "\"features\":[";
    for (size_t i = 0; i < d.features.size(); ++i) {
        const Feature& f = d.features[i];
        if (i) o << ",";
        o << "{\"featureId\":\"" << esc(f.oracle.featureId) << "\",\"class\":\""
          << className(f.oracle.cls) << "\",\"status\":\"" << statusName(f.status) << "\"";
        if (f.status == Status::Split) o << ",\"splitN\":" << f.splitN;
        o << ",\"credit\":" << jnum(f.credit) << ",\"areaMM2\":" << jnum(f.oracle.w)
          << ",\"areaFraction\":" << jnum(f.areaFraction)
          << ",\"triangles\":" << f.oracle.tris.size()
          << ",\"distinctVertices\":" << f.oracle.verts.size()
          << ",\"oraclePieces\":" << f.oracle.oraclePieces
          << ",\"oracleMinSeparation\":" << f.oracle.oracleMinSeparation
          << ",\"maxResidMM\":" << jnum(f.oracle.maxResid)
          << ",\"maxNormalDevRad\":" << jnum(f.oracle.maxNormalDev)
          << ",\"coverage\":" << jnum(f.coverage) << ",\"params\":" << paramsJson(f.oracle)
          << ",\"stepFaces\":[";
        for (size_t k = 0; k < f.stepFaces.size(); ++k) {
            if (k) o << ",";
            o << "{\"entity\":" << f.stepFaces[k].entity
              << ",\"areaMM2\":" << jnum(f.stepFaces[k].areaMM2) << "}";
        }
        o << "]}";
    }
    o << "],\n";
    o << "\"intersections\":[";
    for (size_t i = 0; i < d.intersections.size(); ++i) {
        const Intersection& x = d.intersections[i];
        if (i) o << ",";
        o << "{\"a\":\"" << esc(x.a) << "\",\"b\":\"" << esc(x.b)
          << "\",\"expectedTier\":" << x.expectedTier << ",\"shipped\":\"" << esc(x.shipped)
          << "\",\"verdict\":\"" << esc(x.verdict) << "\",\"edgeEntity\":" << x.edgeEntity << "}";
    }
    o << "],\n";
    o << "\"unmatchedStepFaces\":[";
    for (size_t i = 0; i < d.unmatched.size(); ++i) {
        if (i) o << ",";
        o << "{\"entity\":" << d.unmatched[i].entity << ",\"type\":\"" << esc(d.unmatched[i].type)
          << "\",\"areaMM2\":" << jnum(d.unmatched[i].areaMM2) << "}";
    }
    o << "]\n}\n";
    return o.str();
}

std::string writeMd(const GradeDocument& d) {
    std::ostringstream o;
    const std::string stlRel = relPath(d.stlPath, d.stlPath);
    auto stem = [](std::string p) {
        for (char& c : p)
            if (c == '\\') c = '/';
        const auto sl = p.find_last_of('/');
        if (sl != std::string::npos) p = p.substr(sl + 1);
        return p;
    };
    o << "# grade — " << stem(d.stlPath) << " × " << stem(d.stepPath) << "\n";
    o << "q = " << fmtG(d.mesh.qf.q) << " mm (" << (d.mesh.qf.ascii ? "ascii" : "binary")
      << ", ulp " << fmtG(d.mesh.qf.res) << " @ maxAbs " << fmtG(d.mesh.qf.maxAbs) << ")\n";
    o << "mesh: " << d.mesh.tris.size() << " triangles, " << d.mesh.verts.size()
      << " welded vertices, openEdges " << d.mesh.openEdges << ", nonManifoldEdges "
      << d.mesh.nonManifoldEdges << "\n";
    o << "step: " << d.step.faces.size() << " faces (" << d.step.nFacet << " facet), "
      << (d.valid ? "valid" : "invalid") << ", " << (d.watertight ? "watertight" : "open") << "\n";
    o << "volume: step " << fmtG(d.Vstep) << " mesh " << fmtG(d.Vmesh) << " delta "
      << fmtG(d.volumeDelta) << " chordBudget " << fmtG(d.chordBudget) << " volumeQ "
      << fmtG(d.volumeQ) << "  -> " << (d.volumeWithin ? "within" : "OUTSIDE") << "\n";
    o << "determinism: " << d.determinism
      << (d.determinism == "skipped" ? " (no --engine-bin)" : "") << "\n\n";
    o << "| class    | oracles | area mm^2 | recovered mm^2 | grade |\n";
    o << "|----------|--------:|----------:|---------------:|------:|\n";
    auto row = [&](const char* name, int n, double area, double rec, bool has, double g) {
        o << "| " << name << " | " << n << " | " << fmtG(area) << " | " << fmtG(rec) << " | ";
        if (has) o << fmtG(g);
        else o << "-";
        o << " |\n";
    };
    row("plane    ", d.nPlane, d.areaPlane, d.recPlane, d.hasPlane, d.gradePlane);
    row("cylinder ", d.nCyl, d.areaCyl, d.recCyl, d.hasCyl, d.gradeCyl);
    row("cone     ", d.nCone, d.areaCone, d.recCone, d.hasCone, d.gradeCone);
    row("sphere   ", d.nSphere, d.areaSphere, d.recSphere, d.hasSphere, d.gradeSphere);
    row("torus    ", d.nTorus, d.areaTorus, d.recTorus, d.hasTorus, d.gradeTorus);
    row("OVERALL  ", d.nPlane + d.nCyl + d.nCone + d.nSphere + d.nTorus,
        d.areaPlane + d.areaCyl + d.areaCone + d.areaSphere + d.areaTorus,
        d.recPlane + d.recCyl + d.recCone + d.recSphere + d.recTorus, d.hasOverall,
        d.gradeOverall);
    const double frac = (d.mesh.surfaceArea > 0.0) ? (d.oracle.residueArea / d.mesh.surfaceArea) : 0.0;
    o << "\nresidue: " << d.oracle.residueTris << " triangles, " << fmtG(d.oracle.residueArea)
      << " mm^2 (" << fmtG(frac * 100.0) << " % of surface), " << d.oracle.residue.size()
      << " components, faceted " << fmtG(d.residueFacetedFraction) << "\n";
    o << "unprovable singletons: " << d.oracle.unprovableSingletons << "\n";
    if (!d.hardZero.empty()) {
        o << "hardZero:";
        for (const auto& h : d.hardZero) o << " " << h;
        o << "\n";
    }
    o << "\n## work queue (unrecovered or split, by area desc)\n";
    o << "| # | featureId              | class    | status  | area mm^2 |  frac | tris | pieces | params                          |\n";
    o << "|--:|------------------------|----------|---------|----------:|------:|-----:|-------:|---------------------------------|\n";
    int nq = 0;
    for (const Feature& f : d.features) {
        if (f.status == Status::Recovered) continue;
        ++nq;
        o << "| " << nq << " | " << f.oracle.featureId << " | " << className(f.oracle.cls) << " | "
          << statusPrint(f) << " | " << fmtG(f.oracle.w) << " | " << fmtG(f.areaFraction * 100.0)
          << "% | " << f.oracle.tris.size() << " | " << f.oracle.oraclePieces << " | ";
        switch (f.oracle.cls) {
            case SurfClass::Plane:
                o << "n " << fmtG(f.oracle.S.n.x) << "," << fmtG(f.oracle.S.n.y) << ","
                  << fmtG(f.oracle.S.n.z);
                break;
            case SurfClass::Cylinder:
                o << "R " << fmtG(f.oracle.S.R);
                break;
            case SurfClass::Cone:
                o << "alpha " << fmtG(f.oracle.S.alpha);
                break;
            case SurfClass::Sphere:
                o << "R " << fmtG(f.oracle.S.R);
                break;
            case SurfClass::Torus:
                o << "Rmaj " << fmtG(f.oracle.S.R) << " Rmin " << fmtG(f.oracle.S.r);
                break;
            default:
                break;
        }
        if (f.oracle.cls == SurfClass::Sphere)
            o << " | register: 1.5";
        o << " |\n";
        if (nq >= 5 && false) break;
    }
    o << "\n## intersections\n";
    o << "| a | b | expected | shipped | verdict |\n";
    o << "|---|---|---|---|---|\n";
    for (const Intersection& x : d.intersections) {
        o << "| " << x.a << " | " << x.b << " | tier" << x.expectedTier << " | " << x.shipped
          << " | " << x.verdict << " |\n";
    }
    return o.str();
}

bool writeFile(const std::string& path, const std::string& text) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const size_t n = std::fwrite(text.data(), 1, text.size(), f);
    std::fclose(f);
    return n == text.size();
}

}  // namespace grade
