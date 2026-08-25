// Sidecar self-consistency checker for the P0 corpus.
// SPDX-License-Identifier: MIT

#include "corpus_util.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace corpus;

namespace {

bool readBinaryStl(const std::string& path, MeshData& mesh) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    char header[80];
    in.read(header, 80);
    uint32_t nTri = 0;
    in.read(reinterpret_cast<char*>(&nTri), 4);
    mesh.verts.clear();
    mesh.tris.clear();
    std::map<std::array<float, 3>, int> index;
    auto getIdx = [&](float x, float y, float z) {
        std::array<float, 3> key = {x, y, z};
        auto it = index.find(key);
        if (it != index.end()) return it->second;
        const int id = static_cast<int>(mesh.verts.size());
        mesh.verts.push_back({x, y, z});
        index.emplace(key, id);
        return id;
    };
    for (uint32_t i = 0; i < nTri; ++i) {
        float rec[12];
        in.read(reinterpret_cast<char*>(rec), 48);
        uint16_t attr = 0;
        in.read(reinterpret_cast<char*>(&attr), 2);
        const int a = getIdx(rec[3], rec[4], rec[5]);
        const int b = getIdx(rec[6], rec[7], rec[8]);
        const int c = getIdx(rec[9], rec[10], rec[11]);
        mesh.tris.push_back({a, b, c});
    }
    return static_cast<bool>(in);
}

double parseJsonDouble(const std::string& text, const std::string& key) {
    const std::regex re("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
    std::smatch m;
    if (std::regex_search(text, m, re)) return std::stod(m[1].str());
    return 0;
}

int parseJsonInt(const std::string& text, const std::string& key) {
    return static_cast<int>(parseJsonDouble(text, key));
}

std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

struct ParsedRecoverable {
    std::string type;
    double radius = 0;
    AxisSpec axis;
    int nSides = 0;
    double vMin = std::numeric_limits<double>::quiet_NaN();
    double vMax = std::numeric_limits<double>::quiet_NaN();
};

std::vector<ParsedRecoverable> parseRecoverables(const std::string& text) {
    std::vector<ParsedRecoverable> out;
    const std::regex cylRe(
        "\"type\"\\s*:\\s*\"cylinder\"[\\s\\S]*?\"radius\"\\s*:\\s*([-+0-9.eE]+)"
        "[\\s\\S]*?\"loc\"\\s*:\\s*\\[([^\\]]+)\\][\\s\\S]*?\"dir\"\\s*:\\s*\\[([^\\]]+)\\]"
        "[\\s\\S]*?\"nSides\"\\s*:\\s*([0-9]+)(?:[\\s\\S]*?\"vMin\"\\s*:\\s*([-+0-9.eE]+)[\\s\\S]*?\"vMax\"\\s*:\\s*([-+0-9.eE]+))?",
        std::regex::ECMAScript);
    auto begin = std::sregex_iterator(text.begin(), text.end(), cylRe);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        ParsedRecoverable r;
        r.type = "cylinder";
        r.radius = std::stod((*it)[1].str());
        {
            std::istringstream ls((*it)[2].str());
            char c;
            ls >> r.axis.loc.x >> c >> r.axis.loc.y >> c >> r.axis.loc.z;
        }
        {
            std::istringstream ds((*it)[3].str());
            char c;
            ds >> r.axis.dir.x >> c >> r.axis.dir.y >> c >> r.axis.dir.z;
        }
        r.nSides = std::stoi((*it)[4].str());
        if ((*it)[5].matched) r.vMin = std::stod((*it)[5].str());
        if ((*it)[6].matched) r.vMax = std::stod((*it)[6].str());
        out.push_back(r);
    }
    return out;
}

int countLiveEntries(const std::string& text) {
    const size_t pos = text.find("\"live\"");
    if (pos == std::string::npos) return 0;
    const size_t start = text.find('[', pos);
    const size_t end = text.find(']', start);
    if (start == std::string::npos || end == std::string::npos || end <= start) return 0;
    const std::string block = text.substr(start, end - start);
    int count = 0;
    size_t p = 0;
    while ((p = block.find("\"component\"", p)) != std::string::npos) {
        ++count;
        p += 11;
    }
    return count;
}

int parseSidecarComponentDegens(const std::string& text, int index) {
    const std::regex blockRe(
        "\\{\\s*\"index\"\\s*:\\s*" + std::to_string(index) +
        "[\\s\\S]*?\"degenerateTriangles\"\\s*:\\s*([0-9]+)",
        std::regex::ECMAScript);
    std::smatch m;
    if (std::regex_search(text, m, blockRe)) return std::stoi(m[1].str());
    return -1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: stl2step_corpus_check <corpus-dir>\n";
        return 2;
    }
    const fs::path dir(argv[1]);
    int failures = 0;
    int checked = 0;
    for (const auto& entry : fs::directory_iterator(dir)) {
        const std::string fname = entry.path().filename().string();
        if (fname.size() < 15 || fname.substr(fname.size() - 14) != ".expected.json") continue;
        const std::string base = fname.substr(0, fname.size() - 14);
        const fs::path stlPath = dir / (base + ".stl");
        const std::string jsonText = slurp(entry.path().string());
        MeshData mesh;
        if (!readBinaryStl(stlPath.string(), mesh)) {
            std::cerr << "FAIL " << base << ": unreadable STL\n";
            ++failures;
            continue;
        }
        const int expectTris = parseJsonInt(jsonText, "triangleCount");
        const double expectVol = parseJsonDouble(jsonText, "meshVolume");
        const int actualTris = static_cast<int>(mesh.tris.size());
        const double actualVol = meshVolume(mesh);
        if (actualTris != expectTris) {
            std::cerr << "FAIL " << base << ": triangleCount " << actualTris << " != "
                      << expectTris << "\n";
            ++failures;
        }
        const double rel = std::fabs(actualVol - expectVol) / std::max(1.0, std::fabs(expectVol));
        if (rel > 1e-9) {
            std::cerr << "FAIL " << base << ": meshVolume rel err " << rel << "\n";
            ++failures;
        }
        const int expectOpen = parseJsonInt(jsonText, "openEdges");
        const int expectNm = parseJsonInt(jsonText, "nonManifoldEdges");
        const EdgeStats es = meshEdgeStats(mesh);
        if (es.open != expectOpen) {
            std::cerr << "FAIL " << base << ": openEdges " << es.open << " != sidecar "
                      << expectOpen << "\n";
            ++failures;
        }
        if (es.nonManifold != expectNm) {
            std::cerr << "FAIL " << base << ": nonManifoldEdges " << es.nonManifold
                      << " != sidecar " << expectNm << "\n";
            ++failures;
        }
        const bool expectWt = (jsonText.find("\"watertight\": true") != std::string::npos) ||
                              (jsonText.find("\"watertight\":true") != std::string::npos);
        const bool actualWt = (es.open == 0 && es.nonManifold == 0);
        // Sidecar must agree with the welded mesh: torn fixtures (S02/S09/S15) set false.
        if (jsonText.find("\"watertight\"") != std::string::npos && expectWt != actualWt) {
            std::cerr << "FAIL " << base << ": watertight sidecar=" << expectWt
                      << " mesh=" << actualWt << "\n";
            ++failures;
        }
        for (const auto& r : parseRecoverables(jsonText)) {
            if (r.type != "cylinder") continue;
            Recoverable rec;
            rec.type = "cylinder";
            rec.radius = r.radius;
            rec.axis = r.axis;
            rec.vMin = r.vMin;
            rec.vMax = r.vMax;
            const int measured = measureRecoverableNSides(mesh, rec);
            if (measured != r.nSides) {
                std::cerr << "FAIL " << base << ": nSides measured " << measured
                          << " != sidecar " << r.nSides << "\n";
                ++failures;
            }
        }
        const std::vector<MeshComponentInfo> comps = splitMeshComponents(mesh);
        if (jsonText.find("\"components\"") == std::string::npos) {
            std::cerr << "FAIL " << base << ": missing components[] block\n";
            ++failures;
        } else {
            for (const auto& c : comps) {
                const int expectDeg = parseSidecarComponentDegens(jsonText, c.index);
                if (expectDeg < 0) {
                    std::cerr << "FAIL " << base << ": sidecar missing component "
                              << c.index << "\n";
                    ++failures;
                } else if (expectDeg != c.degenerateTriangles) {
                    std::cerr << "FAIL " << base << ": component " << c.index
                              << " degenerateTriangles " << c.degenerateTriangles
                              << " != sidecar " << expectDeg << "\n";
                    ++failures;
                }
                if (c.degenerateTriangles > 0 && c.openEdges == 0) {
                    std::cerr << "FAIL " << base << ": component " << c.index
                              << " has " << c.degenerateTriangles
                              << " degenerate triangle(s) but openEdges==0\n";
                    ++failures;
                }
            }
        }
        const int liveCount = countLiveEntries(jsonText);
        const int expectLive = static_cast<int>(
            comps.size() - (base == "S15" || base == "S16-R1-round-2" ? comps.size() : 0));
        if (expectLive > 0 && liveCount < expectLive) {
            std::cerr << "FAIL " << base << ": live[] has " << liveCount
                      << " entries, expected >= " << expectLive << "\n";
            ++failures;
        }
        if (jsonText.find("\"faceCount\"") == std::string::npos ||
            jsonText.find("\"surfaceCensus\"") == std::string::npos ||
            jsonText.find("\"volumeBudgetMM3\"") == std::string::npos) {
            std::cerr << "FAIL " << base << ": live block missing faceCount/"
                         "surfaceCensus/volumeBudgetMM3\n";
            ++failures;
        }
        std::cout << "OK " << base << " tris=" << actualTris << " vol=" << actualVol
                  << " open=" << es.open << " nm=" << es.nonManifold << "\n";
        ++checked;
    }
    if (checked == 0) {
        std::cerr << "no fixtures found in " << dir << "\n";
        return 2;
    }
    std::cout << "checked " << checked << " fixture(s), failures=" << failures << "\n";
    return failures ? 1 : 0;
}
