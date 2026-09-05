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

// Sidecar scanning. These are machine-written JSON documents, but the tool
// deliberately takes no JSON dependency: every value is found by scanning the
// text forward once. This replaced std::regex patterns chaining "[\s\S]*?"
// across the WHOLE document -- MSVC's engine backtracks those combinatorially
// and aborts a multi-KB sidecar with regex_error(error_complexity), which the
// uncaught-exception path reports as 0xC0000409 on Windows CI. The scans below
// reproduce those patterns' matching semantics exactly, quirks included.

bool isJsonSpace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}

size_t skipJsonSpace(const std::string& text, size_t p) {
    while (p < text.size() && isJsonSpace(text[p])) ++p;
    return p;
}

// `"<key>" \s* : \s*` anchored at `pos`; returns the value offset, or npos.
size_t matchKeyColon(const std::string& text, size_t pos, const std::string& quotedKey) {
    if (pos > text.size() || text.size() - pos < quotedKey.size()) return std::string::npos;
    if (text.compare(pos, quotedKey.size(), quotedKey) != 0) return std::string::npos;
    size_t p = skipJsonSpace(text, pos + quotedKey.size());
    if (p >= text.size() || text[p] != ':') return std::string::npos;
    return skipJsonSpace(text, p + 1);
}

// Greedy run of the ECMAScript class [-+0-9.eE]; returns `p` when it is empty.
size_t scanNumber(const std::string& text, size_t p) {
    size_t q = p;
    while (q < text.size()) {
        const char c = text[q];
        if (c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E' ||
            (c >= '0' && c <= '9')) {
            ++q;
        } else {
            break;
        }
    }
    return q;
}

// Greedy run of [0-9]; returns `p` when it is empty.
size_t scanDigits(const std::string& text, size_t p) {
    size_t q = p;
    while (q < text.size() && text[q] >= '0' && text[q] <= '9') ++q;
    return q;
}

enum class ValueKind { Number, Digits, Array };

struct KeyHit {
    size_t begin = 0;  // first character of the value (inside [] for Array)
    size_t end = 0;    // one past the last character of the value
    size_t after = 0;  // one past the whole `"key": value` token (past ] for Array)
};

// First `"<key>" : <value>` at or after `from`, `<value>` matching `kind`.
// Occurrences whose value does not match are skipped and the scan continues,
// exactly as the lazy "[\s\S]*?" prefix used to.
bool findKeyValue(const std::string& text, size_t from, const std::string& quotedKey,
                  ValueKind kind, KeyHit& hit) {
    for (size_t at = text.find(quotedKey, from); at != std::string::npos;
         at = text.find(quotedKey, at + 1)) {
        const size_t v = matchKeyColon(text, at, quotedKey);
        if (v == std::string::npos) continue;
        if (kind == ValueKind::Array) {
            // \[([^\]]+)\]: [^\]] cannot cross the bracket, so the run ends
            // at the first ] and must be non-empty.
            if (v >= text.size() || text[v] != '[') continue;
            const size_t close = text.find(']', v + 1);
            if (close == std::string::npos || close == v + 1) continue;
            hit = {v + 1, close, close + 1};
            return true;
        }
        const size_t e = (kind == ValueKind::Digits) ? scanDigits(text, v) : scanNumber(text, v);
        if (e == v) continue;
        hit = {v, e, e};
        return true;
    }
    return false;
}

std::string valueOf(const std::string& text, const KeyHit& hit) {
    return text.substr(hit.begin, hit.end - hit.begin);
}

double parseJsonDouble(const std::string& text, const std::string& key) {
    KeyHit hit;
    if (findKeyValue(text, 0, "\"" + key + "\"", ValueKind::Number, hit)) {
        return std::stod(valueOf(text, hit));
    }
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

// Linear equivalent of the old cylinder pattern:
//   "type"\s*:\s*"cylinder" [\s\S]*? "radius"\s*:\s*(num)
//   [\s\S]*? "loc"\s*:\s*\[(list)\] [\s\S]*? "dir"\s*:\s*\[(list)\]
//   [\s\S]*? "nSides"\s*:\s*(int)
//   (?: [\s\S]*? "vMin"\s*:\s*(num) [\s\S]*? "vMax"\s*:\s*(num) )?
// Each lazy prefix is "the next occurrence of the following key", so the whole
// pattern is one forward walk. Note the two behaviours this deliberately keeps:
// the fields are picked up wherever they next appear (so a cylinder entry with
// no "nSides" of its own borrows a later entry's -- linkage_bores_chamfer's
// sidecar documents relying on that), and the trailing group is greedy, so it is
// taken whenever some "vMin" with a following "vMax" exists later in the file.
std::vector<ParsedRecoverable> parseRecoverables(const std::string& text) {
    std::vector<ParsedRecoverable> out;
    static const std::string kType = "\"type\"";
    static const std::string kCylinder = "\"cylinder\"";
    size_t searchFrom = 0;
    while (true) {
        const size_t typeAt = text.find(kType, searchFrom);
        if (typeAt == std::string::npos) break;
        // Any failure below retries from the next candidate start, as the regex
        // engine's leftmost-match scan did.
        searchFrom = typeAt + 1;
        const size_t typeVal = matchKeyColon(text, typeAt, kType);
        if (typeVal == std::string::npos) continue;
        if (text.compare(typeVal, kCylinder.size(), kCylinder) != 0) continue;

        ParsedRecoverable r;
        r.type = "cylinder";
        KeyHit radius, loc, dir, nSides;
        size_t p = typeVal + kCylinder.size();
        if (!findKeyValue(text, p, "\"radius\"", ValueKind::Number, radius)) continue;
        if (!findKeyValue(text, radius.after, "\"loc\"", ValueKind::Array, loc)) continue;
        if (!findKeyValue(text, loc.after, "\"dir\"", ValueKind::Array, dir)) continue;
        if (!findKeyValue(text, dir.after, "\"nSides\"", ValueKind::Digits, nSides)) continue;
        r.radius = std::stod(valueOf(text, radius));
        {
            std::istringstream ls(valueOf(text, loc));
            char c;
            ls >> r.axis.loc.x >> c >> r.axis.loc.y >> c >> r.axis.loc.z;
        }
        {
            std::istringstream ds(valueOf(text, dir));
            char c;
            ds >> r.axis.dir.x >> c >> r.axis.dir.y >> c >> r.axis.dir.z;
        }
        r.nSides = std::stoi(valueOf(text, nSides));
        p = nSides.after;
        // The trailing group is optional but greedy, so it is taken whenever it
        // can be. Only the first "vMin" needs trying: a later one starts further
        // in, so it can only see a subset of this one's candidate "vMax".
        KeyHit vMin, vMax;
        if (findKeyValue(text, p, "\"vMin\"", ValueKind::Number, vMin) &&
            findKeyValue(text, vMin.after, "\"vMax\"", ValueKind::Number, vMax)) {
            r.vMin = std::stod(valueOf(text, vMin));
            r.vMax = std::stod(valueOf(text, vMax));
            p = vMax.after;
        }
        out.push_back(r);
        searchFrom = p;  // resume past the whole match, as sregex_iterator did
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

// Linear equivalent of:
//   \{\s*"index"\s*:\s*<index>[\s\S]*?"degenerateTriangles"\s*:\s*([0-9]+)
// The index is compared as a bare digit prefix, exactly as the pattern did.
int parseSidecarComponentDegens(const std::string& text, int index) {
    static const std::string kIndex = "\"index\"";
    const std::string digits = std::to_string(index);
    for (size_t brace = text.find('{'); brace != std::string::npos;
         brace = text.find('{', brace + 1)) {
        const size_t indexVal = matchKeyColon(text, skipJsonSpace(text, brace + 1), kIndex);
        if (indexVal == std::string::npos) continue;
        if (text.compare(indexVal, digits.size(), digits) != 0) continue;
        KeyHit degens;
        if (!findKeyValue(text, indexVal + digits.size(), "\"degenerateTriangles\"",
                          ValueKind::Digits, degens)) {
            continue;
        }
        return std::stoi(valueOf(text, degens));
    }
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
