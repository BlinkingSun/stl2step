// DXF R12 reader for tests/tools only. Not linked into src/.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_DXF_READ_HPP
#define STL2STEP_DXF_READ_HPP

#include <string>
#include <vector>

namespace stl2step {
namespace dxfcheck {

enum class EntType { Line, Arc, Circle };

struct DxfEnt {
    EntType     type = EntType::Line;
    std::string layer;
    double x0 = 0, y0 = 0, x1 = 0, y1 = 0;  // LINE ends; ARC/CIRCLE endpoints
    double cx = 0, cy = 0, R = 0;
    double a0deg = 0, a1deg = 0;            // ARC only, DXF degrees
};

struct DxfLoop {
    std::vector<DxfEnt> segs;
    bool closed = false;
    bool outer = false;
};

struct DxfFile {
    std::string path;
    std::string projectName;
    int         declinedHeader = 0;
    int         nLine = 0, nArc = 0, nCircle = 0, nDeclined = 0;
    int         nBanned = 0;
    std::vector<DxfEnt>  ents;
    std::vector<DxfLoop> loops;
    std::vector<double>  radii;
};

bool readDxf(const std::string& path, DxfFile& out, std::string& err);
void reconstructLoops(DxfFile& f, double closeTol);
bool loopsClosed(const DxfFile& f, double closeTol, std::string& err);

}  // namespace dxfcheck
}  // namespace stl2step

#endif
