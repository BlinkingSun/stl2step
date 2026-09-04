// stl2step — the mesh file's own coordinate quantization (D-130-12).
//
// The absorb's admission certificate may not be a part-scale budget and may not
// be a statistic of some other population: it is the resolution of the number
// format the mesh was written in, at the mesh's own coordinate magnitude. That
// is a property of the FILE, not of the loaded double-precision points, so it is
// measured here — before welding, before scaling — and threaded down as
// MeshView::quantFloor.
//
// Pure <cstdio>/<string> — no OCCT, no engine header. Compiled directly into the
// unit test (tests/unit/quant_floor_test.cpp) the same way refit_cone_math.cpp is.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_STL_QUANT_HPP
#define STL2STEP_STL_QUANT_HPP

#include <string>

namespace stl2step {

// The measurement, all of it, so a caller can print what it used.
struct StlQuantFloor {
    bool   ok       = false;  // false => file unreadable / no coordinate found; q is 0
    bool   ascii    = false;  // false => binary STL (the 84 + 50*n layout)
    long long nTri  = 0;      // triangles seen
    double maxAbs   = 0.0;    // max |coordinate| over the mesh, FILE units
    double res      = 0.0;    // the number format's resolution at maxAbs, FILE units
    int    sigDigits = 0;     // ASCII only: significant decimal digits the file carries
    double q        = 0.0;    // res * sqrt(3)/2, FILE units
};

// Measures `path`. Never throws; returns ok=false with q=0 on anything it cannot
// read, so a caller that multiplies by its unit scale gets a certificate of 0 —
// max(bandResid, 0) == bandResid, i.e. the band's own residual and nothing more.
StlQuantFloor stlQuantFloor(const std::string& path);

}  // namespace stl2step

#endif  // STL2STEP_STL_QUANT_HPP
