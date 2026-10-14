// Shared profile certificate (D-fc-1.10). Bodies moved verbatim from refit_torus.cpp.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_REFIT_PROFILE_CERT_HPP
#define STL2STEP_REFIT_PROFILE_CERT_HPP

#include <vector>

namespace stl2step {
namespace refit {

struct ProfPt {
    double rho = 0, z = 0, psi = 0;
    int lv = -1;
};

struct TorusCert {
    bool ok = false;
    bool determinate = false;
    const char* why = "none";
    double Rmaj = 0, z0 = 0, Rmin = 0;
    double maxDev = 0;
    double lineRes = 0;
    int levels = 0, cols = 0, profPts = 0;
    double rmsDev = 0;
};

int countColumns(std::vector<double> psi, const std::vector<double>& rho, double tau);

int countLevels(std::vector<double> z, double tau);

int countProfPts(const std::vector<ProfPt>& P, double tau);

double profileLineMaxResid(const std::vector<ProfPt>& P);

bool fitProfileCircle(const std::vector<ProfPt>& P, double tau, double& Rmaj, double& z0,
                      double& Rmin);

TorusCert certifyTorus(const std::vector<ProfPt>& P, double tau);

}}  // namespace stl2step::refit

#endif  // STL2STEP_REFIT_PROFILE_CERT_HPP
