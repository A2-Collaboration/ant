#pragma once

#include <iostream>

#include "TH1D.h"

#include "base/ParticleType.h"
#include "base/WrapTFile.h"

namespace ant {
namespace analysis {
namespace utils {

/**
 * @brief Corrections to calocluster energies,
 * allows access to True - Rec vs Rec values from an existing TH1 histogram
 *
 * E.g. the TH1 histogram can be created using https://github.com/lheijken/ClusterECorr
 * which takes as input the TH2 histogram created by analysis/physics/tunings/TrueRecCheck_ClusterE
 */

class ClusterECorr_simple {

public:
    ClusterECorr_simple();
    ~ClusterECorr_simple();
    void LoadECorr(const std::string& filename, const std::string &hname);
    double GetECorr(const double CluEin);

private:
    TH1D hECorr;
    bool hECorrSet = false;
};

}}} // namespace ant::analysis::utils
