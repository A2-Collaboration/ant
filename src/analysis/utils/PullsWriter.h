#pragma once

#include "base/WrapTTree.h"
#include "analysis/utils/Fitter.h"

namespace ant {
namespace analysis {

class HistogramFactory;

namespace utils {

class PullsWriter {

public:
    struct PullTree_t : WrapTTree {

        ADD_BRANCH_T(double, TaggW) // for prompt-random subtraction
        ADD_BRANCH_T(double, FitProb)
        ADD_BRANCH_T(double, FittedZVertex)
        ADD_BRANCH_T(unsigned, Multiplicity)

        // for Particle ID debugging (time of flight),
        // and coverage testing
        ADD_BRANCH_T(double, ProtonE)
        ADD_BRANCH_T(double, ProtonTheta)
        ADD_BRANCH_T(double, ProtonTime)

        ADD_BRANCH_T(double, E)
        ADD_BRANCH_T(double, Theta)
        ADD_BRANCH_T(double, Phi)
        ADD_BRANCH_T(double, ShowerDepth)

        // see utils::Fitter::FitParticle for the parametrization
        ADD_BRANCH_T(std::vector<double>, Values)
        ADD_BRANCH_T(std::vector<double>, Sigmas)
        ADD_BRANCH_T(std::vector<double>, Pulls)
    };

protected:
    PullTree_t photons_cb;
    PullTree_t photons_taps;

    PullTree_t proton_cb;
    PullTree_t proton_taps;

    PullTree_t& getPullTree(const ant::analysis::utils::Fitter::FitParticle& particle);

public:

    PullsWriter(ant::analysis::HistogramFactory& histfac);
    ~PullsWriter();

    void Fill(const std::vector<Fitter::FitParticle>& fitParticles,
              double tagger_weight, double fitprob, double fitted_z_vertex);

};

}
}
}
