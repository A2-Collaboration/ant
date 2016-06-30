#include "base/WrapTTree.h"
#include "analysis/utils/Fitter.h"

namespace ant {
namespace analysis {

class HistogramFactory;

namespace utils {

class PullOutput {

public:
    struct PullTree_t : WrapTTree {

        ADD_BRANCH_T(double, E)
        ADD_BRANCH_T(double, Theta)
        ADD_BRANCH_T(double, Phi)

        ADD_BRANCH_T(double, PullE)
        ADD_BRANCH_T(double, PullTheta)
        ADD_BRANCH_T(double, PullPhi)

        ADD_BRANCH_T(unsigned, Multiplicity)
    };

protected:
    PullTree_t photons_cb;
    PullTree_t photons_taps;

    PullTree_t proton_cb;
    PullTree_t proton_taps;

    PullTree_t& getProtonTree(const ant::analysis::utils::Fitter::FitParticle& particle);

public:

    PullOutput(ant::analysis::HistogramFactory& histfac);
    ~PullOutput();

    void Fill(const std::vector<Fitter::FitParticle>& fitParticles);

};

}
}
}
