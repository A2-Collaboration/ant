#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/Fitter.h"
#include "base/WrapTTree.h"

#include "analysis/utils/MCFakeReconstructed.h"
#include "analysis/utils/MCSmear.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class EtapSergey : public Physics {
public:

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(double,   TrueZVertex)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)



        ADD_BRANCH_T(double,   ProtonTime)
        ADD_BRANCH_T(double,   ProtonE)
        ADD_BRANCH_T(double,   ProtonTheta)
        ADD_BRANCH_T(double,   ProtonVetoE)
        ADD_BRANCH_T(double,   ProtonShortE)

        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)
        ADD_BRANCH_T(std::vector<double>,  PhotonsE)
        ADD_BRANCH_T(std::vector<double>,  PhotonsTheta)

        ADD_BRANCH_T(double,   PhotonSum)
        ADD_BRANCH_T(double,   MissingMass)

        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(int,      KinFitIterations)

        ADD_BRANCH_T(double,   FittedProtonE)
        ADD_BRANCH_T(double,   FittedProtonTheta)

        ADD_BRANCH_T(std::vector<double>,  FittedPhotonsE)
        ADD_BRANCH_T(std::vector<double>,  FittedPhotonsTheta)
        ADD_BRANCH_T(double,   FittedPhotonSum)

        ADD_BRANCH_T(double,   FittedZVertex)
    };

protected:

    PromptRandom::Switch promptrandom;

    TH1D* steps;

    Tree_t t;

    std::unique_ptr<utils::KinFitter> fitter;

public:
    EtapSergey(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}
