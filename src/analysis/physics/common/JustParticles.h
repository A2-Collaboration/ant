#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"
#include "plot/PromptRandomHist.h"
#include "utils/Fitter.h"
#include "expconfig/detectors/TAPS.h"
#include "base/WrapTTree.h"

#include "TLorentzVector.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class JustParticles : public Physics {

protected:
    const PiecewiseInterval<unsigned> multiplicities;
    const bool save_events;

    TH1D* steps;

    TTree* tree = nullptr;

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(unsigned, b_nCB)
        ADD_BRANCH_T(unsigned,  b_nTAPS)
        ADD_BRANCH_T(double,    b_CBAvgTime)
        ADD_BRANCH_T(double,    b_CBSumVetoE)

        ADD_BRANCH_T(TLorentzVector, b_PhotonSum)
        ADD_BRANCH_T(TLorentzVector, b_Proton)
        ADD_BRANCH_T(double, b_Proton_vetoE)
        ADD_BRANCH_T(double, b_ProtonCopl)
        ADD_BRANCH_T(double, b_ProtonBeta)
        ADD_BRANCH_T(double, b_ProtonToF)
        ADD_BRANCH_T(double, b_ProtonPSA_R)
        ADD_BRANCH_T(double, b_ProtonPSA_Angle)

        ADD_BRANCH_T(double, b_FitChi2)
        ADD_BRANCH_T(unsigned, b_FitStatus)
        ADD_BRANCH_T(unsigned, b_NFitIterations)
        ADD_BRANCH_T(double, b_TaggW)
        ADD_BRANCH_T(double, b_TaggE)
        ADD_BRANCH_T(double, b_TaggT)
        ADD_BRANCH_T(unsigned, b_TaggCh)
        ADD_BRANCH_T(TLorentzVector, b_Missing)

        ADD_BRANCH_T(double, b_FittedTaggE)
        ADD_BRANCH_T(TLorentzVector, b_FittedProton)
        ADD_BRANCH_T(TLorentzVector, b_FittedPhotonSum)
        ADD_BRANCH_T(double, b_FittedProtonCopl)
    };

    Tree_t t;

    PromptRandom::Switch promptrandom;
    std::vector<std::unique_ptr<utils::KinFitter>> fitters;
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;

public:

    JustParticles(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
