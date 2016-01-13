#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"
#include "plot/PromptRandomHist.h"
#include "utils/KinFitter.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class EtapProton : public Physics {
public:

    TTree* tree = nullptr;

    unsigned  b_nCB = 0;
    unsigned  b_nTAPS = 0;
    double    b_CBAvgTime  = 0.0;
    double    b_CBAvgVetoE  = 0.0;

    TCandidate b_Proton;
    std::vector<TCandidate> b_Photons;
    TLorentzVector b_PhotonSum;
    double b_ProtonCopl;

    PromptRandom::Switch promptrandom;

    std::vector<std::unique_ptr<utils::KinFitter>> fitters;
    double b_BestChi2;
    double b_TaggW;
    double b_TaggE;
    double b_TaggT;
    unsigned b_TaggCh;

    TLorentzVector b_FittedProton;
    std::vector<TLorentzVector> b_FittedPhotons;
    TLorentzVector b_FittedPhotonSum;
    double b_FittedProtonCopl;

    EtapProton(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
