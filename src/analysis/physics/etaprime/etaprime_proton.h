#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"

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
    //std::vector<TLorentzVector> b_Photons;

    TLorentzVector b_PhotonSum;
    double b_ProtonCopl;


    EtapProton(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
