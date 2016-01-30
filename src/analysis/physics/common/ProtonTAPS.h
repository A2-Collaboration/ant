#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class ProtonTAPS : public Physics {
public:

    TTree* tree = nullptr;

    unsigned  b_nCB = 0;
    unsigned  b_nTAPS = 0;
    double    b_CBAvgTime  = 0.0;
    double    b_CBAvgVetoE  = 0.0;

    TCandidate b_Proton;


    ProtonTAPS(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const TEvent& event) override;
    void ShowResult() override;
};

}
}
}
