#pragma once

#include "analysis/physics/Physics.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class ReconstructCheck : public Physics {
protected:

    TH2D* EnergyRec_cb;
    TH2D* EnergyRec_taps;

public:
    ReconstructCheck(PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
