#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class TAPS_ShortEnergy : public Physics {

protected:

    TH2D* h_pedestals = nullptr;
    TH2D* h_rel_gamma = nullptr;

public:

    TAPS_ShortEnergy(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics