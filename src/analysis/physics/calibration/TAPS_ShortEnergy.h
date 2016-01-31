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

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics