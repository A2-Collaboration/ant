#pragma once

#include "analysis/physics/Physics.h"

#include "root-addons/cbtaps_display/TH2CB.h"

namespace ant {
namespace analysis {
namespace physics {

class CB_Energy : public Physics {

protected:
    TH2* ggIM = nullptr;
    TH2CB* h_cbdisplay = nullptr;

public:

    CB_Energy(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics