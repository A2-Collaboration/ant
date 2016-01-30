#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class TAPSVeto_Energy : public Physics {

protected:
    TH2D* h_pedestals = nullptr;
    TH3D* h_bananas = nullptr;
public:

    TAPSVeto_Energy(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const TEvent& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics