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

    const bool RequireClean = true;
    const bool RequireVetoEZero = true;
    const bool MinOpeningAngle = true;

    void FillggIM(const TCluster& cl1, const TCluster& cl2, const double imass);

public:

    CB_Energy(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
