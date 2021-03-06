#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {


class MCSmearing : public Physics {
protected:
    TH3D* energies;
    TH3D* angles;
    TH1D* IM;
    TH1D* IM_2g;


public:
    MCSmearing(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~MCSmearing();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
