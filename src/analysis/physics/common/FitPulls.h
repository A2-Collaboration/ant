#pragma once

#include "physics/Physics.h"
#include "utils/Fitter.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class FitPulls : public Physics {

public:
    FitPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}