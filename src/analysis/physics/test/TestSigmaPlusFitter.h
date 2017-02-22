#pragma once

#include "physics/Physics.h"
#include "utils/fitter/SigmaPlusFitter.h"

namespace ant {
namespace analysis {
namespace physics {

class TestSigmaPlusFitter : public Physics {



public:
    TestSigmaPlusFitter(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}} // end of ant::analysis::physics
