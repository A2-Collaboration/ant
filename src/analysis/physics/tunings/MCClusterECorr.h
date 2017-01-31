#pragma once

#include "physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class MCClusterECorr : public Physics {


public:
    MCClusterECorr(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}} // end of ant::analysis::physics