#pragma once

#include "variables/Variable.h"
#include "variables/TaggerScalers.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {

// keep list of available variables here, to be used by physics classes
struct Variables {
    // holds a list of shared_ptr of the static variables below
    // see SlowControlVariables.cc how this is filled
    // the list is used by the SlowControlManager to check which
    // variables were requested
    static std::list<VariablePtr> All;

    struct AddToAll {
        AddToAll(const VariablePtr& var) {
            All.push_back(var);
        }
    };


    static const std::shared_ptr<variable::TaggerScalers> TaggerScalers;
};

}}}