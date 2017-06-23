#pragma once

#include "variables/Variable.h"
#include "variables/TaggerScalers.h"
#include "variables/Clocks.h"
#include "variables/Beam.h"
#include "variables/Pairspec.h"
#include "variables/Trigger.h"
#include "variables/TaggEff.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {

struct AddToAll;

// keep list of available variables here, to be used by physics classes
struct Variables {

    static const std::shared_ptr<const variable::TaggerScalers> TaggerScalers;
    static const std::shared_ptr<const variable::Clocks>        Clocks;
    static const std::shared_ptr<const variable::Beam>          Beam;
    static const std::shared_ptr<const variable::Pairspec>      Pairspec;
    static const std::shared_ptr<const variable::Trigger>       Trigger;
    static const std::shared_ptr<const variable::TaggEff>       TaggEff;



protected:
    // holds a list of shared_ptr of the static variables below
    // see SlowControlVariables.cc how this is filled using AddToAll
    // the list is used by the SlowControlManager to check which
    // variables were requested
    friend class ant::analysis::SlowControlManager;
    friend struct ant::analysis::slowcontrol::AddToAll;
    static std::list<VariablePtr> All;
};

}}}
