#pragma once

#include "processors/EPT_Scalers.h"
#include "processors/Beampolmon.h"
#include "processors/ExpTrigger.h"

#include <memory>

namespace ant {
namespace analysis {
namespace slowcontrol {

// keep list of available processors here, to be used by SlowControlVariables
struct Processors {

    static const std::shared_ptr<processor::Beampolmon>  Beampolmon;
    static const std::shared_ptr<processor::EPT_Scalers> EPT_Scalers;
    static const std::shared_ptr<processor::EPT_Or>      EPT_Or;
    static const std::shared_ptr<processor::ExpTrigger>  ExpTrigger;

};

}}}
