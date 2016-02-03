#pragma once

#include "Processor.h"

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct EPT_Scalers : Processor {


    virtual return_t ProcessEventData(const TEventData& recon,  physics::manager_t& manager) override
    {
        return return_t::Complete;
    }

};

}}}} // namespace ant::analysis::slowcontrol::processor
