#pragma once

#include "slowcontrol/processors/Processor.h"

#include <memory>
#include <list>

namespace ant {
namespace analysis {

class SlowControlManager;

namespace slowcontrol {

struct Variable {

    // to be called by physics classes
    virtual void Request() {
        requested = true;
    }

    virtual ~Variable() = default;

protected:
    bool requested = false;

    using ProcessorPtr = std::shared_ptr<Processor>;
    virtual std::list<ProcessorPtr> GetNeededProcessors() =0;

    friend class ant::analysis::SlowControlManager;
    std::list<ProcessorPtr> GetProcessors() {
        if(requested)
            return GetNeededProcessors();
        else
            return {};
    }
};

using VariablePtr = std::shared_ptr<Variable>;

}}}