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

    virtual bool HasChanged() {
        for(auto& p : GetNeededProcessors())
            if(p->HasChanged())
               return true;
        return false;
    }

    virtual ~Variable() = default;

protected:
    friend class ant::analysis::SlowControlManager;

    bool requested = false;
    using ProcessorPtr = std::shared_ptr<Processor>;
    virtual std::list<ProcessorPtr> GetNeededProcessors() =0;
};

using VariablePtr = std::shared_ptr<Variable>;

}}}