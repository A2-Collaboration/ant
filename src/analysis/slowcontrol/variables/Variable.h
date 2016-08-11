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
    virtual void Request() const {
        requested = true;
    }

    virtual bool HasChanged() const {
        for(auto& p : GetNeededProcessors())
            if(p->HasChanged())
               return true;
        return false;
    }

    virtual ~Variable() = default;

protected:
    friend class ant::analysis::SlowControlManager;

    mutable bool requested = false;
    virtual void Init() {} // accessing the ExpConfig in the ctor is too early

    using ProcessorPtr = std::shared_ptr<Processor>;
    virtual std::list<ProcessorPtr> GetNeededProcessors() const =0;
};

using VariablePtr = std::shared_ptr<Variable>;

}}}