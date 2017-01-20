#pragma once

#include "event_t.h"
#include "SlowControlProcessors.h"

#include <queue>
#include <memory>

namespace ant {
namespace analysis {

class SlowControlManager {

protected:

    std::queue<slowcontrol::event_t> eventbuffer;

    using ProcessorPtr = std::shared_ptr<slowcontrol::Processor>;

    struct processor_t {
        processor_t(ProcessorPtr proc) : Processor(proc) {}
        ProcessorPtr    Processor;
        std::list<TID>  CompletionPoints;

        enum class type_t {
            Unknown, Backward, Forward
        };
        type_t Type = type_t::Unknown;

        bool IsComplete() const;
    };

    std::vector<processor_t> processors;

    void AddProcessor(ProcessorPtr p);

public:
    SlowControlManager();

    bool ProcessEvent(TEvent event);

    slowcontrol::event_t PopEvent();

    size_t BufferSize() const { return eventbuffer.size(); }

};

}} // namespace ant::analysis
