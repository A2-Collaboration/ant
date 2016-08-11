#pragma once

#include "analysis/physics/manager_t.h"
#include "tree/TEventData.h"

namespace ant {
namespace analysis {
namespace slowcontrol {

struct Processor {
    enum class return_t {
        Process,   // event processable without forward buffering
        Complete,  // this event made the processor complete
        Buffer,    // buffer more events
        Skip,      // skip that whole event (unprocessable)
    };

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) =0;

    virtual void PopQueue() = 0;
    virtual void Reset() = 0;

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    virtual ~Processor() = default;
};

}}}