#pragma once

#include "analysis/physics/manager_t.h"
#include "tree/TEventData.h"

namespace ant {
namespace analysis {
namespace slowcontrol {

struct Processor {
    enum class return_t {
        Buffer,    // buffer more events
        Complete,  // this event made
        Skip       // skip that whole event
    };

    virtual return_t ProcessEventData(const TEventData& recon, physics::manager_t& manager) =0;

    virtual void PopQueue() =0;

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    virtual ~Processor() = default;
};

}}}