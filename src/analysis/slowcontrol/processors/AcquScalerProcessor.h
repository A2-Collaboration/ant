#pragma once

#include "Processor.h"

#include "tree/TEventData.h"

#include <queue>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct AcquScalerVector : Processor {

    AcquScalerVector(const std::string& name) : name(name) {}

    using value_t = decltype(TSlowControl::Payload_Int);

    virtual return_t ProcessEventData(const TEventData& recon,  physics::manager_t& manager) override;

    virtual void PopQueue() override;
    virtual void Reset() override;

    value_t Get() const;


private:
    bool firstScalerSeen = false;
    const std::string name;
    std::queue<value_t> queue;
};


struct AcquScalerScalar : AcquScalerVector {
    using return_t = decltype(value_t::value_type::Value);

    using AcquScalerVector::AcquScalerVector;

    return_t Get() {
        return AcquScalerVector::Get().front().Value;
    }

protected:
    // hide this get method from base class
    using AcquScalerVector::Get;

};



}}}} // namespace ant::analysis::slowcontrol::processor