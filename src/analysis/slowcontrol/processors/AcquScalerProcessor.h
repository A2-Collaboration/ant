#pragma once

#include "Processor.h"

#include "tree/TEventData.h"

#include <queue>

namespace ant {
namespace analysis {
namespace slowcontrol {
namespace processor {

struct AcquScalerVector : Processor {

    using value_t = decltype(TSlowControl::Payload_Int);

    virtual return_t ProcessEventData(const TEventData& recon,  physics::manager_t& manager) override {
        for(const TSlowControl& sc : recon.SlowControls) {
            if(sc.Name == name) {
                if(sc.Validity != TSlowControl::Validity_t::Backward)
                    throw Exception("Encountered AcquScaler with forward validity. That's strange.");

                queue.push(sc.Payload_Int);
                manager.SaveEvent();

                return return_t::Complete;
            }
        }
        return return_t::Buffer;
    }
    virtual void PopQueue() override {
        queue.pop();
    }

    value_t Get() const {
        return queue.front();
    }

    AcquScalerVector(const std::string& name) : name(name) {}

private:
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