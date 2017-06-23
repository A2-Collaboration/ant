#include "AcquScalerProcessor.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::processor;

Processor::return_t AcquScalerVector::ProcessEventData(const TEventData& recon,
                                                       manager_t& manager) {
    // search for the TSlowControl with the type/name
    for(const TSlowControl& sc : recon.SlowControls) {
        if(sc.Type != TSlowControl::Type_t::AcquScaler)
            continue;
        if(sc.Name != name)
            continue;
        if(sc.Validity != TSlowControl::Validity_t::Backward)
            throw Exception("Encountered AcquScaler with forward validity. That's strange.");

        manager.SaveEvent();

        // very first scaler seen?
        if(!firstScalerSeen) {
            firstScalerSeen = true;
            return return_t::Skip;
        }

        queue.push(sc.Payload_Int);
        return return_t::Complete;

    }

    if(firstScalerSeen)
        return return_t::Buffer;
    else
        return return_t::Skip;
}

void AcquScalerVector::PopQueue() {
    queue.pop();
}

AcquScalerVector::value_t AcquScalerVector::Get() const {
    // if this assert fails, probably a physics class forgot
    // to request the slowcontrol variable in its constructor
    // see DebugPhysics how to it properly
    assert(!queue.empty());
    return queue.front();
}
