#include "TaggEff.h"

#include "base/Logger.h"
#include "calibration/modules/TaggEff.h"

#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::processor;

Processor::return_t TaggEff::ProcessEventData(const TEventData& recon,
                                              manager_t& manager) {
    // search for the TSlowControl with the name
    for(const TSlowControl& sc : recon.SlowControls) {
        if(sc.Type != TSlowControl::Type_t::Ant)
            continue;
        /// \bug Here we ignore that there might be more than one tagger present
        if(!std_ext::string_ends_with(sc.Name, calibration::TaggEff::GetModuleNameSuffix()))
            continue;
        if(sc.Validity != TSlowControl::Validity_t::Forward)
            throw Exception("Encountered TaggEff slowcontrol with backward validity. That's strange.");

        manager.SaveEvent();

        stringstream ss(sc.Payload_String.front().Value);
        cereal::BinaryInputArchive ar(ss);

        queue.emplace();
        ar(queue.back());

        return return_t::Complete;
    }

    if(queue.empty())
        return return_t::Skip;
    else
        return return_t::Process;
}

void TaggEff::PopQueue() {
    queue.pop();
}

TaggEff::value_t TaggEff::Get() const {
    // if this assert fails, probably a physics class forgot
    // to request the slowcontrol variable in its constructor
    // see DebugPhysics how to it properly
    assert(!queue.empty());
    return queue.front();
}
