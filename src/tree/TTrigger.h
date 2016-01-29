#pragma once

#include "base/types.h"
#include "base/printable.h"

#include <list>
#include <memory>

namespace ant {

struct TTrigger : printable_traits {

    struct DAQError : printable_traits {
        index_t ModuleID;
        index_t ModuleIndex;
        int     ErrorCode;

        DAQError(index_t id, index_t index, int error_code):
            ModuleID(id),
            ModuleIndex(index),
            ErrorCode(error_code) {}
        DAQError() {}

        virtual ~DAQError() {}

        template<class Archive>
        void serialize(Archive& archive) {
            archive(ModuleID, ModuleIndex, ErrorCode);
        }


        std::ostream& Print(std::ostream& s) const {
            s << "DAQError(Module ID=" << ModuleID << " ModuleIndex="
              << ModuleIndex << " ErrorCode=" << ErrorCode << ")";
            return s;
        }
    };

    mev_t           CBEnergySum;
    unsigned int    ClusterMultiplicity;
    ns_t            CBTiming;

    std::list<DAQError>     Errors;

    TTrigger( mev_t CBESum=0.0,
                 unsigned int multiplicity=0):
        CBEnergySum(CBESum),
        ClusterMultiplicity(multiplicity),
        Errors()
    {}

    virtual ~TTrigger() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CBEnergySum, ClusterMultiplicity, CBTiming, Errors);
    }

    std::ostream& Print(std::ostream& stream) const {
        stream << "Trigger("
               << " CB Energy Sum=" << CBEnergySum << " MeV"
               << " Multipicity=" << ClusterMultiplicity
               << ")";
        for(auto& error : Errors) {
            stream << "\t" << error << "\n";
        }
        return stream;
    }
};

}

