#pragma once

#include "base/types.h"
#include "base/printable.h"

#include <list>
#include <memory>

namespace ant {

struct TDAQError : printable_traits {
    index_t ModuleID;
    index_t ModuleIndex;
    int     ErrorCode;
    std::string ModuleName; // might be empty or unknown

    TDAQError(index_t id, index_t index, int error_code,
              const std::string& modname = ""):
        ModuleID(id),
        ModuleIndex(index),
        ErrorCode(error_code),
        ModuleName(modname)
    {}
    TDAQError() {}

    virtual ~TDAQError() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(ModuleID, ModuleIndex, ErrorCode, ModuleName);
    }


    std::ostream& Print(std::ostream& s) const {
        s << "TDAQError(Module ID=" << ModuleID;
        if(!ModuleName.empty())
            s << " Name=" << ModuleName;
        s  << " ModuleIndex="
          << ModuleIndex << " ErrorCode=" << ErrorCode << ")";
        return s;
    }
};

struct TTrigger : printable_traits {

    mev_t           CBEnergySum;
    unsigned int    ClusterMultiplicity;
    ns_t            CBTiming;

    std::vector<TDAQError> DAQErrors;

    TTrigger( mev_t CBESum=0.0,
                 unsigned int multiplicity=0):
        CBEnergySum(CBESum),
        ClusterMultiplicity(multiplicity),
        DAQErrors()
    {}

    virtual ~TTrigger() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CBEnergySum, ClusterMultiplicity, CBTiming, DAQErrors);
    }

    std::ostream& Print(std::ostream& stream) const {
        stream << "Trigger("
               << " CB Energy Sum=" << CBEnergySum << " MeV"
               << " Multipicity=" << ClusterMultiplicity
               << ")";
        for(auto& error : DAQErrors) {
            stream << "\t" << error << "\n";
        }
        return stream;
    }
};

}

