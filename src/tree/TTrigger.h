#pragma once

#include "base/types.h"
#include "base/printable.h"

#include <list>
#include <limits>
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
    unsigned        DAQEventID; // might wrap-around at 16bit integer...

    std::vector<TDAQError> DAQErrors;

    TTrigger():
        CBEnergySum(std::numeric_limits<double>::quiet_NaN()),
        ClusterMultiplicity(0),
        CBTiming(std::numeric_limits<double>::quiet_NaN()),
        DAQEventID(0),
        DAQErrors()
    {}

    virtual ~TTrigger() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CBEnergySum, ClusterMultiplicity, CBTiming,
                DAQEventID, DAQErrors);
    }

    std::ostream& Print(std::ostream& s) const {
        s << "Trigger"
          << " CBEnergySum=" << CBEnergySum
          << " Multipicity=" << ClusterMultiplicity
          << " CBTiming=" << CBTiming
          << " DAQEventId=0x" << std::hex << DAQEventID << std::dec
          << " nDAQErrors=" << DAQErrors.size();
        return s;
    }
};

}

