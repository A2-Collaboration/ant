#pragma once

#include "base/types.h"
#include "base/printable.h"

#include <vector>
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

    // NB: Those values are measured by the system,
    // in particular, the CBEnergySum/CBTiming will be NaN for most A2 beamtimes,
    // as this is NOT measured by any device in the setup
    mev_t           CBEnergySum;         // usually not provided by A2 setup
    unsigned int    ClusterMultiplicity;
    ns_t            CBTiming;            // usually not provided by A2 setup
    unsigned        DAQEventID;          // typically wraps around at 16bit integer...

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
          << " CBEnergySum(measured)=" << CBEnergySum
          << " Multipicity=" << ClusterMultiplicity
          << " CBTiming(measured)=" << CBTiming
          << " DAQEventId=0x" << std::hex << DAQEventID << std::dec
          << " nDAQErrors=" << DAQErrors.size();
        return s;
    }
};

}

