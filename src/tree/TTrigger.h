#pragma once

#include "base/types.h"

#include <vector>
#include <limits>
#include <memory>

namespace ant {

struct TDAQError {
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
    TDAQError() = default;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(ModuleID, ModuleIndex, ErrorCode, ModuleName);
    }


    friend std::ostream& operator<<(std::ostream& s, const TDAQError& o) {
        s << "TDAQError(Module ID=" << o.ModuleID;
        if(!o.ModuleName.empty())
            s << " Name=" << o.ModuleName;
        s  << " ModuleIndex="
          << o.ModuleIndex << " ErrorCode=" << o.ErrorCode << ")";
        return s;
    }
};

struct TTrigger {

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

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CBEnergySum, ClusterMultiplicity, CBTiming,
                DAQEventID, DAQErrors);
    }

    friend std::ostream& operator<<(std::ostream& s, const TTrigger& o) {
        s << "Trigger"
          << " CBEnergySum(measured)=" << o.CBEnergySum
          << " Multipicity=" << o.ClusterMultiplicity
          << " CBTiming(measured)=" << o.CBTiming
          << " DAQEventId=0x" << std::hex << o.DAQEventID << std::dec
          << " nDAQErrors=" << o.DAQErrors.size();
        return s;
    }
};

}

