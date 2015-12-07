#pragma once

#include "base/types.h"
#include "base/printable.h"
#include "tree/TDataRecord.h"

#include <list>
#include <memory>

namespace ant {
namespace analysis {
namespace data {

struct DAQError : printable_traits {
    index_t ModuleID;
    index_t ModuleIndex;
    int     ErrorCode;

    DAQError(index_t id, index_t index, int error_code):
        ModuleID(id),
        ModuleIndex(index),
        ErrorCode(error_code) {}

    virtual ~DAQError() {}

    std::ostream& Print(std::ostream& stream) const;
};



struct Trigger_t :  printable_traits {
    mev_t           CBEnergySum;
    unsigned int    ClusterMultiplicity;
    TID             EventID;

    std::list<DAQError>     Errors;

    Trigger_t( mev_t CBESum=0.0,
                 unsigned int multiplicity=0,
                 const TID& eventID = TID()):
        CBEnergySum(CBESum),
        ClusterMultiplicity(multiplicity),
        EventID(eventID),
        Errors()
    {}

    virtual ~Trigger_t() {}

    std::ostream& Print(std::ostream& stream) const;
};

}
}
}
