#pragma once

#include "base/types.h"
#include "base/printable.h"
#include "tree/TDataRecord.h"

#include <list>
#include <memory>

namespace ant {
namespace analysis {
namespace data {

class DAQError: public printable_traits {
protected:
    index_t module_id;
    index_t module_index;
    int     error;

public:

    DAQError(index_t id, index_t index, int error_code):
        module_id(id),
        module_index(index),
        error(error_code) {}

    virtual ~DAQError() {}

    index_t ModuleID() const { return module_id; }
    index_t ModuleIndex() const { return module_index; }
    int     ErrorCode() const { return error; }

    std::ostream& Print(std::ostream& stream) const;
};



class TriggerInfo: public printable_traits {
protected:
    mev_t           cb_energy_sum;
    unsigned int    cluster_multiplicity;
    TID             event_id;

    std::list<DAQError>     errors;

public:

    TriggerInfo( mev_t CBESum=0.0,
                 unsigned int multiplicity=0,
                 const TID& event_id_ = TID()):
        cb_energy_sum(CBESum),
        cluster_multiplicity(multiplicity),
        event_id(event_id_),
        errors()
    {}

    virtual ~TriggerInfo() {}

    std::list<DAQError>  Errors() const { return errors; }
    std::list<DAQError>& Errors()       { return errors; }

    mev_t  CBEenergySum() const { return cb_energy_sum; }
    mev_t& CBEenergySum()      { return cb_energy_sum; }

    unsigned int  Multiplicity() const { return cluster_multiplicity; }
    unsigned int& Multiplicity()       { return cluster_multiplicity; }

    TID  EventID() const { return event_id; }
    TID& EventID()       { return event_id; }

    std::ostream& Print(std::ostream& stream) const;
};

}
}
}
