#ifndef TRIGGERINFO_H
#define TRIGGERINFO_H

#include "base/types.h"
#include "base/printable.h"

#include <list>
#include <memory>

namespace ant {

class DAQError: public ant::printable_traits {
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

    const index_t ModuleID() const { return module_id; }
    const index_t ModuleIndex() const { return module_index; }
    const int     ErrorCode() const { return error; }

    std::ostream &Print(std::ostream &stream) const;
};



class TriggerInfo: public ant::printable_traits {
public:
    typedef std::list<ant::DAQError> ErrorList_t;  //TODO: use shared_ptr also for DAQErrors?

protected:
    mev_t           cb_energy_sum;
    unsigned int    cluster_multiplicity;

    ErrorList_t     errors;

public:

    TriggerInfo( mev_t CBESum=0.0, unsigned int multiplicity=0):
        cb_energy_sum(CBESum),
        cluster_multiplicity(multiplicity)
    {}

    virtual ~TriggerInfo() {}

    ErrorList_t& Errors()       { return errors; }
    const ErrorList_t& Errors() const { return errors; }

    mev_t CBEenergySum() const { return cb_energy_sum; }
    mev_t& CBEenergySum()      { return cb_energy_sum; }

    unsigned int Multiplicity() const { return cluster_multiplicity; }
    unsigned int& Multiplicity()       { return cluster_multiplicity; }

    std::ostream &Print(std::ostream &stream) const;
};

using TriggerInfoPtr = std::shared_ptr<TriggerInfo>;

}

#endif
