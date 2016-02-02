#pragma once

#include "analysis/physics/Physics.h"

#include "base/interval.h"



namespace ant {
namespace analysis {
namespace physics {

class EventFilter : public Physics {
protected:

    const interval<size_t> nCands;
    const double CBEsum;

    TH1D* steps;

public:
    EventFilter(const std::string& name, OptionsPtr opts);
    virtual ~EventFilter();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
};

}
}
}
