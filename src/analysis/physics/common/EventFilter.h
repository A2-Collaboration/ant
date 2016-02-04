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
    const double maxCoplAngle;

    TH1D* steps;

    /**
     * @brief check if at least one combination passes "coplanarity" test
     * @param cands list of candidates
     * @param maxangle
     * @return
     */
    static bool checkCoplanarity(const TCandidateList& cands, const double maxangle);

public:
    EventFilter(const std::string& name, OptionsPtr opts);
    virtual ~EventFilter();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
};

}
}
}
