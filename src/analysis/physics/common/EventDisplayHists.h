#pragma once

#include "analysis/physics/Physics.h"

#include <string>
#include <memory>
#include <vector>


namespace ant {
namespace analysis {
namespace physics {

class EventDisplayHists : public Physics {
protected:
    unsigned n = 0;

    double tapsZ = {};

    static constexpr int clusterMarker = 25;  // Empty square
    static constexpr int trueMarker    = 24;  // Empty circle

public:
    EventDisplayHists(const std::string& name, OptionsPtr opts);
    virtual ~EventDisplayHists();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
};

}
}
}
