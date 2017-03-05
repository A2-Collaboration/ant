#pragma once

#include "analysis/physics/Physics.h"
#include "base/ParticleTypeTree.h"

#include <map>

namespace ant {
namespace analysis {
namespace physics {


class MCTrueOverview : public Physics {
protected:

    struct perChannel_t{

    };

    std::map<ParticleTypeTreeDatabase::Channel, perChannel_t> channels;

public:
    MCTrueOverview(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
