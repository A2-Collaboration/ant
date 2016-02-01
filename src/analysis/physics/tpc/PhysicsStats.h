#pragma once

#include "physics/Physics.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief Physics class to investigate angular distributions of charged particles.
 *        Intended for studies of coverage of the CBTPC and usefulness for A2 physics.
 *        Runs on MCTrue Data only.
 */
class TPC_PhysicsStats: public Physics {
protected:
    std::map<std::string, TH1D*> decay_angles;

public:
    TPC_PhysicsStats(const std::string& name, OptionsPtr& opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
