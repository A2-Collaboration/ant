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
    TPC_PhysicsStats(const std::string& name, PhysOptPtr& opts);
    void ProcessEvent(const data::Event& event);
    void Finish();
    void ShowResult();
};

}
}
}
