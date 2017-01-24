#pragma once

#include "physics/Physics.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {

namespace utils { struct MCWeighting_test; }

namespace physics {

class TestMCWeighting : public Physics {

    struct tree_t : WrapTTree {
        ADD_BRANCH_T(double, BeamE)
        ADD_BRANCH_T(double, CosTheta)
    };
    tree_t t;

    std::unique_ptr<utils::MCWeighting_test> mcWeighting;

public:
    TestMCWeighting(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};


}}} // end of ant::analysis::physics