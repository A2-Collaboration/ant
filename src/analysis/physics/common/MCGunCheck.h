#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

#include <string>

class TH1D;
class TH2D;
class TH3D;


namespace ant {


namespace analysis {
namespace physics {

class MCGunCheck : public Physics {
protected:

    struct tree_t : WrapTTree {
        ADD_BRANCH_T(std::vector<double>,      MCopenings)
    };

    tree_t t;


public:
    MCGunCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
