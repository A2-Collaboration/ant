#pragma once

#include "analysis/physics/Physics.h"
#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm_mc : public Physics {
public:
    Omega_EpEm_mc(const std::string& name, OptionsPtr opts);

    void ProcessEvent(const TEvent& event, manager_t&) override;
    void Finish() override;
    void ShowResult() override;

    struct tree_t : WrapTTree {
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(unsigned, nClusters)
//        ADD_BRANCH_T(std::vector<double>, Numbers, 3) // vector Numbers is three items large now
    };

protected:
    //using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;
    //const std::shared_ptr<decaytree_t> channelTree;

private:
    tree_t t;

};

}}} // end of ant::analysis::physics
