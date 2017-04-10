#pragma once

#include "analysis/physics/Physics.h"
#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"
#include "utils/A2GeoAcceptance.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class Omega_EpEm_mc : public Physics {
public:
    Omega_EpEm_mc(const std::string& name, OptionsPtr opts);

    void ProcessEvent(const TEvent& event, manager_t&) override;
    void Finish() override;
    void ShowResult() override;


protected:
    struct tree_t : WrapTTree {
        // charged e
        ADD_BRANCH_T(unsigned,                      nClusters)
        ADD_BRANCH_T(int,                           hitsCB)
        ADD_BRANCH_T(int,                           hitsTAPS)
        ADD_BRANCH_T(int,                           nEcharged)
        ADD_BRANCH_T(int,                           nEminus)
        ADD_BRANCH_T(int,                           nEplus)
        ADD_BRANCH_T(double,                        eeOpenAngle)
        ADD_BRANCH_T(double,                        eeIM)
        ADD_BRANCH_T(double,                        eAngle)
        ADD_BRANCH_T(double,                        eEk)

        // contains all charged particles
        ADD_BRANCH_T(std::vector<TLorentzVector>,   eVector)

        // proton:
        ADD_BRANCH_T(double,                        pEk)
        ADD_BRANCH_T(double,                        pAngle)

        void fillAndReset() // to fill the tree and reset all values which add up
        {
            Tree->Fill();
            hitsCB   = 0;
            hitsTAPS = 0;
            eVector().resize(0);
        }
    };

    ant::analysis::utils::A2SimpleGeometry a2geo;


    //using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;
    //const std::shared_ptr<decaytree_t> channelTree;


private:
    tree_t t;

};

}}} // end of ant::analysis::physics
