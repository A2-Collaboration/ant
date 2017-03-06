#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"
#include "utils/A2GeoAcceptance.h"

#include "TLorentzVector.h"
#include <string>

class TH1D;
class TH2D;
class TH3D;


namespace ant {


namespace analysis {
namespace physics {

class MCPhotonPairCheck : public Physics {
protected:

    struct tree_t : WrapTTree {
        ADD_BRANCH_T(int,                           multiplicity)
        ADD_BRANCH_T(int,                           hitsCB)
        ADD_BRANCH_T(int,                           hitsTAPS)
        ADD_BRANCH_T(std::vector<double>,           openings)
        ADD_BRANCH_T(std::vector<TLorentzVector>,   p)


        void fillAndReset()
        {
            Tree->Fill();
            multiplicity = 0;
            hitsCB   = 0;
            hitsTAPS = 0;
            openings().resize(0);
            p().resize(0);
        }
    };

    tree_t t;

    ant::analysis::utils::A2SimpleGeometry a2geo;

public:
    MCPhotonPairCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    virtual ~MCPhotonPairCheck(){}
};

}
}
}
