#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

#include <string>


class TH1D;

namespace ant {
namespace analysis {
namespace physics {


class scratch_wagner_MCGunChecks : public Physics {
protected:

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(std::vector<std::string>,  names)
        ADD_BRANCH_T(std::vector<double>,       energies)
        ADD_BRANCH_T(std::vector<double>,       energies_true)
        ADD_BRANCH_T(std::vector<double>,       thetas)
        ADD_BRANCH_T(std::vector<double>,       thetas_true)
        ADD_BRANCH_T(std::vector<double>,       phis)
        ADD_BRANCH_T(std::vector<double>,       phis_true)
        ADD_BRANCH_T(unsigned,                  multiplicity)

        void fillAndReset()
        {
            Tree->Fill();
            names().resize(0);
            energies().resize(0);
            energies_true().resize(0);
            thetas().resize(0);
            thetas_true().resize(0);
            phis().resize(0);
            phis_true().resize(0);
            multiplicity = 0;
        }
    };

    Tree_t t;


public:
    scratch_wagner_MCGunChecks(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    virtual ~scratch_wagner_MCGunChecks(){}
};

}
}
}
