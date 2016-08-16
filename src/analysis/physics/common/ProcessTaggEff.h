#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


struct ProcessTaggEff: public Physics {

    unsigned seenEvents = 0;

    TH1D* taggerChannels;


    struct TreeScalarReads : WrapTTree {
        ADD_BRANCH_T(TID,   LastID)
        ADD_BRANCH_T(int,   nEvtsPerRead)

        ADD_BRANCH_T(double,    ExpLifeTime)
        ADD_BRANCH_T(double,    ExpTriggerRate)
        ADD_BRANCH_T(int,       Exp1MHz)

        ADD_BRANCH_T(int,   BeamPolMon1MHz)

        ADD_BRANCH_T(double,              LGRate)
        ADD_BRANCH_T(std::vector<double>, TaggRates)

        ADD_BRANCH_T(std::vector<int>,                  TDCHits)
        ADD_BRANCH_T(std::vector<int>,                  CoincidentTDCHits)
        ADD_BRANCH_T(std::vector<std::vector<double>>,  TaggTimings)

    };

    TreeScalarReads scalarReads;

    ProcessTaggEff(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~ProcessTaggEff();

    virtual void ProcessEvent(const TEvent& ev, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
