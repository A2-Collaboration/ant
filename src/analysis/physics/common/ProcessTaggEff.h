#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


struct ProcessTaggEff: public Physics {

    unsigned seenEvents = 0;
    TH1D* byChannel;


    struct TreeData : WrapTTree {
        ADD_BRANCH_T(TID, StartID)
        ADD_BRANCH_T(std::vector<double>,   AvgTaggFreqs)
        ADD_BRANCH_T(double,   AvgLGFreq)
    };

    TreeData mainTree;

    ProcessTaggEff(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~ProcessTaggEff();

    virtual void ProcessEvent(const TEvent& ev, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
