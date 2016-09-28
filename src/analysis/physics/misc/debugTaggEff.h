#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"


namespace ant {
namespace analysis {
namespace physics {


struct debugTaggEff: public Physics {

    static size_t getNchannels();

    unsigned SeenEvents = 0;

    const size_t nchannels;

    struct taggEffTree_t : WrapTTree {
        ADD_BRANCH_T(double, TaggerOrRate)
        ADD_BRANCH_T(double, P2Rate)
        ADD_BRANCH_T(std::vector<double>, TaggEffs)
        ADD_BRANCH_T(std::vector<double>, TaggEffErrors)
    };

    taggEffTree_t TaggEffTree;

    debugTaggEff(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~debugTaggEff();

    virtual void ProcessEvent(const TEvent&, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
