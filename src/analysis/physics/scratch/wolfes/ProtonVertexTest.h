#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

#include "base/std_ext/string.h"

#include "base/WrapTTree.h"

#include "TLorentzVector.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


class ProtonVertexTest: public Physics{

public:
    //storage
    struct ptree : WrapTTree
    {
        ADD_BRANCH_T(double, theta)
        ADD_BRANCH_T(double, ekin)
        ADD_BRANCH_T(double, coplanarity)
        ADD_BRANCH_T(double, mm_im)
        ADD_BRANCH_T(double, mm_angle)
        ADD_BRANCH_T(double, prob)
        ADD_BRANCH_T(double, zvertex)

        ADD_BRANCH_T(double,                      photonVeto)
        ADD_BRANCH_T(double,                      corrPhotonVeto)

        ADD_BRANCH_T(double, prW)
    };
    ptree tree;

    //hists
    TH1D* hist_steps  = nullptr;
    TH1D* hist_theta  = nullptr;

    // implementation
    ProtonVertexTest(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override {}
    virtual void ShowResult() override;

private:

    // cuts
    const interval<size_t>    NCands     = {7,7};
    const IntervalD           CBESum     = {550,std_ext::inf};
    const IntervalD           ProtonCopl = {-25,25};
    const IntervalD           MM         = {600,1300};
    const IntervalD           MMAngle    = {0,25};
    const IntervalD           EMB_prob   = {0.005,1};

    // prompt random
    ant::analysis::PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    const IntervalD              Range_Prompt  =   { -5,  5};
    const std::vector<IntervalD> Ranges_Random = { {-55,-10},
                                                   { 10, 55}  };

    //kinFit
    const double fitter_ZVertex = 3;
    std::shared_ptr<utils::UncertaintyModel> uncertModel;
    utils::KinFitter kinFitterEMB;

    //helpers
    void FillStep(const std::string& step) {hist_steps->Fill(step.c_str(),1);}



};


//namespaces
}
}
}
