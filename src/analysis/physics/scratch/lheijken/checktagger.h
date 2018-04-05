#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "expconfig/detectors/Tagger.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A playground class for checking tagger stuff
 *
 */
class scratch_lheijken_checktagger: public Physics {
protected:

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;


    TH2D* hTime;
    TH2D* hTimeToTriggerRef;
    TH2D* hTimeZoomed;
    TH2D* hTimeToTagger;
    TH2D* hTimeMultiplicity;
    TH2D* hTimeMultiplicityFromTHits;
    TH2D* hTimeMultTaggVsReadHits;
    TH1D* hTriggerRefTiming;
    TH1D* hNrHitInTDCModPerEv[3];
    TH2D* hTime_cutNrHitInTDCMod0[3];
    TH1D* hTestChanNrMod[3];
    TH1D* hFaradayCupScaler;

   struct TaggHitTree_t : WrapTTree {

        ADD_BRANCH_T(int, EventNumber)
        ADD_BRANCH_T(double, Time)
        ADD_BRANCH_T(int, Channel)
        ADD_BRANCH_T(int, ChannelMult)
        ADD_BRANCH_T(int, TimingInd)
    };
    TaggHitTree_t TaggHitTree;



protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();

public:
    scratch_lheijken_checktagger(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_checktagger() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
