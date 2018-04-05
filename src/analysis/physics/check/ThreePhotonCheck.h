#pragma once

#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class ThreePhotonCheck : public Physics {

    PromptRandom::Switch promptrandom;

    utils::TriggerSimulation triggersimu;
    utils::KinFitter fitter;

    TH1D* h_Steps;
    TH1D* h_photonsIM_fitted;
    TH1D* h_photonsIM_raw;
    TH1D* h_mm;

    template <typename it_type>
    static LorentzVec LVSum(it_type begin, it_type end) {
        LorentzVec v;

        while(begin!=end) {
            v += **begin;
            ++begin;
        }

        return v;
    }

public:
    ThreePhotonCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}}}
