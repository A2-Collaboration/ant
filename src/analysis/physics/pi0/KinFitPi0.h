#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/Fitter.h"
#include "base/ParticleTypeTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class KinFitPi0 : public Physics {
protected:

    struct MultiPi0 {
        MultiPi0(HistogramFactory& histFac, unsigned nPi0);

        void ProcessData(const TEventData& data);
        void ShowResult();

    protected:

        const unsigned multiplicity;

        utils::KinFitter fitter;

        TH1D* steps;

        TH1D* Proton_Coplanarity;
        TH1D* h_taggtime;

        PromptRandom::Switch promptrandom;

        PromptRandom::Hist1 h_missingmass;
        PromptRandom::Hist1 h_fitprobability;
        PromptRandom::Hist1 IM_2g_byFit;
        PromptRandom::Hist1 IM_2g_fitted;

    };

    std::vector<std::unique_ptr<MultiPi0>> multiPi0;

public:
    KinFitPi0(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
