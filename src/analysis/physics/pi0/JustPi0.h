#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/KinFitter.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
protected:

    struct MultiPi0 {
        MultiPi0(SmartHistFactory& histFac, unsigned nPi0);

        void ProcessData(const data::Event::Data& data);
        void ShowResult();

    protected:

        const unsigned multiplicity;

        utils::KinFitter fitter;

        TTree* tree;
        double BestFitProbability;
        TCandidate Proton;
        std::vector<TCandidate> Photons;

        TH1D* steps;
        TH1D* Proton_Coplanarity;
        PromptRandom::Switch promptrandom;

        PromptRandom::Hist1 h_missingmass;
        PromptRandom::Hist1 h_fitprobability;
        PromptRandom::Hist1 IM_2g_byMM;
        PromptRandom::Hist1 IM_2g_byFit;
        PromptRandom::Hist1 IM_2g_fitted;

    };

    std::vector<std::unique_ptr<MultiPi0>> multiPi0;

public:
    JustPi0(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};

}}}
