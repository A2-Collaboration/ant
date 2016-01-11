#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/KinFitter.h"
#include "base/ParticleTypeTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
protected:

    struct MultiPi0 {
        MultiPi0(SmartHistFactory& histFac, unsigned nPi0, bool nofitandnotree = false);

        void ProcessData(const data::Event::Data& data, const data::ParticleTree_t& ptree);
        void ShowResult();

    protected:

        const unsigned multiplicity;
        const bool skipfit;
        ParticleTypeTree directPi0;

        utils::KinFitter fitter;

        TTree* tree;
        double   Tagg_W;  // tagger prompt/random weight for subtraction
        unsigned Tagg_Ch; // tagger channel
        double   Tagg_E;  // tagger energy

        TCandidate Proton;
        TCandidate ProtonMCTrue;
        unsigned ProtonMCTrueMatches;

        std::vector<TCandidate> Photons;


        TH1D* steps;
        TH1D* Proton_Coplanarity;
        TH1D* Proton_Angle_True;
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
