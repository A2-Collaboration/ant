#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/Fitter.h"
#include "base/ParticleTypeTree.h"
#include "TLorentzVector.h"
#include "TVector2.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
protected:

    struct MultiPi0 {
        MultiPi0(HistogramFactory& histFac, unsigned nPi0, bool nofitandnotree = false);

        void ProcessData(const TEventData& data, const TParticleTree_t& ptree);
        void ShowResult();

    protected:

        const unsigned multiplicity;
        const bool skipfit;
        ParticleTypeTree directPi0;

        std::shared_ptr<utils::Fitter::UncertaintyModel> model;
        utils::KinFitter fitter;

        TTree* tree;
        TTree* buffer;

        // branches
        double   Tagg_W;  // tagger prompt/random weight for subtraction
        unsigned Tagg_Ch; // tagger channel
        double   Tagg_E;  // tagger energy

        TLorentzVector b_proton;
        TLorentzVector b_proton_fitted;
        TVector2  b_proton_PSA;
        TVector2  b_proton_PSA1;
        TVector2  b_proton_PSA2;

        double b_proton_vetoE;

        double treefit_chi2dof;
        double treefit_prob;

        double gg_IM;

        unsigned ProtonMCTrueMatches;

        TH1D* steps;
        TH1D* Proton_Coplanarity;
        TH1D* Proton_Angle_True;
        PromptRandom::Switch promptrandom;

        PromptRandom::Hist1 h_missingmass;
        PromptRandom::Hist1 h_fitprobability;
        PromptRandom::Hist1 IM_2g_byMM;
        PromptRandom::Hist1 IM_2g_byFit;
        PromptRandom::Hist1 IM_2g_fitted;

        utils::TreeFitter treefitter;

        static ParticleTypeTree getParticleTree(const unsigned nPi0);

    };

    std::vector<std::unique_ptr<MultiPi0>> multiPi0;

public:
    JustPi0(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}}
