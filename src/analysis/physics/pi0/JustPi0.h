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

        std::vector<std::pair<utils::TreeFitter::tree_t,utils::TreeFitter::tree_t>> pions;

        TTree* tree;

        // branches
        double   b_Tagg_W;  // tagger prompt/random weight for subtraction
        unsigned b_Tagg_Ch; // tagger channel
        double   b_Tagg_E;  // tagger energy


        TVector2  b_proton_PSA;
        TVector2  b_proton_PSA1;
        TVector2  b_proton_PSA2;
        std::vector<double> b_ggIM;

        double b_proton_vetoE;

        double b_treefit_chi2dof;
        double b_treefit_prob;

        unsigned ProtonMCTrueMatches;


        TLorentzVector b_proton_fitted;
        std::vector<TLorentzVector> b_photons_fitted;
        double b_beamE_fitted;
        double b_fit_beamE_pull;

        std::vector<double> b_fit_photons_E_pulls;
        std::vector<double> b_fit_photons_Theta_pulls;
        std::vector<double> b_fit_photons_Phi_pulls;
        double b_fit_proton_E_pull;
        double b_fit_proton_Theta_pull;
        double b_fit_proton_Phi_pull;

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
