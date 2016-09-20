#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/Fitter.h"
#include "base/ParticleTypeTree.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "base/WrapTTree.h"
#include "analysis/utils/PullsWriter.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class JustPi0 : public Physics {
public:

    struct MultiPi0 {
        MultiPi0(HistogramFactory& histFac, unsigned nPi0, utils::UncertaintyModelPtr FitterModel, bool nofitandnotree = false);

        void ProcessData(const TEventData& data, const TParticleTree_t& ptree);
        void ShowResult();

    protected:

        const unsigned multiplicity;
        HistogramFactory HistFac;
        const unsigned nPhotons_expected;
        const bool skipfit;
        ParticleTypeTree directPi0;

        utils::UncertaintyModelPtr model;
        utils::KinFitter fitter;

        std::vector<std::pair<utils::TreeFitter::tree_t,utils::TreeFitter::tree_t>> pions;

        public:
        struct MultiPi0Tree : WrapTTree {

            ADD_BRANCH_T(bool,  isMC)

            ADD_BRANCH_T(double,   Tagg_W)
            ADD_BRANCH_T(unsigned, Tagg_Ch)
            ADD_BRANCH_T(double,   Tagg_E)

            ADD_BRANCH_T(double, CBAvgTime)

            ADD_BRANCH_T(std::vector<double>, ggIM)

            ADD_BRANCH_T(double, kinfit_chi2dof)
            ADD_BRANCH_T(double, kinfit_prob)
            ADD_BRANCH_T(double, treefit_chi2dof)
            ADD_BRANCH_T(double, treefit_prob)

            ADD_BRANCH_T(unsigned, ProtonMCTrueMatches)

            ADD_BRANCH_T(TLorentzVector, proton)
            ADD_BRANCH_T(TLorentzVector, proton_fitted)
            ADD_BRANCH_T(TVector2,       proton_PSA)
            ADD_BRANCH_T(double,         proton_vetoE)
            ADD_BRANCH_T(double,         proton_Time)
            ADD_BRANCH_T(unsigned,       proton_vetoCh)
            ADD_BRANCH_T(unsigned,       proton_det)


            ADD_BRANCH_T(std::vector<TLorentzVector>, photons)
            ADD_BRANCH_T(std::vector<TLorentzVector>, photons_fitted)
            ADD_BRANCH_T(std::vector<TVector2>, photons_PSA)
            ADD_BRANCH_T(std::vector<double>, photons_vetoE)
            ADD_BRANCH_T(std::vector<double>, photons_Time)

            ADD_BRANCH_T(double, Tagg_E_fitted)
            ADD_BRANCH_T(double, fit_Tagg_E_pull)

            ADD_BRANCH_T(std::vector<double>, fit_photons_E_pulls)
            ADD_BRANCH_T(std::vector<double>, fit_photons_Theta_pulls)
            ADD_BRANCH_T(std::vector<double>, fit_photons_Phi_pulls)
            ADD_BRANCH_T(double, fit_proton_E_pull)
            ADD_BRANCH_T(double, fit_proton_Theta_pull)
            ADD_BRANCH_T(double, fit_proton_Phi_pull)
        };

    protected:

        TTree* tree;
        MultiPi0Tree t;

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

protected:

    utils::UncertaintyModelPtr model;
    std::vector<std::unique_ptr<MultiPi0>> multiPi0;

public:
    JustPi0(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};

}}}
