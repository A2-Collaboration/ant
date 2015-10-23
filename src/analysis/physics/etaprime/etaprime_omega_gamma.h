#pragma once

#include "analysis/physics/Physics.h"
#include "utils/particle_tools.h"

#include "base/ParticleTypeTree.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class EtapOmegaG : public Physics {

    struct expected_peak_t {
        double Mean;
        double Sigma;
        expected_peak_t(double mean, double sigma) :
            Mean(mean), Sigma(sigma) {}
        interval<double> makeCutInterval(unsigned nSigma=2) const {
            return {Mean-nSigma*Sigma,Mean+nSigma*Sigma};
        }
    };

    // means/sigma extracted from gg/ggg/gggg histograms for signal channel
    const expected_peak_t Pi0 = {126, 15};
    const expected_peak_t Eta = {515, 18};
    const expected_peak_t Omega = {735, 32};
    const expected_peak_t EtaPrime_sig = {895, 27};
    // extracted from gg histogram for reference channel
    const expected_peak_t EtaPrime_ref = {905, 29};

    ParticleTypeTree treeSignal;
    ParticleTypeTree treeReference;


    SmartHistFactory sig_HistFac;
    SmartHistFactory ref_HistFac;

    TH1D* h_TotalEvents;

    struct histogram_t {
        TH1D* Steps;
        TH1D* MissedBkg;
        histogram_t(SmartHistFactory HistFac);
    };

    histogram_t sig_hists;
    histogram_t ref_hists;

    template<typename T>
    struct perDecayHists_t {
        ParticleTypeTree Tree;
        std::string ShortName;
        std::string DecayString;
        T PerDecayHists;
        perDecayHists_t(const std::string& shortName,
                    const SmartHistFactory& HistFac_parent,
                    ParticleTypeTree tree
                    ) :
            Tree(tree),
            ShortName(shortName),
            DecayString(utils::ParticleTools::GetDecayString(tree)),
            PerDecayHists(SmartHistFactory(ShortName, HistFac_parent, DecayString))
        {}
        perDecayHists_t(const std::string& shortName, const SmartHistFactory& HistFac_parent) :
            Tree(nullptr),
            ShortName(shortName),
            DecayString(shortName),
            PerDecayHists(SmartHistFactory(ShortName, HistFac_parent, DecayString))
        {}
    };

    struct sig_perDecayHists_t {
        TH1D* Steps;

        TH1D* gggg;
        TH1D* ggg;
        TH1D* gg;

        TH2D* IM_gg_gg;
        TH2D* IM_gg_gg_cut;


        TH1D* Proton_Copl;

        TH2D* IM_etap_omega;
        TH1D* IM_pi0;

        TH1D* MM_gggg;
        TH1D* MM_etap;

        TH1D* Chi2_All;
        TH1D* Chi2_Best;

        TH1D* g_EtaPrime_E;

        sig_perDecayHists_t(SmartHistFactory HistFac);
    };

    std::vector<perDecayHists_t<sig_perDecayHists_t>> sig_perDecayHists;

    struct ref_perDecayHists_t {
        TH1D* Steps;

        TH1D* gg;

        TH2D* Proton_ThetaPhi;
        TH1D* Proton_Energy;

        TH1D* IM_etap;

        TH1D* MM_etap;

        TH1D* Proton_Copl;

        ref_perDecayHists_t(SmartHistFactory HistFac);
    };

    std::vector<perDecayHists_t<ref_perDecayHists_t>> ref_perDecayHists;

    struct sig_TTree_t {
        sig_TTree_t(TTree* tree) : Tree(tree) {}

        TTree* Tree;
        int MCTrueIndex = -1;

        utils::ParticleVars Proton;
        utils::ParticleVars ProtonTrue;

        double Chi2;
        utils::ParticleVars g_Pi0_0;
        utils::ParticleVars g_Pi0_1;
        utils::ParticleVars g_Omega;
        utils::ParticleVars g_EtaPrime;
        utils::ParticleVars g_EtaPrime_Boosted;
        utils::ParticleVars Pi0;
        utils::ParticleVars Omega;
        utils::ParticleVars EtaPrime;

        void SetBranches();
    };
    sig_TTree_t sig_TTree;

    struct ref_TTree_t {
        ref_TTree_t(TTree* tree) : Tree(tree) {}
        TTree* Tree;
        int MCTrueIndex = -1;

        utils::ParticleVars Proton;

        void SetBranches();
    };
    ref_TTree_t ref_TTree;


    void ProcessSig(const data::ParticleTree_t& particletree, const data::Event::Data& data);
    void ProcessRef(const data::ParticleTree_t& particletree, const data::Event::Data& data);


public:
    EtapOmegaG(const std::string& name, PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};


}}}
