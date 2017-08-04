#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/utils/uncertainties/Optimized.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/plot/PromptRandomHist.h"
#include "base/WrapTTree.h"
#include "analysis/utils/TriggerSimulation.h"
#include "analysis/utils/ClusterTools.h"


#include "root-addons/cbtaps_display/TH2CB.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class MesonDalitzDecays : public Physics {

public:
    struct Tree_t : WrapTTree {
        Tree_t();

        ADD_BRANCH_T(unsigned,                    nCands)

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_kinfitted, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_treefitted, 3)
        ADD_BRANCH_T(std::vector<double>,         photons_Time, 3)
        ADD_BRANCH_T(std::vector<TVector2>,       photons_PSA, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_detector, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_clusterSize, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_centralElem, 3)
        ADD_BRANCH_T(std::vector<double>,         photons_vetoE, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_vetoChannel, 3)

        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_E_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_theta_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_phi_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_E_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_theta_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_phi_pulls, 3)

        ADD_BRANCH_T(TLorentzVector,              p)
        ADD_BRANCH_T(TLorentzVector,              p_kinfitted)
        ADD_BRANCH_T(TLorentzVector,              p_treefitted)
        ADD_BRANCH_T(double,                      p_Time)
        ADD_BRANCH_T(TVector2,                    p_PSA)
        ADD_BRANCH_T(int,                         p_detector)
        ADD_BRANCH_T(int,                         p_clusterSize, 3)
        ADD_BRANCH_T(int,                         p_centralElem, 3)
        ADD_BRANCH_T(double,                      p_vetoE)
        ADD_BRANCH_T(int,                         p_vetoChannel, 3)

        ADD_BRANCH_T(double,                      p_kinfit_theta_pull)
        ADD_BRANCH_T(double,                      p_kinfit_phi_pull)
        ADD_BRANCH_T(double,                      p_treefit_theta_pull)
        ADD_BRANCH_T(double,                      p_treefit_phi_pull)

        ADD_BRANCH_T(double,                      TaggW)
        ADD_BRANCH_T(double,                      TaggW_wide)
        ADD_BRANCH_T(double,                      TaggE)
        ADD_BRANCH_T(double,                      TaggT)
        ADD_BRANCH_T(unsigned,                    TaggCh)

        ADD_BRANCH_T(double,                      beam_E_kinfitted)
        ADD_BRANCH_T(double,                      beam_E_treefitted)
        ADD_BRANCH_T(double,                      beam_kinfit_E_pull)
        ADD_BRANCH_T(double,                      beam_treefit_E_pull)
        ADD_BRANCH_T(double,                      kinfit_ZVertex)
        ADD_BRANCH_T(double,                      kinfit_ZVertex_pull)
        ADD_BRANCH_T(double,                      treefit_ZVertex)
        ADD_BRANCH_T(double,                      treefit_ZVertex_pull)

        ADD_BRANCH_T(double,                      kinfit_chi2)
        ADD_BRANCH_T(double,                      kinfit_probability)
        ADD_BRANCH_T(unsigned,                    kinfit_iterations)
        ADD_BRANCH_T(unsigned,                    kinfit_DoF)
        ADD_BRANCH_T(double,                      treefit_chi2)
        ADD_BRANCH_T(double,                      treefit_probability)
        ADD_BRANCH_T(unsigned,                    treefit_iterations)
        ADD_BRANCH_T(unsigned,                    treefit_DoF)

        ADD_BRANCH_T(double,                      CBSumE)
        ADD_BRANCH_T(double,                      CBAvgTime)

        ADD_BRANCH_T(unsigned,                    channel)
        ADD_BRANCH_T(bool,                        MCtrue)
        ADD_BRANCH_T(double,                      trueZVertex)

        ADD_BRANCH_T(TLorentzVector,              eta)
        ADD_BRANCH_T(TLorentzVector,              eta_kinfit)
        ADD_BRANCH_T(TLorentzVector,              eta_treefit)
        ADD_BRANCH_T(TLorentzVector,              mm)
        ADD_BRANCH_T(double,                      copl)
    };

protected:
    TH1D* h_tagger_time;
    TH1D* h_tagger_time_CBavg;

    TH2* h_eegPID = nullptr;
    TH2* h_eegPID_proton = nullptr;
    TH2* h_eegPID_combined = nullptr;
    TH1D* h_pTheta = nullptr;
    TH1D* h_protonVeto = nullptr;
    TH1D* h_etaIM_final = nullptr;
    TH2D* h_IM2d = nullptr;
    TH2* h_eta = nullptr;
    TH2* h_proton = nullptr;

    TH1D* h_counts = nullptr;
    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;

    static constexpr unsigned N_FINAL_STATE = 4;
    static constexpr double ETA_IM = 547.853;
    static constexpr double ETA_SIGMA = 50.;
    // cuts
    static constexpr bool PROBABILITY_CUT = false;
    static constexpr double PROBABILITY = .02;
    static constexpr bool ANTI_PI0_CUT = false;
    static constexpr double ANTI_PI0_LOW = 102.;
    static constexpr double ANTI_PI0_HIGH = 170.;
    static constexpr bool IM2D_LINEAR_CUT = false;
    static constexpr bool LEPTON_PI0_CUT = false;
    static constexpr double LEPTON_PI0_THRESH = 130.;
    // which fit should be used to determine best candidate combination?
    static constexpr bool USE_TREEFIT = false;

    struct PerChannel_t {
        std::string title;
        std::string name;
        TH2* eegPID = nullptr;
        TH1D* steps = nullptr;
        TH1D* etaIM = nullptr;
        TH1D* etaIM_kinfit = nullptr;
        TH1D* etaIM_treefit = nullptr;
        TH1D* etaIM_cand = nullptr;
        TH1D* etaIM_final = nullptr;
        TH2D* IM2d = nullptr;
        TH1D* MM = nullptr;
        TH1D* hCopl = nullptr;
        TH1D* hCopl_final = nullptr;
        TH1D* treefitChi2 = nullptr;
        TH1D* treefitProb = nullptr;
        TH1D* treefitIter = nullptr;
        TH1D* kinfitChi2 = nullptr;
        TH1D* kinfitProb = nullptr;
        TH1D* kinfitIter = nullptr;
        TH1D* effect_rad = nullptr;
        TH2D* effect_rad_E = nullptr;
        TH1D* cluster_size = nullptr;
        TH2D* cluster_size_E = nullptr;

        TH2* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Name, const std::string& Title, HistogramFactory& hf);

        void Show();
        void Fill(const TEventData& d);
    };

    std::map<std::string, PerChannel_t> channels;
    std::map<std::string, HistogramFactory&> productions;

    Tree_t t;
    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;
    using uncertainty_model_t = utils::UncertaintyModels::Optimized_Oli1;
    utils::UncertaintyModelPtr model;
    utils::KinFitter kinfit;
    utils::TreeFitter treefitter_eta;

    std::shared_ptr<ant::Detector_t> cb;

    template<typename T>
    void shift_right(std::vector<T>&);

    utils::ClusterTools clustertools;

    double effective_radius(const TCandidatePtr) const;

    ParticleTypeTree base_tree();
    ParticleTypeTree eta_3g();

    double linear_cut(const double) const;

public:

    MesonDalitzDecays(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    bool doFit_checkProb(const TTaggerHit& taggerhit,
                                const TParticlePtr proton,
                                const TParticleList photons,
                                PerChannel_t& h,
                                Tree_t& t,
                                double& best_prob_fit);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    struct ReactionChannel_t {
        std::string name = "";
        std::shared_ptr<decaytree_t> tree = nullptr;
        int color = kBlack;

        ReactionChannel_t() = default;
        ReactionChannel_t(const std::string& n);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const int c);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const std::string& n, const int c);
        ~ReactionChannel_t();
    };

    struct ReactionChannelList_t {
        static const unsigned other_index;
        std::map<unsigned, ReactionChannel_t> channels;
        unsigned identify(const TParticleTree_t &tree) const;
    };

    static ReactionChannelList_t makeChannels();
    static const ReactionChannelList_t reaction_channels;
};

}}} // namespace ant::analysis::physics
