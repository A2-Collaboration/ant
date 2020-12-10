#pragma once

#include "analysis/physics/etaprime/etaprime_dalitz.h"

#include "analysis/utils/A2GeoAcceptance.h"

namespace ant {
namespace analysis {
namespace physics {

class Etap2gMC;

class EtapDalitzMC : public Physics, public EtapDalitzTools {

public:
    using SigTree_t = EtapDalitz::SigTree_t;
    using RefTree_t = EtapDalitz::RefTree_t;

protected:
    struct MCTree_t : WrapTTree {
        ADD_BRANCH_T(std::vector<std::string>,  names)
        ADD_BRANCH_T(std::vector<double>,       energies)
        ADD_BRANCH_T(std::vector<double>,       energies_true)
        ADD_BRANCH_T(std::vector<double>,       thetas)
        ADD_BRANCH_T(std::vector<double>,       thetas_true)
        ADD_BRANCH_T(std::vector<double>,       phis)
        ADD_BRANCH_T(std::vector<double>,       phis_true)
        ADD_BRANCH_T(unsigned,                  multiplicity)
        ADD_BRANCH_T(double,                    imee)
        ADD_BRANCH_T(double,                    opening)
        ADD_BRANCH_T(unsigned,                  nCB)
        ADD_BRANCH_T(unsigned,                  nTAPS)

        void fillAndReset()
        {
            Tree->Fill();
            names().resize(0);
            energies().resize(0);
            energies_true().resize(0);
            thetas().resize(0);
            thetas_true().resize(0);
            phis().resize(0);
            phis_true().resize(0);
            multiplicity = 0;
            imee = std_ext::NaN;
            opening = std_ext::NaN;
            nCB = 0;
            nTAPS = 0;
        }
    };

    MCTree_t mc;


    double imee;

    TH1D* h_pTheta = nullptr;
    TH1D* h_protonVeto = nullptr;
    TH1D* h_etapIM_final = nullptr;
    TH2D* h_IM2d = nullptr;
    TH2* h_etap = nullptr;
    TH2* h_proton = nullptr;
    TH1D* h_subIM_2g = nullptr;
    TH1D* h_subIM_2g_fit = nullptr;
    TH1D* h_IMee_true = nullptr;

    TH2D* h_etapIM_vs_IMee = nullptr;
    TH2D* h_MM_vs_IMee = nullptr;
    TH2D* h_etapIM_fitted_vs_IMee = nullptr;
    TH1D* h_counts = nullptr;
    TH1D* h_nCands = nullptr;
    TH1D* h_cluster_CB = nullptr;
    TH1D* h_cluster_TAPS = nullptr;
    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;

    // energy and theta dependence and deviations of kinematics in MC
    TH1D* h_energy_deviation = nullptr;
    TH2D* h_fsClE_vs_pluto_geant_dE = nullptr;
    TH2D* h_theta_vs_vz = nullptr;
    TH2D* h_theta_vs_pluto_geant_dE = nullptr;
    TH2D* h_vz_vs_pluto_geant_dE = nullptr;
    TH2D* h_delta_vz_vs_pluto_geant_dE = nullptr;
    // differences/resolutions between true and reconstructed/fitted particles
    // focus on energy
    TH1D* h_energy_resolution_g = nullptr;
    TH1D* h_energy_resolution_em = nullptr;
    TH1D* h_energy_resolution_ep = nullptr;
    TH2D* h_energy_resolution_vs_theta_g = nullptr;
    TH2D* h_energy_resolution_vs_theta_em = nullptr;
    TH2D* h_energy_resolution_vs_theta_ep = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g = nullptr;
    TH2D* h_energy_resolution_vs_trueE_em = nullptr;
    TH2D* h_energy_resolution_vs_trueE_ep = nullptr;
    TH1D* h_energy_resolution_g_fit = nullptr;
    TH1D* h_energy_resolution_em_fit = nullptr;
    TH1D* h_energy_resolution_ep_fit = nullptr;
    TH2D* h_energy_resolution_vs_theta_g_fit = nullptr;
    TH2D* h_energy_resolution_vs_theta_em_fit = nullptr;
    TH2D* h_energy_resolution_vs_theta_ep_fit = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g_fit = nullptr;
    TH2D* h_energy_resolution_vs_trueE_em_fit = nullptr;
    TH2D* h_energy_resolution_vs_trueE_ep_fit = nullptr;
    // focus on theta
    TH1D* h_theta_resolution_g = nullptr;
    TH1D* h_theta_resolution_em = nullptr;
    TH1D* h_theta_resolution_ep = nullptr;
    TH2D* h_theta_resolution_vs_energy_g = nullptr;
    TH2D* h_theta_resolution_vs_energy_em = nullptr;
    TH2D* h_theta_resolution_vs_energy_ep = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_em = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_ep = nullptr;
    TH1D* h_theta_resolution_g_fit = nullptr;
    TH1D* h_theta_resolution_em_fit = nullptr;
    TH1D* h_theta_resolution_ep_fit = nullptr;
    TH2D* h_theta_resolution_vs_energy_g_fit = nullptr;
    TH2D* h_theta_resolution_vs_energy_em_fit = nullptr;
    TH2D* h_theta_resolution_vs_energy_ep_fit = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g_fit = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_em_fit = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_ep_fit = nullptr;

    // some more histograms to check certain kinematic conditions,
    // especially in different dilepton mass regions
    TH1D* h_IMee = nullptr;
    TH1D* h_IMee_total = nullptr;
    TH1D* h_IMee_acceptance = nullptr;
    TH1D* h_IMee_Trigger4Cl = nullptr;
    TH1D* h_IMee_fraction3CB1TAPS_total = nullptr;
    TH1D* h_IMee_fraction3CB1TAPS_acceptance = nullptr;
    TH1D* h_IMee_fraction3CB1TAPS_Trigger4Cl = nullptr;
    TH2D* h_nCands_vs_IMee = nullptr;
    TH2D* h_openingAngle_vs_IMee = nullptr;
    TH1D* h_CBEsum = nullptr;
    TH1D* h_CBEsum_true = nullptr;
    TH2D* h_CBEsum_vs_IMee = nullptr;
    TH2D* h_E_vs_IMee_eCharged_true = nullptr;
    TH2D* h_E_vs_IMee_photon_true = nullptr;
    TH2D* h_E_vs_IMee_proton_true = nullptr;
    TH2D* h_E_vs_IMee_eCharged_rec = nullptr;
    TH2D* h_E_vs_IMee_photon_rec = nullptr;
    TH2D* h_E_vs_IMee_proton_rec = nullptr;

    // test histograms for some checks related to radiative corrections
    TH1D* h_radCorr_x = nullptr;
    TH1D* h_radCorr_y = nullptr;
    TH2D* h_radCorr_checkBoundaries = nullptr;
    TH2D* h_radCorr_y_vs_IMee = nullptr;
    TH2D* h_radCorr_y_vs_x = nullptr;

    // test for cone prediction of proton candidate
    TH1D* h_theta_miss_res = nullptr;
    TH1D* h_angle_miss_res = nullptr;
    TH2D* h_theta_vs_phi_miss_res = nullptr;


    using Cuts_t = EtapDalitz::Cuts_t;
    using Settings_t = EtapDalitz::Settings_t;

    Settings_t& settings = Settings_t::get();

    EtapDalitzTools tools;

    struct PerChannel_t {
        std::string title;
        std::string name;
        TH1D* steps = nullptr;
        TH2D* steps_vs_IMee = nullptr;
        TH1D* etapIM = nullptr;
        TH1D* etapIM_kinfit = nullptr;
        TH1D* etapIM_kinfit_freeZ = nullptr;
        TH1D* etapIM_treefit = nullptr;
        TH1D* etapIM_treefit_freeZ = nullptr;
        TH1D* etapIM_cand = nullptr;
        TH1D* etapIM_final = nullptr;
        TH2D* IM2d = nullptr;
        TH1D* MM = nullptr;
        TH1D* hCopl = nullptr;
        TH1D* hCopl_final = nullptr;
        TH1D* trueZVertex = nullptr;
        TH2D* true_rec_ZVertex = nullptr;
        TH1D* treefitChi2 = nullptr;
        TH1D* treefitProb = nullptr;
        TH1D* treefitIter = nullptr;
        TH1D* treefit_freeZ_chi2 = nullptr;
        TH1D* treefit_freeZ_prob = nullptr;
        TH1D* treefit_freeZ_iter = nullptr;
        TH1D* kinfitChi2 = nullptr;
        TH1D* kinfitProb = nullptr;
        TH1D* kinfitIter = nullptr;
        TH1D* kinfit_freeZ_chi2 = nullptr;
        TH1D* kinfit_freeZ_prob = nullptr;
        TH1D* kinfit_freeZ_iter = nullptr;
        TH1D* kinfit_ZVertex = nullptr;
        TH1D* kinfit_freeZ_ZVertex = nullptr;
        TH1D* treefit_ZVertex = nullptr;
        TH1D* treefit_freeZ_ZVertex = nullptr;
        TH1D* antiPionProb = nullptr;
        TH1D* effect_rad = nullptr;
        TH2D* effect_rad_E = nullptr;
        TH1D* cluster_size = nullptr;
        TH2D* cluster_size_E = nullptr;
        TH1D* lat_moment = nullptr;
        TH2D* lat_moment_E = nullptr;

        TH2* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Name, const std::string& Title, HistogramFactory& hf);

        void Show();
        void Fill(const TEventData& d);
    };

    channel_id_t chan_id;
    std::map<std::string, PerChannel_t> channels;
    std::map<std::string, HistogramFactory&> productions;

    SigTree_t sig;
    RefTree_t ref;

    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    utils::A2SimpleGeometry geo;
    mev_t calcEnergySum(const TParticleList&) const;
    TParticleList getGeoAccepted(const TParticleList&) const;
    TParticleList getGeoAcceptedDetector(const TParticleList&, const Detector_t::Type_t) const;
    template <typename T, typename = std::enable_if<std::is_integral<T>::value>>
    T geoAccepted(const TParticleList&) const;
    template <typename T, typename = std::enable_if<std::is_integral<T>::value>>
    T geoAccepted(const TCandidateList&) const;
    template <typename T, typename = std::enable_if<std::is_integral<T>::value>>
    T geoAcceptedDetector(const TParticleList&, const Detector_t::Type_t) const;
    template <typename T, typename = std::enable_if<std::is_integral<T>::value>>
    T geoAcceptedDetector(const TCandidateList&, const Detector_t::Type_t) const;
    size_t geoAccepted(const TParticleList& p) const { return geoAccepted<size_t>(p); }
    size_t geoAcceptedDetector(const TParticleList& p, const Detector_t::Type_t d) const
    { return geoAcceptedDetector<size_t>(p, d); }

    utils::UncertaintyModelPtr model_MC;

    utils::KinFitter kinfit;
    utils::KinFitter kinfit_freeZ;
    utils::TreeFitter treefitter_etap;
    utils::TreeFitter treefitter_etap_freeZ;

    //using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_comb_t = fake_comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    std::shared_ptr<ant::Detector_t> cb;

    PerChannel_t manage_channel_histograms_get_current(const bool, const TEvent&);

    Etap2gMC* etap2g;

public:

    EtapDalitzMC(const std::string& name, OptionsPtr opts);

    bool doFit_checkProb(const TTaggerHit& taggerhit,
                         const particle_comb_t& comb,
                         PerChannel_t& h,
                         SigTree_t& t,
                         double& best_prob_fit);
    double anti_pion_fit(const TTaggerHit& taggerhit, const particle_comb_t& comb);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    using ReactionChannel_t = EtapDalitz::ReactionChannel_t;
    using ReactionChannelList_t = EtapDalitz::ReactionChannelList_t;

    static ReactionChannelList_t makeChannels() { return EtapDalitz::makeChannels(); }
    static const ReactionChannelList_t reaction_channels;
};


class Etap2gMC : public Physics, public EtapDalitzTools {

protected:
    struct MCTree_t : WrapTTree {
        ADD_BRANCH_T(std::vector<std::string>,  names)
        ADD_BRANCH_T(std::vector<double>,       energies)
        ADD_BRANCH_T(std::vector<double>,       energies_true)
        ADD_BRANCH_T(std::vector<double>,       thetas)
        ADD_BRANCH_T(std::vector<double>,       thetas_true)
        ADD_BRANCH_T(std::vector<double>,       phis)
        ADD_BRANCH_T(std::vector<double>,       phis_true)
        ADD_BRANCH_T(unsigned,                  multiplicity)
        ADD_BRANCH_T(double,                    opening)
        ADD_BRANCH_T(unsigned,                  nCB)
        ADD_BRANCH_T(unsigned,                  nTAPS)

        void fillAndReset()
        {
            Tree->Fill();
            names().resize(0);
            energies().resize(0);
            energies_true().resize(0);
            thetas().resize(0);
            thetas_true().resize(0);
            phis().resize(0);
            phis_true().resize(0);
            multiplicity = 0;
            opening = std_ext::NaN;
            nCB = 0;
            nTAPS = 0;
        }
    };

    MCTree_t mc;

    const bool less_plots;


    TH2D* h_taggChannel_vs_trueIM = nullptr;

    // differences/resolutions between true and reconstructed/fitted particles
    // focus on energy
    TH2D* h_energy_resolution_vs_theta_g1 = nullptr;
    TH2D* h_energy_resolution_vs_theta_g2 = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g1 = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g2 = nullptr;
    TH2D* h_energy_resolution_vs_theta_g1_fit = nullptr;
    TH2D* h_energy_resolution_vs_theta_g2_fit = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g1_fit = nullptr;
    TH2D* h_energy_resolution_vs_trueE_g2_fit = nullptr;
    // focus on theta
    TH2D* h_theta_resolution_vs_energy_g1 = nullptr;
    TH2D* h_theta_resolution_vs_energy_g2 = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g1 = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g2 = nullptr;
    TH2D* h_theta_resolution_vs_energy_g1_fit = nullptr;
    TH2D* h_theta_resolution_vs_energy_g2_fit = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g1_fit = nullptr;
    TH2D* h_theta_resolution_vs_trueTheta_g2_fit = nullptr;

    PromptRandom::Switch* promptrandom;
    utils::TriggerSimulation triggersimu;

    utils::UncertaintyModelPtr model_MC;

    utils::KinFitter kinfit;
    utils::TreeFitter treefitter_etap;

    std::shared_ptr<ant::Detector_t> ept;

    using RefTree_t = EtapDalitzMC::RefTree_t;

    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    RefTree_t* t;

    void fill_tree(const APLCON::Result_t&,
                   const APLCON::Result_t&,
                   const TParticlePtr,
                   const TParticleList&);
    bool simple2CB1TAPS(const TCandidateList& cands,
                        TParticlePtr& proton,
                        TParticleList& photons);
    bool doFit_checkProb(const TTaggerHit& taggerhit,
                         const TParticlePtr proton,
                         const TParticleList& photons,
                         double& best_prob_fit);

public:
    Etap2gMC(const std::string& name, OptionsPtr opts);

    void setPromptRandom(PromptRandom::Switch&);
    void linkTree(RefTree_t&);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    void Process(const TEvent& event);
};

}}} // namespace ant::analysis::physics
