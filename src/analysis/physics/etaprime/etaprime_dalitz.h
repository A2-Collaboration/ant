#pragma once

#include <vector>
#include <list>
#include <type_traits>
#include <random>
#include <mutex>

#include "analysis/physics/Physics.h"
#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/utils/ClusterTools.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"
#include "analysis/utils/ProtonPhotonCombs.h"
#include "base/WrapTTree.h"
#include "base/std_ext/vector.h"
#include "tree/TSimpleParticle.h"

#include "root-addons/cbtaps_display/TH2CB.h"

#include "TLorentzVector.h"

namespace ant {

namespace expconfig {
namespace detector {
    struct TAPS;
}}

namespace analysis {
namespace physics {

class Etap2g;

// some tools to access from several classes (as well as nested ones)
struct EtapDalitzTools {

    struct fake_comb_t {
        fake_comb_t() = default;
        TParticleList Photons;
        TParticlePtr  Proton;

        LorentzVec PhotonSum{{0,0,0},0};
        double MissingMass{std_ext::NaN};
        double DiscardedEk{0.};

        void reset();
        void calc_values(const TTaggerHit& taggerhit);
    };

    double effective_radius(const TCandidatePtr) const;
    double lat_moment(const TCandidatePtr) const;

    template <typename iter>
    LorentzVec sumlv(iter start, iter end) {
        LorentzVec s;
        while (start != end) {
            s += **(start);
            ++start;
        }
        return s;
    }

    void count_clusters(const TCandidateList&, size_t&, size_t&);
    bool q2_preselection(const TEventData&, const double) const;

protected:
    utils::ClusterTools clustertools;

    struct channel_id_t {
        // identifier for current channel information
        std::string production;
        std::string decaystring;
        std::string decay_name;
    };

public:
    void channel_id(const TEvent&, channel_id_t&);

    ParticleTypeTree base_tree();
    ParticleTypeTree etap_3g();

    /* helper methods to handle veto energies; implemented in header to avoid linker problems without modifying CMake files */

    template <typename Particle>
    static std::vector<size_t> get_sorted_indices_vetoE(std::vector<Particle> particles)
    {
        return std_ext::get_sorted_indices_desc(get_veto_energies(particles));
    }

    static std::vector<double> get_veto_energies(std::vector<TSimpleParticle> particles)
    {
        std::vector<double> veto_energies;
        for (const auto& p : particles)
            veto_energies.emplace_back(p.VetoE);

        return veto_energies;
    }

    static std::vector<double> get_veto_energies(std::vector<TParticlePtr> particles)
    {
        std::vector<double> veto_energies;
        for (const auto& p : particles)
            veto_energies.emplace_back(p->Candidate->VetoEnergy);

        return veto_energies;
    }
};

class EtapDalitz : public Physics, public EtapDalitzTools {

    friend class EtapDalitzMC;

public:
    struct common_tree : virtual WrapTTree {
        ADD_BRANCH_T(unsigned, beamtime, 0)
        ADD_BRANCH_T(unsigned, nCands)
        ADD_BRANCH_T(unsigned, channel)
        ADD_BRANCH_T(bool,     MCtrue)
        ADD_BRANCH_T(double,   trueZVertex)

        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)

        void init();
        void reset();
    };

    struct proton_tree : virtual WrapTTree {
        ADD_BRANCH_T(TSimpleParticle, p)
        ADD_BRANCH_T(double,          p_PSAangle)
        ADD_BRANCH_T(double,          p_PSAradius)
        ADD_BRANCH_T(int,             p_detector)
        ADD_BRANCH_T(unsigned,        p_centralElem)
        ADD_BRANCH_T(unsigned,        p_vetoChannel)
        ADD_BRANCH_T(double,          p_vetoTime)

        void reset();
        void set_proton_information(const TParticlePtr);
    };

    template <size_t N>
    struct photon_tree : virtual WrapTTree, public EtapDalitzTools {
        ADD_BRANCH_T(std::vector<TSimpleParticle>, photons, N)
        ADD_BRANCH_T(std::vector<int>,             photons_detector, N)
        ADD_BRANCH_T(std::vector<unsigned>,        photons_centralElem, N)
        ADD_BRANCH_T(std::vector<unsigned>,        photons_vetoChannel, N)
        ADD_BRANCH_T(std::vector<double>,          photons_vetoTime, N)
        ADD_BRANCH_T(std::vector<double>,          photons_effect_radius, N)
        ADD_BRANCH_T(std::vector<double>,          photons_lat_moment, N)

        void reset();
        void set_photon_information(const TParticleList&, const bool shower_shape = false);
    };

    template <size_t Nphotons>
    struct fit_tree : virtual WrapTTree {
        ADD_BRANCH_T(double,                       beam_E_kinfitted)
        ADD_BRANCH_T(double,                       beam_E_treefitted)
        ADD_BRANCH_T(double,                       beam_kinfit_E_pull)
        ADD_BRANCH_T(double,                       beam_treefit_E_pull)
        ADD_BRANCH_T(double,                       kinfit_ZVertex)
        ADD_BRANCH_T(double,                       kinfit_ZVertex_pull)
        ADD_BRANCH_T(double,                       treefit_ZVertex)
        ADD_BRANCH_T(double,                       treefit_ZVertex_pull)

        ADD_BRANCH_T(TLorentzVector,               p_kinfitted)
        ADD_BRANCH_T(TLorentzVector,               p_treefitted)

        ADD_BRANCH_T(double,                       p_kinfit_E_pull)
        ADD_BRANCH_T(double,                       p_kinfit_theta_pull)
        ADD_BRANCH_T(double,                       p_kinfit_phi_pull)
        ADD_BRANCH_T(double,                       p_treefit_E_pull)
        ADD_BRANCH_T(double,                       p_treefit_theta_pull)
        ADD_BRANCH_T(double,                       p_treefit_phi_pull)

        ADD_BRANCH_T(std::vector<TLorentzVector>,  photons_kinfitted, Nphotons)
        ADD_BRANCH_T(std::vector<TLorentzVector>,  photons_treefitted, Nphotons)

        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_E_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_theta_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_phi_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_E_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_theta_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_phi_pulls, Nphotons)

        ADD_BRANCH_T(double,                       kinfit_chi2)
        ADD_BRANCH_T(double,                       kinfit_probability)
        ADD_BRANCH_T(int,                          kinfit_iterations)
        ADD_BRANCH_T(int,                          kinfit_DoF)
        ADD_BRANCH_T(double,                       treefit_chi2)
        ADD_BRANCH_T(double,                       treefit_probability)
        ADD_BRANCH_T(int,                          treefit_iterations)
        ADD_BRANCH_T(int,                          treefit_DoF)

        ADD_BRANCH_T(TLorentzVector,               etap_kinfit)
        ADD_BRANCH_T(TLorentzVector,               etap_treefit)

        void reset();
        void set_kinfit_information(const utils::KinFitter&, const APLCON::Result_t&);
        void set_treefit_information(const utils::TreeFitter&, const APLCON::Result_t&);
    };

    template <size_t Nphotons>
    struct fit_freeZ_tree : virtual WrapTTree {
        ADD_BRANCH_T(double,                       beam_E_kinfit_freeZ)
        ADD_BRANCH_T(double,                       beam_E_treefit_freeZ)
        ADD_BRANCH_T(double,                       beam_kinfit_freeZ_E_pull)
        ADD_BRANCH_T(double,                       beam_treefit_freeZ_E_pull)
        ADD_BRANCH_T(double,                       kinfit_freeZ_ZVertex)
        ADD_BRANCH_T(double,                       kinfit_freeZ_ZVertex_pull)
        ADD_BRANCH_T(double,                       treefit_freeZ_ZVertex)
        ADD_BRANCH_T(double,                       treefit_freeZ_ZVertex_pull)

        ADD_BRANCH_T(TLorentzVector,               p_kinfit_freeZ)
        ADD_BRANCH_T(TLorentzVector,               p_treefit_freeZ)

        ADD_BRANCH_T(double,                       p_kinfit_freeZ_theta_pull)
        ADD_BRANCH_T(double,                       p_kinfit_freeZ_phi_pull)
        ADD_BRANCH_T(double,                       p_treefit_freeZ_theta_pull)
        ADD_BRANCH_T(double,                       p_treefit_freeZ_phi_pull)

        ADD_BRANCH_T(std::vector<TLorentzVector>,  photons_kinfit_freeZ, Nphotons)
        ADD_BRANCH_T(std::vector<TLorentzVector>,  photons_treefit_freeZ, Nphotons)

        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_freeZ_E_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_freeZ_theta_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_kinfit_freeZ_phi_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_freeZ_E_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_freeZ_theta_pulls, Nphotons)
        ADD_BRANCH_T(std::vector<double>,          photon_treefit_freeZ_phi_pulls, Nphotons)

        ADD_BRANCH_T(double,                       kinfit_freeZ_chi2)
        ADD_BRANCH_T(double,                       kinfit_freeZ_probability)
        ADD_BRANCH_T(int,                          kinfit_freeZ_iterations)
        ADD_BRANCH_T(int,                          kinfit_freeZ_DoF)
        ADD_BRANCH_T(double,                       treefit_freeZ_chi2)
        ADD_BRANCH_T(double,                       treefit_freeZ_probability)
        ADD_BRANCH_T(int,                          treefit_freeZ_iterations)
        ADD_BRANCH_T(int,                          treefit_freeZ_DoF)

        ADD_BRANCH_T(TLorentzVector,               etap_kinfit_freeZ)
        ADD_BRANCH_T(TLorentzVector,               etap_treefit_freeZ)

        void reset();
        void set_fit_freeZ_results(const utils::KinFitter&, const utils::TreeFitter&,
                                   const APLCON::Result_t&, const APLCON::Result_t&);
    };

    struct SigTree_t : common_tree, proton_tree, photon_tree<3>, fit_tree<3>, fit_freeZ_tree<3> {
        ADD_BRANCH_T(std::vector<double>,          photons_PSAangle, 3)
        ADD_BRANCH_T(std::vector<double>,          photons_PSAradius, 3)

        ADD_BRANCH_T(double,                       p_effect_radius)
        ADD_BRANCH_T(double,                       p_lat_moment)
        ADD_BRANCH_T(double,                       DiscardedEk)

        ADD_BRANCH_T(double,                       prob_antiPionFit)
        ADD_BRANCH_T(double,                       p_predictedEnergy)

        ADD_BRANCH_T(TLorentzVector,               etap)
        ADD_BRANCH_T(double,                       mm)
        ADD_BRANCH_T(double,                       copl)

        void init();
        void reset();
        void set_additional_photon_information(const TParticleList&);
    };

    struct RefTree_t : common_tree, proton_tree, photon_tree<2>, fit_tree<2> {
        ADD_BRANCH_T(TLorentzVector, etap)

        void init();
        void reset();
    };

protected:
    TH1D* h_tagger_time = nullptr;
    TH1D* h_tagger_time_CBavg = nullptr;

    TH1D* h_pTheta = nullptr;
    TH1D* h_protonVeto = nullptr;
    TH1D* h_etapIM_final = nullptr;
    TH2D* h_IM2d = nullptr;
    TH2* h_etap = nullptr;
    TH2* h_proton = nullptr;
    TH1D* h_subIM_2g = nullptr;
    TH1D* h_subIM_2g_fit = nullptr;

    TH1D* h_counts = nullptr;
    TH1D* h_nCands = nullptr;
    TH1D* h_cluster_CB = nullptr;
    TH1D* h_cluster_TAPS = nullptr;
    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;

    struct Cuts_t {
        static constexpr unsigned N_FINAL_STATE = 4;
        // threshold to check if double value should be treated as zero
        static constexpr double EPSILON = 2*std::numeric_limits<double>::epsilon();
        // cuts
        static constexpr double ETAP_SIGMA = 50.;
        static constexpr bool Q2_PRESELECTION = false;
        static constexpr double Q2_MIN_VALUE = 50.;
        static constexpr bool PROBABILITY_CUT = true;
        static constexpr double PROBABILITY = .02;
        static constexpr bool ANTI_PI0_CUT = false;
        static constexpr double ANTI_PI0_LOW = 102.;
        static constexpr double ANTI_PI0_HIGH = 170.;
    };

    // use singleton pattern for settings
    struct Settings_t {
        virtual ~Settings_t() = default;
        Settings_t(Settings_t const&) = delete;             // copy ctr
        Settings_t(Settings_t&&) = delete;                  // move ctr
        Settings_t& operator=(Settings_t const&) = delete;  // copy assign
        Settings_t& operator=(Settings_t&&) = delete;       // move assign

        static Settings_t& get()
        {
            static Settings_t instance;
            return instance;
        }

        void init(const bool _ref = false,
                  const bool _ref_only = false,
                  const bool _lessplots = false) noexcept
        {
            auto& instance = Settings_t::get();
            std::call_once(initialized, [&_ref, &_ref_only, &_lessplots, &instance] () {
                    instance.ref = _ref;
                    instance.ref_only = _ref_only;
                    instance.lessplots = _lessplots;
            });
        }

        bool reference() const { return ref || ref_only; }
        bool reference_only() const { return ref_only; }
        bool less_plots() const { return lessplots; }

        const size_t n_final_state = 4;
        const size_t n_final_state_etap = 3;
        const double max_proton_theta = 25.;
        const double max_discarded_energy = 0.;
        const double mm_window_size = 800.;
        // which fit should be used to determine best candidate combination?
        const bool use_treefit = false;

    private:
        Settings_t() = default;

        static std::once_flag initialized;

        // should the reference channel be analysed?
        bool ref = false;
        bool ref_only = false;
        // produce less plots
        bool lessplots = false;
    };

    Settings_t& settings = Settings_t::get();

    struct PerChannel_t {
        std::string title;
        std::string name;
        TH1D* steps = nullptr;
        TH1D* etapIM = nullptr;
        TH1D* etapIM_final = nullptr;
        TH2D* IM2d = nullptr;
        TH1D* MM = nullptr;
        TH1D* trueZVertex = nullptr;
        TH1D* treefitProb = nullptr;
        TH1D* treefit_freeZ_prob = nullptr;
        TH1D* kinfitProb = nullptr;
        TH1D* kinfit_freeZ_prob = nullptr;
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

    utils::UncertaintyModelPtr model_sergey;
    utils::UncertaintyModelPtr model_data;
    utils::UncertaintyModelPtr model_MC;
    utils::UncertaintyModelPtr model_data_protonMeasured;
    utils::UncertaintyModelPtr model_MC_protonMeasured;

    utils::KinFitter kinfit;
    utils::KinFitter kinfit_freeZ;
    utils::TreeFitter treefitter_etap;
    utils::TreeFitter treefitter_etap_freeZ;

    //using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_comb_t = fake_comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    std::shared_ptr<ant::Detector_t> cb;
    std::shared_ptr<expconfig::detector::TAPS> taps;

    PerChannel_t manage_channel_histograms_get_current(const bool, const TEvent&);

    void set_beamtime(common_tree*);

    Etap2g* etap2g;

public:

    EtapDalitz(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    bool doFit_checkProb(const TTaggerHit& taggerhit,
                         const particle_comb_t& comb,
                         PerChannel_t& h,
                         SigTree_t& t,
                         double& best_prob_fit);
    double anti_pion_fit(const TTaggerHit& taggerhit, const particle_comb_t& comb);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    struct ReactionChannel_t {
        std::string name = "";
        std::shared_ptr<decaytree_t> tree = nullptr;
        short color = kBlack;

        ReactionChannel_t() = default;
        ReactionChannel_t(const std::string& n);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const short c);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const std::string& n, const short c);
        ~ReactionChannel_t() = default;
    };

    struct ReactionChannelList_t {
        static const unsigned other_index;
        std::map<unsigned, ReactionChannel_t> channels;
        unsigned identify(const TParticleTree_t &tree) const;
    };

    static ReactionChannelList_t makeChannels();
    static const ReactionChannelList_t reaction_channels;
};


class Etap2g : public Physics, public EtapDalitzTools {

    friend class Etap2gMC;

protected:
    static constexpr unsigned N_FINAL_STATE = 3;
    // which method should be used to determine the proton?
    static constexpr bool USE_KINFIT = true;
    // which fit should be used to determine best candidate combination?
    // (see definition of EtapDalitz class for consistency)
    static constexpr bool USE_TREEFIT = false;
    // test combinations of proton and photons to find best proton candidate?
    static constexpr bool TEST_COMBS = true;

    PromptRandom::Switch* promptrandom;
    utils::TriggerSimulation triggersimu;

    utils::UncertaintyModelPtr model_data;
    utils::UncertaintyModelPtr model_MC;

    utils::KinFitter kinfit;
    utils::TreeFitter treefitter_etap;

    using RefTree_t = EtapDalitz::RefTree_t;

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
    Etap2g(const std::string& name, OptionsPtr opts);

    void setPromptRandom(PromptRandom::Switch&);
    void linkTree(RefTree_t&);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    void Process(const TEvent& event);
};

}}} // namespace ant::analysis::physics
