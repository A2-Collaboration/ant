#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/A2GeoAcceptance.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "analysis/utils/particle_tools.h"
#include "base/std_ext/math.h"
#include "base/WrapTTree.h"

#include "analysis/utils/Fitter.h"
#include "base/interval.h"
#include "analysis/plot/PromptRandomHist.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Rtypes.h"

#include <map>


class TH1D;
class TH2D;
class TH3D;


namespace ant {

namespace analysis {
namespace physics {

class OmegaMCTruePlots: public Physics {
public:
    struct PerChannel_t {
        std::string title;
        TH2D* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Title, HistogramFactory& hf);

        void Show();
        void Fill(const TEventData& d);
    };

    std::map<std::string,PerChannel_t> channels;

    OmegaMCTruePlots(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

class OmegaBase: public Physics {

public:
    enum class DataMode {
        MCTrue,
        Reconstructed
    };

protected:
    utils::A2SimpleGeometry geo;
    double calcEnergySum(const TParticleList& particles) const;
    TParticleList getGeoAccepted(const TParticleList& p) const;
    unsigned geoAccepted(const TCandidateList& cands) const;

    DataMode mode = DataMode::Reconstructed;

    virtual void Analyse(const TEventData& data, const TEvent& event, manager_t& manager) =0;



public:
    OmegaBase(const std::string &name, OptionsPtr opts);
    virtual ~OmegaBase() = default;

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;


};

class OmegaEtaG: public OmegaBase {

protected:

    TH2D* ggg_gg;
    TH2D* ggg_gg_bg;    // if not from omega decay
    TH2D* ggg_gg_all;
    //TH1D* gg_bg;        // if from omega but not pi0 or eta decay

    TH1D* ggg;
    TH1D* ggg_omega;
    TH1D* ggg_bg;
    TH1D* ggg_omega_pi0oreta;

    TH2D* ggg_gg_omega_eta;
    TH2D* ggg_gg_omega_pi0;

    TH1D* steps;

    IntervalD omega_range = IntervalD(680,780);

    struct perDecayhists_t {
        TH1D* gg = nullptr;
        TH1D* ggg = nullptr;
        TH1D* mm = nullptr;
        TH1D* angle_p;
        TH1D* angle_p_ggg;
        TH1D* p_phi_diff;
        TH2D* calc_proton_energy_theta;
        TH2D* calc_proton_special;
        TH1D* nCand;
    };

    perDecayhists_t makePerDecayHists(const std::string &title="");

    std::map<std::string, perDecayhists_t> gg_decays;

    virtual void Analyse(const TEventData& data, const TEvent& event, manager_t&) override;

    BinSettings imbinning = BinSettings(1000);
    BinSettings mmbinning = BinSettings(1000, 400,1400);

public:
    OmegaEtaG(const std::string& name, OptionsPtr opts);
    virtual ~OmegaEtaG() = default;
    virtual void ShowResult() override;
};



class OmegaMCTree : public Physics {
protected:
    TTree* tree = nullptr;
    TLorentzVector proton_vector;
    TLorentzVector omega_vector;
    TLorentzVector eta_vector;
    TLorentzVector gamma1_vector;
    TLorentzVector gamma2_vector;
    TLorentzVector gamma3_vector;

public:

    OmegaMCTree(const std::string& name, OptionsPtr opts);
    virtual ~OmegaMCTree();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    LorentzVec getGamma1() const;
    void setGamma1(const LorentzVec& value);
};


class OmegaEtaG2 : public OmegaBase {
public:
    struct OmegaTree_t : WrapTTree {
        OmegaTree_t();

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_fitted, 3)
        ADD_BRANCH_T(TLorentzVector,              p)
        ADD_BRANCH_T(double,                      p_Time)
        ADD_BRANCH_T(double,                      p_PSA_Angle)
        ADD_BRANCH_T(double,                      p_PSA_Radius)
        ADD_BRANCH_T(int,                         p_detector)

        ADD_BRANCH_T(TLorentzVector,              p_true)
        ADD_BRANCH_T(TLorentzVector,              p_fitted)

        ADD_BRANCH_T(TLorentzVector,              ggg)
        ADD_BRANCH_T(TLorentzVector,              ggg_fitted)
        ADD_BRANCH_T(TLorentzVector,              mm)
        ADD_BRANCH_T(double,                      copl_angle)
        ADD_BRANCH_T(double,                      p_mm_angle)

        ADD_BRANCH_T(std::vector<double>,         ggIM, 3)
        ADD_BRANCH_T(std::vector<double>,         ggIM_fitted, 3)

        ADD_BRANCH_T(std::vector<double>,         BachelorE, 3)
        ADD_BRANCH_T(std::vector<double>,         BachelorE_fitted, 3)

        ADD_BRANCH_T(double,                      ggIM_real)  // only if Signal/Ref
        ADD_BRANCH_T(std::vector<double>,         ggIM_comb, 2)  // only if Signal/Ref

        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(double,   TaggW_tight)
        ADD_BRANCH_T(double,   TaggE)
        ADD_BRANCH_T(double,   TaggT)
        ADD_BRANCH_T(unsigned, TaggCh)

        ADD_BRANCH_T(double,   KinFitChi2)
        ADD_BRANCH_T(unsigned, KinFitIterations)

        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)

        ADD_BRANCH_T(bool,     p_matched)

        ADD_BRANCH_T(unsigned,      Channel)

        ADD_BRANCH_T(std::vector<double>,       pi0chi2, 3)
        ADD_BRANCH_T(int,                        iBestPi0)
        ADD_BRANCH_T(std::vector<double>,       etachi2, 3)
        ADD_BRANCH_T(int,                        iBestEta)

        ADD_BRANCH_T(int,                        bestHyp)

    };

protected:
    void Analyse(const TEventData &data, const TEvent& event, manager_t&manager) override;


    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;


    //======== Tree ===============================================================

    TTree*  tree = nullptr;
    OmegaTree_t t;

    static const std::vector<std::vector<std::size_t>> combs;




    //======== Settings ===========================================================

    const double cut_ESum;
    const double cut_Copl;
    const interval<double> photon_E_cb;
    const interval<double> photon_E_taps;
    const interval<double> proton_theta;


    ant::analysis::PromptRandom::Switch promptrandom;

    std::shared_ptr<utils::Fitter::UncertaintyModel> model;

    utils::KinFitter fitter;

    struct MyTreeFitter_t {
        utils::TreeFitter treefitter;
        utils::TreeFitter::tree_t fitted_Omega;
        utils::TreeFitter::tree_t fitted_g_Omega;
        utils::TreeFitter::tree_t fitted_X;
        utils::TreeFitter::tree_t fitted_g1_X;
        utils::TreeFitter::tree_t fitted_g2_X;

        MyTreeFitter_t(const ParticleTypeTree& ttree, const ParticleTypeDatabase::Type& mesonT, const std::shared_ptr<const utils::Fitter::UncertaintyModel>& model);
    };
    static int CombIndex(const ant::TParticleList& orig, const MyTreeFitter_t&);

    MyTreeFitter_t fitter_pi0;
    MyTreeFitter_t fitter_eta;

    bool AcceptedPhoton(const TParticlePtr& photon);
    bool AcceptedProton(const TParticlePtr& proton);

    TParticleList FilterPhotons(const TParticleList& list);
    TParticleList FilterProtons(const TParticleList& list);

    bool opt_save_after_kinfit = false;

public:

    OmegaEtaG2(const std::string& name, OptionsPtr opts);
    virtual ~OmegaEtaG2();

    void Finish() override;

    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    struct ReactionChannel_t {
        std::string name="";
        std::shared_ptr<decaytree_t> tree=nullptr;
        int color=kBlack;
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const std::string& n, const int c);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const int c);
        ReactionChannel_t(const std::string& n);
        ReactionChannel_t() = default;
        ~ReactionChannel_t();
    };

    struct ReactionChannelList_t {
        static const unsigned other_index;
        std::map<unsigned, ReactionChannel_t> channels;
        unsigned identify(const TParticleTree_t &tree) const;
    };

    static ReactionChannelList_t makeChannels();
    static const ReactionChannelList_t reaction_channels;

    std::map<unsigned, TH1D*> stephists;

};

}
}
}

std::string to_string(const ant::analysis::physics::OmegaBase::DataMode& m);
