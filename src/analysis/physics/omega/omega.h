#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/A2GeoAcceptance.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/Fitter.h"
#include "analysis/plot/PromptRandomHist.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "base/std_ext/math.h"
#include "base/WrapTTree.h"

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

    static double getTime(const TParticlePtr& p) {
        return p->Candidate != nullptr ? p->Candidate->Time : std_ext::NaN;
    }

    template <typename it_type>
    static LorentzVec LVSum(it_type begin, it_type end) {
        LorentzVec v;

        while(begin!=end) {
            v += **begin;
            ++begin;
        }

        return v;
    }

    template <typename it_type>
    static LorentzVec LVSumL(it_type begin, it_type end) {
        LorentzVec v;

        while(begin!=end) {
            v += *begin;
            ++begin;
        }

        return v;
    }

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

struct TagChMultiplicity {
    TH1D* hTagChMult;
    unsigned nchannels;

    TagChMultiplicity(HistogramFactory& hf);

    void Fill(const std::vector<TTaggerHit>& t);
};

/**
 * @brief Hacked up version of the Kin fitter that allows access to the fit particles
 */
struct AccessibleFitter : utils::KinFitter {

    using KinFitter::KinFitter;

    KinFitter::FitParticle& FitPhotons(std::size_t n) {
        return *(Photons.at(n));
    }

    void SetPhoton(size_t i, const TParticlePtr& p) {
        SetPhotonEkThetaPhi(*Photons.at(i), p);
    }
};

class OmegaEtaG2 : public OmegaBase {
public:
    struct OmegaTree_t : WrapTTree {
        OmegaTree_t();

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_fitted, 3)

//        ADD_BRANCH_T(std::vector<double>,         photon_E_pulls, 3)
//        ADD_BRANCH_T(std::vector<double>,         photon_theta_pulls, 3)
//        ADD_BRANCH_T(std::vector<double>,         photon_phi_pulls, 3)
//        ADD_BRANCH_T(double,                      beam_E_pull)
//        ADD_BRANCH_T(double,                      p_theta_pull)
//        ADD_BRANCH_T(double,                      p_phi_pull)

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
        ADD_BRANCH_T(double,   KinFitProb)
        ADD_BRANCH_T(unsigned, KinFitIterations)

        ADD_BRANCH_T(double,   CBSumE)
        ADD_BRANCH_T(double,   CBAvgTime)
        ADD_BRANCH_T(unsigned, nPhotonsCB)
        ADD_BRANCH_T(unsigned, nPhotonsTAPS)

        ADD_BRANCH_T(bool,     p_matched)

        ADD_BRANCH_T(unsigned,      Channel)

        ADD_BRANCH_T(std::vector<double>,       pi0chi2, 3)
        ADD_BRANCH_T(std::vector<double>,       pi0prob, 3)
        ADD_BRANCH_T(int,                       iBestPi0)
        ADD_BRANCH_T(std::vector<double>,       etachi2, 3)
        ADD_BRANCH_T(std::vector<double>,       etaprob, 3)
        ADD_BRANCH_T(int,                       iBestEta)

        ADD_BRANCH_T(int,                        bestHyp)

        ADD_BRANCH_T(std::vector<double>,       pi0_im,3)
        ADD_BRANCH_T(std::vector<double>,       eta_im,3)
        ADD_BRANCH_T(std::vector<double>,       eta_omega_im,3)
        ADD_BRANCH_T(std::vector<double>,       pi0_omega_im,3)


        ADD_BRANCH_T(double,   Pi0EtaFitChi2)
        ADD_BRANCH_T(double,   Pi0EtaFitProb)
        ADD_BRANCH_T(unsigned, Pi0EtaFitIterations)
        ADD_BRANCH_T(TLorentzVector,            lost_gamma_guess)
        ADD_BRANCH_T(TLorentzVector,            extra_gamma)
        ADD_BRANCH_T(std::vector<TLorentzVector>, bachelor_extra, 3)

    };

protected:
    void Analyse(const TEventData &data, const TEvent& event, manager_t&manager) override;


    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;


    //======== Tree ===============================================================

    TTree*  tree = nullptr;
    OmegaTree_t t;

    using combs_t = std::vector<std::vector<std::size_t>>;
    static const combs_t combs;




    //======== Settings ===========================================================

    const double cut_ESum;
    const double cut_Copl;
    const interval<double> photon_E_cb;
    const interval<double> photon_E_taps;
    const interval<double> proton_theta;


    ant::analysis::PromptRandom::Switch promptrandom;

    std::shared_ptr<utils::Fitter::UncertaintyModel> model;

    utils::KinFitter fitter;
    AccessibleFitter pi0eta_fitter;

    using doubles = std::vector<double>;

    struct MyTreeFitter_t {
        utils::TreeFitter treefitter;
        utils::TreeFitter::tree_t fitted_Omega;
        utils::TreeFitter::tree_t fitted_g_Omega;
        utils::TreeFitter::tree_t fitted_X;
        utils::TreeFitter::tree_t fitted_g1_X;
        utils::TreeFitter::tree_t fitted_g2_X;

        MyTreeFitter_t(const ParticleTypeTree& ttree, const ParticleTypeDatabase::Type& mesonT, const std::shared_ptr<const utils::Fitter::UncertaintyModel>& model);

        void HypTestCombis(const TParticleList& photons, doubles& chi2s, doubles& probs, doubles& ggims, doubles& gggims, int& bestIndex);
    };
    static size_t CombIndex(const ant::TParticleList& orig, const MyTreeFitter_t&);

    MyTreeFitter_t fitter_pi0;
    MyTreeFitter_t fitter_eta;

    bool AcceptedPhoton(const TParticlePtr& photon);
    bool AcceptedProton(const TParticlePtr& proton);

    TParticleList FilterPhotons(const TParticleList& list);
    TParticleList FilterProtons(const TParticleList& list);



    // Analysis options

    const bool   opt_save_after_kinfit = false;
    const double opt_kinfit_chi2cut    = 10.0;

    TagChMultiplicity tagChMult;

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
