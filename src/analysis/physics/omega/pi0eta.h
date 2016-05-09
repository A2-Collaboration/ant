#pragma once
#include "analysis/physics/omega/omega.h"
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



class Pi0Eta : public OmegaBase {
public:
    struct Pi0EtaTree_t : WrapTTree {
        Pi0EtaTree_t();

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons, 4)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_fitted, 4)

        ADD_BRANCH_T(TLorentzVector,              p)
        ADD_BRANCH_T(double,                      p_Time)
        ADD_BRANCH_T(int,                         p_detector)

        ADD_BRANCH_T(TLorentzVector,              p_true)
        ADD_BRANCH_T(TLorentzVector,              p_fitted)

        ADD_BRANCH_T(TLorentzVector,              gggg)
        ADD_BRANCH_T(TLorentzVector,              gggg_fitted)
        ADD_BRANCH_T(TLorentzVector,              mm)
        ADD_BRANCH_T(double,                      copl_angle)
        ADD_BRANCH_T(double,                      p_mm_angle)

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

        ADD_BRANCH_T(unsigned,      Channel)


        ADD_BRANCH_T(double,   TreeFitChi2)
        ADD_BRANCH_T(double,   TreeFitProb)
        ADD_BRANCH_T(TLorentzVector, TreeFitpi0)
        ADD_BRANCH_T(TLorentzVector, TreeFiteta)
        ADD_BRANCH_T(TLorentzVector, TreeFitgggg)

    };

protected:
    void Analyse(const TEventData& data, const TEvent& event, manager_t&manager) override;


    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;


    //======== Tree ===============================================================

    TTree*  tree = nullptr;
    Pi0EtaTree_t t;

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

    using doubles = std::vector<double>;

    struct MyTreeFitter_t {
        utils::TreeFitter treefitter;

        utils::TreeFitter::tree_t fitted_pi0;
        utils::TreeFitter::tree_t fitted_pi0_g1;
        utils::TreeFitter::tree_t fitted_pi0_g2;

        utils::TreeFitter::tree_t fitted_eta;
        utils::TreeFitter::tree_t fitted_eta_g1;
        utils::TreeFitter::tree_t fitted_eta_g2;

        MyTreeFitter_t(const ParticleTypeTree& ttree, const std::shared_ptr<const utils::Fitter::UncertaintyModel>& model);

        void HypTestCombis(const TParticleList& unfitted, const TParticleList& kinfitted, double& chi2,
                           double& prob,
                           TLorentzVector& pi0,
                           TLorentzVector& eta,
                           TLorentzVector& gggg);
    };

    MyTreeFitter_t treefitter;

    bool AcceptedPhoton(const TParticlePtr& photon);
    bool AcceptedProton(const TParticlePtr& proton);

    TParticleList FilterPhotons(const TParticleList& list);
    TParticleList FilterProtons(const TParticleList& list);



    // Analysis options

    const bool   opt_save_after_treefit = false;
    const double opt_treefit_chi2cut    = 10.0;


public:

    Pi0Eta(const std::string& name, OptionsPtr opts);
    virtual ~Pi0Eta();

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
