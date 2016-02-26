#pragma once

#include "analysis/utils/Fitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class Etap3pi0 : public Physics {

protected:
    //geometry
    ant::analysis::utils::A2SimpleGeometry geometry;

    // =======================   constants =====================================================

    const double IM_mean_etaprime = 895.0;
    const double IM_sigma_etaprime = 26.3;

    const double IM_mean_eta = 515.5;
    const double IM_sigma_eta = 19.4;

    const double IM_mean_pi0 = 126.0;
    const double IM_sigma_pi0 = 15.0;


    const double copl_opening_sigma = 7.80;

    ParticleTypeTree signal_tree;
    ParticleTypeTree reference_tree;
    ParticleTypeTree bkg_tree;

    struct settings_t
    {
        std::map<int,std::string> EventTypes= {{0,"signal"},
                                               {1,"reference"},
                                               {2,"background"},
                                               {-1,"other"}};
        const double EMBChi2Cut= std::numeric_limits<double>::infinity();
        const double fourConstrainChi2Cut= 40;
        const double EsumCB= 550;
    };
    settings_t phSettings;

    const std::vector<std::vector<std::pair<unsigned,unsigned>>> combinations =
    {
        { {0, 1}, {2, 3}, {4, 5} },
        { {0, 1}, {2, 4}, {3, 5} },
        { {0, 1}, {2, 5}, {3, 4} },

        { {0, 2}, {1, 3}, {4, 5} },
        { {0, 2}, {1, 4}, {3, 5} },
        { {0, 2}, {1, 5}, {3, 4} },

        { {0, 3}, {1, 2}, {4, 5} },
        { {0, 3}, {1, 4}, {2, 5} },
        { {0, 3}, {1, 5}, {2, 4} },

        { {0, 4}, {1, 2}, {3, 5} },
        { {0, 4}, {1, 3}, {2, 5} },
        { {0, 4}, {1, 5}, {2, 3} },

        { {0, 5}, {1, 2}, {3, 4} },
        { {0, 5}, {1, 3}, {2, 4} },
        { {0, 5}, {1, 4}, {2, 3} }
    };

    utils::TreeFitter fitterSig;
    utils::TreeFitter fitterRef;


    ant::analysis::PromptRandom::Switch promptrandom;
    utils::KinFitter kinFitterEMB;

    TTree* tree;

    struct branches {
        TLorentzVector etaprimeCand= {};

        TLorentzVector proton= {};

        TLorentzVector fittedProton= {};

        TLorentzVector trueProton= {};

        TLorentzVector MM= {};

        double coplanarity= {};
        double EsumCB= {};

        double      taggWeight= {};
        double      taggE= {};
        unsigned    taggCh= {};
        double      taggTime= {};

        double EMB_chi2= std::numeric_limits<double>::infinity();

        std::vector<TLorentzVector> pi0 = std::vector<TLorentzVector>(3);
        double pi0_chi2[3]= {  std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()  };
        double pi0_prob[3]= {};
        int    pi0_iteration[3]= {};
        int    pi0_status[3]= {};

        double event_chi2_ref= std::numeric_limits<double>::infinity();
        double event_chi2_sig= std::numeric_limits<double>::infinity();
        double event_prob= {};
        int    event_iteration= {};
        int    event_status= {};

        int type= -1;
        int truetype= -1;

        std::string decayString= {};

        void SetBranches(TTree* tree);

    };

    branches vars;



    // =======================   datastorage  ==================================================


    //histograms
    std::map<std::string,std::map<std::string,TH1*>> hists;
    void AddHist1D(const std::string& category, const std::string& hname,
                   const std::string& title,
                   const std::string& xlabel, const std::string& ylabel,
                   const BinSettings& bins);
    void AddHist2D(const std::string& category, const std::string& hname,
                   const std::string& title,
                   const std::string& xlabel, const std::string& ylabel,
                   const BinSettings& xbins, const BinSettings& ybins);


    //helpers:
    TLorentzVector MakeLoretzSum(const TParticleList& particles);

    void MakeSignal(const TParticleList& photonLeaves);
    void MakeReference(const TParticleList& photonLeaves);
    bool MakeMCProton(const TEventData& mcdata, TParticlePtr& proton);


    double getEnergyMomentumConservation(double EBeam, const TParticleList& photons, const TParticlePtr& proton);
public:
    Etap3pi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};


}}}
