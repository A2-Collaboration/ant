#pragma once

#include "analysis/utils/KinFitter.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class Etap3pi0 : public Physics {

protected:

    // =======================   constants =====================================================

    const double IM_mean_etaprime = 895.0;
    const double IM_sigma_etaprime = 26.3;

    const double IM_mean_eta = 515.5;
    const double IM_sigma_eta = 19.4;

    const double IM_mean_pi0 = 126.0;
    const double IM_sigma_pi0 = 15.0;


    const double copl_opening_sigma = 7.80;

    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> signal_tree;
    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> reference_tree;
    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> bkg_tree;

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

    ant::analysis::PromptRandom::Switch promptrandom;
    utils::KinFitter fitter;
    TTree* tree;

    struct branches {
        TLorentzVector proton= {};

        TLorentzVector fittedProton= {};

        TLorentzVector trueProton= {};

        TLorentzVector MM= {};

        double coplanarity= {};

        double      taggWeight= {};
        double      taggE= {};
        unsigned    taggCh= {};
        double      taggTime= {};

        std::vector<TLorentzVector> pi0 = std::vector<TLorentzVector>(3);
        double pi0_chi2[3]= {};
        double pi0_prob[3]= {};
        int    pi0_iteration[3]= {};
        int    pi0_status[3]= {};

        TLorentzVector etaprime= {};
        double event_chi2= {};
        double event_prob= {};
        int    event_iteration= {};
        int    event_status= {};

        int type= {};

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







    TLorentzVector MakeLoretzSum(const TParticleList& particles);
public:
    Etap3pi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
private:
    bool MakeMCProton(const TEventData& mcdata, TParticlePtr& proton);
};


}}}
