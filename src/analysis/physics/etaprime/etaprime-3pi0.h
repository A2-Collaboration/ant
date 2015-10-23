#pragma once

#include "analysis/physics/Physics.h"
class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class Etap3pi0 : public Physics {

protected:

   // =======================   constants =====================================================

    const double IM_mean_etaprime = 906.0;
    const double IM_sigma_etaprime = 26.3;

    const double IM_mean_eta = 515.5;
    const double IM_sigma_eta = 19.4;

    const double IM_mean_pi0 = 126.0;
    const double IM_sigma_pi0 = 15.0;

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

   // =======================   structs   =====================================================

    using MesonCandidate = std::pair<data::ParticlePtr,double>;    // <particle,chi2>

    struct result_t {
        double Chi2_intermediate = std::numeric_limits<double>::infinity();
        double Chi2_mother       = std::numeric_limits<double>::infinity();
        double chi2() const { return Chi2_mother + Chi2_intermediate; }

        bool success = false;

        std::vector<data::ParticlePtr> g_final;
        std::vector<MesonCandidate> mesons;

        TLorentzVector mother;

        result_t() : g_final(6), mesons(3), mother(0,0,0,0){}
    };

    struct DalitzVars
    {
        double TMean;

        double s1;
        double s2;
        double s3;

        double x;
        double y;

        double z;

        DalitzVars(result_t r);
    };

   // =======================   datastorage  ==================================================

    std::string dataset;

    // xchecks
    TH1D* hNgamma;
    TH1D* hNgammaMC;
    TH1D* h2g;
    TH1D* h6g;
    TH1D* h6photonEvents;

    TH1D* ch_3pi0_IM_etap;
    TH1D* ch_3pi0_IM_pi0;

    TH1D* ch_eta2pi0_IM_etap;
    TH1D* ch_eta2pi0_IM_pions;
    TH1D* ch_eta2pi0_IM_etas;

    TH1D* mcdalitz_z;
    TH2D* mcdalitz_xy;
    TH1D* dalitz_z;
    TH2D* dalitz_xy;

    void FillCrossChecks(const data::ParticleList& photons, const data::ParticleList& mcphotons);

    Etap3pi0::result_t Make3pi0(const data::ParticleList& photons);
    Etap3pi0::result_t MakeEta2pi0(const data::ParticleList& photons);
    Etap3pi0::result_t MakeMC3pi0(const data::Event::Data &mcEvt);

    void FillIm(const Etap3pi0::result_t& result, const ParticleTypeDatabase::Type& type, TH1D* hist);
    void FillImEtaPrime(const Etap3pi0::result_t& result, TH1D* hist);
public:
    Etap3pi0(const std::string& name, PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void ShowResult() override;
};


}}}
