#pragma once

#include "analysis/physics/Physics.h"

#include <APLCON.hpp>

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

    // =======================   aplcon    =====================================================


    class KinFitter
    {

    private:
        struct kinVector
        {
            const std::string Name;
            const unsigned nPhotons = 6;

            double Ek;
            double Theta;
            double Phi;

            double sEk;
            double sTheta;
            double sPhi;

            double energySmear(const double& E) const;

            std::vector<double*> Adresses()
            {
                return { std::addressof(Ek),
                         std::addressof(Theta),
                         std::addressof(Phi)};
            }
            std::vector<double*> Adresses_Sigma()
            {
                return { std::addressof(sEk),
                         std::addressof(sTheta),
                         std::addressof(sPhi)};
            }

            void SetEkThetaPhi(double ek, double theta, double phi);

            kinVector(const std::string& name): Name(name) {}
        };

        double taggerSmear( const double& E) const;
        TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, const double m) const;

        APLCON aplcon;

        std::pair<double,double> EgammaBeam;
        kinVector ProtonTAPS;
        std::vector<kinVector> Photons;

        const double IM_Mother;

    public:
        const std::string egammaName = "EBEAM";




        KinFitter(const ParticleTypeDatabase::Type& motherParticle = ParticleTypeDatabase::EtaPrime);

        void SetEgammaBeam(const double& ebeam);
        void SetProtonTAPS(const data::ParticlePtr& proton);
        void SetPhotons(const std::vector<data::ParticlePtr>& photons);

        APLCON::Result_t DoFit() { return aplcon.DoFit(); }
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

    KinFitter fitToEtaPrime;
    APLCON::Result_t result_fitToEtaPrime;

    enum class filterType {
        Chi2,
        KinFit,
        ProtonInTaps,
        ProtonTOF
    };

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




    void FillCrossChecks(const data::ParticleList& photons, const data::ParticleList& mcphotons);

    bool MakeMCProton(const data::Event::Data& mcdata, data::ParticlePtr& proton);

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
