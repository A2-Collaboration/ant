#ifndef APLCON_H
#define APLCON_H

#include "AntPhysics.h"

#include <APLCON.hpp>
#include "plot/Histogram.h"


#include <vector>
#include <map>
#include <random>

#include "TH1D.h"
#include "TH2D.h"

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> v = v1;
    v.insert(v.end(),v2.begin(),v2.end());
    return v;
}

namespace ant {
namespace analysis {

class TestAPLCON: public Physics {

protected:


    // choose here what you want to do
    // please also provide GoAT trees with matching MC true information...
    static constexpr bool includeIMconstraint = false;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nPhotons = 6;
//    const double IM = ParticleTypeDatabase::Pi0.Mass();
    const double IM = ParticleTypeDatabase::EtaPrime.Mass();

    const map<short, string> component = {{0, "Energy"}, {1, "Theta"}, {2, "Phi"}};


    TH2D* banana;
    TH1D* particles;
    TH1D* tagger;
    TH1D* ntagged;
    TH1D* cbesum;

    std::map<const ParticleTypeDatabase::Type*, TH1D*> numParticleType;

    TH1D* chisquare;
    TH1D* probability;
    TH1D* iterations;
    std::map<std::string, TH1D*> pulls;


    TH1D* im_true;
    TH1D* im_smeared;
    TH1D* im_fit;

    TH1D* vertex_z_after;
    TH1D* vertex_z_before;



    // lightweight structure for linking to fitter
    struct FitParticle {
        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
        }

        static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
                                           const Double_t m);
        static TLorentzVector Make(const FitParticle& p,
                                   const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
        }

        std::vector<double*> Link() {
            return {std::addressof(Ek),
                    std::addressof(Theta),
                    std::addressof(Phi)};
        }
        std::vector<double*> LinkSigma() {
            return {std::addressof(Ek_Sigma),
                    std::addressof(Theta_Sigma),
                    std::addressof(Phi_Sigma)};
        }

        void Smear();

        double Ek;
        double Ek_Sigma;
        double Theta;
        double Theta_Sigma;
        double Phi;
        double Phi_Sigma;
    private:
        static std::default_random_engine generator;
    };

    void FillIM(TH1D* h, const std::vector<FitParticle>& photons);

    APLCON fitter;
    FitParticle beam;
    std::vector<FitParticle> photons;
    FitParticle proton;


public:
    TestAPLCON(const mev_t energy_scale=1000.0);
    virtual ~TestAPLCON() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
