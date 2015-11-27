#pragma once

#include "data/Particle.h"

#include <string>
#include <vector>
#include <memory>
#include <array>

#include "TLorentzVector.h"

#include "APLCON.hpp"

class APLCON;
class TTree;

class KinFitter
{

private:
    struct kinVector
    {
        double Ek            = 0.0;
        double Theta         = 0.0;
        double Phi           = 0.0;

        double sigmaEk       = 0.0;
        double sigmaTheta    = 0.0;
        double sigmaPhi      = 0.0;

        double pullEk        = 0.0;
        double pullTheta     = 0.0;
        double pullPhi       = 0.0;

        const std::string Name;

        double EnergyResolution(const double& E) const;

        std::vector<double*> Addresses()
        {
            return { std::addressof(Ek),
                     std::addressof(Theta),
                     std::addressof(Phi)};
        }
        std::vector<double*> Addresses_Sigma()
        {
            return { std::addressof(sigmaEk),
                     std::addressof(sigmaTheta),
                     std::addressof(sigmaPhi)};
        }

        std::vector<double*> Addresses_Pulls()
        {
            return { std::addressof(pullEk),
                     std::addressof(pullTheta),
                     std::addressof(pullPhi)};
        }

        kinVector(const std::string& name): Name(name) {}

        void SetupBranches(TTree* tree, const std::string& prefix);
    };

    struct PhotonBeamVector {
        double E     = 0.0;
        double Sigma = 0.0;

        const std::string name;

        PhotonBeamVector(const std::string& Name);

        std::vector<double*> Adresses()
        {
            return { std::addressof(E) };
        }

        std::vector<double*> Adresses_Sigma()
        {
            return { std::addressof(Sigma)};
        }

    };

    static constexpr auto nGamma = 3;
    std::vector<kinVector> Photons = {
         kinVector("Photon1"),
         kinVector("Photon2"),
         kinVector("Photon3")};

    kinVector Proton = kinVector("Proton");

    PhotonBeamVector Beam = PhotonBeamVector("Beam");

    double EnergyResolution(const ant::analysis::data::ParticlePtr& p) const;
    double ThetaResolution(const ant::analysis::data::ParticlePtr& p) const;
    double PhiResolution(const ant::analysis::data::ParticlePtr& p) const;

    static TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, const double m);

    std::unique_ptr<APLCON> aplcon;
    double result_chi2ndof       =  0.0;
    int result_iterations        =  0;
    int result_status            = -1;
    double result_probability    =  0.0;


    static double fct_GlobalEResolution(double E);
    static double fct_TaggerEGausSigma(double E);





public:

    KinFitter();

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const ant::analysis::data::ParticlePtr& proton);
    void SetPhotons(const std::vector<ant::analysis::data::ParticlePtr>& photons);

    APLCON::Result_t DoFit();

    void SetupBranches(TTree* tree);
};
