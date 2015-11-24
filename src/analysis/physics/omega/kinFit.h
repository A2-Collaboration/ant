#pragma once

#include "data/Particle.h"

#include <string>
#include <vector>
#include <memory>
#include <array>

#include "TLorentzVector.h"

#include "APLCON.hpp"

class APLCON;

class KinFitter
{

private:
    struct kinVector
    {
        double Ek;
        double Theta;
        double Phi;

        double sEk;
        double sTheta;
        double sPhi;

        const std::string Name;

        double EnergyResolution(const double& E) const;

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

        kinVector(const std::string& name): Name(name) {}
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


    static double fct_GlobalEResolution(double E);
    static double fct_TaggerEGausSigma(double E);





public:

    KinFitter();

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const ant::analysis::data::ParticlePtr& proton);
    void SetPhotons(const std::vector<ant::analysis::data::ParticlePtr>& photons);

    APLCON::Result_t DoFit();
};
