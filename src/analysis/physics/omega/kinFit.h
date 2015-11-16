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
        const std::string Name;
        const unsigned nPhotons = 6;

        double Ek;
        double Theta;
        double Phi;

        double sEk;
        double sTheta;
        double sPhi;

        double energySmear(const double& E) const;
        double taggerSmear(const double& E) const;

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

    double taggerSmear(const double& E) const;
    static TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, const double m);

    std::unique_ptr<APLCON> aplcon;

    static constexpr auto nGamma = 3;

    std::pair<double,double> EgammaBeam;

    kinVector Proton = kinVector("Proton");

    std::vector<kinVector> Photons = {
        kinVector("Photon1"),
         kinVector("Photon2"),
         kinVector("Photon3")};

public:
    const std::string egammaName = "EBEAM";

    KinFitter();

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const ant::analysis::data::ParticlePtr& proton);
    void SetPhotons(const std::vector<ant::analysis::data::ParticlePtr>& photons);

    APLCON::Result_t DoFit();
};
