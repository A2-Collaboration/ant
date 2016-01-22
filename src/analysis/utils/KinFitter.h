#pragma once

#include "analysis/data/Particle.h"

#include <string>
#include <vector>
#include <memory>
#include <array>

#include "TLorentzVector.h"

#include "APLCON.hpp"

#include <stdexcept>

class TTree;
class TH1D;
class TF1;

namespace ant {

class WrapTFile;

namespace analysis {
namespace utils {

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

    std::vector<kinVector> Photons;

    kinVector Proton = kinVector("Proton");

    PhotonBeamVector Beam = PhotonBeamVector("Beam");

    double EnergyResolution(const data::ParticlePtr& p) const;
    double ThetaResolution(const data::ParticlePtr& p) const;
    double PhiResolution(const data::ParticlePtr& p) const;

    static TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, const double m);

    std::unique_ptr<APLCON> aplcon;
    double result_chi2ndof       =  0.0;
    int result_iterations        =  0;
    int result_status            = -1;
    double result_probability    =  0.0;

    static double fct_TaggerEGausSigma(double E);

public:

    struct angular_sigma {
        using Hist = std::shared_ptr<TH1D>;
        Hist p0 = nullptr;
        Hist p1 = nullptr;
        Hist p2 = nullptr;

        double GetSigma(const unsigned element, const double E) const;

        static double f(const double x, const double p0, const double p1, const double p2) noexcept;
        static double f_root(const double* x, const double* p) noexcept;

        static TF1* GetTF1(const std::string& name="SigmaFit");

        void Load(ant::WrapTFile& f, const std::string& prefix, const int bins);
        Hist LoadHist(ant::WrapTFile& f, const std::string& name, const int bins);

        angular_sigma();
        ~angular_sigma();


    };

    angular_sigma cb_sigma_theta;
    angular_sigma cb_sigma_phi;
    angular_sigma taps_sigma_theta;
    angular_sigma taps_sigma_phi;

    KinFitter(const std::string& name, unsigned numGammas);

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const data::ParticlePtr& proton);
    void SetPhotons(const data::ParticleList& photons);

    data::ParticlePtr GetFittedProton() const;
    data::ParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;

    APLCON::Result_t DoFit();

    void SetupBranches(TTree* tree, std::string branch_prefix="");

    void LoadSigmaData(const std::string& filename);

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
};

}}} // namespace ant::analysis::utils
