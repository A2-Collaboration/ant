#pragma once

#include "tree/TParticle.h"
#include "base/ParticleTypeTree.h"

#include "TLorentzVector.h"

#include "APLCON.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <functional>

class TTree;
class TH1D;
class TF1;

namespace ant {

class WrapTFile;

namespace analysis {
namespace utils {

class Fitter {

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

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    void LoadSigmaData(const std::string& filename);

    APLCON::Result_t DoFit();

protected:

    Fitter(const std::string& fittername);

    std::unique_ptr<APLCON> aplcon;

    void SetupBranches(TTree* tree, std::string branch_prefix="");

    struct FitParticle
    {
        struct Var_t {
            double Value = 0;
            double Sigma = 0;
            double Pull = 0;
            void SetupBranches(TTree* tree, const std::string& prefix);

        };

        Var_t Ek;
        Var_t Theta;
        Var_t Phi;

        const std::string Name;

        std::vector<double*> Addresses()
        {
            return { std::addressof(Ek.Value),
                     std::addressof(Theta.Value),
                     std::addressof(Phi.Value)};
        }
        std::vector<double*> Addresses_Sigma()
        {
            return { std::addressof(Ek.Sigma),
                     std::addressof(Theta.Sigma),
                     std::addressof(Phi.Sigma)};
        }

        std::vector<double*> Addresses_Pulls()
        {
            return { std::addressof(Ek.Pull),
                     std::addressof(Theta.Pull),
                     std::addressof(Phi.Pull)};
        }

        FitParticle(const std::string& name): Name(name) {}

        void SetupBranches(TTree* tree, const std::string& prefix);
    };

    static TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, const double m);

    angular_sigma cb_sigma_theta;
    angular_sigma cb_sigma_phi;
    angular_sigma taps_sigma_theta;
    angular_sigma taps_sigma_phi;

    double EnergyResolution(const TParticlePtr& p) const;
    double ThetaResolution(const TParticlePtr& p) const;
    double PhiResolution(const TParticlePtr& p) const;

    static double fct_TaggerEGausSigma(double E);

private:
    double result_chi2ndof       =  0.0;
    int result_iterations        =  0;
    int result_status            = -1;
    double result_probability    =  0.0;

};

class KinFitter : public Fitter
{
public:

    KinFitter(const std::string& name, unsigned numGammas);

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const TParticlePtr& proton);
    void SetPhotons(const TParticleList& photons);

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;

    void SetupBranches(TTree* tree, std::string branch_prefix="");

protected:

    struct PhotonBeamVector {
        double E     = 0.0;
        double Sigma = 0.0;

        const std::string Name;

        PhotonBeamVector(const std::string& name);

        std::vector<double*> Adresses()
        {
            return { std::addressof(E) };
        }

        std::vector<double*> Adresses_Sigma()
        {
            return { std::addressof(Sigma)};
        }

    };

    std::vector<FitParticle> Photons;

    FitParticle Proton = FitParticle("Proton");

    PhotonBeamVector Beam = PhotonBeamVector("Beam");

};

class TreeFitter : public Fitter
{

public:
    TreeFitter(const std::string& name, ParticleTypeTree ptree);

protected:

    struct Node_t {
        Node_t(ParticleTypeTree ptree) : TypeTree(ptree) {}
        ParticleTypeTree TypeTree;
        TParticlePtr Particle;
        TParticlePtr FittedParticle;
    private:
        friend class TreeFitter;
        std::uint64_t leave_sum = 0;
    };

    using tree_t = std::shared_ptr<Tree<Node_t>>;

    void WalkTree(ParticleTypeTree ptree, tree_t tree);

    static void CalcSum(tree_t tree);
    static void CopyTree(tree_t source, tree_t dest);
    static bool IsEqual(tree_t a, tree_t b);

    tree_t tree;

    std::vector<tree_t> tree_leaves;
};



}}} // namespace ant::analysis::utils
