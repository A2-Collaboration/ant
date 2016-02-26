#pragma once

#include "tree/TParticle.h"
#include "base/ParticleTypeTree.h"

#include "analysis/utils/combinatorics.h"

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

protected:

    Fitter(const std::string& fittername);

    std::unique_ptr<APLCON> aplcon;

    struct FitParticle
    {
        TParticlePtr Particle;

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

        static TLorentzVector GetVector(const std::vector<double>& EkThetaPhi, double m);
    };

    void LinkVariable(FitParticle& particle);

    angular_sigma cb_sigma_theta;
    angular_sigma cb_sigma_phi;
    angular_sigma taps_sigma_theta;
    angular_sigma taps_sigma_phi;

    double EnergyResolution(const TParticlePtr& p) const;
    double ThetaResolution(const TParticlePtr& p) const;
    double PhiResolution(const TParticlePtr& p) const;

    static double fct_TaggerEGausSigma(double E);

};

class KinFitter : public Fitter
{
public:

    KinFitter(const std::string& name, unsigned numGammas);

    void SetEgammaBeam(const double& ebeam);
    void SetProton(const TParticlePtr& proton);
    virtual void SetPhotons(const TParticleList& photons);

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;

    void SetupBranches(TTree* tree, std::string branch_prefix="");

    APLCON::Result_t DoFit();

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

    std::vector<std::shared_ptr<FitParticle>> Photons;

    FitParticle Proton = FitParticle("Proton");

    PhotonBeamVector Beam = PhotonBeamVector("Beam");
private:
    double result_chi2ndof       =  0.0;
    int result_iterations        =  0;
    int result_status            = -1;
    double result_probability    =  0.0;
};

class TreeFitter : public KinFitter
{

public:
    struct nodesetup_t {
        double IM_Sigma;
        double Excluded;
        explicit nodesetup_t(double IM_sigma = 1.0, bool excluded = false) :
            IM_Sigma(IM_sigma),
            Excluded(excluded)
        {}
    };

    // construct TreeFitter with KinFit
    TreeFitter(const std::string& name,
               ParticleTypeTree ptree,
               unsigned kinFitGammas,
               std::function<nodesetup_t(ParticleTypeTree)> nodeSetup = [] (ParticleTypeTree) {return nodesetup_t{};}
              );

    // construct TreeFitter without additional KinFit
    TreeFitter(const std::string& name,
               ParticleTypeTree ptree,
               std::function<nodesetup_t(ParticleTypeTree)> nodeSetup = [] (ParticleTypeTree) {return nodesetup_t{};}
              ) : TreeFitter(name, ptree, 0, nodeSetup)
    {}


    struct node_t {
        node_t(const ParticleTypeTree& ptree) : TypeTree(ptree) {}
        const ParticleTypeTree TypeTree;
        TLorentzVector LVSum;
        std::shared_ptr<FitParticle> Leave;
        bool operator<(const node_t& rhs) const {
            return TypeTree->Get() < rhs.TypeTree->Get();
        }
    };

    using tree_t = Tree<node_t>::node_t;

    void SetLeaves(const TParticleList& photons);

    virtual void SetPhotons(const TParticleList& photons) override {
        SetLeaves(photons);
    }

    tree_t GetTreeNode(const ParticleTypeDatabase::Type& type) const {
        tree_t treenode = nullptr;
        tree->Map_nodes([&treenode, &type] (tree_t t) {
            if(t->Get().TypeTree->Get() == type)
                treenode = t;
        });
        return treenode;
    }

    bool NextFit(APLCON::Result_t& fit_result);

    std::vector<size_t> GetCurrentPermutation() const {
        return *current_perm;
    }

protected:

    static tree_t MakeTree(ParticleTypeTree ptree);

    const tree_t tree;
    using permutations_t = std::vector<std::vector<size_t>>;
    permutations_t permutations;
    permutations_t::const_iterator current_perm;

    // use unique_ptr since KofNvector does not have default ctor
    using current_comb_t = KofNvector<TParticlePtr>;
    std::unique_ptr<current_comb_t> current_comb_ptr;

};



}}} // namespace ant::analysis::utils
