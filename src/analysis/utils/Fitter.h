#pragma once

#include "base/vec/LorentzVec.h"
#include "base/ParticleTypeTree.h"
#include "tree/TParticle.h"
#include "analysis/utils/combinatorics.h"
#include "base/std_ext/math.h"

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
namespace analysis {
namespace utils {

class Fitter {

public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    /**
     * @brief Uncertainties for E, theta, and phi
     */
    struct Uncertainties_t {
        double sigmaE     = {};
        double sigmaTheta = {};
        double sigmaPhi   = {};

        Uncertainties_t() = default;
        Uncertainties_t(const double E, const double Theta, const double Phi) : sigmaE(E), sigmaTheta(Theta), sigmaPhi(Phi) {}

        bool operator==(const Uncertainties_t& other) const noexcept {
            return sigmaE == other.sigmaE && sigmaTheta == other.sigmaTheta && sigmaPhi == other.sigmaPhi;
        }
    };

    /**
     * @brief Virtual base class for different Uncertainty Models for kion fitter.
     *        Derive and implement the GetSigmas() method.
     * @see UncertaintyModels::Constant
     * @see UncertaintyModels::MCExtracted
     */
    class UncertaintyModel {
    public:
        virtual ~UncertaintyModel();
        virtual Uncertainties_t GetSigmas(const TParticle& particle) const =0;
    };
    using UncertaintyModelPtr = std::shared_ptr<const Fitter::UncertaintyModel>;

    static const APLCON::Fit_Settings_t DefaultSettings;

    Fitter(const Fitter&) = delete;
    Fitter& operator=(const Fitter&) = delete;
    virtual ~Fitter();

    struct FitParticle
    {
        TParticlePtr Particle; // pointer to unfitted value

        struct Var_t {
            double Value = 0;
            double Sigma = 0;
            double Pull = 0;
        protected:
            friend struct FitParticle;
            void SetupBranches(TTree* tree, const std::string& prefix);
        };

        Var_t Ek;
        Var_t Theta;
        Var_t Phi;

        FitParticle(const std::string& name): Name(name) {}


    protected:
        friend class Fitter;
        friend class KinFitter;
        friend class TreeFitter;


        const std::string Name;
        void SetupBranches(TTree* tree, const std::string& prefix);
        static ant::LorentzVec GetVector(const std::vector<double>& EkThetaPhi, double m);

        TParticlePtr AsFitted(const ParticleTypeDatabase::Type& type);

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
    };

protected:

    Fitter(const std::string& fittername,
           const APLCON::Fit_Settings_t& settings, std::shared_ptr<const Fitter::UncertaintyModel>& uncertainty_model);

    Fitter(Fitter&&) = default;
    Fitter& operator=(Fitter&&) = default;

    std::shared_ptr<const Fitter::UncertaintyModel> uncertainty;
    std::unique_ptr<APLCON> aplcon;

    void LinkVariable(FitParticle& particle);

    void SetPhotonEkThetaPhi(FitParticle& photon, const TParticlePtr& p) const;

    static double fct_TaggerEGausSigma(double E);

private:
    static APLCON::Fit_Settings_t MakeDefaultSettings();

};


class KinFitter : public Fitter
{
public:

    KinFitter(const std::string& name,
              unsigned numGammas,
              UncertaintyModelPtr Uncertainty_model,
              const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    virtual ~KinFitter();

    KinFitter(const KinFitter&) = delete;
    KinFitter& operator=(const KinFitter&) = delete;
    KinFitter(KinFitter&&) = default;
    KinFitter& operator=(KinFitter&&) = default;

    void SetEgammaBeam(double ebeam);
    void SetProton(const TParticlePtr& proton);
    virtual void SetPhotons(const TParticleList& photons);

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    double GetFittedBeamEPull() const;

    std::vector<FitParticle> GetFitParticles() const;

    void SetupBranches(TTree* tree, std::string branch_prefix="");

    APLCON::Result_t DoFit();

protected:

    struct PhotonBeamVector {
        double E_before = std_ext::NaN;
        double E     = std_ext::NaN;
        double Sigma = std_ext::NaN;
        double Pull  = std_ext::NaN;

        const std::string Name;

        PhotonBeamVector(const std::string& name);

        std::vector<double*> Adresses()
        {
            return { std::addressof(E) };
        }

        std::vector<double*> Adresses_Sigma()
        {
            return { std::addressof(Sigma) };
        }

        std::vector<double*> Adresses_Pulls()
        {
            return { std::addressof(Pull) };
        }

    };

    // it's pretty important that those things are pointers,
    // since the members are linked to APLCON in ctor!
    // A move/copy of those members may not happen, so we just
    // point to their fixed location in memory.
    std::vector<std::shared_ptr<FitParticle>> Photons;
    std::unique_ptr<FitParticle> Proton;
    std::unique_ptr<PhotonBeamVector> Beam;

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
        bool Excluded;
        nodesetup_t(double IM_sigma = 1.0, bool excluded = false) :
            IM_Sigma(IM_sigma),
            Excluded(excluded)
        {}

        // the getter lets the user decide how to setup the node
        struct getter : std::function<nodesetup_t(ParticleTypeTree)> {
            // use base class constructors
            using std::function<nodesetup_t(ParticleTypeTree)>::function;
            // provide default which gets the default nodesetup_t
            getter() :
                std::function<nodesetup_t(ParticleTypeTree)>([] (ParticleTypeTree) {return nodesetup_t{};})
            {}
        };
    };

    // construct TreeFitter always with KinFit
    TreeFitter(const std::string& name,
               ParticleTypeTree ptree,
               UncertaintyModelPtr uncertainty_model,
               nodesetup_t::getter nodeSetup = {},
               const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    TreeFitter(const TreeFitter&) = delete;
    TreeFitter& operator=(const TreeFitter&) = delete;
    TreeFitter(TreeFitter&&) = default;
    TreeFitter& operator=(TreeFitter&&) = default;

    virtual ~TreeFitter();

    struct node_t {
        node_t(const ParticleTypeTree& ptree) : TypeTree(ptree) {}
        const ParticleTypeTree TypeTree;
        LorentzVec LVSum;
        std::shared_ptr<FitParticle> Leave;
        int PhotonLeaveIndex = -1; // according to list given by SetPhotons
        bool operator<(const node_t& rhs) const {
            return TypeTree->Get() < rhs.TypeTree->Get();
        }
    };

    using tree_t = Tree<node_t>::node_t;

    virtual void SetPhotons(const TParticleList& photons) override;

    tree_t GetTreeNode(const ParticleTypeDatabase::Type& type) const {
        tree_t treenode = nullptr;
        tree->Map_nodes([&treenode, &type] (tree_t t) {
            if(t->Get().TypeTree->Get() == type)
                treenode = t;
        });
        return treenode;
    }

    bool NextFit(APLCON::Result_t& fit_result);

    using current_comb_t = KofNvector<TParticlePtr>;
    const current_comb_t& GetCurrentCombination() const {
        return *current_comb_ptr;
    }



protected:

    // use while(NextFit()) {} instead
    // to run all fits
    using KinFitter::DoFit;

    static tree_t MakeTree(ParticleTypeTree ptree);
    static unsigned CountGammas(ParticleTypeTree ptree);


    const tree_t tree;
    std::vector<tree_t> tree_leaves;
    using permutations_t = std::vector<std::vector<size_t>>;
    permutations_t permutations;
    permutations_t::const_iterator current_perm;

    // use unique_ptr since KofNvector does not have default ctor
    std::unique_ptr<current_comb_t> current_comb_ptr;

};


}}} // namespace ant::analysis::utils


