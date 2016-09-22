#pragma once

#include "base/vec/LorentzVec.h"
#include "base/ParticleTypeTree.h"
#include "tree/TParticle.h"
#include "analysis/utils/combinatorics.h"
#include "analysis/utils/Uncertainties.h"
#include "base/std_ext/math.h"

#include "APLCON.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <functional>
#include <type_traits>

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

    static const APLCON::Fit_Settings_t DefaultSettings;

    Fitter(const Fitter&) = delete;
    Fitter& operator=(const Fitter&) = delete;
    virtual ~Fitter();

    struct FitVariable {
        // initialize linked values at zero, important for Z Vertex
        double Value = 0;
        double Value_before = std_ext::NaN;
        double Sigma = 0;
        double Sigma_before = std_ext::NaN;
        double Pull = 0;
        void SetValueSigma(double v, double s) {
            Value = v;
            Value_before = v;
            Sigma = s;
            Sigma_before = s;
        }
    };

    struct FitParticle
    {
        TParticlePtr      Particle; // pointer to unfitted particle
        Detector_t::Any_t Detector; // remember detector type provided by uncertainty model

        std::vector<FitVariable> Vars;

        TParticlePtr AsFitted() const;

        using pulls_t = std::vector<double>; // used also in KinFitter
        std::vector<double> GetSigmas() const;
        pulls_t GetPulls() const;

        FitParticle(const std::string& name,
                    APLCON& aplcon,
                    std::shared_ptr<FitVariable> z_vertex);
        virtual ~FitParticle();

    protected:

        friend class Fitter;
        friend class KinFitter;
        friend class TreeFitter;

        FitVariable& GetEk();
        void Set(const TParticlePtr& p, const UncertaintyModel& uncertainty);

        const std::string Name;
        const std::shared_ptr<const FitVariable> Z_Vertex;

        ant::LorentzVec GetLorentzVec(const std::vector<double>& values,
                                      double z_vertex) const;
    };

protected:

    Fitter(const std::string& fittername,
           const APLCON::Fit_Settings_t& settings,
           UncertaintyModelPtr uncertainty_model);

    Fitter(Fitter&&) = default;
    Fitter& operator=(Fitter&&) = default;

    UncertaintyModelPtr uncertainty;
    std::unique_ptr<APLCON> aplcon;

private:
    static APLCON::Fit_Settings_t MakeDefaultSettings();

};


class KinFitter : public Fitter
{
public:

    KinFitter(const std::string& name,
              unsigned numGammas,
              UncertaintyModelPtr uncertainty_model,
              bool fit_Z_vertex = false,
              const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    virtual ~KinFitter();

    KinFitter(const KinFitter&) = delete;
    KinFitter& operator=(const KinFitter&) = delete;
    KinFitter(KinFitter&&) = default;
    KinFitter& operator=(KinFitter&&) = default;

    void SetEgammaBeam(double ebeam);
    void SetZVertexSigma(double sigma);
    void SetProton(const TParticlePtr& proton);
    virtual void SetPhotons(const TParticleList& photons);

    bool IsZVertexFitEnabled() const noexcept;

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    TParticlePtr GetFittedBeamParticle() const;
    double GetFittedZVertex() const;

    double GetBeamEPull() const;
    double GetZVertexPull() const;

    FitParticle::pulls_t GetProtonPulls() const;
    /**
     * @brief GetPhotonsPulls
     * @return matrix with first index specifying parameter (0...3), second the photons.
     * in congruence with GetProtonPulls
     */
    std::vector<std::vector<double>> GetPhotonsPulls() const;

    /**
     * @brief GetFitParticles returns as first item the proton, then all n photons
     * @return FitParticle contain all info about the fitted state
     */
    std::vector<FitParticle> GetFitParticles() const;

    APLCON::Result_t DoFit();

protected:

    struct BeamE_t : FitVariable {
        const std::string Name = "Beam";
    };

    struct Z_Vertex_t : FitVariable {
        const std::string Name = "ZVertex";
    };

    // it's pretty important that those things are pointers,
    // since the members are linked to APLCON in ctor!
    // A move/copy of those members may not happen, so we just
    // point to their fixed location in memory.
    std::vector<std::shared_ptr<FitParticle>> Photons;
    std::shared_ptr<FitParticle> Proton;
    std::unique_ptr<BeamE_t>    BeamE;
    std::shared_ptr<Z_Vertex_t> Z_Vertex;

    static LorentzVec MakeBeamLorentzVec(double BeamE);

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
               bool fit_Z_vertex = false,
               nodesetup_t::getter nodeSetup = {},
               const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    TreeFitter(const TreeFitter&) = delete;
    TreeFitter& operator=(const TreeFitter&) = delete;
    TreeFitter(TreeFitter&&) = default;
    TreeFitter& operator=(TreeFitter&&) = default;

    virtual ~TreeFitter();

    virtual void SetPhotons(const TParticleList& photons) override;

    using iteration_filter_t = std::function<double()>;
    /**
     * @brief SetIterationFilter
     * @param filter function returning a quality factor for the iteration. factor=0 means skip iteration.
     * @param max only runs the best max number of iterations
     * @note the info about the current state should be managed by captured pointers to tree nodes of interest
     */
    void SetIterationFilter(iteration_filter_t filter, unsigned max = 0);

    /**
     * @brief The node_t struct represents
     */
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
    /**
     * @brief GetTreeNode returns the first pointer to a tree node
     * @param type
     * @return nullptr if not found
     * @see GetTreeNodes to get all matching nodes
     */
    tree_t GetTreeNode(const ParticleTypeDatabase::Type& type) const;
    std::vector<tree_t> GetTreeNodes(const ParticleTypeDatabase::Type& type) const;

    /**
     * @brief NextFit runs the next fit iteration
     * @param fit_result
     * @return true if fit successfully run, false if no fit executed and thus fit_result unchanged
     */
    bool NextFit(APLCON::Result_t& fit_result);

protected:

    // use while(NextFit()) {} instead
    // to run all fit iterations
    using KinFitter::DoFit;

    static tree_t MakeTree(ParticleTypeTree ptree);
    static unsigned CountGammas(ParticleTypeTree ptree);

    const tree_t tree;
    std::vector<tree_t> tree_leaves;

    using permutations_t = std::vector<std::vector<size_t>>;
    permutations_t permutations;

    using sum_daughters_t = std::vector<std::function<void()>>;
    sum_daughters_t sum_daughters;

    struct iteration_t {
        struct particle_t {
            particle_t(const TParticlePtr& p, int leaveIndex = -1) :
                Particle(p), LeaveIndex(leaveIndex)
            {}
            TParticlePtr Particle;
            int LeaveIndex;
        };

        std::vector<particle_t> Particles;
        // given by iterationFilter
        double QualityFactor = std_ext::NaN;
        // list::sort makes highest quality come first
        bool operator<(const iteration_t& o) const {
            return QualityFactor > o.QualityFactor;
        }
    };

    using iterations_t = std::list<iteration_t>;
    iterations_t iterations;

    unsigned           max_iterations = 0; // 0 means no filtering
    iteration_filter_t iteration_filter;

};


}}} // namespace ant::analysis::utils


