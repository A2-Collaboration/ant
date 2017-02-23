#pragma once

#include "KinFitter.h"

#include "base/ParticleTypeTree.h"

namespace ant {
namespace analysis {
namespace utils {

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


    /**
     * @brief TreeFitter creates fitter for the given ParticleTypeTree.
     * Applies IM constraint to each node, can be customized by optional nodeSetup.
     * Supports only photon leaves at the moment.
     *
     * @param name unique name for this fitter
     * @param ptree tree describing the decay to be fitted
     * @param uncertainty_model the uncertainties, see KinFitter
     * @param fit_Z_vertex make z vertex unfixed (not equal to zero)
     * @param nodeSetup fine-grained control over each node in the particle type tree
     * @param settings fit settings for APLCON (iterations, epsilons, ...)
     */
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
        int PhotonLeaveIndex = -1; // according to list given by SetPhotons, or -1 if proton
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

    using permutation_t = std::vector<int>;
    std::vector<permutation_t> permutations;

    int i_leave_offset;

    using sum_daughters_t = std::vector<std::function<void()>>;
    sum_daughters_t sum_daughters;

    struct iteration_t {
        struct photon_t {
            photon_t(const TParticlePtr& p, int leaveIndex) :
                Particle(p), LeaveIndex(leaveIndex)
            {}
            TParticlePtr Particle;
            int LeaveIndex;
        };

        std::vector<photon_t> Photons;
        // given by iterationFilter (if defined by user)
        double QualityFactor = std_ext::NaN;
        // list::sort makes highest quality come first
        bool operator<(const iteration_t& o) const {
            return QualityFactor > o.QualityFactor;
        }
    };

    std::list<iteration_t> iterations;

    unsigned           max_iterations = 0; // 0 means no filtering
    iteration_filter_t iteration_filter;

};

}}} // namespace ant::analysis::utils