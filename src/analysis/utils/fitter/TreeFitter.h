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
        bool Excluded;
        explicit nodesetup_t(bool excluded = false) :
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
     * Supports only photon leaf permutations at the moment (but single extra
     * final state particles are ok).
     *
     * @param ptree tree describing the decay to be fitted
     * @param uncertainty_model the uncertainties, see KinFitter
     * @param fit_Z_vertex make z vertex unfixed (not equal to zero)
     * @param nodeSetup fine-grained control over each node in the particle type tree
     * @param settings fit settings for APLCON (iterations, epsilons, ...)
     */
    TreeFitter(ParticleTypeTree ptree,
               UncertaintyModelPtr uncertainty_model,
               bool fit_Z_vertex = false,
               nodesetup_t::getter nodeSetup = {},
               const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    TreeFitter(const TreeFitter&) = delete;
    TreeFitter& operator=(const TreeFitter&) = delete;
    TreeFitter(TreeFitter&&) = default;
    TreeFitter& operator=(TreeFitter&&) = default;

    void PrepareFits(double ebeam,
                     const TParticlePtr& proton,
                     const TParticleList& photons);

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
        const ParticleTypeTree TypeTree;     // the underlying type tree (from which the tree this node belongs to was created)
        LorentzVec LVSum;                    // if Leave==nullptr, the sum of the daughters
        const FitParticle* Leaf = nullptr;  // observer pointer for assigned fitparticle leaf (is permuted by NextFit)
        int PhotonLeafIndex = -1;           // according to photon list given by PrepareFits, or -1 if proton

        // sorted by our underlying type tree
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

    // force usage of "PrepareFits(...)" and "while(NextFit()) {}" interface
    using KinFitter::DoFit;

    static tree_t MakeTree(ParticleTypeTree ptree);
    static unsigned CountGammas(ParticleTypeTree ptree);

    const tree_t tree;
    std::vector<tree_t> tree_leaves;

    using permutation_t = std::vector<int>;
    std::vector<permutation_t> permutations;

    int i_leaf_offset;

    using sum_daughters_t = std::vector<std::function<void()>>;
    sum_daughters_t sum_daughters;

    using node_constraint_t = std::function<double()>;
    std::vector<node_constraint_t> node_constraints;

    struct iteration_t {
        struct photon_t {
            photon_t(const TParticlePtr& p, int leafIndex) :
                Particle(p), LeafIndex(leafIndex)
            {}
            TParticlePtr Particle;
            int LeafIndex;
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

    void PrepareFit(const iteration_t& it);

    unsigned           max_iterations = 0; // 0 means no filtering
    iteration_filter_t iteration_filter;

    void do_sum_daughters() const;

    // this constraint needs stuff from the class instance
    // it's gonna be wrapped in a lambda for the DoFit call
    std::vector<double> constraintIMatNodes() const;
};

}}} // namespace ant::analysis::utils