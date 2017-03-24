#include "TreeFitter.h"

#include "base/Logger.h"
#include "utils/ParticleTools.h"
#include "base/std_ext/string.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

TreeFitter::TreeFitter(ParticleTypeTree ptree,
                       UncertaintyModelPtr uncertainty_model,
                       bool fit_Z_vertex,
                       nodesetup_t::getter nodeSetup,
                       const APLCON::Fit_Settings_t& settings) :
    KinFitter(uncertainty_model, fit_Z_vertex, settings),
    tree(MakeTree(ptree))
{
    // the tree fitter knows already the number of photons from the tree
    Photons.resize(CountGammas(ptree), Proton); // copy from Proton (will be overriden by Set calls)

    tree->GetUniquePermutations(tree_leaves, permutations, i_leaf_offset);

    // if the leaf offset is not 1, there's more in the tree than the proton
    if(i_leaf_offset > 1)
        throw Exception("Given particle type tree is too complex");
    if(i_leaf_offset == 1 && tree_leaves.front()->Get().TypeTree->Get() != ParticleTypeDatabase::Proton)
        throw Exception("Proton in final state expected");

    LOG(INFO) << "Initialized TreeFitter for " << ParticleTools::GetDecayString(ptree, false)
              << " with " << permutations.size() << " permutations, including KinFit";

    if(i_leaf_offset==1) {
        // check if the first leaf
        node_t& proton_leaf = tree_leaves.front()->Get();

        if(proton_leaf.TypeTree->GetParent()->Get() != ParticleTypeDatabase::BeamTarget) {
            // connect the proton as very first leave
            proton_leaf.Leaf = addressof(Proton);
        }
        else {
            // handle this tree as if the proton isn't present
            // pop the first leaf, the proton (may improve performance)
            i_leaf_offset = 0;
            tree_leaves.erase(tree_leaves.begin());
        }
    }

    for(unsigned i=0;i<Photons.size();i++) {
        node_t& photon_leaf = tree_leaves[i+i_leaf_offset]->Get();
        photon_leaf.Leaf = addressof(Photons[i]); // link via pointer
    }

    // prepare the calculation at each constraint
    // the Map_nodes call can be done once and the
    // calculation is stored in little functions

    tree->Map_nodes([this, nodeSetup] (const tree_t& tnode) {
        // do not include leaves
        if(tnode->IsLeaf())
            return;

        // do not include beamparticles
        if(tnode->Get().TypeTree->Get() == ParticleTypeDatabase::BeamTarget)
            return;

        // always sum up the tree nodes
        sum_daughters.emplace_back([tnode] () {
            node_t& node = tnode->Get();
            node.LVSum = LorentzVec{{0,0,0},0};
            assert(!tnode->Daughters().empty());
            for(const auto& d : tnode->Daughters())
                node.LVSum += d->Get().LVSum;
        });

        const nodesetup_t& setup = nodeSetup(tnode->Get().TypeTree);
        if(setup.Excluded)
            return;

        const auto IM_Sigma = setup.IM_Sigma;
        LOG(INFO) << "IM constraint for " << tnode->Get().TypeTree->Get().Name()
                  << " with sigma=" << IM_Sigma;
        node_constraints.emplace_back([tnode, IM_Sigma] () {
            node_t& node = tnode->Get();
            const double IM_calc = node.LVSum.M();
            const double IM_expected = tnode->Get().TypeTree->Get().Mass();
            return (IM_expected - IM_calc)/IM_Sigma;
        });
    });

    LOG(INFO) << "Have " << node_constraints.size() << " constraints at " << sum_daughters.size() << " nodes";

}

void TreeFitter::PrepareFits(double ebeam,
                             const TParticlePtr& proton,
                             const TParticleList& photons)
{
    // check for photon multiplicity
    if(photons.size() != Photons.size())
        throw Exception(std_ext::formatter()
                        << "TreeFitter: Given number of photons " << photons.size()
                        << " does not match expected " << Photons.size());

    // prepare the underlying kinematic fit
    // this may also set the proton's kinetic energy to missing E
    // do some more checks
    KinFitter::PrepareFit(ebeam, proton, photons);

    // iterations should normally be empty at this point,
    // but the user might call PrepareFits multiple times before running NextFit
    // (for whatever reason...)
    iterations.clear();

    for(const auto& current_perm : permutations)
    {
        iterations.emplace_back();
        iteration_t& it = iterations.back();

        for(unsigned i=0;i<Photons.size();i++) {
            const auto perm_idx = current_perm.at(i);
            it.Photons.emplace_back(photons.at(perm_idx), perm_idx);
        }
    }

    // filter iterations if requested
    if(!iteration_filter)
        return;

    for(auto& it : iterations) {

        PrepareFit(it);

        // after PrepareFit, we can obtain the initial LVSum now from GetLorentzVec
        do_sum_daughters();

        it.QualityFactor = iteration_filter();
    }

    // remove all iterations with factor=0
    iterations.remove_if([] (const iteration_t& it) {
        return it.QualityFactor == 0;
    });

    // if requested, keep only max best iterations
    if(max_iterations>0 && max_iterations<=iterations.size()) {
        iterations.sort();
        iterations.resize(max_iterations);
    }
}

void TreeFitter::PrepareFit(const TreeFitter::iteration_t& it)
{
    // update the current leave index,
    // gather the photons (in the right permuation!)
    // for the KinFitter, which actually sets the FitParticles in this order!
    TParticleList photons;
    for(unsigned i=0; i<Photons.size(); i++) {
        const auto& p = it.Photons.at(i);
        node_t& photon_leaf = tree_leaves[i+i_leaf_offset]->Get();
        photon_leaf.PhotonLeafIndex = p.LeafIndex;
        photons.emplace_back(p.Particle);
    }

    KinFitter::PrepareFit(BeamE.Value_before, Proton.Particle, photons);
}

void TreeFitter::do_sum_daughters() const
{
    // set the leaf LVSum to the leaf's LorentzVec
    for(auto& leave : tree_leaves) {
        node_t& node = leave->Get();
        node.LVSum = node.Leaf->GetLorentzVec();
    }

    // sum the daughters
    for(const auto& f : sum_daughters)
        f();
}

std::vector<double> TreeFitter::constraintIMatNodes() const
{
    do_sum_daughters();

    // calculate the IM constraint by evaluating the pre-defined functions
    vector<double> IM_diff(node_constraints.size());
    for(unsigned i=0;i<node_constraints.size();i++)
        IM_diff[i] = node_constraints[i]();
    return IM_diff;
}

bool TreeFitter::NextFit(APLCON::Result_t& fit_result)
{
    if(iterations.empty())
        return false;
    PrepareFit(iterations.front());

    auto wrap_constraintIMatNodes = [this] (const BeamE_t&, const Proton_t&, const Photons_t&, const Z_Vertex_t&) {
        return this->constraintIMatNodes();
    };

    fit_result = aplcon.DoFit(BeamE, Proton, Photons, Z_Vertex,
                              KinFitter::constraintEnergyMomentum,
                              wrap_constraintIMatNodes
                              );

    iterations.pop_front();
    return true;
}

void TreeFitter::SetIterationFilter(TreeFitter::iteration_filter_t filter, unsigned max)
{
    iteration_filter = filter;
    max_iterations = max;
}

TreeFitter::tree_t TreeFitter::GetTreeNode(const ParticleTypeDatabase::Type& type) const {
    auto nodes = GetTreeNodes(type);
    return nodes.empty() ? nullptr : nodes.front();
}

std::vector<TreeFitter::tree_t> TreeFitter::GetTreeNodes(const ParticleTypeDatabase::Type& type) const {
    std::vector<tree_t> nodes;
    tree->Map_nodes([&nodes, &type] (tree_t t) {
        if(t->Get().TypeTree->Get() == type)
            nodes.emplace_back(t);
    });
    return nodes;
}

TreeFitter::tree_t TreeFitter::MakeTree(ParticleTypeTree ptree)
{
    auto t =  ptree->DeepCopy<node_t>([] (const ParticleTypeTree& n) { return n; });
    t->Sort();
    return t;
}

unsigned TreeFitter::CountGammas(ParticleTypeTree ptree)
{
    unsigned nGammas = 0;
    ptree->Map([&nGammas] (const ParticleTypeDatabase::Type& n) {
        if(n == ParticleTypeDatabase::Photon)
            nGammas++;
    });
    return nGammas;
}


