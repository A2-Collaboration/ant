#include "TreeFitter.h"

#include "base/Logger.h"
#include "utils/particle_tools.h"
#include "utils/combinatorics.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

TreeFitter::TreeFitter(const string& name,
                       ParticleTypeTree ptree,
                       UncertaintyModelPtr uncertainty_model,
                       bool fit_Z_vertex,
                       nodesetup_t::getter nodeSetup,
                       const APLCON::Fit_Settings_t& settings) :
    KinFitter(name, CountGammas(ptree), uncertainty_model, fit_Z_vertex, settings),
    tree(MakeTree(ptree))
{
    tree->GetUniquePermutations(tree_leaves, permutations);

    LOG(INFO) << "Initialized TreeFitter '" << name
              << "' for " << ParticleTools::GetDecayString(ptree, false)
              << " with " << permutations.size() << " permutations, including KinFit";

    // setup fitter variables, collect leave names and particles for constraint
    // similar to KinFitter ctor
    vector<string> variable_names;
    for(unsigned i=0;i<tree_leaves.size();i++) {
        node_t& node = tree_leaves[i]->Get();
        node.Leave = Photons[i]; // link via shared_ptr
        variable_names.emplace_back(node.Leave->Name);
    }

    // the variable is already setup in KinFitter ctor
    if(fit_Z_vertex) {
        variable_names.emplace_back(Z_Vertex->Name);
    }

    // prepare the calculation at each constraint
    // the Map_nodes call can be done once and the
    // calculation is stored in little functions
    using node_constraint_t = function<double()>;

    vector<node_constraint_t> node_constraints;

    tree->Map_nodes([this, &node_constraints, nodeSetup] (const tree_t& tnode) {
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


    // define the constraint
    // use local copies instead of this capture. capturing "this" is dangerous since
    // fitter might be moved around in memory...
    auto tree_leaves_copy   = tree_leaves;
    auto sum_daughters_copy = sum_daughters;
    auto IM_at_nodes = [tree_leaves_copy, sum_daughters_copy,
                       fit_Z_vertex, node_constraints] (const vector<vector<double>>& v) {
        const auto  k = tree_leaves_copy.size();
        // k serves as an offset here
        const auto  z_vertex = fit_Z_vertex ? v[k+0][0] : 0.0;

        // assign values v to leaves' LVSum
        for(unsigned i=0;i<k;i++) {
            node_t& node = tree_leaves_copy[i]->Get();
            node.LVSum = node.Leave->GetLorentzVec(v[i], z_vertex);
        }

        // sum daughters' Particle
        for(const auto& f : sum_daughters_copy)
            f();

        // calculate the IM constraint by evaluating the pre-defined functions
        vector<double> IM_diff(node_constraints.size());
        for(unsigned i=0;i<node_constraints.size();i++)
            IM_diff[i] = node_constraints[i]();
        return IM_diff;
    };

    aplcon->AddConstraint("IM_at_nodes", variable_names, IM_at_nodes);
}

TreeFitter::~TreeFitter()
{}

void TreeFitter::SetPhotons(const TParticleList& photons)
{
    if(photons.size() != Photons.size())
        throw Exception("Given leave particles does not match configured TreeFitter");

    // iterations should normally be empty at this point,
    // but the user might call SetPhotons multiple times before running NextFit
    iterations.clear();

    const auto k = tree_leaves.size();
    const auto n = Photons.size();
    assert(k<=n);

    for(auto current_comb = makeCombination(photons, k);
        !current_comb.Done(); ++current_comb)
    {
        for(const auto& current_perm : permutations)
        {
            iterations.emplace_back();
            iteration_t& it = iterations.back();

            const auto& comb_indices = current_comb.Indices();
            for(unsigned i=0;i<k;i++) {
                const auto perm_idx = current_perm.at(i);
                it.Particles.emplace_back(current_comb.at(perm_idx), comb_indices[perm_idx]);
            }

            auto it_not_comb = current_comb.begin_not();
            assert((unsigned)distance(it_not_comb, current_comb.end_not()) == n - k);

            // and by construction, the non-leaves are from k..n-1
            for(auto i=k;i<n;i++) {
                it.Particles.emplace_back(*it_not_comb);
                ++it_not_comb;
            }
        }
    }

    // filter iterations if requested
    if(!iteration_filter)
        return;

    for(auto& it : iterations) {
        // set the leaves' LVSum and sum up the tree (as in the constraint)
        // the filtering function might use this to calculate the quality factor
        for(unsigned i=0; i<k;i++) {
            const auto& p = it.Particles[i];
            node_t& leave = tree_leaves[i]->Get();
            leave.LVSum = *p.Particle;
            leave.Leave->Particle = p.Particle;
        }

        // sum the daughters
        for(const auto& f : sum_daughters)
            f();

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

void TreeFitter::SetIterationFilter(TreeFitter::iteration_filter_t filter, unsigned max)
{
    iteration_filter = filter;
    max_iterations = max;
}

bool TreeFitter::NextFit(APLCON::Result_t& fit_result)
{
    if(iterations.empty())
        return false;

    const iteration_t& it = iterations.front();

    const auto k = tree_leaves.size();
    const auto n = Photons.size();

    for(unsigned i=0; i<n;i++) {
        const auto& p = it.Particles[i];
        // by construction, the photon leaves are 0..k-1
        if(i<k)
            tree_leaves[i]->Get().PhotonLeaveIndex = p.LeaveIndex;
        Photons[i]->Set(p.Particle, *uncertainty);
    }

    // restore previous values
    SetProton(Proton->Particle);
    SetEgammaBeam(BeamE->Value_before);

    fit_result = KinFitter::DoFit();

    iterations.pop_front();
    return true;
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

    auto t = ptree->DeepCopy<node_t>([] (const ParticleTypeTree& n) { return n; });

    if(t->Get().TypeTree->Get() == ParticleTypeDatabase::BeamTarget) {
        // do not use constref'ed daughter here, see ant::Tree::Unlink
        for(auto daughter : t->Daughters()) {
            if(daughter->Get().TypeTree->Get() == ParticleTypeDatabase::Nucleon) {
                LOG(INFO) << "Removing nucleon from tree";
                daughter->Unlink();
                break;
            }
        }
    }

    t->Sort();
    return t;
}

unsigned TreeFitter::CountGammas(ParticleTypeTree ptree)
{
    unsigned nGammas = 0;
    ptree->Map([&nGammas] (const ParticleTypeDatabase::Type& n) { if(n == ParticleTypeDatabase::Photon) nGammas++; });
    return nGammas;
}
