#include "TreeFitter.h"

#include "base/Logger.h"
#include "utils/particle_tools.h"

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
    tree->GetUniquePermutations(tree_leaves, permutations, i_leave_offset);

    // if the leave offset is not 1, there's more in the tree than the proton
    if(i_leave_offset > 1)
        throw Exception("Given particle type tree is too complex");
    if(i_leave_offset == 1 && tree_leaves.front()->Get().TypeTree->Get() != ParticleTypeDatabase::Proton)
        throw Exception("Proton in final state expected");

    // expect tree_leaves correspond to photons, module possibly present proton
    if(tree_leaves.size()-i_leave_offset != Photons.size())
        throw Exception("Something went wront in underlying KinFitter?");


    LOG(INFO) << "Initialized TreeFitter '" << name
              << "' for " << ParticleTools::GetDecayString(ptree, false)
              << " with " << permutations.size() << " permutations, including KinFit";

    // setup fitter variables, collect leave names and particles for constraint
    // similar to KinFitter ctor
    vector<string> variable_names;

    if(i_leave_offset==1) {
        node_t& proton_node = tree_leaves.front()->Get();
        if(proton_node.TypeTree->GetParent()->Get() != ParticleTypeDatabase::BeamTarget) {
            // connect the proton as very first leave
            // only if the proton is not a daughter of the (pseudo) beam particle
            proton_node.Leave = Proton;
            variable_names.emplace_back(proton_node.Leave->Name);
        }
        else {
            // handle this tree as if the proton isn't present
            i_leave_offset = 0;
            tree_leaves.erase(tree_leaves.begin()); // pop the first leave, the proton
        }
    }

    for(unsigned i=0;i<Photons.size();i++) {
        node_t& photon_node = tree_leaves[i+i_leave_offset]->Get();
        photon_node.Leave = Photons[i]; // link via shared_ptr
        variable_names.emplace_back(photon_node.Leave->Name);
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

void TreeFitter::PrepareFits(double ebeam,
                             const TParticlePtr& proton,
                             const TParticleList& photons)
{
    if(photons.size() != Photons.size())
        throw Exception("Given leave particles does not match configured TreeFitter");

    // prepare the underlying kinematic fit
    PrepareFit(ebeam, proton, photons);

    // iterations should normally be empty at this point,
    // but the user might call SetPhotons multiple times before running NextFit
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
        // set the leaves' LVSum and sum up the tree (as in the constraint)
        // the filtering function might use this then to calculate the quality factor
        for(unsigned i=0; i<Photons.size();i++) {
            const auto& p = it.Photons[i];
            node_t& photon_leave = tree_leaves[i+i_leave_offset]->Get();
            photon_leave.LVSum = *p.Particle;
            photon_leave.Leave->Particle = p.Particle;
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

    TParticleList photons;
    for(unsigned i=0; i<Photons.size(); i++) {
        const auto& p = it.Photons.at(i);
        node_t& photon_leave = tree_leaves[i+i_leave_offset]->Get();
        photon_leave.PhotonLeaveIndex = p.LeaveIndex;
        photons.emplace_back(p.Particle);
    }

    fit_result = KinFitter::DoFit(BeamE->Value_before, Proton->Particle, photons);

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
