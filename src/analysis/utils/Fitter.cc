#include "Fitter.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"

#include "utils/particle_tools.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/ParticleType.h"
#include "base/Logger.h"

#include "APLCON.hpp" // external project

#include "TTree.h"

#include <cassert>
#include <functional>
#include <cmath>
#include <algorithm>
#include <cstdlib>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;

const APLCON::Fit_Settings_t Fitter::Fitter::DefaultSettings = Fitter::MakeDefaultSettings();

Fitter::Fitter(const string& fittername, const APLCON::Fit_Settings_t& settings, utils::UncertaintyModelPtr uncertainty_model):
    uncertainty(uncertainty_model)
{
    aplcon = std_ext::make_unique<APLCON>(fittername, settings);
}

Fitter::~Fitter()
{}


APLCON::Fit_Settings_t Fitter::MakeDefaultSettings()
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 30;
    return settings;
}

void Fitter::LinkVariable(Fitter::FitParticle& particle)
{
    aplcon->LinkVariable(particle.Name,
                         particle.Addresses(),
                         particle.Addresses_Sigma(),
                         particle.Addresses_Pulls());
}

void Fitter::FitParticle::SetEkThetaPhi(const TParticlePtr& p,
                                        const UncertaintyModelPtr& uncertainty)
{
    const auto sigmas = uncertainty->GetSigmas(*p);

    Particle = p;

    Ek.SetValue(p->Ek());
    Theta.SetValue(p->Theta());
    Phi.SetValue(p->Phi());

    Ek.SetSigma(sigmas.sigmaE);
    Theta.SetSigma(sigmas.sigmaTheta);
    Phi.SetSigma(sigmas.sigmaPhi);
}

double Fitter::fct_TaggerEGausSigma(double)
{
    return  3.0/sqrt(12.0);
}


void Fitter::FitParticle::Var_t::SetupBranches(TTree* tree, const string& prefix)
{
    tree->Branch(prefix.c_str(), addressof(Value));
    tree->Branch((prefix+"_pull").c_str(), addressof(Pull));
    tree->Branch((prefix+"_sigma").c_str(), addressof(Sigma));
}

Fitter::FitParticle::FitParticle(const string& name, std::shared_ptr<Fitter::FitVariable> z_vertex) :
    Name(name), Z_Vertex(z_vertex)
{

}

Fitter::FitParticle::~FitParticle()
{

}

TParticlePtr Fitter::FitParticle::AsFitted() const
{
    const auto z_vertex = Z_Vertex ? Z_Vertex->Value : 0;

    auto p = make_shared<TParticle>(Particle->Type(),
                                    GetVector({Ek.Value, Theta.Value, Phi.Value},
                                              z_vertex
                                              ));
    p->Candidate = Particle->Candidate;
    return p;
}

void Fitter::FitParticle::SetupBranches(TTree* tree, const string& prefix)
{
    Ek.SetupBranches(tree, prefix+"_"+Name+"_Ek");
    Theta.SetupBranches(tree, prefix+"_"+Name+"_Theta");
    Phi.SetupBranches(tree, prefix+"_"+Name+"_Phi");
}

LorentzVec Fitter::FitParticle::GetVector(const std::vector<double>& EkThetaPhi,
                                          const double z_vertex) const
{
    const mev_t& Ek = EkThetaPhi[0];
    const radian_t& theta = EkThetaPhi[1];
    const radian_t& phi   = EkThetaPhi[2];


    const mev_t& E = Ek + Particle->Type().Mass();
    const mev_t& p = sqrt( sqr(E) - sqr(Particle->Type().Mass()) );

    if(z_vertex == 0.0)
        return LorentzVec::EPThetaPhi(E, p, theta, phi);


    double theta_corr = std_ext::NaN;

    if(!Particle->Candidate)
        throw Exception("Z Vertex fitting requires particles with candidates");

    auto calocluster = Particle->Candidate->FindCaloCluster();

    if(!calocluster)
        throw Exception("No calo cluster found");

    if(calocluster->DetectorType == Detector_t::Type_t::CB) {
        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        auto elem = cb->GetClusterElement(calocluster->CentralElement);
        const auto R  = cb->GetInnerRadius() + elem->RadiationLength*std::log2(Ek/elem->CriticalE);
        theta_corr = std::acos(( R*std::cos(theta) - z_vertex) / R );
    }
    else if(calocluster->DetectorType == Detector_t::Type_t::TAPS) {
        static auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
        auto elem = taps->GetClusterElement(calocluster->CentralElement);
        const auto Z  =  taps->GetZPosition() + elem->RadiationLength*std::log2(Ek/elem->CriticalE);
        theta_corr = std::atan( Z*std::tan(theta) / (Z - z_vertex));
    }
    else {
        throw Exception("Unknown calo cluster encountered");
    }

    return LorentzVec::EPThetaPhi(E, p, theta_corr, phi);
}



/**
 * @brief KinFitter::KinFitter
 * @param name Name for the fitter
 * @param numGammas number of photons involved in the fit
 * @param Uncertainty_model model to predict uncertainties
 * @param settings
 */
KinFitter::KinFitter(const std::string& name,
                     unsigned numGammas,
                     utils::UncertaintyModelPtr Uncertainty_model,
                     bool fit_Z_vertex,
                     const APLCON::Fit_Settings_t& settings) :
    Fitter(name, settings, Uncertainty_model)
{
    if(numGammas==0)
        throw Exception("No gammas are not allowed");

    if(fit_Z_vertex)
        Z_Vertex = std::make_shared<Z_Vertex_t>();

    for(unsigned i=0; i<numGammas;++i) {
        Photons.emplace_back(make_shared<FitParticle>("Photon"+to_string(i), Z_Vertex));
    }



    Proton = std::make_shared<FitParticle>("Proton", Z_Vertex);
    LinkVariable(*Proton);

    vector<string> variable_names      = {Proton->Name};
    vector<std::shared_ptr<FitParticle>> fit_particles{Proton};

    for ( auto& photon: Photons)
    {
        LinkVariable(*photon);
        variable_names.emplace_back(photon->Name);
        fit_particles.emplace_back(photon);
    }

    BeamE = std_ext::make_unique<BeamE_t>();
    aplcon->LinkVariable(BeamE->Name,
                         {std::addressof(BeamE->Value)},
                         {std::addressof(BeamE->Sigma)},
                         {std::addressof(BeamE->Pull)}
                         );
    variable_names.emplace_back(BeamE->Name);

    if(fit_Z_vertex) {
        aplcon->LinkVariable(Z_Vertex->Name,
                             {std::addressof(Z_Vertex->Value)},
                             {std::addressof(Z_Vertex->Sigma)},
                             {std::addressof(Z_Vertex->Pull)}
                             );
        variable_names.emplace_back(Z_Vertex->Name);
    }

    auto LorentzInvariance = [fit_Z_vertex, fit_particles] (const vector<vector<double>>& values)
    {

        const auto  n = fit_particles.size();
        // n serves as an offset here
        const auto& Ebeam    = values[n+0][0];
        const auto  z_vertex = fit_Z_vertex ? values[n+1][0] : 0.0;

        // Beam-LV:
        // beam    LorentzVec(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
        // target  LorentzVec(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
        const LorentzVec beam(0, 0, Ebeam, Ebeam);
        const LorentzVec target(0,0,0,ParticleTypeDatabase::Proton.Mass());

        LorentzVec constraint = target + beam;

        for(size_t i=0;i<n;i++)
            constraint -= fit_particles[i]->GetVector(values[i], z_vertex);

        return vector<double>(
               { constraint.p.x,
                 constraint.p.y,
                 constraint.p.z,
                 constraint.E }
               );
    };

    aplcon->AddConstraint("LInv", variable_names, LorentzInvariance);

}

KinFitter::~KinFitter()
{

}

void KinFitter::SetEgammaBeam(const double ebeam)
{
    BeamE->Value        = ebeam;
    BeamE->Value_before = ebeam;
    BeamE->Sigma = fct_TaggerEGausSigma(ebeam);
}

void KinFitter::SetZVertexSigma(double sigma)
{
    if(!Z_Vertex)
        throw Exception("Z Vertex fitting not enabled");
    Z_Vertex->Sigma = sigma;
    Z_Vertex->Sigma_before = sigma;
}

void KinFitter::SetProton(const TParticlePtr& proton)
{
    Proton->SetEkThetaPhi(proton, uncertainty);
}

void KinFitter::SetPhotons(const TParticleList& photons)
{
    if(Photons.size() != photons.size())
        throw Exception("Given number of photons does not match configured fitter");

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        Photons[i]->SetEkThetaPhi(photons[i], uncertainty);
    }
}

TParticlePtr KinFitter::GetFittedProton() const
{
    return Proton->AsFitted();
}

TParticleList KinFitter::GetFittedPhotons() const
{
    TParticleList photons;
    for(unsigned i=0;i<Photons.size();i++) {
        photons.emplace_back(Photons[i]->AsFitted());
    }
    return photons;
}

double KinFitter::GetFittedBeamE() const
{
    return BeamE->Value;
}

double KinFitter::GetFittedZVertex() const
{
    if(Z_Vertex)
        return Z_Vertex->Value;
    else
        return std_ext::NaN; // ignore silently
}

double KinFitter::GetBeamEPull() const
{
    return BeamE->Pull;
}

double KinFitter::GetZVertexPull() const
{
    if(Z_Vertex)
        return Z_Vertex->Pull;
    else
        return std_ext::NaN; // ignore silently
}

double KinFitter::GetProtonEPull() const
{
    return Proton->Ek.Pull;
}

double KinFitter::GetProtonThetaPull() const
{
    return Proton->Theta.Pull;

}

double KinFitter::GetProtonPhiPull() const
{
    return Proton->Phi.Pull;
}

std::vector<double> KinFitter::GetPhotonEPulls() const
{
    std::vector<double> pulls;
    for(auto& photon : Photons)
        pulls.push_back(photon->Ek.Pull);
    return pulls;
}

std::vector<double> KinFitter::GetPhotonThetaPulls() const
{
    std::vector<double> pulls;
    for(auto& photon : Photons)
        pulls.push_back(photon->Theta.Pull);
    return pulls;
}

std::vector<double> KinFitter::GetPhotonPhiPulls() const
{
    std::vector<double> pulls;
    for(auto& photon : Photons)
        pulls.push_back(photon->Phi.Pull);
    return pulls;
}

std::vector<Fitter::FitParticle> KinFitter::GetFitParticles() const
{
    std::vector<Fitter::FitParticle> particles{*Proton};
    for(auto& photon : Photons)
        particles.emplace_back(*photon);
    return particles;
}



void KinFitter::SetupBranches(TTree* tree, string branch_prefix)
{
    if(branch_prefix.empty())
        branch_prefix = aplcon->GetName();

    Proton->SetupBranches(tree, branch_prefix);
    for(auto& p : Photons) {
        p->SetupBranches(tree, branch_prefix);
    }

    tree->Branch((branch_prefix+"_chi2dof").c_str(),     &result_chi2ndof);
    tree->Branch((branch_prefix+"_iterations").c_str(),  &result_iterations);
    tree->Branch((branch_prefix+"_status").c_str(),      &result_status);
    tree->Branch((branch_prefix+"_probability").c_str(), &result_probability);
    tree->Branch((branch_prefix+"_EBeam").c_str(),       &BeamE->Value);
    tree->Branch((branch_prefix+"_EBeam_Pull").c_str(),  &BeamE->Pull);
    tree->Branch((branch_prefix+"_EBeam_Sigma").c_str(),  &BeamE->Sigma);
}

APLCON::Result_t KinFitter::DoFit() {
    if(Z_Vertex) {
        if(!std::isfinite(Z_Vertex->Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Z_Vertex->Value = 0;
        Z_Vertex->Sigma = Z_Vertex->Sigma_before;
    }

    const auto res = aplcon->DoFit();

    result_chi2ndof    = res.ChiSquare / res.NDoF;
    result_iterations  = res.NIterations;
    result_status      = static_cast<int>(res.Status);
    result_probability = res.Probability;

    return res;
}

TreeFitter::TreeFitter(const string& name,
                       ParticleTypeTree ptree,
                       utils::UncertaintyModelPtr uncertainty_model,
                       bool fit_Z_vertex,
                       nodesetup_t::getter nodeSetup,
                       const APLCON::Fit_Settings_t& settings) :
    KinFitter(name, CountGammas(ptree), uncertainty_model, fit_Z_vertex, settings),
    tree(MakeTree(ptree))
{
    tree->GetUniquePermutations(tree_leaves, permutations);
    current_perm = permutations.end();

    LOG(INFO) << "Initialized TreeFitter '" << name
              << "' for " << utils::ParticleTools::GetDecayString(ptree, false)
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
    using sum_daughters_t = function<void()>;
    using node_constraint_t = function<double()>;

    vector<sum_daughters_t> sum_daughters;
    vector<node_constraint_t> node_constraints;

    tree->Map_nodes([&sum_daughters, &node_constraints, nodeSetup] (const tree_t& tnode) {
        // do not include leaves
        if(tnode->IsLeaf())
            return;

        // do not include beamparticles
        if(tnode->Get().TypeTree->Get() == ParticleTypeDatabase::BeamTarget)
            return;

        // always sum up the tree nodes
        sum_daughters.emplace_back([tnode] () {
            node_t& node = tnode->Get();
            node.LVSum = LorentzVec{0,0,0,0};
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
            return (IM_calc - IM_expected)/IM_Sigma;
        });
    });

    LOG(INFO) << "Have " << node_constraints.size() << " constraints at " << sum_daughters.size() << " nodes";

    // define the constraint
    auto tree_leaves_copy = tree_leaves;
    auto IM_at_nodes = [fit_Z_vertex, tree_leaves_copy,
                       sum_daughters, node_constraints] (const vector<vector<double>>& v) {
        const auto  n = tree_leaves_copy.size();
        // n serves as an offset here
        const auto  z_vertex = fit_Z_vertex ? v[n+0][0] : 0.0;

        // assign values v to leaves' LVSum
        for(unsigned i=0;i<n;i++) {
            node_t& node = tree_leaves_copy[i]->Get();
            node.LVSum = node.Leave->GetVector(v[i], z_vertex);
        }

        // sum daughters' Particle
        for(const auto f : sum_daughters)
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

    current_perm = permutations.begin();
    current_comb_ptr = std_ext::make_unique<current_comb_t>(photons, current_perm->size());
}

bool TreeFitter::NextFit(APLCON::Result_t& fit_result)
{
    assert(!permutations.empty());

    if(!current_comb_ptr)
        return false;

    auto& current_comb = *current_comb_ptr;

    if(current_perm == permutations.end()) {
        current_perm = permutations.begin();
        ++current_comb;
    }

    if(current_comb.Done())
        return false;

    const auto k = current_perm->size();
    const auto n = Photons.size();

    // by construction, the photon leaves are 0..k-1
    const auto& comb_indices = current_comb.Indices();
    for(unsigned i=0;i<k;i++) {
        const auto perm_idx = current_perm->at(i);
        tree_leaves[i]->Get().PhotonLeaveIndex = comb_indices[perm_idx];
        const TParticlePtr& p = current_comb.at(perm_idx);
        Photons[i]->SetEkThetaPhi(p, uncertainty);
    }

    auto it_not_comb = current_comb.begin_not();
    assert(n>=k);
    assert((unsigned)distance(it_not_comb, current_comb.end_not()) == n - k);

    // and by construction, the non-leaves are from k..n-1
    for(auto i=k;i<n;i++) {
        Photons[i]->SetEkThetaPhi(*it_not_comb, uncertainty);
        ++it_not_comb;
    }

    // restore previous values
    SetProton(Proton->Particle);
    SetEgammaBeam(BeamE->Value_before);

    fit_result = KinFitter::DoFit();

    ++current_perm;

    return true;
}

TreeFitter::tree_t TreeFitter::GetTreeNode(const ParticleTypeDatabase::Type& type) const {
    tree_t treenode = nullptr;
    tree->Map_nodes([&treenode, &type] (tree_t t) {
        if(t->Get().TypeTree->Get() == type)
            treenode = t;
    });
    return treenode;
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
