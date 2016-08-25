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

Fitter::FitParticle::FitParticle(const string& name,
                                 APLCON& aplcon,
                                 std::shared_ptr<Fitter::FitVariable> z_vertex) :
    Detector(Detector_t::Any_t::None),
    Vars(4), // it's a lucky coincidence that any particle is parametrized by 4 values
    Name(name),
    Z_Vertex(z_vertex)
{
    auto vectorize = [] (vector<FitVariable>& vars, double FitVariable::* member) {
        vector<double*> ptrs(vars.size());
        transform(vars.begin(), vars.end(), ptrs.begin(),
                  [member] (FitVariable& v) { return std::addressof(v.*member); });
        return ptrs;
    };

    aplcon.LinkVariable(
                Name,
                vectorize(Vars, &FitVariable::Value),
                vectorize(Vars, &FitVariable::Sigma),
                vectorize(Vars, &FitVariable::Pull)
                );
}

Fitter::FitParticle::~FitParticle()
{

}

TParticlePtr Fitter::FitParticle::AsFitted() const
{
    const auto z_vertex = Z_Vertex ? Z_Vertex->Value : 0;

    auto vectorize = [] (const vector<FitVariable>& vars) {
        vector<double> values(vars.size());
        transform(vars.begin(), vars.end(), values.begin(),
                  [] (const FitVariable& v) { return v.Value; });
        return values;
    };

    auto p = make_shared<TParticle>(Particle->Type(),
                                    GetLorentzVec(vectorize(Vars), z_vertex)
                                    );
    *p *= 1000.0;
    p->Candidate = Particle->Candidate;
    return p;
}

void Fitter::FitParticle::Set(const TParticlePtr& p,
                              const UncertaintyModel& uncertainty)
{
    Particle = p;
    const auto sigmas = uncertainty.GetSigmas(*p);
    Detector = sigmas.Detector;

    const auto inverse_sigmaE = sigmas.sigmaE/p->Ek()/p->Ek();
    Vars[0].SetValueSigma(1000.0/p->Ek(), inverse_sigmaE*1000.0);
    Vars[2].SetValueSigma(p->Phi(),       sigmas.sigmaPhi);

    // the parametrization, and thus the meaning of the linked fitter variables,
    // depends on the calorimeter
    if(Detector & Detector_t::Type_t::CB)
    {
        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        Vars[1].SetValueSigma(p->Theta(), sigmas.sigmaTheta);
        const auto& CB_R = cb->GetInnerRadius() + sigmas.ShowerDepth;
        Vars[3].SetValueSigma(CB_R, sigmas.sigmaCB_R);
    }
    else if(Detector & Detector_t::Type_t::TAPS)
    {
        static auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
        const auto& TAPS_Lz = taps->GetZPosition() + sigmas.ShowerDepth;
        const auto& TAPS_Rxy = std::tan(p->Theta())*TAPS_Lz;
        Vars[1].SetValueSigma(TAPS_Rxy, sigmas.sigmaTAPS_Rxy);
        Vars[3].SetValueSigma(TAPS_Lz,  sigmas.sigmaTAPS_Lz);
    }
    else {
        throw Exception("Unknown/none detector type provided from uncertainty model");
    }
}

LorentzVec Fitter::FitParticle::GetLorentzVec(const std::vector<double>& vars,
                                              const double z_vertex) const
{
    double theta_corr = std_ext::NaN;

    if(Detector & Detector_t::Type_t::CB)
    {
        // for CB, parametrization is (Ek, theta, phi, CB_R)
        const radian_t& theta = vars[1];
        const auto&     CB_R  = vars[3];
        theta_corr = std::acos(( CB_R*std::cos(theta) - z_vertex) / CB_R );

    }
    else if(Detector & Detector_t::Type_t::TAPS)
    {
        // for TAPS, parametrization is (Ek, TAPS_Rxy, phi, TAPS_Lz)
        const auto& TAPS_Rxy = vars[1];
        const auto& TAPS_Lz  = vars[3];
        theta_corr = std::atan(TAPS_Rxy / (TAPS_Lz - z_vertex));
    }
    else {
        throw Exception("Unknown/none detector type provided from uncertainty model");
    }

    const mev_t& Ek       = 1.0/vars[0];
    const radian_t& phi   = vars[2];

    const mev_t& E = Ek + Particle->Type().Mass()/1000.0;
    const mev_t& p = sqrt( sqr(E) - sqr(Particle->Type().Mass()/1000.0) );

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
        Photons.emplace_back(make_shared<FitParticle>("Photon"+to_string(i), *aplcon, Z_Vertex));
    }

    Proton = std::make_shared<FitParticle>("Proton", *aplcon, Z_Vertex);

    vector<string> variable_names      = {Proton->Name};
    vector<std::shared_ptr<FitParticle>> fit_particles{Proton};

    for ( auto& photon: Photons)
    {
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

    auto EnergyMomentum = [fit_Z_vertex, fit_particles] (const vector<vector<double>>& values)
    {

        const auto  n = fit_particles.size();
        // n serves as an offset here
        const auto& BeamE    = values[n+0][0];
        const auto  z_vertex = fit_Z_vertex ? values[n+1][0] : 0.0;

        // start with the incoming particle
        auto diff = MakeBeamLorentzVec(BeamE);

        for(size_t i=0;i<n;i++)
            diff -= fit_particles[i]->GetLorentzVec(values[i], z_vertex); // minus outgoing

        return vector<double>(
               { diff.p.x,
                 diff.p.y,
                 diff.p.z,
                 diff.E }
               );
    };

    aplcon->AddConstraint("E-p", variable_names, EnergyMomentum);

}

KinFitter::~KinFitter()
{

}

void KinFitter::SetEgammaBeam(const double ebeam)
{
    BeamE->SetValueSigma(ebeam/1000.0, uncertainty->GetBeamEnergySigma(ebeam)/1000.0);
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
    Proton->Set(proton, *uncertainty);
}

void KinFitter::SetPhotons(const TParticleList& photons)
{
    if(Photons.size() != photons.size())
        throw Exception("Given number of photons does not match configured fitter");

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        Photons[i]->Set(photons[i], *uncertainty);
    }
}

bool KinFitter::IsZVertexFitEnabled() const noexcept
{
    return Z_Vertex != nullptr;
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
    return BeamE->Value*1000.0;
}

TParticlePtr KinFitter::GetFittedBeamParticle() const
{
    auto p = std::make_shared<TParticle>(ParticleTypeDatabase::BeamProton,
                                         MakeBeamLorentzVec(BeamE->Value));
    *p *= 1000.0;
    return p;
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
//    return Proton->Ek.Pull;
}

double KinFitter::GetProtonThetaPull() const
{
//    return Proton->Theta.Pull;
}

double KinFitter::GetProtonPhiPull() const
{
//    return Proton->Phi.Pull;
}

std::vector<double> KinFitter::GetPhotonEPulls() const
{
    std::vector<double> pulls;
//    for(auto& photon : Photons)
//        pulls.push_back(photon->Ek.Pull);
    return pulls;
}

std::vector<double> KinFitter::GetPhotonThetaPulls() const
{
    std::vector<double> pulls;
//    for(auto& photon : Photons)
//        pulls.push_back(photon->Theta.Pull);
    return pulls;
}

std::vector<double> KinFitter::GetPhotonPhiPulls() const
{
    std::vector<double> pulls;
//    for(auto& photon : Photons)
//        pulls.push_back(photon->Phi.Pull);
    return pulls;
}

std::vector<Fitter::FitParticle> KinFitter::GetFitParticles() const
{
    std::vector<Fitter::FitParticle> particles{*Proton};
    for(auto& photon : Photons)
        particles.emplace_back(*photon);
    return particles;
}


APLCON::Result_t KinFitter::DoFit() {
    if(Z_Vertex) {
        if(!std::isfinite(Z_Vertex->Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Z_Vertex->Value = 0;
        Z_Vertex->Sigma = Z_Vertex->Sigma_before;
    }

//    double missing_E = BeamE->Value;
//    for(auto& photon : Photons) {
//        missing_E -= photon->Particle->Ek()/1000.0;
//    }
//    // asumme that 0th component is Ek
//    Proton->Vars[0].Value = missing_E;
//    Proton->Vars[0].Value_before = missing_E;

    const auto res = aplcon->DoFit();

    result_chi2ndof    = res.ChiSquare / res.NDoF;
    result_iterations  = res.NIterations;
    result_status      = static_cast<int>(res.Status);
    result_probability = res.Probability;

    return res;
}

LorentzVec KinFitter::MakeBeamLorentzVec(double BeamE)
{
    // Beam Lorentz vector:
    // beam    LorentzVec(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
    // target  LorentzVec(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
    const LorentzVec beam({0, 0, BeamE}, BeamE);
    /// \todo Target is always assumed proton...
    const LorentzVec target({0,0,0}, ParticleTypeDatabase::Proton.Mass()/1000.0);

    return target + beam;
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
            const double IM_expected = tnode->Get().TypeTree->Get().Mass()/1000.0;
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
            node.LVSum = node.Leave->GetLorentzVec(v[i], z_vertex);
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
        Photons[i]->Set(p, *uncertainty);
    }

    auto it_not_comb = current_comb.begin_not();
    assert(n>=k);
    assert((unsigned)distance(it_not_comb, current_comb.end_not()) == n - k);

    // and by construction, the non-leaves are from k..n-1
    for(auto i=k;i<n;i++) {
        Photons[i]->Set(*it_not_comb, *uncertainty);
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
