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

Fitter::Fitter(const string& fittername, const APLCON::Fit_Settings_t& settings,
               utils::UncertaintyModelPtr uncertainty_model) :
    uncertainty(uncertainty_model),
    aplcon(std_ext::make_unique<APLCON>(fittername, settings))
{}

Fitter::~Fitter()
{}


APLCON::Fit_Settings_t Fitter::MakeDefaultSettings()
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 30;
    settings.SkipCovariancesInResult = true;
    return settings;
}

vector<double*> vectorize(vector<Fitter::FitVariable>& vars, double Fitter::FitVariable::* member) {
    vector<double*> ptrs(vars.size());
    transform(vars.begin(), vars.end(), ptrs.begin(),
              [member] (Fitter::FitVariable& v) { return std::addressof(v.*member); });
    return ptrs;
}

Fitter::Vertex_t::Vertex_t(APLCON& aplcon) :
    XYZ(3)
{
    aplcon.LinkVariable(
                Name,
                vectorize(XYZ, &FitVariable::Value),
                vectorize(XYZ, &FitVariable::Sigma),
                vectorize(XYZ, &FitVariable::Pull)
                );
}

std::vector<double> Fitter::Vertex_t::GetXYZ() const
{
    return {XYZ[0].Value, XYZ[1].Value, XYZ[2].Value};
}

Fitter::FitParticle::FitParticle(const string& name,
                                 APLCON& aplcon,
                                 std::shared_ptr<const Vertex_t> vertex) :
    Detector(Detector_t::Any_t::None),
    Vars(4), // it's a lucky coincidence that particles in CB/TAPS are both parametrized by 4 values
    Name(name),
    Vertex(vertex)
{
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
    const auto& vertex = Vertex ? Vertex->GetXYZ() : vector<double>{0.0, 0.0, 0.0};

    vector<double> values(Vars.size());
    transform(Vars.begin(), Vars.end(), values.begin(),
                  [] (const FitVariable& v) { return v.Value; });

    auto p = make_shared<TParticle>(Particle->Type(),
                                    GetLorentzVec(values, vertex)
                                    );
    p->Candidate = Particle->Candidate;
    return p;
}

std::vector<double> Fitter::FitParticle::GetSigmas() const
{
    vector<double> sigmas(Vars.size());
    transform(Vars.begin(), Vars.end(), sigmas.begin(),
                  [] (const FitVariable& v) { return v.Sigma_before; });
    return sigmas;
}

Fitter::FitParticle::pulls_t Fitter::FitParticle::GetPulls() const
{
    pulls_t pulls(Vars.size());
    transform(Vars.begin(), Vars.end(), pulls.begin(),
                  [] (const FitVariable& v) { return v.Pull; });
    return pulls;
}

void Fitter::FitParticle::Set(const TParticlePtr& p,
                              const UncertaintyModel& uncertainty)
{
    Particle = p;
    const auto& sigmas = uncertainty.GetSigmas(*p);
    Detector = sigmas.Detector;

    if(!p->Candidate)
        throw Exception("Need particle with candidate for fitting");

    const auto inverseEk = 1000.0/p->Ek();
    const auto sigma_inverseEk = sigmas.sigmaEk*std_ext::sqr(inverseEk)/1000.0;
    Vars[0].SetValueSigma(inverseEk, sigma_inverseEk);

    Vars[2].SetValueSigma(p->Phi(), sigmas.sigmaPhi);

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
        const auto& TAPS_Lz = taps->GetZPosition() + sigmas.ShowerDepth*std::cos(p->Theta());
        auto& pos = p->Candidate->FindCaloCluster()->Position;
        const vec3 TAPS_L{pos.x, pos.y, TAPS_Lz};
        const auto& TAPS_Rxy = std::sin(p->Theta())*TAPS_L.R();

        Vars[1].SetValueSigma(TAPS_Rxy,   sigmas.sigmaTAPS_Rxy);
        Vars[3].SetValueSigma(TAPS_L.R(), sigmas.sigmaTAPS_L);
    }
    else {
        throw Exception("Unknown/none detector type provided from uncertainty model");
    }
}

LorentzVec Fitter::FitParticle::GetLorentzVec(const std::vector<double>& values,
                                              const std::vector<double>& vertex) const
{

    const radian_t& phi   = values[2];

    // start at the lab origin in the frame of the vertex
    vec3 x{-vertex.at(0),-vertex.at(1),-vertex.at(2)};

    if(Detector & Detector_t::Type_t::CB)
    {
        // for CB, parametrization is (Ek, theta, phi, CB_R)
        const radian_t& theta = values[1];
        const auto&     CB_R  = values[3];
        x += vec3::RThetaPhi(CB_R, theta, phi);
    }
    else if(Detector & Detector_t::Type_t::TAPS)
    {
        // for TAPS, parametrization is (Ek, TAPS_Rxy, phi, TAPS_L)
        const auto& TAPS_Rxy = values[1];
        const auto& TAPS_L   = values[3];
        const auto& TAPS_L_z = sqrt(std_ext::sqr(TAPS_L) - std_ext::sqr(TAPS_Rxy));
        x += vec3(vec2::RPhi(TAPS_Rxy, phi), TAPS_L_z);
    }
    else {
        throw Exception("Unknown/none detector type provided from uncertainty model");
    }

    const mev_t& Ek = 1000.0/values[0];
    const mev_t& E = Ek + Particle->Type().Mass();
    const mev_t& p = sqrt( sqr(E) - sqr(Particle->Type().Mass()) );

    const vec3& p_vec = x*p/x.R();
    return {p_vec, E};
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
                     bool fit_vertex,
                     const APLCON::Fit_Settings_t& settings) :
    Fitter(name, settings, Uncertainty_model)
{
    if(numGammas==0)
        throw Exception("No gammas are not allowed");

    BeamE = std_ext::make_unique<BeamE_t>();
    aplcon->LinkVariable(BeamE->Name,
                         {std::addressof(BeamE->Value)},
                         {std::addressof(BeamE->Sigma)},
                         {std::addressof(BeamE->Pull)}
                         );
    vector<string> variable_names;

    if(fit_vertex)
        Vertex = std::make_shared<Vertex_t>(*aplcon);

    vector<std::shared_ptr<FitParticle>> fit_particles;

    for(unsigned i=0; i<numGammas;++i) {
        auto photon = make_shared<FitParticle>("Photon"+to_string(i), *aplcon, Vertex);
        Photons.emplace_back(photon);
        variable_names.emplace_back(photon->Name);
        fit_particles.emplace_back(photon);
    }

    Proton = std::make_shared<FitParticle>("Proton", *aplcon, Vertex);
    variable_names.emplace_back(Proton->Name);
    fit_particles.emplace_back(Proton);

    variable_names.emplace_back(BeamE->Name);
    if(fit_vertex)
        variable_names.emplace_back(Vertex->Name);


    auto EnergyMomentum = [fit_vertex, fit_particles] (const vector<vector<double>>& values)
    {

        const auto  n = fit_particles.size();
        // n serves as an offset here
        const auto& BeamE    = values[n+0][0];
        const auto& vertex = fit_vertex ? values[n+1] : vector<double>{0.0, 0.0, 0.0};

        // start with the incoming particle
        auto diff = MakeBeamLorentzVec(1000.0/BeamE);

        for(size_t i=0;i<n;i++)
            diff -= fit_particles[i]->GetLorentzVec(values[i], vertex); // minus outgoing

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
    const auto invBeam = 1000.0/ebeam;
    const auto sigma_invBeam = uncertainty->GetBeamEnergySigma(ebeam)*std_ext::sqr(invBeam)/1000.0;
    BeamE->SetValueSigma(invBeam, sigma_invBeam);
}

void KinFitter::SetZVertexSigma(double sigma)
{
    if(!Vertex)
        throw Exception("Z Vertex fitting not enabled");
    Vertex->XYZ[2].SetValueSigma(0.0, sigma);
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
    return Vertex != nullptr;
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

TParticlePtr KinFitter::GetFittedBeamParticle() const
{
    return std::make_shared<TParticle>(ParticleTypeDatabase::BeamProton,
                                         MakeBeamLorentzVec(1000.0/BeamE->Value));
}

double KinFitter::GetFittedZVertex() const
{
    if(Vertex)
        return Vertex->XYZ[2].Value;
    else
        return std_ext::NaN; // ignore silently
}

double KinFitter::GetBeamEPull() const
{
    return BeamE->Pull;
}

double KinFitter::GetZVertexPull() const
{
    if(Vertex)
        return Vertex->XYZ[2].Pull;
    else
        return std_ext::NaN; // ignore silently
}

Fitter::FitParticle::pulls_t KinFitter::GetProtonPulls() const
{
    return Proton->GetPulls();
}


std::vector<std::vector<double>> KinFitter::GetPhotonsPulls() const
{
    /// \bug hardcoded number of parameters!
    std::vector<std::vector<double>> pulls(4, vector<double>(Photons.size()));
    for(unsigned i=0;i<Photons.size();i++) {
        auto p = Photons.at(i)->GetPulls();
        for(unsigned j=0;j<pulls.size();j++) {
            pulls.at(j).at(i) = p[j];
        }
    }
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
    if(Vertex) {
        if(!std::isfinite(Vertex->XYZ[2].Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Vertex->XYZ[0].SetValueSigma(0, 0.2);
        Vertex->XYZ[1].SetValueSigma(0, 0.2);
        Vertex->XYZ[2].SetValueSigma(0, Vertex->XYZ[2].Sigma_before);
    }


    // only set Proton Ek to missing energy if unmeasured
    auto& Var_invEk = Proton->Vars[0];
    if(Var_invEk.Sigma == 0) {
        double missing_E = 1.0/BeamE->Value;
        for(auto& photon : Photons) {
            missing_E -= 1.0/photon->Vars[0].Value;
        }
        Var_invEk.SetValueSigma(1.0/missing_E, Var_invEk.Sigma);
    }

    return aplcon->DoFit();
}

LorentzVec KinFitter::MakeBeamLorentzVec(double BeamE)
{
    // Beam Lorentz vector:
    // beam    LorentzVec(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
    // target  LorentzVec(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
    const LorentzVec beam({0, 0, BeamE}, BeamE);
    /// \todo Target is always assumed proton...
    const LorentzVec target({0,0,0}, ParticleTypeDatabase::Proton.Mass());

    return target + beam;
}


TreeFitter::TreeFitter(const string& name,
                       ParticleTypeTree ptree,
                       utils::UncertaintyModelPtr uncertainty_model,
                       bool fit_vertex,
                       nodesetup_t::getter nodeSetup,
                       const APLCON::Fit_Settings_t& settings) :
    KinFitter(name, CountGammas(ptree), uncertainty_model, fit_vertex, settings),
    tree(MakeTree(ptree))
{
    tree->GetUniquePermutations(tree_leaves, permutations);

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
    if(fit_vertex) {
        variable_names.emplace_back(Vertex->Name);
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
            return (IM_calc - IM_expected)/IM_Sigma;
        });
    });

    LOG(INFO) << "Have " << node_constraints.size() << " constraints at " << sum_daughters.size() << " nodes";


    // define the constraint
    // use local copies instead of this capture. capturing "this" is dangerous since
    // fitter might be moved around in memory...
    auto tree_leaves_copy   = tree_leaves;
    auto sum_daughters_copy = sum_daughters;
    auto IM_at_nodes = [tree_leaves_copy, sum_daughters_copy,
                       fit_vertex, node_constraints] (const vector<vector<double>>& v) {
        const auto  k = tree_leaves_copy.size();
        // k serves as an offset here
        const auto& vertex = fit_vertex ? v[k+0] : vector<double>{0.0, 0.0, 0.0};

        // assign values v to leaves' LVSum
        for(unsigned i=0;i<k;i++) {
            node_t& node = tree_leaves_copy[i]->Get();
            node.LVSum = node.Leave->GetLorentzVec(v[i], vertex);
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
