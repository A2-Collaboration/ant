#include "NoProtonFitter.h"

#include "base/Logger.h"
#include "base/std_ext/map.h"
#include "TVector3.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

NoProtonFitter::NoProtonFitter(UncertaintyModelPtr uncertainty_model,
                     bool fit_Z_vertex,
                     const APLCON::Fit_Settings_t& settings) :
    Z_Vertex(fit_Z_vertex),
    aplcon(settings),
    Model(uncertainty_model) // may be nullptr
{

}


void NoProtonFitter::SetZVertexSigma(double sigma)
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    Z_Vertex.Sigma = sigma;
    Z_Vertex.Sigma_before = sigma;
}

bool NoProtonFitter::IsZVertexFitEnabled() const noexcept
{
    return Z_Vertex.IsEnabled;
}

bool NoProtonFitter::IsZVertexUnmeasured() const
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    if(!std::isfinite(Z_Vertex.Sigma_before))
        throw Exception("Z Vertex sigma not set");
    return Z_Vertex.Sigma_before == 0.;
}

void NoProtonFitter::SetTarget(double length, double center)
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    Target.length = length;
    Target.center = center;
}

LorentzVec NoProtonFitter::GetFittedProton() const
{
    return Proton.AsFitted();
}

TParticleList NoProtonFitter::GetFittedPhotons() const
{
    TParticleList photons;
    for(unsigned i=0;i<Photons.size();i++) {
        photons.emplace_back(Photons[i].AsFitted());
    }
    return photons;
}

double NoProtonFitter::GetFittedBeamE() const
{
    return BeamE.Value;
}

TParticlePtr NoProtonFitter::GetFittedBeamParticle() const
{
    return std::make_shared<TParticle>(ParticleTypeDatabase::BeamProton, BeamE.GetLorentzVec());
}

double NoProtonFitter::GetFittedZVertex() const
{
    return Z_Vertex.Value;
}

double NoProtonFitter::GetBeamEPull() const
{
    return BeamE.Pull;
}

double NoProtonFitter::GetZVertexPull() const
{
    return Z_Vertex.Pull;
}

std::vector<Fitter::FitParticle> NoProtonFitter::GetFitParticles() const
{
    std::vector<Fitter::FitParticle> photons;
    for(auto& photon : Photons)
        photons.emplace_back(photon);
    return photons;
}

APLCON::Result_t NoProtonFitter::DoFit(double ebeam, const TParticleList& photons)
{
    PrepareFit(ebeam, photons);

    const auto& r = aplcon.DoFit(BeamE, Proton, Photons, Z_Vertex, constraintEnergyMomentum);

    // tell the particles the z-vertex after fit
    if (r.Status==APLCON::Result_Status_t::Success && !isfinite( Z_Vertex.Value))
        throw Exception("Fitted Z-vertex not finite!");
    for(auto& photon : Photons)
        photon.SetFittedZVertex(Z_Vertex.Value);

    return r;
}

std::array<double, 4> NoProtonFitter::constraintEnergyMomentum(const NoProtonFitter::BeamE_t& beam,
        const NoProtonFitter::Proton_t& proton,
        const NoProtonFitter::Photons_t& photons,
        const Fitter::Z_Vertex_t& z_vertex)
{
    // start with the incoming particle minus outgoing proton
    auto diff = beam.GetLorentzVec() - proton.GetLorentzVec();

    // subtract outgoing photons
    for(const auto& photon : photons)
        diff -= photon.GetLorentzVec(z_vertex.Value);

    return {diff.E, diff.p.x, diff.p.y, diff.p.z};
}

void NoProtonFitter::PrepareFit(double ebeam, const TParticleList& photons)
{
    if(!Model) {
        throw Exception("No uncertainty provided in ctor or set with SetUncertaintyModel");
    }

    BeamE.SetValueSigma(ebeam, Model->GetBeamEnergySigma(ebeam));

    Photons.resize(photons.size());
    LorentzVec photon_sum; // for proton's missing_E calculation later
    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        Photons[i].Set(photons[i], *Model);
        photon_sum += *photons[i];
    }

    const LorentzVec missing = BeamE.GetLorentzVec() - photon_sum;
    Proton.Set(missing);

    if(IsZVertexUnmeasured())
        throw Exception("Z Vertex sigma can't be unmeasured when the proton is");
    if(Z_Vertex.IsEnabled) {
        if(!std::isfinite(Z_Vertex.Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Z_Vertex.Sigma = Z_Vertex.Sigma_before;
        Z_Vertex.Value = 0;
    }

}

LorentzVec NoProtonFitter::BeamE_t::GetLorentzVec() const noexcept
{
    // Beam Lorentz vector:
    // beam    LorentzVec(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
    // target  LorentzVec(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
    const LorentzVec beam({0, 0, Value}, Value);
    /// \todo Target is always assumed proton...
    const LorentzVec target({0,0,0}, ParticleTypeDatabase::Proton.Mass());

    return target + beam;
}

void NoProtonFitter::FitProtonUnMeas::Set(const LorentzVec MissingP)
{

    //double Ek = MissingP.E - ParticleTypeDatabase::Proton.Mass();
    const double M = ParticleTypeDatabase::Proton.Mass();
    using std_ext::sqr;
    const double missing_Ek = sqrt(sqr(MissingP.P()) + sqr(M)) - M;
    double invEk = 1./missing_Ek;
    Vars[0].SetValueSigma(invEk, 0.);
    Vars[1].SetValueSigma(MissingP.Phi(), 0.);
    Vars[2].SetValueSigma(MissingP.Theta(), 0.);

    // remember old vars
    Vars_before = Vars;
}

LorentzVec NoProtonFitter::FitProtonUnMeas::GetLorentzVec() const noexcept
{
    const double E = 1.0/Vars[0] + ParticleTypeDatabase::Proton.Mass();
    const double p = sqrt(pow(E,2) - pow(ParticleTypeDatabase::Proton.Mass(),2));
    const double theta = Vars[1];
    const double phi = Vars[2];

    TVector3 tv;
    tv.SetMagThetaPhi(p,theta,phi);
    vec3 v(tv);

    return {v, E};
}

LorentzVec NoProtonFitter::FitProtonUnMeas::AsFitted() const
{
    return GetLorentzVec();
}
