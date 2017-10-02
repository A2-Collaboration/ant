#include "KinFitter.h"

#include "base/Logger.h"
#include "base/std_ext/map.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

KinFitter::KinFitter(UncertaintyModelPtr uncertainty_model,
                     bool fit_Z_vertex,
                     const APLCON::Fit_Settings_t& settings) :
    Z_Vertex(fit_Z_vertex),
    aplcon(settings),
    Model(uncertainty_model) // may be nullptr
{

}


void KinFitter::SetZVertexSigma(double sigma)
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    Z_Vertex.Sigma = sigma;
    Z_Vertex.Sigma_before = sigma;
}

bool KinFitter::IsZVertexFitEnabled() const noexcept
{
    return Z_Vertex.IsEnabled;
}

bool KinFitter::IsZVertexUnmeasured() const
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    if(!std::isfinite(Z_Vertex.Sigma_before))
        throw Exception("Z Vertex sigma not set");
    return Z_Vertex.Sigma_before == 0.;
}

void KinFitter::SetTarget(double length, double center)
{
    if(!Z_Vertex.IsEnabled)
        throw Exception("Z Vertex fitting not enabled");
    Target.length = length;
    Target.center = center;
}

TParticlePtr KinFitter::GetFittedProton() const
{
    return Proton.AsFitted();
}

TParticleList KinFitter::GetFittedPhotons() const
{
    TParticleList photons;
    for(unsigned i=0;i<Photons.size();i++) {
        photons.emplace_back(Photons[i].AsFitted());
    }
    return photons;
}

double KinFitter::GetFittedBeamE() const
{
    return BeamE.Value;
}

TParticlePtr KinFitter::GetFittedBeamParticle() const
{
    return std::make_shared<TParticle>(ParticleTypeDatabase::BeamProton, BeamE.GetLorentzVec());
}

double KinFitter::GetFittedZVertex() const
{
    return Z_Vertex.Value;
}

double KinFitter::GetBeamEPull() const
{
    return BeamE.Pull;
}

double KinFitter::GetZVertexPull() const
{
    return Z_Vertex.Pull;
}

std::vector<Fitter::FitParticle> KinFitter::GetFitParticles() const
{
    std::vector<Fitter::FitParticle> particles{Proton};
    for(auto& photon : Photons)
        particles.emplace_back(photon);
    return particles;
}

APLCON::Result_t KinFitter::DoFit(double ebeam, const TParticlePtr& proton, const TParticleList& photons)
{
    PrepareFit(ebeam, proton, photons);

    const auto& r = aplcon.DoFit(BeamE, Proton, Photons, Z_Vertex, constraintEnergyMomentum);

    // tell the particles the z-vertex after fit
    if (r.Status==APLCON::Result_Status_t::Success && !isfinite( Z_Vertex.Value))
        throw Exception("Fitted Z-vertex not finite!");
    Proton.SetFittedZVertex(Z_Vertex.Value);
    for(auto& photon : Photons)
        photon.SetFittedZVertex(Z_Vertex.Value);

    return r;
}

std::array<double, 4> KinFitter::constraintEnergyMomentum(
        const KinFitter::BeamE_t& beam,
        const KinFitter::Proton_t& proton,
        const KinFitter::Photons_t& photons,
        const Fitter::Z_Vertex_t& z_vertex)
{
    // start with the incoming particle minus outgoing proton
    auto diff = beam.GetLorentzVec() - proton.GetLorentzVec(z_vertex.Value);

    // subtract outgoing photons
    for(const auto& photon : photons)
        diff -= photon.GetLorentzVec(z_vertex.Value);

    return {diff.E, diff.p.x, diff.p.y, diff.p.z};
}

void KinFitter::PrepareFit(double ebeam, const TParticlePtr& proton, const TParticleList& photons)
{
    if(!Model) {
        throw Exception("No uncertainty provided in ctor or set with SetUncertaintyModel");
    }

    BeamE.SetValueSigma(ebeam, Model->GetBeamEnergySigma(ebeam));
    Proton.Set(proton, *Model);

    Photons.resize(photons.size());
    LorentzVec photon_sum; // for proton's missing_E calculation later
    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        Photons[i].Set(photons[i], *Model);
        photon_sum += *photons[i];
    }

    if(Z_Vertex.IsEnabled) {
        if(!std::isfinite(Z_Vertex.Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Z_Vertex.Sigma = Z_Vertex.Sigma_before;
        Z_Vertex.Value = 0;

        // if target length was set, calculate starting point for z vertex if parameter is unmeasured
        if(std::isfinite(Target.length)) {
            if(IsZVertexUnmeasured())
                Z_Vertex.Value = CalcZVertexStartingPoint();
        }
    }

    // only set Proton Ek to missing energy if unmeasured
    if(Proton.IsEkUnmeasured()) {
        const LorentzVec missing = BeamE.GetLorentzVec() - photon_sum;
        const double M = Proton.Particle->Type().Mass();
        using std_ext::sqr;
        const double missing_Ek = sqrt(sqr(missing.P()) + sqr(M)) - M;
        Proton.SetEk(missing_Ek);
    }

}

double KinFitter::CalcZVertexStartingPoint() const
{
    const double step_width = .5;
    using pos_map = std::map<double, double>;
    pos_map pos;

    // minimize energy constraint to obtain best starting point
    const auto energy_constraint = [=] (const double z_vertex) -> double {
        // start with the incoming particle minus outgoing proton
        auto diff = BeamE.GetLorentzVec() - Proton.GetLorentzVec(z_vertex);

        // subtract outgoing photons
        for(const auto& photon : Photons)
            diff -= photon.GetLorentzVec(z_vertex);

        // return absolute value to simplify finding the minimum
        return std::abs(diff.E);
    };

    for(auto target_pos = Target.start(); target_pos <= Target.end(); target_pos += step_width)
        pos.emplace(std::make_pair(target_pos, energy_constraint(target_pos)));

    auto it = std_ext::min_map_element(pos);

    return it->first;
}



LorentzVec KinFitter::BeamE_t::GetLorentzVec() const noexcept
{
    // Beam Lorentz vector:
    // beam    LorentzVec(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
    // target  LorentzVec(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
    const LorentzVec beam({0, 0, Value}, Value);
    /// \todo Target is always assumed proton...
    const LorentzVec target({0,0,0}, ParticleTypeDatabase::Proton.Mass());

    return target + beam;
}
