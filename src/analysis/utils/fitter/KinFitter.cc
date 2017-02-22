#include "KinFitter.h"

#include "base/std_ext/memory.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

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

        return vector<double>({
                        diff.E,
                        diff.p.x,
                        diff.p.y,
                        diff.p.z,
                    });
    };

    aplcon->AddConstraint("E-p", variable_names, EnergyMomentum);

}

KinFitter::~KinFitter()
{

}

void KinFitter::SetEgammaBeam(const double ebeam)
{
    BeamE->SetValueSigma(ebeam, uncertainty->GetBeamEnergySigma(ebeam));
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
    return BeamE->Value;
}

TParticlePtr KinFitter::GetFittedBeamParticle() const
{
    return std::make_shared<TParticle>(ParticleTypeDatabase::BeamProton,
                                         MakeBeamLorentzVec(BeamE->Value));
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

std::vector<double> KinFitter::GetProtonPulls() const
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
    if(Z_Vertex) {
        if(!std::isfinite(Z_Vertex->Sigma_before))
            throw Exception("Z Vertex sigma not set although enabled");
        Z_Vertex->Value = 0;
        Z_Vertex->Sigma = Z_Vertex->Sigma_before;
    }


    // only set Proton Ek to missing energy if unmeasured
    auto& Var_Ek = Proton->Vars[0];
    if(Var_Ek.Sigma == 0) {
        LorentzVec missing = MakeBeamLorentzVec(BeamE->Value);
        for(const auto& photon : Photons)
            missing -= *photon->Particle;
        const double M = Proton->Particle->Type().Mass();
        using std_ext::sqr;
        const double missing_E = sqrt(sqr(missing.P()) + sqr(M)) - M;
        Var_Ek.SetValueSigma(missing_E, Var_Ek.Sigma);
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