#include "Fitter.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"

#include "APLCON.hpp" // external project

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

const APLCON::Fit_Settings_t Fitter::Fitter::DefaultSettings = Fitter::MakeDefaultSettings();

APLCON::Fit_Settings_t Fitter::MakeDefaultSettings()
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = 30;
    return settings;
}

TParticlePtr Fitter::FitParticle::AsFitted() const
{
    if(!isfinite(Fitted_Z_Vertex))
        throw Exception("Need z vertex to calculate LorentzVec");

    auto p = make_shared<TParticle>(Particle->Type(), GetLorentzVec(Fitted_Z_Vertex));
    p->Candidate = Particle->Candidate; // link Candidate
    return p;
}

std::vector<double> Fitter::FitParticle::GetValues_before() const
{
    std::vector<double> values(Vars.size());
    transform(Vars_before.begin(), Vars_before.end(), values.begin(),
              [] (const V_S_t& v) { return v.Value; });
    return values;
}

std::vector<double> Fitter::FitParticle::GetSigmas_before() const
{
    vector<double> sigmas(Vars.size());
    transform(Vars_before.begin(), Vars_before.end(), sigmas.begin(),
              [] (const V_S_t& v) { return v.Sigma; });
    return sigmas;
}

std::vector<double> Fitter::FitParticle::GetPulls() const
{
    return std::vector<double>(Pulls.begin(), Pulls.end());
}

void Fitter::FitParticle::Set(const TParticlePtr& p, const UncertaintyModel& model)
{
    Particle = p;
    const auto& sigmas = model.GetSigmas(*p);

    if(!p->Candidate)
        throw Exception("Need particle with candidate for fitting");

    const double invEk = 1.0/p->Ek();
    const double sigma_invEk = sigmas.sigmaEk*std_ext::sqr(invEk);

    Vars[0].SetValueSigma(invEk,  sigma_invEk);
    Vars[2].SetValueSigma(p->Phi(), sigmas.sigmaPhi);
    ShowerDepth = sigmas.ShowerDepth;

    // the parametrization, and thus the meaning of the linked fitter variables,
    // depends on the calorimeter
    auto& detector = Particle->Candidate->Detector;
    if(detector & Detector_t::Type_t::CB)
    {
        static auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
        Vars[1].SetValueSigma(p->Theta(), sigmas.sigmaTheta);
        const auto& CB_R = cb->GetInnerRadius() + ShowerDepth;
        Vars[3].SetValueSigma(CB_R, sigmas.sigmaCB_R);
    }
    else if(detector & Detector_t::Type_t::TAPS)
    {
        static auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
        const auto& TAPS_Lz = taps->GetZPosition() + ShowerDepth*std::cos(p->Theta());
        auto& pos = p->Candidate->FindCaloCluster()->Position;
        const vec3 TAPS_L{pos.x, pos.y, TAPS_Lz};
        const auto& TAPS_Rxy = std::sin(p->Theta())*TAPS_L.R();

        Vars[1].SetValueSigma(TAPS_Rxy,   sigmas.sigmaTAPS_Rxy);
        Vars[3].SetValueSigma(TAPS_L.R(), sigmas.sigmaTAPS_L);
    }
    else {
        throw Exception("Unknown/none detector type provided from uncertainty model");
    }

    // remember old vars (copies also Pulls, which is a bit more than needed)
    Vars_before = Vars;

    // clear the vertex, during fitting, it's explicitly passed to GetLorentzVec
    Fitted_Z_Vertex = std_ext::NaN;
}

LorentzVec Fitter::FitParticle::GetLorentzVec(double z_vertex) const noexcept
{
    const radian_t& phi   = Vars[2];

    // start at the lab origin in the frame of the vertex
    vec3 x{0,0,-z_vertex};

    auto& detector = Particle->Candidate->Detector;
    if(detector & Detector_t::Type_t::CB)
    {
        // for CB, parametrization is (Ek, theta, phi, CB_R)
        const radian_t& theta = Vars[1];
        const auto&     CB_R  = Vars[3];
        x += vec3::RThetaPhi(CB_R, theta, phi);
    }
    else if(detector & Detector_t::Type_t::TAPS)
    {
        // for TAPS, parametrization is (Ek, TAPS_Rxy, phi, TAPS_L)
        const auto& TAPS_Rxy = Vars[1];
        const auto& TAPS_L   = Vars[3];
        const auto& TAPS_L_z = sqrt(std_ext::sqr<double>(TAPS_L) - std_ext::sqr<double>(TAPS_Rxy));
        x += vec3(vec2::RPhi(TAPS_Rxy, phi), TAPS_L_z);
    }

    const mev_t& Ek = 1.0/Vars[0];
    const mev_t& E = Ek + Particle->Type().Mass();
    using std_ext::sqr;
    const mev_t& p = sqrt( sqr(E) - sqr(Particle->Type().Mass()) );

    const vec3& p_vec = x*p/x.R();
    return {p_vec, E};
}







