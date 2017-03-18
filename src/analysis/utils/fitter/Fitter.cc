#include "Fitter.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"
#include "base/std_ext/memory.h"

#include "APLCON.hpp" // external project

using namespace std;
using namespace ant;
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

Fitter::FitParticle::FitParticle(const string& name,
                                 APLCON& aplcon,
                                 std::shared_ptr<Fitter::FitVariable> z_vertex) :
    Vars(4), // it's a lucky coincidence that particles in CB/TAPS are both parametrized by 4 values
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
    auto p = make_shared<TParticle>(Particle->Type(), GetLorentzVec());
    p->Candidate = Particle->Candidate; // link Candidate
    return p;
}

double Fitter::FitParticle::GetShowerDepth() const
{
    return ShowerDepth;
}

std::vector<double> Fitter::FitParticle::GetValues_before() const
{
    std::vector<double> values(Vars.size());
    transform(Vars.begin(), Vars.end(), values.begin(),
              [] (const FitVariable& v) { return v.Value_before; });
    return values;
}

std::vector<double> Fitter::FitParticle::GetSigmas_before() const
{
    vector<double> sigmas(Vars.size());
    transform(Vars.begin(), Vars.end(), sigmas.begin(),
              [] (const FitVariable& v) { return v.Sigma_before; });
    return sigmas;
}

std::vector<double> Fitter::FitParticle::GetPulls() const
{
    std::vector<double> pulls(Vars.size());
    transform(Vars.begin(), Vars.end(), pulls.begin(),
              [] (const FitVariable& v) { return v.Pull; });
    return pulls;
}

void Fitter::FitParticle::Set(const TParticlePtr& p,
                              const UncertaintyModel& uncertainty)
{
    Particle = p;
    const auto& sigmas = uncertainty.GetSigmas(*p);

    if(!p->Candidate)
        throw Exception("Need particle with candidate for fitting");

    Vars[0].SetValueSigma(p->Ek(),  sigmas.sigmaEk);
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
}

LorentzVec Fitter::FitParticle::GetLorentzVec() const
{
    const auto z_vertex = Z_Vertex ? Z_Vertex->Value : 0;

    vector<double> values(Vars.size());
    transform(Vars.begin(), Vars.end(), values.begin(),
                  [] (const FitVariable& v) { return v.Value; });

    return GetLorentzVec(values, z_vertex);
}

LorentzVec Fitter::FitParticle::GetLorentzVec(const std::vector<double>& values,
                                              const double z_vertex) const
{

    const radian_t& phi   = values[2];

    // start at the lab origin in the frame of the vertex
    vec3 x{0,0,-z_vertex};

    auto& detector = Particle->Candidate->Detector;
    if(detector & Detector_t::Type_t::CB)
    {
        // for CB, parametrization is (Ek, theta, phi, CB_R)
        const radian_t& theta = values[1];
        const auto&     CB_R  = values[3];
        x += vec3::RThetaPhi(CB_R, theta, phi);
    }
    else if(detector & Detector_t::Type_t::TAPS)
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

    const mev_t& Ek = values[0];
    const mev_t& E = Ek + Particle->Type().Mass();
    using std_ext::sqr;
    const mev_t& p = sqrt( sqr(E) - sqr(Particle->Type().Mass()) );

    const vec3& p_vec = x*p/x.R();
    return {p_vec, E};
}







