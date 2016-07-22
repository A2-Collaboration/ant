#include "MCSmear.h"
#include <memory>
#include "tree/TParticle.h"
#include "analysis/utils/Fitter.h"
#include "base/std_ext/memory.h"

#include "TRandom2.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

MCSmear::MCSmear(utils::UncertaintyModelPtr m):
    model(m), rng(std_ext::make_unique<TRandom2>()) {}

MCSmear::~MCSmear() {}

ant::TParticlePtr MCSmear::Smear(const TParticlePtr& p, Uncertainties_t& sigmas) const
{
    const auto& type = p->Type();

    const bool isBeam = type == ParticleTypeDatabase::BeamTarget;

    sigmas = isBeam ? Uncertainties_t{ model->GetBeamEnergySigma(p->Ek()), std_ext::NaN, std_ext::NaN}
                    : model->GetSigmas(*p);

    const double Ek    = rng->Gaus(p->Ek(), sigmas.sigmaE);
    const double Theta = isBeam ? p->Theta() : rng->Gaus(p->Theta(), sigmas.sigmaTheta);
    const double Phi   = isBeam ? p->Phi()   : rng->Gaus(p->Phi(),   sigmas.sigmaPhi);

    auto sp = make_shared<TParticle>(type, Ek, Theta, Phi);
    sp->Candidate = p->Candidate;
    return sp;
}

TParticlePtr MCSmear::Smear(const TParticlePtr& p) const
{
    Uncertainties_t sigmas;
    return Smear(p, sigmas);
}
