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

    sigmas = model->GetSigmas(*p);

    const double E     = rng->Gaus(p->Ek(),    sigmas.sigmaE);
    const double Theta = rng->Gaus(p->Theta(), sigmas.sigmaTheta);
    const double Phi   = rng->Gaus(p->Phi(),   sigmas.sigmaPhi);

    auto sp = make_shared<TParticle>(type, E, Theta, Phi);
    sp->Candidate = p->Candidate;

    return sp;
}

TParticlePtr MCSmear::Smear(const TParticlePtr& p) const
{
    Uncertainties_t sigmas;
    return Smear(p, sigmas);
}
