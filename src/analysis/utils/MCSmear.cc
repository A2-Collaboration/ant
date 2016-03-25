#include "MCSmear.h"
#include <memory>
#include "tree/TParticle.h"
#include "analysis/utils/Fitter.h"
#include "base/std_ext/memory.h"

#include "TRandom2.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

ant::TParticlePtr MCSmear::Smear(const TParticlePtr& p) const
{
    const auto& type = p->Type();

    const auto sigma = model->GetSigmas(*p);

    const double E     = rng->Gaus(p->Ek(),    sigma.sigmaE);
    const double Theta = rng->Gaus(p->Theta(), sigma.sigmaTheta);
    const double Phi   = rng->Gaus(p->Phi(),   sigma.sigmaPhi);

    auto sp = make_shared<TParticle>(type, E, Theta, Phi);
    sp->Candidate = p->Candidate;

    return sp;
}

MCSmear::MCSmear(const std::shared_ptr<const utils::Fitter::UncertaintyModel>& m):
    model(m), rng(std_ext::make_unique<TRandom2>()) {}

MCSmear::~MCSmear() {}
