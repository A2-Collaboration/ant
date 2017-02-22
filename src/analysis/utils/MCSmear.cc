#include "MCSmear.h"
#include <memory>
#include "tree/TParticle.h"
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

    shared_ptr<TParticle> smeared;

    if(type == ParticleTypeDatabase::BeamTarget) {
        sigmas = { model->GetBeamEnergySigma(p->Ek()), std_ext::NaN, std_ext::NaN};

        // be careful about this composite particle
        // beamparticle = gamma + nucleon (at rest)

        const double Ek = rng->Gaus(p->Ek(), sigmas.sigmaEk); // photon energy

        smeared = make_shared<TParticle>(
                      type,
                      // spatial components are given by photon direction,
                      // time component is sum of resting target and photon energy
                      LorentzVec::EPThetaPhi(Ek + type.Mass(), Ek, p->Theta(), p->Phi())
                      );
    }
    else {
        sigmas = model->GetSigmas(*p);

        const double Ek    = rng->Gaus(p->Ek(),    sigmas.sigmaEk);
        const double Theta = rng->Gaus(p->Theta(), sigmas.sigmaTheta);
        const double Phi   = rng->Gaus(p->Phi(),   sigmas.sigmaPhi);

        smeared = make_shared<TParticle>(type, Ek, Theta, Phi);
        smeared->Candidate = p->Candidate;
    }


    return smeared;
}

TParticlePtr MCSmear::Smear(const TParticlePtr& p) const
{
    Uncertainties_t sigmas;
    return Smear(p, sigmas);
}
