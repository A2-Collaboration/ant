#pragma once
#include <memory>
#include "analysis/utils/Uncertainties.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"

class TRandom;

namespace ant {

namespace analysis {
namespace utils {
class MCSmear {
protected:
    UncertaintyModelPtr model;
    std::unique_ptr<TRandom> rng;

public:

    MCSmear(UncertaintyModelPtr m);

    ~MCSmear();

    ant::TParticlePtr  Smear(const ant::TParticlePtr& p) const;

    ant::TParticlePtr  Smear(const ant::TParticlePtr& p, Uncertainties_t& sigmas) const;


};

}
}
}
