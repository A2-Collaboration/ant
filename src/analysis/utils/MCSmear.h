#pragma once
#include <memory>
#include "analysis/utils/Fitter.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"

class TRandom;

namespace ant {

namespace analysis {
namespace utils {
class MCSmear {
protected:
    std::shared_ptr<const utils::Fitter::UncertaintyModel> model;
    std::unique_ptr<TRandom> rng;

public:

    MCSmear(const std::shared_ptr<const utils::Fitter::UncertaintyModel>& m);

    ~MCSmear();

    ant::TParticlePtr  Smear(const ant::TParticlePtr& p) const;

};

}
}
}
