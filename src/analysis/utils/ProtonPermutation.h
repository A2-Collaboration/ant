#pragma once

#include "tree/TCandidate.h"
#include "tree/TParticle.h"

namespace ant {
namespace analysis {
namespace utils {

class ProtonPermutation {
protected:
    const TCandidateList& cands;
    TCandidateList::const_iterator p_it = cands.begin();

    TParticlePtr  proton;
    TParticleList photons;

    void Fill();

public:
    ProtonPermutation(const TCandidateList& candidates);

    const TParticlePtr   Proton()  const { return proton; }
    const TParticleList& Photons() const { return photons; }

    bool Next();
};

}
}
}
