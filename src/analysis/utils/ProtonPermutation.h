#pragma once

#include "tree/TCandidate.h"
#include "tree/TParticle.h"

namespace ant {
namespace analysis {
namespace utils {

class ProtonPermutation {
protected:
    TCandidatePtrList cands;
    TCandidatePtrList::const_iterator p_it = cands.begin();

    const TCandidatePtr true_proton;

    TParticlePtr  proton;
    TParticleList photons;

    bool trueMatch = false;

    void Fill();

public:
    ProtonPermutation(const TCandidatePtrList& candidates, const TCandidatePtr& true_p=nullptr);

    const TParticlePtr   Proton()  const { return proton; }
    const TParticleList& Photons() const { return photons; }

    void Next();
    bool Good() const { return p_it != cands.cend(); }

    bool isTrueProton() const { return trueMatch; }
};

}
}
}
