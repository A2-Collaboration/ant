#include "ProtonPermutation.h"

#include "base/ParticleType.h"

#include <memory>

using namespace ant;
using namespace ant::analysis::utils;
using namespace std;

void ProtonPermutation::Fill()
{
    photons.clear();

    for(auto i = cands.cbegin(); i!=cands.cend(); ++i) {
        if(i != p_it) {
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, *i));
        } else {
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, *i);
            trueMatch = (*i == true_proton);
        }
    }
}

ProtonPermutation::ProtonPermutation(const TCandidateList& candidates, const TCandidatePtr& true_p):
    cands(candidates),
    true_proton(true_p)
{
    if(!cands.empty())
        photons.reserve(cands.size()-1);

    p_it = cands.begin();

    Fill();
}

void ProtonPermutation::Next()
{
    ++p_it;

    const auto end = (p_it == cands.cend());

    if(!end)
        Fill();
}
