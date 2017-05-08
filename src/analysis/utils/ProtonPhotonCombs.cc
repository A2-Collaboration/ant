#include "ProtonPhotonCombs.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

ProtonPhotonCombs::Combinations_t&
ProtonPhotonCombs::Combinations_t::Observe(const Observer_t& observer, const string& prefix) noexcept
{
    Observer = observer;
    ObserverPrefix = prefix;
    if(Observer && !prefix.empty()) {
        // maybe a bit tedious, but that's how it's counted...
        for(auto i=0u;i<this->size();i++)
            Observer(prefix);
    }
    return *this;
}

ProtonPhotonCombs::Combinations_t&
ProtonPhotonCombs::Combinations_t::FilterMult(unsigned nPhotonsRequired, double maxDiscardedEk) noexcept
{
    auto it = this->begin();
    while(it != this->end()) {
        const auto nPhotons = it->Photons.size();
        if(nPhotons < nPhotonsRequired) {
            it = this->erase(it);
            continue;
        }
        // calc discarded Ek and do cut
        it->DiscardedEk = 0;
        for(auto i=nPhotonsRequired;i<nPhotons;i++) {
            it->DiscardedEk += it->Photons[i]->Ek();
        }
        if(it->DiscardedEk >= maxDiscardedEk) {
            it = this->erase(it);
            continue;
        }
        if(Observer && isfinite(maxDiscardedEk)) {
            Observer(std_ext::formatter() << ObserverPrefix << "DiscEk<" << maxDiscardedEk);
        }
        it->Photons.resize(nPhotonsRequired); // will always shrink, as nPhotons >= nPhotonsRequired
        ++it;
    }
    return *this;
}

ProtonPhotonCombs::Combinations_t&
ProtonPhotonCombs::Combinations_t::FilterIM(const IntervalD& photon_IM_sum_cut) noexcept
{
    auto it = this->begin();
    while(it != this->end()) {
        it->PhotonSum = LorentzVec{{0,0,0}, 0};
        for(const auto& p : it->Photons)
            it->PhotonSum += *p;
        if(!photon_IM_sum_cut.Contains(it->PhotonSum.M())) {
            it = this->erase(it);
        }
        else {
            if(Observer && photon_IM_sum_cut != nocut)
                Observer(ObserverPrefix+photon_IM_sum_cut.AsRangeString("IM(#gamma)"));
            ++it;
        }
    }
    called_FilterIM = true;
    return *this;
}

ProtonPhotonCombs::Combinations_t&
ProtonPhotonCombs::Combinations_t::FilterMM(const TTaggerHit& taggerhit,
                                            const IntervalD& missingmass_cut,
                                            const ParticleTypeDatabase::Type& target) noexcept
{
    // automatically call FilterIM (without any cut) to ensure PhotonSum is calculated
    if(!called_FilterIM)
        FilterIM();

    auto it = this->begin();
    const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(target.Mass());
    while(it != this->end()) {
        // remember hit and cut on missing mass
        it->MissingMass = (beam_target - it->PhotonSum).M();
        if(!missingmass_cut.Contains(it->MissingMass))
            it = this->erase(it);
        else {
            if(Observer && missingmass_cut != nocut)
                // note that in A2's speech is often "missing mass of proton",
                // but it's actually the "missing mass of photons" expected to be close to the
                // rest mass of the proton
                Observer(ObserverPrefix+missingmass_cut.AsRangeString("MM(#gamma)"));
            ++it;
        }
    }
    return *this;
}

ProtonPhotonCombs::Combinations_t&
ProtonPhotonCombs::Combinations_t::FilterCustom(const cut_t& cut, const string& name)
{
    auto it = this->begin();
    while(it != this->end()) {
        if(cut(*it)) {
            it = this->erase(it);
        }
        else {
            if(Observer && name != "")
                Observer(ObserverPrefix+name);
            ++it;
        }
    }
    return *this;
}

ProtonPhotonCombs::Combinations_t
ProtonPhotonCombs::MakeCombinations(const TCandidateList& cands, const combfilter_t& filter) noexcept
{
    TParticleList all_protons;
    TParticleList all_photons;
    for(auto cand : cands.get_iter()) {
        all_protons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, cand));
        all_photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, cand));
    }

    // important for DiscardedEk cut later
    std::sort(all_photons.begin(), all_photons.end(), [] (const TParticlePtr& a, const TParticlePtr& b) {
        // sort by descending kinetic energy
        return a->Ek() > b->Ek();
    });

    Combinations_t combs;
    for(const auto& proton : all_protons) {
        combs.emplace_back(proton);
        auto& comb = combs.back();
        for(auto photon : all_photons) {
            if(photon->Candidate == proton->Candidate)
                continue;
            comb.Photons.emplace_back(photon);
        }
        assert(comb.Photons.size()+1 == all_photons.size());
        // allow custom modification to combinations
        filter(comb);
    }
    // as every candidate is the proton, it's only
    return combs;
}
