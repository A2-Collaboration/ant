#pragma once


#include "Particle.h"
#include "Candidate.h"
#include "TaggerHit.h"
#include "TriggerInfo.h"

#include "base/types.h"
#include "base/printable.h"
#include "base/Tree.h"

#include <vector>
#include <memory>
#include <map>

namespace ant {
namespace analysis {
namespace data {


class Event: public ant::printable_traits {
public:
    class Data: public ant::printable_traits {
    public:

        class PTypeLists: public ant::printable_traits {
        protected:
            ParticleList particles;
            std::map<const ant::ParticleTypeDatabase::Type*, ParticleList> lists;

            static const ParticleList empty;

        public:
            void AddParticle(ParticlePtr&& particle) {
                lists[&particle->Type()].emplace_back(particle);
                particles.emplace_back(particle);
            }

            void AddParticle(ParticlePtr& particle) {
                lists[&particle->Type()].emplace_back(particle);
                particles.emplace_back(particle);
            }


            const ParticleList& GetAll() const { return particles; }

            const ParticleList& Get(const ant::ParticleTypeDatabase::Type& type) const {
                auto entry = lists.find(&type);
                if(entry == lists.end()) {
                    return empty;
                }
                return entry->second;
            }

            std::ostream& Print(std::ostream& stream) const;
        };

    protected:

        PTypeLists   particles;      // final state / reconstructred particles
        PTypeLists   intermediates;  // intermediate particles (if any)
        ParticleTree_t particletree; // particle tree (if available)

        CandidateList  candidates;     // particle candidates (if any)
        TaggerHitList  taggerhits;     // tagger hits
        TriggerInfo    triggerinfo;
        ClusterList    allclusters;


    public:

        const PTypeLists& Particles() const { return particles; }
              PTypeLists& Particles()       { return particles; }

        const PTypeLists& Intermediates() const { return intermediates; }
              PTypeLists& Intermediates()       { return intermediates; }

        const ParticleTree_t& ParticleTree() const { return particletree; }
              ParticleTree_t& ParticleTree()       { return particletree; }


        const CandidateList& Candidates() const { return candidates; }
              CandidateList& Candidates()       { return candidates; }

        const TaggerHitList& TaggerHits() const { return taggerhits; }
              TaggerHitList& TaggerHits()       { return taggerhits; }

        const TriggerInfo& TriggerInfos() const { return triggerinfo; }
              TriggerInfo& TriggerInfos()       { return triggerinfo; }

        const ClusterList& AllClusters() const { return allclusters; }
              ClusterList& AllClusters()       { return allclusters; }


              std::ostream& Print(std::ostream& stream) const;
    }; // class Data

protected:
    Data    reconstructed;
    Data    mctrue;
public:
    const Data& Reconstructed() const { return reconstructed; }
          Data& Reconstructed()       { return reconstructed; }

    const Data& MCTrue() const { return mctrue; }
          Data& MCTrue()       { return mctrue; }

    std::ostream& Print(std::ostream& stream) const;
};
}
}
}
