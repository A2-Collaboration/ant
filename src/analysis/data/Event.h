#pragma once

#include "Particle.h"
#include "TaggerHit.h"
#include "Trigger.h"
#include "Target.h"

#include "base/types.h"
#include "base/printable.h"
#include "base/Tree.h"

#include <vector>
#include <memory>
#include <map>

namespace ant {
namespace analysis {
namespace data {

using ClusterList = std::vector<TCluster>;

struct Event : printable_traits {

    struct Data : printable_traits {

        struct PTypeLists : printable_traits {

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

            virtual ~PTypeLists() = default;
            std::ostream& Print(std::ostream& stream) const;

        protected:
            ParticleList particles;
            std::map<const ant::ParticleTypeDatabase::Type*, ParticleList> lists;

            static const ParticleList empty;
        }; // PTypeLists


        PTypeLists     Particles;      // final state / reconstructred particles
        PTypeLists     Intermediates;  // intermediate particles (if any)
        ParticleTree_t ParticleTree;   // particle tree (if available)

        CandidateList  Candidates;     // particle candidates (if any)
        TaggerHitList  TaggerHits;     // tagger hits
        Trigger_t      Trigger;
        Target_t       Target;
        ClusterList    AllClusters;

        virtual ~Data() = default;
        std::ostream& Print(std::ostream& stream) const;
    }; // Data

    Data    Reconstructed;
    Data    MCTrue;

    virtual ~Event() = default;
    std::ostream& Print(std::ostream& stream) const;
};

}}} // namespace ant::analysis::data
