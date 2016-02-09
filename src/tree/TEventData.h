#pragma once

#include "base/printable.h"

#include "TID.h"
#include "TDetectorReadHit.h"
#include "TSlowControl.h"
#include "TUnpackerMessage.h"

#include "TTaggerHit.h"
#include "TTrigger.h"
#include "TTarget.h"

#include "TCluster.h"
#include "TCandidate.h"
#include "TParticle.h"

namespace ant {

struct TEventData;
using TEventDataPtr = std::unique_ptr<TEventData> ;

struct TEventData : printable_traits
{
    struct PTypeList {

        void Add(TParticlePtr&& particle) {
            lists[std::addressof(particle->Type())].emplace_back(particle);
            all.emplace_back(particle);
        }

        void Add(TParticlePtr& particle) {
            lists[std::addressof(particle->Type())].emplace_back(particle);
            all.emplace_back(particle);
        }

        TParticleList GetAll() const { return all; }

        TParticleList Get(const ant::ParticleTypeDatabase::Type& type) const {
            auto entry = lists.find(std::addressof(type));
            if(entry == lists.end()) {
                return {};
            }
            return entry->second;
        }

        template<class Archive>
        void save(Archive& archive) const {
            archive(all);
        }

        template<class Archive>
        void load(Archive& archive) {
            archive(all);
            for(const auto& p : all)
                lists[std::addressof(p->Type())].emplace_back(p);
        }

    private:
        TParticleList all;
        std::map<const ParticleTypeDatabase::Type*, TParticleList> lists;

    }; // PTypeList

    TEventData(const TID& id);
    TEventData();
    virtual ~TEventData();

    TID ID;
    std::vector<TDetectorReadHit> DetectorReadHits;
    std::vector<TSlowControl>     SlowControls;
    std::vector<TUnpackerMessage> UnpackerMessages;

    std::vector<TTaggerHit>  TaggerHits;
    TTrigger Trigger;
    TTarget  Target;

    TClusterList     Clusters;
    TCandidateList   Candidates;
    PTypeList        Particles;    // MCTrue final state, or identified from reconstructed candidates
    TParticleTree_t  ParticleTree;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(ID,
                DetectorReadHits, SlowControls, UnpackerMessages,
                TaggerHits, Trigger, Target,
                Clusters, Candidates, Particles, ParticleTree);
    }

    virtual std::ostream& Print(std::ostream& s) const override;

};

}