#pragma once

#ifndef __CINT__
#include "TID.h"
#include "TDetectorReadHit.h"
#include "TSlowControl.h"
#include "TUnpackerMessage.h"

#include "TTagger.h"
#include "TTrigger.h"
#include "TTarget.h"

#include "TCluster.h"
#include "TCandidate.h"
#include "TParticle.h"

#include <iomanip>
#include <ctime>
#include <memory>
#endif

namespace ant {

#ifndef __CINT__
struct TEvent;
using TEventPtr = std::unique_ptr<TEvent>;
struct TEvent : printable_traits
#else
struct TEvent
#endif
{


#ifndef __CINT__

    struct Data : printable_traits
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

            const TParticleList& GetAll() const { return all; }

            const TParticleList& Get(const ant::ParticleTypeDatabase::Type& type) const {
                auto entry = lists.find(std::addressof(type));
                if(entry == lists.end()) {
                    return empty;
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
            static const TParticleList empty;
            TParticleList all;
            std::map<const ParticleTypeDatabase::Type*, TParticleList> lists;

        }; // PTypeList

        Data(const TID& id);
        Data();
        virtual ~Data();

        TID ID;
        std::vector<TDetectorReadHit> DetectorReadHits;
        std::vector<TSlowControl>     SlowControls;
        std::vector<TUnpackerMessage> UnpackerMessages;

        TTagger  Tagger;
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
                    Tagger, Trigger, Target,
                    Clusters, Candidates, Particles, ParticleTree);
        }

        virtual std::ostream& Print(std::ostream& s) const override;

    }; // Data

    using DataPtr = std::unique_ptr<Data> ;

    DataPtr Reconstructed;
    DataPtr MCTrue;

    template<class Archive>
    void serialize(Archive archive) {
        archive(Reconstructed, MCTrue);
    }

    virtual std::ostream& Print( std::ostream& s) const override;

    static std::unique_ptr<TEvent> MakeReconstructed(const TID& id);

    // TEvent is moveable
    TEvent(TEvent&&) = default;
    TEvent& operator=(TEvent&&) = default;

#endif

    TEvent();
    virtual ~TEvent();
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

private:
    // prevent ROOTcint from creating copy-constructors
    TEvent(const TEvent&);
    TEvent& operator=(const TEvent&);

};

}
