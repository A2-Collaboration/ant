#ifndef EVENT_H
#define EVENT_H

#include "base/types.h"

#include "Particle.h"
#include "Track.h"
#include "TaggerHit.h"
#include "TriggerInfo.h"
#include "base/printable.h"

#include <vector>
#include <memory>
#include <map>

namespace ant {


class Event: public ant::printable_traits {
public:
    class Data: public ant::printable_traits {
    public:

        class PTypeLists: public ant::printable_traits {
        protected:
            ant::ParticleList particles;
            std::map<const ParticleTypeDatabase::Type*, ant::ParticleList> lists;

            static const ant::ParticleList empty;

        public:
            void AddParticle(ParticlePtr&& particle) {
                lists[&particle->Type()].emplace_back(particle);
                particles.emplace_back(particle);
            }


            const ParticleList& GetAll() const { return particles; }

            const ParticleList& Get(const ParticleTypeDatabase::Type& type) const {
                auto entry = lists.find(&type);
                if(entry == lists.end()) {
                    return empty;
                }
                return entry->second;
            }

            std::ostream &Print(std::ostream &stream) const;
        };


        PTypeLists   particles;      //final state / reconstructred particles
        PTypeLists   intermediates;  //intermediate particles (if any)

        ant::TrackList      tracks;         //detector tracks (if any)
        ant::TaggerHistList taggerhits;     //tagger hits
        ant::TriggerInfo    triggerinfo;

    public:

        const PTypeLists& Particles() const { return particles; }
              PTypeLists& Particles()       { return particles; }

        const PTypeLists& Intermediates() const { return intermediates; }
              PTypeLists& Intermediates()       { return intermediates; }

        const TrackList& Tracks() const { return tracks; }
              TrackList& Tracks()       { return tracks; }

        const TaggerHistList& TaggerHits() const { return taggerhits; }
              TaggerHistList& TaggerHits()       { return taggerhits; }

        const TriggerInfo& TriggerInfos() const { return triggerinfo; }
              TriggerInfo& TriggerInfos()       { return triggerinfo; }

              std::ostream &Print(std::ostream &stream) const;
    };

protected:
    Data    reconstructed;
    Data    mctrue;
public:
    const Data& Reconstructed() const { return reconstructed; }
          Data& Reconstructed()       { return reconstructed; }

    const Data& MCTrue() const { return mctrue; }
          Data& MCTrue()       { return mctrue; }

    std::ostream &Print(std::ostream &stream) const;
};

}

#endif
