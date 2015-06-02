#ifndef ANT_PARTICLE_H
#define ANT_PARTICLE_H

#include "base/printable.h"
#include "base/types.h"
#include "ParticleType.h"
#include "Track.h"

#include "TLorentzVector.h"
#include <vector>
#include <list>
#include <ostream>


namespace ant {

class Particle;

using ParticlePtr  = std::shared_ptr<Particle>;
using ParticleList = std::vector<ParticlePtr>;

/**
 * @brief Base particle class
 */
class Particle: public TLorentzVector, public printable_traits {
protected:
    const ant::ParticleTypeDatabase::Type* type;
    ParticleList parents;
    ParticleList daughters;
    TrackList tracks;

public:

    Particle(const ant::ParticleTypeDatabase::Type& _type, ant::mev_t _Ek, ant::radian_t _theta, ant::radian_t _phi);

    Particle(const ParticleTypeDatabase::Type &_type, const TLorentzVector &_lorentzvector):
        TLorentzVector(_lorentzvector),
        type(&_type)
    {}

    Particle(const ParticleTypeDatabase::Type& _type, ant::TrackPtr track):
        Particle(_type,track->ClusterEnergy(),track->Theta(), track->Phi()) {
        tracks.emplace_back(track);
    }

    virtual ~Particle() {}

    mev_t Ek() const { return E() - type->Mass(); }


    const ParticleTypeDatabase::Type& Type() const { return *type; }
    void ChangeType(const ParticleTypeDatabase::Type& newtype );

    void SetLorentzVector( const TLorentzVector& lv ) { *((TLorentzVector*)this) = lv; }


    bool hasParent() const { return daughters.size() != 0; }
    ParticlePtr Partent() const { return parents.front(); }
    ParticleList Partents() { return parents; }
    const ParticleList Partents() const { return parents; }

    bool hasDaughters() const { return daughters.size() != 0; }
    ParticleList& Daughters() { return daughters; }
    const ParticleList& Daughters() const { return daughters; }

    bool hasTracks() const { return tracks.size() != 0; }
    TrackList& Tracks() { return tracks; }
    const TrackList& Tracks() const { return tracks; }

    virtual std::ostream& Print(std::ostream& stream) const;

    using TLorentzVector::Angle;

    static double calcAngle( const ParticlePtr& p1, const ParticlePtr& p2 ) {
        return p1->Angle(p2->Vect());
    }

};

}

#endif
