#pragma once

#include "base/printable.h"
#include "base/types.h"
#include "analysis/ParticleType.h"
#include "analysis/data/Candidate.h"

#include "TLorentzVector.h"
#include <vector>
#include <list>
#include <iostream>


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
    CandidateList candidates;

public:

    Particle(const ant::ParticleTypeDatabase::Type& _type, ant::mev_t _Ek, ant::radian_t _theta, ant::radian_t _phi);

    Particle(const ParticleTypeDatabase::Type &_type, const TLorentzVector &_lorentzvector):
        TLorentzVector(_lorentzvector),
        type(&_type),
        parents(),
        daughters(),
        candidates()
    {}

    Particle(const ParticleTypeDatabase::Type& _type, ant::CandidatePtr candidate):
        Particle(_type,candidate->ClusterEnergy(),candidate->Theta(), candidate->Phi()) {
        candidates.emplace_back(candidate);
    }

    Particle(const Particle&) = delete;
    Particle& operator= (const Particle&) = delete;
    virtual ~Particle() = default;

    mev_t Ek() const { return E() - type->Mass(); }


    const ParticleTypeDatabase::Type& Type() const { return *type; }
    void ChangeType(const ParticleTypeDatabase::Type& newtype );

    void SetLorentzVector( const TLorentzVector& lv ) { *((TLorentzVector*)this) = lv; }


    bool hasParent() const { return daughters.size() != 0; }
    ParticlePtr Partent() const { return parents.front(); }
    ParticleList Parents() { return parents; }
    const ParticleList Parents() const { return parents; }

    bool hasDaughters() const { return daughters.size() != 0; }
    ParticleList& Daughters() { return daughters; }
    const ParticleList& Daughters() const { return daughters; }

    bool hasCandidates() const { return candidates.size() != 0; }
    CandidateList& Candidates() { return candidates; }
    const CandidateList& Candidates() const { return candidates; }

    virtual std::ostream& Print(std::ostream& stream) const;

    using TLorentzVector::Angle;

    static double calcAngle( const ParticlePtr& p1, const ParticlePtr& p2 ) {
        return p1->Angle(p2->Vect());
    }

    static void RecPrint(const ParticlePtr& p, std::ostream& stream=std::cout);

};

}
