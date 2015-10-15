#pragma once

#include "analysis/data/Candidate.h"

#include "base/printable.h"
#include "base/types.h"
#include "base/ParticleType.h"
#include "base/Tree.h"

#include "TLorentzVector.h"

#include <vector>
#include <list>
#include <iostream>
#include <memory>

namespace ant {
namespace analysis {
namespace data {

class Particle;

using ParticlePtr    = std::shared_ptr<Particle>;
using ParticleList   = std::vector<ParticlePtr>;
using ParticleTree_t = std::shared_ptr<Tree<ParticlePtr>>;


/**
 * @brief Base particle class
 */
class Particle: public TLorentzVector, public printable_traits {
protected:
    const ant::ParticleTypeDatabase::Type* type;
    CandidatePtr candidate;

public:

    Particle(const ant::ParticleTypeDatabase::Type& _type,
             ant::mev_t _Ek, ant::radian_t _theta, ant::radian_t _phi);

    Particle(const ParticleTypeDatabase::Type &_type, const TLorentzVector &_lorentzvector):
        TLorentzVector(_lorentzvector),
        type(&_type)
    {}

    Particle(const ParticleTypeDatabase::Type& _type, const CandidatePtr& candidate_):
        Particle(_type,candidate_->ClusterEnergy(),candidate_->Theta(), candidate_->Phi())
    {
        candidate = candidate_;
    }

    Particle(const Particle&) = delete;
    Particle& operator= (const Particle&) = delete;
    virtual ~Particle() = default;

    mev_t Ek() const { return E() - type->Mass(); }


    const ParticleTypeDatabase::Type& Type() const { return *type; }
    void ChangeType(const ParticleTypeDatabase::Type& newtype );

    void SetLorentzVector( const TLorentzVector& lv ) { *((TLorentzVector*)this) = lv; }

    bool hasCandidate() const noexcept { return !candidate; }
    CandidatePtr& Candidate() noexcept { return candidate; }
    const CandidatePtr& Candidate() const noexcept { return candidate; }

    virtual std::ostream& Print(std::ostream& stream) const;

    using TLorentzVector::Angle;

    static double calcAngle( const ParticlePtr& p1, const ParticlePtr& p2 ) {
        return p1->Angle(p2->Vect());
    }
};

}
}
}
