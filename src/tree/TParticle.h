#pragma once

#include "tree/TCandidate.h"

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

struct TParticle;

using TParticlePtr    = std::shared_ptr<TParticle>;
using TParticleTree_t = std::shared_ptr<Tree<TParticlePtr>>;


/**
 * @brief Base TParticle class
 */
struct TParticle : TLorentzVector, printable_traits {

    TCandidatePtr Candidate;
    const ParticleTypeDatabase::Type& Type() const { return *type; }

    TParticle(const ParticleTypeDatabase::Type& _type,
             ant::mev_t _Ek, ant::radian_t _theta, ant::radian_t _phi);

    TParticle(const ParticleTypeDatabase::Type &_type, const TLorentzVector &_lorentzvector):
        TLorentzVector(_lorentzvector),
        type(std::addressof(_type))
    {}

    TParticle(const ParticleTypeDatabase::Type& _type, const TCandidatePtr& candidate):
        TParticle(_type, candidate->CaloEnergy, candidate->Theta, candidate->Phi)
    {
        Candidate = candidate;
    }


    mev_t Ek() const { return E() - type->Mass(); }

    void ChangeType(const ParticleTypeDatabase::Type& newtype);

    using TLorentzVector::Print;
    virtual std::ostream& Print(std::ostream& stream) const;

    using TLorentzVector::Angle;
    static double CalcAngle( const TParticlePtr& p1, const TParticlePtr& p2 ) {
        return p1->Angle(p2->Vect());
    }

    template<class Archive>
    void save(Archive& archive) const {
        archive(Candidate, type->UID, Px(), Py(), Pz(), E());
    }

    template<class Archive>
    void load(Archive& archive) {
        unsigned uid;
        double px,py,pz,e;
        archive(Candidate, uid, px, py, pz, e);
        SetPxPyPzE(px,py,pz,e);
        type = std::addressof(ParticleTypeDatabase::types.at(uid));
    }

    TParticle(const TParticle&) = delete;
    TParticle& operator= (const TParticle&) = delete;
    TParticle(TParticle&&) = default;
    TParticle& operator= (TParticle&&) = default;

    TParticle() : TLorentzVector(), Candidate(nullptr), type(nullptr) {}
    virtual ~TParticle() = default;

protected:
    const ParticleTypeDatabase::Type* type;
};

}

