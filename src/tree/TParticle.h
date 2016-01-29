#pragma once

#include "tree/TCandidate.h"

#include "base/printable.h"
#include "base/types.h"
#include "base/ParticleType.h"
#include "base/Tree.h"
#include "base/cereal/types/base_class.hpp"

#include "TLorentzVector.h"

#include <vector>
#include <list>
#include <iostream>
#include <memory>

namespace ant {

struct TParticle;

using TParticlePtr    = std::shared_ptr<TParticle>;
using TParticleList   = std::vector<TParticlePtr>;
using TParticleTree_t = std::shared_ptr<Tree<TParticlePtr>>;


/**
 * @brief Base TParticle class
 */
struct TParticle : TLorentzVector, printable_traits {

    TCandidatePtr Candidate;
    const ParticleTypeDatabase::Type& Type() const { return *type; }

    TParticle(const ParticleTypeDatabase::Type& type_,
             ant::mev_t Ek_, ant::radian_t theta_, ant::radian_t phi_);

    TParticle(const ParticleTypeDatabase::Type& type_, const TLorentzVector& lorentzvector_):
        TLorentzVector(lorentzvector_),
        type(std::addressof(type_))
    {}

    TParticle(const ParticleTypeDatabase::Type& type_, const TCandidatePtr& candidate):
        TParticle(type_, candidate->CaloEnergy, candidate->Theta, candidate->Phi)
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
        archive(cereal::base_class<TLorentzVector>(this), Candidate, type->UID);
    }

    template<class Archive>
    void load(Archive& archive) {
        unsigned uid;
        archive(cereal::base_class<TLorentzVector>(this), Candidate, uid);
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

