#pragma once

#include "tree/TCandidate.h"

#include "base/printable.h"
#include "base/types.h"
#include "base/ParticleType.h"
#include "base/Tree.h"

#include "base/vec/LorentzVec.h"

#include <vector>
#include <list>
#include <iostream>
#include <memory>

namespace ant {

struct TParticle;

using TParticlePtr  = std_ext::cc_shared_ptr<TParticle>;
using TParticleList = std::vector<TParticlePtr>;
using TParticleTree_t = Tree<TParticlePtr>::node_t;


/**
 * @brief Base TParticle class
 */
struct TParticle : LorentzVec, printable_traits {

    TCandidatePtr Candidate;
    const ParticleTypeDatabase::Type& Type() const { return *type; }

    TParticle(const ParticleTypeDatabase::Type& type_,
             ant::mev_t Ek_, ant::radian_t theta_, ant::radian_t phi_);

    TParticle(const ParticleTypeDatabase::Type& type_, const LorentzVec& lorentzvector_):
        LorentzVec(lorentzvector_),
        type(std::addressof(type_))
    {}

    TParticle(const ParticleTypeDatabase::Type& type_, const TCandidatePtr& candidate):
        TParticle(type_, candidate->CaloEnergy, candidate->Theta, candidate->Phi)
    {
        Candidate = candidate;
    }


    mev_t Ek() const { return E - type->Mass(); }

    void ChangeType(const ParticleTypeDatabase::Type& newtype);

    virtual std::ostream& Print(std::ostream& stream) const;


    static double CalcAngle( const TParticlePtr& p1, const TParticlePtr& p2 ) {
        return p1->Angle(*p2);
    }

    template<class Archive>
    void save(Archive& archive) const {
        archive(static_cast<const LorentzVec&>(*this), Candidate, type->UID);
    }

    template<class Archive>
    void load(Archive& archive) {
        unsigned uid;
        archive(static_cast<LorentzVec&>(*this), Candidate, uid);
        type = std::addressof(ParticleTypeDatabase::types.at(uid));
    }

    TParticle(const TParticle&) = delete;
    TParticle& operator= (const TParticle&) = delete;
    TParticle(TParticle&&) = default;
    TParticle& operator= (TParticle&&) = default;

    TParticle() : LorentzVec(), Candidate(nullptr), type(nullptr) {}
    virtual ~TParticle() = default;

protected:
    const ParticleTypeDatabase::Type* type;
};

}

