#pragma once

#include "TID.h"
#include "TDetectorReadHit.h"
#include "TSlowControl.h"
#include "TUnpackerMessage.h"

#include "TTagger.h"
#include "TCluster.h"
#include "TCandidate.h"
#include "TParticle.h"

#include <iomanip>
#include <ctime>
#include <memory>

namespace ant {

struct TEventData;

using TEventDataPtr = std::shared_ptr<TEventData>;

struct TEventData : printable_traits
{
    TEventData(const TID& id) : ID(id) {}
    TEventData() {}
    virtual ~TEventData() {}

    TID ID;
    std::vector<TDetectorReadHit> DetectorReadHits;
    std::vector<TSlowControl>     SlowControls;
    std::vector<TUnpackerMessage> UnpackerMessages;

    TTagger                    Tagger;
    std::vector<TClusterPtr>   Clusters;
    std::vector<TCandidatePtr> Candidates;
    std::vector<TParticlePtr>  Particles;    // MCTrue final state, or identified from reconstructed candidates
    TParticleTree_t            ParticleTree;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(ID,
                DetectorReadHits, SlowControls, UnpackerMessages,
                Tagger, Clusters, Candidates, Particles, ParticleTree);
    }

    virtual std::ostream& Print(std::ostream& s) const override {
        return s << "TEventData ID=" << ID;
    }

}; // TEventData


} // namespace ant
