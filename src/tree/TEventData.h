#pragma once

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

struct TEventData
{
    TEventData(const TID& id);
    TEventData();

    TID ID;
    std::vector<TDetectorReadHit> DetectorReadHits;
    std::vector<TSlowControl>     SlowControls;
    std::vector<TUnpackerMessage> UnpackerMessages;

    std::vector<TTaggerHit>  TaggerHits;
    TTrigger Trigger;
    TTarget  Target;

    TClusterList     Clusters;
    TCandidateList   Candidates;
    TParticleTree_t  ParticleTree; // only on MC

    template<class Archive>
    void serialize(Archive& archive) {
        archive(ID,
                DetectorReadHits, SlowControls, UnpackerMessages,
                TaggerHits, Trigger, Target,
                Clusters, Candidates, ParticleTree);
    }

    friend std::ostream& operator<<(std::ostream& s, const TEventData& o);

    void ClearDetectorReadHits();

};

}