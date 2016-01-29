#pragma once

#ifndef __CINT__
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
#endif

namespace ant {

#ifndef __CINT__
struct TEvent : printable_traits
#else
struct TEvent
#endif
{


#ifndef __CINT__

    struct Data : printable_traits
    {
        Data(const TID& id) : ID(id) {}
        Data() {}
        virtual ~Data() {}

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

        virtual std::ostream& Print(std::ostream& s) const override;

    }; // Data

    using DataPtr = std::shared_ptr<Data> ;

    DataPtr Reconstructed;
    DataPtr MCTrue;

    template<class Archive>
    void serialize(Archive archive) {
        archive(Reconstructed, MCTrue);
    }

    virtual std::ostream& Print( std::ostream& s) const override;

#endif

    TEvent() {}
    virtual ~TEvent() {}
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)
};

}
