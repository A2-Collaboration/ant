#pragma once

#include "TID.h"
#include "TDetectorRead.h"
#include "TTagger.h"
#include "TSlowControl.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

#ifndef __CINT__
#include "TCluster.h"
#include "TCandidate.h"
#include "TParticle.h"
#include <memory>
#endif


namespace ant {

struct TDetectorReadHits : std::vector<TDetectorReadHit> {
    TDetectorReadHits() : std::vector<TDetectorReadHit>() {}
    virtual ~TDetectorReadHits() {}
    ClassDef(TDetectorReadHits, ANT_UNPACKER_ROOT_VERSION)
};

#ifndef __CINT__
struct TEventData : printable_traits
#else
struct TEventData
#endif
{
    TEventData(const TID& id) : ID(id) {}

    // we have a custom Streamer method, so mark all
    // members visible to ROOTcint as transient with comment //!
    TID ID;                               //!
    TDetectorReadHits DetectorHits;       //!
    TTagger Tagger;                       //!

#ifndef __CINT__

    std::vector<TClusterPtr>   Clusters;
    std::vector<TCandidatePtr> Candidates;
    std::vector<TParticlePtr>  Particles; // final state, or identified from reconstructed candidates


    template<class Archive>
    void serialize(Archive& archive) {
        archive(Clusters, Candidates, Particles);
    }

    virtual std::ostream& Print(std::ostream& s) const override;
#endif

    TEventData() : ID() {}
    virtual ~TEventData() {}
    ClassDef(TEventData, ANT_UNPACKER_ROOT_VERSION)

}; // TEventData


} // namespace ant
