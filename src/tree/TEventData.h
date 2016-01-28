#pragma once

#include "TID.h"
#include "TDetectorRead.h"
#include "TTagger.h"
#include "TSlowControl.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

#ifndef __CINT__
#include "TCandidate.h"
#include "TCluster.h"
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

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Clusters, Candidates);
    }

    virtual std::ostream& Print(std::ostream& s) const override;
#endif

    TEventData() : ID() {}
    virtual ~TEventData() {}
    ClassDef(TEventData, ANT_UNPACKER_ROOT_VERSION)

}; // TEventData


} // namespace ant
