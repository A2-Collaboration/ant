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

#ifndef __CINT__
struct TEventData : printable_traits
#else
struct TEventData
#endif
{
    TEventData() : ID() {}
    TEventData(const TID& id) : ID(id) {}
    virtual ~TEventData() {}

    // we have a custom Streamer method, so mark all
    // members visible to ROOTcint as transient with comment //!
    TID ID;                               //!
    TTagger Tagger;                       //!
    std::vector<TDetectorReadHit> Hits;   //!

#ifndef __CINT__

    std::vector<std::shared_ptr<TCluster>>   Clusters;
    std::vector<std::shared_ptr<TCandidate>> Candidates;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Clusters);
    }

    virtual std::ostream& Print(std::ostream& s) const override;
#endif

    ClassDef(TEventData, ANT_UNPACKER_ROOT_VERSION)

}; // TEventData


} // namespace ant
