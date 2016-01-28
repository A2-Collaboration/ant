#pragma once

#include "TID.h"
#include "TDetectorRead.h"
#include "TCluster.h"
#include "TCandidate.h"
#include "TTagger.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

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

    TID ID;
    std::vector<TDetectorReadHit> DetectorHits;
    std::vector<TCluster>         Clusters;
    std::vector<TCandidate>       Candidates;

    // Particles, ParticleTree

    TTagger Tagger;


#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override;
#endif

    ClassDef(TEventData, ANT_UNPACKER_ROOT_VERSION)

}; // TEventData


} // namespace ant
