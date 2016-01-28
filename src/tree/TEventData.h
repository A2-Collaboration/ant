#pragma once

#include "TID.h"
#include "TDetectorRead.h"
#include "TCluster.h"
#include "TCandidate.h"
#include "TTagger.h"
#include "TSlowControl.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

#ifndef __CINT__
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

    TID ID;         //!
    TTagger Tagger; //!

#ifndef __CINT__

    std::vector<std::shared_ptr<TCluster>>   Clusters;
    //std::vector<std::shared_ptr<TCandidate>> Candidates;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Clusters);
    }

    virtual std::ostream& Print(std::ostream& s) const override;
#endif

    ClassDef(TEventData, ANT_UNPACKER_ROOT_VERSION)

}; // TEventData


} // namespace ant
