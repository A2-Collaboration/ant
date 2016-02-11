#pragma once

#include "Rtypes.h"

#ifndef __CINT__
#include "base/printable.h"
#include <memory>
#include <stdexcept>
#endif

#define ANT_TEVENT_VERSION 0

namespace ant {

#ifndef __CINT__
struct TID;
struct TEventData;
struct TEvent;
using TEventPtr = std::unique_ptr<TEvent>;
#endif


#ifndef __CINT__
struct TEvent : printable_traits
#else
struct TEvent
#endif
{


#ifndef __CINT__

    std::unique_ptr<TEventData> Reconstructed;
    std::unique_ptr<TEventData> MCTrue;
    // indicates that this event was only saved for SlowControl processing
    bool SavedForSlowControls = false;

    template<class Archive>
    void serialize(Archive archive, const std::uint32_t version) {
        if(version != ANT_TEVENT_VERSION)
            throw std::runtime_error("TEvent version mismatch");
        archive(Reconstructed, MCTrue, SavedForSlowControls);
    }

    virtual std::ostream& Print( std::ostream& s) const override;

    static std::unique_ptr<TEvent> MakeReconstructed(const TID& id);

    void ClearDetectorReadHits();

    // TEvent is moveable
    TEvent(TEvent&&) = default;
    TEvent& operator=(TEvent&&) = default;

#endif

    TEvent();
    virtual ~TEvent();
    ClassDef(TEvent, ANT_TEVENT_VERSION)

private:
    // prevent ROOTcint from creating copy-constructors
    TEvent(const TEvent&);
    TEvent& operator=(const TEvent&);

};

}
