#pragma once

#include "Rtypes.h"

#ifndef __CINT__
#include "base/printable.h"
#include <memory>
#include <stdexcept>
#endif

#define ANT_TEVENT_VERSION 3

namespace ant {

#ifndef __CINT__
struct TID;
struct TEventData;
#endif


#ifndef __CINT__
struct TEvent : printable_traits
#else
struct TEvent
#endif
{

#ifndef __CINT__

    const TEventData& Reconstructed() const { return *reconstructed; }
    TEventData& Reconstructed() { return *reconstructed; }
    const TEventData& MCTrue() const { return *mctrue; }
    TEventData& MCTrue() { return *mctrue; }

    explicit operator bool() const {
        return reconstructed || mctrue;
    }

    // indicates that this event was only saved for SlowControl processing
    bool SavedForSlowControls = false;

    template<class Archive>
    void serialize(Archive& archive, const std::uint32_t version) {
        if(version != ANT_TEVENT_VERSION)
            throw std::runtime_error("TEvent version mismatch");
        archive(reconstructed, mctrue, SavedForSlowControls);
    }

    virtual std::ostream& Print( std::ostream& s) const override;

    explicit TEvent(const TID& id_reconstructed);
    explicit TEvent(const TID& id_reconstructed, const TID& id_mctrue);

    // TEvent is moveable
    TEvent(TEvent&&);
    TEvent& operator=(TEvent&&);

protected:
    std::unique_ptr<TEventData> reconstructed;
    std::unique_ptr<TEventData> mctrue;

#endif

public:

    TEvent();
    virtual ~TEvent();
    ClassDef(TEvent, ANT_TEVENT_VERSION)

private:
    // prevent ROOTcint from creating copy-constructors
    TEvent(const TEvent&);
    TEvent& operator=(const TEvent&);

};

}
