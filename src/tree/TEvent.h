#pragma once

#include "Rtypes.h"

#ifndef __CINT__
#include <memory>
#include <stdexcept>
#endif

#define ANT_TEVENT_VERSION 5

namespace ant {

#ifndef __CINT__
struct TID;
struct TEventData;
#endif


struct TEvent
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

    friend std::ostream& operator<<( std::ostream& s, const TEvent& o);

    explicit TEvent(const TID& id_reconstructed);
    explicit TEvent(const TID& id_reconstructed, const TID& id_mctrue);

    // TEvent is moveable
    TEvent(TEvent&&);
    TEvent& operator=(TEvent&&);

protected:
    // exclamation mark at the beginning of the comment below tells ROOT
    // to exclude the data members from the Streamer (added because of ROOT6)
    std::unique_ptr<TEventData> reconstructed;  //! reconstructed detector information, either Geant or raw data
    std::unique_ptr<TEventData> mctrue;  //! MC true information from event generator

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
