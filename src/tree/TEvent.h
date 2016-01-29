#pragma once

#ifndef __CINT__
#include "TEventData.h"
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TEvent : printable_traits
#else
struct TEvent
#endif
{


#ifndef __CINT__

    virtual std::ostream& Print( std::ostream& s) const override;

    TEventDataPtr Reconstructed;
    TEventDataPtr MCTrue;

    template<class Archive>
    void serialize(Archive archive) {
        archive(Reconstructed, MCTrue);
    }

#endif

    TEvent() {}
    virtual ~TEvent() {}
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)
};

}
