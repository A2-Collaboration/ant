#pragma once

#include "TEventData.h"

#ifndef __CINT__
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
    TEventData Reconstructed;
    TEventData MCTrue;

#ifndef __CINT__

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "=== Reconstructed:\n" << Reconstructed
                 << "=== MCTrue:\n" << MCTrue;
    }
#endif

    TEvent() : Reconstructed(), MCTrue() {}
    virtual ~TEvent() {}
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)
#ifdef __clang__
#pragma clang diagnostic pop
#endif
};

}
