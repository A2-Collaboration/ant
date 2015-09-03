#pragma once

#include "TID.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

namespace ant {

#ifndef __CINT__
struct TDataRecord : printable_traits
#else
struct TDataRecord
#endif
{
    TDataRecord() : ID() {}
    TDataRecord(const TID& id) : ID(id) {}
    virtual ~TDataRecord() {}

    TID ID;

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TDataRecord ID=" << ID;
    }
#endif

    ClassDef(TDataRecord, ANT_UNPACKER_ROOT_VERSION)

}; // TDataRecord


} // namespace ant
#pragma once
