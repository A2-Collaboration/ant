#pragma once

#include "TDataRecord.h"
#include "TTaggerHit.h"

namespace ant {

#ifndef __CINT__
struct TTagger : printable_traits
#else
struct TTagger
#endif
{
    std::vector<ant::TTaggerHit> Hits;
    std::vector< TKeyValue<double> > Scalers;
    std::vector< TKeyValue<double> > Energies;

    void Clear() {
        Hits.resize(0);
        Scalers.resize(0);
        Energies.resize(0);
    }


    /// \todo add fields for Moeller and PairSpec

    TTagger() {}

    virtual ~TTagger() {}

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        return s;
    }
#endif

    ClassDef(TTagger, ANT_UNPACKER_ROOT_VERSION)

};

}
