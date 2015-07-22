#ifndef ANT_TAGGERHIT_H
#define ANT_TAGGERHIT_H

#include "TDataRecord.h"
#include <vector>

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TTaggerHit : printable_traits
#else
struct TTaggerHit
#endif
{
    double PhotonEnergy;
    double Time;
    std::vector< TKeyValue<double> > Electrons; // Key=channel, Value=timing


#ifndef __CINT__
    TTaggerHit(double photonE, const TKeyValue<double>& hit) :
        PhotonEnergy(photonE),
        Time(hit.Value),
        Electrons{hit}
    {}
    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TTaggerHit: Electrons=" << Electrons.size()
                 << " PhotonEnergy=" << PhotonEnergy << " Time=" << Time;
    }
#endif

    TTaggerHit() {}
    virtual ~TTaggerHit() {}
    ClassDef(TTaggerHit, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // ANT_TAGGERHIT_H
