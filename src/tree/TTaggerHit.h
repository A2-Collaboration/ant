#ifndef ANT_TAGGERHIT_H
#define ANT_TAGGERHIT_H

#include "TDataRecord.h"
#include "TCluster.h"
#include <vector>

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TTaggerHit: ant::printable_traits
#else
struct TTaggerHit
#endif
{
    std::vector<ant::TCluster> Clusters;
    double PhotonEnergy;
    double Time;

    TTaggerHit() {}
    TTaggerHit(double photonE, double t): PhotonEnergy(photonE), Time(t) {}

    virtual ~TTaggerHit() {}

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TTaggerHit: " << Clusters.size() << " clusters PhotonEnergy=" << PhotonEnergy << " Time=" << Time;
  }
#endif

    ClassDef(TTaggerHit, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TAGGERHIT_H
