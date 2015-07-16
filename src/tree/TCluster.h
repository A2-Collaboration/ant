#ifndef ANT_TCLUSTER_H
#define ANT_TCLUSTER_H

#include "TDataRecord.h"
#include "TDetectorRead.h"
#include <TVector3.h>
#include <vector>

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#include "base/root_printable.h"
#endif

namespace ant {
#ifndef __CINT__
struct TCluster: public ant::printable_traits
#else
struct TCluster
#endif
{
    TCluster(): Position(), Energy(0.0), DetectorType(0) {}



    virtual ~TCluster() {}

    std::vector<ant::TDetectorReadHit> Hits;
    TVector3 Position;
    double Energy;
    std::uint8_t DetectorType;

#ifndef __CINT__

    TCluster(const TVector3& pos, double E, const ant::Detector_t::Type_t& type):
        Position(pos),
        Energy(E),
        DetectorType(static_cast<std::uint8_t>(type)) {}

    Detector_t::Type_t GetDetectorType() const {
        return static_cast<Detector_t::Type_t>(DetectorType);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TCluster: " << Hits.size() << " hits @" << Position <<", Energy=" << Energy << " Detector=" << Detector_t::ToString(GetDetectorType());
    }
#endif

    ClassDef(TCluster, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TCLUSTER_H
