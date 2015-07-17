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

struct TClusterHitDatum
{
  std::uint8_t Type;
  double Value;
  std::int16_t ValueInt;

#ifndef __CINT__
  TClusterHitDatum(Channel_t::Type_t type, double value, std::int16_t value_int = 0) :
    Type(static_cast<std::uint8_t>(type)),
    Value(value),
    ValueInt(value_int)
  {}
  Channel_t::Type_t GetType() const {
    return static_cast<Channel_t::Type_t>(Type);
  }
#endif

  TClusterHitDatum() {}
  virtual ~TClusterHitDatum() {}
  ClassDef(TClusterHitDatum, ANT_UNPACKER_ROOT_VERSION)
};

#ifndef __CINT__
struct TClusterHit: public ant::printable_traits
#else
struct TClusterHit
#endif
{
  std::uint32_t Channel;
  std::vector<TClusterHitDatum> Data;

#ifndef __CINT__

  TClusterHit(unsigned channel,
              const std::vector<TClusterHitDatum>& data):
    Channel(channel),
    Data(data) {
    static_assert(sizeof(Channel)>=sizeof(channel),
                  "Parameter channel does not fit into TClusterHit::Channel");

  }

  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TClusterHit: ";
    for(auto& datum : Data) {
      s << Channel_t::ToString(datum.GetType()) << "=" << datum.Value << " ";
    }
    return s;
  }
#endif

  TClusterHit() {}
  virtual ~TClusterHit() {}
  ClassDef(TClusterHit, ANT_UNPACKER_ROOT_VERSION)

};


#ifndef __CINT__
struct TCluster: public ant::printable_traits
#else
struct TCluster
#endif
{

  std::vector<ant::TClusterHit> Hits;
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
    return s << "TCluster: " << Hits.size() << " hits @" << Position <<", Energy=" << Energy
             << " Detector=" << Detector_t::ToString(GetDetectorType());
  }
#endif

  TCluster(): Position(), Energy(0.0), DetectorType(0) {}
  virtual ~TCluster() {}
  ClassDef(TCluster, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // ANT_TCLUSTER_H
