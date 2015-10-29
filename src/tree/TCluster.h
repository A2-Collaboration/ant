#pragma once

#include "TDataRecord.h"


#include <TVector3.h>
#include <vector>
#include <cmath>

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#include "base/root_printable.h"
#include "base/Detector_t.h"
#endif

namespace ant {

struct TClusterHitDatum
{
  std::uint8_t Type;
  double Value;

#ifndef __CINT__
  TClusterHitDatum(Channel_t::Type_t type, double value) :
    Type(static_cast<std::uint8_t>(type)),
    Value(value)
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
    s << "TClusterHit Ch=" << Channel << ": ";
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

  TVector3 Position;
  double Energy;
  double Time;
  std::uint32_t Flags;
  std::uint8_t DetectorType;
  std::uint32_t CentralElement;

  std::vector<TClusterHit> Hits;

#ifndef __CINT__

  TCluster(
      const TVector3& pos,
      double E,
      double t,
      const Detector_t::Type_t& type,
      const unsigned central,
      const std::vector<TClusterHit>& hits = {}
      ):
    Position(pos),
    Energy(E),
    Time(t),
    Flags(0),
    DetectorType(static_cast<std::uint8_t>(type)),
    CentralElement(central),
    Hits(hits)
  {}

  Detector_t::Type_t GetDetectorType() const {
    return static_cast<Detector_t::Type_t>(DetectorType);
  }

  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TCluster: " << Hits.size() << " hits @" << Position <<", Energy=" << Energy
             << " Central Element=" << CentralElement
             << " Detector=" << Detector_t::ToString(GetDetectorType());
  }

  // the cluster alorithm may set those flags
  enum class Flags_t : std::uint8_t {
      TouchesHole, Split, Unmatched
  };

  void SetFlag(Flags_t flag, bool set = true) {
      const std::uint32_t mask = 1 << static_cast<std::uint8_t>(flag);
      if(set) {
          Flags |= mask;
      }
      else {
          Flags &= ~mask;
      }
  }
  bool HasFlag(Flags_t flag) const {
      const std::uint32_t mask = 1 << static_cast<std::uint8_t>(flag);
      return (Flags & mask) != 0;
  }


#endif

  bool isSane() const {
      return std::isfinite(Energy) && std::isfinite(Time);
  }

  TCluster(): Position(), Energy(0.0), Time(0.0), Flags(0), DetectorType(0), CentralElement(0) {}
  virtual ~TCluster() {}
  ClassDef(TCluster, ANT_UNPACKER_ROOT_VERSION)
};

}
