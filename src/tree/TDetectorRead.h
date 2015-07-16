#ifndef ANT_TDETECTORREAD_H
#define ANT_TDETECTORREAD_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "expconfig/Detector_t.h"
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TDetectorReadHit  : printable_traits
#else
struct TDetectorReadHit
#endif
{
  std::uint8_t  Detector;
  std::uint8_t  Type;
  std::uint32_t Channel;

  std::vector<std::uint8_t>  RawData;

  std::vector<double>        Values;
  std::vector<bool>          ValueBits;

  const char* GetDetectorAsString() const;
  const char* GetTypeAsString() const;

#ifndef __CINT__
  TDetectorReadHit(const LogicalChannel_t& element,
      const std::vector<std::uint8_t>& rawData) :
    Detector(static_cast<std::uint8_t>(element.Detector)),
    Type(static_cast<std::uint8_t>(element.Type)),
    Channel(element.Channel),
    RawData(rawData)
  {
    static_assert(sizeof(Channel)>=sizeof(decltype(element.Channel)),
                  "LogicalElement_t::Channel does not fit into TDetecorReadHit::Channel");
  }



  virtual std::ostream& Print( std::ostream& s) const override {
    std::ostringstream rawdata;

    rawdata << std::hex << std::uppercase << std::setfill( '0' );
    for(int c : RawData) {
      rawdata << std::setw( 2 ) << c;
    }

    return s << "Hit Detector="
                //<< std::left << std::setfill(' ') << std::setw(4)
             << GetDetectorAsString() << std::right
             << " Channel="
             << std::setw(3)
             << Channel
             << " Type="
             << GetTypeAsString()
             << " RawData=0x" << rawdata.str()
                ;
  }
#endif

  TDetectorReadHit() {}
  virtual ~TDetectorReadHit() {}
  ClassDef(TDetectorReadHit, ANT_UNPACKER_ROOT_VERSION)
};

struct TDetectorRead : TDataRecord
{
  TDetectorRead(const TID& id) :
    TDataRecord(id),
    Hits()
  {}



  std::vector<TDetectorReadHit> Hits;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TDetectorRead ID=" << ID << " Hits=" << Hits.size() << '\n';
    for(size_t i=0;i<Hits.size();i++) {
      s << "  i="
        << std::setw(3)
        << i << " "
        << Hits[i] << '\n';
    }
    return s;
  }
#endif

  TDetectorRead() : TDataRecord() {}
  ClassDef(TDetectorRead, ANT_UNPACKER_ROOT_VERSION)

};

#ifndef __CINT__
inline const char* TDetectorReadHit::GetDetectorAsString() const {
  return Detector_t::ToString(static_cast<Detector_t::Type_t>(Detector));
}
inline const char* TDetectorReadHit::GetTypeAsString() const {
  return Channel_t::ToString(static_cast<Channel_t::Type_t>(Type));
}
#endif

}

#endif // ANT_TDETECTORREAD_H
