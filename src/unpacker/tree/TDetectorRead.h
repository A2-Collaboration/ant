#ifndef TDETECTORREAD_H
#define TDETECTORREAD_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "expconfig/ExpConfig.h"
#endif

namespace ant {

struct TDetectorRead : TDataRecord
{
  TDetectorRead(const TDataRecord::ID_t& id) :
    TDataRecord(id),
    Hits()
  {}

#ifndef __CINT__
  struct Hit  : printable_traits
#else
  struct Hit
#endif
  {
    std::uint8_t  Detector;
    std::uint8_t  Type;
    std::uint32_t Channel;

    std::vector<std::uint8_t>  RawData;

    std::vector<double>        Values;
    std::vector<bool>          ValueBits;

#ifndef __CINT__
    Hit(const LogicalChannel_t& element,
        const std::vector<std::uint8_t>& rawData) :
      Detector(static_cast<std::uint8_t>(element.Detector)),
      Type(static_cast<std::uint8_t>(element.Type)),
      Channel(element.Channel),
      RawData(rawData)
    {
      static_assert(sizeof(Channel)>=sizeof(decltype(element.Channel)),
                    "LogicalElement_t::Channel does not fit into this.Channel");
    }

    virtual std::ostream& Print( std::ostream& s) const override {
      return s << "Hit Detector=" << static_cast<int>(Detector);
    }
#endif

    Hit() {}
    virtual ~Hit() {}
    ClassDef(Hit, ANT_UNPACKER_ROOT_VERSION)
  };

  std::vector<Hit> Hits;




#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TDetectorRead ID=" << ID;
  }
#endif

  TDetectorRead() : TDataRecord() {}
  ClassDef(TDetectorRead, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // TDETECTORREAD_H
