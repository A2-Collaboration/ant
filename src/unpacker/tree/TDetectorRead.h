#ifndef TDETECTORREAD_H
#define TDETECTORREAD_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "expconfig/ExpConfig.h"
#endif

namespace ant {

struct TDetectorRead : TDataRecord
{

#ifndef __CINT__
  struct Hit  : printable_traits
#else
  struct Hit
#endif
  {
    std::uint8_t  Detector;
    std::uint8_t  Kind;
    std::uint32_t LogicalChannel;

    std::vector<std::uint8_t>  RawData;
    std::vector<double>        Values;
    std::vector<bool>          ValueBits;

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
      return s << "Hit " << Detector;
    }
#endif

    virtual ~Hit() {}
    ClassDef(Hit, ANT_UNPACKER_ROOT_VERSION)
  };

  std::vector<Hit> Hits;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TDetectorRead ID=" << ID;
  }
#endif

  ClassDef(TDetectorRead, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // TDETECTORREAD_H
