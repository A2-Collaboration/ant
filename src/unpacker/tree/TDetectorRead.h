#ifndef TDETECTORREAD_H
#define TDETECTORREAD_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "expconfig/ExpConfig.h"
#endif

namespace ant {

struct TDetectorRead : TDataRecord
{


  std::vector<std::uint8_t>  RawData;
  std::vector<double>        Values;
  std::vector<bool>          ValueBits;

  ClassDef(TDetectorRead, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // TDETECTORREAD_H
