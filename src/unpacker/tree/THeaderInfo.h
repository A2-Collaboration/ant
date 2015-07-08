#ifndef THEADERINFO_H
#define THEADERINFO_H

#include "TDataRecord.h"

#ifndef __CINT__
#include <ctime>
#include "base/std_ext.h"
#endif

namespace ant {

struct THeaderInfo : TDataRecord
{

  THeaderInfo(
      const TDataRecord::ID_t& id,
      std::time_t timestamp,
      const std::string& description,
      unsigned runnumber
      )
    :
      TDataRecord(id),
      Timestamp(timestamp),
      RunNumber(runnumber),
      Description(description)
  {
    static_assert(sizeof(decltype(timestamp)) <= sizeof(decltype(Timestamp)),
                  "Bug: type of timestamp too big for THeaderInfo"
                  );
    static_assert(sizeof(runnumber) <= sizeof(decltype(RunNumber)),
                  "Bug: type of runnumber too big for THeaderInfo"
                  );
  }

  std::uint64_t Timestamp;   // unix epoch
  std::uint32_t RunNumber;   // runnumber
  std::string   Description; // full descriptive string

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "THeaderInfo ID=" << ID
             << " Timestamp='" << std_ext::ctime(Timestamp) << "'"
             << " RunNumber=" << RunNumber
             << " Description='" << Description << "'";

  }
#endif

  ClassDef(THeaderInfo, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // THEADERINFO_H
