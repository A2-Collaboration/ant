#ifndef THEADERINFO_H
#define THEADERINFO_H

#include "TDataRecord.h"

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
      Description(description),
      RunNumber(runnumber)
  {
    static_assert(sizeof(decltype(timestamp)) <= sizeof(decltype(Timestamp)),
                  "Bug: type of timestamp too big for THeaderInfo"
                  );
    static_assert(sizeof(runnumber) <= sizeof(decltype(RunNumber)),
                  "Bug: type of runnumber too big for THeaderInfo"
                  );
  }

  std::uint64_t Timestamp;   // unix epoch
  std::string   Description; // full descriptive string
  std::uint32_t RunNumber;   // runnumber

  ClassDef(THeaderInfo, 1)
};

}

#endif // THEADERINFO_H
