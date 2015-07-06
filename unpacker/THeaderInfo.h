#ifndef THEADERINFO_H
#define THEADERINFO_H

#include "TDataRecord.h"

#include <string>
#include <vector>
#include <ctime>

namespace ant {

class THeaderInfo : public TDataRecord
{
public:
  THeaderInfo(
      TDataRecord::UUID_t uuid,
      std::time_t timestamp,
      const std::string& description,
      unsigned runnumber
      )
    :
      TDataRecord(uuid),
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

private:
  uint64_t Timestamp; // unix epoch
  std::string Description; // full descriptive string
  uint32_t RunNumber;
};

}

#endif // THEADERINFO_H
