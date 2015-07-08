#ifndef TUNPACKERMESSAGE_H
#define TUNPACKERMESSAGE_H

#include "TDataRecord.h"
#include "base/Format.h"
#include <string>


namespace ant {

struct TUnpackerMessage : TDataRecord
{

  enum class Level_t : unsigned {
    Info, Warn, SoftError, HardError
  };

  TUnpackerMessage(TDataRecord::ID_t id,
                   Level_t level,
                   const std::string& message) :
    TDataRecord(id),
    Level(level),
    Message(message)
  {}

  Level_t Level;
  std::string Message;
  /// \todo Support formatted message, but still store the arguments?
};

}

#endif // TUNPACKERMESSAGE_H
