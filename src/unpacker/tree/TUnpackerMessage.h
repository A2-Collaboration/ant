#ifndef TUNPACKERMESSAGE_H
#define TUNPACKERMESSAGE_H

#include "TDataRecord.h"

namespace ant {

struct TUnpackerMessage : TDataRecord
{

#ifndef __CINT__
  enum class Level_t : std::uint32_t {
    Info, Warn, Error
  };
#endif

  TUnpackerMessage(TDataRecord::ID_t id,
                   std::uint32_t level,
                   const std::string& message) :
    TDataRecord(id),
    Level(level),
    Message(message)
  {}

  TUnpackerMessage(TDataRecord::ID_t id,
                   Level_t level,
                   const std::string& message) :
    TUnpackerMessage(id, static_cast<std::uint32_t>(level), message)
  {}

  std::uint32_t Level;
  std::string Message;

  /// \todo Support formatted message, but still store the arguments?


  ClassDef(TUnpackerMessage, 1)
};

}

#endif // TUNPACKERMESSAGE_H
