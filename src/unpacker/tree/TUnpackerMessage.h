#ifndef TUNPACKERMESSAGE_H
#define TUNPACKERMESSAGE_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "base/Format.h"
#endif

namespace ant {

struct TUnpackerMessage : TDataRecord
{


  TUnpackerMessage(TDataRecord::ID_t id,
                   std::uint32_t level,
                   const std::string& message) :
    TDataRecord(id),
    Level(level),
    Message(message)
  {}



  std::uint32_t Level;
  std::string Message;
  std::vector<double> Payload;

#ifndef __CINT__
  enum class Level_t : std::uint32_t {
    Info, Warn, DataError, DataDiscard, HardwareError
  };
  TUnpackerMessage(TDataRecord::ID_t id,
                   Level_t level,
                   const std::string& message) :
    TUnpackerMessage(id, static_cast<std::uint32_t>(level), message)
  {}
  std::string FormattedMessage() const {
    return fmt::format_vector(Message, Payload);
  }
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TUnpackerMessage ID=" << ID
             << " Level=" << Level << " Msg='" << FormattedMessage() << "'";
  }
#endif


  ClassDef(TUnpackerMessage, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // TUNPACKERMESSAGE_H
