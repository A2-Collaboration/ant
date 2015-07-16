#ifndef ANT_TUNPACKERMESSAGE_H
#define ANT_TUNPACKERMESSAGE_H

#include "TDataRecord.h"

#ifndef __CINT__
#include "base/Format.h"
#endif

namespace ant {

struct TUnpackerMessage : TDataRecord
{


  TUnpackerMessage(TID id,
                   std::uint32_t level,
                   const std::string& message) :
    TDataRecord(id),
    Level(level),
    Message(message)
  {}



  std::uint8_t Level;
  std::string Message;
  std::vector<double> Payload;

  const char* LevelToString() const;

#ifndef __CINT__
  enum class Level_t : std::uint8_t {
    Info, Warn, DataError, DataDiscard, HardwareError
  };
  static const char* LevelToString(const Level_t&);
  TUnpackerMessage(TID id,
                   Level_t level,
                   const std::string& message) :
    TUnpackerMessage(id, static_cast<decltype(Level)>(level), message)
  {}
  std::string FormattedMessage() const {
    return fmt::format_vector(Message, Payload);
  }
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TUnpackerMessage ID=" << ID
             << " Level=" << LevelToString()
             << " Msg='" << FormattedMessage() << "'";
  }
#endif

  TUnpackerMessage() : TDataRecord() {}
  ClassDef(TUnpackerMessage, ANT_UNPACKER_ROOT_VERSION)
};

#ifndef __CINT__
inline const char* TUnpackerMessage::LevelToString(const Level_t& level) {
  switch(level) {
  case Level_t::Info:
    return "Info";
  case Level_t::Warn:
    return "Warn";
  case Level_t::DataError:
    return "DataError";
  case Level_t::DataDiscard:
    return "DataDiscard";
  case Level_t::HardwareError:
    return "HardwareError";
  }
  throw std::runtime_error("Not implemented");
}
inline const char* TUnpackerMessage::LevelToString() const {
  return LevelToString(static_cast<Level_t>(Level));
}
#endif

} // namespace ant

#endif // ANT_TUNPACKERMESSAGE_H
