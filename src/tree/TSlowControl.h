#ifndef ANT_TSLOWCONTROL_H
#define ANT_TSLOWCONTROL_H

#include "TDataRecord.h"

namespace ant {

// a slow control record
// can carry very different types of payload
struct TSlowControl : TDataRecord
{
  TSlowControl(const TID& id,
               std::uint8_t type,
               std::time_t timestamp,
               const std::string& name,
               const std::string& description
               ) :
    TDataRecord(id),
    Type(type),
    Timestamp(timestamp),
    Name(name),
    Description(description)
  {
    static_assert(sizeof(decltype(timestamp)) <= sizeof(decltype(Timestamp)),
                  "Bug: type of timestamp too big for TSlowControl"
                  );
  }

  std::uint8_t Type;
  std::int64_t Timestamp;   // unix epoch, only meaningful for Epics data
  std::string  Name;
  std::string  Description;

  std::vector< TKeyValue<std::int64_t> > Payload_Int;
  std::vector< TKeyValue<double> >       Payload_Float;
  std::vector< TKeyValue<std::string> >  Payload_String;

  const char* TypeToString() const;

#ifndef __CINT__
  /// \todo Maybe those types are bit too Acqu-like...?!
  enum class Type_t : std::uint8_t {
    AcquScaler, EpicsOneShot, EpicsScaler, EpicsTimer
  };
  static const char* TypeToString(const Type_t&);
  TSlowControl(TID id,
               Type_t type,
               std::time_t timestamp,
               const std::string& name,
               const std::string& description) :
    TSlowControl(id, static_cast<decltype(Type)>(type),
                 timestamp, name, description)
  {}
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TSlowControl ID=" << ID
             << " Type=" << static_cast<int>(Type)
             << " Timestamp='" << std_ext::ctime(Timestamp) << "'"
             << " Name='" << Name << "'"
             << " Description='" << Description << "'";
  }
#endif

  TSlowControl() : TDataRecord() {}
  ClassDef(TSlowControl, ANT_UNPACKER_ROOT_VERSION)
};

#ifndef __CINT__
inline const char* TSlowControl::TypeToString(const TSlowControl::Type_t& level) {
  switch(level) {
  case Type_t::AcquScaler:
    return "AcquScaler";
  case Type_t::EpicsOneShot:
    return "EpicsOneShort";
  case Type_t::EpicsScaler:
    return "EpicsScaler";
  case Type_t::EpicsTimer:
    return "EpicsTimer";
  }
  throw std::runtime_error("Not implemented");
}
inline const char* TSlowControl::TypeToString() const {
  return TypeToString(static_cast<Type_t>(Type));
}
#endif


} // namespace ant

#endif // ANT_TSLOWCONTROL_H
