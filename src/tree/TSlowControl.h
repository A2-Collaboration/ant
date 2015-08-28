#pragma once

#include "TDataRecord.h"

#ifndef __CINT__
#include "base/std_ext/time.h"
#include <tuple>
#endif

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
  enum class Type_t : std::uint8_t {
    AcquScaler, EpicsOneShot, EpicsScaler, EpicsTimer
  };
  Type_t GetType() const {
      return static_cast<Type_t>(Type);
  }
  static const char* type_to_string(Type_t type);

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
             << " Type=" << TypeToString()
             << " Timestamp='" << std_ext::ctime(Timestamp) << "'"
             << " Name='" << Name << "'"
             << " Description='" << Description << "'";
  }

  struct Key  : printable_traits
  {
      Type_t Type;
      std::string Name;
      Key(Type_t type, const std::string& name) :
          Type(type),
          Name(name)
      {}
      Key(const TSlowControl& sc) : Key(sc.GetType(), sc.Name) {}
      bool operator<(const Key& rhs) const {
          return std::tie(Type, Name) < std::tie(rhs.Type, rhs.Name);
      }

      virtual std::ostream& Print( std::ostream& s) const override {
          return s << "TSlowControl::Key Type=" << TSlowControl::type_to_string(Type)
                   << " " << Name;
      }
  };

  Key GetKey() const {
      // use implicit conversion constructor
      return *this;
  }

#endif

  TSlowControl() : TDataRecord() {}
  ClassDef(TSlowControl, ANT_UNPACKER_ROOT_VERSION)
};

#ifndef __CINT__

inline const char* TSlowControl::TypeToString() const {
    return type_to_string(GetType());
}

inline const char* TSlowControl::type_to_string(Type_t type) {
    switch(type) {
    case Type_t::AcquScaler:
      return "AcquScaler";
    case Type_t::EpicsOneShot:
      return "EpicsOneShot";
    case Type_t::EpicsScaler:
      return "EpicsScaler";
    case Type_t::EpicsTimer:
      return "EpicsTimer";
    }
    throw std::runtime_error("Not implemented");
}
#endif


} // namespace ant
