#pragma once

#include "TID.h" // for TKeyValue
#include "base/std_ext/time.h"

#include <cstdint>
#include <tuple>

namespace ant {

struct TSlowControl
{
    enum class Type_t : std::uint8_t {
        AcquScaler, EpicsOneShot, EpicsScaler, EpicsTimer, Ant
    };

    /**
      * @brief The Validity_t enum determines if the slowcontrol item
      * is valid for the upcoming or the previous events
      */
    enum class Validity_t : std::uint8_t {
        Backward, Forward
    };

    Type_t Type;
    Validity_t Validity;
    std::int64_t Timestamp;   // unix epoch, only meaningful for some Type items
    std::string  Name;
    std::string  Description;

    std::vector< TKeyValue<std::int64_t> > Payload_Int;
    std::vector< TKeyValue<double> >       Payload_Float;
    std::vector< TKeyValue<std::string> >  Payload_String;

    TSlowControl(Type_t type,
                 Validity_t validity,
                 std::time_t timestamp,
                 const std::string& name,
                 const std::string& description
                 ) :
        Type(type),
        Validity(validity),
        Timestamp(timestamp),
        Name(name),
        Description(description)
    {
        static_assert(sizeof(decltype(timestamp)) <= sizeof(decltype(Timestamp)),
                      "Bug: type of timestamp too big for TSlowControl"
                      );
    }
    TSlowControl() = default;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Type, Validity, Timestamp, Name, Description,
                Payload_Int, Payload_Float, Payload_String);
    }

    friend std::ostream& operator<<( std::ostream& s, const TSlowControl& o) {
        return s << "TSlowControl"
                 << " Type=" << o.TypeToString()
                 << " Timestamp='" << std_ext::to_iso8601(o.Timestamp) << "'"
                 << " Name='" << o.Name << "'"
                 << " Description='" << o.Description << "'";
    }

    const char* TypeToString() const {
        switch(Type) {
        case Type_t::AcquScaler:
            return "AcquScaler";
        case Type_t::EpicsOneShot:
            return "EpicsOneShot";
        case Type_t::EpicsScaler:
            return "EpicsScaler";
        case Type_t::EpicsTimer:
            return "EpicsTimer";
        case Type_t::Ant:
            return "Ant";
        }
        throw std::runtime_error("Not implemented");
    }

    const char* ValidityToString() const {
        switch(Validity) {
        case Validity_t::Backward:
            return "Backward";
        case Validity_t::Forward:
            return "Forward";
        }
        throw std::runtime_error("Not implemented");
    }
};

} // namespace ant
