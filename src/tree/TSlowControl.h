#pragma once

#include "TID.h" // for TKeyValue
#include "base/std_ext/time.h"
#include "base/printable.h"

#include <cstdint>
#include <tuple>

namespace ant {

struct TSlowControl : printable_traits
{
    enum class Type_t : std::uint8_t {
        AcquScaler, EpicsOneShot, EpicsScaler, EpicsTimer
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
    std::int64_t Timestamp;   // unix epoch, only meaningful for Type==Epics* items
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
    TSlowControl() {}
    virtual ~TSlowControl() {}


    template<class Archive>
    void serialize(Archive& archive) {
        archive(Type, Validity, Timestamp, Name, Description,
                Payload_Int, Payload_Float, Payload_String);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TSlowControl"
                 << " Type=" << TypeToString()
                 << " Timestamp='" << std_ext::to_iso8601(Timestamp) << "'"
                 << " Name='" << Name << "'"
                 << " Description='" << Description << "'";
    }

    const char* TypeToString() const {
        return type_to_string(Type);
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

    static const char* type_to_string(Type_t type) {
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

    /**
     * @brief The Key struct represents a TSlowControl item without payload
     */
    struct Key  : printable_traits
    {
        Type_t Type;
        std::string Name;
        Key(Type_t type, const std::string& name) :
            Type(type),
            Name(name)
        {}
        Key(const TSlowControl& sc) : Key(sc.Type, sc.Name) {}
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

};

} // namespace ant
