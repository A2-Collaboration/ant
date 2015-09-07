#pragma once

#include "Rtypes.h"

#ifdef __CINT__

// simulate minimal cstdint for ROOTCINT
namespace std {
typedef UChar_t    uint8_t;
typedef UInt_t     uint32_t;
typedef ULong_t    uint64_t;
typedef Short_t    int16_t;
typedef Long_t     int64_t;
}

#else

#include "base/printable.h"
#include "base/std_ext/time.h"

#include <type_traits>
#include <tuple>
#include <cstdint>
#include <stdexcept>

static_assert(std::is_same<ULong_t, std::uint64_t>::value, "Type size mismatch");
static_assert(std::is_same<Long_t, std::int64_t>::value, "Type size mismatch");
#endif // __CINT__

#include <iomanip>
#include <ctime>

// if we want to change something
// of the data format defined by the ant::T* classes
#define ANT_UNPACKER_ROOT_VERSION 1

namespace ant {

#ifndef __CINT__
struct TID  : printable_traits
#else
struct TID
#endif
{

    std::uint32_t Flags;
    std::uint32_t Timestamp; // unix epoch
    std::uint32_t Lower;     // usually event counter
    std::uint32_t Reserved;  // unused

#ifndef __CINT__

    // you may append flags, but never remove or change order!
    enum class Flags_t : std::uint8_t {
        Invalid, MC
    };

    // ensure correct init in default constructor
    static_assert(static_cast<std::uint8_t>(Flags_t::Invalid) == 0, "Invalid flag should be first item in enum class");

    TID(
            std::uint32_t timestamp,
            std::uint32_t lower = 0,
            bool isMC = false)
        :
          Flags(0),
          Timestamp(timestamp),
          Lower(lower),
          Reserved(0)
    {
        Flags |= static_cast<decltype(Flags)>(isMC) << static_cast<std::uint8_t>(Flags_t::MC);
    }
    // prevent implicit conversion calls
    TID(std::uint32_t, bool) = delete;
    TID(std::uint32_t, std::uint32_t, int) = delete;

    bool IsInvalid() const {
        return Flags & (1 << static_cast<std::uint8_t>(Flags_t::Invalid));
    }

    std::uint64_t Value() const {
        return (static_cast<std::uint64_t>(Timestamp) << 8*sizeof(std::uint32_t)) + Lower;
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        if(IsInvalid())
            return s << "INVALID";

        s << "(";
        if(Flags)
            s  << "flags=0x" << Flags << ",";
        s << "'" << std_ext::to_iso8601(Timestamp) <<"',";
        s << "0x" << std::hex <<  std::setw(sizeof(decltype(Lower))*2) << std::setfill('0')
          << Lower << std::dec;
        s  << ")" ;
        return s;
    }

    bool operator<(const TID& rhs) const
    {
        if(IsInvalid() || rhs.IsInvalid())
            return false;
        auto make_tuple = [] (const TID& id) {
            return std::tie(id.Flags, id.Timestamp, id.Lower);
        };
        return make_tuple(*this) < make_tuple(rhs);
    }

    bool operator!=(const TID& rhs) const
    {
        if(IsInvalid() || rhs.IsInvalid())
            return true;
        return *this < rhs || rhs < *this;
    }

    bool operator<=(const TID& rhs) const
    {
        return *this < rhs || *this == rhs;
    }

    bool operator>=(const TID& rhs) const
    {
        return *this > rhs || *this == rhs;
    }

    bool operator>(const TID& rhs) const
    {
        return rhs < *this;
    }

    bool operator==(const TID& rhs) const
    {
        return !(*this != rhs);
    }

    TID& operator++() {
        ++Lower;
        return *this;
    }

    TID& operator--() {
        --Lower;
        return *this;
    }

#endif

    TID() : Flags(1), Timestamp(0), Lower(0), Reserved(0) {} // set the invalid flag by default
    virtual ~TID() {}
    ClassDef(TID, ANT_UNPACKER_ROOT_VERSION)

}; // TID

template<typename ValueType>
#ifndef __CINT__
struct TKeyValue  : printable_traits
#else
struct TKeyValue
#endif
{
    std::uint32_t  Key;
    ValueType Value;

#ifndef __CINT__
    TKeyValue(unsigned key, ValueType value) :
        Key(key), Value(std::move(value))
    {}
    virtual std::ostream& Print( std::ostream& s) const override {
        return s << Key << "=" << Value;
    }
#endif
    TKeyValue() : Key(), Value() {}
    virtual ~TKeyValue() {}
    // for some very strange reason,
    // ClassDef MUST NOT be used here
    // otherwise reading it from the tree gives a segfault!
    // Anyhow, reading/writing still works...
    // ClassDef(TKeyValue,0)
}; // TKeyValue<ValueType>

}
