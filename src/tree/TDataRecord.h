#pragma once

#include "Rtypes.h"

#ifdef __CINT__
// simulate minimal cstdint for ROOTCINT
namespace std {
typedef UChar_t    uint8_t;
typedef UInt_t     uint32_t;
typedef ULong_t  uint64_t;
typedef Short_t    int16_t;
typedef Long_t   int64_t;
}
#else
#include <type_traits>
static_assert(std::is_same<ULong_t, std::uint64_t>::value, "Type size mismatch");
static_assert(std::is_same<Long_t, std::int64_t>::value, "Type size mismatch");
#include <cstdint>
#include <stdexcept>
#include "base/printable.h"
#endif // __CINT__

#include <string>
#include <iomanip>
#include <ctime>


// if we want to change something
// of the data format defined by the ant::T* classes
#define ANT_UNPACKER_ROOT_VERSION 1

namespace ant {

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

#ifndef __CINT__
struct TID  : printable_traits
#else
struct TID
#endif
{

    TID() : Value(0), Flags(1) {} // set the invalid flag by default
    virtual ~TID() {}

    std::uint64_t Value;
    std::uint32_t Flags;

#ifndef __CINT__

    // you may append flags, but never remove or change order!
    enum class Flags_t : std::uint8_t {
        Invalid, MC, Acqu, GoAT
    };

    // ensure correct init in default constructor
    static_assert(static_cast<std::uint8_t>(Flags_t::Invalid) == 0, "Invalid flag should be first item in enum class");

    TID(
            UInt_t upper,
            UInt_t lower,
            bool isMC = false
            )
        :
          Value(lower),
          Flags(0)
    {
        Flags |= static_cast<decltype(Flags)>(isMC) << static_cast<std::uint8_t>(Flags_t::MC);
        Value |= static_cast<decltype(Value)>(upper) << sizeof(std::uint32_t)*8;
    }

    bool IsInvalid() const {
        return Flags & (1 << static_cast<std::uint8_t>(Flags_t::Invalid));
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        if(IsInvalid())
            return s << "INVALID";

        return s << std::hex << "(flags=0x" << Flags << ",0x"
                 << std::setw(sizeof(decltype(Value))*2) << std::setfill('0')
                 << Value
                 << ")" << std::dec;
    }
#endif

    bool operator<(const TID& other) const
    {
        if(IsInvalid())
            return false;
        if (Flags == other.Flags)
            return (Value < other.Value);
        return (Flags < other.Flags);
    }

    bool operator!=(const TID& other) const
    {
        if(IsInvalid())
            return true;
        return *this < other || other < *this;
    }

    bool operator<=(const TID& other) const
    {
        return *this < other || *this == other;
    }

    bool operator>=(const TID& other) const
    {
        return *this > other || *this == other;
    }

    bool operator>(const TID& other) const
    {
        return other < *this;
    }

    bool operator==(const TID& other) const
    {
        return !(*this != other);
    }

    TID& operator++() {
        ++Value;
        return *this;
    }

    TID& operator--() {
        --Value;
        return *this;
    }

    ClassDef(TID, ANT_UNPACKER_ROOT_VERSION)

}; // TID

#ifndef __CINT__
struct TDataRecord : printable_traits
#else
struct TDataRecord
#endif
{
    TDataRecord() : ID() {}
    TDataRecord(const TID& id) : ID(id) {}
    virtual ~TDataRecord() {}

    TID ID;

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TDataRecord ID=" << ID;
    }
#endif

    ClassDef(TDataRecord, ANT_UNPACKER_ROOT_VERSION)

}; // TDataRecord


} // namespace ant
#pragma once
