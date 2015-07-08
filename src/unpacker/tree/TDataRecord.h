#ifndef TDATARECORD_H
#define TDATARECORD_H

#include <cstdint>
#include <base/printable.h>
#include <iomanip>

namespace ant {

struct TDataRecord : printable_traits
{
  struct ID_t : printable_traits {
    // you may append flags, but never remove or change order!
    enum class Flags_t : unsigned {
      MC
    };

    ID_t(
        uint32_t upper,
        uint32_t lower,
        bool isMC = false
        ) {
      Flags = 0;
      Flags |= static_cast<decltype(Flags)>(isMC) << static_cast<unsigned>(Flags_t::MC);
      Value = lower;
      Value |= static_cast<decltype(Value)>(upper) << sizeof(uint32_t)*8;
    }

    std::uint64_t Value;
    std::uint32_t Flags;

    virtual std::ostream& Print( std::ostream& s) const {
      return s << std::hex << "(flags=0x" << Flags << ",0x"
               << std::setw(sizeof(std::uint64_t)*2) << std::setfill('0')
               << Value
               << ")" << std::dec;
    }

  }; // UUID_t


  TDataRecord(const ID_t& id) : ID(id) {}
  virtual ~TDataRecord() = default;

  ID_t ID;

  virtual std::ostream& Print( std::ostream& s) const {
    return s << "TDataRecord ID=" << ID;
  }

};

} // namespace ant

#endif // TDATARECORD_H
