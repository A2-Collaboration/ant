#ifndef TDATARECORD_H
#define TDATARECORD_H

#include <cstdint>

namespace ant {

struct TDataRecord
{
  struct ID_t {
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

    uint64_t Value;
    uint32_t Flags;
  }; // UUID_t


  TDataRecord(const ID_t& id) : ID(id) {}
  virtual ~TDataRecord() = default;

  ID_t ID;
};

} // namespace ant

#endif // TDATARECORD_H
