#ifndef TDATARECORD_H
#define TDATARECORD_H

#include <cstdint>

namespace ant {

struct TDataRecord
{
  struct UUID_t {
    UUID_t(
        uint32_t upper,
        uint32_t lower,
        bool isMC = false
        ) {
      flags = 0;
      flags |= isMC << 0;
      value = lower;
      value |= static_cast<decltype(value)>(upper) << sizeof(uint32_t)*8;
    }

  private:
    uint64_t value;
    uint32_t flags;
  }; // UUID_t

  TDataRecord(UUID_t UID) : UID_(UID) {}
  virtual ~TDataRecord() = default;
  UUID_t GetUID() const { return UID_; }
private:
  UUID_t UID_;
};

} // namespace ant

#endif // TDATARECORD_H
