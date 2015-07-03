#ifndef TDATARECORD_H
#define TDATARECORD_H

#include <cstdint>

namespace ant {

struct TDataRecord
{
  TDataRecord(std::uint64_t UID) : UID_(UID) {}
  virtual ~TDataRecord() = default;
  std::uint64_t GetUID() const { return UID_; }
private:
  std::uint64_t UID_;
};

} // namespace ant

#endif // TDATARECORD_H
