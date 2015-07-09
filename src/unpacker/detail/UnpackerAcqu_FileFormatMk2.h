#ifndef UNPACKERACQU_FILEFORMATMK2_H
#define UNPACKERACQU_FILEFORMATMK2_H

#include "UnpackerAcqu_detail.h"

namespace ant {
namespace unpacker {
namespace acqu {

class FileFormatMk2 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<std::uint32_t> &buffer) const override;
  virtual void FillInfo() override;
  virtual void FillFirstDataBuffer(queue_t& queue) override;
  virtual bool UnpackDataBuffer(queue_t &queue) noexcept override;

private:
  bool SearchFirstDataBuffer(queue_t &queue, size_t offset);
  using it_t = std::vector<uint32_t>::const_iterator;
  void HandleEPICSBuffer(queue_t& queue, it_t& it, const it_t& it_end, bool& good) const noexcept;
  void HandleScalerBuffer(queue_t& queue, it_t& it, const it_t& it_end, bool& good) const noexcept;
  void HandleReadError(queue_t& queue, it_t& it, const it_t& it_end, bool& good) const noexcept;
};

}}} // namespace ant::unpacker::acqu


#endif // UNPACKERACQU_FILEFORMATMK2_H
