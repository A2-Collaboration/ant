#ifndef UNPACKERACQU_FILEFORMATMK1_H
#define UNPACKERACQU_FILEFORMATMK1_H

#include "UnpackerAcqu_detail.h"

namespace ant {
namespace unpacker {
namespace acqu {

class FileFormatMk1 : public FileFormatBase {

  // UnpackerAcquFile interface
protected:
  virtual size_t SizeOfHeader() const override;
  virtual bool InspectHeader(const std::vector<std::uint32_t>& buffer) const override;
  virtual void FillInfo() override;
  virtual void FillFirstDataBuffer(queue_t& queue) override;
  virtual bool UnpackDataBuffer(queue_t &queue) noexcept override;

};

}}} // namespace ant::unpacker::acqu


#endif // UNPACKERACQU_FILEFORMATMK1_H
