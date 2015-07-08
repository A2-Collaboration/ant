#ifndef UNPACKERACQU_H
#define UNPACKERACQU_H

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include <memory>
#include <deque>
#include <vector>
#include <cstdint>

namespace ant {

class UnpackerAcquFileFormat; // implemented in detail/UnpackerAcqu.h

class UnpackerAcqu : public Unpacker::Module
{
public:
  UnpackerAcqu();
  virtual bool OpenFile(const std::string& filename) override;
  virtual std::shared_ptr<TDataRecord> NextItem() noexcept override;

  class Exception : public Unpacker::Exception {
    using Unpacker::Exception::Exception; // use base class constructor
  };

private:
  std::deque< std::unique_ptr<TDataRecord> > queue;
  std::unique_ptr<UnpackerAcquFileFormat> file;

};

// we define some methods here which
// we require the configs to implement
class UnpackerAcquConfig : public ExpConfig::Unpacker<UnpackerAcquConfig> {
public:
  struct bits_t {
    std::uint32_t RawChannel;
    std::uint32_t Mask;
  };

  using mapping_t = std::pair< std::vector<bits_t> , LogicalElement_t >;
  virtual void BuildMapping(std::vector<mapping_t>& mapping) = 0;
};

} // namespace ant

#endif // UNPACKERACQU_H
