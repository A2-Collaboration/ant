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
// the configs are required to implement
class UnpackerAcquConfig : public ExpConfig::Unpacker<UnpackerAcquConfig> {
public:

  struct RawChannel_t {
    std::uint32_t RawChannel;
    std::uint32_t Mask;
    // provide some handy constructors
    RawChannel_t(const std::initializer_list<uint32_t>& l);
    RawChannel_t(const uint32_t& ch);
  };

  // this defines how one LogicalChannel is built from
  // the given RawChannels
  struct mapping_t {
    LogicalChannel_t LogicalElement;
    std::vector<RawChannel_t> RawChannels;
  };

  virtual void BuildMappings(std::vector<mapping_t>& mappings) = 0;
};

} // namespace ant

#endif // UNPACKERACQU_H
