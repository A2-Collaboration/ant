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


  template<typename T>
  struct RawChannel_t {
    T RawChannel;
    T Mask;
    // provide some handy constructors
    RawChannel_t(const std::initializer_list<T>& l);
    RawChannel_t(const T& ch);
  };

  // this defines how one LogicalChannel is built from
  // the given RawChannels
  struct mapping_t {
    LogicalChannel_t LogicalElement;
    std::vector< RawChannel_t<std::uint16_t> > RawChannels;
  };

  virtual void BuildMappings(std::vector<mapping_t>& mappings) = 0;

  // scalers in acqu can be handled as additional information
  // for a logical detector channel, or as TSlowControl items


};

} // namespace ant

#endif // UNPACKERACQU_H
