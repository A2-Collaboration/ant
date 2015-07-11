#ifndef UNPACKERACQU_H
#define UNPACKERACQU_H

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include <memory>
#include <deque>
#include <vector>
#include <cstdint>
#include <limits>

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
    static_assert(std::is_unsigned<T>::value, "T must be unsigned");
    T RawChannel;
    T Mask;
    // provide some handy constructors (implemented below)
    // and constants
    RawChannel_t(const std::initializer_list<T>& l);
    RawChannel_t(const T& ch);
    static constexpr T NoMask = std::numeric_limits<T>::max();
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
  struct scaler_mapping_t {
    LogicalChannel_t LogicalElement;
    std::vector< RawChannel_t<std::uint32_t> > RawChannels;
  };

};

// define the templated constructors here to keep the class definition clean
template<typename T>
inline UnpackerAcquConfig::RawChannel_t<T>::RawChannel_t(const std::initializer_list<T> &l) {
  if(l.size()==2) {
    const std::vector<T> v(l);
    RawChannel = v[0];
    Mask = v[1];
  }
  else
    throw std::runtime_error("RawChannel_t can only be initialized with 2 values.");
}

template<typename T>
inline UnpackerAcquConfig::RawChannel_t<T>::RawChannel_t(const T &ch)
{
  RawChannel = ch;
  Mask = RawChannel_t<T>::NoMask;
}


} // namespace ant

#endif // UNPACKERACQU_H
