#ifndef UNPACKERACQU_H
#define UNPACKERACQU_H

#include "Unpacker.h"

#include <memory>
#include <deque>

namespace ant {

class UnpackerAcquFileFormat; // implemented in UnpackerAcqu_detail.h

class UnpackerAcqu : public Unpacker::Module
{
public:
  UnpackerAcqu();
  virtual bool OpenFile(const std::string &filename) override;
  virtual std::shared_ptr<TDataRecord> NextItem() override;

  class Exception : public Unpacker::Exception {
    using Unpacker::Exception::Exception; // use base class constructor
  };

private:
  std::deque< std::unique_ptr<TDataRecord> > queue;
  std::unique_ptr<UnpackerAcquFileFormat> file;

};



} // namespace ant

#endif // UNPACKERACQU_H
