#ifndef UNPACKERACQU_H
#define UNPACKERACQU_H

#include "Unpacker.h"

namespace ant {

class UnpackerAcqu : public Unpacker::Module
{
public:
  UnpackerAcqu();
  virtual bool OpenFile(const std::string &filename);
  virtual std::shared_ptr<TDataRecord> NextItem();
private:
};

} // namespace ant

#endif // UNPACKERACQU_H
