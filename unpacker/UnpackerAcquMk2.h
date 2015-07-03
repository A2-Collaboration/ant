#ifndef UNPACKERACQUMK2_H
#define UNPACKERACQUMK2_H

#include "Unpacker.h"

namespace ant {

class UnpackerAcquMk2 : public Unpacker::Module
{
public:
  UnpackerAcquMk2();
  virtual bool OpenFile(const std::string &filename);
  virtual std::shared_ptr<TDataRecord> NextItem();
};

} // namespace ant

#endif // UNPACKERACQUMK2_H
