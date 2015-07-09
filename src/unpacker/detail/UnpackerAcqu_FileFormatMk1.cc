#include "UnpackerAcqu_FileFormatMk1.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"

using namespace std;
using namespace ant;
using namespace ant::unpacker;

void acqu::FileFormatMk1::FillInfo()
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}

void acqu::FileFormatMk1::FillFirstDataBuffer(queue_t& queue)
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}

bool acqu::FileFormatMk1::UnpackDataBuffer(UnpackerAcquFileFormat::queue_t &queue) noexcept
{
  /// \todo Implement Mk1 unpacking
  return true;
}

size_t acqu::FileFormatMk1::SizeOfHeader() const
{
  return sizeof(AcquExptInfo_t);
}