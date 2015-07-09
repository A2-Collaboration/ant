#include "UnpackerAcqu_FileFormatMk1.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

using namespace std;
using namespace ant;
using namespace ant::unpacker;

size_t acqu::FileFormatMk1::SizeOfHeader() const
{
  return sizeof(AcquExptInfo_t);
}

bool acqu::FileFormatMk1::InspectHeader(const vector<uint32_t>& buffer) const
{
  return inspectHeaderMk1Mk2<AcquExptInfo_t>(buffer);
}

void acqu::FileFormatMk1::FillInfo(reader_t &reader, buffer_t &buffer, Info &info) const
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}

void acqu::FileFormatMk1::FillFirstDataBuffer(queue_t& queue, reader_t &reader, buffer_t &buffer) const
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}

bool acqu::FileFormatMk1::UnpackDataBuffer(queue_t &queue, it_t& it, const it_t& it_endbuffer) noexcept
{
  /// \todo Implement Mk1 unpacking
  return true;
}

