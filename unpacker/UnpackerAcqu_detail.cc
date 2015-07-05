#include "UnpackerAcqu_detail.h"
#include "UnpackerAcqu_legacy.h"


#include "Logger.h"
#include "RawFileReader.h"

#include "stl_helpers.h"

#include <algorithm>
#include <exception>
#include <list>


using namespace std;
using namespace ant;
using namespace ant::unpacker;


unique_ptr<UnpackerAcquFileFormat> UnpackerAcquFileFormat::Get(const string &filename)
{
  // make a list of all available acqu file format classes
  using format_t = unique_ptr<UnpackerAcquFileFormat>;
  list< format_t > formats;
  formats.emplace_back(new acqu::FileFormatMk1());
  formats.emplace_back(new acqu::FileFormatMk2());

  // how many words shall we read for header inspection?
  const auto it_max =
      max_element(formats.cbegin(), formats.cend(),
                  [] (const format_t& f1, const format_t& f2)
  {
    return f1->SizeOfHeader() < f2->SizeOfHeader();
  });

  /// \todo Check for big-endian/little-endian formatted files?
  // convert bytes to words
  const size_t bufferSize = (*it_max)->SizeOfHeader()/sizeof(uint32_t) + 1;

  // now we try to open the file
  unique_ptr<RawFileReader> reader(new RawFileReader());
  reader->open(filename);
  vector<uint32_t> buffer(bufferSize);
  reader->read(buffer.data(), buffer.size());

  // then remove all formats which fail
  // to inspect the header in the given buffer
  formats.remove_if([&buffer] (const format_t& f) { return !f->InspectHeader(buffer); });

  if(formats.size() != 1) {
    LOG_IF(formats.size()>1, WARNING) << "More than one matching Acqu format found for file " << filename;
    return nullptr;
  }

  // we found only one candidate
  // give him the reader and the buffer for further processing
  const format_t& format = formats.back();
  format->Setup(move(reader), move(buffer));

  // return the fully setup UnpackerAcquFormat instance
  return move(formats.back());
}

size_t acqu::FileFormatMk1::SizeOfHeader() const
{
  return sizeof(AcquExptInfo_t);
}

size_t acqu::FileFormatMk2::SizeOfHeader() const
{
  return sizeof(AcquMk2Info_t);
}


template<typename T>
bool checkMk2(const T* h, true_type) {
  return h->fMk2 != acqu::EHeadBuff;
}

template<typename T>
bool checkMk2(const T*, false_type) {
  return false;
}

template<typename T>
bool inspectHeaderMk1Mk2(const vector<uint32_t>& buffer) {
  // ensure the correct T at compile time
  static_assert(is_same<T, acqu::AcquExptInfo_t>::value
                || is_same<T, acqu::AcquMk2Info_t>::value,
                "T can only be either Mk1 or Mk2 header struct");

  if(buffer[0] != acqu::EHeadBuff)
    return false;

  // try to interpret the buffer as some Mk2 header
  const T* h = reinterpret_cast<const T*>(buffer.data()+1);

  // Mk2 has some additional head marker in the struct,
  // but the other type does not have it
  // we use tag-based dispatching here:
  // if T is AcquMk2Info_t, then call the checkMk2 which inspects the struct,
  // if T is something else, then call the empty
  using tag = integral_constant<bool, is_same<T, acqu::AcquMk2Info_t>::value >;
  if(checkMk2(h, tag{}))
    return false;

  if(h->fRun<0)
    return false;

  if(std_ext::string_sanitize(h->fTime).length() != 24)
    return false;

  if(std_ext::string_sanitize(h->fOutFile).empty())
    return false;

  return true;
}

bool acqu::FileFormatMk1::InspectHeader(const vector<uint32_t>& buffer) const
{
  return inspectHeaderMk1Mk2<AcquExptInfo_t>(buffer);
}

bool acqu::FileFormatMk2::InspectHeader(const vector<uint32_t>& buffer) const
{
  return inspectHeaderMk1Mk2<AcquMk2Info_t>(buffer);
}


void acqu::FileFormatMk1::Setup(std::unique_ptr<RawFileReader> &&reader, std::vector<uint32_t> &&buffer)
{

}

void acqu::FileFormatMk2::Setup(std::unique_ptr<RawFileReader> &&reader, std::vector<uint32_t> &&buffer)
{

}
