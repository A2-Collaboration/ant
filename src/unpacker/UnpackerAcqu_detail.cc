#include "UnpackerAcqu_detail.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu.h"

#include "TDataRecord.h"
#include "THeaderInfo.h"

#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "RawFileReader.h"

#include "base/stl_helpers.h"

#include <algorithm>
#include <exception>
#include <list>
#include <ctime>

using namespace std;
using namespace ant;
using namespace ant::unpacker;


unique_ptr<UnpackerAcquFileFormat>
UnpackerAcquFileFormat::Get(const string &filename,
                            deque<unique_ptr<TDataRecord> >& queue)
{
  // make a list of all available acqu file format classes
  using format_t = unique_ptr<UnpackerAcquFileFormat>;
  list< format_t > formats;
  formats.emplace_back(new acqu::FileFormatMk1());
  formats.emplace_back(new acqu::FileFormatMk2());

  // the rather complicated setup of FileFormat implementation
  // is due to the fact that we want to open the file only once,
  // but cannot seek in it (due to compressed files)

  // that means that this factory method needs to know something about

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
  // also fill some header-like events into the queue
  const format_t& format = formats.back();
  format->Setup(move(reader), move(buffer));
  format->FillHeader(queue);

  // return the UnpackerAcquFormat instance
  return move(formats.back());
}

void acqu::FileFormatBase::Setup(std::unique_ptr<RawFileReader> &&reader_, std::vector<uint32_t> &&buffer_) {
  reader = move(reader_);
  buffer = move(buffer_);
}

void acqu::FileFormatMk1::FillEvents(std::deque<std::unique_ptr<TDataRecord> >& queue)
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}

void acqu::FileFormatMk1::FillInfo()
{
  throw UnpackerAcqu::Exception("Mk1 format not implemented yet");
}


void acqu::FileFormatMk2::FillEvents(std::deque<std::unique_ptr<TDataRecord> >& queue)
{

}


void acqu::FileFormatBase::FillHeader(std::deque<std::unique_ptr<TDataRecord> >& queue)
{
  FillInfo();




  auto headerInfo = BuildTHeaderInfo();
  // try to find some config with the headerInfo
  config = ExpConfig::Unpacker<UnpackerAcquConfig>::Get(*headerInfo);
  // then enqueue the header info
  // but upcast the pointer for this
  queue.emplace_back(
        std_ext::static_cast_uptr<THeaderInfo, TDataRecord>(move(headerInfo))
        );

}



unique_ptr<THeaderInfo> acqu::FileFormatBase::BuildTHeaderInfo()
{
  // this unpacker has a constant ID_upper
  // based on the timestamp inside the file
  /// \todo make 100% unique due to daylight saving time,

  const time_t timestamp = mktime(&info.Time); // convert to unix epoch
  ID_upper = static_cast<decltype(ID_upper)>(timestamp);

  // construct the unique ID, header record as lower ID=0
  const TDataRecord::ID_t id(ID_upper, 0);


  // build the genernal description
  stringstream description;
  description << "AcquData "
              << "Number=" << info.RunNumber << " "
              << "OutFile='" << info.OutFile << "' "
              << "Description='"+info.Description+"' "
              << "Note='"+info.RunNote+"' ";


  return unique_ptr<THeaderInfo>(
        new THeaderInfo(id, timestamp, description.str(), info.RunNumber)
        );
}

void acqu::FileFormatMk2::FillInfo()
{
  const acqu::AcquMk2Info_t* h = reinterpret_cast<const acqu::AcquMk2Info_t*>(buffer.data()+1);

  // fill the Info

  // the C++11 istringstream way is not working on older systems
  // since std::get_time is missing :(
  //istringstream ss_time(h->fTime);
  //ss_time >> get_time(&info.Time,"%a %b %d %T %Y");
  strptime(h->fTime, "%a %b %d %T %Y", &info.Time);
  info.Time.tm_isdst = 0; // std::get_time does not set this, but mktime wants to know


  info.Description = std_ext::string_sanitize(h->fDescription);
  info.RunNote = std_ext::string_sanitize(h->fRunNote);
  info.OutFile = std_ext::string_sanitize(h->fOutFile);
  info.RunNumber = static_cast<unsigned>(h->fRun);
  info.RecordLength = static_cast<unsigned>(h->fRecLen);

  // extract the module info
  const unsigned nModules = static_cast<unsigned>(h->fNModule);
  VLOG(9) << "Header says: Have " << nModules << " modules";

  // calculate sizes, ensure that on this platform the extraction works
  using buffer_type = decltype(buffer);
  static_assert(sizeof(acqu::AcquMk2Info_t) % sizeof(buffer_type::value_type) == 0,
                "AcquMk2Info_t is not aligned to buffer");
  static_assert(sizeof(acqu::ModuleInfoMk2_t) % sizeof(buffer_type::value_type) == 0,
                "ModuleInfoMk2_t is not aligned to buffer");
  constexpr size_t infoSize = 1 + // dont't forget the 32bit start marker
      sizeof(acqu::AcquMk2Info_t)/4; // division by four is ok since static_assert'ed
  constexpr size_t moduleSize = sizeof(acqu::ModuleInfoMk2_t)/4;
  const size_t totalSize =
      infoSize +
      nModules*moduleSize;

  // read enough into buffer
  reader->expand_buffer(buffer, totalSize);

  // add modules to lists
  /// \todo Check for overlapping of raw channels?

  unsigned totalADCs = 0;
  unsigned totalScalers = 0;
  for(unsigned i=0;i<nModules;i++) {
    const auto buffer_ptr = &buffer[infoSize+i*moduleSize];
    const acqu::ModuleInfoMk2_t* m =
        reinterpret_cast<const acqu::ModuleInfoMk2_t*>(buffer_ptr);
    auto it = ModuleIDToString.find(m->fModID);
    if(it == ModuleIDToString.end()) {
      LOG(WARNING) << "Skipping unknown module with ID=0x" << hex << m->fModID << dec;
      continue;
    }

    Info::HardwareModule module;
    module.Identifier = it->second;
    module.Index = m->fModIndex;
    module.Bits = m->fBits;
    module.FirstRawChannel = m->fAmin;

    // ADC and Scaler will be added twice!
    if(m->fModType & acqu::EDAQ_ADC) {
      totalADCs += m->fNChannel;
      module.NRawChannels = m->fNChannel;
      info.ADCModules.emplace_back(module);
    }
    if(m->fModType & acqu::EDAQ_Scaler) {
      totalScalers += m->fNScChannel;
      module.NRawChannels = m->fNScChannel;
      info.ScalerModules.emplace_back(module);
    }
  }
  VLOG(9) << "Header says: Have " << info.ADCModules.size() << " ADC modules with "
          << totalADCs << " channels";
  VLOG(9) << "Header says: Have " << info.ScalerModules.size() << " Scaler modules with "
          << totalScalers << " channels";

  // finally search for the Mk2Header

  // first search at 0x8000 bytes
  if(SearchMk2Signature(0x8000))
    return;

  // then search at 10*0x8000
  if(SearchMk2Signature(10*0x8000))
    return;

  // else fail
  throw UnpackerAcqu::Exception("Did not find first data buffer with Mk2 signature");
}

bool acqu::FileFormatMk2::SearchMk2Signature(size_t offset)
{
  VLOG(9) << "Searching first Mk2 buffer at offset 0x" << hex << offset << dec;
  // read full header record, and one additional word
  reader->expand_buffer(buffer, offset/4+1);
  // this word should be the Mk2DataBuff...
  if(buffer.back() != acqu::EMk2DataBuff)
    return false;
  // check header info
  LOG_IF(info.RecordLength != offset, WARNING)
      << "Record length in header 0x" << hex << info.RecordLength
      << " does not match true file record length 0x" << offset << dec;
  // overwrite the RecordLength
  info.RecordLength = offset;
  return true;
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
bool inspectHeaderMk1Mk2_(const T* h, true_type) {
  return h->fMk2 != acqu::EHeadBuff;
}

template<typename T>
bool inspectHeaderMk1Mk2_(const T*, false_type) {
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
  if(inspectHeaderMk1Mk2_(h, tag{}))
    return false;

  /// \todo Improve the checking here

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
