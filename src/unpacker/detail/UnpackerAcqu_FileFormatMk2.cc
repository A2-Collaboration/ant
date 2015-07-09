#include "UnpackerAcqu_FileFormatMk2.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

#include "tree/THeaderInfo.h"
#include "tree/TUnpackerMessage.h"
#include "tree/TSlowControl.h"

#include "RawFileReader.h"

#include "base/std_ext.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::unpacker;

size_t acqu::FileFormatMk2::SizeOfHeader() const
{
  return sizeof(AcquMk2Info_t);
}

bool acqu::FileFormatMk2::InspectHeader(const vector<uint32_t>& buffer) const
{
  return inspectHeaderMk1Mk2<AcquMk2Info_t>(buffer);
}

void acqu::FileFormatMk2::FillInfo(reader_t &reader, buffer_t &buffer, Info &info) const
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
  static_assert(sizeof(acqu::AcquMk2Info_t) % sizeof(buffer_t::value_type) == 0,
                "AcquMk2Info_t is not aligned to buffer");
  static_assert(sizeof(acqu::ModuleInfoMk2_t) % sizeof(buffer_t::value_type) == 0,
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
}

void acqu::FileFormatMk2::FillFirstDataBuffer(queue_t& queue, reader_t& reader, buffer_t& buffer) const
{
  // finally search for the Mk2Header, this also
  // fills the buffer correctly with the first Mk2DataBuffer (if available)

  // first search at 0x8000 bytes
  if(SearchFirstDataBuffer(queue, reader, buffer, 0x8000))
    return;

  // then search at 10*0x8000 bytes
  if(SearchFirstDataBuffer(queue, reader, buffer, 10*0x8000))
    return;

  // else fail
  throw UnpackerAcqu::Exception("Did not find first data buffer with Mk2 signature");
}


bool acqu::FileFormatMk2::SearchFirstDataBuffer(queue_t& queue,
                                                reader_t& reader, buffer_t& buffer,
                                                size_t offset) const
{
  VLOG(9) << "Searching first Mk2 buffer at offset 0x"
          << hex << offset << dec;
  // offset is given in bytes, convert to number of words in uint32_t buffer
  const streamsize nWords = offset/sizeof(uint32_t);
  // read full header record, the file
  // should be at least this long
  reader->expand_buffer(buffer, nWords);
  // if this is a header-only file, the next expand might file
  try {
    reader->expand_buffer(buffer, nWords+1);
  }
  catch(RawFileReader::Exception e) {
    if(reader->eof()) {
      LOG(WARNING) << "File is exactly " << offset
                    << " bytes long, and contains only header.";
      // indicate eof by empty buffer
      buffer.clear();
      return true;
    }
    // else throw e
    throw e;
  }

  // this word should be the Mk2DataBuff...
  if(buffer.back() != acqu::EMk2DataBuff)
    return false;

  // check header info, emit message if there's a problem
  if(info.RecordLength != offset) {
    LogMessage(queue,
          TUnpackerMessage::Level_t::Warn,
          std_ext::formatter()
          << "Record length in header 0x" << hex << info.RecordLength
          << " does not match true file record length 0x" << offset << dec
          );
  }

  VLOG(9) << "Found first Mk2 buffer at offset 0x"
          << hex << offset << dec;

  // we finally prepare the first data buffer
  // buffer is at the moment nWords+1 large, and the last word
  // is the first word of the first data buffer
  buffer[0] = buffer.back();

  reader->read(&buffer[1], nWords-1);
  // ensure that nWords-1 are actually read
  streamsize expectedBytes = (nWords-1)*sizeof(uint32_t);
  if(reader->gcount() != expectedBytes) {
    throw UnpackerAcqu::Exception(
          std_ext::formatter()
          << "Only " << reader->gcount() << " bytes read from file, but "
          << expectedBytes << " required for first data buffer"
          );
  }
  // get rid of duplicate header word at very end
  // now the buffer.size() has exactly the length of one file record
  buffer.resize(nWords);

  return true;
}

bool acqu::FileFormatMk2::UnpackDataBuffer(UnpackerAcquFileFormat::queue_t& queue, it_t& it, const it_t& it_endbuffer) noexcept
{
  // check header word
  if(*it != acqu::EMk2DataBuff) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               std_ext::formatter() <<
               "Buffer starts with unexpected header word 0x" << hex << *it << dec
               );
    return false;
  }
  it++;

  // now loop over buffer contents
  while(*it != acqu::EBufferEnd && it != it_endbuffer) {

    // extract and check serial ID
    const unsigned acquID = *it;
    if(AcquID_last>acquID) {
//      LogMessage(queue,
//                 TUnpackerMessage::Level_t::Info,
//                 std_ext::formatter()
//                 << "Overflow of Acqu EventId detected from "
//                 << AcquID_last << " to " << acquID
//                 );
        VLOG(8) << "Overflow of Acqu EventId detected from "
                << AcquID_last << " to " << acquID;
    }
    AcquID_last = acquID;
    it++;

    // extract and check eventLength
    const unsigned eventLength = *it/sizeof(decltype(*it));

    const auto it_endevent = next(it, eventLength);
    if(it_endevent == it_endbuffer) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 std_ext::formatter() <<
                 "Event with size 0x" << eventLength
                 << " too big to fit in buffer of remaining size "
                 << distance(it, it_endbuffer)
                 );
      return false;
    }
    if(*it_endevent != acqu::EEndEvent) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 std_ext::formatter() <<
                 "At designated end of event, found unexpected word 0x"
                 << hex << *it_endevent << dec
                 );
      return false;
    }
    it++;

    // now work on one event inside buffer
    /// \todo Scan config if there's an ADC channel defined which mimicks those blocks

    bool good;
    while(*it != acqu::EEndEvent && it != it_endbuffer) {
      // note that the Handle* methods move the iterator
      // themselves and set good to true if nothing went wrong
      good = false;

      switch(*it) {
      case acqu::EEPICSBuffer:
        // EPICS buffer
        HandleEPICSBuffer(queue, it, it_endevent, good);
        break;
      case acqu::EScalerBuffer:
        // Scaler read in this event
        HandleScalerBuffer(queue, it, it_endevent, good);
        break;
      case acqu::EReadError:
        // read error block, some hardware-related information
        HandleReadError(queue, it, it_endevent, good);
        break;
      default:
        // unfortunately, normal hits don't have a marker
        // so we hope for the best at this position
        static_assert(sizeof(acqu::AcquBlock_t) <= sizeof(decltype(*it)),
                      "acqu::AcquBlock_t does not fit into word of buffer");
        const acqu::AcquBlock_t* acqu_hit = reinterpret_cast<const acqu::AcquBlock_t*>(addressof(*it));
        //VLOG(9) << "ADC ID=" << acqu_hit->id << " Value=" << acqu_hit->adc;

        good = true;

        it++;
        break;
      }
      if(!good)
        break;
    }

    // assume EEndEvent, skip it
    it++;

    // increment official unique event ID
    ID_lower++;
  }


  return true;
}


void acqu::FileFormatMk2::HandleEPICSBuffer(
    UnpackerAcquFileFormat::queue_t &queue,
    it_t& it, const it_t& it_end,
    bool& good
    ) const noexcept
{
  // ignore EPICS buffer marker
  it++;

  // is there enough space in the event at all?
  constexpr size_t wordbytes = sizeof(decltype(*it));
  static_assert(sizeof(acqu::EpicsHeaderInfo_t) % wordbytes == 0,
                "acqu::EpicsHeaderInfo_t is not word aligned");
  constexpr int headerWordsize = sizeof(acqu::EpicsHeaderInfo_t)/wordbytes;

  if(distance(it, it_end)<headerWordsize) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "acqu::EpicsHeaderInfo_t header not completely present in buffer"
               );
    return;
  }

  // then cast it to data structure
  const acqu::EpicsHeaderInfo_t* hdr =
      reinterpret_cast<const acqu::EpicsHeaderInfo_t*>(addressof(*it));

  // check the given header info
  // hdr->len aka epicsTotalWords
  // is the maximum size of the EPICS data including the info header
  if(hdr->len % wordbytes != 0) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "EPICS data not word aligned"
               );
    return;
  }
  const int epicsTotalWords = hdr->len/wordbytes;
  const string epicsModName = hdr->name;
  const size_t nChannels = hdr->nchan;
  const time_t hdr_timestamp = hdr->time;
  //const string epicsTime = std_ext::ctime(hdr->time);

  if(epicsModName.length()>32) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "acqu::EpicsHeaderInfo_t header has malformed module name"
               );
    return;
  }
  if(distance(it, it_end)<epicsTotalWords) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "EPICS data not completely present in buffer"
               );
    return;
  }

  // unfortunately, the EPICS data is no longer word-aligned,
  // so create some additional uint8_t buffer here
  /// \todo Find some way without copying the data?
  /// \todo correct for endian-ness / machine byte ordering?
  const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(addressof(*it));
  const vector<uint8_t> bytes(byte_ptr, byte_ptr + epicsTotalWords*wordbytes);

  auto it_byte = bytes.cbegin();

  // then skip the epics header info
  advance(it_byte, headerWordsize*wordbytes);

  // the header told us how many EPICS channels there are,
  // so start decoding them

  for(size_t i=0; i<nChannels; i++) {
    constexpr int chHdrBytes = sizeof(acqu::EpicsChannelInfo_t);
    if(distance(it_byte, bytes.cend()) < chHdrBytes) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 "EPICS channel header not completely present in buffer"
                 );
      return;
    }

    const acqu::EpicsChannelInfo_t* ch =
        reinterpret_cast<const acqu::EpicsChannelInfo_t*>(addressof(*it_byte));

    if(distance(it_byte, bytes.cend()) < ch->bytes) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 "EPICS channel payload not completely present in buffer"
                 );
      return;
    }

    auto it_map = map_EpicsTypes.find(ch->type);
    if(it_map == map_EpicsTypes.cend()) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 "EPICS channel type unknown"
                 );
      return;
    }

    const acqu::EpicsDataTypes_t ch_datatype = it_map->second.first;
    const int16_t ch_typesize = it_map->second.second;
    const int16_t ch_nElements = ch->nelem;
    const string  ch_Name = ch->pvname;

    // another size check for the channel payload
    if(ch->bytes != chHdrBytes + ch_nElements * ch_typesize) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 "EPICS channel payload size inconsistent"
                 );
      return;
    }

    // finally we can create the TSlowControl record

    TSlowControl::Type_t record_type = TSlowControl::Type_t::EpicsOneShot;
    stringstream description;
    if(hdr->period<0) {
      record_type = TSlowControl::Type_t::EpicsTimer;
      description << "Period='" << -hdr->period << " ms'";
    }
    else if(hdr->period>0) {
      record_type = TSlowControl::Type_t::EpicsScaler;
      description << "Period='" << hdr->period << " events'";
    }

    auto record = createDataRecord<TSlowControl>(
          TDataRecord::ID_t(ID_upper, ID_lower),
          record_type,
          hdr_timestamp,
          ch_Name,
          description.str()
          );

    // advance to the EPICS channel data (skip channel info header)
    advance(it_byte, chHdrBytes);

    // fill the payload depending on the EPICS data type
    // upcast float to double and short,byte to long

    for(int16_t elem=0;elem<ch_nElements;elem++) {
      switch(ch_datatype) {
      case acqu::EpicsDataTypes_t::BYTE:
        record->Payload_Int.push_back(*it_byte);
        break;
      case acqu::EpicsDataTypes_t::SHORT: {
        const int16_t* value = reinterpret_cast<const int16_t*>(addressof(*it_byte));
        record->Payload_Int.push_back(*value);
        break;
      }
      case acqu::EpicsDataTypes_t::LONG: {
        const int64_t* value = reinterpret_cast<const int64_t*>(addressof(*it_byte));
        record->Payload_Int.push_back(*value);
        break;
      }
      case acqu::EpicsDataTypes_t::FLOAT: {
        static_assert(sizeof(float)==4,"Float should be 4 bytes long");
        const float* value = reinterpret_cast<const float*>(addressof(*it_byte));
        record->Payload_Float.push_back(*value);
        break;
      }
      case acqu::EpicsDataTypes_t::DOUBLE: {
        static_assert(sizeof(double)==8,"Float should be 8 bytes long");
        const double* value = reinterpret_cast<const double*>(addressof(*it_byte));
        record->Payload_Float.push_back(*value);
        break;
      }
      case acqu::EpicsDataTypes_t::STRING: {
        const char* value = reinterpret_cast<const char*>(addressof(*it_byte));
        // interpret as string
        const string value_str(value);
        if((signed)value_str.length()>=ch_typesize) {
          LogMessage(queue,
                     TUnpackerMessage::Level_t::DataError,
                     "EPICS channel string data too long (no terminating \\0?)"
                     );
          return;
        }
        record->Payload_String.push_back(value);
        break;
      }
      default:
        throw UnpackerAcqu::Exception("Not implemented");

      } // end switch

      advance(it_byte, ch_typesize);
    }

    VLOG(9) << *record;

    // enqueue the nicely created EPICS slowcontrol record
    fillQueue(queue, move(record));

  } // end channel loop

  // we successfully parsed the EPICS buffer
  advance(it, epicsTotalWords);

  VLOG(9) << "Successfully parsed EPICS buffer";

  good = true;
}

void acqu::FileFormatMk2::HandleScalerBuffer(UnpackerAcquFileFormat::queue_t &queue,
    it_t& it, const it_t& it_end,
    bool& good) const noexcept
{
  throw UnpackerAcqu::Exception("Not implemented");

  good = true;
}

void acqu::FileFormatMk2::HandleReadError(
    UnpackerAcquFileFormat::queue_t &queue,
    it_t& it, const it_t& it_end,
    bool& good) const noexcept
{
  // is there enough space in the event at all?
  static_assert(sizeof(acqu::ReadErrorMk2_t) % sizeof(decltype(*it)) == 0,
                "acqu::ReadErrorMk2_t is not word aligned");
  constexpr int wordsize = sizeof(acqu::ReadErrorMk2_t)/sizeof(decltype(*it));

  if(distance(it, it_end)<wordsize) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "acqu::ReadErrorMk2_t block not completely present in buffer"
               );
    return;
  }

  // then cast it to data structure
  const acqu::ReadErrorMk2_t* err =
      reinterpret_cast<const acqu::ReadErrorMk2_t*>(addressof(*it));

  // check for trailer word
  if(err->fTrailer != acqu::EReadError) {
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               "Acqu ErrorBlock does not end with expected trailer word"
               );
    return;
  }

  // lookup the module name
  auto it_modname = acqu::ModuleIDToString.find(err->fModID);
  const string& modname = it_modname == acqu::ModuleIDToString.cend()
      ? "UNKNOWN" : it_modname->second;

  // build TUnpackerMessage record from error info
  auto record = createDataRecord<TUnpackerMessage>(
        TDataRecord::ID_t(ID_upper, ID_lower),
        TUnpackerMessage::Level_t::HardwareError,
        std_ext::formatter()
        << "Acqu HardwareError ModuleID={} (" << modname << ") "
        << "Index={} ErrorCode={}"
        );
  record->Payload.push_back(err->fModID);
  record->Payload.push_back(err->fModIndex);
  record->Payload.push_back(err->fErrCode);

  VLOG(9) << *record;

  fillQueue(queue, move(record));
  advance(it, wordsize);
  good = true;
}

