
#include "UnpackerAcqu.h" // for UnpackerAcquConfig
#include "UnpackerAcqu_detail.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"
#include "UnpackerAcqu_FileFormatMk1.h"
#include "UnpackerAcqu_FileFormatMk2.h"


#include "tree/TDataRecord.h"
#include "tree/THeaderInfo.h"
#include "tree/TSlowControl.h"


#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "RawFileReader.h"

#include "base/std_ext.h"

#include <algorithm>
#include <exception>
#include <list>
#include <ctime>
#include <iterator> // for std::next

using namespace std;
using namespace ant;
using namespace ant::unpacker;


unique_ptr<UnpackerAcquFileFormat>
UnpackerAcquFileFormat::Get(const string &filename,
                            queue_t &queue)
{
  // make a list of all available acqu file format classes
  using format_t = unique_ptr<UnpackerAcquFileFormat>;
  list< format_t > formats;
  formats.push_back(std_ext::make_unique<acqu::FileFormatMk1>());
  formats.push_back(std_ext::make_unique<acqu::FileFormatMk2>());

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
  auto reader = std_ext::make_unique<RawFileReader>();
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

void acqu::FileFormatBase::Setup(reader_t &&reader_, buffer_t &&buffer_) {
  reader = move(reader_);
  buffer = move(buffer_);
}

void acqu::FileFormatBase::FillHeader(queue_t& queue)
{
  // let child class fill the info
  FillInfo(reader, buffer, info);

  auto headerInfo = BuildTHeaderInfo();

  // try to find some config with the headerInfo
  auto config = ExpConfig::Unpacker<UnpackerAcquConfig>::Get(*headerInfo);

  // then enqueue the header info
  fillQueue<THeaderInfo>(queue, move(headerInfo));

  // now try to fill the first data buffer
  FillFirstDataBuffer(queue, reader, buffer);
  unpackedBuffers = 0; // not yet unpacked

  // remember the record length size
  trueRecordLength = buffer.size();

  // get the mappings once
  config->BuildMappings(hit_mappings, scaler_mappings);
}

unique_ptr<THeaderInfo> acqu::FileFormatBase::BuildTHeaderInfo()
{
  // this unpacker has a constant ID_upper
  // based on the timestamp inside the file
  /// \todo make ID really unique due to daylight saving time...

  const time_t timestamp = mktime(&info.Time); // convert to unix epoch
  ID_upper = static_cast<decltype(ID_upper)>(timestamp);
  ID_lower = 0;

  // construct the unique ID, header record as lower ID=0
  const TID id(ID_upper, ID_lower);


  // build the genernal description
  stringstream description;
  description << "AcquData "
              << "Number=" << info.RunNumber << " "
              << "OutFile='" << info.OutFile << "' "
              << "Description='"+info.Description+"' "
              << "Note='"+info.RunNote+"' ";


  return std_ext::make_unique<THeaderInfo>(id, timestamp, description.str(), info.RunNumber);
}

void acqu::FileFormatBase::LogMessage(
    UnpackerAcquFileFormat::queue_t &queue,
    TUnpackerMessage::Level_t level,
    const string& msg
    ) const
{
  auto record = std_ext::make_unique<TUnpackerMessage>(
        TID(ID_upper, ID_lower),
        level,
        msg
        );

  const string& text = "[TUnpackerMessage] " + record->Message ;

  switch(level) {
  case TUnpackerMessage::Level_t::Info:
    LOG(INFO) << text;
    break;
  case TUnpackerMessage::Level_t::Warn:
    LOG(WARNING) << text;
    break;
  case TUnpackerMessage::Level_t::DataError:
  case TUnpackerMessage::Level_t::HardwareError:
  case TUnpackerMessage::Level_t::DataDiscard:
    LOG(ERROR) << text;
    break;
  }
  fillQueue(queue, move(record));
}


void acqu::FileFormatBase::FillEvents(queue_t& queue) noexcept
{
  // this method never throws exceptions, but just adds TUnpackerMessage to queue
  // if something strange while unpacking is encountered

  // we use the buffer as some state-variable
  // if the buffer is already empty now, there is nothing more to read
  if(buffer.empty())
    return;

  // start parsing the filled buffer
  // however, we fill a temporary queue first
  auto it = buffer.cbegin();
  queue_t queue_buffer;
  if(!UnpackDataBuffer(queue_buffer, it, buffer.cend())) {
    // handle errors on buffer scale
    VLOG(7) << "Error while unpacking buffer, discarding all unpacked data from buffer.";
    // add the last item in the queue (if any), which should be a TUnpackerMessage instance
    if(!queue_buffer.empty()) {
      const auto& lastItem = queue.back();
      auto ptr = dynamic_cast<const TUnpackerMessage*>(lastItem.get());
      if(ptr != nullptr) {
        queue.splice(queue.cend(), move(queue_buffer),
                     next(queue_buffer.cend(),-1), queue_buffer.cend());
      }
    }
    // always add an datadiscard info record
    auto record = std_ext::make_unique<TUnpackerMessage>(
          TID(ID_upper, ID_lower),
          TUnpackerMessage::Level_t::DataDiscard,
          "Discarded buffer number {}"
          );
    record->Payload.push_back(unpackedBuffers);
    fillQueue(queue, move(record));
  }
  else {
    // successful, so add all to output
    const int unpackedWords = distance(buffer.cbegin(), it);
    VLOG(7) << "Successfully unpacked " << unpackedWords << " words ("
            << 100.0*unpackedWords/buffer.size() << " %) from buffer ";
    queue.splice(queue.cend(), move(queue_buffer));
  }

  unpackedBuffers++;


  // refill the buffer
  try {
    reader->read(buffer.data(), trueRecordLength);
  }
  catch(ant::RawFileReader::Exception e) {
    // clear buffer if there was a problem when reading
    LogMessage(queue,
               TUnpackerMessage::Level_t::DataError,
               std_ext::formatter()
               << "Error while reading input: " << e.what());
    buffer.clear();
  }

  // check if actually enough bytes were read
  if(reader->gcount() != 4*trueRecordLength) {
    // reached the end of file properly?
    if(reader->gcount() == 0 && reader->eof()) {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::Info,
                 std_ext::formatter()
                 << "Found proper end of file");
    }
    else {
      LogMessage(queue,
                 TUnpackerMessage::Level_t::DataError,
                 std_ext::formatter()
                 << "Read only " << reader->gcount()
                 << " bytes, not enough for record length " << 4*trueRecordLength);
    }
    buffer.clear();
  }
}
