
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

UnpackerAcquFileFormat::~UnpackerAcquFileFormat() {}

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

    // and prepare the member variables for fast unpacking of hits
    for(const UnpackerAcquConfig::hit_mapping_t& hit_mapping : hit_mappings) {
        for(const UnpackerAcquConfig::RawChannel_t<uint16_t>& rawChannel : hit_mapping.RawChannels) {
            const uint16_t ch = rawChannel.RawChannel;
            if(hit_mappings_ptr.size()<=ch)
                hit_mappings_ptr.resize(ch+1);
            hit_mappings_ptr[ch].push_back(addressof(hit_mapping));
        }
    }
}

acqu::FileFormatBase::~FileFormatBase()
{

}

unique_ptr<THeaderInfo> acqu::FileFormatBase::BuildTHeaderInfo()
{
    // this unpacker has a constant ID_upper
    // based on the timestamp inside the file
    /// \todo make ID really unique due to daylight saving time...

    unsigned upperbits = GetUpperBitsTID() & 0xf;

    const time_t timestamp = mktime(&info.Time); // convert to unix epoch
    tm offset;
    strptime("Jan 1 00:00:00 2000", "%b %d %T %Y", &offset);
    offset.tm_isdst = 0;
    const time_t timestamp_offset = mktime(&offset);
    if(timestamp_offset > timestamp)
        throw UnpackerAcqu::Exception("File was recorded earlier than 2000");
    const time_t timebits = timestamp - timestamp_offset;
    if(timebits > 0x3fffffff) // 30bits maximum
        throw UnpackerAcqu::Exception("File was recorded later than ~2034");

    // construct the unique ID, header record as lower ID=0
    id = TID(static_cast<std::uint64_t>(timestamp) << sizeof(std::uint32_t));


    // build the genernal description
    stringstream description;
    description << "AcquData "
                << "Number=" << info.RunNumber << " "
                << "OutFile='" << info.OutFile << "' "
                << "Description='"+info.Description+"' "
                << "Note='"+info.RunNote+"' ";


    return std_ext::make_unique<THeaderInfo>(id, timestamp, description.str(), info.RunNumber);
}

unsigned acqu::FileFormatBase::GetUpperBitsTID()
{
    // there's no other way than hardcode the runs
    // recorded when the the MEST->MET transition occurred
    const vector<pair<unsigned, unsigned>> upper_nibble_mest2met = {
        {6592, 0x0},
        {6593, 0x1},
        {6594, 0x2},
        {6595, 0x3},
    };
    // handle run in transition
    if(std_ext::is_mest2met_transition(info.Time)) {
        auto it = upper_nibble_mest2met.begin();
        while(it != upper_nibble_mest2met.end()) {
            if(it->first == info.RunNumber) {
                return it->second;
            }
            ++it;
        }
        // not found...should be added than to list above
        if(it == upper_nibble_mest2met.end()) {
            throw UnpackerAcqu::Exception(
                        std_ext::formatter()
                        << "Cannot unpack file in MEST->MET transition "
                        << "without additional information for RunNumber " << info.RunNumber);
        }
    }

    // normal run

    for(size_t i=0;i<upper_nibble_mest2met.size();i++) {

    }

}

void acqu::FileFormatBase::LogMessage(
        UnpackerAcquFileFormat::queue_t &queue,
        TUnpackerMessage::Level_t level,
        const string& msg
        ) const
{
    auto record = std_ext::make_unique<TUnpackerMessage>(
                id,
                level,
                msg
                );

    stringstream ss_text;

    ss_text << "Buffer n=" << unpackedBuffers
            << " [TUnpackerMessage] " << record->Message ;

    const string& text = ss_text.str();

    switch(level) {
    case TUnpackerMessage::Level_t::Info:
        VLOG(3) << text;
        break;
    case TUnpackerMessage::Level_t::Warn:
        VLOG(2) << text;
        break;
    case TUnpackerMessage::Level_t::DataError:
    case TUnpackerMessage::Level_t::HardwareError:
    case TUnpackerMessage::Level_t::DataDiscard:
        VLOG(1) << text;
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
        LOG(WARNING) << "Error while unpacking buffer n=" << unpackedBuffers
                     << ", discarding all unpacked data from buffer.";
        // add the last item in the queue (if any), which should be a TUnpackerMessage instance
        if(!queue_buffer.empty()) {
            const auto& lastItem = queue.back();
            auto ptr = dynamic_cast<const TUnpackerMessage*>(lastItem.get());
            if(ptr != nullptr) {
                queue.splice(queue.end(), move(queue_buffer),
                             next(queue_buffer.end(),-1), queue_buffer.end());
            }
        }
        // always add an datadiscard info record
        auto record = std_ext::make_unique<TUnpackerMessage>(
                          id,
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
        queue.splice(queue.end(), move(queue_buffer));
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
