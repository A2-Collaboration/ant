#include "UnpackerAcqu_FileFormatMk1.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

#include "RawFileReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"

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

void acqu::FileFormatMk1::FillInfo(reader_t& reader, buffer_t& buffer, Info& info)
{

    const acqu::AcquExptInfo_t* h = reinterpret_cast<const acqu::AcquExptInfo_t*>(buffer.data()+1);

    info.Format = Info::Format_t::Mk1;
    info.Time = std_ext::to_tm(h->fTime, "%a %b %d %T %Y");
    info.Description = std_ext::string_sanitize(h->fDescription);
    info.RunNote = std_ext::string_sanitize(h->fRunNote);
    info.OutFile = std_ext::string_sanitize(h->fOutFile);
    info.RunNumber = static_cast<unsigned>(h->fRun);
    info.RecordLength = static_cast<unsigned>(h->fRecLen);

    nScalers = h->fNscaler;


    /// \todo parse some more stuff from the Mk1 header here,
    /// but don't forget to read enough into the buffer using reader
    (void)reader; // prevent unused variable warning for now...

}

void acqu::FileFormatMk1::FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const
{
    // search at 0x8000 bytes
    if(SearchFirstDataBuffer(reader, buffer, 0x8000))
        return;

    // else fail
    throw UnpackerAcqu::Exception("Did not find first data buffer with Mk1 signature");
}

void acqu::FileFormatMk1::UnpackEvent(queue_t& queue,
                                      it_t& it, const it_t& it_endbuffer,
                                      bool& good) noexcept
{
    // Mk1 format does not tell us about the event length
    // so we scan before we try to unpack anything
    auto it_endevent = it;
    while(it_endevent != it_endbuffer &&
          *it_endevent != acqu::EEndEvent) {
        it_endevent++;
    }
    // check proper EEndEvent
    if(it == it_endbuffer) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   std_ext::formatter() <<
                   "While unpacking event, found premature end of buffer."
                   );
        return;
    }

    /// \todo Scan mappings if there's an ADC channel defined which mimicks those blocks
    queue.emplace_back(std_ext::make_unique<TEvent>(id));
    TEventData& eventdata = queue.back()->Reconstructed();
    eventdata.Trigger.DAQEventID = AcquID_last;

    hit_storage.clear();
    // there might be more than one scaler block in each event, so
    // so collect them first in this map
    scalers_t scalers;
    while(it != it_endevent) {
        // note that the Handle* methods move the iterator
        // themselves and set good to true if nothing went wrong

        good = false;

        switch(*it) {
        case acqu::EScalerBuffer:
            // Scaler read in this event
            HandleScalerBuffer(
                        scalers,
                        it, it_endevent, good);
            break;
        case acqu::EReadError:
            // read error block, some hardware-related information
            HandleDAQError(eventdata.Trigger.DAQErrors, it, it_endevent, good);
            break;
        default:
            // unfortunately, normal hits don't have a marker
            // so we hope for the best at this position in the buffer
            /// \todo Implement better handling of malformed event buffers
            static_assert(sizeof(acqu::AcquBlock_t) <= sizeof(decltype(*it)),
                          "acqu::AcquBlock_t does not fit into word of buffer");
            auto acqu_hit = reinterpret_cast<const acqu::AcquBlock_t*>(addressof(*it));
            // during a buffer, hits can come in any order,
            // and multiple hits with the same ID can happen
            hit_storage.add_item(acqu_hit->id, acqu_hit->adc);
            // decoding hits always works
            good = true;
            it++;
            break;
        }
        // stop immediately in case of problem
        /// \todo Implement more fine-grained unpacking error handling
        if(!good)
            return;
    }



    // hit_storage is member variable for better memory allocation performance
    FillDetectorReadHits(eventdata.DetectorReadHits);
    FillSlowControls(scalers, eventdata.SlowControls);

    it++; // go to start word of next event (if any)
}

void acqu::FileFormatMk1::HandleDAQError(vector<TDAQError>& errors,
                                         it_t& it, const it_t& it_end,
                                         bool& good) const noexcept
{
    // is there enough space in the event at all?
    static_assert(sizeof(acqu::ReadError_t) % sizeof(decltype(*it)) == 0,
                  "acqu::ReadError_t is not word aligned");
    constexpr int wordsize = sizeof(acqu::ReadError_t)/sizeof(decltype(*it));


    if(distance(it, it_end)<wordsize) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "acqu::ReadError_t block not completely present in buffer"
                   );
        return;
    }

    // then cast it to data structure
    const acqu::ReadError_t* err =
            reinterpret_cast<const acqu::ReadError_t*>(addressof(*it));

    // mapping is according to Acqu's Mk1ErrorCheck routine in TAcquRoot.h
    errors.emplace_back(err->fCrate, err->fBus, err->fCode, "Mk1-UNKNOWN");

    VLOG_N_TIMES(1000, 2) << errors.back();

    advance(it, wordsize);
    good = true;
}

void acqu::FileFormatMk1::HandleScalerBuffer(
        scalers_t& scalers,
        it_t& it, const it_t& it_end,
        bool& good
        ) const noexcept
{
    // ignore Scaler buffer marker
    it++;

    if(distance(it, it_end)<nScalers) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Mk1 scaler block not completely present in buffer"
                   );
        return;
    }

    /// \todo check if that scaler block decoding makes any sense
    /// we assume that a scaler block comes last in event
    /// also what the hell is Acqu doing with the SplitScaler stuff?!

    for(uint32_t scalerIndex = 0; scalerIndex < nScalers; scalerIndex++) {
        scalers[scalerIndex].push_back(*it);
        it++;
    }

    good = true;
}
