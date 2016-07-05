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

bool acqu::FileFormatMk1::UnpackDataBuffer(UnpackerAcquFileFormat::queue_t& queue, it_t& it,
                                           const it_t& it_endbuffer) noexcept
{
    // at least 4 words should be in buffer...
    if(distance(it, it_endbuffer)<4)
        return false;

    // check header word
    if(*it != GetDataBufferMarker()) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   std_ext::formatter() <<
                   "Buffer starts with unexpected header word 0x" << hex << *it << dec
                   );
        return false;
    }
    it++;

    // scan for EndEvent markers in between
    std::vector<it_t> endevent_markers;
    {
        // this degenerate choice of markers makes decoding kinda complicated
        // additionally, scaler reads can mimick this marker...sigh
        assert(acqu::EBufferEnd == acqu::EEndEvent);

        // so we search for places were we find a possible EEndEvent candidate
        // which must be followed by the current event number
        // the very first one is given by first word in the data buffer after the marker
        auto event_id = *it;
        for(auto it_ = it; it_ != std::prev(it_endbuffer); ++it_) {
            auto it_next = std::next(it_);

            if(*it_ == acqu::EEndEvent &&
               *it_next == event_id+1) {
                endevent_markers.emplace_back(it_);
                ++event_id;
            }
        }

        if(endevent_markers.empty()) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       "No EEndEvent marker candidates found"
                       );
            return false;
        }

        // the last event marker is a bit harder to find
        // we assume that there are no super-short events with event nr only for example
        if(std::distance(endevent_markers.back(), it_endbuffer)<3) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       std_ext::formatter() <<
                       "Last but one event is too close to end of buffer"
                       );
            return false;
        }


        bool found_threetimes0xffffffff = false;
        // either its end is marked as the first appearance of three times 0xffffffff
        for(auto it_ = endevent_markers.back(); it_ != std::prev(it_endbuffer,2); ++it_) {
            auto it_next = std::next(it_);
            auto it_nextnext = std::next(it_next);
            if(*it_ == acqu::EEndEvent &&
               *it_next == acqu::EBufferEnd &&
               *it_nextnext == acqu::EBufferEnd) {
                endevent_markers.emplace_back(it_);
                found_threetimes0xffffffff = true;
                // important to stop after the first occurence since rest of buffer is filled
                // with garbage from previous buffers (acqu does not clear buffers....sigh)
                break;
            }
        }

        // or the last event just fit into the buffer, then only one (or two?) 0xffffffff
        // at the very end are there
        if(!found_threetimes0xffffffff) {
            auto it_prevprev = std::prev(it_endbuffer, 2);
            auto it_prev = std::prev(it_endbuffer, 1);

            // prefer the second last one as end of event marker
            if(*it_prevprev == acqu::EEndEvent) {
                endevent_markers.emplace_back(it_prevprev);
            }
            else if(*it_prev == acqu::EEndEvent) {
                endevent_markers.emplace_back(it_prev);
            }
            else {
                // unexpected end of event, stop unpacking...
                LogMessage(TUnpackerMessage::Level_t::DataError,
                           std_ext::formatter() <<
                           "Last event did not have event buffer marker as last word but 0x" << hex << *it_prev
                           );
                return false;
            }

            // notify about this special case
            LogMessage(TUnpackerMessage::Level_t::Info,
                       std_ext::formatter()
                       << "Buffer was exactly filled with " << endevent_markers.size()
                       << " events, no buffer endmarker present"
                       );
        }
    }


    for(auto it_end : endevent_markers) {
        // if buffers are discarded in between, this can happen...
        if(id.Lower != *it) {
            VLOG(9) << "Ant TID out-of-sync with Acqu event id " << *it;
        }

        bool good = false;
        UnpackEvent(queue, it, it_end, good);
        if(!good)
            return false;

        it = std::next(it_end);
        ++id;
    }

//    cout << endl;
//    for(auto offset : endevent_markers) {
//        cout << setfill('0') << setw(4) << std::distance(it, offset) << " ";
//    }

//    cout << endl;

//    // now loop over buffer contents, aka single events
//    unsigned nEventsInBuffer = 0;
//    while(it != it_endbuffer && *it != acqu::EBufferEnd) {

//        // extract and check serial ID
//        const unsigned acquID = *it;
//        if(AcquID_last>acquID) {
//            VLOG(8) << "Overflow of Acqu EventId detected from "
//                    << AcquID_last << " to " << acquID;
//        }
//        AcquID_last = acquID;
//        it++;

//        bool good = false;
//        UnpackEvent(queue, it, it_endbuffer, good);
//        if(!good)
//            return false;
//        // append the messages to some successfully unpacked event
//        AppendMessagesToEvent(queue.back());

//        // increment official unique event ID
//        ++id;
//        nEventsInBuffer++;
//    }

//    // we reached the end of buffer before finding the acqu::EBufferEnd
//    if(it == it_endbuffer) {
//        // there's one exception when the sum of the events
//        // inside one buffer fit exactly into the buffer,
//        // then only the end marker for the event is present
//        // but there's no way to detect this properly, since
//        // acqu::EEndEvent == acqu::EBufferEnd (grrrrr)
//        if(*next(it,-1) == acqu::EEndEvent) {
//            LogMessage(TUnpackerMessage::Level_t::Info,
//                       std_ext::formatter()
//                       << "Buffer was exactly filled with " << nEventsInBuffer
//                       << " events, no buffer endmarker present"
//                       );
//            return true;
//        }

//        LogMessage(TUnpackerMessage::Level_t::DataError,
//                   std_ext::formatter()
//                   << "Buffer did not have proper end buffer marker:"
//                   << "  1. lastword=0x" << hex << setw(8) << setfill('0') << *next(it,-1)
//                   << ", 2. lastword=0x" << hex << setw(8) << setfill('0') << *next(it,-2)
//                   << ", buffersize_bytes=0x" << buffersize_bytes
//                   );
//        return false;
//    }

//    queue.emplace_back(id);


    return true;
}

void acqu::FileFormatMk1::UnpackEvent(queue_t& queue, it_t it, const it_t& it_end, bool& good) noexcept
{


    /// \todo Scan mappings if there's an ADC channel defined which mimicks those blocks
    queue.emplace_back(id);
    TEventData& eventdata = queue.back().Reconstructed();

    // expect the first word to be the event id
    eventdata.Trigger.DAQEventID = *it;
    ++it;

    hit_storage.clear();
    // there might be more than one scaler block in each event, so
    // so collect them first in this map
    scalers_t scalers;
    while(it != it_end) {
        // note that the Handle* methods move the iterator
        // themselves and set good to true if nothing went wrong

        good = false;

        switch(*it) {
        case acqu::EScalerBuffer:
            // Scaler read in this event
            HandleScalerBuffer(
                        scalers,
                        it, it_end, good);
            break;
        case acqu::EReadError:
            // read error block, some hardware-related information
            HandleDAQError(eventdata.Trigger.DAQErrors, it, it_end, good);
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

    cout << endl;
    cout << "Scaler block at TID=" << id << endl;

    auto it_endscaler = it_end;
//    if(distance(it_endscaler, it_end)<nScalers) {
//        cout << "WARNING: Scaler block too short, max length=" << distance(it_endscaler, it_end) << endl << endl;
//        it_endscaler = it_end;
//    }
//    else {
//        std::advance(it_endscaler, nScalers);
//    }



    unsigned n = 0;
    while(it != it_endscaler) {
        cout << hex << setw(8) << setfill('0') << *it << " ";
        ++it;
        ++n;
        if(n % 8 == 0)
            cout << endl;
    }

    cout << endl << endl;


    good = true;
    return;

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
