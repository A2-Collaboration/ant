#include "UnpackerAcqu_FileFormatMk2.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "RawFileReader.h"

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

void acqu::FileFormatMk2::FillInfo(reader_t &reader, buffer_t &buffer, Info &info)
{
    const acqu::AcquMk2Info_t* h = reinterpret_cast<const acqu::AcquMk2Info_t*>(buffer.data()+1);

    // fill the Info

    // the C++11 istringstream way is not working on older systems
    // since std::get_time is missing :(
    //istringstream ss_time(h->fTime);
    //ss_time >> get_time(&info.Time,"%a %b %d %T %Y");

    info.Format = Info::Format_t::Mk2;
    info.Time = std_ext::to_tm(h->fTime, "%a %b %d %T %Y");
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

    for(unsigned i=0;i<nModules;i++) {
        const auto buffer_ptr = addressof(buffer[infoSize+i*moduleSize]);
        const acqu::ModuleInfoMk2_t* m =
                reinterpret_cast<const acqu::ModuleInfoMk2_t*>(buffer_ptr);
        auto it = ModuleIDToString.find(m->fModID);
        if(it == ModuleIDToString.end()) {
            LogMessage(TUnpackerMessage::Level_t::Warn,
                       std_ext::formatter()
                       << "Skipping unknown module with ID=0x" << hex << m->fModID << dec
                       );
            continue;
        }

        Info::HardwareModule module;
        module.Identifier = it->second;
        module.Index = m->fModIndex;
        module.Bits = m->fBits;
        module.FirstRawChannel = m->fAmin;

        // ADC and Scaler will be added twice!
        if(m->fModType & acqu::EDAQ_ADC) {
            module.NRawChannels = m->fNChannel;
            info.Modules.emplace_back(module);
        }
        if(m->fModType & acqu::EDAQ_Scaler) {
            module.NRawChannels = m->fNScChannel;
            info.Modules.emplace_back(module);
        }
    }

    VLOG(9) << "Header says: Have " << info.Modules.size() << " modules";

    LogMessage(TUnpackerMessage::Level_t::Info,
               "Acqu Mk2 header successfully unpacked"
               );
}

void acqu::FileFormatMk2::FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const
{
    // finally search for the Mk2Header, this also
    // fills the buffer correctly with the first Mk2DataBuffer (if available)

    // first search at 0x8000 bytes
    if(SearchFirstDataBuffer(reader, buffer, 0x8000))
        return;

    // then search at 10*0x8000 bytes
    if(SearchFirstDataBuffer(reader, buffer, 10*0x8000))
        return;

    // then search at 16*0x8000 bytes
    if(SearchFirstDataBuffer(reader, buffer, 16*0x8000))
        return;

    // if we did not find the MK2DataBuffer at the positions above, try to search it
    // with a maximum multiplier of 32 and assert that it's at a multiplicity of 0x8000
    if(FindFirstDataBuffer(reader, buffer, 32, true))
        return;

    // else fail
    throw UnpackerAcqu::Exception("Did not find first data buffer with Mk2 signature");
}



void acqu::FileFormatMk2::UnpackEvent(
        TEventData& eventdata,
        it_t& it, const it_t& it_endbuffer,
        bool& good
        ) noexcept
{
    // extract and check eventLength
    const unsigned eventLength = *it/sizeof(decltype(*it));
    {
        const auto nRemainingWords = std::distance(it, it_endbuffer);
        if(nRemainingWords < eventLength) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       std_ext::formatter() <<
                       "Event with size 0x" << eventLength
                       << " too big to fit in buffer of remaining size "
                       << nRemainingWords
                       );
            return;
        }
    }
    const auto it_endevent = next(it, eventLength);
    if(*it_endevent != acqu::EEndEvent) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   std_ext::formatter() <<
                   "At designated end of event, found unexpected word 0x"
                   << hex << *it_endevent << dec
                   );
        return;
    }
    it++;

    // now work on one event inside buffer

    /// \todo Scan mappings if there's an ADC channel defined which mimicks those blocks

    hit_storage.clear();
    // there might be more than one scaler block in each event, so
    // so collect them first in this map
    scalers_t scalers;
    while(it != it_endbuffer && *it != acqu::EEndEvent) {
        // note that the Handle* methods move the iterator
        // themselves and set good to true if nothing went wrong

        good = false;

        switch(*it) {
        case acqu::EEPICSBuffer:
            // EPICS buffer
            HandleEPICSBuffer(eventdata.SlowControls, it, it_endevent, good);
            break;
        case acqu::EScalerBuffer:
            // Scaler read in this event
            // there might also be some error buffers in between, args...
            HandleScalerBuffer(
                        scalers,
                        it, it_endevent, good,
                        eventdata.Trigger.DAQErrors);
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

    // check proper EEndEvent
    if(it == it_endbuffer) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   std_ext::formatter() <<
                   "While unpacking event, found premature end of buffer."
                   );
        return;
    }

    // hit_storage is member variable for better memory allocation performance
    FillDetectorReadHits(hit_storage, hit_mappings_ptr, eventdata.DetectorReadHits);
    FillSlowControls(scalers, scaler_mappings, eventdata.SlowControls);

    it++; // go to start word of next event (if any)
}

void acqu::FileFormatMk2::HandleScalerBuffer(
        scalers_t& scalers,
        it_t& it, const it_t& it_end,
        bool& good,
        std::vector<TDAQError>& errors
        ) const noexcept
{
    // ignore Scaler buffer marker
    it++;

    if(it==it_end) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Acqu ScalerBlock only start marker found"
                   );
        return;
    }

    // get the scaler block length in words
    const int scalerLength = *it;
    constexpr int wordsize = sizeof(decltype(*it));
    if(scalerLength % wordsize != 0
       || distance(it,it_end)<scalerLength/wordsize) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Acqu ScalerBlock length invalid"
                   );
        return;
    }

    const auto it_endscaler = next(it, scalerLength/wordsize);
    if(*it_endscaler != acqu::EScalerBuffer) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   std_ext::formatter()
                   << "Acqu ScalerBlock did not have proper end marker: "
                   << hex << "0x" << *it_endscaler
                   );
        return;
    }
    it++; // skip the length word now

    // fill simple scalers map
    while(it != it_endscaler) {
        // within a scaler block, there might be error blocks
        if(*it == acqu::EReadError) {
            HandleDAQError(errors, it, it_end, good);
            if(!good)
                return;
            good = false;
        }

        // check if enough space left
        if(distance(it, it_endscaler) < 2) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       "Acqu ScalerBlock contains malformed scaler read"
                       );
            return;
        }

        const uint32_t index = *it++;
        const uint32_t value = *it++;

        scalers[index].push_back(value);
    }

    // this position is only reached when it==it_endscaler
    // skip the scaler buffer end marker
    // already checked above with it_endscaler
    it++;

    good = true;
}

void acqu::FileFormatMk2::HandleDAQError(std::vector<TDAQError>& errors,
                                         it_t& it, const it_t& it_end,
                                         bool& good) const noexcept
{
    // is there enough space in the event at all?
    static_assert(sizeof(acqu::ReadErrorMk2_t) % sizeof(decltype(*it)) == 0,
                  "acqu::ReadErrorMk2_t is not word aligned");
    constexpr int wordsize = sizeof(acqu::ReadErrorMk2_t)/sizeof(decltype(*it));

    if(distance(it, it_end)<wordsize) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "acqu::ReadErrorMk2_t block not completely present in buffer"
                   );
        return;
    }

    // then cast it to data structure
    const acqu::ReadErrorMk2_t* err =
            reinterpret_cast<const acqu::ReadErrorMk2_t*>(addressof(*it));

    // check for trailer word
    if(err->fTrailer != acqu::EReadError) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Acqu ErrorBlock does not end with expected trailer word"
                   );
        return;
    }

    // lookup the module name
    auto it_modname = acqu::ModuleIDToString.find(err->fModID);
    const string& modname = it_modname == acqu::ModuleIDToString.cend()
                            ? "UNKNOWN" : it_modname->second;

    errors.emplace_back(err->fModID, err->fModIndex, err->fErrCode, modname);

    VLOG_N_TIMES(1000, 2) << errors.back();

    advance(it, wordsize);
    good = true;
}

void acqu::FileFormatMk2::HandleEPICSBuffer(
        std::vector<TSlowControl>& slowcontrols,
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
        LogMessage(TUnpackerMessage::Level_t::DataError,
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
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "EPICS data not word aligned"
                   );
        return;
    }
    const int epicsTotalWords = hdr->len/wordbytes;
    const string epicsModName = hdr->name;
    const size_t nChannels = hdr->nchan;
    const time_t hdr_timestamp = hdr->time;

    if(epicsModName.length()>32) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "acqu::EpicsHeaderInfo_t header has malformed module name"
                   );
        return;
    }
    if(distance(it, it_end)<epicsTotalWords) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
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
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       "EPICS channel header not completely present in buffer"
                       );
            return;
        }

        const acqu::EpicsChannelInfo_t* ch =
                reinterpret_cast<const acqu::EpicsChannelInfo_t*>(addressof(*it_byte));

        if(distance(it_byte, bytes.cend()) < ch->bytes) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       "EPICS channel payload not completely present in buffer"
                       );
            return;
        }

        auto it_map = map_EpicsTypes.find(ch->type);
        if(it_map == map_EpicsTypes.cend()) {
            LogMessage(TUnpackerMessage::Level_t::DataError,
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
            LogMessage(TUnpackerMessage::Level_t::DataError,
                       "EPICS channel payload size inconsistent"
                       );
            return;
        }

        // finally we can create the TSlowControl record

        auto record_type = TSlowControl::Type_t::EpicsOneShot;
        auto validity = TSlowControl::Validity_t::Forward;
        stringstream description;
        if(hdr->period<0) {
            record_type = TSlowControl::Type_t::EpicsTimer;
            validity = TSlowControl::Validity_t::Backward;
            description << "Period='" << -hdr->period << " ms'";
        }
        else if(hdr->period>0) {
            record_type = TSlowControl::Type_t::EpicsScaler;
            validity = TSlowControl::Validity_t::Backward;
            description << "Period='" << hdr->period << " scalers'";
        }

        slowcontrols.emplace_back(
                    record_type,
                    validity,
                    hdr_timestamp,
                    ch_Name,
                    description.str()
                    );
        auto& sc = slowcontrols.back();

        // advance to the EPICS channel data (skip channel info header)
        advance(it_byte, chHdrBytes);

        // fill the payload depending on the EPICS data type
        // upcast float to double and short,byte to long

        for(unsigned elem=0;elem<(unsigned)ch_nElements;elem++) {
            switch(ch_datatype) {
            case acqu::EpicsDataTypes_t::BYTE:
                sc.Payload_Int.emplace_back(elem, *it_byte);
                break;
            case acqu::EpicsDataTypes_t::SHORT: {
                const int16_t* value = reinterpret_cast<const int16_t*>(addressof(*it_byte));
                sc.Payload_Int.emplace_back(elem, *value);
                break;
            }
            case acqu::EpicsDataTypes_t::LONG: {
                const int64_t* value = reinterpret_cast<const int64_t*>(addressof(*it_byte));
                sc.Payload_Int.emplace_back(elem, *value);
                break;
            }
            case acqu::EpicsDataTypes_t::FLOAT: {
                static_assert(sizeof(float)==4,"Float should be 4 bytes long");
                const float* value = reinterpret_cast<const float*>(addressof(*it_byte));
                sc.Payload_Float.emplace_back(elem, *value);
                break;
            }
            case acqu::EpicsDataTypes_t::DOUBLE: {
                static_assert(sizeof(double)==8,"Float should be 8 bytes long");
                const double* value = reinterpret_cast<const double*>(addressof(*it_byte));
                sc.Payload_Float.emplace_back(elem, *value);
                break;
            }
            case acqu::EpicsDataTypes_t::STRING: {
                const char* value = reinterpret_cast<const char*>(addressof(*it_byte));
                // interpret as string
                const string value_str(value);
                if((signed)value_str.length()>=ch_typesize) {
                    LogMessage(TUnpackerMessage::Level_t::DataError,
                               "EPICS channel string data too long (no terminating \\0?)"
                               );
                    return;
                }
                sc.Payload_String.emplace_back(elem, value);
                break;
            }
            default:
                LOG(ERROR) << "Not implemented";

            } // end switch

            advance(it_byte, ch_typesize);
        }

        VLOG(9) << sc;

    } // end channel loop

    // we successfully parsed the EPICS buffer
    advance(it, epicsTotalWords);

    VLOG(9) << "Successfully parsed EPICS buffer";

    good = true;
}
