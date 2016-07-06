#include "UnpackerAcqu_FileFormatMk1.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

#include "RawFileReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"

#include <numeric>

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

    // first 32bit word is header marker, that's why
    const acqu::AcquExptInfo_t* h = reinterpret_cast<const acqu::AcquExptInfo_t*>(buffer.data()+1);

    info.Format = Info::Format_t::Mk1;
    info.Time = std_ext::to_tm(h->fTime, "%a %b %d %T %Y");
    info.Description = std_ext::string_sanitize(h->fDescription);
    info.RunNote = std_ext::string_sanitize(h->fRunNote);
    info.OutFile = std_ext::string_sanitize(h->fOutFile);
    info.RunNumber = static_cast<unsigned>(h->fRun);
    info.RecordLength = static_cast<unsigned>(h->fRecLen);

    const auto maxADCIndex = h->fNspect;
    const auto maxScalerIndex = h->fNscaler;
    const auto nModules = h->fNmodule;

    // see TDAQsupervise::CreateMk1Header( void* buff ) in Acqu
    const size_t minBytes = 4 // dont't forget the 32bit start marker
                            + sizeof(acqu::AcquExptInfo_t) // not multiple of 4!
                            + sizeof(acqu::ADCInfo_t)*(maxADCIndex + maxScalerIndex)
                            + sizeof(acqu::ModuleInfo_t)*nModules;

    // the pointer h is only up to now, since buffer might be relocated
    // since AcquExptInfo_t is not word-aligned (4byte), we need to count in bytes...
    reader->expand_buffer(buffer, minBytes/4+1);
    auto ADCInfo_offset =
            reinterpret_cast<const acqu::ADCInfo_t*>(
                reinterpret_cast<const acqu::AcquExptInfo_t*>(
                    buffer.data()
                    +1 // skip 4byte long header word
                    )
                +1 // skip AcquExptInfo_t
                );

    auto ScalerInfo_offset = ADCInfo_offset + maxADCIndex;
    auto ModuleInfo_offset = reinterpret_cast<const acqu::ModuleInfo_t*>(ScalerInfo_offset + maxScalerIndex);

    /// \todo ADCs could be checked against hit mappings

    // use ADCInfos of scalers to determine split of scaler buffers
    vector<string> scaler_modnames;
    for(unsigned i=0;i<maxScalerIndex;i++) {
        const acqu::ADCInfo_t* scalerinfo = ScalerInfo_offset + i;
        if(scalerinfo->fModIndex >= nModules) {
            throw Exception("Invalid fModIndex encountered");
        }
        const acqu::ModuleInfo_t* m = ModuleInfo_offset + scalerinfo->fModIndex;

        scaler_modnames.emplace_back(m->fName);
    }

    FindScalerBlocks(scaler_modnames);

    for(unsigned i=0;i<nModules;i++) {
        const acqu::ModuleInfo_t* m = ModuleInfo_offset + i;

        Info::HardwareModule module;
        module.Identifier = m->fName;
        module.Index = m->fBusType; // is fIndex according to TDAQmodule::ReadHeader( ModuleInfo_t* mod )
        module.Bits = m->fBits;
        module.FirstRawChannel = m->fAmin;
        module.NRawChannels =  m->fAmax - m->fAmin + 1;
        info.Modules.emplace_back(move(module));
    }

    VLOG(9) << "Header says: Have " << info.Modules.size() << " modules";
}

void acqu::FileFormatMk1::FindScalerBlocks(const std::vector<string>& scaler_modnames)
{
    /// \todo the heuristic here test only with 2007 data, where it produces the meaningful values
    // search for "LRS2551" as indicator of scaler block
    // then calculate scaler block sizes
    vector<unsigned> block_offsets;
    auto it = scaler_modnames.begin();
    while(it != scaler_modnames.end()) {
        if(*it == "LRS2551") {
            block_offsets.emplace_back(std::distance(scaler_modnames.begin(), it));
        }
        do {
            ++it;
        }
        while(it != scaler_modnames.end() && *it == *std::prev(it));
    }

    if(block_offsets.empty()) {
        throw Exception("No scaler blocks beginning with LRS2551 found");
    }

    if(block_offsets.front() != 0) {
        throw Exception("Unexpected first scaler block found");
    }

    // contains at least two elements
    block_offsets.push_back(scaler_modnames.size());
    assert(block_offsets.size()>1);

    ScalerBlockSizes.resize(block_offsets.size());
    std::adjacent_difference(block_offsets.begin(), block_offsets.end(),
                             ScalerBlockSizes.begin());
    ScalerBlockSizes.pop_front();

}

void acqu::FileFormatMk1::FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const
{
    // search at 0x8000 bytes
    if(SearchFirstDataBuffer(reader, buffer, 0x8000))
        return;

    // else fail
    throw UnpackerAcqu::Exception("Did not find first data buffer with Mk1 signature");
}

void acqu::FileFormatMk1::UnpackEvent(TEventData& eventdata, it_t& it, const it_t& it_endbuffer, bool& good) noexcept
{
    /// \todo Scan mappings if there's an ADC channel defined which mimicks those markers

    // expect the first word to be the event id
    eventdata.Trigger.DAQEventID = *it;
    ++it;

    hit_storage.clear();
    // there might be more than one scaler block in each event, so
    // so collect them first in this map
    scalers_t scalers;
    auto it_scalerblock = ScalerBlockSizes.cbegin();
    while(it != it_endbuffer && *it != acqu::EEndEvent) {
        // note that the Handle* methods move the iterator
        // themselves and set good to true if nothing went wrong

        good = false;

        switch(*it) {
        case acqu::EScalerBuffer:
            // Scaler read in this event
            HandleScalerBuffer(it_scalerblock, scalers, it, it_endbuffer, good);
            break;
        case acqu::EReadError:
            // read error block, some hardware-related information
            HandleDAQError(eventdata.Trigger.DAQErrors, it, it_endbuffer, good);
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

    if(!scalers.empty() && it_scalerblock != ScalerBlockSizes.cend()) {
        /// \todo check if this requirement is not too hard for some beam times...
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Not all scaler blocks found in same event"
                   );
        good = false;
        return;
    }

    // hit_storage is member variable for better memory allocation performance
    FillDetectorReadHits(hit_storage, hit_mappings_ptr, eventdata.DetectorReadHits);
    FillSlowControls(scalers, scaler_mappings, eventdata.SlowControls);

    ++it; // go to start word of next event (if any)
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
        ScalerBlockSizes_t::const_iterator& it_scalerblock,
        scalers_t& scalers,
        it_t& it, const it_t& it_end,
        bool& good
        ) const noexcept
{

    // skip Scaler buffer marker
    it++;

    auto nWordsRemaining = std::distance(it, it_end);

    if(nWordsRemaining < *it_scalerblock) {
        LogMessage(TUnpackerMessage::Level_t::DataError,
                   "Expected scaler block not completely present in event"
                   );
        return;
    }

    auto it_endscaler = std::next(it, *it_scalerblock);

    // start with offset
    /// \todo would be better to use header info to make scaler
    /// indices more meaningful, for example blocked by module names...
    auto scalerIndex = std::accumulate(ScalerBlockSizes.cbegin(), it_scalerblock, 0);

    while(it != it_endscaler) {
        scalers[scalerIndex].push_back(*it);
        ++scalerIndex;
        ++it;
    }

    // go past last scaler, and to next scaler block (if any)
    ++it;
    ++it_scalerblock;

    good = true;
}
