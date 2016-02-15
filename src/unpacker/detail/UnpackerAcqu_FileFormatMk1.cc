#include "UnpackerAcqu_FileFormatMk1.h"

#include "UnpackerAcqu.h"
#include "UnpackerAcqu_legacy.h"
#include "UnpackerAcqu_templates.h"

#include "RawFileReader.h"

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

void acqu::FileFormatMk1::FillInfo(reader_t& reader, buffer_t& buffer, Info& info) const
{
    const acqu::AcquExptInfo_t* h = reinterpret_cast<const acqu::AcquExptInfo_t*>(buffer.data()+1);

    info.Format = Info::Format_t::Mk1;
    info.Time = std_ext::to_tm(h->fTime, "%a %b %d %T %Y");
    info.Description = std_ext::string_sanitize(h->fDescription);
    info.RunNote = std_ext::string_sanitize(h->fRunNote);
    info.OutFile = std_ext::string_sanitize(h->fOutFile);
    info.RunNumber = static_cast<unsigned>(h->fRun);
    info.RecordLength = static_cast<unsigned>(h->fRecLen);

    /// \todo parse some more stuff from the Mk1 header here

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
                                      it_t& it,
                                      const it_t& it_endbuffer, bool& good) noexcept
{

}


