#pragma once

#include "UnpackerAcqu_detail.h"

namespace ant {

struct TDAQError;

namespace unpacker {
namespace acqu {

class FileFormatMk1 : public FileFormatBase {


protected:

    unsigned AcquID_last = 0;
    unsigned nScalers;

    virtual size_t SizeOfHeader() const override;
    virtual bool InspectHeader(const std::vector<std::uint32_t>& buffer) const override;
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) override;
    virtual void FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const override;
    virtual bool UnpackDataBuffer(queue_t& queue, it_t& it, const it_t& it_endbuffer) noexcept override;

    void UnpackEvent(queue_t& queue, it_t& it, const it_t& it_endbuffer, bool& good) noexcept;
    void HandleDAQError(std::vector<TDAQError>& errors,
                        it_t& it, const it_t& it_end, bool& good) const noexcept;
    void HandleScalerBuffer(scalers_t& scalers,
                            it_t& it, const it_t& it_end, bool& good) const noexcept;
};

}}} // namespace ant::unpacker::acqu
