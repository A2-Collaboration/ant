#pragma once

#include "UnpackerAcqu_detail.h"

#include "base/mapped_vectors.h"

namespace ant {
namespace unpacker {
namespace acqu {

class FileFormatMk2 : public FileFormatBase {

    // UnpackerAcquFile interface
protected:
    virtual size_t SizeOfHeader() const override;
    virtual bool InspectHeader(const std::vector<std::uint32_t> &buffer) const override;
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) const override;
    virtual void FillFirstDataBuffer(queue_t& queue, reader_t& reader, buffer_t& buffer) const override;

    virtual bool UnpackDataBuffer(queue_t &queue, it_t& it, const it_t& it_endbuffer) noexcept override;

private:
    //using hits_t = std::map<uint16_t, std::vector<uint16_t> >;
    using scalers_t = std::map<uint32_t, std::vector<uint32_t> >;

    using hits_t = std_ext::mapped_vectors<uint16_t, uint16_t>;

    bool SearchFirstDataBuffer(queue_t &queue, reader_t& reader, buffer_t& buffer, size_t offset) const;
    void UnpackEvent(queue_t& queue, it_t& it, const it_t& it_endbuffer, bool& good) noexcept;
    void FillTDetectorRead(queue_t& queue,
                           const hits_t& hits,
                           const scalers_t& scalers) const noexcept;
    void HandleScalerBuffer(queue_t& queue, it_t& it, const it_t& it_end, bool& good,
                            scalers_t &scalers) const noexcept;
    void HandleReadError(queue_t& queue, it_t& it, const it_t& it_end, bool& good) const noexcept;
    void HandleEPICSBuffer(queue_t& queue, it_t& it, const it_t& it_end, bool& good) const noexcept;
};

}}} // namespace ant::unpacker::acqu
