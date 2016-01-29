#pragma once

#include "UnpackerAcqu_detail.h"



namespace ant {

struct TDAQError;
struct TSlowControl;

namespace unpacker {
namespace acqu {

class FileFormatMk2 : public FileFormatBase {

    // UnpackerAcquFile interface
protected:
    virtual size_t SizeOfHeader() const override;
    virtual bool InspectHeader(const std::vector<std::uint32_t> &buffer) const override;
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) const override;
    virtual void FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) override;

    virtual bool UnpackDataBuffer(queue_t &queue, it_t& it, const it_t& it_endbuffer) noexcept override;

private:
    using scalers_t = std::map<uint32_t, std::vector<uint32_t> >;

    bool SearchFirstDataBuffer(reader_t& reader, buffer_t& buffer, size_t offset);
    void UnpackEvent(queue_t& queue, it_t& it, const it_t& it_endbuffer, bool& good) noexcept;
    void FillDetectorReadHits(std::vector<TDetectorReadHit>& hits) const noexcept;
    void HandleScalerBuffer(std::vector<TSlowControl>& slowcontrols,
                            it_t& it, const it_t& it_end, bool& good,
                            std::vector<TDAQError>& errors) noexcept;
    void HandleDAQError(std::vector<TDAQError>& errors,
                        it_t& it, const it_t& it_end, bool& good) noexcept;
    void HandleEPICSBuffer(std::vector<TSlowControl>& slowcontrols, it_t& it, const it_t& it_end, bool& good) noexcept;
};

}}} // namespace ant::unpacker::acqu
