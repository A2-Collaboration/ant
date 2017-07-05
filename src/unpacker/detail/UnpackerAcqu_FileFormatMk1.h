#pragma once

#include "UnpackerAcqu_detail.h"

namespace ant {

struct TDAQError;

namespace unpacker {
namespace acqu {

class FileFormatMk1 : public FileFormatBase {
public:
    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

protected:

    using ScalerBlockSizes_t = std::list<unsigned>;
    ScalerBlockSizes_t ScalerBlockSizes;

    virtual size_t SizeOfHeader() const override;
    virtual bool InspectHeader(const std::vector<std::uint32_t>& buffer) const override;
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) override;
    virtual void FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const override;
    virtual void UnpackEvent(TEventData& eventdata, it_t& it, const it_t& it_endbuffer, bool& good) noexcept override;

    void FindScalerBlocks(const std::vector<std::uint16_t>& scaler_bustype, const std::vector<std::uint16_t>& scaler_bits);

    void HandleDAQError(std::vector<TDAQError>& errors,
                        it_t& it, const it_t& it_end, bool& good) const noexcept;
    void HandleScalerBuffer(ScalerBlockSizes_t::const_iterator& it_scalerblock,
                            scalers_t& scalers, it_t& it, const it_t& it_end, bool& good,
                            std::vector<TDAQError>& errors) const noexcept;
};

}}} // namespace ant::unpacker::acqu
