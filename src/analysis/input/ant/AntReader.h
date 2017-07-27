#pragma once

#include "analysis/input/DataReader.h"

#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct_traits.h"
#include "base/WrapTFile.h"

#include <memory>
#include <string>

namespace ant {
namespace analysis {
namespace input {

namespace detail {
struct AntReaderInternal;
}

class AntReader : public DataReader {

protected:
    std::unique_ptr<detail::AntReaderInternal> reader;
    std::unique_ptr<Reconstruct_traits>        reconstruct;

public:
    AntReader(const std::shared_ptr<WrapTFileInput>& rootfiles,
              std::unique_ptr<Unpacker::Module> unpacker,
              std::unique_ptr<Reconstruct_traits> reconstruct_);
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;

    // DataReader interface
    virtual reader_flags_t GetFlags() const override;
    virtual bool ReadNextEvent(event_t& event) override;

    double PercentDone() const override;
};

}
}
}
