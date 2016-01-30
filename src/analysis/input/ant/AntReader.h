#pragma once

#include "analysis/input/DataReader.h"

#include "unpacker/Unpacker.h"

#include "reconstruct/Reconstruct_traits.h"

#include "base/WrapTFile.h"

#include <memory>
#include <string>


class TTree;

namespace ant {

struct Event;

namespace analysis {

namespace data {
    struct Event;
}

namespace input {


class AntReader : public DataReader {
protected:
    std::unique_ptr<Unpacker::Module>   unpacker;
    std::unique_ptr<Reconstruct_traits> reconstruct;

public:
    AntReader(const std::shared_ptr<WrapTFileInput>& rootfiles,
              std::unique_ptr<Unpacker::Module> unpacker,
              std::unique_ptr<Reconstruct_traits> reconstruct);
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;

    // DataReader interface
    virtual bool IsSource() override { return true; }
    virtual bool ReadNextEvent(TEvent& event) override;

    double PercentDone() const override;
};

}
}
}
