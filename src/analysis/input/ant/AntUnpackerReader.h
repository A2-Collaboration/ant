#pragma once

#include "analysis/input/DataReader.h"

#include "unpacker/Unpacker.h"
#include "tree/UnpackerWriter.h"

#include "reconstruct/Reconstruct_traits.h"

#include "Rtypes.h"

#include <memory>
#include <string>


class TTree;

namespace ant {

class TEvent;
class ReadTFiles;

namespace input {


class AntUnpackerReader : public DataReader {
protected:
    std::unique_ptr<Unpacker::Reader> reader;
    std::unique_ptr<tree::UnpackerWriter> writer;
    std::unique_ptr<Reconstruct_traits> reconstruct;


    long long nEvents;
    bool writeUncalibrated;
    bool writeCalibrated;

public:
    AntUnpackerReader(std::unique_ptr<Unpacker::Reader> unpacker_reader,
                      std::unique_ptr<Reconstruct_traits> reconstruct);
    virtual ~AntUnpackerReader();
    AntUnpackerReader(const AntUnpackerReader&) = delete;
    AntUnpackerReader& operator= (const AntUnpackerReader&) = delete;


    // DataReader interface
    std::shared_ptr<Event> ReadNextEvent();
    virtual bool hasData() const override;
    virtual long long EventsRead() const override;
    virtual long long TotalEvents() const override;
};

}
}
