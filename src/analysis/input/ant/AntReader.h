#pragma once

#include "analysis/input/DataReader.h"

#include "Rtypes.h"

#include <memory>
#include <string>


class TTree;

namespace ant {

class TEvent;
class ReadTFiles;


namespace input {


class AntReader: public FileDataReader {
protected:

    std::unique_ptr<ReadTFiles> files;

    TTree* tree    = nullptr;
    TEvent* buffer = nullptr;

    Long64_t current = -1;

public:
    AntReader();
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;

    void AddInputFile(const std::string& filename) override;
    void Initialize() override;

    /**
     * @brief Get number of events in tree
     * @see TotalEvents()
     * @return number of events total
     */
    Long64_t  GetNEvents() const;

    std::shared_ptr<Event> ReadNextEvent();
    bool hasData() const override;

    long long EventsRead() const override;
    long long TotalEvents() const override;

};

}
}
