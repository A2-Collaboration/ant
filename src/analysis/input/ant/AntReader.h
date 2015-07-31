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


class AntReader : public DataReader {
protected:
    std::shared_ptr<ReadTFiles>    files;


    TTree* tree    = nullptr;
    TEvent* buffer = nullptr;

    Long64_t current = -1;

    /**
     * @brief Get number of events in tree
     * @see TotalEvents()
     * @return number of events total
     */
    Long64_t  GetNEvents() const;

public:
    AntReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;



    std::shared_ptr<Event> ReadNextEvent();
    virtual bool hasData() const override;

    virtual long long EventsRead() const override;
    virtual long long TotalEvents() const override;

};

}
}
