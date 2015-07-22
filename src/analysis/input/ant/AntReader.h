#ifndef ANTREDER_H
#define ANTREDER_H

#include "analysis/input/DataReader.h"

#include "Rtypes.h"

#include <memory>
#include <string>


class TTree;

namespace ant {

class TEvent;

namespace input {

class FileManager;

class AntReader: public DataReader {
protected:

    std::unique_ptr<ant::input::FileManager> files;

    TTree* tree    = nullptr;
    TEvent* buffer = nullptr;

    Long64_t current = -1;

public:
    AntReader();
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;

    void AddInputFile(const std::string& filename);
    void Initialize();

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

#endif
