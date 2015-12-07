#pragma once

#include "analysis/input/DataReader.h"
#include "analysis/data/Event.h"



#include <memory>
#include <string>
#include <list>

class TTree;

namespace ant {

class WrapTFileInput;
struct TID;

namespace analysis {

namespace data {
    class Event;
}

namespace input {

class GeantReader : public DataReader {
protected:

    std::shared_ptr<WrapTFileInput> files; // save pointer to keep extracted TTree pointers valid

    TTree*          tree = nullptr;

    Long64_t    current_entry = 0;


    Float_t         fvertex[3] = {};

public:
    GeantReader(const std::shared_ptr<ant::WrapTFileInput>& rootfiles);
    virtual ~GeantReader();
    GeantReader(const GeantReader&) = delete;
    GeantReader& operator= (const GeantReader&) = delete;

    virtual bool IsSource() override { return false; }

    virtual bool ReadNextEvent(data::Event& event) override;

    double PercentDone() const override;
};

}
}
}
