#pragma once

#include "analysis/input/DataReader.h"
#include "analysis/data/Event.h"

#include "analysis/utils/A2GeoAcceptance.h"

#include <memory>
#include <string>
#include <list>



class PStaticData;
class TClonesArray;
class TTree;
class PParticle;

namespace ant {

class WrapTFileInput;
class TID;

namespace analysis {

namespace data {
    class Event;
}

namespace input {

class TreeManager;

class PlutoReader: public DataReader {
protected:

    utils::A2SimpleGeometry geometry;

    std::shared_ptr<TaggerDetector_t> tagger;

    std::shared_ptr<WrapTFileInput> files; // save pointer to keep extracted TTree pointers valid

    TTree*          tree = nullptr;
    TClonesArray*   PlutoMCTrue = nullptr;
    TID*            tid = nullptr;
    bool tid_from_file = false;

    Long64_t    current_entry = 0;

    void CopyPluto(data::Event& event);

    PStaticData* pluto_database;
    const ParticleTypeDatabase::Type* GetType(const PParticle* p) const;

public:
    PlutoReader(const std::shared_ptr<ant::WrapTFileInput>& rootfiles, const std::shared_ptr<TaggerDetector_t> tagger_);
    virtual ~PlutoReader();
    PlutoReader(const PlutoReader&) = delete;
    PlutoReader& operator= (const PlutoReader&) = delete;

    virtual bool IsSource() override { return false; }

    virtual bool ReadNextEvent(data::Event& event) override;

};

}
}
}
