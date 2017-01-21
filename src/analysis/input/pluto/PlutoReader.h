#pragma once

#include "analysis/input/DataReader.h"

#include "analysis/utils/A2GeoAcceptance.h"

#include "base/ParticleType.h"

#include <memory>
#include <string>
#include <list>



class PStaticData;
class TClonesArray;
class TTree;
class PParticle;

namespace ant {

class WrapTFileInput;
struct TID;
struct TEventData;

namespace analysis {
namespace input {

class PlutoReader: public DataReader {
protected:

    utils::A2SimpleGeometry geometry;

    std::shared_ptr<TaggerDetector_t> tagger;

    std::shared_ptr<WrapTFileInput> files; // save pointer to keep extracted TTree pointers valid

    TTree*          tree = nullptr;
    TClonesArray*   PlutoMCTrue = nullptr;
    TID*            tid = nullptr;
    bool tid_from_file = false;

    long long current_entry = 0;

    void CopyPluto(TEventData& mctrue);

    PStaticData* pluto_database;

public:
    PlutoReader(const std::shared_ptr<ant::WrapTFileInput>& rootfiles);
    virtual ~PlutoReader();
    PlutoReader(const PlutoReader&) = delete;
    PlutoReader& operator= (const PlutoReader&) = delete;

    virtual bool IsSource() override { return false; }

    virtual bool ReadNextEvent(input::event_t& event) override;

    double PercentDone() const override;
};

}
}
}
