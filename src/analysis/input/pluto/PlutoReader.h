#pragma once

#include "analysis/input/DataReader.h"

#include "analysis/utils/A2GeoAcceptance.h"

#include "base/ParticleType.h"
#include "base/WrapTTree.h"

#include <memory>
#include <string>
#include <list>

#include "TClonesArray.h"

class PStaticData;
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

    struct PlutoTree_t : WrapTTree {
        // important to tell TClonesArray what stuff is inside
        ADD_BRANCH_T(TClonesArray, Particles, "PParticle");
    };

    struct TIDTree_t : WrapTTree {
        ADD_BRANCH_T(TID, tid);
    };

    PlutoTree_t plutoTree;
    TIDTree_t tidTree;

    long long current_entry = 0;

    void CopyPluto(TEventData& mctrue);

    PStaticData* pluto_database;

public:
    PlutoReader(const std::shared_ptr<ant::WrapTFileInput>& rootfiles);
    virtual ~PlutoReader();
    PlutoReader(const PlutoReader&) = delete;
    PlutoReader& operator= (const PlutoReader&) = delete;

    virtual reader_flags_t GetFlags() const override { return {}; }
    virtual bool ReadNextEvent(input::event_t& event) override;

    double PercentDone() const override;
};

}
}
}
