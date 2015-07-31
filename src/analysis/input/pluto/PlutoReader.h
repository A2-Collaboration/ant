#pragma once

#include "analysis/input/DataReader.h"
#include "analysis/data/Event.h"

#include <memory>
#include <string>
#include <list>



class PStaticData;
class TClonesArray;
class TTree;
class PParticle;

namespace ant {

class ReadTFiles;

namespace input {

class TreeManager;

class PlutoReader: public DataReader {
protected:

    std::shared_ptr<ReadTFiles>    files;

    TTree*          tree = nullptr;
    TClonesArray*   PlutoMCTrue = nullptr;
    Long64_t        plutoID = 0;
    Long64_t        plutoRandomID = 0;

    using PParticleVector = std::vector<const PParticle*>;
    PParticleVector PlutoParticles;

    Long64_t    current_entry = 0;

    void CopyPluto(Event& event);

    PStaticData* pluto_database;
    const ParticleTypeDatabase::Type* GetType(const PParticle* p) const;

public:
    PlutoReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~PlutoReader();
    PlutoReader(const PlutoReader&) = delete;
    PlutoReader& operator= (const PlutoReader&) = delete;

    virtual bool ReadNextEvent(Event& event, TSlowControl&) override;

};

}
}
