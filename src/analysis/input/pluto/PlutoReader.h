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

    Long64_t    current_entry = -1;

    void CopyPluto(std::shared_ptr<Event>& event);

    PStaticData* pluto_database;
    const ParticleTypeDatabase::Type* GetType(const PParticle* p) const;

public:
    PlutoReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~PlutoReader();
    PlutoReader(const PlutoReader&) = delete;
    PlutoReader& operator= (const PlutoReader&) = delete;

    std::shared_ptr<Event> ReadNextEvent();
    virtual bool hasData() const override;

    virtual long long EventsRead() const override;
    virtual long long TotalEvents() const override;
};

}
}
