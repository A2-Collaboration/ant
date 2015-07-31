#pragma once

#include "analysis/input/DataReader.h"
#include "analysis/data/Event.h"

#include <memory>
#include <string>
#include <list>

#include "detail/InputModule.h"

#include "detail/TriggerInput.h"
#include "detail/TaggerInput.h"
#include "detail/DetectorHitInput.h"
#include "detail/TrackInput.h"
#include "detail/PlutoInput.h"
#include "detail/ParticleInput.h"


class PStaticData;


namespace ant {

class ReadTFiles;

namespace input {

class TreeManager;

class GoatReader: public DataReader {
protected:

    class InputWrapper {

    };

    class ModuleManager: public std::list<BaseInputModule*> {
    public:
        ModuleManager() = default;
        ModuleManager(const std::initializer_list<BaseInputModule*>& initlist):
            std::list<BaseInputModule*>(initlist) {}

        void GetEntry() {
            for(auto& module : *this) {
                module->GetEntry();
            }
        }
    };

    std::unique_ptr<TreeManager>   trees;


    TriggerInput        trigger;
    TaggerInput         tagger;
    TrackInput          tracks;
    DetectorHitInput    detectorhits;
    PlutoInput          pluto;
    ParticleInput       photons   = ParticleInput("photons");
    ParticleInput       protons   = ParticleInput("protons");
    ParticleInput       pichagred = ParticleInput("pions");
    ParticleInput       echarged  = ParticleInput("echarged");
    ParticleInput       neutrons  = ParticleInput("neutrons");

    ModuleManager active_modules = {
        &trigger,
        &tagger,
        &tracks,
        &detectorhits,
        &pluto,
        &photons,
        &protons,
        &pichagred,
        &echarged,
        &neutrons
    };

    Long64_t    current_entry = -1;
    Long64_t    max_entry = std::numeric_limits<Long64_t>::max();

    static clustersize_t MapClusterSize(const int& size);

    void CopyTagger(std::shared_ptr<Event>& event);
    void CopyTrigger(std::shared_ptr<Event>& event);

    /**
     * @brief CopyDetectorHits
     * @param event
     * @todo implement
     */
    void CopyDetectorHits(std::shared_ptr<Event>& event);
    void CopyTracks(std::shared_ptr<Event>& event);
    void CopyPluto(std::shared_ptr<Event>& event);
    void CopyParticles(std::shared_ptr<Event>& event, ParticleInput& input_module, const ParticleTypeDatabase::Type& type);

    PStaticData* pluto_database;
    const ParticleTypeDatabase::Type* GetType(const PParticle* p) const;

    /**
     * @brief Get number of events in tree
     * @see TotalEvents()
     * @return number of events total
     */
    Long64_t  GetNEvents() const;

public:
    GoatReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~GoatReader();
    GoatReader(const GoatReader&) = delete;
    GoatReader& operator= (const GoatReader&) = delete;

    std::shared_ptr<Event> ReadNextEvent();
    virtual bool hasData() const override;

    virtual long long EventsRead() const override;
    virtual long long TotalEvents() const override;

    void SetMaxEntries(const long long max);

    const PlutoInput& GetPlutoInput() { return pluto; }
};

}
}
