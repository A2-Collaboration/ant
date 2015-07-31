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

    std::shared_ptr<ReadTFiles>    files;
    std::unique_ptr<TreeManager>   trees;


    TriggerInput        trigger;
    TaggerInput         tagger;
    TrackInput          tracks;
    DetectorHitInput    detectorhits;
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
        &photons,
        &protons,
        &pichagred,
        &echarged,
        &neutrons
    };

    Long64_t    current_entry = 0;

    static clustersize_t MapClusterSize(const int& size);

    void CopyTagger(Event& event);
    void CopyTrigger(Event& event);
    void CopyTracks(Event& event);
    void CopyParticles(Event& event, ParticleInput& input_module, const ParticleTypeDatabase::Type& type);


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

    virtual bool ReadNextEvent(Event& event,  TSlowControl& slowControl) override;
};

}
}
