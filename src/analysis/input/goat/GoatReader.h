#pragma once

#include "analysis/input/DataReader.h"

#include <memory>
#include <string>
#include <list>

#include "detail/InputModule.h"

#include "detail/TriggerInput.h"
#include "detail/EventParameters.h"
#include "detail/TaggerInput.h"
#include "detail/DetectorHitInput.h"
#include "detail/TrackInput.h"

#include "base/types.h"

namespace ant {

class WrapTFileInput;
struct TEventData;

namespace expconfig {
namespace detector {
struct Trigger;
}}

namespace analysis {
namespace input {

class TreeManager;

class GoatReader: public DataReader {
protected:

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

    std::shared_ptr<WrapTFileInput>    files;
    std::unique_ptr<TreeManager>   trees;

    std::shared_ptr<expconfig::detector::Trigger> det_trigger;

    TriggerInput        trigger;
    EventParameters     eventParameters;
    TaggerInput         tagger;
    TrackInput          tracks;
    DetectorHitInput    detectorhits;

    ModuleManager active_modules = {
        &trigger,
        &eventParameters,
        &tagger,
        &tracks,
        &detectorhits,
    };

    Long64_t    current_entry = 0;

    static clustersize_t MapClusterSize(const int& size);

    void CopyDetectorHits(TEventData& recon);
    void CopyTagger(TEventData& recon);
    void CopyTrigger(TEventData& recon);
    void CopyTracks(TEventData& recon);

    /**
     * @brief Get number of events in tree
     * @see TotalEvents()
     * @return number of events total
     */
    Long64_t  GetNEvents() const;

public:
    GoatReader(const std::shared_ptr<WrapTFileInput>& rootfiles);
    virtual ~GoatReader();
    GoatReader(const GoatReader&) = delete;
    GoatReader& operator= (const GoatReader&) = delete;

    virtual bool IsSource() override;
    virtual bool ReadNextEvent(event_t& event) override;

    double PercentDone() const override;
};

}
}
}
