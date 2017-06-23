#pragma once

#include "calibration/Calibration.h"
#include "base/interval.h"

class TGraph;

namespace ant {

struct TaggerDetector_t;

namespace calibration {

class DataManager;

class TaggEff :
        public Calibration::BaseModule,
        public Updateable_traits,
        public ReconstructHook::EventData
{
public:
    TaggEff(const std::shared_ptr<ant::TaggerDetector_t>&  tagger,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~TaggEff();

    static std::string GetModuleName(Detector_t::Type_t type);

    // Updateable interface
    virtual std::list<Loader_t> GetLoaders() override;

    // ReconstructHook interface
    virtual void ApplyTo(TEventData& reconstructed) override;

protected:
    std::shared_ptr<ant::TaggerDetector_t> Tagger;
    std::shared_ptr<DataManager> CalibrationManager;
    std::vector<TaggerDetector_t::taggeff_t> currentTaggEff;
    TID loadedTaggEff; // timestamp when currentTaggEff was loaded

};

}}
