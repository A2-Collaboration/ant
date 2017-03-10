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
        public Updateable_traits
{
public:
    TaggEff(const std::shared_ptr<ant::TaggerDetector_t>&  tagger,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~TaggEff();

    static std::string GetModuleName(Detector_t::Type_t type);

    virtual std::list<Loader_t> GetLoaders() override;

    virtual void UpdatedTIDFlags(const TID& tid) override;


protected:
    std::shared_ptr<ant::TaggerDetector_t> Tagger;
    std::shared_ptr<DataManager> CalibrationManager;
};

}}
