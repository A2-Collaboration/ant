#pragma once

#include "calibration/Calibration.h"
#include "base/interval.h"

class TGraph;

namespace ant {

struct TaggerDetector_t;

namespace calibration {

class DataManager;

class TaggEff : public Calibration::Module
{
public:
    TaggEff(const std::shared_ptr<ant::TaggerDetector_t>&  tagger,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~TaggEff();

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >&) override;

    static const std::string GetDataName() {return "TaggEff";}
    static const std::string GetDataErrorsName() {return "TaggEffErorrs";}

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

protected:
    std::shared_ptr<ant::TaggerDetector_t> Tagger;
    std::shared_ptr<DataManager> CalibrationManager;
};

}}
