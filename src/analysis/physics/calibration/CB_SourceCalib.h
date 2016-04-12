#include "analysis/physics/Physics.h"
#include "root-addons/cbtaps_display/TH2CB.h"

#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"
#include "expconfig/detectors/CB.h"


namespace ant {
namespace analysis {
namespace physics {

class CB_SourceCalib : public Physics {


protected:
    TH2* HitsADC = nullptr;
    TH2* TimeHits = nullptr;
    TH2* HitsADC_Cluster = nullptr;
    TH2* VerworfeneHits = nullptr;
    TH1* KristallHits = nullptr;
    TH2CB* h_cbdisplay = nullptr;
    calibration::converter::GeSiCa_SADC adc_converter;
    calibration::converter::CATCH_TDC   tdc_converter;

    std::shared_ptr<expconfig::detector::CB> cb_detector;

public:

    CB_SourceCalib(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}
}
} //namespace ant::analysis::physics

