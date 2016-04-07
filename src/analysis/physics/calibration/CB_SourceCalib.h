#include "analysis/physics/Physics.h"
#include "root-addons/cbtaps_display/TH2CB.h"

#include "calibration/converters/GeSiCa_SADC.h"


namespace ant {
namespace analysis {
namespace physics {

class CB_SourceCalib : public Physics {


protected:
    TH2* ggIM = nullptr;
    TH2CB* h_cbdisplay = nullptr;
    calibration::converter::GeSiCa_SADC adc_converter;

public:

    CB_SourceCalib(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};

}
}
} //namespace ant::analysis::physics

