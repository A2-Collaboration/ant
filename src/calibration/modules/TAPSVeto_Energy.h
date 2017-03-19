#pragma once

#include "Energy.h"

class TH1;

namespace ant {

namespace expconfig {
namespace detector {
struct TAPSVeto;
}}

namespace calibration {

class TAPSVeto_Energy : public Energy
{


public:

    using detector_ptr_t = std::shared_ptr<const expconfig::detector::TAPSVeto>;

    TAPSVeto_Energy(const detector_ptr_t& tapsveto,
                    const std::shared_ptr<DataManager>& calmgr,
                    const Calibration::Converter::ptr_t& converter,
                    double defaultPedestal = 100,
                    double defaultGain_BaF2 = 0.010,
                    double defaultGain_PbWO4 = 0.050,
                    double defaultThreshold_MeV = 0.1,
                    double defaultRelativeGain = 1.0);

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;

protected:
    const detector_ptr_t tapsveto_detector;
};

}} // namespace ant::calibration
