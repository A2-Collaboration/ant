#pragma once

#include "calibration/Calibration.h"
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

    TAPSVeto_Energy(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                    std::shared_ptr<DataManager> calmgr,
                    Calibration::Converter::ptr_t converter,
                    double defaultPedestal = 100,
                    double defaultGain = 0.010,
                    double defaultThreshold = 0.1,
                    double defaultRelativeGain = 1.0);

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis) override;

protected:
    std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto_detector;
};

}} // namespace ant::calibration
