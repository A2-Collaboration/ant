#include "TAPSVeto_Energy.h"

#include "calibration/fitfunctions/FitLandau.h"

#include "expconfig/detectors/TAPSVeto.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPSVeto_Energy::TAPSVeto_Energy(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                                 std::shared_ptr<DataManager> calmgr,
                                 Calibration::Converter::ptr_t converter,
                                 double defaultPedestal,
                                 double defaultGain_BaF2,
                                 double defaultGain_PbWO4,
                                 double defaultThreshold,
                                 double defaultRelativeGain):
    Energy(Detector_t::Type_t::TAPSVeto,
           calmgr,
           converter,
           {defaultPedestal},
           {},   /* gains are set below */
           {defaultThreshold},
           {defaultRelativeGain}),
    tapsveto_detector(tapsveto)
{
    Gains.DefaultValues.clear();
    Gains.DefaultValues.resize(tapsveto->GetNChannels(), defaultGain_BaF2);
    for(unsigned ch=0; ch<tapsveto->GetNChannels(); ch++)
    {
        if(tapsveto->IsPbWO4(ch))
        {
            Gains.DefaultValues[ch] = defaultGain_PbWO4;
        }
    }
}

void TAPSVeto_Energy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          tapsveto_detector,
                          make_shared<gui::FitLandau>()
                          ));
    guis.emplace_back(std_ext::make_unique<GUI_Banana>(
                          GetName(),
                          RelativeGains,
                          calibrationManager,
                          tapsveto_detector
                          ));
}
