#include "TAPSVeto_Energy.h"

#include "calibration/fitfunctions/FitLandau.h"

#include "expconfig/detectors/TAPSVeto.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

vector<double> makeGains(const std::shared_ptr<const expconfig::detector::TAPSVeto>& tapsveto,
                         double defaultGain_BaF2,
                         double defaultGain_PbWO4)
{
    vector<double> gains(tapsveto->GetNChannels(), defaultGain_BaF2);
    for(unsigned ch=0; ch<tapsveto->GetNChannels(); ch++)
    {
        if(tapsveto->IsPbWO4(ch))
        {
            gains[ch] = defaultGain_PbWO4;
        }
    }
    return gains;
}

TAPSVeto_Energy::TAPSVeto_Energy(
        const detector_ptr_t& tapsveto,
        const std::shared_ptr<DataManager>& calmgr,
        const Calibration::Converter::ptr_t& converter,
        double defaultPedestal,
        double defaultGain_BaF2,
        double defaultGain_PbWO4,
        double defaultThreshold_Raw,
        double defaultThreshold_MeV,
        double defaultRelativeGain) :
    Energy(tapsveto,
           calmgr,
           converter,
           {defaultPedestal},
           makeGains(tapsveto, defaultGain_BaF2, defaultGain_PbWO4),
           {defaultThreshold_Raw},
           {defaultThreshold_MeV},
           {defaultRelativeGain}),
    tapsveto_detector(tapsveto)
{

}

void TAPSVeto_Energy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr options) {
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          options,
                          Pedestals,
                          calibrationManager,
                          tapsveto_detector,
                          make_shared<gui::FitLandau>()
                          ));
    guis.emplace_back(std_ext::make_unique<GUI_Banana>(
                          GetName(),
                          options,
                          RelativeGains,
                          calibrationManager,
                          tapsveto_detector,
                          interval<double>(600.0,700.0),
                          1.27 // MeV, from 2pi0 MC cocktail, -> same as PID, banana is quite flat there
                          ));
}
