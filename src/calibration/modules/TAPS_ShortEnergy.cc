#include "TAPS_ShortEnergy.h"
#include "TF1.h"
#include "expconfig/detectors/TAPS.h"
#include "calibration/fitfunctions/FitLandau.h"
#include "calibration/fitfunctions/FitGaus.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"
#include "calibration/gui/CalCanvas.h"
#include "tree/TDataRecord.h"
#include "base/Logger.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"

#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPS_ShortEnergy::TAPS_ShortEnergy(std::shared_ptr<expconfig::detector::TAPS> taps,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        double defaultPedestal,
        double defaultGain,
        double defaultThreshold,
        double defaultRelativeGain) :
    Energy(Detector_t::Type_t::TAPS,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain,
           Channel_t::Type_t::IntegralShort),
    taps_detector(taps)
{

}

void TAPS_ShortEnergy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          taps_detector,
                          make_shared<gui::FitLandau>()
                          ));
}


TAPS_ShortEnergy::GUI_Pedestals::GUI_Pedestals(const string& basename,
                                               Energy::CalibType& type,
                                               const std::shared_ptr<DataManager>& calmgr,
                                               const std::shared_ptr<expconfig::detector::TAPS>& taps,
                                               std::shared_ptr<gui::PeakingFitFunction> fitfunction) :
    Energy::GUI_Pedestals(basename, type, calmgr, taps, fitfunction),
    taps_detector(taps)
{
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ShortEnergy::GUI_Pedestals::DoFit(TH1* hist, unsigned channel)
{
    if(taps_detector->IsPbWO4(channel))
        return DoFitReturn_t::Skip;
    return Energy::GUI_Pedestals::DoFit(hist, channel);
}
