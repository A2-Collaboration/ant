#include "TAPSVeto_Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include "calibration/fitfunctions/FitLandau.h"

#include "expconfig/detectors/TAPSVeto.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPSVeto_Energy::TAPSVeto_Energy(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                                 std::shared_ptr<DataManager> calmgr,
                                 Calibration::Converter::ptr_t converter,
                                 double defaultPedestal,
                                 double defaultGain,
                                 double defaultThreshold,
                                 double defaultRelativeGain):
    Energy(Detector_t::Type_t::TAPSVeto,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain),
    tapsveto_detector(tapsveto)
{

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
