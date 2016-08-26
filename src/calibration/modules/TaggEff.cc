#include "TaggEff.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitGaus.h"
#include "tree/TCalibrationData.h"

#include "expconfig/detectors/PID.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TGraph.h"
#include "TFitResult.h"

#include <limits>
#include <cmath>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

TaggEff::TaggEff(
        const std::shared_ptr<expconfig::detector::TaggerDetector_t>& pid,
        const std::shared_ptr<DataManager>& calmgr) :
    Module("TaggEff"),
    pid_detector(pid),
    calibrationManager(calmgr)
{
}

TaggEff::~TaggEff()
{
}

void TaggEff::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& ) {
}

std::list<Updateable_traits::Loader_t> TaggEff::GetLoaders()
{
    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(!calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;
            if(cdata.Data.size() != 1)
                return;
            const TKeyValue<double>& kv = cdata.Data.front();
//            pid_detector->SetPhiOffset(kv.Value);
        }
    };
}


