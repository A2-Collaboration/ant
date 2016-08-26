#include "TaggEff.h"

#include "DataManager.h"
#include "tree/TCalibrationData.h"

#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;


TaggEff::TaggEff(
        const shared_ptr<ant::TaggerDetector_t>& tagger,
        const shared_ptr<DataManager>& calmgr) :
    Module(GetDataName()),
    Tagger(tagger),
    CalibrationManager(calmgr)
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
            TCalibrationData cdataErrrors;
            if(!CalibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;
            if(cdata.Data.size() != 1)
                return;
            // if errors are defined:
            const TKeyValue<double>& kv = cdata.Data.front();
            if(CalibrationManager->GetData(GetDataErrorsName(), currPoint, cdataErrrors,nextChangePoint))
            {
                const TKeyValue<double>& kverr = cdataErrrors.Data.front();
                Tagger->SetTaggEff(kv.Key,{kv.Value,kverr.Value});
                return;
            }
            // if errors not defined:
            Tagger->SetTaggEff(kv.Key,kv.Value); // error on tagg eff still missing
        }
    };
}


