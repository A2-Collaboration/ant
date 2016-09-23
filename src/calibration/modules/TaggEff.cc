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
    BaseModule(GetDataName()),
    Tagger(tagger),
    CalibrationManager(calmgr)
{
}

TaggEff::~TaggEff()
{
}



std::list<Updateable_traits::Loader_t> TaggEff::GetLoaders()
{
    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(!CalibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;
            for ( const auto& data: cdata.Data )
            {
                auto channel = data.Key;
                auto taggEff = data.Value;
                auto taggEffError = std_ext::NaN;

                for ( const auto& fitP: cdata.FitParameters)
                {
                    if (channel == fitP.Key)
                    {
                        taggEffError = fitP.Value.front();
                        break;
                    }
                }

                Tagger->SetTaggEff(channel,{taggEff,taggEffError});

            }
        }
    };
}


