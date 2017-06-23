#include "TaggEff.h"

#include "DataManager.h"
#include "tree/TCalibrationData.h"
#include "tree/TEventData.h"

#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"


using namespace std;
using namespace ant;
using namespace ant::calibration;


TaggEff::TaggEff(
        const shared_ptr<ant::TaggerDetector_t>& tagger,
        const shared_ptr<DataManager>& calmgr) :
    BaseModule(GetModuleName(tagger->Type)),
    Tagger(tagger),
    CalibrationManager(calmgr)
{
}

TaggEff::~TaggEff()
{
}

string TaggEff::GetModuleNameSuffix()
{
    return "TaggEff";
}

string TaggEff::GetModuleName(Detector_t::Type_t type) {
    return std_ext::formatter()
            << Detector_t::ToString(type)
            << "_" << GetModuleNameSuffix();
}



std::list<Updateable_traits::Loader_t> TaggEff::GetLoaders()
{
    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(!CalibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;

            currentTaggEff.resize(0);
            for ( const auto& data: cdata.Data )
            {
                const auto channel = data.Key;
                currentTaggEff.resize(channel+1);
                currentTaggEff.at(channel).Value = data.Value;
            }

            for ( const auto& data : cdata.FitParameters)
            {
                // expect the currentTaggEff to be resized from previous filling
                // might throw index-out-of-bound exception if this assumption is not true
                currentTaggEff.at(data.Key).Value = data.Value.front();
            }

            // flag that we have just loaded something
            loadedTaggEff = currPoint;
        }
    };
}


void ant::calibration::TaggEff::ApplyTo(TEventData& reconstructed)
{
    if(!loadedTaggEff.IsInvalid()) {
        // serialize tagging efficiencies into stringstream buffer
        stringstream ss;
        cereal::BinaryOutputArchive ar(ss);
        ar(currentTaggEff);

        // insert the tagging efficiencies, given from TCalibrationData,
        // as cereal'ized TSlowControl item
        reconstructed.SlowControls.emplace_back(
                    TSlowControl::Type_t::Ant,
                    TSlowControl::Validity_t::Forward,
                    loadedTaggEff.Timestamp,
                    GetName(),
                    "cereal'ized TaggEff"
                    );

        // create the Payload_String inside TSlowControl
        // with key=0 inside Payload_String
        reconstructed.SlowControls.back()
                .Payload_String.emplace_back(0,ss.str());

        // clear flag by invalidating
        loadedTaggEff = TID();
    }
}