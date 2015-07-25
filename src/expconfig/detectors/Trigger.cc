#include "Trigger.h"

#include "tree/THeaderInfo.h"
#include "base/std_ext.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;


void Trigger::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                            std::vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
    // the CATCH TDC reference channels have never changed
    // so map them here in the base class (to be used by derived classes)
    hit_mappings.emplace_back(
                Reference_CATCH_TaggerCrate,
                1400
                );
    hit_mappings.emplace_back(
                Reference_CATCH_CBCrate,
                2000
                );
}


bool Trigger_2014::Matches(const THeaderInfo& headerInfo) const {
    // timepoint roughly to Eta Prime beamtime...
    return std_ext::time_after(headerInfo.Timestamp, "2014-07-01");
}

void Trigger_2014::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                                 std::vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // call base method for mappings which have never changed over the years
    Trigger::BuildMappings(hit_mappings, scaler_mappings);

    std::list<scaler_mapping_t> scalers;

    scalers.emplace_back(
                "Exptrigger_1MHz",
                Type,
                Scaler_Exptrigger_1MHz,
                191
                );

    scalers.emplace_back(
                "Beampolmon_1MHz",
                Type,
                Scaler_Beampolmon_1MHz,
                315
                );

    // the four most important scalers
    scalers.emplace_back(
                "TotalLivetime",
                Type,
                100,
                190
                );
    scalers.emplace_back(
                "FaradayCup",
                Type,
                200,
                313
                );
    scalers.emplace_back(
                "IonChamber",
                Type,
                201,
                312
                );
    scalers.emplace_back(
                "PbGlass",
                Type,
                202,
                311
                );

    // those most important scalers are all added as TSlowControl
    // and as logical channels
    for(const scaler_mapping_t& scaler : scalers) {
        scaler_mappings.push_back(scaler);
        // add again as DetectorRead item
        scaler_mappings.emplace_back(
                    scaler.LogicalChannel.DetectorType,
                    scaler.LogicalChannel.Channel,
                    scaler.RawChannels[0].RawChannel // only one channel assumed
                    );
    }

}
