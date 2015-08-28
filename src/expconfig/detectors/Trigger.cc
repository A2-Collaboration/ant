#include "Trigger.h"

#include "tree/THeaderInfo.h"

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

    std::list<scaler_mapping_t> reference_scalers;

    reference_scalers.emplace_back(
                "Exptrigger_1MHz",
                Type,
                Scaler_Exptrigger_1MHz.Channel,
                191
                );
    reference_scalers.emplace_back(
                "Beampolmon_1MHz",
                Type,
                Scaler_Beampolmon_1MHz.Channel,
                315
                );

    // the reference scalers are added as TSlowControl
    // and as logical channels
    // (might be used by scaler calibrations from other detectors)
    for(const scaler_mapping_t& scaler : reference_scalers) {
        scaler_mappings.push_back(scaler);
        // add again as DetectorRead item
        scaler_mappings.emplace_back(scaler.Entries);
    }

    // some interesting scalers, added as TSlowControl only
    scaler_mappings.emplace_back(
                "TotalLivetime",
                Type,
                100,
                190
                );
    scaler_mappings.emplace_back(
                "FaradayCup",
                Type,
                200,
                313
                );
    scaler_mappings.emplace_back(
                "IonChamber",
                Type,
                201,
                312
                );
    scaler_mappings.emplace_back(
                "PbGlass",
                Type,
                202,
                311
                );
}
