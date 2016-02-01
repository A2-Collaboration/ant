#include "Trigger.h"

#include "tree/TID.h"

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


bool Trigger_2014::Matches(const TID& tid) const {
    // timepoint roughly to Eta Prime beamtime...
    return std_ext::time_after(tid.Timestamp, "2014-07-01");
}

void Trigger_2014::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                                 std::vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // call base method for mappings which have never changed over the years
    Trigger::BuildMappings(hit_mappings, scaler_mappings);

    // The PbWO4 timings also have some reference,
    // but they are quite new...
    /// \todo check when those V1190 modules actually were installed
    hit_mappings.emplace_back(
                Reference_V1190_TAPSPbWO4,
                29192
                );

    scaler_mappings.emplace_back(
                "Exptrigger_1MHz",
                Type,
                Scaler_Exptrigger_1MHz.Channel,
                191
                );
    scaler_mappings.emplace_back(
                "Beampolmon_1MHz",
                Type,
                Scaler_Beampolmon_1MHz.Channel,
                315
                );

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


void ant::expconfig::detector::Trigger::ApplyTo(TEvent::Data& reconstructed)
{

    /// @todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012
    /// so it might be best to implement such algorithms with some nicely designed interface into the
    /// pseudo-detector Trigger in expconfig/detectors

    double Esum = 0.0;
    double TimeEsum = 0.0;
    double TimeE = 0.0;
    for(const TClusterPtr& cluster : reconstructed.Clusters) {
        if(cluster->DetectorType == Detector_t::Type_t::CB) {
            if(isfinite(cluster->Energy)) {
                Esum += cluster->Energy;
                if(isfinite(cluster->Time)) {
                    TimeEsum += cluster->Energy;
                    TimeE += cluster->Energy*cluster->Time;
                }
            }
        }
    }

    auto& triggerinfos = reconstructed.Trigger;
    triggerinfos.CBEnergySum = Esum;
    triggerinfos.CBTiming = TimeE/TimeEsum;

}