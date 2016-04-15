#include "Trigger.h"

#include "tree/TID.h"
#include "tree/TEventData.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

const std::string Trigger::ScalerName::Exptrigger_1MHz = "Exptrigger_1Mhz";
const std::string Trigger::ScalerName::Beampolmon_1MHz = "Beampolmon_1Mhz";

const std::string Trigger::ScalerName::TotalLivetime = "TotalLivetime";
const std::string Trigger::ScalerName::FaradayCup    = "FaradayCup";
const std::string Trigger::ScalerName::IonChamber    = "IonChamber";
const std::string Trigger::ScalerName::PbGlass       = "PbGlass";

// the CATCH TDC reference channels have never changed
// so map them here as const static
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_CATCH_TaggerCrate = {1000, 1400};
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_CATCH_CBCrate = {1001, 2000};
const Trigger::ReferenceTimingHitMapping_t Trigger_2014::Reference_V1190_TAPSPbWO4 = {1002, 29192};

void Trigger::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                            std::vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
    hit_mappings.emplace_back(
                Reference_CATCH_TaggerCrate.LogicalChannel,
                Reference_CATCH_TaggerCrate.AcquRawChannel
                );
    hit_mappings.emplace_back(
                Reference_CATCH_CBCrate.LogicalChannel,
                Reference_CATCH_CBCrate.AcquRawChannel
                );
}

void Trigger::ApplyTo(TEventData& reconstructed)
{

    /// @todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012
    /// so it might be best to implement such algorithms with some nicely designed interface into the
    /// pseudo-detector Trigger in expconfig/detectors

    double Esum = 0.0;
    for(const TDetectorReadHit& dethit : reconstructed.DetectorReadHits) {
        if(dethit.DetectorType != Detector_t::Type_t::CB)
            continue;
        if(dethit.ChannelType != Channel_t::Type_t::Integral)
            continue;
        for(double energy : dethit.Values)
            Esum += energy;
    }

    double TimeEsum = 0.0;
    double TimeE = 0.0;
    for(const TCluster& cluster : reconstructed.Clusters) {
        if(cluster.DetectorType != Detector_t::Type_t::CB)
            continue;
        if(!isfinite(cluster.Energy) || !isfinite(cluster.Time))
            continue;
        TimeEsum += cluster.Energy;
        TimeE += cluster.Energy*cluster.Time;
    }


    auto& triggerinfos = reconstructed.Trigger;
    triggerinfos.CBEnergySum = Esum;
    triggerinfos.CBTiming = TimeE/TimeEsum;

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
    /// \todo check when those V1190 modules actually were installed and
    /// move the emplace_back to base class maybe...
    hit_mappings.emplace_back(
                Reference_V1190_TAPSPbWO4.LogicalChannel,
                Reference_V1190_TAPSPbWO4.AcquRawChannel
                );

    // add the scaler mappings from internal map
    for(const auto& m : scaler_mapping) {
        scaler_mappings.emplace_back(m.first, m.second);
    }
}

string Trigger_2014::GetScalerReference(const string& scalername) const
{
    auto& m = scaler_mapping; // shortcut

    auto item = m.find(scalername);
    if(item == m.end())
        return "";

    if(item->second <= m.at(ScalerName::Exptrigger_1MHz))
        return ScalerName::Exptrigger_1MHz;
    if(item->second <= m.at(ScalerName::Beampolmon_1MHz))
        return ScalerName::Beampolmon_1MHz;

    return "";
}

Trigger_2014::scaler_mapping_t Trigger_2014::MakeScalerMapping()
{
    // build the mapping from human-readable names
    // to Acqu scaler indices
    scaler_mapping_t m;
    m[ScalerName::Exptrigger_1MHz] = 191;
    m[ScalerName::Beampolmon_1MHz] = 315;
    m[ScalerName::TotalLivetime]   = 190;
    m[ScalerName::FaradayCup]      = 313;
    m[ScalerName::IonChamber]      = 312;
    m[ScalerName::PbGlass]         = 311;
    return m;
}




