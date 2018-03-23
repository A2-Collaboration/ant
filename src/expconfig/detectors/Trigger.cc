#include "Trigger.h"

#include "tree/TID.h"
#include "tree/TEventData.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

const std::string Trigger::ScalerName::TotalLivetime   = "TotalLivetime";
const std::string Trigger::ScalerName::Exptrigger_1MHz = "Exptrigger_1Mhz";
const std::string Trigger::ScalerName::ExpTrigger      = "ExpTrigger";
const std::string Trigger::ScalerName::L1Trigger       = "L1Trigger";

const std::string Trigger::ScalerName::EPTReferenceOR  = "EPTReferenceOR";

const std::string Trigger::ScalerName::PairSpecGate      = "PairSpecGate";
const std::string Trigger::ScalerName::TaggerReferenceOR = "TaggerReferenceOR";
const std::string Trigger::ScalerName::PbGlass           = "PbGlass";
const std::string Trigger::ScalerName::Paddle            = "Paddle";
const std::string Trigger::ScalerName::IonChamber        = "IonChamber";
const std::string Trigger::ScalerName::FaradayCup        = "FaradayCup";
const std::string Trigger::ScalerName::Beampolmon_1MHz   = "Beampolmon_1Mhz";

// the CATCH TDC reference channels have never changed
// so map them here as const static
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_CATCH_TaggerCrate = {1000, 1400};
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_CATCH_CBCrate = {1001, 2000};
const Trigger::ReferenceTimingHitMapping_t Trigger_2014::Reference_V1190_TAPSPbWO4 = {1002, 29192};
// these are for the new tagger, starting second half of 2017
/// \todo Possibly create new trigger struct for the new tagger
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_V1190_TaggerTDC1 = {1010, 927};
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_V1190_TaggerTDC2 = {1011, 1055};
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_V1190_TaggerTDC3_1 = {1012, 1151};
const Trigger::ReferenceTimingHitMapping_t Trigger::Reference_V1190_TaggerTDC3_2 = {1013, 1183};

void Trigger::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                            std::vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
    const auto& references = {
            Reference_CATCH_CBCrate,
            Reference_CATCH_TaggerCrate,
            Reference_V1190_TaggerTDC1,
            Reference_V1190_TaggerTDC2,
            Reference_V1190_TaggerTDC3_1,
            Reference_V1190_TaggerTDC3_2
        };

        for(const auto& reference : references) {
            hit_mappings.emplace_back(reference.LogicalChannel, reference.AcquRawChannel);
        }
}

void Trigger_2014::BuildMappings(std::vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                                 std::vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    // call base method for mappings which have never changed over the years
    Trigger::BuildMappings(hit_mappings, scaler_mappings);

    // the patterns are read here from the VUPROM, fBaseIndex is 0 in that ReadIRQ method
    // https://github.com/A2-Collaboration/acqu/blob/daq-master/acqu_core/AcquDAQ/src/TVME_VUPROM.cc
    for(unsigned j=0;j<patterns.size();j++) {
        hit_mappings.emplace_back(
                    LogicalChannel_t{Detector_t::Type_t::Trigger, Channel_t::Type_t::BitPattern, j},
                    j
                    );
    }

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

void ant::expconfig::detector::Trigger_2014::ApplyTo(const ReconstructHook::Base::readhits_t& hits)
{
    std::fill(patterns.begin(), patterns.end(), 0); // clear all patterns
    auto& dethits = hits.get_item(Type);
    for(TDetectorReadHit& readhit : dethits) {
        if(readhit.ChannelType != Channel_t::Type_t::BitPattern)
            continue;
        assert(readhit.Channel < patterns.size());
        assert(readhit.RawData.size() == 2); // expect 2 bytes
        patterns[readhit.Channel] = *reinterpret_cast<const uint16_t*>(addressof(readhit.RawData[0]));
    }
}

std::bitset<16> Trigger_2014::GetL1Pattern() const
{
    return patterns.at(0);
}

std::bitset<16> Trigger_2014::GetL2Pattern() const
{
    return patterns.at(1);
}

std::bitset<64> Trigger_2014::GetMultiplicityPattern() const
{
    /// \todo check if this is the right byte-order for the multiplicty pattern
    return std::bitset<64>(
            patterns.at(5).to_string() +
            patterns.at(4).to_string() +
            patterns.at(3).to_string() +
            patterns.at(2).to_string());
}

uint16_t Trigger_2014::GetMultiplicityValue() const
{
    return patterns.at(6).to_ulong();
}

std::bitset<16> Trigger_2014::GetHelicityPattern() const
{
    return patterns.at(7);
}

std::bitset<16> Trigger_2014::GetTriggerFiredPattern() const
{
    return patterns.at(8);
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
    m[ScalerName::TotalLivetime]   = 190;   // live counter vme-exptrigger
    m[ScalerName::Exptrigger_1MHz] = 191;   // free clock
    m[ScalerName::ExpTrigger]      = 192;
    m[ScalerName::L1Trigger]       = 194;

    m[ScalerName::EPTReferenceOR]  = 307;

    m[ScalerName::PairSpecGate]      = 308;
    m[ScalerName::TaggerReferenceOR] = 309;
    m[ScalerName::PbGlass]           = 310;
    m[ScalerName::Paddle]            = 311;
    m[ScalerName::IonChamber]        = 312;
    m[ScalerName::FaradayCup]        = 313;
                                    // 314 is unused
    m[ScalerName::Beampolmon_1MHz]   = 315; // free clock vme-beampolmon
    return m;
}





string Trigger_2007::GetScalerReference(const string&) const
{
    return "";
}


