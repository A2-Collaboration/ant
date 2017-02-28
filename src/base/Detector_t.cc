#include "Detector_t.h"

#include "base/interval.h"
#include "base/std_ext/string.h"

#include <sstream>
#include <map>

using namespace ant;
using namespace std;

const Detector_t::Any_t Detector_t::Any_t::None(0);
const Detector_t::Any_t Detector_t::Any_t::Tracker(Type_t::MWPC0 | Type_t::MWPC1);
const Detector_t::Any_t Detector_t::Any_t::CB_Apparatus(Detector_t::Any_t::Tracker | Type_t::PID | Type_t::CB );
const Detector_t::Any_t Detector_t::Any_t::TAPS_Apparatus(Type_t::TAPS | Type_t::TAPSVeto);
const Detector_t::Any_t Detector_t::Any_t::Calo(Type_t::CB | Type_t::TAPS);
const Detector_t::Any_t Detector_t::Any_t::Veto(Type_t::PID | Type_t::TAPSVeto);


ostream& ant::Detector_t::Any_t::Print(ostream& stream) const  {
    if(bitfield == 0)
        return stream << "None";

    using i_t = typename underlying_type<Type_t>::type;
    using bitfield_t = decltype(bitfield);

    bitfield_t temp = bitfield;

    i_t i_max = 0;
    while(temp >>= 1)
        i_max++;

    i_t i = 0;
    temp = 1;
    while(temp<=bitfield) {
        if(bitfield & temp) {
            stream << Detector_t::ToString(
                          static_cast<Detector_t::Type_t>(i)
                          );
            if(i<i_max)
                stream << "|";
        }
        ++i;
        temp <<= 1;
    }
    return stream;
}

#define MAKE_DETECTOR_TYPE_ENTRY(det) {Detector_t::Type_t::det, #det}

const map<Detector_t::Type_t, string> detectorTypeMap = {
    MAKE_DETECTOR_TYPE_ENTRY(CB),
    MAKE_DETECTOR_TYPE_ENTRY(Cherenkov),
    MAKE_DETECTOR_TYPE_ENTRY(MWPC0),
    MAKE_DETECTOR_TYPE_ENTRY(MWPC1),
    MAKE_DETECTOR_TYPE_ENTRY(PID),
    MAKE_DETECTOR_TYPE_ENTRY(Tagger),
    MAKE_DETECTOR_TYPE_ENTRY(TaggerMicro),
    MAKE_DETECTOR_TYPE_ENTRY(EPT),
    MAKE_DETECTOR_TYPE_ENTRY(Moeller),
    MAKE_DETECTOR_TYPE_ENTRY(PairSpec),
    MAKE_DETECTOR_TYPE_ENTRY(TAPS),
    MAKE_DETECTOR_TYPE_ENTRY(TAPSVeto),
    MAKE_DETECTOR_TYPE_ENTRY(Trigger),
    MAKE_DETECTOR_TYPE_ENTRY(Raw),
};

const char* ant::Detector_t::ToString(const Type_t& type)
{
    auto it = detectorTypeMap.find(type);
    if(it == detectorTypeMap.end())
        throw runtime_error("Unknown detector type");
    return it->second.c_str();
}

Detector_t::Type_t Detector_t::FromString(const string& str)
{
    using map_t = decltype(detectorTypeMap);
    auto it = find_if(detectorTypeMap.begin(), detectorTypeMap.end(), [str] (const map_t::value_type& v ) {
        return str == v.second;
    });
    if(it == detectorTypeMap.end())
        throw runtime_error("Unknown detector string "+str);
    return it->first;
}

void Detector_t::SetElementFlag(const std::vector<unsigned>& channels, ElementFlag_t flag)
{
    for(auto& ch : channels)
        SetElementFlag(ch, flag);
}


ant::Detector_t::Any_t::operator string() const {
    stringstream s;
    Print(s);
    return s.str();
}


const char* ant::Channel_t::ToString(const Type_t& type)
{
    switch(type) {
    case Channel_t::Type_t::BitPattern:
        return "BitPattern";
    case Channel_t::Type_t::Integral:
        return "Integral";
    case Channel_t::Type_t::IntegralAlternate:
        return "IntegralAlternate";
    case Channel_t::Type_t::IntegralShort:
        return "IntegralShort";
    case Channel_t::Type_t::IntegralShortAlternate:
        return "IntegralShortAlternate";
    case Channel_t::Type_t::Timing:
        return "Timing";
    case Channel_t::Type_t::Raw:
        return "Raw";
    }
    throw runtime_error("Not implemented");
}


double TaggerDetector_t::GetPhotonEnergyWidth(unsigned channel) const
{
    if(channel >= GetNChannels())
        throw out_of_range(std_ext::formatter() << "Tagger channel index out of range: " << channel << " (" << GetNChannels() << ")");

    if(channel == 0) {
        return abs(GetPhotonEnergy(channel+1) - GetPhotonEnergy(channel));
    }
    else if(channel == GetNChannels() -1) {
        return abs(GetPhotonEnergy(channel) - GetPhotonEnergy(channel-1));
    }
    else {
        // =  (E(ch+1) - central)/2 + (central - E(ch-1))/2
        return abs(GetPhotonEnergy(channel+1) - GetPhotonEnergy(channel-1)) / 2.0;
    }
}

bool TaggerDetector_t::TryGetChannelFromPhoton(double photonEnergy, unsigned& channel) const
{
    for(channel=0;channel<GetNChannels();++channel) {
        const auto& i = interval<double>::CenterWidth(
                            GetPhotonEnergy(channel),
                            GetPhotonEnergyWidth(channel)
                            );
        if(i.Contains(photonEnergy))
            return true;
    }
    return false;
}
