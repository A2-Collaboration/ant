#include "TAPS.h"

#include "detail/TAPS_2007_BaF2_elements.h"
#include "detail/TAPS_2009_03_BaF2_elements.h"
#include "detail/TAPS_2009_03_PbWO4_elements.h"
#include "detail/TAPS_2013_11_BaF2_elements.h"
#include "detail/TAPS_2013_11_PbWO4_elements.h"


#include "tree/TID.h"
#include "tree/TCandidate.h"

#include "base/interval.h"

#include <iostream>
#include <cassert>
#include <algorithm>

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;


void TAPS::SetElementFlags(unsigned channel, const ElementFlags_t& flags) {
    if(flags & ElementFlag_t::Missing)
        throw Exception("TAPS does not have any intentionally missing elements");

    clusterelements[channel]->Flags |= flags;
    // set touches hole if element is marked broken
    if(flags & ElementFlag_t::Broken) {
        for(auto neighbour : clusterelements[channel]->Neighbours) {
            clusterelements[neighbour]->TouchesHole = true;
        }
    }
}

double TAPS::GetTimeOfFlight(double clustertime, unsigned channel, double trigger_reftime) const {
    return clustertime - clusterelements.at(channel)->ToFOffset - trigger_reftime;
}

double TAPS::GetBeta(const TCandidate& cand_taps, double trigger_reftime) const {
    const auto taps_cluster = cand_taps.FindCaloCluster();
    if(!taps_cluster)
        return std_ext::NaN;
    const auto dt = GetTimeOfFlight(taps_cluster->Time, taps_cluster->CentralElement,
                                    trigger_reftime);
    const auto s = GetZPosition();
    constexpr auto c = 30.0; // velocity of light in cm/ns
    return s / (s + c * dt * cos(cand_taps.Theta));
}

double TAPS::GetZPosition() const
{
    double z = 145.7;
    if (PizzaInstalled && CherenkovInstalled)
        z = 198;
    else if (CherenkovInstalled)
        z = 174.2;
    else if (PizzaInstalled)
        z = 188;
    return z;
}

double TAPS::GetRadius() const
{
    /// \todo check value?
    return 70.0;
}

void TAPS::BuildMappings(
        vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
        vector<UnpackerAcquConfig::scaler_mapping_t>&) const
{
    for(const BaF2_Element_t& element : BaF2_elements)  {

        // TAC provides timing information

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TAC);

        // the flag UseSensitiveChannels swaps the mapping
        // of LGS/LG to type Integral/IntegralAlternate
        // only Integral is used in the clustering

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  UseSensitiveChannels ? element.LGS : element.LG
                                                         );

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::IntegralAlternate,
                                  element.Channel,
                                  UseSensitiveChannels ? element.LG : element.LGS
                                                         );

        // the same logic applies to the short gate integral

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::IntegralShort,
                                  element.Channel,
                                  UseSensitiveChannels ? element.SGS : element.SG
                                                         );

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::IntegralShortAlternate,
                                  element.Channel,
                                  UseSensitiveChannels ? element.SG : element.SGS
                                                         );
    }

    // the PbWO4 are a bit simpler

    for(const PbWO4_Element_t& element : PbWO4_elements)  {

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC);

        /// \todo maybe use different switch for BaF2 and PbWO4 sensitive?!
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  UseSensitiveChannels ? element.QDCL : element.QDCH
                                                         );

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::IntegralAlternate,
                                  element.Channel,
                                  UseSensitiveChannels ? element.QDCH : element.QDCL
                                                         );
    }

}

unsigned TAPS::GetRing(const unsigned channel) const
{
    const vector< interval<unsigned> > ringRanges = {
        {0, 0}, {1, 2}, {3, 5}, {6, 9}, {10, 14}, {15, 20},
        {21, 27}, {28, 35}, {36, 44}, {45, 54}, {55, 63}
    };
    assert(ringRanges.size() == 11);
    // the GetHexChannel is always between 0 and 383
    // the %64 maps it to the first sector
    constexpr unsigned HexElementsPerSector = NHexElements/NSectors;
    unsigned hexChannelFirstSector = GetHexChannel(channel) % HexElementsPerSector;

    for(size_t i=0;i<ringRanges.size();i++) {
        if(ringRanges[i].Contains(hexChannelFirstSector))
            return i+1; // ring0 is central element (never installed)
    }
    throw Exception("Cannot find ring for TAPS channel "+to_string(channel));
}

unsigned TAPS::GetHexChannel(const unsigned channel) const
{
    constexpr unsigned HexElementsPerSector = NHexElements/NSectors;

    // no PbWO4s at all, so nothing to do
    if(BaF2_elements.size() == NHexElements && PbWO4_elements.size() == 0)
        return channel;

    const unsigned elementsPerSector = (BaF2_elements.size()+PbWO4_elements.size())/NSectors;

    const unsigned channelSector = channel / elementsPerSector;
    const unsigned channelFirstSector = channel % elementsPerSector;

    const unsigned PbWO4_elementsPerSector = PbWO4_elements.size() / NSectors;
    assert(PbWO4_elementsPerSector % PbWO4PerHex == 0);

    // Handle PbWO4 channel first
    if(channelFirstSector < PbWO4_elementsPerSector) {
        // there are 4 PbWO4 elements per hex
        const unsigned PbWO4_hexIndex = channelFirstSector / PbWO4PerHex;
        return PbWO4_hexIndex + channelSector * HexElementsPerSector;
    }
    else {
        const unsigned PbWO4_offset = PbWO4_elementsPerSector / PbWO4PerHex;
        const unsigned BaF2_hexIndex = (channelFirstSector - PbWO4_elementsPerSector) + PbWO4_offset;
        return BaF2_hexIndex + channelSector * HexElementsPerSector;
    }
}

bool TAPS::IsPbWO4(const unsigned channel) const
{
    if(PbWO4_elements.empty())
        return false;

    const unsigned elementsPerSector = (BaF2_elements.size()+PbWO4_elements.size())/NSectors;

    const unsigned channelFirstSector = channel % elementsPerSector;

    const unsigned PbWO4_elementsPerSector = PbWO4_elements.size() / NSectors;
    assert(PbWO4_elementsPerSector % PbWO4PerHex == 0);

    return channelFirstSector < PbWO4_elementsPerSector;
}

std::vector<unsigned> TAPS::GetBaF2Channels() const
{
    vector<unsigned> channels(BaF2_elements.size());

    transform(BaF2_elements.begin(), BaF2_elements.end(), channels.begin(),
              [](const BaF2_Element_t& element){ return element.Channel; });

    return channels;
}

std::vector<unsigned> TAPS::GetPbWO4Channels() const
{
    vector<unsigned> channels(PbWO4_elements.size());

    transform(PbWO4_elements.begin(), PbWO4_elements.end(), channels.begin(),
              [](const PbWO4_Element_t& element){ return element.Channel; });

    return channels;
}

void TAPS::InitClusterElements()
{
    assert(BaF2_elements.size()>0);
    assert(BaF2_elements.size() % NSectors == 0);
    assert(PbWO4_elements.size() % NSectors == 0);

    clusterelements.resize(BaF2_elements.size()+PbWO4_elements.size(), nullptr);



    // apply the z-position depending on Cherenkov
    // and build the clusterelements
    // we assume that the channel elements are consecutive
    const double zpos = GetZPosition();
    for(auto& baf2 : BaF2_elements) {
        baf2.Position.z = zpos;
        clusterelements[baf2.Channel] = addressof(baf2);
    }
    for(auto& pbwo4 : PbWO4_elements) {
        pbwo4.Position.z = zpos;
        clusterelements[pbwo4.Channel] = addressof(pbwo4);
    }

    assert(clusterelements.size()>0);

    // check that every element was actually set
    for(auto elem : clusterelements) {
        (void)elem; // prevent unused variable warning in Release build
        assert(elem != nullptr);
    }

    // set touches hole
    /// \todo think of something better than counting neighbours...

    for(auto& baf2 : BaF2_elements) {
        // non-edge baf2 elements have at least 6 neighbours
        if(baf2.Neighbours.size()<6)
            baf2.TouchesHole = true;
    }

    for(auto& pbwo4 : PbWO4_elements) {
        bool hasBaF2 = false;
        std::for_each(pbwo4.Neighbours.begin(), pbwo4.Neighbours.end(),
                      [this,&hasBaF2] (unsigned n) { hasBaF2 |= IsBaF2(n); });
        if(hasBaF2)
            continue;
        if(pbwo4.Neighbours.size()<8)
            pbwo4.TouchesHole = true;
    }
}



