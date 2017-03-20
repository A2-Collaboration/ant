#include "CB.h"
#include <cassert>
#include "base/std_ext/container.h"

#include "detail/CB_elements.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;


CB::CB() :
    ClusterDetector_t(Detector_t::Type_t::CB),
    elements(elements_init)
{
    std::vector<unsigned> holes;
    std_ext::insertRange(holes,  26,  26);
    std_ext::insertRange(holes,  29,  38);
    std_ext::insertRange(holes,  40,  40);
    std_ext::insertRange(holes, 311, 311);
    std_ext::insertRange(holes, 315, 316);
    std_ext::insertRange(holes, 318, 319);
    std_ext::insertRange(holes, 353, 366);
    std_ext::insertRange(holes, 400, 402);
    std_ext::insertRange(holes, 405, 405);
    std_ext::insertRange(holes, 408, 408);
    std_ext::insertRange(holes, 679, 679);
    std_ext::insertRange(holes, 681, 689);
    std_ext::insertRange(holes, 691, 692);
    for(auto hole : holes) {
        SetElementFlags(hole, ElementFlag_t::Missing);
    }

    // the gap is reflected in the Geant simulation,
    // but not in the element positions from Acqu
    // we use 1.0 cm as some rough estimate here
    constexpr double gap = 1.0; // in cm

    for(auto& element : elements) {
        auto& pos = element.Position;
        if(pos.y<0)
            pos.y -= gap/2;
        else
            pos.y += gap/2;
    }
}

void CB::SetElementFlags(unsigned channel, const ElementFlags_t& flags) {
    elements.at(channel).Flags |= flags;
    if(IsIgnored(channel))
        SetTouchesHoleOfNeighbours(channel);
}

void CB::SetTouchesHoleOfNeighbours(unsigned hole)
{
    for(auto neighbour : elements.at(hole).Neighbours)
        elements.at(neighbour).TouchesHole = true;
}


bool CB::IsHole(unsigned channel) const
{
    return elements[channel].Flags.test(ElementFlag_t::Missing);
}

void CB::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
                       vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
    // CB has only hit_mappings to add,
    // no scalers
    unsigned true_elements = 0;
    for(const Element_t& element : elements)  {
        // exclude holes from mapping
        if(IsHole(element.Channel))
            continue;

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC);

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC);

        true_elements++;
    }

    assert(true_elements == 672);
}




