#include "CB.h"
#include <cassert>
#include "base/std_ext/vector.h"

#include "detail/CB_elements.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

CB::CB() : ClusterDetector_t(Detector_t::Type_t::CB) {
    auto& holes = ignoredChannels;
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

}

void CB::SetIgnored(unsigned channel) {
    ignoredChannels.push_back(channel);
}

bool CB::IsIgnored(unsigned channel) const {
    return std_ext::contains(ignoredChannels, channel);
}

void CB::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
                       vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
    // CB has only hit_mappings to add,
    // no scalers
    unsigned true_elements = 0;
    for(const Element_t& element : elements)  {
        // exclude ignored elements from mapping
        if(IsIgnored(element.Channel))
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

    assert(true_elements <= 672);
}



