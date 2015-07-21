#include "EPT.h"
#include <cassert>
#include "base/std_ext.h"
#include "tree/THeaderInfo.h"

#include "detail/EPT_2014_elements.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;



bool EPT_2014::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp,
                                 "2014-07-01", "2014-12-31");
}

void EPT::BuildMappings(
        vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
        vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    for(const Element_t& element : elements) {
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC
                                  );
        // ADC information are rarely present for the EPT
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC
                                  );


    }
}



