#include "PID.h"

#include "tree/THeaderInfo.h"


#include "detail/PID_2004_elements.h"
#include "detail/PID_2009_05_elements.h"
#include "detail/PID_2009_06_elements.h"
#include "detail/PID_2009_07_elements.h"
#include "detail/PID_2014_elements.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;



void PID::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                        vector<UnpackerAcquConfig::scaler_mapping_t>&) const
{
    for(const Element_t& element : elements)
    {
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC);

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC);
    }
}

bool PID_2004::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp, "2004-04-30", "2009-05-09");
}


bool PID_2009_05::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp, "2009-05-10", "2009-06-29");
}


bool PID_2009_06::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp, "2009-06-30", "2009-07-12");
}


bool PID_2009_07::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp, "2009-07-13", "2014-01-25");
}


bool PID_2014::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_after(headerInfo.Timestamp, "2014-01-26");
}
